/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "scalarGeneralExchange.H"
#include "addToRunTimeSelectionTable.H"
#include "dataExchangeModel.H"
#define ALARGECONCENTRATION 1e32

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(scalarGeneralExchange, 0);

addToRunTimeSelectionTable
(
    forceModel,
    scalarGeneralExchange,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
scalarGeneralExchange::scalarGeneralExchange
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    scalarTransportProperties_                  //this is clumsy, but effective
    (
        IOobject
        (
            "scalarTransportProperties",
            sm.mesh().time().constant(),
            sm.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    generalPropsDict_(scalarTransportProperties_.subDict("generalManualProps")),
    voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),             //common names/data
    velFieldName_(propsDict_.lookup("velFieldName")),
    tempFieldName_(propsDict_.lookup("tempFieldName")),                             //temperature names/data
    partTempName_(propsDict_.lookup("partTempName")),
    partHeatFluxName_(propsDict_.lookupOrDefault<word>(      "partHeatFluxName", "na")),
    partHeatTransCoeffName_(propsDict_.lookupOrDefault<word>("partHeatTransCoeffName", "na")),
    partHeatFluidName_(propsDict_.lookupOrDefault<word>(     "partHeatFluidName", "na")),
    partDat_(NULL),
    partDatFlux_(NULL),
    partDatTransCoeff_(NULL),
    partDatFluid_(NULL),
    partDatTmpExpl_(NULL),
    partDatTmpImpl_(NULL),
    validPartFlux_(false),
    validPartTransCoeff_(false),
    validPartFluid_(false),
    haveTemperatureEqn_(false),
    useLiMason_(false),
    lambda_(readScalar(propsDict_.lookup("lambda"))),
    Prandtl_(readScalar(propsDict_.lookup("Prandtl"))),
    eulerianFieldNames_( generalPropsDict_.lookup("eulerianFields")), 
    partSpeciesNames_(propsDict_.lookup("partSpeciesNames")),                       
    partSpeciesFluxNames_(propsDict_.lookup("partSpeciesFluxNames")),
    partSpeciesTransCoeffNames_(propsDict_.lookup("partSpeciesTransCoeffNames")),
    partSpeciesFluidNames_(propsDict_.lookup("partSpeciesFluidNames")),
    DMolecular_(propsDict_.lookup("DMolecular")),
    parameterVap_(propsDict_.lookup("parameterVap")),
    Rvap_(propsDict_.lookupOrDefault<scalar>("Rvap", 461.5)),
   // alphaImExSplit_(propsDict_.lookupOrDefault<scalar>("alphaImExSplit", 0.5)),
    deltaHEvap_("deltaHEvap", dimLength*dimLength/dimTime/dimTime, 1),
    maxSource_(1e30),
    scaleDia_(1.)
{
    allocateMyArrays(0.0);

    if (propsDict_.found("maxSource"))
    {
        maxSource_=readScalar(propsDict_.lookup ("maxSource"));
        Info << "limiting eulerian source field to: " << maxSource_ << endl;
    }

    if (propsDict_.found("useLiMason"))
    {
        useLiMason_=readBool(propsDict_.lookup ("useLiMason"));
        Info << "setting for useLiMason: " << useLiMason_ << endl;
    }
    if(useLiMason_)
        Nusselt=&scalarGeneralExchange::NusseltLiMason;
    else
        Nusselt=&scalarGeneralExchange::NusseltDeenEtAl;

    if(partHeatFluxName_!="na")        
    {   validPartFlux_=true;
        Info << "Found a valid partHeatFluxName: " << partHeatFluxName_ << endl;
    }
    if(partHeatTransCoeffName_!="na")
    {   validPartTransCoeff_=true;  
        Info << "Found a valid partHeatTransCoeffName: " << partHeatTransCoeffName_ << endl;
    }
    if(partHeatFluidName_!="na")
    {   validPartFluid_=true;
        Info << "Found a valid partHeatFluidName: " << partHeatFluidName_ << endl;
    }

    if( validPartTransCoeff_ && !validPartFluid_ )
        FatalError <<"Transfer coefficient set, but and fluid name missing. Check your entries in the couplingProperties! \n" 
                   << abort(FatalError);    

    if( !validPartTransCoeff_ && validPartFluid_ )
        FatalError <<"Fluid name set, but transfer coefficient missing. Check your entries in the couplingProperties! \n" 
                   << abort(FatalError);    
    
    if(!validPartFlux_ && !(validPartTransCoeff_ && validPartFluid_) )
        FatalError <<"You must set a valid heat flux name, or a valid transfer coefficient and fluid name \n" 
                   << abort(FatalError);

    // init force sub model
    setForceSubModels(propsDict_);

    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(3,true); // activate search for verbose switch
    forceSubM(0).setSwitchesList(4,true); // activate search for interpolate switch
    forceSubM(0).setSwitchesList(8,true); // activate scalarViscosity switch

    // read those switches defined above, if provided in dict
    for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
        forceSubM(iFSub).readSwitches();

    particleCloud_.checkCG(true);

    if (propsDict_.found("scale"))
        scaleDia_=scalar(readScalar(propsDict_.lookup("scale")));

    //check species names
    Info << "scalarGeneralExchange found the following eulerianFieldName: " << eulerianFieldNames_ << endl;
    int numTempEqn=0;
    for(int iEul=0;iEul<eulerianFieldNames_.size(); iEul++)
        if(eulerianFieldNames_[iEul]==tempFieldName_)
        {
            haveTemperatureEqn_ = true;
            numTempEqn = 1;
            Info << "scalarGeneralExchange identified the eulerianField '" << tempFieldName_ 
                 << "' as being the temperature field" << endl;
        }

    //check if enough particle properties have been provided
    if(partSpeciesNames_.size()!=(eulerianFieldNames_.size()-numTempEqn))   
        FatalError <<"Not enough partSpeciesNames specified in the couplingProperties file. \n" 
                   << abort(FatalError);
    else
        Info << "Found valid partSpeciesNames: " << partSpeciesNames_ << endl;
    if(partSpeciesFluxNames_.size()!=(eulerianFieldNames_.size()-numTempEqn))   
        FatalError <<"Not enough partSpeciesFluxNames specified in the couplingProperties file. \n" 
                   << abort(FatalError);
    else
        Info << "Found valid partSpeciesFluxNames: " << partSpeciesFluxNames_ << endl;
    if(partSpeciesTransCoeffNames_.size()!=(eulerianFieldNames_.size()-numTempEqn))   
        FatalError <<"Not enough partSpeciesTransCoeffNames specified in the couplingProperties file. \n" 
                   << abort(FatalError);
    else
        Info << "Found valid partSpeciesTransCoeffNames: " << partSpeciesTransCoeffNames_ << endl;
    if(partSpeciesFluidNames_.size()!=(eulerianFieldNames_.size()-numTempEqn))   
        FatalError <<"Not enough partSpeciesFluidNames specified in the couplingProperties file. \n" 
                   << abort(FatalError);
    else
        Info << "Found valid partSpeciesFluidNames: " << partSpeciesFluidNames_ << endl;
    if(DMolecular_.size()!=(eulerianFieldNames_.size()-numTempEqn))   
        FatalError <<"Not enough DMolecular specified in the couplingProperties file. \n" 
                   << abort(FatalError);

    for(int iPart=0;iPart<partSpeciesNames_.size(); iPart++)
    {
        if(partSpeciesNames_[iPart]=="none")
            particleSpeciesValue_.setSize(particleSpeciesValue_.size()+1, -1); //will not consider this species for coupling
        else if(partSpeciesNames_[iPart]=="zero")
            particleSpeciesValue_.setSize(particleSpeciesValue_.size()+1, 0.0); //will not consider this species for coupling
        else
        {
            particleSpeciesValue_.setSize(particleSpeciesValue_.size()+1, -2*ALARGECONCENTRATION); //set to a very large value to request pull
        }
    }

    if(nrForceSubModels()>1)
        FatalError <<"nrForceSubModels() must be zero or one for scalarGeneralExchange. \n" 
                   << abort(FatalError);

    if (probeIt_ && propsDict_.found("suppressProbe"))
        probeIt_=!Switch(propsDict_.lookup("suppressProbe"));
    if(probeIt_)
    {
        particleCloud_.probeM().initialize(typeName, "scalarGeneralExchange.logDat");
        particleCloud_.probeM().vectorFields_.append("Urel");               //first entry must the be the vector to probe
        particleCloud_.probeM().scalarFields_.append("Nu");                 //other are debug
        particleCloud_.probeM().scalarFields_.append("Rep");                //other are debug
        particleCloud_.probeM().writeHeader();
    }
}

// Construct from components for scalarGeneralExchangePhaseChange
scalarGeneralExchange::scalarGeneralExchange
(
    const dictionary& dict,
    cfdemCloud& sm,
    word        dictName
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(dictName + "Props")),
    scalarTransportProperties_                  //this is clumsy, but effective
    (
        IOobject
        (
            "scalarTransportProperties",
            sm.mesh().time().constant(),
            sm.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    generalPropsDict_(scalarTransportProperties_.subDict("generalManualProps")),
    voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),             //common names/data
    velFieldName_(propsDict_.lookup("velFieldName")),
    tempFieldName_(propsDict_.lookup("tempFieldName")),                             //temperature names/data
    partTempName_(propsDict_.lookup("partTempName")),
    partHeatFluxName_(propsDict_.lookupOrDefault<word>(      "partHeatFluxName", "na")),
    partHeatTransCoeffName_(propsDict_.lookupOrDefault<word>("partHeatTransCoeffName", "na")),
    partHeatFluidName_(propsDict_.lookupOrDefault<word>(     "partHeatFluidName", "na")),
    partDat_(NULL),
    partDatFlux_(NULL),
    partDatTransCoeff_(NULL),
    partDatFluid_(NULL),
    partTemp_(NULL),
    partDatTmpExpl_(NULL),
    partDatTmpImpl_(NULL),
    partDatSaturation_(NULL),
    partCoolingFlux_(NULL),
    validPartFlux_(false),
    validPartTransCoeff_(false),
    validPartFluid_(false),
    haveTemperatureEqn_(false),
    useLiMason_(false),
    lambda_(readScalar(propsDict_.lookup("lambda"))),
    Prandtl_(readScalar(propsDict_.lookup("Prandtl"))),
    eulerianFieldNames_( generalPropsDict_.lookup("eulerianFields")), 
    partSpeciesNames_(propsDict_.lookup("partSpeciesNames")),                       
    partSpeciesFluxNames_(propsDict_.lookup("partSpeciesFluxNames")),
    partSpeciesTransCoeffNames_(propsDict_.lookup("partSpeciesTransCoeffNames")),
    partSpeciesFluidNames_(propsDict_.lookup("partSpeciesFluidNames")),
    DMolecular_(propsDict_.lookup("DMolecular")),
    parameterVap_(propsDict_.lookup("parameterVap")),
    Rvap_(propsDict_.lookupOrDefault<scalar>("Rvap", 461.5)),
   // alphaImExSplit_(propsDict_.lookupOrDefault<scalar>("alphaImExSplit", 0.5)),
    deltaHEvap_("deltaHEvap", dimLength*dimLength/dimTime/dimTime, 1),
    maxSource_(1e30),
    scaleDia_(1.)
{
    allocateMyArrays(0.0);

    if (propsDict_.found("maxSource"))
    {
        maxSource_=readScalar(propsDict_.lookup ("maxSource"));
        Info << "limiting eulerian source field to: " << maxSource_ << endl;
    }

    if (propsDict_.found("useLiMason")) //TODO for Schmidt Number correlation 
    {
        useLiMason_=readBool(propsDict_.lookup ("useLiMason"));
        Info << "setting for useLiMason: " << useLiMason_ << endl;
    }
    if(useLiMason_)
        Nusselt=&scalarGeneralExchange::NusseltLiMason;
    else
        Nusselt=&scalarGeneralExchange::NusseltDeenEtAl;

    if(partHeatFluxName_!="na")        
    {   validPartFlux_=true;
        Info << "Found a valid partHeatFluxName: " << partHeatFluxName_ << endl;
    }
    if(partHeatTransCoeffName_!="na")
    {   validPartTransCoeff_=true;  
        Info << "Found a valid partHeatTransCoeffName: " << partHeatTransCoeffName_ << endl;
    }
    if(partHeatFluidName_!="na")
    {   validPartFluid_=true;
        Info << "Found a valid partHeatFluidName: " << partHeatFluidName_ << endl;
    }

    if( validPartTransCoeff_ && !validPartFluid_ )
        FatalError <<"Transfer coefficient set, but and fluid name missing. Check your entries in the couplingProperties! \n" 
                   << abort(FatalError);    

    if( !validPartTransCoeff_ && validPartFluid_ )
        FatalError <<"Fluid name set, but transfer coefficient missing. Check your entries in the couplingProperties! \n" 
                   << abort(FatalError);    
    
    if(!validPartFlux_ && !(validPartTransCoeff_ && validPartFluid_) )
        FatalError <<"You must set a valid heat flux name, or a valid transfer coefficient and fluid name \n" 
                   << abort(FatalError);

    // init force sub model
    setForceSubModels(propsDict_);

    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(3,true); // activate search for verbose switch
    forceSubM(0).setSwitchesList(4,true); // activate search for interpolate switch
    forceSubM(0).setSwitchesList(8,true); // activate scalarViscosity switch

    // read those switches defined above, if provided in dict
    for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
        forceSubM(iFSub).readSwitches();

    particleCloud_.checkCG(true);

    if (propsDict_.found("scale"))
        scaleDia_=scalar(readScalar(propsDict_.lookup("scale")));

    //check species names
    Info << "scalarGeneralExchange found the following eulerianFieldName: " << eulerianFieldNames_ << endl;
    int numTempEqn=0;
    for(int iEul=0;iEul<eulerianFieldNames_.size(); iEul++)
        if(eulerianFieldNames_[iEul]==tempFieldName_)
        {
            haveTemperatureEqn_ = true;
            numTempEqn = 1;
            Info << "scalarGeneralExchange identified the eulerianField '" << tempFieldName_ 
                 << "' as being the temperature field" << endl;
        }

    //check if enough particle properties have been provided
    if(partSpeciesNames_.size()!=(eulerianFieldNames_.size()-numTempEqn))   
        FatalError <<"Not enough partSpeciesNames specified in the couplingProperties file. \n" 
                   << abort(FatalError);
    else
        Info << "Found valid partSpeciesNames: " << partSpeciesNames_ << endl;
    if(partSpeciesFluxNames_.size()!=(eulerianFieldNames_.size()-numTempEqn))   
        FatalError <<"Not enough partSpeciesFluxNames specified in the couplingProperties file. \n" 
                   << abort(FatalError);
    else
        Info << "Found valid partSpeciesFluxNames: " << partSpeciesFluxNames_ << endl;
    if(partSpeciesTransCoeffNames_.size()!=(eulerianFieldNames_.size()-numTempEqn))   
        FatalError <<"Not enough partSpeciesTransCoeffNames specified in the couplingProperties file. \n" 
                   << abort(FatalError);
    else
        Info << "Found valid partSpeciesTransCoeffNames: " << partSpeciesTransCoeffNames_ << endl;
    if(partSpeciesFluidNames_.size()!=(eulerianFieldNames_.size()-numTempEqn))   
        FatalError <<"Not enough partSpeciesFluidNames specified in the couplingProperties file. \n" 
                   << abort(FatalError);
    else
        Info << "Found valid partSpeciesFluidNames: " << partSpeciesFluidNames_ << endl;
    if(DMolecular_.size()!=(eulerianFieldNames_.size()-numTempEqn))   
        FatalError <<"Not enough DMolecular specified in the couplingProperties file. \n" 
                   << abort(FatalError);

    for(int iPart=0;iPart<partSpeciesNames_.size(); iPart++)
    {
        if(partSpeciesNames_[iPart]=="none")
            particleSpeciesValue_.setSize(particleSpeciesValue_.size()+1, -1); //will not consider this species for coupling
        else if(partSpeciesNames_[iPart]=="zero")
            particleSpeciesValue_.setSize(particleSpeciesValue_.size()+1, 0.0); //will not consider this species for coupling
        else
        {
            particleSpeciesValue_.setSize(particleSpeciesValue_.size()+1, 2*ALARGECONCENTRATION); //set to a very large value to request pull
        }
    }

    if(nrForceSubModels()>1)
        FatalError <<"nrForceSubModels() must be zero or one for scalarGeneralExchange. \n" 
                   << abort(FatalError);

    if (probeIt_ && propsDict_.found("suppressProbe"))
        probeIt_=!Switch(propsDict_.lookup("suppressProbe"));
    if(probeIt_)
    {
        particleCloud_.probeM().initialize(typeName, "scalarGeneralExchange.logDat");
        particleCloud_.probeM().vectorFields_.append("Urel");               //first entry must the be the vector to probe
        particleCloud_.probeM().scalarFields_.append("Nu");                 //other are debug
        particleCloud_.probeM().scalarFields_.append("Rep");                //other are debug
        particleCloud_.probeM().writeHeader();
    }

}

//Todo for schmidt number???


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

scalarGeneralExchange::~scalarGeneralExchange()
{
    particleCloud_.dataExchangeM().destroy(partDat_,1);
    particleCloud_.dataExchangeM().destroy(partDatFlux_,1);

    particleCloud_.dataExchangeM().destroy(partDatTmpExpl_,1);
    particleCloud_.dataExchangeM().destroy(partDatTmpImpl_,1);
    particleCloud_.dataExchangeM().destroy(partDatSaturation_,1);

    if(validPartTransCoeff_)
        particleCloud_.dataExchangeM().destroy(partDatTransCoeff_,1);

    if(validPartFluid_)
        particleCloud_.dataExchangeM().destroy(partDatFluid_,1);

}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
void scalarGeneralExchange::allocateMyArrays(scalar initVal) const
{
    // Heat - get memory for 2d arrays
    particleCloud_.dataExchangeM().allocateArray(partDat_,initVal,1);  // field/initVal/with/lenghtFromLigghts
    particleCloud_.dataExchangeM().allocateArray(partDatFlux_,initVal,1);
    particleCloud_.dataExchangeM().allocateArray(partDatTmpExpl_,initVal,1);
    particleCloud_.dataExchangeM().allocateArray(partDatTmpImpl_,initVal,1);

    if(validPartTransCoeff_)
        particleCloud_.dataExchangeM().allocateArray(partDatTransCoeff_,initVal,1);    
    if(validPartFluid_)
        particleCloud_.dataExchangeM().allocateArray(partDatFluid_,initVal,1);
}
// * * * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * * //

void scalarGeneralExchange::setForce() const
{
    // do nothing
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void scalarGeneralExchange::manipulateScalarField(volScalarField& explicitEulerSource,
                                                  volScalarField& implicitEulerSource,
                                                  int speciesID) const
{

    // reset Scalar field
    explicitEulerSource.internalField() = 0.0;
    implicitEulerSource.internalField() = 0.0;

    if(speciesID>=0 && particleSpeciesValue_[speciesID]<0.0)    //skip if species is not active
        return;

    //Set the names of the exchange fields
    word    fieldName;
    word    partDatName;
    word    partFluxName;
    word    partTransCoeffName;
    word    partFluidName;
    scalar  transportParameter;

    // realloc the arrays to pull particle data


    if(speciesID<0) //this is the temperature - always pull from LIGGGHTS
    {
        fieldName          = tempFieldName_;
        partDatName        = partTempName_;
        partFluxName       = partHeatFluxName_;
        partTransCoeffName = partHeatTransCoeffName_;
        partFluidName      = partHeatFluidName_;
        transportParameter = lambda_;

        allocateMyArrays(0.0);
        particleCloud_.dataExchangeM().getData(partDatName,"scalar-atom", partDat_);
    }
    else
    {
        fieldName          = eulerianFieldNames_[speciesID];
        partDatName        = partSpeciesNames_[speciesID];
        partFluxName       = partSpeciesFluxNames_[speciesID]; 
        partTransCoeffName = partSpeciesTransCoeffNames_[speciesID]; 
        partFluidName      = partSpeciesFluidNames_[speciesID]; 
        transportParameter = DMolecular_[speciesID];

        allocateMyArrays(0.0);
        if(particleSpeciesValue_[speciesID]>ALARGECONCENTRATION)   
            particleCloud_.dataExchangeM().getData(partDatName,"scalar-atom", partDat_);
    }

    if (scaleDia_ > 1)
        Info << typeName << " using scale = " << scaleDia_ << endl;
    else if (particleCloud_.cg() > 1)
    {
        scaleDia_=particleCloud_.cg();
        Info << typeName << " using scale from liggghts cg = " << scaleDia_ << endl;
    }

    //==============================
    // get references
    const volScalarField& voidfraction_(particleCloud_.mesh().lookupObject<volScalarField> (voidfractionFieldName_));    // ref to voidfraction field
    const volVectorField& U_(particleCloud_.mesh().lookupObject<volVectorField> (velFieldName_));
    const volScalarField& fluidScalarField_(particleCloud_.mesh().lookupObject<volScalarField> (fieldName));            // ref to scalar field
    const volScalarField& nufField = forceSubM(0).nuField();
    //==============================

    // calc La based heat flux
    vector position(0,0,0);
    scalar voidfraction(1);
    vector Ufluid(0,0,0);
    scalar fluidValue(0);
    label  cellI=0;
    vector Us(0,0,0);
    vector Ur(0,0,0);
    scalar dscaled(0);
    scalar nuf(0);
    scalar magUr(0);
    scalar As(0);
    scalar Rep(0);
    scalar Pr(0);

    scalar sDth(scaleDia_*scaleDia_*scaleDia_);

    interpolationCellPoint<scalar> voidfractionInterpolator_(voidfraction_);
    interpolationCellPoint<vector> UInterpolator_(U_);
    interpolationCellPoint<scalar> fluidScalarFieldInterpolator_(fluidScalarField_);

    #include "setupProbeModel.H"

    for(int index = 0;index < particleCloud_.numberOfParticles(); ++index)
    {
            cellI = particleCloud_.cellIDs()[index][0];
            if(cellI >= 0)
            {
                if(forceSubM(0).interpolation())
                {
	                position = particleCloud_.position(index);
                    voidfraction = voidfractionInterpolator_.interpolate(position,cellI);
                    Ufluid = UInterpolator_.interpolate(position,cellI);
                    fluidValue = fluidScalarFieldInterpolator_.interpolate(position,cellI);
                }else
                {
					voidfraction = voidfraction_[cellI];
                    Ufluid       = U_[cellI];
                    fluidValue   = fluidScalarField_[cellI];
                }

                // calc relative velocity
                Us      = particleCloud_.velocity(index);
                Ur      = Ufluid-Us;
                magUr   = mag(Ur);
                dscaled = 2*particleCloud_.radius(index)/scaleDia_;
                As      = dscaled*dscaled*M_PI*sDth;
                nuf     = nufField[cellI];
                Rep     = dscaled*magUr/nuf;
                if(speciesID<0) //have temperature
                    Pr      = Prandtl_; 
                else
                    Pr      = max(SMALL,nuf/transportParameter); //This is Sc for species

                scalar alpha = transportParameter*(this->*Nusselt)(Rep,Pr,voidfraction)/(dscaled);

                // calc convective heat flux [W]
                scalar areaTimesTransferCoefficient = alpha * As;
                scalar tmpPartFlux     =  areaTimesTransferCoefficient 
                                       * (fluidValue - partDat_[index][0]);
                partDatFlux_[index][0] = tmpPartFlux;

                // split implicit/explicit contribution
                forceSubM(0).explicitCorrScalar( partDatTmpImpl_[index][0], 
                                                 partDatTmpExpl_[index][0],
                                                 areaTimesTransferCoefficient,
                                                 fluidValue,
                                                 fluidScalarField_[cellI],
                                                 partDat_[index][0],
                                                 forceSubM(0).verbose()
                                               );

                if(validPartTransCoeff_)
                    partDatTransCoeff_[index][0] = alpha;

                if(validPartFluid_)
                    partDatFluid_[index][0]      = fluidValue;


                if( forceSubM(0).verbose())
                {
                    Pout << "fieldName = " << fieldName << endl;
                    Pout << "partTransCoeffName = " << partTransCoeffName << endl;
                    Pout << "index    = " <<index << endl;
                    Pout << "partFlux = " << tmpPartFlux << endl;
                    Pout << "magUr = " << magUr << endl;
                    Pout << "As = " << As << endl;
                    Pout << "r = " << particleCloud_.radius(index) << endl;
                    Pout << "dscaled = " << dscaled << endl;
                    Pout << "nuf = " << nuf << endl;
                    Pout << "Rep = " << Rep << endl;
                    Pout << "Pr/Sc = " << Pr << endl;
                    Pout << "Nup/Shp = " << (this->*Nusselt)(Rep,Pr,voidfraction) << endl;
                    Pout << "partDatTransCoeff: " <<  partDatTransCoeff_[index][0] << endl;
                    Pout << "voidfraction = " << voidfraction << endl;
                    Pout << "partDat_[index][0] = " << partDat_[index][0] << endl  ;
                    Pout << "fluidValue = " << fluidValue << endl  ;
                }
                
                //Set value fields and write the probe
                if(probeIt_)
                {
                    #include "setupProbeModelfields.H"
                    // Note: for other than ext one could use vValues.append(x)
                    // instead of setSize
                    vValues.setSize(vValues.size()+1, Ur);           //first entry must the be the force
                    sValues.setSize(sValues.size()+1, (this->*Nusselt)(Rep,Pr,voidfraction)); 
                    sValues.setSize(sValues.size()+1, Rep);
                    particleCloud_.probeM().writeProbe(index, sValues, vValues);
                }
            }
    }

    //Handle explicit and implicit source terms on the Euler side
    //these are simple summations!
    particleCloud_.averagingM().setScalarSum
    (
        explicitEulerSource,
        partDatTmpExpl_,
        particleCloud_.particleWeights(),
        NULL
    );

    particleCloud_.averagingM().setScalarSum
    (
        implicitEulerSource,
        partDatTmpImpl_,
        particleCloud_.particleWeights(),
        NULL
    );

    // scale with the cell volume to get (total) volume-specific source 
    explicitEulerSource.internalField() /= -explicitEulerSource.mesh().V();
    implicitEulerSource.internalField() /= -implicitEulerSource.mesh().V();

    // limit explicit source term
    scalar explicitEulerSourceInCell;
    forAll(explicitEulerSource,cellI)
    {
        explicitEulerSourceInCell = explicitEulerSource[cellI];

        if(mag(explicitEulerSourceInCell) > maxSource_ )
        {
             explicitEulerSource[cellI] = sign(explicitEulerSourceInCell) * maxSource_;
        }
    }

    if(speciesID<0) //have temperature
        Info << "total convective particle-fluid heat flux [W] (Eulerian) = " 
             << gSum((explicitEulerSource+implicitEulerSource*fluidScalarField_)*explicitEulerSource.mesh().V()) 
             << endl;
    else
        Info << "speciesID: " << speciesID 
             << ": total convective particle-fluid species flux [kmol/s] (Eulerian) = " 
             << gSum((explicitEulerSource+implicitEulerSource*fluidScalarField_)*explicitEulerSource.mesh().V()) 
             << endl;

    // give DEM data
    if(validPartFlux_)
        particleCloud_.dataExchangeM().giveData(partFluxName ,
                                                "scalar-atom", 
                                                partDatFlux_);
    if(validPartTransCoeff_)
        particleCloud_.dataExchangeM().giveData(partTransCoeffName,
                                                "scalar-atom", 
                                                partDatTransCoeff_);
                
    if(validPartFluid_)
        particleCloud_.dataExchangeM().giveData(partFluidName,
                                                "scalar-atom", 
                                                partDatFluid_);

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
double scalarGeneralExchange::NusseltLiMason(double Rep, double Pr, double voidfraction) const
{
    double h1(0);
    double h2(0);
    double Nup(0);
    if (Rep < 200.0)
    {
         Nup = 2.0
             + 0.6
             * voidfraction*voidfraction*voidfraction*sqrt(voidfraction) //voidfraction^3.5
             * sqrt(Rep)
             * pow(Pr,0.33333333333); //This is Sh for species
    }
    else if (Rep < 1500.0)
    {
          h1  = voidfraction*voidfraction*voidfraction*sqrt(voidfraction); //voidfraction^3.5
          h2  = pow(Pr,0.3333333333);
          Nup = 2.0
              + 0.5 *h1*sqrt(Rep)*h2
              + 0.02*h1*pow(Rep,0.8)*h2;
    }
    else
    {
          Nup = 2.0
              + 0.000045
              * voidfraction*voidfraction*voidfraction*sqrt(voidfraction) //voidfraction^3.5
              * pow(Rep,1.8);
    }

    return Nup;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
double scalarGeneralExchange::NusseltDeenEtAl(double Rep, double Pr, double voidfraction) const
{
    //WARNING: This function is fitted for Reynolds numbers between 0 and 100!!!
    double Nup(0);
    double PrPowOneThird  = pow(Pr,0.3333333333) ;
    Nup = 
       ( 7.0 - 10.0 * voidfraction + 5 * voidfraction * voidfraction ) 
     *
       (1.0 + 
        0.17 
       * pow(Rep,0.2) 
       * PrPowOneThird 
       )
     + 
     ( 1.33 - 2.31 * voidfraction + 1.16 * voidfraction * voidfraction ) 
     * pow(Rep,0.7)  
     * PrPowOneThird ;

    return Nup;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
