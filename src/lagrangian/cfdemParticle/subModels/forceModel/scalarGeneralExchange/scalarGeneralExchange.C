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
    partHeatFluxName_(propsDict_.lookupOrDefault<word>(      "partHeatFluxName", "none")),
    partHeatTransCoeffName_(propsDict_.lookupOrDefault<word>("partHeatTransCoeffName", "none")),
    partHeatFluidName_(propsDict_.lookupOrDefault<word>(     "partHeatFluidName", "none")),
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
    useGeneralCorrelation_(false),
    generalCorrelationParameters_(propsDict_.lookupOrDefault<scalarList>("generalCorrelationParameters",scalarList(0))),
    lambda_(readScalar(propsDict_.lookup("lambda"))),
    Prandtl_(readScalar(propsDict_.lookup("Prandtl"))),
    eulerianFieldNames_( generalPropsDict_.lookup("eulerianFields")), 
    partSpeciesNames_(propsDict_.lookup("partSpeciesNames")),                       
    partSpeciesFluxNames_(propsDict_.lookup("partSpeciesFluxNames")),
    partSpeciesTransCoeffNames_(propsDict_.lookup("partSpeciesTransCoeffNames")),
    partSpeciesFluidNames_(propsDict_.lookup("partSpeciesFluidNames")),
    DMolecular_(propsDict_.lookup("DMolecular")),
    maxSource_(1e30),
    scaleDia_(1.)
{
    setForceSubModels(propsDict_);
    setupModel();
    if (probeIt_ && propsDict_.found("suppressProbe"))
        probeIt_=!Switch(propsDict_.lookup("suppressProbe"));
    if(probeIt_)
    {
      forAll(eulerianFieldNames_, fieldIt)
      {
        particleCloud_.probeM().initialize(typeName, typeName + "_" + eulerianFieldNames_[fieldIt] + ".logDat");
        particleCloud_.probeM().vectorFields_.append("Urel");               //first entry must the be the vector to probe
        if(eulerianFieldNames_[fieldIt]==tempFieldName_) //this is the temperature
        {
            particleCloud_.probeM().scalarFields_.append("Rep");            
            particleCloud_.probeM().scalarFields_.append("Nu"); 
        }
        else                
            particleCloud_.probeM().scalarFields_.append("Sh"); 
        particleCloud_.probeM().scalarFields_.append("exchangeRate");
        particleCloud_.probeM().writeHeader();
      }
    }

    for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
        forceSubM(iFSub).constructorCalls(typeName);

    particleCloud_.setAllowCFDsubTimestep(false);
}

// ***********************************************************
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
    partHeatFluxName_(propsDict_.lookupOrDefault<word>(      "partHeatFluxName", "none")),
    partHeatTransCoeffName_(propsDict_.lookupOrDefault<word>("partHeatTransCoeffName", "none")),
    partHeatFluidName_(propsDict_.lookupOrDefault<word>(     "partHeatFluidName", "none")),
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
    useGeneralCorrelation_(false),
    lambda_(readScalar(propsDict_.lookup("lambda"))),
    Prandtl_(readScalar(propsDict_.lookup("Prandtl"))),
    eulerianFieldNames_( generalPropsDict_.lookup("eulerianFields")), 
    partSpeciesNames_(propsDict_.lookup("partSpeciesNames")),                       
    partSpeciesFluxNames_(propsDict_.lookup("partSpeciesFluxNames")),
    partSpeciesTransCoeffNames_(propsDict_.lookup("partSpeciesTransCoeffNames")),
    partSpeciesFluidNames_(propsDict_.lookup("partSpeciesFluidNames")),
    DMolecular_(propsDict_.lookup("DMolecular")),
    partHeatFluxPositionInRegister_(-1),
    partHeatTransCoeffPositionInRegister_(-1),
    partHeatFluidPositionInRegister_(-1),
    maxSource_(1e30),
    scaleDia_(1.)

{
    setForceSubModels(propsDict_);
    setupModel();
    if (probeIt_ && propsDict_.found("suppressProbe"))
        probeIt_=!Switch(propsDict_.lookup("suppressProbe"));
    if(probeIt_)
    {
      forAll(eulerianFieldNames_, fieldIt)
      {
        particleCloud_.probeM().initialize(dictName, dictName + "_" + eulerianFieldNames_[fieldIt] + ".logDat");
        particleCloud_.probeM().vectorFields_.append("Urel");               //first entry must the be the vector to probe
        if(eulerianFieldNames_[fieldIt]==tempFieldName_) //this is the temperature
        {
            particleCloud_.probeM().scalarFields_.append("Rep");            
            particleCloud_.probeM().scalarFields_.append("Nu");                 //other are debug
        }
        else                
            particleCloud_.probeM().scalarFields_.append("Sh"); 
        particleCloud_.probeM().scalarFields_.append("exchangeRate");
        particleCloud_.probeM().writeHeader();
      }
    }

    for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
        forceSubM(iFSub).constructorCalls(dictName);

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

scalarGeneralExchange::~scalarGeneralExchange()
{
    particleCloud_.dataExchangeM().destroy(partDat_,1);
    particleCloud_.dataExchangeM().destroy(partDatTmpExpl_,1);
    particleCloud_.dataExchangeM().destroy(partDatTmpImpl_,1);

    //external particle data (e.g., fluxes, transCoeff, or fluid data) will be destroyed in cloud, 
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
void scalarGeneralExchange::allocateMyArrays(scalar initVal) const
{
    //Only allocate temporary Lagrangian arrays, correct length will be provided by ExternalCode
    if(particleCloud_.numberOfParticlesChanged())
    {
        particleCloud_.dataExchangeM().allocateArray(partDat_,initVal,1);  
        particleCloud_.dataExchangeM().allocateArray(partDatTmpExpl_,initVal,1);
        particleCloud_.dataExchangeM().allocateArray(partDatTmpImpl_,initVal,1);
    }
    //external particle data (e.g., fluxes, transCoeff, or fluid data) will be allocated in cloud, 
    //just need to set pointers before access

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

    // reset Scalar field (== means hard reset)
    explicitEulerSource == dimensionedScalar("zero", explicitEulerSource.dimensions(), 0.);
    implicitEulerSource == dimensionedScalar("zero", implicitEulerSource.dimensions(), 0.);

    if(speciesID>=0 && particleSpeciesValue_[speciesID]<0.0)    //skip if species is not active
        return;

    //Set the names of the exchange fields
    word    fieldName;
    word    partDatName;
    scalar  transportParameter;

    if(speciesID<0) //this is the temperature 
    {
        fieldName          = tempFieldName_;
        partDatName        = partTempName_;
        transportParameter = lambda_;

        setPointersToExternalArrays( partHeatFluxName_,       partHeatFluxPositionInRegister_,
                                     partHeatTransCoeffName_, partHeatTransCoeffPositionInRegister_,
                                     partHeatFluidName_,      partHeatFluidPositionInRegister_
                                   );
        if(probeIt_)
            particleCloud_.probeM().setOutputFile(typeName+"_"+tempFieldName_+".logDat");
    }
    else
    {
        fieldName          = eulerianFieldNames_[speciesID];
        partDatName        = partSpeciesNames_[speciesID];
        transportParameter = DMolecular_[speciesID];

        setPointersToExternalArrays( partSpeciesFluxNames_[speciesID],       partSpeciesFluxPositionInRegister_[speciesID],
                                     partSpeciesTransCoeffNames_[speciesID], partSpeciesTransCoeffPositionInRegister_[speciesID],
                                     partSpeciesFluidNames_[speciesID],      partSpeciesFluidPositionInRegister_[speciesID]
                                   );
        if(probeIt_)
            particleCloud_.probeM().setOutputFile(typeName + "_" + fieldName + ".logDat");
    }


    //==============================
    // get references
    const volScalarField& voidfraction_(particleCloud_.mesh().lookupObject<volScalarField> (voidfractionFieldName_));    // ref to voidfraction field
    const volVectorField& U_(particleCloud_.mesh().lookupObject<volVectorField> (velFieldName_));
    const volScalarField& fluidScalarField_(particleCloud_.mesh().lookupObject<volScalarField> (fieldName));            // ref to scalar field
    const volScalarField& nufField = forceSubM(0).nuField();
    //==============================

    if (scaleDia_ > 1)
        Info << typeName << " using scale = " << scaleDia_ << endl;
    else if (particleCloud_.cg() > 1)
    {
        scaleDia_=particleCloud_.cg();
        Info << typeName << " using scale from liggghts cg = " << scaleDia_ << endl;
    }

    // realloc the arrays and get data
    allocateMyArrays(0.0);
    if(speciesID<0)
        particleCloud_.dataExchangeM().getData(partDatName,"scalar-atom", partDat_);
    else 
      if(particleSpeciesValue_[speciesID]>ALARGECONCENTRATION)   //only pull if needed
            particleCloud_.dataExchangeM().getData(partDatName,"scalar-atom", partDat_);

    // calc La based heat flux
    vector position(0,0,0);
    scalar voidfraction(1);
    vector Ufluid(0,0,0);
    scalar fluidValue(0);
    label  cellI=0;
    vector Us(0,0,0);
    vector Ur(0,0,0);
    scalar dscaled(0);
    scalar dparcel(0);
    scalar numberParticlesInParcel(1);
    scalar nuf(0);
    scalar magUr(0);
    scalar As(0);
    scalar Rep(0);
    scalar Pr(0);
    scalar vCell(0);
    scalar cellDpRatio(1);

    #include "resetVoidfractionInterpolator.H"
    #include "resetUInterpolator.H"
    #include "resetFluidScalarFieldInterpolator.H"

    for(int index = 0;index < particleCloud_.numberOfParticles(); ++index)
    {
            cellI = particleCloud_.cellIDs()[index][0];
            if(cellI >= 0)
            {
                if(forceSubM(0).interpolation())
                {
	                position = particleCloud_.position(index);
                    voidfraction = voidfractionInterpolator_().interpolate(position,cellI);
                    Ufluid = UInterpolator_().interpolate(position,cellI);
                    fluidValue = fluidScalarFieldInterpolator_().interpolate(position,cellI);
                }else
                {
					voidfraction = voidfraction_[cellI];
                    Ufluid       = U_[cellI];
                    fluidValue   = fluidScalarField_[cellI];
                }

                dscaled = 2*particleCloud_.radius(index);
                dparcel = dscaled;
                if(forceSubM(0).useCorrectedVoidage())
                {
                    vCell = U_.mesh().V()[cellI];
                    cellDpRatio =  pow(vCell,0.3333333)/(dparcel);

                    for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
                         voidfraction = forceSubM(iFSub).calculateCorrectedVoidage(voidfraction,cellDpRatio); 
                }

                // calc relative velocity
                Us      = particleCloud_.velocity(index);
                Ur      = Ufluid-Us;
                magUr   = mag(Ur);

                forceSubM(0).scaleDia(dscaled,index); //caution: this fct will scale ds!
                numberParticlesInParcel    = dparcel/dscaled;
                numberParticlesInParcel   *= numberParticlesInParcel*numberParticlesInParcel;
                As      = dscaled*dscaled*M_PI*numberParticlesInParcel;
                nuf     = nufField[cellI];
                Rep     = dscaled*magUr*voidfraction/nuf; //MUST use superficial velocity here!
                if(speciesID<0) //have temperature
                    Pr      = Prandtl_; 
                else
                    Pr      = max(SMALL,nuf/transportParameter); //This is Sc for species

                scalar alpha = transportParameter*(this->*Nusselt)(Rep,Pr,voidfraction)/dscaled;

                // calc convective heat flux [W]
                scalar areaTimesTransferCoefficient = alpha * As;
                scalar tmpPartFlux     =  areaTimesTransferCoefficient 
                                       * (fluidValue - partDat_[index][0]);

                // split implicit/explicit contribution
                forceSubM(0).explicitCorrScalar( partDatTmpImpl_[index][0], 
                                                 partDatTmpExpl_[index][0],
                                                 areaTimesTransferCoefficient,
                                                 fluidValue,
                                                 fluidScalarField_[cellI],
                                                 partDat_[index][0],
                                                 forceSubM(0).verbose()
                                               );
                if(validPartFlux_)
                    partDatFlux_[index][0]+= tmpPartFlux; //MUST ADD total source for ALL particles in parcel

                if(validPartTransCoeff_)
                    partDatTransCoeff_[index][0]+= alpha; //MUST ADD total source coefficient here

                if(validPartFluid_)
                    partDatFluid_[index][0]      = fluidValue; //MUST NOT add here, since this is factor


                if( forceSubM(0).verbose())
                {
                    Pout << "fieldName = " << fieldName << endl;
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
                    Pout << "voidfraction = " << voidfraction << endl;
                    Pout << "fluidValue = " << fluidValue << endl  ;
                    Pout << "partDat_[index][0] = " << partDat_[index][0] << endl;
                    if(validPartFlux_)
                        Pout << "partDatFlux: "       <<  partDatFlux_[index][0] << endl;
                    if(validPartTransCoeff_)
                        Pout << "partDatTransCoeff: " <<  partDatTransCoeff_[index][0] << endl;
                }
                
                //Set value fields and write the probe
                if(probeIt_)
                {
                    #include "setupProbeModelfields.H"
                    // Note: for other than ext one could use vValues.append(x)
                    // instead of setSize
                    vValues.setSize(vValues.size()+1, Ur);
                    if(speciesID<0) //this is the temperature, then also report Rep 
                        sValues.setSize(sValues.size()+1, Rep);
                    sValues.setSize(sValues.size()+1, (this->*Nusselt)(Rep,Pr,voidfraction));
                    sValues.setSize(sValues.size()+1, tmpPartFlux);
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
    particleCloud_.makeSpecific(explicitEulerSource);
    explicitEulerSource*=-1;
    particleCloud_.makeSpecific(implicitEulerSource);
    implicitEulerSource*=-1;

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


    //Reporting of integral quantities
    //TODO: write to different file for speciesId>0
    Field<scalar> writeValues; bool writeDiskNow=forceSubM(0).verboseToDisk(); //must call 'verboseToDisk()' only once since this function is incremeting a counter!
    writeValues.clear();
    if( forceSubM(0).verbose() || writeDiskNow)
    {
	    scalar exchangeRate = gSum(-(explicitEulerSource
				                    +implicitEulerSource*fluidScalarField_)
			                        *explicitEulerSource.mesh().V()
					              );
        // Note: for other than ext one could use writeValues.append(x)
        // instead of setSize
        writeValues.setSize(writeValues.size()+1, exchangeRate);

        if(forceSubM(0).verbose())
        {
          if(speciesID<0) //have temperature
            Info << "total convective particle-fluid heat flux [W] (Eulerian) = " 
                << exchangeRate
                << endl;
          else
            Info << "speciesID: " << speciesID 
                 << ": total convective particle-fluid species flux [kmol/s] (or [kg/s]) (Eulerian) = " 
                 << exchangeRate
                 << endl;
        }
    }
    if( writeDiskNow )
      for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
        forceSubM(iFSub).verboseToDiskWrite(writeValues);
    
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
double scalarGeneralExchange::NusseltGeneralCorrelation(double Rep, double Pr, double voidfraction) const
{
    //WARNING: This function my be fitted to data for a limited range of Reynolds number!!
    double Nup(0);
    double PrPowOneThird  = pow(Pr,0.3333333333) ;
    Nup = 
       (  generalCorrelationParameters_[0] 
        + generalCorrelationParameters_[1] * voidfraction 
        + generalCorrelationParameters_[2] * voidfraction * voidfraction ) 
     *
       (  generalCorrelationParameters_[3] 
        + generalCorrelationParameters_[4]
       * pow(Rep,0.2) 
       * PrPowOneThird 
       )
     + 
     (  generalCorrelationParameters_[5]
      + generalCorrelationParameters_[6] * voidfraction 
      + generalCorrelationParameters_[7] * voidfraction * voidfraction ) 
     * pow(Rep,0.7)  
     * PrPowOneThird ;

    return Nup;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void scalarGeneralExchange::setupModel() const
{
    //just allocate arrays for internal use
    allocateMyArrays(0.0);

    //analyze dict to decide on allocation of external arrays
    //FLUXES (explicit coupling strategy)
    if(partHeatFluxName_!="none")
    {
        validPartFlux_=true;
        particleCloud_.registerNamesFieldsUserCFDEMToExt(propsDict_.lookup("partHeatFluxName"), partHeatFluxPositionInRegister_);
    }
    bool validspeciesFlux = 
        particleCloud_.checkAndregisterNamesFieldsUserCFDEMToExt(partSpeciesFluxNames_,
                                                                 partSpeciesFluxPositionInRegister_);

    if( (validPartFlux_ && !validspeciesFlux && partSpeciesFluxNames_.size()>0) || (!validPartFlux_ && validspeciesFlux ))    
        FatalError <<  "scalarGeneralExchange::setupModel: you have set a valid species flux name, but a non-valid heatflux name (or vice versa). This will mess up memory allocation. Please use both valid or non-valid flux names" << abort(FatalError);

    //TRANSCOEFF and FLUID (implicit coupling strategy)
    if(partHeatTransCoeffName_!="none")    {   
        validPartTransCoeff_=true;  
        Info << "Found a valid partHeatTransCoeffName: " << partHeatTransCoeffName_ << endl;
        particleCloud_.registerNamesFieldsUserCFDEMToExt(partHeatTransCoeffName_, partHeatTransCoeffPositionInRegister_);
    }
    if(partHeatFluidName_!="none")    {   
        validPartFluid_=true;
        Info << "Found a valid partHeatFluidName: " << partHeatFluidName_ << endl;
        particleCloud_.registerNamesFieldsUserCFDEMToExt(partHeatFluidName_, partHeatFluidPositionInRegister_);
    }
    if( validPartTransCoeff_ && !validPartFluid_ )
        FatalError <<"Transfer coefficient set, but and fluid name missing. Check your entries in the couplingProperties! \n" 
                   << abort(FatalError);    
    if( !validPartTransCoeff_ && validPartFluid_ )
        FatalError <<"Fluid name set, but transfer coefficient missing. Check your entries in the couplingProperties! \n" 
                   << abort(FatalError);    
    
    bool validspeciesTransCoeff = 
        particleCloud_.checkAndregisterNamesFieldsUserCFDEMToExt(partSpeciesTransCoeffNames_,
                                                                 partSpeciesTransCoeffPositionInRegister_);
    bool validspeciesFluid = 
        particleCloud_.checkAndregisterNamesFieldsUserCFDEMToExt(partSpeciesFluidNames_,
                                                                 partSpeciesFluidPositionInRegister_);

    if( (validPartTransCoeff_ && !validspeciesTransCoeff && partSpeciesTransCoeffNames_.size()>0 ) || (!validPartTransCoeff_ && validspeciesTransCoeff ))    
        FatalError <<  "scalarGeneralExchange::setupModel: you have set a valid species transCoeff name, but a non-valid heat trans coeff name (or vice versa). This will mess up memory allocation. Please use both valid or non-valid transCoeff names. You might want to use 'none' in partSpeciesTransCoeffNames to de-activate the species transCoeff name!" << abort(FatalError);

    if( (validPartFluid_ && !validspeciesFluid && partSpeciesFluidNames_.size()>0) || (!validPartFluid_ && validspeciesFluid ))    
        FatalError <<  "scalarGeneralExchange::setupModel: you have set a valid species fluid name, but a non-valid heat fluid name (or vice versa). This will mess up memory allocation. Please use both valid or non-valid fluid names. You might want to use 'none' in partSpeciesFluidNames to de-activate the species fluid name!" << abort(FatalError);


    //FINALY check and Info
    if(validPartFlux_ && validspeciesFlux)    {
        Info << "Found a valid partHeatFluxName and partSpeciesFluxName: " << partHeatFluxName_ << endl;
        Info << "scalarGeneralExchange (or derived model) will now proceed with EXPLICIT flux coupling " << endl;
    }
    else if(validPartTransCoeff_ && validPartFluid_)
    {
        Info << "Found a valid partHeatTransCoeffName and partHeatFluidfName (& corresponding species names)." << endl;
        Info << "scalarGeneralExchange (or derived model) will now proceed with IMPLICIT flux coupling " << endl;
    }
    else if(validPartFlux_ )    {
        Info << "Found a valid partHeatFluxName: " << partHeatFluxName_ << endl;
        Info << "scalarGeneralExchange (or derived model) will now proceed with EXPLICIT flux coupling for heat. " << endl;
    }
    else
        FatalError <<  "scalarGeneralExchange::setupModel: you did not specify a valid flux or transCoeff name set. Please either specify valid flux names, or valid transCoeff and fluid names." << abort(FatalError);


    //Make internal settings
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
    if (propsDict_.found("useGeneralCorrelation"))
    {
        useGeneralCorrelation_=readBool(propsDict_.lookup ("useGeneralCorrelation"));
        Info << "setting for useGeneralCorrelation: " << useGeneralCorrelation_ << endl;
    }
    
    if(useLiMason_ && useGeneralCorrelation_)
        FatalError <<"You cannot set 'useLiMason' AND 'useGeneralCorrelation' to true. Just set one seeting to true.  \n" 
                   << abort(FatalError);

    if(useLiMason_)
        Nusselt=&scalarGeneralExchange::NusseltLiMason;
    else if(useGeneralCorrelation_)
    {
    	//generalCorrelationParameters_ = propsDict_.lookup("generalCorrelationParameters");
        if(generalCorrelationParameters_.size()<8 || generalCorrelationParameters_.size()>8)
            FatalError <<"The data array specified as 'generalCorrelationParameters' is too short or too long. Must specify exactly 8 values. You specified: \n" 
                       << generalCorrelationParameters_ << endl
                       << abort(FatalError);
        Nusselt=&scalarGeneralExchange::NusseltGeneralCorrelation;
    }
    else
        Nusselt=&scalarGeneralExchange::NusseltDeenEtAl;

    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(3,true); // activate search for verbose switch
    forceSubM(0).setSwitchesList(4,true); // activate search for interpolate switch
    forceSubM(0).setSwitchesList(8,true); // activate scalarViscosity switch
    forceSubM(0).setSwitchesList(9,true); // activate verboseToDisk switch
    forceSubM(0).setSwitchesList(10,true); // activate correctedVoidage switch

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
            particleSpeciesValue_.setSize(particleSpeciesValue_.size()+1, 0.0);//will set this species to zero for coupling
        else
            particleSpeciesValue_.setSize(particleSpeciesValue_.size()+1, 2*ALARGECONCENTRATION);//set to a very large value to request pull
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void scalarGeneralExchange::setPointersToExternalArrays(    word nameFlux,          int positionFlux,
                                                            word nameTransCoeff,    int positionTransCoeff,
                                                            word nameFluid,         int positionFluid
                                                       ) const
{
    if(particleCloud_.particleDatFieldsUserCFDEMToExt.size() !=
               particleCloud_.namesFieldsUserCFDEMToExt.size()
      )
      FatalError <<  "\n\n****CATASTROPHIC ERROR MOST LIKELY CAUSED BY USER!!! \n" 
                 <<  "particleCloud_.particleDatFieldsUserCFDEMToExt.size() is NOT EQUAL to particleCloud_.namesFieldsUserCFDEMToExt.size()." 
                 <<  "This may be caused by an incorrect time step, or coupling interval, resulting in the fact that an array inside particleCloud_ was not allocated. "
                 <<  "Please check your time step and coupling interval settings, such that fluid-particle coupling is done EVERY fluid time step. \n\n" 
                 << abort(FatalError);
	
    if(validPartFlux_)  //EXPLICIT coupling strategy for Lagrangian part
    {
            partDatFlux_           = particleCloud_.particleDatFieldsUserCFDEMToExt[positionFlux];
            if( forceSubM(0).verbose() )    {
                Info <<"scalarParticleFilter will push particle fluxes to '" << nameFlux 
                     << "' (register:" << positionFlux <<  "). \n" << endl;
            }
    }

    if(validPartTransCoeff_) //IMPLICIT coupling strategy for Lagrangian part
    {
            partDatTransCoeff_     = particleCloud_.particleDatFieldsUserCFDEMToExt[positionTransCoeff];
            partDatFluid_          = particleCloud_.particleDatFieldsUserCFDEMToExt[positionFluid];
            if( forceSubM(0).verbose() )    {
                Info <<"scalarParticleFilter will push transfer coefficients to '" << nameTransCoeff 
                     << "' (register:" << positionTransCoeff <<  "). \n" << endl;
                Info <<"scalarParticleFilter will push fluid values to '" << nameFluid 
                     << "' (register:" << positionFluid <<  "). \n" << endl;
            }
    }

}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
