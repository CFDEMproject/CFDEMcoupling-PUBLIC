/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
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

#include "LaEuScalarTemp.H"
#include "addToRunTimeSelectionTable.H"
#include "dataExchangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(LaEuScalarTemp, 0);

addToRunTimeSelectionTable
(
    forceModel,
    LaEuScalarTemp,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
LaEuScalarTemp::LaEuScalarTemp
(
    const dictionary& dict,
    cfdemCloud& sm,
    word name
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(name == "" ? typeName + "Props" : name + "Props")),
    tempFieldName_(propsDict_.lookup("tempFieldName")),
    T_(sm.mesh().lookupObject<volScalarField> (tempFieldName_)),
    voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    maxSource_(1e30),
    velFieldName_(propsDict_.lookup("velFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    partTempName_(propsDict_.lookup("partTempName")),
    partTemp_(NULL),
    partHeatFluxName_(propsDict_.lookupOrDefault<word>(      "partHeatFluxName", "na")),
    partHeatFlux_(NULL),
    validPartHeatFlux_(false),
    partHeatTransCoeffName_(propsDict_.lookupOrDefault<word>("partHeatTransCoeffName", "na")),
    partHeatTransCoeff_(NULL),
    validPartHeatTransCoeff_(false),
    partHeatFluidName_(propsDict_.lookupOrDefault<word>(     "partHeatFluidName", "na")),
    partHeatFluid_(NULL),
    validPartHeatFluid_(false),
    lambda_(readScalar(propsDict_.lookup("lambda"))),
    Cp_(readScalar(propsDict_.lookup("Cp"))),
    compressible_(propsDict_.lookupOrDefault<Switch>("compressible",false)),
    doExchange_(propsDict_.lookupOrDefault<Switch>("doExchange",true))
{
    if (lambda_<SMALL)
    {
        Warning << "you tried using an invalid value for lambda - using SMALL instead!" << endl;
        lambda_=SMALL;
    }

    allocateMyArrays();

    if (propsDict_.found("maxSource"))
    {
        maxSource_=readScalar(propsDict_.lookup ("maxSource"));
        Info << "limiting eulerian source field to: " << maxSource_ << endl;
    }
    
    if(compressible_) Info << "Compressibility mode enabled." << endl;

    if(partHeatFluxName_!="na") validPartHeatFlux_=true;
    if(partHeatTransCoeffName_!="na") validPartHeatTransCoeff_=true;
    if(partHeatFluidName_!="na") validPartHeatFluid_=true;

    if( validPartHeatTransCoeff_ && !validPartHeatFluid_ )
        FatalError <<"Transfer coefficient set, but and fluid name missing. Check your entries in the couplingProperties! \n" 
                   << abort(FatalError);    

    if( !validPartHeatTransCoeff_ && validPartHeatFluid_ )
        FatalError <<"Fluid name set, but transfer coefficient missing. Check your entries in the couplingProperties! \n" 
                   << abort(FatalError);    
    
    if(!validPartHeatFlux_ && !(validPartHeatTransCoeff_ && validPartHeatFluid_) )
        FatalError <<"You must set a valid heat flux name, or a valid transfer coefficient and fluid name \n" 
                   << abort(FatalError);

    // init force sub model
    setForceSubModels(propsDict_);

    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(3,true); // activate search for verbose switch
    forceSubM(0).setSwitchesList(4,true); // activate search for interpolate switch
    forceSubM(0).setSwitchesList(8,true); // activate scalarViscosity switch

    // read those switches defined above, if provided in dict
    forceSubM(0).readSwitches();
    //for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
    //    forceSubM(iFSub).readSwitches();

    particleCloud_.checkCG(true);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

LaEuScalarTemp::~LaEuScalarTemp()
{
    delete partTemp_;
    delete partHeatFlux_;

    if(validPartHeatTransCoeff_)
       delete partHeatTransCoeff_;

    if(validPartHeatFluid_)
        delete partHeatFluid_;

}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
void LaEuScalarTemp::allocateMyArrays() const
{
    // get memory for 2d arrays
    double initVal=0.0;
    particleCloud_.dataExchangeM().allocateArray(partTemp_,initVal,1);  // field/initVal/with/lenghtFromLigghts
    particleCloud_.dataExchangeM().allocateArray(partHeatFlux_,initVal,1);

    if(validPartHeatTransCoeff_)
        particleCloud_.dataExchangeM().allocateArray(partHeatTransCoeff_,initVal,1);    

    if(validPartHeatFluid_)
        particleCloud_.dataExchangeM().allocateArray(partHeatFluid_,initVal,1);
}
// * * * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * * //

void LaEuScalarTemp::setForce() const
{
    // do nothing
}

void LaEuScalarTemp::manipulateScalarField(volScalarField& EuField) const
{
    // realloc the arrays
    allocateMyArrays();
    if(doExchange_)
    {
        // reset Scalar field (== means hard reset)
        EuField == dimensionedScalar("zero", EuField.dimensions(), 0.);

        // get DEM data
        particleCloud_.dataExchangeM().getData(partTempName_,"scalar-atom",partTemp_);

        const volScalarField& nufField = forceSubM(0).nuField();
        const volScalarField& rhoField = forceSubM(0).rhoField();

        // calc La based heat flux
        vector position(0,0,0);
        scalar voidfraction(1);
        vector Ufluid(0,0,0);
        scalar Tfluid(0);
        label cellI=0;
        vector Us(0,0,0);
        scalar ds(0);
        scalar dparcel(0);
        scalar numberParticlesInParcel(1);
        scalar nuf(0);
        scalar rho(0);
        scalar magUr(0);
        scalar As(0);
        scalar Rep(0);
        scalar Pr(0);
        scalar Nup(0);
        scalar n = 3.5; // model parameter

        #include "resetVoidfractionInterpolator.H"
        #include "resetUInterpolator.H"
        #include "resetTInterpolator.H"

        for(int index = 0;index < particleCloud_.numberOfParticles(); ++index)
        {
            //if(particleCloud_.regionM().inRegion()[index][0])
            //{
                cellI = particleCloud_.cellIDs()[index][0];
                if(cellI >= 0)
                {
                    if(forceSubM(0).interpolation())
                    {
                        position = particleCloud_.position(index);
                        voidfraction = voidfractionInterpolator_().interpolate(position,cellI);
                        Ufluid = UInterpolator_().interpolate(position,cellI);
                        Tfluid = TInterpolator_().interpolate(position,cellI);
                    }else
                    {
                        voidfraction = voidfraction_[cellI];
                        Ufluid = U_[cellI];
                        Tfluid = T_[cellI];
                    }

                    // calc relative velocity
                    Us = particleCloud_.velocity(index);
                    ds = 2*particleCloud_.radius(index);
                    dparcel = ds;
                    forceSubM(0).scaleDia(ds,index); //caution: this fct will scale ds!
                    numberParticlesInParcel    = dparcel/ds;
                    numberParticlesInParcel   *= numberParticlesInParcel*numberParticlesInParcel;
                    nuf = nufField[cellI];
                    rho = rhoField[cellI];

                    //Update any scalar or vector quantity
                    for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
                          forceSubM(iFSub).update(  index, 
                                                    cellI,
                                                    ds,
                                                    Ufluid, 
                                                    Us, 
                                                    nuf,
                                                    rho,
                                                    forceSubM(0).verbose()
                                                 );

                    magUr = mag(Ufluid-Us);
                    As = ds*ds*M_PI*numberParticlesInParcel;

                    Rep = ds*magUr/nuf;
                    Pr = max(SMALL,Cp_*nuf*rho/lambda_);

                    if (Rep < 200)
                    {
                        Nup = 2+0.6*pow(voidfraction,n)*sqrt(Rep)*pow(Pr,0.33);
                    }
                    else if (Rep < 1500)
                    {
                        Nup = 2. + (0.5 * sqrt(Rep) + 0.02 * pow(Rep,0.8)) * pow(voidfraction,n) * pow(Pr,0.33);
                    }
                    else
                    {
                        Nup = 2+0.000045*pow(voidfraction,n)*pow(Rep,1.8);
                    }
                    scalar h = lambda_*Nup/(ds);

                    // calc convective heat flux [W]
                    scalar partHeatFlux     = h * As * (Tfluid - partTemp_[index][0]);
                    partHeatFlux_[index][0] = partHeatFlux;

                    if(validPartHeatTransCoeff_)
                        partHeatTransCoeff_[index][0] = h;

                    if(validPartHeatFluid_)
                        partHeatFluid_[index][0]      = Tfluid;


                    if(forceSubM(0).verbose() && index >=0 && index <2)
                    {
                        Pout << "magUr = " << magUr << endl;
                        Pout << "As = " << As << endl;
                        Pout << "r = " << particleCloud_.radius(index) << endl;
                        Pout << "dprim = " << ds << endl;
                        Pout << "nuf = " << nuf << endl;
                        Pout << "Rep = " << Rep << endl;
                        Pout << "Pr = " << Pr << endl;
                        Pout << "Nup = " << Nup << endl;
                        Pout << "voidfraction = " << voidfraction << endl;
                        Pout << "partTemp_[index][0] = " << partTemp_[index][0] << endl;
                        Pout << "partHeatFlux_[index][0] = " << partHeatFlux_[index][0] << endl;
                        Pout << "Tfluid = " << Tfluid << endl  ;
                    }
                }
            //}
        }

        particleCloud_.averagingM().setScalarSum
        (
            EuField,
            partHeatFlux_,
            particleCloud_.particleWeights(),
            NULL
        );

        // scale with -1./(Vcell*rho*Cp)
        
        if(compressible_)
        {
            forAll(EuField,cellI)
                EuField[cellI] /= -EuField.mesh().V()[cellI];
        } else
        {
            forAll(EuField,cellI)
                EuField[cellI] /= -rhoField[cellI]*Cp_*EuField.mesh().V()[cellI];
        }

        // limit source term
        scalar EuFieldInCell;
        forAll(EuField,cellI)
        {
            EuFieldInCell = EuField[cellI];

            if(mag(EuFieldInCell) > maxSource_ )
            {
                 EuField[cellI] = sign(EuFieldInCell) * maxSource_;
            }
        }

        if(compressible_)
        {
            Info << "total convective particle-fluid heat flux [W] (Eulerian) = " 
                 << gSum(EuField*1*EuField.mesh().V()) 
                 << endl;
        } else
        {
            Info << "total convective particle-fluid heat flux [W] (Eulerian) = " 
                 << gSum(EuField*rhoField*Cp_*EuField.mesh().V()) 
                 << endl;
        }
    }
}

void LaEuScalarTemp::commToDEM() const
{
    if(doExchange_)
    {
        // give DEM data
        if(validPartHeatFlux_)
            particleCloud_.dataExchangeM().giveData(partHeatFluxName_,
                                                    "scalar-atom", 
                                                    partHeatFlux_);
            
        if(validPartHeatTransCoeff_)
            particleCloud_.dataExchangeM().giveData(partHeatTransCoeffName_,
                                                    "scalar-atom", 
                                                    partHeatTransCoeff_);
                    
        if(validPartHeatFluid_)
            particleCloud_.dataExchangeM().giveData(partHeatFluidName_,
                                                    "scalar-atom", 
                                                    partHeatFluid_);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
