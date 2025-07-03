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
    compressible_(propsDict_.lookupOrDefault<Switch>("compressible",false)),
    EuFieldName_(propsDict_.lookupOrDefault<word> ("EuFieldName",(compressible_ == false) ? "Tsource" : "Qsource")),
    tempFieldName_(propsDict_.lookup("tempFieldName")),
    T_(sm.mesh().lookupObject<volScalarField> (tempFieldName_)),
    voidfractionFieldName_(propsDict_.lookupOrDefault<word>("voidfractionFieldName","voidfraction")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    maxSource_(1e30),
    velFieldName_(propsDict_.lookupOrDefault<word>("velFieldName","U")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    lambda_(propsDict_.lookupOrDefault<scalar>("lambda", -1.0)),
    Cp_(propsDict_.lookupOrDefault<scalar>("Cp", -1.0)),
    CpField_
    (
        IOobject
        (
            "CpField",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 2, -2, -1, 0, 0, 0), Cp_)
    ),
    LambdaField_
    (
        IOobject
        (
            "LambdaField",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(1, 1, -3, -1, 0, 0, 0), lambda_)
    ),
    NuCorrelation_(propsDict_.lookupOrDefault<word>("NuCorrelation", "LiMason"))
{
    /*if (lambda_<0)
    {
        Warning << "You tried using an invalid value (<0) for lambda - using SMALL instead!" << endl;
        lambda_=SMALL;
    } else if (lambda_==0)
        Warning << "Ignoring fluid-particle heat transfer as lambda was is to 0!" << endl;*/

    if (propsDict_.found("maxSource"))
    {
        maxSource_=readScalar(propsDict_.lookup ("maxSource"));
        Info << "limiting eulerian source field to: " << maxSource_ << endl;
    }

    if(compressible_) Info << "Compressibility mode enabled." << endl;

    if (NuCorrelation_ == "LiMason")
    {
        Info << "LaEuScalarTemp: using LiMason correlation for Nusselt number." << endl;
        Nusselt=&LaEuScalarTemp::NuLiMason;
    }
    else if (NuCorrelation_ == "Deen")
    {
        Info << "LaEuScalarTemp: using Deen correlation for Nusselt number." << endl;
        Nusselt = &LaEuScalarTemp::NuDeen;
    }
    else if (NuCorrelation_ == "Gunn")
    {
        Info << "LaEuScalarTemp: using Gunn correlation for Nusselt number." << endl;
        Nusselt = &LaEuScalarTemp::NuGunn;
    }
    else
        FatalError << "Unknown name for Nusselt number correlation: '" << NuCorrelation_
            << "'. Must be 'LiMason' or 'Deen' or 'Gunn'."
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

    //set default switches (hard-coded default = false)
    forceSubM(0).setSwitches(32,true); // pullConvectiveHeatFlux
    forceSubM(0).setSwitches(33,true); // pullTemp

    // re-setup required communication
    for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
        forceSubM(iFSub).setupCommunication();

    particleCloud_.checkCG(true);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

LaEuScalarTemp::~LaEuScalarTemp()
{
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * * //

void LaEuScalarTemp::setForce() const
{
        // build reference
        volScalarField& EuField(particleCloud_.mesh().lookupObjectRef<volScalarField> (EuFieldName_));

        // reset Scalar field
        forAll(EuField,cellI) EuField[cellI] = 0.;

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

        #include "resetVoidfractionInterpolator.H"
        #include "resetUInterpolator.H"
        #include "resetTInterpolator.H"

        if(Cp_ < 0) CpField_ = EuField.mesh().lookupObject<volScalarField> ("Cp");
        if(lambda_<0)
        {
            const volScalarField& alpha = EuField.mesh().lookupObject<volScalarField> ("thermo:alpha");
            LambdaField_ = alpha*CpField_;
        }

        for(int index = 0;index < particleCloud_.numberOfParticles(); ++index)
        {
            //if(particleCloud_.regionM().inRegion()[index][0])
            //{
                cellI = particleCloud_.cfdemCloud::cellIDs()[index][0];
                if(cellI >= 0)
                {
                    if(forceSubM(0).interpolation())
                    {
                        position = particleCloud_.cfdemCloud::position(index);
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
                    Us = particleCloud_.cfdemCloud::velocity(index);
                    ds = 2*particleCloud_.radius(index);
                    dparcel = ds;
                    forceSubM(0).scaleDia(ds,index); //caution: this fct will scale ds!
                    numberParticlesInParcel    = dparcel/ds;
                    numberParticlesInParcel   *= numberParticlesInParcel*numberParticlesInParcel;
                    nuf = nufField[cellI];
                    rho = rhoField[cellI];
                    if (LambdaField_[cellI]<0)
                    {
                        Warning << "You tried using an invalid value (<0) for lambda - using SMALL instead!" << endl;
                        LambdaField_[cellI]=SMALL;
                    }
                    if (CpField_[cellI]<0)
                    {
                        Warning << "You tried using an invalid value (<0) for Cp - using SMALL instead!" << endl;
                        CpField_[cellI]=SMALL;
                    }

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
                    scalar h(0.);
                    Nup = 0.;
                    Pr = 0.;
                    if(LambdaField_[cellI]>0)
                    {
                        Pr = max(SMALL,CpField_[cellI]*nuf*rho/LambdaField_[cellI]);

                        Nup = (this->*Nusselt)(Rep, Pr, voidfraction);

                        h = LambdaField_[cellI]*Nup/(ds);
                    }

                    // calc convective heat flux [W]
                    scalar partHeatFlux     = h * As * (Tfluid - particleCloud_.fieldsToDEM[particleCloud_.idTemp()][index][0]);
                    particleCloud_.fieldsToDEM[particleCloud_.idConvectiveHeatFlux()][index][0] = partHeatFlux;

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
                        Pout << "particleCloud_.fieldsToDEM[particleCloud_.idTemp()][index][0] = " << particleCloud_.fieldsToDEM[particleCloud_.idTemp()][index][0] << endl;
                        Pout << "particleCloud_.fieldsToDEM[particleCloud_.idConvectiveHeatFlux()][index][0] = " << particleCloud_.fieldsToDEM[particleCloud_.idConvectiveHeatFlux()][index][0] << endl;
                        Pout << "Tfluid = " << Tfluid << endl  ;
                    }
                }
            //}
        }

        particleCloud_.averagingM().setScalarSum
        (
            EuField,
            particleCloud_.fieldsToDEM[particleCloud_.idConvectiveHeatFlux()],
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
                EuField[cellI] /= -rhoField[cellI]*CpField_[cellI]*EuField.mesh().V()[cellI];
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

        EuField.correctBoundaryConditions();

        if(compressible_)
        {
            Info << "total convective particle-fluid heat flux [W] (Eulerian) = "
                 << gSum(EuField*1*EuField.mesh().V())
                 << endl;
        } else
        {
            Info << "total convective particle-fluid heat flux [W] (Eulerian) = "
                 << gSum(EuField*rhoField*CpField_*EuField.mesh().V())
                 << endl;
        }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
