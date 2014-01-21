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

#include "implicitCouple.H"
#include "addToRunTimeSelectionTable.H"
#include "forceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(implicitCouple, 0);

addToRunTimeSelectionTable
(
    momCoupleModel,
    implicitCouple,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
implicitCouple::implicitCouple
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    momCoupleModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    velFieldName_(propsDict_.lookup("velFieldName")),
    granVelFieldName_(propsDict_.lookup("granVelFieldName")),
    voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    Us_(sm.mesh().lookupObject<volVectorField> (granVelFieldName_)),
    alpha_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    KslLimit_(1e10),
    KslPrev_
    (   IOobject
        (
            "KslPrev",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,//MUST_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh().lookupObject<volScalarField> ("Ksl")
        //sm.mesh(),
        //dimensionedScalar("zero", dimensionSet(1,-3,-1,0,0), 0)  // N/m3 / m/s
    ),
    KslNext_
    (   IOobject
        (
            "KslNext",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,//MUST_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh().lookupObject<volScalarField> ("Ksl")
        //sm.mesh(),
        //dimensionedScalar("zero", dimensionSet(1,-3,-1,0,0), 0)  // N/m3 / m/s
    )
{
    Info << "" << endl;

    if (propsDict_.found("KslLimit"))
    {
        KslLimit_=readScalar(propsDict_.lookup ("KslLimit"));
        Info << "implicit momentum exchange field is limited to : " << KslLimit_ << endl;
    }

    if (propsDict_.found("minAlphaP"))
        maxAlpha_ = 1-readScalar(propsDict_.lookup ("minAlphaP"));

    Info << "implicit momentum exchange field calculate if alphaP larger than : " <<  maxAlpha_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

implicitCouple::~implicitCouple()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> implicitCouple::impMomSource() const
{
    tmp<volScalarField> tsource
    (
        new volScalarField
        (
            IOobject
            (
                "Ksl_implicitCouple",
                particleCloud_.mesh().time().timeName(),
                particleCloud_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            particleCloud_.mesh(),
            dimensionedScalar
            (
                "zero",
                dimensionSet(1, -3, -1, 0, 0), // N/m3 / m/s
                0
            )
        )
    );

    scalar tsf = particleCloud_.dataExchangeM().timeStepFraction();

    // calc Ksl
    scalar Ur;

    if(1-tsf < 1e-4) //tsf==1
    {
        forAll(KslNext_,cellI)
        {
            Ur = mag(U_[cellI] - Us_[cellI]);

            if(Ur > SMALL && alpha_[cellI] < maxAlpha_) //momentum exchange switched off if alpha too big
            {
                KslNext_[cellI] = mag(particleCloud_.forceM(0).impParticleForces()[cellI])
                            / Ur
                            / particleCloud_.mesh().V()[cellI];
            }
            else KslNext_[cellI] = 0;

            // limiter
            if (KslNext_[cellI] > KslLimit_) KslNext_[cellI] = KslLimit_;
        }
        tsource() = KslPrev_;
    }else
    {
        tsource() = (1 - tsf) * KslPrev_ + tsf * KslNext_;
    }

    return tsource;
}

void Foam::implicitCouple::resetMomSourceField() const
{
    KslPrev_.internalField() = KslNext_.internalField();
    KslNext_.internalField() = 0;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
