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
    and OpenFOAM. Note: this code is not part of OpenFOAM (see DISCLAIMER).
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "explicitCouple.H"
#include "addToRunTimeSelectionTable.H"
#include "forceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(explicitCouple, 0);

addToRunTimeSelectionTable
(
    momCoupleModel,
    explicitCouple,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
explicitCouple::explicitCouple
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    momCoupleModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    fPrev_
    (   IOobject
        (
            "fPrev",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,//MUST_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh().lookupObject<volVectorField> ("f")
        //sm.mesh(),
        //dimensionedVector("zero", dimensionSet(1,-2,-2,0,0), vector(0,0,0)) // N/m3
    ),
    fNext_
    (   IOobject
        (
            "fNext",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,//MUST_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh().lookupObject<volVectorField> ("f")
        //sm.mesh(),
        //dimensionedVector("zero", dimensionSet(1,-2,-2,0,0), vector(0,0,0)) // N/m3
    ),
    fLimit_(1e10,1e10,1e10)
{
    if (propsDict_.found("fLimit"))
    {
        fLimit_=vector(propsDict_.lookup ("fLimit"));
        Info << "explicit momentum exchange field is limited to : " << fLimit_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

explicitCouple::~explicitCouple()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
tmp<volVectorField> explicitCouple::expMomSource() const
{
    tmp<volVectorField> tsource
    (
        new volVectorField
        (
            IOobject
            (
                "f_explicitCouple",
                particleCloud_.mesh().time().timeName(),
                particleCloud_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            particleCloud_.mesh(),
            dimensionedVector
            (
                "zero",
                dimensionSet(1, -2, -2, 0, 0), // N/m3
                vector::zero
            )
        )
    );

    // calc fNext
    forAll(fNext_,cellI)
    {
        fNext_[cellI] = particleCloud_.forceM(0).expParticleForces()[cellI] / particleCloud_.mesh().V()[cellI];

        // limiter
        for (int i=0;i<3;i++)
        {
            if (fNext_[cellI][i] > fLimit_[i]) fNext_[cellI][i] = fLimit_[i];
        }
    }

    // underrelaxation of f
    if (particleCloud_.dataExchangeM().couplingStep() > 1)
    {
        tsource() = (1 - particleCloud_.dataExchangeM().timeStepFraction()) * fPrev_
                    + particleCloud_.dataExchangeM().timeStepFraction() * fNext_;
    }
    else
    {
        tsource() = fNext_;
    }

    return tsource;
}

void Foam::explicitCouple::resetMomSourceField() const
{
    fPrev_.internalField() = fNext_.internalField();
    fNext_.internalField() = vector::zero;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
