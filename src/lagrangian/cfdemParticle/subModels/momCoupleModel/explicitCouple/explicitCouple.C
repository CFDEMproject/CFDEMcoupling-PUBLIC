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
    fLimit_(1e10,1e10,1e10),
    sourceField_
    (   IOobject
        (
            "sourceField",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        sm.mesh(),
        dimensionedVector("zero", dimensionSet(1,-2,-2,0,0), vector(0,0,0)) // N/m3
    )
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

void Foam::explicitCouple::applyDebugSettings(bool debug) const
{
    if(!debug)
    {
        fPrev_.writeOpt() = IOobject::NO_WRITE;
        fNext_.writeOpt() = IOobject::NO_WRITE;
    }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dimensionedVector explicitCouple::returnIntegralSourceField() const
{
    //Sum is over cells of PROCESOR ONLY!!
    dimensionedVector intSource = dimensionedVector("0", dimensionSet(1, 1, -2, 0, 0), vector::zero);
    forAll(sourceField_,cellI)
        intSource.value() += sourceField_.internalField()[cellI] * particleCloud_.mesh().V()[cellI];

    return intSource;  // in Newton
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
tmp<volVectorField> explicitCouple::expMomSource() const
{
    //This INCLUDES the source field that is set separately!
    scalar tsf = particleCloud_.dataExchangeM().timeStepFraction();

    // update fNext at last subTS
    //if(1-tsf < 1e-4) //tsf==1

    // update KslNext in first subTS
    // NOTE: without following if we could update f every subTS (based on current values) and use this value
    if(tsf < particleCloud_.mesh().time().deltaT().value()/particleCloud_.dataExchangeM().couplingTime() + 0.000001 )
    {
        // calc fNext
        forAll(fNext_,cellI)
        {
            fNext_[cellI] = arrayToField(cellI);
    
            // limiter
            for (int i=0;i<3;i++)
            {
                if (mag(fNext_[cellI][i]) > fLimit_[i]) fNext_[cellI][i] = fLimit_[i];
            }
        }
    }

    /*if(1-tsf < 1e-4) //tsf==1
    {
        return tmp<volVectorField>
        (
            new volVectorField("f_explicitCouple", fPrev_)
        );
    }else*/
    {
        return tmp<volVectorField>
        (
            new volVectorField("f_explicitCouple", (1 - tsf) * fPrev_ + tsf * fNext_)
        );
    }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void Foam::explicitCouple::resetMomSourceField() const
{
    fPrev_ == fNext_;
    fNext_ == dimensionedVector("zero", fNext_.dimensions(), vector::zero);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
inline vector Foam::explicitCouple::arrayToField(label cellI) const
{
    //This INCLUDES the source field that is set separately!
    return particleCloud_.forceM(0).expParticleForces()[cellI] / particleCloud_.mesh().V()[cellI] 
         + sourceField_[cellI];
}

void Foam::explicitCouple::setSourceField(volVectorField & field) const
{
    sourceField_ = field;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
