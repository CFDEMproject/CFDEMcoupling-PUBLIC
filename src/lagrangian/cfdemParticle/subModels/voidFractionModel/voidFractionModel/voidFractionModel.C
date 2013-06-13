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
#include "voidFractionModel.H"
#include "dataExchangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(voidFractionModel, 0);

defineRunTimeSelectionTable(voidFractionModel, dictionary);

// * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
voidFractionModel::voidFractionModel
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    dict_(dict),
    particleCloud_(sm),
    voidfractionPrev_
    (   IOobject
        (
            "voidfractionPrev",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,//MUST_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh().lookupObject<volScalarField> ("voidfraction")
        /*sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0), 1)*/
    ),
    voidfractionNext_
    (   IOobject
        (
            "voidfractionNext",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,//MUST_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh().lookupObject<volScalarField> ("voidfraction")
        /*sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0), 1)*/
    ),
    cellsPerParticle_(NULL),
    maxCellsPerParticle_(1),
    weight_(1.)
{
    particleCloud_.dataExchangeM().allocateArray(cellsPerParticle_,1,maxCellsPerParticle_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

voidFractionModel::~voidFractionModel()
{
    particleCloud_.dataExchangeM().destroy(cellsPerParticle_,1);
}

// * * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * * //
tmp<volScalarField> Foam::voidFractionModel::voidFractionInterp() const
{
    tmp<volScalarField> tsource
    (
        new volScalarField
        (
            IOobject
            (
                "alpha_voidFractionModel",
                particleCloud_.mesh().time().timeName(),
                particleCloud_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            particleCloud_.mesh(),
            dimensionedScalar
            (
                "zero",
                dimensionSet(0, 0, 0, 0, 0),
                0
            )
        )
    );

    scalar tsf = particleCloud_.dataExchangeM().timeStepFraction();
    if(1-tsf < 1e-4) //tsf==1
    {
        tsource() = voidfractionPrev_;
    }
    else
    {
        tsource() = (1 - tsf) * voidfractionPrev_ + tsf * voidfractionNext_;
    }
    return tsource;
}

void Foam::voidFractionModel::resetVoidFractions() const
{
    voidfractionPrev_.internalField() = voidfractionNext_.internalField();
    voidfractionNext_.internalField() = 1;
}

/*void Foam::voidFractionModel::undoVoidFractions(double**const& mask) const
{
    voidfractionPrev_.internalField() = voidfractionNext_.internalField();

    for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {
        if(mask[index][0])
        {
            // undo voidfraction cause by particle
            label cellI = particleCloud_.cellIDs()[index][0];
            scalar cellVolume=voidfractionNext_.mesh().V()[cellI];
            voidfractionNext_[cellI] += particleCloud_.particleVolumes()[index][0]/cellVolume;
        }
    }
}*/

double** const& Foam::voidFractionModel::cellsPerParticle() const
{
    return cellsPerParticle_;
}

int Foam::voidFractionModel::maxCellsPerParticle() const
{
    return maxCellsPerParticle_;
}

void Foam::voidFractionModel::reAllocArrays() const
{
    if(particleCloud_.numberOfParticlesChanged())
    {
        // get arrays of new length
        particleCloud_.dataExchangeM().allocateArray(cellsPerParticle_,1,1);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
