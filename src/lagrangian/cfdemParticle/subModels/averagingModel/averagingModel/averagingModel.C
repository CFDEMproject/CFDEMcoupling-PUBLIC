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
#include "averagingModel.H"
#include "voidFractionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(averagingModel, 0);

defineRunTimeSelectionTable(averagingModel, dictionary);

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void averagingModel::undoVectorAverage
(
    volVectorField& fieldPrev,
    volVectorField& fieldNext,
    volScalarField& weightField,
    double** const& value,
    double** const& weight,
    double**const& mask,
    bool single
) const
{
// WARNING - not sure if this is valid for dilute model!!!

    if(!single) fieldPrev.internalField() = fieldNext.internalField();

    label cellI;
    vector valueVec;
    scalar weightP;

    for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {
        //if(mask[index][0])
        //{
            for(int subCell=0;subCell<particleCloud_.voidFractionM().cellsPerParticle()[index][0];subCell++)
            {
                //Info << "subCell=" << subCell << endl;
                cellI = particleCloud_.cellIDs()[index][subCell];

                if (cellI >= 0)
                {
                    for(int i=0;i<3;i++) valueVec[i] = value[index][i];
                    weightP = weight[index][subCell];

                    if(weightField[cellI] == weightP)
                    {
                        fieldNext[cellI] = vector(0,0,0);
                    }else
                    {
                        fieldNext[cellI] = (fieldNext[cellI]*weightField[cellI]-valueVec*weightP)/(weightField[cellI]-weightP);
                    }
                }
            }
        //}
    }

    // correct cell values to patches
    fieldNext.correctBoundaryConditions();
}

void averagingModel::undoVectorSum
(
    volVectorField& field,
    double** const& value,
    double** const& weight,
    double**const& mask
) const
{
    label cellI;
    vector valueVec;
    scalar weightP;

    for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {
        //if(mask[index][0])
        //{
            for(int subCell=0;subCell<particleCloud_.voidFractionM().cellsPerParticle()[index][0];subCell++)
            {
                //Info << "subCell=" << subCell << endl;
                cellI = particleCloud_.cellIDs()[index][subCell];

                if (cellI >= 0)
                {
                    for(int i=0;i<3;i++) valueVec[i] = value[index][i];
                    weightP = weight[index][subCell];

                    field[cellI] -= valueVec*weightP;
                }
            }//forAllSubPoints
        //}
    }

    // correct cell values to patches
    field.correctBoundaryConditions();
}

void averagingModel::setVectorSum
(
    volVectorField& field,
    double**& value,
    double**& weight,
    double**const& mask
) const
{
    label cellI;
    vector valueVec;
    scalar weightP;

    for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {
        //if(mask[index][0])
        //{
            for(int subCell=0;subCell<particleCloud_.voidFractionM().cellsPerParticle()[index][0];subCell++)
            {
                cellI = particleCloud_.cellIDs()[index][subCell];

                if (cellI >= 0)
                {
                    for(int i=0;i<3;i++) valueVec[i] = value[index][i];
                    weightP = weight[index][subCell];
                    field[cellI] += valueVec*weightP;
                }
            }//forAllSubPoints
        //}
    }

    // correct cell values to patches
    field.correctBoundaryConditions();
}

void averagingModel::setScalarSum
(
    volScalarField& field,
    double**& value,
    double**const& weight,
    double**const& mask
) const
{
    label cellI;
    scalar valueScal;
    scalar weightP;

    for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {
        //if(mask[index][0])
        //{
            for(int subCell=0;subCell<particleCloud_.voidFractionM().cellsPerParticle()[index][0];subCell++)
            {
                //Info << "subCell=" << subCell << endl;
                cellI = particleCloud_.cellIDs()[index][subCell];

                if (cellI >= 0)
                {
                    valueScal = value[index][0];
                    weightP = weight[index][subCell];
                    field[cellI] += valueScal*weightP;
                }
            }//forAllSubPoints
        //}
    }

    // correct cell values to patches
    field.correctBoundaryConditions();
}

void averagingModel::resetVectorAverage(volVectorField& prev,volVectorField& next,bool single) const
{
    if(!single) prev.internalField() = next.internalField();
    next.internalField() = vector::zero;
}

void averagingModel::resetWeightFields() const
{
    UsWeightField_.internalField() = 0;
}


void Foam::averagingModel::undoWeightFields(double**const& mask) const
{
    for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {
        //if(mask[index][0])
        //{
            // undo voidfraction cause by particle
            label cellI = particleCloud_.cellIDs()[index][0];
            UsWeightField_[cellI] -= particleCloud_.particleWeights()[index][0];
        //}
    }
}

tmp<volVectorField> Foam::averagingModel::UsInterp() const
{
   tmp<volVectorField> tsource
    (
        new volVectorField
        (
            IOobject
            (
                "Us_averagingModel",
                particleCloud_.mesh().time().timeName(),
                particleCloud_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            particleCloud_.mesh(),
            dimensionedVector
            (
                "zero",
                dimensionSet(0, 1, -1, 0, 0),
                vector::zero
            )
        )
    );

    if (particleCloud_.dataExchangeM().couplingStep() > 1)
    {
        tsource() = (1 - particleCloud_.dataExchangeM().timeStepFraction()) * UsPrev_
                    + particleCloud_.dataExchangeM().timeStepFraction() * UsNext_;
    }
    else
    {
        tsource() = UsNext_;
    }

    return tsource;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
averagingModel::averagingModel
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    dict_(dict),
    particleCloud_(sm),
    UsWeightField_
    (
        IOobject
        (
            "UsWeightField_",
            particleCloud_.mesh().time().timeName(),
            particleCloud_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        particleCloud_.mesh(),
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0), 0.0)
    ),
    UsPrev_
    (   IOobject
        (
            "UsPrev",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,//MUST_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh().lookupObject<volVectorField> ("Us")
        /*sm.mesh(),
        dimensionedVector("zero", dimensionSet(0,1,-1,0,0),vector::zero)*/
    ),
    UsNext_
    (   IOobject
        (
            "UsNext",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,//MUST_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh().lookupObject<volVectorField> ("Us")
        /*sm.mesh(),
        dimensionedVector("zero", dimensionSet(0,1,-1,0,0),vector::zero)*/
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

averagingModel::~averagingModel()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
