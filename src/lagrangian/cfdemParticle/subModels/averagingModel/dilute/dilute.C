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

#include "dilute.H"
#include "addToRunTimeSelectionTable.H"
#include "voidFractionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dilute, 0);

addToRunTimeSelectionTable
(
    averagingModel,
    dilute,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dilute::dilute
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    averagingModel(dict,sm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dilute::~dilute()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dilute::setScalarAverage
(
    volScalarField& field,
    double**& value,
    double**& weight,
    volScalarField& weightField,
    double**const& mask,
    double**const& weight2,       //allows the specification of a 2nd weight field
    bool      weightWithWeight2   //switch to activate 2nd weight field
) const
{
    label cellI;
    scalar valueScal;
    scalar weightP;

    if(weightWithWeight2) 
        FatalError << "dilute::setScalarAverage: attempt to weight with weight2, which is not implemented" << abort(FatalError);

    for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {
        for(int subCell=0;subCell<particleCloud_.cellsPerParticle()[index][0];subCell++)
        {
            //Info << "subCell=" << subCell << endl;
            cellI = particleCloud_.cfdemCloud::cellIDs()[index][subCell];

            if (cellI >= 0)
            {
                valueScal = value[index][0];
                weightP = weight[index][0];
                weightField[cellI] += weightP;

                field[cellI] = valueScal/weightP;
            }
        }
    }

    // correct cell values to patches
    field.correctBoundaryConditions();
}

void dilute::setVectorAverage
(
    volVectorField& field,
    double**& value,
    double**& weight,
    volScalarField& weightField,
    double**const& mask,
    double**const& weight2,         //allows the specification of a 2nd weight field
    bool weightWithWeight2    //switch to activate 2nd weight field
) const
{
    label cellI;
    vector valueVec;
    scalar weightP;

    if(weightWithWeight2) //use weight2, e.g., mass-averaged - has no effect, just weight is DIFFERENT!
    {
        for(int index=0; index< particleCloud_.numberOfParticles(); index++)
        {
            for(int subCell=0;subCell<particleCloud_.cellsPerParticle()[index][0];subCell++)
            {
                cellI = particleCloud_.cfdemCloud::cellIDs()[index][subCell];
                if (cellI >= 0)
                {
                    for(int i=0;i<3;i++) valueVec[i] = value[index][i];
                    weightP = weight[index][subCell]*weight2[index][subCell];
                    weightField[cellI] += weightP;
                    if(weightP > 0) field[cellI] = valueVec; //field[cellI] = valueVec/weightP;
                }
            }
        }
    }
    else //standard, i.e., volume-averaged - has no effect, just weight is DIFFERENT!
    {
        for(int index=0; index< particleCloud_.numberOfParticles(); index++)
        {
            for(int subCell=0;subCell<particleCloud_.cellsPerParticle()[index][0];subCell++)
            {
                //Info << "subCell=" << subCell << endl;
                cellI = particleCloud_.cfdemCloud::cellIDs()[index][subCell];

                if (cellI >= 0)
                {
                    for(int i=0;i<3;i++) valueVec[i] = value[index][i];
                    weightP = weight[index][subCell];
                    weightField[cellI] += weightP;
                    if(weightP > 0) field[cellI] = valueVec; //field[cellI] = valueVec/weightP;
                    //else Warning << "!!! W A R N I N G --- weightP <= 0" << endl;
                }
            }
        }
    }

    // correct cell values to patches
    field.correctBoundaryConditions();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
