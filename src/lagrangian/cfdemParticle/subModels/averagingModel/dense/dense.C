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

#include "dense.H"
#include "addToRunTimeSelectionTable.H"
#include "voidFractionModel.H"

//#include "mpi.h"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dense, 0);

addToRunTimeSelectionTable
(
    averagingModel,
    dense,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dense::dense
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    averagingModel(dict,sm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dense::~dense()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dense::setScalarAverage
(
    volScalarField& field,
    double**& value,
    double**& weight,
    volScalarField& weightField,
    double**const& mask,
    double**const& weight2,          //allows the specification of a 2nd weight field
    bool      weightWithWeight2      //switch to activate 2nd weight field
) const
{
    label cellI;
    scalar valueScal;
    scalar weightP;

    if(weightWithWeight2) 
        FatalError << "dense::setScalarAverage: attempt to weight with weight2, which is not implemented" << abort(FatalError);

    for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {
        if(!checkParticleType(index)) continue; //skip this particle if not correct type

        for(int subCell=0;subCell<particleCloud_.cellsPerParticle()[index][0];subCell++)
        {
            //Info << "subCell=" << subCell << endl;
            cellI = particleCloud_.cellIDs()[index][subCell];

            if (cellI >= 0)
            {
                valueScal = value[index][0];
                weightP = weight[index][0];

                // first entry in this cell
                if(weightField[cellI] == 0)
                {
                    field[cellI] = valueScal;
                    weightField[cellI] = weightP;
                }
                else
                {
                    field[cellI] = (field[cellI]*weightField[cellI]+valueScal*weightP)/(weightField[cellI]+weightP);
                    weightField[cellI] += weightP;
                }
            }
        }
    }

    // correct cell values to patches
    field.correctBoundaryConditions();
}

void dense::setVectorAverage
(
    volVectorField& field,
    double**& value,
    double**& weight,
    volScalarField& weightField,
    double**const& mask,
    double**const& weight2,                  //allows the specification of a 2nd weight field
    bool      weightWithWeight2   //switch to activate 2nd weight field
) const
{
    label cellI;
    vector valueVec;
    scalar weightP;

    //====================
    // debug
    /*for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {
        Pout << "--index=" << index << endl;
        for(int subCell=0;subCell<particleCloud_.cellsPerParticle()[index][0];subCell++)
        {
            cellI = particleCloud_.cellIDs()[index][subCell];

            if (cellI >= 0)
            {
                Pout << "---subCell=" << subCell << endl;
                Pout << "---weight=" << weight[index][subCell] << endl;
                if(weightWithWeight2)    
                    Pout << "---weight2=" << weight2[index][subCell] << endl;
            }
            else
                Pout << "---subCell not found" << endl;
        }
    }*/
    //====================
    if(weightWithWeight2) //use weight2, e.g., mass-averaged
    for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {
        if(!checkParticleType(index)) continue; //skip this particle if not correct type

        for(int subCell=0;subCell<particleCloud_.cellsPerParticle()[index][0];subCell++)
        {
            cellI = particleCloud_.cellIDs()[index][subCell];

            if (cellI >= 0)
            {
                for(int i=0;i<3;i++) valueVec[i] = value[index][i];
                weightP = weight[index][subCell]*weight2[index][subCell];

                // first entry in this cell
                if(weightField[cellI] == 0)
                {
                    field[cellI] = valueVec;
                    weightField[cellI] = weightP;
                }
                else
                {
                    field[cellI] = (field[cellI]*weightField[cellI]+valueVec*weightP)/(weightField[cellI]+weightP); //Running average calc. style
                    weightField[cellI] += weightP;
                }
            }
        }//forAllSubPoints
    }
    else //standard, i.e., volume-averaged
    for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {
        if(!checkParticleType(index)) continue; //skip this particle if not correct type

        for(int subCell=0;subCell<particleCloud_.cellsPerParticle()[index][0];subCell++)
        {
            cellI = particleCloud_.cellIDs()[index][subCell];

            if (cellI >= 0)
            {
                for(int i=0;i<3;i++) valueVec[i] = value[index][i];
                weightP = weight[index][subCell];

                // first entry in this cell
                if(weightField[cellI] == 0)
                {
                    field[cellI] = valueVec;
                    weightField[cellI] = weightP;
                }
                else
                {
                    field[cellI] = (field[cellI]*weightField[cellI]+valueVec*weightP)/(weightField[cellI]+weightP); //Running average calc. style
                    weightField[cellI] += weightP;
                }
            }
        }//forAllSubPoints
    }

    // correct cell values to patches
    field.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
