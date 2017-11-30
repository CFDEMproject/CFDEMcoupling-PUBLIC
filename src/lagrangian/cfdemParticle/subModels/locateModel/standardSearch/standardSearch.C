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
#include "standardSearch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(standardSearch, 0);

addToRunTimeSelectionTable
(
    locateModel,
    standardSearch,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
standardSearch::standardSearch
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    locateModel(dict,sm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

standardSearch::~standardSearch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

label standardSearch::findCell
(
    double** const& mask,
    double**& positions,
    double**& cellIDs,
    int size,
    bool checkRad
) const
{
    vector position;
    for(int index = 0;index < size; ++index)
    {
        //if(mask[index][0] && particleCloud_.radius(index) > SMALL)
        if(!checkRad || particleCloud_.radius(index) > SMALL)
        {
            // create pos vector
            for(int i=0;i<3;i++) position[i] = positions[index][i];

            // find cell
            #if defined(version30)
                cellIDs[index][0] = particleCloud_.mesh().findCell(position, polyMesh::FACE_PLANES);
            #elif defined(version21)
                cellIDs[index][0] = particleCloud_.mesh().findCell(position, polyMesh::FACEPLANES);
            #elif defined(version16ext)
                cellIDs[index][0] = particleCloud_.mesh().findCell(position);
            #endif
        }
        else cellIDs[index][0]=-1;
    }

    return 1;
}

label standardSearch::findSingleCell
(
    vector& position,
    label& oldCellID
) const
{
    // find cell
    #if defined(version30)
        return particleCloud_.mesh().findCell(position, polyMesh::FACE_PLANES);
    #elif defined(version21)
        return particleCloud_.mesh().findCell(position, polyMesh::FACEPLANES);
    #elif defined(version16ext)
        return particleCloud_.mesh().findCell(position);
    #endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
