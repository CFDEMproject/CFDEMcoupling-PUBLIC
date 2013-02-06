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

#include "engineSearch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(engineSearch, 0);

addToRunTimeSelectionTable
(
    locateModel,
    engineSearch,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
engineSearch::engineSearch
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    locateModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    //faceDecomp_(propsDict_.lookup("faceDecomp")),
    treeSearch_(propsDict_.lookup("treeSearch")),
    #ifdef version16ext
        searchEngine_(particleCloud_.mesh(),false) //(particleCloud_.mesh(),faceDecomp_)
    #elif defined(version21)
        searchEngine_(particleCloud_.mesh(),polyMesh::FACEPLANES) // FACEPLANES or FACECENTRETETS; FACEDIAGTETS not stable
    #endif
    //searchEngine_(particleCloud_.mesh(),faceDecomp_) // only 2.0.x
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

engineSearch::~engineSearch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

label engineSearch::findCell
(
    double** const& mask,
    double**& positions,
    double**& cellIDs,
    int size
) const
{
    vector position;
    for(int index = 0;index < size; ++index)
    {
        cellIDs[index][0]=-1;

        //if(mask[index][0] && particleCloud_.radius(index) > SMALL)
        if(particleCloud_.radius(index) > SMALL)
        {

            // create pos vector
            for(int i=0;i<3;i++) position[i] = positions[index][i];
            // find cell
            cellIDs[index][0] =searchEngine_.findCell(position,cellIDs[index][0],treeSearch_);
        }
    }
    return 1;
}

label engineSearch::findSingleCell
(
    vector& position,
    label& oldCellID
) const
{
    // find cell
    return searchEngine_.findCell(position,oldCellID,treeSearch_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
