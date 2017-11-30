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
    cfdemCloud& sm,
    word name
)
:
    locateModel(dict,sm),
    propsDict_(dict.subDict(name == "" ? typeName + "Props" : name + "Props")),
    treeSearch_(propsDict_.lookupOrDefault<Switch>("treeSearch", true)),
    #if defined(version30)
        searchEngine_(particleCloud_.mesh(),polyMesh::FACE_PLANES)
    #elif defined(version21)
        searchEngine_(particleCloud_.mesh(),polyMesh::FACEPLANES) // FACEPLANES or FACECENTRETETS; FACEDIAGTETS not stable
    #elif defined(version16ext)
        searchEngine_(particleCloud_.mesh(),false) //(particleCloud_.mesh(),faceDecomp_)
    #endif
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
    int size,
    bool checkRad
) const
{
    vector position;
    for(int index = 0;index < size; ++index)
    {
        if(!checkRad || particleCloud_.radius(index) > SMALL)
        {
            // create pos vector
            for(int i=0;i<3;i++) position[i] = positions[index][i];

            cellIDs[index][0] = searchEngine_.findCell(position,-1,treeSearch_);

            // give it another try with CELL_TETS (very expensive)
            //if(cellIDs[index][0]==-1) cellIDs[index][0] = particleCloud_.mesh().findCell(position,polyMesh::CELL_TETS);
        }
        else cellIDs[index][0] = -1;
    }
    return 1;
}

label engineSearch::findSingleCell
(
    vector& position,
    label& oldCellID
) const
{
    label cellI = searchEngine_.findCell(position,oldCellID,treeSearch_);

    // give it another try with CELL_TETS (very expensive)
    //if(cellI==-1) cellI = particleCloud_.mesh().findCell(position,polyMesh::CELL_TETS);

    return cellI;
}

label engineSearch::intersection
(
    const point& pStart,
    const point& pEnd
) const
{
    // find intersection with boundary
    label face = searchEngine_.findNearestBoundaryFace(pEnd);

    // try alternative
    if (face==-1)
    {
        face = searchEngine_.intersection(pStart,pEnd).index();

        if (face==-1 && mag(pStart-point(0,0,0))<SMALL)
        {
            point pStart2 = pEnd+0.0001*(pStart-pEnd)/mag(pStart-pEnd);
            face = searchEngine_.intersection(pStart2,pEnd).index();
        }
    }
    return face;
}

label engineSearch::intersections
(
    const point& pStart,
    const point& pEnd
) const
{
    return searchEngine_.intersection(pEnd,pStart).index();
}

label engineSearch::findNearestCell
(
    const point& pStart
) const
{
    return searchEngine_.findNearestCell(pStart,-1,true);
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
