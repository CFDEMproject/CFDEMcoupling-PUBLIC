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

Class
    engineSearch

SourceFiles
    engineSearch.C

\*---------------------------------------------------------------------------*/

#ifndef engineSearch_H
#define engineSearch_H

#include "locateModel.H"

# include "meshSearch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class noDrag Declaration
\*---------------------------------------------------------------------------*/

class engineSearch
:
    public locateModel
{

protected:

    dictionary propsDict_;

    Switch treeSearch_;

    mutable meshSearch searchEngine_;

public:

    //- Runtime type information
    TypeName("engine");


    // Constructors

        //- Construct from components
        engineSearch
        (
            const dictionary& dict,
            cfdemCloud& sm,
            word name = word("")
        );

    // Destructor

        ~engineSearch();


    // Member Functions
        virtual label findCell
        (
            double** const& mask,
            double**& positions,
            double**& cellIDs,
            int size,
            bool checkRad=true
        ) const;

        label findSingleCell
        (
            vector& position,
            label& oldCellID
        ) const;

        label intersection
        (
                const point& pStart,
                const point& pEnd
        ) const;

        label intersections
        (
                const point& pStart,
                const point& pEnd
        ) const;

        label findNearestCell
        (
                const point& pStart
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
