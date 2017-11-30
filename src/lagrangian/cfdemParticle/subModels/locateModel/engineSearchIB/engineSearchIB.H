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
    engineSearchIB

SourceFiles
    engineSearchIB.C

\*---------------------------------------------------------------------------*/

#ifndef engineSearchIB_H
#define engineSearchIB_H

#include "engineSearch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class noDrag Declaration
\*---------------------------------------------------------------------------*/

class engineSearchIB
:
    public engineSearch
{

protected:
    const label zSplit_;

    const label xySplit_;

    scalar coef_;

    mutable autoPtr<boundBox> bbPtr_;

    bool verbose_;

    const label numberOfSatellitePoints_;

    std::vector<vector> satellitePoints_;

public:

    //- Runtime type information
    TypeName("engineIB");


    // Constructors

    //- Construct from components
    engineSearchIB
    (
        const dictionary& dict,
        cfdemCloud& sm
    );

    // Destructor

    ~engineSearchIB();


    // Member Functions
    label findCell
    (
        double** const& mask,
        double**& positions,
        double**& cellIDs,
        int size,
        bool checkRad=true
    ) const;

    vector generateSatellitePoint
    (
        int countPoints
    ) const;
    virtual vector getSatellitePoint
    (
        int index,
        int countPoints
    ) const;
    bool isInsideRectangularDomain
    (
        vector point,
        scalar skin
    ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
