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
    setDEMGravity

SourceFiles
    setDEMGravity.C

\*---------------------------------------------------------------------------*/

#ifndef setDEMGravity_H
#define setDEMGravity_H

#include "liggghtsCommandModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class setDEMGravity Declaration
\*---------------------------------------------------------------------------*/

class setDEMGravity
:
    public liggghtsCommandModel
{
private:

    dictionary propsDict_;

    int nrModel_;

    int insertionNr_;

    wordList command_;

    wordList filePathList_;

    scalarList scalarList_;

    labelList labelList_;

    #ifdef version21
        const uniformDimensionedVectorField& g_;    // ref to gravity
    #elif defined(version16ext) || defined(version15)
        const dimensionedVector& g_;    // ref to gravity
    #endif

    Switch unfix_;

public:

    //- Runtime type information
    TypeName("setDEMGravity");


    // Constructors

        //- Construct from components
        setDEMGravity
        (
            const dictionary& dict,
            cfdemCloud& sm,
            int i
        );

    // Destructor

        ~setDEMGravity();


    // Member Functions

    word name(){return typeName;};

    const char* command(int);

    bool runCommand(int);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
