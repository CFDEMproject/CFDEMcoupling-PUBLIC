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

    execute the liggghts command run $nrSteps

Class
    runLiggghts

SourceFiles
    runLiggghts.C

\*---------------------------------------------------------------------------*/

#ifndef runLiggghts_H
#define runLiggghts_H

#include "liggghtsCommandModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class runLiggghts Declaration
\*---------------------------------------------------------------------------*/

class runLiggghts
:
    public liggghtsCommandModel
{
private:
    dictionary propsDict_;

    word command_;

    bool preNo_;

    int stdInterval_;

public:

    //- Runtime type information
    TypeName("runLiggghts");


    // Constructors

        //- Construct from components
        runLiggghts
        (
            const dictionary& dict,
            cfdemCloud& sm,
            int i
        );

    // Destructor

        ~runLiggghts();


    // Member Functions

    const char* command(int);

    string createCommand(word="",int=0,word="",word="",word="",word="");

    bool runCommand(int);

    word name(){return "runLiggghts";};

    void set(int);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
