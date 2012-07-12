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

    execute the liggghts command write_restart

Class
    writeLiggghts

SourceFiles
    writeLiggghts.C

\*---------------------------------------------------------------------------*/

#ifndef writeLiggghts_H
#define writeLiggghts_H

#include "liggghtsCommandModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class writeLiggghts Declaration
\*---------------------------------------------------------------------------*/

class writeLiggghts
:
    public liggghtsCommandModel
{
private:

    dictionary propsDict_;

    word command_;

    fileName path_;

    fileName writeName_;

    Switch writeLast_;

    Switch overwrite_;

public:

    //- Runtime type information
    TypeName("writeLiggghts");


    // Constructors

        //- Construct from components
        writeLiggghts
        (
            const dictionary& dict,
            cfdemCloud& sm,
            int i
        );

    // Destructor

        ~writeLiggghts();


    // Member Functions

    const char* command();

    bool runCommand(int);

    word name(){return "writeLiggghts";};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
