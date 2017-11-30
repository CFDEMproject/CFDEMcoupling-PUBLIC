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

#include "trackIO.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(trackIO, 0);

addToRunTimeSelectionTable
(
    IOModel,
    trackIO,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
trackIO::trackIO
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    sophIO(dict,sm),
    partID_(NULL)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

trackIO::~trackIO()
{
        particleCloud_.dataExchangeM().destroy(partID_,-1);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Public Member Functions

int trackIO::dumpDEMdata() const
{
    int npProcs(-1);

    if (dumpNow())
    {
        npProcs=sophIO::dumpDEMdata();


        // get id data from liggghts
        particleCloud_.dataExchangeM().allocateArray(partID_,0.,1);
        particleCloud_.dataExchangeM().getData("id","scalar-atom",partID_);

        // stream data to file
        streamDataToPath(lagPath_, partID_,npProcs,"id","scalar","scalarField","");
    }
    return npProcs;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Private Member Functions

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
