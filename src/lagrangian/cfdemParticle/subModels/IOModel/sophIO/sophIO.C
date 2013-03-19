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

#include "sophIO.H"
#include "addToRunTimeSelectionTable.H"
#include "forceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(sophIO, 0);

addToRunTimeSelectionTable
(
    IOModel,
    sophIO,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
sophIO::sophIO
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    basicIO(dict,sm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

sophIO::~sophIO()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Public Member Functions

void sophIO::dumpDEMdata() const
{
    if (time_.outputTime())
    {
        basicIO::dumpDEMdata();

        // stream data to file
        streamDataToPath(lagPath_, particleCloud_.voidfractions(), "voidfractions","scalar","scalarField","");
        streamDataToPath(lagPath_, particleCloud_.impForces(), "impForces","vector","vectorField","");
        streamDataToPath(lagPath_, particleCloud_.expForces(), "expForces","vector","vectorField","");
        streamDataToPath(lagPath_, particleCloud_.DEMForces(), "DEMForces","vector","vectorField","");
    }
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Private Member Functions

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
