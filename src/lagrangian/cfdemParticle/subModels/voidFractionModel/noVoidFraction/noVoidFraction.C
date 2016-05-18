/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
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

#include "noVoidFraction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(noVoidFraction, 0);

addToRunTimeSelectionTable
(
    voidFractionModel,
    noVoidFraction,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
noVoidFraction::noVoidFraction
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    voidFractionModel(dict,sm)
{
    maxCellsPerParticle_ = 29; // cloud uses maxCellsPerParticle_ from this voidfraction model
                               // here we assume voidFractionModelMS is dividedMS
                               // better solution: maxCellsPerParticle_ lives in cloud

    Warning <<"Please use 'noVoidFraction' in combination with a voidFractionModelMS dividedMS.\n" << endl;

    // avoid that this model dumps Prev and Next fields
    voidfractionPrev_.writeOpt() = IOobject::NO_WRITE;
    voidfractionNext_.writeOpt() = IOobject::NO_WRITE;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

noVoidFraction::~noVoidFraction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void noVoidFraction::setvoidFraction(double** const& mask,double**& voidfractions,double**& particleWeights,double**& particleVolumes,double**& particleV) const
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
