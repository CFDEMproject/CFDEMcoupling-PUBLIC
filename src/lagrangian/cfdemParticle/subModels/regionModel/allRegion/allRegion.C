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

#include "allRegion.H"
#include "addToRunTimeSelectionTable.H"
#include "momCoupleModel.H"
#include "voidFractionModel.H"
#include "averagingModel.H"
#include "forceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(allRegion, 0);

addToRunTimeSelectionTable
(
    regionModel,
    allRegion,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
allRegion::allRegion
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    regionModel(dict,sm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

allRegion::~allRegion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void allRegion::defineRegion() const
{
    reAllocArrays();
    // do nothing
}

void allRegion::expandRegion(volVectorField& U) const
{
    // do nothing
}

void allRegion::resetVolFields(volVectorField& Us) const
{
    // reset voidFraction, particle velocity, weightField, partForcesNext
    particleCloud_.averagingM().resetVectorAverage(particleCloud_.averagingM().UsPrev()
                                                  ,particleCloud_.averagingM().UsNext());
    particleCloud_.voidFractionM().resetVoidFractions();
    particleCloud_.averagingM().resetVectorAverage(particleCloud_.forceM(0).impParticleForces()
                                                  ,particleCloud_.forceM(0).impParticleForces()
                                                  ,true);
    particleCloud_.averagingM().resetVectorAverage(particleCloud_.forceM(0).expParticleForces()
                                                  ,particleCloud_.forceM(0).expParticleForces()
                                                  ,true);
    particleCloud_.averagingM().resetWeightFields();
    particleCloud_.momCoupleM(0).resetMomSourceField();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
