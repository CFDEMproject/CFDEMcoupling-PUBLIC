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

#include "checkCouplingInterval.H"
#include "addToRunTimeSelectionTable.H"

//#include "mpi.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(checkCouplingInterval, 0);

addToRunTimeSelectionTable
(
    forceModel,
    checkCouplingInterval,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
checkCouplingInterval::checkCouplingInterval
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    densityFieldName_(propsDict_.lookup("densityFieldName")),
    rho_(sm.mesh().lookupObject<volScalarField> (densityFieldName_)),
    rhoP_(readScalar(propsDict_.lookup("rhoP")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

checkCouplingInterval::~checkCouplingInterval()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void checkCouplingInterval::setForce() const
{
    if(particleCloud_.mesh().time().write())
    {
        // get viscosity field
        #ifdef comp
            const volScalarField nufField = particleCloud_.turbulence().mu() / rho_;
        #else
            const volScalarField& nufField = particleCloud_.turbulence().nu();
        #endif

        // find min particle relaxation time
        scalar minTauP = 1000;
        scalar tauP = -1;
        scalar scaledRad = -1.;
        for(int index = 0;index <  particleCloud_.numberOfParticles(); index++)
        {
            label cellI = particleCloud_.cellIDs()[index][0];

            if (cellI > -1) // particle Found
            {
                scaledRad = particleCloud_.radius(index)/particleCloud_.cg();
                tauP = rhoP_*4*scaledRad*scaledRad/
                        (18 * nufField[cellI] * rho_[cellI]);
                minTauP = min(minTauP,tauP);
            }
        }

        // calc max couplingTime/particle relaxation time ratio
        scalar DEMtime = particleCloud_.dataExchangeM().DEMts()
                        *particleCloud_.dataExchangeM().couplingInterval();
        double accNr = DEMtime/minTauP;
        double accNrAll=-1.;

        MPI_Allreduce(&accNr, &accNrAll, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        Info << "max. occurring acceleration nr: " << accNrAll << endl;
        if(accNrAll > 0.1) Warning << "you should use a smaller coupling interval!" << endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
