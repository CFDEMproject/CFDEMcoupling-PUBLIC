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

    check the coupling interval,
    by comparing couplingTime/particle relaxation time

Class
    checkCouplingInterval

SourceFiles
    checkCouplingInterval.C

\*---------------------------------------------------------------------------*/

#ifndef checkCouplingInterval_H
#define checkCouplingInterval_H

#include "forceModel.H"
#include "dataExchangeModel.H"
#include "mpi.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class checkCouplingInterval Declaration
\*---------------------------------------------------------------------------*/

class checkCouplingInterval
:
    public forceModel
{
private:
    dictionary propsDict_;

    Switch warnOnly_;

    word velocityFieldName_;

    const volVectorField& U_;

    const scalar rhoP_;

    scalar maxCFL_;

    scalar maxPCFL_;

    scalar maxAccNr_;

    scalar UmaxExpected_;

    scalar minAllowedVcellByVparcel_;

    mutable int nextRun_;

    scalar timeInterval_;

    int couplingStepInterval_;

    inline bool doCheck() const;

public:

    //- Runtime type information
    TypeName("checkCouplingInterval");


    // Constructors

        //- Construct from components
        checkCouplingInterval
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~checkCouplingInterval();


    // Member Functions
        void setForce() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
