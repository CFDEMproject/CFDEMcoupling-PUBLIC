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

    -grad(p)*Vp * rho  ... if p is normalized with pressure
    -grad(p)*Vp        ... if p is real pressure
    pressure gradient force
    including interpolation of the velocity to the exact position

Class
    gradPForce

SourceFiles
    gradPForce.C

\*---------------------------------------------------------------------------*/

#ifndef gradPForce_H
#define gradPForce_H

#include "forceModel.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class gradPForce Declaration
\*---------------------------------------------------------------------------*/

class gradPForce
:
    public forceModel
{
private:
    dictionary propsDict_;

    word pFieldName_;

    const volScalarField& p_;

    word velocityFieldName_;

    const volVectorField& U_;

    bool useRho_;

    bool useU_;          // if false: substitution p=0.5*rho*U^2

    mutable double addedMassCoeff_; //added mass coefficient

public:

    //- Runtime type information
    TypeName("gradPForce");


    // Constructors

        //- Construct from components
        gradPForce
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~gradPForce();


    // Member Functions
        void setForce() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
