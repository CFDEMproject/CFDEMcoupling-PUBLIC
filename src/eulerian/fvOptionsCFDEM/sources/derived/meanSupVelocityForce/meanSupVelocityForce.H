/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright (C) 2011 OpenFOAM Foundation
                                Copyright (C) 2012-     DCS Computing GmbH,Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling.  If not, see <http://www.gnu.org/licenses/>.

Description
    This code is based on meanSupVelocity model of OpenFOAM 3.0.x
    It is extended to account for the presence of a particle phase.

    Calculates and applies the force necessary to maintain the specified mean
    velocity.

    Note: Currently only handles kinematic pressure (incompressible solvers).

    \heading Source usage
    Example usage:
    \verbatim
    meanSupVelocityForceCoeffs
    {
        selectionMode   all;                    // Apply force to all cells
        fieldNames      (U);                    // Name of velocity field
        Ubar            (10.0 0 0);             // Desired mean velocity
        relaxation      0.2;                    // Optional relaxation factor
    }
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef meanSupVelocityForce_H
#define meanSupVelocityForce_H

#include "meanVelocityForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
               Class meanSupVelocityForce Declaration
\*---------------------------------------------------------------------------*/

class meanSupVelocityForce
:
    public meanVelocityForce
{
protected:

    // Protected data

        const bool twoPhase_;

        volScalarField alpha_;

        const bool coupled_;

	    volScalarField voidfraction_;

        const word modelName_;

        const scalar alphaMin_;

    // Protected Member Functions

        virtual scalar magUbarAve(const volVectorField& U) const;

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        meanSupVelocityForce(const meanSupVelocityForce&);

        //- Disallow default bitwise assignment
        void operator=(const meanSupVelocityForce&);

public:

    //- Runtime type information
    TypeName("meanSupVelocityForce");


    // Constructors

        //- Construct from explicit source name and mesh
        meanSupVelocityForce
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Member Functions

        // Evaluate

        virtual void correct(volVectorField& U);

        virtual void addSup
        (
            fvMatrix<vector>& eqn,
            const label fieldI
        );

        // IO

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
