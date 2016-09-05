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

    Gidaspow drag law
    - only valid for low-Reynolds number systems (Re_p<1000)
    - including interpolation of the velocity to the exact position
    - splits off explicit drag component due to fluctuation in fluid and particle
      velocity (optional via forceSubModel "ImExCorr")

Class
    GidaspowDrag

SourceFiles
    GidaspowDrag.C

\*---------------------------------------------------------------------------*/

#ifndef GidaspowDrag_H
#define GidaspowDrag_H

#include "forceModel.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class GidaspowDrag Declaration
\*---------------------------------------------------------------------------*/

class GidaspowDrag
:
    public forceModel
{
private:
    dictionary propsDict_;

    word velFieldName_;

    const volVectorField& U_;

    word voidfractionFieldName_;

    const volScalarField& voidfraction_;

    const scalar phi_;

    word UsFieldName_;   

    const volVectorField& UsField_;  // the average particle velocity field

    mutable scalar switchingVoidfraction_; //voidfraction above which dilute formulation will be used

public:

    //- Runtime type information
    TypeName("GidaspowDrag");


    // Constructors

        //- Construct from components
        GidaspowDrag
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~GidaspowDrag();


    // Member Functions
        void setForce() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
