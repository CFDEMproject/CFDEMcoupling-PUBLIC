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

Class
    interface

SourceFiles
    interface.C

\*---------------------------------------------------------------------------*/

#ifndef interface_H
#define interface_H

#include "forceModel.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class interface Declaration
\*---------------------------------------------------------------------------*/

class interface
:
    public forceModel
{
private:
    dictionary propsDict_;

    word VOFvoidfractionFieldName_;

    const volScalarField& alpha_;    // ref to VOF voidfraction

    word gradAlphaName_;

    const volVectorField& gradAlpha_;    // ref to fvc::grad(VOF voidfraction)

    const scalar sigma_;

    const scalar theta_;

    const scalar alphaThreshold_;

    const scalar  deltaAlphaIn_;

    const scalar  deltaAlphaOut_;

    scalar C_; // model constant

    bool interpolation_; // use interpolated field values (does not work for fine grids???)

public:

    //- Runtime type information
    TypeName("interface");


    // Constructors

        //- Construct from components
        interface
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~interface();


    // Member Functions
        void setForce() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
