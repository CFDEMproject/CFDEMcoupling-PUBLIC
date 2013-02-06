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
    dividedVoidFraction

SourceFiles
    dividedVoidFraction.C

\*---------------------------------------------------------------------------*/

#ifndef dividedVoidFraction_H
#define dividedVoidFraction_H

#include "voidFractionModel.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class noDrag Declaration
\*---------------------------------------------------------------------------*/

class dividedVoidFraction
:
    public voidFractionModel
{

private:
        dictionary propsDict_;

        const scalar alphaMin_;          //NP min value of voidFraction

	    mutable bool alphaLimited_;

        mutable scalar tooMuch_;         //NP particle volume which is lost due to voidFraction limitation

        const scalar scaleUpVol_;        //NP scaling radius, keeping volume of particle

        bool interpolation_;

public:

    //- Runtime type information
    TypeName("divided");


    // Constructors

        //- Construct from components
        dividedVoidFraction
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~dividedVoidFraction();


    // Member Functions
        void setvoidFraction(double** const& ,double**&, double**&, double**&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
