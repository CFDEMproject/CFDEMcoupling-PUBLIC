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
    explicitCouple

SourceFiles
    explicitCouple.C

\*---------------------------------------------------------------------------*/

#ifndef explicitCouple_H
#define explicitCouple_H

#include "momCoupleModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class noDrag Declaration
\*---------------------------------------------------------------------------*/

class explicitCouple
:
    public momCoupleModel
{
private:

    dictionary propsDict_;

    mutable volVectorField fPrev_;

    mutable volVectorField fNext_;

    vector fLimit_;     // limit for for exchange field

    virtual inline vector arrayToField(label) const;

    mutable volVectorField sourceField_; // N/m3

public:

    //- Runtime type information
    TypeName("explicitCouple");


    // Constructors

        //- Construct from components
        explicitCouple
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~explicitCouple();


    // Member Functions

        void applyDebugSettings(bool) const;

        tmp<volVectorField> expMomSource() const;

        void resetMomSourceField() const;

        void setSourceField(volVectorField &) const;

        dimensionedVector returnIntegralSourceField() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
