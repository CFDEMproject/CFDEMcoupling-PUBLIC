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
    smoothingModel

SourceFiles
    smoothingModel.C

\*---------------------------------------------------------------------------*/

#ifndef smoothingModel_H
#define smoothingModel_H

#include "fvCFD.H"
#include "cfdemCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class smoothingModel Declaration
\*---------------------------------------------------------------------------*/

class smoothingModel
{
protected:

    // Protected data
    const dictionary& dict_;
    cfdemCloud& particleCloud_;

public:

    //- Runtime type information
    TypeName("smoothingModel");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            smoothingModel,
            dictionary,
            (
                const dictionary& dict,
                cfdemCloud& sm
            ),
            (dict,sm)
        );

    // Constructors

        //- Construct from components
        smoothingModel
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        virtual ~smoothingModel();

    // Selector

        static autoPtr<smoothingModel> New
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Member Functions

        virtual bool doSmoothing() const;

        // the particle diameter for smoothing
        virtual void dSmoothing(volScalarField&) const;

        virtual void smoothen(volScalarField&) const=0;

        virtual void smoothen(volVectorField&) const=0;

        virtual void smoothenReferenceField(volVectorField&) const=0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
