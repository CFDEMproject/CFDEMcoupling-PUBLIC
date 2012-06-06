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
    averagingModel

SourceFiles
    averagingModel.C

\*---------------------------------------------------------------------------*/

#ifndef averagingModel_H
#define averagingModel_H

#include "fvCFD.H"
#include "cfdemCloud.H"
#include "dataExchangeModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class averagingModel Declaration
\*---------------------------------------------------------------------------*/

class averagingModel
{

protected:

    // Protected data
        const dictionary& dict_;

        cfdemCloud& particleCloud_;

        mutable volScalarField UsWeightField_;

        mutable volVectorField UsPrev_;

        mutable volVectorField UsNext_;

    // Protected member functions

public:

    //- Runtime type information
    TypeName("averagingModel");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            averagingModel,
            dictionary,
            (
                const dictionary& dict,
                cfdemCloud& sm
            ),
            (dict,sm)
        );


    // Constructors

        //- Construct from components
        averagingModel
        (
            const dictionary& dict,
            cfdemCloud& sm
        );


    // Destructor

        virtual ~averagingModel();


    // Selector

        static autoPtr<averagingModel> New
        (
            const dictionary& dict,
            cfdemCloud& sm
        );


    // Member Function
        virtual void setScalarAverage
        (
            volScalarField& field,
            double**& value,
            double**& weight,
            volScalarField& weightField,
            double**const& mask
        ) const = 0;

        virtual void setVectorAverage
        (
            volVectorField& field,
            double**& value,
            double**& weight,
            volScalarField& weightField,
            double**const& mask
        ) const = 0;

        void undoVectorAverage
        (
            volVectorField& fieldPrev,
            volVectorField& fieldNext,
            volScalarField& weightField,
            double** const& value,
            double** const& weight,
            double**const& mask,
            bool single=false
        ) const;

        void undoVectorSum
        (
            volVectorField& field,
            double** const& value,
            double** const& weight,
            double**const& mask
        ) const;


        void setVectorSum
        (
            volVectorField& field,
            double**& value,
            double**& weight,
            double**const& mask
        ) const;

        void setScalarSum
        (
            volScalarField& field,
            double**& value,
            double**const& weight,
            double**const& mask
        ) const;

        void resetVectorAverage(volVectorField& prev,volVectorField& next,bool single=false) const;

        void resetWeightFields() const;

        void undoWeightFields(double**const&) const;

        tmp<volVectorField> UsInterp() const;

    // Access
        inline volVectorField& UsPrev() const {return UsPrev_;};

        inline volVectorField& UsNext() const {return UsNext_;};

        inline volScalarField& UsWeightField() const {return UsWeightField_;};

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
