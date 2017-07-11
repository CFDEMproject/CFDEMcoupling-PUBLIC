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

    Gaussian distributed particle voidfraction model
    contribution from RQ

Class
    voidFractionModel

SourceFiles
    voidFractionModel.C

\*---------------------------------------------------------------------------*/

#ifndef voidFractionModel_H
#define voidFractionModel_H

#include "fvCFD.H"
#include "cfdemCloud.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class voidFractionModel Declaration
\*---------------------------------------------------------------------------*/

class voidFractionModel
{

protected:

    // Protected data
        const dictionary& dict_;

        cfdemCloud& particleCloud_;

        mutable volScalarField voidfractionPrev_;

        mutable volScalarField voidfractionNext_;

        mutable double ** cellsPerParticle_;

        int maxCellsPerParticle_;

        mutable scalar weight_;

        mutable scalar porosity_;

        bool requiresSuperquadric_;

    // Protected member functions

public:

    //- Runtime type information
    TypeName("voidFractionModel");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            voidFractionModel,
            dictionary,
            (
                const dictionary& dict,
                cfdemCloud& sm
            ),
            (dict,sm)
        );


    // Constructors

        //- Construct from components
        voidFractionModel
        (
            const dictionary& dict,
            cfdemCloud& sm
        );


    // Destructor

        virtual ~voidFractionModel();


    // Selector

        static autoPtr<voidFractionModel> New
        (
            const dictionary& dict,
            cfdemCloud& sm
        );


    // public member functions
        void applyDebugSettings(bool) const;

        virtual void setvoidFraction(double** const&,double**&,double**&,double**&,double**&) const = 0;

        tmp<volScalarField> voidFractionInterp() const;

        inline volScalarField& voidFractionPrev()const { return voidfractionPrev_; }

        inline volScalarField& voidFractionNext()const { return voidfractionNext_; }

        inline scalar weight()const { return weight_; }

        inline scalar porosity()const { return porosity_; }

        inline void checkWeightNporosity(dictionary& propsDict) const
        {
            if (propsDict.found("weight")) weight_ = readScalar(propsDict.lookup("weight"));
            if (propsDict.found("porosity")) porosity_ = readScalar(propsDict.lookup("porosity"));
        };

        void resetVoidFractions() const;

        //void undoVoidFractions(double**const&) const;

        double** const& cellsPerParticle() const;

        int maxCellsPerParticle() const;

        void reAllocArrays() const;

        void reAllocArrays(int nP) const;  //force number of particles during reallocation, for CFD offline-use

        virtual void setParticleType(label type) const {};

        virtual bool checkParticleType(label) const {return true;}; //consider all particles by default

        bool requiresSuperquadric() const;

        virtual double minPeriodicDistance
        (
            int index,
            vector    cellCentrePosition,
            vector    positionCenter,
            boundBox  globalBb,
            vector&   minPeriodicPos,
            vector    dirCheckRange = vector(1,1,1)
        ) const;
        virtual double pointInParticle(int index, vector positionCenter, vector point, double scale) const;
        virtual double pointInParticle(int index, vector positionCenter, vector point) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
