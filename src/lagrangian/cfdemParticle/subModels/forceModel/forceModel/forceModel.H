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
    forceModel

SourceFiles
    forceModel.C

\*---------------------------------------------------------------------------*/

#ifndef forceModel_H
#define forceModel_H

#include "fvCFD.H"
#include "cfdemCloud.H"
#include "probeModel.H"
#include "forceSubModel.H"
#include "interpolationCellPointFace.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class forceModel Declaration
\*---------------------------------------------------------------------------*/

class forceModel
{

protected:

    // Protected data
        const dictionary& dict_;

        cfdemCloud& particleCloud_;

        //Switch treatExplicit_;                      // marker to treat force in implicit way (otherwise explicit)
        //Switch treatDEM_;                           // marker to use the force only on DEM side
        //Switch implDEM_;                           // marker to use the implicit force on DEM side

        mutable volVectorField impParticleForces_;  // sum of implicit particle forces [N]

        mutable volVectorField expParticleForces_;  // sum of explicit particle forces [N]

        const word modelType_;

        bool probeIt_;

        bool particleBased_;

        wordList forceSubModels_;

        autoPtr<forceSubModel>* forceSubModel_;

        mutable autoPtr<interpolation<scalar> > GInterpolator_;

        mutable autoPtr<interpolation<scalar> > voidfractionInterpolator_;

        mutable autoPtr<interpolation<vector> > UInterpolator_;

        mutable autoPtr<interpolation<vector> > vorticityInterpolator_;

        mutable autoPtr<interpolation<vector> > gradPInterpolator_;

        mutable autoPtr<interpolation<vector> > gradUInterpolator_;

        mutable autoPtr<interpolation<vector> > gradVoidfractionInterpolator_;

        mutable autoPtr<interpolation<vector> > Up1Interpolator_;

        mutable autoPtr<interpolation<vector> > Up2Interpolator_;

        mutable autoPtr<interpolation<scalar> > dSauterInterpolator_;

        mutable autoPtr<interpolation<scalar> > phiP1Interpolator_;

        mutable autoPtr<interpolation<scalar> > phiP2Interpolator_;

        mutable autoPtr<interpolation<scalar> > alphaInterpolator_;

        mutable autoPtr<interpolation<vector> > gradAlphaInterpolator_;

        mutable autoPtr<interpolation<scalar> > TInterpolator_;

        mutable autoPtr<interpolation<vector> > UsInterpolator_;

        mutable autoPtr<interpolation<scalar> > fluidScalarFieldInterpolator_;

        mutable autoPtr<interpolation<vector> > gradPsolidInterpolator_;

        mutable autoPtr<interpolation<symmTensor> > shearRateInterpolator_;

        mutable autoPtr<interpolation<vector> > DDtUInterpolator_;

        mutable autoPtr<interpolation<vector> > divTauInterpolator_;

        mutable autoPtr<interpolation<scalar> > RhoInterpolator_;

        mutable autoPtr<interpolation<scalar> > kInterpolator_;

        mutable autoPtr<interpolation<scalar> > epsilonInterpolator_;

public:

    //- Runtime type information
    TypeName("forceModel");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            forceModel,
            dictionary,
            (
                const dictionary& dict,
                cfdemCloud& sm
            ),
            (dict,sm)
        );


    // Constructors

        //- Construct from components
        forceModel
        (
            const dictionary& dict,
            cfdemCloud& sm
        );


    // Destructor

        virtual ~forceModel();


    // Selector

        static autoPtr<forceModel> New
        (
            const dictionary& dict,
            cfdemCloud& sm,
            word forceType
        );


    // Member Functions

        void applyDebugSettings(bool) const;

        virtual void setForce() const = 0;

        virtual void multiphaseTurbulence(volScalarField&, bool) const;

        //tmp<volScalarField> provideScalarField();

    // Access
        word modelType(){ return modelType_; };

        bool modelIsParticleBased() {return particleBased_;};

        inline volVectorField& impParticleForces() const { return impParticleForces_;};

        inline volVectorField& expParticleForces() const { return expParticleForces_;};

//         virtual inline double ** impForces() const { return particleCloud_.impForces_;};
//
//         virtual inline double ** expForces() const { return particleCloud_.expForces_;};

        inline double ** omegaFluid() const { return particleCloud_.omegaFluid();};

        virtual inline bool requiresQuaternion() { return forceSubM(0).useQuat();};

        virtual inline bool requiresSuperquadric() { return forceSubM(0).sq();};

        void repartitionImExForces() const; //Repartition Implixit/Explicit forces

        void treatVoidCells() const;

        inline const wordList& forceSubModels(){ return forceSubModels_; };

        inline const forceSubModel& forceSubM(int i) const { return forceSubModel_[i]; };

        inline int nrForceSubModels() const { return forceSubModels_.size(); };

        void setForceSubModels(dictionary& dict);

    // MS member funtions
        void readDHcorr(dictionary& dict);

        void readArea(dictionary& dict);

        virtual void MSinit() {}; // extra init function for MS stuff

    // MS access
        virtual inline double ** impForces() const
        {
            if (forceSubM(0).ms())
                return particleCloud_.impForcesCM();
            else
                return particleCloud_.impForces();
        };

        virtual inline double ** expForces() const
        {
            if (forceSubM(0).ms())
                return particleCloud_.expForcesCM();
            else
                return particleCloud_.expForces();
        };

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
