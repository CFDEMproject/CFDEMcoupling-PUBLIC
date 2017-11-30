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
        
        bool requiresEx_;           //requires an orientation vector

        bool requiresShape_;        //requires shape information 

        bool requiresQuaternion_;

        bool requiresSuperquadric_;
        
        bool pullPushRotation_;     //will exchange rotation information

        bool implicitDrag_;             //will exchange Ksl 
        
        bool implicitAnisotropicDrag_;  //will exchange anisotropic Ksl 
        
        bool implicitRotation_;         //will exchange rotational Ksl value

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

        virtual void setForcesOnParticle() const {};
        
        virtual void multiphaseTurbulence(volScalarField&, bool) const;

        //tmp<volScalarField> provideScalarField();

        virtual void manipulateScalarField(volScalarField&) const;

        virtual void manipulateScalarField(volScalarField&, volScalarField&, int) const;

        virtual void commToDEM() const;

    // Access
        word modelType(){ return modelType_; };

        inline volVectorField& impParticleForces() const { return impParticleForces_;};

        inline volVectorField& expParticleForces() const { return expParticleForces_;};

        inline double ** impForces() const { return particleCloud_.impForces_;};

        inline double ** expForces() const { return particleCloud_.expForces_;};

        inline double ** DEMForces() const { return particleCloud_.DEMForces_;};
        
        inline double ** DEMTorques() const { return particleCloud_.DEMTorques();};

        inline double ** Cds() const { return particleCloud_.Cds_;};

        inline double ** CdsExtra() const { return particleCloud_.CdsExtra();};

        inline double ** CdsRotation() const { return particleCloud_.CdsRotation();};

        inline double ** omegaFluid() const { return particleCloud_.omegaFluid();};

        inline double ** fluidVel() const { return particleCloud_.fluidVel_;};

        virtual inline bool& requiresEx() { return requiresEx_;};

        virtual inline bool& requiresQuaternion() { return requiresQuaternion_;};

        virtual inline bool& requiresSuperquadric() { return requiresSuperquadric_;};

        virtual inline bool& requiresShape() { return requiresShape_;};

        virtual inline bool& pullPushRotation() {return pullPushRotation_;};

        virtual inline bool& implicitDrag() {return implicitDrag_;};
        
        virtual inline bool& implicitAnisotropicDrag() {return implicitAnisotropicDrag_;};
        
        virtual inline bool& implicitRotation() {return implicitRotation_;};
        
        void repartitionImExForces() const; //Repartition Implixit/Explicit forces

        void treatVoidCells() const;

        inline const wordList& forceSubModels(){ return forceSubModels_; };

        inline const forceSubModel& forceSubM(int i) const { return forceSubModel_[i]; };

        inline int nrForceSubModels() const { return forceSubModels_.size(); };

        void setForceSubModels(dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
