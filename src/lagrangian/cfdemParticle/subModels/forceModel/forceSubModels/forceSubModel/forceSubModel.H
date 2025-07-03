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
    forceSubModel

SourceFiles
    forceSubModel.C

\*---------------------------------------------------------------------------*/

#ifndef forceSubModel_H
#define forceSubModel_H

#include "fvCFD.H"
#include "cfdemCloud.H"
#include "probeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class forceSubModel Declaration
\*---------------------------------------------------------------------------*/

class forceSubModel
{
private:


protected:

    // Protected data
        const dictionary& dict_;

        cfdemCloud& particleCloud_;

        forceModel& forceModel_;

        label nrDefaultSwitches_;                   // nr of switches defined in mother class

        wordList switchesNameList_;                 // names of switches available

        mutable List<Switch> switchesList_;         // switches which are requested in dict

        mutable List<Switch> switches_;

        mutable volScalarField nu_;

        //mutable volScalarField mu_;

        mutable volVectorField divTau_;

        mutable volVectorField IBDragPerV_;

        word densityFieldName_;

        const   volScalarField& rho_;

        mutable int     verboseDiskIntervall_;

        mutable int     verboseDiskCounter_;

        scalar cg_;

        scalar scaleDrag_;

        scalar scaleDragLocal_;

        scalar scaleTorque_;

        scalar scaleDH_;

public:

    //- Runtime type information
    TypeName("forceSubModel");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            forceSubModel,
            dictionary,
            (
                const dictionary& dict,
                cfdemCloud& sm,
                forceModel& fm
            ),
            (dict,sm,fm)
        );


    // Constructors

        //- Construct from components
        forceSubModel
        (
            const dictionary& dict,
            cfdemCloud& sm,
            forceModel& fm
        );


    // Destructor

        virtual ~forceSubModel();


    // Selector

        static autoPtr<forceSubModel> New
        (
            const dictionary& dict,
            cfdemCloud& sm,
            forceModel& fm,
            word forceType
        );


    // Member Functions
        void partToArray(const label&, const vector&, const vector&, const vector& Ufluid=vector::zero, scalar Cd=scalar(0),const vector& CdExtra=vector::zero) const;

        //function to set anisotropic drag components for DEM side only.
        void partToArrayAnisotropic(const label&, const vector&, const vector& dragEx=vector::zero) const;

        //function to set anisotropic torque components for DEM side only.
        void partToArrayAnisotropicTorque(const label&, const vector&, const vector& torqueTotal=vector::zero) const;

        virtual void constructorCalls(word typeName) const {}

        virtual void preParticleLoop(bool verbose=false) const {}

        virtual void postParticleLoop(bool verbose=false) const {}

        virtual void explicitCorr(vector&, vector&, scalar&, vector&, const vector&, vector&, const vector&, bool,label index=100) const;

        virtual void explicitCorrScalar(scalar&, scalar&, scalar&,
                                        scalar&, const scalar&, scalar&, bool, label index=100) const;

        virtual void update(label particleI,        label cellI,            scalar& d,
                            scalar& scalToUpdate1,  scalar&  scalToUpdate2, bool  verbose) const;

        virtual void update(label particleI,        label cellI,            scalar& d,
                            vector& vecToUpdate1,   vector&  vecToUpdate2,
                            scalar& scalToUpdate1,  scalar&  scalToUpdate2,
                            bool  verbose) const;

        virtual void scaleDia(scalar& d, int index=0) const;

        virtual void scaleForce(vector& force, scalar& d, int index=0) const;

        virtual void scaleCoeff(scalar& coeff, scalar& d, scalar& Rep, int index=0) const;

        virtual void scaleForceLocal(vector& force) const;

        virtual void scaleTorque(vector& torque) const;

        virtual void explicitLimit(vector&, vector&, scalar&) const;

        virtual void verboseToDiskWrite(Field<scalar>& writeValues) const {};//implement writing to disk

        virtual void calculateCorrectedVoidage(
            scalar& voidfraction, vector& Ufluid, scalar cellVolume, scalar diameter
        ) const {};

        virtual void  updateField(volScalarField* scalField, volVectorField* vecField) const {};

        virtual void  getField(vector position, label cellI, scalar& scalVal, vector& vecVal) const {};

    // Access
        inline bool treatForceExplicit() const { return switches_[0]; }

        inline bool treatForceDEM() const { return switches_[1]; }

        inline bool implForceDEM() const { return switches_[2]; }

        inline bool verbose() const      { return switches_[3]; }

        inline bool interpolation() const { return switches_[4]; }

        inline bool useFilteredDragModel() const { return switches_[5]; }

        inline bool useParcelSizeDependentFilteredDrag() const { return switches_[6]; }

        inline bool verboseToDisk() const
        {
            verboseDiskCounter_++;
            if(verboseDiskCounter_>=verboseDiskIntervall_)
            {
                verboseDiskCounter_=0;
                return switches_[9];
            }
            else
                return false;
        }

        virtual word myType() const=0;

        inline bool ms() const { return switches_[11]; }
        inline bool sq() const { return switches_[25]; }
        inline bool convex() const { return switches_[29]; }

        inline bool useTorque() const { return switches_[12]; }

        inline bool anisotropicDrag() const { return switches_[13]; }

        inline bool pullRotation() const { return switches_[14]; }

        inline bool pullOrientation() const { return switches_[15]; }

        inline bool pullShape() const { return switches_[16]; }

        inline bool voidageFunctionDiFelice() const { return switches_[17]; }

        inline bool voidageFunctionRong() const { return switches_[18]; }
        inline bool voidageFunctionTang() const { return switches_[23]; }

        inline bool useUf() const { return switches_[19]; }

        inline bool useFhydro() const { return switches_[20]; }

        inline bool implTorqueDEM() const { return switches_[21]; }

        inline bool useVisc() const { return switches_[22]; }

        inline bool useMpData() const { return switches_[24]; }

        inline bool useQuat() const { return switches_[26]; }

        inline bool pushTurbulence() const { return switches_[27]; }

        inline bool particleCG() const { return switches_[28]; }

        inline bool pullType() const { return switches_[30]; }

        inline bool pullDensity() const { return switches_[31]; }

        inline bool pushConvectiveHeatFlux() const { return switches_[32]; }

        inline bool pullTemp() const { return switches_[33]; }

        inline forceModel& myForceM() const { return forceModel_; }

        inline const List<Switch>& switches() const { return switches_; }

        inline const wordList& switchesNameList() const { return switchesNameList_; }

        void setSwitchesList(label i, bool v, bool defaultValue=false) const { switchesList_[i] = v; switches_[i]=defaultValue; }

        void setSwitches(label i, Switch v) const { switches_[i] = v; }

        virtual void readSwitches() const;

        void setupCommunication() const;

        const label& nrDefaultSwitches() const { return nrDefaultSwitches_; }

        const volScalarField& nuField() const;

        const volScalarField& muField() const;

        const volScalarField& rhoField() const;

        const volVectorField& divTauField(const volVectorField&) const;

        const volVectorField& IBDragPerV(const volVectorField&,const volScalarField&) const;

        inline scalar getCG(int index=0) const { return particleCloud_.cg(index); };

        const scalar scaleDrag(int index=0) const;

        const scalar& scaleDragLocal() const {return scaleDragLocal_;};

        const scalar& scaleTorque() const {return scaleTorque_;};

        virtual void calcXi(const scalar& ds, const scalar& voidfaction, const scalar& magUr, const scalar& nuf, scalar& Xi ) const=0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
