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

    cloud class managing DEM data for CFD-DEM coupling

Class
    Foam::cfdemCloud

SourceFiles
    cfdemCloud.C
    cfdemCloudIO.C

\*---------------------------------------------------------------------------*/

#ifndef cfdemCloud_H
#define cfdemCloud_H

// choose version
#include "OFversion.H"
#include <vector>

#include "fvCFD.H"
#include "IFstream.H"

#if defined(version21) || defined(version16ext)
    #include "turbulenceModel.H"
#elif defined(version15)
    #include "RASModel.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// forward declarations
class forceModel;
class locateModel;
class voidFractionModel;
class dataExchangeModel;
class IOModel;
class probeModel;
class registryModel;
class averagingModel;
class clockModel;
class smoothingModel;
class momCoupleModel;
class meshMotionModel;
class liggghtsCommandModel;

/*---------------------------------------------------------------------------*\
                           Class cfdemCloud Declaration
\*---------------------------------------------------------------------------*/

class cfdemCloud
{

// protected data
protected:
    const fvMesh& mesh_;

    IOdictionary couplingProperties_;

    IOdictionary liggghtsCommandDict_;

    Switch allowAdjustTimeStep_;

    Switch solveFlow_;

    Switch solveScalarTransport_;

    Switch verbose_;

    Switch expCorrDeltaUError_;

    bool debug_;

    bool allowCFDsubTimestep_;

    bool ignore_;

    bool writeTimePassed_;

    bool resetWriteTimePassed_;

    const word modelType_;

    mutable double **impForces_;

    mutable double **expForces_;

    mutable double **voidfractions_;

    mutable double **cellIDs_;

    mutable double **particleWeights_;

    mutable double **particleVolumes_;

    mutable double **particleV_;

    mutable double **dragPrev_;

    int numberOfParticles_;

    //mutable int maxCellsPerParticle_; // lives now in voidFractionModel.H

    scalar d32_;

    bool numberOfParticlesChanged_;

    mutable bool arraysReallocated_;

    const wordList forceModels_;

    const wordList momCoupleModels_;

    const wordList liggghtsCommandModelList_;

    const word turbulenceModelType_;

    mutable bool isLES_;

    mutable scalar cg_;

    mutable std::vector<double> cgTypeSpecific_;

    bool cgOK_;

    bool impDEMdrag_;

    bool impDEMdragAcc_;

    mutable scalar imExSplitFactor_;

    mutable bool treatVoidCellsAsExplicitForce_; //will treat the coupling force in cells with no Us data explicitly

    mutable word useDDTvoidfraction_;

    mutable volScalarField ddtVoidfraction_;

    //Variable used to de-activate mirroring across periodic boundary conditions
    mutable Switch checkPeriodicCells_;
    //de-activation and tolerance variables
    //if set to  (for a specific direction), the periodic check will NOT be done
    //default = (1,1,1), i.e., periodic checks will be done
    //important for probing ambient points. Only read-in in case checkPeriodicCells is active
    mutable vector wall_periodicityCheckRange_;
    scalar         wall_periodicityCheckTolerance_;

    bool meshHasUpdated_;

    #if defined(version24Dev)
        const turbulenceModel& turbulence_;
    #elif defined(version21) || defined(version16ext)
        #ifdef compre
            const compressible::turbulenceModel& turbulence_;
        #else
            const incompressible::turbulenceModel& turbulence_;
        #endif
    #elif defined(version15)
        const incompressible::RASModel& turbulence_;
    #endif

    //Multiphase Turbulence (e.g., slip-induced turbulence)
    mutable volScalarField turbulenceMultiphase_;

    autoPtr<forceModel>* forceModel_;

    autoPtr<locateModel> locateModel_;

    autoPtr<momCoupleModel>* momCoupleModel_;

    autoPtr<dataExchangeModel> dataExchangeModel_;

    autoPtr<IOModel> IOModel_;

    autoPtr<probeModel> probeModel_;

    autoPtr<registryModel> registryModel_;

    autoPtr<voidFractionModel> voidFractionModel_;

    autoPtr<averagingModel> averagingModel_;

    autoPtr<clockModel> clockModel_;

    autoPtr<smoothingModel> smoothingModel_;

    autoPtr<meshMotionModel> meshMotionModel_;

    autoPtr<liggghtsCommandModel>* liggghtsCommand_;

// Protected member functions
    virtual void getDEMdata();

    virtual void giveDEMdata();

    virtual void setNumberOfParticles(int);

    virtual void setNumberOfClumps(int);

    virtual void setPositionsCM(label,double*);

    virtual void setCellIDsCM(label,int*);

    virtual void findCells();

    virtual void setForces();

    virtual void setVoidFraction();

    virtual void resetVoidFraction();

    virtual void setAlpha(volScalarField&);

    virtual void setParticleForceField();

    virtual void setVectorAverages();

public:

    friend class dataExchangeModel;
    friend class voidFractionModel;
    friend class forceModel;
    friend class forceSubModel;
    friend class generalManual;
    friend class temperatureModel;

// Constructors

    //- Construct from mesh and a list of particles
    cfdemCloud
    (
            const fvMesh& mesh
    );

    //- Destructor
    virtual ~cfdemCloud();

// public Member Functions

    // Access
       const turbulenceModel& turbulence() {return turbulence_;};

       bool isLES() const {return isLES_; };

       bool debugMode(){ return debug_;};

       bool allowCFDsubTimestep(){ return allowCFDsubTimestep_;};

       void setAllowCFDsubTimestep(bool b){allowCFDsubTimestep_=b;};

       bool expCorrDeltaUError(){ return expCorrDeltaUError_;};

       void setExpCorrDeltaUError(){expCorrDeltaUError_=couplingProperties_.lookupOrDefault<Switch>("expCorrDeltaUError", false);};

       void checkCG(bool);

       void setPos(double **&);

       word modelType(){ return modelType_; };

       label particleCell(int);

       virtual inline vector position(int, bool ovr=false, int mp=-1);

       virtual inline vector velocity(int, bool ovr=false, int mp=-1);

       vector fluidVel(int);

       virtual const forceModel& forceM(int);

       virtual int nrForceModels();

       virtual double** cellsPerParticle();

       scalar voidfraction(int);

       label liggghtsCommandModelIndex(word);

       inline void setCG(double) const;

       inline void setCGTypeSpecific(int, double);

       mutable bool cgParticleSpecific_;

       mutable bool cgTypeSpecificDifferent_;

       inline const scalar& cg() const;

       inline scalar cg(int index) const;

       inline const bool& impDEMdrag() const;

       inline const bool& impDEMdragAcc() const;

       inline const scalar& imExSplitFactor() const;

       inline const bool& treatVoidCellsAsExplicitForce() const;

       inline const bool& ignore() const;

       inline const bool& writeTimePassed() const;

       inline void resetWriteTimePassed();

       inline const fvMesh& mesh() const;

       inline bool allowAdjustTimeStep() const;

       inline bool solveFlow() const;

       inline bool solveScalarTransport() const;

       inline bool verbose() const;

       inline const IOdictionary& couplingProperties() const;

       inline double ** positions() const;

       inline double ** velocities() const;

       inline double ** fluidVels() const;

       inline double ** fAccs() const;

       inline double ** impForces() const;

       inline double ** expForces() const;

       inline double ** Cds() const;

       inline double ** radii() const;

       inline double ** voidfractions() const;

       inline void get_radii(double**&) const;

       virtual inline double ** cellIDs(bool ovr=false) const;

       virtual inline label cellID(int, bool ovr=false, int mp=-1) const;

       inline void get_cellIDs(double**&) const;

       inline double ** particleWeights() const;

       virtual inline label body(int);

       virtual inline double particleVolume(int);

       inline scalar radius(int, bool ovr=false);

       virtual inline double diameter(int, bool ovr=false, int mp=-1);

       inline scalar d32(bool recalc=true);
       label particleType(label index) const
       {
            if(idPartTypes()>-1) return fieldsToDEM[idPartTypes()][index][0];
            return 1;
       };

       //access to the particle's rotation and torque data
       virtual inline double ** DEMTorques() const {return NULL;};
       virtual inline double ** omegaArray() const {return NULL;};
       virtual vector omega(int) const {return Foam::vector(0,0,0);};

//        //access to the particles' orientation information
       virtual vector ex(int) const {
         FatalError << "cfdemCloud::ex(): you are attempting to get an access to ex array for spheroids/superquadrics. This is not allowed for spheres. Recompile the code with #define anisotropicRotation in OFversion.H or Are you using the correct solver?" << abort(FatalError);
         return Foam::vector(0,0,0);
       };

       virtual vector ey(int) const {
         FatalError << "cfdemCloud::ey(): you are attempting to get an access to ey array for spheroids/superquadrics. This is not allowed for spheres. Recompile the code with #define anisotropicRotation in OFversion.H or Are you using the correct solver?" << abort(FatalError);
         return Foam::vector(0,0,0);
       };

       virtual vector shape(int) const {
         FatalError << "cfdemCloud::shapeArray(): you are attempting to get an access to shape array for spheroids/superquadrics. This is not allowed for spheres. Recompile the code with #define anisotropicRotation in OFversion.H or Are you using the correct solver?" << abort(FatalError);
         return Foam::vector(0,0,0);
       };

       // access to complex shape information
       virtual scalar volume(int) const {
         FatalError << "cfdemCloud::volume(): you are attempting to get an access to volume array for superquadrics. This is not allowed for spheres. Are you using the correct solver?" << abort(FatalError);
         return 0;
       };

       virtual scalar area(int) const {
         FatalError << "cfdemCloud::area(): you are attempting to get an access to area array for superquadrics. This is not allowed for spheres. Are you using the correct solver?" << abort(FatalError);
         return 0;
       };
       virtual inline scalarList& area() { static scalarList foo(0); return foo; }

       virtual quaternion quat(int) const {
         FatalError << "cfdemCloud::quat(): you are attempting to get an access to quaternion array for superquadrics. This is not allowed for spheres. Are you using the correct solver?" << abort(FatalError);
         return quaternion(1);
       };

       virtual vector2D blockiness(int) const {
         FatalError << "cfdemCloud::blockiness(): you are attempting to get an access to blockiness array for superquadrics. This is not allowed for spheres. Are you using the correct solver?" << abort(FatalError);
         return Foam::vector2D(0,0);
       };

       // MS access functions
#define ErrorMsg FatalError<<"cfdemCloud: You are attempting to use a MS interface routine. However, your solver is not compiled to use multispheres."<<abort(FatalError)
       virtual inline int nClumpTypes() const {ErrorMsg; return 0;};
       virtual inline int clumpType(int i) const {ErrorMsg; return 0;};
       virtual inline double **& clumpVol() const {ErrorMsg; static double **foo = nullptr; return foo;};
       virtual inline double **& clumpArea() const {ErrorMsg; static double **foo = nullptr; return foo;};
       virtual inline double **& clumpDH() const {ErrorMsg; static double **foo = nullptr; return foo;};
       virtual inline label clumpIndexOfParticle(int) {ErrorMsg; return 0;};
       virtual inline labelList& particlesOfClumpIndex(int) {ErrorMsg; static labelList foo(0); return foo;};
       virtual inline bool useforcePerClump() const {ErrorMsg; return false;};
       virtual inline label nrigid(int) {ErrorMsg; return 0;};
       virtual inline double **& impForcesCM() const {ErrorMsg; static double **foo = nullptr; return foo;};
       virtual inline double **& expForcesCM() const {ErrorMsg; static double **foo = nullptr; return foo;};

       virtual inline scalarList& DHc() {ErrorMsg; static scalarList foo(0); return foo;};
       virtual inline void setDHc(scalarList) {FatalError << "setDHc in cfdemCloud" << abort(FatalError);};
       virtual inline void setArea(scalarList) {FatalError << "setArea in cfdemCloud" << abort(FatalError);};

       virtual inline vector resolution() {FatalError << "resolution() in cfdemCloud" << abort(FatalError); return vector(0,0,0);};
       virtual inline void initMpData() {FatalError << "setMpData in cfdemCloud" << abort(FatalError);};
       virtual inline void setMpData(int, int, int,  vector, vector, scalar) {/*do nothing*/};
       virtual inline void setMpDataPoint(int, int, int, int, double) {FatalError << "setMpDataPoint in cfdemCloud" << abort(FatalError);};;
       virtual inline word shapeTypeName() {return "sphere";};

       // wrapper function for volume calculation based on particle radius
       virtual inline scalar particleVolumeScaled(int index, const scalarList& scaleVol)
       {
           return 4./3.*M_PI * radius(index)*radius(index)*radius(index) * scaleVol[0];
       };

       //Optional Extra arrays for implicit coupling with anisotropic drag
       virtual inline double ** omegaFluid() const { return NULL;};     //the fluid rotation rate
       double **dragPrev();

       //Detector if SRF module is enable or not
       virtual inline bool SRFOn(){return false;}

       //allow handshake with refinement field
       virtual void setRefinementField(volScalarField* refine_) {};

       virtual inline int numberOfObjects(bool ovr=false) const;

       inline int numberOfParticles() const;

       inline bool numberOfParticlesChanged() const;

       virtual inline int numberOfClumps() const;

       //inline void setMaxCellsPerParticle(int) const;

       //inline int maxCellsPerParticle() const;

       inline bool arraysReallocated() const;

       inline const wordList& forceModels();

       inline const voidFractionModel& voidFractionM() const;

       inline const locateModel& locateM() const;

       inline const momCoupleModel& momCoupleM(int) const;

       inline const dataExchangeModel& dataExchangeM() const;

       inline const IOModel& IOM() const;

       inline const probeModel& probeM() const;

       inline const registryModel& registryM() const;

       inline const averagingModel& averagingM() const;

       inline const clockModel& clockM() const;

       inline const smoothingModel& smoothingM() const;

       inline const meshMotionModel& meshMotionM() const;

       inline const wordList& liggghtsCommandModelList() const;

       inline autoPtr<liggghtsCommandModel>* liggghtsCommand() const;

        inline void makeSpecific(volScalarField&);

        inline void makeSpecific(volVectorField&);

        inline void scaleWithVcell(volScalarField&);

        inline void scaleWithVcell(volVectorField&);

    // Write

      // write cfdemCloud internal data
        virtual bool evolve(volScalarField&,volVectorField&,volVectorField&);

        virtual bool reAllocArrays() const;

        virtual bool reAllocArrays(int nP, bool forceRealloc) const;  //force number of particles during reallocation


    // IO
        void writeScalarFieldToTerminal(double**&);

        void writeVectorFieldToTerminal(double**&);

    // functions
        tmp<fvVectorMatrix> divVoidfractionTau(volVectorField& ,volScalarField&) const;

        tmp<fvVectorMatrix> divVoidfractionTau(volVectorField& ,volScalarField& ,volScalarField&) const;

        tmp<volScalarField> ddtVoidfraction() const;

        void calcDdtVoidfraction(volScalarField& voidfraction, volVectorField& Us) const;

        void calcMultiphaseTurbulence();

        volScalarField& turbulenceMultiphase() const {return turbulenceMultiphase_;};

        //tmp<fvVectorMatrix> ddtVoidfractionU(volVectorField& ,volScalarField&) const;

        tmp<volScalarField> voidfractionNuEff(volScalarField&) const;

        void resetArray(double**&,int,int,double resetVal=0.);

        std::vector< std::vector<double*> >* getVprobe();

        std::vector< std::vector<double> >* getSprobe();

        //=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
        // EXPERIMENAL CODE
        //SHARED GLOBAL PARTICLE Arrays (public containers for data management)
        //1-USER-defined particle fields for CFDEM-->ExternalCode transfer (e.g., to handle fluxes)
        mutable std::vector<word>       namesfieldsToDEM;                    // names of field which are communicated to DEM
        mutable std::vector<word>       typesfieldsToDEM;                    // types of field which are communicated to DEM
        mutable std::vector<bool>       pullfieldsToDEM;                     // types of field which are pulled (only) from DEM
        mutable std::vector<double**>   fieldsToDEM;                         // fields which are communicated to DEM

        void registerFieldsToDEM(word name, word type, int& id, bool pull=false);
        bool checkAndRegisterFieldsToDEM(const wordList names,const word types,std::vector<int> & id);
        virtual void allocateFieldsToDEM(word shapeType="atom") const;
        virtual void reAllocateFieldsToDEM(word shapeType="atom") const;
        virtual void zeroizeFieldsToDEM(word shapeType="atom");
        void accessFieldsToDEM(word name, double **& field);
        int existsFieldsToDEM(word name);

        mutable int idRadius_; // id of particle radius in DEM in modular data exchange framework
        mutable int idPos_; // id of particle position in DEM in modular data exchange framework
        mutable int idVel_; // id of particle velocity in DEM in modular data exchange framework
        mutable int idFacc_; // id of accumulated force on particle in DEM in modular data exchange framework
        mutable int idPartTypes_; // id of particle type in DEM in modular data exchange framework
        mutable int idDragExp_; // id of explicit drag in DEM in modular data exchange framework
        mutable int idKsl_; // id of implicit drag coefficient in DEM in modular data exchange framework
        mutable int idKslExtra_; // id of implicit anisotropic drag coefficient in DEM in modular data exchange framework
        mutable int idUf_; // id of fluid vel in DEM in modular data exchange framework
        mutable int idTorqueExp_; // id of explicit torque in DEM in modular data exchange framework
        mutable int idKslRotation_; // id of implicit torque coefficient in DEM in modular data exchange framework
        mutable int idPullRotation_; // id of rotation omega in DEM in modular data exchange framework
        mutable int idPullOrientation_; // id of orientation (ex) in DEM in modular data exchange framework
        mutable int idPullOrientation1_; // id of orientation (ey) in DEM in modular data exchange framework
        mutable int idPullShape_; // id of orientation in DEM in modular data exchange framework
        mutable int idDragExpCM_; // id of explicit drag for multisphere in DEM in modular data exchange framework
        mutable int idKslCM_; // id of implicit anisotropic drag for multisphere in DEM in modular data exchange framework
        mutable int idKslExtraCM_; // id of implicit drag for multisphere in DEM in modular data exchange framework
        mutable int idUfCM_; // id of fluid vel for multisphere in DEM in modular data exchange framework
        mutable int idTorqueExpCM_; // id of explicit torque for multisphere in DEM in modular data exchange framework
        mutable int idKslRotationCM_; // id of implicit torque coefficient for multisphere in DEM in modular data exchange framework
        mutable int idFhydro_; // id of void fraction in DEM in modular data exchange framework
        mutable int idVisc_; // id of fluid vel in DEM in modular data exchange framework
        mutable int idBlockiness_; // id of sq blockiness in DEM in modular data exchange framework
        mutable int idArea_; // id of sq area in DEM in modular data exchange framework
        mutable int idVol_; // id of sq volume in DEM in modular data exchange framework
        mutable int idQuat_; // id of quaternion in DEM in modular data exchange framework
        mutable int idK_; // id of turbulence k in DEM in modular data exchange framework
        mutable int idEpsilon_; // id of turbulence epsilon in DEM in modular data exchange framework
        mutable int idParticleCG_; // id of particle specifig CG factor in DEM in modular data exchange framework
        mutable int idMass_; // id of convex volume in DEM in modular data exchange framework
        mutable int idDensity_; // id of convex volume in DEM in modular data exchange framework
        mutable int idType_; // id of particle type in DEM in modular data exchange framework
        mutable int idConvectiveHeatFlux_; // id of particle convectiveHeatFlux in DEM in modular data exchange framework
        mutable int idTemp_; // id of particle convectiveHeatFlux in DEM in modular data exchange framework

        int& idRadius() const {return idRadius_;};
        int& idPos() const {return idPos_;};
        int& idVel() const {return idVel_;};
        int& idFacc() const {return idFacc_;};
        int& idPartTypes() const {return idPartTypes_;};
        int& idDragExp() const {return idDragExp_;};
        int& idKsl() const {return idKsl_;};
        int& idKslExtra() const {return idKslExtra_;};
        int& idUf() const {return idUf_;};
        int& idTorqueExp() const {return idTorqueExp_;};
        int& idKslRotation() const {return idKslRotation_;};
        int& idPullRotation() const {return idPullRotation_;};
        int& idPullOrientation() const {return idPullOrientation_;};
        int& idPullOrientation1() const {return idPullOrientation1_;};
        int& idPullShape() const {return idPullShape_;};
        int& idDragExpCM() const {return idDragExpCM_;};
        int& idKslCM() const {return idKslCM_;};
        int& idKslExtraCM() const {return idKslExtraCM_;};
        int& idUfCM() const {return idUfCM_;};
        int& idTorqueExpCM() const {return idTorqueExpCM_;};
        int& idKslRotationCM() const {return idKslRotationCM_;};
        int& idFhydro() const {return idFhydro_;};
        int& idVisc() const {return idVisc_;};
        int& idBlockiness() const {return idBlockiness_;};
        int& idArea() const {return idArea_;};
        int& idVol() const {return idVol_;};
        int& idQuat() const {return idQuat_;};
        int& idK() const {return idK_;};
        int& idEpsilon() const {return idEpsilon_;};
        int& idParticleCG() const {return idParticleCG_;};
        int& idMass() const {return idMass_;};
        int& idDensity() const {return idDensity_;};
        int& idType() const {return idType_;};
        int& idConvectiveHeatFlux() const {return idConvectiveHeatFlux_;};
        int& idTemp() const {return idTemp_;};
        //=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
        bool checkPeriodicCells() { return checkPeriodicCells_; }
        int  wall_periodicityCheckRange(int dir) {return static_cast<int>(wall_periodicityCheckRange_[dir]);}
        vector  wall_periodicityCheckRange() {return wall_periodicityCheckRange_;}
        bool meshHasUpdated() { return meshHasUpdated_; }
        void setMeshHasUpdatedFlag(bool hasUpdated) { meshHasUpdated_ = hasUpdated; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cfdemCloudI.H"

#endif

// ************************************************************************* //
