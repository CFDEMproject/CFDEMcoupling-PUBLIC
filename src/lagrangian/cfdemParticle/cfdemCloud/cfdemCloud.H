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

    Switch solveFlow_;

    bool verbose_;

    bool ignore_;

    const word modelType_;

    mutable double **positions_;

    mutable double **velocities_;

    mutable double **fluidVel_;

    mutable double **impForces_;

    mutable double **expForces_;

    mutable double **DEMForces_;

    mutable double **Cds_;

    mutable double **radii_;

    mutable double **voidfractions_;

    mutable double **cellIDs_;

    mutable double **particleWeights_;

    mutable double **particleVolumes_;

    int numberOfParticles_;

    bool numberOfParticlesChanged_;

    mutable bool arraysReallocated_;

    const wordList forceModels_;

    const wordList momCoupleModels_;

    const wordList liggghtsCommandModelList_;

    const word turbulenceModelType_;

    mutable scalar cg_;

    bool cgOK_;

    bool impDEMdrag_;

    mutable scalar imExSplitFactor_;

    bool useDDTvoidfraction_;

    mutable volScalarField ddtVoidfraction_;

    #if defined(version21) || defined(version16ext)
        #ifdef comp
            const compressible::turbulenceModel& turbulence_;
        #else
            const incompressible::turbulenceModel& turbulence_;
        #endif
    #elif defined(version15)
        const incompressible::RASModel& turbulence_;
    #endif

    autoPtr<forceModel>* forceModel_;

    autoPtr<locateModel> locateModel_;

    autoPtr<momCoupleModel>* momCoupleModel_;

    autoPtr<dataExchangeModel> dataExchangeModel_;

    autoPtr<IOModel> IOModel_;

    autoPtr<probeModel> probeModel_;

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

    virtual void findCells();

    virtual void setForces();

    virtual void setParticleForceField();

    virtual void setVectorAverages();

public:

    friend class dataExchangeModel;
    friend class voidFractionModel;
    friend class forceModel;

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
       void checkCG(bool);

       void setPos(double **&);

       word modelType(){ return modelType_; };

       label particleCell(int);

       double d(int);

       vector position(int);

       vector velocity(int);

       vector fluidVel(int);

       virtual const forceModel& forceM(int);

       virtual int nrForceModels();

       scalar radius(int);

       scalar voidfraction(int);

       label liggghtsCommandModelIndex(word);

       inline void setCG(double) const;

       inline const scalar& cg() const;

       inline const bool& impDEMdrag() const;

       inline const scalar& imExSplitFactor() const;

       inline const bool& ignore() const;

       inline const fvMesh& mesh() const;

       inline bool solveFlow() const;

       inline bool verbose() const;

       inline const IOdictionary& couplingProperties() const;

       inline double ** positions() const;

       inline double ** velocities() const;

       inline double ** fluidVels() const;

       inline double ** impForces() const;

       inline double ** expForces() const;

       inline double ** DEMForces() const;

       inline double ** Cds() const;

       inline double ** radii() const;

       inline double ** voidfractions() const;

       inline void get_radii(double**&) const;

       inline double ** cellIDs() const;

       inline void get_cellIDs(double**&) const;

       inline double ** particleWeights() const;

       inline double ** particleVolumes() const;

       virtual inline double dMin() {return -1;}
       virtual inline double dMax() {return -1;}
       virtual inline int minType() {return -1;}
       virtual inline int maxType() {return -1;}
       virtual inline bool multipleTypesDMax() {return false;}
       virtual inline bool multipleTypesDMin() {return false;}
       virtual inline double ** particleDensity() const {return NULL;};
       virtual inline int ** particleTypes() const {return NULL;};
       virtual label particleType(label index) const {return -1;};

       //access to the particles' orientation information
       virtual inline double ** exArray() const {return NULL;};
       virtual vector ex(int) const {return Foam::vector(0,0,0);};

       inline int numberOfParticles() const;

       inline bool numberOfParticlesChanged() const;

       inline int numberOfClumps() const;

       inline bool arraysReallocated() const;

       inline const wordList& forceModels();

       inline const voidFractionModel& voidFractionM() const;

       inline const locateModel& locateM() const;

       inline const momCoupleModel& momCoupleM(int) const;

       inline const dataExchangeModel& dataExchangeM() const;

       inline const IOModel& IOM() const;

       inline const probeModel& probeM() const;

       inline const averagingModel& averagingM() const;

       inline const clockModel& clockM() const;

       inline const smoothingModel& smoothingM() const;

       inline const meshMotionModel& meshMotionM() const;

       inline const wordList& liggghtsCommandModelList() const;

       inline autoPtr<liggghtsCommandModel>* liggghtsCommand() const;

       #if defined(version21) || defined(version16ext)
            #ifdef comp
                inline const compressible::turbulenceModel& turbulence() const;
            #else
                inline const incompressible::turbulenceModel& turbulence() const;
            #endif
       #elif defined(version15)
           inline const incompressible::RASModel& turbulence() const;
       #endif

    // Write

      // write cfdemCloud internal data
        virtual bool evolve(volScalarField&,volVectorField&,volVectorField&);

        virtual bool reAllocArrays() const;

    // IO
        void writeScalarFieldToTerminal(double**&);

        void writeVectorFieldToTerminal(double**&);

    // functions
        tmp<fvVectorMatrix> divVoidfractionTau(volVectorField& ,volScalarField&) const;

        tmp<volScalarField> ddtVoidfraction() const;

        void calcDdtVoidfraction(volScalarField& voidfraction) const;

        //tmp<fvVectorMatrix> ddtVoidfractionU(volVectorField& ,volScalarField&) const;

        tmp<volScalarField> voidfractionNuEff(volScalarField&) const;

        void resetArray(double**&,int,int,double resetVal=0.);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cfdemCloudI.H"

#endif

// ************************************************************************* //
