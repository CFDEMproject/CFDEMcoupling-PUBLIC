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
\*---------------------------------------------------------------------------*/

#include "fileName.H"
#include "cfdemCloud.H"
#include "global.H"
#include "forceModel.H"
#include "locateModel.H"
#include "momCoupleModel.H"
#include "meshMotionModel.H"
#include "voidFractionModel.H"
#include "dataExchangeModel.H"
#include "IOModel.H"
#include "probeModel.H"
#include "registryModel.H"
#include "averagingModel.H"
#include "clockModel.H"
#include "smoothingModel.H"
#include "liggghtsCommandModel.H"
#include "cyclicPolyPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::cfdemCloud::cfdemCloud
(
    const fvMesh& mesh
)
:
    mesh_(mesh),
    couplingProperties_
    (
        IOobject
        (
            "couplingProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    liggghtsCommandDict_
    (
        IOobject
        (
            "liggghtsCommands",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    allowAdjustTimeStep_(couplingProperties_.lookupOrDefault<Switch>("allowAdjustTimeStep", false)),
    solveFlow_(true),
    solveScalarTransport_(true),
    verbose_(couplingProperties_.lookupOrDefault<Switch>("verbose", false)),
    expCorrDeltaUError_(false),
    debug_(false),
    allowCFDsubTimestep_(true),
    ignore_(false),
    writeTimePassed_(false),
    modelType_(couplingProperties_.lookup("modelType")),
    positions_(NULL),
    velocities_(NULL),
    fluidVel_(NULL),
    fAcc_(NULL),
    impForces_(NULL),
    expForces_(NULL),
    DEMForces_(NULL),
    Cds_(NULL),
    radii_(NULL),
    voidfractions_(NULL),
    cellIDs_(NULL),
    particleWeights_(NULL),
    particleVolumes_(NULL),
    particleV_(NULL),
    dragPrev_(NULL),
    numberOfParticles_(0),
    d32_(-1),
    numberOfParticlesChanged_(false),
    arraysReallocated_(false),
    forceModels_(couplingProperties_.lookup("forceModels")),
    momCoupleModels_(couplingProperties_.lookup("momCoupleModels")),
    liggghtsCommandModelList_(liggghtsCommandDict_.lookup("liggghtsCommandModels")),
    turbulenceModelType_(couplingProperties_.lookup("turbulenceModelType")),
    isLES_(false),
    cg_(1.),
    cgOK_(true),
    impDEMdrag_(false),
    impDEMdragAcc_(false),
    imExSplitFactor_(1.0),
    treatVoidCellsAsExplicitForce_(false),
    useDDTvoidfraction_("off"),
    ddtVoidfraction_
    (   
        IOobject
        (
            "ddtVoidfraction",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0,0,-1,0,0), 0)  // 1/s
    ),
    checkPeriodicCells_(false),
    wall_periodicityCheckRange_(vector(1,1,1)),
    wall_periodicityCheckTolerance_(1e-07),
    meshHasUpdated_(false),
    turbulence_
    (
        #if defined(version24Dev)
            mesh.lookupObject<turbulenceModel>
        #elif defined(version21) || defined(version16ext)
            #ifdef compre
                mesh.lookupObject<compressible::turbulenceModel>
            #else
                mesh.lookupObject<incompressible::turbulenceModel>
            #endif
        #elif defined(version15)
            mesh.lookupObject<incompressible::RASModel>
        #endif
        (
            turbulenceModelType_
        )
    ),
    turbulenceMultiphase_
    (   
        IOobject
        (
            "turbulenceMultiphase",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
#ifdef compre
        dimensionedScalar("zero", dimensionSet(1,-1,-1,0,0), 0)  // kg/m/s
#else
        dimensionedScalar("zero", dimensionSet(0,2,-1,0,0), 0)  // m²/s
#endif
    ),    
    locateModel_
    (
        locateModel::New
        (
            couplingProperties_,
            *this
        )
    ),
    /*momCoupleModel_
    (
        momCoupleModel::New
        (
            couplingProperties_,
            *this
        )
    ),*/
    dataExchangeModel_
    (
        dataExchangeModel::New
        (
            couplingProperties_,
            *this
        )
    ),
    IOModel_
    (
        IOModel::New
        (
            couplingProperties_,
            *this
        )
    ),
    probeModel_
    (
        probeModel::New
        (
            couplingProperties_,
            *this,
            const_cast<char *>("none"),
            const_cast<char *>("none")
        )
    ),
    registryModel_
    (
        registryModel::New
        (
            couplingProperties_,
            *this
        )
    ),
    voidFractionModel_
    (
        voidFractionModel::New
        (
            couplingProperties_,
            *this
        )
    ),
    averagingModel_
    (
        averagingModel::New
        (
            couplingProperties_,
            *this
        )
    ),
    clockModel_
    (
        clockModel::New
        (
            couplingProperties_,
            *this
        )
    ),
    smoothingModel_
    (
        smoothingModel::New
        (
            couplingProperties_,
            *this
        )
    ),
    meshMotionModel_
    (
        meshMotionModel::New
        (
            couplingProperties_,
            *this
        )
    )
{
    #include "versionInfo.H"
    global buildInfo(couplingProperties_,*this);
    buildInfo.info();

    //-- apply debug Mode to sub models

    // set debug flag according to env
    debug_ = buildInfo.debugMode();

    // overwrite debug flag if found in dict
    if (couplingProperties_.found("debug"))
        debug_=Switch(couplingProperties_.lookup("debug"));

    // apply flag
    if(!debugMode()) ddtVoidfraction_.writeOpt() = IOobject::NO_WRITE;
    if(!debugMode()) turbulenceMultiphase_.writeOpt() = IOobject::NO_WRITE;
    voidFractionM().applyDebugSettings(debugMode());
    averagingM().applyDebugSettings(debugMode());
    //--

    //push dummy to type-specific cg factor since types start with 1
    cgTypeSpecific_.push_back(-1);
    cgTypeSpecificDifferent=false;
    dataExchangeM().setCG();

    Info << "If BC are important, please provide volScalarFields -imp/expParticleForces-" << endl;
    if (couplingProperties_.found("solveFlow"))
        solveFlow_=Switch(couplingProperties_.lookup("solveFlow"));
    if (couplingProperties_.found("solveScalarTransport"))
        solveScalarTransport_=Switch(couplingProperties_.lookup("solveScalarTransport"));
    if (couplingProperties_.found("imExSplitFactor"))
        imExSplitFactor_ = readScalar(couplingProperties_.lookup("imExSplitFactor"));

    if(imExSplitFactor_>1.0)
            FatalError  << "You have set imExSplitFactor > 1 in your couplingProperties. Must be <=1."
                       << abort(FatalError);
    if(imExSplitFactor_<0.0)
            FatalError  << "You have set imExSplitFactor < 0 in your couplingProperties. Must be >=0"
                       << abort(FatalError);

    if (couplingProperties_.found("treatVoidCellsAsExplicitForce"))
        treatVoidCellsAsExplicitForce_ = readBool(couplingProperties_.lookup("treatVoidCellsAsExplicitForce"));
    if (couplingProperties_.found("ignore")) ignore_=true;
    if (turbulenceModelType_=="LESProperties")
    {
        isLES_ = true;
        Info << "WARNING - LES functionality not yet tested!" << endl;
    }

    if (couplingProperties_.found("useDDTvoidfraction"))
    {
        useDDTvoidfraction_=word(couplingProperties_.lookup("useDDTvoidfraction"));

        if(useDDTvoidfraction_==word("a") || 
           useDDTvoidfraction_==word("b") ||
           useDDTvoidfraction_==word("off")
          )
            Info << "choice for ddt(voidfraction) = " << useDDTvoidfraction_ << endl;
        else
            FatalError << "Model " << useDDTvoidfraction_ 
                       << " is not a valid choice for ddt(voidfraction). Choose a or b or off."
                       << abort(FatalError);
    }
    else        
        Info << "ignoring ddt(voidfraction)" << endl;

    momCoupleModel_ = new autoPtr<momCoupleModel>[momCoupleModels_.size()];
    for (int i=0;i<momCoupleModels_.size();i++)
    {
        momCoupleModel_[i] = momCoupleModel::New
        (
            couplingProperties_,
            *this,
            momCoupleModels_[i]
        );
        momCoupleModel_[i]().applyDebugSettings(debugMode());
        registryM().addProperty(momCoupleModel_[i]().type()+"_index",i);
    }

    forceModel_ = new autoPtr<forceModel>[nrForceModels()];
    for (int i=0;i<nrForceModels();i++)
    {
        forceModel_[i] = forceModel::New
        (
            couplingProperties_,
            *this,
            forceModels_[i]
        );
        forceModel_[i]().applyDebugSettings(debugMode());
    }

    if (nrForceModels()<SMALL)
        FatalError  << "Please use at least one forceModel ! "
                    << "(e.g. noDrag) \n"
                    << abort(FatalError);

    // run liggghts commands from cfdem
    liggghtsCommand_ = new autoPtr<liggghtsCommandModel>[liggghtsCommandModelList_.size()];
    for (int i=0;i<liggghtsCommandModelList_.size();i++)
    {
        liggghtsCommand_[i] = liggghtsCommandModel::New
        (
            liggghtsCommandDict_,
            *this,
            liggghtsCommandModelList_[i],
            i
        );
    }

    Switch cgWarnOnly_(couplingProperties_.lookupOrDefault<Switch>("cgWarnOnly", false));
    if (!cgOK_ && cg_ > 1)
    {
        if(cgWarnOnly_)
            Warning<< "at least one of your models is not fit for cg !!!"<< endl; 
        else
            FatalError<< "at least one of your models is not fit for cg !!!"<< abort(FatalError); 
    }

    // check if sim is fully peridic box
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    int nPatchesCyclic(0);
    int nPatchesNonCyclic(0);
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        #if defined(versionExt32)
        if (isA<cyclicPolyPatch>(pp))
            nPatchesCyclic++;
        else if (!isA<processorPolyPatch>(pp))
            nPatchesNonCyclic++;
        #else
        if (isA<cyclicPolyPatch>(pp) || isA<cyclicAMIPolyPatch>(pp))
            nPatchesCyclic++;
        else if (!isA<processorPolyPatch>(pp))
            nPatchesNonCyclic++;
        #endif
    }
    if(nPatchesNonCyclic==0)
        checkPeriodicCells_=true;

    //hard set checkperiodic cells if wished
    if(this->couplingProperties().found("checkPeriodicCells"))
        checkPeriodicCells_ = couplingProperties().lookupOrDefault<Switch>("checkPeriodicCells", checkPeriodicCells_);

    if(nPatchesCyclic>0 && nPatchesNonCyclic>0)
    {
        if(verbose_) Info << "nPatchesNonCyclic=" << nPatchesNonCyclic << ", nPatchesCyclic=" << nPatchesCyclic << endl;
        Warning << "Periodic handing is disabled because the domain is not fully periodic!\n" << endl;
    }

    //Check if user attempts to change fluid time step
    if( mesh_.time().controlDict().lookupOrDefault<Switch>("adjustTimeStep", false) && !allowAdjustTimeStep_ )
    {
        FatalError << "cfdemCloud:: you want to adjustTimeStep in controlDict. This is not allowed in this version of CFDEM."
                   << abort(FatalError);
    }

    //Check if mesh is empty (decomposition error)
    if(mesh_.cells().size() < 1) 
    {
        int nprocs;
        MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
        if(nprocs > 2)  FatalError << endl << "cfdemCloud:: local mesh has zero cells. Please check the mesh and the decomposition!" << abort(FatalError);
        
        Pout << "WARNING: cfdemCloud:: local mesh has zero cells. Please check the mesh and the decomposition!" << endl;
    }
}

// * * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * //
Foam::cfdemCloud::~cfdemCloud()
{
    clockM().evalPar();
    clockM().normHist();
    dataExchangeM().destroy(positions_,3);
    dataExchangeM().destroy(velocities_,3);
    dataExchangeM().destroy(fluidVel_,3);
    dataExchangeM().destroy(fAcc_,3);
    dataExchangeM().destroy(impForces_,3);
    dataExchangeM().destroy(expForces_,3);
    dataExchangeM().destroy(DEMForces_,3);
    dataExchangeM().destroy(Cds_,1);
    dataExchangeM().destroy(radii_,1);
    dataExchangeM().destroy(voidfractions_,1);
    dataExchangeM().destroy(cellIDs_,1);
    dataExchangeM().destroy(particleWeights_,1);
    dataExchangeM().destroy(particleVolumes_,1);
    dataExchangeM().destroy(particleV_,1);

    //=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
    int iUser=0;
    for( std::vector<double**>::iterator 
         it  = particleDatFieldsUserCFDEMToExt.begin(); 
         it != particleDatFieldsUserCFDEMToExt.end(); 
       ++it)
    {
        Info << "cfdemCloud destroys UserCFDEM data: " << namesFieldsUserCFDEMToExt[iUser++] << endl;
        dataExchangeM().destroy((*it),1);
    }
    //=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
}
// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
void Foam::cfdemCloud::getDEMdata()
{
    if(verbose_) Info << "Foam::cfdemCloud::getDEMdata()" << endl;
    dataExchangeM().getData("radius","scalar-atom",radii_);
    dataExchangeM().getData("x","vector-atom",positions_);
    dataExchangeM().getData("v","vector-atom",velocities_);

    if(impDEMdragAcc_)
        dataExchangeM().getData("dragAcc","vector-atom",fAcc_); // array is used twice - might be necessary to clean it first

    if(verbose_) Info << "Foam::cfdemCloud::getDEMdata() - done." << endl;
}

void Foam::cfdemCloud::giveDEMdata()
{

    dataExchangeM().giveData("dragforce","vector-atom",DEMForces_);

    if(impDEMdrag_)
    {
        if(verbose_) Info << "sending Ksl and uf" << endl;
        dataExchangeM().giveData("Ksl","scalar-atom",Cds_);
        dataExchangeM().giveData("uf","vector-atom",fluidVel_);
    }
    if(verbose_) Info << "giveDEMdata done." << endl;
}

//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
// * * *   write top level fields   * * * //
void Foam::cfdemCloud::giveUSERdata()
{
    //Handover USER-defined data 
    for(std::vector<word>::iterator it = namesFieldsUserCFDEMToExt.begin(); it != namesFieldsUserCFDEMToExt.end(); ++it) 
    {
        int positionInRegister = std::distance(namesFieldsUserCFDEMToExt.begin(), it);
        dataExchangeM().giveData(namesFieldsUserCFDEMToExt[positionInRegister],"scalar-atom",
                                 particleDatFieldsUserCFDEMToExt[positionInRegister]
                                );
        Info << "giveData field with name '" << *it << "' at position: " << positionInRegister << endl;
    }
    if(verbose_) Info << "giveUSERdata done." << endl;
}

// * * *   write top level fields   * * * //
//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=

// * * * * * * * * * * * * * * * protected Member Functions  * * * * * * * * * * * * * //

void Foam::cfdemCloud::setNumberOfParticles(int nP)
{
    if(nP != numberOfParticles())
    {
        numberOfParticlesChanged_ = true;
        numberOfParticles_ = nP;
    }
}

void Foam::cfdemCloud::setNumberOfClumps(int nC)
{
    //Info << "Foam::cfdemCloud::setNumberOfClumps(int nC) ... do nothing" << endl;
}

void Foam::cfdemCloud::setPositionsCM(label n,double* pos)
{
    //Info << "Foam::cfdemCloud::setPositionsCM(int nC) ... do nothing" << endl;
}

void Foam::cfdemCloud::setCellIDsCM(label n,int* ID)
{
    //Info << "Foam::cfdemCloud::setCellIDsCM(int nC) ... do nothing" << endl;
}

void Foam::cfdemCloud::findCells()
{
    locateM().findCell(NULL,positions_,cellIDs_,numberOfParticles());
}

void Foam::cfdemCloud::setForces()
{
    resetArray(fluidVel_,numberOfParticles(),3);
    resetArray(impForces_,numberOfParticles(),3);
    resetArray(expForces_,numberOfParticles(),3);
    resetArray(DEMForces_,numberOfParticles(),3);
    resetArray(Cds_,numberOfParticles(),1);

    //=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
    //reset all USER-defined particle fields
    zeroizeParticleDatFieldsUserCFDEMToExt();
    //=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=

    for (int i=0;i<cfdemCloud::nrForceModels();i++) cfdemCloud::forceM(i).setForce();
}

void Foam::cfdemCloud::setVoidFraction()
{
    cfdemCloud::voidFractionM().setvoidFraction(NULL,voidfractions_,particleWeights_,particleVolumes_,particleV_);
}

void Foam::cfdemCloud::resetVoidFraction()
{
    cfdemCloud::voidFractionM().resetVoidFractions();
}

void Foam::cfdemCloud::setAlpha(volScalarField& alpha)
{
    alpha = cfdemCloud::voidFractionM().voidFractionInterp();
}

void Foam::cfdemCloud::setParticleForceField()
{
    averagingM().setVectorSum
    (
        forceM(0).impParticleForces(),
        impForces_,
        particleWeights_,
        NULL //mask
    );
    averagingM().setVectorSum
    (
        forceM(0).expParticleForces(),
        expForces_,
        particleWeights_,
        NULL //mask
    );
}

void Foam::cfdemCloud::setVectorAverages()
{
    if(verbose_) Info << "- setVectorAverage(Us,velocities_,weights_)" << endl;
    averagingM().setVectorAverage
    (
        averagingM().UsNext(),
        velocities_,
        particleWeights_,
        averagingM().UsWeightField(),
        NULL, //mask
        NULL,
        false
    );
    if(verbose_) Info << "setVectorAverage done." << endl;
}
// * * * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * * //
void Foam::cfdemCloud::checkCG(bool ok)
{
    if(!cgOK_) return;
    if(!ok) cgOK_ = ok;
}
void Foam::cfdemCloud::setPos(double**& pos)
{
    for(int index = 0;index <  numberOfParticles(); ++index)
    {
        for(int i=0;i<3;i++){
            positions_[index][i] = pos[index][i];
        }
    }
}
// * * * * * * * * * * * * * * * ACCESS  * * * * * * * * * * * * * //

label Foam::cfdemCloud::particleCell(int index)
{
    label cellI = cellIDs()[index][0];
    return cellI;
}

vector Foam::cfdemCloud::position(int index)
{
    vector pos;
    for(int i=0;i<3;i++) pos[i] = positions()[index][i];
    return pos;
}

vector Foam::cfdemCloud::velocity(int index)
{
    vector vel;
    for(int i=0;i<3;i++) vel[i] = velocities()[index][i];
    return vel;
}

vector Foam::cfdemCloud::fluidVel(int index)
{
    vector vel;
    for(int i=0;i<3;i++) vel[i] = fluidVels()[index][i];
    return vel;
}

const forceModel& Foam::cfdemCloud::forceM(int i)
{
    return forceModel_[i];
}

int Foam::cfdemCloud::nrForceModels()
{
    return forceModels_.size();
}

double** Foam::cfdemCloud::cellsPerParticle()
{
    return voidFractionModel_().cellsPerParticle();
}

scalar Foam::cfdemCloud::voidfraction(int index)
{
    return voidfractions()[index][0];
}

label Foam::cfdemCloud::liggghtsCommandModelIndex(word name)
{
    int index=-1;
    forAll(liggghtsCommandModelList_,i)
    {
        if(liggghtsCommand()[i]().name() == name)
        {
            index = i;
            break;
        }
    }
    return index;
}

std::vector< std::vector<double*> >* Foam::cfdemCloud::getVprobe()
{
 return probeModel_->getVprobe();
}

std::vector< std::vector<double> >* Foam::cfdemCloud::getSprobe()
{
 return probeModel_->getSprobe();
}

// * * * * * * * * * * * * * * * WRITE  * * * * * * * * * * * * * //

// * * *   write cfdemCloud internal data   * * * //

bool Foam::cfdemCloud::evolve
(
    volScalarField& alpha,
    volVectorField& Us,
    volVectorField& U
)
{
    numberOfParticlesChanged_ = false;
    arraysReallocated_=false;
    bool doCouple=false;

    if(!ignore())
    {
        if(!writeTimePassed_ && mesh_.time().outputTime()) writeTimePassed_=true;
        if (dataExchangeM().doCoupleNow())
        {
            Info << "\n Coupling..." << endl;
            dataExchangeM().couple(0);
            doCouple=true;

            // reset vol Fields
            clockM().start(16,"resetVolFields");
            if(verbose_)
            {
                Info << "couplingStep:" << dataExchangeM().couplingStep() 
                     << "\n- resetVolFields()" << endl;
            }
            averagingM().resetVectorAverage(averagingM().UsPrev(),averagingM().UsNext(),false);
            resetVoidFraction();
            averagingM().resetVectorAverage(forceM(0).impParticleForces(),forceM(0).impParticleForces(),true);
            averagingM().resetVectorAverage(forceM(0).expParticleForces(),forceM(0).expParticleForces(),true);
            averagingM().resetWeightFields();
            for (int i=0;i<momCoupleModels_.size(); i++)
                momCoupleM(i).resetMomSourceField();
            if(verbose_) Info << "resetVolFields done." << endl;
            clockM().stop("resetVolFields");

            if(verbose_) Info << "- getDEMdata()" << endl;
            clockM().start(17,"getDEMdata");
            getDEMdata();
            clockM().stop("getDEMdata");
            if(verbose_) Info << "- getDEMdata done." << endl;

            // search cellID of particles
            clockM().start(18,"findCell");
            if(verbose_) Info << "- findCell()" << endl;
            findCells();
            if(verbose_) Info << "findCell done." << endl;
            clockM().stop("findCell");

            // set void fraction field
            clockM().start(19,"setvoidFraction");
            if(verbose_) Info << "- setvoidFraction()" << endl;
            setVoidFraction();
            if(verbose_) Info << "setvoidFraction done." << endl;
            clockM().stop("setvoidFraction");

            // set average particles velocity field
            clockM().start(20,"setVectorAverage");
            setVectorAverages();


            //Smoothen "next" fields            
            smoothingM().dSmoothing();
            smoothingM().smoothen(voidFractionM().voidFractionNext());

            //only smoothen if we use implicit force coupling in cells void of particles
            //because we need unsmoothened Us field to detect cells for explicit 
            //force coupling
            if(!treatVoidCellsAsExplicitForce())
                smoothingM().smoothenReferenceField(averagingM().UsNext());
            
            clockM().stop("setVectorAverage");
        }
        
        //============================================
        //CHECK JUST TIME-INTERPOATE ALREADY SMOOTHENED VOIDFRACTIONNEXT AND UsNEXT FIELD 
        //      IMPLICIT FORCE CONTRIBUTION AND SOLVER USE EXACTLY THE SAME AVERAGED
        //      QUANTITIES AT THE GRID!
        Info << "\n timeStepFraction() = " << dataExchangeM().timeStepFraction() << endl;
        if( dataExchangeM().timeStepFraction() > 1.001 || dataExchangeM().timeStepFraction() < -0.001 )
        {
            FatalError << "cfdemCloud::dataExchangeM().timeStepFraction() = "<< dataExchangeM().timeStepFraction() <<" must be >1 or <0 : This might be due to the fact that you used a adjustable CFD time step. Please use a fixed CFD time step." << abort(FatalError);
        }
        clockM().start(24,"interpolateEulerFields");

        // update voidFractionField
        setAlpha(alpha);
        if(dataExchangeM().couplingStep() < 2)
        {
            alpha.oldTime() = alpha; // supress volume src
            alpha.oldTime().correctBoundaryConditions();
        }
        alpha.correctBoundaryConditions();

        // calc ddt(voidfraction)
        calcDdtVoidfraction(alpha,Us);

        // update mean particle velocity Field
        Us = averagingM().UsInterp();
        Us.correctBoundaryConditions();

        clockM().stop("interpolateEulerFields");
        //============================================

        if(doCouple)
        {
            // set particles forces
            clockM().start(21,"setForce");
            if(verbose_) Info << "- setForce(forces_)" << endl;
            setForces();
            if(verbose_) Info << "setForce done." << endl;
            calcMultiphaseTurbulence();
            if(verbose_) Info << "calcMultiphaseTurbulence done." << endl;
            clockM().stop("setForce");

            // get next force field
            clockM().start(22,"setParticleForceField");
            if(verbose_) Info << "- setParticleForceField()" << endl;
            setParticleForceField();
            if(verbose_) Info << "- setParticleForceField done." << endl;
            clockM().stop("setParticleForceField");

            // write DEM data
            if(verbose_) Info << " -giveDEMdata()" << endl;
            clockM().start(23,"giveDEMdata");
            giveDEMdata();
            clockM().stop("giveDEMdata");

            dataExchangeM().couple(1);
        }//end dataExchangeM().couple()


        if(verbose_){
            #include "debugInfo.H"
        }

        clockM().start(25,"dumpDEMdata");
        // do particle IO
        IOM().dumpDEMdata();
        clockM().stop("dumpDEMdata");

    }//end ignore
    return doCouple;
}

bool Foam::cfdemCloud::reAllocArrays() const
{
    if(numberOfParticlesChanged_ && !arraysReallocated_)
    {
        // get arrays of new length
        dataExchangeM().allocateArray(positions_,0.,3);
        dataExchangeM().allocateArray(velocities_,0.,3);
        dataExchangeM().allocateArray(fluidVel_,0.,3);
        dataExchangeM().allocateArray(fAcc_,0.,3);
        dataExchangeM().allocateArray(impForces_,0.,3);
        dataExchangeM().allocateArray(expForces_,0.,3);
        dataExchangeM().allocateArray(DEMForces_,0.,3);
        dataExchangeM().allocateArray(Cds_,0.,1);
        dataExchangeM().allocateArray(radii_,0.,1);
        dataExchangeM().allocateArray(voidfractions_,1.,voidFractionM().maxCellsPerParticle());
        dataExchangeM().allocateArray(cellIDs_,-1.,voidFractionM().maxCellsPerParticle());
        dataExchangeM().allocateArray(particleWeights_,0.,voidFractionM().maxCellsPerParticle());
        dataExchangeM().allocateArray(particleVolumes_,0.,voidFractionM().maxCellsPerParticle());
        dataExchangeM().allocateArray(particleV_,0.,1);
        
        //=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
        if(namesFieldsUserCFDEMToExt.size()!=particleDatFieldsUserCFDEMToExt.size())
            allocateParticleDatFieldsUserCFDEMToExt();
        else
            reAllocateParticleDatFieldsUserCFDEMToExt();
        //=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=

        arraysReallocated_ = true;
        return true;
    }
    return false;
}

bool Foam::cfdemCloud::reAllocArrays(int nP, bool forceRealloc) const
{
    if( (numberOfParticlesChanged_ && !arraysReallocated_) || forceRealloc)
    {
        // get arrays of new length
        dataExchangeM().allocateArray(positions_,0.,3,nP);
        dataExchangeM().allocateArray(velocities_,0.,3,nP);
        dataExchangeM().allocateArray(fluidVel_,0.,3,nP);
        dataExchangeM().allocateArray(impForces_,0.,3,nP);
        dataExchangeM().allocateArray(expForces_,0.,3,nP);
        dataExchangeM().allocateArray(DEMForces_,0.,3,nP);
        dataExchangeM().allocateArray(Cds_,0.,1,nP);
        dataExchangeM().allocateArray(radii_,0.,1,nP);
        dataExchangeM().allocateArray(voidfractions_,1.,voidFractionM().maxCellsPerParticle(),nP);
        dataExchangeM().allocateArray(cellIDs_,0.,voidFractionM().maxCellsPerParticle(),nP);
        dataExchangeM().allocateArray(particleWeights_,0.,voidFractionM().maxCellsPerParticle(),nP);
        dataExchangeM().allocateArray(particleVolumes_,0.,voidFractionM().maxCellsPerParticle(),nP);

        //=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
        if(namesFieldsUserCFDEMToExt.size()!=particleDatFieldsUserCFDEMToExt.size())
            allocateParticleDatFieldsUserCFDEMToExt();
        else
            reAllocateParticleDatFieldsUserCFDEMToExt();
        //=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
        arraysReallocated_ = true;
        return true;
    }
    return false;
}

tmp<fvVectorMatrix> cfdemCloud::divVoidfractionTau(volVectorField& U,volScalarField& voidfraction) const
{
    return
    (
      - fvm::laplacian(voidfractionNuEff(voidfraction), U)
      - fvc::div(voidfractionNuEff(voidfraction)*dev2(fvc::grad(U)().T()))
    );
}

tmp<volScalarField> cfdemCloud::ddtVoidfraction() const
{
    if (useDDTvoidfraction_==word("off"))
    {
        return tmp<volScalarField> (ddtVoidfraction_ * 0.);
        if(verbose_)
            Info << "suppressing ddt(voidfraction)" << endl;
    }
    return tmp<volScalarField> (ddtVoidfraction_ * 1.) ;
}

void cfdemCloud::calcDdtVoidfraction(volScalarField& voidfraction, volVectorField& Us) const
{
    if (useDDTvoidfraction_==word("a"))
    {
        // Calculation of new ddtVoidfraction by using divergence of particle fluxes instead
        ddtVoidfraction_=fvc::div(Us*(1.-voidfraction));
    }else // "b" or "off"
    {
        ddtVoidfraction_ = fvc::ddt(voidfraction);
    }
}

//****************************************
void cfdemCloud::calcMultiphaseTurbulence()
{
    //Temporary field for collecting the sources
    volScalarField tmpSource
    (   IOobject
        (
            "tmpSource",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        turbulenceMultiphase_*0.0
    );

    //reset sources, and compute the sources due to the cloud
    //Will accumulate all sources for all force models
    for (int iModel=0; iModel<nrForceModels(); iModel++) //
    {
        #ifdef compre
            forceM(iModel).multiphaseTurbulence(tmpSource, true);
        #else
            forceM(iModel).multiphaseTurbulence(tmpSource, false);
        #endif
        if(iModel==0)
            turbulenceMultiphase_    = tmpSource;
        else
            turbulenceMultiphase_   += tmpSource;
    }

}

/*tmp<fvVectorMatrix> cfdemCloud::ddtVoidfractionU(volVectorField& U,volScalarField& voidfraction) const
{
    if (dataExchangeM().couplingStep() <= 2) return fvm::ddt(U);
    
    return fvm::ddt(voidfraction,U);
}*/

tmp<volScalarField> cfdemCloud::voidfractionNuEff(volScalarField& voidfraction) const
{
    if (modelType_=="B" || modelType_=="Bfull")
    {
        return tmp<volScalarField>
        (
            #ifdef compre
                new volScalarField("viscousTerm", (  turbulence_.mut()  
                                                   + turbulence_.mu() 
                                                   + turbulenceMultiphase_
                                                  )
                                  )
            #else
                new volScalarField("viscousTerm", (  turbulence_.nut() 
                                                   + turbulence_.nu()
                                                   + turbulenceMultiphase_
                                                  )
                                  )
            #endif
        );
    }
    else
    {
        return tmp<volScalarField>
        (
            #ifdef compre
                new volScalarField("viscousTerm", voidfraction*(  turbulence_.mut() 
                                                                + turbulence_.mu()
                                                                + turbulenceMultiphase_
                                                               )
                                  )
            #else
                new volScalarField("viscousTerm", voidfraction*(  turbulence_.nut() 
                                                                + turbulence_.nu()
                                                                + turbulenceMultiphase_
                                                               )
                                  )
            #endif
        );
    }
}

void cfdemCloud::resetArray(double**& array,int length,int width,double resetVal)
{
    for(int index = 0;index < length; ++index){
        for(int i=0;i<width;i++){
            array[index][i] = resetVal;
        }
    }
}

double **cfdemCloud::dragPrev()
{
    return dragPrev_;
}
// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //
//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
void cfdemCloud::registerNamesFieldsUserCFDEMToExt(word fieldToRegister, int& positionInRegister)
{
    //check if field is available
    Info << "cfdemCloud is registering field '" << fieldToRegister <<"'" << endl;
    std::vector<word>::iterator it;
    it = std::find(namesFieldsUserCFDEMToExt.begin(), namesFieldsUserCFDEMToExt.end(), fieldToRegister);
    if ( it != namesFieldsUserCFDEMToExt.end() )
    {
        positionInRegister = std::distance(namesFieldsUserCFDEMToExt.begin(), it);
        Info << "cfdemCloud found field '" << fieldToRegister << "' at position: " << positionInRegister << endl;
    }
    else
    {
        //if not, add to list of names
        Info << "cfdemCloud COULD NOT find field '" << fieldToRegister <<"', will push to end." << endl;
        namesFieldsUserCFDEMToExt.push_back(fieldToRegister);
        positionInRegister = namesFieldsUserCFDEMToExt.size()-1;
    }
}

//****************************************
bool cfdemCloud::checkAndregisterNamesFieldsUserCFDEMToExt(const wordList names, std::vector<int> & positionInRegister)
{
    bool validFieldName=false;
    forAll(names,i)    {
        int tempPosition=-1; //by default use -1 to indicate invalid field
        if(names[i]!="none")
        {
            validFieldName = true;
            registerNamesFieldsUserCFDEMToExt(names[i], tempPosition);
        }
        positionInRegister.push_back(tempPosition);
    }
    return validFieldName;
}

//****************************************
void cfdemCloud::allocateParticleDatFieldsUserCFDEMToExt() const
{
    if(particleDatFieldsUserCFDEMToExt.size()>0)
        FatalError << "cfdemCloud::allocateParticleDatFieldsUserCFDEMToExt(): you are attempting to allocate fields in a container that already contains elements. This is not allowed, please clear container." << abort(FatalError);
    //Go through list and allocate
    for(std::vector<word>::const_iterator it = namesFieldsUserCFDEMToExt.begin(); it != namesFieldsUserCFDEMToExt.end(); ++it) 
    {
        Info << "allocating field with name '" << *it << "'" << endl;
        particleDatFieldsUserCFDEMToExt.push_back(NULL); //Must be NULL, otherwise this might confuse external code
        dataExchangeM().allocateArray(particleDatFieldsUserCFDEMToExt.back(),0.0,1);
    }
}

//****************************************
void cfdemCloud::reAllocateParticleDatFieldsUserCFDEMToExt() const
{
    //Go through list and allocate
    for(std::vector<word>::iterator it = namesFieldsUserCFDEMToExt.begin(); it != namesFieldsUserCFDEMToExt.end(); ++it) 
    {
        int positionInRegister = std::distance(namesFieldsUserCFDEMToExt.begin(), it);
        if(verbose_)
            Info << "reAllocating field with name '" << *it << "' at position: " << positionInRegister << endl;
        dataExchangeM().allocateArray(particleDatFieldsUserCFDEMToExt[positionInRegister],0.0,1);
    }
}

//****************************************
void cfdemCloud::zeroizeParticleDatFieldsUserCFDEMToExt()
{
    //Go through list and set zero
    for(std::vector<word>::iterator it = namesFieldsUserCFDEMToExt.begin(); it != namesFieldsUserCFDEMToExt.end(); ++it) 
    {
        int positionInRegister = std::distance(namesFieldsUserCFDEMToExt.begin(), it);
        if(verbose_)
            Info << "Zeroizing field with name '" << *it << "' at position: " << positionInRegister << endl;
        resetArray(particleDatFieldsUserCFDEMToExt[positionInRegister],numberOfParticles(),1);
    }

}

//****************************************
void cfdemCloud::accessParticleDatFieldsUserCFDEMToExt(word fieldToAccess, double **& fieldData)
{
    //set pointer to correct location in the memory
    if(verbose_)
        Info << "cfdemCloud is accessing field '" << fieldToAccess << "'" << endl;
    std::vector<word>::iterator it;
    it = std::find(namesFieldsUserCFDEMToExt.begin(), namesFieldsUserCFDEMToExt.end(), fieldToAccess);
    if ( it != namesFieldsUserCFDEMToExt.end() )
    {
        int positionInRegister = std::distance(namesFieldsUserCFDEMToExt.begin(), it);
        if(verbose_)
            Info << "cfdemCloud found field '" << fieldToAccess << "' at position: " << positionInRegister << endl;
        fieldData = particleDatFieldsUserCFDEMToExt[positionInRegister];
    }
    else
            FatalError << "field " << fieldToAccess 
                       << " not found."
                       << abort(FatalError);
}
//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "cfdemCloudIO.C"

// ************************************************************************* //
