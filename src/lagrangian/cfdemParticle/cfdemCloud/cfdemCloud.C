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
    cgParticleSpecific_(false),
    cgTypeSpecificDifferent_(false),
    allowAdjustTimeStep_(couplingProperties_.lookupOrDefault<Switch>("allowAdjustTimeStep", false)),
    solveFlow_(true),
    solveScalarTransport_(true),
    verbose_(couplingProperties_.lookupOrDefault<Switch>("verbose", false)),
    expCorrDeltaUError_(false),
    debug_(false),
    allowCFDsubTimestep_(true),
    ignore_(false),
    writeTimePassed_(false),
    resetWriteTimePassed_(false),
    modelType_(couplingProperties_.lookup("modelType")),
    impForces_(NULL),
    expForces_(NULL),
    voidfractions_(NULL),
    cellIDs_(NULL),
    particleWeights_(NULL),
    particleVolumes_(NULL),
    particleV_(NULL),
    dragPrev_(NULL),
    numberOfParticles_(-1),
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
    checkPeriodicCells_(couplingProperties_.lookupOrDefault<Switch>("checkPeriodicCells",false)),
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
    ),
    idRadius_(-1),
    idPos_(-1),
    idVel_(-1),
    idFacc_(-1),
    idPartTypes_(-1),
    idDragExp_(-1),
    idKsl_(-1),
    idKslExtra_(-1),
    idUf_(-1),
    idTorqueExp_(-1),
    idKslRotation_(-1),
    idPullRotation_(-1),
    idPullOrientation_(-1),
    idPullOrientation1_(-1),
    idPullShape_(-1),
    idDragExpCM_(-1),
    idKslCM_(-1),
    idKslExtraCM_(-1),
    idUfCM_(-1),
    idTorqueExpCM_(-1),
    idKslRotationCM_(-1),
    idFhydro_(-1),
    idVisc_(-1),
    idBlockiness_(-1),
    idArea_(-1),
    idVol_(-1),
    idQuat_(-1),
    idK_(-1),
    idEpsilon_(-1),
    idParticleCG_(-1),
    idMass_(-1),
    idDensity_(-1),
    idType_(-1),
    idConvectiveHeatFlux_(-1),
    idTemp_(-1)
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
    cgTypeSpecificDifferent_ = false;
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
        registryM().addProperty(forceModel_[i]().type()+"_index",i);
    }

    IOM().allocFieldsToDEM();

    registryM().addProperty(voidFractionM().type()+"_index", 0);
    registryM().addProperty(dataExchangeM().type()+"_index", 0);

    // Note: this check for getProperty("multisphere") does not work as the cloud is
    // built before cloudMS which sets the property
    //if (nrForceModels()<SMALL && registryM().getProperty("multisphere")<1)
    //    FatalError  << "Please use at least one forceModel ! "
    //                << "(e.g. noDrag) \n"
    //                << abort(FatalError);

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
        registryM().addProperty(liggghtsCommand_[i]().type()+"_index",i);
    }

    #include "sanityChecks/level0_Cloud.H"
}

// * * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * //
Foam::cfdemCloud::~cfdemCloud()
{
    clockM().evalPar();
    clockM().normHist();
    dataExchangeM().destroy(impForces_,3);
    dataExchangeM().destroy(expForces_,3);
    dataExchangeM().destroy(voidfractions_,1);
    dataExchangeM().destroy(cellIDs_,1);
    dataExchangeM().destroy(particleWeights_,1);
    dataExchangeM().destroy(particleVolumes_,1);
    dataExchangeM().destroy(particleV_,1);

    //=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
    int iUser=0;
    for( std::vector<double**>::iterator
         it  = fieldsToDEM.begin();
         it != fieldsToDEM.end();
       ++it)
    {
        Info << "cfdemCloud destroys UserCFDEM data: " << namesfieldsToDEM[iUser++] << endl;
        int id = std::distance(fieldsToDEM.begin(), it);
        if(typesfieldsToDEM[id]=="scalar-atom" || typesfieldsToDEM[id]=="scalar-multisphere")
            dataExchangeM().destroy((*it),1);
        else if(typesfieldsToDEM[id]=="vector-atom" || typesfieldsToDEM[id]=="vector-multisphere")
            dataExchangeM().destroy((*it),3);
        else if(typesfieldsToDEM[id]=="vector2D-atom" || typesfieldsToDEM[id]=="vector2D-multisphere")
            dataExchangeM().destroy((*it),2);
        else if(typesfieldsToDEM[id]=="quaternion-atom" || typesfieldsToDEM[id]=="quaternion-multisphere")
            dataExchangeM().destroy((*it),4);
        else
            FatalError << "cfdemCloud::~cfdemCloud(): unknown data type "<<  typesfieldsToDEM[id] << " to destroy." << abort(FatalError);
    }
    //=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
}
// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
void Foam::cfdemCloud::getDEMdata()
{
    if(verbose_) Info << "Foam::cfdemCloud::getDEMdata()" << endl;
    for(std::vector<word>::iterator it = namesfieldsToDEM.begin(); it != namesfieldsToDEM.end(); ++it)
    {
        int id = std::distance(namesfieldsToDEM.begin(), it);
        if(pullfieldsToDEM[id])
        {
            Info << "get field with name '" << *it << "' at position: " << id << endl;
            if(namesfieldsToDEM[id]=="shapetype" || namesfieldsToDEM[id]=="type" || namesfieldsToDEM[id]=="body" || namesfieldsToDEM[id]=="id") //type is of type int and we do not yet have fieldsToDEM for int
            {
                int** h=NULL;
                dataExchangeM().allocateArray(h,0.,1);
                dataExchangeM().getData(namesfieldsToDEM[id],typesfieldsToDEM[id],h);
                for(int i=0;i<numberOfParticles_;i++)
                {
                    fieldsToDEM[id][i][0]=h[i][0];
                    if(namesfieldsToDEM[id]=="shapetype") fieldsToDEM[id][i][0] += 1; // conversion to clumpType convention
                }
            }
            else if(namesfieldsToDEM[id]=="clumptype" || namesfieldsToDEM[id]=="nrigid" || namesfieldsToDEM[id]=="id_multisphere")
            {
                int** h=NULL;
                dataExchangeM().allocateArray(h,0.,1);
                dataExchangeM().getData(namesfieldsToDEM[id],typesfieldsToDEM[id],h);
                for(int i=0;i<numberOfClumps();i++)
                {
                    fieldsToDEM[id][i][0]=h[i][0];
                }
            }
            else
                dataExchangeM().getData(namesfieldsToDEM[id],typesfieldsToDEM[id],fieldsToDEM[id]);
        }
    }
    if(verbose_) Info << "Foam::cfdemCloud::getDEMdata() - done." << endl;
}

void Foam::cfdemCloud::giveDEMdata()
{
    for(std::vector<word>::iterator it = namesfieldsToDEM.begin(); it != namesfieldsToDEM.end(); ++it)
    {
        int id = std::distance(namesfieldsToDEM.begin(), it);
        if(!pullfieldsToDEM[id])
        {
            Info << "giveData field with name '" << *it << "' at position: " << id << endl;
            dataExchangeM().giveData(namesfieldsToDEM[id],typesfieldsToDEM[id],fieldsToDEM[id]);
        }
    }
    if(verbose_) Info << "giveDEMdata done." << endl;
}

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
    locateM().findCell(NULL,fieldsToDEM[idPos()],cellIDs_,numberOfParticles());
}

void Foam::cfdemCloud::setForces()
{
    resetArray(impForces_,numberOfParticles(),3);
    resetArray(expForces_,numberOfParticles(),3);

    cfdemCloud::zeroizeFieldsToDEM();

    for (int i=0;i<cfdemCloud::nrForceModels();i++)
        cfdemCloud::forceM(i).setForce();
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
    if(verbose_) Info << "- setVectorAverage(Us,fieldsToDEM[idVel()],weights_)" << endl;
    averagingM().setVectorAverage
    (
        averagingM().UsNext(),
        fieldsToDEM[idVel()],
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
            fieldsToDEM[idPos()][index][i] = pos[index][i];
        }
    }
}
// * * * * * * * * * * * * * * * ACCESS  * * * * * * * * * * * * * //

label Foam::cfdemCloud::particleCell(int index)
{
    label cellI = cellIDs()[index][0];
    return cellI;
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

        if(resetWriteTimePassed_)
        {
            writeTimePassed_=false;
            resetWriteTimePassed_=false;
        }

    }//end ignore
    return doCouple;
}

bool Foam::cfdemCloud::reAllocArrays() const
{
    if(verbose_) Info << "cfdemCloud::reAllocArrays()" << endl;
    if(numberOfParticlesChanged_ && !arraysReallocated_)
    {
        // get arrays of new length
        dataExchangeM().allocateArray(impForces_,0.,3);
        dataExchangeM().allocateArray(expForces_,0.,3);
        dataExchangeM().allocateArray(voidfractions_,1.,voidFractionM().maxCellsPerParticle());
        dataExchangeM().allocateArray(cellIDs_,-1.,voidFractionM().maxCellsPerParticle());
        dataExchangeM().allocateArray(particleWeights_,0.,voidFractionM().maxCellsPerParticle());
        dataExchangeM().allocateArray(particleVolumes_,0.,voidFractionM().maxCellsPerParticle());
        dataExchangeM().allocateArray(particleV_,0.,1);

        //=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
        if(namesfieldsToDEM.size()!=fieldsToDEM.size())
            cfdemCloud::allocateFieldsToDEM();
        else
            cfdemCloud::reAllocateFieldsToDEM();
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
        dataExchangeM().allocateArray(impForces_,0.,3,nP);
        dataExchangeM().allocateArray(expForces_,0.,3,nP);
        dataExchangeM().allocateArray(voidfractions_,1.,voidFractionM().maxCellsPerParticle(),nP);
        dataExchangeM().allocateArray(cellIDs_,0.,voidFractionM().maxCellsPerParticle(),nP);
        dataExchangeM().allocateArray(particleWeights_,0.,voidFractionM().maxCellsPerParticle(),nP);
        dataExchangeM().allocateArray(particleVolumes_,0.,voidFractionM().maxCellsPerParticle(),nP);

        //=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
        if(namesfieldsToDEM.size()!=fieldsToDEM.size())
            cfdemCloud::allocateFieldsToDEM();
        else
            cfdemCloud::reAllocateFieldsToDEM();
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

tmp<fvVectorMatrix> cfdemCloud::divVoidfractionTau(volVectorField& U,volScalarField& voidfraction,volScalarField& rho) const
{
    return
    (
      - fvm::laplacian(rho*voidfractionNuEff(voidfraction), U)
      - fvc::div(rho*voidfractionNuEff(voidfraction)*dev2(fvc::grad(U)().T()))
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
void cfdemCloud::registerFieldsToDEM(word name, word type, int& id, bool pull)
{
    //check if field is available
    Info << "cfdemCloud is registering field '" << name <<"', and type '" << type << "' ..." << endl;
    std::vector<word>::iterator it;
    it = std::find(namesfieldsToDEM.begin(), namesfieldsToDEM.end(), name);
    if ( it != namesfieldsToDEM.end() )
    {
        id = std::distance(namesfieldsToDEM.begin(), it);
        Info << "cfdemCloud already found field '" << name << "' at position: " << id << endl;
    }
    else
    {
        Info << "cfdemCloud adds field '" << name <<"' at end of list." << endl;
        namesfieldsToDEM.push_back(name);
        typesfieldsToDEM.push_back(type);
        pullfieldsToDEM.push_back(pull);
        id = namesfieldsToDEM.size()-1;
    }
}

//****************************************
bool cfdemCloud::checkAndRegisterFieldsToDEM(const wordList names, const word type, std::vector<int> & id)
{
    bool validFieldName=false;
    forAll(names,i)    {
        int tempPosition=-1; //by default use -1 to indicate invalid field
        if(names[i]!="none")
        {
            validFieldName = true;
            registerFieldsToDEM(names[i], type, tempPosition);
        }
        id.push_back(tempPosition);
    }
    return validFieldName;
}

//****************************************
void cfdemCloud::allocateFieldsToDEM(word shapeType) const
{
    if(verbose_) Info << "cfdemCloud::allocateFieldsToDEM()" << endl;
    if(fieldsToDEM.size()>0)
        FatalError << "cfdemCloud::allocateFieldsToDEM(): you are attempting to allocate fields in a container that already contains elements. This is not allowed, please clear container." << abort(FatalError);
    //Go through list and allocate
    for(std::vector<word>::iterator it = namesfieldsToDEM.begin(); it != namesfieldsToDEM.end(); ++it)
    {
        int id = std::distance(namesfieldsToDEM.begin(), it);

        fieldsToDEM.push_back(NULL); //Must be NULL, otherwise this might confuse external code

        word type(typesfieldsToDEM[id]);
        int len(-1);
        if(type=="scalar-" + shapeType)
            len=1;
        else if(type=="vector-" + shapeType)
            len=3;
        else if(type=="vector2D-" + shapeType)
            len=2;
        else if(type=="quaternion-" + shapeType)
            len=4;
//         else
//             FatalError << "cfdemCloud::allocateFieldsToDEM -- unknown array shape: "
//                 << type << abort(FatalError);

        if(len>0)
        {
            if(verbose_) Info << "cfdemCloud::allocateFieldsToDEM() allocating field with name '" << *it << "'" << endl;
            dataExchangeM().allocateArray(fieldsToDEM.back(),0.0,len);
        }
    }
}

//****************************************
void cfdemCloud::reAllocateFieldsToDEM(word shapeType) const
{
    if(verbose_) Info << "cfdemCloud::reAllocateFieldsToDEM()" << endl;
    //Go through list and allocate
    for(std::vector<word>::iterator it = namesfieldsToDEM.begin(); it != namesfieldsToDEM.end(); ++it)
    {
        int id = std::distance(namesfieldsToDEM.begin(), it);

        word type(typesfieldsToDEM[id]);
        int len(-1);
        if(type=="scalar-" + shapeType)
            len=1;
        else if(type=="vector-" + shapeType)
            len=3;
        else if(type=="vector2D-" + shapeType)
            len=2;
        else if(type=="quaternion-" + shapeType)
            len=4;
//         else
//             FatalError << "cfdemCloud::reAllocateFieldsToDEM -- unknown array shape: "
//                 << type << abort(FatalError);


        if(len>0)
        {
            if(verbose_) Info << "cfdemCloud::reAllocateFieldsToDEM() reAllocating field with name '" << *it << "' at position: " << id << endl;
            dataExchangeM().allocateArray(fieldsToDEM[id],0.0,len);
        }
    }
}

//****************************************
void cfdemCloud::zeroizeFieldsToDEM(word shapeType)
{
    //Go through list and set zero
    for(std::vector<word>::iterator it = namesfieldsToDEM.begin(); it != namesfieldsToDEM.end(); ++it)
    {
        int id = std::distance(namesfieldsToDEM.begin(), it);
        if(!pullfieldsToDEM[id])
        {
            if(typesfieldsToDEM[id]=="scalar-" + shapeType)
            {
                //if(verbose_)
                Info << "Zeroizing field with name '" << *it
                     << "' at position: " << id
                     << " with type " << typesfieldsToDEM[id] << endl;
                resetArray(fieldsToDEM[id],numberOfParticles(),1);
            }
            else if(typesfieldsToDEM[id]=="vector-" + shapeType)
            {
                //if(verbose_)
                Info << "Zeroizing field with name '" << *it
                     << "' at position: " << id
                     << " with type " << typesfieldsToDEM[id] << endl;
                resetArray(fieldsToDEM[id],numberOfParticles(),3);
            }
        }
    }
}

//****************************************
void cfdemCloud::accessFieldsToDEM(word name, double **& field)
{
    //set pointer to correct location in the memory
    if(verbose_)
        Info << "cfdemCloud is accessing field '" << name << "'" << endl;
    std::vector<word>::iterator it;
    it = std::find(namesfieldsToDEM.begin(), namesfieldsToDEM.end(), name);
    if ( it != namesfieldsToDEM.end() )
    {
        int id = std::distance(namesfieldsToDEM.begin(), it);
        if(verbose_)
            Info << "cfdemCloud found field '" << name << "' at position: " << id << endl;
        field = fieldsToDEM[id];
    }
    else
            FatalError << "field " << name
                       << " not found."
                       << abort(FatalError);
}
//****************************************
int cfdemCloud::existsFieldsToDEM(word name)
{
    //set pointer to correct location in the memory
    if(verbose_)
        Info << "cfdemCloud is testing existence field '" << name << "'" << endl;

    int id(-1);
    std::vector<word>::iterator it;
    it = std::find(namesfieldsToDEM.begin(), namesfieldsToDEM.end(), name);
    if ( it != namesfieldsToDEM.end() )
    {
        id = std::distance(namesfieldsToDEM.begin(), it);
        if(verbose_)
        {
            if(pullfieldsToDEM[id]) Info << "  cfdemCloud found the pull field '" << name << "' at position: " << id << endl;
            else  Info << "  cfdemCloud found the push field '" << name << "' at position: " << id << endl;
        }
    }
    else
    {
        if(verbose_)
            Info << "  cfdemCloud could not fine field '" << name << "'" << endl;
    }
    return id;
}
//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=

// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "cfdemCloudIO.C"

// ************************************************************************* //
