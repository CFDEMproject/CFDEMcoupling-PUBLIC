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
#include "forceModel.H"
#include "locateModel.H"
#include "momCoupleModel.H"
#include "meshMotionModel.H"
#include "voidFractionModel.H"
#include "dataExchangeModel.H"
#include "IOModel.H"
#include "averagingModel.H"
#include "clockModel.H"
#include "smoothingModel.H"
#include "liggghtsCommandModel.H"

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
    verbose_(false),
    ignore_(false),
    modelType_(couplingProperties_.lookup("modelType")),
    positions_(NULL),
    velocities_(NULL),
    fluidVel_(NULL),
    impForces_(NULL),
    expForces_(NULL),
    DEMForces_(NULL),
    Cds_(NULL),
    radii_(NULL),
    voidfractions_(NULL),
    cellIDs_(NULL),
    particleWeights_(NULL),
    particleVolumes_(NULL),
    numberOfParticles_(0),
    numberOfParticlesChanged_(false),
    arraysReallocated_(false),
    forceModels_(couplingProperties_.lookup("forceModels")),
    momCoupleModels_(couplingProperties_.lookup("momCoupleModels")),
    liggghtsCommandModelList_(liggghtsCommandDict_.lookup("liggghtsCommandModels")),
    turbulenceModelType_(couplingProperties_.lookup("turbulenceModelType")),
    cgOK_(true),
    impDEMdrag_(false),
    useDDTvoidfraction_(false),
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
    turbulence_
    (
        #if defined(version21) || defined(version16ext)
            #ifdef comp
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

    Info << "If BC are important, please provide volScalarFields -imp/expParticleForces-" << endl;

    if (couplingProperties_.found("verbose")) verbose_=true;
    if (couplingProperties_.found("ignore")) ignore_=true;
    if (turbulenceModelType_=="LESProperties")
        Info << "WARNING - LES functionality not yet tested!" << endl;

    if (couplingProperties_.found("useDDTvoidfraction"))
        useDDTvoidfraction_=true;
    else        
        Info << "ignoring ddt(voidfraction)" << endl;

    forceModel_ = new autoPtr<forceModel>[nrForceModels()];
    for (int i=0;i<nrForceModels();i++)
    {
        forceModel_[i] = forceModel::New
        (
            couplingProperties_,
            *this,
            forceModels_[i]
        );
    }

    momCoupleModel_ = new autoPtr<momCoupleModel>[momCoupleModels_.size()];
    for (int i=0;i<momCoupleModels_.size();i++)
    {
        momCoupleModel_[i] = momCoupleModel::New
        (
            couplingProperties_,
            *this,
            momCoupleModels_[i]
        );
    }

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

    dataExchangeM().setCG();
    if (!cgOK_ && forceM(0).cg() > 1) FatalError<< "at least one of your models is not fit for cg !!!"<< abort(FatalError); 
}

// * * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * //
Foam::cfdemCloud::~cfdemCloud()
{
    clockM().evalPar();
    clockM().normHist();
    dataExchangeM().destroy(positions_,3);
    dataExchangeM().destroy(velocities_,3);
    dataExchangeM().destroy(fluidVel_,3);
    dataExchangeM().destroy(impForces_,3);
    dataExchangeM().destroy(expForces_,3);
    dataExchangeM().destroy(DEMForces_,3);
    dataExchangeM().destroy(Cds_,1);
    dataExchangeM().destroy(radii_,1);
    dataExchangeM().destroy(voidfractions_,1);
    dataExchangeM().destroy(cellIDs_,1);
    dataExchangeM().destroy(particleWeights_,1);
    dataExchangeM().destroy(particleVolumes_,1);
}
// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
void Foam::cfdemCloud::getDEMdata()
{
    dataExchangeM().getData("radius","scalar-atom",radii_);
    dataExchangeM().getData("x","vector-atom",positions_);
    dataExchangeM().getData("v","vector-atom",velocities_);
}

void Foam::cfdemCloud::giveDEMdata()
{
    if(forceM(0).coupleForce())
    {
        dataExchangeM().giveData("dragforce","vector-atom",DEMForces_);

        if(impDEMdrag_)
        {
            dataExchangeM().giveData("Ksl","scalar-atom",Cds_);
            dataExchangeM().giveData("uf","vector-atom",fluidVel_);
        }
    }
    if(verbose_) Info << "giveDEMdata done." << endl;
}

// * * *   write top level fields   * * * //

// * * * * * * * * * * * * * * * protected Member Functions  * * * * * * * * * * * * * //

void Foam::cfdemCloud::setNumberOfParticles(int nP)
{
    if(nP != numberOfParticles())
    {
        numberOfParticlesChanged_ = true;
        numberOfParticles_ = nP;
    }
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
    for (int i=0;i<cfdemCloud::nrForceModels();i++) cfdemCloud::forceM(i).setForce();
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

double Foam::cfdemCloud::d(int index)
{
    return 2*radii()[index][0];
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

scalar Foam::cfdemCloud::radius(int index)
{
    scalar r = radii()[index][0];
    return r;
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
        if (dataExchangeM().couple())
        {
            Info << "\n Coupling..." << endl;
            doCouple=true;

            // reset vol Fields
            clockM().start(16,"resetVolFields");
            if(verbose_){
                Info << "couplingStep:" << dataExchangeM().couplingStep() 
                     << "\n- resetVolFields()" << endl;
            }
            averagingM().resetVectorAverage(averagingM().UsPrev(),averagingM().UsNext());
            voidFractionM().resetVoidFractions();
            averagingM().resetVectorAverage(forceM(0).impParticleForces(),forceM(0).impParticleForces(),true);
            averagingM().resetVectorAverage(forceM(0).expParticleForces(),forceM(0).expParticleForces(),true);
            averagingM().resetWeightFields();
            momCoupleM(0).resetMomSourceField();
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
            voidFractionM().setvoidFraction(NULL,voidfractions_,particleWeights_,particleVolumes_);
            if(verbose_) Info << "setvoidFraction done." << endl;
            clockM().stop("setvoidFraction");

            // set particles velocity field
            clockM().start(20,"setVectorAverage");
            if(verbose_) Info << "- setVectorAverage(Us,velocities_,weights_)" << endl;
            averagingM().setVectorAverage
            (
                averagingM().UsNext(),
                velocities_,
                particleWeights_,
                averagingM().UsWeightField(),
                NULL //mask
            );
            if(verbose_) Info << "setVectorAverage done." << endl;
            clockM().stop("setVectorAverage");

            // set particles forces
            clockM().start(21,"setForce");
            if(verbose_) Info << "- setForce(forces_)" << endl;
            setForces();
            if(verbose_) Info << "setForce done." << endl;
            clockM().stop("setForce");

            // get next force field
            clockM().start(22,"setParticleForceField");
            if(verbose_) Info << "- setParticleForceField()" << endl;
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
            if(verbose_) Info << "- setParticleForceField done." << endl;
            clockM().stop("setParticleForceField");

            // write DEM data
            if(verbose_) Info << " -giveDEMdata()" << endl;
            clockM().start(23,"giveDEMdata");
            giveDEMdata();
            clockM().stop("giveDEMdata");
        }//end dataExchangeM().couple()
        Info << "\n timeStepFraction() = " << dataExchangeM().timeStepFraction() << endl;

        clockM().start(24,"interpolateEulerFields");
        // update voidFractionField
        alpha.oldTime().internalField() = voidFractionM().voidFractionInterp();
        
        // smoothing exchange field
        smoothingM().smoothen(alpha);
        alpha.correctBoundaryConditions();

        // calc ddt(voidfraction)
        if (doCouple) calcDdtVoidfraction(voidFractionM().voidFractionNext());
        //calcDdtVoidfraction(alpha); // alternative with scale=1! (does not see change in alpha?)

        // update particle velocity Field
        Us.internalField() = averagingM().UsInterp();
        Us.correctBoundaryConditions();
        clockM().stop("interpolateEulerFields");

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
        dataExchangeM().allocateArray(impForces_,0.,3);
        dataExchangeM().allocateArray(expForces_,0.,3);
        dataExchangeM().allocateArray(DEMForces_,0.,3);
        dataExchangeM().allocateArray(Cds_,0.,1);
        dataExchangeM().allocateArray(radii_,0.,1);
        dataExchangeM().allocateArray(voidfractions_,1.,voidFractionM().maxCellsPerParticle());
        dataExchangeM().allocateArray(cellIDs_,0.,voidFractionM().maxCellsPerParticle());
        dataExchangeM().allocateArray(particleWeights_,0.,voidFractionM().maxCellsPerParticle());
        dataExchangeM().allocateArray(particleVolumes_,0.,voidFractionM().maxCellsPerParticle());
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
      - fvc::div(voidfractionNuEff(voidfraction)*dev(fvc::grad(U)().T()))
    );
}

tmp<volScalarField> cfdemCloud::ddtVoidfraction() const
{
    if (dataExchangeM().couplingStep() <= 2 || !useDDTvoidfraction_)
    {
        Info << "suppressing ddt(voidfraction)" << endl;
        return tmp<volScalarField> (ddtVoidfraction_ * 0.);
    }
    return tmp<volScalarField> (ddtVoidfraction_) ;
}

void cfdemCloud::calcDdtVoidfraction(volScalarField& voidfraction) const
{
    Info << "calculating ddt(voidfraction) based on couplingTime" << endl;
    scalar scale=mesh().time().deltaT().value()/dataExchangeM().couplingTime();
    ddtVoidfraction_ = fvc::ddt(voidfraction) * scale;
}

tmp<volScalarField> cfdemCloud::voidfractionNuEff(volScalarField& voidfraction) const
{
    if (modelType_=="A")
    {
        return tmp<volScalarField>
        (
            #ifdef comp
                new volScalarField("viscousTerm", (turbulence_.mut() + turbulence_.mu()))
            #else
                new volScalarField("viscousTerm", (turbulence_.nut() + turbulence_.nu()))
            #endif
        );
    }
    else
    {
        return tmp<volScalarField>
        (
            #ifdef comp
                new volScalarField("viscousTerm", voidfraction*(turbulence_.mut() + turbulence_.mu()))
            #else
                new volScalarField("viscousTerm", voidfraction*(turbulence_.nut() + turbulence_.nu()))
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
// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "cfdemCloudIO.C"

// ************************************************************************* //
