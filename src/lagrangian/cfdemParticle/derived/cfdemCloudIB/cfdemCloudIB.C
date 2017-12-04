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
#include "cfdemCloudIB.H"
#include "voidFractionModel.H"
#include "forceModel.H"
#include "locateModel.H"
#include "dataExchangeModel.H"
#include "IOModel.H"
#include "mpi.h"
#include "IOmanip.H"
#include "OFversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
cfdemCloudIB::cfdemCloudIB
(
    const fvMesh& mesh
)
:
    cfdemCloud(mesh),
    angularVelocities_(NULL),
    DEMTorques_(NULL),
    pRefCell_(readLabel(mesh.solutionDict().subDict("PISO").lookup("pRefCell"))),
    pRefValue_(readScalar(mesh.solutionDict().subDict("PISO").lookup("pRefValue"))),
    haveEvolvedOnce_(false),
    skipLagrangeToEulerMapping_(false),
    skipAfter_(false),
    timeStepsToSkip_(0),
    calculateTortuosity_(false),
    frontMeshRefine_(false)
{
    if(this->couplingProperties().found("skipLagrangeToEulerMapping"))
    {
        Info << "Will skip lagrange-to-Euler mapping..." << endl;
        skipLagrangeToEulerMapping_=true;
    }
    if(this->couplingProperties().found("timeStepsBeforeSkipping"))
    {
        skipAfter_=true;
        timeStepsToSkip_ =  readScalar
        (
            this->couplingProperties().lookup("timeStepsBeforeSkipping")
        );
        Info << "Will skip LagrangeToEuler mapping after " << timeStepsToSkip_ << " time steps" <<  endl;
    }
    if(this->couplingProperties().found("tortuosity"))
    {
        calculateTortuosity_ = true;
        flowDir_ = this->couplingProperties().subDict("tortuosity").lookup("flowDirection");
        flowDir_ = flowDir_ / mag(flowDir_);
        Info << "Will calculate tortuosity in the mean flow direction ("<<flowDir_[0]<<" "<<flowDir_[1]<<" "<<flowDir_[2]<<")"<< endl;
    }

     //Must check for walls in case of checkPeriodicCells
     //periodic check will mirror particles and probing points to ensure proper behavior near processor bounds
     if(checkPeriodicCells_)
     {
        //Enforce reading of the blocking for periodic checks
        if(readBool(this->couplingProperties().subDict("wall_blockPeriodicityCheck").lookup("x")))
            wall_periodicityCheckRange_[0] = 0;
        if(readBool(this->couplingProperties().subDict("wall_blockPeriodicityCheck").lookup("y")))
            wall_periodicityCheckRange_[1] = 0;
        if(readBool(this->couplingProperties().subDict("wall_blockPeriodicityCheck").lookup("z")))
            wall_periodicityCheckRange_[2] = 0;

        if(this->couplingProperties().found("wall_periodicityCheckTolerance"))
            wall_periodicityCheckTolerance_ = readScalar (this->couplingProperties().lookup("wall_periodicityCheckTolerance"));
     }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cfdemCloudIB::~cfdemCloudIB()
{
    dataExchangeM().destroy(angularVelocities_,3);
    dataExchangeM().destroy(dragPrev_,3);
    dataExchangeM().destroy(DEMTorques_,3);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::cfdemCloudIB::getDEMdata()
{
    cfdemCloud::getDEMdata();
    dataExchangeM().getData("omega","vector-atom",angularVelocities_);
}

bool Foam::cfdemCloudIB::reAllocArrays() const
{
    if(cfdemCloud::reAllocArrays())
    {
        Info <<"Foam::cfdemCloudIB::reAllocArrays()"<<endl;
        dataExchangeM().allocateArray(angularVelocities_,0,3);
        dataExchangeM().allocateArray(dragPrev_,0,3);
        dataExchangeM().allocateArray(DEMTorques_,0,3);
        return true;
    }
    return false;
}

void Foam::cfdemCloudIB::giveDEMdata()
{

    cfdemCloud::giveDEMdata();
    dataExchangeM().giveData("hdtorque","vector-atom", DEMTorques_);
}

inline double ** Foam::cfdemCloudIB::DEMTorques() const
{
    return DEMTorques_;
}


bool Foam::cfdemCloudIB::evolve
(
    volScalarField& alpha,
    volScalarField& interFace
)
{
    numberOfParticlesChanged_ = false;
    arraysReallocated_=false;
    bool doCouple=false;

    if(skipAfter_) {
      if(timeStepsToSkip_<1)
        skipLagrangeToEulerMapping_=true;
    }

    if(!writeTimePassed_ && mesh_.time().outputTime()) writeTimePassed_=true;
    if (dataExchangeM().doCoupleNow())
    {
        Info << "\n timeStepFraction() = " << dataExchangeM().timeStepFraction() << endl;
        dataExchangeM().couple(0);
        doCouple=true;

//        Info << "skipLagrangeToEulerMapping_: " << skipLagrangeToEulerMapping_ 
//             << " haveEvolvedOnce_: " << haveEvolvedOnce_ << endl;
        if(!skipLagrangeToEulerMapping_ || !haveEvolvedOnce_)
        {
            if(verbose_) Info << "- getDEMdata()" << endl;
            getDEMdata();
            Info << "nr particles = " << numberOfParticles() << endl;

            // search cellID of particles
            if(verbose_) Info << "- findCell()" << endl;
            locateM().findCell(NULL,positions_,cellIDs_,numberOfParticles());
            if(verbose_) Info << "findCell done." << endl;

            // set void fraction field
            if(verbose_) Info << "- setvoidFraction()" << endl;
            voidFractionM().setvoidFraction(NULL,voidfractions_,particleWeights_,particleVolumes_,particleV_);
            if(verbose_) Info << "setvoidFraction done." << endl;

            if(verbose_) Info <<"setInterFace"<<endl;
            setInterFace(interFace);
            if(verbose_) Info <<"setInterFace done"<<endl;
        }

        // update voidFractionField
        alpha == voidFractionM().voidFractionNext(); // there might be a better approach, see cfdemCloud.C
        alpha.correctBoundaryConditions();

        // set particles forces
        if(verbose_) Info << "- setForce(forces_)" << endl;
        for(int index = 0;index <  numberOfParticles_; ++index){
            for(int i=0;i<3;i++){
                impForces_[index][i] = 0;
                expForces_[index][i] = 0;
                DEMForces_[index][i] = 0;
            }
        }
        for (int i=0;i<nrForceModels();i++) forceM(i).setForce();
        if(verbose_) Info << "setForce done." << endl;

        // write DEM data
        if(verbose_) Info << " -giveDEMdata()" << endl;
        giveDEMdata();

        dataExchangeM().couple(1);
        
        haveEvolvedOnce_=true;
    }
    Info << "evolve done." << endl;

    //if(verbose_)    #include "debugInfo.H";

    // do particle IO
    IOM().dumpDEMdata();
    if(skipAfter_)
    {
        timeStepsToSkip_--;
        Info << "Will skip LagrangeToEuler mapping after " << timeStepsToSkip_ << " time steps" <<  endl;
    }

    return doCouple;
}

//defines the mesh refinement zone around a particle
//twice the particle size in each direction
void Foam::cfdemCloudIB::setInterFace
(
    volScalarField& interFace
)
{
    interFace == dimensionedScalar("zero", interFace.dimensions(), 0.);
    for(int par=0; par< numberOfParticles(); par++)
    {
        vector ParPos(positions()[par][0],positions()[par][1],positions()[par][2]);
        const boundBox& globalBb = mesh().bounds();
        double skin = 2.0;
        forAll(mesh_.C(),cellI)
        {
            vector posC = mesh_.C()[cellI];
            if(checkPeriodicCells_)
            {
                // Some cells may be located on the other side of a periodic boundary.
                // In this case, the particle center has to be mirrored in order to correctly
                // evaluate the interpolation points.
                vector minPeriodicParticlePos=ParPos;
                voidFractionM().minPeriodicDistance(par,posC, ParPos, globalBb, 
                                                    minPeriodicParticlePos, 
                                                    wall_periodicityCheckRange());

                ParPos = minPeriodicParticlePos;
            }
            double value = voidFractionM().pointInParticle(par, ParPos, posC, skin);
            if(value <= 0.0)
            {
                interFace[cellI] = value + 1.0;
            }
        }
    }
}

void Foam::cfdemCloudIB::calcVelocityCorrection
(
    volScalarField& p,
    volVectorField& U,
    volScalarField& phiIB,
    volScalarField& voidfraction
)
{
    setParticleVelocity(U);

    // make field divergence free - set reference value in case it is needed
    fvScalarMatrix phiIBEqn
    (
        fvm::laplacian(phiIB) == fvc::div(U) + fvc::ddt(voidfraction)
    );
    if(phiIB.needReference()) 
    {
         phiIBEqn.setReference(pRefCell_, pRefValue_);
    }
    
    phiIBEqn.solve();

    U=U-fvc::grad(phiIB);
    U.correctBoundaryConditions();

    // correct the pressure as well
    p=p+phiIB/U.mesh().time().deltaT();  // do we have to  account for rho here?
    p.correctBoundaryConditions();

    if (couplingProperties_.found("checkinterface"))
    {
          Info << "checking no-slip on interface..." << endl;
//          #include "checkInterfaceVelocity.H" //TODO: check carefully!
    }
}

void Foam::cfdemCloudIB::setParticleVelocity
(
    volVectorField& U
)
{
    label cellI=0;
    vector uParticle(0,0,0);
    vector rVec(0,0,0);
    vector velRot(0,0,0);
    vector angVel(0,0,0);
    
    for(int index=0; index < numberOfParticles(); index++)
    {
        for(int subCell=0;subCell<cellsPerParticle()[index][0];subCell++)
        {
            //Info << "subCell=" << subCell << endl;
            cellI = cellIDs()[index][subCell];

            if (cellI >= 0)
            {
                // calc particle velocity
                for(int i=0;i<3;i++) rVec[i]=U.mesh().C()[cellI][i]-position(index)[i];
                for(int i=0;i<3;i++) angVel[i]=angularVelocities()[index][i];
                velRot=angVel^rVec;
                for(int i=0;i<3;i++) uParticle[i] = velocities()[index][i]+velRot[i];

                // impose field velocity
                U[cellI]=(1-voidfractions_[index][subCell])*uParticle+voidfractions_[index][subCell]*U[cellI];
            }
        }
    }
    U.correctBoundaryConditions();
}

vector Foam::cfdemCloudIB::angularVelocity(int index)
{
    vector vel;
    for(int i=0;i<3;i++) vel[i] = angularVelocities_[index][i]; 
    return vel;
}

double Foam::cfdemCloudIB::getTortuosity(vector dir)
{
    volVectorField U = mesh_.lookupObject<volVectorField>("U");
    volScalarField voidfraction = mesh_.lookupObject<volScalarField>("voidfraction");
    double ux = 0.0;
    double umag = 0.0;
    forAll(mesh_.V(),cellI)
    {
        if(voidfraction[cellI] > 0.5)
        {
            double V = mesh_.V()[cellI];
            ux += ((U[cellI] & dir))*V;
            umag += mag(U[cellI])*V;
        }
    }
    //double ux_reduced = 0.0;
    //double umag_reduced = 0.0;
    //MPI_Allreduce(&ux, &ux_reduced, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    //MPI_Allreduce(&umag, &umag_reduced, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    reduce(umag, sumOp<scalar>());
    reduce(ux, sumOp<scalar>());
    double tortuosity = ux == 0.0 ? 1.0 : umag / ux;
    return tortuosity;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void cfdemCloudIB::setRefinementField(volScalarField* refine_)
{
 //Function to allow for setting and activating special refinement operations
 frontMeshRefineField_ = refine_;
 frontMeshRefine_ = true;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
