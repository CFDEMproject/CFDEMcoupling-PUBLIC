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
    remoteScalarInterp_(NULL),
    remoteVectorInterp_(NULL),
    displs_(NULL),
    pRefCell_(readLabel(mesh_.solutionDict().subDict("PISO").lookup("pRefCell"))),
    pRefValue_(readScalar(mesh_.solutionDict().subDict("PISO").lookup("pRefValue"))),
    haveEvolvedOnce_(false),
    skipLagrangeToEulerMapping_(false),
    useHFDIBM_(false),
    checkPeriodicCells_(false)
{

    if(this->couplingProperties().found("skipLagrangeToEulerMapping"))
    {
        Info << "Will skip lagrange-to-Euler mapping..." << endl;
        skipLagrangeToEulerMapping_=true;
    }
    
    if(this->couplingProperties().found("useHFDIBM"))
    {
        Info << "Will use Hybrid Fictitious Domain / Immerse Boundary Method" << endl;
        useHFDIBM_=true;
        HFDIBMinterpDict_=this->couplingProperties().subDict("HFDIBMProps").subDict("interpFunctions");
        if(this->couplingProperties().subDict("HFDIBMProps").found("checkPeriodicCells"))
         checkPeriodicCells_=true;
    }
    
    //get MPI info for communication
    MPI_Comm_rank(MPI_COMM_WORLD,&me);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    
    dataExchangeM().allocateArray(displs_, 0, nprocs);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cfdemCloudIB::~cfdemCloudIB()
{
    dataExchangeM().destroy(angularVelocities_,1);
    dataExchangeM().destroy(remoteScalarInterp_);
    dataExchangeM().destroy(remoteVectorInterp_);
    dataExchangeM().destroy(displs_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::cfdemCloudIB::getDEMdata()
{
    cfdemCloud::getDEMdata();
    Info << "=== cfdemCloudIB::getDEMdata() === particle rotation not considered in CFD" << endl;
    //dataExchangeM().getData("omega","vector-atom",angularVelocities_);
}

bool Foam::cfdemCloudIB::reAllocArrays() const
{
    if(cfdemCloud::reAllocArrays())
    {
        dataExchangeM().allocateArray(angularVelocities_,0,3);
    }
    return true;
}

bool Foam::cfdemCloudIB::evolve
(
    volScalarField& alpha
)
{
    numberOfParticlesChanged_ = false;
    arraysReallocated_=false;
    bool doCouple=false;

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
            
            if(useHFDIBM_)
            {
              // set void fraction field
             if(verbose_) Info << "- setInterpolationPoints()" << endl;
             setInterpolationPoints();
             if(verbose_) Info << "setInterpolationPoints done." << endl;            
            }
        }

        // update voidFractionField
        alpha.internalField() = voidFractionM().voidFractionNext().internalField(); // there might be a better approach, see cfdemCloud.C
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

    return doCouple;
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

    for(int index=0; index< numberOfParticles(); index++)
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

               //That is misleading! The function is named setParticleVelocity and thus should return
               //ONLY the particle velocity field. Then the user can change it as he likes. 
                // impose field velocity
              if(useHFDIBM_)
              {
               U[cellI]=uParticle;
              }
              else
              {
               U[cellI]=(1-voidfractions_[index][subCell])*uParticle+voidfractions_[index][subCell]*U[cellI];
              }
              
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
// * * * * * * * * * * * * * * *  HFDIBM Functions  * * * * * * * * * * * * * //
void Foam::cfdemCloudIB::setInterpolationPoints()
{
  
    if(verbose_) Info << "Reallocating arrays for HFDIBM" << endl; 
    reallocateHFDIBMarrays();
    if(verbose_) Info << "-Reallocation done"<<endl;
    
    if(verbose_) Pout << "Number of particles: " << numberOfParticles()<<endl;
    
    std::vector<double>     localRemoteInterpolationPoints_;  
    
    for(int par=0; par< numberOfParticles(); par++)
    {
       
       scalar radius = radii()[0][par];
        
        for(int subCell=0;subCell<cellsPerParticle()[par][0];subCell++)
        {
           label cellI = cellIDs()[par][subCell];
           vector ParPos(positions()[par][0],positions()[par][1],positions()[par][2]);
          
         // Loop over all the particle cells and collect the ones belonging
         // to the "surface" (i.e 0.001 < phi < 0.999)
           
          if(cellI < 0) continue;
                 
          if( voidFractionM().voidFractionNext()[cellI] < 1e-03       || 
              voidFractionM().voidFractionNext()[cellI] > ( 1 - 1e-03)
            ) continue;          //Take only boundary cells
          
           vector posC = mesh_.C()[cellI];
     
           double res_ = 1.2*std::pow( mesh_.V()[cellI] , 0.3333 ); //distance between interpPoints
         
                      
         if(checkPeriodicCells_) 
         {
          // Some cells may be located on the other side of a periodic boundary.
          // In this case, the particle center has to be mirrored in order to correctly
          // evaluate the interpolation points. 
          // 
          // Notice that a tolerance of 20% is used which means that the domain should be relatively
          // larger ( >140% of the particle size) in order to trigger the algorithm OR the mesh should
          // not be coarser that 5 cells per particle diameter.
         
           for(int dir=0;dir<3;dir++)
           {
            if( std::abs(ParPos[dir] - posC[dir]) > 1.2*radius)
            {
             
              ParPos[dir] -= mesh().bounds().max()[dir] - mesh().bounds().min()[dir];
              if( std::abs(ParPos[dir] - posC[dir]) > 1.2*radius)
              {
                ParPos[dir] += 2*(mesh().bounds().max()[dir] - mesh().bounds().min()[dir]);
              }
            
            }
           
           }
          }
           
           posC = posC - ParPos;
           scalar rC = mag(posC) + 1e-08;
           scalar theta_ = std::acos(posC[2]/rC);
           scalar phi_   = std::atan2(posC[1],posC[0] );
           
           // Surface position and interpolation points are initially calculated relatively 
           // to the particle center
     
           vector posS(radius * std::sin(theta_) * std::cos(phi_),
                       radius * std::sin(theta_) * std::sin(phi_),
                       radius * std::cos(theta_)                    
                      ); 
                
     
           vector posP1((radius + res_) * std::sin(theta_) * std::cos(phi_),
                        (radius + res_) * std::sin(theta_) * std::sin(phi_),
                        (radius + res_) * std::cos(theta_)                    
                       ); 
     
     
          
    
           vector posP2((radius + 2*res_) * std::sin(theta_) * std::cos(phi_),
                        (radius + 2*res_) * std::sin(theta_) * std::sin(phi_),
                        (radius + 2*res_) * std::cos(theta_)                    
                       );
           
           // Get absolute position in mesh
                
           posS  = posS  + ParPos;
           posP1 = posP1 + ParPos;
           posP2 = posP2 + ParPos; 
           
          if(checkPeriodicCells_)
          {
           // If P1 or P2 are crossing the periodic boundary, they should be mirrored
           
           // Check P1
           for(int dir=0;dir<3;dir++)
           {
            if( std::abs( posP1[dir] ) > mesh().bounds().max()[dir])
            {
             
              posP1[dir] -= mesh().bounds().max()[dir] - mesh().bounds().min()[dir];
             }
            if( std::abs(posP1[dir] ) < mesh().bounds().min()[dir])
            {
                posP1[dir] += (mesh().bounds().max()[dir] - mesh().bounds().min()[dir]);
            }
            
           }
            
           
           // Check P2
           for(int dir=0;dir<3;dir++)
           {
            if( std::abs(posP2[dir] ) > mesh().bounds().max()[dir])
            {
             
              posP2[dir] -= mesh().bounds().max()[dir] - mesh().bounds().min()[dir];
             }
            if( std::abs(posP2[dir] ) < mesh().bounds().min()[dir])
            {
                posP2[dir] += (mesh().bounds().max()[dir] - mesh().bounds().min()[dir]);
            }
            
           }
           
          }
           
           
          // Now we look for the cells containing the interpolation points.
          // If they are not available (i.e. cellID=-1) they are added to the
          // list of points to interpolate remotely.
           
           label cellP1 = locateModel_->findSingleCell(posP1,cellI);
           
           //if cannot find, it is probably in a different processor
           if(cellP1==-1)
           {
             localRemoteInterpolationPoints_.push_back(posP1.component(0));
             localRemoteInterpolationPoints_.push_back(posP1.component(1));
             localRemoteInterpolationPoints_.push_back(posP1.component(2));
           }
           else
           {
            //Check if inside particle
           if( voidFractionM().voidFractionNext()[cellP1] < 1- 1e-03 ) 
               cellP1 = -2;   
           
           }
    
           label cellP2 = locateModel_->findSingleCell(posP2,cellI);
   
           //if cannot find, it is probably in a different processor 
           if(cellP2==-1)
           {
             localRemoteInterpolationPoints_.push_back(posP2.component(0));
             localRemoteInterpolationPoints_.push_back(posP2.component(1));
             localRemoteInterpolationPoints_.push_back(posP2.component(2));
           }
           else
           {
            
           //Check if inside particle
           if( voidFractionM().voidFractionNext()[cellP2] < 1- 1e-03 ) 
               cellP2 = -2;   
          
           }
           
           
           
          
           // The cell ID is pushed back even when -1.
           // In this way, the order in wich -1 appears in the 
           // interpCells_ vectors is the same in the localRemoteInterpolationPoints_
           // (i.e. they refer to the same interpolation point).
           
           interpP_[0][par].push_back(posP1);
           interpP_[1][par].push_back(posP2);
           
           interpCells_[0][par].push_back(cellP1);
           interpCells_[1][par].push_back(cellP2);
           
           surfaceCells_[par].push_back(cellI);
           
         
         
        } 
    }
    
    communicateRemoteInterpolationPoints(localRemoteInterpolationPoints_);
  
}
//-------------------------------------------------//
void Foam::cfdemCloudIB::reallocateHFDIBMarrays()
{
  std::vector< label  >  labelTmp_;
  std::vector< vector >  doubleTmp_;
  
  //delete old vectors and reserve new space
  
   interpP_[0].clear();
   interpP_[0].reserve(numberOfParticles());
   
   interpP_[1].clear();
   interpP_[1].reserve(numberOfParticles());
  
  
  interpCells_[0].clear();
  interpCells_[1].clear();
  
  interpCells_[0].reserve(numberOfParticles());
  interpCells_[1].reserve(numberOfParticles());
  
  
  surfaceCells_.clear();
  
  surfaceCells_.reserve(numberOfParticles()); 
  
  
  //Reserve space
  for(int par=0; par< numberOfParticles(); par++)
  {
    //create vectors for particle and reserve space
    
    //Provide an estimation of the space to reserve, i.e NofCells^(2/3) to estimate discrete particle surface
    int spaceToReserve_ =  int( std::pow( ( cellsPerParticle()[par][0] ) , 0.6666 ) );
    
    interpP_[0].push_back(doubleTmp_);
    interpP_[0][par].reserve(spaceToReserve_);
    interpP_[1].push_back(doubleTmp_);
    interpP_[1][par].reserve(spaceToReserve_);
        
        
    interpCells_[0].push_back(labelTmp_);
    interpCells_[1].push_back(labelTmp_);
    
    interpCells_[0][par].reserve(spaceToReserve_);
    interpCells_[1][par].reserve(spaceToReserve_);
    
        
    surfaceCells_.push_back(labelTmp_);    
    surfaceCells_[par].reserve(spaceToReserve_);


  }
 
}

//-------------------------------------------------//
void Foam::cfdemCloudIB::communicateRemoteInterpolationPoints(std::vector<double> localRemoteInterpolationPoints_)
{
  
  if(nprocs == 1) return;
  
  int sizeP = localRemoteInterpolationPoints_.size();
  int nOfPoints = localRemoteInterpolationPoints_.size()/3;
  
  if(verbose_) 
   Pout << "Processor " << me << " found " << nOfPoints << " remote interpolation points" << endl;
   
  int numfrags[nprocs];
  
  // The localRemoteInterpolationPoints_ vectors are gathered by
  // every processor since (also because of periodicity) the "not found"
  // interpolation points could be located in any processor. 
  
  //Get sizeP from every processor
   MPI_Allgather( &sizeP, 1, MPI_INT, numfrags, 1, MPI_INT, MPI_COMM_WORLD);

  //create displacement array
  int totSize = 0;
  for(int p=0;p<nprocs;p++)
  {
   displs_[p] = totSize;
   
   totSize += numfrags[p];
  
  }
  
  //Allocate space for remoteInterpolationPoints
  //and initialise empty vector
  double tmpInterpoPoints[totSize];
   
  MPI_Allgatherv( &localRemoteInterpolationPoints_[0],
                                          numfrags[me],
                                            MPI_DOUBLE,
                                  &tmpInterpoPoints[0],
                                              numfrags,
                                                displs_,
                                            MPI_DOUBLE,
                                        MPI_COMM_WORLD
                );
   
  
  remoteInterpolationPoints.clear();
  remoteInterpolationPoints.assign(tmpInterpoPoints,tmpInterpoPoints+totSize);
  
  // Now that the communication is complete, we perform some more operations 
  // to improve the speed during interpolation.
  
  //Create legend & list of cells for fast access (avoid using findCell() during interpolation)
  int totNofRemotePoints = remoteInterpolationPoints.size()/3;
  
  //initialize vector of velid remote points
  double remoteFDTmp_[totNofRemotePoints];
  for( int id_=0; id_<totSize;id_++)
   remoteFDTmp_[id_]=0;
   
  switchRemoteToFD_.assign(totNofRemotePoints,0);
  
  for(int point_=0;point_<totNofRemotePoints;point_++)
  {
   
   int currentId = point_*3;
   vector P(  remoteInterpolationPoints[currentId    ],
              remoteInterpolationPoints[currentId + 1],
              remoteInterpolationPoints[currentId + 2]
           );
   
   label cell0 = 0;
   label cellP = locateModel_->findSingleCell(P,cell0);
   
   if(cellP != -1)
   {
     
     //Check if inside particle
    if( voidFractionM().voidFractionNext()[cellP] < 1- 1e-03 ) 
    {
      remoteFDTmp_[point_] = 1; 
      continue;
    }
    remoteInterpLegend.push_back(currentId);
    remoteInterpCells.push_back(cellP);
    
      
   
   }
  
  
  }
  
  dataExchangeM().allocateArray(remoteScalarInterp_, 0, totNofRemotePoints);
  dataExchangeM().allocateArray(remoteVectorInterp_, 0, totSize);

  //now everyone knows if is valid or should switch to FD
  MPI_Allreduce(remoteFDTmp_, &switchRemoteToFD_[0], totNofRemotePoints ,MPI_INT, MPI_SUM, MPI_COMM_WORLD);
 
}
//--------------------------------------------------------------------------------------------//
void Foam::cfdemCloudIB::vectorInterpolateRemote( volVectorField& V, volVectorField& Vs)
{
 
 autoPtr<interpolation<vector> > interpV_ = interpolation<vector>::New(HFDIBMinterpDict_, V);

 
 int totSize=remoteInterpolationPoints.size();
 
 double remoteVectorInterpTmp_[totSize];
 //Set old array to zero
 for( int id_=0; id_<totSize;id_++)
 {
   remoteVectorInterp_[id_]=0.0;
   remoteVectorInterpTmp_[id_]=0.0;
 } 

 // Following the remoteInterpLegend, each processor interpolates the
 // remote points it owns.
 
 for(unsigned int point_=0; point_< remoteInterpLegend.size();point_++)
 {
   //Get cellId and position of the interpolation point
   label cellI = remoteInterpCells[point_];
   int currentId = remoteInterpLegend[point_];
   vector P(  remoteInterpolationPoints[currentId    ],
              remoteInterpolationPoints[currentId + 1],
              remoteInterpolationPoints[currentId + 2]
           );
           
   //interpolate
   vector VP =  interpV_->interpolate( P, cellI ) ;
   
   //assign to array
   remoteVectorInterpTmp_[currentId    ] = VP.component(0);
   remoteVectorInterpTmp_[currentId + 1] = VP.component(1);
   remoteVectorInterpTmp_[currentId + 2] = VP.component(2);
 
 }
 
 
 // The use of MPI_Allreduce allows to keep the original order.
 // Thus, the reduced vector is still synchronized with the local
 // interpCells_ vector.
 
 MPI_Allreduce(remoteVectorInterpTmp_, remoteVectorInterp_, totSize ,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

}

//--------------------------------------------------------------------------------------------//
void Foam::cfdemCloudIB::scalarInterpolateRemote( volScalarField& S,volScalarField& Ss)
{
 autoPtr<interpolation<scalar> > interpS_ = interpolation<scalar>::New(HFDIBMinterpDict_, S);
 int totSize=remoteInterpolationPoints.size()/3;
 //Set old array to zero
  double remoteScalarInterpTmp_[totSize];
 for( int id_=0; id_<totSize;id_++)
 {
   remoteScalarInterp_[id_]=0.0;
   remoteScalarInterpTmp_[id_]=0.0;
 }
 
 // Following the remoteInterpLegend, each processor interpolates the
 // remote points it owns.
 
 for(unsigned int point_=0; point_< remoteInterpLegend.size();point_++)
 {
  
   //Get cellId and position of the interpolation point
   label cellI = remoteInterpCells[point_];
   int currentId = remoteInterpLegend[point_];
   vector P(  remoteInterpolationPoints[currentId    ],
              remoteInterpolationPoints[currentId + 1],
              remoteInterpolationPoints[currentId + 2]
           );
           
   //interpolate
   scalar SP =  interpS_->interpolate( P, cellI ) ;
   
   //assign to array
   remoteScalarInterpTmp_[currentId/3] = SP;
   
   //Ss[cellI] = 9999999999999*2;
   
 }
 
 
 // The use of MPI_Allreduce allows to keep the original order.
 // Thus, the reduced vector is still synchronized with the local
 // interpCells_ vector.
 
 MPI_Allreduce(remoteScalarInterpTmp_, remoteScalarInterp_, totSize ,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

}
//--------------------------------------------------------------------------------------------//

void Foam::cfdemCloudIB::interpolateVectorField( volVectorField                        &V,   //main field
                                                 volVectorField                        &Vs   //imposed field
                                               )
{
  if(!useHFDIBM_) return;
  
  //create Interpolator
  autoPtr<interpolation<vector> > interpV_ = interpolation<vector>::New(HFDIBMinterpDict_, V);
  
  
  //First of all, interpolate remote points
  vectorInterpolateRemote(V, Vs);
  
  int remoteCount = 0;
  
   for(int par=0; par< numberOfParticles(); par++)
   {
    
  
   
    scalar radius = radii()[0][par];
    
    int surfCls_ = surfaceCells_[par].size();
    
    for(int surf_=0; surf_<surfCls_; surf_++)
    {
     
      label cellI = surfaceCells_[par][surf_];
      
      vector posC = mesh_.C()[cellI];
      
      vector ParPos(positions()[par][0],positions()[par][1],positions()[par][2]);

      bool useFD=false;
      
      if(checkPeriodicCells_)
      {      
       // Again, the particle position could need a mirroring
       // in order to correctly calculate the distance between the 
       // surface cell and the particle center.
       
           for(int dir=0;dir<3;dir++)
           {
            if( std::abs(ParPos[dir] - posC[dir]) > 1.2*radius)
            {
             
              ParPos[dir] -= mesh().bounds().max()[dir] - mesh().bounds().min()[dir];
              if( std::abs(ParPos[dir] - posC[dir]) > 1.2*radius)
              {
                ParPos[dir] += 2*(mesh().bounds().max()[dir] - mesh().bounds().min()[dir]);
              }
            
            }
           
           }
      
      }
      
      double rC = mag(posC - ParPos);
      
      vector VP1(0.0,0.0,0.0);
      vector VP2(0.0,0.0,0.0);
      
      // If the cellID is -1, look in remoteVectorInterp_ to find the
      // remotely intrpolated values.
      
      if(interpCells_[0][par][surf_] > -1)
      {
        VP1 =  interpV_->interpolate(      interpP_[0][par][surf_],
                                           interpCells_[0][par][surf_]
                         
                                     )  - Vs[cellI];
      }
      else if(interpCells_[0][par][surf_] != -2)
      {
       
        if(switchRemoteToFD_[displs_[me]/3 + remoteCount]==1)
        {
         useFD=true;
          
        
        }
        else
        {
         VP1.component(0) = remoteVectorInterp_[displs_[me] + remoteCount*3   ];
         VP1.component(1) = remoteVectorInterp_[displs_[me] + remoteCount*3 +1];
         VP1.component(2) = remoteVectorInterp_[displs_[me] + remoteCount*3 +2];
        
         VP1 = VP1 - Vs[cellI];
        }
         remoteCount++;
      }
      else
      {
        useFD=true;
      } 
      
      
      if(interpCells_[1][par][surf_] > -1)
      {
        VP2 =  interpV_->interpolate(      interpP_[1][par][surf_],
                                           interpCells_[1][par][surf_]
                         
                                     )  - Vs[cellI];
      }
      else if(interpCells_[1][par][surf_] != -2)
      {
       
        if(switchRemoteToFD_[displs_[me]/3 + remoteCount]==1)
        {
         useFD=true;
        }
        else
        {
         VP2.component(0) = remoteVectorInterp_[displs_[me] + remoteCount*3   ];
         VP2.component(1) = remoteVectorInterp_[displs_[me] + remoteCount*3 +1];
         VP2.component(2) = remoteVectorInterp_[displs_[me] + remoteCount*3 +2];
         
         VP2 = VP2 - Vs[cellI];
       
        }
         remoteCount++;
      }                         
                                                   
      
     if(useFD)
     {
      //Vs[cellI] =  voidFractionM().voidFractionNext()[cellI]*V[cellI] + (1- voidFractionM().voidFractionNext()[cellI])*Vs[cellI];
      continue;
     }
      
      double res_ = 1.2*std::pow( mesh_.V()[cellI] , 0.3333 ); //distance between interpPoints
      
      vector quadCoeff = 1/(res_*res_) * ( VP2/2 - VP1 );
      vector linCoeff  = 1/(2*res_) * ( 4*VP1 - VP2 );
   
    //Correct imposed field
      Vs[cellI] = quadCoeff*(rC-radius)*(rC-radius) + linCoeff * (rC-radius) + Vs[cellI]  ;
    
    }
   }
}
//------------------------------------------------------------------------------------------//
void Foam::cfdemCloudIB::interpolateScalarField( volScalarField                       &S,  //main field
                                                 volScalarField                       &Ss  //imposed field
                                               )
{
  if(!useHFDIBM_) return;
  
  autoPtr<interpolation<scalar> > interpS_ = interpolation<scalar>::New(HFDIBMinterpDict_, S);
  
  //First of all, interpolate remote points
  scalarInterpolateRemote(S, Ss);
  
  int remoteCount = 0;
   for(int par=0; par< numberOfParticles(); par++)
   {
    
   
    scalar radius = radii()[0][par];
    
    int surfCls_ = surfaceCells_[par].size();
    
    for(int surf_=0; surf_<surfCls_; surf_++)
    {
      
      bool useFD=false;
      label cellI = surfaceCells_[par][surf_];
      vector posC = mesh_.C()[cellI];
      
        vector ParPos(positions()[par][0],positions()[par][1],positions()[par][2]);
      
      
      if(checkPeriodicCells_)
      {     
       // Again, the particle position could need a mirroring
       // in order to correctly calculate the distance between the 
       // surface cell and the particle center.
       
           for(int dir=0;dir<3;dir++)
           {
            if( std::abs(ParPos[dir] - posC[dir]) > 1.2*radius)
            {
             
              ParPos[dir] -= mesh().bounds().max()[dir] - mesh().bounds().min()[dir];
              if( std::abs(ParPos[dir] - posC[dir]) > 1.2*radius)
              {
                ParPos[dir] += 2*(mesh().bounds().max()[dir] - mesh().bounds().min()[dir]);
              }
            
            }
           
           }
      
       }
      
      double rC = mag(posC - ParPos);
      
      scalar SP1=0.0;
      scalar SP2=0.0;
      
      // If the cellID is -1, look in remoteVectorInterp_ to find the
      // remotely intrpolated values.
      
      if(interpCells_[0][par][surf_] > -1)
      {
        SP1 =  interpS_->interpolate(      interpP_[0][par][surf_],
                                            interpCells_[0][par][surf_]
                                           
                                            ) - Ss[cellI];
      }
      else if(interpCells_[0][par][surf_] != -2)
      {
       
        if(switchRemoteToFD_[displs_[me]/3 + remoteCount]==1)
        {
         useFD=true;
          
        
        }
        else
        {
         SP1 = remoteScalarInterp_[displs_[me]/3 + remoteCount] -Ss[cellI];
        }
        
        remoteCount++;
      } 
      else
      {
       useFD=true;
      }
      
      if(interpCells_[1][par][surf_] > -1)
      {
        SP2 =  interpS_->interpolate(      interpP_[1][par][surf_],
                                            interpCells_[1][par][surf_]
                                           
                                            ) - Ss[cellI];
      }
       else if(interpCells_[1][par][surf_] != -2)
      {
       
        if(switchRemoteToFD_[displs_[me]/3 + remoteCount]==1)
        {
         useFD=true;
        }
        else
        {
         SP2 = remoteScalarInterp_[displs_[me]/3 + remoteCount] -Ss[cellI];
        }
         remoteCount++;
      } 
      else
      {
       useFD=true;
      }
      
      
     if(useFD)
     {
     // Ss[cellI] =  voidFractionM().voidFractionNext()[cellI]*S[cellI] + (1- voidFractionM().voidFractionNext()[cellI])*Ss[cellI];
      continue;
     }
      
      double res_ = 1.2*std::pow( mesh_.V()[cellI] , 0.3333 ); //distance between interpPoints
      
      scalar quadCoeff = 1/(res_*res_) * ( SP2/2 - SP1 );
      scalar linCoeff  = 1/(2*res_) * ( 4*SP1 - SP2 );
   
      //Correct imposed field
    //  Ss[interpCells_[0][par][surf_]] = -999999999;
    //  Ss[interpCells_[1][par][surf_]] = 999999999;
      Ss[cellI] = quadCoeff*(rC-radius)*(rC-radius) + linCoeff * (rC-radius) + Ss[cellI]  ;
       }
   }
}
//---------------------------------------------------------------------------//
void Foam::cfdemCloudIB::checkInterfaceFlowRate()
{
 
 //This function should just be used to check the solver and for debugging
 
 //Each particle is discretized in several lagrangian points where the velocity field is interpolated
 
 if(!this->couplingProperties().subDict("HFDIBMProps").found("checkInterfaceFlowRate"))
 return;
 
 int thetaDiscr_ = 180;
 int phiDiscr_   = 180;
 
 // Loop over all the particles
 //
 // Interpolate at the discrete surface 
 // point 
 
  double totalFlux_=0.0;
  double pi_=3.14159265359;
  double totSurf=0.0;
  
  volVectorField U = mesh_.lookupObject<volVectorField>("U");
  
  //create Interpolator
  autoPtr<interpolation<vector> > interpV_ = interpolation<vector>::New(HFDIBMinterpDict_, U);
  
  for(int par=0; par< numberOfParticles(); par++)
  {
   
    scalar radius = radii()[0][par];
    
    //Advance of two to correctly consider the angle at the poles
    for(int thetaId_=0;thetaId_<thetaDiscr_/2;thetaId_++)
    {
     
      double theta = (1+thetaId_*2)*pi_/thetaDiscr_;
      vector ParPos(positions()[par][0],positions()[par][1],positions()[par][2]);
     
     for(int phiId_=0;phiId_<phiDiscr_;phiId_++)
     {
     
      double phi=phiId_*phiId_*2*pi_/phiDiscr_;
      
      double surf_= radius*radius*(2*pi_/phiDiscr_)*( cos(theta - pi_/(thetaDiscr_) ) - cos(theta + pi_/(thetaDiscr_) ) );
      totSurf+=surf_;
      
        vector posS(radius * std::sin(theta) * std::cos(phi),
                    radius * std::sin(theta) * std::sin(phi),
                    radius * std::cos(theta)                    
                   );
       posS=posS+ParPos;
       label cell0 = 0; 
       label cellS = locateModel_->findSingleCell(posS,cell0);
       
       if(cellS==-1) continue;
       
       
       vector US =  interpV_->interpolate( posS, cellS);
                                           
       posS = (posS-ParPos)/radius;
       
       //flow normal to the surface
       
       totalFlux_ += std::abs( US&posS*surf_ );  
     }
    }
  }
  
  Pout << "Total flow rate normal to the immersed surfaces : " << totalFlux_ << endl;
  
  Pout << "Total immersed surfaces: " << totSurf << endl; //Should be pi for one particle with dp=1


}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
