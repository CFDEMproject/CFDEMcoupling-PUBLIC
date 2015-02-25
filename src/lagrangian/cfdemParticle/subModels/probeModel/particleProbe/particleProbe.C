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

#include "error.H"

#include "particleProbe.H"
#include "addToRunTimeSelectionTable.H"
#include "mpi.h"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(particleProbe, 0);

addToRunTimeSelectionTable
(
    probeModel,
    particleProbe,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
particleProbe::particleProbe
(
    const dictionary& dict,
    cfdemCloud& sm,
    word   typeName,
    char*  logFileName
)
:
    probeModel(dict,sm,typeName,logFileName),
    propsDict_(dict.subDict(typeName + "Props")),
    name_(typeName),
    particleCloud_(sm),
    verbose_(false),
    verboseToFile_(false),
    writePrecision_(3),
    dirName_("particleProbes"),
    rank_(-1),
    sPtr(NULL),
    printEvery_(1),
    sampleAll_(false),
    probeDebug_(false),
    includePosition_(false),
    particleIDsToSample_(propsDict_.lookup("particleIDsToSample")),
    itemsToSample_(NULL),
    sPtrList_(NULL),
    itemCounter_(0),
    currItemId_(0),
    printCounter_(0),
    printNow_(false)
{
    if (propsDict_.found("verbose")) verbose_=true;
    if (propsDict_.found("verboseToFile")) verboseToFile_=true;

    if (propsDict_.found("printEvery")) printEvery_= readScalar(propsDict_.lookup("printEvery"));
    if (propsDict_.found("sampleAll")) sampleAll_=true;
    if (propsDict_.found("probeDebug")) probeDebug_=true;
    if (propsDict_.found("includePosition")) includePosition_=true;

    if (propsDict_.found("writePrecision")) writePrecision_= readScalar(propsDict_.lookup("writePrecision"));

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

particleProbe::~particleProbe()
{
 clearProbes();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void particleProbe::setOutputFile() const
{
    //set the current item ID
    if(currItemId_== itemCounter_)
        currItemId_=1;
    else
        currItemId_+=1;
    sPtr = sPtrList_[currItemId_-1]; //set the pointer to the output file from list
    probeIndex_=currItemId_-1;
}


void particleProbe::initialize(word typeName, word  logFileName) const
{
  //update the list of items to be sampled
  itemCounter_ += 1; 
  itemsToSample_.append(logFileName);

  // init environment
  //propsDict_ = particleCloud_.couplingProperties().subDict(typeName + "Props");
  name_ = typeName;
  const char* fileNameOut_ = wordToChar(logFileName);
  
  if(verboseToFile_)
  {

    Info << "Will sample these particle IDs: " << particleIDsToSample_ << " every " <<  printEvery_ <<  endl;

    //initialize the output files
    int myrank_(-1);
    int numprocs_(-1);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank_);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs_);
    rank_=myrank_;

    //open a separate file for each processor
    char* filecurrent_;
    filecurrent_= new char[strlen(fileNameOut_) + 4]; //reserve 4 chars for processor name
    if (myrank_ < numprocs_)  
             sprintf(filecurrent_,"%s%s%d", fileNameOut_, ".", myrank_);
    else  //open one file for proc 0
            sprintf(filecurrent_,"%s", fileNameOut_); 
     word file_(filecurrent_);
     Info << "particleProbe for model " <<  name_ << " will write to file " << file_ << endl;

      //generate the file streams
      fileName probeSubDir = dirName_;
      if (particleCloud_.mesh().name() != polyMesh::defaultRegion)
      {
                probeSubDir = probeSubDir/particleCloud_.mesh().name();
       }
      probeSubDir = probeSubDir/particleCloud_.mesh().time().timeName();

       fileName probeDir_;
       if (Pstream::parRun())
       {
                // Put in undecomposed case
                // (Note: gives problems for distributed data running)
               probeDir_ = particleCloud_.mesh().time().path()/".."/probeSubDir;
        }
        else
        {
              probeDir_ = particleCloud_.mesh().time().path()/probeSubDir;
        }

    //manage files and OFstreams
    mkDir(probeDir_);
    sPtr = new OFstream(probeDir_/file_);
    sPtrList_.append(sPtr);

    //Clear the containers for the fields to be probed
    scalarFields_.clear();
    vectorFields_.clear();
  }
  return;

}
void particleProbe::writeHeader() const
{

     if(verboseToFile_ )
     {
         *sPtr<<"#processor: " << rank_ << endl;
         *sPtr<<"#index   time   "  << "   ";  


        *sPtr<<"||  vectorData:  "  << "   "; 
    
         forAll(vectorFields_, iter)
         {
             if(!probeDebug_ && iter>0) break;
             *sPtr << vectorFields_(iter) << "   "; 
          }

         if(probeDebug_)
         {
            *sPtr<<"||   scalarData:  "  << "   ";  
            forAll(scalarFields_, iter)
            {
                 *sPtr << scalarFields_(iter)  << "   ";
            }
         }

         if(includePosition_) *sPtr<<" ||  position" << endl; 
         else *sPtr << endl; 
     }

}

void particleProbe::clearProbes() const
{
  for (unsigned int i=0; i<vProbes_.size(); i++)
    vProbes_[i].clear();
  
  for (unsigned int j=0; j<sProbes_.size(); j++)
    sProbes_[j].clear();
      
  sProbes_.clear();
  vProbes_.clear();
}

void particleProbe::updateProbes(int index, Field<scalar> sValues, Field<vector> vValues) const
{
  int vSize_=vProbes_.size();
  int sSize_=sProbes_.size();
 
 //check if the particle already has an allocated vector. If not, create it. It should be only called at the beginning. 
 while(index >= vSize_)
 {
   std::vector<double*> particleVector_;
   vProbes_.push_back(particleVector_);
   vSize_=vProbes_.size();
 }
  
  while(index >= sSize_)
  {
   std::vector<double> particleScalar_;
   sProbes_.push_back(particleScalar_);
   sSize_=sProbes_.size();
  }
  
 //register vector probes on the corresponding vector
  forAll(vValues, iter)
  {
   int ProbeSize_=vProbes_[index].size();
   
   if(probeIndex_<ProbeSize_) //The corresponding probe for this particle already exists, values are overwritten.
   {
    vProbes_[index][probeIndex_][0]=vValues[iter][0];
    vProbes_[index][probeIndex_][1]=vValues[iter][1];
    vProbes_[index][probeIndex_][2]=vValues[iter][2];
    
   }
   else //The corresponding probe for this particle has to be created
   {
    double * probe_= new double[3];
   
    probe_[0]=vValues[iter][0];
    probe_[1]=vValues[iter][1];
    probe_[2]=vValues[iter][2];
              
    vProbes_[index].push_back(probe_); 
   }
  }
  
  //register scalar probes on the corresponding vector
  forAll(sValues, iter)
  {
   int ProbeSize_=sProbes_[index].size();
   
   if(probeIndex_<ProbeSize_) //The corresponding probe for this particle already exists, values are overwritten.
   {
    sProbes_[index][probeIndex_]=sValues[iter];
   }
   else //The corresponding probe for this particle has to be created
   {           
    sProbes_[index].push_back(sValues[iter]); 
   } 
  
  }
  
  
}

void particleProbe::writeProbe(int index, Field<scalar> sValues, Field<vector> vValues) const
{
    updateProbes(index,sValues,vValues); //update probe vectors
    
   
    if(printNow_ && checkIDForPrint(index) &&  verboseToFile_) 
    {

        //index and time
       *sPtr <<    setprecision(IOstream::defaultPrecision()+7) ;
       *sPtr << index  << tab 
                << particleCloud_.mesh().time().value()  << "   " ;
        *sPtr << "||   ";
        
        //vectorFields
        *sPtr <<    setprecision(writePrecision_) ;
         forAll(vValues, iter)
         {
             // if(!probeDebug_ && iter>0) break;
             *sPtr << vValues[iter][0] << "   ";
             *sPtr << vValues[iter][1] << "   ";
             *sPtr << vValues[iter][2] << "   "; 
   
          }

        //scalarFields
        if(probeDebug_)
        {
          *sPtr << "||   ";
           forAll(sValues, iter)
           {
               *sPtr << sValues[iter] << "   "; 
            }
        }

        if(includePosition_)
        {
            *sPtr << "||   ";
             *sPtr <<     particleCloud_.position(index)[0] << "   " 
                      <<     particleCloud_.position(index)[1] << "   "
                      <<     particleCloud_.position(index)[2]
                      << endl;
        }
        else *sPtr << endl;

    }
    
    return;
    
}

bool particleProbe::checkIDForPrint(int index) const
{
  
      bool sampleThisId_ = false;
      if(sampleAll_) sampleThisId_ = true;
      else
      {
            forAll(particleIDsToSample_, iSample)
            {
                  if(index==particleIDsToSample_[iSample]) sampleThisId_ = true;
            }
      }
      return sampleThisId_;
}

void particleProbe::setCounter() const
{

    //reset or increment counter for printing to file 
    //Do only if called by first item in the list of items!
    if(currItemId_==1)
    {
        printCounter_++;
        if(  printCounter_ >= printEvery_ )
        {
            printCounter_=0;
            printNow_ = true;
        }
        else printNow_ = false;

    }
    return;    

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
