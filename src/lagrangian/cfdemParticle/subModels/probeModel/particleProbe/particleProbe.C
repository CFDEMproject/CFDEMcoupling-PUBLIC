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
    const word& typeName,
    const char* logFileName
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
    itemCounter_(0),
    currItemId_(0),
    printCounter_(0),
    printNow_(false),
    printOnlyAt_(-1)
{
    if (propsDict_.found("verbose")) verbose_=true;
    if (propsDict_.found("verboseToFile")) verboseToFile_=true;

    if (propsDict_.found("printEvery")) printEvery_= readScalar(propsDict_.lookup("printEvery"));
    if (propsDict_.found("printOnlyAtStep")) printOnlyAt_= readScalar(propsDict_.lookup("printOnlyAtStep"));
    if (propsDict_.found("sampleAll")) sampleAll_=true;
    if (propsDict_.found("probeDebug")) probeDebug_=true;
    if (propsDict_.found("includePosition")) includePosition_=true;

    if (propsDict_.found("writePrecision")) writePrecision_= readScalar(propsDict_.lookup("writePrecision"));

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

particleProbe::~particleProbe()
{
 clearProbes();
 forAll(sPtrList_, i)
     delete sPtrList_[i];
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void particleProbe::setOutputFile(const word& logFileName) const
{
    bool foundFile = false;
    if(itemCounter_>0 && verboseToFile_ )
    {
        forAll(itemsToSample_, iterator)
        {
            if(itemsToSample_[iterator] == logFileName)
            {
                probeIndex_ = iterator ;
                foundFile = true;
            }
        }

        if(!foundFile)
             FatalError <<  "particleProbe::setOutputFile for logFileName " <<  logFileName << " : "<< "File not found" << abort(FatalError);
        currItemId_ = probeIndex_ + 1;
        setCounter(); //set counter, will auto-check if first item
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
void particleProbe::initialize(const word& modelName,const word& logFileName) const
{
  //update the list of items to be sampled
  itemCounter_ += 1; 
  // Note: for other than ext one could use vValues.append(x)
  // instead of setSize
  itemsToSample_.setSize(itemsToSample_.size()+1, logFileName);

  // init environment
  name_ = modelName;
  const char* fileNameOut_ = logFileName.c_str();

  if(verboseToFile_)
  {

   if(printOnlyAt_==-1)
   {
    Info << "Will sample these particle IDs: " << particleIDsToSample_ << " every " <<  printEvery_ <<  endl;
   }
   else if(printOnlyAt_>-1)
   {
     Info << "Will sample these particle IDs: " << particleIDsToSample_ << " at step " <<  printOnlyAt_ <<  endl;
   }
   else
   {
    FatalError <<  "particleProbe for model " <<  name_ << " : "<< "printOnlyAtStep cannot be negative or zero" << abort(FatalError);
   }

    //initialize the output files
    MPI_Comm_rank(MPI_COMM_WORLD,&rank_);

    //open a separate file for each processor
    char* filecurrent_ = new char[strlen(fileNameOut_) + 5]; //reserve 4 chars for processor name
    sprintf(filecurrent_,"%s%s%d", fileNameOut_, ".", rank_);
    Info << "particleProbe for model " <<  name_ << " will write to file " << filecurrent_ 
         << ". This is item " << itemCounter_ << " to probe." << endl;

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
    sPtr = new OFstream(probeDir_/filecurrent_);
    // Note: for other than ext one could use xx.append(x)
    // instead of setSize
    sPtrList_.setSize(sPtrList_.size()+1, sPtr);

    delete [] filecurrent_;
    //Clear the containers for the fields to be probed
    scalarFields_.clear();
    vectorFields_.clear();
  }

  Info << "particleProbe::initialize done!" << endl;
  return;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
void particleProbe::clearProbes() const
{

 for(std::vector< std::vector<double*> >::iterator 
       it  = vProbes_.begin(); 
       it != vProbes_.end(); 
     ++it) 
 {
    for(std::vector<double*>::iterator 
           jt  = (*it).begin(); 
           jt != (*it).end(); 
         ++jt) 
    {
        delete[] (*jt);
    }
    (*it).clear();
    
 }

  for (unsigned int j=0; j<sProbes_.size(); j++)
    sProbes_[j].clear();
      
  sProbes_.clear();
  vProbes_.clear();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
void particleProbe::updateProbes(int index, Field<scalar> sValues, Field<vector> vValues) const
{
  //check if the particle already has an allocated vector. If not, create it. It should be only called at the beginning. 
  while(index >= static_cast<int>(vProbes_.size()))
  {
    std::vector<double*> particleVector_;
    vProbes_.push_back(particleVector_);
  }
  
  while(index >= static_cast<int>(sProbes_.size()))
  {
    std::vector<double> particleScalar_;
    sProbes_.push_back(particleScalar_);
  }
  
  //register vector probes on the corresponding vector
  forAll(vValues, iter)
  {
   if(probeIndex_<static_cast<int>(vProbes_[index].size())) //The corresponding probe for this particle already exists, values are overwritten.
   {
    vProbes_[index][probeIndex_][0]=vValues[iter][0];
    vProbes_[index][probeIndex_][1]=vValues[iter][1];
    vProbes_[index][probeIndex_][2]=vValues[iter][2];
   }
   else //The corresponding probe for this particle has to be created
   {
    vProbes_[index].push_back(new double[3]); //this pointer MUST be freed in the destructor
    vProbes_[index].back()[0]=vValues[iter][0];
    vProbes_[index].back()[1]=vValues[iter][1];
    vProbes_[index].back()[2]=vValues[iter][2];
   }
  }
  
  //register scalar probes on the corresponding vector
  forAll(sValues, iter)
  {
   if(probeIndex_<static_cast<int>(sProbes_[index].size())) //The corresponding probe for this particle already exists, values are overwritten.
    sProbes_[index][probeIndex_]=sValues[iter];
   else //The corresponding probe for this particle has to be created
    sProbes_[index].push_back(sValues[iter]); 
  }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
void particleProbe::writeProbe(int index, Field<scalar> sValues, Field<vector> vValues) const
{
    updateProbes(index,sValues,vValues); //update probe vectors
   
    if(printNow_ && checkIDForPrint(index) &&  verboseToFile_) 
    {
        sPtr = sPtrList_[probeIndex_]; //set the pointer to the output file from list

        //index and time
        *sPtr <<    setprecision(IOstream::defaultPrecision()+7) ;
        *sPtr << index  << "   " 
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
void particleProbe::setCounter() const
{

    //reset or increment counter for printing to file 
    //Do only if called by first item in the list of items!
    if(currItemId_==1)
    {
       printCounter_++;
       if(printOnlyAt_==-1)
       {
        if(  printCounter_ >= printEvery_ )
        {
            printCounter_=0;
            printNow_ = true;
        }
        else printNow_ = false;
       }
       else if(printCounter_==printOnlyAt_)
       {
        printNow_ = true;
       }
       else 
       {
        printNow_ = false;
       }
    
    }
    return;    

}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
