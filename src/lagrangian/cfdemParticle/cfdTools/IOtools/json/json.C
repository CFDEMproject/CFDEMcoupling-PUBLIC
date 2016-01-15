/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
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
    This code implements a simple JSON file writer
    
Contributing Authors

    copyright:  Federico Municchi, TU Graz, 2015
                Stefan Radl, TU Graz, 2015


\*---------------------------------------------------------------------------*/

#include "json.H"
#include "error.H"
#include "IOmanip.H"
#include <typeinfo>
#include <stdlib.h>
#include <sys/stat.h>

using namespace Foam;

jsonObject::jsonObject()
:
lastObject(false)
{
};

jsonObject::~jsonObject()
{
};
/*-----------------Auxiliary Functions------------------------------------------*/
void indentation(int ind, OFstream* file)
{
 
 for(int i=0;i<ind;i++)
  *file << " ";

}

/*-----------------Add Functions-----------------------------------------------*/
void jsonObject::addjObject(std::string elemName, jsonObject* data) const
{
 //Register data name
 jObjectNames_.push_back(elemName);
 
 //Register data
 jObjects_.push_back(data);
 data->setIndent(indent+2);
 data->setPtr(sPtr);
  
};

/*-------------------------------------------------------------------------------*/
void jsonObject::addjVector(std::string elemName, std::vector<double>* data) const
{
 //Register data name
 jVectorNames_.push_back(elemName);
 
 //Register data
 jVectors_.push_back(data);
  
};

/*-------------------------------------------------------------------------------*/
void jsonObject::addjScalar(std::string elemName, double* data) const
{
 //Register data name
 jScalarNames_.push_back(elemName);
 
 //Register data
 jScalars_.push_back(data);
  
};

/*-------------------------------------------------------------------------------*/
void jsonObject::addjString(std::string elemName, std::string* data) const
{
 //Register data name
 jStringNames_.push_back(elemName);
 
 //Register data
 jStrings_.push_back(data);
  
};

/*-------------------------------------------------------------------------------*/
void jsonObject::addjBool(std::string elemName, bool* data) const
{
 //Register data name
 jBoolNames_.push_back(elemName);
 
 //Register data
 jBools_.push_back(data);
  
};

/*-----------------Write Functions-----------------------------------------------*/
void jsonObject::write() const
{
 
  //start writing
  indentation(indent-2,sPtr);
  *sPtr<<"{" << endl;
  
  bool comma=false;
  
  //Write Objects
  for(unsigned int i=0;i<jObjects_.size();i++)
  {
   if(i!=0)
    *sPtr<<","<<endl<<endl;
   else
    *sPtr<<endl<<endl; 
    
   indentation(indent,sPtr);
   *sPtr<< jObjectNames_[i] << ":" <<endl;
   
   jObjects_[i]->write();
 
  }
  
  
  if(jObjects_.size()!=0)
   comma=true;
  
  //Write vectors
  for(unsigned int i=0;i<jVectors_.size();i++)
  {
    if(i!=0)
     *sPtr<<","<<endl<<endl;
    else
    {
     if(comma) *sPtr<<",";
     *sPtr<<endl<<endl;  
    }
   indentation(indent,sPtr);
   *sPtr<< jVectorNames_[i] << ": [";
   
   for(unsigned int s=0;s<jVectors_[i]->size();s++)
   {
    if(s!=0)
      *sPtr<<","<<endl;
    else
     *sPtr<<endl;
      
    indentation(indent,sPtr);
     *sPtr<< (*jVectors_[i])[s];
   }
   
    *sPtr<<endl;
    indentation(indent,sPtr);
   *sPtr<<  "]" ;
   
  }
  
   if(jVectors_.size()!=0)
   comma=true;
  
  //Write strings
  for(unsigned int i=0;i<jStrings_.size();i++)
  {
    if(i!=0)
     *sPtr<<","<<endl<<endl;
    else
    {
     if(comma) *sPtr<<",";
     *sPtr<<endl<<endl;  
    }
   indentation(indent,sPtr);
   *sPtr<< jStringNames_[i] << ": " << *jStrings_[i];
  
   
  }
  
   if(jStrings_.size()!=0)
   comma=true;
  
  //Write strings
  for(unsigned int i=0;i<jScalars_.size();i++)
  {
    if(i!=0)
     *sPtr<<","<<endl<<endl;
    else
    {
     if(comma) *sPtr<<",";
     *sPtr<<endl<<endl;  
    }
   indentation(indent,sPtr);
   *sPtr<< jScalarNames_[i] << ": " << *jScalars_[i];
  
   
  }
  
   if(jScalars_.size()!=0)
   comma=true;
  
  //Write strings
  for(unsigned int i=0;i<jBools_.size();i++)
  {
    if(i!=0)
     *sPtr<<","<<endl<<endl;
    else
    {
     if(comma) *sPtr<<",";
     *sPtr<<endl<<endl;  
    }
   indentation(indent,sPtr);
   *sPtr<< jBoolNames_[i] << ": " << *jBools_[i];
  
   
  }
  
  
  *sPtr<<endl<<endl;
  indentation(indent-2,sPtr);
   *sPtr<<"}"<<endl;
}

/*-----------------------------------------------------------------------------*/
/*-----------------jsonFile Class----------------------------------------------*/
/*-----------------------------------------------------------------------------*/
jsonFile::jsonFile(std::string _dir, std::string _file)
:
dirName_(_dir),
fileName_(_file)
{
    //create directory
    int myProcID(Pstream::myProcNo());
    if(myProcID==0)
    {
        char buf[128];
        sprintf(buf,"mkdir %s",dirName_.c_str());
        system(buf);
    }
    sPtr = new OFstream(dirName_ + "/" + fileName_);
    mainObject.setIndent(2);
    mainObject.setPtr(sPtr);
}

/*-----------------------------------------------------------------------------*/
jsonFile::~jsonFile()
{
    delete sPtr;
    for(unsigned int n=0;n<jObjectNames_.size();n++)
        delete objList_[n];
}

/*-----------------------------------------------------------------------------*/
int jsonFile::lookupObject(std::string name) const
{
    for(unsigned int n=0;n<jObjectNames_.size();n++)
    {
        if(name==jObjectNames_[n])
        return n;
    }

    FatalErrorIn("ERROR") << "\n ERROR: the object named " << name << " could not be found in json file " << fileName_;

    return 1;
}

/*-----------------------------------------------------------------------------*/
void jsonFile::setFileName(std::string newname) const
{
    fileName_.assign(newname);

    // WARNING: rm: cannot remove ‘AutoDraGStatus.json’: No such file or directory might stem from here
    delete sPtr;

    sPtr = new OFstream(dirName_ + "/" + fileName_);
}

/*-----------------------------------------------------------------------------*/
void jsonFile::newObject(std::string objName, std::string rootName) const
{
    jsonObject* ob = new jsonObject;

    if(rootName=="mainObject")
    {
        mainObject.addjObject(objName,ob);
    }
    else
    {
        int id=this->lookupObject(rootName);

        objList_[id]->addjObject(objName,ob);
    }
    objList_.push_back(ob);
    jObjectNames_.push_back(objName);
}

/*-----------------------------------------------------------------------------*/
void jsonFile::addjVector(std::vector<double>* data , std::string dataname, std::string rootName) const
{
    if(rootName=="mainObject")
    {
        mainObject.addjVector(dataname,data);
    }
    else
    {
        int id=this->lookupObject(rootName);

        objList_[id]->addjVector(dataname, data);
    }
}

/*-----------------------------------------------------------------------------*/
void jsonFile::addjScalar(double* data , std::string dataname, std::string rootName) const
{
    if(rootName=="mainObject")
    {
        mainObject.addjScalar(dataname,data);
    }
    else
    {
        int id=this->lookupObject(rootName);

        objList_[id]->addjScalar(dataname, data);
    }
}

/*-----------------------------------------------------------------------------*/
void jsonFile::addjString(std::string* data , std::string dataname, std::string rootName) const
{
    if(rootName=="mainObject")
    {
        mainObject.addjString(dataname,data);
    }
    else
    {
        int id=this->lookupObject(rootName);

        objList_[id]->addjString(dataname, data);
    }
}

/*-----------------------------------------------------------------------------*/
void jsonFile::addjBool(bool* data , std::string dataname, std::string rootName) const
{
    if(rootName=="mainObject")
    {
        mainObject.addjBool(dataname,data);
    }
    else
    {
        int id=this->lookupObject(rootName);

        objList_[id]->addjBool(dataname, data);
    }   
}

/*-----------------------------------------------------------------------------*/
void jsonFile::write()
{
    //remove old file
    int myProcID(Pstream::myProcNo());
    if(myProcID==0)
    {
        char buf[128];
        sprintf(buf,"rm %s",fileName_.c_str());
        system(buf);
    }

    setFileName(fileName_);
    mainObject.write();
}
// ************************************************************************* //
