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

#ifndef OF_JSON_H
#define OF_JSON_H

#include <string>
#include <vector>
#include <fstream>
#include "OFstream.H"
#include "Pstream.H"


namespace Foam
{

 class jsonObject
 {
   friend class jsonFile;
  private: 
   
   mutable OFstream*                                           sPtr;
   mutable int                                               indent;
   mutable bool                                          lastObject;
   mutable std::vector<jsonObject*>                        jObjects_;
   mutable std::vector<std::string>                   jObjectNames_;
    
   mutable std::vector<std::vector<double>*>               jVectors_;
   mutable std::vector<std::string>                   jVectorNames_;
    
   mutable std::vector<double*>                            jScalars_;
   mutable std::vector<std::string>                   jScalarNames_;
    
   mutable std::vector<std::string*>                       jStrings_;
   mutable std::vector<std::string>                   jStringNames_;
       
   mutable std::vector<bool*>                                jBools_;
   mutable std::vector<std::string>                     jBoolNames_;
   
  
   jsonObject();
   ~jsonObject();
    
   void addjObject(std::string elemName, jsonObject* data) const;
   
   void addjVector(std::string elemName, std::vector<double>* data) const;
  void addjScalar(std::string elemName, double* data) const;
  void addjString(std::string elemName, std::string* data) const;
  void addjBool(std::string elemName, bool* data) const;
  
  void setIndent(int ind) const {indent=ind;};
  void setPtr(OFstream* ptr) const {sPtr=ptr;};
         
    
  void write() const;
  
  void lastObj() const {lastObject=true;};
 
  public:
     
   
 };
 
 class jsonFile
 {
 
  
  private:
  
   mutable std::string                   dirName_;
   mutable std::string                   fileName_;
   jsonObject                           mainObject;
   mutable OFstream*                          sPtr;
   
   mutable std::vector<jsonObject*>       objList_;
   mutable std::vector<std::string>  jObjectNames_;  
 
   int lookupObject(std::string name) const;
  public:
   
   jsonFile(std::string, std::string);
   ~jsonFile();
   
   void setFileName(std::string newname) const;
   
   
   jsonObject* getMainObject() {return &mainObject;};
   
   void write();
   
   void newObject(std::string objName, std::string rootName="mainObject") const;
   
   void addjVector( std::vector<double>* data , std::string dataname, std::string objName="mainObject") const;
   void addjScalar( double* data, std::string dataname, std::string objName="mainObject") const;
   void addjString(std::string* data, std::string dataname, std::string objName="mainObject") const;
   void addjBool(bool* data, std::string dataname, std::string objName="mainObject") const;
      
   
 };



}

/*--------USAGE EXAMPLE---------------------------

    jsonFile file("newDir","test.json");
   
    std::string teststring1("Is the");
    std::string teststring2("output working?");
    
    bool ans=true;
    
    double scal=9.9999;
     
    std::vector<double> vec2_;
    
    for(int i=0;i<7;i++)
     vec2_.push_back(i);
     
    file.addjScalar(&scal,"smaller than 10");
    
    file.newObject("test_object");
    
    file.addjVector(&vec2_,"nested_vector","test_object");
    
    file.newObject("nested_object","test_object");
    
    file.addjString(&teststring1,"teststring1","nested_object");
    file.addjString(&teststring2,"teststring2","nested_object");
    file.addjBool(&ans,"answer","nested_object");
    
    
    file.write();
    
-------THIS CODE WILL PRODUCE A JSON FILE-----------
{


  "test_object":
  {


    "nested_object":
    {


      "teststring1": "Is the",

      "teststring2": "output working?",

      "answer": 1

    },

    "nested_vector": [
    0,
    1,
    2,
    3,
    4,
    5,
    6
    ]

  },

  "smaller than 10": 9.9999

}
    



*/
#endif
