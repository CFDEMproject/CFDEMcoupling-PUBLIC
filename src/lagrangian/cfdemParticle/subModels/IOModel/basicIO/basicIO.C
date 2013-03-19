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

#include "basicIO.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(basicIO, 0);

addToRunTimeSelectionTable
(
    IOModel,
    basicIO,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
basicIO::basicIO
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    IOModel(dict,sm),
    //propsDict_(dict.subDict(typeName + "Props")),
    dirName_("lagrangian"),
    path_("dev/null"),
    parOutput_(true),
    nPProc_(-1),
    lagPath_("dev/null")
{
	if (
            particleCloud_.dataExchangeM().myType()=="oneWayVTK" ||
            dict_.found("serialOutput")
       )
    {
        parOutput_=false;
        Warning << "IO model is in serial write mode, only data on proc 0 is written" << endl;
    }

    //if (propsDict_.found("dirName")) dirName_=word(propsDict_.lookup("dirName"));
    path_ = buildFilePath(dirName_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

basicIO::~basicIO()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Public Member Functions

void basicIO::dumpDEMdata() const
{
    if (time_.outputTime())
    {
        // make time directory
        if (parOutput_) lagPath_=buildFilePath(dirName_);
        else
        {
            Info << "createTimeDir(path_), path="<<path_ << endl;
            Info << "lagPath_=createTimeDir(fileName(lagPath_/lagrangian)), lagPath="<<path_ << endl;
        	lagPath_=createTimeDir(path_);
        	lagPath_=createTimeDir(fileName(lagPath_/"lagrangian"));
        }
        // calc the number of particles on proc
        int count(0);
        for(int index = 0;index < particleCloud_.numberOfParticles(); ++index)
            if (particleCloud_.cellIDs()[index][0] > -1) count++;
        nPProc_=count;
        
        // stream data to file
        streamDataToPath(lagPath_, particleCloud_.positions(), "positions","vector","Cloud<passiveParticle>","0");
        streamDataToPath(lagPath_, particleCloud_.velocities(), "v","vector","vectorField","");
        streamDataToPath(lagPath_, particleCloud_.radii(), "r","scalar","scalarField","");
    }
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Private Member Functions

fileName basicIO::buildFilePath(word dirName) const
{
    // create file structure
	fileName path("");
    if(parOutput_)
    {
    	path=fileName(particleCloud_.mesh().time().path()/particleCloud_.mesh().time().timeName()/dirName/"particleCloud");
    	mkDir(path,0777);
    } else
    {
		path=fileName("."/dirName);
    	mkDir(path,0777);
    	mkDir(fileName(path/"constant"),0777);
    	OFstream* stubFile = new OFstream(fileName(path/"particles.foam"));
    	delete stubFile;
        }
    return path;
}

void basicIO::streamDataToPath(fileName path, double** array,word name,word type,word className,word finaliser) const
{
    vector vec;
    OFstream* fileStream = new OFstream(fileName(path/name));
    *fileStream << "FoamFile\n";
    *fileStream << "{version 2.0; format ascii;class "<< className << "; location 0;object  "<< name <<";}\n";
    *fileStream << nPProc_ <<"\n";
    //*fileStream << "(\n";

    if(type!="origProcId")*fileStream << "(\n";
    else if(type=="origProcId")
    {
        if(nPProc_>0) *fileStream <<"{0}"<< "\n";
        else *fileStream <<"{}"<< "\n";
    }

    for(int index = 0;index < particleCloud_.numberOfParticles(); ++index)
    {
        if (particleCloud_.cellIDs()[index][0] > -1) // particle Found
        {
            if (type=="scalar"){
                *fileStream << array[index][0] << " \n";
            }else if (type=="position" || type=="vector"){
                for(int i=0;i<3;i++) vec[i] = array[index][i];
                *fileStream <<"( "<< vec[0] <<" "<<vec[1]<<" "<<vec[2]<<" ) "<< finaliser << " \n";
            }else if (type=="label"){
                *fileStream << index << finaliser << " \n";
            }
        }
    }
    //*fileStream << ")\n";
    if(type!="origProcId")*fileStream << ")\n";
    delete fileStream;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
