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

#include "trackIO.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(trackIO, 0);

addToRunTimeSelectionTable
(
    IOModel,
    trackIO,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
trackIO::trackIO
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    IOModel(dict,sm),
    //propsDict_(dict.subDict(typeName + "Props")),
    dirName_(""),
    path_("dev/null"),
    lagPath_("dev/null")
{
    //if (propsDict_.found("dirName")) dirName_=word(propsDict_.lookup("dirName"));
    path_ = buildFilePath(dirName_);

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

trackIO::~trackIO()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Public Member Functions

void trackIO::dumpDEMdata() const
{
    if (time_.outputTime())
    {
        // make time directory
        lagPath_=createTimeDir(path_);
        lagPath_=createLagrangianDir(fileName(lagPath_));

        // stream data to file
        streamDataToPath(lagPath_, particleCloud_.positions(), particleCloud_.numberOfParticles(), "positions","position","Cloud<passiveParticle>","0");
        streamDataToPath(lagPath_, particleCloud_.velocities(), particleCloud_.numberOfParticles(), "v","vector","vectorField","");
        streamDataToPath(lagPath_, particleCloud_.velocities(), particleCloud_.numberOfParticles(), "origId","label","labelField","");
        streamDataToPath(lagPath_, particleCloud_.velocities(), particleCloud_.numberOfParticles(), "origProcId","origProcId","labelField","");
        streamDataToPath(lagPath_, particleCloud_.radii(), particleCloud_.numberOfParticles(), "r","scalar","scalarField","");
    }
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Private Member Functions

fileName trackIO::buildFilePath(word dirName) const
{
    // create file structure
    fileName path("."/dirName);
    OFstream* stubFile = new OFstream(fileName(path/"particles.foam"));
    delete stubFile;
    return path;
}

void trackIO::streamDataToPath(fileName path, double** array,int n,word name,word type,word className,word finaliser) const
{
    vector vec;
    OFstream* fileStream = new OFstream(fileName(path/name));
    *fileStream << "FoamFile\n";
    *fileStream << "{version 2.0; format ascii;class "<< className << "; location 0;object  "<< name <<";}\n";
    *fileStream << n <<"\n";
    if(type!="origProcId")*fileStream << "(\n";
    else if(type=="origProcId")*fileStream <<"{0}"<< "\n";

    for(int index = 0;index < n; ++index)
    {
        if (type=="scalar"){
            *fileStream << array[index][0] << finaliser << " \n";
        }else if (type=="position"){
            for(int i=0;i<3;i++) vec[i] = array[index][i];
//          You may need to use these two lines if you have cyclics
//		    if(vec[0]<0)vec[0]+=0.12;if(vec[0]>0.12)vec[0]-=0.12;
//		    if(vec[1]<0)vec[1]+=0.06;if(vec[1]>0.06)vec[1]-=0.06;
            *fileStream <<"( "<< vec[0] <<" "<<vec[1]<<" "<<vec[2]<<" ) "<< finaliser << " \n";
        }else if (type=="label"){
            *fileStream << index << finaliser << " \n";
        }else  if (type=="vector"){
            for(int i=0;i<3;i++)  vec[i] = array[index][i]; 
            *fileStream <<"( "<< vec[0] <<" "<<vec[1]<<" "<<vec[2]<<" ) " << finaliser << " \n";
        }
    }
    if(type!="origProcId")*fileStream << ")\n";
    delete fileStream;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
