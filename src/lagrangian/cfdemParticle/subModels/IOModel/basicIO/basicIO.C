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
    and OpenFOAM. Note: this code is not part of OpenFOAM (see DISCLAIMER).
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
    dirName_("particles"),
    path_("dev/null")
{
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
        fileName lagPath=createTimeDir(path_);
        lagPath=createTimeDir(fileName(lagPath/"lagrangian"));

        // stream data to file
        streamDataToPath(lagPath, particleCloud_.positions(), particleCloud_.numberOfParticles(), "positions","vector","Cloud<passiveParticle>","0");
        streamDataToPath(lagPath, particleCloud_.velocities(), particleCloud_.numberOfParticles(), "v","vector","vectorField","");
        streamDataToPath(lagPath, particleCloud_.radii(), particleCloud_.numberOfParticles(), "r","scalar","scalarField","");
    }
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Private Member Functions

fileName basicIO::buildFilePath(word dirName) const
{
    // create file structure
    fileName path("."/dirName);
    mkDir(path,0777);
    mkDir(fileName(path/"constant"),0777);
    OFstream* stubFile = new OFstream(fileName(path/"particles.foam"));
    delete stubFile;
    return path;
}

void basicIO::streamDataToPath(fileName path, double** array,int n,word name,word type,word className,word finaliser) const
{
    vector position;
    OFstream* fileStream = new OFstream(fileName(path/name));
    *fileStream << "FoamFile\n";
    *fileStream << "{version 2.0; format ascii;class "<< className << "; location 0;object  "<< name <<";}\n";
    *fileStream << n <<"\n";
    *fileStream << "(\n";

    for(int index = 0;index < n; ++index)
    {
        if (type=="scalar"){
            *fileStream << array[index][0] << " \n";
        }else {
            for(int i=0;i<3;i++) position[i] = array[index][i];
            *fileStream <<"( "<< position[0] <<" "<<position[1]<<" "<<position[2]<<" ) "<< finaliser << " \n";
        }
    }
    *fileStream << ")\n";
    delete fileStream;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
