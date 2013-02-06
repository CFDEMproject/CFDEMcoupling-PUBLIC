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

#include "twoWayFiles.H"
#include "addToRunTimeSelectionTable.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(twoWayFiles, 0);

addToRunTimeSelectionTable
(
    dataExchangeModel,
    twoWayFiles,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
twoWayFiles::twoWayFiles
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    dataExchangeModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props"))
{
    readDEMtsfromDict(propsDict_);

    // set max nr of particles from dict
    maxNumberOfParticles_ = readScalar(propsDict_.lookup("maxNumberOfParticles"));

    // give max nr of particles to cloud (corrected later)
    setNumberOfParticles(maxNumberOfParticles_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

twoWayFiles::~twoWayFiles()
{}


// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //

const char* twoWayFiles::wordToChar(word& inWord) const
{
    string HH = string(inWord);
    return HH.c_str();
}

const char* twoWayFiles::fileNameToChar(fileName& inWord) const
{
    string HH = string(inWord);
    return HH.c_str();
}

fileName twoWayFiles::getFilePath(word& name, bool in) const
{
    const char* charName = wordToChar(name);
    char timeStep[40];

    // file touched by DEM
    strcpy(timeStep, charName);
    strcat(timeStep,"1");
    fileName particleFilePathOld(particleCloud_.mesh().time().path()/"couplingFiles"/timeStep);

    //NP no waiting when writing out at first time
    if (couplingStep() > 1 || in)
    {
        Info << "wait for file " << particleFilePathOld << endl;
        struct stat st;
        while (stat(fileNameToChar(particleFilePathOld),&st)) sleep(0.03);
    }

    return particleFilePathOld;
}

void twoWayFiles::renameFilePath(fileName& particleFilePathOld,word& name) const
{
    const char* charName = wordToChar(name);
    char timeStep[40];

    // file touched by CFD
    strcpy(timeStep, charName);
    strcat(timeStep,"0");
    fileName particleFilePath(particleCloud_.mesh().time().path()/"couplingFiles"/timeStep);

    // rename old file
    rename(fileNameToChar(particleFilePathOld),fileNameToChar(particleFilePath));
}

// * * * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * * //
void twoWayFiles::getData
(
    word name,
    word type,
    double ** const& field,
    label step
) const
{
    // get input path
    fileName particleFilePath = getFilePath(name,true);
    Info << "reading from file: " << particleFilePath << endl;

    // set file pointer
    IFstream* inputPtr = new IFstream(particleFilePath);

    // write data to variable
    int numberOfParticles;
    /*if(name != "outRegion1" && name != "inRegion1")*/ *inputPtr >> numberOfParticles;

    // give nr of particles to cloud
    setNumberOfParticles(numberOfParticles);

    // re-allocate arrays of cloud
    particleCloud_.reAllocArrays();

    for(int index = 0;index < numberOfParticles; ++index)
    {
        if (type == "scalar-atom")
        {
            *inputPtr >> field[index][0];
        }
        else if (type == "vector-atom")
        {
            for(int i=0;i<3;i++) *inputPtr >> field[index][i];
        }
        else
        {
            FatalError<<"unknown type in twoWayFiles::getData!!!\n" << abort(FatalError);
        }
    }

    // clean up inputStream
    delete inputPtr;

    // rename file
    renameFilePath(particleFilePath,name);
}

void twoWayFiles::giveData
(
    word name,
    word type,
    double ** const& field,
    const char* datatype
) const
{
    // get output path
    fileName particleFilePath = getFilePath(name,false);
    Info << "writing to file: " << particleFilePath << endl;

    // set file pointer
    OFstream* outputPtr = new OFstream(particleFilePath);

    // write data to file
    int numberOfParticles = particleCloud_.numberOfParticles();
    *outputPtr << numberOfParticles << endl;

    for(int index = 0;index < numberOfParticles; ++index)
    {
        if (type == "scalar-atom")
        {
            *outputPtr << field[index][0] << endl;
        }
        else if (type == "vector-atom")
        {
            for(int i=0;i<3;i++) *outputPtr << field[index][i] << " ";
            *outputPtr << endl;
        }
        else
        {
            FatalError<<"unknown type in twoWayFiles::giveData!!!\n" << abort(FatalError);
        }
    }

    // clean up outputStream
    delete outputPtr;

    // rename file
    renameFilePath(particleFilePath,name);

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
