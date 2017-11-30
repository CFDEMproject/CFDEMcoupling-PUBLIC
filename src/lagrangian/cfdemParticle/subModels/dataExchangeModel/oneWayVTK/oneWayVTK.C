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

#include "oneWayVTK.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(oneWayVTK, 0);

addToRunTimeSelectionTable
(
    dataExchangeModel,
    oneWayVTK,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
oneWayVTK::oneWayVTK
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    dataExchangeModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    filename_(propsDict_.lookup("couplingFilename")),
    relativePath_(propsDict_.lookup("relativePath"))
{
    readDEMtsfromDict(propsDict_);

    // set max nr of particles from dict
    maxNumberOfParticles_ = readScalar(propsDict_.lookup("maxNumberOfParticles"));
    setNumberOfParticles(maxNumberOfParticles_);

    // make a const char* from word
    //string HH=string(filename_);
    //charFilename_=const_cast<char*>(HH.c_str());
    charFilename_ = wordToChar(filename_);

    Info << "relativePath_" << relativePath_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

oneWayVTK::~oneWayVTK()
{}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
char* oneWayVTK::wordToChar(word& inWord) const
{
    return const_cast<char*>(inWord.c_str());
}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void oneWayVTK::getData
(
    word name,
    word type,
    double ** const& field,
    int step
) const
{
    if (type == "scalar-atom")
    {
        // get path to particle VTK files
        char index[100];
        sprintf(index, charFilename_, step);
        //fileName H(particleCloud_.mesh().time().path()/".."/"DEM"/"post"/index);
        fileName H(particleCloud_.mesh().time().path()/relativePath_/index);
        Info << "opening file: " << H << endl;

        // set file pointer
        string HH=string(H);
        const char * paricleFilePath=HH.c_str();
        ifstream* inputPtr;
        inputPtr = new ifstream(paricleFilePath);
        if(!*inputPtr) FatalError << "File not found!, " << H << "\n" << abort(FatalError);

        if (name == "radius")
        {
            // read data
            string just_read = " ";
            while(just_read.compare(name) != 0)  *inputPtr >> just_read;   //read until we read "name"
            *inputPtr >> just_read;                                    // skip text for dataType
            *inputPtr >> just_read;                                    // skip text for "1"
            *inputPtr >> just_read;                                    // skip text for "LookUp"
            *inputPtr >> just_read;                                    // skip text for "default"
            for(int index = 0;index <  particleCloud_.numberOfParticles(); ++index)
            {
                *inputPtr >> field[index][0];
            }
        }
        else
        {
            // read data
            string just_read = " ";
            while(just_read.compare(name) != 0)  *inputPtr >> just_read;   //read until we read "name"
            *inputPtr >> just_read;                                    // skip text for dataType
            for(int index = 0;index <  particleCloud_.numberOfParticles(); ++index)
            {
                *inputPtr >> field[index][0];
            }
        }

        // clean up inputStream
        delete inputPtr;
    } else if (type == "vector-atom")
    {
        // get path to particle VTK files
        char index[100];
        sprintf(index, charFilename_, step);
        Info << "debug: index is " << index << endl; //JOKER
        //fileName H(particleCloud_.mesh().time().path()/".."/"DEM"/"post"/index);
        fileName H(particleCloud_.mesh().time().path()/relativePath_/index);
        Info << "opening file: " << H << endl;

        // set file pointer
        string HH=string(H);
        const char * paricleFilePath=HH.c_str();
        ifstream* inputPtr;
        inputPtr = new ifstream(paricleFilePath);
        if(!*inputPtr) FatalError << "File not found!, " << H << "\n" << abort(FatalError);

        // read position data from VTK file
        //NP: secial case as position data has no "name" in the vtk file
        if (name == "x")
        {
            int numberOfParticles;  // remove this?

            string just_read = " ";
//            if(!*inputPtr) FatalIOError << "File not found!, " << H << "\n" << abort(FatalError);
            while(just_read.compare("POINTS") != 0)  *inputPtr >> just_read;   //read until we read "POINTS"
            *inputPtr >> numberOfParticles;                           //this is now the number of points in the file
            *inputPtr >> just_read;                                    // skip text for dataType

            // give nr of particles to cloud
            setNumberOfParticles(numberOfParticles);

            // re-allocate arrays of cloud
            particleCloud_.reAllocArrays();

            for(int index = 0;index <  numberOfParticles; ++index)
            {
                *inputPtr >> field[index][0] >> field[index][1] >> field[index][2];
            }
        }
        else
        {
            string just_read = " ";
            while(just_read.compare(name) != 0)  *inputPtr >> just_read;   //read until we read "name"
            *inputPtr >> just_read;                                    // skip "3"
            *inputPtr >> just_read;                                    // skip nr entries
            *inputPtr >> just_read;                                    // skip text for dataType
            for(int index = 0;index <  particleCloud_.numberOfParticles(); ++index)
            {
                *inputPtr >> field[index][0] >> field[index][1] >> field[index][2];
            }
        }

        // clean up inputStream
        delete inputPtr;
    }
    else
    {
        Info << "unknown type in getData!!!" << endl;
    }
}

void oneWayVTK::giveData
(
    word name,
    word type,
    double ** const& field,
    const char* datatype
) const
{
    // do nothing
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
