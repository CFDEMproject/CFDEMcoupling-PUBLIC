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

#include "readLiggghtsData.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(readLiggghtsData, 0);

addToRunTimeSelectionTable
(
    liggghtsCommandModel,
    readLiggghtsData,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
readLiggghtsData::readLiggghtsData
(
    const dictionary& dict,
    cfdemCloud& sm,
    int i
)
:
    liggghtsCommandModel(dict,sm,i),
    nrModel_(i),
    insertionNr_(0.),
    command_("read_data"),
    myName_("notYetGiven"),
    propsDict_(dict),
    filePathList_(0)
{
    // define dictionary
    char h[80];
    sprintf(h,"%d",nrModel_);
    myName_=word(typeName + "Props" + h);
    propsDict_=dictionary(dict.subDict(myName_));

    if (propsDict_.found("exactTiming"))
        exactTiming_=true;
    Info << "exactTiming==" << exactTiming_ << endl;

    if (propsDict_.found("verbose")) verbose_=true;

    // read first index of data file to be injected
    insertionNr_=readScalar(propsDict_.lookup("startIndex"));

    // read command from dict
    filePathList_ = wordList(propsDict_.lookup("filePath"));

    command_ += " ";
    forAll(filePathList_,i)
    {
        fileName add = filePathList_[i];

        //- handle symbols
        if (add=="dot")    add = ".";
        else if (add=="dotdot") add = "..";
        else if (add=="slash")  add = "/";

        // compose command
        command_ += add;
    }

    checkTimeMode(propsDict_);

    checkTimeSettings(propsDict_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

readLiggghtsData::~readLiggghtsData()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
const char* readLiggghtsData::command(int commandLine)
{
    char h[50];
    sprintf(h,"_%d",insertionNr_);
    word add = h;
    insertionNr_++;
    strCommand_=string(command_ + add + " add");

    return strCommand_.c_str();
}

bool readLiggghtsData::runCommand(int couplingStep)
{
    checkTimeSettings(propsDict_);
    return runThisCommand(couplingStep);
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
