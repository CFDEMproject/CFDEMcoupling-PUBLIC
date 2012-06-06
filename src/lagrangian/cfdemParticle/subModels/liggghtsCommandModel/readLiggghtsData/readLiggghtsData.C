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

    // check if being run only at first coupling step
    runFirst_=Switch(propsDict_.lookup("runFirst"));

    if(!runFirst_)
    {
        // check if being run every coupling step
        runEveryCouplingStep_=Switch(propsDict_.lookup("runEveryCouplingStep"));

        scalar DEMts = particleCloud_.dataExchangeM().DEMts();
        scalar couplingInterval = particleCloud_.dataExchangeM().couplingInterval();

        if (!runEveryCouplingStep_)
        {
            // read time options
            startTime_ = readScalar(propsDict_.lookup("startTime"));
            endTime_ = readScalar(propsDict_.lookup("endTime"));
            timeInterval_ = readScalar(propsDict_.lookup("timeInterval"));

            // calculate coupling times
            firstCouplingStep_ = floor(startTime_/DEMts/couplingInterval);
            lastCouplingStep_ = floor(endTime_/DEMts/couplingInterval);
            couplingStepInterval_ = floor(timeInterval_/DEMts/couplingInterval);
        }
        else
        {
            firstCouplingStep_ =1;
            lastCouplingStep_ =1e9;
            couplingStepInterval_ =1;
        }
    }
    else
    {
            firstCouplingStep_ =1;
            lastCouplingStep_ =1;
            couplingStepInterval_ =1;
    }

    nextRun_ = firstCouplingStep_;

    Info << "firstCouplingStep = " << firstCouplingStep_ << endl;
    Info << "lastCouplingStep = " << lastCouplingStep_ << endl;
    Info << "couplingStepInterval = " << couplingStepInterval_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

readLiggghtsData::~readLiggghtsData()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
const char* readLiggghtsData::command()
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

    return runCommand(couplingStep);
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
