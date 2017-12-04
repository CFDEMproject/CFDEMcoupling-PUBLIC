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
#include <sys/stat.h>
#include "liggghtsCommandModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(liggghtsCommandModel, 0);

defineRunTimeSelectionTable(liggghtsCommandModel, dictionary);

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
liggghtsCommandModel::liggghtsCommandModel
(
    const dictionary& dict,
    cfdemCloud& sm,
    int i
)
:
    dict_(dict),
    particleCloud_(sm),
    strCommand_("notDefined"),
    nextRun_(-1),
    lastRun_(-1),
    runFirst_(false),
    runLast_(false),
    runEveryCouplingStep_(false),
    runEveryWriteStep_(false),
    startTime_(-1.),
    endTime_(-1.),
    timeInterval_(0.),
    firstCouplingStep_(-1),
    lastCouplingStep_(-1),
    couplingStepInterval_(0),
    exactTiming_(false),
    commandLines_(1),
    verbose_(false),
    timeStamp_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

liggghtsCommandModel::~liggghtsCommandModel()
{}

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //

void liggghtsCommandModel::checkTimeMode(dictionary& propsDict)
{
    runFirst_=Switch(propsDict.lookup("runFirst"));

    if(!runFirst_)
    {
        // check if being run only at last coupling step
        runLast_=Switch(propsDict.lookup("runLast"));

        if(!runLast_)
        {
            // check if being run every coupling step
            runEveryCouplingStep_=Switch(propsDict.lookup("runEveryCouplingStep"));

            if(!runEveryCouplingStep_)
            {
                runEveryWriteStep_=Switch(propsDict.lookup("runEveryWriteStep"));
            }
        }
    }
    if(verbose_){
        Info << "runFirst = " << runFirst_ << endl;
        Info << "runLast = " << runLast_ << endl;
        Info << "runEveryCouplingStep = " << runEveryCouplingStep_ << endl;
        Info << "runEveryWriteStep = " << runEveryWriteStep_ << endl;
    }
}

void liggghtsCommandModel::checkTimeSettings(const dictionary& propsDict)
{
    if(!runFirst_) //lastRun or runEveryCouplingStep or every n steps or every writeStep
    {
        scalar DEMts = particleCloud_.dataExchangeM().DEMts();
        scalar couplingInterval = particleCloud_.dataExchangeM().couplingInterval();
        scalar simStartTime = particleCloud_.mesh().time().startTime().value();

        if(runLast_) // last run
        {
            // read time options from subdict
            endTime_ = particleCloud_.mesh().time().endTime().value()-simStartTime;
            startTime_ = endTime_;
            timeInterval_ = -1;

            // calculate coupling times
            firstCouplingStep_ = floor((startTime_+SMALL)/DEMts/couplingInterval);
            lastCouplingStep_ = floor((endTime_+SMALL)/DEMts/couplingInterval);
            couplingStepInterval_ = -1;
        }
        else         //runEveryCouplingStep of every n steps or every writeStep
        {
            if (!runEveryCouplingStep_ && !runEveryWriteStep_) // every n steps
            {
                // read time options from subdict
                startTime_ = readScalar(propsDict.lookup("startTime"));
                endTime_ = readScalar(propsDict.lookup("endTime"));
                timeInterval_ = readScalar(propsDict.lookup("timeInterval"));

                // calculate coupling times
                firstCouplingStep_ = floor((startTime_+SMALL-simStartTime)/DEMts/couplingInterval);
                lastCouplingStep_ = floor((endTime_+SMALL-simStartTime)/DEMts/couplingInterval);
                couplingStepInterval_ = floor(timeInterval_+SMALL/DEMts/couplingInterval);
            }
            else      //runEveryCouplingStep  or writeStep
            {
                firstCouplingStep_ =1;
                lastCouplingStep_ =1e9;
                couplingStepInterval_ =1;
            }
        }
    }
    else // runFirst
    {
            firstCouplingStep_ =1;
            lastCouplingStep_ =1;
            couplingStepInterval_ =-1;
    }

    if(verbose_){
        Info << "firstCouplingStep = " << firstCouplingStep_ << endl;
        Info << "lastCouplingStep = " << lastCouplingStep_ << endl;
        Info << "couplingStepInterval = " << couplingStepInterval_ << endl;
        Info << "nextRun = " << nextRun_ << endl;
    }
}

bool liggghtsCommandModel::runThisCommand(int couplingStep)
{
    if(verbose_) Info << "couplingStep = " << couplingStep << endl;
    bool runIt=false;
    if(
       (!runEveryWriteStep_ && firstCouplingStep_  <= couplingStep && lastCouplingStep_  >= couplingStep)  ||
       (runEveryWriteStep_  && particleCloud_.writeTimePassed())
      )
    {
        if(couplingStep >= nextRun_)
        {
            runIt=true;
            nextRun_=couplingStep + couplingStepInterval_;
            lastRun_=couplingStep;
        }
    }
    return runIt;
}

string liggghtsCommandModel::addTimeStamp(word command)
{
    char add[100];
    sprintf(add,"%f",particleCloud_.mesh().time().value());

    return string(command+add);
}

DynamicList<scalar> liggghtsCommandModel::executionsWithinPeriod(scalar TSstart,scalar TSend)
{
    Info << "liggghtsCommandModel::executionsWithinPeriod TSstart" << TSstart << endl;
    Info << "liggghtsCommandModel::executionsWithinPeriod TSend" << TSend << endl;
    Info << "startTime =" << startTime_ << endl;
    Info << "endTime =" << endTime_ << endl;

    // init exec times array
    DynamicList<scalar> executions(0);

    // current TS within active period
    if(startTime_+SMALL<TSend && endTime_>TSstart-SMALL )
    {
        Info << "working time within this TS" << endl;

        // find first execution within TS (better routine with modulo)
        int startNr = 0;
        scalar t = startTime_ + startNr * timeInterval_;

        if(timeInterval_ > SMALL)
        {
            while (TSend - t > SMALL)
            {
                t = startTime_ + startNr * timeInterval_;
                startNr++;
            }
            t -= timeInterval_;
        }
        // check if first exec found within TS
        if(TSstart < t + SMALL && t +SMALL < TSend)
        {
            // check for more executions
            while (t < endTime_ + SMALL && TSend - t > SMALL)
            {
                executions.append(t);
                t += timeInterval_;
            }
        }
        //else
        //    Info << "liggghtsCommandModel::executionsWithinPeriod error???" << endl;     

        // debug
        Info << "liggghtsCommandModel::executionsWithinPeriod executions=" << executions << endl;
    }

    // return dummy

    return executions;
}

bool liggghtsCommandModel::checkPath(fileName path)
{
    string strPath = string(path);
    struct stat buffer;   
    return (stat (strPath.c_str(), &buffer) == 0); 
}


void liggghtsCommandModel::parseCommandList(wordList& commandList,labelList& labelList,
                                            scalarList& scalarList,
                                            word& command, dictionary& propsDict, 
                                            bool& timeStamp)
{
    bool addBlank = true;  // std no blanks after each word
    fileName add;
    label numberCount=0;   // nr of scalars inserted to command
    label labelCount=0;   // nr of labels inserted to command

    forAll(commandList,i)
    {
        add = word(commandList[i]);

        //- handle symbols
        if (add == "$couplingInterval")
        {
            char h[50];
            sprintf(h,"%d",particleCloud_.dataExchangeM().couplingInterval());
            add = h;
        }
        else if (add=="dot")    add = ".";
        else if (add=="dotdot") add = "..";
        else if (add=="slash")  add = "/";
        else if (add=="noBlanks")  // no blanks after the following words
        {
            add = "";
            addBlank = false;
        }else if (add=="blanks") // add a blank here and after the following words
        {
            add = "";
            addBlank = true;
        }else if (add=="timeStamp") // next command will be a number read from labelList
        {
            add = "";
            timeStamp=true;
        }else if (add=="number") // next command will be a number read from labelList
        {
            /*if (!propsDict.found("scalars"))
            {
                FatalError<<"you want to use a number in the command\n - specify a scalar list with all numbers"
                << abort(FatalError);
            }*/
            char h[50];
            sprintf(h,"%f",scalarList[numberCount]);
            add = h;
            numberCount ++;
        }else if (add=="label") // next command will be a number read from labelList
        {
            /*if (!propsDict.found("labels"))
            {
                FatalError<<"you want to use a label in the command\n - specify a label list with all numbers"
                << abort(FatalError);
            }*/
            char h[50];
            sprintf(h,"%d",labelList[labelCount]);
            add = h;
            labelCount ++;
        }

        // compose command
        if (addBlank)
        {
            command += add + " ";
        }else
        {
            command += add;
        }
    }
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
