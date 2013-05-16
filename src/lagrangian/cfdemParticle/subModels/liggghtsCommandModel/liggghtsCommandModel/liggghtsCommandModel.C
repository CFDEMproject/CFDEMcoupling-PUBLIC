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
    couplingStepInterval_(0)
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
}

void liggghtsCommandModel::checkTimeSettings(const dictionary& propsDict)
{
    if(!runFirst_) //lastRun or runEveryCouplingStep or every n steps or every writeStep
    {
        scalar DEMts = particleCloud_.dataExchangeM().DEMts();
        scalar couplingInterval = particleCloud_.dataExchangeM().couplingInterval();

        if(runLast_) // last run
        {
            // read time options from subdict
            endTime_ = particleCloud_.mesh().time().endTime().value();
            startTime_ = endTime_;
            timeInterval_ = 1;

            // calculate coupling times
            firstCouplingStep_ = floor(startTime_/DEMts/couplingInterval);
            lastCouplingStep_ = floor(endTime_/DEMts/couplingInterval);
            couplingStepInterval_ = floor(timeInterval_/DEMts/couplingInterval);
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
                firstCouplingStep_ = floor(startTime_/DEMts/couplingInterval)+1;
                lastCouplingStep_ = floor(endTime_/DEMts/couplingInterval)+1;
                couplingStepInterval_ = floor(timeInterval_/DEMts/couplingInterval)+1;
            }
            else      //runEveryCouplingStep  or writeStep
            {
                firstCouplingStep_ =1;
                lastCouplingStep_ =1e9;
                couplingStepInterval_ =1;
            }
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

bool liggghtsCommandModel::runThisCommand(int couplingStep)
{
    bool runIt=false;
    if(
       (!runEveryWriteStep_ && firstCouplingStep_  <= couplingStep && lastCouplingStep_  >= couplingStep)  ||
       (runEveryWriteStep_  && particleCloud_.mesh().time().outputTime())
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
