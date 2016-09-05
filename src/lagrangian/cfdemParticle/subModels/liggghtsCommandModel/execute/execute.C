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

#include "execute.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(execute, 0);

addToRunTimeSelectionTable
(
    liggghtsCommandModel,
    execute,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
execute::execute
(
    const dictionary& dict,
    cfdemCloud& sm,
    int i
)
:
    liggghtsCommandModel(dict,sm,i),
    nrModel_(i),
    myName_("notYetGiven"),
    propsDict_(dict),
    commandList_(0),
    command_(""),
    scalarList_(0),
    labelList_(0)
{
    // define dictionary
    char h[80];
    sprintf(h,"%d",nrModel_);
    myName_=word(typeName + "Props" + h);
    propsDict_=dictionary(dict.subDict(myName_));

    // read command from dict
    commandList_ = wordList(propsDict_.lookup("command"));

    // read list of scalars
    if(propsDict_.found("scalars")) scalarList_ = scalarList(propsDict_.lookup("scalars"));

    // read list of labels
    if(propsDict_.found("labels")) labelList_ = labelList(propsDict_.lookup("labels"));

    // check if verbose
    if (propsDict_.found("verbose")) verbose_=true;

    parseCommandList(commandList_, labelList_, scalarList_, command_, propsDict_, timeStamp_);
    Info << "liggghtsCommand " << command_ << endl;

    strCommand_=string(command_);

    checkTimeMode(propsDict_);

    checkTimeSettings(propsDict_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

execute::~execute()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const char* execute::command(int commandLine)
{
    return strCommand_.c_str();
}

bool execute::runCommand(int couplingStep)
{
    if(timeStamp_) strCommand_=addTimeStamp(command_);
    checkTimeSettings(propsDict_);

    return runThisCommand(couplingStep);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
