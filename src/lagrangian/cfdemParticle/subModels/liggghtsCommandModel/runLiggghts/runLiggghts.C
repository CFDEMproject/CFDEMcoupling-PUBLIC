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

#include "runLiggghts.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(runLiggghts, 0);

addToRunTimeSelectionTable
(
    liggghtsCommandModel,
    runLiggghts,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
runLiggghts::runLiggghts
(
    const dictionary& dict,
    cfdemCloud& sm,
    int i
)
:
    liggghtsCommandModel(dict,sm,i),
    propsDict_(dict),
    command_("run"),
    preNo_(true)
{
    word myName=word(typeName + "Props");
    if (dict.found(myName))    
    {
        propsDict_=dictionary(dict.subDict(myName));
        preNo_=Switch(propsDict_.lookup("preNo"));
    }

    runEveryCouplingStep_=true;

    strCommand_=createCommand(command_);

    checkTimeSettings(dict_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

runLiggghts::~runLiggghts()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const char* runLiggghts::command()
{
    return strCommand_.c_str();
}

string runLiggghts::createCommand(word command, word appendix, word appendix2, word appendix3, word appendix4)
{
    fileName add;
    char h[50];
    sprintf(h,"%d",particleCloud_.dataExchangeM().couplingInterval());
    add = h;
    command += " " + add + " " + appendix + " " + appendix2 + " " + appendix3 + " " + appendix4;

    return string(command);
}

bool runLiggghts::runCommand(int couplingStep)
{
    //change command to  "run xxx pre no"
    if (preNo_ && (couplingStep > firstCouplingStep_))
        strCommand_=createCommand(command_,"pre","no","post","no");

    return runThisCommand(couplingStep);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
