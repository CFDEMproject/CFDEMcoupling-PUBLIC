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

#include "writeLiggghts.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(writeLiggghts, 0);

addToRunTimeSelectionTable
(
    liggghtsCommandModel,
    writeLiggghts,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
writeLiggghts::writeLiggghts
(
    const dictionary& dict,
    cfdemCloud& sm,
    int i
)
:
    liggghtsCommandModel(dict,sm,i),
    propsDict_(dict),
    command_("write_restart"),
    path_(word("..")/word("DEM")),
    writeName_("liggghts.restartCFDEM"),
    writeLast_(true),
    overwrite_(true)
{
    if (dict.found(typeName + "Props"))    
    {
        propsDict_=dictionary(dict.subDict(typeName + "Props"));

        if(propsDict_.found("writeLast"))
        {
            writeLast_=Switch(propsDict_.lookup("writeLast"));
        }

        if (!writeLast_ && propsDict_.found("overwrite"))
        {
            overwrite_=Switch(propsDict_.lookup("overwrite"));
        }
    }
    if(writeLast_)
        runLast_=true;
    else
    {
        //Warning << "Using invalid options of writeLiggghts, please use 'writeLast' option." << endl;
        runEveryWriteStep_=true;
    }    


    command_ += " " + path_ + "/" + writeName_;
    if(overwrite_)
        strCommand_=string(command_);
    else
        command_ += "_";

    Info << "writeLiggghts: A restart file writeName_"<< command_ <<" is written." << endl;

    checkTimeSettings(dict_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

writeLiggghts::~writeLiggghts()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const char* writeLiggghts::command()
{
    return strCommand_.c_str();
}

bool writeLiggghts::runCommand(int couplingStep)
{
    if(!overwrite_) strCommand_=addTimeStamp(command_);
    return runThisCommand(couplingStep);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
