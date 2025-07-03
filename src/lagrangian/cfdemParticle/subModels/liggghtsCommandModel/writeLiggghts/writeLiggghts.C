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
    writeLastOnly_(false),
    overwrite_(false),
    writeEvery_(1),
    counter_(-1)
{
    if (dict.found(typeName + "Props"))    
    {
        Info << "   Found a dictionary " << word(typeName + "Props") << ", reading user defined values (if specified):" << endl;
        propsDict_=dictionary(dict.subDict(typeName + "Props"));

        // check if verbose
        if (propsDict_.found("verbose")) verbose_=true;

        if(propsDict_.found("writeLastOnly"))
        {
            writeLastOnly_=Switch(propsDict_.lookup("writeLastOnly"));
        }
        if(propsDict_.found("writeLast")) FatalError<<"You are using the old style variable name <writeLast> - this keyword has changed to <writeLastOnly>." << abort(FatalError);

        if (propsDict_.found("path"))
        {
            path_=fileName(propsDict_.lookup("path"));
            if (!checkPath(path_))
                FatalError<<"The path you provided in writeLiggghtsProps is incorrect: " << path_ << "\n" << abort(FatalError);
            else
                Info << "   Using user defined path to write LIGGGHTS restart file: " << path_ << endl;
        }

        if(propsDict_.found("writeName"))
        {
            propsDict_.lookup("writeName") >> writeName_;
        }

        if(propsDict_.found("writeEvery"))
        {
            propsDict_.lookup("writeEvery") >> writeEvery_;
            if(writeEvery_<1) FatalError << "writeEvery must be chosen > 0" << abort(FatalError);
        }

        overwrite_=Switch(propsDict_.lookup("overwrite"));
    }
    else
    {
        Info << "   Did not find a dictionary " << word(typeName + "Props") << ", using default values:" << endl;
    }

    if(writeLastOnly_)
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

    // output for user
    Info << "     overwrite = " << Switch(overwrite_) << endl;
    Info << "     path = " << path_ << endl;
    if(overwrite_) Info << "     writeName = " << writeName_ << " (Restart file is written to " << command_ << " ." << endl;
    else Info << "     writeName = " << writeName_ << " (Restart file is written to " << command_ << "<timeStamp> ." << endl;
    Info << "     writeLastOnly = " << writeLastOnly_ << endl;
    Info << "     writeEvery = " << writeEvery_ << " (Restart file is written every " << writeEvery_ 
         << " CFD's writeInterval.)" << endl;
    Info << "     verbose = " << Switch(verbose_) << endl;

    checkTimeSettings(dict_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

writeLiggghts::~writeLiggghts()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const char* writeLiggghts::command(int commandLine)
{
    return strCommand_.c_str();
}

bool writeLiggghts::runCommand(int couplingStep)
{
    if(!overwrite_) strCommand_=addTimeStamp(command_);
    if(particleCloud_.writeTimePassed()) counter_ = counter_+1;
    label writeInterval = counter_ % writeEvery_;
    bool a=runThisCommand(couplingStep); // needs to be called as it triggers resetWriteTimePassed()
    if(writeInterval == 0) return a;
    else return 0;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
