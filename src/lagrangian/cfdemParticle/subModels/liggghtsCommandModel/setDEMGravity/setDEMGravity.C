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

#include "setDEMGravity.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(setDEMGravity, 0);

addToRunTimeSelectionTable
(
    liggghtsCommandModel,
    setDEMGravity,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
setDEMGravity::setDEMGravity
(
    const dictionary& dict,
    cfdemCloud& sm,
    int i
)
:
    liggghtsCommandModel(dict,sm,i),
    propsDict_(dict),
    nrModel_(i),
    insertionNr_(0.),
    command_(0),
    filePathList_(100),
    scalarList_(40),
    labelList_(1),
    #if defined(version21) || defined(version16ext)
        g_(sm.mesh().lookupObject<uniformDimensionedVectorField> ("g")),
    #elif defined(version15)
        g_(dimensionedVector(sm.mesh().lookupObject<IOdictionary>("environmentalProperties").lookup("g")).value()),
    #endif
    unfix_(false)
{
    // define dictionary
    char h[80];
    sprintf(h,"%d",nrModel_);
    word myName_=word(typeName + "Props" + h);
    propsDict_=dictionary(dict.subDict(myName_));

    // read unfix flag
    unfix_=Switch(propsDict_.lookupOrDefault("unfix",false));

    if (propsDict_.found("exactTiming"))
        exactTiming_=true;
    Info << "exactTiming==" << exactTiming_ << endl;

    if (propsDict_.found("verbose")) verbose_=true;

    scalar gNorm=mag(g_.value());
    vector g=g_.value()/mag(g_.value());


    //=======================================
    // create command

    if(unfix_)
    {
        filePathList_[0] = "unfix";
        filePathList_[1] = "gravity";
    }
    else
    {
        filePathList_[0] = "fix";
        filePathList_[1] = "gravity";
        filePathList_[2] = "all";
        filePathList_[3] = "gravity";
        filePathList_[4] = "number";
        filePathList_[5] = "vector";
        filePathList_[6] = "number";
        filePathList_[7] = "number";
        filePathList_[8] = "number";
        scalarList_[0]=gNorm;
        scalarList_[1]=g[0];
        scalarList_[2]=g[1];
        scalarList_[3]=g[2];
    }

    command_=wordList(1);

    parseCommandList(filePathList_, labelList_, scalarList_, command_[0], propsDict_, timeStamp_);
    Info << "liggghtsCommand " << command_[0] << endl;

    checkTimeSettings(propsDict_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

setDEMGravity::~setDEMGravity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
const char* setDEMGravity::command(int commandLine)
{
Info << "commandLine=" << commandLine << endl;
    strCommand_=string(command_[commandLine]);
    
    return strCommand_.c_str();
}

bool setDEMGravity::runCommand(int couplingStep)
{
    checkTimeSettings(propsDict_);
    return runThisCommand(couplingStep);
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
