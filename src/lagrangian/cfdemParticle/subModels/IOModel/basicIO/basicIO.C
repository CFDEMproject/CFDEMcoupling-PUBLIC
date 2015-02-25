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

#include "basicIO.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(basicIO, 0);

addToRunTimeSelectionTable
(
    IOModel,
    basicIO,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
basicIO::basicIO
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    IOModel(dict,sm),
    //propsDict_(dict.subDict(typeName + "Props")),
    dirName_("lagrangian"),
    path_("dev/null"),
    nPProc_(-1),
    lagPath_("dev/null")
{
    //if (propsDict_.found("dirName")) dirName_=word(propsDict_.lookup("dirName"));
    path_ = buildFilePath(dirName_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

basicIO::~basicIO()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Public Member Functions

int basicIO::dumpDEMdata() const
{
    if (dumpNow())
    {
        // make time directory
        if (parOutput_) lagPath_=buildFilePath(dirName_);
        else
        {
            Info << "createTimeDir(path_), path="<<path_ << endl;
            Info << "lagPath_=createTimeDir(fileName(lagPath_/lagrangian)), lagPath="<<path_ << endl;
        	lagPath_=createTimeDir(path_);
        	lagPath_=createTimeDir(fileName(lagPath_/"lagrangian"));
        }
        // calc the number of particles on proc
        int count(0);
        for(int index = 0;index < particleCloud_.numberOfParticles(); ++index)
            if (particleCloud_.cellIDs()[index][0] > -1) count++;
        nPProc_=count;
        
        // stream data to file
        streamDataToPath(lagPath_, particleCloud_.positions(),nPProc_,"positions","vector","Cloud<passiveParticle>","0");
        streamDataToPath(lagPath_, particleCloud_.velocities(),nPProc_,"v","vector","vectorField","");
        streamDataToPath(lagPath_, particleCloud_.radii(),nPProc_,"r","scalar","scalarField","");
    }
    return nPProc_;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Private Member Functions


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
