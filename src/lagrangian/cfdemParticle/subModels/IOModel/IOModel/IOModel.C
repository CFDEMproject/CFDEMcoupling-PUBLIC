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
#include "IOModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(IOModel, 0);

defineRunTimeSelectionTable(IOModel, dictionary);

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void IOModel::dumpDEMdata() const
{}

fileName IOModel::createTimeDir(fileName path) const
{
    fileName timeDirPath(path/time_.timeName());
    mkDir(timeDirPath,0777);
    return timeDirPath;
}

fileName IOModel::createLagrangianDir(fileName path) const
{
    fileName lagrangianDirPath(path/"lagrangian");
    mkDir(lagrangianDirPath,0777);
    fileName cfdemCloudDirPath(lagrangianDirPath/"cfdemCloud1");
    mkDir(cfdemCloudDirPath,0777);
    return cfdemCloudDirPath;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
IOModel::IOModel
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    dict_(dict),
    particleCloud_(sm),
    time_(sm.mesh().time())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

IOModel::~IOModel()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
