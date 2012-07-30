/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright (C) 1991-2009 OpenCFD Ltd.
                                Copyright (C) 2009-2012 JKU, Linz
                                Copyright (C) 2012-     DCS Computing GmbH,Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling.  If not, see <http://www.gnu.org/licenses/>.

Application
    writeUfluidwriteUfluid

Description
    Writes the the cell center fluid velocity to particles in the lagrangian 
    time directory.
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Cloud.H"
#include "IOdictionary.H"
#include "fvCFD.H"
#include "fvMesh.H"
#include "Time.H"
#include "timeSelector.H"
#include "OFstream.H"
#include "passiveParticleCloud.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "addRegionOption.H"

    #include "setRootCase.H"

    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int nParticle=0;
 forAll(timeDirs, timeI)
 {
    runTime.setTime(timeDirs[timeI], timeI);
    Info<< "Time = " << runTime.timeName() << endl;
    IOobject UHeader
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );
    if (UHeader.headerOk())
    {
        volVectorField U(UHeader,mesh);
        passiveParticleCloud myCloud(mesh, cloudName);
		myCloud.write();
        nParticle = myCloud.size();
	    IOField<vector> Ufluid(myCloud.fieldIOobject("Ufluid",IOobject::NO_READ),nParticle);
        label i = 0;
        forAllConstIter(passiveParticleCloud, myCloud, iter)
        {
		    Ufluid[i]=U[iter().cell()];
            i++;
        }
        Ufluid.write();
    }
    else
    {
        Info << "velocity field not found" << endl;
    }
 }
    return 0;
}


// ************************************************************************* //
