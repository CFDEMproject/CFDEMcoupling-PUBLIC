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
    cfdemPostproc

Description
    Tool for DEM->CFD (Lagrange->Euler) mapping to calculate local voidfraction
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"

#include "cfdemCloud.H"
#include "dataExchangeModel.H"
#include "voidFractionModel.H"
#include "locateModel.H"
#include "averagingModel.H"
#include "momCoupleModel.H"
#include "forceModel.H"
#include "IOModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // create cfdemCloud
    cfdemCloud particleCloud(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    int count=0;
    int DEM_dump_Interval=1000;
    particleCloud.reAllocArrays();

    double **positions_;
    double **velocities_;
    double **radii_;
    double **voidfractions_;
    double **particleWeights_;
    double **particleVolumes_;
    double **cellIDs_;
    
    particleCloud.dataExchangeM().allocateArray(positions_,0.,3);
    particleCloud.dataExchangeM().allocateArray(velocities_,0.,3);
    particleCloud.get_radii(radii_);  // get ref to radii
    //particleCloud.dataExchangeM().allocateArray(radii_,0.,1);
    particleCloud.dataExchangeM().allocateArray(voidfractions_,0.,1);
    particleCloud.dataExchangeM().allocateArray(particleWeights_,0.,1);
    particleCloud.dataExchangeM().allocateArray(particleVolumes_,0.,1);
    particleCloud.get_cellIDs(cellIDs_);  // get ref to cellIDs
    //particleCloud.dataExchangeM().allocateArray(cellIDs_,0.,1);
    

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        count+=DEM_dump_Interval;// proceed loading new data

        particleCloud.averagingM().resetVectorAverage(particleCloud.averagingM().UsPrev(),particleCloud.averagingM().UsNext());
        particleCloud.voidFractionM().resetVoidFractions();
        particleCloud.averagingM().resetVectorAverage(particleCloud.forceM(0).impParticleForces(),particleCloud.forceM(0).impParticleForces(),true);
        particleCloud.averagingM().resetVectorAverage(particleCloud.forceM(0).expParticleForces(),particleCloud.forceM(0).expParticleForces(),true);
        particleCloud.averagingM().resetWeightFields();
        particleCloud.momCoupleM(0).resetMomSourceField();

        particleCloud.dataExchangeM().couple();

        particleCloud.dataExchangeM().getData("x","vector-atom",positions_,count);
        particleCloud.dataExchangeM().getData("v","vector-atom",velocities_,count);
        particleCloud.dataExchangeM().getData("radius","scalar-atom",radii_,count);

        particleCloud.locateM().findCell(NULL,positions_,cellIDs_,particleCloud.numberOfParticles());
        particleCloud.setPos(positions_);

        particleCloud.voidFractionM().setvoidFraction(NULL,voidfractions_,particleWeights_,particleVolumes_);

        voidfraction.internalField() = particleCloud.voidFractionM().voidFractionInterp();
        voidfraction.correctBoundaryConditions();

        particleCloud.averagingM().setVectorAverage
        (
            particleCloud.averagingM().UsNext(),
            velocities_,
            particleWeights_,
            particleCloud.averagingM().UsWeightField(),
            NULL
        );

        runTime.write();

        particleCloud.IOM().dumpDEMdata();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

    }

    particleCloud.dataExchangeM().destroy(positions_,3);
    particleCloud.dataExchangeM().destroy(velocities_,3);
    //particleCloud.dataExchangeM().destroy(radii_); // destroyed in cloud
    particleCloud.dataExchangeM().destroy(voidfractions_,1);
    particleCloud.dataExchangeM().destroy(particleWeights_,1);
    particleCloud.dataExchangeM().destroy(particleVolumes_,1);
    //particleCloud.dataExchangeM().destroy(cellIDs_); // destroyed in cloud

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
