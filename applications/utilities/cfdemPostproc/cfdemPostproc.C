/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

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
#include "regionModel.H"
#include "locateModel.H"
#include "averagingModel.H"

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

    double **positions_;
    double **velocities_;
    double **radii_;
    double **voidfractions_;
    double **particleWeights_;
    double **particleVolumes_;
    double **cellIDs_;
    
    particleCloud.dataExchangeM().allocateArray(positions_,0.,3);
    particleCloud.dataExchangeM().allocateArray(velocities_,0.,3);
    particleCloud.dataExchangeM().allocateArray(radii_,0.,1);
    particleCloud.dataExchangeM().allocateArray(voidfractions_,0.,1);
    particleCloud.dataExchangeM().allocateArray(particleWeights_,0.,1);
    particleCloud.dataExchangeM().allocateArray(particleVolumes_,0.,1);
    particleCloud.dataExchangeM().allocateArray(cellIDs_,0.,1);

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        particleCloud.regionM().resetVolFields(Us);

        particleCloud.dataExchangeM().couple();

        particleCloud.dataExchangeM().getData("x","vector-atom",positions_,count);
        particleCloud.dataExchangeM().getData("v","vector-atom",velocities_,count);
        particleCloud.dataExchangeM().getData("radius","scalar-atom",radii_,count);

        particleCloud.set_radii(radii_);

        particleCloud.locateM().findCell(particleCloud.regionM().inRegion(),positions_,cellIDs_,particleCloud.numberOfParticles());

        particleCloud.set_cellIDs(cellIDs_);

        particleCloud.voidFractionM().setvoidFraction
        (
            particleCloud.regionM().inRegion(),voidfractions_,particleWeights_,particleVolumes_
        );

        voidfraction.internalField() = particleCloud.voidFractionM().voidFractionInterp();
        voidfraction.correctBoundaryConditions();

        particleCloud.averagingM().setVectorAverage
        (
            particleCloud.averagingM().UsNext(),
            velocities_,
            particleWeights_,
            particleCloud.averagingM().UsWeightField(),
            particleCloud.regionM().inRegion()
        );

        runTime.write();
        count++;  // proceed loading new data

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

    }

    delete positions_;
    delete velocities_;
    delete radii_;
    delete voidfractions_;
    delete particleWeights_;
    delete particleVolumes_;
    delete cellIDs_;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
