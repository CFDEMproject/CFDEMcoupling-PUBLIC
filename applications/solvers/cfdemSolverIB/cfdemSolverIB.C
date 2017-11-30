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
    cfdemSolverIB

Description
    Transient solver for incompressible flow.
    The code is an evolution of the solver pisoFoam in OpenFOAM(R) 1.6, 
    where additional functionality for CFD-DEM coupling using immersed body
    (fictitious domain) method is added.
Contributions
    Alice Hager
\*---------------------------------------------------------------------------*/


#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "OFversion.H"

#if defined(version30)
    #include "turbulentTransportModel.H"
    #include "pisoControl.H"
#else
    #include "turbulenceModel.H"
#endif
#include "cfdemCloudIB.H"
#if defined(superquadrics_flag)
#include "cfdemCloudIBSuperquadric.H"
#endif
#include "implicitCouple.H"

#include "averagingModel.H"
#include "voidFractionModel.H"

#include "dynamicFvMesh.H"

#include "cellSet.H"

#if defined(version22)
    #include "meshToMeshNew.H"
    #include "fvIOoptionList.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"

    #include "createDynamicFvMesh.H"

    #if defined(version30)
        pisoControl piso(mesh);
        #include "createTimeControls.H"
    #endif

    #include "createFields.H"

    #include "initContinuityErrs.H"

    #if defined(version22)
        #include "createFvOptions.H"
    #endif

    // create cfdemCloud
    #include "readGravitationalAcceleration.H"
    #if defined(superquadrics_flag)
        cfdemCloudIBSuperquadric particleCloud(mesh);
    #else
        cfdemCloudIB particleCloud(mesh);
    #endif

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        //=== dyM ===================
        interFace = mag(mesh.lookupObject<volScalarField>("voidfractionNext"));
        particleCloud.setMeshHasUpdatedFlag(mesh.update()); //dyM

        #if defined(version30)
            #include "readTimeControls.H"
            #include "CourantNo.H"
            #include "setDeltaT.H"
        #else
            #include "readPISOControls.H"
            #include "CourantNo.H"
        #endif

        // do particle stuff
        Info << "- evolve()" << endl;
        particleCloud.evolve(voidfraction, interFace);

        // Pressure-velocity PISO corrector
        if(particleCloud.solveFlow())
        {
            // Momentum predictor

            fvVectorMatrix UEqn
            (
                fvm::ddt(U) //fvm::ddt(voidfraction,U)
              + fvm::div(phi, U)
              + turbulence->divDevReff(U)
                #if defined(version22)
                ==
                fvOptions(U)
                #endif
            );

            UEqn.relax();

            #if defined(version22)
            fvOptions.constrain(UEqn);
            #endif

            #if defined(version30)
                if (piso.momentumPredictor())
            #else
                if (momentumPredictor)
            #endif
            {
                solve(UEqn == -fvc::grad(p));
            }

            // --- PISO loop
            #if defined(version30)
                while (piso.correct())
            #else
                for (int corr=0; corr<nCorr; corr++)
            #endif
            {
                volScalarField rUA = 1.0/UEqn.A();
                surfaceScalarField rUAf(fvc::interpolate(rUA));

                U = rUA*UEqn.H();
                #ifdef version23
                phi = (fvc::interpolate(U) & mesh.Sf())
                    + rUAf*fvc::ddtCorr(U, phi);
                #else
                phi = (fvc::interpolate(U) & mesh.Sf())
                    + fvc::ddtPhiCorr(rUA, U, phi);
                #endif
                adjustPhi(phi, U, p);

                #if defined(version22)
                fvOptions.relativeFlux(phi);
                #endif

                // Non-orthogonal pressure corrector loop
                #if defined(version30)
                    while (piso.correctNonOrthogonal())
                #else
                    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
                #endif
                {
                    // Pressure corrector

                    fvScalarMatrix pEqn
                    (
                        fvm::laplacian(rUA, p) == fvc::div(phi) + particleCloud.ddtVoidfraction()
                    );

                    pEqn.setReference(pRefCell, pRefValue);

                    #if defined(version30)
                        pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));
                        if (piso.finalNonOrthogonalIter())
                            phi -= pEqn.flux();
                    #else
                        if( corr == nCorr-1 && nonOrth == nNonOrthCorr )
                            #if defined(versionExt32)
                                pEqn.solve(mesh.solutionDict().solver("pFinal"));
                            #else
                                pEqn.solve(mesh.solver("pFinal"));
                            #endif
                        else
                            pEqn.solve();

                        if (nonOrth == nNonOrthCorr)
                            phi -= pEqn.flux();
                    #endif
                }

                #include "continuityErrs.H"

                U -= rUA*fvc::grad(p);
                U.correctBoundaryConditions();
            } 
        } //end solveFlow

        laminarTransport.correct();
        turbulence->correct();

        Info << "particleCloud.calcVelocityCorrection() " << endl;
        volScalarField voidfractionNext=mesh.lookupObject<volScalarField>("voidfractionNext");
        particleCloud.calcVelocityCorrection(p,U,phiIB,voidfractionNext);

        #if defined(version22)
        fvOptions.correct(U);
        #endif

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
