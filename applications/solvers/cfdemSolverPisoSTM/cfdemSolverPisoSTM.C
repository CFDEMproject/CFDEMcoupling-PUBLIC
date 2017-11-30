/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright (C) 1991-2009 OpenCFD Ltd.
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
    cfdemSolverPisoSTM

Description
    Transient solver for incompressible flow.
    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.
    The code is an evolution of the solver pisoFoam in OpenFOAM(R) 1.6,
    where additional functionality for CFD-DEM coupling is added.
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
#if defined(versionv1606plus) || defined(version40)
    #include "fvOptions.H"
#else
    #include "fvIOoptionList.H"
#endif
#include "fixedFluxPressureFvPatchScalarField.H"
#ifdef MS
    #include "cfdemCloudMS.H"
#else
    #include "cfdemCloud.H"
#endif
#if defined(anisotropicRotation)
    #include "cfdemCloudRotation.H"
#endif
#include "implicitCouple.H"
#include "clockModel.H"
#include "smoothingModel.H"
#include "forceModel.H"

#include "scalarTransportModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #if defined(version30)
        pisoControl piso(mesh);
        #include "createTimeControls.H"
    #endif
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    // create cfdemCloud
    #include "readGravitationalAcceleration.H"
    #include "checkImCoupleM.H"
    #if defined(anisotropicRotation)
        cfdemCloudRotation particleCloud(mesh);
    #else
        #ifdef MS
            cfdemCloudMS particleCloud(mesh);
        #else
            cfdemCloud particleCloud(mesh);
        #endif
    #endif
    #include "checkModelType.H"

    // create a scalarTransportModel
    autoPtr<scalarTransportModel> stm
    (
        scalarTransportModel::New(particleCloud.couplingProperties(),particleCloud)
    );

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #if defined(version30)
            #include "readTimeControls.H"
            #include "CourantNo.H"
            #include "setDeltaT.H"
        #else
            #include "readPISOControls.H"
            #include "CourantNo.H"
        #endif

        // do particle stuff
        particleCloud.clockM().start(1,"Global");
        particleCloud.clockM().start(2,"Coupling");
        bool hasEvolved = particleCloud.evolve(voidfraction,Us,U);

        if(hasEvolved)
        {
            particleCloud.smoothingM().smoothenAbsolutField(particleCloud.forceM(0).impParticleForces());
        }
    
        Ksl = particleCloud.momCoupleM(particleCloud.registryM().getProperty("implicitCouple_index")).impMomSource();
        Ksl.correctBoundaryConditions();

        surfaceScalarField voidfractionf = fvc::interpolate(voidfraction);
        phi = voidfractionf*phiByVoidfraction;

        //Force Checks
        #include "forceCheckIm.H"

        #include "solverDebugInfo.H"
        particleCloud.clockM().stop("Coupling");

        particleCloud.clockM().start(26,"Flow");

        //Scalar transport if desired. Use "none" (noTransport) if no scalar transport is desired
        stm().update();

        if(particleCloud.solveFlow())
        {
            // Pressure-velocity PISO corrector
            {
                // Momentum predictor
                fvVectorMatrix UEqn
                (
                    fvm::ddt(voidfraction,U) - fvm::Sp(fvc::ddt(voidfraction),U)
                  + fvm::div(phi,U) - fvm::Sp(fvc::div(phi),U)
//                + turbulence->divDevReff(U)
                  + particleCloud.divVoidfractionTau(U, voidfraction)
                  ==
                  - fvm::Sp(Ksl/rho,U)
                  + fvOptions(U)
                );

                UEqn.relax();
                fvOptions.constrain(UEqn);

                #if defined(version30)
                    if (piso.momentumPredictor())
                #else
                    if (momentumPredictor)
                #endif
                {
                    if (modelType=="B" || modelType=="Bfull")
                        solve(UEqn == - fvc::grad(p) + Ksl/rho*Us);
                    else
                        solve(UEqn == - voidfraction*fvc::grad(p) + Ksl/rho*Us);

                    fvOptions.correct(U);
                }

                // --- PISO loop
                #if defined(version30)
                    while (piso.correct())
                #else
                    for (int corr=0; corr<nCorr; corr++)
                #endif
                {
                    volScalarField rUA = 1.0/UEqn.A();

                    surfaceScalarField rUAf("(1|A(U))", fvc::interpolate(rUA));
                    volScalarField rUAvoidfraction("(voidfraction2|A(U))",rUA*voidfraction);
                    surfaceScalarField rUAfvoidfraction("(voidfraction2|A(U)F)", fvc::interpolate(rUAvoidfraction));

                    U = rUA*UEqn.H();

                    #ifdef version23
                        phi = ( fvc::interpolate(U) & mesh.Sf() )
                            + rUAfvoidfraction*fvc::ddtCorr(U, phiByVoidfraction);
                    #else
                        phi = ( fvc::interpolate(U) & mesh.Sf() )
                            + fvc::ddtPhiCorr(rUAvoidfraction, U, phiByVoidfraction);
                    #endif
                    surfaceScalarField phiS(fvc::interpolate(Us) & mesh.Sf());
                    phi += rUAf*(fvc::interpolate(Ksl/rho) * phiS);

                    if (modelType=="A")
                        rUAvoidfraction = volScalarField("(voidfraction2|A(U))",rUA*voidfraction*voidfraction);

                    // Update the fixedFluxPressure BCs to ensure flux consistency
                    #include "fixedFluxPressureHandling.H"                    

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
                            fvm::laplacian(rUAvoidfraction, p) == fvc::div(voidfractionf*phi) + particleCloud.ddtVoidfraction()
                        );
                        pEqn.setReference(pRefCell, pRefValue);

                        #if defined(version30)
                            pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));
                            if (piso.finalNonOrthogonalIter())
                            {
                                phiByVoidfraction = phi - pEqn.flux()/voidfractionf;
                            }
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
                            {
                                phiByVoidfraction = phi - pEqn.flux()/voidfractionf;
                            }
                        #endif

                    } // end non-orthogonal corrector loop

                    phi = voidfractionf*phiByVoidfraction;
                    #include "continuityErrorPhiPU.H"

                    if (modelType=="B" || modelType=="Bfull")
                        U -= rUA*fvc::grad(p) - Ksl/rho*Us*rUA;
                    else
                        U -= voidfraction*rUA*fvc::grad(p) - Ksl/rho*Us*rUA;

                    U.correctBoundaryConditions();
                    fvOptions.correct(U);

                } // end piso loop
            }

            laminarTransport.correct();
            turbulence->correct();
        }// end solveFlow
        else
        {
            Info << "skipping flow solution." << endl;
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        particleCloud.clockM().stop("Flow");
        particleCloud.clockM().stop("Global");
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
