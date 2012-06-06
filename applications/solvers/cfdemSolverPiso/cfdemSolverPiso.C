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
    cfdemSolverPiso

Description
    Transient solver for incompressible flow.
    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.
    The code is an evolution of the solver pisoFoam in OpenFOAM 1.6,
    where additional functionality for CFD-DEM coupling is added.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"

#include "cfdemCloud.H"
#include "implicitCouple.H"
#include "clockModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    #include "initContinuityErrs.H"

    // create cfdemCloud
    #include "readGravitationalAcceleration.H"
    cfdemCloud particleCloud(mesh);

    #include "checkModelType.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    particleCloud.clockM().start(1,"Global");

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "readPISOControls.H"
        #include "CourantNo.H"

        // do particle stuff
        particleCloud.clockM().start(2,"Coupling");
        Info << "- evolve()" << endl;
        particleCloud.evolve(voidfraction,Us,U);

        Info << "update Ksl.internalField()" << endl;
        Ksl.internalField() = particleCloud.momCoupleM(0).impMomSource();
        Ksl.correctBoundaryConditions();

        #include "solverDebugInfo.H"
        particleCloud.clockM().stop("Coupling");

        particleCloud.clockM().start(10,"Flow");
        // Pressure-velocity PISO corrector
        {
            // Momentum predictor
            fvVectorMatrix UEqn
            (
                fvm::ddt(voidfraction,U)
              + fvm::div(phi, U)
//              + turbulence->divDevReff(U)
              + particleCloud.divVoidfractionTau(U, voidfraction)
             ==
              - fvm::Sp(Ksl/rho,U)
            );

            UEqn.relax();

            if (momentumPredictor)
            {
                //solve UEqn
                if (modelType=="B")
                    solve(UEqn == - fvc::grad(p) + Ksl/rho*Us);
                else
                    solve(UEqn == - voidfraction*fvc::grad(p) + Ksl/rho*Us);
            }

            // --- PISO loop

            //for (int corr=0; corr<nCorr; corr++)
            int nCorrSoph = nCorr + 5 * pow((1-particleCloud.dataExchangeM().timeStepFraction()),1);

            Info << "nCorrSoph = " << nCorrSoph << endl;
            Info << "particleCloud.dataExchangeM().timeStepFraction() = " << particleCloud.dataExchangeM().timeStepFraction() << endl;

            for (int corr=0; corr<nCorrSoph; corr++)
            {
                volScalarField rUA = 1.0/UEqn.A();
                surfaceScalarField rUAf("(1|A(U))", fvc::interpolate(rUA));

                U = rUA*UEqn.H();

                phi = fvc::interpolate(U*voidfraction) & mesh.Sf();
                                      //+ fvc::ddtPhiCorr(rUA, U, phi)
                surfaceScalarField phiS(fvc::interpolate(Us*voidfraction) & mesh.Sf());
                surfaceScalarField phiGes = phi + rUAf*(fvc::interpolate(Ksl/rho) * phiS);

                volScalarField rUAvoidfraction("(voidfraction2|A(U))",rUA*voidfraction);
                if (modelType=="A")
                    rUAvoidfraction = volScalarField("(voidfraction2|A(U))",rUA*voidfraction*voidfraction);

                // Non-orthogonal pressure corrector loop
                for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
                {
                    // Pressure corrector
                    fvScalarMatrix pEqn
                    (
                        fvm::laplacian(rUAvoidfraction, p) == fvc::div(phiGes) + fvc::ddt(voidfraction)
                    );
                    pEqn.setReference(pRefCell, pRefValue);

                    if
                    (
                        corr == nCorr-1
                     && nonOrth == nNonOrthCorr
                    )
                    {
                        pEqn.solve(mesh.solver("pFinal"));
                    }
                    else
                    {
                        pEqn.solve();
                    }

                    if (nonOrth == nNonOrthCorr)
                    {
                        phiGes -= pEqn.flux();
                    }

                } // end non-orthogonal corrector loop

                #include "continuityErrs.H"

                if (modelType=="B")
                    U -= rUA*fvc::grad(p) - Ksl/rho*Us*rUA;
                else
                    U -= voidfraction*rUA*fvc::grad(p) - Ksl/rho*Us*rUA;

                U.correctBoundaryConditions();

            } // end piso loop
        }

        turbulence->correct();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        particleCloud.clockM().stop("Flow");
    }

    Info<< "End\n" << endl;

    particleCloud.clockM().stop("Global");

    return 0;
}


// ************************************************************************* //
