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

#include "cfdemCloudIB.H"
#include "voidFractionModel.H"
#include "forceModel.H"
#include "regionModel.H"
#include "locateModel.H"
#include "dataExchangeModel.H"
#include "IOModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
cfdemCloudIB::cfdemCloudIB
(
    const fvMesh& mesh
)
:
    cfdemCloud(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cfdemCloudIB::~cfdemCloudIB()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::cfdemCloudIB::evolve()
{
    numberOfParticlesChanged_ = false;
    arraysReallocated_=false;
    bool doCouple=false;

    if (dataExchangeM().couple())
    {
        Info << "\n timeStepFraction() = " << dataExchangeM().timeStepFraction() << endl;
        doCouple=true;

        if(verbose_) Info << "- defineRegion()" << endl;
        regionM().defineRegion();
        if(verbose_) Info << "defineRegion done." << endl;

        if(verbose_) Info << "- getDEMdata()" << endl;
        getDEMdata();
        Info << "nr particles = " << numberOfParticles() << endl;

        // search cellID of particles
        if(verbose_) Info << "- findCell()" << endl;
        locateM().findCell(regionM().inRegion(),positions_,cellIDs_,numberOfParticles());
        if(verbose_) Info << "findCell done." << endl;

        // set void fraction field
        if(verbose_) Info << "- setvoidFraction()" << endl;
        voidFractionM().setvoidFraction(regionM().inRegion(),voidfractions_,particleWeights_,particleVolumes_);
        if(verbose_) Info << "setvoidFraction done." << endl;

        // set particles forces
        if(verbose_) Info << "- setForce(forces_)" << endl;
        for(int index = 0;index <  numberOfParticles_; ++index){
            for(int i=0;i<3;i++){
                impForces_[index][i] = 0;
                expForces_[index][i] = 0;
            }
        }
        for (int i=0;i<nrForceModels();i++) forceM(i).setForce(regionM().inRegion(),impForces_,expForces_,expForces_);
        if(verbose_) Info << "setForce done." << endl;

        // write DEM data
        if(verbose_) Info << " -giveDEMdata()" << endl;
        giveDEMdata();
    }
    Info << "evolve done." << endl;

    //if(verbose_)    #include "debugInfo.H";

    // do particle IO
    IOM().dumpDEMdata();

    return doCouple;
}

void Foam::cfdemCloudIB::calcVelocityCorrection
(
    volScalarField& p,
    volVectorField& U,
    volScalarField& phiIB
)
{
    label cellI=0;
    vector uParticle(0,0,0);
    for(int index=0; index< numberOfParticles(); index++)
    {
        if(regionM().inRegion()[index][0])
        {
            for(int subCell=0;subCell<voidFractionM().cellsPerParticle()[index][0];subCell++)
            {
                //Info << "subCell=" << subCell << endl;
                cellI = cellIDs()[index][subCell];

                if (cellI >= 0)
                {
                    for(int i=0;i<3;i++) uParticle[i] = velocities()[index][i];
                    U[cellI] = uParticle;
                }
            }
        }
    }

    // make field divergence free
    solve(fvm::laplacian(phiIB) == fvc::div(U));

    U=U-fvc::grad(phiIB);
    U.correctBoundaryConditions();

    // correct the pressure as well
    p=p+phiIB/U.mesh().time().deltaT();  // do we have to  account for rho here?
    p.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
