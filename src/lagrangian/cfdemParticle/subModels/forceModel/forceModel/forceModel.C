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
#include "forceModel.H"
#include "mathExtra.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(forceModel, 0);

defineRunTimeSelectionTable(forceModel, dictionary);

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
forceModel::forceModel
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    dict_(dict),
    particleCloud_(sm),
    //treatExplicit_(false),
    //treatDEM_(false),
    //implDEM_(false),
    impParticleForces_
    (   IOobject
        (
            "impParticleForces",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        sm.mesh(),
        dimensionedVector("zero", dimensionSet(1,1,-2,0,0), vector(0,0,0)) // N
    ),
    expParticleForces_
    (   IOobject
        (
            "expParticleForces",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        sm.mesh(),
        dimensionedVector("zero", dimensionSet(1,1,-2,0,0), vector(0,0,0)) // N
    ),
    coupleForce_(true),
    modelType_(sm.modelType()),
    probeIt_(sm.probeM().active()),
    requiresEx_(false),
    forceSubModels_(wordList(0)),
    forceSubModel_(new autoPtr<forceSubModel>[nrForceSubModels()])
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

forceModel::~forceModel()
{}

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //
/*tmp<volScalarField> forceModel::provideScalarField()
{
Info << "now providing a scalar field" << endl;
    tmp<volScalarField> tsource
    (
        new volScalarField
        (
            IOobject
            (
                "xxx",
                particleCloud_.mesh().time().timeName(),
                particleCloud_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            particleCloud_.mesh(),
            dimensionedScalar
            (
                "zero",
                dimensionSet(0,0,0,0,0),
                0.0
            )
        )
    );

    manipulateScalarField(tsource());
    return tsource;
};*/

void forceModel::manipulateScalarField(volScalarField& field) const
{
    Info << "no scalar manipulation done" << endl;
    // do nothing
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void forceModel::repartitionImExForces() const
{
  if(particleCloud_.imExSplitFactor()<1.0)
  {
    Info << "Will re-partition split of implicit and explicit forces: alpha = " 
         << particleCloud_.imExSplitFactor() << endl;
    // Update implicit particle
    expParticleForces_ += (1.0-particleCloud_.imExSplitFactor())*impParticleForces_;
    impParticleForces_ *= particleCloud_.imExSplitFactor();
  }
}

void forceModel::treatVoidCells() const
{
  //force coupling force in cells where there are no particles to be explicit force
  if(particleCloud_.treatVoidCellsAsExplicitForce())
  {
        int counter(0);
        volVectorField& Us = particleCloud_.averagingM().UsNext();
        forAll(Us,cellI)
        {
            if ( mag(Us[cellI]) == 0.0)  // cell is void of particles
            {
                expParticleForces_[cellI] += impParticleForces_[cellI];
                impParticleForces_[cellI] *= 0.0;
                counter +=1;
            }
        }
        Info << "Re-partitioned "<<  counter << " cells void of particles" << endl;
  }
}

void forceModel::setForceSubModels(dictionary& dict)
{
    if (dict.found("forceSubModels"))
        forceSubModels_ = wordList(dict.lookup("forceSubModels"));
    else{
        forceSubModels_ = wordList(1);
        forceSubModels_[0] = "ImEx";
    }

    delete[] forceSubModel_;
    forceSubModel_ = new autoPtr<forceSubModel>[nrForceSubModels()];
    Info << "nrForceSubModels()=" << nrForceSubModels() << endl;
    for (int i=0;i<nrForceSubModels();i++)
    {
        forceSubModel_[i] = forceSubModel::New
        (
            dict,
            particleCloud_,
            *this,
            forceSubModels_[i]
        );
    }
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
