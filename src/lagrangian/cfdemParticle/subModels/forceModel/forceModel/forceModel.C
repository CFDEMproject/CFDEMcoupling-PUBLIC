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
    modelType_(sm.modelType()),
    probeIt_(sm.probeM().active()),
    requiresEx_(false),
    requiresShape_(false),
    requiresQuaternion_(false),
    requiresSuperquadric_(false),
    pullPushRotation_(false),
    implicitDrag_(false),
    implicitAnisotropicDrag_(false),
    implicitRotation_(false),
    forceSubModels_(0),
    forceSubModel_(new autoPtr<forceSubModel>[nrForceSubModels()]),
    voidfractionInterpolator_(NULL),
    UInterpolator_(NULL),
    vorticityInterpolator_(NULL),
    gradPInterpolator_(NULL),
    gradUInterpolator_(NULL),
    gradVoidfractionInterpolator_(NULL),
    Up1Interpolator_(NULL),
    Up2Interpolator_(NULL),
    dSauterInterpolator_(NULL),
    phiP1Interpolator_(NULL),
    phiP2Interpolator_(NULL),
    alphaInterpolator_(NULL),
    gradAlphaInterpolator_(NULL),
    TInterpolator_(NULL),
    UsInterpolator_(NULL),
    fluidScalarFieldInterpolator_(NULL),
    gradPsolidInterpolator_(NULL),
    shearRateInterpolator_(NULL),
    DDtUInterpolator_(NULL),
    divTauInterpolator_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

forceModel::~forceModel()
{}

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //

void Foam::forceModel::applyDebugSettings(bool debug) const
{
    if(!debug)
    {
        impParticleForces_.writeOpt() = IOobject::NO_WRITE;
        expParticleForces_.writeOpt() = IOobject::NO_WRITE;
    }
}

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

//Function for to add turbulence due to multiphase interaction
void forceModel::multiphaseTurbulence(volScalarField& field, bool) const
{
    // just return zero
    field     *= 0.0;
}

//Function for simple explicit treatment of coupling terms, only for temperature
void forceModel::manipulateScalarField(volScalarField& field) const
{
    // just return zero
    field     *= 0.0;
}

//Function for explicit or implicit treatment of coupling terms, for heat and species balance equations
void forceModel::manipulateScalarField(volScalarField& field, volScalarField& fieldImpl, int speciesID) const
{
    // just return zero
    field     *= 0.0;
    fieldImpl *= 0.0;
}

//Function for return particle based data to DEM
void forceModel::commToDEM() const
{
    // do noting.
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
    {
        forceSubModels_ = wordList(dict.lookup("forceSubModels"));
    }
    else if (dict.found("forceSubModel"))
    {
        FatalError << "Did you mean the forceSubModels keyword? " << abort(FatalError);
    }
    else
    {
        forceSubModels_.setSize(1, "ImEx");
    }

    delete[] forceSubModel_;
    forceSubModel_ = new autoPtr<forceSubModel>[nrForceSubModels()];
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
