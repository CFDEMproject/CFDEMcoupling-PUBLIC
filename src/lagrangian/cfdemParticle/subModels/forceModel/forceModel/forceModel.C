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
    particleBased_(dict.lookupOrDefault<Switch>("particleBased", false)),
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
    divTauInterpolator_(NULL),
    RhoInterpolator_(NULL),
    kInterpolator_(NULL),
    epsilonInterpolator_(NULL)
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

//Function for to add turbulence due to multiphase interaction
void forceModel::multiphaseTurbulence(volScalarField& field, bool) const
{
    // just return zero
    field     *= 0.0;
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
        if(forceSubModels_.size()==0) // empty list found
        {
            Info << " Found empty list of forceSubModels - setting to default forceSubModel ImEx" << endl;
            forceSubModels_.setSize(1, "ImEx");
        }
        else
        {
            // check if an ImEx model is there
            bool foundImEx=false;
            forAll(forceSubModels_,i)
            {
                if(forceSubModels_[i]=="ImEx" || forceSubModels_[i]=="ImExCorr" || forceSubModels_[i]=="ImExDipole" || forceSubModels_[i]=="ImExFibre") foundImEx=true;
            }

            if(!foundImEx) // add ImEx if none found
            {
                Info << " Adding the default forceSubModel ImEx on top of your selection (as first model)." << endl;
                wordList h(0);
                h.setSize(1, "ImEx");
                h.append(forceSubModels_);
                forceSubModels_=h;
            }

            if(forceSubModels_.size()>1)
            {
               Warning << " You are using more than one forceSubModel: "
                        << forceSubModels_ <<"Make sure all operations can be superposed.\n" << endl;
            }
        }
    }
    else if (dict.found("forceSubModel"))
    {
        FatalError << "Did you mean the forceSubModels keyword? " << abort(FatalError);
    }
    else // use ImEx as default forceSubModel if nothing is specifed
    {
        Info << " No forceSubModel specified - setting to default forceSubModel ImEx." << endl;
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
        particleCloud_.registryM().addProperty(type()+"_"+forceSubModel_[i]().type()+"_index",i);
    }
}

void forceModel::readDHcorr(dictionary& dict)
{
    scalarList DHc(dict.lookup("DHc"));
    particleCloud_.setDHc(DHc);
}

void forceModel::readArea(dictionary& dict)
{
    scalarList area(dict.lookup("area"));
    particleCloud_.setArea(area);

    if(forceSubM(0).verbose() && area.size() > 0 &&  forceSubM(0).getCG()>1)
    {
        Warning << "\n\n==============================================================================\n"
         << "       ! ! !  W A R N I N G  ! ! !\n"
         << "  You specified area and coarse graining is active - take care you specify\n"
         << "  area for the coarse grained particles!\n" 
         << "       ! ! !  W A R N I N G  ! ! !\n"
         << "==============================================================================\n\n" << endl;
    }
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
