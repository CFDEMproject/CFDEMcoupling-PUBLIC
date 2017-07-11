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
#include "voidFractionModel.H"
#include "dataExchangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(voidFractionModel, 0);

defineRunTimeSelectionTable(voidFractionModel, dictionary);

// * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
voidFractionModel::voidFractionModel
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    dict_(dict),
    particleCloud_(sm),
    voidfractionPrev_
    (   IOobject
        (
            "voidfractionPrev",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,//MUST_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh().lookupObject<volScalarField> ("voidfraction")
        /*sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0), 1)*/
    ),
    voidfractionNext_
    (   IOobject
        (
            "voidfractionNext",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,//MUST_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh().lookupObject<volScalarField> ("voidfraction")
        /*sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0), 1)*/
    ),
    cellsPerParticle_(NULL),
    maxCellsPerParticle_(1),
    weight_(1.),
    porosity_(1.),
    requiresSuperquadric_(false)
{
    particleCloud_.dataExchangeM().allocateArray(cellsPerParticle_,1,maxCellsPerParticle_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

voidFractionModel::~voidFractionModel()
{
    particleCloud_.dataExchangeM().destroy(cellsPerParticle_,1);
}

// * * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * * //

void Foam::voidFractionModel::applyDebugSettings(bool debug) const
{
    if(!debug)
    {
        voidfractionPrev_.writeOpt() = IOobject::NO_WRITE;
        voidfractionNext_.writeOpt() = IOobject::NO_WRITE;
    }
}

tmp<volScalarField> Foam::voidFractionModel::voidFractionInterp() const
{
    scalar tsf = particleCloud_.dataExchangeM().timeStepFraction();
    /*if(1-tsf < 1e-4 && particleCloud_.dataExchangeM().couplingStep() > 1)   // if no subTS &&  !firstTS
    {
        //Info << "using voidfractionPrev" << endl;
        return tmp<volScalarField>
        (
            new volScalarField("alpha_voidFractionModel", voidfractionPrev_)
        );
    }
    else*/                                                                    // if subTS || firstTS
    {
        //Info << "using voidfraction blend, tsf=" << tsf << endl;
        return tmp<volScalarField>
        (
            new volScalarField("alpha_voidFractionModel", (1 - tsf) * voidfractionPrev_ + tsf * voidfractionNext_)
        );
    }
}

void Foam::voidFractionModel::resetVoidFractions() const
{
    voidfractionPrev_ == voidfractionNext_;
    voidfractionNext_ == dimensionedScalar("one", voidfractionNext_.dimensions(), 1.);
}

/*void Foam::voidFractionModel::undoVoidFractions(double**const& mask) const
{
    voidfractionPrev_.internalField() = voidfractionNext_.internalField();

    for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {
        if(mask[index][0])
        {
            // undo voidfraction cause by particle
            label cellI = particleCloud_.cellIDs()[index][0];
            scalar cellVolume=voidfractionNext_.mesh().V()[cellI];
            voidfractionNext_[cellI] += particleCloud_.particleVolumes()[index][0]/cellVolume;
        }
    }
}*/

double** const& Foam::voidFractionModel::cellsPerParticle() const
{
    return cellsPerParticle_;
}

int Foam::voidFractionModel::maxCellsPerParticle() const
{
    return maxCellsPerParticle_;
}

void Foam::voidFractionModel::reAllocArrays() const
{
    if(particleCloud_.numberOfParticlesChanged())
    {
        // get arrays of new length
        particleCloud_.dataExchangeM().allocateArray(cellsPerParticle_,1,1);
    }
}

void Foam::voidFractionModel::reAllocArrays(int nP) const
{
    if(particleCloud_.numberOfParticlesChanged())
    {
        // get arrays of new length
        particleCloud_.dataExchangeM().allocateArray(cellsPerParticle_,1,1,nP);
    }
}

bool Foam::voidFractionModel::requiresSuperquadric() const
{
  return requiresSuperquadric_;
}

double Foam::voidFractionModel::pointInParticle(int index, vector positionCenter, vector point, double scale) const
{
    scalar radius =  particleCloud_.radius(index);
    //Pout << "radius=" << radius << endl;
    if(radius>SMALL)
    {
        scalar pointDistSq = magSqr(point - positionCenter);
        return pointDistSq / (scale*scale*radius*radius) - 1.0;
    }
    else
    {
        return 0.;
    }
}

double Foam::voidFractionModel::pointInParticle(int index, vector positionCenter, vector point) const {
  return pointInParticle(index, positionCenter, point, 1.0);
}

//Function to determine minimal distance of point
//to one of the periodic images of a particle
double Foam::voidFractionModel::minPeriodicDistance(int index,
                                                    vector    cellCentrePosition,
                                                    vector    positionCenter,
                                                    boundBox  globalBb,
                                                    vector&   minPeriodicPos,
                                                    vector    dirCheckRange
                                                   )const
{
    double f=999e32;
    vector positionCenterPeriodic;

    for(  int xDir=-static_cast<int>(dirCheckRange[0]); 
              xDir<=static_cast<int>(dirCheckRange[0]); 
              xDir++)
    {
        positionCenterPeriodic[0] =  positionCenter[0]
                                  + static_cast<double>(xDir)
                                  * (globalBb.max()[0]-globalBb.min()[0]);
        for(int yDir=-static_cast<int>(dirCheckRange[1]);
                yDir<=static_cast<int>(dirCheckRange[1]); 
                yDir++)
        {
            positionCenterPeriodic[1] =  positionCenter[1]
                                      + static_cast<double>(yDir)
                                      * (globalBb.max()[1]-globalBb.min()[1]);
            for(int zDir=-static_cast<int>(dirCheckRange[2]); 
                    zDir<=static_cast<int>(dirCheckRange[2]); 
                    zDir++)
            {
                positionCenterPeriodic[2] =  positionCenter[2]
                                          + static_cast<double>(zDir)
                                          * (globalBb.max()[2]-globalBb.min()[2]);
                //if( mag(cellCentrePosition-positionCenterPeriodic)<f)
                if(pointInParticle(index, positionCenterPeriodic, cellCentrePosition) < f)
                {
                    f = pointInParticle(index, positionCenterPeriodic, cellCentrePosition);
                    minPeriodicPos = positionCenterPeriodic;
                }
            }
        }
    }

    return f;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
