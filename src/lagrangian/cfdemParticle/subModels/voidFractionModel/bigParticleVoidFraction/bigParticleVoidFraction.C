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

#include "bigParticleVoidFraction.H"
#include "addToRunTimeSelectionTable.H"
#include "locateModel.H"
#include "dataExchangeModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(bigParticleVoidFraction, 0);

addToRunTimeSelectionTable
(
    voidFractionModel,
    bigParticleVoidFraction,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
bigParticleVoidFraction::bigParticleVoidFraction
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    voidFractionModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    alphaMin_(propsDict_.lookupOrDefault<scalar>("alphaMin",0.1)),
    alphaLimited_(0)
{
    Info << "\n\n W A R N I N G - do not use in combination with differentialRegion model! \n\n" << endl;
    Info << "\n\n W A R N I N G - this model does not yet work properly! \n\n" << endl;
    //reading maxCellsPerParticle from dictionary
    maxCellsPerParticle_=readLabel(propsDict_.lookup("maxCellsPerParticle"));
    //particleCloud_.setMaxCellsPerParticle(readLabel(propsDict_.lookup("maxCellsPerParticle"))); // alternative to line above
    if(alphaMin_ > 1 || alphaMin_ < 0.01){ FatalError<< "alphaMin shloud be > 1 and < 0.01." << abort(FatalError); }

    checkWeightNporosity(propsDict_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

bigParticleVoidFraction::~bigParticleVoidFraction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void bigParticleVoidFraction::setvoidFraction(double** const& mask,double**& voidfractions,double**& particleWeights,double**& particleVolumes,double**& particleV) const
{
    reAllocArrays();

    voidfractionNext_ == dimensionedScalar("one", voidfractionNext_.dimensions(), 1.);

    scalar radius(-1);
    scalar volume(0);
    scalar scaleVol= weight();
    scalar scaleRadius = pow(porosity(),1./3.);

    for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {
        //if(mask[index][0])
        //{
            //reset
            for(int subcell=0;subcell<cellsPerParticle_[index][0];subcell++)
            {
                particleWeights[index][subcell]=0;
                particleVolumes[index][subcell]=0;
            }
            cellsPerParticle_[index][0]=1.0;
            particleV[index][0]=0;

            //collecting data
            label particleCenterCellID=particleCloud_.cfdemCloud::cellIDs()[index][0];
            radius = particleCloud_.radius(index);
            volume = 4./3.*M_PI*radius*radius*radius*scaleVol;
            radius *= scaleRadius;
            vector positionCenter=particleCloud_.cfdemCloud::position(index);

            if (particleCenterCellID >= 0)
            {
                labelHashSet hashSett;

                //determining label and degree of coveredness of cells covered by the particle
                buildLabelHashSet(radius, positionCenter, particleCenterCellID, hashSett);
                  //Info << "completeSize=" << hashSett.size() << ", completeList =\n" << endl;
                  //for(label i=0;i<hashSett.size();i++) Info << " ," << hashSett.toc()[i] << endl;


                //generating list with cell and subcells
                scalar hashSetLength = hashSett.size();
                if (hashSetLength > maxCellsPerParticle_)
                {
                    FatalError<< "big particle algo found more cells ("<< hashSetLength
                              <<") than storage is prepared ("<<maxCellsPerParticle_<<")" << abort(FatalError);
                }
                else if (hashSetLength > 0)
                {
                    cellsPerParticle_[index][0]=hashSetLength;

                    //making sure that the cell containing the center is the first subcell
                    particleCloud_.cfdemCloud::cellIDs()[index][0]=particleCenterCellID;
                    //deleting the cell containing the center of the particle
                    hashSett.erase(particleCenterCellID);

                    //==========================//
                    //setting the voidfractions

                    // volume occupied in every covered cell
                    scalar occupiedVolume = volume/hashSetLength;

                    // correct volumefraction of centre
                    voidfractionNext_[particleCenterCellID] -=occupiedVolume/particleCloud_.mesh().V()[particleCenterCellID];

                    particleWeights[index][0] += 1.0/hashSetLength;
                    particleVolumes[index][0] += occupiedVolume;
                    particleV[index][0] += occupiedVolume;

                      //Info << "Centre:set voidfraction in cellI=" << particleCenterCellID
                      //     << ", voidfraction =" << voidfractionNext_[particleCenterCellID] << endl;

                    // correct volumefraction of sub-cells
                    for(label i=0;i<hashSetLength-1;i++)
                    {
                        label cellI=hashSett.toc()[i];
                        particleCloud_.cfdemCloud::cellIDs()[index][i+1]=cellI; //adding subcell represenation

                        voidfractionNext_[cellI] -=occupiedVolume/particleCloud_.mesh().V()[cellI];
                        particleWeights[index][i+1] += 1.0/hashSetLength;
                        particleVolumes[index][i+1] += occupiedVolume;
                        particleV[index][0] += occupiedVolume;

                          //Info << "AFTER:set voidfraction in cellI=" << cellI
                          //     << ", voidfraction =" << voidfractionNext_[cellI] << endl;

                    }

                    // debug
                    if(index==0)
                    {
                        Info << "particle 0 is represented by " << hashSetLength << "cells" << endl;
                    }
                    //==========================//
                }//end cells found on this proc
            }// end found cells
        //}// end if masked
    }// end loop all particles

    //bringing eulerian field to particle array
    for(label index=0; index< particleCloud_.numberOfParticles(); index++)
    {
        for(label subcell=0;subcell<cellsPerParticle_[index][0];subcell++)
        {
            label cellID = particleCloud_.cfdemCloud::cellIDs()[index][subcell];

            if(cellID >= 0)
            {
                 // limiting voidfraction
                 if (voidfractionNext_[cellID] < alphaMin_) voidfractionNext_[cellID]=alphaMin_;

                 // set particle based voidfraction
                 voidfractions[index][subcell] = voidfractionNext_[cellID];
                 //Info<<"setting the voidfraction, index = "<<index<<endl;
            }
            else
            {
                voidfractions[index][subcell] = -1.;
            }
        }
    }
}

void bigParticleVoidFraction::buildLabelHashSet
(
    const scalar radius,
    const vector position,
    const label cellID,
    labelHashSet& hashSett
)const
{
    hashSett.insert(cellID);
    //Info<<"cell inserted"<<cellID<<endl;
    const labelList& nc = particleCloud_.mesh().cellCells()[cellID];
    forAll(nc,i){
        label neighbor=nc[i];
        if(!hashSett.found(neighbor) && mag(position-particleCloud_.mesh().C()[neighbor])<radius){
            buildLabelHashSet(radius,position,neighbor,hashSett);
        }
    }

}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
