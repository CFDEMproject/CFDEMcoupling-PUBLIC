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

#include "centreVoidFraction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(centreVoidFraction, 0);

addToRunTimeSelectionTable
(
    voidFractionModel,
    centreVoidFraction,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
centreVoidFraction::centreVoidFraction
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    voidFractionModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    alphaMin_(readScalar(propsDict_.lookup("alphaMin"))),
    alphaLimited_(0)
{
    checkWeightNporosity(propsDict_);
    if(porosity()!=1) FatalError << "porosity not used in centreVoidFraction" << abort(FatalError);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

centreVoidFraction::~centreVoidFraction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void centreVoidFraction::setvoidFraction(double** const& mask,double**& voidfractions,double**& particleWeights,double**& particleVolumes,double**& particleV) const
{
    reAllocArrays();

    scalar radius(-1);
    scalar volume(0);
    scalar cellVol(0);
    scalar scaleVol= weight();

    for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {
        //if(mask[index][0])
        //{
            // reset
            particleWeights[index][0]=0;
            cellsPerParticle_[index][0]=1;

            label cellI = particleCloud_.cellIDs()[index][0];

            if (cellI >= 0)  // particel centre is in domain
            {
                cellVol = voidfractionNext_.mesh().V()[cellI];
                radius = particleCloud_.radius(index);
                volume = 4.188790205*radius*radius*radius*scaleVol;

                // store volume for each particle
                particleVolumes[index][0] = volume;
                particleV[index][0] = volume;

                voidfractionNext_[cellI] -= volume/cellVol;

                if(voidfractionNext_[cellI] < alphaMin_ )
                {
                    voidfractionNext_[cellI] = alphaMin_;
                    alphaLimited_ = 1;
                }

                if(index==0 && alphaLimited_) Info<<"alpha limited to" <<alphaMin_<<endl;

                // store voidFraction for each particle
                voidfractions[index][0] = voidfractionNext_[cellI];

                // store cellweight for each particle  - this should not live here
                particleWeights[index][0] = 1;

                /*//OUTPUT
                if (index==0)
                {
                    Info << "centre cellI = " << cellI << endl;
                    Info << "cellsPerParticle_=" << cellsPerParticle_[index][0] << endl;

                    for(int i=0;i<cellsPerParticle_[index][0];i++)
                    {
                       if(i==0)Info << "cellids, voidfractions, particleWeights, : \n";
                       Info << particleCloud_.cellIDs()[index][i] << " ," << endl;
                       Info << voidfractions[index][i] << " ," << endl;
                       Info << particleWeights[index][i] << " ," << endl;
                     }
                }*/
            }
        //}
    }
    voidfractionNext_.correctBoundaryConditions();

    // bring voidfraction from Eulerian Field to particle array
    for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {
        label cellID = particleCloud_.cellIDs()[index][0];

        if(cellID >= 0)
        {
            voidfractions[index][0] = voidfractionNext_[cellID];
        }
        else
        {
            voidfractions[index][0] = -1.;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
