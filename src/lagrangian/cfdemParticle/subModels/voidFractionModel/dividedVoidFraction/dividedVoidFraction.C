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

#include "dividedVoidFraction.H"
#include "addToRunTimeSelectionTable.H"
#include "locateModel.H"
#include "dataExchangeModel.H"

//#include "mpi.h"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dividedVoidFraction, 0);

addToRunTimeSelectionTable
(
    voidFractionModel,
    dividedVoidFraction,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dividedVoidFraction::dividedVoidFraction
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    voidFractionModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    alphaMin_(readScalar(propsDict_.lookup("alphaMin"))),
    alphaLimited_(0),
    tooMuch_(0.0),
    scaleUpVol_(readScalar(propsDict_.lookup("scaleUpVol"))),
    interpolation_(false)
{
    maxCellsPerParticle_ = 29;

    if(scaleUpVol_ > 1.3 || scaleUpVol_ < 1){ FatalError<< "scaleUpVol shloud be > 1 and < 1.3 !!!" << abort(FatalError); }
    if(alphaMin_ > 1 || alphaMin_ < 0.01){ FatalError<< "alphaMin shloud be > 1 and < 0.01 !!!" << abort(FatalError); }
    if (propsDict_.found("interpolation")){
        interpolation_=true;
        Warning << "interpolation for dividedVoidFraction does not yet work correctly!" << endl;
        Info << "Using interpolated voidfraction field - do not use this in combination with interpolation in drag model!"<< endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dividedVoidFraction::~dividedVoidFraction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dividedVoidFraction::setvoidFraction(double** const& mask,double**& voidfractions,double**& particleWeights,double**& particleVolumes) const
{
    reAllocArrays();

    scalar pi = M_PI;
    vector position(0,0,0);
    label cellID=-1;
    scalar radius(-1);
    scalar cellVol(0);

    for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {
        //if(mask[index][0])
        //{
            // reset
            for(int subcell=0;subcell<cellsPerParticle_[index][0];subcell++){
                particleWeights[index][subcell]=0;
                particleVolumes[index][subcell]=0;
            }

            cellsPerParticle_[index][0]=1;
            position = particleCloud_.position(index);
            cellID = particleCloud_.cellIDs()[index][0];
            radius = particleCloud_.radii()[index][0];
            scalar volume =  4./3.*radius*radius*radius*3.1415;
            radius = radius*pow(scaleUpVol_,1/3);
            cellVol=0;

            //--variables for sub-search
            int nPoints = 29;
            int nNotFound=0,nUnEqual=0,nTotal=0;
            vector offset(0,0,0);
            int cellsSet = 0;

            if (cellID >= 0)  // particel centre is in domain
            {
                cellVol = particleCloud_.mesh().V()[cellID];

                //NP for 2 different radii
                for(scalar r = 0.623926*radius;r < radius;r+=0.293976*radius)
                {
                    //NP try 8 subpoint derived from spherical coordinates
	            for (scalar zeta=pi/4.;zeta<(2.*pi);zeta+=(pi/2.))
	            {
                        for (scalar theta=(pi/4.);theta<pi;theta+=(pi/2.))
	                {
	                    offset[0]=double(r)*Foam::sin(theta)*Foam::cos(zeta);
	                    offset[1]=double(r)*Foam::sin(theta)*Foam::sin(zeta);
	                    offset[2]=double(r)*Foam::cos(theta);
                            #include "setWeightedSource.H"   // set source terms at position+offset
	                }
                    }
	            //NP try 2 more subpoints for each coordinate direction (6 total)
	            for (int j=-1;j<=1;j+=2)
	            {
	    	        offset[0]=double(r)*(double(j));
	                offset[1]=double(0.);offset[2]=double(0.);
                        #include "setWeightedSource.H"   //NP set source terms at position+offset
	                offset[1]=double(r)*(double(j));
	                offset[0]=double(0.);offset[2]=double(0.);
                        #include "setWeightedSource.H"   //NP set source terms at position+offset

	                offset[2]=double(r)*(double(j));
	                offset[0]=double(0.);offset[1]=double(0.);

                        #include "setWeightedSource.H"   //NP set source terms at position+offset
	            }
                }// end loop radiivoidfractions

	        if(cellsSet>29 || cellsSet<0)
                {
	            Info << "ERROR  cellsSet =" << cellsSet << endl;
        	}

                //NP set source for particle center; source 1/nPts+weight of all subpoints that have not been found
                scalar centreWeight = 1./nPoints*(nPoints-cellsSet);

                // update voidfraction for each particle read
                scalar newAlpha = voidfractionNext_[cellID]- volume*centreWeight/cellVol;
                if(newAlpha > alphaMin_) voidfractionNext_[cellID] = newAlpha;
                else
                {
                    voidfractionNext_[cellID] = alphaMin_;
                    tooMuch_ += (alphaMin_-newAlpha) * cellVol;
                }

                // store cellweight for each particle --- this should be done for subpoints as well!!
                particleWeights[index][0] += centreWeight;

                // store particleVolume for each particle
                particleVolumes[index][0] += volume*centreWeight;

                /*//OUTPUT
                if (index==0)
                {
                    Info << "centre cellID = " << cellID << endl;
                    Info << "cellsPerParticle_=" << cellsPerParticle_[index][0] << endl;

                    for(int i=0;i<cellsPerParticle_[index][0];i++)
                    {
                       if(i==0)Info << "cellids, voidfractions, particleWeights, : \n";
                       Info << particleCloud_.cellIDs()[index][i] << " ," << endl;
                       Info << voidfractions[index][i] << " ," << endl;
                       Info << particleWeights[index][i] << " ," << endl;
                     }
                }*/

            }// end if in cell
        //}// end if in mask
        //NP reset counter of lost volume
        if(index == particleCloud_.numberOfParticles()-1) Info << "Total particle volume neglected: " << tooMuch_<< endl;
    }// end loop all particles

    // bring voidfraction from Eulerian Field to particle arra
//    interpolationCellPoint<scalar> voidfractionInterpolator_(voidfractionNext_);
    scalar voidfractionAtPos(0);
    for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {
/*        if(interpolation_)
        {
            label cellI = particleCloud_.cellIDs()[index][0];
            if(cellI >= 0)
            {
                position = particleCloud_.position(index);
                voidfractionAtPos=voidfractionInterpolator_.interpolate(position,cellI);
            }else{
                voidfractionAtPos=-1;
            }
    
            for(int subcell=0;subcell<cellsPerParticle_[index][0];subcell++)
            {
                label cellID = particleCloud_.cellIDs()[index][subcell];

                if(cellID >= 0)
                {
                    if(voidfractionAtPos > 0)
                        voidfractions[index][subcell] = voidfractionAtPos;
                    else
                        voidfractions[index][subcell] = voidfractionNext_[cellID];
                } 
                else
                {
                    voidfractions[index][subcell] = -1.;
                }
            }
        }
        else*/
        {
            for(int subcell=0;subcell<cellsPerParticle_[index][0];subcell++)
            {
                label cellID = particleCloud_.cellIDs()[index][subcell];

                if(cellID >= 0)
                {
                    voidfractions[index][subcell] = voidfractionNext_[cellID];
                } 
                else
                {
                    voidfractions[index][subcell] = -1.;
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
