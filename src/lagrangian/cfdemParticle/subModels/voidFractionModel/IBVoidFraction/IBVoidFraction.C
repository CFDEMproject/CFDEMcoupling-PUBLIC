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

#include "IBVoidFraction.H"
#include "addToRunTimeSelectionTable.H"
#include "locateModel.H"
#include "dataExchangeModel.H"

#include "mpi.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(IBVoidFraction, 0);

addToRunTimeSelectionTable
(
    voidFractionModel,
    IBVoidFraction,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
IBVoidFraction::IBVoidFraction
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    voidFractionModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    alphaMin_(readScalar(propsDict_.lookup("alphaMin"))),
    alphaLimited_(0),
    scaleUpVol_(readScalar(propsDict_.lookup("scaleUpVol"))),
    checkPeriodicCells_(false)
{
    Info << "\n\n W A R N I N G - do not use in combination with differentialRegion model! \n\n" << endl;
    //Info << "\n\n W A R N I N G - this model does not yet work properly! \n\n" << endl;
    maxCellsPerParticle_=readLabel(propsDict_.lookup("maxCellsPerParticle"));
    //particleCloud_.setMaxCellsPerParticle(readLabel(propsDict_.lookup("maxCellsPerParticle"))); // alternative to line above

    if(scaleUpVol_ < 1){ FatalError<< "scaleUpVol shloud be > 1."<< abort(FatalError); }
    if(alphaMin_ > 1 || alphaMin_ < 0.01){ FatalError<< "alphaMin shloud be > 1 and < 0.01." << abort(FatalError); }
    
    if(propsDict_.found("checkPeriodicCells")) checkPeriodicCells_=true;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

IBVoidFraction::~IBVoidFraction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void IBVoidFraction::setvoidFraction(double** const& mask,double**& voidfractions,double**& particleWeights,double**& particleVolumes,double**& particleV) const
{

    int numprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    const boundBox& globalBb = particleCloud_.mesh().bounds();

    reAllocArrays();

    voidfractionNext_.internalField()=1;

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
            label particleCenterCellID=particleCloud_.cellIDs()[index][0];
            scalar radius =  particleCloud_.radius(index);
            vector positionCenter=particleCloud_.position(index);

            if (particleCenterCellID >= 0)
            {
                labelHashSet hashSett;

                //compute the voidfraction for the cell "particleCentreCellID
                vector cellCentrePosition = particleCloud_.mesh().C()[particleCenterCellID];
                scalar centreDist=mag(cellCentrePosition-positionCenter);

                vector minPeriodicParticlePos=positionCenter;
                if(checkPeriodicCells_) //consider minimal distance to all periodic images of this particle
                {
                    centreDist = minPeriodicDistance(cellCentrePosition, positionCenter, globalBb,          
                                                     minPeriodicParticlePos);
                }

                if(centreDist + 0.5*sqrt(3.0)*pow(particleCloud_.mesh().V()[particleCenterCellID],0.33333) < radius)
                {
                    voidfractionNext_[particleCenterCellID] = 0;
                }
                else
                {
                	const labelList& vertices = particleCloud_.mesh().cellPoints()[particleCenterCellID];
                	forAll(vertices, i)
                    {
                		vector vertexPosition = particleCloud_.mesh().points()[vertices[i]];
                		scalar centreVertexDist = mag(vertexPosition-positionCenter);
                        if(checkPeriodicCells_) //consider minimal distance to all periodic images of this particle
                        {                		
                		    centreVertexDist = minPeriodicDistance(vertexPosition, positionCenter, globalBb,
                		                                           minPeriodicParticlePos);
                		}
                		
                		if(centreDist<radius &&  centreVertexDist<radius)
                		{
                			voidfractionNext_[particleCenterCellID]-=0.125;
                		}
                		else if(centreDist<radius && centreVertexDist>radius)
                		{
                			//compute lambda
                			scalar a = (vertexPosition - cellCentrePosition)
                			         & (vertexPosition - cellCentrePosition);
                			scalar b = 2. * (vertexPosition - cellCentrePosition)
                			              & (cellCentrePosition-minPeriodicParticlePos);
                			scalar c = ((cellCentrePosition-minPeriodicParticlePos)
                			         &  (cellCentrePosition-minPeriodicParticlePos)
                			           )
                			         - radius*radius;

                			scalar lambda = 0.;

                			if (b*b-4*a*c>=0)  lambda =  (-b+sqrt(b*b-4*a*c))/(2*a);
                			if (lambda > 0 && lambda <=1) voidfractionNext_[particleCenterCellID]-=lambda*.125;
                			else 
                			{
                			    lambda =  (-b-sqrt(b*b-4*a*c))/(2*a);
                			    if (lambda > 0 && lambda <=1) voidfractionNext_[particleCenterCellID]-=lambda*.125;
                			}
                		}
                		else if(centreDist>radius && centreVertexDist<radius)
                		{
                		    //compute another lambda too
                            scalar a = (vertexPosition - cellCentrePosition)
                                     & (vertexPosition - cellCentrePosition);
                            scalar b = 2.* (vertexPosition - cellCentrePosition)
                                     &     (cellCentrePosition-minPeriodicParticlePos);
                            scalar c = ( (cellCentrePosition-minPeriodicParticlePos)
                                        &(cellCentrePosition-minPeriodicParticlePos)
                                       )
                                     - radius*radius;
                            scalar lambda = 0.;

                            if(b*b-4*a*c>=0)  lambda =  (-b+sqrt(b*b-4*a*c))/(2*a);
                            if(lambda > 0 && lambda <=1) voidfractionNext_[particleCenterCellID]-=(1-lambda)*0.125;
                            else 
                            {
                                lambda =  (-b-sqrt(b*b-4*a*c))/(2*a);
                                if (lambda > 0 && lambda <=1) voidfractionNext_[particleCenterCellID]-=(1-lambda)*0.125;
                            }
                		}
                	}
                } //end particle partially overlapping with cell

                //generating list with cell and subcells
                buildLabelHashSet(radius, minPeriodicParticlePos, particleCenterCellID, hashSett, true);

                //Add cells of periodic particle images on same processor
      			if(checkPeriodicCells_) 
               	{
               	    int doPeriodicImage[3];
               	    for(int iDir=0;iDir<3;iDir++)
               	    {
              	      doPeriodicImage[iDir]= 0;
                      if( (minPeriodicParticlePos[iDir]+radius)>globalBb.max()[iDir] ) 
                      {
                         doPeriodicImage[iDir] =-1;
                      }
                      if( (minPeriodicParticlePos[iDir]-radius)<globalBb.min()[iDir] )
                      {
                         doPeriodicImage[iDir] = 1;
                      }
                   	}
               	    
               	    //scan directions and map particles
               	    List<vector> particlePosList;         //List of particle center position
              	    List<label>  particleLabelList;

              	    vector       nearestPosInMesh=vector(0.0,0.0,0.0);
               	    int copyCounter=0;
                    // Note: for other than ext one could use xx.append(x)
                    // instead of setSize
                    particlePosList.setSize(particlePosList.size()+1, minPeriodicParticlePos);
               	    
               	    //x-direction
               	    if(doPeriodicImage[0]!=0) 
               	    {
               	        particlePosList.setSize(particlePosList.size()+1, particlePosList[copyCounter]
               	                              + vector(
                                                               static_cast<double>(doPeriodicImage[0])
               	                                              *(globalBb.max()[0]-globalBb.min()[0]),
               	                                              0.0,
               	                                              0.0)
               	                               );
               	        copyCounter++;
               	    }
               	    //y-direction
               	    int currCopyCounter=copyCounter;
               	    if(doPeriodicImage[1]!=0) 
               	    {
               	       for(int yDirCop=0; yDirCop<=currCopyCounter; yDirCop++)
               	       {
               	        particlePosList.setSize(particlePosList.size()+1, particlePosList[yDirCop]
               	                              + vector(
               	                                              0.0,
                                                               static_cast<double>(doPeriodicImage[1])
               	                                              *(globalBb.max()[1]-globalBb.min()[1]),
               	                                              0.0)
               	                               );
               	        copyCounter++;
               	       }
               	    }
               	    //z-direction
               	    currCopyCounter=copyCounter;
               	    if(doPeriodicImage[2]!=0) 
               	    {
               	       for(int zDirCop=0; zDirCop<=currCopyCounter; zDirCop++)
               	       {
               	        particlePosList.setSize(particlePosList.size()+1, particlePosList[zDirCop]
               	                              + vector(
               	                                              0.0,
               	                                              0.0,
                                                               static_cast<double>(doPeriodicImage[2])
               	                                              *(globalBb.max()[2]-globalBb.min()[2])
               	                                              )
               	                               );
               	        copyCounter++;
               	       }
               	    }               	    

                    //add the nearest cell labels
                    particleLabelList.setSize(particleLabelList.size()+1,particleCenterCellID);
                    for(int iPeriodicImage=1;iPeriodicImage<=copyCounter; iPeriodicImage++)
                    {
                        label copyCellID=-1;                                        
                        label partCellId = 

                        particleCloud_.mesh().findNearestCell(particlePosList[iPeriodicImage]);
                        particleLabelList.setSize(particleLabelList.size()+1,partCellId);

                        buildLabelHashSet(radius, particlePosList[iPeriodicImage], particleLabelList[iPeriodicImage], hashSett, false);
                        
                    }

               	} //end checkPeriodicCells_


                scalar hashSetLength = hashSett.size();
                if (hashSetLength > maxCellsPerParticle_)
                {
                    FatalError<< "big particle algo found more cells ("<< hashSetLength 
                              <<") than storage is prepared ("<<maxCellsPerParticle_<<")" << abort(FatalError);
                }
                else if (hashSetLength > 0)
                {
                    cellsPerParticle_[index][0]=hashSetLength;
                    hashSett.erase(particleCenterCellID);

                    for(label i=0;i<hashSetLength-1;i++)
                    {
                        label cellI=hashSett.toc()[i];
                        particleCloud_.cellIDs()[index][i+1]=cellI; //adding subcell represenation
                    }
                }//end cells found on this proc
            }// end found cells
        //}// end if masked
    }// end loop all particles

    for(label index=0; index< particleCloud_.numberOfParticles(); index++)
    {
        for(label subcell=0;subcell<cellsPerParticle_[index][0];subcell++)
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

void IBVoidFraction::buildLabelHashSet
(
    const scalar radius,
    const vector position,
    const label cellID,
    labelHashSet& hashSett, 
    bool initialInsert //initial insertion of own cell
)const
{   

    int me;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);

    if(initialInsert)  hashSett.insert(cellID);
    
    const labelList& nc = particleCloud_.mesh().cellCells()[cellID];
    forAll(nc,i){
        label neighbor=nc[i];
        vector cellCentrePosition = particleCloud_.mesh().C()[neighbor];
        scalar centreDist = mag(cellCentrePosition-position);
        
        if(!hashSett.found(neighbor) && centreDist + 0.5*sqrt(3.0)*pow(particleCloud_.mesh().V()[neighbor],0.33333) < radius){
            voidfractionNext_[neighbor] = 0;
            buildLabelHashSet(radius,position,neighbor,hashSett,true);
        }
        else if(!hashSett.found(neighbor) && centreDist < radius + sqrt(3.0)*pow(particleCloud_.mesh().V()[neighbor],0.33333)){
            scalar scale = 1;
            const labelList& vertexPoints = particleCloud_.mesh().cellPoints()[neighbor];

            forAll(vertexPoints, j){
                vector vertexPosition = particleCloud_.mesh().points()[vertexPoints[j]];
                scalar vertexDist = mag(vertexPosition - position);

                if (centreDist < radius){
                    if (vertexDist < radius) scale -= 0.125;
                    else {
                        scalar a = (vertexPosition - cellCentrePosition)&(vertexPosition - cellCentrePosition);
                        scalar b = 2.* (vertexPosition - cellCentrePosition)&(cellCentrePosition-position);
                        scalar c = ((cellCentrePosition-position)&(cellCentrePosition-position))-radius*radius;
                        scalar lambda = 0.;

                        if(b*b-4*a*c>=0)  lambda =  (-b+sqrt(b*b-4*a*c))/(2*a);
                        if (lambda > 0 && lambda <=1) scale -=lambda * 0.125;
                        else {
                            lambda =  (-b-sqrt(b*b-4*a*c))/(2*a);
                            if (lambda > 0 && lambda <=1) scale -=lambda * 0.125;
                        }
                    }
                }
                else if (vertexDist < radius){
                    scalar a = (vertexPosition - cellCentrePosition)&(vertexPosition - cellCentrePosition);
                    scalar b = 2.* (vertexPosition - cellCentrePosition)&(cellCentrePosition-position);
                    scalar c = ((cellCentrePosition-position)&(cellCentrePosition-position))-radius*radius;
                    scalar lambda = 0.;

                    if(b*b-4*a*c>=0)  lambda =  (-b+sqrt(b*b-4*a*c))/(2*a);
                    if (lambda > 0 && lambda <=1) scale -=(1-lambda) * 0.125;
                    else {
                        lambda =  (-b-sqrt(b*b-4*a*c))/(2*a);
                        if (lambda > 0 && lambda <=1) scale -=(1-lambda) * 0.125;
                    }
                }
            }

            if(voidfractionNext_[neighbor]==1) voidfractionNext_[neighbor] = scale;
            else {
            	voidfractionNext_[neighbor] -= (1-scale);
            	if(voidfractionNext_[neighbor]<0) voidfractionNext_[neighbor] = 0;
            }
            if(!(scale == 1))  buildLabelHashSet(radius,position,neighbor,hashSett, true);
        }
    }
}

//Function to determine minimal distance of point
//to one of the periodic images of a particle
double IBVoidFraction::minPeriodicDistance(vector    cellCentrePosition, 
                                           vector    positionCenter, 
                                           boundBox  globalBb,
                                           vector&   minPeriodicPos)const
{
    double centreDist=999e32;
    vector positionCenterPeriodic;
    
    for(int xDir=-1; xDir<=1; xDir++)
    {
        positionCenterPeriodic[0] =  positionCenter[0]
                                  + static_cast<double>(xDir)
                                  * (globalBb.max()[0]-globalBb.min()[0]);
        for(int yDir=-1; yDir<=1; yDir++)
        {
            positionCenterPeriodic[1] =  positionCenter[1]
                                      + static_cast<double>(yDir)
                                      * (globalBb.max()[1]-globalBb.min()[1]);                        
            for(int zDir=-1; zDir<=1; zDir++)
            {
                positionCenterPeriodic[2] =  positionCenter[2]
                                          + static_cast<double>(zDir)
                                          * (globalBb.max()[2]-globalBb.min()[2]);
                if( mag(cellCentrePosition-positionCenterPeriodic)<centreDist)
                {
                    centreDist     = mag(cellCentrePosition-positionCenterPeriodic);
                    minPeriodicPos = positionCenterPeriodic;
                }
            }
        }
    }

    return centreDist;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
