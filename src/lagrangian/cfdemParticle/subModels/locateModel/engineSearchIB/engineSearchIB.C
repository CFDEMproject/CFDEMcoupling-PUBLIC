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


#include "engineSearchIB.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

#include "mpi.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(engineSearchIB, 0);

addToRunTimeSelectionTable
(
    engineSearch,
    engineSearchIB,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
engineSearchIB::engineSearchIB
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    engineSearch(dict,sm,typeName),
    zSplit_(readLabel(propsDict_.lookup("zSplit"))),
    xySplit_(readLabel(propsDict_.lookup("xySplit"))),
    coef_(2.0),
    verbose_(propsDict_.lookupOrDefault<Switch>("verbose", false)),
    numberOfSatellitePoints_((zSplit_-1)*xySplit_ + 2)
{
    bbPtr_.reset(new boundBox(particleCloud_.mesh().points(), false));
    if(verbose_)
    {
        Pout<<"MinBounds (x,y,z): "<<bbPtr_().min()<<endl;
        Pout<<"MaxBounds (x,y,z): "<<bbPtr_().max()<<endl;
    }
    for(int countPoints = 0; countPoints < numberOfSatellitePoints_; ++countPoints)
    {
        satellitePoints_.push_back(generateSatellitePoint(countPoints));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

engineSearchIB::~engineSearchIB()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


label engineSearchIB::findCell
(
    double** const& mask,
    double**& positions,
    double**& cellIDs,
    int size,
    bool checkRad
) const
{
    bool checkPeriodicCells(particleCloud_.checkPeriodicCells());
    const boundBox& globalBb = particleCloud_.mesh().bounds();

    if(particleCloud_.meshHasUpdated())
    {
        searchEngine_.correct();
        bbPtr_.reset(new boundBox(particleCloud_.mesh().points(), false));
        if(verbose_)
        {
            Pout<<"MinBounds (x,y,z): "<<bbPtr_().min()<<endl;
            Pout<<"MaxBounds (x,y,z): "<<bbPtr_().max()<<endl;
        }
    }

    vector position;
    for(int index = 0;index < size; ++index)
    {
        cellIDs[index][0]=-1;
        double radius=particleCloud_.radius(index);
        //if(mask[index][0] && radius > SMALL)
        if(radius > SMALL)
        {
            // create pos vector
            for(int i=0;i<3;i++) position[i] = positions[index][i];

            bool isInside = isInsideRectangularDomain(position, coef_*radius);

            if(!isInside && checkPeriodicCells)
            {
                vector positionCenterPeriodic;
                for(int xDir=-1*particleCloud_.wall_periodicityCheckRange(0); 
                        xDir<=particleCloud_.wall_periodicityCheckRange(0); 
                        xDir++)
                {
                    positionCenterPeriodic[0] =  position[0]
                                              + static_cast<double>(xDir)
                                              * (globalBb.max()[0]-globalBb.min()[0]);
                    for(int yDir=-1*particleCloud_.wall_periodicityCheckRange(1); 
                            yDir<=particleCloud_.wall_periodicityCheckRange(1); 
                            yDir++)
                    {
                        positionCenterPeriodic[1] =  position[1]
                                                  + static_cast<double>(yDir)
                                                  * (globalBb.max()[1]-globalBb.min()[1]);
                        for(int zDir=-1*particleCloud_.wall_periodicityCheckRange(2); 
                                zDir<=particleCloud_.wall_periodicityCheckRange(2); 
                                zDir++)
                        {
                            positionCenterPeriodic[2] =  position[2]
                                                      + static_cast<double>(zDir)
                                                      * (globalBb.max()[2]-globalBb.min()[2]);
                            isInside = isInsideRectangularDomain(positionCenterPeriodic, coef_*radius);
                            if(isInside) break;
                        }
                        if(isInside) break;
                    }
                    if(isInside) break;
                }
            }


            if(isInside)
            {
              // find cell
              label oldID = cellIDs[index][0];
              cellIDs[index][0] = findSingleCell(position,oldID);

              if(cellIDs[index][0] < 0)
              {
                  label altStartPos = -1;

                  for(unsigned int countPoints = 0; countPoints < satellitePoints_.size(); ++countPoints)
                  {
                    vector pos = getSatellitePoint(index, countPoints);
                    isInside = isInsideRectangularDomain(pos, SMALL);

                    if(isInside)
                        altStartPos = findSingleCell(pos,oldID);

                    //check for periodic domains, only do if check range is larger than 0
                    if(checkPeriodicCells)
                    {
                        for(int iDir=0;iDir<3;iDir++)
                        {
                            if( pos[iDir] > globalBb.max()[iDir] && particleCloud_.wall_periodicityCheckRange(iDir)>0 ) 
                                pos[iDir]-=globalBb.max()[iDir]-globalBb.min()[iDir];
                            else if( pos[iDir] < globalBb.min()[iDir] && particleCloud_.wall_periodicityCheckRange(iDir)>0 )
                                pos[iDir]+=globalBb.max()[iDir]-globalBb.min()[iDir];
                        }
                        isInside = isInsideRectangularDomain(pos, SMALL);

                        if(isInside)
                          altStartPos=findSingleCell(pos,oldID); //particleCloud_.mesh().findCell(pos);//
                    }

                    if(altStartPos >= 0) // found position, we're done
                    {
                        cellIDs[index][0] = altStartPos;
                        break;
                    }
                }
              }
            }
        }
    }
    return 1;
}

bool engineSearchIB::isInsideRectangularDomain(vector centre, scalar skin) const
{
    vector offset(skin, skin, skin);
    boundBox bb(bbPtr_().min()-offset, bbPtr_().max() + offset);
    return bb.contains(centre);
}

vector engineSearchIB::generateSatellitePoint(int countPoints) const
{
    scalar theta, phi;
    const scalar thetaSize = 180./zSplit_, phiSize = 360./xySplit_;
    const scalar deg2rad = M_PI/180.;
    vector pos(0.0, 0.0, 0.0);
    // 1 point at bottom, 1 point at top
    if(countPoints == 0)
    {
        pos[2] = 1.0;
    } else if(countPoints == 1)
    {
        pos[2] = -1.0;
    } else {
        scalar thetaLevel = (countPoints - 2) / xySplit_;
        theta = deg2rad * thetaSize * (thetaLevel+1);
        phi = deg2rad * phiSize * (countPoints - 2 - thetaLevel*xySplit_);
        pos[0] = sin(theta) * cos(phi);
        pos[1] = sin(theta) * sin(phi);
        pos[2] = cos(theta);
    }
    return pos;
}

vector engineSearchIB::getSatellitePoint(int index, int countPoints) const
{
    double radius=particleCloud_.radius(index);
    vector position = particleCloud_.position(index);
    vector pos = radius*satellitePoints_[countPoints] + position;
    return pos;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
