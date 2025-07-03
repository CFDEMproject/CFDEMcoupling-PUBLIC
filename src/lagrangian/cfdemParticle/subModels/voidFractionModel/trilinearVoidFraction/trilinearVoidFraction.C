/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
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

#include "trilinearVoidFraction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(trilinearVoidFraction, 0);

addToRunTimeSelectionTable
(
    voidFractionModel,
    trilinearVoidFraction,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
trilinearVoidFraction::trilinearVoidFraction
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    voidFractionModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    alphaMin_(propsDict_.lookupOrDefault<scalar>("alphaMin",0.1)),
    bb_(particleCloud_.mesh().points(),false)

{
    maxCellsPerParticle_=8;
    if(alphaMin_ > 1 || alphaMin_ < 0.01){ FatalError<< "alphaMin should be < 1 and > 0.01 !!!" << abort(FatalError); }
    checkWeightNporosity(propsDict_);
    if(porosity()!=1) FatalError << "porosity not used in trilinearVoidFraction" << abort(FatalError);

    //Warning << "trilinearVoidFraction model is ot yet complete and does not work near boundaries" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

trilinearVoidFraction::~trilinearVoidFraction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void trilinearVoidFraction::setvoidFraction(double** const& mask,double**& voidfractions,double**& particleWeights,double**& particleVolumes,double**& particleV) const
{
    reAllocArrays();

    scalar radius(-1);
    scalar volume(0);
    scalar scaleVol=weight();

    vector partPos(0,0,0);
    vector cellCenter(0,0,0);
    vector pt(0,0,0);
    vector posShift(0,0,0);
    vector offsetCell(0,0,0);
    vector offsetOrigin(0,0,0);

    label i000(-1);
    label i100(-1);
    label i110(-1);
    label i101(-1);
    label i111(-1);
    label i010(-1);
    label i011(-1);
    label i001(-1);



    //bool alphaLimited(false);

    for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {

        i000 = i100 = i010 = i001 = i101 = i011 = i110 = i111 = -1;
        
        scalar C000(0);
        scalar C100(0);
        scalar C110(0);
        scalar C101(0);
        scalar C111(0);
        scalar C010(0);
        scalar C011(0);
        scalar C001(0);

        scalar x(0);
        scalar y(0);
        scalar z(0);

        scalar a(0);
        scalar b(0);
        scalar c(0);

        label cellI = particleCloud_.cfdemCloud::cellIDs()[index][0];

        if (cellI >= 0)
        {
            // get all neighbors of the cell (n1):
            labelList n1 = particleCloud_.mesh().cellCells()[cellI];
            // get all neighbors from the neighbors: 
            labelList allNeighs;
            forAll(n1,neigh)
            {
                allNeighs.append(particleCloud_.mesh().cellCells()[n1[neigh]]);
            }

            // Find all neighs that occur more than once (part. for 3D case this has to be done twice): 
            labelList finalNeighs;
            forAll(allNeighs,neigh)
            {
                label elem = allNeighs[neigh];
                if (findIndices(allNeighs,elem).size()>1 && findIndices(finalNeighs,elem).size()==0 && elem != cellI)
                    finalNeighs.append(elem);         
            }
            labelList allNeighs3D;
            forAll(finalNeighs,neigh)
            {
                allNeighs3D.append(particleCloud_.mesh().cellCells()[finalNeighs[neigh]]);
            }
            labelList finalNeighs3D;
            forAll(allNeighs3D,neigh)
            {
                label elem = allNeighs3D[neigh];
                if (findIndices(allNeighs3D,elem).size()>1 && findIndices(finalNeighs,elem).size()==0  && findIndices(finalNeighs3D,elem).size()==0 && findIndices(n1,elem).size() == 0 && elem != cellI)
                    finalNeighs3D.append(elem);  
            }
            // collect all cells that could be in the marked area (completeList)
            labelList completeList;
            completeList.append(cellI);
            completeList.append(n1);
            completeList.append(finalNeighs);
            completeList.append(finalNeighs3D);

            // now find out which ones are actually required:
            partPos = particleCloud_.cfdemCloud::position(index);
            cellCenter = particleCloud_.mesh().C()[cellI];
            offsetCell = partPos-cellCenter;

            a = offsetCell[0];
            b = offsetCell[1];
            c = offsetCell[2];

            label keeper = cellI; 

            // search for the starting cell
            if(a+SMALL>=0 && b+SMALL>=0 && c+SMALL>=0)
                i000 = cellI;
            else
            {
                // candidates
                labelList candidates;
                for (int i=0; i<completeList.size(); i++)
                {
                    label elem = completeList[i];
                    vector ccelem = particleCloud_.mesh().C()[elem];
                    if(ccelem[0]<=cellCenter[0]+SMALL && ccelem[1]<=cellCenter[1]+SMALL && ccelem[2]<=cellCenter[2]+SMALL)
                        candidates.append(elem);
                }

                labelList possibleX;
                if (a < 0)
                {
                    for (int i=0; i<candidates.size(); i++)
                    {
                        label elem = candidates[i];
                        if (particleCloud_.mesh().C()[elem][0] < cellCenter[0])
                            possibleX.append(elem);
                    }
                }
                else
                {
                    for (int i = 0; i < candidates.size(); i++)
                    {
                        label elem = candidates[i];
                        if (particleCloud_.mesh().C()[elem][0]+SMALL >= cellCenter[0])
                            possibleX.append(elem);
                    }
                }
                if(possibleX.size()==1)
                    keeper = possibleX[0];

                labelList possibleY;
                if (b < 0)
                {
                    for (int i = 0; i < possibleX.size(); i++)
                    {
                        label elem = possibleX[i];
                        if (particleCloud_.mesh().C()[elem][1] < cellCenter[1])
                            possibleY.append(elem);
                    }
                }
                else
                {
                    for (int i = 0; i < possibleX.size(); i++)
                    {
                        label elem = possibleX[i];
                        if (particleCloud_.mesh().C()[elem][1]+SMALL >= cellCenter[1])
                            possibleY.append(elem);
                    }
                }
                if(possibleY.size()==1)
                    keeper=possibleY[0];

                if (c < 0)
                {
                    for (int i = 0; i < possibleY.size(); i++)
                    {
                        label elem = possibleY[i];
                        if (particleCloud_.mesh().C()[elem][2] < cellCenter[2])
                            i000 = elem;
                    }
                }
                else
                {
                    for (int i = 0; i < possibleY.size(); i++)
                    {
                        label elem = possibleY[i];
                        if (particleCloud_.mesh().C()[elem][2]+SMALL >= cellCenter[2])
                            i000 = elem;
                    }
                }
            }

            // if for some reason the starting cell is not in the investigated domain set to current cell
            //Info << "i000: " << i000 << endl;
            if (i000 == -1)
                i000 = keeper; //cellI;

            // i100, i010 and i001 are amongst the neigbors of i000 if they can be found:
            vector centerI000 = particleCloud_.mesh().C()[i000];
            labelList neighsI000 = particleCloud_.mesh().cellCells()[i000];
            for (int i=0; i<neighsI000.size();i++)
            {
                label elem = neighsI000[i];
                vector ccelem = particleCloud_.mesh().C()[elem];
                if (ccelem[0] >= centerI000[0]+SMALL )
                    i100 = elem;
                if (ccelem[1] >= centerI000[1]+SMALL)
                    i010 = elem;
                if (ccelem[2] >= centerI000[2]+SMALL)
                    i001 = elem;
            }

            // i110, i101, i011 are shared neighbors of two of the three previous
            // get i110: 
            if (i100 != -1 && i010 != -1)
            {
                labelList allNeighs = particleCloud_.mesh().cellCells()[i100];
                allNeighs.append(particleCloud_.mesh().cellCells()[i010]);
                for(int i=0;i<allNeighs.size();i++)
                {
                    label elem = allNeighs[i];
                    if (findIndices(allNeighs,elem).size()>1 && elem != i000)
                    i110 = elem;         
                }
            }

            // get i101
            if (i100 != -1 && i001 != -1)
            {
                labelList allNeighs = particleCloud_.mesh().cellCells()[i100];
                allNeighs.append(particleCloud_.mesh().cellCells()[i001]);
                for(int i=0;i<allNeighs.size();i++)
                {
                    label elem = allNeighs[i];
                    if (findIndices(allNeighs,elem).size()>1 && elem != i000)
                    {
                        i101 = elem;         
                    }
                }
            }

            // get i011
            if (i010 != -1 && i001 != -1)
            {
                labelList allNeighs = particleCloud_.mesh().cellCells()[i010];
                allNeighs.append(particleCloud_.mesh().cellCells()[i001]);
                for(int i=0;i<allNeighs.size();i++)
                {
                    label elem = allNeighs[i];
                    if (findIndices(allNeighs,elem).size()>1 && elem != i000)
                        i011 = elem;         
                }
            }

            // i111 is a shared neighbor of i110, i101 and i011 
            // get i111
            //if (i010 != -1 && i001 != -1)
            if (i110 != -1 && i101 != -1 && i011 != -1)
            {
                labelList allNeighs = particleCloud_.mesh().cellCells()[i101];
                allNeighs.append(particleCloud_.mesh().cellCells()[i011]);
                allNeighs.append(particleCloud_.mesh().cellCells()[i110]);
                for(int i=0;i<allNeighs.size();i++)
                {
                    label elem = allNeighs[i];
                    if (findIndices(allNeighs,elem).size()>2 && elem != i000)
                        i111 = elem;         
                }
            }

            offsetOrigin = particleCloud_.mesh().C()[i000] - partPos;

            if (i100 != -1)
                x = mag(offsetOrigin[0])/mag(particleCloud_.mesh().C()[i000][0]-particleCloud_.mesh().C()[i100][0]);
            if (i010 != -1)
                y = mag(offsetOrigin[1])/mag(particleCloud_.mesh().C()[i000][1]-particleCloud_.mesh().C()[i010][1]);
            if (i001 != -1)
                z = mag(offsetOrigin[2])/mag(particleCloud_.mesh().C()[i000][2]-particleCloud_.mesh().C()[i001][2]);

            // calculate the mapping coeffs - distribute to neighboring cells 
            C000=(1-x)*(1-y)*(1-z);
            C100=x*(1-y)*(1-z);
            if (i110 == -1)
            {
                C100 += 0.5*x*y*(1-z);
                C010 += 0.5*x*y*(1-z);
            }
            else
                C110=x*y*(1-z);
            C010=(1-x)*y*(1-z);
            C001=(1-x)*(1-y)*z;
            if (i101 == -1)
            {
                C100 += 0.5*x*(1-y)*z;
                C001 += 0.5*x*(1-y)*z;
            }
            else
                C101=x*(1-y)*z;
            if (i111 == -1)
            {
                C101 +=0.33*x*y*z;
                C011 +=0.33*x*y*z;
                C110 +=0.34*x*y*z;
            }
            else
                C111=x*y*z;
            if (i011 == -1)
            {
                C010 += 0.5*(1-x)*y*z;
                C001 += 0.5*(1-x)*y*z;
            }
            else 
                C011=(1-x)*y*z;

            // set weights
            particleWeights[index][0]=C000;
            particleWeights[index][1]=C100;
            particleWeights[index][2]=C110;
            particleWeights[index][3]=C010;
            particleWeights[index][4]=C001;
            particleWeights[index][5]=C101;
            particleWeights[index][6]=C111;
            particleWeights[index][7]=C011;

            // set cellIDs
            particleCloud_.cfdemCloud::cellIDs()[index][0]=i000;
            particleCloud_.cfdemCloud::cellIDs()[index][1]=i100;
            particleCloud_.cfdemCloud::cellIDs()[index][2]=i110;
            particleCloud_.cfdemCloud::cellIDs()[index][3]=i010;
            particleCloud_.cfdemCloud::cellIDs()[index][4]=i001;
            particleCloud_.cfdemCloud::cellIDs()[index][5]=i101;
            particleCloud_.cfdemCloud::cellIDs()[index][6]=i111;
            particleCloud_.cfdemCloud::cellIDs()[index][7]=i011;

            radius = particleCloud_.radius(index);
            volume = 4./3.*M_PI*radius*radius*radius*scaleVol;

            // change voidfraction according to covered volume
            if (i000 > -1 && alphaMin_ < voidfractionNext_[i000]-volume*C000/particleCloud_.mesh().V()[i000]) 
                voidfractionNext_[i000]-=volume*C000/particleCloud_.mesh().V()[i000];
            if (i100 > -1 && alphaMin_ < voidfractionNext_[i100]-volume*C100/particleCloud_.mesh().V()[i100])
                voidfractionNext_[i100]-=volume*C100/particleCloud_.mesh().V()[i100];
            if (i010 > -1 && alphaMin_ < voidfractionNext_[i010]-volume*C010/particleCloud_.mesh().V()[i010])
                voidfractionNext_[i010]-=volume*C010/particleCloud_.mesh().V()[i010];
            if (i001 > -1 && alphaMin_ < voidfractionNext_[i001]-volume*C001/particleCloud_.mesh().V()[i001])
                voidfractionNext_[i001]-=volume*C001/particleCloud_.mesh().V()[i001];
            if (i101 > -1 && alphaMin_ < voidfractionNext_[i101]-volume*C101/particleCloud_.mesh().V()[i101])
                voidfractionNext_[i101]-=volume*C101/particleCloud_.mesh().V()[i101];
            if (i011 > -1 && alphaMin_ < voidfractionNext_[i011]-volume*C011/particleCloud_.mesh().V()[i011])
                voidfractionNext_[i011]-=volume*C011/particleCloud_.mesh().V()[i011];
            if (i110 > -1 && alphaMin_ < voidfractionNext_[i110]-volume*C110/particleCloud_.mesh().V()[i110])
                voidfractionNext_[i110]-=volume*C110/particleCloud_.mesh().V()[i110];
            if (i111 > -1 && alphaMin_ < voidfractionNext_[i111]-volume*C111/particleCloud_.mesh().V()[i111])
                voidfractionNext_[i111]-=volume*C111/particleCloud_.mesh().V()[i111];

            // debugging
            /*Info << "partPos: " << partPos << endl;
            Info << "cellI=" << cellI << endl;
            Info << "a=" << a << endl;
            Info << "b=" << b << endl;
            Info << "c=" << c << endl;
            Info << "x=" << x << endl;
            Info << "y=" << y << endl;
            Info << "z=" << z << endl;
            Info << "i000=" << i000 << endl;
            Info << "i100=" << i100 << endl;
            Info << "i010=" << i010 << endl;
            Info << "i001=" << i001 << endl;
            Info << "i101=" << i101 << endl;
            Info << "i011=" << i011 << endl;
            Info << "i110=" << i110 << endl;
            Info << "i111=" << i111 << endl;

            Info << "C000=" << C000 << endl;
            Info << "C100=" << C100 << endl;
            Info << "C010=" << C010 << endl;
            Info << "C001=" << C001 << endl;
            Info << "C101=" << C101 << endl;
            Info << "C011=" << C011 << endl;
            Info << "C110=" << C110 << endl;
            Info << "C111=" << C111 << endl;
            Info << "sum(Cijk)=" << C000+C100+C010+C001+C101+C011+C110+C111 << endl;*/

            // store voidfraction for each particle
            voidfractions[index][0] = voidfractionNext_[cellI];

            // store cellweight for each particle  - this should not live here
            particleWeights[index][0] = 1;
        }
    }

    voidfractionNext_.correctBoundaryConditions();

    // bring voidfraction from Eulerian Field to particle array
    for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {
        label cellID = particleCloud_.cfdemCloud::cellIDs()[index][0];

        if(cellID >= 0)
        {
            if (i000 != -1) voidfractions[index][0] = voidfractionNext_[i000];
            if (i100 != -1) voidfractions[index][1] = voidfractionNext_[i100];
            if (i110 != -1) voidfractions[index][2] = voidfractionNext_[i110];
            if (i010 != -1) voidfractions[index][3] = voidfractionNext_[i010];
            if (i001 != -1) voidfractions[index][4] = voidfractionNext_[i001];
            if (i101 != -1) voidfractions[index][5] = voidfractionNext_[i101];
            if (i111 != -1) voidfractions[index][6] = voidfractionNext_[i111];
            if (i011 != -1) voidfractions[index][7] = voidfractionNext_[i011];
        }
        else
        {
            for(int i=0;i<8;i++)
                voidfractions[index][i] = -1.;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
