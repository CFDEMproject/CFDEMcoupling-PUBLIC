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
    alphaMin_(readScalar(propsDict_.lookup("alphaMin"))),
    bb_(particleCloud_.mesh().points(),false),
    cellVol_(particleCloud_.mesh().V()[0]),
    cellLength_(pow(cellVol_,1./3.)),
    nCellXYZ_
    (
        round((bb_.max()[0]-bb_.min()[0])/cellLength_),
        round((bb_.max()[1]-bb_.min()[1])/cellLength_),
        round((bb_.max()[2]-bb_.min()[2])/cellLength_)
    )
{
    maxCellsPerParticle_=8;
    checkWeightNporosity(propsDict_);
    if(porosity()!=1) FatalError << "porosity not used in trilinearVoidFraction" << abort(FatalError);

    Warning << "trilinearVoidFraction model is ot yet complete and does not work near boundaries" << endl;
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
    vector pt(0,0,0);
    vector posShift(0,0,0);
    vector offsetCell(0,0,0);
    vector offsetOrigin(0,0,0);

    label i000(0);
    label i100(0);
    label i110(0);
    label i101(0);
    label i111(0);
    label i010(0);
    label i011(0);
    label i001(0);

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

    //bool alphaLimited(false);

    for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {
        // reset
        cellsPerParticle_[index][0]=8;
        //TODO do we need to set particleVolumes, particleV?
        // ===

        label cellI = particleCloud_.cellIDs()[index][0];

        if (cellI >= 0)  // particel centre is in domain
        {
            radius = particleCloud_.radius(index);
            volume = 4.188790205*radius*radius*radius*scaleVol;

            // store volume for each particle
            particleVolumes[index][0] = volume;
            particleV[index][0] = volume;

            // find a,b,c
            partPos=particleCloud_.position(index);
            offsetCell=partPos-particleCloud_.mesh().C()[cellI];
            a=offsetCell[0];
            b=offsetCell[1];
            c=offsetCell[2];

            // find "origin" index for mapping
            if(a>0){
                if(b>0){
                    if(c>0) //FNE
                        i000=cellI;
                    else    //BNE
                        i000=cellI-nCellXYZ_[0]*nCellXYZ_[1];
                }else{
                    if(c>0) //FSE
                        i000=cellI-nCellXYZ_[0];
                    else    //BSE
                        i000=cellI-nCellXYZ_[0]-nCellXYZ_[0]*nCellXYZ_[1];
                }
            }else{
                if(b>0){
                    if(c>0) //FNW
                        i000=cellI-1;
                    else    //BNW
                        i000=cellI-1-nCellXYZ_[0]*nCellXYZ_[1];
                }else{
                    if(c>0) //FSW
                        i000=cellI-1-nCellXYZ_[0];                     
                    else    //BSW
                        i000=cellI-1-nCellXYZ_[0]-nCellXYZ_[0]*nCellXYZ_[1];
                }
            }

            // check boundaries
            // TODO different handling for periodic and processor boundaries
            pt=particleCloud_.mesh().C()[cellI];
            posShift=vector(0,0,0);
            if(a>0){
                pt+=vector(cellLength_,0,0);
                if(pt[0]>bb_.max()[0]){
                    i000-=1;
                    posShift[0]=-a;
                }                
            }else{
                pt-=vector(cellLength_,0,0);
                if(pt[0]<bb_.min()[0])
                {
                    i000+=1;
                    posShift[0]=-a;
                }
            }
            if(b>0){
                pt+=vector(0,cellLength_,0);
                if(pt[1]>bb_.max()[1]){
                    i000-=nCellXYZ_[0];
                    posShift[1]=-b;
                }                
            }else{
                pt-=vector(0,cellLength_,0);
                if(pt[1]<bb_.min()[1]){
                    i000+=nCellXYZ_[0];                
                    posShift[1]=-b;
                }
            }
            if(c>0){
                pt+=vector(0,0,cellLength_);
                if(pt[2]>bb_.max()[2]){
                    i000-=nCellXYZ_[0]*nCellXYZ_[1];
                    posShift[2]=-c;
                }                
            }else{
                pt-=vector(0,0,cellLength_);
                if(pt[2]<bb_.min()[2]){
                    i000+=nCellXYZ_[0]*nCellXYZ_[1];
                    posShift[2]=-c;
                }
            }

            // define other 7 indices
            i100=i000+1;
            i110=i100+nCellXYZ_[0];
            i010=i000+nCellXYZ_[0];
            i001=i000+nCellXYZ_[0]*nCellXYZ_[1];
            i101=i100+nCellXYZ_[0]*nCellXYZ_[1];
            i111=i110+nCellXYZ_[0]*nCellXYZ_[1];
            i011=i010+nCellXYZ_[0]*nCellXYZ_[1];

            // find x,y,z
            // TODO changes needed here when generalized for quader cells
            offsetOrigin=particleCloud_.mesh().C()[i000]-(partPos+posShift);
            x= mag(offsetOrigin[0])/cellLength_;
            y= mag(offsetOrigin[1])/cellLength_;
            z= mag(offsetOrigin[2])/cellLength_;

            // calculate the mapping coeffs
            C000=(1-x)*(1-y)*(1-z);
            C100=x*(1-y)*(1-z);
            C110=x*y*(1-z);
            C010=(1-x)*y*(1-z);
            C001=(1-x)*(1-y)*z;
            C101=x*(1-y)*z;
            C111=x*y*z;
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
            particleCloud_.cellIDs()[index][0]=i000;
            particleCloud_.cellIDs()[index][1]=i100;
            particleCloud_.cellIDs()[index][2]=i110;
            particleCloud_.cellIDs()[index][3]=i010;
            particleCloud_.cellIDs()[index][4]=i001;
            particleCloud_.cellIDs()[index][5]=i101;
            particleCloud_.cellIDs()[index][6]=i111;
            particleCloud_.cellIDs()[index][7]=i011;

            //distribute volume
            // TODO use different cell volume when generalized for quader cells
            voidfractionNext_[i000]-=volume*C000/cellVol_;
            voidfractionNext_[i100]-=volume*C100/cellVol_;
            voidfractionNext_[i010]-=volume*C010/cellVol_;
            voidfractionNext_[i001]-=volume*C001/cellVol_;
            voidfractionNext_[i101]-=volume*C101/cellVol_;
            voidfractionNext_[i011]-=volume*C011/cellVol_;
            voidfractionNext_[i110]-=volume*C110/cellVol_;
            voidfractionNext_[i111]-=volume*C111/cellVol_;

            // debugging
            /*Pout << "cellI=" << cellI << endl;
            Pout << "a=" << a << endl;
            Pout << "b=" << b << endl;
            Pout << "c=" << c << endl;
            Pout << "x=" << x << endl;
            Pout << "y=" << y << endl;
            Pout << "z=" << z << endl;
            Pout << "i000=" << i000 << endl;
            Pout << "i100=" << i100 << endl;
            Pout << "i010=" << i010 << endl;
            Pout << "i001=" << i001 << endl;
            Pout << "i101=" << i101 << endl;
            Pout << "i011=" << i011 << endl;
            Pout << "i110=" << i110 << endl;
            Pout << "i111=" << i111 << endl;

            Pout << "C000=" << C000 << endl;
            Pout << "C100=" << C100 << endl;
            Pout << "C010=" << C010 << endl;
            Pout << "C001=" << C001 << endl;
            Pout << "C101=" << C101 << endl;
            Pout << "C011=" << C011 << endl;
            Pout << "C110=" << C110 << endl;
            Pout << "C111=" << C111 << endl;
            Pout << "sum(Cijk)=" << C000+C100+C010+C001+C101+C011+C110+C111 << endl;*/

            /*voidfractionNext_[i000]=0.999;
            voidfractionNext_[i100]=0.999;
            voidfractionNext_[i010]=0.999;
            voidfractionNext_[i001]=0.999;
            voidfractionNext_[i101]=0.999;
            voidfractionNext_[i011]=0.999;
            voidfractionNext_[i110]=0.999;
            voidfractionNext_[i111]=0.999;*/

            // limit volumefraction
            // TODO implement limiter for all 8 indices
            /*if(voidfractionNext_[cellI] < alphaMin_ )
            {
                voidfractionNext_[cellI] = alphaMin_;
                alphaLimited = true;
            }
            if(index==0 && alphaLimited) Info<<"alpha limited to" <<alphaMin_<<endl;*/

            // store voidFraction for each particle
            voidfractions[index][0] = voidfractionNext_[cellI];

            // store cellweight for each particle  - this should not live here
            particleWeights[index][0] = 1;
        }
    }
    voidfractionNext_.correctBoundaryConditions();

    // bring voidfraction from Eulerian Field to particle array
    for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {
        label cellID = particleCloud_.cellIDs()[index][0];

        if(cellID >= 0)
        {
            voidfractions[index][0] = voidfractionNext_[i000];
            voidfractions[index][1] = voidfractionNext_[i100];
            voidfractions[index][2] = voidfractionNext_[i110];
            voidfractions[index][3] = voidfractionNext_[i010];
            voidfractions[index][4] = voidfractionNext_[i001];
            voidfractions[index][5] = voidfractionNext_[i101];
            voidfractions[index][6] = voidfractionNext_[i111];
            voidfractions[index][7] = voidfractionNext_[i011];
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
