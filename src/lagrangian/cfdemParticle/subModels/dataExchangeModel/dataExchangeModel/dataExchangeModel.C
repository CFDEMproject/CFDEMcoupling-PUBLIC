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
#include "dataExchangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dataExchangeModel, 0);

defineRunTimeSelectionTable(dataExchangeModel, dictionary);

// * * * * * * * * * * * * * * protected Member Functions  * * * * * * * * * * * * * //

void Foam::dataExchangeModel::setNumberOfParticles(int numberOfParticles) const
{
    particleCloud_.setNumberOfParticles(numberOfParticles);
}

// * * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * * //

//====
// double **

void Foam::dataExchangeModel::allocateArray
(
    double**& array,
    double initVal,
    int width,
    int length
) const
{
    // this model should only be used for VTK (and File exchange model)
    //if(particleCloud_.dataExchangeM().type() != "oneWayVTK")
    //    FatalError<< "dataExchangeModel::allocateArray should not be used with your dataExchangeModel" << abort(FatalError);
    
    // allocate and init double array
    array = new double*[length];
    for (int i=0; i<length; i++)
    {
        array[i] = new double [width];

        // init array
        for(int j=0; j<width;j++) array[i][j] = initVal;
    }
}

void Foam::dataExchangeModel::allocateArray
(
    double**& array,
    double initVal,
    int width,
    const char* length
) const
{
    int len=0;
    if(strcmp(length,"nparticles")==0) len = particleCloud_.numberOfParticles();
    else if (strcmp(length,"nbodies")==0) len = particleCloud_.numberOfClumps();
    else FatalError<<"call allocateArray with length, nparticles or nbodies!\n" << abort(FatalError);
    allocateArray(array,initVal,width,len);
}
void Foam::dataExchangeModel::destroy(double** array,int len) const
{
    if (array == NULL) return;

    for ( int i = 0; i < len; i++ )
        delete [] array[i];  

    delete [] array;
}

//====
// int **
void Foam::dataExchangeModel::allocateArray
(
    int**& array,
    int initVal,
    int width,
    int length
) const
{
    // this model should only be used for VTK (and File exchange model)
    //if(particleCloud_.dataExchangeM().type() != "oneWayVTK")
    //    FatalError<< "dataExchangeModel::allocateArray should not be used with your dataExchangeModel" << abort(FatalError);

    // allocate and init double array
    array = new int*[length];
    for (int i=0; i<length; i++)
    {
        array[i] = new int [width];

        // init array
        for(int j=0; j<width;j++) array[i][j] = initVal;
    }
}

void Foam::dataExchangeModel::allocateArray
(
    int**& array,
    int initVal,
    int width,
    const char* length
) const
{
    int len=0;
    if(strcmp(length,"nparticles")==0) len = particleCloud_.numberOfParticles();
    else if (strcmp(length,"nbodies")==0) len = nClumpTypes_;
    else FatalError<<"call allocateArray with length, nparticles or nbodies!\n" << abort(FatalError);
    allocateArray(array,initVal,width,len);
}
void Foam::dataExchangeModel::destroy(int** array,int len) const
{
    if (array == NULL) return;

    for ( int i = 0; i < len; i++ )
        delete [] array[i];  

    delete [] array;
}
//====

//====
// int *
void Foam::dataExchangeModel::allocateArray
(
    int*& array,
    int initVal,
    int length
) const
{
    // this model should only be used for VTK (and File exchange model)
    //if(particleCloud_.dataExchangeM().type() != "oneWayVTK")
    //    FatalError<< "dataExchangeModel::allocateArray should not be used with your dataExchangeModel" << abort(FatalError);

    // allocate and init int array
    array = new int[length];
    for (int i=0; i<length; i++)
        array[i] = initVal;
}
void Foam::dataExchangeModel::destroy(int* array) const
{
    if (array == NULL) return;
    delete [] array;
}
//====

//====
// double *
void Foam::dataExchangeModel::allocateArray
(
    double*& array,
    double initVal,
    int length
) const
{
    // this model should only be used for VTK (and File exchange model)
    //if(particleCloud_.dataExchangeM().type() != "oneWayVTK")
    //    FatalError<< "dataExchangeModel::allocateArray should not be used with your dataExchangeModel" << abort(FatalError);

    // allocate and init double array
    array = new double[length];
    for (int i=0; i<length; i++)
        array[i] = initVal;
}
void Foam::dataExchangeModel::destroy(double* array) const
{
    if (array == NULL) return;
    delete [] array;
}
//====


bool Foam::dataExchangeModel::couple() const
{
    bool coupleNow = false;
    if (doCoupleNow())
    {
        couplingStep_++;
        coupleNow = true;
    }
    return coupleNow;
}

scalar Foam::dataExchangeModel::timeStepFraction() const
{
    //return fraction between previous coupling TS and actual TS
    //scalar DEMtime = DEMts_ * couplingInterval_;
    //scalar frac = ( ( particleCloud_.mesh().time().value()-particleCloud_.mesh().time().startTime().value() ) - (couplingStep_) * DEMtime) / DEMtime; //Chr 05.03.2013
    scalar frac = ( particleCloud_.mesh().time().value()-particleCloud_.mesh().time().startTime().value() - couplingStep_ * couplingTime() ) / couplingTime();
    if (frac<1e-4) frac = 1;

    return frac;
}
int Foam::dataExchangeModel::getNumberOfParticles() const
{
    Warning << "ask for nr of clumps - which is not supported for this dataExchange model" << endl;
    return -1;
}

int Foam::dataExchangeModel::getNumberOfClumps() const
{
    Warning << "ask for nr of clumps - which is not supported for this dataExchange model" << endl;
    return -1;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dataExchangeModel::dataExchangeModel
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    dict_(dict),
    particleCloud_(sm),
    maxNumberOfParticles_(0),
    nClumpTypes_(1),
    couplingStep_(0),
    DEMts_(-1.),
    couplingInterval_(readScalar(dict_.lookup("couplingInterval")))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dataExchangeModel::~dataExchangeModel()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
