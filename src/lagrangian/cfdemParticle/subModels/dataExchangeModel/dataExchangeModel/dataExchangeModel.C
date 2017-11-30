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

// * * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * * //

void Foam::dataExchangeModel::setNumberOfParticles(int numberOfParticles) const
{
    particleCloud_.setNumberOfParticles(numberOfParticles);
}

void Foam::dataExchangeModel::setNumberOfClumps(int numberOfClumps) const
{
    particleCloud_.setNumberOfClumps(numberOfClumps);
}

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
    // allocate and init double array
    destroy(array, -1);
    double *data = new double[width*length];
    std::fill_n(data, width*length, initVal);
    array = new double*[length];

    int n = 0;
    for (int i=0; i<length; i++)
    {
        array[i] = &data[n];
        n += width;
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
    if (strcmp(length,"nparticles")==0) len = particleCloud_.numberOfParticles();
    else if (strcmp(length,"nbodies")==0) len = particleCloud_.numberOfClumps();
    else FatalError<<"call allocateArray with length, nparticles or nbodies!\n" << abort(FatalError);
    allocateArray(array,initVal,width,len);
}

void Foam::dataExchangeModel::destroy(double** array,int /*len*/) const
{
    if (array == NULL) return;

    delete [] array[0];
    delete [] array;
    array = NULL;
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
    // allocate and init int array
    destroy(array, -1);
    int *data = new int[width*length];
    std::fill_n(data, width*length, initVal);
    array = new int*[length];

    int n = 0;
    for (int i=0; i<length; i++)
    {
        array[i] = &data[n];
        n += width;
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
    if (strcmp(length,"nparticles")==0) len = particleCloud_.numberOfParticles();
    else if (strcmp(length,"nbodies")==0) len = particleCloud_.numberOfClumps();
    else FatalError<<"call allocateArray with length, nparticles or nbodies!\n" << abort(FatalError);
    allocateArray(array,initVal,width,len);
}

void Foam::dataExchangeModel::destroy(int** array,int /*len*/) const
{
    if (array == NULL) return;

    delete [] array[0];
    delete [] array;
    array = NULL;
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
    dataExchangeModel::destroy(array);

    // allocate and init int array
    array = new int[length];
    std::fill_n(array, length, initVal);
}

void Foam::dataExchangeModel::destroy(int* array) const
{
    delete [] array;
    array = NULL;
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
    dataExchangeModel::destroy(array);

    // allocate and init double array
    array = new double[length];
    std::fill_n(array, length, initVal);
}

void Foam::dataExchangeModel::destroy(double* array) const
{
    delete [] array;
    array = NULL;
}
//====


bool Foam::dataExchangeModel::couple(int i) const
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
    scalar frac = ( (particleCloud_.mesh().time().timeIndex()-timeIndexOffset_)*particleCloud_.mesh().time().deltaT().value() - (couplingStep_-1) * couplingTime() ) / couplingTime();
    return frac;
}

int Foam::dataExchangeModel::getNumberOfParticles() const
{
    FatalError << "ask for nr of particles - which is not supported for this dataExchange model" << abort(FatalError);
    return -1;
}

int Foam::dataExchangeModel::getNumberOfClumps() const
{
    FatalError << "ask for nr of clumps - which is not supported for this dataExchange model" << abort(FatalError);
    return -1;
}

int Foam::dataExchangeModel::getNumberOfTypes() const
{
    FatalError << "ask for nr of types - which is not supported for this dataExchange model" << abort(FatalError);
    return -1;
}

double* Foam::dataExchangeModel::getTypeVol() const
{
    FatalError << "ask for type volume - which is not supported for this dataExchange model" << abort(FatalError);
    return NULL;
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
    couplingStep_(0),
    DEMts_(-1.),
    couplingInterval_(readScalar(dict_.lookup("couplingInterval"))),
    timeIndexOffset_(particleCloud_.mesh().time().timeIndex())
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dataExchangeModel::~dataExchangeModel()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
