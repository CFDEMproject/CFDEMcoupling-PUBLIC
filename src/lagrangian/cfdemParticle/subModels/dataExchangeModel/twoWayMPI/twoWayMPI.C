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
#include "twoWayMPI.H"
#include "addToRunTimeSelectionTable.H"
#include "clockModel.H"
#include "pair.h"
#include "force.h"
#include "forceModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(twoWayMPI, 0);

addToRunTimeSelectionTable
(
    dataExchangeModel,
    twoWayMPI,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
twoWayMPI::twoWayMPI
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    dataExchangeModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    lmp(NULL)
{
    Info<<"Starting up LIGGGHTS for first time execution"<<endl;

    MPI_Comm_dup(MPI_COMM_WORLD, &comm_liggghts);

    // read path from dictionary
    const fileName liggghtsPath(propsDict_.lookup("liggghtsPath"));

    // open LIGGGHTS input script
    Info<<"Executing input script '"<< liggghtsPath.c_str() <<"'"<<endl;
    lmp = new LAMMPS_NS::LAMMPS(0,NULL,comm_liggghts);
    lmp->input->file(liggghtsPath.c_str());

    // get DEM time step size
    DEMts_ = lmp->update->dt;
    checkTSsize();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

twoWayMPI::~twoWayMPI()
{
    delete lmp;
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
char* twoWayMPI::wordToChar(word& inWord) const
{
    return const_cast<char*>(inWord.c_str());
}


// * * * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * * //

void twoWayMPI::getData
(
    word name,
    word type,
    double ** const& field,
    label step
) const
{
    char* charName = wordToChar(name);
    char* charType = wordToChar(type);
    data_liggghts_to_of(charName,charType, lmp, (void*&) field, (char *)"double");
}

void twoWayMPI::getData
(
    word name,
    word type,
    int ** const& field,
    label step
) const
{
    char* charName = wordToChar(name);
    char* charType = wordToChar(type);
    data_liggghts_to_of(charName,charType, lmp, (void*&) field, (char *)"int");
}

void twoWayMPI::giveData
(
    word name,
    word type,
    double ** const& field,
    const char* datatype
) const
{
    char* charName = wordToChar(name);
    char* charType = wordToChar(type);
    char* charDatatype= const_cast<char*> (datatype);
    data_of_to_liggghts(charName,charType, lmp, (void*) field,charDatatype);
}
//============
// double **
void Foam::twoWayMPI::allocateArray
(
    double**& array,
    double initVal,
    int width,
    int length
) const
{
    //if(length==-1) then LIGGGHTS uses own length data
    allocate_external_double(array, width,length,initVal,lmp);
}

void Foam::twoWayMPI::allocateArray
(
    double**& array,
    double initVal,
    int width,
    const char* length
) const
{
    //if(length==-1) then LIGGGHTS uses own length data
    char* charLength= const_cast<char*> (length);
    allocate_external_double(array, width,charLength,initVal,lmp);
}
void Foam::twoWayMPI::destroy(double** array,int len) const
{
    if (array == NULL) return;

    //for ( int i = 0; i < len; i++ ) // does not work
    for ( int i = 0; i < 1; i++ )
        free(array[i]); 

    free(array);
}
//============
// int **
void Foam::twoWayMPI::allocateArray
(
    int**& array,
    int initVal,
    int width,
    int length
) const
{
    //if(length==-1) then LIGGGHTS uses own length data
    allocate_external_int(array, width,length,initVal,lmp);
}

void Foam::twoWayMPI::allocateArray
(
    int**& array,
    int initVal,
    int width,
    const char* length
) const
{
    //if(length==-1) then LIGGGHTS uses own length data
    char* charLength= const_cast<char*> (length);
    allocate_external_int(array, width,charLength,initVal,lmp);
}
void Foam::twoWayMPI::destroy(int** array,int len) const
{
    if (array == NULL) return;

    //for ( int i = 0; i < len; i++ ) // does not work
    for ( int i = 0; i < 1; i++ )
        free(array[i]); 

    free(array);
}

bool Foam::twoWayMPI::couple(int i) const
{
    bool coupleNow = false;
    if (i==0)
    {
        couplingStep_++;
        coupleNow = true;

        // start liggghts
        {
            // run commands from liggghtsCommands dict
            Info<<"Starting up LIGGGHTS" << endl;
            particleCloud_.clockM().start(3,"LIGGGHTS");

            // check if liggghtsCommandModels with exaxt timing are being run
            bool exactTiming(false);
            int runComNr = -10;
            DynamicList<scalar> interruptTimes(0);
            DynamicList<int> DEMstepsToInterrupt(0);
            DynamicList<int> lcModel(0);

            forAll(particleCloud_.liggghtsCommandModelList(),i)
            {
                // Check if exact timing is needed
                // get time for execution
                // store time for execution in list
                if(particleCloud_.liggghtsCommand()[i]().exactTiming())
                {
                    exactTiming = true;
                    DynamicList<scalar> h = particleCloud_.liggghtsCommand()[i]().executionsWithinPeriod(TSstart(),TSend());

                    forAll(h,j)
                    {
                        // save interrupt times (is this necessary)
                        interruptTimes.append(h[j]);

                        // calc stepsToInterrupt
                        DEMstepsToInterrupt.append(DEMstepsTillT(h[j]));

                        // remember which liggghtsCommandModel to run
                        lcModel.append(i);
                    }

                    // make cumulative
                    label len = DEMstepsToInterrupt.size();
                    label ind(0);
                    forAll(DEMstepsToInterrupt,i)
                    {
                        ind = len-i-1;
                        if(ind>0)
                            DEMstepsToInterrupt[ind] -= DEMstepsToInterrupt[ind-1];
                    }                    

                    Info << "Foam::twoWayMPI::couple(i): interruptTimes=" << interruptTimes << endl;
                    Info << "Foam::twoWayMPI::couple(i): DEMstepsToInterrupt=" << DEMstepsToInterrupt << endl;
                    Info << "Foam::twoWayMPI::couple(i): lcModel=" << lcModel << endl;
                }

                if(particleCloud_.liggghtsCommand()[i]().type()=="runLiggghts")
                    runComNr=i;
            }

            // models with exact timing exists
            label commandLines(0);
            if(exactTiming)
            {
                // extension for more liggghtsCommands active the same time:
                //    sort interrupt list within this run period
                //    keep track of corresponding liggghtsCommand
                int DEMstepsRun(0);
                
                forAll(interruptTimes,j)
                {                  
                    // set run command till interrupt
                    DEMstepsRun += DEMstepsToInterrupt[j];          
                    particleCloud_.liggghtsCommand()[runComNr]().set(DEMstepsToInterrupt[j]);
                    const char* command = particleCloud_.liggghtsCommand()[runComNr]().command(0);
                    Info << "Executing run command: '"<< command <<"'"<< endl;
                    lmp->input->one(command);
    
                    // run liggghts command with exact timing
                    command = particleCloud_.liggghtsCommand()[lcModel[j]]().command(0);
                    Info << "Executing command: '"<< command <<"'"<< endl;
                    lmp->input->one(command);
                }

                // do the run
                if(particleCloud_.liggghtsCommand()[runComNr]().runCommand(couplingStep()))
                {
                    particleCloud_.liggghtsCommand()[runComNr]().set(couplingInterval() - DEMstepsRun);
                    const char* command = particleCloud_.liggghtsCommand()[runComNr]().command(0);
                    Info << "Executing run command: '"<< command <<"'"<< endl;
                    lmp->input->one(command);
                }

                // do the other non exact timing models
                forAll(particleCloud_.liggghtsCommandModelList(),i)
                {
                    if
                    (
                      ! particleCloud_.liggghtsCommand()[i]().exactTiming() &&
                        particleCloud_.liggghtsCommand()[i]().runCommand(couplingStep())
                    )
                    {
                        commandLines=particleCloud_.liggghtsCommand()[i]().commandLines();
                        for(int j=0;j<commandLines;j++)
                        {
                            const char* command = particleCloud_.liggghtsCommand()[i]().command(j);
                            Info << "Executing command: '"<< command <<"'"<< endl;
                            lmp->input->one(command);
                        }
                    }
                }
            }
            // no model with exact timing exists
            else
            {
                forAll(particleCloud_.liggghtsCommandModelList(),i)
                {
                    if(particleCloud_.liggghtsCommand()[i]().runCommand(couplingStep()))
                    {
                        commandLines=particleCloud_.liggghtsCommand()[i]().commandLines();
                        for(int j=0;j<commandLines;j++)
                        {
                            const char* command = particleCloud_.liggghtsCommand()[i]().command(j);
                            Info << "Executing command: '"<< command <<"'"<< endl;
                            lmp->input->one(command);
                        }
                    }
                }
            }

            particleCloud_.clockM().stop("LIGGGHTS");
            Info<<"LIGGGHTS finished"<<endl;
        }

        // give nr of particles to cloud
        double newNpart = liggghts_get_maxtag(lmp);
        setNumberOfParticles(newNpart);
        setNumberOfClumps(-2);

        // re-allocate arrays of cloud
        particleCloud_.clockM().start(4,"LIGGGHTS_reallocArrays");
        particleCloud_.reAllocArrays();
        particleCloud_.clockM().stop("LIGGGHTS_reallocArrays");
    }

    return coupleNow;
}

int Foam::twoWayMPI::getNumberOfParticles() const
{
    return liggghts_get_maxtag(lmp);
}

int Foam::twoWayMPI::getNumberOfClumps() const
{
    #ifdef multisphere
        return liggghts_get_maxtag_ms(lmp);
    #endif

    Warning << "liggghts_get_maxtag_ms(lmp) is not available here!" << endl;
    return -1;       
}

int Foam::twoWayMPI::getNumberOfTypes() const
{
    #ifdef multisphere
        return liggghts_get_ntypes_ms(lmp);
    #endif
    Warning << "liggghts_get_maxtag_ms(lmp) is not available here!" << endl;
    return -1;
}

double* Foam::twoWayMPI::getTypeVol() const
{
    #ifdef multisphere
        return liggghts_get_vclump_ms(lmp);
    #endif

    Warning << "liggghts_get_vclump_ms(lmp) is not available here!" << endl;
    return NULL;       
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
