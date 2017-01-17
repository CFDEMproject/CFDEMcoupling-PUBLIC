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
#include "mpi.h"
#include "clockModel.H"
#include <unistd.h>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(clockModel, 0);

defineRunTimeSelectionTable(clockModel, dictionary);

// * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * //

void Foam::clockModel::start(int pos) const
{
    start(pos,"");
    return;
}

void Foam::clockModel::start(int pos,const std::string& ident) const
{
    if(particleCloud_.mesh().time().value() > startTime_)
    {
        if (pos >= n_) // alternatively one fixed size?
        {
            n_ = 2*n_;
            deltaT_.resize(n_,0);
            identifier_.resize(n_,"");
            nOfRuns_.resize(n_,0);
            level_.resize(n_,-1);
            parent_.resize(n_,-2);
        }
        identifier_[pos]=ident;
        level_[pos] = curLev_;
        curLev_ += 1;
        parent_[pos]=curParent_;
        curParent_ = pos;
        nOfRuns_[pos] += 1;
        deltaT_[pos]-=std::clock();
    }
    return;
}

void Foam::clockModel::stop() const
{
    if(particleCloud_.mesh().time().value() > startTime_)
    {
        deltaT_[curParent_]+=std::clock();
        curLev_ -= 1;
        if (curParent_ >= 0)
        {
            curParent_ = parent_[curParent_];
        }
        else
        {
            curParent_ = -1;
        }
    }
    return;
}

void Foam::clockModel::stop(const std::string& ident) const
{
    if(particleCloud_.mesh().time().value() > startTime_)
    {
        deltaT_[curParent_] += std::clock();
        if (curParent_ > 0 && identifier_[curParent_].compare(ident)!=0)
        {
            Pout<<"Warning: stop identifier did not equal start identifier! "<<ident<<" & "<<identifier_[curParent_]<<nl;
        }
        curLev_ -= 1;
        if (curParent_ >= 0)
        {
            curParent_ = parent_[curParent_];
        }
        else
        {
            curParent_ = -1;
        }
    }
    return;
}

std::string Foam::clockModel::eval() const
{
    std::ostringstream strs("Measurements in CPU-seconds:\n");
    strs << "Name\tdeltaT\tnOfRuns\tlevel\tparentNr\tparentName\n";
    strs.setf(std::ios_base::scientific);
    std::vector<int> shifts = calcShift();

    for (int i=0; i<n_; i++)
    {
        if (parent_[i] != -2)
        {
            strs << identifier_[i] << "\t";
            strs << static_cast<double>(deltaT_[i])/(CLOCKS_PER_SEC) << "\t";
            strs << nOfRuns_[i] << "\t";
            strs << level_[i] << "\t";

            if (parent_[i] >= 0)
            {
                strs << (shifts[parent_[i]]) << "\t";
                strs << identifier_[parent_[i]] << "\n";
            }
            else
            {
                strs << parent_[i] << "\t";
                strs << "none\n";
            }
        }
    }
    return strs.str();
}

void Foam::clockModel::evalFile() const
{
    std::ofstream outFile;
    std::string fileName(path_/"timeEval.txt");
    outFile.open(fileName.c_str(), ios_base::trunc);
    outFile << "Time Evaluation" << nl;
    outFile << eval();
    outFile.close();
}

void Foam::clockModel::evalPar() const
{
    int myrank, numprocs;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    std::ofstream outFile;
    std::ostringstream strs;
    strs.setf(std::ios_base::scientific);

    std::string fileName(path_/"timeEval_");
    strs << myrank << ".txt";
    fileName.append(strs.str());

    outFile.open(fileName.c_str(), ios_base::trunc);
    outFile << "Time Evaluation for Processor Nr." << myrank << nl;
    outFile << eval();
    outFile.close();

    // MPI_REDUCE SUM NODES
    MPI_Barrier(MPI_COMM_WORLD);
    strs.str("Parallel Measurements in CPU-seconds of all Processors (starting after first t.s.):\n");
    strs << "Name\tavgdeltaT\tmaxdeltaT\tnOfRuns\tlevel\tparentNr\tparentName\n";
    double buffOut = 0.;
    double buffIn = 0.;
    std::vector<int> shifts = calcShift();

    for (int i=0; i<n_; i++)
    {
        if (parent_[i] != -2)
        {
            strs << identifier_[i] << "\t";

            buffIn = static_cast<double>(deltaT_[i])/(CLOCKS_PER_SEC);

            MPI_Reduce(&buffIn, &buffOut, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);

            strs << buffOut/numprocs << "\t";

            MPI_Reduce(&buffIn, &buffOut, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);

            strs << buffOut << "\t";
            strs << nOfRuns_[i] << "\t";
            strs << level_[i] << "\t";

            if (parent_[i] >= 0)
            {
                strs << (shifts[parent_[i]]) << "\t";
                strs << identifier_[parent_[i]] << "\n";
            }
            else
            {
                strs << parent_[i] << "\t";
                strs << "none\n";
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (myrank == 0)
    {
        std::string fileName(path_/"timeEvalFull.txt");
        outFile.open(fileName.c_str(),ios_base::trunc);
        outFile << strs.str();
        outFile.close();
    }
}


void Foam::clockModel::initElems()
{
    //init elems
    for (int i = 0;i < n_; i++)
    {
        deltaT_[i] = 0;
        identifier_[i] = "";
        nOfRuns_[i] = 0;
        level_[i] = -1;
        parent_[i] = -2;
    }
}

std::vector<int> Foam::clockModel::calcShift() const
{
    std::vector<int> shifts = std::vector<int> (n_);
    shifts[0]=0;
    for (int i=1;i<n_;i++)
    {
        if (parent_[i] == -2)
        {
            shifts[i] = shifts[i-1];
        }
        else
        {
            shifts[i] = shifts[i-1]+1;
        }
    }
    return shifts;
}

void Foam::clockModel::normHist() const
{
    int myrank=-10;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    int numprocs=-10;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    double buffOut=0.;
    double buffIn=0.;

    Info << "==========================" << endl;
    Info << " PROCESSOR LOAD HISTOGRAM" << endl;

    //Global = 1
    buffIn = double(deltaT_[1]);
    MPI_Allreduce(&buffIn, &buffOut, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if(buffOut>SMALL) buffIn = buffIn*double(numprocs)/buffOut;
    plotHist(buffIn,identifier_[1],numprocs,myrank);

    //LIGGGHTS = 3
    buffIn = double(deltaT_[3]);
    MPI_Allreduce(&buffIn, &buffOut, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if(buffOut>SMALL) buffIn = buffIn*double(numprocs)/buffOut;
    plotHist(buffIn,identifier_[3],numprocs,myrank);

    //Coupling - LIGGGHTS = 2 - 3
    buffIn = double(deltaT_[2]) - buffIn;
    MPI_Allreduce(&buffIn, &buffOut, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if(buffOut>SMALL) buffIn = buffIn*double(numprocs)/buffOut;
    plotHist(buffIn,"Coupling (routines)",numprocs,myrank);

    //Flow = 26
    buffIn = double(deltaT_[26]);
    MPI_Allreduce(&buffIn, &buffOut, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if(buffOut>SMALL) buffIn = buffIn*double(numprocs)/buffOut;
    plotHist(buffIn,identifier_[26],numprocs,myrank);
    Info << "===========================" << endl;

    getRAMUsage();
    return;
}

void Foam::clockModel::plotHist(double buffIn,const std::string& identifier,int numprocs,int myrank) const
{
    double* globalTime_all = NULL;
    if (myrank == 0) globalTime_all = new double[numprocs];
    MPI_Gather(&buffIn, 1, MPI_DOUBLE, globalTime_all, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (myrank == 0)
        for (int j=0; j<numprocs; j++)
            printf("%4f  ",globalTime_all[j]);

    Info << "\t" << identifier << endl;

    delete [] globalTime_all;
}

void Foam::clockModel::Hist() const
{
    int myrank=-10;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

    //Global = 1 / Coupling = 2 / LIGGGHTS = 3 /Flow = 26

    //Global = 1
    Pout << "[" << myrank << "]: " << identifier_[1] << " " << (deltaT_[1]/CLOCKS_PER_SEC) << '\n';
    //LIGGGHTS = 3
    Pout << "[" << myrank << "]: " << identifier_[3] << " " << (deltaT_[3]/CLOCKS_PER_SEC) << '\n';
    //Coupling - LIGGGHTS = 2 - 3
    Pout << "[" << myrank << "]: " << "Coupling - LIGGGHTS" << " " << ((deltaT_[2]-deltaT_[3])/CLOCKS_PER_SEC) << '\n';
    //Flow = 26
    Pout << "[" << myrank << "]: " << identifier_[26] << " " << (deltaT_[26]/CLOCKS_PER_SEC) << '\n';

    return;
}

void Foam::clockModel::getRAMUsage() const
{
    int myrank=-10;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    int numprocs=-10;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    pid_t myPID = getpid();    //get PID of running process
    //Pout << myPID << "\n";

    std::string fileName = "/proc/"; //build path to /proc/PID/smaps and open file
    std::stringstream strs;
    strs << myPID;

    fileName.append(strs.str());
    fileName.append("/smaps");
    std::ifstream inFile;
    inFile.open(fileName.data(),ios_base::in);

    std::string line;
    int RssMem = 0;
    int SwapMem = 0;
    int temp = 0;
    strs.str("");
    if (inFile.is_open())    //search in File smaps for Rss and Swap entries
    {
        while(inFile.good())
        {
            getline(inFile,line);
            strs.str("");
            if (line.substr(0,4).compare("Rss:") == 0)
            {
                strs << line;
                strs >> line >> temp;
                RssMem = RssMem + temp;
                //Pout << temp << " ";
            }
            else if (line.substr(0,5).compare("Swap:") == 0)
            {
                strs << line;
                strs >> line >> temp;
                SwapMem = SwapMem + temp;
                //Pout << strs.str() << " ";
            }
        }
    }
    double SwapMB = double(SwapMem)/1024.0; //kB -> MB
    double RssMB = double(RssMem)/1024.0;

    inFile.close();

    // set up communication between Procs and plot Stuff
    Info << " RAM USAGE HISTOGRAM in MB" << endl;
    plotHist(RssMB,"RSS memory used",numprocs,myrank);
    if (SwapMem > 0)
    {
        plotHist(SwapMB,"WARNING: Swap",numprocs,myrank);
    }
    Info << "===========================" << endl;

    //Pout << "SWAP Memory used: " << SwapMem <<"MB\n";
    //Pout << "Rss Memory used: " << RssMem <<"MB\n";


    return;
}
// * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::clockModel::clockModel
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    dict_(dict),
    particleCloud_(sm),
    path_("clockData"),
    startTime_(sm.mesh().time().startTime().value()+2.*sm.dataExchangeM().couplingTime()), // delay start of measurement by 2*tCouple
    //startTime_(0),                                //no delay
    n_(30),
    deltaT_(n_),
    identifier_(n_),
    nOfRuns_(n_),
    level_(n_),
    curLev_(0),
    parent_(n_),
    curParent_(0)
{

    Info << "start clock measurement at t >"  << startTime_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::clockModel::~clockModel()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
