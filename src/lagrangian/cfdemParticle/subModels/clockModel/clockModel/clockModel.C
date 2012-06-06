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
    and OpenFOAM. Note: this code is not part of OpenFOAM (see DISCLAIMER).
\*---------------------------------------------------------------------------*/

#include "error.H"
#include "mpi.h"
#include "clockModel.H"

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

void Foam::clockModel::start(int pos,std::string ident) const
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
	return;
}

void Foam::clockModel::stop() const
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
	return;
}

void Foam::clockModel::stop(std::string ident) const
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
	return;
}

std::string Foam::clockModel::eval() const
{
	std::string msg = "Measurements in CPU-seconds:";
	msg.append("\n");
	msg.append("Name \t deltaT \t nOfRuns \t level \t parentNr \t parentName \n");
	std::ostringstream strs;
	strs.setf(std::ios_base::scientific);
	for (int i=0;i<n_;i++)
	{
		if (parent_[i] != -2)
		{
			msg.append(identifier_[i]);
			msg.append("\t");
			strs << (double(deltaT_[i])/(CLOCKS_PER_SEC));
			msg.append(strs.str());
			msg.append("\t");
			strs.str("");

			strs << nOfRuns_[i];
			msg.append(strs.str());
			msg.append("\t");
			strs.str("");

			strs << level_[i];
			msg.append(strs.str());
			msg.append("\t");
			strs.str("");

			strs << parent_[i];
			msg.append(strs.str());
			msg.append("\t");
			strs.str("");

			if (parent_[i] >= 0)
			{
				strs << identifier_[parent_[i]];
			}
			else
			{
				strs << "none";
			}

			msg.append(strs.str());
			msg.append("\n");
			strs.str("");
		}
	}
	return msg;
}

void Foam::clockModel::evalFile() const
{
	std::ofstream outFile;
	std::string fileName(path_/"timeEval.txt");
	outFile.open(fileName.data(),ios_base::app);
	outFile << "Time Evaluation"<<nl;
	outFile << eval();
	outFile.close();

}

void Foam::clockModel::evalPar() const
{
	int myrank=-10;
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
	int numprocs=-10;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	std::ofstream outFile;
	std::ostringstream strs;
	strs.setf(std::ios_base::scientific);
	std::string fileName(path_/"timeEval_");
	strs << myrank;
	fileName.append(strs.str());
	fileName.append(".txt");
	outFile.open(fileName.data(),ios_base::app);
	outFile << "Time Evaluation for Processor Nr." << myrank <<nl;
	outFile << eval();
	outFile.close();

	// MPI_REDUCE SUM NODES
	MPI_Barrier(MPI_COMM_WORLD);
	strs.str("");
	std::string msg = "Parallel Measurements in CPU-seconds of all Processors:";
	msg.append("\n");
	msg.append("Name \t avgdeltaT \t maxdeltaT \t nOfRuns \t level \t parentNr \t parentName \n");
	double buffOut, buffIn;

	for (int i=0;i<n_;i++)
	{
		if (parent_[i]!=-2)
		{
			msg.append(identifier_[i]);
			msg.append("\t");

			buffIn = double(deltaT_[i])/(CLOCKS_PER_SEC);
			MPI_Reduce(&buffIn, &buffOut, 1, MPI_DOUBLE, MPI_SUM,0 , MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
			strs << (buffOut)/(numprocs);
			msg.append(strs.str());
			msg.append("\t");
			strs.str("");

			//buffIn = double(deltaT[i])/(CLOCKS_PER_SEC);
			MPI_Reduce(&buffIn, &buffOut, 1, MPI_DOUBLE, MPI_MAX,0 , MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
			strs << (buffOut);
			msg.append(strs.str());
			msg.append("\t");
			strs.str("");

			strs << nOfRuns_[i];
			msg.append(strs.str());
			msg.append("\t");
			strs.str("");

			strs << level_[i];
			msg.append(strs.str());
			msg.append("\t");
			strs.str("");

			strs << parent_[i];
			msg.append(strs.str());
			msg.append("\t");
			strs.str("");

			if (parent_[i] >= 0)
			{
				strs << identifier_[parent_[i]];
			}
			else
			{
				strs << "none";
			}

			msg.append(strs.str());
			msg.append("\n");
			strs.str("");
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	if(myrank == 0)
	{
        std::string fileName(path_/"timeEvalFull.txt");
		outFile.open(fileName.data(),ios_base::app);
		outFile << msg;
		outFile.close();
	}

	return;
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
// * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
clockModel::clockModel
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    dict_(dict),
    particleCloud_(sm),
    path_("clockData"),
    n_(30),
    deltaT_(std::vector<clock_t> (n_)),
    identifier_(std::vector<std::string> (n_)),
    nOfRuns_(std::vector<int> (n_)),
    level_(std::vector<short> (n_)),
    curLev_(0),
    parent_(std::vector<int> (n_)),
    curParent_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

clockModel::~clockModel()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
