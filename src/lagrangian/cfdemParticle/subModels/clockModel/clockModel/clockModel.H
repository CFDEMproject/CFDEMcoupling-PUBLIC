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

Class
    clockModel

SourceFiles
    clockModel.C

\*---------------------------------------------------------------------------*/

#ifndef clockModel_H
#define clockModel_H

#define START(x) start(__COUNTER__,x)

#include "fvCFD.H"
#include "cfdemCloud.H"
#include "dataExchangeModel.H"

#include <vector>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class clockModel Declaration
\*---------------------------------------------------------------------------*/

class clockModel
{
protected:

    // Protected data
    const dictionary& dict_;
    cfdemCloud& particleCloud_;
    fileName path_;
    scalar startTime_;

    mutable int n_;
    mutable std::vector<clock_t> deltaT_;
    mutable std::vector<std::string> identifier_;
    mutable std::vector<int> nOfRuns_;
    mutable std::vector<short> level_;
    mutable short curLev_;
    mutable std::vector<int> parent_;
    mutable int curParent_;

public:

    //- Runtime type information
    TypeName("clockModel");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            clockModel,
            dictionary,
            (
                const dictionary& dict,
                cfdemCloud& sm
            ),
            (dict,sm)
        );


    // Constructors

        //- Construct from components
        clockModel
        (
            const dictionary& dict,
            cfdemCloud& sm
        );


    // Destructor

        virtual ~clockModel();


    // Selector

        static autoPtr<clockModel> New
        (
            const dictionary& dict,
            cfdemCloud& sm
        );


    // Member Functions

	    virtual void start(int pos) const;	                    //start measurement with custom string identifier
	    virtual void start(int pos,std::string identifier) const;	//start measurement with custom string identifier
	    virtual void stop() const;				                //stop last started measurement
	    virtual void stop(std::string identifier) const;	        //stop last started measurement with check if identifier is equal
	    virtual std::string eval() const;
	    virtual void evalFile() const;
	    virtual void evalPar() const;
	    void initElems();
	    std::vector<int> calcShift() const;			//detects empty indices in vector, when times are evaluated
	    void Hist() const;					//calc Histogram
	    virtual void normHist() const;				//calc normalized Histogram
        void plotHist(double,std::string,int,int) const;    //plot histogramm to terminal
    	void getRAMUsage() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
