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
    liggghtsCommandModel

SourceFiles
    liggghtsCommandModel.C

\*---------------------------------------------------------------------------*/

#ifndef liggghtsCommandModel_H
#define liggghtsCommandModel_H

#include "fvCFD.H"
#include "cfdemCloud.H"
#include "dataExchangeModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class liggghtsCommandModel Declaration
\*---------------------------------------------------------------------------*/

class liggghtsCommandModel
{

protected:

    // Protected data
    const dictionary& dict_;

    cfdemCloud& particleCloud_;

    string strCommand_;

    int nextRun_;

    int lastRun_;

    Switch runFirst_;

    Switch runLast_;

    Switch runEveryCouplingStep_;

    Switch runEveryWriteStep_;

    scalar startTime_;

    scalar endTime_;

    scalar timeInterval_;

    int firstCouplingStep_;

    int lastCouplingStep_;

    int couplingStepInterval_;

public:

    //- Runtime type information
    TypeName("liggghtsCommandModel");

    // Declare runtime constructor selection table


        declareRunTimeSelectionTable
        (
            autoPtr,
            liggghtsCommandModel,
            dictionary,
            (
                const dictionary& dict,
                cfdemCloud& sm,
                int i
            ),
            (dict,sm,i)
        );


    // Constructors

        //- Construct from components
        liggghtsCommandModel
        (
            const dictionary& dict,
            cfdemCloud& sm,
            int i
        );


    // Destructor

        virtual ~liggghtsCommandModel();


    // Selector

        static autoPtr<liggghtsCommandModel> New
        (
            const dictionary& dict,
            cfdemCloud& sm,
            word liggghtsCommandType,
            int i
        );


    // Member Functions

        virtual const char* command()=0;

        void checkTimeMode(dictionary&);

        void checkTimeSettings(const dictionary&);

        virtual bool runCommand(int)=0;

        bool runThisCommand(int);

        string addTimeStamp(word);

    // Access

        int nextRun(){return nextRun_;};

        int lastRun(){return lastRun_;};

        virtual word name()=0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
