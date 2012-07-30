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
    IOModel

SourceFiles
    IOModel.C

\*---------------------------------------------------------------------------*/

#ifndef IOModel_H
#define IOModel_H

#include "fvCFD.H"
#include "cfdemCloud.H"
#include "OFstream.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class IOModel Declaration
\*---------------------------------------------------------------------------*/

class IOModel
{

protected:

    // Protected data
        const dictionary& dict_;

        cfdemCloud& particleCloud_;

        const Time& time_;

public:

    //- Runtime type information
    TypeName("IOModel");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            IOModel,
            dictionary,
            (
                const dictionary& dict,
                cfdemCloud& sm
            ),
            (dict,sm)
        );


    // Constructors

        //- Construct from components
        IOModel
        (
            const dictionary& dict,
            cfdemCloud& sm
        );


    // Destructor

        virtual ~IOModel();


    // Selector

        static autoPtr<IOModel> New
        (
            const dictionary& dict,
            cfdemCloud& sm
        );


    // Member Functions

    virtual void dumpDEMdata() const;

    fileName createTimeDir(fileName) const;

    fileName createLagrangianDir(fileName) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
