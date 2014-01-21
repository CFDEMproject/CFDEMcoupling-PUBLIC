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
    probeModel

SourceFiles
    probeModel.C

\*---------------------------------------------------------------------------*/

#ifndef probeModel_H
#define probeModel_H

#include "fvCFD.H"
#include "cfdemCloud.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class probeModel Declaration
\*---------------------------------------------------------------------------*/

class probeModel
{

protected:

    // Protected data
        const dictionary& dict_;

        cfdemCloud& particleCloud_;

        //- Class used for grouping field types
        template<class Type>
        class fieldGroup
        :
            public DynamicList<word>
        {
        public:
            //- Construct null
            fieldGroup()
            :
                DynamicList<word>(0)
            {}

        };


public:

    // Variables
    //- Categorized scalar/vector/tensor vol fields
    mutable fieldGroup<scalar> scalarFields_;
    mutable fieldGroup<vector> vectorFields_;

    //- Runtime type information
    TypeName("probeModel");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            probeModel,
            dictionary,
            (
                const dictionary& dict,
                cfdemCloud& sm,
                word   typeName,
                char*  logFileName
            ),
            (dict,sm,typeName,logFileName)
        );


    // Constructors

        //- Construct from components
        probeModel
        (
            const dictionary& dict,
            cfdemCloud& sm,
            word   typeName,
            char*  logFileName
        );


    // Destructor

        virtual ~probeModel();


    // Selector

        static autoPtr<probeModel> New
        (
            const dictionary& dict,
            cfdemCloud& sm,
            word   typeName,
            char*  logFileName
        );

    // Member Functions

    virtual void initialize(word typeName, word  logFileName) const {};
    virtual void setOutputFile() const {};
    virtual void writeHeader() const {};
    virtual void writeProbe(int index, Field<scalar> sValues, Field<vector> vValues) const {};
    virtual bool checkIDForPrint(int) const {return false;};
    virtual void setCounter() const {};
    virtual bool active() const {return true;}; 
    const char* wordToChar(word&) const;

    // Access

};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
