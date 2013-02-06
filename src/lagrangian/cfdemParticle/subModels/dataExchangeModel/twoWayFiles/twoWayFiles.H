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

    two way DEM-CFD coupling via files

Class
    twoWayFiles

SourceFiles
    twoWayFiles.C

\*---------------------------------------------------------------------------*/

#ifndef twoWayFiles_H
#define twoWayFiles_H

#include "dataExchangeModel.H"
#include "OFstream.H"
#include "sys/stat.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class twoWayFiles Declaration
\*---------------------------------------------------------------------------*/

class twoWayFiles
:
    public dataExchangeModel
{
private:

  // private data
    dictionary propsDict_;

  // private member functions

    const char* wordToChar(word&) const;

    const char* fileNameToChar(fileName&) const;

    fileName getFilePath(word&, bool) const;

    void renameFilePath(fileName&,word&) const;

public:

    //- Runtime type information
    TypeName("twoWayFiles");


    // Constructors

        //- Construct from components
        twoWayFiles
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~twoWayFiles();


    // Member Functions
        void getData
        (
            word name,
            word type,
            double ** const& field,
            label step
        ) const;

        void getData
        (
            word name,
            word type,
            int ** const& field,
            label step
        ) const {FatalError<<"function not implemented !!! twoWayFiles::getData!!!\n" << abort(FatalError);};

        void giveData
        (
            word name,
            word type,
            double ** const& field,
            const char* datatype = ""
        ) const;

        word myType() const{return typeName; };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
