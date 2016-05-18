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

    calc time average of scalar or vector field

Class
    fieldTimeAverage

SourceFiles
    fieldTimeAverage.C

\*---------------------------------------------------------------------------*/

#ifndef fieldTimeAverage_H
#define fieldTimeAverage_H

#include "forceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fieldTimeAverage Declaration
\*---------------------------------------------------------------------------*/

class fieldTimeAverage
:
    public forceModel
{
protected:

    dictionary propsDict_;

    const fvMesh& mesh_;

    scalar startTime_;

    wordList scalarFieldNames_;

    wordList vectorFieldNames_;

    mutable PtrList<volScalarField> scalarFields_;

    mutable PtrList<volVectorField> vectorFields_;

    mutable double nrAverages_;

    void init(word fieldPrefix = word("timeAverage_"));

public:

    //- Runtime type information
    TypeName("fieldTimeAverage");


    // Constructors

        //- Construct from components
        fieldTimeAverage
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

        //- Construct from components
        fieldTimeAverage
        (
            const dictionary& dict,
            const dictionary& mydict,
            const word& fieldPrefix,
            cfdemCloud& sm
        );

    // Destructor

        ~fieldTimeAverage();


    // Member Functions
        void setForce() const;

        virtual inline void averaging(const volScalarField&,int&,double&) const;

        virtual inline void averaging(const volVectorField&,int&,double&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
