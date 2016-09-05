/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
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

    transport model for scalar properties

Class
    temperatureModel

SourceFiles
    temperatureModel.C

\*---------------------------------------------------------------------------*/

#ifndef temperatureModel_H
#define temperatureModel_H

#include "scalarTransportModel.H"
#include "forceModel.H"
#include "eulerianScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class temperatureModel Declaration
\*---------------------------------------------------------------------------*/

class temperatureModel
:
    public scalarTransportModel
{

    dictionary                      propsDict_;

    autoPtr<eulerianScalarField>    temperatureField_;

    scalar                          PrT_; // - Turbulent Prandtl number (optional, default = 0.7)

public:

    //- Runtime type information
    TypeName("temperatureModel");


    // Constructors

    //- Construct from components
    temperatureModel
    (
            const dictionary& dict,
            cfdemCloud& sm
    );

    // Destructor

    	~temperatureModel();


    // Member Functions

    void createFields();

    void update();

    const volScalarField& sourceField();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
