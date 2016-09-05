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
\*---------------------------------------------------------------------------*/

#include "scalarTransportModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(scalarTransportModel, 0);

defineRunTimeSelectionTable(scalarTransportModel, dictionary);

// * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * //
void scalarTransportModel::createFields()
{}

void scalarTransportModel::update()
{}

const volScalarField& scalarTransportModel::sourceField()
{
    // dummy source field
    tmp<volScalarField> tsource
    (
        volScalarField
        (
            IOobject
            (
                "tsource",
                particleCloud_.mesh().time().timeName(),
                particleCloud_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            particleCloud_.mesh(),
            dimensionedScalar
            (
                "zero",
                dimensionSet(0,0,0,0,0),
                0
            )
        )
    );
    return tsource();
}

// * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from Components
scalarTransportModel::scalarTransportModel
(
    const dictionary& dict, //not used at the moment
    cfdemCloud& sm
)
:
    particleCloud_(sm),
    mesh_(sm.mesh()),
    scalarTransportProperties_
    (
        IOobject
        (
            "scalarTransportProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    verbose_(false),
    ignore_(false)
{
    createFields();
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

scalarTransportModel::~scalarTransportModel()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
