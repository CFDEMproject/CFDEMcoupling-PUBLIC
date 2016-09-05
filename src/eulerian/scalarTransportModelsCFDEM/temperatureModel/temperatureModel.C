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

#include "temperatureModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(temperatureModel, 0);

addToRunTimeSelectionTable
(
	scalarTransportModel,
	temperatureModel,
	dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
// Construct from components
temperatureModel::temperatureModel
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    scalarTransportModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    PrT_(0.7)
{
    propsDict_.readIfPresent("PrT", PrT_);

    Info << "Using PrT " << PrT_ << endl;

    temperatureField_ = eulerianScalarField::New
    (
         propsDict_,
         sm,
         "T",
         0
    );
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
temperatureModel::~temperatureModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void temperatureModel::createFields()
{}

// ************************************************************
void temperatureModel::update()
{

    //==============================
    // get references
    const surfaceScalarField& phi(particleCloud_.mesh().lookupObject<surfaceScalarField> ("phi"));
    const volScalarField& voidfraction(particleCloud_.mesh().lookupObject<volScalarField> ("voidfraction"));
    //==============================
    temperatureField_->update(phi, voidfraction, particleCloud_.turbulence().nuEff(), PrT_);

}

// ************************************************************
const volScalarField& temperatureModel::sourceField()
{
    return temperatureField_->mSource();
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
