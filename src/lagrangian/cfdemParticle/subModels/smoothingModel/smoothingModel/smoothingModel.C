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
\*---------------------------------------------------------------------------*/

#include "error.H"
#include "smoothingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(smoothingModel, 0);

defineRunTimeSelectionTable(smoothingModel, dictionary);

// * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * //

// * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::smoothingModel::smoothingModel
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    dict_(dict),
    particleCloud_(sm),
    vSmoothField_
    (   
        IOobject
        (
            "vSmoothField",
            particleCloud_.mesh().time().timeName(),
            particleCloud_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        particleCloud_.mesh(),
        dimensionedVector("zero", dimensionSet(0,0,0,0,0), vector::zero)
    ),
    sSmoothField_
    (   
        IOobject
        (
            "sSmoothField",
            particleCloud_.mesh().time().timeName(),
            particleCloud_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        particleCloud_.mesh(),
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0), 0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::smoothingModel::~smoothingModel()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void smoothingModel::checkFields(volScalarField& sSmoothField_) const
{
    // currently it is detected if field was auto generated or defined
    // improvement would be changing the type here automatically
    forAll(sSmoothField_.boundaryField(),patchI)
        if(sSmoothField_.boundaryField()[patchI].type()=="calculated")
            FatalError <<"Scalar field:"<< sSmoothField_.name() << " must be defined.\n" << abort(FatalError);

    sSmoothField_.writeOpt() = IOobject::AUTO_WRITE;
}

void smoothingModel::checkFields(volVectorField& vSmoothField_) const
{
    // currently it is detected if field was auto generated or defined
    // improvement would be changing the type here automatically
    forAll(vSmoothField_.boundaryField(),patchI)      
        if(vSmoothField_.boundaryField()[patchI].type()=="calculated")
            FatalError <<"Vector field:"<< vSmoothField_.name() << " must be defined.\n" << abort(FatalError);

    vSmoothField_.writeOpt() = IOobject::AUTO_WRITE;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
bool smoothingModel::doSmoothing() const
{
    return false;
}

void smoothingModel::dSmoothing() const
{
    // do nothing
    //dSmooth *= 0.0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void smoothingModel::smoothenAbsolutField(volScalarField& scalField) const
{

    //1 - First make the field volume-specific
    particleCloud_.makeSpecific(scalField);

    //2 - smoothen now the volume-specific field (the volume integral of this field will be conserved!)
    smoothen(scalField);

    //3 - Finally, make field absolute again
    particleCloud_.scaleWithVcell(scalField);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void smoothingModel::smoothenAbsolutField(volVectorField& vecField) const
{

    //1 - First make the field volume-specific
    particleCloud_.makeSpecific(vecField);

    //2 - smoothen now the volume-specific field (the volume integral of this field will be conserved!)
    smoothen(vecField);

    //3 - Finally, make field absolute again
    particleCloud_.scaleWithVcell(vecField);
}
} // End namespace Foam

// ************************************************************************* //
