/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
                                Copyright (C) 2013-     Graz University of  
                                                        Technology, IPPT
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

#include "constDiffSmoothing.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(constDiffSmoothing, 0);

addToRunTimeSelectionTable
(
    smoothingModel,
    constDiffSmoothing,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
constDiffSmoothing::constDiffSmoothing
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    smoothingModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    lowerLimit_(readScalar(propsDict_.lookup("lowerLimit"))),
    upperLimit_(readScalar(propsDict_.lookup("upperLimit"))),
    smoothingLength_(dimensionedScalar("smoothingLength",dimensionSet(0,1,0,0,0,0,0), readScalar(propsDict_.lookup("smoothingLength")))),
    smoothingLengthReferenceField_(dimensionedScalar("smoothingLengthReferenceField",dimensionSet(0,1,0,0,0,0,0), readScalar(propsDict_.lookup("smoothingLength")))),
    DT_("DT", dimensionSet(0,2,-1,0,0), 0.),
    verbose_(false)
{

    if(propsDict_.found("verbose"))  
        verbose_ = true;

    if(propsDict_.found("smoothingLengthReferenceField"))  
       smoothingLengthReferenceField_.value() = double(readScalar(propsDict_.lookup("smoothingLengthReferenceField")));

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

constDiffSmoothing::~constDiffSmoothing()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
bool constDiffSmoothing::doSmoothing() const
{
    return true;
}


void Foam::constDiffSmoothing::smoothen(volScalarField& fieldSrc) const
{
    // transfer data to working field to not mess up ddt
    volScalarField field=fieldSrc;
    field.correctBoundaryConditions();
    field.oldTime()=fieldSrc;
    field.oldTime().correctBoundaryConditions();

    double deltaT = field.mesh().time().deltaTValue();
    DT_.value() = smoothingLength_.value() * smoothingLength_.value() / deltaT;

    // do smoothing
    solve
    (
        fvm::ddt(field)
       -fvm::laplacian(DT_, field)
    );

    // bound field
    forAll(field,cellI)
    {
        field[cellI]=max(lowerLimit_,min(upperLimit_,field[cellI]));
    }  

    // get data from working field - will copy only values at new time
    fieldSrc=field;
    fieldSrc.correctBoundaryConditions(); 

    if(verbose_)
    {
        Info << "min/max(fieldoldTime) (unsmoothed): " << min(field.oldTime()) << tab << max(field.oldTime()) << endl;
        Info << "min/max(fieldSrc): " << min(fieldSrc) << tab << max(fieldSrc) << endl;
        Info << "min/max(fieldSrc.oldTime): " << min(fieldSrc.oldTime()) << tab << max(fieldSrc.oldTime()) << endl;
    }
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void Foam::constDiffSmoothing::smoothen(volVectorField& fieldSrc) const
{
    // transfer data to working field to not mess up ddt
    volVectorField field=fieldSrc;
    field.correctBoundaryConditions();
    field.oldTime()=fieldSrc;
    field.oldTime().correctBoundaryConditions();

    double deltaT = field.mesh().time().deltaTValue();
    DT_.value() = smoothingLength_.value() * smoothingLength_.value() / deltaT;

    // do smoothing
    solve
    (
        fvm::ddt(field)
       -fvm::laplacian(DT_, field)
    );

    // get data from working field
    fieldSrc=field;
    fieldSrc.correctBoundaryConditions(); 

    if(verbose_)
    {
        Info << "min/max(fieldoldTime) (unsmoothed): " << min(field.oldTime()) << tab << max(field.oldTime()) << endl;
        Info << "min/max(fieldSrc): " << min(fieldSrc) << tab << max(fieldSrc) << endl;
        Info << "min/max(fieldSrc.oldTime): " << min(fieldSrc.oldTime()) << tab << max(fieldSrc.oldTime()) << endl;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void Foam::constDiffSmoothing::smoothenReferenceField(volVectorField& fieldSrc) const
{
    // transfer data to working field to not mess up ddt
    volVectorField field=fieldSrc;
    field.correctBoundaryConditions();
    field.oldTime()=fieldSrc;
    field.oldTime().correctBoundaryConditions();

    double sourceStrength = 1e5; //large number to keep reference values constant

    dimensionedScalar deltaT = field.mesh().time().deltaT();
    DT_.value() = smoothingLengthReferenceField_.value() 
                         * smoothingLengthReferenceField_.value() / deltaT.value();

    tmp<volScalarField> NLarge
    (
        new volScalarField
        (
            IOobject
            (
                "xxx",
                particleCloud_.mesh().time().timeName(),
                particleCloud_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            particleCloud_.mesh(),
            0.0
        )
    );


    //loop over particles and map max particle diameter to Euler Grid
    forAll(field,cellI)
    {
        if ( mag(field.oldTime().internalField()[cellI]) > 0.0f)  // have a vector in the OLD field, so keep it!
            NLarge()[cellI] = sourceStrength;
    }

    // do the smoothing
    solve
    (
        fvm::ddt(field)
       -fvm::laplacian( DT_, field)
       == 
        NLarge() / deltaT * field.oldTime()  //add source to keep cell values constant
       -fvm::Sp( NLarge() / deltaT, field)   //add sink to keep cell values constant
    );

    // get data from working field
    fieldSrc=field;
    fieldSrc.correctBoundaryConditions(); 

    if(verbose_)
    {
        Info << "min/max(fieldoldTime) (unsmoothed): " << min(field.oldTime()) << tab << max(field.oldTime()) << endl;
        Info << "min/max(fieldSrc): " << min(fieldSrc) << tab << max(fieldSrc) << endl;
        Info << "min/max(fieldSrc.oldTime): " << min(fieldSrc.oldTime()) << tab << max(fieldSrc.oldTime()) << endl;
    }

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
