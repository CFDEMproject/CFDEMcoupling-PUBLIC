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

#include "fieldStore.H"
#include "addToRunTimeSelectionTable.H"
#include "dataExchangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(fieldStore, 0);

addToRunTimeSelectionTable
(
    forceModel,
    fieldStore,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
fieldStore::fieldStore
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    mesh_(particleCloud_.mesh()),
    scalarFieldNames_(propsDict_.lookup("scalarFieldNames")),
    vectorFieldNames_(propsDict_.lookup("vectorFieldNames"))
{
    // create time average scalar fields
    scalarFields_.setSize(scalarFieldNames_.size());

    for (int i=0;i < scalarFieldNames_.size(); i++)
    {
        word fieldName = "stored_" + scalarFieldNames_[i];

        Info<< "Creating field " << fieldName << endl;
        scalarFields_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    fieldName,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("0", mesh_.lookupObject<volScalarField>(scalarFieldNames_[i]).dimensions(), 0)
            )
        );
    }

    // create time average vector fields
    vectorFields_.setSize(vectorFieldNames_.size());

    for (int i=0;i < vectorFieldNames_.size(); i++)
    {
        word fieldName = "stored_" + vectorFieldNames_[i];

        Info<< "Creating field " << fieldName << endl;
        vectorFields_.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    fieldName,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedVector("0", mesh_.lookupObject<volVectorField>(vectorFieldNames_[i]).dimensions(), vector::zero)
            )
        );
    }

    particleCloud_.checkCG(true);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

fieldStore::~fieldStore()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fieldStore::setForce() const
{
    if(particleCloud_.mesh().time().outputTime())
    {
        if(particleCloud_.verbose()) Info << "fieldStore.C - setForce()" << endl;

        for (int i=0;i < scalarFieldNames_.size(); i++)
        {
            // get reference to actual field
            const volScalarField& field = mesh_.lookupObject<volScalarField>(scalarFieldNames_[i]);

            // save field
            scalarFields_[i] = field;
        }

        for (int i=0;i < vectorFieldNames_.size(); i++)
        {
            // get reference to actual field
            const volVectorField& field = mesh_.lookupObject<volVectorField>(vectorFieldNames_[i]);

            // save field
            vectorFields_[i] = field;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
