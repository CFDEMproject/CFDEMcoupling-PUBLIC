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

#include "fieldTimeAverage.H"
#include "addToRunTimeSelectionTable.H"
#include "dataExchangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(fieldTimeAverage, 0);

addToRunTimeSelectionTable
(
    forceModel,
    fieldTimeAverage,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
fieldTimeAverage::fieldTimeAverage
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    mesh_(particleCloud_.mesh())
{
    init();
    particleCloud_.checkCG(true);
}

// Construct from components
fieldTimeAverage::fieldTimeAverage
(
    const dictionary& dict,
    const dictionary& mydict,
    const word& fieldPrefix,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(mydict),
    mesh_(particleCloud_.mesh())
{
    init(fieldPrefix);
    particleCloud_.checkCG(true);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

fieldTimeAverage::~fieldTimeAverage()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void fieldTimeAverage::init(word fieldPrefix)
{
    startTime_=0.;
    scalarFieldNames_=wordList(propsDict_.lookup("scalarFieldNames"));
    vectorFieldNames_=wordList(propsDict_.lookup("vectorFieldNames"));
    nrAverages_=0.0;

    // create time average scalar fields
    scalarFields_.setSize(scalarFieldNames_.size());

    if (propsDict_.found("startTime"))
        startTime_=readScalar(propsDict_.lookup("startTime"));

    for (int i=0;i < scalarFieldNames_.size(); i++)
    {
        word fieldName = fieldPrefix + scalarFieldNames_[i];

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
                    IOobject::READ_IF_PRESENT,
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
        word fieldName = fieldPrefix + vectorFieldNames_[i];

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
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedVector("0", mesh_.lookupObject<volVectorField>(vectorFieldNames_[i]).dimensions(), vector::zero)
            )
        );
    }
}

void fieldTimeAverage::setForce() const
{
    if(mesh_.time().value() >= startTime_)
    {
        if(particleCloud_.verbose()) Info << "fieldTimeAverage.C - setForce()" << endl;

        for (int i=0;i < scalarFieldNames_.size(); i++)
        {
            // get reference to actual field
            const volScalarField& field = mesh_.lookupObject<volScalarField>(scalarFieldNames_[i]);
            averaging(field,i,nrAverages_);
        }

        for (int i=0;i < vectorFieldNames_.size(); i++)
        {
            // get reference to actual field
            const volVectorField& field = mesh_.lookupObject<volVectorField>(vectorFieldNames_[i]);
            averaging(field,i,nrAverages_);
        }
        nrAverages_++;
    }// end if time >= startTime_
}

void fieldTimeAverage::averaging(const volScalarField& field,int& i,double& nrAverages_) const
{
    // first entry in this field
    if(nrAverages_ == 0)
        scalarFields_[i] = field;
    else
        scalarFields_[i] = (scalarFields_[i]*nrAverages_+field*1)/(nrAverages_+1);
}

void fieldTimeAverage::averaging(const volVectorField& field,int& i,double& nrAverages_) const
{
    // first entry in this field
    if(nrAverages_ == 0)
        vectorFields_[i] = field;
    else
        vectorFields_[i] = (vectorFields_[i]*nrAverages_+field*1)/(nrAverages_+1);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
