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

#include "particleCellVolume.H"
#include "addToRunTimeSelectionTable.H"
#include "dataExchangeModel.H"
#include "mpi.h"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(particleCellVolume, 0);

addToRunTimeSelectionTable
(
    forceModel,
    particleCellVolume,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
particleCellVolume::particleCellVolume
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    mesh_(particleCloud_.mesh()),
    startTime_(0.),
    scalarFieldName_("voidfraction"),
    scalarField_
    (   
        IOobject
        (
            "particleCellVolume",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0), 0)
    ),
    upperThreshold_(readScalar(propsDict_.lookup("upperThreshold"))),
    lowerThreshold_(readScalar(propsDict_.lookup("lowerThreshold"))),
    verbose_(false)
{
    if (propsDict_.found("startTime")){
        startTime_=readScalar(propsDict_.lookup("startTime"));
    }

    if (propsDict_.found("verbose")){
        verbose_ = true;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

particleCellVolume::~particleCellVolume()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void particleCellVolume::setForce() const
{
    if(mesh_.time().value() >= startTime_)
    {
        if(verbose_) Info << "particleCellVolume.C - setForce()" << endl;

        scalarField_.internalField()=0.;

        // get reference to actual field
        volScalarField& field = (volScalarField&) mesh_.lookupObject<volScalarField>(scalarFieldName_);

        scalar fieldValue=-1;
        scalar cellVol=-1;

        forAll(field,cellI)
        {
            fieldValue = field[cellI];
            if(fieldValue < upperThreshold_ && fieldValue > lowerThreshold_)
            {
                cellVol = mesh_.V()[cellI];
                scalarField_[cellI] = (1-fieldValue) * cellVol;
            }
            else
            {
                scalarField_[cellI] = 0.;
            }
        }
        scalarField_.internalField() = gSum(scalarField_);

        if(verbose_)
        {
           Info << "calculated integral of field: " << scalarFieldName_
                << " = " << scalarField_[0]
                << ",\n considering cells where the field < " << upperThreshold_
                << ", and > " << lowerThreshold_ << endl;
        }
    }// end if time >= startTime_
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
