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
    startTime_(propsDict_.lookupOrDefault<scalar>("startTime",0.)),
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
    scalarField2_
    (   
        IOobject
        (
            "cellVolume",
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
    path_("postProcessing/particleCellVolume"),
    sPtr_(NULL),
    writeToFile_(propsDict_.lookupOrDefault<Switch>("writeToFile",false))
{
    // init force sub model
    setForceSubModels(propsDict_);

    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(3,true); // activate search for verbose switch

    // read those switches defined above, if provided in dict
    forceSubM(0).readSwitches();

    particleCloud_.checkCG(true);

    // create the path and generate output file
    if(writeToFile_)
    {
        path_=particleCloud_.IOM().createTimeDir(path_);
	    sPtr_ = new OFstream(path_/"particleCellVolume.txt");
        *sPtr_ << "# time | total particle volume in cells | total volume of cells with particles | average volume fraction | min(voidfraction) | max(voidfraction)" << endl;
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
        if(forceSubM(0).verbose()) Info << "particleCellVolume.C - setForce()" << endl;

        scalarField_ == dimensionedScalar("zero", scalarField_.dimensions(), 0.);

        // get reference to actual field
        const volScalarField& field = mesh_.lookupObject<volScalarField>(scalarFieldName_);

        scalar fieldValue=-1;
        scalar cellVol=-1;
        scalar minFieldVal=1e18;
        scalar maxFieldVal=-1e18;

        forAll(field,cellI)
        {
            fieldValue = field[cellI];
            if(fieldValue < upperThreshold_ && fieldValue > lowerThreshold_)
            {
                cellVol = mesh_.V()[cellI];
                scalarField_[cellI] = (1-fieldValue) * cellVol;
                scalarField2_[cellI] = cellVol;
                minFieldVal=min(minFieldVal,fieldValue);
                maxFieldVal=max(maxFieldVal,fieldValue);
            }
            else
            {
                scalarField_[cellI] = 0.;
                scalarField2_[cellI] = 0.;
            }
        }
        scalarField_ == dimensionedScalar("zero", scalarField_.dimensions(), gSum(scalarField_));
        scalarField2_ == dimensionedScalar("zero", scalarField_.dimensions(), gSum(scalarField2_));
        reduce(minFieldVal, minOp<scalar>());
        reduce(maxFieldVal, maxOp<scalar>());

        if(forceSubM(0).verbose())
        {
           Info << "calculated integral particle volume "
                << " = " << scalarField_[0]
                << ",\n considering cells where the field < " << upperThreshold_
                << ", and > " << lowerThreshold_
                << ",\n the total volume of cells holding particles = " << scalarField2_[0]
                << ",\n this results in an average volume fraction of:" << scalarField_[0]/(scalarField2_[0]+SMALL)
                << ",\n the min occurring " << scalarFieldName_ << " is:" << minFieldVal
                << ",\n the max occurring " << scalarFieldName_ << " is:" << maxFieldVal
                << endl;
        }

        // write to file
        if(writeToFile_)
        {
            *sPtr_<< mesh_.time().value() << " " ;
            *sPtr_<< scalarField_[0] << " " ;
            *sPtr_<< scalarField2_[0] << " " ;
            *sPtr_<< scalarField_[0]/(scalarField2_[0]+SMALL) << " " ;
            *sPtr_<< minFieldVal << " " ;
            *sPtr_<< maxFieldVal << endl;
        }
    }// end if time >= startTime_
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
