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

#include "volWeightedAverage.H"
#include "addToRunTimeSelectionTable.H"
#include "dataExchangeModel.H"
#include "mpi.h"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(volWeightedAverage, 0);

addToRunTimeSelectionTable
(
    forceModel,
    volWeightedAverage,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
volWeightedAverage::volWeightedAverage
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    mesh_(particleCloud_.mesh()),
    startTime_(propsDict_.lookupOrDefault<scalar>("startTime",0.)),
    scalarFieldNames_(propsDict_.lookup("scalarFieldNames")),
    vectorFieldNames_(propsDict_.lookup("vectorFieldNames")),
    volumeFractionName_(propsDict_.lookupOrDefault<word>("volumeFractionName","voidfraction")),
    volumeFraction_(particleCloud_.mesh().lookupObject<volScalarField>(volumeFractionName_)),
    upperThreshold_(readScalar(propsDict_.lookup("upperThreshold"))),
    lowerThreshold_(readScalar(propsDict_.lookup("lowerThreshold"))),
    useVolumeFraction_(propsDict_.lookupOrDefault<Switch>("useVolumeFraction",false)),
    verbose_(propsDict_.lookupOrDefault<Switch>("verbose",false)),
    path_("postProcessing/volWeightedAverage"),
    sPtr_(NULL),
    writeToFile_(propsDict_.lookupOrDefault<Switch>("writeToFile",false))
{
    // create vol weighted average scalar fields
    scalarFields_.setSize(scalarFieldNames_.size());

    for (int i=0;i < scalarFieldNames_.size(); i++)
    {
        word fieldName = "volAverage_" + scalarFieldNames_[i];

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

    // create vol weighted average vector fields
    vectorFields_.setSize(vectorFieldNames_.size());

    for (int i=0;i < vectorFieldNames_.size(); i++)
    {
        word fieldName = "volAverage_" + vectorFieldNames_[i];

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

    // create the path and generate output file
    if(writeToFile_)
    {
        path_=particleCloud_.IOM().createTimeDir(path_);
	    sPtr_ = new OFstream(path_/"volWeightedAverage.txt");
        //*sPtr_ << "time | vol av. scalar field i | ... | vol av. vector field i" << nl;
    }

    particleCloud_.checkCG(true);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

volWeightedAverage::~volWeightedAverage()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void volWeightedAverage::setForce() const
{
    if(mesh_.time().value() >= startTime_)
    {
        if(verbose_) Info << "volWeightedAverage.C - setForce()" << endl;

        // write to file
        if(writeToFile_) *sPtr_<< mesh_.time().value() << " " ;

        for (int i=0;i < scalarFieldNames_.size(); i++)
        {
            // get reference to actual field
            const volScalarField& field = mesh_.lookupObject<volScalarField>(scalarFieldNames_[i]);

            scalar fieldValue=-1.;
            scalar volWeightedAverage=-1.;
            scalar cellVol=-1.;
            scalar totVol=0.;
            scalar totVol_all=0.;
            scalar integralValue=0.;
            scalar volumeFraction(1.);

            forAll(field,cellI)
            {
                fieldValue = field[cellI];
                if(fieldValue < upperThreshold_ && fieldValue > lowerThreshold_)
                {
                    cellVol = mesh_.V()[cellI];
                    if(useVolumeFraction_)
                        volumeFraction = volumeFraction_[cellI];
                    scalarFields_[i][cellI] = fieldValue * cellVol * volumeFraction;
                    totVol += cellVol;
                }
                else
                {
                    scalarFields_[i][cellI] = 0.;
                }
            }

            MPI_Allreduce(&totVol, &totVol_all, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            integralValue = gSum(scalarFields_[i]);
            volWeightedAverage = integralValue / (totVol_all+SMALL);
            scalarFields_[i] == dimensionedScalar("value", scalarFields_[i].dimensions(), volWeightedAverage);

            if(verbose_)
            {
                Info << "calculated vol. weighted average of field: " << scalarFieldNames_[i]
                     << " = " << volWeightedAverage
                     << ",\n integral value "
                     << " = " << integralValue
                     << ",\n volume accounted for "
                     << " = " << totVol_all
                     << ",\n considering cells where the field < " << upperThreshold_
                     << ", and > " << lowerThreshold_ << endl;
            }

            // write to file
            if(writeToFile_) *sPtr_<< volWeightedAverage << " " ;
        }

        for (int i=0;i < vectorFieldNames_.size(); i++)
        {
            // get reference to actual field
            const volVectorField& field = mesh_.lookupObject<volVectorField>(vectorFieldNames_[i]);

            vector fieldValue(-1.,-1.,-1.);
            vector volWeightedAverage(-1.,-1.,-1.);
            scalar magvolWeightedAverage=-1.;
            scalar cellVol=-1.;
            scalar totVol=0.;
            scalar totVol_all=0.;
            scalar volumeFraction(1.);

            forAll(field,cellI)
            {
                fieldValue = field[cellI];
                magvolWeightedAverage = mag(fieldValue);
                if(magvolWeightedAverage < upperThreshold_ && magvolWeightedAverage > lowerThreshold_)
                {
                    cellVol = mesh_.V()[cellI];
                    if(useVolumeFraction_)
                        volumeFraction = volumeFraction_[cellI];
                    vectorFields_[i][cellI] = fieldValue * cellVol * volumeFraction;
                    totVol += cellVol;
                }
                else
                {
                    vectorFields_[i][cellI] = vector::zero;
                }
            }

            MPI_Allreduce(&totVol, &totVol_all, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            volWeightedAverage = gSum(vectorFields_[i]) / (totVol_all+SMALL);
            vectorFields_[i] == volWeightedAverage;

            if(verbose_)
            {
                Info << "calculated vol. weighted average of field: " << vectorFieldNames_[i]
                     << " = " << volWeightedAverage
                     << ",\n considering cells where the mag(field) < " << upperThreshold_
                     << ", and > " << lowerThreshold_ << endl;
            }

            // write to file
            if(writeToFile_) *sPtr_<< volWeightedAverage << " ";
        }

        // write to file
        if(writeToFile_) *sPtr_<< endl;

    }// end if time >= startTime_
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
