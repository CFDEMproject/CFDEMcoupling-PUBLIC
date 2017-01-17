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

#include "error.H"

#include "particleVolume.H"
#include "addToRunTimeSelectionTable.H"
#include "dataExchangeModel.H"
#include "mpi.h"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(particleVolume, 0);

addToRunTimeSelectionTable
(
    forceModel,
    particleVolume,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
particleVolume::particleVolume
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    startTime_(propsDict_.lookupOrDefault<scalar>("startTime",0.)),
    scaleDia_(1.),
    path_("postProcessing/particleVolume"),
    sPtr_(NULL),
    writeToFile_(propsDict_.lookupOrDefault<Switch>("writeToFile",false))
{
    // init force sub model
    setForceSubModels(propsDict_);

    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(3,true); // activate search for verbose switch

    // read those switches defined above, if provided in dict
    forceSubM(0).readSwitches();

    if (propsDict_.found("scale"))
        scaleDia_=scalar(readScalar(propsDict_.lookup("scale")));

    particleCloud_.checkCG(true);

    // create the path and generate output file
    if(writeToFile_)
    {
        path_=particleCloud_.IOM().createTimeDir(path_);
	    sPtr_ = new OFstream(path_/"particleVolume.txt");
        //*sPtr_ << "time | total particle volume" << nl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

particleVolume::~particleVolume()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void particleVolume::setForce() const
{
    if(particleCloud_.mesh().time().value() >= startTime_)
    {
        if (scaleDia_ > 1)
            Info << typeName << " using scale = " << scaleDia_ << endl;
        else if (particleCloud_.cg() > 1){
            scaleDia_=particleCloud_.cg();
            Info << typeName << " using scale from liggghts cg = " << scaleDia_ << endl;
        }

        if(forceSubM(0).verbose()) Info << "particleVolume.C - setForce()" << endl;

        scalar ds(0);
        scalar VsTot(0);
        scalar fpth=4.188790204786390525;//4./3.*pi;
        
        for(int index = 0;index <  particleCloud_.numberOfParticles(); ++index)
        {
            if(particleCloud_.cellIDs()[index][0] >= 0)
            {
                ds = 2*particleCloud_.radius(index)/scaleDia_;
                VsTot += ds*ds*ds*fpth;
            }
        }

        reduce(VsTot, sumOp<scalar>());
        if(forceSubM(0).verbose()) Info << "Total particle volume (located in domain) = " << VsTot << endl;

        // write to file
        if(writeToFile_)
        {
            *sPtr_<< particleCloud_.mesh().time().value() << " " ;
            *sPtr_<< VsTot << endl;
        }
    }// end if time >= startTime_
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
