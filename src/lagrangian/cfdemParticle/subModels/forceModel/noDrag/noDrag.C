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

#include "noDrag.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(noDrag, 0);

addToRunTimeSelectionTable
(
    forceModel,
    noDrag,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
noDrag::noDrag
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    noDEMForce_(propsDict_.lookupOrDefault("noDEMForce",false)),
    keepCFDForce_(propsDict_.lookupOrDefault("keepCFDForce",false))
{
    // init force sub model
    setForceSubModels(propsDict_);

    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(3,true); // activate search for verbose switch

    // read those switches defined above, if provided in dict
    forceSubM(0).readSwitches();
    particleCloud_.checkCG(true);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

noDrag::~noDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void noDrag::setForce() const
{
    if(forceSubM(0).verbose())
    {
        Info << "noDrag::setForce:" << endl;
        Info << "noDEMForce=" << noDEMForce_ << endl;
        Info << "keepCFDForce=" << keepCFDForce_ << endl;
    }

    label cellI=0;
    for(int index = 0;index <  particleCloud_.numberOfParticles(); ++index)
    {
        cellI = particleCloud_.cellIDs()[index][0];
        if (cellI > -1) // particle Found
        {
            //==========================
            // set force on particle (proposed new code)
            // write particle based data to global array
            //forceSubM(0).partToArray(index,drag,dragExplicit);
            //==========================
            // set force on particle (old code)
            if(!keepCFDForce_)
            {
                for(int j=0;j<3;j++)
                {
                    expForces()[index][j] = 0.;
                    impForces()[index][j] = 0.;
                }
            }
            if(noDEMForce_)
            {
                for(int j=0;j<3;j++) DEMForces()[index][j] = 0.;
                if(particleCloud_.impDEMdrag())
                {
                    Cds()[index][0] = 0.;
                    for(int j=0;j<3;j++) fluidVel()[index][j] = 0.;
                }
            }
            //==========================
        }        
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
