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
    noDEMForce_(propsDict_.lookupOrDefault<Switch>("noDEMForce",false)),
    keepCFDForce_(propsDict_.lookupOrDefault<Switch>("keepCFDForce",false))
{
    Info << "noDragTestMe: " << noDEMForce_ << " " << keepCFDForce_ << endl;

    particleCloud_.checkCG(true);

    // init force sub model
    setForceSubModels(propsDict_);

    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(3,true); // activate search for verbose switch

    //set default switches (hard-coded default = false)
    //forceSubM(0).setSwitches(XXX,true);

    // read those switches defined above, if provided in dict
    for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
        forceSubM(iFSub).readSwitches();

    // setup required communication
    for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
        forceSubM(iFSub).setupCommunication();
}

void noDrag::MSinit()
{
    // NOTE: all MS related operations need to be performed during init of MS cloud
    if (particleCloud_.shapeTypeName() == "multisphere" && !particleBased_)
    {
        Info << type() << ": activating multisphere mode..." << endl;
        forceSubM(0).setSwitches(11,true); // this is a MS model

        // re-setup required communication for MS
        for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
            forceSubM(iFSub).setupCommunication();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

noDrag::~noDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void noDrag::setForce() const
{
    if (forceSubM(0).verbose())
        Info << "noDrag::setForce:"
             << "  noDEMForce=" << noDEMForce_
             << "  , keepCFDForce=" << keepCFDForce_
             << endl;

    label cellI=0;
    int idDragExp=0;
    int idKsl=0;
    int idUf=0;
    if (forceSubM(0).ms())
    {
        idDragExp=particleCloud_.idDragExpCM();
        idKsl=particleCloud_.idKslCM();
        idUf=particleCloud_.idUfCM();
    }
    else
    {
        idDragExp=particleCloud_.idDragExp();
        idKsl=particleCloud_.idKsl();
        idUf=particleCloud_.idUf();
    }


    for (int index = 0;index < particleCloud_.numberOfObjects(particleBased_); index++)
    {
        cellI = particleCloud_.cellIDs(particleBased_)[index][0];
        if (cellI > -1) // particle Found
        {
            //==========================
            // set force on particle (proposed new code)
            // write particle based data to global array
            //forceSubM(0).partToArray(index,drag,dragExplicit);
            //==========================
            // set force on particle (old code)
            if (!keepCFDForce_)
            {
                for (int j=0;j<3;j++)
                {
                    expForces()[index][j] = 0.;
                    impForces()[index][j] = 0.;
                }
            }
            if (noDEMForce_)
            {
                for (int j=0;j<3;j++)
                    particleCloud_.fieldsToDEM[idDragExp][index][j] = 0.;

                if (particleCloud_.impDEMdrag())
                {
                    particleCloud_.fieldsToDEM[idKsl][index][0] = 0.;
                    for (int j=0;j<3;j++)
                        particleCloud_.fieldsToDEM[idUf][index][j] = 0.;
                }
            }
            //==========================
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
