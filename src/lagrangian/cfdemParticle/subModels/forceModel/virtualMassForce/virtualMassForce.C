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

#include "virtualMassForce.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(virtualMassForce, 0);

addToRunTimeSelectionTable
(
    forceModel,
    virtualMassForce,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
virtualMassForce::virtualMassForce
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    verbose_(false),
    velFieldName_(propsDict_.lookup("velFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    densityFieldName_(propsDict_.lookup("densityFieldName")),
    rho_(sm.mesh().lookupObject<volScalarField> (densityFieldName_)),
    UrelOld_(NULL)
{
    if (propsDict_.found("verbose")) verbose_=true;

    if (particleCloud_.dataExchangeM().maxNumberOfParticles() > 0)
    {
        // get memory for 2d array
        particleCloud_.dataExchangeM().allocateArray(UrelOld_,0.,3);
    }
    if (propsDict_.found("treatExplicit")) treatExplicit_=true;
    particleCloud_.checkCG(true);

    //Append the field names to be probed
    particleCloud_.probeM().initialize(typeName, "virtualMass.logDat");
    particleCloud_.probeM().vectorFields_.append("virtualMassForce"); //first entry must the be the force
    particleCloud_.probeM().vectorFields_.append("Urel");
    particleCloud_.probeM().vectorFields_.append("UrelOld");
    particleCloud_.probeM().scalarFields_.append("Vs");
    particleCloud_.probeM().scalarFields_.append("rho");
    particleCloud_.probeM().writeHeader();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

virtualMassForce::~virtualMassForce()
{
    delete UrelOld_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void virtualMassForce::setForce() const
{
    reAllocArrays();

    scalar dt = U_.mesh().time().deltaT().value();

    #include "setupProbeModel.H"

    for(int index = 0;index <  particleCloud_.numberOfParticles(); index++)
    {
        //if(mask[index][0])
        //{
            vector virtualMassForce(0,0,0);
            label cellI = particleCloud_.cellIDs()[index][0];

            if (cellI > -1) // particle Found
            {
                vector Us = particleCloud_.velocity(index);
                vector Ur = U_[cellI]-Us;
                vector UrelOld;
                for(int j=0;j<3;j++)
                {
                    UrelOld[j] = UrelOld_[index][j];
                    UrelOld_[index][j] = Ur[j];
                }

                vector ddtUrel = (Ur-UrelOld)/dt;
                scalar ds = 2*particleCloud_.radius(index);
                scalar Vs = ds*ds*ds*M_PI/6;
                scalar rho = rho_[cellI];

                virtualMassForce = 0.5 * rho * Vs * ddtUrel;

                //Set value fields and write the probe
                if(probeIt_)
                {
                    #include "setupProbeModelfields.H"
                    vValues.append(virtualMassForce);           //first entry must the be the force
                    vValues.append(Ur);
                    vValues.append(UrelOld);
                    sValues.append(Vs);
                    sValues.append(rho);
                    particleCloud_.probeM().writeProbe(index, sValues, vValues);
                }
            }
            // set force on particle
            if(treatExplicit_) for(int j=0;j<3;j++) expForces()[index][j] += virtualMassForce[j];
            else  for(int j=0;j<3;j++) impForces()[index][j] += virtualMassForce[j];
            for(int j=0;j<3;j++) DEMForces()[index][j] += virtualMassForce[j];
        //}
    }

}

void Foam::virtualMassForce::reAllocArrays() const
{
    if(particleCloud_.numberOfParticlesChanged())
    {
        // get arrays of new length
        particleCloud_.dataExchangeM().allocateArray(UrelOld_,1.,1);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
