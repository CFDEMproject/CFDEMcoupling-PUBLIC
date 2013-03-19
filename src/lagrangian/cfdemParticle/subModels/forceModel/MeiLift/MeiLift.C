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

#include "MeiLift.H"
#include "addToRunTimeSelectionTable.H"

//#include "mpi.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(MeiLift, 0);

addToRunTimeSelectionTable
(
    forceModel,
    MeiLift,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
MeiLift::MeiLift
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
    rho_(sm.mesh().lookupObject<volScalarField> (densityFieldName_))/*,
    vorticityFieldName_(propsDict_.lookup("vorticityFieldName")),
    vorticity_(sm.mesh().lookupObject<volVectorField> (vorticityFieldName_))*/
{
    if (propsDict_.found("verbose")) verbose_=true;
    if (propsDict_.found("treatExplicit")) treatExplicit_=true;
    particleCloud_.checkCG(false);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

MeiLift::~MeiLift()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void MeiLift::setForce() const
{
    // get viscosity field
    #ifdef comp
        const volScalarField nufField = particleCloud_.turbulence().mu() / rho_;
    #else
        const volScalarField& nufField = particleCloud_.turbulence().nu();
    #endif

    vector lift(0,0,0);
    vector Us(0,0,0);
    vector Ur(0,0,0);
    scalar magUr(0);
    scalar ds(0);
    scalar nuf(0);
    scalar rho(0);
    scalar voidfraction(1);
    scalar Rep(0);
    scalar Rew(0);
    scalar Cl(0);
    scalar Cl_star(0);
    scalar J_star(0);
    scalar Omega_eq(0);
    scalar omega_star(0);
    vector vorticity(0,0,0);
    volVectorField vorticityField = fvc::curl(U_);


    for(int index = 0;index <  particleCloud_.numberOfParticles(); index++)
    {
        //if(mask[index][0])
        //{
            lift=vector::zero;
            label cellI = particleCloud_.cellIDs()[index][0];

            if (cellI > -1) // particle Found
            {
                //NP note: one could add pointInterpolated values instead of cell centered
                Us = particleCloud_.velocity(index);
                Ur = U_[cellI]-Us;
                magUr = mag(Ur);
                vorticity=vorticityField[cellI];

                if (magUr > 0 && mag(vorticity) > 0)
                {
                    ds = 2*particleCloud_.radius(index);
                    nuf = nufField[cellI];
                    rho = rho_[cellI];
                    voidfraction = particleCloud_.voidfraction(index);
                    omega_star=mag(vorticity)*ds/magUr;

                    // calc particle Re Nr
                    Rep = ds*magUr/nuf;
		    Rew = mag(vorticity)*ds*ds/nuf;

                    Omega_eq = omega_star/2.0*(1.0-0.0075*Rew)*(1.0-0.062*sqrt(Rep)-0.001*Rep);
                    J_star = 0.3*(1.0+tanh(2.5*(log10(sqrt(omega_star/Rep))+0.191)))
                             *(2.0/3.0+tanh(6.0*sqrt(omega_star/Rep)-1.92));
                    Cl_star=1.0-(0.675+0.15*(1.0+tanh(0.28*(omega_star/2.0-2.0))))*tanh(0.18*sqrt(Rep));
                    Cl=J_star*12.92/M_PI*sqrt(omega_star/Rep)+Omega_eq*Cl_star;
                    lift = 0.125*rho*M_PI*Cl*magUr*Ur^vorticity/mag(vorticity)*ds*ds;

                    if (modelType_=="B")
                        lift /= voidfraction;
                }

                if(verbose_ && index >100 && index <102)
                {
                    Pout << "index = " << index << endl;
                    Pout << "Us = " << Us << endl;
                    Pout << "Ur = " << Ur << endl;
                    Pout << "ds = " << ds << endl;
                    Pout << "rho = " << rho << endl;
                    Pout << "nuf = " << nuf << endl;
                    Pout << "Rep = " << Rep << endl;
                    Pout << "lift = " << lift << endl;
                }
            }
            // set force on particle
            if(!treatDEM_){
                if(!treatExplicit_) for(int j=0;j<3;j++) impForces()[index][j] += lift[j];
                else  for(int j=0;j<3;j++) expForces()[index][j] += lift[j];
            }
            for(int j=0;j<3;j++) DEMForces()[index][j] += lift[j];
        //}
    }

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
