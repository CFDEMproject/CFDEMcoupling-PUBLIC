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

#include "SchillerNaumannDrag.H"
#include "addToRunTimeSelectionTable.H"

//#include "mpi.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(SchillerNaumannDrag, 0);

addToRunTimeSelectionTable
(
    forceModel,
    SchillerNaumannDrag,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
SchillerNaumannDrag::SchillerNaumannDrag
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
    rho_(sm.mesh().lookupObject<volScalarField> (densityFieldName_))
{
    if (propsDict_.found("verbose")) verbose_=true;
    if (propsDict_.found("treatExplicit")) treatExplicit_=true;
    particleCloud_.checkCG(false);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

SchillerNaumannDrag::~SchillerNaumannDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void SchillerNaumannDrag::setForce() const
{
    // get viscosity field
    #ifdef comp
        const volScalarField nufField = particleCloud_.turbulence().mu() / rho_;
    #else
        const volScalarField& nufField = particleCloud_.turbulence().nu();
    #endif

    for(int index = 0;index <  particleCloud_.numberOfParticles(); index++)
    {
        //if(mask[index][0])
        //{
            vector drag(0,0,0);
            label cellI = particleCloud_.cellIDs()[index][0];

            if (cellI > -1) // particle Found
            {
                //NP note: one could add pointInterpolated values instead of cell centered
                vector Us = particleCloud_.velocity(index);
                vector Ur = U_[cellI]-Us;
                scalar ds = 2*particleCloud_.radius(index);
                scalar nuf = nufField[cellI];
                scalar rho = rho_[cellI];
                scalar voidfraction = particleCloud_.voidfraction(index);
                scalar magUr = mag(Ur);
                scalar Rep = 0;
                scalar Cd = 0;

                if (magUr > 0)
                {
                   // calc particle Re Nr
                    Rep = ds*magUr/nuf;

                    // calc fluid drag Coeff
                    Cd = max(0.44,24.0/Rep*(1.0+0.15*pow(Rep,0.687)));

                    // calc particle's drag
                    drag = 0.125*Cd*rho*M_PI*ds*ds*magUr*Ur;

                    if (modelType_=="B")
                        drag /= voidfraction;
                }

                if(verbose_ && index >100 && index <102)
                {
                    Info << "index = " << index << endl;
                    Info << "Us = " << Us << endl;
                    Info << "Ur = " << Ur << endl;
                    Info << "ds = " << ds << endl;
                    Info << "rho = " << rho << endl;
                    Info << "nuf = " << nuf << endl;
                    Info << "voidfraction = " << voidfraction << endl;
                    Info << "Rep = " << Rep << endl;
                    Info << "Cd = " << Cd << endl;
                    Info << "drag = " << drag << endl;
                }
            }
            // set force on particle
            if(treatExplicit_) for(int j=0;j<3;j++) expForces()[index][j] += drag[j];
            else  for(int j=0;j<3;j++) impForces()[index][j] += drag[j];
            for(int j=0;j<3;j++) DEMForces()[index][j] += drag[j];
        //}
    }

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
