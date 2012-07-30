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

#include "KochHillDrag.H"
#include "addToRunTimeSelectionTable.H"

#include "mpi.h"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(KochHillDrag, 0);

addToRunTimeSelectionTable
(
    forceModel,
    KochHillDrag,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
KochHillDrag::KochHillDrag
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
    voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    interpolation_(false)
{
    if (propsDict_.found("verbose")) verbose_=true;
    if (propsDict_.found("treatExplicit")) treatExplicit_=true;
    if (propsDict_.found("interpolation")) interpolation_=true;
    if (propsDict_.found("scale"))
        scale_=scalar(readScalar(propsDict_.lookup("scale")));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

KochHillDrag::~KochHillDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void KochHillDrag::setForce
(
    double** const& mask,
    double**& impForces,
    double**& expForces,
    double**& DEMForces
) const
{
    // get viscosity field
    #ifdef comp
        const volScalarField nufField = particleCloud_.turbulence().mu()/rho_;
    #else
        const volScalarField& nufField = particleCloud_.turbulence().nu();
    #endif

    vector position(0,0,0);
    scalar voidfraction(1);
    vector Ufluid(0,0,0);
    vector drag(0,0,0);
    label cellI=0;

    vector Us(0,0,0);
    vector Ur(0,0,0);
    scalar ds(0);
    scalar nuf(0);
    scalar rho(0);
    scalar magUr(0);
    scalar Rep(0);
	scalar Vs(0);
	scalar volumefraction(0);

    interpolationCellPoint<scalar> voidfractionInterpolator_(voidfraction_);
    interpolationCellPoint<vector> UInterpolator_(U_);

    for(int index = 0;index <  particleCloud_.numberOfParticles(); index++)
    {
        if(mask[index][0])
        {

            cellI = particleCloud_.cellIDs()[index][0];
            drag = vector(0,0,0);

            if (cellI > -1) // particle Found
            {
                if(interpolation_)
                {
	                position = particleCloud_.position(index);
                    voidfraction = voidfractionInterpolator_.interpolate(position,cellI);
                    Ufluid = UInterpolator_.interpolate(position,cellI);
                }else
                {
					voidfraction = particleCloud_.voidfraction(index);
                    Ufluid = U_[cellI];
                }

                Us = particleCloud_.velocity(index);
                Ur = Ufluid-Us;
                ds = 2*particleCloud_.radius(index);
                nuf = nufField[cellI];
                rho = rho_[cellI];
                magUr = mag(Ur);
				Rep = 0;
                Vs = ds*ds*ds*M_PI/6;
                volumefraction = 1-voidfraction+SMALL;

                if (magUr > 0)
                {
                    // calc particle Re Nr
                    Rep = ds/scale_*voidfraction*magUr/(nuf+SMALL);

                    // calc model coefficient F0
                    scalar F0=0.;
                    if(volumefraction < 0.4)
                    {
                        F0 = (1+3*sqrt((volumefraction)/2)+135/64*volumefraction*log(volumefraction)+16.14*volumefraction)/
                             (1+0.681*volumefraction-8.48*sqr(volumefraction)+8.16*volumefraction*volumefraction*volumefraction);
                    } else {
                        F0 = 10*volumefraction/(voidfraction*voidfraction*voidfraction);
                    }

                    // calc model coefficient F3
                    scalar F3 = 0.0673+0.212*volumefraction+0.0232/pow(voidfraction,5);

                    // calc model coefficient beta
                    scalar beta = 18*nuf*rho*voidfraction*voidfraction*volumefraction/(ds/scale_*ds/scale_)*
                                  (F0 + 0.5*F3*Rep);

                    // calc particle's drag
                    drag = Vs*beta/volumefraction*Ur;

                    if (modelType_=="B")
                        drag /= voidfraction;
                }

                if(verbose_ && index >=0 && index <2)
                {
                    Info << "index = " << index << endl;
                    Info << "Us = " << Us << endl;
                    Info << "Ur = " << Ur << endl;
                    Info << "ds = " << ds << endl;
                    Info << "ds/scale = " << ds/scale_ << endl;
                    Info << "rho = " << rho << endl;
                    Info << "nuf = " << nuf << endl;
                    Info << "voidfraction = " << voidfraction << endl;
                    Info << "Rep = " << Rep << endl;
                    Info << "drag = " << drag << endl;
                }
            }
            // set force on particle
            if(treatExplicit_) for(int j=0;j<3;j++) expForces[index][j] += drag[j];
            else  for(int j=0;j<3;j++) impForces[index][j] += drag[j];
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
