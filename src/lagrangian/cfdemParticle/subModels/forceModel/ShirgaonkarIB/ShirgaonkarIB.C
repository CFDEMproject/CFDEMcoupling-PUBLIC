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

#include "ShirgaonkarIB.H"
#include "addToRunTimeSelectionTable.H"
#include "voidFractionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(ShirgaonkarIB, 0);

addToRunTimeSelectionTable
(
    forceModel,
    ShirgaonkarIB,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
ShirgaonkarIB::ShirgaonkarIB
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    verbose_(false),
    twoDimensional_(false),
    depth_(1),
    velFieldName_(propsDict_.lookup("velFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    densityFieldName_(propsDict_.lookup("densityFieldName")),
    rho_(sm.mesh().lookupObject<volScalarField> (densityFieldName_)),
    pressureFieldName_(propsDict_.lookup("pressureFieldName")),
    p_(sm.mesh().lookupObject<volScalarField> (pressureFieldName_))
{
    if (propsDict_.found("verbose")) verbose_=true;
    if (propsDict_.found("twoDimensional"))
    {
        twoDimensional_=true;
        depth_ = propsDict_.lookup("depth");
        Info << "2-dimensional simulation - make sure DEM side is 2D" << endl;
        Info << "depth of domain is assumed to be :" << depth_ << endl;
    }
    if (propsDict_.found("treatExplicit")) treatExplicit_=true;
    particleCloud_.checkCG(false);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

ShirgaonkarIB::~ShirgaonkarIB()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ShirgaonkarIB::setForce() const
{

    label cellI;
    vector drag;

    #ifdef comp
        // get viscosity field
        const volScalarField& mufField = particleCloud_.turbulence().mu();
        volVectorField h = (mufField*fvc::laplacian(U_)-fvc::grad(p_));
    #else
        // get viscosity field
        const volScalarField& nufField = particleCloud_.turbulence().nu();
        volVectorField h = rho_*(nufField*fvc::laplacian(U_)-fvc::grad(p_));
    #endif

    for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {
        //if(mask[index][0])
        //{
            drag=vector::zero;

            for(int subCell=0;subCell<particleCloud_.voidFractionM().cellsPerParticle()[index][0];subCell++)
            {
                //Info << "subCell=" << subCell << endl;
                cellI = particleCloud_.cellIDs()[index][subCell];

                if (cellI > -1) // particle Found
                {
                    drag += h[cellI]*h.mesh().V()[cellI];
                }

            }

            // set force on particle
            if(twoDimensional_) drag /= depth_;
            if(treatExplicit_) for(int j=0;j<3;j++) expForces()[index][j] += drag[j];
            else  for(int j=0;j<3;j++) impForces()[index][j] += drag[j];
            for(int j=0;j<3;j++) DEMForces()[index][j] += drag[j];

            if(verbose_) Info << "impForces = " << impForces()[index][0]<<","<<impForces()[index][1]<<","<<impForces()[index][2] << endl;
        //}
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
