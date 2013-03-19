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

#include "ArchimedesIB.H"
#include "addToRunTimeSelectionTable.H"
#include "voidFractionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(ArchimedesIB, 0);

addToRunTimeSelectionTable
(
    forceModel,
    ArchimedesIB,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
ArchimedesIB::ArchimedesIB
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    twoDimensional_(false),
    densityFieldName_(propsDict_.lookup("densityFieldName")),
    rho_(sm.mesh().lookupObject<volScalarField> (densityFieldName_)),
    voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")), //mod by alice
    voidfractions_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),//mod by alice
    gravityFieldName_(propsDict_.lookup("gravityFieldName")),
    #if defined(version21) || defined(version16ext)
        g_(sm.mesh().lookupObject<uniformDimensionedVectorField> (gravityFieldName_))
    #elif  defined(version15)
        g_(dimensionedVector(sm.mesh().lookupObject<IOdictionary>("environmentalProperties").lookup(gravityFieldName_)).value())
    #endif
{
    if (propsDict_.found("twoDimensional"))
    {
        twoDimensional_=true;
        Info << "2-dimensional simulation - make sure DEM side is 2D" << endl;
    }

    if (propsDict_.found("treatExplicit")) treatExplicit_=true;
    treatDEM_=true;
    Info << "accounting for Archimedes only on DEM side!" << endl;
    particleCloud_.checkCG(true);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

ArchimedesIB::~ArchimedesIB()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ArchimedesIB::setForce() const
{
    vector force;
    for(int index = 0;index <  particleCloud_.numberOfParticles(); ++index)
    {
        //if(mask[index][0])
        //{
            force=vector::zero;
            for(int subCell=0;subCell<particleCloud_.voidFractionM().cellsPerParticle()[index][0];subCell++)
            {
                label cellI = particleCloud_.cellIDs()[index][subCell];
                if (cellI > -1) // particle Found
                {
                    //force += -g_.value()*rho_[cellI]*rho_.mesh().V()[cellI]*(1-particleCloud_.voidfractions()[index][subCell]);//mod by alice
                	force += -g_.value()*rho_[cellI]*rho_.mesh().V()[cellI]*(1-voidfractions_[cellI]);//mod by alice
        	    }
            }
            // set force on particle
            if(twoDimensional_) Warning<<"ArchimedesIB model doesn't work for 2D right now!!\n"<< endl;
            if(!treatDEM_)
            {
                if(treatExplicit_) for(int j=0;j<3;j++) expForces()[index][j] += force[j];
                else for(int j=0;j<3;j++) impForces()[index][j] += force[j];
            }
            for(int j=0;j<3;j++) DEMForces()[index][j] += force[j];
        //}
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
