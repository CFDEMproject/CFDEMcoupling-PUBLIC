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

#include "interface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(interface, 0);

addToRunTimeSelectionTable
(
    forceModel,
    interface,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
interface::interface
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    VOFvoidfractionFieldName_(propsDict_.lookup("VOFvoidfractionFieldName")),
    alpha_(sm.mesh().lookupObject<volScalarField> (VOFvoidfractionFieldName_)),
    gradAlphaName_(propsDict_.lookup("gradAlphaName")),
    gradAlpha_(sm.mesh().lookupObject<volVectorField> (gradAlphaName_)),
    sigma_(readScalar(propsDict_.lookup("sigma"))),
    theta_(readScalar(propsDict_.lookup("theta"))),
    alphaThreshold_(readScalar(propsDict_.lookup("alphaThreshold"))),
    deltaAlphaIn_(readScalar(propsDict_.lookup("deltaAlphaIn"))),
    deltaAlphaOut_(readScalar(propsDict_.lookup("deltaAlphaOut"))),
    C_(1.0),
    interpolation_(false),
    alphaInterpolator_(interpolation<scalar>::New("cellPoint", alpha_)),
    gradAlphaInterpolator_(interpolation<vector>::New("cellPoint", gradAlpha_))
{
    if (propsDict_.found("C")) C_=readScalar(propsDict_.lookup("C"));
    if (propsDict_.found("interpolation")) interpolation_=true;
    if (propsDict_.found("treatExplicit")) treatExplicit_=true;

    Info << "check if interpolation really works - use directly interpolationCellPoint<vector> ???" << endl;
    particleCloud_.checkCG(false);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

interface::~interface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void interface::setForce() const
{
Info << "interface::setForce" << endl;
    for(int index = 0;index <  particleCloud_.numberOfParticles(); ++index)
    {
        //if(mask[index][0])
        //{
            // definition of spherical particle
            scalar dp = 2*particleCloud_.radius(index);
            vector position = particleCloud_.position(index);
            label cellI = particleCloud_.cellIDs()[index][0];

            if(cellI >-1.0) // particle found on proc domain
            {
                scalar alphap;
                vector magGradAlphap;

                if(interpolation_) // use intepolated values for alpha (normally off!!!)
                {
                    // make interpolation object for alpha
                    alphap = alphaInterpolator_().interpolate(position,cellI);

                    // make interpolation object for grad(alpha)/|grad(alpha)|
                    vector gradAlphap = gradAlphaInterpolator_().interpolate(position,cellI);
                    magGradAlphap = gradAlphap/max(mag(gradAlphap),SMALL);
                }
                else // use cell centered values for alpha
                {
                    //// for any reason fvc::grad(alpha_) cannot be executed here!?
                    //volVectorField gradAlpha=fvc::grad(alpha_);
                    //volVectorField a = gradAlpha/
                    //                   max(mag(gradAlpha),dimensionedScalar("a",dimensionSet(0,-1,0,0,0), SMALL));
                    //magGradAlphap = a[cellI];

                    alphap = alpha_[cellI];
                    volVectorField a = gradAlpha_/
                                       max(mag(gradAlpha_),dimensionedScalar("a",dimensionSet(0,-1,0,0,0), SMALL));
                    magGradAlphap = a[cellI];
                }

                // Initialize an interfaceForce vector
                vector interfaceForce = Foam::vector(0,0,0);

                // Calculate the interfaceForce (range of alphap needed for stability)

                if ((alphaThreshold_-deltaAlphaIn_) < alphap && alphap < (alphaThreshold_+deltaAlphaOut_))
                {
                    Info << "within threshold limits" << endl;
                    // Calculate estimate attachment force as
                    // |6*sigma*sin(pi-theta/2)*sin(pi+theta/2)|*2*pi*dp
                    scalar Fatt =   mag(
                                   6
                                 * sigma_
                                 * sin(M_PI - theta_/2)
                                 * sin(M_PI + theta_/2)
                            )
                      * M_PI
                      * dp;

                    interfaceForce = - magGradAlphap
                         * tanh(alphap-alphaThreshold_)
                         * Fatt
                         * C_;
                }

                if(true && mag(interfaceForce) > 0)
                {
                Info << "dp = " << dp << endl;
                Info << "position = " << position << endl;
                Info << "cellI = " << cellI << endl;
                Info << "alpha cell = " << alpha_[cellI] << endl;
                Info << "alphap = " << alphap << endl;
                Info << "magGradAlphap = " << magGradAlphap << endl;
                Info << "interfaceForce = " << interfaceForce << endl;
                Info << "mag(interfaceForce) = " << mag(interfaceForce) << endl;
                }

                // limit interface force
                /*scalar rhoP=3000;
                scalar mP=dp*dp*dp*3.1415/4*rhoP;
                scalar fMax=5*mP*9.81;
                if(mag(interfaceForce)>fMax){
                    interfaceForce /= mag(interfaceForce)/fMax;
                    Info << "interface force is limited to " << interfaceForce << endl;
                }*/

               if(treatExplicit_) for(int j=0;j<3;j++) expForces()[index][j] += interfaceForce[j];
               else  for(int j=0;j<3;j++) impForces()[index][j] += interfaceForce[j];
               for(int j=0;j<3;j++) DEMForces()[index][j] += interfaceForce[j];
            } // end if particle found on proc domain
        //}// end if in mask
    }// end loop particles
Info << "interface::setForce - done" << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
