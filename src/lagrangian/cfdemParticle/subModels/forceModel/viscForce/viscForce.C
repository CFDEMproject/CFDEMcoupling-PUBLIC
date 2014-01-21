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

#include "viscForce.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(viscForce, 0);

addToRunTimeSelectionTable
(
    forceModel,
    viscForce,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
viscForce::viscForce
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    verbose_(false),
    velocityFieldName_(propsDict_.lookup("velocityFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velocityFieldName_)),
    densityFieldName_(propsDict_.lookup("densityFieldName")),
    rho_(sm.mesh().lookupObject<volScalarField> (densityFieldName_)),
    interpolation_(false)
{
    if (modelType_ == "B")
    {
        FatalError <<"using  model viscForce with model type B is not valid\n" << abort(FatalError);
    }else
    {
        treatDEM_=true;
        Info << "viscForce is applied only to DEM side" << endl;
    }
    if (propsDict_.found("verbose")) verbose_=true;
    if (propsDict_.found("treatExplicit")) treatExplicit_=true;
    if (propsDict_.found("interpolation"))
    {
        Info << "using interpolated value of pressure gradient." << endl;
        interpolation_=true;
    }
    particleCloud_.checkCG(true);

    //Append the field names to be probed
    particleCloud_.probeM().initialize(typeName, "visc.logDat");
    particleCloud_.probeM().vectorFields_.append("viscForce"); //first entry must the be the force
    particleCloud_.probeM().scalarFields_.append("Vs");
    particleCloud_.probeM().writeHeader();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

viscForce::~viscForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void viscForce::setForce() const
{

    // get viscosity field
    #ifdef comp
        const volScalarField& mufField = particleCloud_.turbulence().mu();

        // calc div(Tau)
        volVectorField divTauField =
        - fvc::laplacian(mufField, U_)
        - fvc::div(mufField*dev(fvc::grad(U_)().T()));
    #else
        const volScalarField& nufField = particleCloud_.turbulence().nu();

        // calc div(Tau)
        volVectorField divTauField =
        - fvc::laplacian(nufField*rho_, U_)
        - fvc::div(nufField*rho_*dev(fvc::grad(U_)().T()));
    #endif

    vector divTau;
    scalar ds;
    scalar Vs;
    vector position;
    vector force;
    label cellI;

    interpolationCellPoint<vector> divTauInterpolator_(divTauField);

    #include "setupProbeModel.H"

    for(int index = 0;index <  particleCloud_.numberOfParticles(); index++)
    {
        //if(mask[index][0])
        //{
            force=vector(0,0,0);
            cellI = particleCloud_.cellIDs()[index][0];

            if (cellI > -1) // particle Found
            {

                position = particleCloud_.position(index);

                if(interpolation_) // use intepolated values for alpha (normally off!!!)
                {
                    divTau = divTauInterpolator_.interpolate(position,cellI);
                }else
                {
                    divTau = divTauField[cellI];
                }

                ds = 2*particleCloud_.radius(index);
                Vs = ds*ds*ds*M_PI/6;

                // calc the contribution of the deviatoric stress 
                // to the generalized buoyancy force
                force = -Vs*divTau;

                if(verbose_ && index >0 && index <2)
                {
                    Info << "index = " << index << endl;
                    Info << "gradP = " << divTau << endl;
                    Info << "force = " << force << endl;
                }

                //Set value fields and write the probe
                if(probeIt_)
                {
                    #include "setupProbeModelfields.H"
                    vValues.append(force);  //first entry must the be the force
                    sValues.append(Vs);
                    particleCloud_.probeM().writeProbe(index, sValues, vValues);
                }
            }

            // set force on particle
            if(!treatDEM_){
                if(!treatExplicit_) for(int j=0;j<3;j++) impForces()[index][j] += force[j];
                else  for(int j=0;j<3;j++) expForces()[index][j] += force[j];
            }
             for(int j=0;j<3;j++) DEMForces()[index][j] += force[j];

        //}
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
