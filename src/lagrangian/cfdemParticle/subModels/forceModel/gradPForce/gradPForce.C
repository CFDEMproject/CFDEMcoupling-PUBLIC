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

#include "gradPForce.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(gradPForce, 0);

addToRunTimeSelectionTable
(
    forceModel,
    gradPForce,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
gradPForce::gradPForce
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    pFieldName_(propsDict_.lookup("pFieldName")),
    p_(sm.mesh().lookupObject<volScalarField> (pFieldName_)),
    velocityFieldName_(propsDict_.lookup("velocityFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velocityFieldName_)),
    useRho_(false),
    useU_(false),
    addedMassCoeff_(0.0)
{
    // init force sub model
    setForceSubModels(propsDict_);

    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(0,true); // activate treatExplicit switch
    forceSubM(0).setSwitchesList(1,true); // activate treatForceDEM switch
    forceSubM(0).setSwitchesList(4,true); // activate search for interpolate switch

    // read those switches defined above, if provided in dict
    forceSubM(0).readSwitches();

    if (modelType_ == "B")
    {
        FatalError <<"using  model gradPForce with model type B is not valid\n" << abort(FatalError);
    }else if (modelType_ == "Bfull")
    {
        if(forceSubM(0).switches()[1])
        {
            Info << "Using treatForceDEM false!" << endl;
            forceSubM(0).setSwitches(1,false); // treatForceDEM = false
        }
    }else // modelType_=="A"
    {
        if(!forceSubM(0).switches()[1])
        {
            Info << "Using treatForceDEM true!" << endl;
            forceSubM(0).setSwitches(1,true); // treatForceDEM = true
        }
    }

    if (propsDict_.found("useU")) useU_=true;
    if (propsDict_.found("useAddedMass")) 
    {
        addedMassCoeff_ =  readScalar(propsDict_.lookup("useAddedMass"));
        Info << "gradP will also include added mass with coefficient: " << addedMassCoeff_ << endl;
        Info << "WARNING: use fix nve/sphere/addedMass in LIGGGHTS input script to correctly account for added mass effects!" << endl;
    }

    if(p_.dimensions()==dimensionSet(0,2,-2,0,0))
        useRho_ = true;

    particleCloud_.checkCG(true);

    particleCloud_.probeM().initialize(typeName, "gradP.logDat");
    particleCloud_.probeM().vectorFields_.append("gradPForce"); //first entry must the be the force
    particleCloud_.probeM().scalarFields_.append("Vs");
    particleCloud_.probeM().scalarFields_.append("rho");
    particleCloud_.probeM().writeHeader();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

gradPForce::~gradPForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void gradPForce::setForce() const
{
    volVectorField gradPField = fvc::grad(p_);
    /*if (useU_)
    {
        // const volScalarField& voidfraction_ = particleCloud_.mesh().lookupObject<volScalarField> ("voidfraction");
        volScalarField U2 = U_&U_;// *voidfraction_*voidfraction_;
        if (useRho_)
            gradPField = fvc::grad(0.5*U2);
        else
            gradPField = fvc::grad(0.5*forceSubM(0).rhoField()*U2);
    }*/
    vector gradP;
    scalar Vs;
    scalar rho;
    vector position;
    vector force;
    label cellI;

    interpolationCellPoint<vector> gradPInterpolator_(gradPField);

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

                if(forceSubM(0).interpolation()) // use intepolated values for alpha (normally off!!!)
                {
                    gradP = gradPInterpolator_.interpolate(position,cellI);
                }else
                {
                    gradP = gradPField[cellI];
                }

                Vs = particleCloud_.particleVolume(index);
                rho = forceSubM(0).rhoField()[cellI];

                // calc particle's pressure gradient force
                if (useRho_)
                    force = -Vs*gradP*rho*(1.0+addedMassCoeff_);
                else
                    force = -Vs*gradP*(1.0+addedMassCoeff_);

                if(forceSubM(0).verbose() && index >=0 && index <2)
                {
                    Info << "index = " << index << endl;
                    Info << "gradP = " << gradP << endl;
                    Info << "force = " << force << endl;
                }

                //Set value fields and write the probe
                if(probeIt_)
                {
                    #include "setupProbeModelfields.H"
                    vValues.append(force);           //first entry must the be the force
                    sValues.append(Vs);
                    sValues.append(rho);
                    particleCloud_.probeM().writeProbe(index, sValues, vValues);
                }
            }

            // write particle based data to global array
            forceSubM(0).partToArray(index,force,vector::zero);

        //}
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
