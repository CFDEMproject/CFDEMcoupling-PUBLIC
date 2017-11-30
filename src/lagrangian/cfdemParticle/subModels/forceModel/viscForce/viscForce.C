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
    velocityFieldName_(propsDict_.lookup("velocityFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velocityFieldName_)),
    addedMassCoeff_(0.0)
{
    // block viscForceModel for type B
    if (modelType_ == "B") FatalError <<"using  model viscForce with model type B is not valid\n" << abort(FatalError);

    // init force sub model
    setForceSubModels(propsDict_);

    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(0,true); // activate search for treatForceExplicit switch
    forceSubM(0).setSwitchesList(3,true); // activate search for verbose switch
    forceSubM(0).setSwitchesList(4,true); // activate search for interpolate switch
    forceSubM(0).setSwitchesList(8,true); // activate scalarViscosity switch

    //set default switches (hard-coded default = false)
    forceSubM(0).setSwitches(0,true);       // make treatForceExplicit=true the default (is desired, otherwise this force would be implicit in slip vel!)
    if (modelType_ == "Bfull")              // type Bfull
        forceSubM(0).setSwitches(1,false);  // treatForceDEM = false
    else                                    // type A
        forceSubM(0).setSwitches(1,true);   // treatForceDEM = true

    // read user defined switches
    for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
        forceSubM(iFSub).readSwitches();

    if (propsDict_.found("useAddedMass")) 
    {
        addedMassCoeff_ =  readScalar(propsDict_.lookup("useAddedMass"));
        Info << "viscForce will also include added mass with coefficient: " << addedMassCoeff_ << endl;
        Info << "WARNING: use fix nve/sphere/addedMass in LIGGGHTS input script to correctly account for added mass effects!" << endl;
    }

    particleCloud_.checkCG(true);

    //Append the field names to be probed
    // suppress particle probe
    if (probeIt_ && propsDict_.found("suppressProbe"))
        probeIt_=!Switch(propsDict_.lookup("suppressProbe"));

    if(probeIt_)
    {
        particleCloud_.probeM().initialize(typeName, typeName+".logDat");
        particleCloud_.probeM().vectorFields_.append("viscForce"); //first entry must the be the force
        particleCloud_.probeM().scalarFields_.append("Vs");
        particleCloud_.probeM().writeHeader();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

viscForce::~viscForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void viscForce::setForce() const
{
    const volVectorField& divTau_ = forceSubM(0).divTauField(U_);

    vector divTau;
    scalar Vs;
    vector position;
    vector force;
    label cellI;

    #include "resetDivTauInterpolator.H"
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
                    divTau = divTauInterpolator_().interpolate(position,cellI);
                }else
                {
                    divTau = divTau_[cellI];
                }

                Vs = particleCloud_.particleVolume(index);

                // calc the contribution of the deviatoric stress 
                // to the generalized buoyancy force
                force = -Vs*divTau*(1.0+addedMassCoeff_);

                if(forceSubM(0).verbose() && index >0 && index <2)
                {
                    Info << "index = " << index << endl;
                    Info << "gradP = " << divTau << endl;
                    Info << "force = " << force << endl;
                }

                //Set value fields and write the probe
                if(probeIt_)
                {
                    #include "setupProbeModelfields.H"
                    // Note: for other than ext one could use vValues.append(x)
                    // instead of setSize
                    vValues.setSize(vValues.size()+1, force);           //first entry must the be the force
                    sValues.setSize(sValues.size()+1, Vs);
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
