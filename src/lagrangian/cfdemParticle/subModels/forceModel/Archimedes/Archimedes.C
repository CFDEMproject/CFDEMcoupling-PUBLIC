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

#include "Archimedes.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Archimedes, 0);

addToRunTimeSelectionTable
(
    forceModel,
    Archimedes,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Archimedes::Archimedes
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    twoDimensional_(false),
    gravityFieldName_(propsDict_.lookupOrDefault<word>("gravityFieldName","g")),
    #if defined(version21) || defined(version16ext)
        g_(sm.mesh().lookupObject<uniformDimensionedVectorField> (gravityFieldName_)),
    #elif defined(version15)
        g_(dimensionedVector(sm.mesh().lookupObject<IOdictionary>("environmentalProperties").lookup(gravityFieldName_)).value()),
    #endif
    rhoFieldName_(propsDict_.lookupOrDefault<word>("densityFieldName","rho")),
    rho_(sm.mesh().lookupObject<volScalarField> (rhoFieldName_))
{
    // suppress particle probe
    if (probeIt_ && propsDict_.found("suppressProbe"))
        probeIt_=!Switch(propsDict_.lookup("suppressProbe"));
    if(probeIt_)
    {
        particleCloud_.probeM().initialize(typeName, typeName+".logDat");
        particleCloud_.probeM().vectorFields_.append("archimedesForce");  //first entry must the be the force
        particleCloud_.probeM().scalarFields_.append("Vp");
        particleCloud_.probeM().writeHeader();
    }

    particleCloud_.checkCG(true);

    if (propsDict_.found("twoDimensional"))
    {
        twoDimensional_=true;
        Info << "2-dimensional simulation - make sure DEM side is 2D" << endl;
    }

    // init force sub model
    setForceSubModels(propsDict_);

    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(3,true); // activate search for verbose switch
    forceSubM(0).setSwitchesList(4,true); // activate search for interpolate switch

    //set default switches (hard-coded default = false)
    forceSubM(0).setSwitches(1,true);       // Archimedes only on DEM side (treatForceDEM=true)

    // read those switches defined above, if provided in dict
    for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
        forceSubM(iFSub).readSwitches();

    // setup required communication
    for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
        forceSubM(iFSub).setupCommunication();

    // sanity check
    if(!forceSubM(0).treatForceDEM()) FatalError << "You use 'treatForceDEM=false' for Archimedes, which is not correct." << abort(FatalError);
}

void Archimedes::MSinit()
{
    // NOTE: all MS related operations need to be performed during init of MS cloud
    if (particleCloud_.shapeTypeName() == "multisphere" && !particleBased_)
    {
        Info << type() << ": activating multisphere mode..." << endl;
        forceSubM(0).setSwitches(11,true); // this is a MS model

        // re-setup required communication for MS
        for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
            forceSubM(iFSub).setupCommunication();
    }
    else if((particleCloud_.shapeTypeName() == "superquadric"))
    {
        //set default switches (hard-coded default = false)
        forceSubM(0).setSwitches(16,true); // pullShape=true
        forceSubM(0).setSwitches(25,true); // superquadric=true
        forceSubM(0).setSwitches(26,true); // useQuaternion=true

        // read those switches defined above, if provided in dict  // NOT NECESSARY AS WE DO NOT ADD THINGS TO READABLE SWITCHES LIST
        //for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
        //    forceSubM(iFSub).readSwitches();

        // re-setup required communication for SQ
        for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
            forceSubM(iFSub).setupCommunication();
    }
    else if((particleCloud_.shapeTypeName() == "convex"))
    {
        forceSubM(0).setSwitches(29,true); // convex=true

        // read those switches defined above, if provided in dict  // NOT NECESSARY AS WE DO NOT ADD THINGS TO READABLE SWITCHES LIST
        //for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
        //    forceSubM(iFSub).readSwitches();

        // re-setup required communication for SQ
        for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
            forceSubM(iFSub).setupCommunication();
    }
    else
        // TODO NOTE: this model supposedly cannot run in a particle based mode
        // --> disallow for now
        FatalError << type() << ": You are trying to run this model in particle based operation mode. This has not been tested. Aborting..." << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Archimedes::~Archimedes()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Archimedes::setForce() const
{
    vector force(0,0,0);
    scalar piBySix(M_PI/6.);
    scalar Vs(0.);
    scalar ds(0.);
    vector position(0,0,0);
    label cellI=0;
    scalar RhoInt(0.);

    #include "resetRhoInterpolator.H"
    #include "setupProbeModel.H"

    for(int index = 0;index <  particleCloud_.numberOfObjects(particleBased_); index++)
    {
        cellI = particleCloud_.cfdemCloud::cellIDs(particleBased_)[index][0];
        force=vector::zero;

        if (cellI > -1) // particle Found
        {
            if(twoDimensional_)
            {
                scalar r = particleCloud_.radius(index);
                force = -g_.value()*forceSubM(0).rhoField()[cellI]*r*r*M_PI; // circle area
                Warning << "Archimedes::setForce() : this functionality is not tested!" << endl;
            }
            else
            {
                if (forceSubM(0).interpolation())
                {
                    position = particleCloud_.cfdemCloud::position(index, particleBased_);
                    RhoInt = RhoInterpolator_().interpolate(position, cellI);
                }
                else
                    RhoInt = forceSubM(0).rhoField()[cellI];

                if (forceSubM(0).sq() || forceSubM(0).convex())
                {
                    Vs = particleCloud_.volume(index);
                }
                else
                {
                    ds = particleCloud_.diameter(index, particleBased_);
                    scalar dParcel = ds;
                    forceSubM(0).scaleDia(dParcel,index); //caution: this fct will scale ds!
                    Vs = dParcel*dParcel*dParcel*piBySix;
                }

                force = -g_.value() * RhoInt * Vs;

                forceSubM(0).scaleForce(force,ds,index);
            }

            if(forceSubM(0).verbose() && index==0)
            {
                Pout << "g = " << g_.value() << endl;
                Pout << "cellI = " << cellI << endl;
                Pout << "index = " << index << endl;
                if (forceSubM(0).interpolation())
                    Pout << "Rho(" << position <<") = "<< RhoInt << endl;
                else
                    Pout << "forceSubM(0).rhoField()[cellI] = " << forceSubM(0).rhoField()[cellI] << endl;
                Pout << "Vs = " << Vs << endl;
                Pout << "force = " << force << endl;

                if (forceSubM(0).ms())
                {
                    label ind = particleCloud_.clumpIndexOfParticle(index);
                    Pout << "clumpType = " << particleCloud_.clumpType(ind) << endl;
                }
            }

            //Set value fields and write the probe
            if(probeIt_)
            {
                #include "setupProbeModelfields.H"
                // Note: for other than ext one could use vValues.append(x)
                // instead of setSize
                vValues.setSize(vValues.size()+1, force);           //first entry must the be the force
                sValues.setSize(sValues.size()+1, particleCloud_.cfdemCloud::particleVolume(index)); //change this to Vs after TH is clean
                particleCloud_.probeM().writeProbe(index, sValues, vValues);
            }
        }

        // write particle based data to global array
        forceSubM(0).partToArray(index,force,vector::zero);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
