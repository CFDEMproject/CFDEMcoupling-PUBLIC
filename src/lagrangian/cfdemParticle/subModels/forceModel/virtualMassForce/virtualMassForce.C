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

#include "virtualMassForce.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

#define NOTONCPU 9999

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(virtualMassForce, 0);

addToRunTimeSelectionTable
(
    forceModel,
    virtualMassForce,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
virtualMassForce::virtualMassForce
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    velFieldName_(propsDict_.lookup("velFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    phiFieldName_(propsDict_.lookup("phiFieldName")),
    phi_(sm.mesh().lookupObject<surfaceScalarField> (phiFieldName_)),
    UrelOld_(NULL),
    splitUrelCalculation_(false),
    Cadd_(0.5)
{

    if (particleCloud_.dataExchangeM().maxNumberOfParticles() > 0)
    {
        // get memory for 2d array
        particleCloud_.dataExchangeM().allocateArray(UrelOld_,NOTONCPU,3);
    }

    // init force sub model
    setForceSubModels(propsDict_);
    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(0,true); // activate treatExplicit switch
    forceSubM(0).setSwitchesList(4,true); // activate search for interpolate switch
    forceSubM(0).readSwitches();

    //Extra switches/settings
    if(propsDict_.found("splitUrelCalculation"))
    {
        splitUrelCalculation_ = readBool(propsDict_.lookup("splitUrelCalculation")); 
        if(splitUrelCalculation_)
            Info << "Virtual mass model: will split the Urel calculation\n";
            Info << "WARNING: be sure that LIGGGHTS integration takes ddtv_p implicitly into account! \n";
    }
    if(propsDict_.found("Cadd"))
    {
        Cadd_ = readScalar(propsDict_.lookup("Cadd")); 
        Info << "Virtual mass model: using non-standard Cadd = " << Cadd_ << endl;
    }

    particleCloud_.checkCG(true);

    //Append the field names to be probed
    particleCloud_.probeM().initialize(typeName, "virtualMass.logDat");
    particleCloud_.probeM().vectorFields_.append("virtualMassForce"); //first entry must the be the force
    particleCloud_.probeM().vectorFields_.append("Urel");
    particleCloud_.probeM().vectorFields_.append("UrelOld");
    particleCloud_.probeM().vectorFields_.append("ddtUrel");
    particleCloud_.probeM().scalarFields_.append("Vs");
    particleCloud_.probeM().scalarFields_.append("rho");
    particleCloud_.probeM().writeHeader();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

virtualMassForce::~virtualMassForce()
{
    delete UrelOld_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void virtualMassForce::setForce() const
{
    reAllocArrays();

    scalar dt = U_.mesh().time().deltaT().value();

    vector position(0,0,0);
    vector Ufluid(0,0,0);
    vector Ur(0,0,0);
    vector DDtU(0,0,0);

    //Compute extra vfields in case it is needed
    volVectorField DDtU_(0.0*U_/U_.mesh().time().deltaT());
    if(splitUrelCalculation_)
        DDtU_ = fvc::ddt(U_) + fvc::div(phi_, U_); //Total Derivative of fluid velocity

    interpolationCellPoint<vector> UInterpolator_(   U_);
    interpolationCellPoint<vector> DDtUInterpolator_(DDtU_);

    #include "setupProbeModel.H"

    bool haveUrelOld_(false); 

    for(int index = 0;index <  particleCloud_.numberOfParticles(); index++)
    {
            vector virtualMassForce(0,0,0);
            label cellI = particleCloud_.cellIDs()[index][0];

            if (cellI > -1) // particle Found
            {

                if(forceSubM(0).interpolation()) 
                {
	                position     = particleCloud_.position(index);
                    Ufluid       = UInterpolator_.interpolate(position,cellI);
                }
                else
                {
                    Ufluid = U_[cellI];
                }
           

                if(splitUrelCalculation_)  //if split, just use total derivative of fluid velocity
                if(forceSubM(0).interpolation()) 
                {
                    DDtU = DDtUInterpolator_.interpolate(position,cellI);
                }
                else
                {
                    DDtU = DDtU_[cellI];
                }
                else
                {
                    vector Us = particleCloud_.velocity(index);
                    Ur = Ufluid - Us;
                }

                
                //Check of particle was on this CPU the last step
                if(UrelOld_[index][0]==NOTONCPU) //use 1. element to indicate that particle was on this CPU the last time step
                    haveUrelOld_ = false;
                else
                    haveUrelOld_ = true;

                vector UrelOld(0.,0.,0.);
                vector ddtUrel(0.,0.,0.);
                for(int j=0;j<3;j++)
                {
                    UrelOld[j]         = UrelOld_[index][j];
                    UrelOld_[index][j] = Ur[j];
                }
                if(haveUrelOld_ ) //only compute force if we have old (relative) velocity
                    ddtUrel = (Ur-UrelOld)/dt;

                if(splitUrelCalculation_) //we can always compute the total derivative in case we split
                    ddtUrel = DDtU;

                scalar ds = 2*particleCloud_.radius(index);
                scalar Vs = ds*ds*ds*M_PI/6;
                scalar rho = forceSubM(0).rhoField()[cellI];

                virtualMassForce = Cadd_ * rho * Vs * ddtUrel;

                if( forceSubM(0).verbose() ) //&& index>100 && index < 105)
                {
                    Pout << "index / cellI = " << index << tab << cellI << endl;
                    Pout << "position = " << particleCloud_.position(index) << endl;
                }

                //Set value fields and write the probe
                if(probeIt_)
                {
                    #include "setupProbeModelfields.H"
                    vValues.append(virtualMassForce);           //first entry must the be the force
                    vValues.append(Ur);
                    vValues.append(UrelOld);
                    vValues.append(ddtUrel);
                    sValues.append(Vs);
                    sValues.append(rho);
                    particleCloud_.probeM().writeProbe(index, sValues, vValues);
                }
            }
            else    //particle not on this CPU
                UrelOld_[index][0]=NOTONCPU;
    
            // write particle based data to global array
            forceSubM(0).partToArray(index,virtualMassForce,vector::zero);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void Foam::virtualMassForce::reAllocArrays() const
{
    if(particleCloud_.numberOfParticlesChanged())
    {
        Pout << "virtualMassForce::reAllocArrays..." << endl;
        particleCloud_.dataExchangeM().allocateArray(UrelOld_,NOTONCPU,3);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
