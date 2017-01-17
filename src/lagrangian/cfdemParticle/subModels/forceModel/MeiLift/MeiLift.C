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

#include "MeiLift.H"
#include "addToRunTimeSelectionTable.H"

//#include "mpi.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(MeiLift, 0);

addToRunTimeSelectionTable
(
    forceModel,
    MeiLift,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
MeiLift::MeiLift
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    velFieldName_(propsDict_.lookup("velFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    useSecondOrderTerms_(false)
{
    if (propsDict_.found("useSecondOrderTerms")) useSecondOrderTerms_=true;

    // init force sub model
    setForceSubModels(propsDict_);
    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(0,true); // activate treatExplicit switch
    forceSubM(0).setSwitchesList(3,true); // activate search for verbose switch
    forceSubM(0).setSwitchesList(4,true); // activate search for interpolate switch
    forceSubM(0).setSwitchesList(8,true); // activate scalarViscosity switch

    //set default switches (hard-coded default = false)
    forceSubM(0).setSwitches(0,true);  // enable treatExplicit, otherwise this force would be implicit in slip vel! - IMPORTANT!

    for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
        forceSubM(iFSub).readSwitches();

    particleCloud_.checkCG(false);

    //Append the field names to be probed
    particleCloud_.probeM().initialize(typeName, typeName+".logDat");
    particleCloud_.probeM().vectorFields_.append("liftForce"); //first entry must the be the force
    particleCloud_.probeM().vectorFields_.append("Urel");        //other are debug
    particleCloud_.probeM().vectorFields_.append("vorticity");  //other are debug
    particleCloud_.probeM().scalarFields_.append("Rep");          //other are debug
    particleCloud_.probeM().scalarFields_.append("Rew");          //other are debug
    particleCloud_.probeM().scalarFields_.append("J_star");       //other are debug
    particleCloud_.probeM().writeHeader();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

MeiLift::~MeiLift()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void MeiLift::setForce() const
{
    const volScalarField& nufField = forceSubM(0).nuField();
    const volScalarField& rhoField = forceSubM(0).rhoField();

    vector position(0,0,0);
    vector lift(0,0,0);
    vector Us(0,0,0);
    vector Ur(0,0,0);
    scalar magUr(0);
    scalar magVorticity(0);
    scalar ds(0);
    scalar dParcel(0);
    scalar nuf(0);
    scalar rho(0);
    scalar voidfraction(1);
    scalar Rep(0);
    scalar Rew(0);
    scalar Cl(0);
    scalar Cl_star(0);
    scalar J_star(0);
    scalar Omega_eq(0);
    scalar alphaStar(0);
    scalar epsilon(0);
    scalar omega_star(0);
    vector vorticity(0,0,0);
    volVectorField vorticity_ = fvc::curl(U_);

    #include "resetVorticityInterpolator.H"
    #include "resetUInterpolator.H"

    #include "setupProbeModel.H"

    for(int index = 0;index <  particleCloud_.numberOfParticles(); index++)
    {
        //if(mask[index][0])
        //{
            lift           = vector::zero;
            label cellI = particleCloud_.cellIDs()[index][0];

            if (cellI > -1) // particle Found
            {
                Us = particleCloud_.velocity(index);

                if( forceSubM(0).interpolation() )
                {
	                position       = particleCloud_.position(index);
                    Ur               = UInterpolator_().interpolate(position,cellI) 
                                        - Us;
                    vorticity       = vorticityInterpolator_().interpolate(position,cellI);
                }
                else
                {
                    Ur =  U_[cellI]
                          - Us;
                    vorticity=vorticity_[cellI];
                }

                magUr           = mag(Ur);
                magVorticity = mag(vorticity);

                if (magUr > 0 && magVorticity > 0)
                {
                    ds  = 2*particleCloud_.radius(index);
                    dParcel = ds;
                    forceSubM(0).scaleDia(ds,index); //caution: this fct will scale ds!
                    nuf = nufField[cellI];
                    rho = rhoField[cellI];

                    //Update any scalar or vector quantity
                    for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
                          forceSubM(iFSub).update(  index, 
                                                    cellI,
                                                    ds,
                                                    nuf,
                                                    rho,
                                                    forceSubM(0).verbose()
                                                 );


                    // calc dimensionless properties
                    Rep = ds*magUr/nuf;
		            Rew = magVorticity*ds*ds/nuf;

                    alphaStar   = magVorticity*ds/magUr/2.0;
                    epsilon     = sqrt(2.0*alphaStar /Rep );
                    omega_star=2.0*alphaStar;

                    //Basic model for the correction to the Saffman lift
                    //Based on McLaughlin (1991)
                    if(epsilon < 0.1)
                    {
                        J_star = -140 *epsilon*epsilon*epsilon*epsilon*epsilon 
                                             *log( 1./(epsilon*epsilon+SMALL) );
                    }
                    else if(epsilon > 20)
                    {
                      J_star = 1.0-0.287/(epsilon*epsilon+SMALL);
                    }
                    else
                    {
                     J_star = 0.3
                                *(     1.0
                                      +tanh(  2.5 * log10(epsilon+0.191)  )
                                 )
                                *(    0.667
                                     +tanh(  6.0 * (epsilon-0.32)  )
                                  );
                    }
                    Cl=J_star*4.11*epsilon; //multiply McLaughlin's correction to the basic Saffman model

                    //Second order terms given by Loth and Dorgan 2009 
                    if(useSecondOrderTerms_)
                    {   
                        Omega_eq = omega_star/2.0*(1.0-0.0075*Rew)*(1.0-0.062*sqrt(Rep)-0.001*Rep);
                        Cl_star=1.0-(0.675+0.15*(1.0+tanh(0.28*(omega_star/2.0-2.0))))*tanh(0.18*sqrt(Rep));
                        Cl += Omega_eq*Cl_star;
                    }

                    lift =  0.125*M_PI
                           *rho
                           *Cl  
                           *magUr*Ur^vorticity/magVorticity
                           *ds*ds; //total force on all particles in parcel

                    forceSubM(0).scaleForce(lift,dParcel,index);

                    if (modelType_=="B")
                    {
                        voidfraction = particleCloud_.voidfraction(index);
                        lift /= voidfraction;
                    }
                }

                //**********************************        
                //SAMPLING AND VERBOSE OUTOUT
                if( forceSubM(0).verbose() )
                {   
                    Pout << "index = " << index << endl;
                    Pout << "Us = " << Us << endl;
                    Pout << "Ur = " << Ur << endl;
                    Pout << "vorticity = " << vorticity << endl;
                    Pout << "dprim = " << ds << endl;
                    Pout << "rho = " << rho << endl;
                    Pout << "nuf = " << nuf << endl;
                    Pout << "Rep = " << Rep << endl;
                    Pout << "Rew = " << Rew << endl;
                    Pout << "alphaStar = " << alphaStar << endl;
                    Pout << "epsilon = " << epsilon << endl;
                    Pout << "J_star = " << J_star << endl;
                    Pout << "lift = " << lift << endl;
                }

                //Set value fields and write the probe
                if(probeIt_)
                {
                    #include "setupProbeModelfields.H"
                    // Note: for other than ext one could use vValues.append(x)
                    // instead of setSize
                    vValues.setSize(vValues.size()+1, lift);           //first entry must the be the force
                    vValues.setSize(vValues.size()+1, Ur);
                    vValues.setSize(vValues.size()+1, vorticity); 
                    sValues.setSize(sValues.size()+1, Rep);
                    sValues.setSize(sValues.size()+1, Rew);
                    sValues.setSize(sValues.size()+1, J_star);
                    particleCloud_.probeM().writeProbe(index, sValues, vValues);
                }
                // END OF SAMPLING AND VERBOSE OUTOUT
                //**********************************        

            }
            // write particle based data to global array
            forceSubM(0).partToArray(index,lift,vector::zero);
        //}
    }

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
