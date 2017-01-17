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

#include "SchillerNaumannDrag.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(SchillerNaumannDrag, 0);

addToRunTimeSelectionTable
(
    forceModel,
    SchillerNaumannDrag,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
SchillerNaumannDrag::SchillerNaumannDrag
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    velFieldName_(propsDict_.lookup("velFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    UsFieldName_(propsDict_.lookup("granVelFieldName")),
    UsField_(sm.mesh().lookupObject<volVectorField> (UsFieldName_))
{
    // suppress particle probe
    if (probeIt_ && propsDict_.found("suppressProbe"))
        probeIt_=!Switch(propsDict_.lookup("suppressProbe"));
    if(probeIt_)
    {
        particleCloud_.probeM().initialize(typeName, typeName+".logDat");
        particleCloud_.probeM().vectorFields_.append("dragForce"); //first entry must the be the force
        particleCloud_.probeM().vectorFields_.append("Urel");      //other are debug
        particleCloud_.probeM().scalarFields_.append("Rep");       //other are debug
        particleCloud_.probeM().scalarFields_.append("Cd");        //other are debug
        particleCloud_.probeM().writeHeader();
    }

    // init force sub model
    setForceSubModels(propsDict_);

    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(0,true); // activate treatExplicit switch
    forceSubM(0).setSwitchesList(2,true); // activate implDEM switch
    forceSubM(0).setSwitchesList(3,true); // activate search for verbose switch
    forceSubM(0).setSwitchesList(4,true); // activate search for interpolate switch
    forceSubM(0).setSwitchesList(8,true); // activate scalarViscosity switch

    // read those switches defined above, if provided in dict
    for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
        forceSubM(iFSub).readSwitches();

    particleCloud_.checkCG(false);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

SchillerNaumannDrag::~SchillerNaumannDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void SchillerNaumannDrag::setForce() const
{
    const volScalarField& nufField = forceSubM(0).nuField();
    const volScalarField& rhoField = forceSubM(0).rhoField();

    //update force submodels to prepare for loop
    for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
        forceSubM(iFSub).preParticleLoop(forceSubM(iFSub).verbose());

    vector position(0,0,0);
    scalar voidfraction(1);
    vector Ufluid(0,0,0);
    vector drag(0,0,0);
    vector dragExplicit(0,0,0);
  	scalar dragCoefficient(0);
    label cellI=0;
    vector Us(0,0,0);
    vector Ur(0,0,0);
    scalar ds(0);
    scalar dParcel(0);
    scalar nuf(0);
    scalar rho(0);
    scalar magUr(0);
    scalar Rep(0);
    scalar Cd(0);

    #include "resetVoidfractionInterpolator.H"
    #include "resetUInterpolator.H"
    #include "setupProbeModel.H"

    for(int index = 0;index <  particleCloud_.numberOfParticles(); index++)
    {
        drag = vector(0,0,0);
        dragExplicit = vector(0,0,0);
        Ufluid = vector(0,0,0);
        cellI = particleCloud_.cellIDs()[index][0];

        if (cellI > -1) // particle Found
        {
            if(forceSubM(0).interpolation())
            {
                position = particleCloud_.position(index);
                voidfraction = voidfractionInterpolator_().interpolate(position,cellI);
                Ufluid = UInterpolator_().interpolate(position,cellI);
            }else
            {
                voidfraction = voidfraction_[cellI];
                Ufluid = U_[cellI];
            }

            Us = particleCloud_.velocity(index);
            ds = 2*particleCloud_.radius(index);
            dParcel = ds;
            forceSubM(0).scaleDia(ds,index); //caution: this fct will scale ds!
            nuf = nufField[cellI];
            rho = rhoField[cellI];

            //Update any scalar or vector quantity
            for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
                  forceSubM(iFSub).update(  index, 
                                            cellI,
                                            ds,
                                            Ufluid, 
                                            Us, 
                                            nuf,
                                            rho,
                                            forceSubM(0).verbose()
                                         );

            Ur = Ufluid-Us;
            magUr = mag(Ur);
            Rep = 0;
            Cd = 0;

            if (magUr > 0)
            {
               // calc particle Re Nr
                Rep = ds*magUr/nuf;

                // calc fluid drag Coeff
                Cd = max(0.44,24.0/Rep*(1.0+0.15*pow(Rep,0.687)));

                // calc particle's drag
                dragCoefficient = 0.125*Cd*rho
                                  *M_PI
                                  *ds*ds
                                  *magUr;

                if (modelType_=="B")
                    dragCoefficient /= voidfraction;

                // calc particle's drag
                forceSubM(0).scaleCoeff(dragCoefficient,dParcel,index);
                drag = dragCoefficient*Ur;

                // explicitCorr
                for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
                    forceSubM(iFSub).explicitCorr( drag, 
                                                   dragExplicit,
                                                   dragCoefficient,
                                                   Ufluid, U_[cellI], Us, UsField_[cellI],
                                                   forceSubM(iFSub).verbose()
                                                 );
            }

            if(forceSubM(0).verbose() && index >-1 && index <102)
            {
                Pout << "index = " << index << endl;
                Pout << "Us = " << Us << endl;
                Pout << "Ur = " << Ur << endl;
                Pout << "dprim = " << ds << endl;
                Pout << "rho = " << rho << endl;
                Pout << "nuf = " << nuf << endl;
                Pout << "voidfraction = " << voidfraction << endl;
                Pout << "Rep = " << Rep << endl;
                Pout << "Cd = " << Cd << endl;
                Pout << "drag = " << drag << endl;
            }

            //Set value fields and write the probe
            if(probeIt_)
            {
                #include "setupProbeModelfields.H"

                // Note: for other than ext one could use vValues.append(x)
                // instead of setSize
                vValues.setSize(vValues.size()+1, drag);           //first entry must the be the force
                vValues.setSize(vValues.size()+1, Ur);
                sValues.setSize(sValues.size()+1, Rep);
                sValues.setSize(sValues.size()+1, Cd);
                particleCloud_.probeM().writeProbe(index, sValues, vValues);
            }
        }

        // write particle based data to global array
        forceSubM(0).partToArray(index,drag,dragExplicit,Ufluid,dragCoefficient);
    }

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
