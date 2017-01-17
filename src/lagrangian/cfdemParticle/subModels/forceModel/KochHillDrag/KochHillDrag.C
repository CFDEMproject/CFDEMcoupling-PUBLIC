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

#include "KochHillDrag.H"
#include "addToRunTimeSelectionTable.H"
#include "dataExchangeModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(KochHillDrag, 0);

addToRunTimeSelectionTable
(
    forceModel,
    KochHillDrag,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
KochHillDrag::KochHillDrag
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
    UsFieldName_(propsDict_.lookupOrDefault("granVelFieldName",word("Us"))),
    UsField_(sm.mesh().lookupObject<volVectorField> (UsFieldName_))
{
    // suppress particle probe
    if (probeIt_ && propsDict_.found("suppressProbe"))
        probeIt_=!Switch(propsDict_.lookup("suppressProbe"));
    if(probeIt_)
    {
        particleCloud_.probeM().initialize(typeName, typeName+".logDat");
        particleCloud_.probeM().vectorFields_.append("dragForce"); //first entry must the be the force
        particleCloud_.probeM().vectorFields_.append("Urel");        //other are debug
        particleCloud_.probeM().scalarFields_.append("Rep");          //other are debug
        particleCloud_.probeM().scalarFields_.append("beta");                 //other are debug
        particleCloud_.probeM().scalarFields_.append("voidfraction");       //other are debug
        particleCloud_.probeM().writeHeader();
    }

    // init force sub model
    setForceSubModels(propsDict_);

    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(0,true); // activate search for treatExplicit switch
    forceSubM(0).setSwitchesList(2,true); // activate search for implDEM switch
    forceSubM(0).setSwitchesList(3,true); // activate search for verbose switch
    forceSubM(0).setSwitchesList(4,true); // activate search for interpolate switch
    forceSubM(0).setSwitchesList(7,true); // activate implForceDEMacc switch
    forceSubM(0).setSwitchesList(8,true); // activate scalarViscosity switch

    // read those switches defined above, if provided in dict
    for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
        forceSubM(iFSub).readSwitches();

    particleCloud_.checkCG(true);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

KochHillDrag::~KochHillDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void KochHillDrag::setForce() const
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
	scalar Vs(0);
	scalar volumefraction(0);
    scalar betaP(0);

    scalar piBySix(M_PI/6);


    int couplingInterval(particleCloud_.dataExchangeM().couplingInterval());

    #include "resetVoidfractionInterpolator.H"
    #include "resetUInterpolator.H"
    #include "setupProbeModel.H"

    for(int index = 0;index <  particleCloud_.numberOfParticles(); index++)
    {
            cellI = particleCloud_.cellIDs()[index][0];
            drag = vector(0,0,0);
            dragExplicit = vector(0,0,0);
            dragCoefficient=0;
            betaP = 0;
            Vs = 0;
            Ufluid =vector(0,0,0);
            voidfraction=0;

            if (cellI > -1) // particle Found
            {
                if(forceSubM(0).interpolation())
                {
	                position = particleCloud_.position(index);
                    voidfraction = voidfractionInterpolator_().interpolate(position,cellI);
                    Ufluid = UInterpolator_().interpolate(position,cellI);

                    //Ensure interpolated void fraction to be meaningful
                    // Info << " --> voidfraction: " << voidfraction << endl;
                    if(voidfraction>1.00) voidfraction = 1.00;
                    if(voidfraction<0.40) voidfraction = 0.40;
                }else
                {
					voidfraction = voidfraction_[cellI];
                    Ufluid = U_[cellI];
                }

                ds = particleCloud_.d(index);
                dParcel = ds;
                forceSubM(0).scaleDia(ds,index); //caution: this fct will scale ds!
                nuf = nufField[cellI];
                rho = rhoField[cellI];

                Us = particleCloud_.velocity(index);

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

                Vs = ds*ds*ds*piBySix;

                volumefraction = max(SMALL,min(1-SMALL,1-voidfraction));

                if (magUr > 0)
                {
                    // calc particle Re Nr
                    Rep = ds*voidfraction*magUr/(nuf+SMALL);

                    // calc model coefficient F0
                    scalar F0=0.;
                    if(volumefraction < 0.4)
                    {
                        F0 = (1. + 3.*sqrt((volumefraction)/2.) + (135./64.)*volumefraction*log(volumefraction)
                              + 16.14*volumefraction
                             )/
                             (1+0.681*volumefraction-8.48*sqr(volumefraction)
                              +8.16*volumefraction*volumefraction*volumefraction
                             );
                    } else {
                        F0 = 10*volumefraction/(voidfraction*voidfraction*voidfraction);
                    }

                    // calc model coefficient F3
                    scalar F3 = 0.0673+0.212*volumefraction+0.0232/pow(voidfraction,5);

                    //Calculate F (the factor 0.5 is introduced, since Koch and Hill, ARFM 33:619–47, use the radius
                    //to define Rep, and we use the particle diameter, see vanBuijtenen et al., CES 66:2368–2376.
                    scalar F = voidfraction * (F0 + 0.5*F3*Rep);

                    // calc drag model coefficient betaP
                    betaP = 18.*nuf*rho/(ds*ds)*voidfraction*F;

                    // calc particle's drag
                    dragCoefficient = Vs*betaP;
                    if (modelType_=="B")
                        dragCoefficient /= voidfraction;

                    forceSubM(0).scaleCoeff(dragCoefficient,dParcel,index);

                    if(forceSubM(0).switches()[7]) // implForceDEMaccumulated=true
                    {
		                //get drag from the particle itself
		                for (int j=0 ; j<3 ; j++) drag[j] = particleCloud_.fAccs()[index][j]/couplingInterval;
                    }else
                    {
                        drag = dragCoefficient * Ur;

                        // explicitCorr
                        for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
                            forceSubM(iFSub).explicitCorr( drag, 
                                                           dragExplicit,
                                                           dragCoefficient,
                                                           Ufluid, U_[cellI], Us, UsField_[cellI],
                                                           forceSubM(iFSub).verbose()
                                                         );
                    }
                }

                if(forceSubM(0).verbose() && index >=0 && index <2)
                {
                    Pout << "cellI = " << cellI << endl;
                    Pout << "index = " << index << endl;
                    Pout << "Us = " << Us << endl;
                    Pout << "Ur = " << Ur << endl;
                    Pout << "dprim = " << ds << endl;
                    Pout << "rho = " << rho << endl;
                    Pout << "nuf = " << nuf << endl;
                    Pout << "voidfraction = " << voidfraction << endl;
                    Pout << "Rep = " << Rep << endl;
                    Pout << "betaP = " << betaP << endl;
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
                    sValues.setSize(sValues.size()+1, betaP);
                    sValues.setSize(sValues.size()+1, voidfraction);
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
