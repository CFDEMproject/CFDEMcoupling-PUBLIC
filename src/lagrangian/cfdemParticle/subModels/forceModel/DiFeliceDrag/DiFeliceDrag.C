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

#include "DiFeliceDrag.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(DiFeliceDrag, 0);

addToRunTimeSelectionTable
(
    forceModel,
    DiFeliceDrag,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
DiFeliceDrag::DiFeliceDrag
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    velFieldName_(propsDict_.lookupOrDefault<word>("velFieldName","U")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    voidfractionFieldName_(propsDict_.lookupOrDefault<word>("voidfractionFieldName","voidfraction")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    UsFieldName_(propsDict_.lookupOrDefault<word>("granVelFieldName","Us")),
    UsField_(sm.mesh().lookupObject<volVectorField> (UsFieldName_))
{
    // suppress particle probe
    if (probeIt_ && propsDict_.found("suppressProbe"))
        probeIt_=!Switch(propsDict_.lookup("suppressProbe"));
    if (probeIt_)
    {
        particleCloud_.probeM().initialize(typeName, typeName+".logDat");
        particleCloud_.probeM().vectorFields_.append("dragForce"); //first entry must be the force
        particleCloud_.probeM().vectorFields_.append("Urel");      //other are debug
        particleCloud_.probeM().scalarFields_.append("Rep");       //other are debug
        particleCloud_.probeM().scalarFields_.append("Cd");        //other are debug
        particleCloud_.probeM().scalarFields_.append("voidfraction");       //other are debug
        particleCloud_.probeM().writeHeader();
    }

    particleCloud_.checkCG(true);

    // init force sub model
    setForceSubModels(propsDict_);

    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(0,true); // activate treatExplicit switch
    forceSubM(0).setSwitchesList(2,true,true); // activate implForceDEM switch and set default to true
    forceSubM(0).setSwitchesList(3,true); // activate search for verbose switch
    forceSubM(0).setSwitchesList(4,true); // activate search for interpolate switch
    forceSubM(0).setSwitchesList(8,true); // activate scalarViscosity switch
    forceSubM(0).setSwitchesList(17,true); // activate DiFelice drag correction switch
    forceSubM(0).setSwitchesList(18,true); // activate Rong drag correction switch
    forceSubM(0).setSwitchesList(23,true); // activate Tang drag correction switch
    forceSubM(0).setSwitchesList(28,true); // allow particle specific coarsegraining

    //set default switches (hard-coded default = false)
    forceSubM(0).setSwitches(17,true);  // DiFelice drag correction as default
    forceSubM(0).setSwitches(18,false); // Rong drag correction NOT default
    forceSubM(0).setSwitches(23,false); // Tang drag correction NOT default

    // read those switches defined above, if provided in dict
    for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
        forceSubM(iFSub).readSwitches();

    // sanity check
    if (forceSubM(0).voidageFunctionRong()) // not default
    {
        Info << "Using non-default model voidageFunctionRong!" << endl;
        forceSubM(0).setSwitches(17,false);
        forceSubM(0).setSwitches(23,false);
    }
    else if(forceSubM(0).voidageFunctionTang()) // not default
    {
        Info << "Using non-default model voidageFunctionTang!" << endl;
        forceSubM(0).setSwitches(17,false);
        forceSubM(0).setSwitches(18,false);
    }

    // setup required communication
    for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
        forceSubM(iFSub).setupCommunication();
}

void DiFeliceDrag::MSinit()
{
    // NOTE: all MS related operations need to be performed during init of MS cloud
    if (particleCloud_.shapeTypeName() == "multisphere" && !particleBased_)
    {
        Info << type() << ": activating multisphere mode..." << endl;
        forceSubM(0).setSwitches(11,true); // this is a MS model
        particleCloud_.checkCG(false);

        // read hydraulic diameter correctors
        readDHcorr(propsDict_);

        // re-setup required communication for MS
        for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
            forceSubM(iFSub).setupCommunication();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

DiFeliceDrag::~DiFeliceDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void DiFeliceDrag::setForce() const
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

    // loop all objects
    for (int index = 0;index <  particleCloud_.numberOfObjects(particleBased_); index++)
    {
        cellI = particleCloud_.cellIDs(particleBased_)[index][0];

        drag=vector::zero;
        dragExplicit=vector::zero;
        Ufluid=vector::zero;
        dragCoefficient = 0;

        if (cellI > -1) // particle Found
        {
            if (forceSubM(0).interpolation())
            {
                position = particleCloud_.position(index, particleBased_);

                voidfraction = voidfractionInterpolator_().interpolate(position,cellI);
                Ufluid = UInterpolator_().interpolate(position,cellI);

                //Ensure interpolated void fraction to be meaningful
                // Info << " --> voidfraction: " << voidfraction << endl;
                if(voidfraction>1.00) voidfraction = 1.00;
                if(voidfraction<0.30) voidfraction = 0.30;
            }
            else
            {
                voidfraction = voidfraction_[cellI];
                Ufluid = U_[cellI];
            }

            // correct voidfraction
            ds = particleCloud_.diameter(index, particleBased_);
//             for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
//                 forceSubM(iFSub).calculateCorrectedVoidage(
//                     voidfraction, Ufluid, U_.mesh().V()[cellI], ds );

            Us = particleCloud_.velocity(index, particleBased_);

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

            if (magUr > SMALL && ds > SMALL)
            {

                // calc particle Re Nr
                Rep = ds*voidfraction*magUr/(nuf+SMALL);

                // calc fluid drag Coeff
                Cd = sqr(0.63 + 4.8/sqrt(Rep));

                // calc model coefficient Xi
                scalar Xi = 0;
                for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
                {
                    forceSubM(iFSub).calcXi(ds,voidfraction,magUr,nuf,Xi);
                }

                // calc particle's drag
                dragCoefficient = 0.125*Cd*rho
                                 *M_PI
                                 *ds*ds
                                 *pow(voidfraction,(2-Xi))*magUr;

                if (modelType_=="B")
                    dragCoefficient /= voidfraction;

                for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
                    forceSubM(iFSub).scaleCoeff(dragCoefficient,dParcel,Rep,index);
                drag = dragCoefficient*Ur; //total drag force!

                // explicitCorr
                for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
                    forceSubM(iFSub).explicitCorr( drag,
                                                   dragExplicit,
                                                   dragCoefficient,
                                                   Ufluid, U_[cellI], Us, UsField_[cellI],
                                                   forceSubM(iFSub).verbose()
                                                 );
            }

            if (forceSubM(0).verbose() && index==0)
            {
                Pout << "cellI = " << cellI << endl;
                Pout << "index = " << index << endl;
                Pout << "Us = " << Us << endl;
                Pout << "Ur = " << Ur << endl;
                Pout << "dprim = " << ds << endl;
                Pout << "dParcel = " << dParcel << endl;
                Pout << "rho = " << rho << endl;
                Pout << "nuf = " << nuf << endl;
                Pout << "voidfraction = " << voidfraction << endl;
                Pout << "Rep = " << Rep << endl;
                Pout << "Cd = " << Cd << endl;
                Pout << "dragCoefficient = " << dragCoefficient << endl;
                Pout << "drag (total) = " << drag << endl;

                if (forceSubM(0).ms())
                {
                    label ind = particleCloud_.clumpIndexOfParticle(index);
                    Pout << "ds/scale = " << ds/forceSubM(0).getCG() << endl;
                    Pout << "scale = " << forceSubM(0).getCG() << endl;
                    Pout << "scaleDrag_ = " << forceSubM(0).scaleDrag() << endl;
                    Pout << "clumpType = "
                         << particleCloud_.clumpType(ind) << endl;
                }
            }

            //Set value fields and write the probe
            if (probeIt_)
            {
                #include "setupProbeModelfields.H"
                // Note: for other than ext one could use vValues.append(x)
                // instead of setSize
                vValues.setSize(vValues.size()+1, drag);           //first entry must the be the force
                vValues.setSize(vValues.size()+1, Ur);
                sValues.setSize(sValues.size()+1, Rep);
                sValues.setSize(sValues.size()+1, Cd);
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
