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

#include "GidaspowDrag.H"
#include "addToRunTimeSelectionTable.H"
#include "averagingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(GidaspowDrag, 0);

addToRunTimeSelectionTable
(
    forceModel,
    GidaspowDrag,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
GidaspowDrag::GidaspowDrag
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    verbose_(false),
    velFieldName_(propsDict_.lookup("velFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    densityFieldName_(propsDict_.lookup("densityFieldName")),
    rho_(sm.mesh().lookupObject<volScalarField> (densityFieldName_)),
    voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    phi_(readScalar(propsDict_.lookup("phi"))),
    interpolation_(false),
    scaleDia_(1.),
    scaleDrag_(1.)
{
    //Append the field names to be probed
    particleCloud_.probeM().initialize(typeName, "gidaspowDrag.logDat");
    particleCloud_.probeM().vectorFields_.append("dragForce"); //first entry must  be the force
    particleCloud_.probeM().vectorFields_.append("Urel");
    particleCloud_.probeM().scalarFields_.append("Rep");
    particleCloud_.probeM().scalarFields_.append("betaP");
    particleCloud_.probeM().scalarFields_.append("voidfraction");
    particleCloud_.probeM().writeHeader();

    if (propsDict_.found("verbose")) verbose_=true;
    if (propsDict_.found("treatExplicit")) treatExplicit_=true;
    if (propsDict_.found("interpolation")) interpolation_=true;
    if (propsDict_.found("implDEM"))
    {
        treatExplicit_=false;
        implDEM_=true;
        setImpDEMdrag();
        Info << "Using implicit DEM drag formulation." << endl;
    }
    particleCloud_.checkCG(true);
    if (propsDict_.found("scale"))
        scaleDia_=scalar(readScalar(propsDict_.lookup("scale")));
    if (propsDict_.found("scaleDrag"))
        scaleDrag_=scalar(readScalar(propsDict_.lookup("scaleDrag")));

    Info << "Gidaspow - interpolation switch: " << interpolation_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

GidaspowDrag::~GidaspowDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void GidaspowDrag::setForce() const
{
    if (scaleDia_ > 1)
        Info << "Gidaspow using scale = " << scaleDia_ << endl;
    else if (particleCloud_.cg() > 1){
        scaleDia_=particleCloud_.cg();
        Info << "Gidaspow using scale from liggghts cg = " << scaleDia_ << endl;
    }

    // get viscosity field
    #ifdef comp
        const volScalarField nufField = particleCloud_.turbulence().mu() / rho_;
    #else
        const volScalarField& nufField = particleCloud_.turbulence().nu();
    #endif

    vector position(0,0,0);
    scalar voidfraction(1);
    vector Ufluid(0,0,0);
    vector drag(0,0,0);
    label cellI=0;

    vector Us(0,0,0);
    vector Ur(0,0,0);
    scalar ds(0);
    scalar nuf(0);
    scalar rho(0);
    scalar magUr(0);
    scalar Rep(0);
    scalar Vs(0);
    scalar localPhiP(0);

    scalar CdMagUrLag(0);       //Cd of the very particle
    scalar KslLag(0);           //momentum exchange of the very particle (per unit volume)
    scalar betaP(0);             //momentum exchange of the very particle

    interpolationCellPoint<scalar> voidfractionInterpolator_(voidfraction_);
    interpolationCellPoint<vector> UInterpolator_(U_);

    #include "setupProbeModel.H"

    for(int index = 0;index <  particleCloud_.numberOfParticles(); ++index)
    {
        //if(mask[index][0])
        //{
            cellI = particleCloud_.cellIDs()[index][0];
            drag = vector(0,0,0);
            betaP = 0;
            Vs = 0;
            Ufluid =vector(0,0,0);
            voidfraction=0;

            if (cellI > -1) // particle Found
            {

                if(interpolation_)
                {
	                position     = particleCloud_.position(index);
                    voidfraction = voidfractionInterpolator_.interpolate(position,cellI);
                    Ufluid       = UInterpolator_.interpolate(position,cellI);
                    //Ensure interpolated void fraction to be meaningful
                    // Info << " --> voidfraction: " << voidfraction << endl;
                    if(voidfraction>1.00) voidfraction = 1.0;
                    if(voidfraction<0.10) voidfraction = 0.10;
                }
                else
                {
					voidfraction = voidfraction_[cellI];
                    Ufluid = U_[cellI];
                }

                Us = particleCloud_.velocity(index);
                Ur = Ufluid-Us;
                magUr = mag(Ur);
                ds = 2*particleCloud_.radius(index)*phi_;
                rho = rho_[cellI];
                nuf = nufField[cellI];

                Rep=0.0;
                localPhiP = 1.0f-voidfraction+SMALL;
                Vs = ds*ds*ds*M_PI/6;

                //Compute specific drag coefficient (i.e., Force per unit slip velocity and per m³ SUSPENSION)
                if(voidfraction > 0.8) //dilute
                {
                    Rep=ds/scaleDia_*voidfraction*magUr/nuf;
                    CdMagUrLag = (24.0*nuf/(ds/scaleDia_*voidfraction)) //1/magUr missing here, but compensated in expression for KslLag!
                                 *(scalar(1)+0.15*Foam::pow(Rep, 0.687));

                    KslLag = 0.75*(
                                            rho*localPhiP*voidfraction*CdMagUrLag
                                          /
                                            (ds/scaleDia_*Foam::pow(voidfraction,2.65))
                                          );
                }
                else  //dense
                {
                    KslLag = (150*Foam::pow(localPhiP,2)*nuf*rho)/
                             (voidfraction*ds/scaleDia_*ds/scaleDia_+SMALL)
                            +
                             (1.75*(localPhiP) * magUr * rho)/
                             ((ds/scaleDia_));
                }

                // calc particle's drag coefficient (i.e., Force per unit slip velocity and per m³ PARTICLE)
                betaP = KslLag / localPhiP;

                // calc particle's drag
                drag = Vs * betaP * Ur * scaleDrag_;

                if (modelType_=="B")
                    drag /= voidfraction;

                if(verbose_ && index >=0 && index <2)
                {
                    Pout << "cellI = " << cellI << endl;
                    Pout << "index = " << index << endl;
                    Pout << "Us = " << Us << endl;
                    Pout << "Ur = " << Ur << endl;
                    Pout << "ds = " << ds << endl;
                    Pout << "ds/scale = " << ds/scaleDia_ << endl;
                    Pout << "phi = " << phi_ << endl;
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
                    vValues.append(drag);   //first entry must the be the force
                    vValues.append(Ur);
                    sValues.append(Rep);
                    sValues.append(betaP);
                    sValues.append(voidfraction);
                    particleCloud_.probeM().writeProbe(index, sValues, vValues);
                }
            }

            // set force on particle
            if(treatExplicit_) for(int j=0;j<3;j++) expForces()[index][j] += drag[j];
            else  for(int j=0;j<3;j++) impForces()[index][j] += drag[j];

            // set Cd
            if(implDEM_)
            {
                for(int j=0;j<3;j++) fluidVel()[index][j]=Ufluid[j];

                if (modelType_=="B" && cellI > -1)
                    Cds()[index][0] = Vs*betaP/voidfraction*scaleDrag_;
                else
                    Cds()[index][0] = Vs*betaP*scaleDrag_;

            }else{
                for(int j=0;j<3;j++) DEMForces()[index][j] += drag[j];
            }

        //}// end if mask
    }// end loop particles
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
