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
    phi_(readScalar(propsDict_.lookup("phi")))
{
    //Append the field names to be probed
    particleCloud_.probeM().initialize(typeName, "gidaspowDrag.logDat");
    particleCloud_.probeM().vectorFields_.append("dragForce"); //first entry must the be the force
    particleCloud_.probeM().vectorFields_.append("Urel");        //other are debug
    particleCloud_.probeM().scalarFields_.append("KslLag");                 //other are debug
    particleCloud_.probeM().scalarFields_.append("voidfraction");       //other are debug
    particleCloud_.probeM().writeHeader();

    if (propsDict_.found("verbose")) verbose_=true;
    if (propsDict_.found("treatExplicit")) treatExplicit_=true;
    particleCloud_.checkCG(false);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

GidaspowDrag::~GidaspowDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void GidaspowDrag::setForce() const
{
    // get viscosity field
    #ifdef comp
        const volScalarField nufField = particleCloud_.turbulence().mu() / rho_;
    #else
        const volScalarField& nufField = particleCloud_.turbulence().nu();
    #endif

    #include "setupProbeModel.H"

    for(int index = 0;index <  particleCloud_.numberOfParticles(); ++index)
    {
        //if(mask[index][0])
        //{
            vector drag(0,0,0);
            label cellI = particleCloud_.cellIDs()[index][0];

            if (cellI > -1) // particle Found
            {
                vector Us = particleCloud_.velocity(index);
                vector Ur = U_[cellI]-Us;
                scalar magUr = mag(Ur);
                scalar ds = 2*particleCloud_.radius(index);
                scalar voidfraction = particleCloud_.voidfraction(index);
                scalar rho = rho_[cellI];
                scalar nuf = nufField[cellI];
                scalar CdMagUrLag=0;//Cd of the very particle
                scalar KslLag;      //momentum exchange of the very particle

                if(voidfraction > 0.8) //dilute
                {
                    CdMagUrLag = (24.0*nuf/(ds*voidfraction))
                                 *(scalar(1)+0.15*Foam::pow(ds*voidfraction*magUr/nuf, 0.687));
                    KslLag = 0.75*((1-voidfraction)*rho*voidfraction*CdMagUrLag
                                /(ds*Foam::pow(voidfraction,2.65)));
                }
                else  //dense
                {
                    KslLag = (150*Foam::pow(1-voidfraction,2)*nuf*rho)/
                             (voidfraction*ds*ds+SMALL)
                            +
                             (1.75*(1-voidfraction) * magUr * rho)/
                             ((ds));
                }

                //divide by number of particles per unit volume - Enwald (Int J Multiphase Flow, 22, 21-61, pp39
                KslLag /= (particleCloud_.averagingM().UsWeightField()[cellI]/particleCloud_.mesh().V()[cellI]);

                drag = KslLag*Ur*phi_;

                if (modelType_=="B")
                    drag /= voidfraction;


                //Set value fields and write the probe
                if(probeIt_)
                {
                    #include "setupProbeModelfields.H"
                    vValues.append(drag);   //first entry must the be the force
                    vValues.append(Ur);
                    sValues.append(KslLag);
                    sValues.append(voidfraction);
                    particleCloud_.probeM().writeProbe(index, sValues, vValues);
                }
            }

            // set force on particle
            if(treatExplicit_) for(int j=0;j<3;j++) expForces()[index][j] += drag[j];
            else  for(int j=0;j<3;j++) impForces()[index][j] += drag[j];
            for(int j=0;j<3;j++) DEMForces()[index][j] += drag[j];
        //}// end if mask
    }// end loop particles
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
