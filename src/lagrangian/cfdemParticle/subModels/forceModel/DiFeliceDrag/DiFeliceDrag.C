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

//#include "mpi.h"

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
    verbose_(false),
    velFieldName_(propsDict_.lookup("velFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    densityFieldName_(propsDict_.lookup("densityFieldName")),
    rho_(sm.mesh().lookupObject<volScalarField> (densityFieldName_)),
    voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    interpolation_(false),
    splitImplicitExplicit_(false),
    UsFieldName_(propsDict_.lookup("granVelFieldName")),
    UsField_(sm.mesh().lookupObject<volVectorField> (UsFieldName_)),
    scaleDia_(1.),
    scaleDrag_(1.)
{
    //Append the field names to be probed
    particleCloud_.probeM().initialize(typeName, "diFeliceDrag.logDat");
    particleCloud_.probeM().vectorFields_.append("dragForce"); //first entry must the be the force
    particleCloud_.probeM().vectorFields_.append("Urel");        //other are debug
    particleCloud_.probeM().scalarFields_.append("Rep");          //other are debug
    particleCloud_.probeM().scalarFields_.append("Cd");                 //other are debug
    particleCloud_.probeM().scalarFields_.append("voidfraction");       //other are debug
    particleCloud_.probeM().writeHeader();

    if (propsDict_.found("verbose")) verbose_=true;
    if (propsDict_.found("treatExplicit")) treatExplicit_=true;
    if (propsDict_.found("interpolation"))
    {
        Info << "using interpolated value of U." << endl;
        interpolation_=true;
    }
    if (propsDict_.found("splitImplicitExplicit"))
    {
        Info << "will split implicit / explicit force contributions." << endl;
        splitImplicitExplicit_ = true;
        if(!interpolation_) 
            Info << "WARNING: will only consider fluctuating particle velocity in implicit / explicit force split!" << endl;
    }
    particleCloud_.checkCG(true);
    if (propsDict_.found("scale"))
        scaleDia_=scalar(readScalar(propsDict_.lookup("scale")));
    if (propsDict_.found("scaleDrag"))
        scaleDrag_=scalar(readScalar(propsDict_.lookup("scaleDrag")));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

DiFeliceDrag::~DiFeliceDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void DiFeliceDrag::setForce() const
{
    if (scaleDia_ > 1)
        Info << "DiFeliceDrag using scale = " << scaleDia_ << endl;
    else if (particleCloud_.cg() > 1){
        scaleDia_=particleCloud_.cg();
        Info << "DiFeliceDrag using scale from liggghts cg = " << scaleDia_ << endl;
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
    scalar Cd(0);

	vector UfluidFluct(0,0,0);
    vector UsFluct(0,0,0);
    vector dragExplicit(0,0,0);
  	scalar dragCoefficient(0);

    interpolationCellPoint<scalar> voidfractionInterpolator_(voidfraction_);
    interpolationCellPoint<vector> UInterpolator_(U_);

    #include "setupProbeModel.H"

    for(int index = 0;index <  particleCloud_.numberOfParticles(); index++)
    {
        //if(mask[index][0])
        //{

            cellI = particleCloud_.cellIDs()[index][0];
            drag = vector(0,0,0);

            if (cellI > -1) // particle Found
            {
                if(interpolation_)
                {
                    position = particleCloud_.position(index);
                    voidfraction = voidfractionInterpolator_.interpolate(position,cellI);
                    Ufluid = UInterpolator_.interpolate(position,cellI);
                }else
                {
                    voidfraction = voidfraction_[cellI];
                    Ufluid = U_[cellI];
                }

                Us = particleCloud_.velocity(index);
                Ur = Ufluid-Us;
                ds = 2*particleCloud_.radius(index);
                nuf = nufField[cellI];
                rho = rho_[cellI];
                magUr = mag(Ur);
                Rep = 0;
                Cd = 0;
                dragCoefficient = 0;

                if (magUr > 0)
                {

                    // calc particle Re Nr
                    Rep = ds/scaleDia_*voidfraction*magUr/(nuf+SMALL);

                    // calc fluid drag Coeff
                    Cd = sqr(0.63 + 4.8/sqrt(Rep));

                    // calc model coefficient Xi
                    scalar Xi = 3.7 - 0.65 * exp(-sqr(1.5-log10(Rep))/2);

                    // calc particle's drag
                    dragCoefficient = 0.125*Cd*rho
                                     *M_PI
                                     *ds*ds     
                                     *scaleDia_ 
                                     *pow(voidfraction,(2-Xi))*magUr
                                     *scaleDrag_;
                    if (modelType_=="B")
                        dragCoefficient /= voidfraction;

                    drag = dragCoefficient*Ur; //total drag force!

                    //Split forces
                    if(splitImplicitExplicit_)
                    {
                        UfluidFluct  = Ufluid - U_[cellI];
                        UsFluct      = Us     - UsField_[cellI];
                        dragExplicit = dragCoefficient*(UfluidFluct - UsFluct); //explicit part of force
                    }
                }

                if(verbose_ && index >-1 && index <102)
                {
                    Pout << "index = " << index << endl;
                    Pout << "Us = " << Us << endl;
                    Pout << "Ur = " << Ur << endl;
                    Pout << "ds/scale = " << ds/scaleDia_ << endl;
                    Pout << "rho = " << rho << endl;
                    Pout << "nuf = " << nuf << endl;
                    Pout << "voidfraction = " << voidfraction << endl;
                    Pout << "Rep = " << Rep << endl;
                    Pout << "Cd = " << Cd << endl;
                    Pout << "drag (total) = " << drag << endl;
                    if(splitImplicitExplicit_)
                    {
                        Pout << "UfluidFluct = " << UfluidFluct << endl;
                        Pout << "UsFluct = " << UsFluct << endl;
                        Pout << "dragExplicit = " << dragExplicit << endl;
                    }
                }

                //Set value fields and write the probe
                if(probeIt_)
                {
                    #include "setupProbeModelfields.H"
                    vValues.append(drag);   //first entry must the be the force
                    vValues.append(Ur);
                    sValues.append(Rep);
                    sValues.append(Cd);
                    sValues.append(voidfraction);
                    particleCloud_.probeM().writeProbe(index, sValues, vValues);
                }
            }
            // set force on particle
            if(treatExplicit_) for(int j=0;j<3;j++) expForces()[index][j] += drag[j];
            else   //implicit treatment, taking explicit force contribution into account
            {
               for(int j=0;j<3;j++) 
               { 
                    impForces()[index][j] += drag[j] - dragExplicit[j]; //only consider implicit part!
                    expForces()[index][j] += dragExplicit[j];
               }
            }
            
            for(int j=0;j<3;j++) DEMForces()[index][j] += drag[j];
        }
    //}
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
