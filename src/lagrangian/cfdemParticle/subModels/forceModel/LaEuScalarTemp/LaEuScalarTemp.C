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

#include "LaEuScalarTemp.H"
#include "addToRunTimeSelectionTable.H"
#include "dataExchangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(LaEuScalarTemp, 0);

addToRunTimeSelectionTable
(
    forceModel,
    LaEuScalarTemp,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
LaEuScalarTemp::LaEuScalarTemp
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    tempFieldName_(propsDict_.lookup("tempFieldName")),
    tempField_(sm.mesh().lookupObject<volScalarField> (tempFieldName_)),
    voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),
    voidfractionField_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    maxSource_(1e30),
    velFieldName_(propsDict_.lookup("velFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    partTempName_(propsDict_.lookup("partTempName")),
    partTemp_(NULL),
    partHeatFluxName_(propsDict_.lookup("partHeatFluxName")),
    partHeatFlux_(NULL),
    lambda_(readScalar(propsDict_.lookup("lambda"))),
    Cp_(readScalar(propsDict_.lookup("Cp"))),
    densityFieldName_(propsDict_.lookup("densityFieldName")),
    rho_(sm.mesh().lookupObject<volScalarField> (densityFieldName_))
{
    allocateMyArrays();

    if (propsDict_.found("maxSource"))
    {
        maxSource_=readScalar(propsDict_.lookup ("maxSource"));
        Info << "limiting eulerian source field to: " << maxSource_ << endl;
    }
    particleCloud_.checkCG(false);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

LaEuScalarTemp::~LaEuScalarTemp()
{
    delete partTemp_;
    delete partHeatFlux_;
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
void LaEuScalarTemp::allocateMyArrays() const
{
    // get memory for 2d arrays
    double initVal=0.0;
    particleCloud_.dataExchangeM().allocateArray(partTemp_,initVal,1);  // field/initVal/with/lenghtFromLigghts
    particleCloud_.dataExchangeM().allocateArray(partHeatFlux_,initVal,1);
}
// * * * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * * //

void LaEuScalarTemp::setForce() const
{
    // do nothing
}

void LaEuScalarTemp::manipulateScalarField(volScalarField& EuField) const
{
    // realloc the arrays
    allocateMyArrays();

    // reset Scalar field
    EuField.internalField() = 0.0;

    // get DEM data
    particleCloud_.dataExchangeM().getData(partTempName_,"scalar-atom",partTemp_);

    // get viscosity field
    #ifdef comp
        const volScalarField& nufField = particleCloud_.turbulence().mu() / rho_;
    #else
        const volScalarField& nufField = particleCloud_.turbulence().nu();
    #endif

    // calc La based heat flux
    vector Us;
    scalar magUr;
    scalar alpha;
    scalar rs;
    scalar As;
    scalar nuf;
    scalar Rep;
    scalar Pr;
    scalar n = 3.5; // model parameter
    scalar Nup;

    for(int index = 0;index < particleCloud_.numberOfParticles(); ++index)
    {
        //if(particleCloud_.regionM().inRegion()[index][0])
        //{
            label cellI = particleCloud_.cellIDs()[index][0];

            if(cellI >= 0)
            {
                // calc relative velocity
                Us = particleCloud_.velocity(index);
                magUr = mag(U_[cellI]-Us);
                alpha = voidfractionField_[cellI];
                rs = particleCloud_.radius(index);
                As = 4*rs*rs*M_PI;
                nuf = nufField[cellI];
                Rep = 2*rs*magUr/nuf;
                Pr = Cp_*nuf*rho_[cellI]/lambda_;

                if (Rep < 200)
                {
                    Nup = 2+0.6*pow(alpha,n)*sqrt(Rep)*pow(Pr,0.33);
                }
                else if (Rep < 1500)
                {
                    Nup = 2+0.5*pow(alpha,n)*sqrt(Rep)*pow(Pr,0.33)
                                 +0.02*pow(alpha,n)*pow(Rep,0.8)*pow(Pr,0.33);
                }
                else
                {
                    Nup = 2+0.000045*pow(alpha,n)*pow(Rep,1.8);
                }
                scalar h = lambda_*Nup/(2*rs);

                // calc convective heat flux [W]
                scalar partHeatFlux = h * As * (tempField_[cellI] - partTemp_[index][0]);
                partHeatFlux_[index][0] = partHeatFlux;


                /*if(index == 101)
                {
                    Info << "partHeatFlux = " << partHeatFlux << endl;
                    Info << "magUr = " << magUr << endl;
                    Info << "As = " << As << endl;
                    Info << "nuf = " << nuf << endl;
                    Info << "Rep = " << Rep << endl;
                    Info << "Pr = " << Pr << endl;
                    Info << "Nup = " << Nup << endl;
                    Info << "alpha = " << alpha << endl;
                    Info << "partTemp_[index][0] = " << partTemp_[index][0] << endl  ;
                    Info << "ptempField_[cellI] = " << tempField_[cellI] << endl  ;
                }*/
            }
        //}
    }

    particleCloud_.averagingM().setScalarSum
    (
        EuField,
        partHeatFlux_,
        particleCloud_.particleWeights(),
        NULL
    );

    // scale with -1/(Vcell*rho*Cp)
    EuField.internalField() /= -rho_.internalField()*Cp_*EuField.mesh().V();

    // limit source term
    scalar EuFieldInCell;
    forAll(EuField,cellI)
    {
        EuFieldInCell = EuField[cellI];

        if(mag(EuFieldInCell) > maxSource_ )
        {
             EuField[cellI] = sign(EuFieldInCell) * maxSource_;
        }
    }

    Info << "total convective particle-fluid heat flux [W] (Eulerian) = " << gSum(EuField*rho_*Cp_*EuField.mesh().V()) << endl;

    // give DEM data
    particleCloud_.dataExchangeM().giveData(partHeatFluxName_,"scalar-atom", partHeatFlux_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
