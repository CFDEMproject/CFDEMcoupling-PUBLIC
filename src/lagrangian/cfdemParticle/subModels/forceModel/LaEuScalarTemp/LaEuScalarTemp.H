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

    two way LaEu Scalar Exchange between DEM and CFD
    convective heat transfer model following
    Li and Mason (2000), A computational investigation of transient heat
    transfer in  pneumatic transport of granular particles, Pow.Tech 112

    This model will put all the source into the explicit coupling term
    for the fluid-side integration (might be unstable)

Class
    LaEuScalarTemp

SourceFiles
    LaEuScalarTemp.C

\*---------------------------------------------------------------------------*/

#ifndef LaEuScalarTemp_H
#define LaEuScalarTemp_H

#include "forceModel.H"
#include "averagingModel.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class LaEuScalarTemp Declaration
\*---------------------------------------------------------------------------*/

class LaEuScalarTemp
:
    public forceModel
{
protected:

    dictionary propsDict_;

    bool compressible_;

    const word EuFieldName_;

    word tempFieldName_;

    const volScalarField& T_;            // ref to temperature field

    word voidfractionFieldName_;

    const volScalarField& voidfraction_;    // ref to voidfraction field

    scalar maxSource_;                           // max (limited) value of src field

    word velFieldName_;

    const volVectorField& U_;

    scalar lambda_;                              // fluid thermal conductivity [W/(m*K)]

    scalar Cp_;                                  // specific heat capacity [W*s/(kg*K)]

    mutable volScalarField CpField_;

    mutable volScalarField LambdaField_;

    word NuCorrelation_;
    double (LaEuScalarTemp::*Nusselt)(double Re, double Pr, double vf) const;

public:

    //- Runtime type information
    TypeName("LaEuScalarTemp");

    // Constructors

        //- Construct from components
        LaEuScalarTemp
        (
            const dictionary& dict,
            cfdemCloud& sm,
            word name = word("")
        );

    // Destructor

        ~LaEuScalarTemp();


    // Member Functions
        void setForce() const;

    //--- Nusselt correlations ---

// Nusselt number correlation from
//    Li and Mason (2000), A computational investigation of transient heat
//    transfer in  pneumatic transport of granular particles, Pow.Tech 112
        double NuLiMason(double Re, double Pr, double vf) const
        {
            double Nup(1.);
            double n(3.5);

            if (Re < 200)
                Nup = 2. + 0.6 * pow(vf,n) * sqrt(Re) * pow(Pr,0.33);
            else if (Re < 1500)
                Nup = 2. + (0.5 * sqrt(Re) + 0.02 * pow(Re,0.8)) * pow(vf,n) * pow(Pr,0.33);
            else
                Nup = 2. + 0.000045 * pow(vf,n) * pow(Re,1.8);

            return Nup;
        }

// Nusselt number correlation from
//    Deen, N.G. et al., Review of direct numerical simulation of
//    fluid–particle mass, momentum and heat transfer in dense gas–solid flows.
//    Chemical Engineering Science 116 (2014) 710–724
        double NuDeen(double Re, double Pr, double vf) const
        {
            //NOTE: This function is fitted for Reynolds numbers between 0 and 100!!!
            double Nup(1.);
            double PrPowOneThird = pow(Pr, 0.3333333333);

            Nup = (7.00 - 10.0*vf + 5*vf*vf)
                * (1.00 + 0.17 * pow(Re, 0.2) * PrPowOneThird)
                + (1.33 - 2.31*vf + 1.16*vf*vf) * pow(Re, 0.7) * PrPowOneThird ;

            return Nup;
        }

// Nusselt number correlation from
//    Gunn, D. (1978). Transfer of heat or mass to particles in fixed and fluidised beds.
//    International Journal of Heat and Mass Transfer, 21(4):467–476.
        double NuGunn(double Re, double Pr, double vf) const
        {
            //NOTE: This function is fitted for Reynolds numbers between 0 and 100!!!
            double Nup(1.);
            double PrPowOneThird = pow(Pr, 0.3333333333);

            Nup = (7.00 - 10.0*vf + 5*vf*vf)
                * (1.00 + 0.70 * pow(Re, 0.2) * PrPowOneThird)
                + (1.33 - 2.40*vf + 1.20*vf*vf) * pow(Re, 0.7) * PrPowOneThird ;

            return Nup;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
