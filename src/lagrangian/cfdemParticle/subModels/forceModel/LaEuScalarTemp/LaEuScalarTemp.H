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

    word tempFieldName_;

    const volScalarField& T_;            // ref to temperature field

    word voidfractionFieldName_;

    const volScalarField& voidfraction_;    // ref to voidfraction field

    scalar maxSource_;                           // max (limited) value of src field

    word velFieldName_;

    const volVectorField& U_;

    word partTempName_;

    mutable double **partTemp_;                  // Lagrangian array

    word partHeatFluxName_;

    mutable double **partHeatFlux_;              // Lagrangian array
    
    bool           validPartHeatFlux_;           //indicator if found

    word partHeatTransCoeffName_;

    mutable double **partHeatTransCoeff_;        // Lagrangian array
    
    bool           validPartHeatTransCoeff_;     //indicator if found
    
    word partHeatFluidName_;

    mutable double **partHeatFluid_;             // Lagrangian array
    
    bool           validPartHeatFluid_;           //indicator if found

    scalar lambda_;                              // fluid thermal conductivity [W/(m*K)]

    scalar Cp_;                                  // specific heat capacity [W*s/(kg*K)]
    
    bool compressible_;
    
    bool doExchange_;

    void allocateMyArrays() const;

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

        void manipulateScalarField(volScalarField&) const;

        void commToDEM() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
