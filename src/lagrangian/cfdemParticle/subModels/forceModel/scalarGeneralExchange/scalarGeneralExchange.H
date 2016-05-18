/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
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

    Two way general scalar exchange between DEM and CFD
    convective heat and species transfer model. The standard model is that of
    Deen, N.G. et al., Review of direct numerical simulation of 
    fluid–particle mass, momentum and heat transfer in dense gas–solid flows.
    Chemical Engineering Science 116 (2014) 710–724
    This correlation is based on that of Gunn (1978).

    Alternatively, the correclation of
    Li and Mason (2000), A computational investigation of transient heat
    transfer in  pneumatic transport of granular particles, Pow.Tech 112
    can be activated. However, this correlation is not suitable for 
    dense granular flows.

    This model allows an implicit/explicit split of the coupling term.
    The implicit/explicit splitting is realized in a force sub-model

    WARNING:
    This model REQUIRES the 'generalManual' speciesTransportModel

Class
    scalarGeneralExchange

SourceFiles
    scalarGeneralExchange.C

Contributing author and copyright holder of this model/file
    Copyright, 2015     Stefan Radl, TU Graz (radl@tugraz.at)

\*---------------------------------------------------------------------------*/

#ifndef scalarGeneralExchange_H
#define scalarGeneralExchange_H

#include "forceModel.H"
#include "averagingModel.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class scalarGeneralExchange Declaration
\*---------------------------------------------------------------------------*/

class scalarGeneralExchange
:
    public forceModel
{
protected:

    mutable dictionary  propsDict_;

    const IOdictionary  scalarTransportProperties_;

    dictionary          generalPropsDict_;

    word                voidfractionFieldName_;

    word                velFieldName_;

    word                tempFieldName_;               
    word                partTempName_;
    word                partHeatFluxName_;
    word                partHeatTransCoeffName_;
    word                partHeatFluidName_;

    mutable double      **partDat_;                   // Lagrangian array ExtCode-->CFDEM

    mutable double      **partDatFlux_;               // Lagrangian array CFDEM-->ExtCode (explicit contribution)
    mutable double      **partDatTransCoeff_;         // Lagrangian array CFDEM-->ExtCode (implicit contribution)
    mutable double      **partDatFluid_;              // Lagrangian array CFDEM-->ExtCode (implicit contribution)

    mutable double       **partTemp_;

    mutable double      **partDatTmpExpl_;            // Lagrangian array - for explicit source to fluid eqns
    mutable double      **partDatTmpImpl_;            // Lagrangian array - for implicit source to fluid eqns

    mutable double      **partDatSaturation_;                   // Lagrangian array
    mutable double      **partCoolingFlux_;
    
    mutable bool        validPartFlux_;               //indicator if found, based on heat, but also used for species 
                                                      //also indicates EXplicit Coupling
    mutable bool        validPartTransCoeff_;         //indicator if found, based on heat, but also used for species
                                                      //also indicates IMplicit Coupling
    mutable bool        validPartFluid_;              //indicator if found, based on heat, but also used for species

    mutable bool        haveTemperatureEqn_;          //indicator for temperature field or not

    mutable bool        useLiMason_;                  //switch to activate calculation using Li-Mason
    mutable bool        useGeneralCorrelation_;       //switch to activate calculation using a generalized correlation
    mutable scalarList  generalCorrelationParameters_; //parameter for general correlation

    scalar              lambda_;                      // fluid thermal conductivity [W/(m*K)]
    scalar              Prandtl_;                     // Prandtl number

    //Species Word Lists
    mutable wordList    eulerianFieldNames_;           //List with Eulerian fields to exchange (handed over)

    mutable scalarList  particleSpeciesValue_;         //list with scalar to indicate particle property is available
    const wordList      partSpeciesNames_;
    const wordList      partSpeciesFluxNames_;
    const wordList      partSpeciesTransCoeffNames_;
    const wordList      partSpeciesFluidNames_;

    scalarList          DMolecular_; 

    //information related to external register
    mutable int                 partHeatFluxPositionInRegister_;
    mutable std::vector<int>    partSpeciesFluxPositionInRegister_;
    mutable int                 partHeatTransCoeffPositionInRegister_;
    mutable std::vector<int>    partSpeciesTransCoeffPositionInRegister_;
    mutable int                 partHeatFluidPositionInRegister_;
    mutable std::vector<int>    partSpeciesFluidPositionInRegister_;

    //Scalar properties
    mutable scalar      maxSource_;                   // max (limited) value of src field

    mutable scalar      scaleDia_;

    void allocateMyArrays(scalar initialValue) const;
    void setupModel() const;
    void setPointersToExternalArrays(   word nameFlux,          int positionFlux,
                                        word nameTransCoeff,    int positionTransCoeff,
                                        word nameFluid,         int positionFluid
                                    ) const;

    mutable double (scalarGeneralExchange::*Nusselt)(double Re, double Pr, double voidfraction) const;
    double NusseltLiMason(double Re, double Pr, double voidfraction) const;
    double NusseltDeenEtAl(double Re, double Pr, double voidfraction)  const;
    double NusseltGeneralCorrelation(double Re, double Pr, double voidfraction)  const;

public:

    //- Runtime type information
    TypeName("scalarGeneralExchange");

    // Constructors

        //- Construct from components
        scalarGeneralExchange
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

        //- Construct from components
        scalarGeneralExchange
        (
            const dictionary& dict,
            cfdemCloud& sm,
            word        dictName
        );

    // Destructor

        ~scalarGeneralExchange();


    // Member Functions
        void setForce() const;

        void manipulateScalarField(volScalarField&, volScalarField&, int) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
