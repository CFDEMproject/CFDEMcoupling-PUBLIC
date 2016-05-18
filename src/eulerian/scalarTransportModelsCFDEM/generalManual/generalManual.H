/*---------------------------------------------------------------------------*\
License

    This is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This code is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with this code.  If not, see <http://www.gnu.org/licenses/>.

    Copyright (C) 2014- Stefan Radl, TU Graz, Austria

    Description
    Evolves an arbitrary number of scalar quantities, exchanges sources
    with the granular phase

\*---------------------------------------------------------------------------*/

#ifndef generalManual_H
#define generalManual_H

#include "scalarTransportModel.H"
#include "forceModel.H"
#include "eulerianScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class generalManual Declaration
\*---------------------------------------------------------------------------*/

class generalManual
:
    public scalarTransportModel
{

protected:

    dictionary                      propsDict_;

    word                            phiFieldName_;

    word                            voidfractionFieldName_;

    const wordList                  eulerianFieldList_;

    autoPtr<eulerianScalarField>*   eulerianFields_;

    scalar                          ScT_; // - Turbulent Schmidt number (optional, default = 0.7)

    scalar                          PrT_; // - Turbulent Prandtl number (optional, default = 0.7)

    int                             idTemp_;    //id with temperature field

    //Extra fields to represent a mixture in the void space
    bool                            updateMixtureProperties_;   // switch to indicate whether phase properties will be updated or not
    volScalarField                  rhoMix_;                  //field (if needed) to hold mixture density
    volScalarField                  cpRho_;                  //field (if needed) to hold volumetric cp-value
             


public:

    //- Runtime type information
    TypeName("generalManual");


    // Constructors

    //- Construct from components
    generalManual
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

    ~generalManual();


    // Member Functions

    const                eulerianScalarField& eulerianTemperatureF();

    const                eulerianScalarField& eulerianScalarF(int);

    void                 createFields();

    void                 update();

    void                 setSources();

    void                 evolveFields();

    volScalarField&      sourceField(int i);



};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
