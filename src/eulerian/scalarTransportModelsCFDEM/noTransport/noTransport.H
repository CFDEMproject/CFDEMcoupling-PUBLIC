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

\*---------------------------------------------------------------------------*/

#ifndef noTransport_H
#define noTransport_H

#include "scalarTransportModel.H"
#include "forceModel.H"
#include "eulerianScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class noTransport Declaration
\*---------------------------------------------------------------------------*/

class noTransport
:
    public scalarTransportModel
{

public:

    //- Runtime type information
    TypeName("none");


    // Constructors

    //- Construct from components
    noTransport
    (
            const dictionary& dict,
            cfdemCloud& sm
   );

    // Destructor

    	~noTransport();


    // Member Functions

    void createFields();

    void update();

    const volScalarField& sourceField();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
