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

#include "noTransport.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(noTransport, 0);

addToRunTimeSelectionTable
(
	scalarTransportModel,
	noTransport,
	dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
// Construct from components
noTransport::noTransport
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    scalarTransportModel(dict,sm)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
noTransport::~noTransport()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void noTransport::createFields()
{}

// ************************************************************
void noTransport::update()
{}

// ************************************************************
const volScalarField& noTransport::sourceField()
{
    FatalErrorIn("const volScalarField& noTransport::sourceField() ")
        << "this source field is NOT implemented, and hence MUST NOT be called!" << abort(FatalError);
    return volScalarField::null();
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
