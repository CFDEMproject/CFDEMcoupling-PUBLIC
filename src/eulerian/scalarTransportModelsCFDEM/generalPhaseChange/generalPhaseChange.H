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

    Copyright (C) 2015 -Stefan Radl, TU Graz, Austria

    Description
    Derived from generalManual transport Model
    Evolves an arbitrary number of scalar quantities, exchanges sources
    with the granular phase, and implements phase-change physics 
    between the Eulerian phases

\*---------------------------------------------------------------------------*/

#ifndef generalPhaseChange_H
#define generalPhaseChange_H

#include "generalManual.H"
#include "forceModel.H"
#include "phaseChangeModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class generalPhaseChange Declaration
\*---------------------------------------------------------------------------*/

class generalPhaseChange
:
    public generalManual
{

protected:

    dictionary                      phaseChangeDict_;

    const wordList                  phaseChangeModelList_;

    autoPtr<phaseChangeModel>*      phaseChangeModels_;

public:

    //- Runtime type information
    TypeName("generalPhaseChange");


    // Constructors
    //- Construct from components
    generalPhaseChange
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor
    ~generalPhaseChange();

    //Access function
    const phaseChangeModel& phaseChangeModelRef(int); //must be const to satisfy compiler!

    // Member Functions

    void                 createFields();

    void                 update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
