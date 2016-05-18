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

#include "generalPhaseChange.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(generalPhaseChange, 0);

addToRunTimeSelectionTable
(
	scalarTransportModel,
	generalPhaseChange,
	dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
generalPhaseChange::generalPhaseChange
(
    const dictionary& dict,
    cfdemCloud&       sm
)
:
    generalManual(dict,sm),
    phaseChangeDict_(propsDict_.subDict("PhaseChangeParameters")),
    phaseChangeModelList_(phaseChangeDict_.lookup("phaseChangeModels"))
{

    phaseChangeModels_ = new autoPtr<phaseChangeModel>[phaseChangeModelList_.size()];
    for (int i=0;i<phaseChangeModelList_.size();i++)
    {
        phaseChangeModels_[i] = phaseChangeModel::New
        (
            phaseChangeDict_,
            sm,
            phaseChangeModelList_[i],
            i
        );
    }

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
generalPhaseChange::~generalPhaseChange()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void generalPhaseChange::createFields()
{}

// ************************************************************
void generalPhaseChange::update()
{
    //Re-set the sources due to particle-fluid interactions
    generalManual::setSources();

    //Apply the phaseChange operation (loop through list of models)
    //phaseChangeModels must ADD any sources to the eulerianScalarFields 
    //(since there might be sources due to particle-fluid interactions)
    const volScalarField&     voidfraction(particleCloud_.mesh().lookupObject<volScalarField> (voidfractionFieldName_));
    for (int i=0;i<phaseChangeModelList_.size();i++)
    {
        int idFieldFrom = phaseChangeModelRef(i).fromID();
        int idFieldTo   = phaseChangeModelRef(i).toID();
        phaseChangeModelRef(i).update(voidfraction,                 eulerianTemperatureF().m(),
                                      eulerianScalarF(idFieldFrom), eulerianScalarF(idFieldTo));

        phaseChangeModelRef(i).setEnthalpySource(eulerianTemperatureF());
    }

    //to stuff for standard scalar transport
    generalManual::evolveFields();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
const phaseChangeModel& generalPhaseChange::phaseChangeModelRef(int i)
{
     return phaseChangeModels_[i];
}

} // End namespace Foam

// ************************************************************************* //
