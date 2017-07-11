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
\*---------------------------------------------------------------------------*/

#include "error.H"
#include "registryModel.H"
#include "dataExchangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(registryModel, 0);

defineRunTimeSelectionTable(registryModel, dictionary);

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void registryModel::addProperty(word name, scalar value) const
{
    //Pout << "trying to add the property" << name << " to registry, with a value of " << value << endl;
    int index(checkIfExists(name,value));
    if(index > -1)
    {
        scalar newValue(max(value,valuesScalar_[index]));
        valuesScalar_[index] = newValue;
        //Pout << "   the property " << name << " was found in registry, the value is updated (if bigger than old) to " << newValue << endl;
    }
    else
    {
        namesScalar_.append(name);
        valuesScalar_.append(value);
        //Pout << "   the property " << name << " is (NEW!!!) put to the registry, the value is set to " << value << endl;
    }
}

scalar registryModel::getProperty(word name) const
{
    //Pout << "trying to receive the property" << name << " from the registry..." << endl;
    int index(checkIfExists(name,scalar(1.)));
    if(index > -1)
    {
        //Pout << "   the property" << name << " was found in registry, the value is " << valuesScalar_[index] << endl;
        return valuesScalar_[index];
    }
    else
    {
        //Pout << "   the property" << name << " is not found in the registry!" << endl;
        return -1.; 
    }
}

inline int registryModel::checkIfExists(word name, scalar value) const
{
    for(int i=0;i<namesScalar_.size();i++)
        if(namesScalar_[i]==name) return i;
    return -1;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
registryModel::registryModel
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    dict_(dict),
    particleCloud_(sm),
    namesScalar_(0),
    valuesScalar_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

registryModel::~registryModel()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
