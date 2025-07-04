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

#include "scalarTransportModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

autoPtr<scalarTransportModel> scalarTransportModel::New
(
    const dictionary& dict, //Not used at the moment
    cfdemCloud& sm
)
{
    IOdictionary tempDict
    (
        IOobject
        (
            "scalarTransportProperties",
            sm.mesh().time().constant(),
            sm.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    word scalarTransportModelType
    (
        tempDict.lookup("scalarTransportModel")
    );

    Info<< "Selecting scalarTransportModel "
		 << scalarTransportModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(scalarTransportModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "scalarTransportModel::New(const dictionary&, const spray&) : "
            << endl
            << "    unknown scalarTransportModelType type "
            << scalarTransportModelType
            << ", constructor not in hash table" << endl << endl
            << "    Valid scalarTransportModel types are :"
            << endl;
        Info<< dictionaryConstructorTablePtr_->toc()
            << abort(FatalError);
    }

    autoPtr<scalarTransportModel> h(cstrIter()(tempDict,sm));

    sm.registryM().addProperty("scalarTransportModel_"+h().type()+"_index", 0);

    return h;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
