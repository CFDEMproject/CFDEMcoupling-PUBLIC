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

    transport model for scalar properties

Class
    Foam::scalarTransportModel

SourceFiles
    scalarTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef scalarTransportModel_H
#define scalarTransportModel_H


#include "fvCFD.H"
#include "cfdemCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class scalarTransportModel Declaration
\*---------------------------------------------------------------------------*/

class scalarTransportModel
{
protected:

    // Protected data
    cfdemCloud&          particleCloud_;

    const fvMesh&        mesh_;

    const IOdictionary   scalarTransportProperties_;

    bool                 verbose_;      //switch for debug output to screen

    bool                 ignore_;       //ignore updating the model, i.e., deactivate

public:

    //- Runtime type information
    TypeName("scalarTransportModel");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            scalarTransportModel,
            dictionary,
            (
                const dictionary& dict,
                cfdemCloud& sm
            ),
            (dict,sm)
        );


    // Constructors

        //- Construct from components
        scalarTransportModel
        (
            const dictionary& dict,
            cfdemCloud& sm
        );


    // Destructor
    	virtual ~scalarTransportModel();

    // Selector

        static autoPtr<scalarTransportModel> New
        (
            const dictionary& dict,
            cfdemCloud& sm
        );


	// Member Functions

        virtual void createFields();

        virtual void update();

        virtual const volScalarField& sourceField();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
