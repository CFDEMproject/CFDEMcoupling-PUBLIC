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

    cloud class managing DEM data for CFD-DEM coupling and IB representation

Class
    Foam::cfdemCloudIB derived from cfdemCloud

SourceFiles
    cfdemCloudIB.C

\*---------------------------------------------------------------------------*/

#ifndef cfdemCloudIB_H
#define cfdemCloudIB_H

#include "cfdemCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cfdemCloudIB Declaration
\*---------------------------------------------------------------------------*/

class cfdemCloudIB
:
    public cfdemCloud
{
protected:

    mutable double **angularVelocities_;
    label pRefCell_;
    scalar pRefValue_;
    
    mutable double **DEMTorques_;

    mutable bool haveEvolvedOnce_;
    mutable bool skipLagrangeToEulerMapping_;
    mutable bool skipAfter_;
    mutable int  timeStepsToSkip_;
    mutable bool calculateTortuosity_;

    vector flowDir_;

    //for refinement - needed to ensure correct probing around the interface
    volScalarField*     frontMeshRefineField_;
    bool                frontMeshRefine_;

public:

    // Constructors

        //- Construct from components
        cfdemCloudIB
        (
            const fvMesh& mesh
        );

    // Destructor

        ~cfdemCloudIB();


    // Member Functions
        void giveDEMdata();

        void getDEMdata();

        bool reAllocArrays() const;

        inline double ** DEMTorques() const;
    
        bool evolve(volScalarField&,volScalarField&);

        void calcVelocityCorrection(volScalarField&,volVectorField&,volScalarField&,volScalarField&);

        virtual void setParticleVelocity(volVectorField&);

      // Access
        vector angularVelocity(int);

        inline double ** angularVelocities() const
        {
            return angularVelocities_;
        };
        void setInterFace(volScalarField& interFace);   //this field can be used for refinement outside of the cloud

        double getTortuosity(vector dir);

       void setRefinementField(volScalarField* refine_); //this pointer is meant for a field that is used refinement outside of the cloud. can point to the same field as interFace

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
