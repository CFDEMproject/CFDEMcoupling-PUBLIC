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
#include  <vector>
#include "interpolation.H"


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
    
    mutable bool haveEvolvedOnce_;
    mutable bool skipLagrangeToEulerMapping_;

    //===
    mutable bool useHFDIBM_;
    mutable bool checkPeriodicCells_;
    
    //HFDIBM members    
    std::vector< std::vector< Foam::vector > > interpP_[2];
    std::vector< std::vector< label > > interpCells_[2];
    std::vector< std::vector< label > > surfaceCells_;
    std::vector<double> remoteInterpolationPoints;
    double* remoteScalarInterp_;
    double* remoteVectorInterp_;
    std::vector<int> remoteInterpLegend;
    std::vector<label> remoteInterpCells;
    dictionary  HFDIBMinterpDict_;
    int me;
    int nprocs;
    int* displs_;
    std::vector<int> switchRemoteToFD_;
    void setInterpolationPoints();
    void reallocateHFDIBMarrays();
    void communicateRemoteInterpolationPoints(std::vector<double> localRemoteInterpolationPoints_);
    void vectorInterpolateRemote( volVectorField   &V, volVectorField  &Vs);
    void scalarInterpolateRemote(  volScalarField  &S, volScalarField  &Ss);
    //===

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
        void getDEMdata();

        bool reAllocArrays() const;
    
        bool evolve(volScalarField&);

        void calcVelocityCorrection(volScalarField&,volVectorField&,volScalarField&,volScalarField&);

        void setParticleVelocity(volVectorField&);

      // Access
        vector angularVelocity(int);

        inline double ** angularVelocities() const
        {
            return angularVelocities_;
        };


       void interpolateVectorField( volVectorField   &V,  //main field
                                    volVectorField  &Vs   //imposed field
                                  );

       void interpolateScalarField( volScalarField   &S,  //main field
                                    volScalarField  &Ss   //imposed field
                                  );

       void checkInterfaceFlowRate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
