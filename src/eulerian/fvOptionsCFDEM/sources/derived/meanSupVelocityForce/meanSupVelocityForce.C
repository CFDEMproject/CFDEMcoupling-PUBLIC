/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright (C) 2011 OpenFOAM Foundation
                                Copyright (C) 2012-     DCS Computing GmbH,Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "meanSupVelocityForce.H"
#include "addToRunTimeSelectionTable.H"

#include "DimensionedField.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(meanSupVelocityForce, 0);

    addToRunTimeSelectionTable
    (
        option,
        meanSupVelocityForce,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::meanSupVelocityForce::meanSupVelocityForce
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    meanVelocityForce(sourceName, modelType, dict, mesh),
    twoPhase_( coeffs_.lookupOrDefault("twoPhase",false) ),
    alpha_
    (   
        IOobject
        (
            "voidfractionPrev",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,//MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("ones", dimensionSet(0,0,0,0,0), 1)
    ),
    coupled_( coeffs_.lookupOrDefault("coupled",false) ),
    voidfraction_
    (   
        IOobject
        (
            "voidfractionPrev",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,//MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("ones", dimensionSet(0,0,0,0,0), 1)
    ),
    modelName_(modelType),
    alphaMin_(coeffs_.lookupOrDefault("alphaMin",0.0))

{
    Warning << "THE FVOPTION meanSupVelocityForce has not been tested/validated!!! " << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::fv::meanSupVelocityForce::correct(volVectorField& U)
{
    if (twoPhase_)
    {
        word alphaName = coeffs_.lookup("alphaField");
        alpha_ = mesh_.lookupObject<volScalarField>(alphaName);
    }

    if (coupled_)
    {
        word voidfractionName = coeffs_.lookup("voidfractionField");
        voidfraction_ = mesh_.lookupObject<volScalarField>(voidfractionName);
    }
    
    // making sure that alpha is only set almost-only fluid regions if desired
    forAll(alpha_,cellI)
    {
        if(alpha_[cellI]<alphaMin_)
            alpha_[cellI]=0;
    }

    const scalarField& rAU = rAPtr_().internalField();

    // Integrate flow variables over cell set
    scalar rAUave = 0.0;
    const scalarField& cv = mesh_.V();
    scalar totV = 0.0;
    forAll(cells_,i)
    {
        label cellI = cells_[i];
        scalar volCell = cv[cellI];
        totV += volCell*alpha_[cellI]; 
        rAUave += rAU[cellI]*volCell*alpha_[cellI]; 
    }

    // Collect accross all processors
    reduce(rAUave, sumOp<scalar>());
    reduce(totV, sumOp<scalar>());
    V_=totV;

    // Volume averages
    rAUave /= V_;

    scalar magUbarAve = this->magUbarAve(U);

    // Calculate the pressure gradient increment needed to adjust the average
    // flow-rate to the desired value
    dGradP_ = relaxation_*(mag(Ubar_) - magUbarAve)/rAUave;

    // Apply correction to velocity field
    if (modelName_=="B" || modelName_=="Bfull")
    {
        forAll(cells_, i)
        {
            label cellI = cells_[i];
            U[cellI] += flowDir_*rAU[cellI]*dGradP_*alpha_[cellI];
        }
    }
    else
    {
        forAll(cells_, i)
        {
            label cellI = cells_[i];
            U[cellI] += voidfraction_[cellI]*flowDir_*rAU[cellI]*dGradP_*alpha_[cellI];
        }
    }

    scalar gradP = gradP0_ + dGradP_;

    Info<< "Pressure gradient source: uncorrected Ubar = " << magUbarAve
        << ", pressure gradient = " << gradP << endl;

    writeProps(gradP);

    Warning << "Pressure gradient force is neglected in this model!!" << endl; 

    // The following lines would compensate the error that occurs due to the splitting 
    // of the pressure gradient. However, the particleCloud_ object is currently not
    // accessible here. 
    /*scalar ds(0.0);
    scalar Vs(0.0);
    label cellI=0;
    for(int index = 0;index <  particleCloud_.numberOfParticles(); index++)
    {
        //if(mask[index][0])
        //{
            cellI = particleCloud_.cellIDs()[index][0];
            
            if (cellI > -1) // particle found on this processor
            {
                //Calc the particle volume
                ds = 2*particleCloud_.radius(index);
                Vs = ds*ds*ds*M_PI/6;
                
                // set force on particle
                for(int j=0;j<3;j++) 
                {
                    // calc particle's static pressure gradient force
                    particleCloud_.DEMForces()[index][j] -= Vs*gradP*flowDir_[j];
                }
            }
    }*/


}

void Foam::fv::meanSupVelocityForce::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    DimensionedField<vector, volMesh> Su
    (
        IOobject
        (
            name_ + fieldNames_[fieldI] + "Sup",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", eqn.dimensions()/dimVolume, vector::zero)
    );

    scalar gradP = gradP0_ + dGradP_;

    if (modelName_=="B" || modelName_=="Bfull")
    {
        UIndirectList<vector>(Su, cells_) = flowDir_*gradP;
    }
    else
    {
        //UIndirectList<vector>(Su, cells_) = voidfraction_*flowDir_*gradP*alpha_; // org version does not work for zones

        UIndirectList<vector>(Su, cells_) = vector::zero;
        label cellI(-1);
        forAll(cells_,i)
        {
            cellI = cells_[i];
            Su[cellI] = voidfraction_[cellI]*flowDir_*gradP*alpha_[cellI];
        }
    }

    eqn += Su;

}

Foam::scalar Foam::fv::meanSupVelocityForce::magUbarAve
(
    const volVectorField& U
) const
{
    scalar magUbarAve = 0.0;

    const scalarField& cv = mesh_.V();
    forAll(cells_, i)
    {
        label cellI = cells_[i];
        scalar volCell = cv[cellI];
        magUbarAve += (flowDir_ & U[cellI])*volCell*alpha_[cellI]*voidfraction_[cellI];
    }

    reduce(magUbarAve, sumOp<scalar>());

    magUbarAve /= V_;

    return magUbarAve;
}




// ************************************************************************* //
