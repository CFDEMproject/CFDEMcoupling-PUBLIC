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

#include "error.H"
#include "eulerianScalarField.H"
#include "OFversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(eulerianScalarField, 0);

defineRunTimeSelectionTable(eulerianScalarField, dictionary);

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
eulerianScalarField::eulerianScalarField
(
    const dictionary&   dict,
    cfdemCloud&         sm,
    word                modelType,
    int                 modelID
)
:
    dict_(dict),
    particleCloud_(sm),
    fieldName_(modelType),
    cpVolumetricFieldName_(dict_.lookupOrDefault<word>("cpVolumetricFieldName", "na")),
    cpVolumetric_(dict_.lookupOrDefault<scalar>("cpVolumetric", 0.0)),
    updateMixtureProperties_(dict_.lookupOrDefault<bool>("updateMixtureProperties", false)),
    rho_(dict_.lookupOrDefault<scalar>("rho"+fieldName_, -1)),
    rhoCarrier_(dict_.lookupOrDefault<scalar>("rhoCarrier", -1)),
    cp_(dict_.lookupOrDefault<scalar>("cp"+fieldName_, -1)),
    cpCarrier_(dict_.lookupOrDefault<scalar>("cpCarrier", -1)),
    m_
    (   IOobject
        (
            fieldName_,
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh()
    ),
    mSource_
    (   IOobject
        (
            fieldName_+"Source",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh()
    ),
    mSourceKImpl_
    (   IOobject
        (
            fieldName_+"SourceKImpl",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0.0*mSource_ /( m_ + dimensionedScalar("dummy", m_.dimensions(), 1e-32) ) //initi with zero
    ),
    fieldType_("undefined")
    #ifndef versionExt32
    ,fvOptions_(sm.mesh())
    #endif
{


    if ( m_.dimensions() == dimensionSet(0, 0, 0, 1, 0) )
    {
        speciesID_ = -1;
        fieldType_ = "temperature";
        Info << "eulerianScalarField:: found a Temperature field! " << endl;
    }
    else
    {
        speciesID_ = modelID;
        fieldType_ = "species";
        Info << "eulerianScalarField:: found a species field, will assign speciesID: " 
             << speciesID_
             << endl;
    }

    #ifndef versionExt32
    fvOptions_.reset(dict.subDict("fvOptions"+fieldName_));
    #endif

    if( (cpVolumetricFieldName_=="na"||!updateMixtureProperties_) && cpVolumetric_<=0.0)
        FatalError <<"You did not specify a cpVolumetricFieldName (or you do not updateMixtureProperties) and also cpVolumetric is zero (or negative)! Either provide the field name, or set cpVolumetric to a reasonable value. \n" 
                   << abort(FatalError);    

    if(speciesID_>-1 && updateMixtureProperties_ && (rho_<=0 || cp_<=0) )
        FatalError <<"You like to update the phase properties, but density and cp of the eulerianScalarField with name '" 
                   << fieldName_
                   <<"' are not specified or zero. \n" 
                   << abort(FatalError);    

    if(speciesID_>-1 && updateMixtureProperties_ && (rhoCarrier_<=0 || cpCarrier_<=0) )
        FatalError <<"You like to update the phase properties, but density and cp of the carrier phase are not specified or zero \n" 
                   << abort(FatalError);    
                   
    //Report options for cp 
    if(fieldType_=="temperature")
    {
        if(cpVolumetric_!=0.0 && cpVolumetricFieldName_!="na")
        FatalError <<"eulerianScalarField:: You have specified 'cpVolumetric' and 'cpVolumetricFieldName' in a dictionary in '/constant'. This might be confusing. Please unset one of these two inputs to avoid confusion. \n" 
                   << abort(FatalError);    

        if(cpVolumetricFieldName_=="na" || !updateMixtureProperties_) //use also if mixture properties are not updated
            Info << "eulerianScalarField:: will use the following FIXED VOLUMETRIC HEAT CAPACITY: " 
                 << cpVolumetric_ << " [J/K/m³]" << endl;
        else
            Info << "eulerianScalarField:: will use the a SPATIALLY-VARAIBLE VOLUMETRIC HEAT CAPACITY with name: " << cpVolumetricFieldName_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

eulerianScalarField::~eulerianScalarField()
{}

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //
void eulerianScalarField::pullCloudFields() const 
{
    //Temporary field for collecting the sources
    volScalarField tmpSource
    (   IOobject
        (
            "tmpSource",
            particleCloud_.mesh().time().timeName(),
            particleCloud_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mSource_*0.0
    );

    volScalarField tmpSourceImpl
    (   IOobject
        (
            "tmpSourceImpl",
            particleCloud_.mesh().time().timeName(),
            particleCloud_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mSourceKImpl_*0.0
    );

    //reset sources, and compute the sources due to the cloud
    //Will accumulate all sources for all force models
    for (int iModel=0; iModel<particleCloud_.nrForceModels(); iModel++)
    {
        particleCloud_.forceM(iModel).manipulateScalarField(tmpSource, tmpSourceImpl, speciesID_);
        if(iModel==0)
        {
            mSource_        = tmpSource;
            mSourceKImpl_   = tmpSourceImpl;
        }
        else
        {
            mSource_       += tmpSource;
            mSourceKImpl_  += tmpSourceImpl;
        }
    }

}

// ************************************************************
void eulerianScalarField::update(surfaceScalarField phi, volScalarField voidfraction, volScalarField nuEff, scalar Sc, bool limitDiffusion) const 
{
    scalar oneByCpVolumetric = 1./(cpVolumetric_+SMALL);
    //Normalize source in case we have a temperature field
    if(fieldType_=="temperature")
    {
        if(cpVolumetricFieldName_=="na" || !updateMixtureProperties_) //use also if mixture properties are not updated
        {
            mSource_ *= oneByCpVolumetric;
            mSourceKImpl_ *= oneByCpVolumetric;
        }
        else
        {
            const volScalarField& cpVolumetricField_(particleCloud_.mesh().lookupObject<volScalarField> (cpVolumetricFieldName_));

            #if defined(version40) || defined(versionv1612plus)
            mSource_.primitiveFieldRef() /= cpVolumetricField_.primitiveField()+SMALL;
            mSourceKImpl_.primitiveFieldRef() /= cpVolumetricField_.primitiveField()+SMALL;
            #else
            mSource_.internalField() /= cpVolumetricField_.internalField()+SMALL;
            mSourceKImpl_.internalField() /= cpVolumetricField_.internalField()+SMALL;
            #endif
        }
    }

    word divScheme("div(phi,m)");
    word laplacianScheme("laplacian(D,m)");

    if( speciesID_== -1)
    {
        divScheme       = "div(phi,T)";
        laplacianScheme = "laplacian((DT*voidfraction),T)";
    }

    // calc diffusionCoeff field
    volScalarField a=nuEff/Sc*voidfraction;
    if(limitDiffusion)
    {
        forAll(a,cellI)
            if(voidfraction[cellI] <2*SMALL)
                a[cellI]=0.;
    }

    // solve scalar transport equation
    fvScalarMatrix mEqn
    (

       fvm::ddt(voidfraction, m_)               //This is the material derivative in a modified form
     - fvm::Sp(fvc::ddt(voidfraction), m_)      //Needed since phi is (U_face * voidfraction)!
     + fvm::div(phi, m_, divScheme)             //This phi must be SUPERFICIAL! (i.e., U_face * voidfraction)!
     - fvm::Sp(fvc::div(phi), m_)

     ==

       fvm::laplacian(a, m_, laplacianScheme) 
     + mSource_
     + fvm::Sp(mSourceKImpl_, m_)
     #ifndef versionExt32
     + fvOptions_(m_)
     #endif
    );

    mEqn.relax();
    #ifndef versionExt32
    fvOptions_.constrain(mEqn);
    #endif
    mEqn.solve();


}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
