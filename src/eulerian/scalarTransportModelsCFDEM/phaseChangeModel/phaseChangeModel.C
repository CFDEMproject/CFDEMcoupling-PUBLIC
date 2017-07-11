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
#include "phaseChangeModel.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(phaseChangeModel, 0);

defineRunTimeSelectionTable(phaseChangeModel, dictionary);

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
phaseChangeModel::phaseChangeModel
(
    const dictionary&   dict,
    cfdemCloud&         sm,
    word                modelType,
    int                 modelID
)
:
    dict_(dict),
    particleCloud_(sm),
    modelName_(modelType),
    mSaturation_
    (   IOobject
        (
            modelName_+"Saturation",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("dummy", dimensionSet(0,0,0,0,0), -1)
    ),
    mSource_
    (   IOobject
        (
            modelName_+"Rate",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("dummy", dimensionSet(0,0,-1,0,0), -1)
    ),
    speciesIDFrom_(0),
    speciesIDTo_(1),
    parameterVap_(dict_.lookup("parameterVap")),
    Rvap_(dict_.lookupOrDefault<scalar>("Rvap", 461.5)),
    alphaImExSplit_(dict_.lookupOrDefault<scalar>("alphaImExSplit", 0.5)),
    cpFromField_(0.0),
    cpToField_(0.0),
    deltaHEvap_(dict_.lookup("deltaHEvap")),
    tEvap_(dict_.lookup("tEvap")),
    verboseDiskIntervall_(-1), //zero or negative--> deaktivate output to disk
    verboseDiskCounter_(0),
    sPtr_(NULL)
{
    if(parameterVap_.size()<5)
        FatalError <<"phaseChangeModel: parameterVap_.size()<5! Provide more parameters to this model. \n" 
                   << abort(FatalError);  

    if(alphaImExSplit_<0 || alphaImExSplit_>1)
        FatalError <<"alphaImExSplit must be between 0 and 1. \n" 
                   << abort(FatalError);  
    if(dict_.found("verboseDiskIntervall"))
        verboseDiskIntervall_=readScalar(dict_.lookup("verboseDiskIntervall"));

    if (verboseDiskIntervall_>0)
    {
        Info << "phaseChangeModel will report to disk with intervall " << verboseDiskIntervall_ << endl;
        initialzeSummation(typeName, "phaseChange.logDat");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

phaseChangeModel::~phaseChangeModel()
{}

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //
void phaseChangeModel::update(const volScalarField&      voidfraction,   //this is 1-particleFraction field
                              const volScalarField&      temp,           //this is temperature field of the fluid phase
                              const eulerianScalarField& fromField,      const eulerianScalarField& toField) const 
{
    //To implement phase change model updates that directly affect "fromField", and "toField"
    //MUST ADD to sources (not reset sources!)
    //MUST be per m³ TOTAL volume, since scalar transport solver is based on this
    
    //update the saturation field and cp quantities
    cpFromField_ = fromField.cpCarrier();
    cpToField_   = toField.cpCarrier();

    forAll(mSaturation_, iter)
        mSaturation_[iter] = pVapor( temp[iter] )  / temp[iter] / Rvap_; 

    //update the reference quantities

    volScalarField tempF =      voidfraction 
	                              / ( 
                                              fromField.m()
                                            + toField.m() * fromField.rho() / toField.rho()
                                            +               fromField.rho() / fromField.rhoCarrier()
                                        ); //phi_liquid ... (global) liquid volume fraction (per m³ total!)
                                           //divided by fromField.m() 


    //leaving mass rate - implicit/explicit term (divided by fromFiel.m())
#if defined(version40) || defined(versionv1612plus)
    fromField.mSource()    -= (1-alphaImExSplit_) * tempF * fromField.m()
#else
    fromField.mSource().internalField()     -= (1-alphaImExSplit_) * tempF.internalField()  * fromField.m()
#endif
                                              / tEvap_.value()  //characteristic evaporation time
                                            * (   
                                                 mSaturation_ 
                                               - toField.m()
                                                *fromField.rhoCarrier()
                                              );

#if defined(version40) || defined(versionv1612plus)
    fromField.mSourceKImpl() -= alphaImExSplit_ * tempF
#else
    fromField.mSourceKImpl().internalField()  -= alphaImExSplit_ * tempF.internalField() 
#endif
                                              / tEvap_.value()  //characteristic evaporation time
                                            * (   
                                                 mSaturation_ 
                                               - toField.m()
                                                *fromField.rhoCarrier()
                                              );

    tempF *= fromField.m() / tEvap_.value(); //phi_liquid ... (global) liquid volume fraction
                                             //divided by evaporation time scale


    //entering mass rate - explicit & implicit term

#if defined(version40) || defined(versionv1612plus)
        toField.mSource()     += tempF * mSaturation_;
        toField.mSourceKImpl()-= tempF * toField.rhoCarrier();
#else
        toField.mSource().internalField()       += tempF.internalField()   * mSaturation_;
        toField.mSourceKImpl().internalField()  -= tempF.internalField()   * toField.rhoCarrier();
#endif


    //set the rate
#if defined(version40) || defined(versionv1612plus)
    mSource_  = tempF 
#else
    mSource_.internalField()  = tempF   .internalField()
#endif
                              * (   
                                    mSaturation_
                                  - toField.m()
                                  * toField.rhoCarrier() 
                                );

    if(verboseToDisk())
        computeIntegral(mSource_);
}

//************************************************
void phaseChangeModel::setEnthalpySource(const eulerianScalarField& Temperature) const
{

    //update the heat source
#if defined(version40) || defined(versionv1612plus)
       Temperature.mSource() -= mSource_
                                           * (   deltaHEvap_.value() 
                                               - Temperature.m() * (cpFromField_ - cpToField_)
                                             );
#else
       Temperature.mSource().internalField() -= mSource_.internalField()
                                           * (   deltaHEvap_.value() 
                                               - Temperature.m().internalField() * (cpFromField_ - cpToField_)
                                             );
#endif


}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void phaseChangeModel::initialzeSummation(word typeName, word  logFileName) const
{
    if (Pstream::master())
    {
        fileName file_ =logFileName;
        fileName probeDir;
        fileName probeSubDir =  typeName;

        Info << "Integral quantity for model " <<  typeName << " will write to file " << file_ << endl;

        if (particleCloud_.mesh().name() != polyMesh::defaultRegion)
        {
            probeSubDir = probeSubDir/particleCloud_.mesh().name();
        }
        probeSubDir = "postProcessing"/probeSubDir/particleCloud_.mesh().time().timeName();

        if (Pstream::parRun())
        {
            // Put in undecomposed case
            // (Note: gives problems for distributed data running)
            probeDir = particleCloud_.mesh().time().path()/".."/probeSubDir;
        }
        else
        {
            probeDir = particleCloud_.mesh().time().path()/probeSubDir;
        }


        // Create directory if does not exist.
        mkDir(probeDir);

        sPtr_ = new OFstream(probeDir+"/"+file_);

        *sPtr_ << '#' 
              << "Time" << "  " 
              << "sourceValue" << endl;
    }
}

//*******************************************************************
void phaseChangeModel::computeIntegral(volScalarField& explicitEulerSource) const
{
    volScalarField expEulerSrc=explicitEulerSource;
    particleCloud_.scaleWithVcell(expEulerSrc);
    scalar integralValue = gSum(expEulerSrc);

    if (Pstream::master() )
    {
      *sPtr_ << setprecision(IOstream::defaultPrecision()) ;
      *sPtr_ << particleCloud_.mesh().time().value() 
            << "  " //setw(IOstream::defaultPrecision() + 6)
            << integralValue
            << endl;
    }
}

} // End namespace Foam

// ************************************************************************* //
