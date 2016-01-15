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
    deltaHEvap_("deltaHEvap", dimLength*dimLength/dimTime/dimTime, 1),
    tEvap_("tEvap", dimTime, 1)
{

    deltaHEvap_  = dict_.lookup("deltaHEvap");
    tEvap_       = dict_.lookup("tEvap");

    if(parameterVap_.size()<5)
        FatalError <<"phaseChangeModel: parameterVap_.size()<5! Provide more parameters to this model. \n" 
                   << abort(FatalError);  

    if(alphaImExSplit_<0 || alphaImExSplit_>1)
        FatalError <<"alphaImExSplit must be between 0 and 1. \n" 
                   << abort(FatalError);  
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
    
    //update the saturation field and cp quantities
    cpFromField_ = fromField.cpCarrier();
    cpToField_   = toField.cpCarrier();

    forAll(mSaturation_.internalField(), iter)
        mSaturation_.internalField()[iter] = pVapor( temp.internalField()[iter] )  / temp[iter] / Rvap_; 

    //update the reference quantities
    volScalarField tempF =      1.0   / ( 
                                              fromField.m()
                                            + toField.m() * fromField.rho() / toField.rho()
                                            +               fromField.rho() / fromField.rhoCarrier()
                                        ); //phi_liquid ... (global) liquid volume fraction
                                           //divided by fromField.m()
                                        

    //leaving mass rate - implicit/explicit term (divided by fromFiel.m())
    fromField.mSource().internalField()    -= (1-alphaImExSplit_) * tempF.internalField() * fromField.m().internalField()
                                              / tEvap_.value()  //characteristic evaporation time
                                            * (   
                                                 mSaturation_.internalField() 
                                               - toField.m().internalField()
                                                *fromField.rhoCarrier()
                                              );

    fromField.mSourceKImpl().internalField() -= alphaImExSplit_ * tempF.internalField() 
                                              / tEvap_.value()  //characteristic evaporation time
                                            * (   
                                                 mSaturation_.internalField() 
                                               - toField.m().internalField()
                                                *fromField.rhoCarrier()
                                              );

    tempF *= fromField.m() / tEvap_.value(); //phi_liquid ... (global) liquid volume fraction
                                             //divided by evaporation time scale

    //entering mass rate - explicit & implicit term
    toField.mSource().internalField()     += tempF * mSaturation_.internalField();
    toField.mSourceKImpl().internalField()-= tempF * toField.rhoCarrier();

    //set the rate
    mSource_.internalField()  = tempF 
                              * (   
                                    mSaturation_.internalField() 
                                  - toField.m().internalField()
                                  * toField.rhoCarrier() 
                                );
}

void phaseChangeModel::setEnthalpySource(const eulerianScalarField& Temperature) const
{
    //update the heat source
    Temperature.mSource().internalField() -= mSource_.internalField()
                                           * (   deltaHEvap_.value() 
                                               - Temperature.m().internalField() * (cpFromField_ - cpToField_)
                                             );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
