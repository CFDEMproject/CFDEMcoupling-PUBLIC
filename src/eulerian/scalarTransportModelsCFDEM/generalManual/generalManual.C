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

#include "generalManual.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(generalManual, 0);

addToRunTimeSelectionTable
(
	scalarTransportModel,
	generalManual,
	dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
generalManual::generalManual
(
    const dictionary& dict,
    cfdemCloud&       sm
)
:
    scalarTransportModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    phiFieldName_(propsDict_.lookupOrDefault<word>("phiFieldName", "phi")),
    voidfractionFieldName_(propsDict_.lookupOrDefault<word>("voidfractionFieldName", "voidfraction")),
    eulerianFieldList_(propsDict_.lookup("eulerianFields")),
    ScT_(0.7),
    PrT_(0.7),
    idTemp_(-1),
    updateMixtureProperties_(propsDict_.lookupOrDefault<bool>("updateMixtureProperties", true)),
    rhoMix_
    (   IOobject
        (
            propsDict_.lookupOrDefault<word>("rhoMixFieldName","rhoMixDefault"),
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("dummy", dimensionSet(1,-3,0,0,0), -1)
    ),
    cpRho_
    (   IOobject
        (
            propsDict_.lookupOrDefault<word>("cpVolumetricFieldName","cpRhoDefault"),
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("dummy", dimensionSet(1,-1,-2,-1,0), -1)
    )
{

    propsDict_.readIfPresent("ScT", ScT_);
    propsDict_.readIfPresent("PrT", PrT_);
    
    Info << "generalManual:: Using the following turbulent dispersion coefficients: ScT = " << ScT_ << " and PrT " << PrT_ << endl;

    eulerianFields_ = new autoPtr<eulerianScalarField>[eulerianFieldList_.size()];
    for (int i=0;i<eulerianFieldList_.size();i++)
    {
        if(eulerianFieldList_[i]=="T")
            idTemp_= i;

        eulerianFields_[i] = eulerianScalarField::New
        (
            propsDict_,
            sm,
            eulerianFieldList_[i],
            i
        );
    }

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
generalManual::~generalManual()
{
    delete [] eulerianFields_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void generalManual::createFields()
{}

// ************************************************************
void generalManual::setSources()
{
    //Loop through all Eulerian fields and and set fields
    for (int i=0;i<eulerianFieldList_.size();i++)
            eulerianScalarF(i).pullCloudFields();

    //Communicate exchange terms to DEM
    for (int iModel=0; iModel<particleCloud_.nrForceModels(); iModel++)
        particleCloud_.forceM(iModel).commToDEM();

    //Send Sources to External Code (i.e., Lagrangian arrays handled by LIGGGHTS)
    particleCloud_.giveUSERdata();
}

// ************************************************************
void generalManual::evolveFields()
{
    //update the fields
    if(updateMixtureProperties_)
    {
        if(eulerianScalarF(0).fieldType()=="Temperature" )
            FatalError <<"generalManual: first eulerianField is temperature, but we need a species. Please re-order your eulerianFields in the input dict. \n" 
                       << abort(FatalError);  

        forAll(rhoMix_, iter)
        {
            double denominator              = 1./eulerianScalarF(0).rhoCarrier();
            rhoMix_[iter]   = 1.0;
            cpRho_[iter]    = eulerianScalarF(0).cpCarrier();
            for (int i=0;i<eulerianFieldList_.size();i++)
            {
              if(eulerianScalarF(i).fieldType()!="temperature")
              {
                denominator += eulerianScalarF(i).m()[iter]
                             / eulerianScalarF(i).rho();

                rhoMix_[iter] += eulerianScalarF(i).m()[iter];
                cpRho_[iter]  += eulerianScalarF(i).m()[iter]
                                                *eulerianScalarF(i).cp();
              }
            }
            rhoMix_[iter]  /= denominator;
            cpRho_[iter]   /= denominator;
        }

    }

    //==============================
    // get references
    const surfaceScalarField& phi(particleCloud_.mesh().lookupObject<surfaceScalarField> (phiFieldName_));
    const volScalarField&     voidfraction(particleCloud_.mesh().lookupObject<volScalarField> (voidfractionFieldName_));
    //==============================

    //Loop through all eulerian fields and update them
    for (int i=0;i<eulerianFieldList_.size();i++)
    {
            if(eulerianScalarF(i).fieldType()=="Temperature")
                eulerianScalarF(i).update(phi, voidfraction, particleCloud_.turbulence().nuEff(), PrT_);
            else 
                eulerianScalarF(i).update(phi, voidfraction, particleCloud_.turbulence().nuEff(), ScT_);
    }
}

// ************************************************************
void generalManual::update()
{
    setSources();
    evolveFields();
}

// ************************************************************
volScalarField& generalManual::sourceField(int i)
{
    return eulerianScalarF(i).mSource();
}

// ************************************************************
const eulerianScalarField& generalManual::eulerianScalarF(int i)
{
     return eulerianFields_[i];
}

// ************************************************************
const eulerianScalarField& generalManual::eulerianTemperatureF()
{
    if(idTemp_<0)
        FatalError <<"You did not specify a temperature field with name 'T', but you are requesting this field in generalManual::eulerianTemperatureF. \n" 
                   << abort(FatalError);  

    return eulerianFields_[idTemp_];
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace Foam

// ************************************************************************* //
