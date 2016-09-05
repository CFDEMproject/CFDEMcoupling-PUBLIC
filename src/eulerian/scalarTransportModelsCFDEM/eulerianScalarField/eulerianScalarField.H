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

Application / Class
    eulerianScalarField
    
Description
    Base class for moment transport model & FVM-based solver.
    WARNING:    the surfaceScalarField 'phi' must be based on U*voidfraction,
                and not on U!

\*---------------------------------------------------------------------------*/

#ifndef eulerianScalarField_H
#define eulerianScalarField_H

#include "cfdemCloud.H"
#include "forceModel.H"
#ifndef versionExt32
#include "fvOptionList.H"
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class eulerianScalarField Declaration
\*---------------------------------------------------------------------------*/

class eulerianScalarField
{

protected:

    // Protected data
        const dictionary&       dict_;

        cfdemCloud&             particleCloud_;

        const word              fieldName_;

        word                    cpVolumetricFieldName_;   // name of the field holding the volumetric heat capacity
                                                          // if not provided, must specify a constant holding this value

        scalar                  cpVolumetric_;            // the (constant) volumetric heat capacity (optional)

        //species properties (if required, these properties will be set, default = -1
        bool                    updateMixtureProperties_;   // switch to indicate whether phase properties will be updated or not
        scalar                  rho_;                     // the density of this field
        scalar                  rhoCarrier_;              // the density of the carrier field
        scalar                  cp_;                      // the heat capacity of this field
        scalar                  cpCarrier_;               // the heat capacity of the carrier field

        mutable volScalarField  m_;                       // transported quantity

        mutable volScalarField  mSource_;   // source of transported quantity; 
                                            // might involve part of implicit term

        mutable volScalarField  mSourceKImpl_;   // pre-factor for implicit source
                                                 // will be added with negative sign in transport eqn.
                                                 // will shadowed

        mutable word            fieldType_;

        #ifndef versionExt32
        mutable fv::optionList  fvOptions_;
        #endif

        int                     speciesID_;


public:

    //- Runtime type information
    TypeName("eulerianScalarField");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            eulerianScalarField,
            dictionary,
            (
                const dictionary&   dict,
                cfdemCloud&         sm,
                word                modelType,
                int                 modelID
            ),
            (dict,sm,modelType,modelID)
        );


    // Constructors

        //- Construct from components
        eulerianScalarField
        (
            const dictionary&   dict,
            cfdemCloud&         sm,
            word                modelType,
            int                 modelID
        );


    // Destructor

        virtual ~eulerianScalarField();


    // Selector

        static autoPtr<eulerianScalarField> New
        (
            const dictionary&   dict,
            cfdemCloud&         sm,
            word                modelType,
            int                 modelID
        );


    // Member Functions
    void pullCloudFields() const ;

    void update(surfaceScalarField, volScalarField, volScalarField, scalar Sc, bool limitDiffusion=false) const ;

    void bound(autoPtr<eulerianScalarField>*) const {};

    // Access
    word fieldType() const                 { return fieldType_; };

    inline bool updateMixtureProperties() const   { return updateMixtureProperties_; };
    inline scalar rho()          const            { return rho_; };
    inline scalar rhoCarrier()   const            { return rhoCarrier_; };
    inline scalar cp()           const            { return cp_; };
    inline scalar cpCarrier()    const            { return cpCarrier_; };

    inline volScalarField& m()            const {return m_;}; //returns the transported quantity

    inline volScalarField& mSource()      const {return mSource_;}; //returns the source
    inline volScalarField& mSourceKImpl() const {return mSourceKImpl_;}; //returns the implicit prefactor of the source

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
