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

    Copyright (C) 2015 - Stefan Radl, TU Graz, Austria

Application / Class
    phaseChangeModel
    
Description
    Base class for phase change modeling (e.g., evaporation)
    TODO: currently handles only evaporation!

\*---------------------------------------------------------------------------*/

#ifndef phaseChangeModel_H
#define phaseChangeModel_H

#include "cfdemCloud.H"
#include "forceModel.H"
#ifndef versionExt32
#include "fvOptionList.H"
#endif
#include "eulerianScalarField.H"
#include "generalManual.H"
#include "OFstream.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class phaseChangeModel Declaration
\*---------------------------------------------------------------------------*/

class phaseChangeModel
{

protected:

    // Protected data
        const dictionary&       dict_;

        cfdemCloud&             particleCloud_;

        const word              modelName_;     

        mutable volScalarField  mSaturation_;   // saturation concentration field

        mutable volScalarField  mSource_;       // volumetric source of transported quantity (for the leaving phase)
    
        int                     speciesIDFrom_; // id of leaving phase

        int                     speciesIDTo_;   // id of receiving phase

        scalarList              parameterVap_;  // parameters for vapor pressure

        scalar                  Rvap_;          // Gas constant of vapor (in J/kg/K)

        scalar                  alphaImExSplit_; // 0...1 (indicate split into implicit and explicit part)

        // Key Evaporation properties
        mutable scalar          cpFromField_;
        mutable scalar          cpToField_;
        dimensionedScalar       deltaHEvap_;    //units: K; deltaHEvap = deltaH_v / cp_g
        dimensionedScalar       tEvap_;         //units: s; tEvap_ = d_d^2 / (6*Sh*D_Vapor)

        //Reporting
        mutable int         verboseDiskIntervall_;

        mutable int         verboseDiskCounter_;

        mutable OFstream*   sPtr_; 


        inline bool verboseToDisk() const
        { 
            if(verboseDiskIntervall_<=0) return false;

            verboseDiskCounter_++;
            if(verboseDiskCounter_>=verboseDiskIntervall_)
            {
                verboseDiskCounter_=0;
                return true; 
            }
            else
                return false;
        };

        //Private/protected member functions
        inline double pVapor(double T) const
        {   //T in K, p in Pa, from http://ddbonline.ddbst.de/AntoineCalculation/AntoineCalculationCGI.exe?component=Water
            const double A = parameterVap_[0];
            const double B = parameterVap_[1];
            const double C = parameterVap_[2];
            const double Tmin = parameterVap_[3];
            const double Tmax = parameterVap_[4];
    
            double result = 0;
            if( T>Tmin )
                result = 133.322368 * pow(10.0, A - B / ( C + min(Tmax,T) ) );
            return result;
        }

public:

    //- Runtime type information
    TypeName("phaseChangeModel");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            phaseChangeModel,
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
        phaseChangeModel
        (
            const dictionary&   dict,
            cfdemCloud&         sm,
            word                modelType,
            int                 modelID
        );


    // Destructor

        virtual ~phaseChangeModel();


    // Selector

        static autoPtr<phaseChangeModel> New
        (
            const dictionary&   dict,
            cfdemCloud&         sm,
            word                modelType,
            int                 modelID
        );


    // Member Functions
    int  fromID() const {return speciesIDFrom_;};
    int  toID()   const {return speciesIDTo_;};

    void update(const volScalarField&, const volScalarField&, const eulerianScalarField&, const eulerianScalarField&) const ;

    void setEnthalpySource(const eulerianScalarField&) const;

    void bound(autoPtr<phaseChangeModel>*) const {};
    void initialzeSummation(word typeName, word  logFileName) const;
    void computeIntegral (volScalarField& explicitEulerSource) const;

    // Access
    inline volScalarField& mSource() const {return mSource_;}; //returns the source

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
