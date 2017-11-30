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

Class
    dataExchangeModel

SourceFiles
    dataExchangeModel.C

\*---------------------------------------------------------------------------*/

#ifndef dataExchangeModel_H
#define dataExchangeModel_H

#include "fvCFD.H"
#include "cfdemCloud.H"

#include "mpi.h"
#include "math.h"
#include <lammps.h>
#include <pair.h>
#include <force.h>
#include <atom.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dataExchangeModel Declaration
\*---------------------------------------------------------------------------*/

class dataExchangeModel
{

protected:

    // Protected data
        const dictionary& dict_;

        cfdemCloud& particleCloud_;

        int maxNumberOfParticles_;

        mutable int couplingStep_;

        scalar DEMts_;

        int couplingInterval_;

        const int timeIndexOffset_;

    // Protected member functions

public:

    //- Runtime type information
    TypeName("dataExchangeModel");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            dataExchangeModel,
            dictionary,
            (
                const dictionary& dict,
                cfdemCloud& sm
            ),
            (dict,sm)
        );


    // Constructors

        //- Construct from components
        dataExchangeModel
        (
            const dictionary& dict,
            cfdemCloud& sm
        );


    // Destructor

        virtual ~dataExchangeModel();


    // Selector

        static autoPtr<dataExchangeModel> New
        (
            const dictionary& dict,
            cfdemCloud& sm
        );


    // Member Function
        void setNumberOfParticles(int) const;

        void setNumberOfClumps(int) const;

        inline const int& maxNumberOfParticles() const {return maxNumberOfParticles_;};

        template <typename T>
        void getData
        (
            word name,
            word type,
            T ** const& field
        ) const { getData(name,type,field,couplingStep_-1); }

        virtual void getData
        (
            word name,
            word type,
            double ** const& field,
            label step
        ) const = 0;

        virtual void getData
        (
            word name,
            word type,
            int ** const& field,
            label step
        ) const=0;

        virtual void giveData
        (
            word name,
            word type,
            double ** const& field,
            const char* datatype="double"
        ) const = 0;

        //====
        // double **
        virtual void allocateArray(double**&, double, int, int) const;
        virtual void allocateArray(double**&, double, int, const char* ="nparticles") const;
        virtual void destroy(double**,int) const;

        //====
        // int **
        virtual void allocateArray(int**&, int, int, int) const;
        virtual void allocateArray(int**&, int, int, const char* ="nparticles") const;
        virtual void destroy(int**,int) const;
        //====

        //====
        // int *
        virtual void allocateArray(int*&, int, int) const;
        virtual void destroy(int*) const;
        //====

        //====
        // double *
        virtual void allocateArray(double*&, double, int) const;
        virtual void destroy(double*) const;
        //====

        virtual bool couple(int) const;

        virtual scalar timeStepFraction() const;

        inline int couplingStep() const {return couplingStep_;};

        inline const scalar& DEMts() const {return DEMts_;};

        inline int couplingInterval() const {return couplingInterval_;};

        inline scalar couplingTime() const {return couplingInterval_*DEMts_;};

        inline scalar TSstart() const { return particleCloud_.mesh().time().startTime().value()
                                        + (couplingStep_-1) * couplingTime();};

        inline scalar TSend() const {return particleCloud_.mesh().time().startTime().value()
                                        + (couplingStep_) * couplingTime();};

        inline int DEMstepsTillT(scalar t) const {return (t - (particleCloud_.mesh().time().value() - couplingTime())  + SMALL) / DEMts_;};

        inline void checkTSsize() const
        {
            if(particleCloud_.mesh().time().deltaT().value() >  couplingInterval_ * DEMts_ + SMALL)
            {
                Info << "particleCloud_.mesh().time().deltaT().value() = " << particleCloud_.mesh().time().deltaT().value() << endl;
                Info << "couplingInterval_ = " << couplingInterval_ << endl;
                Info << "DEMts_ = " << DEMts_ << endl;
                FatalError<<"\nError - CFD time-step bigger than coupling time (= DEM time step * coupling interval)!\n"<< abort(FatalError);
            }
            if ( std::fabs(( round(couplingTime()/particleCloud_.mesh().time().deltaT().value()) * particleCloud_.mesh().time().deltaT().value() )-couplingTime()) > SMALL )
            {
                Info << "particleCloud_.mesh().time().deltaT().value() = " << particleCloud_.mesh().time().deltaT().value() << endl;
                Info << "couplingInterval_ = " << couplingInterval_ << endl;
                Info << "DEMts_ = " << DEMts_ << endl;
                FatalError <<"\nWarning - Coupling time (= DEM time step * coupling interval) is not a multiple of  CFD time-step!\n"<< abort(FatalError);
            }
            if(!particleCloud_.allowCFDsubTimestep())
                if(particleCloud_.mesh().time().deltaT().value() < couplingInterval_ * DEMts_ + SMALL)
                    FatalError<<"\nYour models require: CFD time-step = coupling interval (= DEM time step * coupling interval)! \n"<< abort(FatalError);

            // warn if sub-TS
            if (particleCloud_.mesh().time().deltaT().value() < couplingTime() - SMALL)
                Warning << "You are using sub-time-steps (i.e. CFD TS < coupling time)! Check your settings properly." << endl;
        }

        /*inline bool checkExactTiming() const
        {
            return false;
        }*/

        //void checkNClumpTypes() const {};

        inline void readDEMtsfromDict(dictionary& propsDict)
        {
            DEMts_ = readScalar(propsDict.lookup("DEMts"));
            checkTSsize();
        }

        inline bool doCoupleNow() const
        {
            // couple at first sub TS
            if ((particleCloud_.mesh().time().timeIndex()-timeIndexOffset_)*particleCloud_.mesh().time().deltaT().value()-SMALL
                > couplingStep_*DEMts_*couplingInterval_)
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        virtual int getNumberOfParticles() const;
        virtual int getNumberOfClumps() const;
        virtual int getNumberOfTypes() const;
        virtual double* getTypeVol() const;

        inline void setPositions(label n,double* pos) const
        {
            for (int i=0;i<n;i++)
                for (int j=0;j<3;j++)
                    particleCloud_.positions_[i][j]=pos[i*3+j];
        };
        inline void setCellIDs(label n,int* ID) const
        {
            for (int i=0;i<n;i++)
                particleCloud_.cellIDs_[i][0]=ID[i];
        };

        inline void setPositionsCM(label n,double* pos) const
        {
            particleCloud_.setPositionsCM(n,pos);
        };
        inline void setCellIDsCM(label n,int* ID) const
        {
            particleCloud_.setCellIDsCM(n,ID);
        };

        virtual word myType() const=0;

        virtual void setCG() const { FatalError<<"\nsetCG() not executed correctly!\n"<< abort(FatalError); };

        void setCGs(LAMMPS_NS::LAMMPS *lmp) const 
        { 
            particleCloud_.setCG(lmp->force->cg(int(1e14)));
            
            Info << "LIGGGHTS hosts " << lmp->atom->ntypes << " atom types." << endl;

            #if defined(typeSecificCG)
            Warning << "\nYou are using untested code (type specific coarse graining)!!!\n" 
                    << "\nThis will most probably give wrong results for most models!!!\n"<< endl;
            for(int iTyp=1; iTyp<=lmp->atom->ntypes; iTyp++)
                particleCloud_.setCGTypeSpecific(iTyp,lmp->force->cg(iTyp));
            #endif
        };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
