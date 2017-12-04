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

\*---------------------------------------------------------------------------*/

#include "momCoupleModel.H"
#include "smoothingModel.H"
#include "meshMotionModel.H"
#include "averagingModel.H"
#include "clockModel.H"
#include "IOModel.H"
#include "voidFractionModel.H"
#include "locateModel.H"
#include "probeModel.H"
#include "registryModel.H"

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
inline void cfdemCloud::setCG(double cg) const
{ 
    cg_ = cg;
    Info << "cg is set to: " << cg_ << endl;
};

inline void cfdemCloud::setCGTypeSpecific(int type, double cg) const
{
    if(int(cgTypeSpecific_.size())!=type)   
        FatalError  << "setCGTypeSpecific attempts to set at location "
                    << int(cgTypeSpecific_.size())
                    << ", but caller requests to set at "
                    << type
                    << ". This is fatal. "
                    << abort(FatalError);
    cgTypeSpecific_.push_back(cg);
    if(cg<cg_*0.99999 || cg>cg_*1.00001) cgTypeSpecificDifferent = true;
    
    Info << "type specific cg for type " << type << " is set to: " << cg << endl;
}

inline const bool& cfdemCloud::impDEMdrag() const
{ 
    return impDEMdrag_;
};

inline const bool& cfdemCloud::impDEMdragAcc() const
{ 
    return impDEMdragAcc_;
};

inline const scalar& cfdemCloud::imExSplitFactor() const
{
	    return imExSplitFactor_;
};

inline const bool& cfdemCloud::treatVoidCellsAsExplicitForce() const
{
    return treatVoidCellsAsExplicitForce_;
}

inline const scalar& cfdemCloud::cg() const
{ 
    return cg_;
};

inline const scalar& cfdemCloud::cg(int type) const
{ 
    return cgTypeSpecific_[type];
};

inline const bool& cfdemCloud::ignore() const
{
    return ignore_;
}

inline const bool& cfdemCloud::writeTimePassed() const
{
    return writeTimePassed_;
}

inline const fvMesh& cfdemCloud::mesh() const
{
    return mesh_;
}

inline bool cfdemCloud::solveFlow() const
{
    return bool(solveFlow_);
}

inline bool cfdemCloud::solveScalarTransport() const
{
    return bool(solveScalarTransport_);
}


inline bool cfdemCloud::verbose() const
{
    return verbose_;
}

inline const IOdictionary& cfdemCloud::couplingProperties() const
{
    return couplingProperties_;
}

inline double ** cfdemCloud::positions() const
{
    return positions_;
}

inline double ** cfdemCloud::velocities() const
{
    return velocities_;
}

inline double ** cfdemCloud::fluidVels() const
{
    return fluidVel_;
}

inline double ** cfdemCloud::fAccs() const
{
    return fAcc_;
}

inline double ** cfdemCloud::impForces() const
{
    return impForces_;
}

inline double ** cfdemCloud::expForces() const
{
    return expForces_;
}

inline double ** cfdemCloud::DEMForces() const
{
    return DEMForces_;
}

inline double ** cfdemCloud::Cds() const
{
    return Cds_;
}

inline double ** cfdemCloud::radii() const
{
    return radii_;
}

inline double ** cfdemCloud::voidfractions() const
{
    return voidfractions_;
}

inline void cfdemCloud::get_radii(double **& values) const
{
//    Info << "set_radii level=" << numberOfParticles_ << endl;
    // make a copy of the array entries
//    for (int i=0;i<numberOfParticles_;i++)
//        radii_[0][i]=values[0][i];
    values=radii_;
}

inline double ** cfdemCloud::cellIDs() const
{
    return cellIDs_;
}

inline void cfdemCloud::get_cellIDs(double **& values) const
{
//    // make a copy of the array entries
//    for (int i=0;i<numberOfParticles_;i++)
//        cellIDs_[0][i]=values[0][i];
    values=cellIDs_;
}

inline double ** cfdemCloud::particleWeights() const
{
    return particleWeights_;
}

inline label Foam::cfdemCloud::body(int index)
{
    return index;
}

inline double cfdemCloud::particleVolume(int index)
{
    return particleV_[index][0];
}

inline scalar cfdemCloud::radius(int index)
{
    return radii_[index][0];
}

inline double cfdemCloud::d(int index)
{
    return 2*radii_[index][0];
}

inline double cfdemCloud::d32(bool recalc)
{
    if(d32_<0 || recalc)
    {
        scalar Ntot(0);
        scalar Dtot(0);
        scalar r(0);
        for(int index = 0;index <  numberOfParticles(); ++index)
        {
            r=radii_[index][0];
            Ntot+=2*r*r*r;
            Dtot+=r*r;
        }
        if(Ntot>SMALL) d32_ = Ntot/Dtot;
        else d32_ = 0.;
    }

    return d32_;
}

inline int cfdemCloud::numberOfParticles() const
{
    return numberOfParticles_;
}

inline bool cfdemCloud::numberOfParticlesChanged() const
{
    return numberOfParticlesChanged_;
}

inline int cfdemCloud::numberOfClumps() const
{
    Warning << "cfdemCloud::numberOfClumps() is not called correctly!" << endl;
    return -1;
}

/*inline void cfdemCloud::setMaxCellsPerParticle(int maxCellsPerParticle) const
{ 
    maxCellsPerParticle_ = maxCellsPerParticle;
};

inline int cfdemCloud::maxCellsPerParticle() const
{
    return maxCellsPerParticle_;
}*/

inline bool cfdemCloud::arraysReallocated() const
{
    return arraysReallocated_;
}

inline const wordList& cfdemCloud::forceModels()
{
    return forceModels_;
}

inline const locateModel& cfdemCloud::locateM() const
{
    return locateModel_;
}

inline const momCoupleModel& cfdemCloud::momCoupleM(int i) const
{
    return momCoupleModel_[i];
}

inline const dataExchangeModel& cfdemCloud::dataExchangeM() const
{
    return dataExchangeModel_;
}

inline const IOModel& cfdemCloud::IOM() const
{
    return IOModel_;
}

inline const probeModel& cfdemCloud::probeM() const
{
    return probeModel_;
}

inline const registryModel& cfdemCloud::registryM() const
{
    return registryModel_;
}

inline const voidFractionModel& cfdemCloud::voidFractionM() const
{
    return voidFractionModel_;
}

inline const averagingModel& cfdemCloud::averagingM() const
{
    return averagingModel_;
}

inline const clockModel& cfdemCloud::clockM() const
{
    return clockModel_;
}

inline const smoothingModel& cfdemCloud::smoothingM() const
{
    return smoothingModel_;
}

inline const meshMotionModel& cfdemCloud::meshMotionM() const
{
    return meshMotionModel_;
}

inline const wordList& cfdemCloud::liggghtsCommandModelList() const
{
    return liggghtsCommandModelList_;
}

inline autoPtr<liggghtsCommandModel>* cfdemCloud::liggghtsCommand() const
{
    return liggghtsCommand_;
}

#if defined(version24Dev)
    inline const turbulenceModel& cfdemCloud::turbulence() const
#elif defined(version21) || defined(version16ext)
    #ifdef compre
        inline const compressible::turbulenceModel& cfdemCloud::turbulence() const
    #else
        inline const incompressible::turbulenceModel& cfdemCloud::turbulence() const
    #endif
#elif defined(version15)
    inline const incompressible::RASModel& cfdemCloud::turbulence() const
#endif
{
    return turbulence_;
}

inline void cfdemCloud::makeSpecific(volScalarField& field)
{
    forAll(field,cellI)
        field[cellI] /= mesh_.V()[cellI];
}

inline void cfdemCloud::makeSpecific(volVectorField& field)
{
    forAll(field,cellI)
        field[cellI] /= mesh_.V()[cellI];
}

inline void cfdemCloud::scaleWithVcell(volScalarField& field)
{
    forAll(field,cellI)
        field[cellI] *= mesh_.V()[cellI];
}

inline void cfdemCloud::scaleWithVcell(volVectorField& field)
{
    forAll(field,cellI)
        field[cellI] *= mesh_.V()[cellI];
}

}
// ************************************************************************* //
