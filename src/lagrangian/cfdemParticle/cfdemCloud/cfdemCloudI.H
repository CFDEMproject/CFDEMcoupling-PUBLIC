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

vector Foam::cfdemCloud::position(int index, bool ovr, int mp)
{
    vector pos;
    for(int i=0;i<3;i++)
        pos[i] = positions()[index][i];
    return pos;
}

vector Foam::cfdemCloud::velocity(int index, bool ovr, int mp)
{
    vector vel;
    for(int i=0;i<3;i++)
        vel[i] = velocities()[index][i];
    return vel;
}

inline void cfdemCloud::setCGTypeSpecific(int type, double cg)
{
    if(int(cgTypeSpecific_.size())!=type)
        FatalError  << "setCGTypeSpecific attempts to set at location "
                    << int(cgTypeSpecific_.size())
                    << ", but caller requests to set at "
                    << type
                    << ". This is fatal. "
                    << abort(FatalError);
    cgTypeSpecific_.push_back(cg);
    if(cg<cg_*0.99999 || cg>cg_*1.00001)
    {
        cgTypeSpecificDifferent_ = true;
        registerFieldsToDEM("type","scalar-atom",idPartTypes(),true);
    }
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

inline scalar cfdemCloud::cg(int index) const
{
    if (cgTypeSpecificDifferent_)
        return cgTypeSpecific_[particleType(index)];
    else if (cgParticleSpecific_)
        // idParticleCG constains the char. diameter -> convert to cg factor
        // on top, there may be a global cg -> apply that too
        return 2*fieldsToDEM[idRadius()][index][0]/fieldsToDEM[idParticleCG()][index][0]*cg_;
    else
        return cg_;
};

inline const bool& cfdemCloud::ignore() const
{
    return ignore_;
}

inline const bool& cfdemCloud::writeTimePassed() const
{
    return writeTimePassed_;
}

inline void cfdemCloud::resetWriteTimePassed()
{
    resetWriteTimePassed_=true;
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
    return fieldsToDEM[idPos()];
}

inline double ** cfdemCloud::velocities() const
{
    return fieldsToDEM[idVel()];
}

inline double ** cfdemCloud::fAccs() const
{
    return fieldsToDEM[idFacc()];
}

inline double ** cfdemCloud::impForces() const
{
    return impForces_;
}

inline double ** cfdemCloud::expForces() const
{
    return expForces_;
}

inline double ** cfdemCloud::radii() const
{
    return fieldsToDEM[idRadius()];
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
//        fieldsToDEM[idRadius()][0][i]=values[0][i];
    values=fieldsToDEM[idRadius()];
}

inline label Foam::cfdemCloud::cellID(int index, bool ovr, int mp) const
{
    return cellIDs_[index][0];
}

inline double ** cfdemCloud::cellIDs(bool ovr) const
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

inline scalar cfdemCloud::radius(int index, bool ovr)
{
    return fieldsToDEM[idRadius()][index][0];
}

inline double cfdemCloud::diameter(int index, bool ovr, int mp)
{
    return 2*fieldsToDEM[idRadius()][index][0];
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
            r=fieldsToDEM[idRadius()][index][0];
            Ntot+=2*r*r*r;
            Dtot+=r*r;
        }
        if(Ntot>SMALL) d32_ = Ntot/Dtot;
        else d32_ = 0.;
    }

    return d32_;
}

inline int cfdemCloud::numberOfObjects(bool ovr) const
{
    return numberOfParticles_;
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
    FatalError << "cfdemCloud::numberOfClumps() is not called correctly!" << abort(FatalError);
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
