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

#include "error.H"

#include "checkCouplingInterval.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(checkCouplingInterval, 0);

addToRunTimeSelectionTable
(
    forceModel,
    checkCouplingInterval,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
checkCouplingInterval::checkCouplingInterval
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    warnOnly_(propsDict_.lookupOrDefault<Switch>("warnOnly", true)),
    velocityFieldName_(propsDict_.lookupOrDefault<word>("velocityFieldName", "U")),
    U_(sm.mesh().lookupObject<volVectorField> (velocityFieldName_)),
    rhoP_(readScalar(propsDict_.lookup("rhoP"))),
    maxCFL_(propsDict_.lookupOrDefault<scalar>("maxCFL", 1.)),
    maxPCFL_(propsDict_.lookupOrDefault<scalar>("maxPCFL", 1.)),
    maxAccNr_(propsDict_.lookupOrDefault<scalar>("maxAccNr", 0.005)),
    UmaxExpected_(readScalar(propsDict_.lookup("UmaxExpected"))),
    minAllowedVcellByVparcel_(propsDict_.lookupOrDefault<scalar>("minAllowedVcellByVparcel", 3.)),
    nextRun_(0),
    timeInterval_(propsDict_.lookupOrDefault<scalar>("timeInterval", sm.dataExchangeM().couplingTime())),
    couplingStepInterval_(floor(timeInterval_/sm.dataExchangeM().couplingTime()))
{
    // init force sub model
    setForceSubModels(propsDict_);

    // read those switches defined above, if provided in dict
    forceSubM(0).readSwitches();

    /*if (probeIt_ && propsDict_.found("suppressProbe"))
        probeIt_=!Switch(propsDict_.lookup("suppressProbe"));
    if(probeIt_)
    {
        particleCloud_.probeM().initialize(typeName, typeName+".logDat");
        particleCloud_.probeM().scalarFields_.append("maxCFL");
        particleCloud_.probeM().scalarFields_.append("maxCFLexpected");
        particleCloud_.probeM().scalarFields_.append("minTauP");
        particleCloud_.probeM().scalarFields_.append("minVcellByVparcel");
        particleCloud_.probeM().scalarFields_.append("minStokes");
        particleCloud_.probeM().scalarFields_.append("maxStokes");
        particleCloud_.probeM().writeHeader();
    }*/

    Info << "running checkCouplingInterval every " << couplingStepInterval_ << " coupling steps,"
         << "which is every " << timeInterval_ << " seconds." << endl;

    particleCloud_.checkCG(true);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

checkCouplingInterval::~checkCouplingInterval()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
inline bool checkCouplingInterval::doCheck() const
{
    label couplingStep = particleCloud_.dataExchangeM().couplingStep();
    if(couplingStep >= nextRun_)
    {
        nextRun_ += couplingStepInterval_;
        return true;
    }
    return false;
}

void checkCouplingInterval::setForce() const
{
    #include "setupProbeModel.H"

    if(doCheck())
    {
        Info << "========================================================================" << endl;
        Info << "Check Coupling Interval" << endl;
        Info << "========================================================================" << endl;

        //==============================================================================================
        // calc the maximum CFL number
        scalar CFL = 0.0;
        scalar CFLexpected = 0.0;
        scalar ts = particleCloud_.mesh().time().deltaT().value();

        if (particleCloud_.mesh().nInternalFaces())
        {
            surfaceScalarField phi(particleCloud_.mesh().lookupObject<surfaceScalarField> ("phi"));
            scalarField sumPhi(fvc::surfaceSum(mag(phi))().internalField());

            CFL = 0.5*gMax(sumPhi/particleCloud_.mesh().V().field())*ts;

            scalar minVol = gMin(particleCloud_.mesh().V());
            scalar minLen = pow(minVol,1./3.);
            CFLexpected = UmaxExpected_ * ts / minLen;
        }

        if (CFL > maxCFL_)
        {
            if(warnOnly_)
                Warning << "CFL exceeds maximum allowed value:" << maxCFL_
                           << ", and reached the value:" << CFL
                           << "\ncheck your settings!" << endl;
            else
                FatalError << "CFL exceeds maximum allowed value:" << maxCFL_
                           << ", and reached the value:" << CFL
                           << "\nThe simulation is stopped, check your settings." << abort(FatalError);
        }
        else if (CFLexpected > maxCFL_)
            if(warnOnly_)
                Warning << "Expected CFL (based on UmaxExpected) exceeds maximum allowed value:" << maxCFL_
                           << ", and reached the value:" << CFLexpected
                           << "\ncheck your settings." << endl;
            else
                FatalError << "Expected CFL (based on UmaxExpected) exceeds maximum allowed value:" << maxCFL_
                           << ", and reached the value:" << CFLexpected
                           << "\nThe simulation is stopped, check your settings." << abort(FatalError);
        else
        {
            Info << "max. CFL = " << CFL << endl;
            Info << "max. expected CFL (based on UmaxExpected) = " << CFLexpected << endl;
        }
        //==============================================================================================


        //==============================================================================================
        // particle relaxation time tau,
        // Stokes nr,
        // cell/particle volume ratio

        const volScalarField& nufField = forceSubM(0).nuField();
        const volScalarField& rhoField = forceSubM(0).rhoField();

        // find min particle relaxation time
        scalar minTauP = 1e10;
        scalar minTauPAll = 1e10;
        scalar tauP = -1;
        scalar rad = -1;
        scalar scaledRad = -1.;

        // find min/max particle stokes Nr
        scalar St = -1;
        scalar minSt = 1e10;
        scalar minStAll = 1e10;
        scalar maxSt = -1;
        scalar maxStAll = -1;

        // max particle CFL Nr
        scalar maxVel = -1;
        scalar maxVelAll = -1;

        // find max parcel diameter
        scalar maxDparcel = -1;
        scalar maxDparcelAll = -1;

        // find min cell vol
        scalar minVcell = 1e10;
        scalar minVcellAll = 1e10;

        // get info on scaleDrag being used.
        // we do not call getProperty in constructor to be independent of order of forceModels
        scalar scaleDrag=particleCloud_.registryM().getProperty("scaleDrag");
        Info << "checkCouplingInterval considers scaleDrag = "<< scaleDrag << endl;

        for(int index = 0;index <  particleCloud_.numberOfParticles(); index++)
        {
            label cellI = particleCloud_.cellIDs()[index][0];

            if (cellI > -1) // particle Found
            {
                rad = particleCloud_.radius(index);
                forceSubM(0).scaleDia(rad,index);
                tauP = rhoP_*4*scaledRad*scaledRad/
                        (18 * nufField[cellI] * rhoField[cellI])/scaleDrag;
                minTauP = min(minTauP,tauP);

                maxDparcel = max(maxDparcel,2*particleCloud_.radius(index));

                minVcell = min(minVcell,particleCloud_.mesh().V()[cellI]);

                // StokesNr
                St = tauP*mag(U_[cellI]-particleCloud_.velocity(index))/(2*rad);
                minSt = min(minSt,St);
                maxSt = max(maxSt,St);

                maxVel = max(maxVel,mag(particleCloud_.velocity(index)));
            }
        }

        // allreduce results
        MPI_Allreduce(&minTauP, &minTauPAll, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&minSt, &minStAll, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&maxSt, &maxStAll, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&maxDparcel, &maxDparcelAll, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&minVcell, &minVcellAll, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&maxVel, &maxVelAll, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        // calc based on reduced results

        // calc max couplingTime/particle relaxation time ratio
        scalar DEMtime = particleCloud_.dataExchangeM().DEMts()
                        *particleCloud_.dataExchangeM().couplingInterval();
        scalar accNr = DEMtime/minTauPAll;

        // calc min cell/parcel volume ratio
        scalar maxVparcelAll = maxDparcelAll*maxDparcelAll*maxDparcelAll*3.1416/6.0;
        scalar minVcellByVparcel = minVcellAll/maxVparcelAll;

        //particle CFL number
        scalar pCFL = maxVelAll*particleCloud_.dataExchangeM().couplingTime()/pow(minVcell,1./3.);

        if(accNr > maxAccNr_)
        {
            if(warnOnly_)
                Warning << "The max. acceleration nr (coupling time / particle relaxation time) exceeds the maximum allowed value " << maxAccNr_
                           << ", and reached the value:" << accNr
                           << "\ncheck your settings." << endl;
            else
                FatalError << "The max. acceleration nr (coupling time / particle relaxation time) exceeds the maximum allowed value " << maxAccNr_
                           << ", and reached the value:" << accNr
                           << "\nThe simulation is stopped, check your settings." << abort(FatalError);
        }
        else
        {
            Info << "min. occurring particle relaxation time [s]: " << minTauPAll << endl;
            Info << "coupling interval [s]: " << DEMtime << endl;
            Info << "max. occurring acceleration nr (coupling time [s] / min. particle relaxation time [s]): " << accNr << endl;
        }

        if(minVcellByVparcel < minAllowedVcellByVparcel_)
        {
            if(warnOnly_)
                Warning << "The min. ratio of Vcell / Vparcel is below the minimum allowed value " << minAllowedVcellByVparcel_
                           << ", and reached the value:" << minVcellByVparcel
                           << "\ncheck your settings." << endl;
            else
                FatalError << "The min. ratio of Vcell / Vparcel is below the minimum allowed value " << minAllowedVcellByVparcel_
                           << ", and reached the value:" << minVcellByVparcel
                           << "\nThe simulation is stopped, check your settings." << abort(FatalError);
        }
        else
            Info << "min. ratio of Vcell / Vparcel is " << minVcellByVparcel << endl;

        Info << "min./max. Stokes Nr: " << minStAll << " / " << maxStAll << endl;

        if(pCFL > maxPCFL_)
        {
            if(warnOnly_)
                Warning << "The max. particle coupling CFL number is aboce the maximum allowed value " << maxPCFL_
                           << ", and reached the value:" << pCFL
                           << "\ncheck your settings." << endl;
            else
                FatalError << "The max. particle coupling CFL number is aboce the maximum allowed value " << maxPCFL_
                           << ", and reached the value:" << pCFL
                           << "\ncheck your settings." << abort(FatalError);
        }
        else
        Info << "max. particle coupling CFL Nr: " << pCFL << endl;

        Info << "========================================================================" << endl;

        /*//Set value fields and write the probe
        if(probeIt_)
        {
            #include "setupProbeModelfields.H"
            // Note: for other than ext one could use vValues.append(x)
            // instead of setSize
            sValues.setSize(sValues.size()+1, CFL);
            sValues.setSize(sValues.size()+1, CFLexpected);
            sValues.setSize(sValues.size()+1, minTauPAll);
            sValues.setSize(sValues.size()+1, minVcellByVparcel);
            sValues.setSize(sValues.size()+1, minStAll);
            sValues.setSize(sValues.size()+1, maxStAll);
            particleCloud_.probeM().writeProbe(0, sValues, vValues);
        }*/
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
