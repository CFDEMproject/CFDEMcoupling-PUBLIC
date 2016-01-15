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

//#include "mpi.h"

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
    U_(sm.mesh().lookupObject<volVectorField> ("U")),
    rhoP_(readScalar(propsDict_.lookup("rhoP"))),
    maxCFL_(50)
{
    if (propsDict_.found("maxCFL"))
        maxCFL_=readScalar(propsDict_.lookup("maxCFL"));

    // init force sub model
    setForceSubModels(propsDict_);

    // read those switches defined above, if provided in dict
    forceSubM(0).readSwitches();

    if (probeIt_ && propsDict_.found("suppressProbe"))
        probeIt_=!Switch(propsDict_.lookup("suppressProbe"));
    if(probeIt_)
    {
        particleCloud_.probeM().initialize(typeName, "checkCouplingInterval.logDat");
        particleCloud_.probeM().scalarFields_.append("minTauP");
        particleCloud_.probeM().scalarFields_.append("minVcellByVparcel");
        particleCloud_.probeM().scalarFields_.append("minStokes");
        particleCloud_.probeM().scalarFields_.append("maxStokes");
        particleCloud_.probeM().writeHeader();
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

checkCouplingInterval::~checkCouplingInterval()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void checkCouplingInterval::setForce() const
{
    #include "setupProbeModel.H"

    if(particleCloud_.mesh().time().write())
    {

        // calc the maximum CFL number
        scalar CFL = 0.0;

        if (particleCloud_.mesh().nInternalFaces())
        {
            surfaceScalarField phi(particleCloud_.mesh().lookupObject<surfaceScalarField> ("phi"));
            scalarField sumPhi(fvc::surfaceSum(mag(phi))().internalField());

            CFL = 0.5*gMax(sumPhi/particleCloud_.mesh().V().field())*particleCloud_.mesh().time().deltaT().value();
        }
        if (CFL > maxCFL_)
            FatalError << "CFL exceeds maximum allowed value:" << maxCFL_
                       << ", and reached the value:" << CFL
                       << "\nThe simulation is stopped, check your settings." << abort(FatalError);


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

        // find max parcel diameter
        scalar maxDparcel = -1;
        scalar maxDparcelAll = -1;

        // find min cell vol
        scalar minVcell = 1e10;
        scalar minVcellAll = 1e10;

        for(int index = 0;index <  particleCloud_.numberOfParticles(); index++)
        {
            label cellI = particleCloud_.cellIDs()[index][0];

            if (cellI > -1) // particle Found
            {
                rad = particleCloud_.radius(index);
                scaledRad = rad/particleCloud_.cg();
                tauP = rhoP_*4*scaledRad*scaledRad/
                        (18 * nufField[cellI] * rhoField[cellI]);
                minTauP = min(minTauP,tauP);

                maxDparcel = max(maxDparcel,2*particleCloud_.radius(index));

                minVcell = min(minVcell,particleCloud_.mesh().V()[cellI]);

                // StokesNr
                St = tauP*mag(U_[cellI]-particleCloud_.velocity(index))/(2*rad);
                minSt = min(minSt,St);
                maxSt = max(maxSt,St);
            }
        }

        // allreduce results
        MPI_Allreduce(&minTauP, &minTauPAll, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&minSt, &minStAll, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&maxSt, &maxStAll, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&maxDparcel, &maxDparcelAll, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&minVcell, &minVcellAll, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

        // calc based on reduced results

        // calc max couplingTime/particle relaxation time ratio
        scalar DEMtime = particleCloud_.dataExchangeM().DEMts()
                        *particleCloud_.dataExchangeM().couplingInterval();
        scalar accNr = DEMtime/minTauPAll;

        // calc min cell/parcel volume ratio
        scalar maxVparcelAll = maxDparcelAll*maxDparcelAll*maxDparcelAll*3.1416/6.0;
        scalar minVcellByVparcel = minVcellAll/maxVparcelAll;


        Info << "min. occurring particle relaxation time [s]: " << minTauPAll << endl;
        Info << "coupling interval [s]: " << DEMtime << endl;
        Info << "max. occurring acceleration nr: " << accNr << endl;
        if(accNr > 0.1) Warning << "you should use a smaller coupling interval!" << endl;

        Info << "min. occurring cell/parcel volume ratio: " << minVcellByVparcel << endl;
        if(minVcellByVparcel < 3) Warning << "you should use bigger cells!" << endl;

        Info << "min. Stokes Nr: " << minStAll << endl;
        Info << "max. Stokes Nr: " << maxStAll << endl;

        //Set value fields and write the probe
        if(probeIt_)
        {
            #include "setupProbeModelfields.H"
            // Note: for other than ext one could use vValues.append(x)
            // instead of setSize
            sValues.setSize(sValues.size()+1, minTauPAll);
            sValues.setSize(sValues.size()+1, minVcellByVparcel);
            sValues.setSize(sValues.size()+1, minStAll);
            sValues.setSize(sValues.size()+1, maxStAll);
            particleCloud_.probeM().writeProbe(0, sValues, vValues);
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
