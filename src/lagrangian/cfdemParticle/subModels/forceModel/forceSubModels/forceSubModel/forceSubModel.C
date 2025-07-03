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
#include "forceSubModel.H"
#include "forceModel.H"
#include "mathExtra.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(forceSubModel, 0);

defineRunTimeSelectionTable(forceSubModel, dictionary);

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
forceSubModel::forceSubModel
(
    const dictionary& dict,
    cfdemCloud& sm,
    forceModel& fm
)
:
    dict_(dict),
    particleCloud_(sm),
    forceModel_(fm),
    nrDefaultSwitches_(34),                                          // !!!
    switchesNameList_(wordList(nrDefaultSwitches_)),
    switchesList_(List<Switch>(nrDefaultSwitches_)),
    switches_(List<Switch>(nrDefaultSwitches_)),
    nu_
    (
        IOobject
        (
            "scalarViscosity",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("nu0", dimensionSet(0, 2, -1, 0, 0), 1.)
    ),
    /*mu_
    (
        IOobject
        (
            "scalarViscosity",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("mu0", dimensionSet(1, -1, -1, 0, 0), 1.)
    ),*/
    divTau_
    (
        IOobject
        (
            "divTau",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedVector("divTau", dimensionSet(1, -2, -2, 0, 0), vector::zero)
    ),
    IBDragPerV_
    (
        IOobject
        (
            "IBDragPerV",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedVector("IBDragPerV", dimensionSet(1, -2, -2, 0, 0), vector::zero)
    ),
    densityFieldName_(dict_.lookupOrDefault<word>("densityFieldName","rho")),
    rho_(sm.mesh().lookupObject<volScalarField> (densityFieldName_)),
    verboseDiskIntervall_(1),
    verboseDiskCounter_(0),
    cg_(dict_.lookupOrDefault<scalar>("scale",1.)),
    scaleDrag_(dict_.lookupOrDefault<scalar>("scaleDrag",1.)),
    scaleDragLocal_(dict_.lookupOrDefault<scalar>("scaleDragLocal",1.)),
    scaleTorque_(dict_.lookupOrDefault<scalar>("scaleTorque",1.)),
    scaleDH_(dict_.lookupOrDefault<scalar>("scaleDH",1.))
{
    // init standard switch list
    int iCounter(0);
    switchesNameList_[iCounter]="treatForceExplicit"; iCounter++;   //0 - will treat force explicity (based on slip velocity)
    switchesNameList_[iCounter]="treatForceDEM";iCounter++;         //1 - will treat forces on DEM side only
    switchesNameList_[iCounter]="implForceDEM";iCounter++;          //2
    switchesNameList_[iCounter]="verbose";iCounter++;               //3
    switchesNameList_[iCounter]="interpolation";iCounter++;         //4
    switchesNameList_[iCounter]="useFilteredDragModel";iCounter++;  //5
    switchesNameList_[iCounter]="useParcelSizeDependentFilteredDrag";iCounter++;  //6
    switchesNameList_[iCounter]="implForceDEMaccumulated";iCounter++;             //7
    switchesNameList_[iCounter]="scalarViscosity";iCounter++;                     //8
    switchesNameList_[iCounter]="verboseToDisk";iCounter++;                       //9
    switchesNameList_[iCounter]="useCorrectedVoidage";iCounter++;                 //10
    switchesNameList_[iCounter]="multisphere";iCounter++;                         //11
    switchesNameList_[iCounter]="useTorque";iCounter++;                           //12
    switchesNameList_[iCounter]="anisotropicDrag";iCounter++;                     //13
    switchesNameList_[iCounter]="pullRotation";iCounter++;                        //14
    switchesNameList_[iCounter]="pullOrientation";iCounter++;                     //15
    switchesNameList_[iCounter]="pullShape";iCounter++;                           //16
    switchesNameList_[iCounter]="voidageFunctionDiFelice";iCounter++;             //17
    switchesNameList_[iCounter]="voidageFunctionRong";iCounter++;                 //18
    switchesNameList_[iCounter]="useUf";iCounter++;                               //19
    switchesNameList_[iCounter]="useFhydro";iCounter++;                           //20
    switchesNameList_[iCounter]="implTorqueDEM";iCounter++;                       //21
    switchesNameList_[iCounter]="useVisc";iCounter++;                             //22
    switchesNameList_[iCounter]="voidageFunctionTang";iCounter++;                 //23
    switchesNameList_[iCounter]="useMpData";iCounter++;                           //24
    switchesNameList_[iCounter]="superquadric";iCounter++;                        //25
    switchesNameList_[iCounter]="useQuaternion";iCounter++;                       //26
    switchesNameList_[iCounter]="pushTurbulence";iCounter++;                      //27
    switchesNameList_[iCounter]="particleSpecificCG";iCounter++;                  //28
    switchesNameList_[iCounter]="convex";iCounter++;                              //29
    switchesNameList_[iCounter]="pullType";iCounter++;                            //30
    switchesNameList_[iCounter]="pullDensity";iCounter++;                         //31
    switchesNameList_[iCounter]="pushConvectiveHeatFlux";iCounter++;              //32
    switchesNameList_[iCounter]="pullTemp";iCounter++;                            //33

    // should be done by default
    //for(int i=0;i<switchesList_.size();i++)
    //{
    //    switchesList_[i]=false;
    //    switches_[i]=false;
    //}

    // sanity check of what is defined above
    if(switchesNameList_.size() != nrDefaultSwitches_)
        FatalError<< "please check the nr of switches defined in forceSubModel class." << abort(FatalError);

    // info about scaleDia being used
    if (cg_ != 1)
        Info << "using scale = " << cg_ << endl;
    else if (particleCloud_.cg() != 1)
    {
        cg_=particleCloud_.cg();
        Info << "using scale from liggghts cg = " << cg_ << endl;
    }

    // info about scaleDrag being used
    if (scaleDrag_ != 1.)
        Info << "using scaleDrag = " << scaleDrag_ << endl;

    if (scaleDrag_ < SMALL)
       FatalError<< "scaleDrag > 0 required" << abort(FatalError);

    particleCloud_.registryM().addProperty("scaleDrag",scaleDrag_);

    // info about scaleDrag being used
    if (scaleDragLocal_ != 1.)
        Info << "using scaleDragLocal = " << scaleDragLocal_ << endl;

    if (scaleDragLocal_ < SMALL)
       FatalError<< "scaleDragLocal > 0 required" << abort(FatalError);

    // info about scaleTorque being used
    if (scaleTorque_ != 1.)
        Info << "using scaleTorque = " << scaleTorque_ << endl;

    if (scaleTorque_ < SMALL)
       FatalError<< "scaleTorque > 0 required" << abort(FatalError);

    // info about scaleDH being used
    if (scaleDH_ != 1)
        Info << "using scaleDH = " << scaleDH_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

forceSubModel::~forceSubModel()
{}

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //
void forceSubModel::partToArray
(
    const label& index,
    const vector& dragTot,
    const vector& dragEx,
    const vector& Ufluid,
    scalar Cd,
    const vector& CdExtra
) const
{
    // forces for CFD
    if(!treatForceDEM())
    {
        if(treatForceExplicit())
        {
            for(int j=0;j<3;j++)
                myForceM().expForces()[index][j] += dragTot[j];
        }
        else   //implicit treatment, taking explicit force contribution into account
        {
            for(int j=0;j<3;j++)
            {
                myForceM().impForces()[index][j] += dragTot[j] - dragEx[j];
                myForceM().expForces()[index][j] += dragEx[j];
            }
        }
    }

    // forces for DEM
    if(implForceDEM())
    {
        if(ms())
        {
            for(int j=0;j<3;j++)
                particleCloud_.fieldsToDEM[particleCloud_.idUfCM()][index][j] = Ufluid[j];

            if(anisotropicDrag())
            {
                for(int j=0;j<3;j++)
                    particleCloud_.fieldsToDEM[particleCloud_.idKslExtraCM()][index][j] = CdExtra[j];

                for(int j=0;j<3;j++)
                    particleCloud_.fieldsToDEM[particleCloud_.idDragExpCM()][index][j] += dragEx[j];
            }
            else
                particleCloud_.fieldsToDEM[particleCloud_.idKslCM()][index][0] = Cd;
        }
        else
        {
            for(int j=0;j<3;j++)
                particleCloud_.fieldsToDEM[particleCloud_.idUf()][index][j] = Ufluid[j];

            if(anisotropicDrag())
            {
                for(int j=0;j<3;j++)
                    particleCloud_.fieldsToDEM[particleCloud_.idKslExtra()][index][j] = CdExtra[j];

                for(int j=0;j<3;j++)
                    particleCloud_.fieldsToDEM[particleCloud_.idDragExp()][index][j] += dragEx[j];
            }
            else
                particleCloud_.fieldsToDEM[particleCloud_.idKsl()][index][0] = Cd;
        }
    }
    else
    {
        if(ms())
        {
            for(int j=0;j<3;j++)
                particleCloud_.fieldsToDEM[particleCloud_.idDragExpCM()][index][j] += dragTot[j];
        }
        else
        {
            for(int j=0;j<3;j++)
                particleCloud_.fieldsToDEM[particleCloud_.idDragExp()][index][j] += dragTot[j];
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void forceSubModel::partToArrayAnisotropicTorque
(
    const label&        index,
    const vector&       CdTorque,
    const vector& torqueTotal  //this is the total torque
) const
{

    if(useTorque())
    {
        if(implTorqueDEM())
        {
            if(ms())
            {
                for(int iDir=0;iDir<3;iDir++)
                {
                    particleCloud_.fieldsToDEM[particleCloud_.idKslRotationCM()][index][iDir] = CdTorque[iDir];
                    particleCloud_.fieldsToDEM[particleCloud_.idTorqueExpCM()][index][iDir] += torqueTotal[iDir];
                }
                /*if(verbose()) Pout << "*********forceSubModel::partToArrayAnisotropic: CdTorque: "
                                      << CdTorque
                                      << ", torqueTotal: " << torqueTotal << endl;*/
            }
            else
            {
                for(int iDir=0;iDir<3;iDir++)
                {
                    particleCloud_.fieldsToDEM[particleCloud_.idKslRotation()][index][iDir] = CdTorque[iDir];
                    particleCloud_.fieldsToDEM[particleCloud_.idTorqueExp()][index][iDir] += torqueTotal[iDir];
                }
            }
        }
        else
        {
            for(int iDir=0;iDir<3;iDir++)
                particleCloud_.fieldsToDEM[particleCloud_.idTorqueExpCM()][index][iDir] += torqueTotal[iDir];
        }

    }
    else
    {
        if(index==0) Info << "forceSubModel::partToArrayAnisotropicTorque, useTorque=false" << endl;
    }

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void forceSubModel::explicitCorr
(
    vector& dragImplicit,
    vector& dragExplicit,
    scalar& dragCoefficient,
    vector& Ufluid,
    const vector& Ucell,
    vector& Us,
    const vector& UsCell,
    bool verbose,
    label index
) const
{
    dragExplicit=vector::zero;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void forceSubModel::explicitCorrScalar(scalar& sourceKImplicit,
                                       scalar& sourceExplicit,
                                       scalar& areaTimesTransferCoefficient,
                                       scalar& fluidProperty,
                                       const   scalar& fluidPropertyCell,
                                       scalar& particleProperty,
                                       bool    verbose,
                                       label   index) const
{

    //everything is explicit, no verbose
    sourceExplicit  = areaTimesTransferCoefficient * (fluidProperty - particleProperty);
    sourceKImplicit = 0.0;

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void forceSubModel::update( label    particleI,
                            label    cellI,
                            scalar&  d,
                            scalar&  scalToUpdate1,
                            scalar&  scalToUpdate2,
                            bool     verbose
                          ) const
{
    //no action
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void forceSubModel::update( label    particleI,
                            label    cellI,
                            scalar&  d,
                            vector&  vecToUpdate1,
                            vector&  vecToUpdate2,
                            scalar&  scalToUpdate1,
                            scalar&  scalToUpdate2,
                            bool     verbose
                          ) const
{
    //no action
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void forceSubModel::scaleDia(scalar& d, int index) const
{
    if (particleCG() || particleCloud_.cgTypeSpecificDifferent_)
        d /= particleCloud_.cg(index) / scaleDH_;
    else
        d /= cg_ / scaleDH_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void forceSubModel::scaleForce(vector& force, scalar& d, int index) const
{
    if (particleCG() || particleCloud_.cgTypeSpecificDifferent_)
        force *= getCG(index)*getCG(index)*getCG(index);
    else
        force *= cg_*cg_*cg_;

    force *= scaleDrag_;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void forceSubModel::scaleCoeff(scalar& coeff, scalar& d, scalar& Rep, int index) const
{
    if (particleCG() || particleCloud_.cgTypeSpecificDifferent_)
        coeff *= getCG(index)*getCG(index)*getCG(index);
    else
        coeff *= cg_*cg_*cg_;

    coeff *= scaleDrag_;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void forceSubModel::scaleForceLocal(vector& force) const
{
    force *= scaleDragLocal_;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void forceSubModel::scaleTorque(vector& torque) const
{
    torque *= scaleTorque_;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
const scalar forceSubModel::scaleDrag(int index) const
{
    return particleCloud_.cg(index);
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void forceSubModel::explicitLimit
(
    vector& dragImplicit,
    vector& dragExplicit,
    scalar& d
) const
{
    dragExplicit=vector::zero;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void forceSubModel::readSwitches() const
{
    bool first(true);
    forAll(switchesNameList_,i)
    {
        if(switchesList_[i] > 0+SMALL) //check if switch is read from dict
        {
            if(first)
            {
                Info << " reading switches:" << endl;
                first = false;
            }
            Info << "  looking for " << switchesNameList_[i] << " ... ";
            if (dict_.found(switchesNameList_[i]))
            {
                Info << " found in dict. " ;
                switches_[i]=Switch(dict_.lookup(switchesNameList_[i]));
            }else
            {
                Info << " not found in dict, using default. " ;
                if(i==2) Warning << "\n You are using the default value for 'implForceDEM' - beware that in CFDEMcoupling versions newer than 3.8.0 this default value has changed from (previously) false to (now) true!" << endl;
            }

            // user set treatForceExplicit to true && no explicitCouple model is defined && treatForceDEM=false(i.e. it will go to f) && modelType != "none"
            if(i==0 && switches_[0] > 0+SMALL && particleCloud_.registryM().getProperty("explicitCouple_index") < 0 && switches_[1] < 0+SMALL && particleCloud_.modelType() != "none")
                FatalError <<  "You are using treatForceExplicit=true here, this requres having an explicit momentum couple model!" << abort(FatalError);
            else
                Info << switchesNameList_[i] << " = " << switches_[i] << endl;
        }
    }
    Info << endl;

    /*// sanity check
    if(switchesList_[0] < 0+SMALL) //check if switch is not read from dict
    {
        // (treatForceExplicit is set to true) && no explicitCouple model is defined && treatForceDEM=false
        if(switches_[0] > 0+SMALL && particleCloud_.registryM().getProperty("explicitCouple_index") < 0 && switches_[1] < 0+SMALL)
            FatalError <<  "treatForceExplicit = true. This requres having an explicit momentum couple model!" << abort(FatalError);
    }


    // debug info
    Info << "other switches for forceSubModel which are auot-set:" << myType() << endl;
    forAll(switchesNameList_,i)
    {
        if(switchesList_[i] < 0+SMALL) //check if switch is NOT read from dict but exists
        {
            Info << "\t" << switchesNameList_[i] << " = " << switches_[i] << endl;
        }
    }
    Info << endl;*/

    if(implForceDEM()) // implForceDEM=true
    {
        // communicate implForceDEM to particleCloud
        particleCloud_.impDEMdrag_=true;

        // do sanity check
        // This can work if the accumulator is used, but is explicitely applied on the CFD side
        // Sanity check is therefore not necessary here
        /*
        if(switches_[0]) // treatExplicit=true
        {
            FatalError << "Please check your settings, treatExplicit together with implForceDEM does not work!."
                       << abort(FatalError);
        }
        */
    }

    if(switches_[7]) // implForceDEMaccumulated=true
    {
        // sanity check for implForceDEMaccumulated
        if(!switches_[2]) //implForceDEM=false
        {
            Warning<< "please check your settings, implForceDEMaccumulated without implForceDEM does not work! (using implForceDEMaccumulated=false)" << endl;
            switches_[3]=false;
        }else
        {
            particleCloud_.impDEMdragAcc_=true;
        }
    }

    if(switches_[8]) // scalarViscosity=true
    {
        Info << "Using a constant viscosity for this force model." << endl;
        dimensionedScalar  nu0_("nu", dimensionSet(0, 2, -1, 0, 0), dict_.lookup("nu"));
        nu_=volScalarField
        (
            IOobject
            (
                "scalarViscosity",
                particleCloud_.mesh().time().timeName(),
                particleCloud_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            particleCloud_.mesh(),
            nu0_
        );

        /*if (dict_.found("mu"))
        {
            Info << "Using a constant viscosity for this force model." << endl;
            dimensionedScalar  mu0_("mu", dimensionSet(1, -1, -1, 0, 0), dict_.lookup("mu"));
            mu_=volScalarField
            (
                IOobject
                (
                    "scalarViscosity",
                    particleCloud_.mesh().time().timeName(),
                    particleCloud_.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                particleCloud_.mesh(),
                mu0_
            );
        }*/
    }

    if(anisotropicDrag() && !implForceDEM())
        FatalError<< "You have set 'implForceDEM' to false, but 'anisotropicDrag' to true."<< abort   (FatalError);

    // read extra variables
    dict_.readIfPresent("verboseDiskIntervall", verboseDiskIntervall_);

    // look for old nomenclature
    if (dict_.found("treatExplicit") || dict_.found("treatDEM") || dict_.found("implDEM"))
        FatalError<< "You are using an old nomenclature for force model settings, please have a look at the forceSubModel doc." << abort(FatalError);

    // look for old nomenclature
    if (dict_.found("verbose"))
        Warning<< "Please make sure you use the new nomenclature for verbose force model settings, please have a look at the forceSubModel doc." << endl;
}

void forceSubModel::setupCommunication() const
{
    particleCloud_.registerFieldsToDEM("radius","scalar-atom",particleCloud_.idRadius(),true);
    particleCloud_.registerFieldsToDEM("x","vector-atom",particleCloud_.idPos(),true);
    particleCloud_.registerFieldsToDEM("v","vector-atom",particleCloud_.idVel(),true);

    if(particleCloud_.impDEMdragAcc())
        particleCloud_.registerFieldsToDEM("dragAcc","vector-atom",particleCloud_.idFacc(),true);

    if(ms())
    {
        if(implForceDEM())
        {
            if(anisotropicDrag())
                particleCloud_.registerFieldsToDEM("Ksl_Extra_cm","vector-multisphere",particleCloud_.idKslExtraCM());
            else
                particleCloud_.registerFieldsToDEM("Ksl_cm","scalar-multisphere",particleCloud_.idKslCM());

            particleCloud_.registerFieldsToDEM("uf_cm","vector-multisphere",particleCloud_.idUfCM());
        }
        else
        {
            particleCloud_.registerFieldsToDEM("dragforce_cm","vector-multisphere",particleCloud_.idDragExpCM());
        }
        if(useTorque())
        {
            if(implTorqueDEM())
            {
                particleCloud_.registerFieldsToDEM("KslRotation","vector-multisphere",particleCloud_.idKslRotationCM()); // TODO: this should have other name on DEM side (is ambigous with MS version)
            }
            particleCloud_.registerFieldsToDEM("hdtorque_cm","vector-multisphere",particleCloud_.idTorqueExpCM());
        }
    }
    else
    {
        if(implForceDEM())
        {
            if(anisotropicDrag())
                particleCloud_.registerFieldsToDEM("KslExtra","vector-atom",particleCloud_.idKslExtra());
            else
                particleCloud_.registerFieldsToDEM("Ksl","scalar-atom",particleCloud_.idKsl());

            particleCloud_.registerFieldsToDEM("uf","vector-atom",particleCloud_.idUf());
        }
        else
        {
            particleCloud_.registerFieldsToDEM("dragforce","vector-atom",particleCloud_.idDragExp());
        }
        if(useTorque())
        {
            if(implTorqueDEM())
            {
                particleCloud_.registerFieldsToDEM("KslRotation","vector-atom",particleCloud_.idKslRotation());
            }
            particleCloud_.registerFieldsToDEM("hdtorque","vector-atom",particleCloud_.idTorqueExp());
        }
    }


    if(pullRotation())
        particleCloud_.registerFieldsToDEM("omega","vector-atom",particleCloud_.idPullRotation(),true);

    if(pullOrientation())
    {
        if(ms())
        {
            particleCloud_.registerFieldsToDEM("ex_space","vector-multisphere",particleCloud_.idPullOrientation(),true);
            particleCloud_.registerFieldsToDEM("ey_space","vector-multisphere",particleCloud_.idPullOrientation1(),true);
        }
        else
        {
            particleCloud_.registerFieldsToDEM("ex","vector-atom",particleCloud_.idPullOrientation(),true);
        }
    }
    if(pullShape())
        particleCloud_.registerFieldsToDEM("shape","vector-atom",particleCloud_.idPullShape(),true);
    if(useUf())
        particleCloud_.registerFieldsToDEM("uf","vector-atom",particleCloud_.idUf());
    if(useFhydro())
        particleCloud_.registerFieldsToDEM("Fhydro","vector-atom",particleCloud_.idFhydro(),true);
    if(useVisc())
        particleCloud_.registerFieldsToDEM("muf","scalar-atom",particleCloud_.idVisc());

    if(sq())
    {
        particleCloud_.registerFieldsToDEM("area","scalar-atom",particleCloud_.idArea(),true);
        particleCloud_.registerFieldsToDEM("volume","scalar-atom",particleCloud_.idVol(),true);
        particleCloud_.registerFieldsToDEM("blockiness","vector2D-atom",particleCloud_.idBlockiness(),true);
    }
    if(convex())
    {
        particleCloud_.registerFieldsToDEM("rmass","scalar-atom",particleCloud_.idMass(),true);
        particleCloud_.registerFieldsToDEM("density","scalar-atom",particleCloud_.idDensity(),true);
        particleCloud_.registerFieldsToDEM("shapetype","scalar-atom",particleCloud_.idType(),true);
    }
    if(pullType())
    {
        word nameForType;
        if(particleCloud_.shapeTypeName()=="convex")
            nameForType=word("shapetype");
        else
        {
            nameForType=word("type");
            Warning <<"\n  Using (material) type to distinguish types (e.g. scaleVol).\n"
                    <<"  You might use separate material types if you want to scale them separately." << endl;
        }
        //TODO use e.g. group to identify different sq types and use groups for SQ templates
        //if(particleCloud_.shapeTypeName()=="superquadric")
        particleCloud_.registerFieldsToDEM(nameForType,"scalar-atom",particleCloud_.idType(),true);
    }
    if(pullDensity())
        particleCloud_.registerFieldsToDEM("density","scalar-atom",particleCloud_.idDensity(),true);
    if(pushConvectiveHeatFlux())
        particleCloud_.registerFieldsToDEM("convectiveHeatFlux","scalar-atom",particleCloud_.idConvectiveHeatFlux());
    if(pullTemp())
        particleCloud_.registerFieldsToDEM("Temp","scalar-atom",particleCloud_.idTemp(),true);
    if(useQuat())
        particleCloud_.registerFieldsToDEM("quaternion","quaternion-atom",particleCloud_.idQuat(),true);

    if(pushTurbulence())
    {
        particleCloud_.registerFieldsToDEM("k","scalar-atom",particleCloud_.idK());
        particleCloud_.registerFieldsToDEM("epsilon","scalar-atom",particleCloud_.idEpsilon());
    }

    if (particleCG())
    {
        particleCloud_.registerFieldsToDEM("dSauter", "scalar-atom", particleCloud_.idParticleCG(), true);
        particleCloud_.cgParticleSpecific_ = true;
    }
}

const volScalarField& forceSubModel::nuField() const
{
    #ifdef compre
        nu_=particleCloud_.turbulence_.mu() / rho_;
        return nu_;
    #else
        if(switches_[8]) // scalarViscosity=true
            return nu_;
        else
            return particleCloud_.turbulence_.nu();
    #endif
}

const volScalarField& forceSubModel::muField() const
{
    #ifdef compre
        return particleCloud_.turbulence_.mu();
    #else
        if(switches_[8]) // scalarViscosity=true
        {
            // usage of constant mu_ is still commented, as not tested
            // particleCloud_.turbulence_.nu()*rho_ does not work properly
            FatalError<< "implementation not complete!" << abort(FatalError);
            //return mu_; // to be used with above code to set mu_ in readSwitches()

            return particleCloud_.turbulence_.nu()*rho_;// for now just to have a return
        }else
            return particleCloud_.turbulence_.nu()*rho_;
    #endif
}

const volScalarField& forceSubModel::rhoField() const
{
    return rho_;
}

const volVectorField& forceSubModel::divTauField(const volVectorField& U) const
{
    // calc div(Tau)
    #ifdef compre
        const volScalarField& mu_ = muField();
        divTau_ = -fvc::laplacian(mu_, U) - fvc::div(mu_*dev(fvc::grad(U)().T()));
        return divTau_;
    #else
        const volScalarField& nu_ = nuField();
        const volScalarField& rho_ = rhoField();
        divTau_ = -fvc::laplacian(nu_*rho_, U)- fvc::div(nu_*rho_*dev(fvc::grad(U)().T()));
        return divTau_;
    #endif
}

const volVectorField& forceSubModel::IBDragPerV(const volVectorField& U,const volScalarField& p) const
{
    #ifdef compre
        IBDragPerV_ = muField()*fvc::laplacian(U)-fvc::grad(p);
    #else
        IBDragPerV_ = rhoField()*(nuField()*fvc::laplacian(U)-fvc::grad(p));
    #endif
    return IBDragPerV_;
}

void forceSubModel::calcXi(const scalar& ds, const scalar& voidfraction, const scalar& magUr, const scalar& nuf, scalar& Xi ) const
{
    if(!(voidageFunctionDiFelice() || voidageFunctionRong() || voidageFunctionTang()))
    {
        Xi = 2;
    } else if(!(voidageFunctionDiFelice() ^ voidageFunctionRong() ^ voidageFunctionTang() ))
    { //triple XOR, "There can be only one!"
        FatalError<< "Only one drag correction function permitted, please select either voidageFunctionDiFelice OR voidageFunctionRong!" << abort(FatalError);
    } else if (voidageFunctionDiFelice())
    {
        // calc DiFelice drag correction

        // calc particle Re number
        scalar Rep = ds*voidfraction*magUr/(nuf+SMALL);

        // calc Xi
        if(Rep < SMALL)
            Xi = 3.7;
        else
            Xi = 3.7 - 0.65 * exp(-sqr(1.5-log10(Rep))/2);

    } else if (voidageFunctionRong())
    {
        // calculate Rong drag correction

        scalar Rep = ds*voidfraction*magUr/(nuf+SMALL);
        if(Rep < SMALL)
            Xi = 2.65 * (voidfraction + 1.0);
        else
            Xi = 2.65 * (voidfraction + 1.0) - (5.3-3.5*voidfraction)*voidfraction*voidfraction*exp(-sqr(1.5-Foam::log10(Rep))/2.0);
    } else if (voidageFunctionTang())
    {
        // calculate Tang drag correction
        scalar Rep = ds*voidfraction*magUr/(nuf+SMALL);
        if(Rep < SMALL || 1 - voidfraction < SMALL)
            Xi = 2;
        else
            Xi = 2 -
            Foam::log10(((1.5*sqrt(1-voidfraction)+1)*voidfraction*voidfraction - (10*(voidfraction-1))/(voidfraction*voidfraction) +
             (0.0644*pow(voidfraction,-4)+0.169*voidfraction)*pow(Rep,0.657) - (0.00456*Rep)*pow(voidfraction,-4) + 0.11*(voidfraction-2)*(voidfraction-1)*Rep)
             / (0.2334*pow(Rep,0.657) - 0.00456*Rep + 1))
            / Foam::log10(voidfraction);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace Foam

// ************************************************************************* //
