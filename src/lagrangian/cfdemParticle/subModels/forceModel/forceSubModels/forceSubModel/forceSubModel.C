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
    nrDefaultSwitches_(11),                                          // !!!
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
    scaleDia_(dict_.lookupOrDefault<scalar>("scale",1.)),
    scaleDrag_(dict_.lookupOrDefault<scalar>("scaleDrag",1.)),
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
    if (scaleDia_ != 1)
        Info << "using scale = " << scaleDia_ << endl;
    else if (particleCloud_.cg() != 1)
    {
        scaleDia_=particleCloud_.cg();
        Info << "using scale from liggghts cg = " << scaleDia_ << endl;
    }

    // info about scaleDrag being used
    if (scaleDrag_ != 1.)
        Info << "using scaleDrag = " << scaleDrag_ << endl;

    if (scaleDrag_ < SMALL)
       FatalError<< "scaleDrag > 0 required" << abort(FatalError);

    particleCloud_.registryM().addProperty("scaleDrag",scaleDrag_);

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
    scalar Cd
) const
{
    // forces for CFD
    if(!switches_[1])// !treatForceDEM
    {
        if(switches_[0]) // treatForceExplicit
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
    if(switches_[2]) // implForceDEM
    {
        for(int j=0;j<3;j++)
            myForceM().fluidVel()[index][j]=Ufluid[j];

        myForceM().Cds()[index][0]=Cd;
    }
    else
    {
        for(int j=0;j<3;j++) 
            myForceM().DEMForces()[index][j] += dragTot[j];
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void forceSubModel::partToArrayAnisotropic
(
    const label& index,
    const vector& CdExtra,
    const vector& dragEx //this is the remaining explictit drag
) const
{

    if(switches_[2]) // implForceDEM
    {
        for(int iDir=0;iDir<3;iDir++)
        {
            myForceM().CdsExtra()[index][iDir]   = CdExtra[iDir];
            myForceM().DEMForces()[index][iDir] += dragEx[iDir];
        }
        if(switches_[3]) Pout << "*********forceSubModel::partToArrayAnisotropic: CdExtra: " 
                              << CdExtra 
                              << ", dragEx: " << dragEx << endl;
    }
    else
         FatalError<< "You are not using 'implForceDEM', but you want to update anisotropic drag data. This will not work." << abort(FatalError);   
    
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void forceSubModel::partToArrayAnisotropicTorque
(
    const label&        index,
    const vector&       CdTorque,
    const vector& torqueTotal  //this is the total torque
) const
{

    if(switches_[2]) // implForceDEM - also implies implicit torque!
    {
        for(int iDir=0;iDir<3;iDir++)
        {
            myForceM().CdsRotation()[index][iDir]   = CdTorque[iDir];
            myForceM().DEMTorques()[index][iDir]   += torqueTotal[iDir];
           
        }
        if(switches_[3]) Pout << "*********forceSubModel::partToArrayAnisotropic: CdTorque: " 
                              << CdTorque 
                              << ", torqueTotal: " << torqueTotal << endl;
    }
    else
         FatalError<< "You are not using 'implForceDEM', but you want to update anisotropic torque data. This will not work." << abort(FatalError);   
    
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
    if(particleCloud_.cgTypeSpecificDifferent)
        d /= particleCloud_.cg(particleCloud_.particleType(index));
    else
        d /= scaleDia_ / scaleDH_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void forceSubModel::scaleForce(vector& force, scalar& d, int index) const
{
    if(particleCloud_.cgTypeSpecificDifferent)
    {
        double cgCurr = particleCloud_.cg(particleCloud_.particleType(index));
        force *= cgCurr*cgCurr*cgCurr;
    }
    else
        force *= scaleDia_*scaleDia_*scaleDia_;
           
    force *= scaleDrag_;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void forceSubModel::scaleCoeff(scalar& coeff, scalar& d, int index) const
{
    if(particleCloud_.cgTypeSpecificDifferent)
    {
        double cgCurr = particleCloud_.cg(particleCloud_.particleType(index));
        coeff *= cgCurr*cgCurr*cgCurr;
    }
    else
        coeff *= scaleDia_*scaleDia_*scaleDia_;
        
    coeff *= scaleDrag_;
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
            }else Info << " not found in dict, using default. " ;
                
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

    if(switches_[2]) // implForceDEM=true
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

    // read extra variables
    dict_.readIfPresent("verboseDiskIntervall", verboseDiskIntervall_);

    // look for old nomenclature
    if (dict_.found("treatExplicit") || dict_.found("treatDEM") || dict_.found("implDEM"))
        FatalError<< "You are using an old nomenclature for force model settings, please have a look at the forceSubModel doc." << abort(FatalError);
        
    // look for old nomenclature
    if (dict_.found("verbose"))
        Warning<< "Please make sure you use the new nomenclature for verbose force model settings, please have a look at the forceSubModel doc." << endl;

    //if (dict_.found("interpolation"))
    //    FatalError<< "Please make sure you use the new nomenclature for interpolation in force model settings, please have a look at the forceSubModel doc." << endl;
}

const volScalarField& forceSubModel::nuField() const
{
    #ifdef compre
        nu_=particleCloud_.turbulence().mu() / rho_;
        return nu_;
    #else
        if(switches_[8]) // scalarViscosity=true
            return nu_;
        else
            return particleCloud_.turbulence().nu();
    #endif
}

const volScalarField& forceSubModel::muField() const
{
    #ifdef compre
        return particleCloud_.turbulence().mu();
    #else
        if(switches_[8]) // scalarViscosity=true
        {
            // usage of constant mu_ is still commented, as not tested
            // particleCloud_.turbulence().nu()*rho_ does not work properly
            FatalError<< "implementation not complete!" << abort(FatalError);
            //return mu_; // to be used with above code to set mu_ in readSwitches()

            return particleCloud_.turbulence().nu()*rho_;// for now just to have a return
        }else
            return particleCloud_.turbulence().nu()*rho_;
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace Foam

// ************************************************************************* //
