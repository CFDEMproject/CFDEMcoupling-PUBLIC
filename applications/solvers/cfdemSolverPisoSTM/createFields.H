    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading physical velocity field U" << endl;
    Info<< "Note: only if voidfraction at boundary is 1, U is superficial velocity!!!\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "createCFDEMFields.H"

//========================
// scalar field modelling
//========================
   /* Info<< "\nCreating dummy density field rho = 1\n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh//,
        //dimensionedScalar("0", dimensionSet(0, 0, -1, 1, 0), 273.15)
    );*/

    Info<< "\nCreating fluid-particle heat flux field\n" << endl;
    volScalarField Tsource
    (
        IOobject
        (
            "Tsource",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(0, 0, -1, 1, 0), 0.0)
    );

    /*IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar DT
    (
        transportProperties.lookup("DT")
    );*/

//========================

//#   include "createPhi.H"
#ifndef createPhi_H
#define createPhi_H
Info<< "Reading/calculating face flux field phi\n" << endl;
surfaceScalarField phi
 (
     IOobject
     (
        "phi",
        runTime.timeName(),
         mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
     linearInterpolate(U*voidfraction) & mesh.Sf()
 );
#endif

Info<< "Generating interstitial face flux field phiByVoidfraction\n" << endl;
surfaceScalarField phiByVoidfraction
(
 IOobject
 (
    "phiByVoidfraction",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
 ),
 linearInterpolate(U) & mesh.Sf()
);

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);


    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );
