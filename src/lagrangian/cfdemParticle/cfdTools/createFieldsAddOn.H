    Info<< "\nConstructing explicit momentum exchange field f.\n" << endl;
    volVectorField f
    (
        IOobject
        (
            "f",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ, // zeroGradient here + zeroGradient in p can cause troubles --> force calculated by not reading the field 
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimensionSet(1, -2, -2, 0, 0), vector::zero)
    );
    f.storePrevIter();

    Info<< "Constructing coupling force error field fError.\n" << endl;
    volVectorField fError
    (
        IOobject
        (
            "fError",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimensionSet(1, -2, -2, 0, 0), vector::zero)
    );
