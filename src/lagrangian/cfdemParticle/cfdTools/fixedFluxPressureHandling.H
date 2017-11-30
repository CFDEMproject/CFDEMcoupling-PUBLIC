#ifdef version50
    if (modelType=="A")
    {
        volScalarField rUsed = rUA*voidfraction;
        constrainPressure(p, U, phi, rUsed);
    }
    else constrainPressure(p, U, phi, rUA);
#endif
#ifndef versionExt32
#ifndef version40
    if (modelType=="A")
    {
        setSnGrad<fixedFluxPressureFvPatchScalarField>
        (
            #ifdef versionv1612plus
                p.boundaryFieldRef(),
            #else
                p.boundaryField(),
            #endif
            (
                phi.boundaryField()
              - (mesh.Sf().boundaryField() & U.boundaryField())
            )/(mesh.magSf().boundaryField()*rUAf.boundaryField()*voidfractionf.boundaryField())
        );
    }else
    {
        setSnGrad<fixedFluxPressureFvPatchScalarField>
        (
            #ifdef versionv1612plus
                p.boundaryFieldRef(),
            #else
                p.boundaryField(),
            #endif
            (
                phi.boundaryField()
              - (mesh.Sf().boundaryField() & U.boundaryField())
            )/(mesh.magSf().boundaryField()*rUAf.boundaryField())
        );
    }
#endif
#endif
