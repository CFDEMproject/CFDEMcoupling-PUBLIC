// NOTE:
// this file is typically included at a position of the solver so that
// it reports the momentum exchange based on Ksl and Us of the current time step
// but U is the solution of the previous time step!!!

dimensionedVector fImpTotal;

if(particleCloud.verbose())
{
    Info << "\nSolver level total Eulerian momentum exchange:"<< endl;

    //calc. total implicit force
    volVectorField fImp(Ksl*(Us-U));
    particleCloud.scaleWithVcell(fImp);
    fImpTotal = gSum(fImp);
    Info << "  TotalForceImp [N]:  " << fImpTotal.value() << endl;
    Info << "  Warning, these values are based on latest Ksl and Us but prev. iteration U!\n" << endl;
}
