
if(particleCloud.verbose())
{
    Info <<"=============================================================================" << endl;
    Info << "totaldragforceEuler calculus" << endl;
    vector totaldragforceEuler(0,0,0);
    forAll(Ksl,cellI)
    {
        totaldragforceEuler += Ksl[cellI]*(Us[cellI]-U[cellI])/rho[cellI] * Ksl.mesh().V()[cellI];
    }
    Pout <<"totaldragforceEuler = "<< mag(totaldragforceEuler) << endl;
    Pout << "dv/dt =" << sum(fvc::ddt(voidfraction)) << endl;
    Info <<"=============================================================================" << endl;
}
