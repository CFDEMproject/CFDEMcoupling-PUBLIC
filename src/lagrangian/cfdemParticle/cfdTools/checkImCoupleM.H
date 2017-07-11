// check if NO explicitCouple model is declared - otherwise error!
// include this file before the cloud class is created

{
    // read dict
    IOdictionary    dict(
                            IOobject
                            (
                                "couplingProperties",
                                mesh.time().constant(),
                                mesh,
                                IOobject::MUST_READ,
                                IOobject::NO_WRITE
                            )
                        );

    // read dict
    wordList a(dict.lookup("momCoupleModels"));

    // look for explicit and throw error if found
    for(int i=0; i<a.size(); i++)
        if(a[i] == "explicitCouple")
            FatalError << "You specified an explicit momemtum coupling model, which is not suitable for this solver." << abort(FatalError);
}
