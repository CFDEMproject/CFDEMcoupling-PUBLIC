//============================================================//
// check model type
// referring to: Zhou et al. 2010,JFM
word modelType = particleCloud.modelType();

//Warning << "model type not being checked" << endl;
if (modelType=="Bfull"){
    Info << "solving volume averaged Navier Stokes equations of type Bfull\n"<< endl;

    // check if gradPForce is used
    bool found=false;
    forAll(particleCloud.forceModels(),i)
    {
        if(particleCloud.forceModels()[i]=="gradPForce")
            found=true;
    }
    if(!found)
        FatalError <<"gradPForce model not found! Consider selecting gradPForce model or switching to model type B.\n" << abort(FatalError);

    // check if viscForce is used
    found=false;
    forAll(particleCloud.forceModels(),i)
    {
        if(particleCloud.forceModels()[i]=="viscForce")
            found=true;
    }
    if(!found)
        FatalError <<"viscForce model not found! Consider selecting viscForce model or switching to model type B.\n" << abort(FatalError);

}else if(modelType=="B"){
    Info << "solving volume averaged Navier Stokes equations of type B\n"<< endl;

    // check if gradP and viscForce are used
    bool found=false;
    forAll(particleCloud.forceModels(),i)
    {
        if(particleCloud.forceModels()[i]=="gradPForce" || particleCloud.forceModels()[i]=="viscForce")
            found=true;
    }
    if(found)
        FatalError <<"do not use gradPForce and viscForce with model type B!\n" << abort(FatalError);

}else if (modelType=="A"){
    Info << "solving volume averaged Navier Stokes equations of type A"<< endl;

    // check if gradP is used
    bool found=false;
    forAll(particleCloud.forceModels(),i)
    {
        if(particleCloud.forceModels()[i]=="gradPForce")
            found=true;
    }
    if(!found)
        FatalError <<"gradPForce model not found! Consider selecting gradPForce model or switching to model type B.\n" << abort(FatalError);

    // check if viscForce is used
    found=false;
    forAll(particleCloud.forceModels(),i)
    {
        if(particleCloud.forceModels()[i]=="viscForce")
            found=true;
    }
    if(!found)
    {
        FatalError <<"viscForce model not found! Consider selecting viscForce model or switching to model type B.\n" << abort(FatalError);
    }
}else if (modelType=="none")
{
    Warning << "You chose model type -none- you might get erroneous results!" << endl;
}
else
{
    FatalError <<"no suitable model type specified:" << modelType << "\n" << abort(FatalError);
}
//============================================================//
