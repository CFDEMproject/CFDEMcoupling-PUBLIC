// check if both implicit and explicitCouple model is declared - otherwise error!
// include this file after the cloud class is created

if(particleCloud.registryM().getProperty("implicitCouple_index")+particleCloud.registryM().getProperty("explicitCouple_index") < 1) FatalError << "You need to specify implicit and explicit momemtum coupling model for this solver." << abort(FatalError);
