
//create PISO object
pisoControl piso(mesh);

//create lables
label nCorr             = piso.nCorrPISO();
label nNonOrthCorr      = piso.nNonOrthCorr();
