 //set probeModel parameters for this force model
 //CANNOT be used for scalarGeneralExchange force model, since this force model handles multiple species
 if(probeIt_)
     particleCloud_.probeM().setOutputFile(typeName+".logDat");

