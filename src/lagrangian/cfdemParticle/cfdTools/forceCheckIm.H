//calc. total implicit force
volVectorField fImp(Ksl*(Us-U));
particleCloud.scaleWithVcell(fImp);
dimensionedVector fImpTotal = gSum(fImp);
//                                      mesh.V()
//                                    * Ksl.internalField()
//                                    * ( Us.internalField() - U.internalField() )
//                                  );
Info << "TotalForceImp:  " << fImpTotal.value() << endl;
