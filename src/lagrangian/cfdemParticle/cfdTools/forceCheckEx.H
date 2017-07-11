// NOTE:
// this file is typically included at a position of the solver so that
// it reports the momentum exchange based on f of the current time step

//calc. total source - this NOT a relevant explicit force!
dimensionedVector fExpSource = -1 * particleCloud.momCoupleM(particleCloud.registryM().getProperty("explicitCouple_index")).returnIntegralSourceField(); //negative because we want the force on the fluid
reduce(fExpSource, sumOp<dimensionedVector>()); //must reduce, since sum in function is not the global sum!

//calc. total explicit force
dimensionedVector fExpTotal = -gSum( 
                                    mesh.V()
                                #if defined(version40)
                                    * f.primitiveField()
                                #else
                                    * f.internalField()
                                #endif
                                    //
                                  ); //negative because we want the force on the fluid

fExpTotal.value() -= fExpSource.value();    //Subtract the explicit source, since this is NOT a coupling force!!

Info << "  TotalSourceExp: " << fExpSource.value() << endl;
Info << "  TotalForceExp:  " << fExpTotal.value()  << endl;
Info << "  These values are based on latest f and fExpSource.\n" << endl;
