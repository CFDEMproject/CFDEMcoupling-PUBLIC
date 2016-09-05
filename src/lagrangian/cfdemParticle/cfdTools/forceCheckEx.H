//calc. total source - this NOT a relevant explicit force!
dimensionedVector fExpSource = -1 * particleCloud.momCoupleM(1).returnIntegralSourceField(); //negative because we want the force on the fluid
reduce(fExpSource, sumOp<dimensionedVector>()); //must reduce, since sum in function is not the global sum!

//calc. total explicit force
dimensionedVector fExpTotal = -gSum( 
                                    mesh.V()
                                   * f.internalField()
                                  ); //negative because we want the force on the fluid

fExpTotal.value() -= fExpSource.value();    //Subtract the explicit source, since this is NOT a coupling force!!

Info << "TotalSourceExp: " << fExpSource.value() << endl;
Info << "TotalForceExp:  " << fExpTotal.value()  << endl;
