//calc. total explicit force
vector fExpTotal =-sum(mesh.V()*f.internalField()); //negative because we want the force on the fluid
reduce(fExpTotal, sumOp<vector>());
Info << "TotalForceExp: " << fExpTotal << endl;
