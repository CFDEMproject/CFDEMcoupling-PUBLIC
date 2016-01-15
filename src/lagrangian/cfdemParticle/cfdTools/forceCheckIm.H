//calc. total implicit force
vector fImpTotal = sum(mesh.V()*Ksl.internalField()*(Us.internalField()-U.internalField()));
reduce(fImpTotal, sumOp<vector>());
Info << "TotalForceImp: " << fImpTotal << endl;
