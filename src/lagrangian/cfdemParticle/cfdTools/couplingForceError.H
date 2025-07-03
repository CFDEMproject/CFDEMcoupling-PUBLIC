//Compute error due to  implicit part of the 
//coupling force acting on particle cloud (using fluid velocity AFTER solving NS)
#if defined(version40)
     fError.ref() -= U.ref() * dimensionedScalar("1", dimensionSet(1, -3, -1, 0, 0), 1.0);
     fError.ref() *= Ksl.oldTime().ref() * dimensionedScalar("1", dimensionSet(-1, 3, 1, 0, 0), 1.0); 
#else
    fError.internalField() -= U.internalField();
    fError.internalField() *= Ksl.oldTime().internalField(); 
#endif

//Remove the old error of coupling force from last fluid time step
//to f to conserve the overall momentum
//Minus here, because -f in UEqn
if(particleCloud.expCorrDeltaUError())
    f += fError;

//Integrate error
vector fErrorSum(0.0,0.0,0.0);
#if defined(version40)
  fErrorSum = sum(mesh.V()*fError.primitiveField());
#else
  fErrorSum = sum(mesh.V()*fError.internalField());
#endif
reduce(fErrorSum, sumOp<vector>());

Info << "TotalError(dU): " <<  fErrorSum;
if(particleCloud.verbose()) // in verbose mode we have data for fImpTotal and fExpTotal to normalize fErrorSum
{
    Info                  << " = " 
                          <<  mag(fErrorSum) / (mag(fImpTotal.value()+fExpTotal.value())+SMALL)*100.0 << "% of TotalForceImp";
}

if(particleCloud.expCorrDeltaUError())
    Info << " (corrected)     " << endl;
else
    Info << " (uncorrected)   " << endl;

//Save U BEFORE solving NS for computing error later
#if defined(version40)
  fError.ref() = U.ref() * dimensionedScalar("1", dimensionSet(1, -3, -1, 0, 0), 1.0);
#else
  fError.internalField() = U.internalField();
#endif

