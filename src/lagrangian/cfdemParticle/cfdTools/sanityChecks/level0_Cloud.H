//======================================================================//
//Check if user attempts to change fluid time step
if( mesh_.time().controlDict().lookupOrDefault<Switch>("adjustTimeStep", false) && !allowAdjustTimeStep_ )
{
    FatalError << "cfdemCloud:: you want to adjustTimeStep in controlDict. This is not allowed in this version of CFDEM."
               << abort(FatalError);
}

//======================================================================//
//Check if mesh is empty (decomposition error)
if(mesh_.cells().size() < 1)
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    if(nprocs > 2)  FatalError << endl << "cfdemCloud:: local mesh has zero cells. Please check the mesh and the decomposition!" << abort(FatalError);

    Pout << "WARNING: cfdemCloud:: local mesh has zero cells. Please check the mesh and the decomposition!" << endl;
}

