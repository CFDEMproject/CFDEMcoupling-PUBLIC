#!/bin/bash

#===================================================================#
# allrun script for testcase as part of test routine 
# run ErgunTestMPI_cgs (working in CGS unit system)
# Christoph Goniva - March 2013
#===================================================================#

#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"

# check if mesh was built
if [ -f "$casePath/CFD/constant/polyMesh/points" ]; then
    echo "mesh was built before - using old mesh"
else
    echo "mesh needs to be built"
    cd $casePath/CFD
    blockMesh
fi

if [ -f "$casePath/DEM/post/restart/liggghts.restart" ];  then
    echo "LIGGGHTS init was run before - using existing restart file"
else
    #- run DEM
    $casePath/parDEMrun.sh
fi

#- run parallel CFD-DEM in new terminal
#gnome-terminal --title='cfdemSolverPiso ErgunTestMPI_cgs CFD'  -e "bash $casePath/parCFDDEMrun.sh"
. $casePath/parCFDDEMrun.sh
