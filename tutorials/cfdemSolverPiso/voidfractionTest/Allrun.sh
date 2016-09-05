#!/bin/bash

#===================================================================#
# allrun script for testcase as part of test routine 
# run voidfractionTest
# Christoph Goniva - Jan. 2016
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

#- run parallel CFD-DEM in new terminal
#gnome-terminal --title='cfdemSolverPiso voidfractionTest CFD'  -e "bash $casePath/parCFDDEMrun.sh"
. $casePath/parCFDDEMrun.sh
