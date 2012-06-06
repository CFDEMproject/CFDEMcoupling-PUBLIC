#!/bin/bash

#===================================================================#
# allrun script for testcase as part of test routine 
# run settlingTest
# Christoph Goniva - July. 2011, mod by Alice Hager - July 2011
#===================================================================#

source $CFDEM_SRC_DIR/etc/functions.sh

#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
pizzaPath="$CFDEM_PIZZA_DIR"
postproc="false"

echo $casePath

# check if mesh was built
if [ -d "$casePath/CFD/constant/polyMesh/boundary" ]; then
    echo "mesh was built before - using old mesh"
else
    echo "mesh needs to be built"
    cd $casePath/CFD
    blockMesh
fi


#gnome-terminal --title='cfdemSolverIB two settling disks CFD' -e "CFDrun()"
gnome-terminal --title='cfdemSolverIB twoSpheresGlowinskiMPI CFD' -e "bash $casePath/parCFDDEMrun.sh"
#gnome-terminal --title='cfdemSolverIB two settling disks DEM' -e "DEMrun()"
#gnome-terminal --title='cfdemSolverIB twoSpheresGlowinskiMPI DEM' -e "bash $casePath/DEMrun.sh"



if [ $postproc == "true" ]
  then
    #- get VTK data from liggghts dump file
    cd $casePath/DEM
    python -i $pizzaPath/pizza.py -f pizzaScriptInit

    #- get VTK data from CFD sim
    cd $casePath/CFD
    foamToVTK

    #- start paraview
    paraview

    #- clean up case
    echo "deleting data at: $casePath"    
    rm -r $casePath/CFD/0.*
    rm -r $casePath/CFD/VTK
    rm -r $casePath/CFD/couplingFiles/*
    rm -r $casePath/DEM/post/*
    rm -r $casePath/DEM/log.*
    rm -r $casePath/log_*
    echo "done"
fi




