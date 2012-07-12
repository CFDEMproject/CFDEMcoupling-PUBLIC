#!/bin/bash

#===================================================================#
# allrun script for cfdemPostproc
# Christoph Goniva - Nov. 2011
#===================================================================#

#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
pizzaPath="/home/cfdem/LIGGGHTS/PIZZA/gran_pizza_17Aug10/src"

liggghtsSim="true"
cfdemPostProc="true"
postproc="true"

# check if mesh was built
if [ -d "$casePath/CFD/constant/polyMesh/boundary" ]; then
    echo "mesh was built before - using old mesh"
else
    echo "mesh needs to be built"
    cd $casePath/CFD
    blockMesh
fi

if [ $liggghtsSim == "true" ]
  then
    cd $casePath/DEM
    liggghts < in.liggghts_init

    #- get VTK data from liggghts dump file
    #cd $casePath/DEM
    #python $pizzaPath/pizza.py -f pizzaScriptInit

    cd $casePath/DEM/post
    python -i $CFDEM_LPP_DIR/lpp.py  dump.liggghts_init

fi

if [ $cfdemPostProc == "true" ]
  then
    cd $casePath/CFD
    cfdemPostproc
fi

#echo "now you can run foamToSurface from 2.0.x to generate *.stl or *.inp files.(press enter)"
#read

if [ $postproc == "true" ]
  then

    #- get VTK data from CFD sim
    #cd $casePath/CFD
    foamToVTK
    
    #- start paraview
    paraview
fi

#- keep terminal open (if started in new terminal)
#echo "...press enter to clean up case"
#echo "press Ctr+C to keep data"
#read

#- clean up case
echo "deleting data at: $casePath"
rm -r $casePath/CFD/0.*
rm -r $casePath/CFD/particles
rm -r $casePath/CFD/VTK
rm -r $casePath/DEM/post/*
rm -r $casePath/DEM/log.*
rm -r $casePath/log*
echo "done"

#- preserve post directory
echo "dummyfile" >> $casePath/DEM/post/dummy
