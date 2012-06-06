#!/bin/bash

#===================================================================#
# allrun script for testcase as part of test routine 
# run settlingTest
# Christoph Goniva - Sept. 2010
#===================================================================#

#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
pizzaPath="$CFDEM_PIZZA_DIR"
postproc="false"
nrPostProcProcessors=2

# check if mesh was built
if [ -d "$casePath/CFD/constant/polyMesh/boundary" ]; then
    echo "mesh was built before - using old mesh"
else
    echo "mesh needs to be built"
    cd $casePath/CFD
    blockMesh
fi

#-------------------------------------------------------#
# adapt settings for init run
cp $casePath/CFD/constant/liggghtsCommands_init $casePath/CFD/constant/liggghtsCommands
cp $casePath/CFD/constant/couplingProperties_init $casePath/CFD/constant/couplingProperties
cp $casePath/CFD/system/controlDict_init $casePath/CFD/system/controlDict
#-------------------------------------------------------#

#- run parallel CFD-DEM in new terminal
gnome-terminal --title='cfdemSolverPiso ErgunTestMPI_restart CFD'  -e "bash $casePath/parCFDDEMrun.sh" 

#- wait until sim has finished then run octave
echo "simulation finished? ...press enter to proceed"
read


#-------------------------------------------------------#
# adapt settings for init or restart run
cp $casePath/CFD/constant/liggghtsCommands_restart $casePath/CFD/constant/liggghtsCommands
cp $casePath/CFD/constant/couplingProperties_restart $casePath/CFD/constant/couplingProperties
cp $casePath/CFD/system/controlDict_restart $casePath/CFD/system/controlDict

#- run parallel CFD-DEM in new terminal
gnome-terminal --title='cfdemSolverPiso ErgunTestMPI_restart CFD'  -e "bash $casePath/parCFDDEMrun.sh" 

#-------------------------------------------------------#

if [ $postproc == "true" ]
  then

    #- keep terminal open (if started in new terminal)
    echo "simulation finisehd? ...press enter to proceed"
    read

    #- get VTK data from liggghts dump file
    cd $casePath/DEM
    python $pizzaPath/pizza.py -f pizzaScriptRestart

    #- get VTK data from CFD sim
    cd $casePath/CFD
    #foamToVTK                                                   #- serial run of foamToVTK
    source $CFDEM_SRC_DIR/etc/functions.sh                       #- include functions
    pseudoParallelRun "foamToVTK" $nrPostProcProcessors          #- pseudo parallel run of foamToVTK

    #- start paraview
    paraview

    #- keep terminal open (if started in new terminal)
    echo "...press enter to clean up case"
    echo "press Ctr+C to keep data"
    read

    #- clean up case
    echo "deleting data at: $casePath"
    rm -r $casePath/CFD/0.*
    rm -r $casePath/CFD/clockData
    rm -r $casePath/CFD/VTK
    rm -r $casePath/CFD/couplingFiles/*
    rm -r $casePath/CFD/particles/*
    rm -r $casePath/DEM/post/*
    rm -r $casePath/DEM/log.*
    rm -r $casePath/log_*
    echo "done"

fi

#- wait until sim has finished then run octave
echo "simulation finished? ...press enter to proceed"
read

#------------------------------#
#  octave

#- change path
cd $casePath/CFD/octave

#- rmove old graph
rm cfdemSolverPiso_ErgunTestMPI.eps

#- run octave
octave totalPressureDrop.m

#- show plot 
evince cfdemSolverPiso_ErgunTestMPI.eps
#------------------------------#


#- copy log file to test harness
cp $casePath/../$logfileName $testHarnessPath
cp cfdemSolverPiso_ErgunTestMPI.eps $testHarnessPath

#- clean up case
rm -rf $casePath/CFD/0.*
rm -r $casePath/CFD/clockData
rm -rf $casePath/CFD/processor*
rm -rf $casePath/CFD/patchAverage_pressureDrop
rm -rf $casePath/CFD/probes
rm -rf $casePath/CFD/particles
rm $casePath/log.liggghts
rm $casePath/DEM/liggghts.restartCFDEM
rm $casePath/DEM/post/dump.*


#- preserve post directory
echo "dummyfile" >> $casePath/DEM/post/dummy
