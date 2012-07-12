#!/bin/bash

#===================================================================#
# allrun script for testcase as part of test routine 
# run settlingTest CFD part
# Christoph Goniva - Feb. 2011
#===================================================================#

#- source CFDEM env vars
. ~/.bashrc

#- include functions
source $CFDEM_SRC_DIR/etc/functions.sh

#--------------------------------------------------------------------------------#
#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
logpath=$casePath
headerText="run_parallel_cfdemSolverPiso_ErgunTestMPI_CFDDEM"
logfileName="log_$headerText"
solverName="cfdemSolverPiso"
nrProcs="8"
machineFileName="none"   # yourMachinefileName | none
debugMode="off"          # on | off
testHarnessPath="$CFDEM_TEST_HARNESS_PATH"
runOctave="true"
postproc="false"
#--------------------------------------------------------------------------------#

#- call function to run a parallel CFD-DEM case
parCFDDEMrun $logpath $logfileName $casePath $headerText $solverName $nrProcs $machineFileName $debugMode


if [ $runOctave == "true" ]
    then
        #------------------------------#
        #  octave

        #- change path
        cd octave

        #- rmove old graph
        rm cfdemSolverPiso_ErgunTestMPI.eps

        #- run octave
        octave totalPressureDrop.m

        #- show plot 
        evince cfdemSolverPiso_ErgunTestMPI.eps

        #- copy log file to test harness
        cp ../../$logfileName $testHarnessPath
        cp cfdemSolverPiso_ErgunTestMPI.eps $testHarnessPath
fi

if [ $postproc == "true" ]
  then

    #- keep terminal open (if started in new terminal)
    echo "simulation finisehd? ...press enter to proceed"
    read

    #- get VTK data from liggghts dump file
    cd $casePath/DEM/post
    python -i $CFDEM_LPP_DIR/lpp.py dump*.liggghts_restart

    #- get VTK data from CFD sim
    cd $casePath/CFD
    foamToVTK                                                   #- serial run of foamToVTK
    #source $CFDEM_SRC_DIR/etc/functions.sh                       #- include functions
    #pseudoParallelRun "foamToVTK" $nrPostProcProcessors          #- pseudo parallel run of foamToVTK

    #- start paraview
    paraview

    #- keep terminal open (if started in new terminal)
    echo "...press enter to clean up case"
    echo "press Ctr+C to keep data"
    read

fi

#- clean up case
echo "deleting data at: $casePath :\n"
rm -r $casePath/CFD/0.*
rm -r $casePath/CFD/log.*
rm -r $casePath/CFD/octave/octave-core
rm -r $casePath/CFD/VTK
rm -r $casePath/CFD/processor*
rm -r $casePath/CFD/couplingFiles/*
rm -r $casePath/DEM/post/*
rm -r $casePath/DEM/log.*
rm -r $casePath/DEM/liggghts.restartCFDEM*
rm -r $casePath/CFD/probes
rm -r $casePath/CFD/particles
rm -r $casePath/CFD/clockData
echo "done"

#- preserve post directory
echo "dummyfile" >> $casePath/DEM/post/dummy
