#!/bin/bash

#===================================================================#
# allrun script for testcase as part of test routine 
# run settlingTest CFD part
# Christoph Goniva - Feb. 2011
#===================================================================#

#- source CFDEM env vars
. ~/.bashrc

#- include functions
source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh

#--------------------------------------------------------------------------------#
#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
logpath=$casePath
headerText="run_parallel_cfdemSolverPiso_settlingTestMPI_CFDDEM"
logfileName="log_$headerText"
solverName="cfdemSolverPiso"
nrProcs="2"
machineFileName="none"   # yourMachinefileName | none
debugMode="off"          # on | off| strict
testHarnessPath="$CFDEM_TEST_HARNESS_PATH"
runOctave="true"
cleanUp="true"
postproc="false"
#--------------------------------------------------------------------------------#

#- call function to run a parallel CFD-DEM case
parCFDDEMrun $logpath $logfileName $casePath $headerText $solverName $nrProcs $machineFileName $debugMode "true"

if [ $runOctave == "true" ]
    then
        #- change path
        cd octave

        #- rmove old graph
        rm cfdemSolverPiso_settlingTestMPI.png

        #- run octave
        octave --no-gui settlingVelocity.m

        #- show plot 
        eog cfdemSolverPiso_settlingTestMPI.png

        #- copy log file to test harness
        cp ../../$logfileName $testHarnessPath
        cp cfdemSolverPiso_settlingTestMPI.png $testHarnessPath
fi

if [ $postproc == "true" ]
  then

    #- keep terminal open (if started in new terminal)
    echo "simulation finished? ...press enter to proceed"
    read

    #- get VTK data from liggghts dump file
    cd $casePath/DEM/post
    python -i $CFDEM_LPP_DIR/lpp.py  dump.liggghts_run

    #- get VTK data from CFD sim
    cd $casePath/CFD
    foamToVTK                                                   #- serial run of foamToVTK
    #source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh                       #- include functions
    #pseudoParallelRun "foamToVTK" $nrPostProcProcessors          #- pseudo parallel run of foamToVTK

    #- start paraview
    paraview

    #- keep terminal open (if started in new terminal)
    echo "...press enter to clean up case"
    echo "press Ctr+C to keep data"
    read
fi

#- clean up case
if [ $cleanUp == "true" ]
  then
    keepDEMrestart="false"
    cleanCFDEMcase $casePath/CFD $keepDEMrestart
fi
