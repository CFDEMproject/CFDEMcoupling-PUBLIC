#!/bin/bash

#===================================================================#
# allrun script for testcase as part of test routine 
# run settlingTest CFD part
# Christoph Goniva - May. 2011
#===================================================================#

#- source CFDEM env vars
. ~/.bashrc

#- include functions
source $CFDEM_SRC_DIR/etc/functions.sh

#--------------------------------------------------------------------------------#
#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
logpath=$casePath
headerText="run_parallel_cfdemSolverPisoScalar_packedBedTemp_CFDDEM"
logfileName="log_$headerText"
solverName="cfdemSolverPisoScalar"
nrProcs="2"
machineFileName="none"   # yourMachinefileName | none
debugMode="off"          # on | off
testHarnessPath="$CFDEM_TEST_HARNESS_PATH"
runOctave="true"
#--------------------------------------------------------------------------------#

#- call function to run a parallel CFD-DEM case
parCFDDEMrun $logpath $logfileName $casePath $headerText $solverName $nrProcs $machineFileName $debugMode

#------------------------------#
if [ $runOctave == "true" ]
  then

    #- change path
    cd octave

    #- rmove old graph
    rm *.eps

    #- run octave
    octave totalPressureDropAndNusselt.m

    #- show plots 
    evince cfdemSolverPisoScalar_Nusselt.eps &
    evince cfdemSolverPisoScalar_pressureDrop.eps
    #------------------------------#

    #- copy log file to test harness
    cp ../../$logfileName $testHarnessPath
    cp cfdemSolverPisoScalar_Nusselt.eps $testHarnessPath
    cp cfdemSolverPisoScalar_pressureDrop.eps $testHarnessPath

    #- clean up case
    echo "deleting data at: $casePath : ???\n"
    rm -r $casePath/CFD/0.*
    rm -r $casePath/CFD/octave/*.eps
    rm -r $casePath/CFD/octave/octave-core
    rm -r $casePath/CFD/VTK
    rm -r $casePath/CFD/processor*
    rm -r $casePath/DEM/post/*
    rm -r $casePath/DEM/log.*
    rm -r $casePath/CFD/log.*
    rm -r $casePath/CFD/probes
    rm -r $casePath/log_*
    echo "done"
fi

#- preserve post directory
echo "dummyfile" >> $casePath/DEM/post/dummy


