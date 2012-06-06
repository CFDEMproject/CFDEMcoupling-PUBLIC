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
headerText="run_parallel_cfdemSolverPiso_settlingTestMPI_CFDDEM"
logfileName="log_$headerText"
solverName="cfdemSolverPiso"
nrProcs="2"
machineFileName="none"   # yourMachinefileName | none
debugMode="off"          # on | off
testHarnessPath="$CFDEM_TEST_HARNESS_PATH"
#--------------------------------------------------------------------------------#

#- call function to run a parallel CFD-DEM case
parCFDDEMrun $logpath $logfileName $casePath $headerText $solverName $nrProcs $machineFileName $debugMode

#------------------------------#
#  octave

#- change path
cd octave

#- rmove old graph
rm cfdemSolverPiso_settlingTestMPI.eps

#- run octave
octave settlingVelocity.m

#- show plot 
evince cfdemSolverPiso_settlingTestMPI.eps
#------------------------------#

#- copy log file to test harness
cp ../../$logfileName $testHarnessPath
cp cfdemSolverPiso_settlingTestMPI.eps $testHarnessPath

    #- clean up case
    echo "deleting data at: $casePath :\n"
    rm -r $casePath/CFD/0.*
    rm -r $casePath/CFD/log.*
    rm -r $casePath/log_*
    rm -r $casePath/CFD/octave/octave-core
    rm -r $casePath/CFD/VTK
    rm -r $casePath/CFD/clockData
    rm -r $casePath/CFD/processor*
    rm -r $casePath/CFD/particles
    rm -r $casePath/CFD/couplingFiles/*
    rm -r $casePath/DEM/post/*
    rm -r $casePath/DEM/log.*
    rm -r $casePath/CFD/probes
    echo "done"

#- preserve post directory
echo "dummyfile" >> $casePath/DEM/post/dummy

