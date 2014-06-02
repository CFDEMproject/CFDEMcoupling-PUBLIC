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
headerText="run_parallel_cfdemSolverPiso_ErgunTestMPI_CFDDEM"
logfileName="log_$headerText"
solverName="cfdemSolverPiso"
nrProcs="4"
machineFileName="none"   # yourMachinefileName | none
debugMode="off"          # on | off| strict
reconstuctCase="true"    # true | false
testHarnessPath="$CFDEM_TEST_HARNESS_PATH"
#--------------------------------------------------------------------------------#

#- call function to run a parallel CFD-DEM case
parCFDDEMrun $logpath $logfileName $casePath $headerText $solverName $nrProcs $machineFileName $debugMode $reconstuctCase

