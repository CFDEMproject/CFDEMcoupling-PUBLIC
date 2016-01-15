#!/bin/bash

#===================================================================#
# allrun script for testcase as part of test routine 
# run settlingTest CFD part
# Christoph Goniva - May. 2011
#===================================================================#

#- source CFDEM env vars
. ~/.bashrc

#- include functions
source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh

#--------------------------------------------------------------------------------#
#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
logpath=$casePath
headerText="run_parallel_cfdemSolverPisoScalar_packedBedTemp_CFDDEM"
logfileName="log_$headerText"
solverName="cfdemSolverPisoScalar"
nrProcs="2"
machineFileName="none"   # yourMachinefileName | none
debugMode="off"          # on | off| strict
testHarnessPath="$CFDEM_TEST_HARNESS_PATH"
runOctave="true"
postproc="false"
cleanup="true"
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

fi

#-------------------------------------------------------#
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
if [ $cleanup == "true" ]
  then
    echo "deleting data at: $casePath : \n"
    source $WM_PROJECT_DIR/bin/tools/CleanFunctions
    cd $casePath/CFD
    cleanCase
    rm -r $casePath/CFD/clockData
    rm $casePath/DEM/post/*.*
    #rm -r $casePath/DEM/post/restart/*.*
    echo "done"

    #- preserve post directory
    touch $casePath/DEM/post/.gitignore
    touch $casePath/DEM/post/restart/.gitignore
fi


