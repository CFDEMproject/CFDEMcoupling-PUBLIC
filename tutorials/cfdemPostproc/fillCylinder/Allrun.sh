#!/bin/bash

#===================================================================#
# allrun script for cfdemPostproc
# Christoph Goniva - Nov. 2011
#===================================================================#

#- source CFDEM env vars
. ~/.bashrc

#- include functions
source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh

#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"

liggghtsSim="true"
cfdemPostProc="true"
postproc="true"

# check if mesh was built
if [ -f "$casePath/CFD/constant/polyMesh/points" ]; then
    echo "mesh was built before - using old mesh"
else
    echo "mesh needs to be built"
    cd $casePath/CFD
    blockMesh
fi

if [ $liggghtsSim == "true" ]
  then
    #--------------------------------------------------------------------------------#
    #- define variables
    logpath="$casePath"
    headerText="run_liggghts_fillCylinder_DEM"
    logfileName="log_$headerText"
    solverName="in.liggghts_init"
    #--------------------------------------------------------------------------------#

    #- clean up case
    rm -r $casePath/DEM/post/*

    #- call function to run DEM case
    DEMrun $logpath $logfileName $casePath $headerText $solverName


    #- generate VTK data (no longer needed as we directly write vtk)
    #cd $casePath/DEM/post
    #python $CFDEM_LPP_DIR/lpp.py  dump.liggghts_init

fi

if [ $cfdemPostProc == "true" ]
  then
    #--------------------------------------------------------------------------------#
    #- define variables
    logpath="$casePath"
    headerText="run_cfdemPostproc_fillCylinder_CFD"
    logfileName="log_$headerText"
    solverName="cfdemPostproc"
    debugMode="off"          # on | off | strict
    #--------------------------------------------------------------------------------#

    #- clean up case
    rm -r $casePath/CFD/0.*

    #- call function to run CFD cas
    CFDrun $logpath $logfileName $casePath $headerText $solverName $debugMode
fi

if [ $postproc == "true" ]
  then

    #- get VTK data from CFD sim
    #foamToVTK
    
    #- start paraview
    echo ""
    echo "trying to start paraview..."
    #paraview4 # use your start command for paraview here
    $HOME/software/ParaView-4.3.1-Linux-64bit/bin/paraview
    read
fi

#- keep terminal open (if started in new terminal)
#echo "...press enter to clean up case"
#echo "press Ctr+C to keep data"
#read

#- clean up case
echo "deleting data at: $casePath ?\n"
read
source $WM_PROJECT_DIR/bin/tools/CleanFunctions
cd $casePath/CFD
cleanCase
rm -r $casePath/CFD/postProcessing
rm -r $casePath/CFD/lagrangian
rm -r $casePath/CFD/clockData
rm $casePath/CFD/octave/octave-workspace
rm -r $casePath/CFD/hpctoolkit-*
rm $casePath/DEM/post/*.*
touch $casePath/DEM/post/.gitignore
#rm $casePath/DEM/post/restart/*.*
rm $casePath/DEM/post/restart/liggghts.restartCFDEM*
touch $casePath/DEM/post/restart/.gitignore
echo "done"
