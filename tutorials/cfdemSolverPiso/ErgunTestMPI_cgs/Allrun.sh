#!/bin/bash

#===================================================================#
# allrun script for testcase as part of test routine 
# run ErgunTestMPI_cgs (working in CGS unit system)
# Christoph Goniva - March 2013
#===================================================================#

#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"

# check if mesh was built
if [ -d "$casePath/CFD/constant/polyMesh/boundary" ]; then
    echo "mesh was built before - using old mesh"
else
    echo "mesh needs to be built"
    cd $casePath/CFD
    blockMesh
fi

echo "WARNING:copying a CGS based controlDict to $HOME/.OpenFOAM/$WM_PROJECT_VERSION"
echo "this will make your simulations use CGS unit system"
read
echo "Make sure $HOME/.OpenFOAM/$WM_PROJECT_VERSION/controlDict is removed after this simulation."

mkdir -p $FOAM_INST_DIR/.OpenFOAM//$WM_PROJECT_VERSION
cp $CFDEM_SRC_DIR/etc/controlDict_cgs_$WM_PROJECT_VERSION $FOAM_INST_DIR/.OpenFOAM/$WM_PROJECT_VERSION/controlDict

#- run parallel CFD-DEM in new terminal
gnome-terminal --title='cfdemSolverPiso ErgunTestMPI CFD'  -e "bash $casePath/parCFDDEMrun.sh" 

echo "removing $FOAM_INST_DIR/.OpenFOAM/$WM_PROJECT_VERSION/controlDict?"
read
rm -r $FOAM_INST_DIR/.OpenFOAM/$WM_PROJECT_VERSION/controlDict*
