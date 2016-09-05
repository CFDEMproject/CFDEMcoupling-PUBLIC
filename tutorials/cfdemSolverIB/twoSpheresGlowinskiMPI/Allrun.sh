#!/bin/bash

#===================================================================#
# allrun script for testcase as part of test routine 
# run settlingTest
# Christoph Goniva - July. 2011, mod by Alice Hager - July 2011
#===================================================================#

source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh

#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"

echo $casePath

# check if mesh was built
if [ -f "$casePath/CFD/constant/polyMesh/points" ]; then
    echo "mesh was built before - using old mesh"
else
    echo "mesh needs to be built"
    cd $casePath/CFD
    blockMesh
fi


#gnome-terminal --title='cfdemSolverIB twoSpheresGlowinskiMPI CFD' -e "bash $casePath/parCFDDEMrun.sh"
. $casePath/parCFDDEMrun.sh
