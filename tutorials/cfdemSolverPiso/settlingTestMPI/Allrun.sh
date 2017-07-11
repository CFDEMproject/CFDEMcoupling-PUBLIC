#!/bin/bash

#===================================================================#
# allrun script for testcase as part of test routine 
# run settlingTest
# Christoph Goniva - Feb. 2011
#===================================================================#

#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"

# check if mesh was built
if [ -f "$casePath/CFD/constant/polyMesh/points" ]; then
    echo "mesh was built before - using old mesh"
else
    echo "mesh needs to be built"
    cd $casePath/CFD
    blockMesh
fi

# keep old couplingProperties
cp $casePath/CFD/constant/couplingProperties $casePath/CFD/constant/couplingProperties_backup
cp $casePath/CFD/system/controlDict $casePath/CFD/system/controlDict_backup

# change to M2M
#changeDictionary -constant -dict changeDicts/changeDictionaryDict_1 -case $casePath/CFD

# change to MPI
#changeDictionary -constant -dict changeDicts/changeDictionaryDict_2 -case $casePath/CFD

# change to subTS
#changeDictionary -constant -dict changeDicts/changeDictionaryDict_3 -case $casePath/CFD

#- run parallel CFD-DEM in new terminal
#gnome-terminal --title='cfdemSolverPiso settlingTest CFD'  -e "bash $casePath/parCFDDEMrun.sh"
. $casePath/parCFDDEMrun.sh

# restore old couplingProperties
mv $casePath/CFD/constant/couplingProperties_backup $casePath/CFD/constant/couplingProperties
mv $casePath/CFD/system/controlDict_backup $casePath/CFD/system/controlDict
