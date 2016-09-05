#!/bin/bash

#===================================================================#
# sytsem settings test routine for cfdem project 
# Christoph Goniva - May. 2011, DCS Computing GmbH
#===================================================================#

#- include functions
source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh

#- show gcc settings
checkGPP="true"

#- sys check for add on
checkAddOn="true"

#- system settings
printHeader

echo "*********************************"
echo "CFDEM(R)coupling system settings:"
echo "*********************************"

echo "CFDEM_VERSION=$CFDEM_VERSION"
echo "couple to OF_VERSION=$WM_PROJECT_VERSION"
echo "compile option=$WM_COMPILE_OPTION"


echo
echo "check if paths are set correctly"
checkDirComment "$CFDEM_PROJECT_DIR" '$CFDEM_PROJECT_DIR' "yes"
checkDirComment "$CFDEM_PROJECT_USER_DIR" '$CFDEM_PROJECT_USER_DIR' "no"
checkDirComment "$CFDEM_SRC_DIR" '$CFDEM_SRC_DIR' "yes"
checkDirComment "$CFDEM_SOLVER_DIR" '$CFDEM_SOLVER_DIR' "yes"
checkDirComment "$CFDEM_TUT_DIR" '$CFDEM_TUT_DIR' "yes"
checkDirComment "$CFDEM_LIGGGHTS_SRC_DIR" '$CFDEM_LIGGGHTS_SRC_DIR' "yes"
checkDirComment "$CFDEM_LPP_DIR" '$CFDEM_LPP_DIR' "yes"
checkDirComment "$CFDEM_ADD_LIBS_DIR" '$CFDEM_ADD_LIBS_DIR' "yes"
checkDirComment "$CFDEM_PIZZA_DIR" '$CFDEM_PIZZA_DIR' "no"
checkDirComment "$CFDEM_TEST_HARNESS_PATH" '$CFDEM_TEST_HARNESS_PATH' "no"
checkDirComment "$C3PO_SRC_DIR" '$C3PO_SRC_DIR' "no"
echo ""

echo "library names"
echo '$CFDEM_LIGGGHTS_LIB_NAME = '"$CFDEM_LIGGGHTS_LIB_NAME"
echo '$CFDEM_LIB_NAME = '"$CFDEM_LIB_NAME"
echo '$LD_LIBRARY_PATH  = '"$LD_LIBRARY_PATH"
echo '$WM_NCOMPPROCS  = '"$WM_NCOMPPROCS"
echo '$WM_LABEL_SIZE = '"$WM_LABEL_SIZE"
if [ $WM_LABEL_SIZE != 32 ]; then 
    echo "!!!! Warning: WM_LABEL_SIZE must be 32!!!!! (Please correct in $FOAM_ETC/bashrc.)"
fi

echo "*******************"


if [ $checkGPP == "true" ]
  then
    echo "g++:"
    which g++
    g++ --version

    echo "gcc:"
    which gcc
    gcc --version

    echo "mpic++:"
    which mpic++
    mpic++ --version

    echo "mpirun:"
    which mpirun
    mpirun --version
fi

echo "**********************"
echo "additional packages..."

if [ $checkAddOn == "true" ]
  then
    packageName=c3po
    filePath=$CFDEM_SRC_DIR/$packageName
    if [ $(checkDir $filePath) == "true" ]; then
        source $filePath/etc/$packageName"SystemTest.sh"
    else
        echo "$packageName does not exist." 
    fi
fi
