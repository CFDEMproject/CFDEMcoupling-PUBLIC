#!/bin/bash

#===================================================================#
# compile routine for LIGGGHTS, part of CFDEMproject 
# Christoph Goniva - May. 2012, DCS Computing GmbH
#===================================================================

#- include functions
source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh

NOW="$(date +"%Y-%m-%d-%H:%M")"
logDir="log"

cd $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc
mkdir -p $logDir

#================================================================================#
# copy LIGGGHTS patch files if available
#================================================================================#
echo "copying patch files for LIGGGHTS if available"
cp $CFDEM_SRC_DIR/LIGGGHTSpatch/* $CFDEM_LIGGGHTS_SRC_DIR

#--------------------------------------------------------------------------------#
#- define variables
logpath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")/$logDir"
logfileName="log_compileLIGGGHTS" #alternative: logfileName="log_compileLIGGGHTS_$NOW"
headerText="$logfileName""-$NOW"
#--------------------------------------------------------------------------------#

#================================================================================#
# compile LIGGGHTS libraries
#================================================================================#
#bash $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/compileLIGGGHTS_lib.sh

#================================================================================#
# compile LIGGGHTS src
#================================================================================#

compileLIGGGHTS $logpath $logfileName $headerText "false"


