#!/bin/bash

#===================================================================#
# pull routine for LIGGGHTS, part of CFDEMproject 
# Christoph Goniva - August. 2013, DCS Computing GmbH
#===================================================================

#- include functions
source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh

NOW="$(date +"%Y-%m-%d-%H:%M")"
logDir="log"

cd $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc
mkdir -p $logDir

#--------------------------------------------------------------------------------#
#- define variables
logpath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")/$logDir"
logfileName="log_pullLIGGGHTS"
casePath="$CFDEM_LIGGGHTS_SRC_DIR"
headerText="$logfileName""-$NOW"
#--------------------------------------------------------------------------------#

pullRepo $logpath $logfileName $casePath $headerText
