#!/bin/bash

#===================================================================#
# compile routine for LIGGGHTS, part of CFDEMproject 
# Christoph Goniva - May. 2012, DCS Computing GmbH
#===================================================================

#- include functions
source $CFDEM_SRC_DIR/etc/functions.sh

NOW="$(date +"%Y-%m-%d-%H:%M")"
logDir="log"

cd $CFDEM_SRC_DIR/etc
mkdir -p $logDir

#--------------------------------------------------------------------------------#
#- define variables
logpath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")/$logDir"
logfileName="log_compileLIGGGHTS" #alternative: logfileName="log_compileLIGGGHTS_$NOW"
headerText="$logfileName""-$NOW"
#--------------------------------------------------------------------------------#

compileLIGGGHTS $logpath $logfileName $headerText
