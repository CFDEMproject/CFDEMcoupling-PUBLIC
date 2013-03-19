#!/bin/bash

#===================================================================#
# compile routine for CFDEMcoupling source, part of CFDEMproject 
# Christoph Goniva - May. 2012, DCS Computing GmbH
#===================================================================#

#- include functions
source $CFDEM_SRC_DIR/etc/functions.sh

NOW="$(date +"%Y-%m-%d-%H:%M")"
logDir="log"


cd $CFDEM_SRC_DIR/etc
mkdir -p $logDir


#================================================================================#
# compile src
#================================================================================#

#--------------------------------------------------------------------------------#
#- define variables
logpath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")/$logDir"
logfileName="log_compileCFDEMcoupling_src" #alternative: logfileName="log_compileLIGGGHTS_$NOW"
casePath="$CFDEM_SRC_DIR"
headerText="$logfileName""-$NOW"
#--------------------------------------------------------------------------------#
compileLib $logpath $logfileName $casePath $headerText
