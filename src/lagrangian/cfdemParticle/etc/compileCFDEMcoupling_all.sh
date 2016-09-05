#!/bin/bash

#===================================================================#
# compile routine for CFDEMcoupling + LIGGGHTS, part of CFDEMproject 
# Christoph Goniva - May. 2012, DCS Computing GmbH
#                    update March 2014
#===================================================================#

#- include functions
source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh

NOW="$(date +"%Y-%m-%d-%H:%M")"
logDir="log"


cd $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc
mkdir -p $logDir

#================================================================================#
# compile LIGGGHTS src + libraries
#================================================================================#
. $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/compileLIGGGHTS.sh

#================================================================================#
# compile CFDEMcoupling
#================================================================================#
. $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/compileCFDEMcoupling.sh
