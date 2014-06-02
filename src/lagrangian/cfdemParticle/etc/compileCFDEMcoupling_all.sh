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
# compile LIGGGHTS src
#================================================================================#
bash $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/compileLIGGGHTS.sh

#================================================================================#
# compile LIGGGHTS libraries
#================================================================================#
bash $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/compileLIGGGHTS_lib.sh

#================================================================================#
# compile CFDEMcoupling
#================================================================================#
bash $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/compileCFDEMcoupling.sh
