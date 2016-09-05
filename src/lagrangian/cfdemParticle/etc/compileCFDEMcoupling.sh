#!/bin/bash

#===================================================================#
# compile routine for CFDEMcoupling, part of CFDEMproject 
# Christoph Goniva - May. 2012, DCS Computing GmbH
#===================================================================#

#- include functions
source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh

NOW="$(date +"%Y-%m-%d-%H:%M")"
logDir="log"


cd $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc
mkdir -p $logDir

#================================================================================#
# compile src
#================================================================================#
. $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/compileCFDEMcoupling_src.sh

#================================================================================#
# compile solvers
#================================================================================#
. $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/compileCFDEMcoupling_sol.sh

#================================================================================#
# compile utilities
#================================================================================#
. $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/compileCFDEMcoupling_uti.sh
