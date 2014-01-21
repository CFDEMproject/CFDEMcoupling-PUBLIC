#!/bin/bash

#===================================================================#
# compile routine for CFDEMcoupling solvers, part of CFDEMproject 
# Christoph Goniva - May. 2012, DCS Computing GmbH
#===================================================================#

#- include functions
source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh

NOW="$(date +"%Y-%m-%d-%H:%M")"
logDir="log"


cd $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc
mkdir -p $logDir

#================================================================================#
# compile utilities
#================================================================================#


for utName in "cfdemPostproc"
do
    #--------------------------------------------------------------------------------#
    #- define variables
    logpath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")/$logDir"
    logfileName="log_compileCFDEMcoupling""_$utName" 
    casePath="$CFDEM_UT_DIR/$utName"
    headerText="$logfileName""_$utName""-$NOW"
    #--------------------------------------------------------------------------------#
    compileSolver $logpath $logfileName $casePath $headerText
done

echo "Note: the list of utilities compiled might be incomplete."
echo "please check $CFDEM_UT_DIR for more utilities available"
