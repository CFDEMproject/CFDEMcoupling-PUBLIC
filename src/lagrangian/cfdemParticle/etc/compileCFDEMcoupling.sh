#!/bin/bash

#===================================================================#
# compile routine for CFDEMcoupling, part of CFDEMproject 
# Christoph Goniva - May. 2012, DCS Computing GmbH
#===================================================================#

#- include functions
source $CFDEM_SRC_DIR/etc/functions.sh

NOW="$(date +"%Y-%m-%d-%H:%M")"
logDir="log"


cd $CFDEM_SRC_DIR/etc
mkdir $logDir


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

#================================================================================#
# compile solvers
#================================================================================#


for solverName in "cfdemSolverPiso" "cfdemSolverIB" "cfdemSolverPisoScalar"
do
    #--------------------------------------------------------------------------------#
    #- define variables
    logfileName="log_compileCFDEMcoupling""_$solverName" 
    casePath="$CFDEM_SOLVER_DIR/$solverName"
    headerText="$logfileName""_$solverName""-$NOW"
    #--------------------------------------------------------------------------------#
    compileSolver $logpath $logfileName $casePath $headerText
done

echo "Note: the list of solvers compiled might be incomplete."
echo "please check $CFDEM_SOLVER_DIR for more solvers available"
