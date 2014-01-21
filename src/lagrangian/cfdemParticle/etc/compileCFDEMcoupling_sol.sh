#!/bin/bash

#===================================================================#
# compile routine for CFDEMcoupling solvers, part of CFDEMproject 
# Christoph Goniva - May. 2012, DCS Computing GmbH
#===================================================================#

whitelist="solver-list.txt"

#- include functions
source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh
logDir="log"
cd $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc
mkdir -p $logDir

CWD="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
NOW="$(date +"%Y-%m-%d-%H:%M")"

echo ""
echo "This routine will compile the solvers specified in solver-list.txt"
echo ""
#echo "Are the variables CFDEM_SOLVER_DIR=$CFDEM_SOLVER_DIR"
#echo "and CFDEM_SRC_DIR=$CFDEM_SRC_DIR/lagrangian/cfdemParticle correct? (y/n)"
#read YN
#if [ "$YN" != "y" ];then
#  	echo "Aborted by user."
#  	exit 1
#fi

echo ""
echo "Please provide the solvers to be compiled in the $CWD/$whitelist file."
echo "structure:"
echo "path  to provide the path relative to CFDEM_SOLVER_DIR"
echo ""
echo "example:"
echo "cfdemSolverPiso/dir"
echo ""

if [ ! -f "$CWD/$whitelist" ];then
    echo "$whitelist does not exist in $CWD"
else
    NLINES=`wc -l < $CWD/$whitelist`
    COUNT=0

    for masterLogFile in "$masterLogName"
    do

        while [ $COUNT -lt $NLINES ]
        do
            let COUNT++  
            LINE=`head -n $COUNT $CWD/$whitelist | tail -1`
  
            # white lines
            if [[ "$LINE" == "" ]]; then
                continue
            # comments
            elif [[ "$LINE" == \#* ]]; then
                continue
             # paths
            elif [[ "$LINE" == */dir ]]; then
                echo "change path"
                LINE=$(echo "${LINE%????}")
                path="$CFDEM_SOLVER_DIR/$LINE"
                cd $path
                #continue
            fi

            #- execute tutorial
            echo "running testcase $path"
            #bash Allrun.sh

            #--------------------------------------------------------------------------------#
            #- define variables
            logpath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")/$logDir"
            logfileName="log_compileCFDEMcoupling""_$LINE" 
            casePath="$CFDEM_SOLVER_DIR/$LINE"
            headerText="$logfileName""_$LINE""-$NOW"
            #--------------------------------------------------------------------------------#
            compileSolver $logpath $logfileName $casePath $headerText

            #echo "did the solvers compile correcty? - press enter to proceed."
            #read
        done
    done
fi
