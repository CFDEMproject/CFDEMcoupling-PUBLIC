#!/bin/bash

#===================================================================#
# script to run the basic examples
# Christoph Goniva - June 2012, DCS Computing GmbH
#===================================================================#

whitelist="tutorial-list.txt"

CWD="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
NOW="$(date +"%Y-%m-%d-%H:%M")"

echo ""
echo "This routine will execute the CFDEMcoupling tutorial cases specified in tutorial-list.txt"
echo ""
echo "Are the variables CFDEM_TUT_DIR=$CFDEM_TUT_DIR"
echo "and CFDEM_SRC_DIR=$CFDEM_SRC_DIR correct? (y/n)"
read YN
if [ "$YN" != "y" ];then
  	echo "Aborted by user."
  	exit 1
fi

echo ""
echo "Please provide the examples to be calculated in the $CWD/$whitelist file."
echo "structure:"
echo "path  to provide the path relative to CFDEM_TUT_DIR"
echo ""
echo "example:"
echo "cfdemSolverPiso/settlingTestMPI/dir"
echo ""

if [ ! -f "$CWD/$whitelist" ];then
    echo "$whitelist does not exist in $CWD"
else
    NLINES=`wc -l < $CWD/$whitelist`
    COUNT=0

    for masterLogFile in "$masterLogName" #"$masterLogName""_valgrind" 
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
                path="$CFDEM_TUT_DIR/$LINE"
                cd $path
                #continue
            fi

            #- execute tutorial
            echo "running testcase $path"
            bash Allrun.sh

            echo "did the case run correcty? - press enter to proceed."
            read
        done
    done
fi

## run pvg tool on logfile
#cd $CWD
#grep "==" "$masterLogName""_valgrind" > parallel_"$masterLogName""_valgrind"

## sort by first arg (+0 -0) and disable last resort comparison (-s)
## so sorted by first arg only
#sort +0 -0 -s parallel_"$masterLogName""_valgrind" > tmp
#mv tmp parallel_"$masterLogName""_valgrind"

