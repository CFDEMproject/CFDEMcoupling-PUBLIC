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

#- remove old success/fail logs
rm $logDir/log_compile_results_success
rm $logDir/log_compile_results_fail

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
    njobs=`wc -l < $CWD/$whitelist` 
    echo ""
    echo "running compilation in pseudo-parallel mode of $njobs solvers"

    #--------------------------------------------------------------------------------#
    logpath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")/$logDir"

    ##number of solvers compiled at a time
    nsteps=$WM_NCOMPPROCS
    echo "do compilation on $nsteps procs"
    nchunk=`echo $njobs/$nsteps+1 | bc`
    if [[ $WM_NCOMPPROCS == "" ]]; then
        echo "do compilation in serial"
        nsteps=1
        nchunk=1
    else
        echo "do compilation on $nsteps procs in $nchunk chunks"      
    fi

    counter=0
    for i in `seq $nchunk`
    do

        #wait until prev. compilation is finished
        echo "waiting..."
        until [ `ps -C make | wc -l` -eq 1 ]; 
        do 
            sleep 2
        done

        for j in `seq $nsteps`
        do
            let solNr=($i-1)*$nsteps+$j
            LINE=`head -n $solNr $CWD/$whitelist | tail -1`

            # white lines
            if [[ "$LINE" == "" ]]; then
                continue
            # comments
            elif [[ "$LINE" == \#* ]]; then
                continue
            # paths
            elif [[ "$LINE" == */dir ]]; then
            #echo "change path"
                LINE=$(echo "${LINE%????}")
                path="$CFDEM_SOLVER_DIR/$LINE"
                #cd $path
                let solNr++
            fi

            if [[ "$counter" -lt "$njobs" ]]; then
                #--------------------------------------------------------------------------------#
                #- define variables
                #logpath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")/$logDir"
                logfileName="log_compileCFDEMcoupling""_$LINE" 
                casePath="$CFDEM_SOLVER_DIR/$LINE"
                headerText="$logfileName""_$LINE""-$NOW"
                parallel="true"
                #--------------------------------------------------------------------------------#

                #echo "compiling $LINE"
                compileSolver $logpath $logfileName $casePath $headerText $parallel
                let counter++
            fi
        done
    done

    echo "compilation done."
fi

#--------------------------------------------------------------------------------#
# loop all solvers and collect the logs
#--------------------------------------------------------------------------------#

#wait until prev. compilation is finished
echo "waiting..."
until [ `ps -C make | wc -l` -eq 1 ]; 
do 
    sleep 2
done

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
                LINE=$(echo "${LINE%????}")
                path="$CFDEM_SOLVER_DIR/$LINE"
                #continue
            fi

            #- execute tutorial
            echo "collecting log of $path"

            #--------------------------------------------------------------------------------#
            #- define variables
            #logpath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")/$logDir"
            logfileName="log_compileCFDEMcoupling""_$LINE" 
            casePath="$CFDEM_SOLVER_DIR/$LINE"
            #--------------------------------------------------------------------------------#            
            collectLogCFDEMcoupling_sol $logpath $logfileName $casePath
        done
    done
fi

