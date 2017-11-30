#!/bin/bash

#===================================================================#
# compile routine for CFDEMcoupling utilities, part of CFDEMproject 
# Christoph Goniva - May. 2012, DCS Computing GmbH
#===================================================================#

whitelist="utilities-list.txt"

#- include functions
source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh

logDir="log"
cd $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc
mkdir -p $logDir

CWD="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
NOW="$(date +"%Y-%m-%d-%H:%M")"

echo ""
echo "This routine will compile the utilities specified in utilities-list.txt"
echo ""
#echo "Are the variables CFDEM_UT_DIR=$CFDEM_UT_DIR"
#echo "and CFDEM_SRC_DIR=$CFDEM_SRC_DIR/lagrangian/cfdemParticle correct? (y/n)"
#read YN
#if [ "$YN" != "y" ];then
#  	echo "Aborted by user."
#  	exit 1
#fi

echo ""
echo "Please provide the utilities to be compiled in the $CWD/$whitelist file."
echo "structure:"
echo "path  to provide the path relative to CFDEM_UT_DIR"
echo ""
echo "example:"
echo "cfdemPostproc/dir"
echo ""

#- create a tmp file and delete comments in it - work with tmp file then.
cp $whitelist "tmpFile.txt"
sed -i '/^#/d' "tmpFile.txt"
sed -i '/^$/d' "tmpFile.txt"
whitelist="tmpFile.txt"

if [ ! -f "$CWD/$whitelist" ];then
    echo "$whitelist does not exist in $CWD"
else
    njobs=`wc -l < $CWD/$whitelist` 
    echo ""
    echo "running compilation in pseudo-parallel mode of $njobs utilities"

    #--------------------------------------------------------------------------------#
    logpath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")/$logDir"

    ##number of utilities compiled at a time

    if [[ $WM_NCOMPPROCS == "" ]] || [ $WM_NCOMPPROCS -eq 1 ]; then
        nsteps=1
        let nchunk=$njobs+1 # +1, to wait for the last compilation too
        echo "do compilation in serial"
    else    
        nsteps=$WM_NCOMPPROCS
        let nchunk=$njobs/$nsteps+1
        echo "do compilation on $nsteps procs in $nchunk chunks" 
        let nchunk++ # +1, to wait for the last compilation too     
    fi

    counter=0
    for i in `seq $nchunk`
    do

        #wait until prev. compilation is finished
        echo "waiting..."
        #until [ `ps -a | grep make | wc -l` -eq 0 ]; 
        until [ `ls -a | grep $logpath/*.tempXYZ | wc -l` -eq 0 ];
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
                path="$CFDEM_UT_DIR/$LINE"
                #cd $path
                let solNr++
            fi

            if [[ "$counter" -lt "$njobs" ]]; then
                #--------------------------------------------------------------------------------#
                #- define variables
                #logpath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")/$logDir"
                logfileName="log_compileCFDEMcoupling""_$LINE" 
                casePath="$CFDEM_UT_DIR/$LINE"
                headerText="$logfileName""_$LINE""-$NOW"
                parallel="true"
                #--------------------------------------------------------------------------------#

                echo "compiling $LINE"
                compileSolver $logpath $logfileName $casePath $headerText $parallel
                let counter++
            fi
        done

        #sleep 1 # wait a second until compilation starts
    done

    echo "compilation done."
fi


#wait until prev. compilation is finished
echo "waiting..."
#until [ `ps -a | grep make | wc -l` -eq 0 ]; 
until [ `ls -a | grep $logpath/*.tempXYZ | wc -l` -eq 0 ];
do 
    sleep 2
done

rm "$CWD/tmpFile.txt"
