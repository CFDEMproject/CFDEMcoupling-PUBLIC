#!/bin/bash

#===================================================================#
# compile routine for CFDEMcoupling source, part of CFDEMproject 
# will create all lnInclude directories before compilation in order
# to avoid missing headers in foreign libraries
# Christoph Goniva - May. 2012, DCS Computing GmbH
# update: Stefan Radl (TU Graz, April 2016)
#===================================================================#
 
#- include functions
source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh

NOW="$(date +"%Y-%m-%d-%H:%M")"
logDir="log"


cd $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc
mkdir -p $logDir

rm $logDir/log_compile_results_src_success
rm $logDir/log_compile_results_src_fail

##================================================================================#
## Must compile (but not clean) LIGGGHTS libraries, since it could have been 
## compiled before with the compileLIGGGHTS command
## Then, check successful compilation
##================================================================================#
#bash $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/compileLIGGGHTS_lib.sh noClean
#echo "...now checking if LIGGGHTS libraries are compiled that are needed for CFDEM's src packages."
#bash $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/compileLIGGGHTS_lib.sh false

#================================================================================#
# compile src
#================================================================================#
whitelist="$CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/library-list.txt"
echo ""
echo "Please provide the libraries to be compiled in the $whitelist file."

if [ ! -f "$whitelist" ];then
    echo "$whitelist does not exist in $CWD. Nothing will be done."
    NLINES=0
    COUNT=0
else
    NLINES=`wc -l < $whitelist`
    COUNT=0
fi

#Generate lnIncludes, only for paths
while [ $COUNT -lt $NLINES ]
do
        let COUNT++  
        LINE=`head -n $COUNT $whitelist | tail -1`

        # white lines
        if [[ "$LINE" == "" ]]; then
            echo "compile $LINE"
            continue
        # comments
        elif [[ "$LINE" == \#* ]]; then
            continue
         # paths
        elif [[ "$LINE" == */dir ]]; then
            echo "will change path and create lnInclude..."
            LINE=$(echo "${LINE%????}")
            path="$CFDEM_SRC_DIR/$LINE"
            cd $path
            #continue
        fi
        wmakeLnInclude .
done
COUNT=0

echo
echo
echo "\n Creation of lnInclude directories finished!"
echo
echo

while [ $COUNT -lt $NLINES ]
do
        let COUNT++  
        LINE=`head -n $COUNT $whitelist | tail -1`

        # white lines
        if [[ "$LINE" == "" ]]; then
            echo "compile $LINE"
            continue
        # comments
        elif [[ "$LINE" == \#* ]]; then
            continue
         # paths
        elif [[ "$LINE" == */dir ]]; then
            echo "will change path..."
            LINE=$(echo "${LINE%????}")
            path="$CFDEM_SRC_DIR/$LINE"
            cd $path
            #continue
        fi
        #--------------------------------------------------------------------------------#
        #- define variables
        logpath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")/$logDir"
        logfileName="log_compileCFDEMcoupling_"$(basename $LINE)""
        casePath="$path"
        headerText="$logfileName""-$NOW"
        #--------------------------------------------------------------------------------#
        # remove old log file
        rm "$logpath/$logfileName"*
        compileLib $logpath $logfileName $casePath $headerText

        if [ ${PIPESTATUS[0]} -ne 0 ]; then         
            exit 1
        fi

        collectLogCFDEMcoupling_src $logpath $logfileName $casePath
        

done

