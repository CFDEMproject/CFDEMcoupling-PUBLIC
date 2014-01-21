#!/bin/bash

#===================================================================#
# compile routine for CFDEMcoupling source, part of CFDEMproject 
# Christoph Goniva - May. 2012, DCS Computing GmbH
# update: Stefan Radl (TU Graz, Jan 2014)
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
    whitelist="$CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/library-list.txt"
    echo ""
    echo "Please provide the libraries to be compiled in the $CWD/$whitelist file."

    if [ ! -f "$CWD/$whitelist" ];then
        echo "$whitelist does not exist in $CWD. Nothing will be done."
        NLINES=0
        COUNT=0
    else
        NLINES=`wc -l < $CWD/$whitelist`
        COUNT=0
    fi

    while [ $COUNT -lt $NLINES ]
    do
            let COUNT++  
            LINE=`head -n $COUNT $CWD/$whitelist | tail -1`
  
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
            compileLib $logpath $logfileName $casePath $headerText
    done

