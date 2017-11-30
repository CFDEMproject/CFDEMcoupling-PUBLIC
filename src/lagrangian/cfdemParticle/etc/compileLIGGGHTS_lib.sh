#!/bin/bash

#===================================================================#
# compile routine for LIGGGHTS libraries, part of CFDEMproject 
# Christoph Goniva - March. 2014, DCS Computing GmbH
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
    whitelist="$CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/package-undo-liggghts-list.txt"
    echo ""
    echo "==========================================="
    echo "deactivating all possible packages of LIGGGHTS now..."
    echo "Please provide the packages to be compiled in the $CWD/$whitelist file."
    echo "Packages must be in: $CFDEM_LAMMPS_LIB_DIR."

    if [ ! -f "$CWD/$whitelist" ];then
        echo "$whitelist does not exist in $CWD. Nothing will be done."
        NLINES=0
        COUNT=0
    else
        NLINES=`wc -l < $CWD/$whitelist`
        COUNT=0
    fi

    logpath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")/$logDir"

    # resetting Makefile.package
    cp $CFDEM_LIGGGHTS_SRC_DIR/Makefile.package.empty $CFDEM_LIGGGHTS_SRC_DIR/Makefile.package
    cp $CFDEM_LIGGGHTS_SRC_DIR/Makefile.package.settings.empty $CFDEM_LIGGGHTS_SRC_DIR/Makefile.package.settings

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
                path="$CFDEM_LIGGGHTS_SRC_DIR"
                cd $path
                echo $PWD
                #continue
            fi

            #--------------------------------------------------------------------------------#
            #- define variables
            logfileName="log_compile$LINE""lib"
            #--------------------------------------------------------------------------------#
            rm $logfileName 2>&1 | tee -a $logpath/$logfileName

            make no-$LINE 2>&1 | tee -a $logpath/$logfileName
    done

    #===================================================================================
    whitelist="$CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/package-liggghts-list.txt"
    echo ""
    echo "==========================================="
    echo "activating packages of LIGGGHTS now..."
    echo "Please provide the packages to be compiled in the $CWD/$whitelist file."
    echo "Packages must be in: $CFDEM_LAMMPS_LIB_DIR."

    if [ ! -f "$CWD/$whitelist" ];then
        echo "$whitelist does not exist in $CWD. Nothing will be done."
        NLINES=0
        COUNT=0
    else
        NLINES=`wc -l < $CWD/$whitelist`
        COUNT=0
    fi

    logpath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")/$logDir"

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
                path="$CFDEM_LIGGGHTS_SRC_DIR"
                cd $path
                echo $PWD
                #continue
            fi

            #--------------------------------------------------------------------------------#
            #- define variables
            logfileName="log_compile$LINE""lib"
            #--------------------------------------------------------------------------------#

            rm $logpath/$logfileName
            make yes-$LINE 2>&1 | tee -a $logpath/$logfileName

            # assuming we need the poems lib if the package POEMS is activated
            if [[ "$LINE" == "POEMS" ]]; then
                echo "compile $LINE"
                cd $CFDEM_POEMSLIB_PATH
                make -f Makefile.g++ clean 2>&1 | tee -a $logpath/$logfileName
                make -j $nProc -f Makefile.g++ lib 2>&1 | tee -a $logpath/$logfileName
                cd $path
            fi

            # special handling of PARSCALE
            if [[ "$LINE" == "PASCAL" ]]; then
                echo "compile $LINE"
                #. $PASCAL_SRC_DIR/refresh 2>&1 | tee -a $logpath/$logfileName
                . $PASCAL_SRC_DIR/refreshLibrary.sh 2>&1 | tee -a $logpath/$logfileName
            fi
    done
