#!/bin/bash

#===================================================================#
# test routine for cfdem project 
# defining functions used by the shell scripts
# Christoph Goniva - June. 2011, DCS Computing GmbH
#===================================================================#

#==================================#
#- function to pull from a repo

pullRepo()
{
    #--------------------------------------------------------------------------------#
    #- define variables
    logpath="$1"
    logfileName="$2"
    casePath="$3"
    headerText="$4"
    #--------------------------------------------------------------------------------#

    #- clean up old log file
    rm $logpath/$logfileName

    #- change path
    cd $casePath

    #- header
    echo 2>&1 | tee -a $logpath/$logfileName
    echo 2>&1 | tee -a $logpath/$logfileName
    echo "//=== $headerText ===//" 2>&1 | tee -a $logpath/$logfileName
    echo 2>&1 | tee -a $logpath/$logfileName

    #- write path
    pwd 2>&1 | tee -a $logpath/$logfileName
    echo 2>&1 | tee -a $logpath/$logfileName

    #- pull
    git pull 2>&1 | tee -a $logpath/$logfileName
}
#==================================#

#==================================#
#- function to compile a cfdem library

compileLib()
{
    #--------------------------------------------------------------------------------#
    #- define variables
    logpath="$1"
    logfileName="$2"
    casePath="$3"
    headerText="$4"
    #doClean="$5"
    #--------------------------------------------------------------------------------#

    #- clean up old log file
    rm $logpath/$logfileName

    #- change path
    cd $casePath

    #- header
    echo 2>&1 | tee -a $logpath/$logfileName
    echo "//   $headerText   //" 2>&1 | tee -a $logpath/$logfileName
    echo 2>&1 | tee -a $logpath/$logfileName

    #- write path
    pwd 2>&1 | tee -a $logpath/$logfileName
    echo 2>&1 | tee -a $logpath/$logfileName

    #- wclean and wmake
    #if [ $doClean != "noClean" ]; then
        # check library to compile is compressible
        str=$casePath
        i=$((${#str}-4))
        ending=${str:$i:4}
        if [[ $ending == "Comp" ]]; then
                echo "Compiling a compressible library - so doing an rmdepall of incomp library first."
                echo "Please make sure to have the incompressible libraries first in the library-list.txt!"
                cd $CFDEM_SRC_DIR/lagrangian/cfdemParticle
                echo "changing to $PWD"
                if [[ $WM_PROJECT_VERSION == "dev" ]]; then
                    wrmdep 2>&1 | tee -a $logpath/$logfileName
                else
                    rmdepall 2>&1 | tee -a $logpath/$logfileName
                fi
                cd $casePath
                echo "changing to $PWD"
            else
                echo "Compiling a incompressible library."
        fi
        if [[ $WM_PROJECT_VERSION == "dev" ]]; then
            wrmdep 2>&1 | tee -a $logpath/$logfileName
        else
            rmdepall 2>&1 | tee -a $logpath/$logfileName
        fi
        wclean 2>&1 | tee -a $logpath/$logfileName
    #fi
    wmake libso 2>&1 | tee -a $logpath/$logfileName

    #- keep terminal open
    #read
}
#==================================#

#==================================#
#- function to compile a cfdem solver

compileSolver()
{
    #--------------------------------------------------------------------------------#
    #- define variables
    logpath="$1"
    logfileName="$2"
    casePath="$3"
    headerText="$4"
    #doClean="$5"
    parallel="$5"
    #--------------------------------------------------------------------------------#

    #- clean up old log file
    rm $logpath/$logfileName

    #- change path
    cd $casePath

    #- header
    echo 2>&1 | tee -a $logpath/$logfileName
    echo "//   $headerText   //" 2>&1 | tee -a $logpath/$logfileName
    echo 2>&1 | tee -a $logpath/$logfileName

    #- write path
    pwd 2>&1 | tee -a $logpath/$logfileName
    echo 2>&1 | tee -a $logpath/$logfileName

    # check if there is an Allwmake
    if [ -f "Allwmake" ]; then
        echo "doing Allwclean and Allwmake for the solver"
        bash Allwclean 2>&1 | tee -a $logpath/$logfileName
        bash Allwmake 2>&1 | tee -a $logpath/$logfileName
    else
        #- wclean and wmake
        #if [ $doClean != "noClean" ]; then
            rmdepall 2>&1 | tee -a $logpath/$logfileName
            wclean 2>&1 | tee -a $logpath/$logfileName
        #fi
    fi
    
    # compile parallel?
    if [[ $parallel == "true" ]]; then
        wmake 2>&1 | tee -a $logpath/$logfileName &
    else
        wmake 2>&1 | tee -a $logpath/$logfileName
    fi

    #- keep terminal open
    #read
}
#==================================#

#==================================#
#- function to compile LIGGGHTS

compileLIGGGHTS()
{
    #--------------------------------------------------------------------------------#
    #- define variables
    logpath="$1"
    logfileName="$2"
    headerText="$3"
    clean="$4"
    #--------------------------------------------------------------------------------#

    #- clean up old log file
    rm $logpath/$logfileName

    #- change path
    cd $CFDEM_LIGGGHTS_SRC_DIR

    #- header
    echo 2>&1 | tee -a $logpath/$logfileName
    echo "//   $headerText   //" 2>&1 | tee -a $logpath/$logfileName
    echo 2>&1 | tee -a $logpath/$logfileName

    #- write path
    pwd 2>&1 | tee -a $logpath/$logfileName
    echo 2>&1 | tee -a $logpath/$logfileName

    #- wclean and wmake
    if [[ $clean == "false" ]]; then
        echo "not cleaning LIGGGHTS"
    else
        rm $CFDEM_LIGGGHTS_SRC_DIR/"lmp_"$CFDEM_LIGGGHTS_MAKEFILE_NAME
        rm $CFDEM_LIGGGHTS_SRC_DIR/"lib"$CFDEM_LIGGGHTS_LIB_NAME".a"
        make clean-$CFDEM_LIGGGHTS_MAKEFILE_NAME 2>&1 | tee -a $logpath/$logfileName
        echo "cleaning LIGGGHTS"
    fi
    if [[ $WM_NCOMPPROCS == "" ]]; then
        echo "compiling LIGGGHTS on one CPU"
        make $CFDEM_LIGGGHTS_MAKEFILE_NAME 2>&1 | tee -a $logpath/$logfileName
    else
        echo "compiling LIGGGHTS on $WM_NCOMPPROCS CPUs"
        make $CFDEM_LIGGGHTS_MAKEFILE_NAME -j $WM_NCOMPPROCS 2>&1 | tee -a $logpath/$logfileName
        #make $CFDEM_LIGGGHTS_MAKEFILE_NAME -j $WM_NCOMPPROCS yes-XYZ 2>&1 | tee -a $logpath/$logfileName
    fi
    make makelib 2>&1 | tee -a $logpath/$logfileName
    make -f Makefile.lib $CFDEM_LIGGGHTS_MAKEFILE_NAME 2>&1 | tee -a $logpath/$logfileName
}

#==================================#
#- function to compile a lammps lib

compileLMPlib()
{
    #--------------------------------------------------------------------------------#
    #- define variables
    logpath="$1"
    logfileName="$2"
    headerText="$3"
    makeFileName="$4"
    libraryPath="$5"
    #--------------------------------------------------------------------------------#

    #- clean up old log file
    rm $logpath/$logfileName

    #- change path
    if [ -d "$libraryPath" ]; then
        cd $libraryPath
    else
        echo ""
        echo "lib path $libraryPath does not exist - check settings in .../etc/bashrc."
        read
    fi

    #- header
    echo 2>&1 | tee -a $logpath/$logfileName
    echo "//   $headerText   //" 2>&1 | tee -a $logpath/$logfileName
    echo 2>&1 | tee -a $logpath/$logfileName

    #- write path
    pwd 2>&1 | tee -a $logpath/$logfileName
    echo 2>&1 | tee -a $logpath/$logfileName

    if [[ $makeFileName == "Makefile.Install" ]]; then
        echo "using Install.sh"
        bash Install.sh 0 2>&1 | tee -a $logpath/$logfileName
        bash Install.sh 1 2>&1 | tee -a $logpath/$logfileName
    else
        #- clean up
        echo "make clean" 2>&1 | tee -a $logpath/$logfileName
        echo 2>&1 | tee -a $logpath/$logfileName
        make -f $makeFileName clean 2>&1 | tee -a $logpath/$logfileName

        #- compile
        echo "make" 2>&1 | tee -a $logpath/$logfileName
        echo 2>&1 | tee -a $logpath/$logfileName
        make -f $makeFileName 2>&1 | tee -a $logpath/$logfileName
    fi
}
#==================================#

#==================================#
#- function to clean CFDEMcoupling solvers and src

cleanCFDEM()
{
    echo "do you really want to clean CFDEM src?"
    echo "if not, abort with ctrl-C"
    read

    #**********************************************
    #cleaning libraries
    whitelist="$CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/library-list.txt"
    echo ""
    echo "Please provide the libraries to be cleaned in the $CWD/$whitelist file."

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
                continue
            # comments
            elif [[ "$LINE" == \#* ]]; then
                continue
             # paths
            elif [[ "$LINE" == */dir ]]; then
                echo "will change path..."
                LINE=$(echo "${LINE%????}")
                path="$CFDEM_PROJECT_DIR/src/$LINE"
                cd $path
                #continue
            fi

            cd  $path
            echo "cleaning library $PWD"
            rmdepall
            wclean    
            rm -r ./Make/linux*
            rm -r ./lnInclude
    done


    #**********************************************
    #cleaning utilities
    echo "removing object files in"
    echo "   $CFDEM_UT_DIR"
    rm -r $CFDEM_UT_DIR/*/Make/linux*
    rm -r $CFDEM_UT_DIR/*/Make/linux*
    rm -r $CFDEM_UT_DIR/*/*.dep



    #**********************************************
    #cleaning solvers
    whitelist="$CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/solver-list.txt"
    echo ""
    echo "Please provide the solvers to be cleaned in the $CWD/$whitelist file."

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
                continue
            # comments
            elif [[ "$LINE" == \#* ]]; then
                continue
             # paths
            elif [[ "$LINE" == */dir ]]; then
                echo "will change path..."
                LINE=$(echo "${LINE%????}")
                path="$CFDEM_SOLVER_DIR/$LINE"
                cd $path
                #continue
            fi

            cd  $path            
            echo "cleaning solver $PWD"
            rmdepall
            wclean    
    done
}
#==================================#

#==================================#
#- function to clean CFDEMcoupling case

cleanCFDEMcase()
{
    #--------------------------------------------------------------------------------#
    #- define variables
    casepath="$1"
    #--------------------------------------------------------------------------------#

    echo "deleting data at: $casePath :\n"
    source $WM_PROJECT_DIR/bin/tools/CleanFunctions
    cd $casePath/CFD
    cleanCase
    rm -r $casePath/DEM/post/*
    echo "dummyfile" >> $casePath/DEM/post/dummy
    cd $casePath
    echo "done"
}

#==================================#
#- function to run a DEM case

DEMrun()
{
    #--------------------------------------------------------------------------------#
    #- define variables
    logpath="$1"
    logfileName="$2"
    casePath="$3"
    headerText="$4"
    solverName="$5"
    debugMode="$6"
    #--------------------------------------------------------------------------------#

    if [ $debugMode == "on" ]; then
        debugMode="valgrind"
    elif [ $debugMode == "strict" ]; then
        #debugMode="valgrind --leak-check=full -v --trace-children=yes --track-origins=yes" 
        debugMode="valgrind --tool=memcheck --track-origins=yes --leak-check=yes --show-reachable=yes --num-callers=20 --track-fds=yes"  
    else
        debugMode=""
    fi

    #- clean up old log file
    rm $logpath/$logfileName

    #- change path
    cd $casePath/DEM

    #- header
    echo 2>&1 | tee -a $logpath/$logfileName
    echo "//   $headerText   //" 2>&1 | tee -a $logpath/$logfileName
    echo 2>&1 | tee -a $logpath/$logfileName

    #- write path
    pwd 2>&1 | tee -a $logpath/$logfileName
    echo 2>&1 | tee -a $logpath/$logfileName

    #- run applictaion
    $debugMode $CFDEM_LIGGGHTS_SRC_DIR/$CFDEM_LIGGGHTS_LIB_NAME < $solverName 2>&1 | tee -a $logpath/$logfileName

    #- keep terminal open (if started in new terminal)
    #read
}
#==================================#

#==================================#
#- function to run a DEM case in parallel

parDEMrun()
{
    #--------------------------------------------------------------------------------#
    #- define variables
    logpath="$1"
    logfileName="$2"
    casePath="$3"
    headerText="$4"
    solverName="$5"
    nrProcs="$6"
    machineFileName="$7"
    debugMode="$8"
    #--------------------------------------------------------------------------------#

    if [ $debugMode == "on" ]; then
        debugMode="valgrind"
    elif [ $debugMode == "strict" ]; then
        #debugMode="valgrind --leak-check=full -v --trace-children=yes --track-origins=yes" 
        debugMode="valgrind --tool=memcheck --leak-check=yes --show-reachable=yes --num-callers=20 --track-fds=yes"  
    else
        debugMode=""
    fi

    #- clean up old log file
    rm $logpath/$logfileName

    #- change path
    cd $casePath/DEM

    #- header
    echo 2>&1 | tee -a $logpath/$logfileName
    echo "//   $headerText   //" 2>&1 | tee -a $logpath/$logfileName
    echo 2>&1 | tee -a $logpath/$logfileName

    #- write path
    pwd 2>&1 | tee -a $logpath/$logfileName
    echo 2>&1 | tee -a $logpath/$logfileName

    #- run applictaion
    if [ $machineFileName == "none" ]; then
        mpirun -np $nrProcs $debugMode $CFDEM_LIGGGHTS_SRC_DIR/$CFDEM_LIGGGHTS_LIB_NAME < $solverName 2>&1 | tee -a $logpath/$logfileName
    else
        mpirun -machinefile $machineFileName -np $nrProcs $debugMode $CFDEM_LIGGGHTS_SRC_DIR/$CFDEM_LIGGGHTS_LIB_NAME < $solverName 2>&1 | tee -a $logpath/$logfileName
    fi

    #- keep terminal open (if started in new terminal)
    #read
}
#==================================#

#==================================#
#- function to run a CFD case

CFDrun()
{
    #--------------------------------------------------------------------------------#
    #- define variables
    logpath="$1"
    logfileName="$2"
    casePath="$3"
    headerText="$4"
    solverName="$5"
    debugMode="$6"
    #--------------------------------------------------------------------------------#

    if [ $debugMode == "on" ]; then
        debugMode="valgrind"
    else
        debugMode=""
    fi

    #- clean up old log file
    rm $logpath/$logfileName

    #- change path
    cd $casePath/CFD

    #- header
    echo 2>&1 | tee -a /$logpath/$logfileName
    echo "//   $headerText   //" 2>&1 | tee -a $logpath/$logfileName
    echo 2>&1 | tee -a $logpath/$logfileName

    #- write path
    pwd 2>&1 | tee -a $logpath/$logfileName
    echo 2>&1 | tee -a $logpath/$logfileName

    #- clean up case
    #rm couplingFiles/*

    #- run applictaion
    $debugMode $solverName 2>&1 | tee -a $logpath/$logfileName

    #- keep terminal open (if started in new terminal)
    #read
}
#==================================#

#==================================#
#- function to run a CFD case in parallel   !!!NOT DEBUGGED!!!
parCFDrun()
{
    #--------------------------------------------------------------------------------#
    #- define variables
    logpath="$1"
    logfileName="$2"
    casePath="$3"
    headerText="$4"
    solverName="$5"
    nrProcs="$6"
    machineFileName="$7"
    debugMode="$8"
    #--------------------------------------------------------------------------------#

    if [ $debugMode == "on" ]; then
        debugMode="valgrind"
    else
        debugMode=""
    fi

    #- clean up old log file
    rm $logpath/$logfileName

    #- change path
    cd $casePath

    #- remove old data
    rm -rf processor*

    #- decompose case
    decomposePar

    #- make proc dirs visible
    count=0
    for i in `seq $nrProcs`
    do
        let count=$i-1
        (cd $casePath/processor$count && touch file.foam)
    done

    #- header
    echo 2>&1 | tee -a /$logpath/$logfileName
    echo "//   $headerText   //" 2>&1 | tee -a $logpath/$logfileName
    echo 2>&1 | tee -a $logpath/$logfileName

    #- write path
    pwd 2>&1 | tee -a $logpath/$logfileName
    echo 2>&1 | tee -a $logpath/$logfileName

    #- run applictaion
    if [ $machineFileName == "none" ]; then
        mpirun -np $nrProcs $debugMode $solverName -parallel 2>&1 | tee -a $logpath/$logfileName
    else
        mpirun -machinefile $machineFileName -np $nrProcs $debugMode $solverName -parallel 2>&1 | tee -a $logpath/$logfileName
    fi

    #- keep terminal open (if started in new terminal)
    #read
}
#==================================#

#==================================#
#- function to run a parallel CFD-DEM case

parCFDDEMrun()
{
    #--------------------------------------------------------------------------------#
    #- define variables
    logpath="$1"
    logfileName="$2"
    casePath="$3"
    headerText="$4"
    solverName="$5"
    nrProcs="$6"
    machineFileName="$7"
    debugMode="$8"
    reconstuctCase="$9"
    cleanCase="$10"
    #--------------------------------------------------------------------------------#

    if [ $debugMode == "on" ]; then
        debugMode="valgrind"
    elif [ $debugMode == "strict" ]; then
        #debugMode="valgrind --leak-check=full -v --trace-children=yes --track-origins=yes" 
        debugMode="valgrind --tool=memcheck --leak-check=yes --show-reachable=yes --num-callers=20 --track-fds=yes"
    elif [ $debugMode == "profile" ]; then
        debugMode="hpcrun"
        rm -r $casePath/CFD/hpctoolkit-$solverName-measurements
    else
        debugMode=""
    fi

    #- clean up old log file
    rm $logpath/$logfileName

    #- change path
    cd $casePath/CFD

    #- remove old data
    rm -rf processor*

    #- decompose case
    decomposePar

    #- make proc dirs visible
    count=0
    for i in `seq $nrProcs`
    do
        let count=$i-1
        (cd $casePath/CFD/processor$count && touch file.foam)
    done

    #- header
    echo 2>&1 | tee -a /$logpath/$logfileName
    echo "//   $headerText   //" 2>&1 | tee -a $logpath/$logfileName
    echo 2>&1 | tee -a $logpath/$logfileName

    #- write path
    pwd 2>&1 | tee -a $logpath/$logfileName
    echo 2>&1 | tee -a $logpath/$logfileName

    #- clean up case
    rm couplingFiles/*

    #- run applictaion
    if [[ $machineFileName == "none" ]]; then
        mpirun -np $nrProcs $debugMode $solverName -parallel 2>&1 | tee -a $logpath/$logfileName

        #- reconstruct case
        if [[ $reconstuctCase == "true" ]]; then   
            #pseudoParallelRun "reconstructPar" $nrProcs
            reconstructPar
        fi
    else
        mpirun -machinefile $machineFileName -np $nrProcs $debugMode $solverName -parallel 2>&1 | tee -a $logpath/$logfileName

        #- reconstruct case
        if [[ $reconstuctCase == "true" ]]; then   
            #pseudoParallelRun "reconstructPar" $nrProcs
            reconstructPar
        fi
    fi

    if [[ $debugMode == "hpcrun" ]]; then
        rm hpctoolkit-$solverName-database*
        hpcprof -S $CFDEM_APP_DIR/$solverName.hpcstruct -I ./'*' hpctoolkit-$solverName-measurements
        echo "huhu"      
    fi

    #- keep terminal open (if started in new terminal)
    #read
}
#==================================#


#==================================#
#- function to collect results from 
#- logfiles to one log file

collectLog()
{
    #--------------------------------------------------------------------------------#
    #- define variables
    logpath="$1"
    logfileName="$2"
    casePath="$3"
    headerText="$4"
    nrOfLines="$5"
    #--------------------------------------------------------------------------------#

    echo  2>&1 | tee -a $logpath/$logfileName
    echo "//   $headerText   //" 2>&1 | tee -a $logpath/$logfileName
    tail --lines=$nrOfLines $casePath |cut -d " " -f1- 2>&1 | tee -a $logpath/$logfileName
}
#==================================#

#==================================#
#- function to collect results from 
#- logfiles to one log file

collectLogCFDEMcoupling_sol()
{
    #--------------------------------------------------------------------------------#
    #- define variables
    logpath="$1"
    logfileName="$2"
    casePath="$3"
    #--------------------------------------------------------------------------------#
    # read name of solver
    SOLVERNAME=$(basename $casePath)
    
    # read last line of log
    LASTLINE=`tac $logpath/$logfileName | egrep -m 1 .`
    LASTSTRING=`echo ${LASTLINE##* }`
    LASTWORD=$(basename $LASTSTRING)

    # log if compilation was success  
    if [[ $LASTWORD == $SOLVERNAME || $LASTWORD == "date." ]]; then
        echo "$SOLVERNAME" >> $logpath/log_compile_results_success
    else
        echo "$SOLVERNAME" >> $logpath/log_compile_results_fail
    fi
}
#==================================#

#==================================#
#- function to replace a line in a file where text consecutive
#  the old line must look like: oldWord
#  and will be replaced by: newWord

replaceLineInFile()
{
    #--------------------------------------------------------------------------------#
    #- define variables
    logpath="$1"
    logfileName="$2"
    filePath="$3"
    oldWord="$4"       # give text only until first "/"
    newWord="$5"
    #--------------------------------------------------------------------------------#

    #- adapt /etc/bashrc
    echo "replaceLineInFile $filePath" 2>&1 | tee -a $logpath/$logfileName
    sed "/$oldWord/ c\ $newWord" $filePath > $filePath"2"  2>&1 | tee -a $logpath/$logfileName
    cp $filePath"2" $filePath  2>&1 | tee -a $logpath/$logfileName
    rm $filePath"2"  2>&1 | tee -a $logpath/$logfileName
}
#==================================#

#==================================#
#- function to replace a line in a file where text is separated by one blank
#  the old line must look like: oldWord1 oldWord2
#  and will be replaced by: newWord1 newWord2

replaceSeparatedLineInFile()
{
    #--------------------------------------------------------------------------------#
    #- define variables
    logpath="$1"
    logfileName="$2"
    filePath="$3"
    oldWord1="$4"
    oldWord2="$5"
    newWord1="$6"
    newWord2="$7"
    #--------------------------------------------------------------------------------#

    #- adapt /etc/bashrc
    echo "replaceLineInFile $filePath" 2>&1 | tee -a $logpath/$logfileName
    sed "/$oldWord1 $oldWord2/ c\ $newWord1 $newWord2" $filePath > $filePath"2"  2>&1 | tee -a $logpath/$logfileName
    cp $filePath"2" $filePath  2>&1 | tee -a $logpath/$logfileName
    rm $filePath"2"  2>&1 | tee -a $logpath/$logfileName
}
#==================================#

#=======================================#
#- script to run a function in pseudo-parallel
#  several runs of the function are started 
#  simultanously. Only makes sense on shared memory sys
#  based on script by K. Wardle 6/22/09
#  published at CFDonline forum

pseudoParallelRun()
{
    #--------------------------------------------------------------------------------#
    #- define variables
    appname=$1
    njobs=$2
    outputfile="log.$appname"
    #--------------------------------------------------------------------------------#

if [ -z "$njobs" ]; then 
   echo ""
   echo "  K. Wardle 6/22/09"
   echo "  bash script to run reconstructPar (or other fct) in pseudo-parallel mode"
   echo "  by breaking time directories into multiple ranges"
   echo "" 
   echo "  USAGE: appname <NP> <outputfile>   [output file is optional] "
   echo ""
   exit
fi


appflag="-noZero"

#let njobs1=$njobs+1  
echo "running $appname $appflag in pseudo-parallel mode on $njobs processors"

#count the number of time directories
nsteps=`ls -d processor0/[0-9]*/ | wc -l`
echo "do $appname on $nsteps time directories"
##nsteps=`ls processor0 | wc -l`
#echo "nsteps= $nsteps"
#let nsteps=$nsteps1-1

nchunk=`echo $nsteps/$njobs+1 | bc`
#echo "nchunk = $nchunk"

#find max time
tmin=`ls processor0 | sort -nr | tail -1`
#echo "tmin = $tmin"
tmax=`ls processor0 | sort -nr | head -1`
#echo "tmax = $tmax"

echo "making temp dir"
tempdir="temp.par$appname"
mkdir $tempdir

tstart=$tmin
p=p

for i in `seq $njobs`
do
  nn=`echo $i*$nchunk | bc`
  tstop=`ls processor0 | sed -n $nn$p`
  if [ $i == $njobs ] 
	then
	tstop=$tmax
  fi
  
  echo "Starting Job $i - $appname time = $tstart through $tstop"
  `$appname $appflag -time "$tstart:$tstop" > $tempdir/output-$tstop &`
  
  let nn=$nn+1
  tstart=`ls processor0 | sed -n $nn$p`  
done

#sleep until jobs finish
#if number of jobs > njobs, hold loop until job finishes
nmore_old=`echo 0`
until [ `ps -C $appname | wc -l` -eq 1 ]; 
  do 
    sleep 10
##    n1=`ps -C $appname | wc -l`
##    n2=`echo $n1-1 | bc`
    nnow=`ls -d [0-9]*/ | wc -l` ##count time directories in case root dir, this will include 0
    nmore=`echo $nsteps-$nnow+1 | bc` ##calculate number left to reconstruct and subtract 0 dir
    if [ $nmore != $nmore_old ]
      then
      echo "$nmore directories remaining..."
    fi
    nmore_old=$nmore
  done

#combine and cleanup
if [ -n "$outputfile" ] 
  then
#check if output file already exists
  if [ -e "$outputfile" ] 
  then
    echo "output file $outputfile exists, moving to $outputfile.bak"
    mv $outputfile $outputfile.bak
  fi

  echo "cleaning up temp files"
  for i in `ls $tempdir`
  do
    cat $tempdir/$i >> $outputfile
  done
fi

rm -rf $tempdir

echo "finished"

}
#==================================#
#- function make a tar.gz copy with date tag from a directory
#  Remove the original directory

backupRemoveDir()
{
    #--------------------------------------------------------------------------------#
    #- define variables
    filePath="$1"
    backupPath="$2"
    date="$3"
    #--------------------------------------------------------------------------------#
    echo "creating dirctory :\n $filePath""-until""$date\n"
    mkdir $filePath"-until"$date
    echo "move : $filePath/* to\n $filePath""-until""$date\n and tar.gz"
    mv $filePath/* $filePath"-until"$date
    tar czvf  $filePath"-until"$date.tar.gz $filePath"-until"$date
    rm -r $filePath"-until"$date
    mv $filePath"-until"$date.tar.gz $backupPath
}

#==================================#
#- function make a tar.gz copy with date tag from a directory
#  Keep the original directory

backupDir()
{
    #--------------------------------------------------------------------------------#
    #- define variables
    filePath="$1"
    backupPath="$2"
    date="$3"
    #--------------------------------------------------------------------------------#
    echo "creating dirctory :\n $filePath""-until""$date\n"
    mkdir $filePath"-until"$date
    echo "move : $filePath/* to\n $filePath""-until""$date\n and tar.gz"
    cp -r $filePath/* $filePath"-until"$date
    tar czvf  $filePath"-until"$date.tar.gz $filePath"-until"$date
    rm -r $filePath"-until"$date
    mv $filePath"-until"$date.tar.gz $backupPath
}

#========================================#
#- function to check if a directory exits
checkDir()
{
    #--------------------------------------------------------------------------------#
    #- define variables
    filePath="$1"
    #--------------------------------------------------------------------------------#
    if [ -d "$filePath" ]; then
         echo "true" 
    else
        echo "false" 
    fi
}

#========================================#
#- function to check if a directory exits
checkDirComment()
{
    #--------------------------------------------------------------------------------#
    #- define variables
    filePath="$1"
    varName="$2"
    critical="$3"
    #--------------------------------------------------------------------------------#
    if [ $(checkDir $filePath) == "true" ]; then
         echo "valid:yes critical:$critical - $varName = $filePath" 
    else
        echo "valid:NO  critical:$critical - $varName = $filePath does not exist" 
    fi
}

#========================================#
#- function to check if a variable exits
checkEnv()
{
    #--------------------------------------------------------------------------------#
    #- define variables
    var="$1"
    #--------------------------------------------------------------------------------#
    if [[ $var == "" ]]; then
        echo "false"
    else
        echo "true"
    fi
}

#========================================#
#- function to check if a variable exits
checkEnvComment()
{
    #--------------------------------------------------------------------------------#
    #- define variables
    var="$1"
    varName="$2"
    critical="$3"
    #--------------------------------------------------------------------------------#
    if [ $(checkEnv $var) == "true" ]; then
         echo "valid:yes critical:$critical - $varName = $var" 
    else
        echo "valid:NO  critical:$critical - $varName = $var variable not set!" 
    fi
}

#========================================#
#- function to print a header to terminal
printHeader()
{
    echo ""
    echo "*********************************************"
    echo "* C F D E M (R) c o u p l i n g             *"
    echo "*                                           *"
    echo "* by DCS Computing GmbH                     *"
    echo "* www.dcs-computing.com                     *"
    echo "*********************************************"
    echo "" 
}

#========================================#
#- track memory usage
trackMem()
{

    #--------------------------------------------------------------------------------#
    #- define variables
    appname="$1"
    fileName="$2"
    #--------------------------------------------------------------------------------#

    rm $fileName

    echo "please use only the the first 15  strings of the command !!!"

    /usr/bin/printf "%-6s %-9s %s\n" "PID" "Total" "Command" >> $fileName
    /usr/bin/printf "%-6s %-9s %s\n" "---" "-----" "-------" >> $fileName

    for PID in $(/bin/ps -e | /usr/bin/awk '$1 ~ /[0-9]+/ { print $1 }')
    do
        CMD=$(/bin/ps -o comm -p $PID | /usr/bin/tail -1)

        if [ $CMD == $appname ]
            then

            TOTAL=$(/usr/bin/pmap $PID 2>/dev/null | /usr/bin/tail -1 | /usr/bin/awk '{ print $2 }')
            [ -n "$TOTAL" ] && /usr/bin/printf "%-6s %-9s %s\n" "$PID" "$TOTAL" "$CMD"
        fi
    done | /usr/bin/sort -n -k2 >> $fileName
}

#========================================#
#- remove brackets from file
removeBracketsFromFile()
{
    #--------------------------------------------------------------------------------#
    #- define variables
    oldFileName="$1"
    newFileName="$oldFileName""_noBrackets"
    #--------------------------------------------------------------------------------#

    sed -e 's/[(,)]//g' $oldFileName > $newFileName
}



