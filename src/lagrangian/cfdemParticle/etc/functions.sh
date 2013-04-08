#!/bin/bash

#===================================================================#
# test routine for cfdem project 
# defining functions used by the shell scripts
# Christoph Goniva - June. 2011, DCS Computing GmbH
#===================================================================#


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
        rmdepall 2>&1 | tee -a $logpath/$logfileName
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
        rmdepall 2>&1 | tee -a $logpath/$logfileName
        wclean 2>&1 | tee -a $logpath/$logfileName
    #fi
    wmake 2>&1 | tee -a $logpath/$logfileName

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
    make clean-all 2>&1 | tee -a $logpath/$logfileName
    make $CFDEM_LIGGGHTS_MAKEFILE_NAME -j $WM_NCOMPPROCS  2>&1 | tee -a $logpath/$logfileName
    make makelib 2>&1 | tee -a $logpath/$logfileName
    make -f Makefile.lib $CFDEM_LIGGGHTS_MAKEFILE_NAME 2>&1 | tee -a $logpath/$logfileName
}
#==================================#

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
    cd $libraryPath

    #- header
    echo 2>&1 | tee -a $logpath/$logfileName
    echo "//   $headerText   //" 2>&1 | tee -a $logpath/$logfileName
    echo 2>&1 | tee -a $logpath/$logfileName

    #- write path
    pwd 2>&1 | tee -a $logpath/$logfileName
    echo 2>&1 | tee -a $logpath/$logfileName

    #- clean up
    echo "make clean" 2>&1 | tee -a $logpath/$logfileName
    echo 2>&1 | tee -a $logpath/$logfileName
    make -f $makeFileName clean 2>&1 | tee -a $logpath/$logfileName

    #- compile
    echo "make" 2>&1 | tee -a $logpath/$logfileName
    echo 2>&1 | tee -a $logpath/$logfileName
    make -f $makeFileName 2>&1 | tee -a $logpath/$logfileName
}
#==================================#

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
    rm couplingFiles/*

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
    cd $casePath/CFD

    #- header
    echo 2>&1 | tee -a /$logpath/$logfileName
    echo "//   $headerText   //" 2>&1 | tee -a $logpath/$logfileName
    echo 2>&1 | tee -a $logpath/$logfileName

    #- write path
    pwd 2>&1 | tee -a $logpath/$logfileName
    echo 2>&1 | tee -a $logpath/$logfileName

    #- run applictaion
    if [ $machineFileName == "none" ]; then
        mpirun -np $nrProcs $solverName 2>&1 | tee -a $logpath/$logfileName
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
    cd $casePath/CFD

    #- remove old data
    rm -rf processor*

    #- decompose case
    decomposePar

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
    if [ $machineFileName == "none" ]; then
        mpirun -np $nrProcs $debugMode $solverName -parallel 2>&1 | tee -a $logpath/$logfileName

        #- reconstruct case
        #pseudoParallelRun "reconstructPar" $nrProcs
        reconstructPar
    else
        mpirun -machinefile $machineFileName -np $nrProcs $debugMode $solverName -parallel 2>&1 | tee -a $logpath/$logfileName

        #- reconstruct case
        reconstructPar
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



