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
                if [[ $WM_PROJECT_VERSION == dev || $WM_PROJECT_VERSION == 3.0.* || $WM_PROJECT_VERSION == 4.*  || $WM_PROJECT_VERSION == 5.* ]]; then
                    wrmdep 2>&1 | tee -a $logpath/$logfileName
                else
                    rmdepall 2>&1 | tee -a $logpath/$logfileName
                fi
                cd $casePath
                echo "changing to $PWD"
            else
                echo "Compiling a incompressible library."
        fi
        if [[ $WM_PROJECT_VERSION == dev || $WM_PROJECT_VERSION == 3.0.* || $WM_PROJECT_VERSION == 4.*  || $WM_PROJECT_VERSION == 5.* ]]; then
            wrmdep 2>&1 | tee -a $logpath/$logfileName
        else
            rmdepall 2>&1 | tee -a $logpath/$logfileName
        fi
        wclean 2>&1 | tee -a $logpath/$logfileName
        rm -r $casePath/Make/cfdemParticle
    #fi
    wmake libso 2>&1 | tee -a $logpath/$logfileName

    if [ ${PIPESTATUS[0]} -ne 0 ]; then 
        return 1
    fi

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
            if [[ $WM_PROJECT_VERSION == dev || $WM_PROJECT_VERSION == 3.0.* || $WM_PROJECT_VERSION == 4.*  || $WM_PROJECT_VERSION == 5.* ]]; then
                wrmdep 2>&1 | tee -a $logpath/$logfileName
            else
                rmdepall 2>&1 | tee -a $logpath/$logfileName
            fi
            wclean 2>&1 | tee -a $logpath/$logfileName
        #fi
    fi
    
    # compile parallel?
    if [[ $parallel == "true" ]]; then
        touch $logpath/$logfileName.tempXYZ
        wmake 2>&1 | tee -a $logpath/$logfileName #&
        rm $logpath/$logfileName.tempXYZ
    else
        touch $logpath/$logfileName.tempXYZ
        wmake 2>&1 | tee -a $logpath/$logfileName
        rm $logpath/$logfileName.tempXYZ

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

    # LIGGGHTS compilation flags
    if [[ $WM_NCOMPPROCS == "" ]]; then
         LIG_COMP_FLAGS=" "
    else
         LIG_COMP_FLAGS="-j $WM_NCOMPPROCS "
    fi

    if [[ ${WM_COMPILE_OPTION} == "Debug" ]] && [[ $CFDEM_LIGGGHTS_MAKEFILE_NAME == "auto" ]]; then
        LIG_COMP_FLAGS="${LIG_COMP_FLAGS} debug=FULL"
    fi

    #- clean and make
    if [[ $clean == "false" ]]; then
        echo "not cleaning LIGGGHTS"
    else
        rm $CFDEM_LIGGGHTS_SRC_DIR/$CFDEM_LIGGGHTS_LIB_NAME
        make clean-$CFDEM_LIGGGHTS_MAKEFILE_NAME postfix=$CFDEM_LIGGGHTS_MAKEFILE_POSTFIX ${LIG_COMP_FLAGS} 2>&1 | tee -a $logpath/$logfileName
        rm $CFDEM_LIGGGHTS_SRC_DIR/"lib"$CFDEM_LIGGGHTS_LIB_NAME".a"
        echo "cleaning LIGGGHTS"
    fi
    
    echo "compiling LIGGGHTS with ${LIG_COMP_FLAGS}"
    make $CFDEM_LIGGGHTS_MAKEFILE_NAME postfix=$CFDEM_LIGGGHTS_MAKEFILE_POSTFIX ${LIG_COMP_FLAGS} 2>&1 | tee -a $logpath/$logfileName
    make makeshlib 2>&1 | tee -a $logpath/$logfileName
    make -f Makefile.shlib $CFDEM_LIGGGHTS_MAKEFILE_NAME postfix=$CFDEM_LIGGGHTS_MAKEFILE_POSTFIX ${LIG_COMP_FLAGS} 2>&1 | tee -a $logpath/$logfileName

    if [ ${PIPESTATUS[0]} -ne 0 ]; then 
        return 1
    fi
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
    compilationModeSwitch="$6"
    #--------------------------------------------------------------------------------#

    #- change path
    if [ -d "$libraryPath" ]; then
        cd $libraryPath
    else
        echo ""
        echo "lib path $libraryPath does not exist - check settings in .../etc/bashrc."
        read
    fi

    # LIGGGHTS compilation flags
    if [[ $WM_NCOMPPROCS == "" ]]; then
         LIG_COMP_FLAGS=" "
    else
         LIG_COMP_FLAGS="-j $WM_NCOMPPROCS "
    fi

    if [[ ${WM_COMPILE_OPTION} == "Debug" ]] && [[ $CFDEM_LIGGGHTS_MAKEFILE_NAME == "auto" ]]; then
        LIG_COMP_FLAGS="${LIG_COMP_FLAGS} debug=FULL"
    fi

    # Fallback if $makeFileName does not exist
    if [[ ! -f "${makeFileName}" ]]; then
      echo "WARNING: userdefined Makefile ${makeFileName} cannot be found in ${libraryPath}, using default Makefile.mpi"
      makeFileName="Makefile.mpi"
    fi

    #Just check if library is there and and abort if not
    if [[ $compilationModeSwitch == "false" ]]; then
          if [ -d "$libraryPath" ]; then
                    echo "lib path $libraryPath EXISTS!"
                    libraries=$(ls ${libraryPath}/*.a 2> /dev/null | wc -l)
                    if [[ $libraries != "0" ]]; then
                        echo "... and contains the following libraries: "
                        ls $libraryPath/*.a
                        echo "Congratulations! Check passed! "
                    else
                        echo ""
                        echo "ERROR!!"
                        echo "... could not find Library!"
                        echo "... it should contain a *.a file to be linked to the final application."
                        echo "You need to ensure all libaries in this path are compiled."
                        echo "$libraryPath"
                        echo "are compiled. Therefore, you may want to use an appropriate script that compiles these libraries first. (e.g. cfdemCompLIGlibs)"
                        read
                    fi
          fi
    else
        if [[ $compilationModeSwitch == "noClean" ]]; then
            echo "compileLMPlib will skip the cleaning step!"
            echo ""
        else
            #Clean and then compile library
            #- clean up old log file
            rm $logpath/$logfileName
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
            if [[ $compilationModeSwitch == "noClean" ]]; then
                echo "compileLMPlib will skip the cleaning step!"
                echo ""
            else
                #- clean up
                echo "make clean" 2>&1 | tee -a $logpath/$logfileName
                echo 2>&1 | tee -a $logpath/$logfileName
                make -f $makeFileName clean 2>&1 | tee -a $logpath/$logfileName
            fi

            #- compile
            echo "compiling LIB with ${LIG_COMP_FLAGS}"
            echo "make" 2>&1 | tee -a $logpath/$logfileName
            echo 2>&1 | tee -a $logpath/$logfileName
            make -f $makeFileName ${LIG_COMP_FLAGS} 2>&1 | tee -a $logpath/$logfileName

            if [ ${PIPESTATUS[0]} -ne 0 ]; then 
                return 1
            fi
        fi
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
    echo "Please provide the libraries to be cleaned in the $whitelist file."

    if [ ! -f "$whitelist" ];then
        echo "$whitelist does not exist in $CWD. Nothing will be done."
        NLINES=0
        COUNT=0
    else
        NLINES=`wc -l < $whitelist`
        COUNT=0
    fi

    while [ $COUNT -lt $NLINES ]
    do
            let COUNT++  
            LINE=`head -n $COUNT $whitelist | tail -1`
  
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
            if [[ $WM_PROJECT_VERSION == dev || $WM_PROJECT_VERSION == 3.0.* || $WM_PROJECT_VERSION == 4.* || $WM_PROJECT_VERSION == 5.* ]]; then
                wrmdep
            else
                rmdepall
            fi
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
    echo "Please provide the solvers to be cleaned in the $whitelist file."

    if [ ! -f "$whitelist" ];then
        echo "$whitelist does not exist in $CWD. Nothing will be done."
        NLINES=0
        COUNT=0
    else
        NLINES=`wc -l < $whitelist`
        COUNT=0
    fi

    while [ $COUNT -lt $NLINES ]
    do
            let COUNT++  
            LINE=`head -n $COUNT $whitelist | tail -1`
  
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
            if [[ $WM_PROJECT_VERSION == dev || $WM_PROJECT_VERSION == 3.0.* || $WM_PROJECT_VERSION == 4.* || $WM_PROJECT_VERSION == 5.* ]]; then
                wrmdep
            else
                rmdepall
            fi
            wclean    
    done

    #**********************************************
    #cleaning logs
    rm $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/log/log_*
}
#==================================#

#==================================#
#- function to clean CFDEMcoupling case

cleanCFDEMcase()
{
    #--------------------------------------------------------------------------------#
    #- define variables
    casepath="$1"
    keepDEMrestart="$2"
    keepCFDmesh="$3"
    #--------------------------------------------------------------------------------#

    echo "deleting data at: $casePath ? otherwise press Ctrl-C:\n"
    read
    source $WM_PROJECT_DIR/bin/tools/CleanFunctions
    #CFD
    cd $casePath/CFD
    if [[ $keepCFDmesh == true ]]; then
        echo "keeping CFD mesh files"
        cp -r constant/polyMesh constant/polyMesh_backup
        cleanCase
        mv constant/polyMesh_backup/* constant/polyMesh
        rm -r constant/polyMesh_backup
    else
        echo "deleting CFD mesh files"
        cleanCase
    fi

    #CFDEM
    rm -r $casePath/CFD/clockData
    rm -r $casePath/CFD/particleProbes
    rm -r $casePath/CFD/averageProps/
    rm -r $casePath/CFD/octave/octave-core
    rm -r $casePath/CFD/octave/octave-workspace
    rm -r $casePath/remotePlace
    rm -r $casePath/CFD/oldProcDirs
    rm -r $casePath/CFD/tmp.balance
    rm $casePath/CFD/callgrind.out.*
    rm -r $casePath/CFD/hpctoolkit-*
    rm  $casePath/log_*
    #DEM
    rm $casePath/DEM/post/*
    touch $casePath/DEM/post/.gitignore
    if [[ $keepDEMrestart == true ]]; then
        echo "keeping DEM restart files"
    else
    rm  $casePath/DEM/post/restart/*
    fi
    touch $casePath/DEM/post/restart/.gitignore
    rm  $casePath/DEM/tmp.lammps.variable
    rm  $casePath/DEM/log*
    #ParScale
    rm $casePath/CFD/*.dat
    rm $casePath/CFD/*.pascal
    rm $casePath/CFD/*.profile
    rm -r $casePath/CFD/pascal/0.*
    rm -r $casePath/CFD/pascal/1
    rm -r $casePath/CFD/pascal/1.*
    rm -r $casePath/CFD/pascal/2
    rm -r $casePath/CFD/pascal/2.*
    rm -r $casePath/CFD/pascal/3
    rm -r $casePath/CFD/pascal/3.*
    rm -r $casePath/CFD/pascal/4
    rm -r $casePath/CFD/pascal/4.*
    rm -r $casePath/CFD/pascal/5
    rm -r $casePath/CFD/pascal/5.*
    rm -r $casePath/CFD/pascal/6
    rm -r $casePath/CFD/pascal/6.*
    rm -r $casePath/CFD/pascal/7
    rm -r $casePath/CFD/pascal/7.*
    rm -r $casePath/CFD/pascal/8
    rm -r $casePath/CFD/pascal/8.*
    rm -r $casePath/CFD/pascal/9
    rm -r $casePath/CFD/pascal/9.*
    rm -r $casePath/CFD/pascal/10
    rm -r $casePath/CFD/pascal/10.*
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

    if [[ $debugMode == "on" ]]; then
        debugMode="valgrind"
    elif [[ $debugMode == "strict" ]]; then
        #debugMode="valgrind --leak-check=full -v --trace-children=yes --track-origins=yes" 
        debugMode="valgrind --tool=memcheck --track-origins=yes --leak-check=yes --show-reachable=yes --num-callers=20 --track-fds=yes"  
    elif [[ $debugMode == "strictXML" ]]; then
        #debugMode="valgrind --leak-check=full -v --trace-children=yes --track-origins=yes" 
        debugMode="valgrind --tool=memcheck --leak-check=yes --leak-check-heuristics=stdstring,length64,newarray,multipleinheritance --show-reachable=no --num-callers=20 --track-fds=yes --xml=yes --xml-file=valgrind.xml"  
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
    $debugMode $CFDEM_LIGGGHTS_EXEC -in $solverName 2>&1 | tee -a $logpath/$logfileName

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

    if [[ $debugMode == "on" ]]; then
        debugMode="valgrind"
    elif [[ $debugMode == "strict" ]]; then
        #debugMode="valgrind --leak-check=full -v --trace-children=yes --track-origins=yes" 
        debugMode="valgrind --tool=memcheck --track-origins=yes --leak-check=yes --show-reachable=yes --num-callers=20 --track-fds=yes"  
    elif [[ $debugMode == "strictXML" ]]; then
        #debugMode="valgrind --leak-check=full -v --trace-children=yes --track-origins=yes" 
        debugMode="valgrind --tool=memcheck --leak-check=yes --leak-check-heuristics=stdstring,length64,newarray,multipleinheritance --show-reachable=no --num-callers=20 --track-fds=yes --xml=yes --xml-file=valgrind.xml"  
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
    if [[ $machineFileName == "none" ]]; then
        mpirun -np $nrProcs $debugMode $CFDEM_LIGGGHTS_EXEC -in $solverName 2>&1 | tee -a $logpath/$logfileName
    else
        mpirun -machinefile $machineFileName -np $nrProcs $debugMode $CFDEM_LIGGGHTS_EXEC -in $solverName 2>&1 | tee -a $logpath/$logfileName
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

    if [[ $debugMode == "on" ]]; then
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

    if [[ $debugMode == "on" ]]; then
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
    if [[ $machineFileName == "none" ]]; then
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
    decomposeCase=${10}
    remoteStorageLocation=${11}
    #--------------------------------------------------------------------------------#

    if [[ $debugMode == "on" ]]; then
        debugMode="valgrind"
    elif [[ $debugMode == "strict" ]]; then
        #debugMode="valgrind --leak-check=full -v --trace-children=yes --track-origins=yes" 
        debugMode="valgrind --tool=memcheck --track-origins=yes --leak-check=yes --show-reachable=yes --num-callers=20 --track-fds=yes"
    elif [[ $debugMode == "strictXML" ]]; then
        #debugMode="valgrind --leak-check=full -v --trace-children=yes --track-origins=yes" 
        debugMode="valgrind --tool=memcheck --leak-check=yes --leak-check-heuristics=stdstring,length64,newarray,multipleinheritance --show-reachable=no --num-callers=20 --track-fds=yes --xml=yes --xml-file=valgrind.xml"  
    elif [[ $debugMode == "profile" ]]; then
        if [[ $WM_COMPILE_OPTION == "Debug" ]]; then
            # make sure you did hpcstruct before
            debugMode="hpcrun"
            rm -r $casePath/CFD/hpctoolkit-$solverName-measurements
        else
            echo ""
            echo "WARNING - you do not use proper Debug flags! (for hpcrun you need Debug)"
            echo "using debugMode off now."
            debugMode=""
            read
        fi
    elif [[ $debugMode == "callgrind" ]]; then
        if [[ $WM_COMPILE_OPTION == "Debug" ]]; then
            debugMode="valgrind --tool=callgrind"
        else
            echo ""
            echo "WARNING - you do not use proper Debug flags! Only when using debug, cachegrind will be able to look through you code."
            debugMode="valgrind --tool=callgrind"
            read
        fi

    else
        debugMode=""
    fi

    #- clean up old log file
    rm $logpath/$logfileName

    #- change path
    cd $casePath/CFD

    #- decompose case
    if [[ $decomposeCase == "false" ]]; then   
        echo "Not decomposing case."
    else
        echo "Decomposing case."
        decomposePar -force

        if [[ $remoteStorageLocation == "" ]]; then
            echo "do nothing."
        else
            echo "do links"
            linkProcDirs $remoteStorageLocation
        fi
    fi

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
            reconstructPar -noLagrangian
        fi
    else
        mpirun -machinefile $machineFileName -np $nrProcs $debugMode $solverName -parallel 2>&1 | tee -a $logpath/$logfileName

        #- reconstruct case
        if [[ $reconstuctCase == "true" ]]; then   
            #pseudoParallelRun "reconstructPar" $nrProcs
            reconstructPar -noLagrangian
        fi
    fi

    if [[ $debugMode == "hpcrun" ]]; then
        if [ -f $CFDEM_APP_DIR/$solverName.hpcstruct ]; then
            rm -r hpctoolkit-$solverName-database*
            hpcprof -S $CFDEM_APP_DIR/$solverName.hpcstruct -I ./'*' hpctoolkit-$solverName-measurements
        else
            echo "you need to run hpcstruct first for your app!"
            read
        fi   
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
        echo "$SOLVERNAME" >> $logpath/log_compile_results_sol_success
    else
        echo "$SOLVERNAME" >> $logpath/log_compile_results_sol_fail
    fi
}

collectLogCFDEMcoupling_src()
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
    if [[ $LASTWORD == $SOLVERNAME || ${LASTWORD: -3} == ".so" ]]; then
        echo "$SOLVERNAME" >> $logpath/log_compile_results_src_success
    else
        echo "$SOLVERNAME" >> $logpath/log_compile_results_src_fail
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

if [[ -z "$njobs" ]]; then 
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

let nchunk=$nsteps/$njobs+1
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
  let nn=$i*$nchunk
  tstop=`ls processor0 | sed -n $nn$p`
  if [[ $i == $njobs ]] 
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
##    let n2=$n1-1
    nnow=`ls -d [0-9]*/ | wc -l` ##count time directories in case root dir, this will include 0
    let nmore=$nsteps-$nnow+1 ##calculate number left to reconstruct and subtract 0 dir
    if [[ $nmore != $nmore_old ]]
      then
      echo "$nmore directories remaining..."
    fi
    nmore_old=$nmore
  done

#combine and cleanup
if [[ -n "$outputfile" ]] 
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

        if [[ $CMD == $appname ]]
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
    overwrite="$2"
    #--------------------------------------------------------------------------------#

    if [[ $overwrite == "true" ]]
    then
        sed -i 's/[(,)]//g' $oldFileName
    else
        newFileName="$oldFileName""_noBrackets"
        sed -e 's/[(,)]//g' $oldFileName > $newFileName
    fi
}

#========================================#
#- remove brackets from file
linkProcDirs()
{
    #--------------------------------------------------------------------------------#
    #- define variables
    remoteDir="$1"
    #--------------------------------------------------------------------------------#

    # check remote directory exist
    if [[ -d "$remoteDir" ]]
    then
        echo "  $remoteDir exists - **check**"
    else
        Fatal "ERROR: remote directory does not exist: $remoteDir"
    fi

    # check case is decomposed
    if [[ -d "processor0" ]]
    then
        echo "  case is decomposed - **check**"
    else
        Fatal "ERROR: case is NOT decomposed - please do so and run again."
    fi

    # check CFD dir exists at remote location
    if [[ -d $remoteDir/CFD ]]
    then
        Fatal "ERROR: there is (!) a dir called CFD at $remoteDir - please rename it and run again."
    else
        echo "  there is no dir called CFD at $remoteDir - **check**"
    fi

    # create a dir CFD at remoteDir
    mkdir $remoteDir/CFD

    # create a backup dir oldProcDirs
    mkdir oldProcDirs

    # check if oldProcDirs is empty
    if [[ -d oldProcDirs/processor0 ]]
    then
        Fatal "ERROR: ./oldProcDirs is not empty - please clean up and run again."
    else
        echo "  ./oldProcDirs is empty - **check**"
    fi

    # moving proc dirs to oldProcsDir (backup)
    cp -r processor* $remoteDir/CFD
    mv processor* oldProcDirs

    # create a link to remote proc dirs
    counter=0
    for procDir in $remoteDir/CFD/*/; do
        ln -s $remoteDir/CFD/processor$counter processor$counter
        counter=$[counter + 1]
    done

    # create a link to postProcessing dir
    mkdir postProcessing
    mv postProcessing $remoteDir/CFD
    ln -s $remoteDir/CFD/postProcessing postProcessing

    # success
    echo "linking was successful"
}

#========================================#
#- function to enable a LIGGGHTS package
enableLiggghtsPackage()
{
    #--------------------------------------------------------------------------------#
    #- define variables
    pkgName="$1"
    whitelist="$CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/package-liggghts-list.txt"
    #--------------------------------------------------------------------------------#

    if [ ! -f "$CWD/$whitelist" ];then
        echo "$whitelist does not exist in $CWD. Nothing will be done."
        NLINES=0
        COUNT=0
    else
        NLINES=`wc -l < $CWD/$whitelist`
        COUNT=0
    fi

		found=false
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
								echo "$LINE"
								if [[ "$LINE" == "$pkgName" ]]; then
										echo "Package $pkgName already enabled"
										found=true
										break
								fi
            fi
    done

		if [ $found = false ]; then
				echo "Package $pkgName not found - add it to list"
				echo $pkgName"/dir" >> $CWD/$whitelist
		fi
}

