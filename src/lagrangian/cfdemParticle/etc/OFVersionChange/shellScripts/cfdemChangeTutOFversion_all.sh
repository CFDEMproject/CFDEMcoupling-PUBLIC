#! /bin/bash

##########################################################################
# Shellscript:	update tutorial cases in the whole tree according to OFV
# Author     :	Andreas Aigner <andreas.aigner@dcs-computing.com>
# Date       :	2016-04-15
# Category   :	Admin Utilities
##########################################################################
# Description
#    update tutorial cases in the whole tree according to OFV
##########################################################################

PN=`basename "$0"` # Program name
BASEDIR=`dirname "$0"`
VER='0.1'

Usage () {
    echo >&2 "$PN - update all tutorials in a whole tree, $VER
usage: $PN [-f] path
    -f: force to ignore security check;
        allows other directories than 'cfdemTut' and sub-directories
    -v: verbose output
    -h: show this help

    path: path to the root directory of all tutorials

Example:
    $PN $HOME/Repositories"
    exit 1
}

Msg () { echo >&2 "$*"; } # removed "$PN:" to shorten output
Fatal () { Msg "Error: $@"; exit 1; }

# defaults
force_flag=false
verbose_flag=false

while getopts hfv opt
do
    case "$opt" in
      f)  force_flag=true;;
      v)  verbose_flag=true;;
      h)  Usage;;
      \?) Usage;;
    esac
done
shift `expr $OPTIND - 1`

[ $# -lt 1 ] && Msg "No path defined" && Usage # check for path

# check if readlink exists (we need it for the savety check)
command -v readlink >/dev/null 2>&1 || { Fatal "readlink is required but it's not installed. Aborting."; }

# global paths
curRunDir=$(pwd)
cfdemTutDir=$(readlink -m $CFDEM_TUT_DIR)

#set -x
for dir in "$@"
do
    Msg "INFO: Checking $dir"

    # is directory?
    if [ ! -d "$dir" ]; then
      Msg "$dir is not a directory"
      continue
    fi

    # clean paths from multiple separators (pure bash)
    #dir=$(echo "$dir" | sed s#//*#/#g)
    #cfdemTutDir=$(echo "$CFDEM_TUT_DIR" | sed s#//*#/#g)

    # resolve absolute and relative paths (with readlink)
    # clean paths from multiple separators and resolve ./..
    if [ "${dir:0:1}" = "/" ]; then # absolute path
      curDir=$(readlink -m $dir)
    else # relative path
      curDir=$(readlink -m $curRunDir/$dir)
    fi

    # savety check: check if dir is within $CFDEM_TUT_DIR; skip if forced
    if [[ $curDir != $cfdemTutDir ]] && [ $force_flag != true ]; then
        Fatal "$curDir is not within $cfdemTutDir. Change the directory or use the option -f."
    fi

    # performe change steps
    ldir=$(find "$dir" -name 'CFD' -type d)
    for cfd in $ldir
    do
      casedir=$(dirname "$cfd") ## get parent directory of the CFD directory
      Msg "INFO: Updating $casedir"
      #echo "do sth here..."
      if [ $verbose_flag = true ]; then
        (cd $casedir && bash $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/OFVersionChange/shellScripts/cfdemChangeTutOFversion.sh)
      else
        (cd $casedir && bash $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/OFVersionChange/shellScripts/cfdemChangeTutOFversion.sh) &>/dev/null
      fi
    done
done
