#!/bin/bash
source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/bashrc
#shopt -s expand_aliases

ETCpath=$CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/OFVersionChange/shellScripts

$ETCpath/unComment.sh $1 CFD/constant/turbulenceProperties
$ETCpath/unComment.sh $1 CFD/constant/transportProperties
$ETCpath/unComment.sh $1 CFD/constant/couplingProperties
$ETCpath/unComment.sh $1 CFD/constant/couplingProperties_run
$ETCpath/unComment.sh $1 CFD/constant/couplingProperties_restart
$ETCpath/unComment.sh $1 CFD/system/fvOptions
