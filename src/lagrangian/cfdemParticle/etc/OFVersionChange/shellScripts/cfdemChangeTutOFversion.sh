#!/bin/bash

source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/bashrc
#shopt -s expand_aliases

# find my OF version
if [ "$WM_PROJECT_VERSION" == "3.0.x" ]; then
    echo 'You are using OpenFOAM 3.0.x. Will activate dicts now...'
    OFV=OFversion30x
else
    echo 'You are using OpenFOAM 2.4.x or lower. Will activate dicts now...'
    OFV=OFversion24x
fi
ETCpath=$CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/OFVersionChange/shellScripts


# remove old comment // in every line ending with OFversionXYZ
# doing this three times to make sure multiple comments are removed
bash $ETCpath/activateVersion.sh OFversion
bash $ETCpath/activateVersion.sh OFversion
bash $ETCpath/activateVersion.sh OFversion

# adding comment // to all lines ending with OFversionXYZ
bash $ETCpath/commentOut.sh OFversion CFD/constant/turbulenceProperties
bash $ETCpath/commentOut.sh OFversion CFD/constant/transportProperties
bash $ETCpath/commentOut.sh OFversion CFD/constant/couplingProperties
bash $ETCpath/commentOut.sh OFversion CFD/constant/couplingProperties_run
bash $ETCpath/commentOut.sh OFversion CFD/constant/couplingProperties_restart
bash $ETCpath/commentOut.sh OFversion CFD/system/fvOptions

# removing the comment // in every line ending with $OFV
bash $ETCpath/activateVersion.sh $OFV
