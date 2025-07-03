#!/bin/bash

#===================================================================#
# script which should support changing the syntax in the liggghts
# input files to new integrators and force coupling.
# WARNING: Use this script at YOUR OWN RISK!!!
# Please use with caution and prepare a copy of your files
# BEFORE (!!!) you run this script.
# DCS Computing GmbH - June 2018
#===================================================================#
# SYNTAX CHANGES FOR SPHERES:
#A) nve/sphere-->nve/cfd_cn/sphere
#B) couple/cfd/force/implicit-->couple/cfd/force
#C) CrankNicolson entry goes from couple/cfd/force/implicit to nve/cfd_cn/sphere
#D) CAddRhoFluid is only used with couple/cfd/force and entries are split up to ${Cadd} and ${rhoFluid}
#E) transfer_type yes is now handled via transfer_property type
#F) fix multisphere is no longer an integrator and needs a separate integrator fix integrator all nve/cfd_cn/nonspherical
#   the fix multisphere and fix integrator command must be before the fix couple/cfd/force/multisphere
#G) IB cases need for the fix couple/cfd/force the extra keyword torque explicit to communicate hdtorque
#===================================================================#

echo "WARNING: Use this script at YOUR OWN RISK (!!!). Please use with caution and prepare a copy of your files. BEFORE (!!!) you run this script."
echo "In doubt press Ctrl-C!"
read


echo "apply rule A) nve/sphere-->nve/cfd_cn/sphere"
echo "to all in.* (files named differently will be ignored!)"
read
find . -iname "in.*" -exec sed -i 's/nve\/sphere/nve\/cfd_cn\/sphere/g' {} +

echo "apply rule B) couple/cfd/force/implicit-->couple/cfd/force"
echo "to all in.* (files named differently will be ignored!)"
read
find . -iname "in.*" -exec sed -i 's/couple\/cfd\/force\/implicit/couple\/cfd\/force/g' {} +

echo "apply rule B)2) couple/cfd/force/integrateImp-->couple/cfd/force && nve/cfd_cn/sphere"
echo "to all in.* (files named differently will be ignored!)"
echo "please add the command fix integr all nve/cfd_cn/sphere to the files being opend in gedit"
echo "Then please save the files and close the gedit window. Then press enter."
read
find . -iname "in.*" -exec grep -rl "integrateImp" {} + | xargs gedit -s;
read
find . -iname "in.*" -exec sed -i 's/couple\/cfd\/force\/integrateImp/couple\/cfd\/force/g' {} +

echo "Please manually apply rule C): CrankNicolson entry goes from couple/cfd/force/implicit to nve/cfd_cn/sphere (the files will open in gedit.)"
echo "Then please save the files and close the gedit window. Then press enter."
read
find . -iname "in.*" -exec grep -rl "CrankNicolson" {} + | xargs gedit -s;
read

echo "Please manually apply rule D): CAddRhoFluid is only used with couple/cfd/force and entries are split up to \${Cadd} and \${rhoFluid}"
echo "Then please save the files and close the gedit window. Then press enter."
read
find . -iname "in.*" -exec grep -rl "CAddRhoFluid" {} + | xargs gedit -s;

echo "Please manually apply rule E): transfer_type yes is now handled via couple/cfd (not force!!!) transfer_property name type vector_atom"
echo "Then please save the files and close the gedit window. Then press enter."
read
find . -iname "in.*" -exec grep -rl "transfer_type" {} + | xargs gedit -s;

echo "Please manually apply rule F): fix multisphere is no longer an integrator and needs a separate integrator fix integrator all nve/cfd_cn/nonspherical."
echo "The fix multisphere and fix integrator command must be before the fix couple/cfd/force/multisphere!"
echo "Then please save the files and close the gedit window. Then press enter."
read
find . -iname "in.*" -exec grep -rl "multisphere" {} + | xargs gedit -s;

