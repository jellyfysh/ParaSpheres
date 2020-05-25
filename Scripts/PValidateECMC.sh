#!/bin/bash
SystemDimension=127
LogKmax=10
Nrun=10 #increase to 1000 for a longer run

set -e # exit on any errors
cd ../Run

mkdir PValidateECMC.$SystemDimension.$$
rm -f Data.run
ln -s ../SetupData/Data.$SystemDimension Data.run
cd PValidateECMC.$SystemDimension.$$

export PATH=$PATH:../../Python/SingleThreadLocalTimeECMC

printf "\n\e[91m                Python validation script\e[0m\n\n"
ValidateDriverPy.sh $LogKmax $Nrun

exit 0
