#!/bin/bash
SystemDimension=127
LogKmax=10
Nthreadmax=10
Nrun=10

set -e
cd ../Run

mkdir CValidateECMC.$SystemDimension.$$
rm -f Data.run
ln -s ../SetupData/Data.$SystemDimension Data.run
cd CValidateECMC.$SystemDimension.$$


printf "\n\e[91m                Build parallel validation program\e[0m\n\n"
(cd ../../CPP/MultiThreadECMC/; make clean ; make -j 8)
export PATH=$PATH:../../CPP/MultiThreadECMC/

printf "\n\e[91m                C++ Validation script\e[0m\n\n"
ValidateDriverC.sh $LogKmax $Nthreadmax $Nrun

exit 0
