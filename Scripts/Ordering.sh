#!/bin/bash
#These parameters they run a specific configuration in SetupData/Data.4/ using the general purpose C++ code.
TrialDimension=4
LogK=1
Nthread=2
Ntest=1000

set -e # exit on any errors
cd ../Run

mkdir Ordering.$TrialDimension.$$
rm -f Data.run
ln -s ../SetupData/Data.$TrialDimension Data.run
cd Ordering.$TrialDimension.$$


printf "\n\e[91m                Build parallel validation program\e[0m\n\n"
(cd ../../CPP/MultiThreadECMC/; make clean ; make -j 8)
export PATH=$PATH:../../CPP/MultiThreadECMC/

printf "\n\e[91m                C++ Validation script small, hand created configuration \e[0m\n\n"
ValidateDriverC.sh  $LogK $Nthread $Ntest

exit 0
