#!/bin/bash
TrialDimension=5

set -e # exit on any errors
cd ../Run

mkdir SequentialC5.$TrialDimension.$$
rm -f Data.run
ln -s ../SetupData/Data.$TrialDimension Data.run
cd SequentialC5.$TrialDimension.$$


export PATH=$PATH:../../Python/SingleThreadLocalTimeECMC

printf "\n\e[91m                Markov chain analysis, system size=5\e[0m\n\n"
python3 ../../Python/SequentialMultiThreadECMC/SequentialMultiThreadECMC.py

exit 0
