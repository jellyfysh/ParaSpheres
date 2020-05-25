#!/bin/bash
TrialDimension=4

set -e # exit on any errors

cd ../Run

mkdir SequentialC4.$TrialDimension.$$
rm -f Data.run
ln -s ../SetupData/Data.$TrialDimension Data.run
cd SequentialC4.$TrialDimension.$$


export PATH=$PATH:../../Python/SingleThreadLocalTimeECMC

printf "\n\e[91m                Markov chain analysis system size = 4\e[0m\n\n"
python3 ../../Python/SequentialMultiThreadECMC/SequentialMultiThreadECMC.py

exit 0
