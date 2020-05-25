#!/bin/bash
SystemDimension=127
set -e

cd ../Run

mkdir ValidateG.$SystemDimension.$$
rm -f Data.run
ln -s ../SetupData/Data.$SystemDimension Data.run
cd ValidateG.$SystemDimension.$$

set -e # exit on any errors

printf "\n\e[91m                Build graph checker based on Cell simulation\e[0m\n\n"
(cd ../../CPP/GraphValidateCellECMC; make clean; make -j 8)

printf "\n\e[91m                Run graph checker\e[0m\n\n"

../../CPP/GraphValidateCellECMC/GraphValidateCellECMC

exit 0
