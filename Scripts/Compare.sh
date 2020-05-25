#!/bin/bash
SystemDimension=127

set -e # exit on any errors
cd ../Run
rm -f Data.run
ln -s ../SetupData/Data.$SystemDimension Data.run

mkdir Compare.$SystemDimension.$$
cd Compare.$SystemDimension.$$

#cp ../CPP/ParallelBenchmark/para.cc .

printf "\n\e[91m                Build program\e[0m\n\n"
(cd ../../CPP/MultiThreadECMC/; make clean ; make -j 8)
export PATH=$PATH:../../CPP/MultiThreadECMC/

printf "\n\e[91m                Compare script between threaded/unthreaded \e[0m\n\n"
bash  CompareDriver.sh
