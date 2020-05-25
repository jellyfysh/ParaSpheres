#!/bin/bash
SystemDimension=127
k=40
Nthreadmax=10
Nrun=10

set -e # exit on any errors
cd ../Run
rm -f Data.run
ln -s ../SetupData/Data.$SystemDimension Data.run

mkdir Bench.$SystemDimension.$$
cd Bench.$SystemDimension.$$

#cp ../CPP/ParallelBenchmark/para.cc .

printf "\n\e[91m                Build parallel C++ program\e[0m\n\n"
(cd ../../CPP/MultiThreadECMC/; make clean ; make -j 8)
export PATH=$PATH:../../CPP/MultiThreadECMC/

printf "\n\e[91m                Benchmark script\e[0m\n\n"
bash  BenchDriver.sh $k $Nthreadmax $Nrun
python3 ../../CPP/MultiThreadECMC/speed_up_simplify.py $Nthreadmax $Nrun
