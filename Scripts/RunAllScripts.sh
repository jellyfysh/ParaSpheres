#!/bin/bash 
set -e #exit on error

#create a working configuration
bash  Setup.sh

#tests based on the generated configuration
for i in ValidateG.sh CValidateECMC.sh Ordering.sh BenchmarkECMC.sh PValidateECMC.sh  Compare.sh
do
	bash $i
done

# tests with the small 4/5 particle configurations
for i in   SequentialC4.sh SequentialC5.sh 
do
	bash $i
done
