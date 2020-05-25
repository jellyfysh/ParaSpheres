#!/bin/bash
set -e # exit on errors
rm -f toto2.dat
logkmax=$1
Nthreadmax=$2
Nrun=$3
export PATH=$PATH:.
for k in $(seq 1 $logkmax)
do
    printf "log_2 k: %d\n"  $k
    for j in $(seq 2 $Nthreadmax)
    do
	if [[ $(( j - 1 )) -eq $(( 2 ** k)) ]]; then
	    break
	fi
	printf "\tthread %d\n"  $j
	for i in $(seq 1 $Nrun)
	do
	    OMP_NUM_THREADS=$j MultiThreadECMC validation $(( 2 ** k )) >>  toto2.dat
	    #OMP_NUM_THREADS=2 simul validation >>  toto2.dat
	done
    done
done
echo "runs: "
cat toto2.dat |grep OK |wc -l
printf "\e[33mBugs  \n"
cat toto2.dat |grep 'Not OK' |wc -l
printf "\e[0m "
