#!/bin/bash
set -e
k=$1
Nthreadmax=$2
Nrun=$3
export PATH=$PATH:.
rm -f EPH.dat
for i in $(seq 1 $Nthreadmax)
do
    printf "\e[33momp_num_threads \e[0m %d\n"  $i
    for j in $(seq 1 $Nrun)
    do
	OMP_NUM_THREADS=$i OMP_PLACEs="{0:20:1}{40:20:1}" OMP_PROC_BIND=false MultiThreadECMC timing $k | grep EPH >>  EPH.dat
	sleep .5
    done
done
