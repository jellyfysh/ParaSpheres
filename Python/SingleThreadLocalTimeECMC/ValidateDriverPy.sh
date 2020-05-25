#!/bin/bash
set -e # exit on errors
rm -f toto2.dat
logkmax=$1
imax=$2
export PATH=$PATH:.
for k in $(seq 1 $logkmax)
do
	printf "log_2 k: %d\n"  $k
	for i in $(seq 1 $imax)
	do
		pypy ../../Python/SingleThreadLocalTimeECMC/SingleThreadLocalTimeECMC.py $(( 2 ** k)) >>  toto2.dat
	done
done
echo "runs: "
cat toto2.dat |grep OK |wc -l
echo "Bugs: "
cat toto2.dat |grep 'Not OK' |wc -l
