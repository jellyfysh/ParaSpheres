#!/bin/bash 
set -e # exit script on error at any stage, in particular script exits if the data directory already exists

SystemDimension=127
Density=0.708
LogKmax=10
ChainLength=4.0


printf "\n\e[94m 		Data files will be generated in a directory SetupData/Data.xxx \e[0m\n\n"
cd ../SetupData
if [ ! -d  "Data.$SystemDimension" ]; then
	mkdir Data.$SystemDimension
fi
cd Data.$SystemDimension

if [ ! -f "info.dat" ]; then
	printf "\n\e[94m 		Substitute SystemDimension and Density in Fortran code\e[0m\n\n"
	(cd ../../Fortran90/CellECMC/ ; cat phys_parameters.sed | sed "s/SIZE/$SystemDimension/g" |sed "s/RHO/$Density/g" > phys_parameters.f90)
	printf "\n\e[91m 		Compile fortran code\e[0m\n\n"
	( cd ../../Fortran90/CellECMC/;make clean; make )
	export PATH=$PATH:../../Fortran90/CellECMC

	printf "\n\e[91m 		Create starting configuration by running Fortran code\e[0m\n\n"
	CellECMC create_init


	printf "\n\e[91m 		Check configuration by rerunning the fortran\e[0m\n\n"
	CellECMC benchmark

	printf "\n\e[91m 		Make constraint graph\e[0m\n\n"
	mv init.dat configuration.dat
	echo ChainLength $ChainLength >> info.dat
	echo LogKmax $LogKmax >> info.dat
else
	echo Existing initial configuration, skip configuration creation
fi

if [ -f "constraints.dat" ]; then
	echo Existing constraint graph, exit
	exit 0
fi
pypy ../../Python/CreateG/GenerateG3.py

printf "\n\e[91m 		Prune graph\e[0m\n\n"
pypy ../../Python/CreateG/PruneG.py
ln -s G_pruned4s.dat constraints.dat


printf "\n\e[91m 		Generate reference Liftings\e[0m\n\n"

for ((LogNo=0;LogNo<=LogKmax;LogNo++))
do
    pypy ../../Python/GlobalTimeECMC/GlobalTimeECMC.py $(( 2 ** LogNo))
done

exit 0
