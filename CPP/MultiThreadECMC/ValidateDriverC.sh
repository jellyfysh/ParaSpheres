#!/bin/bash
# JeLLyFysh/ParaSpheres - Multithreaded event-chain Monte Carlo with local times
# - https://github.com/jellyfysh/paraspheres
# Copyright (C) 2020 The JeLLyFysh organization
# (see the AUTHORS file on jellyfysh/paraspheres for the full list of authors)
#
# This file is part of JeLLyFysh/ParaSpheres.
#
# JeLLyFysh/ParaSpheres is free software: you can redistribute it and/or modify
# it under the terms of the GNU General
# Public License as published by the Free Software Foundation, either > version
# 3 of the License, or (at your option)
# any later version.
#
# JeLLyFysh/ParaSpheres is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even 
# the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# JeLLyFysh/ParaSpheres in the LICENSE file. If not, see
# <https://www.gnu.org/licenses/>.
#
# If you use JeLLyFysh/ParaSpheres in published work, please cite the following
# reference
# (see [Li2020] in References.bib):
# Botao Li, Synge Todo, A. C. Maggs, Werner Krauth
# Multithreaded event-chain Monte Carlo with local times,
# arXiv e-prints: 2004.11040 (2020), https://arxiv.org/abs/2004.11040

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
