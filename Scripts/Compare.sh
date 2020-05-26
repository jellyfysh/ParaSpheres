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
