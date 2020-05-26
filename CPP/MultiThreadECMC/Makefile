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
#

UNAME := $(shell uname)

DEBUG= -fsanitize=address -fno-omit-frame-pointer  -Wpadded 
VECTOR= -ftree-vectorize -fopt-info-loop-optimized 

#CXXFLAGS=-Wall -O2 -g -Wextra -Wshadow -std=c++2a -fopenmp  #  $(DEBUG)

ifeq ($(UNAME), Linux)
    CXX=g++
    CXXFLAGS=-Wall -Wextra -O2 -g -Wall -fopenmp -std=c++17 -Wshadow	
endif
ifeq ($(UNAME), Darwin)
    CXX=g++-mp-10	
    CXXFLAGS=-Wall -Wextra -O2 -g -Wall -fopenmp -std=c++20 -Wshadow # $(DEBUG)
endif


all: MultiThreadECMC

MultiThreadECMC: simul.o main.o inout.o MultiThreadECMC.o singlethread.o
	$(CXX) $(CXXFLAGS) simul.o main.o inout.o  MultiThreadECMC.o singlethread.o -o $@ -lpthread
clean:
	$(RM) -rf *.o *.html html simul j.ps *~ MultiThreadECMC

run:constraint_graph
	 ./constraint_graph

enscript:
	enscript -1rG --line-numbers --highlight=cpp  --color -pj.ps *cc *h  

lint:
	/opt/local/libexec/llvm-9.0/bin/scan-build  -o html make MultiThreadECMC

speed:
	python3 speed_up.py
clobber:clean
	$(RM) *.dat
lint2:
	cppcheck --enable=all *.cc

# DO NOT DELETE

inout.o: Simul.h Param.h
main.o: Simul.h Param.h
main1.o: Simul.h Param.h
MultiThreadECMC.o: Simul.h Param.h
simul.o: Simul.h Param.h
singlethread.o: Simul.h Param.h
