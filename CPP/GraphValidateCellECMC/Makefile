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

ifeq ($(UNAME), Linux)
    CXX=g++
endif
ifeq ($(UNAME), Darwin)
    CXX=g++-mp-9
endif

CXXFLAGS=-Wall  -g -O  -Wextra  -std=c++17 -I/opt/local/include -Wshadow

OBJ= cell.o init.o refresh.o main.o inout.o GraphValidateCellECMC.o checker.o 

all:GraphValidateCellECMC

GraphValidateCellECMC: $(OBJ)
	$(CXX)  $(CXXFLAGS) $(OBJ)  -o $@  -lm

clean:
	$(RM) *.o benchmark Checker.tar.gz simul bench *~ checker GraphValidateCellECMC

clobber:clean
	$(RM)  unfoundConstraints.dat foundConstraints.dat active_history.dat nzero.dat count.dat checker

# DO NOT DELETE

cell.o: Checker.h Param.h
checker.o: Checker.h Param.h
init.o: Checker.h Param.h
inout.o: Checker.h Param.h
main.o: Checker.h Param.h
refresh.o: Checker.h Param.h
GraphValidateCellECMC.o: Checker.h Param.h
Checker.o: Param.h
