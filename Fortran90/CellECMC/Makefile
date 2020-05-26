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
    FC=gfortran
endif
ifeq ($(UNAME), Darwin)
    FC=gfortran-mp-9
endif

FFLAGS=-g -O 

all: CellECMC

%.o: %.f90
	$(FC) $(FFLAGS) -c $(switch) $<

%.mod: %.f90
	$(FC) $(FFLAGS) -c $(switch) $<

FILES=phys_parameters.o cell_parameters.o others.o  benchmark.o alg_parameters.o  CellECMC.o expl_cell_x.o expl_cell_y.o init_cell.o init_pos.o   pick_move.o refresh_cell_x.o refresh_cell_y.o write_pos.o  resetseed.o

CellECMC:$(FILES)
	$(FC) $(FFLAGS) $(FILES) -o $@

clean:
	$(RM) -rf *.o *.mod CellECMC latex html run.dSYM *.dat
