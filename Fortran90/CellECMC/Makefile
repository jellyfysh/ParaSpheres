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
