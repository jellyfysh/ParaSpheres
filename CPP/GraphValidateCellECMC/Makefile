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
