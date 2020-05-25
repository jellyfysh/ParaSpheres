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
