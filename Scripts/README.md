# Scripts

This directory provides a number of shell scripts, which assemble all the required elements to run the `Fortran90`, `C++`, and `Python3` programs described of [Li2020](https://arxiv.org/abs/2004.11040) (see [References.bib](References.bib)).

The shell scripts must be run from within the present directory. Configurations and data files are then found by the use of relative path names.

Execution assumes a Unix environment, and has been tested under both Linux, and MacOS. Use of OpenMP requires the use of a non-system compiler on MacOS, either gcc, or clang. The variables compiler variables FC, and CXX may need to be set within the Makefiles: `CPP/GraphValidateCellECMC/Makefile`, `CPP/MultiThreadECMC/Makefile` , `Fortran90/CellECMC/Makefile`.

Generated data files are placed in `SetupData/` by `Setup.sh`, all other scripts create directories under `Run/`.

## Run-time Parameters
The scripts are given in a form such that they can be run on a desktop machine, in about ten minutes of clock time. We thus use system sizes and threading parameters which are *smaller* than those used in [Li2020](https://arxiv.org/abs/2004.11040). In particular, the scripts use `SystemDimension=127`, `LogKmax=10`, `Nrun=10`, `Nthreadmax=10`.

To reproduce the results of [Li2020](https://arxiv.org/abs/2004.11040), these variables should be increased to: `SystemDimension=256`, `LogKmax=13`, `Nrun>1000`, `Nthreadmax=40`. The benchmark performance will depend strongly on the number of hardware threads on the computer.


## `Setup.sh`
Calls `Fortran90/CellECMC/CellECMC`, `Python/CreateG/GenerateG3.py`, `Python/CreateG/PruneG.py`, and  `Python/GlobalTimeECMC/GlobalTimeECMC.py`.

Is used to:
* create a configuration
* calculate the pruned constraint graph
* generate reference liftings

### Requirements
`Make`, `Python3`, as well as a `Fortran90` compiler.

### Usage
* Run `./Setup.sh`
* Parameters in `Setup.sh`:
    * `SystemDimension`: Square root of number of spheres
    * `Density`: Packing density of hard sphere (Only used when creating new configurations)
    * `LogKmax`: log_2(maximal number of active spheres for creating reference liftings)
    * `ChainLength`: Total length of chain (all active spheres combined) when doing validation

### Input
* No input is needed for creating a new configuration.
* `configuration.dat`, `info.dat`, and possibly `active*.dat` in `SetupData/Data.$SystemDimension` if creating constraint graph and reference liftings for existing configurations.

### Output
* `configuration.dat`, `info.dat`, `constraints.dat`, `active*.dat` ,and `liftings*.dat` in `SetupData/Data.$SystemDimension` when creating a new configuration.
* `constraints.dat`, `active*.dat` ,and `liftings.dat` in `SetupData/Data.$SystemDimension` when creating a new configuration. (`constraints.dat` and `liftings.dat` when `active*.dat` are provided.)

## `ValidateG.sh`
Performs a simulation of the hard sphere system using a cell based algorithm. Each lifting in the simulation is then compared to a constraint graph, to check the validity of the graph. The program counts the number of constraints visited as a function of the number of liftings.

Calls `CPP/GraphValidateCellECMC/GraphValidateCellECMC`

### Requirements
`Make`, `C++` (together with the boost library for array creation, https://www.boost.org/)

### Usage
* Run `./ValidateG.sh`
* Parameters in `ValidateG.sh`:
    * `SystemDimension`: square root of number of spheres

### Input
`configuration.dat`, `info.dat`, and `constraints.dat` in `SetupData/Data.$SystemDimension`.

### Output
Outputs to `CPP/GraphValidateCellECMC/GraphValidateCellECMC` in `Run/ValidateG.$SystemDimension.$$` (where `$$` is the process id).

## `PValidateECMC.sh`
Calls `Python/SingleThreadLocalTimeECMC/validation_py.sh`, which calls `Python/SingleThreadLocalTimeECMC/SingleThreadLocalTimeECMC.py`.

Validates the Python implementation of Algorithm 4 in [Li2020](https://arxiv.org/abs/2004.11040).

### Requirements
`Python90` or `PyPy`

### Usage
* Run `./PValidateECMC.sh`
* Parameters in `PValidateECMC.sh`:
    * `SystemDimension`: Square root of number of spheres
    * `LogKmax`: log_2(maximal number of active spheres during the run)
    * `Nrun`: Number of runs for each number of active spheres

### Input
All the files in `SetupData/Data.$SystemDimension`

### Output
Outputs of `Python/SingleThreadLocalTimeECMC/SingleThreadLocalTimeECMC.py` in `Run/PValidateECMC.$SystemDimension.$$` (where `$$` is the process id). Prints the number of total runs and the number of problematic runs.

## `CValidateECMC.sh`
Calls `CPP/MultiThreadECMC/MultiThreadECMC`.

Validates the Python3 implementation of Algorithm 6 in [Li2020](https://arxiv.org/abs/2004.11040).

### Requirements
`Make`, `C++`

### Usage
* Run `./CValidateECMC.sh`
* Parameters in `CValidateECMC.sh`:
    * `SystemDimension`: Square root of number of spheres
    * `LogKmax`: log_2(maximal number of active spheres during the run)
    * `Nthreadmax`: Maximal number of threads
    * `Nrun`: Number of runs for each number of active spheres

### Input
All the files in `SetupData/Data.$SystemDimension`

### Output
Outputs of `CPP/MultiThreadECMC/MultiThreadECMC` in `Run/CValidateECMC.$SystemDimension.$$` (where `$$` is the process id). Prints the number of total runs and the number of problematic runs.

## `Ordering.sh`
Calls `CPP/MultiThreadECMC/MultiThreadECMC`.

Validates the Python3 implementation of Algorithm 6 in [Li2020](https://arxiv.org/abs/2004.11040) in small systems.

### Requirements
`Make`, `C++`

### Usage
* Run by `./Ordering.sh`
* Parameters in `Ordering.sh`:
    * `TrialDimension`: Number of spheres (4 or 5 for provided small configurations)
    * `LogK`: log_2(maximal number of active spheres during the run)
    * `Nthread`: Maximal number of threads
    * `Ntest`: Number of runs for each number of active spheres

### Input
All the files in `SetupData/Data.$SystemDimension`

### Output
Outputs of `CPP/MultiThreadECMC/MultiThreadECMC` in `Run/Ordering.$SystemDimension.$$` (where `$$` is the process id). Prints the number of total runs and the number of problematic runs.

## `SeqentialC4.sh`
Calls `Python/SequentialMultiThreadECMC/SequentialMultiThreadECMC.py`.

Validates the Python implementation of Algorithm 5 of [Li2020](https://arxiv.org/abs/2004.11040) in small systems.

### Requirements
`Python3`

### Usage
* Run by `./SeqentialC4.sh`
* `SeqentialC4.sh` is designed for the 4-sphere configuration.


### Input
`configuration.dat` and `info.dat` in `SetupData/Data.4`

### Output
Output of `Python/SequentialMultiThreadECMC/SequentialMultiThreadECMC.py` in `Run/SequentialC4.$SystemDimension.$$` (where `$$` is the process id).

## `SeqentialC5.sh`
Calls `Python/SequentialMultiThreadECMC/SequentialMultiThreadECMC.py`.

Validates the Python implementation of Algorithm 5 of [Li2020](https://arxiv.org/abs/2004.11040) in small systems.

### Requirements
`Python3`

### Usage
* Run by `./SeqentialC5.sh`
* `SeqentialC5.sh` is designed for the 5-sphere configuration.


### Input
`configuration.dat` and `info.dat` in `SetupData/Data.5`

### Output
Output of `Python/SequentialMultiThreadECMC/SequentialMultiThreadECMC.py` in `Run/SequentialC5.$SystemDimension.$$` (where `$$` is the process id. It serves as random-number generator here.).

## `Benchmark.sh`
Calls `CPP/MultiThreadECMC/BenchDriver.sh`, which calls `CPP/MultiThreadECMC/MultiThreadECMC`.

Measures the performance of this program.

### Requirements
`Make`, `C++`

### Usage
* Run by `./Benchmark.sh $name_of_the_run`, where `$name_of_the_run` is the name of the run
* Parameters in `Benchmark.sh`:
    * `SystemDimension`: Square root of number of spheres
    * `k`: Number of active spheres
    * `Nthreadmax`: Maximal number of threads
    * `Nrun`: Number of runs for each number of active spheres
* Environmental variables related to OpenMP are in `CPP/MultiThreadECMC/time.sh`

### Input
`configuration.dat`, `constraints.dat` and `info.dat` in `SetupData/Data.$SystemDimension`

### Output
* Numbers of events per hour in file `EPH.dat` in `Run/Bench.$SystemDimension.$name_of_the_run`
* A box plot (similar to Figure 7 [Li2020](https://arxiv.org/abs/2004.11040)) is created as `EPH.pdf` in `Run/Bench.$SystemDimension.$name_of_the_run`.

## `RunAllScripts.sh`
Loops over all other scripts, to perform a complete test of all code.
