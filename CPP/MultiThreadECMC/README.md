## Multithreaded ECMC (`MultiThreadECMC.cc`)
This `C++` program implements  Algorithm 6 of [Li2020](https://arxiv.org/abs/2004.11040) (see [References.bib](References.bib)).

The main work is performed in `MultiThreadECMC.cc`.

Uses OpenMP to control the multithreading.

### Requirements
`Make`, `C++`. Note on MacOS, the system compiler does not manage `OpenMP`. In this case, the user has to install `gcc`, or `clang`.

### Usage
* Called by `ValidateDriverC.sh`, `BenchDriver.sh`, and `CompareDriver.sh` in this directory
* `ValidateDriverC.sh` is itself called by both `Scripts/Ordering.sh` and `Scripts/CValidateECMC.sh`. 
* `BenchDriver.sh` called by `Scripts/Benchmark.sh`. 
* `CompareDriver.sh` called by `Scripts/Benchmark.sh`.


### Input
* All the files in the test suit (test suit specified by scripts).
* `ValidateDriverC.sh` needs 3 arguments. They are "log_2(maximal number of active spheres during the run)"", "maximal number of threads", and "number of runs for each case".
* `BenchDriver.sh` needs 3 arguments. They are "number of active spheres", "maximal number of threads", and "number of runs" for each case.
* `CompareDriver.sh` needs no argument.
* `MultiThreadECMC` needs one argument: the "number of active spheres".

### Output
* When called by `Scripts/CValidateECMC.sh` or `Scripts/Ordering.sh`
    * Print 'OK' if the first horizon violation occurs before liftings and reference liftings  start to differ.
    * Print 'Not OK' if the first horizon violation occurs after liftings and reference liftings start to differ.
    * `ValidateDriverC.sh` counts the number of 'OK' and 'Not OK', and output the number of bugs
* When called by `Scripts/Benchmark.sh`
    * Print number of events per hour (EPH)
    * `BenchDriver.sh` writes the EPHs [events per hour] into `EPH.dat`
    * `speed_up_simplify.py` makes a plot similar to Figure 7 of [Li2020](https://arxiv.org/abs/2004.11040), using `EPH.dat`, writing `EPH.pdf`.
* When called by `Scripts/Compare.sh`, `MultiThreadECMC` check the impact of atomic operation on single-thread performance. It prints EPH when having and not having atomic variables, respectively.
