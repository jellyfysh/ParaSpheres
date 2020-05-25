# Generation of initial state for the simulations

The `Fortran90` code in this directory slightly modifies the historic `Fortran90` code (in the `Fortran90/Historic` directory) (see [Bernard2011](https://doi.org/10.1103/PhysRevLett.107.155704) in [References.bib](References.bib)).

The code implements hard-sphere ECMC with spheres stored in a local cell system (spheres are not globally numbered).

### Requirements
`Fortran90` and `Make`.

### Usage
The code is configured to be run from `Scripts/Setup.sh`. In this case the final configuration files will be stored in the directory `SetupData/`.

* To compile the program, set the `FC` variable in the `Makefile` to the name of a `Fortran90` compiler and type `make`.
* Run by `./CellECMC argument`, where `argument` is a string specifying what this program will do.
    * If `argument` is `create_init`, this program creates a initial configuration using the cell-based ECMC algorithm
    * If `argument` is `benchmark`, this program checks that the created initial configuration does not induce three-body collision when using the cell-based ECMC algorithm.
* The snapshots should be reproducible between runs due to the initialization of the random number generator
in `resetseed.f90`. To generate different configurations between runs, the call to `resetseed` in `CellECMC.f90` (line 111) should be commented out.

### Input
The necessary parameters are written as variables in the `.f90` files. There are two variables that may be modified:

* `N_b`, number of spheres (`phys_parameters.f90`, line 4).
* `eta`, volume fraction (`phys_parameters.f90`, line 5).

The lateral box size is always `1`. The radius of spheres is calculated from `eta`. When running this program with argument 

* `create_init`, no input file is needed. When run from the `Scripts/` directory, the variables are automatically set by editing files using the scripted editor `sed`.
* `benchmark`, `init.dat`, which is the generated initial configuration, is needed.

### Output

* When called by `Scripts/Setup.sh` with `create_init`, two files will be created in `SetupData/Data.$SystemDimension`:
    * `init.dat`, the generated initial configuration.
    * `info.dat`, a file containing system size and size of the spheres. 
* With `benchmark`, a warning message is printed to standard output when a thee-body collision is detected.
