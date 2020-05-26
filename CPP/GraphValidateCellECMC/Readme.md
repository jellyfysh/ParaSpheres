## C++ cell-based simulation for 2D hard spheres (` GraphValidateCellECMC.cc`).

This code is based on a direct translation of the historic
Fortran90 code in `Fortran90/Historic/` into `C++` (see
Section 3 of [Li2020](https://arxiv.org/abs/2004.11040) in
[References.bib](References.bib)).

Particles move in a single direction, and explore the constraining polytope.

A given constraint graph is read into the program. The program checks that every collision in the 
cell-based simulation is compatible with the constraints contained in the constraint graph.

Physical parameters are set in the file `Param.h`.

Multi-dimensional arrays are managed using `boost::multi_array`, see [here](https://www.boost.org/).


### Requirements
`Make`, `C++`, `boost`

### Usage
* Compiled and called by `Scripts/ValidateG.sh`

### Input
* `configuration.dat`, `constraints.dat` and `info.dat` from `SetupData/`

### Output
* Terminate and print an error message if a lifting takes place which is not described by the constraint graph.
* Data to plot Figure 6(c) of [Li2020](https://arxiv.org/abs/2004.11040).
