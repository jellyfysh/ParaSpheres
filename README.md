[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v1.4%20adopted-ff69b4.svg)](CODE_OF_CONDUCT.md)


# JellyFysh/ParaSpheres code

This code project is associated with the paper "Multithreaded event-chain Monte
Carlo with local times", Botao Li, Synge Todo, A. C. Maggs, Werner Krauth
(see [Li2020](https://arxiv.org/abs/2004.11040) in [References.bib](References.bib)).

## Installing

The code project is written in Fortran 90, C++, Python, as well as shell.
The following steps are required for running and testing the code, namely:

* Creating initial configurations of hard spheres.
* Creating the constraint graph that is used in the simulations.
* Running reference simulations with variable numbers of threads in order to
validate the multithreaded algorithm.
* Running the multithreaded code in Validate and Benchmark modes.

All programs can be compiled and run individually. For ease of use, driver scripts are provided
in the `Scripts/` directory. Use of scripts ensures that output is placed
in a new directory, in the `Run/` directory tree. The scripts also suggest consistent options
for all the programs. The scripts also create symbolic links to the starting
data and initialization files.

## Running tests

The user should switch to the `Scripts` directory (`cd Scripts` on Linux) and
run from within this directory, then examine the results in subdirectories of
the  `Run/` directory. Each script creates a unique subdirectory.

## Using JellyFysh/ParaSpheres

Use of the programs requires firstly running `Setup.sh`, in order to create the
starting data files, which are placed in `SetupData/`. Other scripts can be run
in any order. The `SystemDimension` variable must be set to the same value
in each script, the default value is `SystemDimension=256` which is the standard value
used in the paper.


The directory `Fortran90/Historic/`, contains code
authored by Etienne P. Bernard, who kindly contributed the
programs used in the paper "Two-Step Melting in Two
Dimensions: First-Order Liquid-Hexatic Transition", Etienne P. Bernard and
Werner Krauth (see [Bernard2011](https://doi.org/10.1103/PhysRevLett.107.155704)
or [ArXiv](https://arxiv.org/abs/1102.4094)).

## Contributing

As an open-source project, the JeLLyFysh organization solicits contributions
from the community. Please read 
the [contribution guideline](CONTRIBUTING.md) for details.

If you find a bug, please raise an Issue here on GitHub to let us know.

Please note that this project is released with the Contributor Covenant [code of
conduct](CODE_OF_CONDUCT.md). By 
participating in this project you agree to abide by its terms. Report
unacceptable behavior to 
[werner.krauth@ens.fr](mailto:werner.krauth@ens.fr).

## License

This project is licensed under the GNU General Public License, version 3 (see
the [LICENSE](LICENSE) file).

## Contact

If you have questions regarding JeLLyFysh/ParaSpheres, just raise an issue
here on GitHub. We are happy to help you!

## Citation

If you use JeLLyFysh/ParaSpheres in published work, please cite the following reference (see
[Li2020](https://arxiv.org/abs/2004.11040) in [References.bib](References.bib)):

Botao Li, Synge Todo, A. C. Maggs, Werner Krauth,
"Multithreaded event-chain Monte Carlo with local times"
[arXiv e-print:  2004.11040 (2020)](https://arxiv.org/abs/2004.11040).

If you use the historic Fortran90 code of JeLLyFysh/ParaSpheres in published work, please cite the
following reference (see [Bernard2011](https://doi.org/10.1103/PhysRevLett.107.155704) in [References.bib](References.bib)):

Etienne P. Bernard and Werner Krauth, "Two-Step Melting in Two
Dimensions: First-Order Liquid-Hexatic Transition"
[Phys. Rev. Lett. 107, 155704 (2011)](https://doi.org/10.1103/PhysRevLett.107.155704).

