## `GlobalTimeECMC.py`
This Python3 program implements Algorithm 2 (see [Li2020](https://arxiv.org/abs/2004.11040) in [References.bib](References.bib)).

### Requirements
Python (`Python3` or `PyPy`).

### Usage
* Called by `Setup.sh` in the `Scripts` directory.

### Input
* `configuration.dat`, `constraints.dat`, `info.dat` ,and potentially `active*.dat` in the test suite (test suite specified by scripts).
* `active*.dat` will be read in if exists.

### Output
* write reference liftings with `x` active spheres into `liftingx.dat`
* If `activex.dat` doesn't exist, a list of `x` active spheres will be generated and written into `activex.dat`
