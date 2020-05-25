## `SingleThreadLocalTimeECMC.py`
This Python3 program implements Algorithm 4 (see [Li2020](https://arxiv.org/abs/2004.11040) in [References.bib](References.bib)).


### Requirements
Python (`Python3` or `PyPy`)

### Usage
* Called by `ValidateDriverPy.sh` in the present directory.
* `ValidateDriverPy.sh` is called by `PValidateECMC.sh`  in the `Scripts` directory.

### Input
* All the files in the test suite (test suite specified by scripts).
* `ValidateDriverPy.sh` needs the following two arguments:
   1 `log_2`(number of active spheres).
   2 number of runs for each number of active spheres.
* `PValidateECMC.sh` needs one argument. It's the number of active spheres.

### Output
* Print 'OK' if the first horizon condition violation takes place before the liftings and reference liftings start to differ.
* Print 'Not OK' if the first horizon condition violation takes place after the liftings and reference liftings start to differ.
* `ValidateDriverPy.sh` counts the number of 'OK' and 'Not OK', and output the number of 'Not OK'.
