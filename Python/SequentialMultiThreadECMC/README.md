## `SequentialMultiThreadECMC.py`
Program for algorithm 5 in [arXiv:2004.11040](https://arxiv.org/abs/2004.11040), written in python.

This Python3 program implements Algorithm 5 (see [Li2020](https://arxiv.org/abs/2004.11040) in [References.bib](References.bib))).


### Requirements
`Python3` (with `Numpy`, `Matplotlib`).

### Usage
Called by `SequentialC4.sh` and `SequentialC5.sh` in the `Scripts` directory.

### Input
`configuration.dat` and `info.dat` in the test suite (test suite specified by scripts).

### Output
Writes to standard output:

* terminate state
* number of terminate states
* buffers that lead to abort
* number of transient states

Writes a simplified figure 4 (see [Li2020](https://arxiv.org/abs/2004.11040)) to `degenerate_table.pdf`.
