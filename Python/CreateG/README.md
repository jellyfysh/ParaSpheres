## `GenerateG3.py`, `PruneG.py`
These Python3 programs create and prune the constraint graphs, as described in Section 3.1 of [Li2020](https://arxiv.org/abs/2004.11040) in [References.bib](References.bib).


### Requirements
Python (`Python3` or `PyPy`).

### Usage
* Called by `Setup.sh` in the `Scripts` directory.

### Input
* `configuration.dat` and `info.dat` for `GenerateG3.py`.
* `G3_1.dat` and `G3_-1.dat` as well as `configuration.dat` and `info.dat` for `PruneG.py`.

### Output
* `GenerateG3.py` writes G^(3), for forward and backward moving spheres, into `G3_1.dat` and `G3_-1.dat` respectively.
* `PruneG.py` writes the pruned constraint graph into `G_pruned4s.dat` with the default setting in this program. In `G_pruned4s.dat`,
`4` means 4th order of pruning and `s` means symmetrized.
* `Setup.sh` creates a shortcut called `constraints.dat` which points towards `G_pruned4s.dat`. If the settings in `PruneG.py`
have been changed, the file name `G_pruned4s.dat` needs to be changed correspondingly in `Setup.sh`.
