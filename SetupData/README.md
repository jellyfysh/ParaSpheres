## Data files
This directory contains data and configurations required by the programs that run from `Scripts/`.

### Provided test suites
* The `Data.4` directory contains a 4-sphere configuration, a list of active spheres, constraint graphs, and reference liftings. This test suite is used by `Scripts/Ordering.sh` and `Scripts/SequentialC4.sh`.
* The `Data.5` directory contains a 5-sphere configuration, a list of active spheres, constraint graphs, and reference liftings. This test suite is used by `Scripts/Ordering.sh` and `Scripts/SequentialC5.sh`.

### Generated test suits
Each run of  `Scripts/Setup.sh` creates a new subdirectory in the present directory. This subdirectory contains the initial configuration, constraint graphs and reference liftings required for the running of scripts.

The data for generated test suites is placed in `Data.$SystemDimension`, where `SystemDimension` is a parameter in `Script/Setup.sh`.
