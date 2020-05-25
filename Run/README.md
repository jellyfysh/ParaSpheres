## Run-time files

When running programs from  `Scripts/` extra run-time directories are created here, containing multiple data files.
Generated data is separated from the code.
The `Makefile` is intended to remove ALL run time data, destructively: `make clean` removes all the directories and data.

Directories have the names `Bench.*.*/`, `CValidateECMC.*.*/`, `Ordering.*.*/`, `PValidateECMC.*.*/`, `SequentialC4.*.*/`, `SequentialC5.*.*/`, `ValidateG.*.*/`, `Compare.*.*/`, `Data.run`.
