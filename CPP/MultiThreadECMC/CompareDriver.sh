#!/bin/bash
set -e
export PATH=$PATH:.
OMP_NUM_THREADS=1 MultiThreadECMC comparison 1
