#!/bin/bash

# if any command inside script returns error, exit and return that error 
set -e

# magic line to ensure that we're always inside the root of our application,
# no matter from which directory we'll run script
# thanks to it we can just enter `./scripts/run-tests.bash`
cd "${0%/*}"

# set paths
TEST_SCRIPT_PATH=~/QCHEM/TESTS/run_tests.py
GAMMCOR_PATH=$PWD/gammcor

# run tests
echo Running $TEST_SCRIPT_PATH for $GAMMCOR_PATH
# $1 is a testing level
$TEST_SCRIPT_PATH $GAMMCOR_PATH $1
#srun --time=30:00 --mem=5G $TEST_SCRIPT_PATH $GAMMCOR_PATH