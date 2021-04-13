#!/bin/bash

echo "Running post-merge hook"
make clean
make -j 8
gammcor_verify.bash 'full'

# $? stores exit value of the last command
if [ $? -ne 0 ]; then
 echo "Tests failed!"
 exit 1
fi