#!/bin/bash

echo "Running pre-push hook"
gammcor_verify.bash 'full'

# $? stores exit value of the last command
if [ $? -ne 0 ]; then
 echo "Tests must pass before push!"
 exit 1
fi

git push --no-verify