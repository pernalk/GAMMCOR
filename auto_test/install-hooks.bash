#!/bin/bash

GIT_DIR=$(git rev-parse --git-dir)

echo "Installing hooks..."
# this command creates symlink to our scripts
ln -s ../../auto_test/pre-commit.bash $GIT_DIR/hooks/pre-commit
ln -s ../../auto_test/pre-push.bash $GIT_DIR/hooks/pre-push
ln -s ../../auto_test/pre-push.bash $GIT_DIR/hooks/post-merge
echo "Done!"