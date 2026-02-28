#!/bin/bash

GIT_DIR=$(git rev-parse --git-dir)

echo "Uninstalling hooks..."
rm $GIT_DIR/hooks/pre-commit
rm $GIT_DIR/hooks/pre-push
rm $GIT_DIR/hooks/post-merge
echo "Done!"