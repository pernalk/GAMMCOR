#!/bin/bash
# Run GAMMCOR tests
# Usage: ./test.sh [FILTER]
#   FILTER  optional substring to run only matching tests (e.g. "AC0" or "INCORE")
#
# Examples:
#   ./test.sh              — run all tests
#   ./test.sh AC0/TEST2    — run only AC0/TEST2 tests
#   ./test.sh INCORE       — run only INCORE tests
#   ./test.sh SAPT         — run only SAPT tests

cd "$(dirname "$0")"
source env.sh

python run_tests_timed.py "$@"
