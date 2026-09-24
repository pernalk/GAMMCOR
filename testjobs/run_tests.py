#! /usr/bin/env python
import os,sys,glob, string,math
from pathlib import Path
import numpy as np
from scripts import test_list, test_utils
from scripts.colors import *

########################################################################
# PREPARE TESTS
########################################################################

### Enter to test directory
SCRIPT_PATH = os.path.dirname(os.path.abspath(__file__)) + "/"
os.chdir(SCRIPT_PATH)

pyscf_basis_path = os.path.normpath(os.path.join(SCRIPT_PATH, "..", "basis-for-pyscf")) + "/"
basis_path = os.path.normpath(os.path.join(SCRIPT_PATH, "..", "basis")) + "/"

INPUTS = sorted(glob.glob("*/*/input.inp"))
PYSCF_INPUTS = [f for f in INPUTS if test_utils.get_interface(f) == "pyscf"]
BASIS_INPUTS = [f for f in INPUTS if test_utils.get_interface(f) != "pyscf"]

def clean():
    patterns = ["FOFO", "FFOO", "fock.dat", "fort*", "cholvecs", "AC0BLK*", "A0BLK*", "ure_casno.dat", "XY0"]
    for test in test_list.get():
        for path in test['units'] :
            rdm2 = test_utils.get_interface(os.path.join(path, "input.inp")) in ["molpro", "orca"]
            for pattern in patterns + (["rdm2.dat"] if rdm2 else []):
                files = glob.glob(os.path.join(path, pattern))
                for file in files:
                    os.remove(file)

if "--clean" in sys.argv[1:]:
    clean()
    exit(0)

### Set executable path
if len(sys.argv) < 2:
    raise Exception("Executable path not specified!")

local_path_gammcor = os.path.abspath(sys.argv[1])

if local_path_gammcor.endswith("/gammcor"):
    RUN_GAMMCOR = local_path_gammcor 
else:
    if local_path_gammcor.endswith("/"):
        RUN_GAMMCOR = local_path_gammcor + 'gammcor'
    else:
        RUN_GAMMCOR = local_path_gammcor + '/gammcor'
#print(RUN_GAMMCOR)
if not Path(RUN_GAMMCOR).is_file():
    raise Exception(f"{RUN_GAMMCOR} not found!")
### Set test level
TEST_LEVEL = 'short'
#
# other options:
# long
# full
#
if len(sys.argv) > 2:
    TEST_LEVEL = sys.argv[2]
print(f"{OKCYAN}TEST LEVEL: {TEST_LEVEL}{ENDC}")


########################################################################
# RUN TESTS
########################################################################
no_passed = no_failed = no_total = 0
tests_failed = []

test_utils.set_basis_path(PYSCF_INPUTS, pyscf_basis_path)
test_utils.set_basis_path(BASIS_INPUTS, basis_path)
try:
    for test in test_list.get():
        units = test_utils.prepare_test(test, TEST_LEVEL)
        if len(units) > 0:
            result = test_utils.start_test(test, units, RUN_GAMMCOR)
            no_total += result.get('passed') + result.get('failed')
            no_passed += result.get('passed')
            no_failed += result.get('failed')
            tests_failed += result.get('failed_tests')
            # if result.get('failed'):
            #     print(f"{FAIL}FAILED{ENDC}")
            # else:
            #     print(f"{OKGREEN}PASSED{ENDC}")
            # print("\n---------------------------------------------\n")
finally:
    test_utils.clear_basis_path(INPUTS)

########################################################################
# SUMMARY
########################################################################
if no_failed == 0:
    print(OKGREEN)
else:
    print(FAIL)

print("===============================================")
print(f"\n{no_passed} / {no_failed} / {no_total} (passed / failed / total)\n")

if no_failed == 0:
    print(f"*** SUCCESS ***")
    print(f"\n===============================================")
    # DELETE UNNECESSARY FILES
    clean()

    exit(0)
else:
    print(f"*** FAILURE ***")
    print("\nFailed tests:")
    for name in tests_failed:
        print(f"* {name}")
    print(f"\n===============================================")    
    exit(1)
