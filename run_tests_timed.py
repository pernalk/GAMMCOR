#!/usr/bin/env python
import os, sys, time, subprocess

SCRIPT_PATH = os.path.dirname(os.path.abspath(__file__))
os.chdir(SCRIPT_PATH)
GAMMCOR = os.path.join(SCRIPT_PATH, "gammcor")

# Collect all test dirs that have input.inp
test_dirs = []
for root, dirs, files in os.walk("TESTS"):
    if "input.inp" in files:
        test_dirs.append(root)
test_dirs.sort()

# Optional filter from command line
if len(sys.argv) > 1:
    filt = sys.argv[1]
    test_dirs = [td for td in test_dirs if filt in td]

print(f"{'Test':<50} {'Time (s)':>10} {'Status'}")
print("-" * 75)

total_time = 0
passed = 0
failed = 0
crashed = 0

for td in test_dirs:
    cwd = os.getcwd()
    os.chdir(td)
    t0 = time.time()
    try:
        result = subprocess.run(
            [GAMMCOR],
            capture_output=True, timeout=300
        )
        elapsed = time.time() - t0
        total_time += elapsed
        if result.returncode == 0:
            status = "OK"
            passed += 1
        else:
            status = "CRASH"
            crashed += 1
    except subprocess.TimeoutExpired:
        elapsed = time.time() - t0
        total_time += elapsed
        status = "TIMEOUT"
        crashed += 1
    except Exception as e:
        elapsed = time.time() - t0
        total_time += elapsed
        status = f"ERROR: {e}"
        crashed += 1
    os.chdir(cwd)
    print(f"{td:<50} {elapsed:>10.2f} {status}")

print("-" * 75)
print(f"{'TOTAL':<50} {total_time:>10.2f}")
print(f"\nPassed: {passed}, Crashed: {crashed}, Total: {len(test_dirs)}")
