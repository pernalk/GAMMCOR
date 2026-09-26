"""Run gammcor on one testjob in a temporary directory and check the energy.

    python run_test.py <testjob dir>

Used by the tests in conda/recipe.yaml. The reference is the one in
testjobs/scripts/test_list.py (AC0/ethene, atol 1e-5).
"""
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REFERENCE = {
    # testjob: (line label in gammcor output, reference AC0-CASSCF energy)
    "ethene": (" ECASSCF+ENuc, AC0-Corr, AC0-CASSCF", -77.45880256),
}
ATOL = 1e-5

job = Path(sys.argv[1]).resolve()
label, ref = REFERENCE[job.name]

with tempfile.TemporaryDirectory() as tmp:
    work = Path(tmp) / job.name
    shutil.copytree(job, work, ignore=shutil.ignore_patterns("*.out"))
    run = subprocess.run(["gammcor"], cwd=work, capture_output=True, text=True)
    out = run.stdout
    if run.returncode != 0:
        print(out[-3000:], run.stderr[-3000:], sep="\n")
        sys.exit(f"gammcor exited with {run.returncode}")

line = next((l for l in out.splitlines() if l.startswith(label)), None)
if line is None:
    print(out[-3000:])
    sys.exit(f"'{label.strip()}' not found in the gammcor output")

energy = float(line.split()[-1])
print(f"{job.name}: {energy:.8f}  reference {ref:.8f}  diff {energy - ref:.1e}")
if abs(energy - ref) > ATOL:
    sys.exit(f"energy differs from the reference by more than {ATOL}")
