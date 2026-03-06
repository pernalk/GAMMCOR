#!/usr/bin/env python3
"""
Generate a test report table comparing computed energies against references.

Usage:
    python test_report.py [FILTER]    # run tests and generate report
    python test_report.py --no-run    # report from existing gammcor_test.out files
    python test_report.py AC0         # only tests matching "AC0"

Output: test_report.md (markdown table) + terminal summary
"""

import os
import sys
import subprocess
import re
import time
from pathlib import Path

SCRIPT_DIR = Path(os.path.dirname(os.path.abspath(__file__)))
TESTS_DIR = SCRIPT_DIR / "TESTS"
GAMMCOR = SCRIPT_DIR / "gammcor"
TIMEOUT = 300
REPORT_FILE = SCRIPT_DIR / "test_report.md"
TEST_OUTPUT = "gammcor_test.out"  # written in each test dir after run

# Energy extraction patterns: (grep_string, value_position, label)
ENERGY_PATTERNS = [
    ("ECASSCF+ENuc, AC0-Corr, AC0-CASSCF",       -1, "AC0-CASSCF"),
    ("ECASSCF+ENuc, AC-Corr, AC-ERPA-CASSCF",     -1, "AC-ERPA-CASSCF"),
    ("ECASSCF+ENuc, AC1-Corr, ERPA-CASSCF",       -1, "AC1-ERPA-CASSCF"),
    ("ECASSCF+ENuc, ACn-Corr, ACn-CASSCF",        -1, "ACn-CASSCF"),
    ("EGVB+ENuc, 0th+1st-order ECorr, AC0-GVB",   -1, "AC0-GVB"),
    ("EGVB+ENuc, Corr, AC-ERPA-GVB",              -1, "AC-ERPA-GVB"),
    ("EGVB+ENuc, Corr, ERPA-GVB",                 -1, "AC1-ERPA-GVB"),
    ("EGVB + ENuc + 1,2-body",                    -1, "EERPA-GVB"),
    ("PiDFT Correlation",                         -1, "PiDFT"),
    ("Dexcitation correction for AC0",            -1, "AC0D-Dexcit"),
    ("CASSCF+ENuc, AC0-CBS[H], Total",             -1, "AC0-CBS[H]"),
    ("CASSCF+ENuc, AC0-CBS[DFT], Total",           -1, "AC0-CBS[DFT]"),
    ("ECASSCF+ENuc, AC1-Corr, AC1-CASSCF",         -1, "AC1-CASSCF"),
    ("lrCASSCF+srDF+ENuc, lrAC0-Corr, Total",      -1, "lrAC0CAS"),
]

# SAPT component tags
SAPT_TAGS = [
    "E1elst", "E1exch(S2)", "E2ind", "E2disp",
    "E2exch-ind", "E2exch-disp", "Eint(SAPT2)", "Eint(SAPT0)",
    "E2ind(unc)", "E2disp(unc)", "E2exch-ind(unc)", "E2exch-disp(unc)",
    "E2disp(CAS)", "E2disp(SCALED)", "E2disp(sp)", "E2disp(sc)",
]


def extract_energies_from_text(text):
    """Extract energy values from gammcor output text."""
    results = {}
    lines = text.splitlines()

    # Standard energy patterns
    for pattern, pos, label in ENERGY_PATTERNS:
        for line in lines:
            if pattern in line:
                nums = re.findall(r'-?\d+\.\d+', line)
                if nums:
                    results[label] = float(nums[pos])
                break

    # SAPT components — match longest tag first to avoid partial matches
    sorted_tags = sorted(SAPT_TAGS, key=len, reverse=True)
    for tag in sorted_tags:
        if tag in results:
            continue
        for line in lines:
            # Match tag followed by spaces and '=', but NOT as prefix of longer tag
            # e.g. "E2disp" should not match "E2disp(unc)" or "E2disp(SCALED)"
            escaped = re.escape(tag)
            # Negative lookahead: tag must not be followed by '(' (unless tag already has it)
            if '(' in tag:
                pattern = rf'^\s*{escaped}\s*='
            else:
                pattern = rf'^\s*{escaped}(?!\()\s*='
            if re.search(pattern, line):
                parts = line.split('=')
                if len(parts) >= 2:
                    nums = re.findall(r'-?\d+\.\d+', parts[-1])
                    if nums:
                        results[tag] = float(nums[0])
                        break

    return results


def extract_energies(filepath):
    """Extract energy values from a gammcor output file."""
    if not filepath.exists():
        return {}
    return extract_energies_from_text(filepath.read_text(errors='replace'))


def find_test_dirs(filter_str=None):
    """Find all test directories containing input.inp."""
    dirs = []
    for root, _, files in os.walk(TESTS_DIR):
        if "input.inp" in files:
            dirs.append(Path(root))
    dirs.sort()
    if filter_str:
        dirs = [d for d in dirs if filter_str in str(d)]
    return dirs


def run_test(test_dir):
    """Run gammcor in a test directory. Returns (status, elapsed, stdout_text)."""
    t0 = time.time()
    try:
        result = subprocess.run(
            [str(GAMMCOR)],
            cwd=test_dir,
            timeout=TIMEOUT,
            capture_output=True,
            text=True,
        )
        elapsed = time.time() - t0
        stdout = result.stdout or ""

        # Save output for --no-run reuse
        out_path = test_dir / TEST_OUTPUT
        out_path.write_text(stdout)

        if result.returncode == 0:
            return "OK", elapsed, stdout
        else:
            return "CRASH", elapsed, stdout
    except subprocess.TimeoutExpired:
        elapsed = time.time() - t0
        return "TIMEOUT", elapsed, ""
    except Exception as e:
        elapsed = time.time() - t0
        return "ERROR", elapsed, ""


def compare_energies(ref_energies, test_energies):
    """Compare two energy dicts, return list of (label, ref, test, diff, ok)."""
    comparisons = []
    for label in ref_energies:
        ref_val = ref_energies[label]
        if label in test_energies:
            test_val = test_energies[label]
            diff = abs(ref_val - test_val)
            tol = 1e-5 if label in SAPT_TAGS else 1e-7
            ok = diff <= tol
            comparisons.append((label, ref_val, test_val, diff, ok))
        else:
            comparisons.append((label, ref_val, None, None, False))
    return comparisons


def generate_report(test_dirs, run_tests=True):
    """Generate the full test report."""
    rows = []

    for i, test_dir in enumerate(test_dirs):
        rel = test_dir.relative_to(TESTS_DIR)
        ref_file = test_dir / "gammcor.out"
        test_out_file = test_dir / TEST_OUTPUT

        # Run test or use cached output
        if run_tests:
            status, elapsed, stdout = run_test(test_dir)
            test_energies = extract_energies_from_text(stdout)
            # Progress
            print(f"  [{i+1}/{len(test_dirs)}] {str(rel):<50} {elapsed:>8.3f}s  {status}")
        else:
            elapsed = None
            if test_out_file.exists():
                status = "cached"
                test_energies = extract_energies(test_out_file)
            else:
                status = "no_output"
                test_energies = {}

        # Extract reference energies
        ref_energies = extract_energies(ref_file) if ref_file.exists() else {}

        # Compare
        if ref_energies and test_energies:
            comparisons = compare_energies(ref_energies, test_energies)
        elif ref_energies and not test_energies:
            comparisons = [(l, v, None, None, False) for l, v in ref_energies.items()]
        else:
            comparisons = []

        # Determine energy status
        if not comparisons:
            energy_status = "—"
        elif all(c[4] for c in comparisons):
            energy_status = "MATCH"
        elif any(c[3] is not None and not c[4] for c in comparisons):
            energy_status = "DIFF"  # real value differences
        else:
            energy_status = "MISSING"  # ref exists but no computed value

        # Max diff
        diffs = [c[3] for c in comparisons if c[3] is not None]
        max_diff = max(diffs) if diffs else None

        # Primary energy label and value
        primary = comparisons[0] if comparisons else None

        rows.append({
            'test': str(rel),
            'status': status,
            'elapsed': elapsed,
            'energy_status': energy_status,
            'max_diff': max_diff,
            'primary': primary,
            'comparisons': comparisons,
            'has_ref': ref_file.exists(),
        })

    return rows


def format_diff(diff):
    """Format energy difference in appropriate units."""
    if diff is None:
        return "—"
    if diff < 1e-9:
        return "0.0"
    elif diff < 1e-6:
        return f"{diff*1e6:.3f} µHa"
    elif diff < 1e-3:
        return f"{diff*1e3:.4f} mHa"
    else:
        return f"{diff:.6f} Ha"


def write_report(rows, filename):
    """Write markdown report."""
    try:
        branch = subprocess.check_output(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"],
            text=True).strip()
        commit = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            text=True).strip()
        commit_msg = subprocess.check_output(
            ["git", "log", "-1", "--format=%s"],
            text=True).strip()
    except Exception:
        branch, commit, commit_msg = "?", "?", "?"

    from datetime import datetime
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")

    n_ok = sum(1 for r in rows if r['status'] == 'OK')
    n_crash = sum(1 for r in rows if r['status'] in ('CRASH', 'TIMEOUT', 'ERROR'))
    n_match = sum(1 for r in rows if r['energy_status'] == 'MATCH')
    n_diff = sum(1 for r in rows if r['energy_status'] == 'DIFF')
    n_missing = sum(1 for r in rows if r['energy_status'] == 'MISSING')
    n_no_ref = sum(1 for r in rows if r['energy_status'] == '—')
    n_total = len(rows)
    total_time = sum(r['elapsed'] for r in rows if r['elapsed'] is not None)

    with open(filename, 'w') as f:
        f.write(f"# GAMMCOR Test Report\n\n")
        f.write(f"- **Date**: {timestamp}\n")
        f.write(f"- **Branch**: `{branch}` @ `{commit}` — {commit_msg}\n")
        f.write(f"- **Run**: {n_ok} OK, {n_crash} CRASH / {n_total} total")
        if total_time > 0:
            f.write(f" ({total_time:.1f}s)")
        f.write(f"\n")
        f.write(f"- **Energy**: {n_match} MATCH, {n_diff} DIFF, "
                f"{n_missing} MISSING, {n_no_ref} no ref\n\n")

        # Main table
        f.write("| Test | Run | Time | Energy | Max Diff | Primary Value |\n")
        f.write("|------|-----|------|--------|----------|---------------|\n")

        for r in rows:
            test = r['test']
            status = r['status']
            e_status = r['energy_status']
            diff_str = format_diff(r['max_diff'])
            time_str = f"{r['elapsed']:.3f}s" if r['elapsed'] is not None else "—"

            if r['primary']:
                label, ref_val, test_val, _, ok = r['primary']
                if test_val is not None:
                    primary_str = f"{label}: {test_val:.8f}"
                else:
                    primary_str = f"{label}: —"
            else:
                primary_str = "—"

            f.write(f"| {test} | {status} | {time_str} | {e_status} | {diff_str} "
                    f"| {primary_str} |\n")

        # Detailed diffs (only real energy differences)
        real_mismatches = [r for r in rows if r['energy_status'] == 'DIFF']
        if real_mismatches:
            f.write(f"\n## Energy Differences (values differ beyond tolerance)\n\n")
            for r in real_mismatches:
                f.write(f"### {r['test']}\n\n")
                f.write("| Component | Reference | Computed | Diff |\n")
                f.write("|-----------|-----------|----------|------|\n")
                for label, ref_val, test_val, diff, ok in r['comparisons']:
                    tv = f"{test_val:.8f}" if test_val is not None else "—"
                    d = format_diff(diff)
                    flag = "" if ok else " **!**"
                    f.write(f"| {label} | {ref_val:.8f} | {tv} | {d}{flag} |\n")
                f.write("\n")

        # Crashes
        crashes = [r for r in rows if r['status'] in ('CRASH', 'TIMEOUT', 'ERROR')]
        if crashes:
            f.write(f"\n## Crashed/Failed Tests\n\n")
            for r in crashes:
                f.write(f"- `{r['test']}` — {r['status']}\n")
            f.write("\n")

    print(f"\nReport written to {filename}")


def main():
    os.chdir(SCRIPT_DIR)

    args = sys.argv[1:]
    run = True
    filter_str = None

    for arg in args:
        if arg == '--no-run':
            run = False
        else:
            filter_str = arg

    test_dirs = find_test_dirs(filter_str)
    if not test_dirs:
        print("No test directories found.")
        sys.exit(1)

    label = f" matching '{filter_str}'" if filter_str else ""
    print(f"Found {len(test_dirs)} tests{label}")

    if run:
        print("Running tests...\n")
    else:
        print("Using cached gammcor_test.out files\n")

    rows = generate_report(test_dirs, run_tests=run)
    write_report(rows, str(REPORT_FILE))

    # Terminal summary
    n_ok = sum(1 for r in rows if r['status'] == 'OK')
    n_crash = sum(1 for r in rows if r['status'] in ('CRASH', 'TIMEOUT', 'ERROR'))
    n_match = sum(1 for r in rows if r['energy_status'] == 'MATCH')
    n_diff = sum(1 for r in rows if r['energy_status'] == 'DIFF')
    n_missing = sum(1 for r in rows if r['energy_status'] == 'MISSING')
    print(f"Run: {n_ok} OK, {n_crash} CRASH")
    print(f"Energy: {n_match} MATCH, {n_diff} DIFF, {n_missing} MISSING")


if __name__ == "__main__":
    main()
