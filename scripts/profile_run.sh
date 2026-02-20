#!/bin/bash
set -uo pipefail

# === LD_LIBRARY_PATH dla MKL ===
export LD_LIBRARY_PATH="/home/hemik/miniconda3/lib:${LD_LIBRARY_PATH:-}"

# === KONFIGURACJA ===
GAMMCOR_BIN="${GAMMCOR_BIN:-./gammcor}"
WORKDIR="${1:-.}"         # katalog z input.inp i integrałami
RESULTS="${2:-results}"   # katalog na wyniki
FLAMEGRAPH_DIR="${FLAMEGRAPH_DIR:-$HOME/tools/FlameGraph}"
REPEAT=3

mkdir -p "$RESULTS"
cd "$WORKDIR"

echo "=== Profilowanie GAMMCOR ==="
echo "Binary: $GAMMCOR_BIN"
echo "Workdir: $(pwd)"
echo "Results: $RESULTS"
echo ""

# === SYSTEM INFO ===
echo "Zbieranie informacji o systemie..."
{
  echo "=== CPU ==="
  lscpu | grep -E 'Model name|Socket|Core|Thread|MHz|Cache|Architecture'
  echo ""
  echo "=== Memory ==="
  free -h
  echo ""
  echo "=== Kernel ==="
  uname -r
  echo ""
  echo "=== Compiler ==="
  gfortran --version 2>/dev/null || ifx --version 2>/dev/null || echo "unknown"
  echo ""
  echo "=== BLAS/LAPACK ==="
  ldd "$GAMMCOR_BIN" | grep -iE 'mkl|blas|lapack' || echo "not found in ldd"
  echo ""
  echo "=== GPU ==="
  nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv 2>/dev/null || echo "no NVIDIA GPU"
} > "$RESULTS/system_info.txt"

# === KROK 1: perf stat — podstawowe metryki ===
echo "[1/5] perf stat: basic metrics..."
perf stat -e cycles,instructions,cache-references,cache-misses,\
branches,branch-misses,L1-dcache-loads,L1-dcache-load-misses,\
LLC-loads,LLC-load-misses,task-clock,context-switches,cpu-migrations,page-faults \
    "$GAMMCOR_BIN" 2> "$RESULTS/perf_stat_basic.txt" || true

# === KROK 2: perf stat — FP metrics (Intel Skylake+) ===
echo "[2/5] perf stat: FP metrics..."
if perf list 2>/dev/null | grep -q fp_arith_inst_retired; then
  perf stat -e fp_arith_inst_retired.scalar_double,\
fp_arith_inst_retired.128b_packed_double,\
fp_arith_inst_retired.256b_packed_double,\
fp_arith_inst_retired.512b_packed_double \
      "$GAMMCOR_BIN" 2> "$RESULTS/perf_stat_fp.txt" || true
else
  echo "FP events not available on this CPU" > "$RESULTS/perf_stat_fp.txt"
fi

# === KROK 3: perf stat — memory metrics ===
echo "[3/5] perf stat: memory metrics..."
perf stat -e dTLB-loads,dTLB-load-misses,\
iTLB-loads,iTLB-load-misses,\
bus-cycles,ref-cycles \
    "$GAMMCOR_BIN" 2> "$RESULTS/perf_stat_memory.txt" || true

# === KROK 4: perf stat — repeatability ===
echo "[4/5] perf stat: repeatability ($REPEAT runs)..."
perf stat -r "$REPEAT" -e cycles,instructions,cache-references,cache-misses,\
branches,branch-misses,LLC-loads,LLC-load-misses,task-clock \
    "$GAMMCOR_BIN" 2> "$RESULTS/perf_stat_repeated.txt" || true

# === KROK 5: perf record + flame graph ===
echo "[5/5] perf record + flame graph..."
perf record -g --call-graph dwarf -F 99 -o "$RESULTS/perf.data" -- "$GAMMCOR_BIN" || true

# Generuj flame graph (tylko jeśli perf.data istnieje)
if [ -f "$RESULTS/perf.data" ]; then
  perf script -i "$RESULTS/perf.data" > "$RESULTS/perf_script.out"
  "$FLAMEGRAPH_DIR/stackcollapse-perf.pl" "$RESULTS/perf_script.out" > "$RESULTS/collapsed.txt"
  "$FLAMEGRAPH_DIR/flamegraph.pl" "$RESULTS/collapsed.txt" > "$RESULTS/flamegraph.svg"

  # Generuj raport tekstowy
  perf report -i "$RESULTS/perf.data" --stdio --no-children > "$RESULTS/perf_report.txt" 2>&1
else
  echo "perf.data not generated — skipping flame graph" >&2
fi

echo ""
echo "=== GOTOWE ==="
echo "Wyniki w: $RESULTS/"
ls -lh "$RESULTS/"
