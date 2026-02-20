#!/bin/bash
# Przygotowanie środowiska do profilowania GAMMCOR
# Wymaga sudo dla ustawień systemowych

set -euo pipefail

echo "=== Przygotowanie środowiska do profilowania ==="
echo ""

echo "Ustawiam perf_event_paranoid..."
sudo sysctl kernel.perf_event_paranoid=1

echo "Ustawiam CPU governor na performance..."
sudo cpupower frequency-set -g performance 2>/dev/null || echo "cpupower niedostępny"

echo "Wyłączam turbo boost (Intel)..."
echo 1 | sudo tee /sys/devices/system/cpu/intel_pstate/no_turbo 2>/dev/null || echo "intel_pstate niedostępny"

echo "Ustawiam MKL_NUM_THREADS=1..."
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

echo ""
echo "Środowisko gotowe. Uruchom profiling:"
echo "  MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 bash scripts/profile_run.sh <workdir> <results>"
