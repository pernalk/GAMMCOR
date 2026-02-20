# Plan: Ustawienie pełnego pipeline'u profilowania CPU dla GAMMCOR

## Kontekst

Celem jest przygotowanie infrastruktury do profilowania programu GAMMCOR (Fortran, quantum chemistry) za pomocą `perf`. Chcemy zebrać hardware counter metrics (`perf stat`), zidentyfikować hotspoty (`perf record`) i wygenerować flame graph — żeby ocenić bottlenecki i potencjał GPU acceleration. Na razie nie mamy danych testowych — skupiamy się na setupie narzędzi i skryptów.

### Stan systemu
- **OS:** Fedora 43, kernel 6.18.12-200.fc43.x86_64
- **CPU:** Intel Xeon W-3245M @ 3.20GHz (Skylake, 16C/32T, AVX-512)
- **perf:** NIE zainstalowany
- **perf_event_paranoid:** 2 (za wysoko, potrzeba ≤1)
- **gfortran:** NIE zainstalowany (pakiet `gcc-gfortran` dostępny w dnf)
- **ifort/ifx:** NIE zainstalowany
- **MKL:** dostępny w conda (`/home/hemik/miniconda3/lib/libmkl_*`)
- **Istniejący binary:** `gammcor` z debug_info, nie stripped, ale brakuje mu bibliotek MKL w runtime (ldd → not found)

---

## Krok 1: Instalacja narzędzi systemowych

```bash
# perf
sudo dnf install perf

# gfortran
sudo dnf install gcc-gfortran

# FlameGraph (do generowania SVG)
git clone https://github.com/brendangregg/FlameGraph.git ~/tools/FlameGraph
```

Weryfikacja:
```bash
perf --version
gfortran --version
ls ~/tools/FlameGraph/stackcollapse-perf.pl
```

## Krok 2: Ustawienie uprawnień perf

```bash
sudo sysctl kernel.perf_event_paranoid=1
```

Weryfikacja: `cat /proc/sys/kernel/perf_event_paranoid` → powinno wypisać `1`.

Opcjonalnie na stałe: `echo 'kernel.perf_event_paranoid=1' | sudo tee -a /etc/sysctl.conf`

## Krok 3: Konfiguracja budowania — dwa warianty kompilatora

### 3a: Nowy Makefile do profilowania (`Makefile.profile`)

Tworzymy nowy Makefile bazujący na `Makefile.gcc`, z flagami profilowymi:

**Pliki do modyfikacji:** Nowy plik `Makefile.profile`

```makefile
MKL_ROOT = /home/hemik/miniconda3/

FCC = gfortran
FFLAGS = -O2 -g -fno-omit-frame-pointer -march=native -fopenmp -I xcfun/fortran

MKL_LIB = -L$(MKL_ROOT)lib/ -lmkl_gf_lp64 -lmkl_sequential -lmkl_core
XCFUN_LIB = -L./xcfun/lib/ -lxcfun

LIBS = $(MKL_LIB) $(XCFUN_LIB)

include Makefile.common
```

Kluczowe różnice vs `Makefile.gcc`:
| Flaga | Zwykły build | Profiling build |
|---|---|---|
| Optymalizacja | `-O3` | `-O2` (lepsze mapowanie symboli) |
| Frame pointer | brak | `-fno-omit-frame-pointer` (perf stack walking) |
| Debug | `-g` | `-g` (zachowane) |
| Arch | `-march=skylake-avx512` | `-march=native` (auto-detect) |
| MKL_ROOT | `/opt/intel/mkl/` | `/home/hemik/miniconda3/` |

### 3b: Przygotowanie do Intel oneAPI (do zainstalowania osobno)

Po zainstalowaniu Intel oneAPI (`ifx`), analogiczny Makefile.profile.intel:

```makefile
MKL_ROOT = /home/hemik/miniconda3/

FCC = ifx
FFLAGS = -O2 -g -fno-omit-frame-pointer -xHost -qopenmp -I xcfun/fortran

MKL_LIB = -L$(MKL_ROOT)lib/ -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core
XCFUN_LIB = -L./xcfun/lib/ -lxcfun

LIBS = $(MKL_LIB) $(XCFUN_LIB)

include Makefile.common
```

## Krok 4: Build profilowanego executable

```bash
cd /home/hemik/Hobby/GAMMCOR
mkdir -p OBJ

# Najpierw xcfun (jeśli nie zbudowany)
cd xcfun && make -f Makefile.gcc && cd ..

# Wyczyść stary build i zbuduj z profiling flags
make -f Makefile.profile clean 2>/dev/null; make -f Makefile.profile
```

Weryfikacja: `file gammcor` → powinno zawierać `with debug_info, not stripped`.

## Krok 5: Skrypt profilujący (`scripts/profile_run.sh`)

Tworzymy skrypt automatyzujący cały pipeline. Nowy plik: `scripts/profile_run.sh`

```bash
#!/bin/bash
set -euo pipefail

# === KONFIGURACJA ===
GAMMCOR_BIN="${GAMMCOR_BIN:-./gammcor}"
WORKDIR="${1:-.}"         # katalog z input.inp i integrałami
RESULTS="${2:-results}"   # katalog na wyniki
FLAMEGRAPH_DIR="${FLAMEGRAPH_DIR:-$HOME/tools/FlameGraph}"
REPEAT=3

export LD_LIBRARY_PATH="/home/hemik/miniconda3/lib:${LD_LIBRARY_PATH:-}"

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
    "$GAMMCOR_BIN" 2> "$RESULTS/perf_stat_basic.txt"

# === KROK 2: perf stat — FP metrics (Intel Skylake+) ===
echo "[2/5] perf stat: FP metrics..."
if perf list 2>/dev/null | grep -q fp_arith_inst_retired; then
  perf stat -e fp_arith_inst_retired.scalar_double,\
fp_arith_inst_retired.128b_packed_double,\
fp_arith_inst_retired.256b_packed_double,\
fp_arith_inst_retired.512b_packed_double \
      "$GAMMCOR_BIN" 2> "$RESULTS/perf_stat_fp.txt"
else
  echo "FP events not available on this CPU" > "$RESULTS/perf_stat_fp.txt"
fi

# === KROK 3: perf stat — memory metrics ===
echo "[3/5] perf stat: memory metrics..."
perf stat -e dTLB-loads,dTLB-load-misses,\
iTLB-loads,iTLB-load-misses,\
bus-cycles,ref-cycles \
    "$GAMMCOR_BIN" 2> "$RESULTS/perf_stat_memory.txt"

# === KROK 4: perf stat — repeatability ===
echo "[4/5] perf stat: repeatability ($REPEAT runs)..."
perf stat -r "$REPEAT" -e cycles,instructions,cache-references,cache-misses,\
branches,branch-misses,LLC-loads,LLC-load-misses,task-clock \
    "$GAMMCOR_BIN" 2> "$RESULTS/perf_stat_repeated.txt"

# === KROK 5: perf record + flame graph ===
echo "[5/5] perf record + flame graph..."
perf record -g --call-graph dwarf -F 99 -o "$RESULTS/perf.data" -- "$GAMMCOR_BIN"

# Generuj flame graph
perf script -i "$RESULTS/perf.data" > "$RESULTS/perf_script.out"
"$FLAMEGRAPH_DIR/stackcollapse-perf.pl" "$RESULTS/perf_script.out" > "$RESULTS/collapsed.txt"
"$FLAMEGRAPH_DIR/flamegraph.pl" "$RESULTS/collapsed.txt" > "$RESULTS/flamegraph.svg"

# Generuj raport tekstowy
perf report -i "$RESULTS/perf.data" --stdio --no-children > "$RESULTS/perf_report.txt" 2>&1

echo ""
echo "=== GOTOWE ==="
echo "Wyniki w: $RESULTS/"
ls -lh "$RESULTS/"
```

## Krok 6: Skrypt przygotowania środowiska (`scripts/profile_env.sh`)

Nowy plik: `scripts/profile_env.sh` — uruchamiany przed profilowaniem.

```bash
#!/bin/bash
# Przygotowanie środowiska do profilowania
# Wymaga sudo

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
```

## Struktura plików do utworzenia

```
GAMMCOR/
├── Makefile.profile          # gfortran z flagami profilowymi
├── Makefile.profile.intel    # ifx z flagami profilowymi (po instalacji oneAPI)
├── scripts/
│   ├── profile_env.sh        # przygotowanie środowiska (sudo)
│   └── profile_run.sh        # pełny pipeline profilowania
```

## Weryfikacja

1. **Instalacja narzędzi:** `perf --version`, `gfortran --version`, `ls ~/tools/FlameGraph/flamegraph.pl`
2. **Build:** `make -f Makefile.profile` → binary z debug_info
3. **Dry run (bez danych):** `bash scripts/profile_run.sh` — powinien się uruchomić i zakończyć z błędem GAMMCOR (brak input.inp), ale perf powinien zaraportować choćby krótki run
4. **Pełny run (z danymi):** Gdy pojawią się dane testowe, uruchomić: `MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 bash scripts/profile_run.sh /path/to/test/dir results/test1`
5. **Flame graph:** Otworzyć `results/flamegraph.svg` w przeglądarce
