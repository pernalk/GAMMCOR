# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

GAMMCOR is a Fortran-based computational chemistry program for advanced electronic structure calculations including SAPT (Symmetry-Adapted Perturbation Theory), AC/AC0/AC1 (Adiabatic Connection), ERPA, EERPA, and CASPiDFT methods. The internal program name is PRDMFT (Pair and Reduced Density Matrix Functional Theory).

## Build Commands

**Prerequisites:** Intel MKL, and either Intel Fortran (`ifort`) or GNU Fortran (`gfortran`). The `env.sh` script configures `MKL_ROOT`, `LD_LIBRARY_PATH`, `PKG_CONFIG_PATH`, and thread settings. Edit `env.sh` to match your MKL installation path before building.

```bash
# Environment setup (edit env.sh first to set your MKL_ROOT path)
source env.sh

# First-time setup
mkdir OBJ

# Build XCFun library (do this first)
cd xcfun && make           # Intel compilers
cd xcfun && make -f Makefile.gcc  # GCC

# Build GAMMCOR
make                       # Intel (uses Makefile -> Makefile.common)
make -f Makefile.gcc       # GCC (uses Makefile.gcc -> Makefile.common)

# Clean build artifacts
make clean
```

Output executable: `./gammcor`

### Convenience Scripts

```bash
./build.sh              # build with GCC (sources env.sh automatically)
./build.sh -c           # clean + build
```

## Running Tests

### Quick run with timing (recommended)

```bash
./test.sh                  # run all tests (sources env.sh automatically)
./test.sh AC0/TEST2        # run only matching tests
./test.sh FOFO_RAM         # run only FOFO_RAM variants
./test.sh SAPT             # run only SAPT tests
```

This invokes `run_tests_timed.py`, which walks `TESTS/` for directories containing `input.inp`, runs gammcor in each, and reports pass/crash/timeout status with execution times. Timeout per test: 300 seconds.

### Full verification with energy comparison

```bash
source env.sh
python test_report.py              # run all tests, compare energies, generate test_report.md
python test_report.py CHOLESKY     # only CHOLESKY tests
python test_report.py --no-run     # report from cached gammcor_test.out files
```

`test_report.py` runs gammcor in each test directory, captures stdout to `gammcor_test.out`, compares computed energies against reference `gammcor.out`, and generates `test_report.md` with a full table of results. Statuses: MATCH (within tolerance), DIFF (real value difference), MISSING (no computed energy).

**After every code change:** always run `test_report.py`, record execution times (to 1ms precision) and energy values, and compare against previous results. This ensures no regressions are introduced. Include the timing/energy comparison table in commit messages or notes when relevant.

### Test variant convention

Each base test directory (e.g., `TESTS/AC0/TEST2/`) contains the default configuration. Variants are added as subdirectories with their own `input.inp` and `gammcor.out`. Integral data files are symlinked from the parent. **Old tests are never modified; new variants are added as subdirectories.**

Common variant subdirectories:
- `INCORE_INTEG/` — `TwoMoInt INCORE` (all integrals in memory)
- `CHOLESKY/` — `Cholesky .true.` (Cholesky decomposition)
- `FOFO_RAM/` — `FOFO_RAM TRUE` (FOFO integrals in RAM)
- `FOFO_DISK/` — `FOFO_RAM FALSE` (explicit disk-based FOFO)
- `FROZEN_CORE/` — frozen core approximation

## Architecture

### Source Organization

All source code is in `SOURCE/`. The codebase mixes fixed-format Fortran 77 (`.f`) with free-format Fortran 90 (`.f90`). Global state is shared via common blocks defined in `SOURCE/commons.inc`.

### Module Dependency Chain

```
types.f90  (data types: InputData, FlagsData, SystemBlock, SaptData)
  ├── inputfill.f90  (parses input.inp files)
  ├── systemdef.f90  (initializes molecular system; uses fofo_data)
  ├── sorter.f90, tran.f90  (integral sorting/transformation)
  ├── abmats.f90, abfofo.f90  (AB matrix construction, depend on tran)
  ├── ab0fofo.f90  (zero-order AB FOFO; uses chol_data, fofo_data)
  ├── diis.f90  (DIIS convergence acceleration)
  ├── sapt_utils.f90 → sapt_pol.f90, sapt_exch.f90 → sapt_main.f90
  └── exmisc.f90 → exappr.f90, exdpino.f90, srefex.f90

fofo_data.f90  (in-memory FOFO/FFOO: IntsFOFO, IntsFFOO, IFOFO_ram flag)
  └── used by: initia.f, systemdef.f90, abfofo.f90, ab0fofo.f90, ac_iter.f

chol_data.f90  (in-memory Cholesky: CholVecsFF, NCholesky_stored)
  └── used by: initia.f, abchol.f90, ab0fofo.f90, ac_fofo.f90, ac_fofo_min.f90, ac_iter.f

batch_dgemm.F90  (batched DGEMM wrapper: MKL / OpenMP / CUDA backends)
  └── used by: tran.f90 (tran4_gen, tran4_gen_incore)
```

### Program Flow (mainp.f)

1. `read_Input` → `check_Calc` → `fill_Flags` → `create_System` (input parsing)
2. Routes to SAPT (`sapt_driver`) or loads integrals (`LdInteg`/`ReadDAL`)
3. Dispatches to calculation method: DMSCF, CASPiDFT, CASPiDFTOPT, or VV10

### Key Subsystems

- **SAPT** (`sapt_main.f90`, `sapt_exch.f90`, `sapt_pol.f90`, `sapt_utils.f90`): Intermolecular interaction energy decomposition. Supports SAPT0 and SAPT2 levels.
- **AC methods** (`acfd.f`, `accas.f`, `ac_fofo.f90`, `ac_iter.f`): Adiabatic connection correlation energy. AC0 is linearized (MP2-like), AC uses full numerical integration.
- **ERPA** (`erpa.f`, `interpa.f`): Extended RPA response equations.
- **Integral handling** (`abmats.f90`, `abfofo.f90`, `ab0fofo.f90`, `tran.f90`, `sorter.f90`): Two-electron integrals stored as (FF|OO) and (FO|FO) in chemists' notation. Three storage modes via `TwoMoInt`:
  - `INCORE` (default): Full 4-index transformation in memory via `tran4_gen_incore`.
  - `FOFO`: Half-transformed integrals as FOFO/FFOO. With `FOFO_RAM TRUE`, arrays from `fofo_data` module are in memory; with `FOFO_RAM FALSE` (default), read from/written to disk.
  - `Cholesky`: Cholesky decomposition vectors in `chol_data` module. Accuracy: DEFAULT/TIGHT/LUDICROUS.
- **Data modules** (`fofo_data.f90`, `chol_data.f90`): Allocatable arrays for in-memory integral storage.
- **Batched BLAS** (`batch_dgemm.F90`): Portable wrapper for batched DGEMM. Uses MKL `dgemm_batch_strided` when `-DUSE_MKL_BATCH` is defined, otherwise falls back to OpenMP loop over standard `dgemm`. CUDA skeleton via `USE_CUDA_BATCH` for future GPU support.
- **DFT/XCFun** (`xcfun.f90`, `dftgrid.f`, `caspidft.f`): Short-range DFT functionals and CASPiDFT. XCFun library in `xcfun/` provides exchange-correlation functionals.
- **Exchange contributions** (`exmisc.f90`, `exdpino.f90`, `exi.f90`, `exappr.f90`): Exchange energy terms for SAPT and DMFT functionals.

### External Interfaces

The program reads integrals from external quantum chemistry codes configured via `InterfaceType`:
- **MOLPRO** (default): reads `AOONEINT.mol`, `AOTWOINT.mol`, `2RDM`, `MOLPRO.MOPUN`
- **DALTON**: reads `SIRIUS.RST`, `AOONEINT`, `AOTWOINT`
- **ORCA**: alternative integral source

### Conventions

- `IFun` flag selects the long-range DMFT functional (0=none/DFT, 1-13=various functionals)
- `IFunSR` flag selects the short-range DFT functional (0=none, 1=srLDA, 2=srPBE, 5-7=CASPiDFT/VV10)
- `TwoMoInt` flag selects integral storage: `INCORE` (1, default), `FFFF` (2), `FOFO` (3)
- `FOFO_RAM` input keyword (default `.false.`): when true, FOFO/FFOO integrals held in memory instead of disk files
- `Cholesky` input keyword (default `.false.`): enables Cholesky decomposition of two-electron integrals. Vectors stored in-memory (`CholVecsFF`) and accessed via pointer (`MatFF => CholVecsFF`) — no disk I/O for `cholvecs` file in non-SAPT paths. Known limitations: incompatible with frozen core (`NCoreOrb`) and LR-ERF.
- `Cholesky_Accuracy` input keyword: `DEFAULT` (TraceError 1e-2), `TIGHT` (1e-3, ~10× more accurate, +30% cost), `LUDICROUS` (1e-4, ~30× more accurate, +70% cost)
- `NBasis` = number of basis functions; `NELE` = number of electrons
- Integer flags use Fortran convention: 0=false/off, 1=true/on
- Constants are in `SOURCE/constants.h`; common blocks in `SOURCE/commons.inc`
- Build produces `.o` files in `OBJ/` and `.mod` files in the project root
- `env.sh` must be sourced before building or running tests
- Test convention: old tests are never modified; new variants are added as subdirectories with symlinks to parent data files
- `-DUSE_MKL_BATCH` compile flag: enables MKL batched DGEMM in integral transformations. Without it, falls back to OpenMP loop (works with any BLAS)
- `-DUSE_CUDA_BATCH` compile flag (future): enables GPU batched DGEMM via cuBLAS. Requires `-lcublas -lcudart`
- `.F90` files (uppercase) use Fortran preprocessor; `.f90` files do not

## Notatki o trudnościach (Issue Tracking)

Podczas każdej sesji pracy prowadź notatki o **wszystkich** napotkanych trudnościach — błędach kompilacji, segfaultach, regresji testów, problemach z konfiguracją, niezgodnościach interfejsów, problemach z zależnościami, itp. Notatki zapisuj w katalogu `memory/issues/`.

Każda notatka powinna zawierać:

1. **Opis problemu** — co się dzieje, jaki błąd, w jakim kontekście
2. **Diagnostyka** — co sprawdzono, jakie hipotezy testowano, co wykluczone
3. **Status** — aktywny / rozwiązany
4. **Rozwiązanie** (gdy znalezione) — co było przyczyną i jak naprawiono

### Workflow
- Gdy napotkasz nową trudność → utwórz plik w `memory/issues/nazwa-problemu.md`
- Aktualizuj notatkę w miarę postępów diagnostyki (dodawaj nowe obserwacje, wyniki testów)
- Gdy problem zostanie rozwiązany → przenieś plik do `memory/resolved/`
- W MEMORY.md utrzymuj krótką listę aktywnych issues z linkami

### Format notatki

Nazwa pliku: krótka-nazwa-problemu.md (np. `cholvecs-write-bug.md`, `mkl-link-error.md`)

Zawartość:
- **Status:** aktywny | rozwiązany
- **Data:** kiedy napotkano
- **Pliki:** lista dotkniętych plików źródłowych
- **Objaw:** krótki opis błędu (komunikat, backtrace)
- **Przyczyna:** (gdy znana) root cause
- **Rozwiązanie:** (gdy znalezione) co zmieniono i dlaczego to działa
- **Diagnostyka:** chronologiczna lista kroków diagnostycznych z wynikami

## TODO

- Rozważyć dodanie `results/` do `.gitignore` — katalog zawiera generowane dane profilowania (perf.data, flamegraph.svg itp.) i nie powinien być śledzony w repo na dłuższą metę
