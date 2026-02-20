# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

GAMMCOR is a Fortran-based computational chemistry program for advanced electronic structure calculations including SAPT (Symmetry-Adapted Perturbation Theory), AC/AC0/AC1 (Adiabatic Connection), ERPA, EERPA, and CASPiDFT methods. The internal program name is PRDMFT (Pair and Reduced Density Matrix Functional Theory).

## Build Commands

**Prerequisites:** Intel MKL, and either Intel Fortran (`ifort`) or GNU Fortran (`gfortran`). Set `MKL_ROOT` in the appropriate Makefile before building.

```bash
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

## Running Tests

```bash
python gammcor_verify.py
```

This requires a `TESTS/` directory (gitignored) containing test cases with reference `gammcor.out` files. The script runs gammcor in each test directory and compares energy values against references. Tolerance: 1e-7 Ha for AC/ERPA methods, 1e-5 mHa for SAPT components.

## Architecture

### Source Organization

All source code is in `SOURCE/`. The codebase mixes fixed-format Fortran 77 (`.f`) with free-format Fortran 90 (`.f90`). Global state is shared via common blocks defined in `SOURCE/commons.inc`.

### Module Dependency Chain

```
types.f90  (data types: InputData, FlagsData, SystemBlock, SaptData)
  ├── inputfill.f90  (parses input.inp files)
  ├── systemdef.f90  (initializes molecular system)
  ├── sorter.f90, tran.f90  (integral sorting/transformation)
  ├── abmats.f90, abfofo.f90  (AB matrix construction, depend on tran)
  ├── diis.f90  (DIIS convergence acceleration)
  ├── sapt_utils.f90 → sapt_pol.f90, sapt_exch.f90 → sapt_main.f90
  └── exmisc.f90 → exappr.f90, exdpino.f90, srefex.f90
```

### Program Flow (mainp.f)

1. `read_Input` → `check_Calc` → `fill_Flags` → `create_System` (input parsing)
2. Routes to SAPT (`sapt_driver`) or loads integrals (`LdInteg`/`ReadDAL`)
3. Dispatches to calculation method: DMSCF, CASPiDFT, CASPiDFTOPT, or VV10

### Key Subsystems

- **SAPT** (`sapt_main.f90`, `sapt_exch.f90`, `sapt_pol.f90`, `sapt_utils.f90`): Intermolecular interaction energy decomposition. Supports SAPT0 and SAPT2 levels.
- **AC methods** (`acfd.f`, `accas.f`, `ac_fofo.f90`, `ac_iter.f`): Adiabatic connection correlation energy. AC0 is linearized (MP2-like), AC uses full numerical integration.
- **ERPA** (`erpa.f`, `interpa.f`): Extended RPA response equations.
- **Integral handling** (`abmats.f90`, `abfofo.f90`, `tran.f90`, `sorter.f90`): Two-electron integrals stored as (FF|OO) and (FO|FO) in chemists' notation. In-core vs disk (FOFO) variants.
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
- `NBasis` = number of basis functions; `NELE` = number of electrons
- Integer flags use Fortran convention: 0=false/off, 1=true/on
- Constants are in `SOURCE/constants.h`; common blocks in `SOURCE/commons.inc`
- Build produces `.o` files in `OBJ/` and `.mod` files in the project root

## TODO

- Rozważyć dodanie `results/` do `.gitignore` — katalog zawiera generowane dane profilowania (perf.data, flamegraph.svg itp.) i nie powinien być śledzony w repo na dłuższą metę
