# Input

## Contents

- [Running `gammcor`](#running-gammcor)
  - [Meson release build](#1-meson-release-build)
  - [Meson debug build](#2-meson-debug-build)
  - [Make build](#3-make-build)
- [General Input Structure](#general-input-structure)
- [Calculation block](#calculation-block)
  - [Interface and required files](#interface)
  - [JobType](#jobtype)
  - [RDMType](#rdmtype)
  - [TwoMoInt](#twomoint)
  - [Other keywords](#other-keywords)
- [CholeskyBlock](#choleskyblock)
- [System block](#system-block)
- [Basis set](#basis-set)
  - [Different basis sets for different atoms](#different-basis-sets-for-different-atoms)
- [Geometry](#geometry)

## Running `gammcor`

GAMMCOR always reads the file `input.inp` from the current directory. Go to
the directory containing the input and the files from the reference
calculation, then run the appropriate executable.

### 1. Meson release build

```bash
cd my_job
/path/to/gammcor/bin/gammcor
```

### 2. Meson debug build

```bash
cd my_job
/path/to/gammcor/bin/gammcor-debug
```

### 3. Make build

```bash
cd my_job
/path/to/gammcor/gammcor
```

<!-- TODO: threads. GAMMCOR uses OpenMP and threaded MKL. Confirm and describe:
export OMP_NUM_THREADS=16  (= number of physical cores reserved for the job) -->

## General Input Structure

The input consists of blocks. Each block starts with its name and ends with
`end`. Keywords and block names are case-insensitive. Lines starting with `!`
are comments.

Below is a mockup illustrating the structure:

```text
Calculation
 Interface  {program}
 JobType    {method}
 RDMType    {reference}
 TwoMoInt   {integral_mode}
 BasisPath  {basis_directory}
 Basis      {basis_file}
end

CholeskyBlock
 Cholesky          {cholesky_mode}
 Cholesky_Accuracy {accuracy}
end

System
 ZNucl     {sum_of_nuclear_charges}
 Charge    {molecular_charge}
 ThrSelAct {threshold}
end

xyz
 {number_of_atoms}
 {atom_symbol} {x} {y} {z}
 ...
end
```

Which blocks are needed depends on the interface and on how the two-electron
integrals are supplied. With Molpro and Dalton, the geometry and basis set come
from the reference files. With PySCF, `BasisPath`, `Basis`, and `xyz` are needed
when GAMMCOR constructs the integrals on the fly, including THC. They are not
needed when PySCF has already written the integrals and the GAMMCOR input uses
`TwoMoInt INCORE`.

<!-- TODO: confirm which blocks are required for each interface. -->

## Calculation block

### Interface

Selects the program that computed the reference wave function.

| Value | Program |
|---|---|
| `PYSCF` | PySCF, via `pyscf_data.h5` written by `gammcor_tools` (see [PySCF interface](03_pyscf_interface.md)) |
| `MOLPRO` | Molpro |
| `DALTON` | Dalton |
| `ORCA` | ORCA |

Put the interface files in the same directory as `input.inp`. File names are
fixed and case-sensitive.

#### PySCF files

- `pyscf_data.h5` contains the molecular-orbital data, energies, density
  matrices, and, when requested, two-electron integrals.

#### Molpro files

The Molpro examples in `testjobs/ACn/AR_CAS` and `testjobs/ACn/AR_HF` use:

- `AOONEINT.mol` -- one-electron AO integrals and molecular information;
- `AOTWOINT.mol` -- two-electron AO integrals or Cholesky vectors;
- `MOLPRO.MOPUN` -- molecular-orbital coefficients;
- `2RDM` -- the two-particle reduced density matrix;
- `DIP` -- dipole integrals, needed only for calculations that use them.

#### Dalton files

The Dalton AC0 example in `testjobs/AC0/ethene` uses:

- `SIRIUS.RST` -- orbital and symmetry information;
- `SIRIFC` -- the Dalton interface data;
- `AOONEINT` -- one-electron AO integrals;
- `AOTWOINT` -- two-electron AO integrals;
- `rdm2.dat` -- the two-particle reduced density matrix.

#### ORCA files

The ORCA examples in `testjobs/ACn/H2O_MOLMPS` and
`testjobs/RDMCORR/H2O_MOLMPS` use:

- `G1.bin` and `G2.bin` -- the one- and two-particle reduced density matrices;
- `FACT.bin` -- one-electron integral data;
- `NUC_REP.bin` -- the nuclear-repulsion energy;
- `C0.bin` -- molecular-orbital coefficients, used when the integrals are
  constructed on the fly;
- `IPQRS.bin` -- additional ORCA interface data included with the test case.

The supplied ORCA inputs use on-the-fly integrals and therefore need `C0.bin`.
The test directories keep the complete exported file set, including
`IPQRS.bin`.

### JobType

Selects the correlation method.

| Value | Description |
|---|---|
| `AC0` | linearized adiabatic connection <!-- TODO: ref --> |
| `AC` | adiabatic connection <!-- TODO: ref --> |
| `AC0PP` | linearized particle-particle adiabatic connection (ppAC0) |
| `ACPP` | particle-particle adiabatic connection (ppAC) |
| `ACFREQNTH` | ACn: adiabatic connection expanded up to order `Max_Cn` <!-- TODO: check wording --> |
| `ERPA` | extended random-phase approximation |
| `SAPT` | symmetry-adapted perturbation theory, for interaction energies of dimers |
| `MP2` | second-order Møller–Plesset perturbation theory |
| `TDDFT` | time-dependent DFT |
| `PDFT`, `CASPIDFT` | <!-- TODO --> |

<!-- TODO: describe or drop the remaining JobType values found in SOURCE/inputfill.f90:
AC1, EERPA (ERPA-2), EERPA-OLD (ERPA-1), CASPIDFTOPT, AC0D, AC0DNOSYMM, AC0DP,
ACFREQ, AC1FREQNTH, RESPONSE, NLOCCORR, SRAC0, SRMP2, SAPT-OS, MSAC0, PPERPA,
PPERPADMP, HHERPADMP, DUCC, AC_REL, AC0_REL.
Mark the ones meant for testing only, as in beyond-rpa's "Additional methods" table. -->

### RDMType

Type of the reference wave function (its one- and two-particle reduced density matrices).

| Value | Reference |
|---|---|
| `CAS` (or `CASSCF`) | complete active space SCF |
| `GVB` | generalized valence bond |
| `APSG` | antisymmetrized product of strongly orthogonal geminals |
| `DMRG` | density matrix renormalization group |
| `HF` | Hartree–Fock |
| `UHF` | unrestricted Hartree–Fock |
| `RKS`, `UKS` | restricted / unrestricted Kohn–Sham |

### TwoMoInt

How the two-electron integrals in the molecular-orbital basis are stored.

| Value | Description |
|---|---|
| `FOFO` | <!-- TODO --> |
| `FULL` (or `FFFF`) | <!-- TODO --> |
| `INCORE` | <!-- TODO --> |

### Other keywords

| Keyword | Description |
|---|---|
| `JobTitle` | free text, printed in the output |
| `Units` | units of the `xyz` block: `angstrom` (default) or `bohr` |
| `BasisPath` | directory with the basis-set files |
| `Basis` | basis-set file inside `BasisPath`, e.g. `cc-pvdz.txt` |
| `Max_Cn` | highest order for `ACFREQNTH` (default 3) |

## CholeskyBlock

Controls the Cholesky decomposition of the two-electron integrals.

| Keyword | Values |
|---|---|
| `Cholesky` | `BIN` (from binary files), `OTF` (computed on the fly), `THC` (tensor hypercontraction) |
| `Cholesky_Accuracy` | `DEFAULT`, `TIGHT`, `LUDICROUS` (short: `D`, `T`, `L`) |

<!-- TODO: when to use BIN vs OTF; recommended accuracy. The PySCF tests use
OTF with LUDICROUS. -->

## System block

| Keyword | Description |
|---|---|
| `ZNucl` | sum of nuclear charges of all atoms (e.g. 14 for N₂) |
| `Charge` | total charge of the molecule |
| `ThrSelAct` | <!-- TODO --> |
| `ThrQVirt` | <!-- TODO --> |
| `ThrQInact` | <!-- TODO --> |
| `NActive` | <!-- TODO --> |
| `NCore` | <!-- TODO --> |
| `Monomer` | `A`, `B` or `AB`, for SAPT <!-- TODO --> |
| `ThrPP` | threshold for selecting orbital pairs in `AC0PP` and `ACPP` calculations (default `1.d-8`) |

## Basis set

The classic way: one basis set for all atoms, given as a file inside `BasisPath`.

```text
Calculation
 BasisPath /path/to/gammcor/basis
 Basis     cc-pVDZ
end
```

The basis-set files are in the GAMESS (US) format, which can be downloaded from
the [Basis Set Exchange](https://www.basissetexchange.org). Ready-made files:

| Directory | Files |
|---|---|
| `basis/` | `cc-pVDZ`, `cc-pVTZ`, `cc-pVQZ`, `aug-cc-pVDZ`, `aug-cc-pVQZ`, and `-OPT` / `_DALTON` variants |
| `basis-for-pyscf/` | pairs `name.txt` (for GAMMCOR) + `name.nw` (for PySCF), e.g. `cc-pvdz.txt` + `cc-pvdz.nw` |

With the PySCF interface, both files of a pair must describe the same basis:
PySCF reads the `.nw` file, GAMMCOR the `.txt` file.

<!-- TODO: what the -OPT and _DALTON variants are. -->

### Different basis sets for different atoms

The `BasisAssignement` block assigns different basis sets to elements, single
atoms, or ranges of atoms. The standard input with one basis set for all atoms
remains supported.

`BasisPath` is still required in the `Calculation` block. If
`BasisAssignement` is present, the `Basis` keyword is ignored:

```text
Calculation
 BasisPath /path/to/gammcor/basis-for-pyscf
 Basis     cc-pvdz.txt
end

BasisAssignement
 * cc-pvdz.txt
end
```

Each line in `BasisAssignement` has the form:

```text
<selector>  <basis-file>
```

`<basis-file>` must be the name of a file inside `BasisPath`, including its
extension when the file has one.

#### Selectors

```text
BasisAssignement
 *    cc-pvdz.txt
 Fe   cc-pvqz.txt
 10   cc-pvtz.txt
 6-8  cc-pvtz.txt
end
```

| Selector | Applies to |
|---|---|
| `*` | every atom (fallback) |
| element symbol, e.g. `Fe` | all atoms of that element |
| atom number, e.g. `10` | one atom, counted from 1 in the `xyz` block |
| atom range, e.g. `6-8` | consecutive atoms in the `xyz` block |

#### Priority

When several rules match an atom, the most specific rule wins:

1. **Atom number or range** (highest priority)
2. **Element symbol**
3. **Fallback `*`** (lowest priority)

The order of the lines does not matter.

#### Examples

One common basis set for all atoms:

```text
BasisAssignement
 * cc-pvdz.txt
end
```

A default basis set with different basis sets for selected elements:

```text
BasisAssignement
 *  cc-pvdz.txt
 Fe cc-pvqz.txt
 S  cc-pvtz.txt
end
```

A basis set assigned to one specific atom:

```text
BasisAssignement
 *  cc-pvdz.txt
 Fe cc-pvqz.txt
 10 cc-pvtz.txt
end
```

Here, `10` refers to the tenth atom in the `xyz` block. The atom-specific rule
takes priority over both the element rule and the default rule.

Use a range when consecutive atoms in the `xyz` block should use the same
basis set:

```text
BasisAssignement
 *    cc-pvdz.txt
 1    cc-pvqz.txt
 6-10 cc-pvtz.txt
end
```

In this example, atom 1 uses `cc-pvqz.txt`, atoms 6 through 10 use
`cc-pvtz.txt`, and all remaining atoms use `cc-pvdz.txt`.

## Geometry

The `xyz` block gives the number of atoms, then one line per atom.
Coordinates are in Ångström unless `Units bohr` is set in `Calculation`.

```text
xyz
2
N    0.000000000    0.000000000    0.000000000
N    0.000000000    0.000000000    1.098000000
end
```

The geometry can also be read from a standard `.xyz` file (as written by
Avogadro, VMD, Molden, OpenBabel, ORCA, …). Put the file name on the `xyz` line:

```text
xyz geom.xyz
end
```

The path may be relative or absolute.

<!-- TODO: describe how to produce the Molpro, Dalton, and ORCA interface files. -->

---

<div align="center">

| Previous | | Next |
|:--|:--:|--:|
| [← Setup](01_setup.md) | [↑ Back to README](../README.md) | [PySCF interface →](03_pyscf_interface.md) |

</div>
