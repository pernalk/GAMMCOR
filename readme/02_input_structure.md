# Input

## Running `gammcor`

GAMMCOR always reads the file `input.inp` from the current directory. Go to
the directory with the input and the files from the reference calculation,
and start the program:

```bash
cd my_job
gammcor > gammcor.out
```

If you built GAMMCOR from source, use the full path, e.g.
`/path/to/gammcor/bin/gammcor > gammcor.out`.

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

Which blocks are needed depends on the interface. With Molpro and Dalton, the
geometry and basis set come from the files of the reference calculation, so
`xyz` and the basis keywords are not needed. With PySCF they are.

<!-- TODO: confirm which blocks are required for each interface. -->

## Calculation block

### Interface

Selects the program that computed the reference wave function.

| Value | Program |
|---|---|
| `PYSCF` | PySCF, via `pyscf_data.h5` written by `gammcor_tools` (see [PySCF interface](#pyscf-interface)) |
| `MOLPRO` | Molpro |
| `DALTON` | Dalton |
| `ORCA` | ORCA |
| `CHRONUSQ` | ChronusQ |
| `OWN`, `NONE` | <!-- TODO --> |

### JobType

Selects the correlation method.

| Value | Description |
|---|---|
| `AC0` | linearized adiabatic connection <!-- TODO: ref --> |
| `AC` | adiabatic connection <!-- TODO: ref --> |
| `ACFREQNTH` | ACn: adiabatic connection expanded up to order `Max_Cn` <!-- TODO: check wording --> |
| `ERPA` | extended random-phase approximation |
| `SAPT` | symmetry-adapted perturbation theory, for interaction energies of dimers |
| `MP2` | second-order Møller–Plesset perturbation theory |
| `TDDFT` | time-dependent DFT |
| `PDFT`, `CASPIDFT` | <!-- TODO --> |

<!-- TODO: describe or drop the remaining JobType values found in SOURCE/inputfill.f90:
AC1, EERPA (ERPA-2), EERPA-OLD (ERPA-1), CASPIDFTOPT, AC0D, AC0DNOSYMM, AC0DP,
ACFREQ, AC1FREQNTH, RESPONSE, NLOCCORR, SRAC0, SRMP2, SAPT-OS, MSAC0, PPERPA,
AC0PP, ACPP, PPERPADMP, HHERPADMP, DUCC, AC_REL, AC0_REL.
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

The `BasisAssignement` block assigns basis sets to elements, single atoms, or
ranges of atoms. `BasisPath` is still needed; `Basis` is ignored when this block
is present.

```text
BasisAssignement
 *    cc-pVDZ
 Fe   cc-pVQZ
 10   cc-pVTZ
 6-8  cc-pVTZ
end
```

| Selector | Applies to |
|---|---|
| `*` | every atom (fallback) |
| element symbol, e.g. `Fe` | all atoms of that element |
| number, e.g. `10` | one atom, counted from 1 in the `xyz` block |
| range, e.g. `6-8` | consecutive atoms |

**Prioritization rules.** When several lines match an atom, the most specific wins:
1. **Atom number or range** (highest priority)
2. **Element symbol**
3. **Fallback `*`** (lowest priority)

The order of the lines does not matter. More examples: [README-BasisAssignement.md](../README-BasisAssignement.md).

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

## PySCF interface

With PySCF, the reference calculation is a Python script. At its end, call
`get_data_for_gammcor`, which writes `pyscf_data.h5` for GAMMCOR:

```python
from gammcor_tools.Pyscf2Gammcor import get_data_for_gammcor

# ... mol, myhf, mycas: your PySCF molecule, SCF and CASSCF objects ...

get_data_for_gammcor(mol, myhf, mycas)
```

The `gammcor_tools` package is installed by `install.sh` and by the conda
package. It also provides command-line tools:

| Command | Purpose |
|---|---|
| `gammcor_basis script.py -o basis.txt` | extract the basis set from a PySCF script, in GAMMCOR format |
| `gammcor_h5info pyscf_data.h5` | show the contents of a GAMMCOR HDF5 file |
| `gammcor_2h5 JOBDIR -o pyscf_data.h5` | convert legacy `.bin`/`.txt` input files to HDF5 |
| `gammcor_2bin pyscf_data.h5 -o JOBDIR` | convert HDF5 back to the legacy files |

Run any command with `-h` to see its options. A complete example is in [Examples](03_examples.md).

<!-- TODO: sections for the Molpro, Dalton and ORCA interfaces: which files
GAMMCOR needs from each program and how to produce them. -->
