# Input keywords

This file explains the input keywords one topic at a time: what each keyword
does, what it changes in the calculation, and when to use it. For the layout of
the input file (blocks, `end`, comments) see
[Input structure](02_input_structure.md).

## Contents

- [Frozen core](#frozen-core)
  - [How to switch it on](#how-to-switch-it-on)
  - [What it does](#what-it-does)
  - [Where it is available](#where-it-is-available)
- [Multistate AC0 (MS-AC0)](#multistate-ac0-ms-ac0)
  - [How to switch it on](#how-to-switch-it-on-1)
  - [What it does](#what-it-does-1)
  - [Where it is available](#where-it-is-available-1)
- [Saving and reusing the THC factorization](#saving-and-reusing-the-thc-factorization)

- [MP2 and MP2-CBS[H]](#mp2-and-mp2-cbsh)
  - [How to switch it on](#how-to-switch-it-on-2)
  - [What it does](#what-it-does-2)
  - [Where it is available](#where-it-is-available-2)

## Frozen core

### How to switch it on

Add `Core 0` to the `Calculation` block:

```text
Calculation
 .
 .
 Core 0
end
```

This freezes **all** inactive orbitals. To freeze only the lowest `n` of
them, also give `NCore n` in the `System` block:

```text
Calculation
 .
 .
 Core 0
end

System
 .
 .
 NCore 2
end
```

| Keyword | Block | Values |
|---|---|---|
| `Core` | `Calculation` | `1` = no frozen core (default); `0` = frozen core |
| `NCore` | `System` | number of frozen orbitals, counted from the lowest inactive orbital. Used only with `Core 0`. If it is omitted or `0`, all inactive orbitals are frozen. `NCoreOrb` is accepted as an alias. |

`NCore` must not be larger than the number of inactive orbitals, which means
the frozen core can never reach into the active space. If it is larger, the
program stops with an error. The output reports the number of frozen orbitals
in the line `Frozen core orbitals (ph pairs)`.

Example for N₂ with a CAS(6,6) reference, freezing the two 1s orbitals:

```text
Calculation
 Interface PYSCF
 JobType   AC0
 RDMType   CAS
 TwoMoInt  FOFO
 Core      0
end

CholeskyBlock
 Cholesky THC
end

System
 ZNucl  14
 Charge 0
 NCore  2
end
```

### What it does

The frozen orbitals stay **doubly occupied in the reference**. They are
excluded only from the correlation treatment.

- **The reference is unchanged.** The frozen electrons still create the Coulomb
  and exchange field that all other electrons feel. The CASSCF energy, the
  Fock matrix and the orbital energies are the same with and without frozen
  core.
- **Excitations out of frozen orbitals are removed.** AC0 obtains the
  correlation energy from the response (ERPA) of the reference to single
  excitations, i.e. orbital pairs `(p,q)` with `p` above `q`. With frozen core,
  every pair whose lower orbital is frozen is dropped:
  - inactive → virtual pairs with a frozen inactive orbital;
  - inactive → active pairs with a frozen inactive orbital.
- **What the correlation energy loses.** The energy loses exactly the
  contributions that involve those pairs, which is the core correlation of the
  frozen electrons. Contributions that involve only active and virtual orbitals
  (active → virtual, active → active) do not change.

**When to use it.** Core correlation is large in absolute terms, but nearly
constant along a reaction path or between electronic states, so it mostly
cancels in energy differences. Standard basis sets such as cc-pVXZ are not
designed to describe it. Freezing the core therefore usually costs little
accuracy for valence properties, and it makes the calculation smaller and
faster, because fewer orbital pairs enter the response problem. For core
properties, or with core-valence basis sets (cc-pCVXZ), keep `Core 1`.

### Where it is available

| Method | Frozen core |
|---|---|
| `JobType AC0` with THC integrals (`TwoMoInt FOFO`, `Cholesky THC`) | `Core 0`, all inactive orbitals or the lowest `NCore` |
| `JobType AC0`, other integral modes | `Core 0` freezes all inactive orbitals; `NCore` is not used |
| `JobType AC0PP`, `ACPP` (particle-particle AC) | **not yet implemented**: `Core` and `NCore` are ignored |
| `JobType MP2` | always on; only `NCore` counts, and `NCore 0` freezes nothing (see [MP2 and MP2-CBS[H]](#mp2-and-mp2-cbsh)) |

<!-- TODO: confirm which other AC0 routes honour Core 0 (the spin-resolved
acph_spinres_driver does; the classic AC0CAS path does not appear to). Remove
this note once frozen core in the ph THC path has been tested against FOFO. -->

## Multistate AC0 (MS-AC0)

### How to switch it on

Set `JobType MSAC0` and give **one** state in `State`: the reference state.
`State n.m` is state `n` of symmetry `m`. The states it couples to are found
automatically. One run gives the couplings of one reference state, so **run
once per state** (`State 1.1`, then `State 2.1`, ...).

| Route | `Interface` | `TwoMoInt` | `CholeskyBlock` | Coupled states |
|---|---|---|---|---|
| Molpro | `MOLPRO` | `INCORE` (default) or `FOFO` | - | all states in the Molpro `2RDM` file |
| PySCF | `PYSCF` | `FOFO` | `Cholesky THC` | all states of symmetry `m` in `pyscf_data.h5` |

**Molpro.** Run a state-averaged CASSCF over the states and write the `2RDM`
file for GammCor. Input, here for state 2 (add `TwoMoInt FOFO` for the FOFO
route):

```text
Calculation
 Interface MOLPRO
 JobType   MSAC0
 RDMType   CAS
end

System
 ZNucl  12
 Charge 0
 State  2.1
end
```

**PySCF.** Run a state-averaged CASSCF with `mycas.natorb = False` (and
`mycas.fix_spin_(ss=0)` for singlets) and call `get_data_for_gammcor`. Input:

```text
Calculation
 Interface PYSCF
 JobType   MSAC0
 RDMType   CAS
 TwoMoInt  FOFO
end

System
 ZNucl  12
 Charge 0
 State  1.1
end

CholeskyBlock
 Cholesky THC
end
```

A complete example (LiF, two states, a bond-length scan) is in
`testjobs/pyscf-tests/test-MSAC0/LiF`.

### What it does

MS-AC0 builds a small effective Hamiltonian `H` over the CASSCF states and
diagonalizes it, so that states close in energy can mix after the AC0
correlation is added.

The run with reference `I` prints its elements of `H` as

```text
 H_mn  I  I  value       the AC0 energy of state I, without nuclear repulsion
 H_mn  J  I  value       the coupling of state J to state I
```

`H(I,J)` and `H(J,I)` come from different runs and are not equal. The tool
`gammcor_heff` (from `gammcor_tools`) collects them from the outputs of all
runs of one geometry, prints `H` and the eigenvalues of `(H + H^T)/2`, the
MS-AC0 energies (without nuclear repulsion):

```text
gammcor_heff gammcor_1.out gammcor_2.out
```

### Where it is available

| Route | Limits |
|---|---|
| Molpro | at most 10 states; all must have the **same symmetry** (not checked) |
| PySCF, THC | state-averaged CASSCF only; only states of the reference's symmetry and spin |

## Saving and reusing the THC factorization

In the `CholeskyBlock`, next to the other THC settings:

```text
CholeskyBlock
 Cholesky THC
 THCDump  thc.h5        compute Xgp/Zgk, write them, stop
end

CholeskyBlock
 Cholesky THC
 THCLoad  thc.h5        read Xgp/Zgk instead of computing them
end
```

The THC factorization (Cholesky decomposition and the THC vectors `Xgp`,
`Zgk`) depends only on the molecule, its geometry, the basis set and the
Cholesky/THC thresholds, not on the orbitals or the state. `THCDump`
computes it once and stops; later runs with `THCLoad` read it and skip the
factorization, which is usually most of the run time. This is useful for
MS-AC0 (one run per reference state at the same geometry) and for restarting
a calculation.

- The file name is optional; the default is `thc.h5` in the working directory.
- `THCDump` and `THCLoad` together are an input error; both need `Cholesky THC`.
- The file belongs to one geometry: in a scan keep one file per point.

MS-AC0 workflow with it:

```text
input_thc.inp   CholeskyBlock: THCDump thc.h5      (any State; stops after the dump)
input_1.inp     State 1.1, THCLoad thc.h5
input_2.inp     State 2.1, THCLoad thc.h5
gammcor_heff gammcor_1.out gammcor_2.out
```

<!-- TODO: implemented (steps 1-5 of .secret_sauce/.implementation_plans/thc/plan-thc-dump-load.md)
but not yet tested (step 7). -->

## MP2 and MP2-CBS[H]

### How to switch it on

MP2 runs on a Hartree-Fock reference from PySCF. In the PySCF script, run
HF (`scf.RHF`; for `spin > 0` it is ROHF) and write the dump without a CASSCF
object:

```python
from gammcor_tools.Pyscf2Gammcor import get_data_for_gammcor
myhf = scf.RHF(mol).run()
mymp = mp.MP2(myhf).set_frozen().run()      # optional, only for checks
get_data_for_gammcor(mol, myhf, mymp=mymp)
```

Then use `JobType MP2`. The basis-set correction is chosen with `CBS`,
exactly as for AC0:

```text
Calculation
 Interface  PYSCF
 JobType    MP2
 RDMType    CAS
 TwoMoInt   FOFO
 Basis      aug-cc-pvdz.txt
 CBS        2
end

CholeskyBlock
 Cholesky          THC
 Cholesky_Accuracy L
end

System
 ZNucl    14
 Charge   0
 NCore    2
end
```

| Keyword | Block | Values |
|---|---|---|
| `JobType MP2` | `Calculation` | MP2 on the HF dump |
| `CBS` | `Calculation` | absent or `0`: plain MP2; `1` (or `DBBSC`): MP2 + CBS[DFT]; `2` (or `H`, `CBS[H]`): MP2-CBS[H] |
| `NCore` (`NCoreOrb`) | `System` | number of frozen core orbitals; `0` or absent: none frozen |
| `Grid` | `Calculation` | grid for mu(r) with `CBS 1` or `2`; default `medium` |

There is no separate `SRMP2` job type: SR-MP2 is `JobType MP2` with `CBS 2`.
`TwoMoInt FOFO` with `Cholesky THC` is required.

### What it does

- **MP2** (`CBS 0`). Closed-shell MP2 with THC integrals. For an open-shell
  (ROHF, high-spin) dump it runs ROMP2: semicanonical alpha and beta
  orbitals, same-spin and opposite-spin doubles, and singles. ROMP2 is not
  what PySCF prints as `mp.MP2(rohf)` (that is UMP2 on the unrotated ROHF
  orbitals).
- **MP2 + CBS[DFT]** (`CBS 1`). MP2 plus the density-based basis-set
  correction of Giner et al., computed with a local range-separation
  parameter mu(r) from the HF pair density. The integrals are not changed.
  For ROHF the PBE correlation is evaluated for the total density
  (zeta = 0), as in the AC0 FOFO code; the spin enters through the on-top
  pair density.
- **MP2-CBS[H]** (`CBS 2`). The two-electron integrals in the MP2 numerator
  are dressed with the short-range interaction weighted by mu(r); the orbital
  energies are not changed. The output also prints plain MP2 and the
  difference. For ROHF the doubles are dressed and the singles are not, as in
  AC0-CBS[H].

**Frozen core.** `NCore n` freezes the lowest `n` orbitals in the
correlation energy and in mu(r). The frozen orbitals stay in the reference
(Fock operator, HF energy). Unlike AC0, `Core 0` is not needed, and
`NCore 0` means no frozen orbitals. `NCore` must be the core of this
molecule (N atom 1, N2 2, SO2 7); a value copied from another molecule can
freeze valence orbitals.

The last output line:

| Job | Line |
|---|---|
| MP2, closed shell | `HF+ENuc, MP2-Corr     , Total` |
| MP2, open shell | `HF+ENuc, ROMP2-Corr   , Total` |
| `CBS 1`, closed shell | `HF+ENuc, MP2-CBS[DFT], Total` |
| `CBS 1`, open shell | `HF+ENuc, ROMP2-CBS[DFT], Total` |
| `CBS 2`, closed shell | `HF+ENuc, MP2-CBS[H]-Corr, Total` |
| `CBS 2`, open shell | `HF+ENuc, ROMP2-CBS[H]-Corr, Total` |

### Where it is available

| Reference | `CBS 0` | `CBS 1` | `CBS 2` |
|---|---|---|---|
| RHF (closed shell) | yes | yes | yes |
| ROHF (high spin) | yes (ROMP2) | yes (ROMP2 + CBS[DFT]) | yes (ROMP2-CBS[H]) |

Only the PySCF interface with THC integrals. A CASSCF dump cannot be used for
MP2, and the HF dump cannot be used for AC0.

<!-- TODO: CBS 1 (closed and open shell), ROMP2-CBS[H] and the removal of
JobType SRMP2 are written but not yet compiled or tested
(.secret_sauce/.skills/mp2/code-mp2.md Sec. 7). -->
