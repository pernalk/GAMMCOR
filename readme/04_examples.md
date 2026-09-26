# Examples

Run each calculation in its own directory. The directory must contain
`input.inp` and the reference data required by the selected interface.

## 1. Exporting a reference from PySCF

There is no need to copy a complete PySCF calculation into the GAMMCOR input.
Import the exporter in your PySCF script and call it after the SCF and CASSCF
calculations have converged:

```python
from gammcor_tools.Pyscf2Gammcor import get_data_for_gammcor

# mol, myhf, and mycas are the converged PySCF objects.
get_data_for_gammcor(mol, myhf, mycas)
```

This writes `pyscf_data.h5` in the current directory. All arguments and export
options are described in [PySCF interface](03_pyscf_interface.md).

The important choice is how GAMMCOR will obtain the two-electron integrals:

- **Integrals stored by PySCF:** call the exporter with `dump_eri=True` and use
  `TwoMoInt INCORE`. GAMMCOR then reads the integrals from `pyscf_data.h5`; no
  GAMMCOR basis file or `xyz` block is needed. This is the non-THC route in the
  ppAC0 and ppAC examples below.
- **Integrals constructed by GAMMCOR:** leave `dump_eri=False`, provide the
  geometry and a GAMMCOR basis file, and select the required Cholesky mode.
  The THC example below uses this route.

Use Python booleans, for example `dump_eri=True`, not `dump_eri='True'`.

## 2. ppAC0 and ppAC without THC

For a calculation without THC, export the two-electron integrals from PySCF:

```python
from gammcor_tools.Pyscf2Gammcor import get_data_for_gammcor

get_data_for_gammcor(mol, myhf, mycas, dump_eri=True)
```

Place the resulting `pyscf_data.h5` next to `input.inp`. A ppAC0 input is:

```text
Calculation
 JobTitle  ppAC0
 Interface PYSCF
 JobType   AC0PP
 RDMType   CAS
 TwoMoInt  INCORE
 Algorithm ppsimple
end

System
 ZNucl     7
 Charge    0
 ThrSelAct 1.d-4
 ThrQVirt  1.d-4
 ThrPP     1.d-12
end
```

There is intentionally no `BasisPath`, `Basis`, `CholeskyBlock`, or `xyz`
block: the required integrals have already been stored in `pyscf_data.h5`.
You can confirm that the integral dataset is present with:

```bash
gammcor_h5info pyscf_data.h5 MOINTS/ERI
```

For the full ppAC calculation, use the same exported data and change the two
job lines to:

```text
 JobTitle  ppAC
 JobType   ACPP
```

These settings follow the `N-ppac0`, `N2-ppac0`, and `N-ppac` tests in
`/home/aleksandra.tucholska/gamtest/pyscf-tests/test-ppac0` on Dragon.

## 3. ppAC0 with THC

With THC, GAMMCOR constructs the integrals. The calculation therefore needs a
geometry and a pair of basis files describing exactly the same basis:

- PySCF reads `cc-pvdz.nw`;
- GAMMCOR reads `cc-pvdz.txt`.

For example, the PySCF script can select its half of the pair as follows:

```python
from gammcor_tools.Pyscf2Gammcor import get_data_for_gammcor

mol.basis = "/path/to/gammcor/basis-for-pyscf/cc-pvdz.nw"

# Build mol and run the SCF and CASSCF calculations before exporting.
get_data_for_gammcor(mol, myhf, mycas)
```

The matching `input.inp` is:

```text
Calculation
 JobTitle  ppAC0-THC
 Interface PYSCF
 JobType   AC0PP
 RDMType   CAS
 TwoMoInt  FOFO
 BasisPath /path/to/gammcor/basis-for-pyscf
 Basis     cc-pvdz.txt
 Units     angstrom
end

CholeskyBlock
 Cholesky          THC
 Cholesky_Accuracy LUDICROUS
end

System
 ZNucl     7
 Charge    0
 ThrSelAct 1.d-4
 ThrQVirt  1.d-4
 ThrPP     1.d-12
end

xyz
1
N  0.0  0.0  0.0
end
```

The geometry, units, and basis must match the PySCF calculation. Unlike the
non-THC case, `dump_eri=True` is not required. This setup follows the
`N-ppac0-thc` and `N2-ppac0-thc` tests on Dragon.

## 4. Molpro example from `testjobs`

`testjobs/ACn/AR_CAS` contains an ACn calculation for Ar. Its working
directory contains `AOONEINT.mol`, `AOTWOINT.mol`, `MOLPRO.MOPUN`, `2RDM`, and
`DIP`, together with this GAMMCOR input:

```text
Calculation
 JobTitle  ACn
 Interface MOLPRO
 JobType   ACFREQNTH
 Max_Cn    10
 RDMType   CAS
 TwoMoInt  FOFO
end

CholeskyBlock
 Cholesky          BIN
 Cholesky_Accuracy DEFAULT
end

System
 ZNucl     18
 Charge    0
 ThrSelAct 1.d-2
 ThrQVirt  1.d-1
end
```

`testjobs/ACn/AR_HF` is the corresponding Hartree--Fock reference example.

## 5. Dalton example from `testjobs`

`testjobs/AC0/ethene` contains an AC0 calculation for ethene. The reference
files are `SIRIUS.RST`, `SIRIFC`, `AOONEINT`, `AOTWOINT`, and `rdm2.dat`:

```text
Calculation
 JobTitle  AC0 dalton job
 Interface DALTON
 JobType   AC0
 RDMType   CAS
 TwoMoInt  FOFO
end

System
 ZNucl     16
 Charge    0
 ThrSelAct 1.d-4
end
```

The script `testjobs/AC0/ethene/run_cas_dalton.py` shows the Dalton reference
calculation and collects the files under the names expected by GAMMCOR.

## 6. ORCA example from `testjobs`

`testjobs/ACn/H2O_MOLMPS` contains an ACn calculation for water. Its ORCA data
are in `G1.bin`, `G2.bin`, `FACT.bin`, `NUC_REP.bin`, `C0.bin`, and
`IPQRS.bin`. The supplied input uses on-the-fly integrals and reads the orbital
coefficients from `C0.bin`; `IPQRS.bin` is retained as part of the complete
ORCA test fixture. The test runner supplies `BasisPath`; when running the job
by hand, replace the path below with the location of the GAMMCOR `basis`
directory.

```text
Calculation
 JobTitle  ACn ORCA example
 Interface ORCA
 JobType   ACFREQNTH
 Max_Cn    10
 NBasis    24
 TwoMoInt  FOFO
 BasisPath /path/to/gammcor/basis
 Basis     cc-pVDZ-OPT
end

CholeskyBlock
 Cholesky          OTF
 Cholesky_Accuracy LUDICROUS
end

xyz
3
O  0.000000000  -0.0669721196   0.0057544380
H  0.000000000   0.5668126641   0.7292018977
H  0.000000000   0.4962602719  -0.8205441913
end

System
 ZNucl     10
 Charge    0
 ThrSelAct 1.d-4
end
```

An ORCA AC0 example with relaxed one-particle density matrix options is in
`testjobs/RDMCORR/H2O_MOLMPS`.

## 7. Further PySCF examples in `testjobs`

The repository also contains ready-to-inspect PySCF reference data and
GAMMCOR inputs for:

- `testjobs/pyscf-tests/test-ac0` -- AC0;
- `testjobs/pyscf-tests/test-ac` -- AC;
- `testjobs/pyscf-tests/test-mcpdft` -- MC-PDFT;
- `testjobs/pyscf-tests/test-sr-pbe` -- AC0 with a short-range PBE functional.

Each directory contains `input.inp` and `pyscf_data.h5`, so it also provides a
small real file on which to try `gammcor_h5info`.

## Running an example

Run GAMMCOR from the job directory, using the executable produced by your
build. For example, for a Meson release build:

```bash
cd my_job
/path/to/gammcor/bin/gammcor
```

The Meson debug and Make executables are listed in [Running
`gammcor`](02_input_structure.md#running-gammcor).
