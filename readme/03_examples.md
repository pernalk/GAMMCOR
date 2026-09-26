# Examples

All examples are taken from `testjobs/`, where you also find the files of the
reference calculations.

### 1. AC0 correlation energy of N₂ from a PySCF CASSCF(6,6) reference

**Step 1.** Run the reference calculation in PySCF (`n2_pyscf.py`). The last
line writes `pyscf_data.h5`:

```python
from pyscf import gto, scf, mcscf
from gammcor_tools.Pyscf2Gammcor import get_data_for_gammcor

mol = gto.Mole()
mol.basis = '/path/to/gammcor/basis-for-pyscf/cc-pvdz.nw'
mol.unit = 'bohr'
mol.symmetry = True
mol.symmetry_subgroup = 'D2h'
mol.build(atom='''N    0.000000000    0.000000000    0.000000000
N    0.000000000    0.000000000    2.075000000
''')

myhf = scf.RHF(mol)
myhf.kernel()

ncas  = {'Ag': 1, 'B2g': 1, 'B3g': 1, 'B1u': 1, 'B2u': 1, 'B3u': 1}
ncore = {'Ag': 2, 'B1u': 2}

mycas = mcscf.CASSCF(myhf, 6, 6)
mcscf.addons.state_specific_(mycas, state=0)
mycas.conv_tol = 1e-9
mycas.conv_tol_grad = 1e-5
mycas.natorb = True
mo = mcscf.sort_mo_by_irrep(mycas, myhf.mo_coeff, ncas, ncore)
mycas.fix_spin_(ss=0.0)
mycas.fcisolver.wfnsym = 'Ag'
mycas.kernel(mo)

get_data_for_gammcor(mol, myhf, mycas)
```

```bash
python3 n2_pyscf.py
```

**Step 2.** Write `input.inp` in the same directory. The basis set and geometry
must match the PySCF script:

```text
Calculation
 JobTitle   AC0
 Interface  PYSCF
 JobType    AC0
 RDMType    CAS
 TwoMoInt   FOFO
 BasisPath  /path/to/gammcor/basis-for-pyscf
 Basis      cc-pvdz.txt
 Units      bohr
end

CholeskyBlock
 Cholesky          OTF
 Cholesky_Accuracy L
end

xyz
2
N    0.000000000    0.000000000    0.000000000
N    0.000000000    0.000000000    2.075000000
end

System
 ZNucl      14
 Charge     0
 ThrSelAct  1.d-4
 ThrQInact  1.d-2
end
```

**Step 3.** Run GAMMCOR:

```bash
gammcor > gammcor.out
```

<!-- TODO: show the lines of gammcor.out with the final energies and explain them. -->

### 2. ACn correlation energy of the Ar atom from a Molpro CASSCF reference

The geometry and the basis set are read from the Molpro files, so `input.inp`
has no `xyz` block. Input from `testjobs/ACn/AR_HF`:

```text
Calculation
 JobTitle   ACn
 Interface  MOLPRO
 JobType    ACFREQNTH
 Max_Cn     10
 RDMType    CAS
 TwoMoInt   FOFO
end

CholeskyBlock
 Cholesky          BIN
 Cholesky_Accuracy DEFAULT
end

System
 ZNucl      18
 Charge     0
 ThrSelAct  1.d-4
 ThrQVirt   1.d-1
end
```

<!-- TODO: how to produce the Molpro files (the Molpro input ar.inp is in the test directory). -->

### 3. AC0 energy of ethene from a Dalton CASSCF reference

Input from `testjobs/AC0/ethene` (the Dalton run is in `run_cas_dalton.py`):

```text
Calculation
 JobTitle   AC0 dalton job
 Interface  DALTON
 JobType    AC0
 RDMType    CAS
 TwoMoInt   FOFO
end

System
 ZNucl      16
 Charge     0
 ThrSelAct  1.d-4
end
```

<!-- TODO: more examples:
 - SAPT interaction energy of a dimer
 - DMRG reference
 - mixed basis sets with BasisAssignement
 - ERPA excitation energies -->
