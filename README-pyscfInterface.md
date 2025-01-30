# PYSCF - GAMMCOR Interface Description

PYSCF works in INCORE or simple THC mode. For THC, you need to set Cholesky OTF in input.

## 1. Setup and Example

**1.1Example Input (eth.py)**

```python
import numpy as np
from pyscf import gto, scf, mcscf, ao2mo, symm, tools, fci, mrpt
import sys
from copy import deepcopy
sys.path.append('/home/aleksandra.tucholska/test/gammcor/')
from Pyscf2Gammcor import get_data_for_gammcor
from Pyscf2Gammcor import irrep_analyze, analysis_of_mo

if __name__ == "__main__":

    basis_file = '/home/aleksandra.tucholska/test/gammcor/bazy-do-pyscf/def2-tzvp.nw'

    mol = gto.Mole()
    mol.basis = basis_file
    mol.unit = 'angs'
    mol.charge = 0
    mol.spin = 0
    mol.cart = False
    mol.symmetry = True
    mol.build(atom='''H     0.0000000     0.9232740     1.2382890                                                             
H     0.0000000     -0.9232740     1.2382890                                                       
H     0.0000000     0.9232740     -1.2382890                                                       
H     0.0000000     -0.9232740     -1.2382890                                                      
C     0.0000000     0.0000000     0.6681880                                                                                                   
C     0.0000000     0.0000000     -0.6681880                                                                                                                                                                                                                                                                                  
''')

    myhf = scf.RHF(mol)
    myhf.kernel()

    nelec = mol.nelectron
    print(f"Total number of electrons: {nelec}")
    irrep_analyze(mol, myhf)

    ncas = {'B3g': 1, 'B2u': 1}
    ncore ={'Ag': 3, 'B3u': 1, 'B1u': 2, 'B3g': 1}


    mycas = mcscf.CASSCF(myhf, 2, 2)
    mcscf.addons.state_specific_(mycas, state=0)


    mycas.conv_tol = 1e-9
    mycas.conv_tol_grad = 1e-5
    print('mycas.conv_tol', mycas.conv_tol)
    print('mycas.conv_tol_grad', mycas.conv_tol_grad)

    mycas.natorb = True
    mo = mcscf.sort_mo_by_irrep(mycas, myhf.mo_coeff, ncas, ncore)
    mycas.fix_spin_(ss=0.0)
    mycas.fcisolver.wfnsym = 'Ag'
    mycas.frozen = 0
    mycas.canonicalization = True # this is true by default                                                                                                                                                                                                                                                                   

    mycas.kernel(mo)

    eorb = mycas.mo_energy

    get_data_for_gammcor(mol, myhf, mycas, dump_eri=True)

```

**1.2 Required Setup**

Before running any calculations:

1. Replace the following line with your path to gammcor where the Pyscf2Gammcor.py file is:
```python
sys.path.append('/home/aleksandra.tucholska/test/gammcor/')
sys.path.append('/path/to/your/gammcor/')
```

2. Replace the following line with  the correct path to basis set files:
```python
basis_file = '/home/aleksandra.tucholska/test/gammcor/bazy-do-pyscf/def2-tzvp.nw'
basis_file = '/path/to/your/gammcor/basis-for-pyscf/basis.nw'
```

### 2 Running Calculations

**2.1 Basic CASSCF Setup:**

- Always use `natorb=True` for state-specific calculations
- Canonicalization is on by default (recommended)
- For state-averaged calculations naturalization is done in Gammcor

**2.2 Dumping Files for GAMMCOR:**

```python
# After CASSCF calculation, add:
get_data_for_gammcor(mol, myhf, mycas, dump_eri=True)
```
   - Set `dump_eri=True` if you want to dump two-electron integrals
   - Use `dump_eri=False` or simply  `get_data_for_gammcor(mol, myhf, mycas)` for direct (non-incore) calculations

## 3 File Reading Logic - this is done automatically. Below is simply the description.

**3.1 File Naming Patterns**

1. State-specific calculations use simple names (e.g., `auxdata.txt`, `rdm2_aaaa.bin`)
2. State-averaged calculations use numbered suffixes in Molpro style (e.g., `auxdata_1.1.txt`, `rdm2_aaaa_1.1.bin`)

**3.2 Binary Files**

- `CAONO.bin`: Orbital coefficients matrix (float64)
  
  - When `natural = 0`: Contains MO coefficients
  - When `natural = 1`: Contains NO coefficients
  
- `HCore.bin`: One-electron integrals matrix (float64)
  
  - In AO basis
  
- `rdm2_aaaa.bin`: Alpha-Alpha reduced density matrix (active space)

- `rdm2_abab.bin`: Alpha-Beta reduced density matrix (active space)

- `rdm2_full.bin`: Full 2-RDM (2.0 * (dm_aaaa + dm_abab))
  
- all rdm.bin files
  
  - When `natural = 0`: In MO basis
  - When `natural = 1`: In NO basis
  
- `TWOEl.bin`: Two-electron integrals (optional, if dump_eri=True)
  
  - In AO basis
  
    

**3.3 Text Files**

- `auxdata.txt`: Contains basic calculation parameters (one per line):
  1. Number of basis functions (nbasis)
  2. Number of inactive orbitals (NI)
  3. Number of active orbitals (NA)
  4. Number of virtual orbitals (NV)
  5. CASSCF total energy
  6. Nuclear repulsion energy
  7. Number of electrons
  8. Natural orbitals flag (1 if used, 0 if not)

- `rdm2.dat`: Full 2-RDM in text format (always in NO basis)

**3.4 Molden Files**

- `moldenhf.inp`: Hartree-Fock orbitals in AO basis
- `molden_mcscf.inp`: CASSCF orbitals

### Reading Cases 

1. **No State Specified in Gammcor input** (when `nst(1) < 0`):
   - If only `auxdata.txt` exists: Uses it, assumes state 1.1
   - If only `auxdata_x.y.txt` exists: Uses it and set nst to x.y
   - If `auxdata.txt` and `auxdata_x.y.txt` files exist: Uses `auxdata.txt`, assumes state 1.1
   - If multiple `auxdata_x.y.txt` files exist and no `auxdata.txt`: Exits with error.

2. **State Specified in Gammcor input** (when `nst(1) > 0`):
   - Looks for either `auxdata.txt` or `auxdata_x.y.txt` matching the specified state
   - Exits if both or neither exist
   - Uses the specified state numbering for reading other files

3. **State Number Setup (anst)**:

When no state is specified (nst(1) < 0):

   - Single `auxdata.txt`: Sets anst = [1,1]
   - Single state-specific file (e.g., `auxdata_2.3.txt`): Sets anst to those numbers [2,3]
   - Multiple state files: Errors unless `auxdata.txt` exists (then uses [1,1])

When state is specified in input to x.y (nst(1) > 0):

   - Sets anst to the specified state numbers
   - Must find exactly one matching file (either `auxdata.txt` or `auxdata_x.y.txt`)
