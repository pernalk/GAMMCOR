import numpy as np
from pyscf import gto, scf, mcscf, ao2mo, symm, tools, fci, mrpt, dft
import sys
from copy import deepcopy

# sys.path.append('/path/to/gammcor/')

try:
    from Pyscf2Gammcor import get_data_for_gammcor
except ImportError:
    print("""
    ====================================================================
    !! ACTION REQUIRED: SCRIPT CONFIGURATION !!

    Could not find the 'Pyscf2Gammcor' module. You need to configure
    the paths in this script before running it.

    Please open this file in a text editor and follow these steps:

    1. Find the line below (around line 7):
       #sys.path.append('/path/to/gammcor/')

       a) Uncomment this line by removing the '#' at the beginning.
       b) Replace the example path with the correct, full path to your
          'gammcor' directory.

    2. Find the line inside the main part of the script:
       basis_file = '/path/to/gammcor/bazy-do-pyscf/cc-pvdz.nw'

       a) Replace this example path with the correct, full path to your
          basis file.

    After editing and saving the file, run the script again.
    ====================================================================
    """)
    sys.exit(1) 
    

if __name__ == "__main__":


    basis_file = '/path/to/gammcor/bazy-do-pyscf/cc-pvdz.nw'
    
    mol = gto.Mole()
    mol.basis = basis_file
    mol.unit = 'bohr'
    mol.charge = 0
    mol.spin = 0
    mol.cart = False
    mol.symmetry = True
    mol.symmetry_subgroup = 'C2v'

    mol.build(atom='''B    0.000000000    0.000000000   -0.426444880
    H    0.000000000    0.000000000    4.573555120
    ''')

    myhf = scf.RHF(mol)
    myhf.kernel()

    nelec = mol.nelectron


    ncas = {'A1': 3, 'B1': 1,  'B2': 1,  'A2': 0}
    ncore = {'A1': 1}

    mycas = mcscf.CASSCF(myhf, 5, 4)
    mcscf.addons.state_specific_(mycas, state=0)
    
    mycas.conv_tol = 1e-9
    mycas.conv_tol_grad = 1e-5
    
    mycas.natorb = True
    mo = mcscf.sort_mo_by_irrep(mycas, myhf.mo_coeff, ncas, ncore)
    mycas.fix_spin_(ss=0.0)
    mycas.fcisolver.wfnsym = 'A1'
    mycas.kernel(mo)

    get_data_for_gammcor(mol, myhf, mycas)

