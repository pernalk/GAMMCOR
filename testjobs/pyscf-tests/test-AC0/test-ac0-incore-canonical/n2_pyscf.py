# N2 CAS(6,6), D2h, natorb = False: NATORB 0, ERI in AO, GAMMCOR does MO->NO (Trans2NO_AO).
# All four test-ac0-incore-* cases must give the same AC0 energy.
import numpy as np
from pyscf import gto, scf, mcscf, ao2mo, symm, tools, fci, mrpt
import sys
from copy import deepcopy
# sys.path.append('/path/to/gammcor/gammcor_tools/')

from Pyscf2Gammcor import get_data_for_gammcor

if __name__ == "__main__":


    basis_file = '/home/aleksandra.tucholska/test/gammcor/bazy-do-pyscf/cc-pvdz.nw'
    
    mol = gto.Mole()
    mol.basis = basis_file
    mol.unit = 'bohr'
    mol.charge = 0
    mol.spin = 0
    mol.cart = False
    mol.symmetry = True
    mol.symmetry_subgroup = 'D2h'

    mol.build(atom='''N    0.000000000    0.000000000   0.0000000000
N    0.000000000    0.000000000    2.07500
''')

    myhf = scf.RHF(mol)
    myhf.kernel()

    nelec = mol.nelectron

    ncas = {'Ag': 1, 'B2g': 1,  'B3g': 1,  'B1u': 1,  'B2u': 1,  'B3u': 1}
    ncore = {'Ag': 2, 'B1u': 2}

    mycas = mcscf.CASSCF(myhf, 6, 6)
    mcscf.addons.state_specific_(mycas, state=0)
    
    mycas.conv_tol = 1e-9
    mycas.conv_tol_grad = 1e-5
    
    mycas.natorb = False
    mo = mcscf.sort_mo_by_irrep(mycas, myhf.mo_coeff, ncas, ncore)
    mycas.fix_spin_(ss=0.0)
    mycas.fcisolver.wfnsym = 'Ag'
    mycas.kernel(mo)

    get_data_for_gammcor(mol, myhf, mycas, dump_eri=True)

