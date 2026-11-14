# N2 CAS(6,6), no symmetry, CASCI on CASSCF orbitals with a randomly rotated active block: NATORB 0, ERI in AO, active 1-RDM non-diagonal, so Trans2NO_AO does a real MO->NO rotation.
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
    mol.symmetry = False

    mol.build(atom='''N    0.000000000    0.000000000   0.0000000000
N    0.000000000    0.000000000    2.07500
''')

    myhf = scf.RHF(mol)
    myhf.kernel()

    nelec = mol.nelectron

    mycas = mcscf.CASSCF(myhf, 6, 6)
    mcscf.addons.state_specific_(mycas, state=0)
    
    mycas.conv_tol = 1e-9
    mycas.conv_tol_grad = 1e-5
    
    mycas.natorb = False
    mo = myhf.mo_coeff
    mycas.fix_spin_(ss=0.0)
    mycas.kernel(mo)


    # Rotate the active block by a fixed random orthogonal matrix and redo the
    # CAS as a CASCI without canonicalization: same energy, but the active
    # 1-RDM is far from diagonal, so GAMMCOR must do a real MO -> NO rotation.
    rng = np.random.default_rng(7)
    U, _ = np.linalg.qr(rng.standard_normal((6, 6)))
    mo = mycas.mo_coeff.copy()
    mo[:, mycas.ncore:mycas.ncore+6] = mo[:, mycas.ncore:mycas.ncore+6] @ U
    mycasci = mcscf.CASCI(myhf, 6, 6)
    mycasci.fix_spin_(ss=0.0)
    mycasci.natorb = False
    mycasci.canonicalization = False
    mycasci.kernel(mo)

    get_data_for_gammcor(mol, myhf, mycasci, dump_eri=True)

