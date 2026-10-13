import numpy as np
import pyscf
from pyscf import gto, scf, mcscf, symm, fci
from Pyscf2Gammcor import get_data_for_gammcor

if __name__ == "__main__":
    print("PySCF version:", pyscf.__version__, flush=True)

    basis_file = '../def2-svp.nw'
    
    mol = gto.Mole()
    mol.basis = basis_file
    mol.unit = 'angs'
    mol.charge = 0
    mol.spin = 0
    mol.cart = False
    mol.symmetry = True
    mol.build(atom='''H     0.0000000     2.1367320     -0.1124450
H     0.0000000     -2.1367320     -0.1124450
H     -0.8813340     1.3337330     -1.4438420
H     0.8813340     -1.3337330     -1.4438420
H     -0.8813340     -1.3337330     -1.4438420
H     0.8813340     1.3337330     -1.4438420
C     0.0000000     0.0000000     0.0000000
C     0.0000000     1.2872530     -0.7959020
C     0.0000000     -1.2872530     -0.7959020
O     0.0000000     0.0000000     1.2276000
''')

    myhf = scf.RHF(mol)
    myhf.conv_tol = 1e-10
    myhf.max_cycle = 200
    myhf.kernel()
    if not myhf.converged:
        myhf = myhf.newton()
        myhf.kernel()
    if not myhf.converged:
        raise RuntimeError("RHF did not converge")
    print("Point groups:", mol.topgroup, mol.groupname)
    labels = symm.label_orb_symm(mol, mol.irrep_name, mol.symm_orb, myhf.mo_coeff)
    for i, (energy, occ, label) in enumerate(zip(myhf.mo_energy, myhf.mo_occ, labels), 1):
        print(i, label, occ, energy)

    nelec = mol.nelectron
    print(f"Total number of electrons: {nelec}")
    
    ncas = {'A1': 2, 'B1': 1, 'B2': 2}
    ncore ={'A1': 6, 'B1': 4, 'B2': 1, 'A2': 1}

    n_states = 4
    weights = np.ones(n_states)/n_states

    solver1 = fci.direct_spin0_symm.FCI(mol)
    solver1.wfnsym= 'A1'
    solver1.spin = 0
    solver1 = fci.addons.fix_spin(solver1, shift=.2, ss=0)
    solver1.nroots = n_states


    mycas = mcscf.CASSCF(myhf, 5, 8)
    mcscf.state_average_mix_(mycas, [solver1], weights)
    
    mycas.conv_tol = 1e-9
    mycas.conv_tol_grad = 1e-5
    print('mycas.conv_tol', mycas.conv_tol)
    print('mycas.conv_tol_grad', mycas.conv_tol_grad)

    mo = mcscf.sort_mo_by_irrep(mycas, myhf.mo_coeff, ncas, ncore)
    mycas.frozen = 0
    mycas.max_cycle_macro = 150
    mycas.kernel(mo)
    print("CASSCF converged:", mycas.converged)
    if not mycas.converged:
        raise RuntimeError("CASSCF did not converge; no GAMMCOR data exported")

    eorb = mycas.mo_energy

    get_data_for_gammcor(mol, myhf, mycas, dump_eri=True, simple=True)
