"""LiF: CAS(6,6), A1, SA-CASSCF over 2 states.
Source: Li-F distance 1.56 angstrom (near equilibrium); SS for root 2 does not converge at 3.51
"""
import numpy as np
import pyscf
from pyscf import gto, scf, mcscf, symm, fci
from gammcor_tools.Pyscf2Gammcor import get_data_for_gammcor

if __name__ == "__main__":
    print("PySCF version:", pyscf.__version__, flush=True)

    mol = gto.Mole()
    mol.basis = '../cc-pvdz.nw'
    mol.unit = 'angs'
    mol.charge = 0
    mol.spin = 0
    mol.cart = False
    mol.symmetry = True
    mol.symmetry_subgroup = 'C2v'
    mol.build(atom='Li 0.0 0.0 0.0\nF 0.0 0.0 1.56')

    myhf = scf.RHF(mol)
    myhf.conv_tol = 1e-10
    myhf.max_cycle = 200
    myhf.kernel()
    if not myhf.converged:
        myhf = myhf.newton()
        myhf.kernel()
    if not myhf.converged:
        raise RuntimeError("RHF did not converge")
    print("Total number of electrons:", mol.nelectron)
    print("Point groups:", mol.topgroup, mol.groupname)
    labels = symm.label_orb_symm(mol, mol.irrep_name, mol.symm_orb, myhf.mo_coeff)
    for i, (energy, occ, label) in enumerate(zip(myhf.mo_energy, myhf.mo_occ, labels), 1):
        print(i, label, occ, energy)

    ncas = {'A1': 2, 'B1': 2, 'B2': 2}
    ncore = {'A1': 3}

    n_states = 2
    weights = np.ones(n_states)/n_states

    solver1 = fci.direct_spin0_symm.FCI(mol)
    solver1.wfnsym = 'A1'
    solver1.spin = 0
    solver1 = fci.addons.fix_spin(solver1, shift=.2, ss=0)
    solver1.nroots = n_states

    mycas = mcscf.CASSCF(myhf, 6, 6)
    mcscf.state_average_mix_(mycas, [solver1], weights)

    mycas.conv_tol = 1e-9
    mycas.conv_tol_grad = 1e-5
    mycas.max_cycle_macro = 150
    mycas.frozen = 0
    mo = mcscf.sort_mo_by_irrep(mycas, myhf.mo_coeff, ncas, ncore)
    mycas.kernel(mo)
    print("CASSCF converged:", mycas.converged)
    if not mycas.converged:
        raise RuntimeError("CASSCF did not converge; no GAMMCOR data exported")

    get_data_for_gammcor(mol, myhf, mycas, dump_eri=True, simple=True)
