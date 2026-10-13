"""LiF: CAS(6,6), A1, SA-CASSCF over 2 states.
Source: Li-F distance 3.51 angstrom
"""
import numpy as np
from pyscf import gto, scf, mcscf, symm
from Pyscf2Gammcor import get_data_for_gammcor

if __name__ == "__main__":
    mol = gto.Mole()
    mol.basis = '../../../basis-for-pyscf/cc-pvdz.nw'
    mol.unit = 'angs'
    mol.charge = 0
    mol.spin = 0
    mol.cart = False
    mol.symmetry = True
    mol.symmetry_subgroup = 'C2v'
    mol.build(atom='Li 0.0 0.0 0.0\nF 0.0 0.0 3.51')

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
    mycas = mcscf.CASSCF(myhf, 6, 6)
    mycas.fix_spin_(ss=0.0)
    mycas.fcisolver.wfnsym = 'A1'
    mycas = mycas.state_average_([0.5, 0.5])
    mycas.natorb = False

    mycas.conv_tol = 1e-9
    mycas.conv_tol_grad = 1e-5
    mycas.max_cycle_macro = 150
    mycas.frozen = 0
    mo = mcscf.sort_mo_by_irrep(mycas, myhf.mo_coeff, ncas, ncore)
    mycas.kernel(mo)
    print("CASSCF converged:", mycas.converged)
    if not mycas.converged:
        raise RuntimeError("CASSCF did not converge; no GAMMCOR data exported")
    get_data_for_gammcor(mol, myhf, mycas, dump_eri=False, simple=True, dump_hdf5=True)
