import h5py
import pyscf
from pyscf import gto, scf, mp, mcscf
from gammcor_tools.Pyscf2Gammcor import get_data_for_gammcor

if __name__ == "__main__":
    print("PySCF version:", pyscf.__version__, flush=True)

    basis_file = '../../aug-cc-pvdz.nw'

    mol = gto.Mole()
    mol.basis = basis_file
    mol.unit = 'angs'
    mol.charge = 0
    mol.spin = 3
    mol.cart = False
    mol.symmetry = True
    mol.symmetry_subgroup = "D2h"
    mol.build(atom='''N     0.0000000     0.0000000     0.00000
''')

    myhf = scf.RHF(mol)
    myhf.conv_tol = 1e-10
    myhf.max_cycle = 200
    myhf.kernel()
    if not myhf.converged:
        raise RuntimeError("RHF did not converge")

    # frozen core: REF/FROZEN = 1
    mymp = mp.MP2(myhf).set_frozen().run()
    print("PySCF MP2 frozen:", mymp.frozen, " E_corr:", mymp.e_corr)

    # fake CAS: 3 singly occupied 2p orbitals with 3 electrons in CAS(3e,3o), E_CASSCF = E_HF
    mycas = mcscf.CASSCF(myhf, 3, 3)
    mycas.natorb = True
    mycas.conv_tol = 1e-9
    mycas.conv_tol_grad = 1e-5
    mycas.kernel()
    if not mycas.converged:
        raise RuntimeError("CASSCF did not converge")

    # CASSCF dump, written to pyscf_data.h5 in this folder
    get_data_for_gammcor(mol, myhf, mycas, mymp)

    with h5py.File('pyscf_data.h5') as f:
        print('CAS', f.attrs['METHOD'], 'NI', f['REF/NI'][0], 'NA', f['REF/NA'][0],
              'FROZEN', f['REF/FROZEN'][0], 'ECAS', f['POSTHF/STATES/1.1/ECAS'][0])
