import pyscf
from pyscf import gto, scf, mp
from gammcor_tools.Pyscf2Gammcor import get_data_for_gammcor

if __name__ == "__main__":
    print("PySCF version:", pyscf.__version__, flush=True)

    basis_file = '../../cc-pvdz.nw'

    mol = gto.Mole()
    mol.basis = basis_file
    mol.unit = 'angs'
    mol.charge = 0
    mol.spin = 0
    mol.cart = False
    mol.symmetry = True
    mol.build(atom='''N     0.0000000     0.0000000     -0.5570000
N     0.0000000     0.0000000      0.5570000
''')

    myhf = scf.RHF(mol)
    myhf.conv_tol = 1e-10
    myhf.max_cycle = 200
    myhf.kernel()
    if not myhf.converged:
        raise RuntimeError("RHF did not converge")

    # optional: only for checks in GAMMCOR (REF/FROZEN vs NCoreOrb, PySCF E_corr)
    mymp = mp.MP2(myhf).set_frozen().run()
    print("PySCF MP2 frozen:", mymp.frozen, " E_corr:", mymp.e_corr)

    # no mycas -> HF dump for JobType MP2
    get_data_for_gammcor(mol, myhf, mymp=mymp)
