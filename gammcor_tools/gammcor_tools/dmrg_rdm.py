import os
from types import MethodType

import numpy as np
from pyscf.dmrgscf.dmrgci import DMRGCI


def add_make_rdm12s(solver):
    def make_rdm12s(self, state, norb, nelec, link_index=None, **kwargs):
        name = '2pdm.npy' if self.nroots == 1 else f'2pdm-{state}-{state}.npy'
        dm2s = np.load(os.path.join(self.scratchDirectory, 'node0', name))
        if dm2s.shape != (3, norb, norb, norb, norb):
            raise ValueError(f'Unexpected 2pdm shape: {dm2s.shape}')
        dm1s = DMRGCI.make_rdm1s(self, state, norb, nelec, link_index, **kwargs)
        return dm1s, tuple(dm2s.transpose(0, 1, 4, 2, 3))

    solver.make_rdm12s = MethodType(make_rdm12s, solver)
    return solver
