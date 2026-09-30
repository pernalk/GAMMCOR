"""Attach an AO embedding potential to a PySCF SCF object.

The potential is attached to one SCF object only; PySCF itself is not
modified::

    from gammcor_tools.embedding_potential import add_embedding_potential

    mf = dft.RKS(mol, xc='pbe').density_fit(auxbasis='def2tzvppjkfit')
    mf = add_embedding_potential(mf, 'emb_ints.pyscf')

Everything that takes its one-electron Hamiltonian from ``mf.get_hcore()``
sees the potential: the SCF itself, CASSCF/DMRGSCF built on ``mf`` (via
``mc._scf``), and ``get_data_for_gammcor`` (INTS/CORE_HAMILTONIAN).  Call it
after ``density_fit()`` / ``newton()`` and on every new SCF object.
"""

import os
from types import MethodType

import numpy as np


def load_embedding_potential(filename, nao):
    """Read and validate a real, symmetric nao x nao text matrix."""
    path = os.path.abspath(os.path.expanduser(filename))
    if not os.path.isfile(path):
        raise FileNotFoundError(f'Embedding potential not found: {path}')

    vemb = np.loadtxt(path)
    if vemb.shape != (nao, nao):
        raise ValueError(f'Embedding matrix in {path} has shape {vemb.shape}; '
                         f'expected {(nao, nao)} for this molecule')
    if not np.isfinite(vemb).all():
        raise ValueError(f'Embedding matrix in {path} contains NaN or inf')
    if not np.allclose(vemb, vemb.T, rtol=1e-10, atol=1e-12):
        raise ValueError(f'Embedding matrix in {path} is not symmetric')
    return path, vemb


def add_embedding_potential(mf, filename='emb_ints.pyscf'):
    """Make *mf* add the matrix in *filename* to its core Hamiltonian.

    The matrix must be in the AO ordering of ``mf.mol`` (PySCF ordering).
    The same SCF object is returned, with ``mf.embedding_potential`` and
    ``mf.embedding_potential_file`` set.
    """
    path, vemb = load_embedding_potential(filename, mf.mol.nao_nr())

    if getattr(mf, 'embedding_potential', None) is not None:
        raise RuntimeError('mf already has an embedding potential '
                           f'({mf.embedding_potential_file})')

    original_get_hcore = mf.get_hcore

    def get_hcore_with_embedding(self, mol=None):
        hcore = np.asarray(original_get_hcore(mol))
        if hcore.shape != vemb.shape:
            raise ValueError(f'Core Hamiltonian has shape {hcore.shape}, but '
                             f'the embedding matrix has shape {vemb.shape}')
        return hcore + vemb

    mf.get_hcore = MethodType(get_hcore_with_embedding, mf)
    mf.embedding_potential_file = path
    mf.embedding_potential = vemb
    return mf
