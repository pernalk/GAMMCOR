"""Compatibility shim -- the real module lives in ``gammcor_tools``.

Lets PySCF input scripts import it the same way as ``Pyscf2Gammcor``::

    from embedding_potential import add_embedding_potential

New scripts should import the package directly::

    from gammcor_tools.embedding_potential import add_embedding_potential
"""

from gammcor_tools.embedding_potential import *  # noqa: F401,F403
from gammcor_tools.embedding_potential import (  # noqa: F401
    add_embedding_potential,
    load_embedding_potential,
)
