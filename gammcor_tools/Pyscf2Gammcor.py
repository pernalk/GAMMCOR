"""Compatibility shim -- the real module now lives in ``gammcor_tools``.

Existing PySCF input scripts do::

    sys.path.append('/path/to/gammcor/')
    from Pyscf2Gammcor import get_data_for_gammcor

and keep working through this file, both from the repository root and after
``pip install``.  New scripts should import the package directly::

    from gammcor_tools.Pyscf2Gammcor import get_data_for_gammcor

The previous, larger interface is preserved as
``gammcor_tools/Pyscf2Gammcor_legacy.py``.
"""

from gammcor_tools.Pyscf2Gammcor import *  # noqa: F401,F403
from gammcor_tools.Pyscf2Gammcor import (  # noqa: F401
    get_data_for_gammcor,
    get_rohf_for_gammcor,
)
