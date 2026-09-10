# gammcor_tools

Python tools for preparing and inspecting GammCor input data.

## Install

From the root of the GammCor repository:

```bash
pip install -e gammcor_tools
```

or, equivalently, `pip install -e .` from inside `gammcor_tools/`. Add the
`pyscf` extra -- `pip install -e 'gammcor_tools[pyscf]'` -- if PySCF is not
already installed; only `gammcor_basis` and the `Pyscf2Gammcor` interface
need it.
The `-e` option installs the package in editable mode, so changes in the source code and `git pull` are used immediately without reinstalling.

## Commands

```bash
gammcor_basis script.py -o basis.txt
```

Extract the basis set from a PySCF input script in format compatible with gammcor (GAMESS_US).

```bash
gammcor_2h5 JOBDIR -o pyscf_data.h5
```

Convert legacy `.bin` / `.txt` GammCor input files to HDF5.

```bash
gammcor_2bin pyscf_data.h5 -o JOBDIR
```

Convert an HDF5 file back to the legacy file format.

```bash
gammcor_h5info pyscf_data.h5
```

Show the contents of a GammCor HDF5 file.

Use `-h` with any command to see available options.

## Python interface

Recommended import:

```python
from gammcor_tools.Pyscf2Gammcor import get_data_for_gammcor
```

Older scripts using:

```python
sys.path.append('/path/to/gammcor/gammcor_tools/')
from Pyscf2Gammcor import get_data_for_gammcor
```

are still supported.

`Pyscf2Gammcor_legacy.py` contains the previous implementation and is kept only for reference.