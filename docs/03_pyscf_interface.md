# PySCF interface

## Contents

- [`get_data_for_gammcor` options](#get_data_for_gammcor-options)
- [Command-line tools](#command-line-tools)

With PySCF, the reference calculation is a Python script. At its end, call
`get_data_for_gammcor`, which writes `pyscf_data.h5` for GAMMCOR:

```python
from gammcor_tools.Pyscf2Gammcor import get_data_for_gammcor

# ... mol, myhf, mycas: your PySCF molecule, SCF and CASSCF objects ...

get_data_for_gammcor(mol, myhf, mycas)
```

HDF5 output is enabled by default, so there is no need to pass
`dump_hdf5=True`.

## `get_data_for_gammcor` options

The complete function signature is:

```python
get_data_for_gammcor(
    mol,
    myhf,
    mycas,
    mymp=None,
    dump_eri=False,
    simple=True,
    dump_hdf5=True,
    dump_dft=False,
    mp2_scale=None,
    dump_grid=False,
    tds=None,
    dump_ab=True,
)
```

The arguments are:

- `mol`: the PySCF `Mole` object.
- `myhf`: a converged SCF object. Use a converged RKS object when
  `dump_dft=True`.
- `mycas`: a converged CASSCF object. It is ignored when `dump_dft=True`, so
  pass `None` in that mode if no CASSCF object exists.
- `mymp`: an optional converged MP2 object. For a CASSCF export it supplies the
  frozen-orbital count. For a DFT export it also supplies the MP2 reference
  energy.
- `dump_eri=False`: do not store the full two-electron integral array. Set it
  to `True` for a calculation with `TwoMoInt INCORE`; this is the route used by
  ppAC0 and ppAC examples without THC.
- `simple=True`: for a CASSCF export, write only the files needed by GAMMCOR.
  Set it to `False` to also write Molden orbital files for inspection.
- `dump_hdf5=True`: for a CASSCF export, write the recommended single file
  `pyscf_data.h5`. This is the default and does not need to be stated in a
  normal call. Set it to `False` only when the legacy collection of `.bin` and
  `.txt` files is required. The DFT export always uses HDF5.
- `dump_dft=False`: export CASSCF reference data. Set it to `True` to export
  RKS/MP2/TDDFT data; in this mode `mycas` is ignored.
- `mp2_scale=None`: in a DFT export, an optional scale factor for the MP2
  contribution. This is used together with `mymp` for a double-hybrid
  reference.
- `dump_grid=False`: in a DFT export, set it to `True` to store the numerical
  grid, weights, AO values and derivatives, density, and exchange-correlation
  quantities.
- `tds=None`: in a DFT export, an optional iterable of converged TDDFT objects
  (for example, a list containing singlet and triplet calculations).
- `dump_ab=True`: in a DFT export with `tds`, also store the full A and B
  response matrices. Set it to `False` to omit them.

Use real Python booleans (`True` and `False`), not the strings `'True'` and
`'False'`.

## Command-line tools

The `gammcor_tools` package is installed by `install.sh` and by the conda
package. It also provides command-line tools:

| Command | Purpose |
|---|---|
| `gammcor_basis script.py -o basis.txt` | extract the basis set from a PySCF script, in GAMMCOR format |
| `gammcor_h5info pyscf_data.h5` | show the contents of a GAMMCOR HDF5 file |
| `gammcor_2h5 JOBDIR -o pyscf_data.h5` | convert legacy `.bin`/`.txt` input files to HDF5 |
| `gammcor_2bin pyscf_data.h5 -o JOBDIR` | convert HDF5 back to the legacy files |

`gammcor_h5info` can inspect the complete file, a group, or an individual
dataset:

```bash
# Print the complete HDF5 structure.
gammcor_h5info pyscf_data.h5

# Print the contents of the REF group.
gammcor_h5info pyscf_data.h5 REF

# Print only the ENUC dataset inside REF.
gammcor_h5info pyscf_data.h5 REF/ENUC
```

Here `REF` is the name of an HDF5 group and `REF/ENUC` is the full path to one
dataset in that group. Run any command with `-h` to see its remaining options.
Practical PySCF inputs are shown in [Examples](04_examples.md).

---

<div align="center">

| Previous | | Next |
|:--|:--:|--:|
| [← Input](02_input_structure.md) | [↑ Back to README](../README.md) | [Examples →](04_examples.md) |

</div>
