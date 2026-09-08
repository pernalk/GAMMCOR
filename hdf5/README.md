# HDF5 for GammCor

Needed by `SOURCE/h5_reader.f90`, which reads `pyscf_data.h5`.

The compiled libraries are committed, so normally you do nothing.

Rebuild only if the compiler says `cannot open module file hdf5.mod`
(it means your `ifx` differs from the one these were built with):

    cd hdf5
    ./build_hdf5.sh
