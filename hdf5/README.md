# HDF5 for GammCor

Needed by `SOURCE/h5_reader.f90`, which reads `pyscf_data.h5`.

HDF5 must be built once, with the same `ifx` that compiles GammCor
(its Fortran `.mod` files only work with that compiler version):

    cd hdf5
    ./build_hdf5.sh

This downloads HDF5 1.14.5 and installs it into `hdf5/include` and `hdf5/lib`.
These are not in git. Rebuild after changing the `ifx` version
(error `cannot open module file hdf5.mod`).

Optional environment variables: `JOBS` (parallel make jobs, default `nproc`),
`ZLIB_ROOT` (zlib prefix, default `/usr`, else `$CONDA_PREFIX`),
`HDF5_TARBALL` (use this tarball instead of downloading).
