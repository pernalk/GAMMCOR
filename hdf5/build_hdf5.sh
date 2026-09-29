#!/bin/bash
set -e

# Builds static HDF5 (with Fortran) into the directory holding this script.
# Position-independent (--with-pic), so that it can be linked into the shared
# libtrexio.so that install.sh builds against it.
#
# Environment variables (all optional):
#   JOBS          parallel make jobs                  (default: nproc)
#   ZLIB_ROOT     zlib install prefix                 (default: /usr if it has
#                 zlib.h, else $CONDA_PREFIX). install.sh sets it to
#                 $CONDA_PREFIX in conda mode, so the default order only
#                 applies in modules mode.
#   HDF5_TARBALL  path to hdf5-$VERSION.tar.gz; no download (conda recipe)

VERSION=1.14.5

# install prefix = the directory holding this script
PREFIX="$(cd "$(dirname "$0")" && pwd)"

JOBS="${JOBS:-$(nproc)}"

if [ -z "$ZLIB_ROOT" ]; then
    if [ -f /usr/include/zlib.h ]; then
        ZLIB_ROOT=/usr
    elif [ -n "$CONDA_PREFIX" ] && [ -f "$CONDA_PREFIX/include/zlib.h" ]; then
        ZLIB_ROOT="$CONDA_PREFIX"
    else
        echo "ERROR: zlib.h not found." >&2
        echo "Install zlib (zlib1g-dev on Ubuntu, or 'conda install zlib')," >&2
        echo "or set ZLIB_ROOT to its install prefix." >&2
        exit 1
    fi
fi
echo "build_hdf5: HDF5 $VERSION, zlib from $ZLIB_ROOT, $JOBS jobs"

mkdir -p "$PREFIX/src"
cd "$PREFIX/src"

if [ -n "$HDF5_TARBALL" ]; then
    TARBALL="$HDF5_TARBALL"
else
    TARBALL="hdf5-$VERSION.tar.gz"
    URL="https://github.com/HDFGroup/hdf5/releases/download/hdf5_$VERSION/$TARBALL"
    if command -v wget > /dev/null; then
        wget -c "$URL"
    else
        curl -fL -C - -o "$TARBALL" "$URL"
    fi
fi
tar xzf "$TARBALL"

# fresh build directory: objects and .mod files from an older ifx must not be reused
rm -rf hdf5-build
mkdir hdf5-build
cd hdf5-build

../hdf5-$VERSION/configure --prefix="$PREFIX" \
    CC=icx FC=ifx CFLAGS=-O2 FCFLAGS=-O2 \
    --enable-fortran --enable-hl \
    --enable-static --disable-shared --with-pic \
    --disable-tests --disable-tools \
    --with-zlib="$ZLIB_ROOT"

make -j"$JOBS"
make install
