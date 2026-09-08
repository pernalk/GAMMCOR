#!/bin/bash
set -e

# install prefix = the directory holding this script
PREFIX="$(cd "$(dirname "$0")" && pwd)"

mkdir -p "$PREFIX/src"
cd "$PREFIX/src"

wget -c https://github.com/HDFGroup/hdf5/releases/download/hdf5_1.14.5/hdf5-1.14.5.tar.gz
tar xzf hdf5-1.14.5.tar.gz

mkdir -p hdf5-build
cd hdf5-build

../hdf5-1.14.5/configure --prefix="$PREFIX" \
    CC=icx FC=ifx CFLAGS=-O2 FCFLAGS=-O2 \
    --enable-fortran --enable-hl \
    --enable-static --disable-shared \
    --disable-tests --disable-tools \
    --with-zlib=/usr

make -j8
make install
