# Setup

## Contents

1. [For developers](#1-for-developers)
   - [With conda (recommended)](#11-with-conda-recommended)
   - [With an existing toolchain](#12-with-an-existing-toolchain)
   - [Testing](#13-testing)
   - [Recompiling after a source change](#14-recompiling-after-a-source-change)
   - [Building with Make](#15-building-with-make)
   - [Debug build](#16-debug-build)
   - [Executables](#17-executables)
2. [For non-developers](#2-for-non-developers)

## 1. For developers

GAMMCOR supports Linux on x86-64 Intel and AMD processors with AVX2.

Clone the repository with its submodules:

```bash
git clone --recursive https://gitlab.com/qchem/gammcor.git
cd gammcor
```

### 1.1. With conda (recommended)

Conda installs the compilers, libraries, and build tools listed in
`environment.yml`. Compiler or MKL modules should not be loaded in the same
shell.

```bash
conda env create -f environment.yml
conda activate gammcor-dev
./install.sh
```

The current conda build has been tested with `ifx` 2026.1.1, `icx` 2026.1.1,
oneMKL 2026.1.0, Python 3.13.15, Meson 1.12.1, and Ninja 1.13.2. These are
tested versions, not minimum requirements unless stated below.

### 1.2. With an existing toolchain

Prerequisites:

- Intel Fortran Compiler (`ifx`; tested with 2025.0.4 and 2026.1.1)
- Intel C/C++ Compiler (`icx`; tested with 2025.0.4 and 2026.1.1)
- Intel oneMKL (tested with 2025.0 and 2026.1.0)
- GCC C++ compiler (`g++`)
- Python 3.10 or newer, with `pip`, NumPy, and h5py
- Meson 1.1 or newer
- Ninja (tested with 1.13.2)
- zlib, including the `zlib.h` header

Build GAMMCOR and all necessary libraries:

```bash
./install.sh
```

The script builds GAMMCOR and its bundled dependencies. Run
`./install.sh --help` to see the available build options.

### 1.3. Testing

Run the test suite from the `testjobs` directory:

```bash
cd testjobs
./run_tests.py gammcor
```

### 1.4. Recompiling after a source change

After editing files in `SOURCE/`, rebuild GAMMCOR with:

```bash
meson compile -C build
```

If you changed `gammcor-integrals`, run the installation script again. The
build is incremental, so unchanged files are not recompiled:

```bash
./install.sh
```

Run `./install.sh --clean` after changing a public module interface (arguments,
types, or module variables), changing the compiler version, or when a full
rebuild is required.

### 1.5. Building with Make

Once the required libraries have been built, GAMMCOR can also be compiled with
the top-level Makefile. The following command runs one build job per
available CPU core:

```bash
make -j $(nproc)
```

### 1.6. Debug build

The repository includes a separate Meson profile for debug builds:

```bash
meson setup build-debug --native-file .meson/profiles/ifx-debug.ini
meson compile -C build-debug
```

### 1.7. Executables

The compiled executables are written to:

- Meson release build: `<path-to-gammcor>/bin/gammcor`
- Meson debug build: `<path-to-gammcor>/bin/gammcor-debug`
- Make build: `<path-to-gammcor>/gammcor`

---

## 2. For non-developers

> **To be added.** The conda package has not been published yet. Replace
> `<channel>` with the GAMMCOR channel once it becomes available.

Install the precompiled package in your active conda environment:

```bash
conda install -c <channel> -c conda-forge gammcor
```

This installs GAMMCOR and its command-line tools. A compiler and the source
code are not required.

---

<div align="center">

| Previous | | Next |
|:--|:--:|--:|
|  | [↑ Back to README](../README.md) | [Input →](02_input_structure.md) |

</div>
