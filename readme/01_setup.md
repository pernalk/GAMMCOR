# Setup

## 1. For non-developers

Install the precompiled package in your active conda environment:

```bash
conda install -c <channel> -c conda-forge gammcor
```

This installs GAMMCOR and its command-line tools. A compiler and the source
code are not required.

> The conda package has not been published yet. Replace `<channel>` with the
> GAMMCOR channel once it becomes available.

## 2. For developers

GAMMCOR supports Linux on x86-64 Intel and AMD processors with AVX2.

Clone the repository with its submodules:

```bash
git clone --recursive https://gitlab.com/qchem/gammcor.git
cd gammcor
```

### 2.1. With conda (recommended)

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

### 2.2. Without conda

The source build requires the following tools to be available in the current
environment:

- Intel oneAPI compilers (`ifx` and `icx`)
- Intel oneMKL
- Python 3.10 or newer with `pip`
- Meson 1.1 or newer and Ninja
- Make and Git
- `wget` or `curl`
- zlib development files

On a computing cluster, first check which Python modules are available and
load one providing Python 3.10 or newer. The exact module name depends on the
cluster:

```bash
module avail python
module load <python-module>
```

Create a private Python environment and install Meson and Ninja in it:

```bash
python3 -m venv ~/.venvs/gammcor
source ~/.venvs/gammcor/bin/activate
python -m pip install --upgrade pip
python -m pip install "meson>=1.1" ninja
```

Then build GAMMCOR:

```bash
./install.sh
```

In a new terminal, load the same compiler, MKL, and Python modules and activate
`~/.venvs/gammcor` before rebuilding GAMMCOR.

The script builds GAMMCOR and its bundled dependencies. Run
`./install.sh --help` to see the available build options.

### 2.3. Recompiling after a source change

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

### 2.4. Building with Make

Once the required libraries have been built, GAMMCOR can also be compiled with
the top-level Makefile. For example, to use eight parallel jobs:

```bash
make -j 8
```

Replace `8` with the number of CPU cores you want to use for compilation.

### 2.5. Debug build

The repository includes a separate Meson profile for debug builds:

```bash
meson setup build-debug --native-file .meson/profiles/ifx-debug.ini
meson compile -C build-debug
```

### 2.6. Executables

The compiled executables are written to:

- Meson release build: `<path-to-gammcor>/bin/gammcor`
- Meson debug build: `<path-to-gammcor>/bin/gammcor-debug`
- Make build: `<path-to-gammcor>/gammcor`
