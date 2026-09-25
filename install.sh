#!/bin/bash
# Builds GAMMCOR from source: gammcor-integrals, XCFun, HDF5, gammcor_tools,
# gammcor. Steps that are already done are skipped, so running it again is cheap.
#
#   ./install.sh [--prefix DIR] [--cpu TARGET] [--jobs N] [--clean]
#
#   --prefix DIR   also install gammcor into DIR/bin (conda recipe: --prefix "$PREFIX")
#   --cpu TARGET   auto (default) | native | avx512 | avx512-amd | avx2 | dispatch
#   --jobs N       parallel build jobs (default: nproc)
#   --clean        rebuild everything
#
# After the first build, recompile with:  meson compile -C build
set -eE -o pipefail

cd "$(dirname "$0")"

PREFIX=""
CPU=auto
JOBS=$(nproc)
CLEAN=0

while [ $# -gt 0 ]; do
    case "$1" in
        --prefix) PREFIX="$(realpath -m "$2")"; shift 2 ;;
        --cpu)    CPU="$2";  shift 2 ;;
        --jobs)   JOBS="$2"; shift 2 ;;
        --clean)  CLEAN=1;   shift ;;
        -h|--help) sed -n '2,12p' "$0"; exit 0 ;;
        *) echo "Unknown option: $1 (see ./install.sh --help)" >&2; exit 1 ;;
    esac
done
export JOBS

STEP=""
step() { STEP="$1"; echo; echo "=== $1"; }
die()  { echo; echo "ERROR: $*" >&2; exit 1; }
trap 'echo; echo "ERROR: step \"$STEP\" failed." >&2' ERR

# ---------------------------------------------------------------------------
# 1. Environment checks
# ---------------------------------------------------------------------------
step "Checking the environment"

command -v ifx > /dev/null || die "ifx not found.

GAMMCOR is compiled with the Intel oneAPI compilers. Load Intel oneAPI into
your environment.

On a workstation this may be:

  source /opt/intel/oneapi/setvars.sh

On an HPC cluster, load the appropriate Intel compiler module
(e.g. \`module avail intel\`)."

command -v icx > /dev/null || die "icx not found (Intel oneAPI C/C++ compiler), see above for ifx."

if [ -z "$MKLROOT" ] && [ ! -f "$CONDA_PREFIX/lib/libmkl_rt.so" ]; then
    die "MKL not found: \$MKLROOT is not set and there is no \$CONDA_PREFIX/lib/libmkl_rt.so."
fi

for tool in python3 meson ninja make; do
    command -v $tool > /dev/null || die "$tool not found."
done
python3 -m pip --version > /dev/null 2>&1 || die "pip not found for $(command -v python3)."

# ---------------------------------------------------------------------------
# 2. gammcor-integrals submodule
# ---------------------------------------------------------------------------
step "gammcor-integrals submodule"

if [ ! -f gammcor-integrals/meson.build ]; then
    git submodule update --init || die "gammcor-integrals is missing. Clone with: git clone --recursive"
fi

# ---------------------------------------------------------------------------
# 3. CPU target (one target for gammcor-integrals and gammcor)
# ---------------------------------------------------------------------------
step "CPU target"

if [ "$CPU" = auto ]; then
    CPU=$(python3 gammcor-integrals/.meson/scripts/detect_cpu.py | head -n 1)
fi
echo "cpu: $CPU"

# ---------------------------------------------------------------------------
# --clean: remove everything built by this script
# ---------------------------------------------------------------------------
if [ $CLEAN = 1 ]; then
    step "Cleaning"
    rm -rf build gammcor-integrals/build hdf5/lib hdf5/include
    make -C xcfun clean
fi

# Configure a Meson build directory, or update its options if it exists.
#   meson_setup <build dir> <source dir> [options...]
meson_setup() {
    local dir=$1 src=$2; shift 2
    if [ -d "$dir/meson-private" ]; then
        meson configure "$dir" "$@"
    else
        meson setup "$dir" "$src" --native-file "$src/.meson/profiles/ifx.ini" "$@"
    fi
}

# ---------------------------------------------------------------------------
# 4. gammcor-integrals  ->  gammcor-integrals/lib/cholesky.a
# ---------------------------------------------------------------------------
step "Building gammcor-integrals"

meson_setup gammcor-integrals/build gammcor-integrals -Dcpu="$CPU"
meson compile -C gammcor-integrals/build -j "$JOBS"

# ---------------------------------------------------------------------------
# 5. XCFun  ->  xcfun/lib/libxcfun.a
# ---------------------------------------------------------------------------
step "Building XCFun"

make -C xcfun intel

# ---------------------------------------------------------------------------
# 6. HDF5  ->  hdf5/lib, hdf5/include
# ---------------------------------------------------------------------------
step "Building HDF5"

if [ -f hdf5/lib/libhdf5_fortran.a ]; then
    echo "already built (rebuild: rm -rf hdf5/lib hdf5/include, or --clean)"
else
    ./hdf5/build_hdf5.sh
fi

# ---------------------------------------------------------------------------
# 7. gammcor_tools (Python), into the current python3
# ---------------------------------------------------------------------------
step "Installing gammcor_tools"

if [ -z "$PREFIX" ]; then
    PIP_ARGS=(-e ./gammcor_tools)
else
    PIP_ARGS=(./gammcor_tools --no-deps --no-build-isolation)
fi
python3 -m pip install -q "${PIP_ARGS[@]}" || die "pip install failed.
If pip says 'externally-managed-environment', use a venv, a conda env,
or run: python3 -m pip install --user ${PIP_ARGS[*]}"

# ---------------------------------------------------------------------------
# 8. gammcor  ->  build/gammcor, bin/gammcor
# ---------------------------------------------------------------------------
step "Building gammcor"

MESON_ARGS=(-Dcpu="$CPU")
if [ -n "$PREFIX" ]; then MESON_ARGS+=(-Dprefix="$PREFIX"); fi

meson_setup build . "${MESON_ARGS[@]}"
meson compile -C build -j "$JOBS"

# ---------------------------------------------------------------------------
# 9. Install (only with --prefix)
# ---------------------------------------------------------------------------
if [ -n "$PREFIX" ]; then
    step "Installing into $PREFIX"
    meson install -C build
fi

echo
echo "=== Done. Executable: $PWD/bin/gammcor${PREFIX:+, $PREFIX/bin/gammcor}"
