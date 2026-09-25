#!/bin/bash
# Builds GAMMCOR from source: gammcor-integrals, XCFun, HDF5, gammcor_tools,
# gammcor. Steps that are already done are skipped, so running it again is cheap.
#
#   ./install.sh [--prefix DIR] [--cpu TARGET] [--jobs N] [--clean] [--verbose]
#
#   --prefix DIR   also install gammcor into DIR/bin (conda recipe: --prefix "$PREFIX")
#   --cpu TARGET   auto (default) | native | avx512 | avx512-amd | avx2 | dispatch
#   --jobs N       parallel build jobs (default: nproc)
#   --clean        rebuild everything
#   --verbose      show compiler output on screen (it always goes to install.log)
#
# After the first build, recompile with:  meson compile -C build
set -eE -o pipefail

cd "$(dirname "$0")"

ARGS="$*"
PREFIX=""
CPU=auto
JOBS=$(nproc)
CLEAN=0
VERBOSE=0

while [ $# -gt 0 ]; do
    case "$1" in
        --prefix) PREFIX="$(realpath -m "$2")"; shift 2 ;;
        --cpu)    CPU="$2";  shift 2 ;;
        --jobs)   JOBS="$2"; shift 2 ;;
        --clean)  CLEAN=1;   shift ;;
        -v|--verbose) VERBOSE=1; shift ;;
        -h|--help) sed -n '2,13p' "$0"; exit 0 ;;
        *) echo "Unknown option: $1 (see ./install.sh --help)" >&2; exit 1 ;;
    esac
done
export JOBS

LOG="$PWD/install.log"
{
    echo "# ./install.sh $ARGS"
    echo "# $(date)  on $(hostname)"
} > "$LOG"

# ---------------------------------------------------------------------------
# Terminal UI
#
# On a terminal (without --verbose) the step list is redrawn in place, with a
# spinner on the running step. Otherwise (--verbose, or output redirected)
# each step prints one line when it finishes. Build tool output always goes
# to install.log, and with --verbose also to the screen.
# ---------------------------------------------------------------------------
if [ -t 1 ] && [ -z "$NO_COLOR" ]; then
    GREEN=$'\e[32m' YELLOW=$'\e[33m' RED=$'\e[31m' CYAN=$'\e[36m'
    DIM=$'\e[2m' BOLD=$'\e[1m' RESET=$'\e[0m'
else
    GREEN="" YELLOW="" RED="" CYAN="" DIM="" BOLD="" RESET=""
fi
if [ -t 1 ] && [ $VERBOSE = 0 ]; then LIVE=1; else LIVE=0; fi
COLS=$(tput cols 2> /dev/null || echo 80)
RULE=$(printf '─%.0s' $(seq 40))
SPIN=(⠋ ⠙ ⠹ ⠸ ⠼ ⠴ ⠦ ⠧ ⠇ ⠏)

# Steps: name, state (pending|running|ok|fail), one-line info, extra lines
NAMES=() STATE=() INFO=() EXTRA=()
add_step() {   # add_step <index variable> <name>
    printf -v "$1" %d ${#NAMES[@]}
    NAMES+=("$2"); STATE+=(pending); INFO+=(""); EXTRA+=("")
}

CUR=-1          # running step
RUNLABEL=""     # its label, e.g. "Building XCFun"
FOOTER=""
DRAWN=0         # lines of the step list currently on screen
SPIN_PID=""

# fit <text> <width>: shorten from the left, so that a line never wraps
fit() {
    local s=$1 w=$2
    if [ ${#s} -gt $w ]; then s="…${s: -$((w - 1))}"; fi
    printf '%s' "$s"
}

step_line() {
    local i=$1 line
    local name=${NAMES[$i]} info
    info=$(fit "${INFO[$i]}" $((COLS - 25)))
    case ${STATE[$i]} in
        pending) printf ' %s·  %s%s\n' "$DIM" "$name" "$RESET" ;;
        running) printf ' %s%s%s  %s\n' "$YELLOW" "${SPIN[0]}" "$RESET" "$RUNLABEL" ;;
        ok)      printf ' %s✓%s  %-18s %s%s%s\n' "$GREEN" "$RESET" "$name" "$CYAN" "$info" "$RESET" ;;
        fail)    printf ' %s✗%s  %-18s %s%s%s\n' "$RED" "$RESET" "$name" "$RED" "$info" "$RESET" ;;
    esac
    if [ -n "${EXTRA[$i]}" ]; then
        while IFS= read -r line; do
            printf '    %s%s%s\n' "$CYAN" "$(fit "$line" $((COLS - 6)))" "$RESET"
        done <<< "${EXTRA[$i]}"
        echo
    fi
}

# number of screen lines above the line of step <i>
lines_before() {
    local i n=0
    for ((i = 0; i < $1; i++)); do
        n=$((n + 1))
        if [ -n "${EXTRA[$i]}" ]; then
            n=$((n + $(printf '%s\n' "${EXTRA[$i]}" | wc -l) + 1))
        fi
    done
    echo $n
}

footer_text() {
    local s=$((SECONDS - STEP_SECONDS))
    printf 'Step %d of %d  %s%d:%02d%s' $((CUR + 1)) ${#NAMES[@]} "$DIM" $((s / 60)) $((s % 60)) "$RESET"
}

# Redraw the step list and footer over the previous drawing.
render() {
    local out i n
    out=$(
        for i in "${!NAMES[@]}"; do step_line $i; done
        echo
        echo "$DIM$RULE$RESET"
        echo "$FOOTER"
    )
    n=$(printf '%s\n' "$out" | wc -l)
    if [ $DRAWN -gt 0 ]; then printf '\e[%dA\r' $DRAWN; fi
    printf '%s\n' "$out" | sed $'s/$/\e[K/'
    if [ $n -lt $DRAWN ]; then printf '\e[J'; fi
    DRAWN=$n
}

# Background spinner: animates the running step's icon and the footer clock.
spin_start() {
    if [ $LIVE = 0 ]; then return 0; fi
    local up=$((DRAWN - $(lines_before $CUR)))
    (
        trap - ERR
        set +e
        k=0
        while kill -0 $$ 2> /dev/null; do
            k=$(((k + 1) % ${#SPIN[@]}))
            # up to the step line, icon; down to the footer, clock; back to the bottom
            printf '\e[%dA\r %s%s%s\e[%dB\r%s\e[K\e[1B\r' \
                $up "$YELLOW" "${SPIN[$k]}" "$RESET" $((up - 1)) "$(footer_text)"
            sleep 0.1
        done
    ) &
    SPIN_PID=$!
}

spin_stop() {
    if [ -n "$SPIN_PID" ]; then
        kill $SPIN_PID 2> /dev/null || true
        wait $SPIN_PID 2> /dev/null || true
        SPIN_PID=""
    fi
}

begin() {   # begin <step index> <label shown while running>
    CUR=$1
    RUNLABEL=$2
    STATE[$CUR]=running
    LOG_MARK=$(wc -l < "$LOG")
    T_STEP=$(now)
    STEP_SECONDS=$SECONDS
    printf '\n=== %s\n' "$RUNLABEL" >> "$LOG"
    if [ $LIVE = 1 ]; then
        FOOTER=$(footer_text)
        render
        spin_start
    elif [ $VERBOSE = 1 ]; then
        printf '\n %s▸%s  %s\n' "$YELLOW" "$RESET" "$RUNLABEL"
    fi
}

ok() {   # ok <info> [<extra lines>]
    spin_stop
    STATE[$CUR]=ok
    INFO[$CUR]=$1
    EXTRA[$CUR]=${2:-}
    printf '=== %s: ok %s\n' "${NAMES[$CUR]}" "$1" >> "$LOG"
    if [ $LIVE = 1 ]; then render; else step_line $CUR; fi
}

fail() {
    spin_stop
    if [ $CUR -lt 0 ]; then return 0; fi
    STATE[$CUR]=fail
    INFO[$CUR]="failed"
    FOOTER="${RED}Failed at step $((CUR + 1)) of ${#NAMES[@]}${RESET}"
    printf '=== %s: FAILED\n' "${NAMES[$CUR]}" >> "$LOG"
    if [ $LIVE = 1 ]; then render; else step_line $CUR; fi
}

die() {
    trap - ERR
    fail
    printf '\n%s%s%s\n' "$RED" "$*" "$RESET" >&2
    printf '%s\n' "$*" >> "$LOG"
    printf '\nFull log: %s\n' "$LOG" >&2
    exit 1
}

on_error() {
    trap - ERR
    fail
    if [ $VERBOSE = 0 ]; then
        printf '\n%sLast 30 lines of install.log:%s\n' "$DIM" "$RESET" >&2
        tail -n 30 "$LOG" | sed 's/^/  /' >&2
    fi
    printf '\nFull log: %s\n' "$LOG" >&2
}
trap on_error ERR
trap 'spin_stop; if [ $LIVE = 1 ]; then printf "\e[?25h"; fi' EXIT

# Run a build command: output to install.log (and to the screen with --verbose).
run() {
    if [ $VERBOSE = 1 ]; then
        "$@" 2>&1 | tee -a "$LOG"
    else
        "$@" >> "$LOG" 2>&1
    fi
}

# ", N warnings" for the log lines written by the current step (empty if none)
warnings() {
    local n
    n=$(tail -n +$((LOG_MARK + 1)) "$LOG" | grep -ci 'warning' || true)
    if [ "$n" = 1 ]; then printf ', 1 warning'
    elif [ "$n" -gt 1 ]; then printf ', %d warnings' "$n"; fi
}

# Timing: t0=$(now); ...; elapsed "$t0"  ->  "1 min 02.3 s" or "4.2 s"
now() { date +%s.%N; }
elapsed() {
    awk -v a="$1" -v b="$(now)" 'BEGIN {
        t = b - a; m = int(t / 60); s = t - 60 * m
        if (m) printf "%d min %04.1f s", m, s; else printf "%.1f s", s }'
}

# ---------------------------------------------------------------------------
# Steps
# ---------------------------------------------------------------------------
add_step S_ENV     "Environment"
add_step S_SRC     "Sources"
add_step S_CPU     "CPU target"
if [ $CLEAN = 1 ]; then add_step S_CLEAN "Clean"; fi
add_step S_INT     "gammcor-integrals"
add_step S_XCFUN   "XCFun"
add_step S_HDF5    "HDF5"
add_step S_TOOLS   "Python tools"
add_step S_GAMMCOR "GAMMCOR"
if [ -n "$PREFIX" ]; then add_step S_INSTALL "Install"; fi

T_TOTAL=$(now)

printf '\n%sGAMMCOR installation%s\n%s%s%s\n' "$BOLD" "$RESET" "$DIM" "$RULE" "$RESET"
if [ $VERBOSE = 0 ]; then
    printf '%sbuild output: install.log  (--verbose to show it here)%s\n' "$DIM" "$RESET"
fi
echo
if [ $LIVE = 1 ]; then printf '\e[?25l'; fi   # hide the cursor while redrawing

# ---------------------------------------------------------------------------
# 1. Environment checks
# ---------------------------------------------------------------------------
begin $S_ENV "Checking the environment"

command -v ifx > /dev/null || die "ifx not found.

GAMMCOR is compiled with the Intel oneAPI compilers. Load Intel oneAPI into
your environment.

On a workstation this may be:

  source /opt/intel/oneapi/setvars.sh

On an HPC cluster, load the appropriate Intel compiler module
(e.g. \`module avail intel\`)."

command -v icx > /dev/null || die "icx not found (Intel oneAPI C/C++ compiler), see above for ifx."

if [ -n "$MKLROOT" ]; then
    MKL="$MKLROOT"
elif [ -f "$CONDA_PREFIX/lib/libmkl_rt.so" ]; then
    MKL="$CONDA_PREFIX/lib/libmkl_rt.so"
else
    die "MKL not found: \$MKLROOT is not set and there is no \$CONDA_PREFIX/lib/libmkl_rt.so."
fi

for tool in python3 meson ninja make; do
    command -v $tool > /dev/null || die "$tool not found."
done
python3 -m pip --version > /dev/null 2>&1 || die "pip not found for $(command -v python3)."

ok "" "ifx:    $(ifx --version 2>&1 | sed -n 1p | awk '{print $3}')
MKL:    $MKL
meson:  $(meson --version)
ninja:  $(ninja --version)"

# ---------------------------------------------------------------------------
# 2. gammcor-integrals submodule
# ---------------------------------------------------------------------------
begin $S_SRC "Checking the gammcor-integrals submodule"

if [ ! -f gammcor-integrals/meson.build ]; then
    run git submodule update --init || die "gammcor-integrals is missing. Clone with: git clone --recursive"
fi

ok "gammcor-integrals $(git -C gammcor-integrals rev-parse --short HEAD 2> /dev/null || true)"

# ---------------------------------------------------------------------------
# 3. CPU target (one target for gammcor-integrals and gammcor)
# ---------------------------------------------------------------------------
begin $S_CPU "Detecting the CPU"

if [ "$CPU" = auto ]; then
    CPU=$(python3 gammcor-integrals/.meson/scripts/detect_cpu.py | sed -n 1p)
fi

ok "$CPU, $JOBS jobs"

# ---------------------------------------------------------------------------
# --clean: remove everything built by this script
# ---------------------------------------------------------------------------
if [ $CLEAN = 1 ]; then
    begin $S_CLEAN "Removing previous builds"
    run rm -rf build gammcor-integrals/build hdf5/lib hdf5/include
    run make -C xcfun clean
    ok "previous builds removed"
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
begin $S_INT "Building gammcor-integrals"

run meson_setup gammcor-integrals/build gammcor-integrals -Dcpu="$CPU"
run meson compile -C gammcor-integrals/build -j "$JOBS"

ok "$(elapsed "$T_STEP")$(warnings)"

# ---------------------------------------------------------------------------
# 5. XCFun  ->  xcfun/lib/libxcfun.a
# ---------------------------------------------------------------------------
begin $S_XCFUN "Building XCFun"

run make -C xcfun intel

ok "$(elapsed "$T_STEP")$(warnings)"

# ---------------------------------------------------------------------------
# 6. HDF5  ->  hdf5/lib, hdf5/include
# ---------------------------------------------------------------------------
begin $S_HDF5 "Building HDF5"

if [ -f hdf5/lib/libhdf5_fortran.a ]; then
    ok "already built"
else
    run ./hdf5/build_hdf5.sh
    ok "$(elapsed "$T_STEP")$(warnings)"
fi

# ---------------------------------------------------------------------------
# 7. gammcor_tools (Python), into the current python3
# ---------------------------------------------------------------------------
begin $S_TOOLS "Installing gammcor_tools"

if [ -z "$PREFIX" ]; then
    PIP_ARGS=(-e ./gammcor_tools)
else
    PIP_ARGS=(./gammcor_tools --no-deps --no-build-isolation)
fi
run python3 -m pip install "${PIP_ARGS[@]}" || die "pip install failed (details in install.log).
If pip says 'externally-managed-environment', use a venv, a conda env,
or run: python3 -m pip install --user ${PIP_ARGS[*]}"

ok "gammcor_tools, $(elapsed "$T_STEP")"

# ---------------------------------------------------------------------------
# 8. gammcor  ->  build/gammcor, bin/gammcor
# ---------------------------------------------------------------------------
begin $S_GAMMCOR "Building GAMMCOR"

MESON_ARGS=(-Dcpu="$CPU")
if [ -n "$PREFIX" ]; then MESON_ARGS+=(-Dprefix="$PREFIX"); fi

run meson_setup build . "${MESON_ARGS[@]}"
run meson compile -C build -j "$JOBS"

ok "$(elapsed "$T_STEP")$(warnings)"

# ---------------------------------------------------------------------------
# 9. Install (only with --prefix)
# ---------------------------------------------------------------------------
if [ -n "$PREFIX" ]; then
    begin $S_INSTALL "Installing into $PREFIX"
    run meson install -C build
    ok "$PREFIX/bin/gammcor"
fi

# ---------------------------------------------------------------------------
# Done
# ---------------------------------------------------------------------------
FOOTER="${GREEN}Done${RESET} in ${CYAN}$(elapsed "$T_TOTAL")${RESET}"
if [ $LIVE = 1 ]; then
    render
else
    printf '\n%s%s%s\n%s\n' "$DIM" "$RULE" "$RESET" "$FOOTER"
fi
printf '\n  gammcor  %s%s%s\n' "$CYAN" "$PWD/bin/gammcor" "$RESET"
if [ -n "$PREFIX" ]; then printf '           %s%s%s\n' "$CYAN" "$PREFIX/bin/gammcor" "$RESET"; fi
printf '  log      %s%s%s\n\n' "$DIM" "$LOG" "$RESET"
