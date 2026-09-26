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
# Keep the terminal on fds 3 and 4: the ERR trap can fire inside run(), where
# stdout and stderr go to install.log.
exec 3>&1 4>&2

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

# Ask the terminal for its background colour (OSC 11); return 0 if it is light.
bg_is_light() {
    [[ -t 0 && -t 1 ]] || return 1
    local resp r g b
    printf '\033]11;?\a' > /dev/tty
    # -t: a terminal that does not answer must not hang the script
    IFS= read -r -s -t 0.3 -d $'\a' resp < /dev/tty || true
    [[ $resp =~ rgb:([0-9a-fA-F]{2})[0-9a-fA-F]*/([0-9a-fA-F]{2})[0-9a-fA-F]*/([0-9a-fA-F]{2}) ]] || return 1
    r=$((16#${BASH_REMATCH[1]})); g=$((16#${BASH_REMATCH[2]})); b=$((16#${BASH_REMATCH[3]}))
    (( (299 * r + 587 * g + 114 * b) / 1000 > 128 ))
}

# Logo colours: lavender (darker on a light background, or GAMMCOR_THEME=light) and red
if [ -z "$RESET" ]; then
    LAV="" LOGO_RED=""
elif [ "${GAMMCOR_THEME:-}" = light ] || { [ -z "${GAMMCOR_THEME:-}" ] && bg_is_light; }; then
    LAV=$'\e[38;5;60m' LOGO_RED=$'\e[38;5;203m'
else
    LAV=$'\e[38;5;146m' LOGO_RED=$'\e[38;5;203m'
fi

print_logo() {
    local g=("╭───╴" "│ ──╮" "╰───╯")
    local a=("  ●  " " ╱ ╲ " "●───●")
    local m=("│╲  ╱│" "│ ╲╱ │" "│    │")
    local c=("╭───╴" "│    " "╰───╴")
    local o=("╭───╮" "│   │" "╰───╯")
    local r=("╭───╮" "├───╯" "│  ╲ ")
    local s=("█████" "█████" "█████")
    local i
    for i in 0 1 2; do
        printf '  %s%s %s%s %s%s %s   %s %s %s   %s%s%s\n' \
            "$LAV" "${g[i]}" "$LOGO_RED" "${a[i]}" "$LAV" "${m[i]}" "${m[i]}" \
            "${c[i]}" "${o[i]}" "${r[i]}" "$LOGO_RED" "${s[i]}" "$RESET"
    done
}

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
    exec 1>&3 2>&4
    fail
    printf '\n%s%s%s\n' "$RED" "$*" "$RESET" >&2
    printf '%s\n' "$*" >> "$LOG"
    printf '\nFull log: %s\n' "$LOG" >&2
    exit 1
}

on_error() {
    trap - ERR
    exec 1>&3 2>&4
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

echo
print_logo
printf '%s%s%s\n' "$DIM" "$RULE" "$RESET"
if [ $VERBOSE = 0 ]; then
    printf '%sbuild output: install.log%s\n' "$DIM" "$RESET"
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
# conda-build/rattler-build set $PYTHON to the host env's python; python3 on
# PATH can be the build env's one there, and gammcor_tools must not go in it.
PYTHON="${PYTHON:-python3}"
"$PYTHON" -m pip --version > /dev/null 2>&1 || die "pip not found for $(command -v "$PYTHON")."

ok "" "ifx:    $(ifx --version 2>&1 | sed -n 1p | awk '{print $3}')
MKL:    $MKL
meson:  $(meson --version)
ninja:  $(ninja --version)"

# ---------------------------------------------------------------------------
# 2. gammcor-integrals submodule
# ---------------------------------------------------------------------------
begin $S_SRC "Checking the gammcor-integrals submodule"

# A clone without --recursive leaves gammcor-integrals/ empty; fetch it here.
# An archive download (zip/tar.gz from GitLab) has no .git and cannot fetch it.
GAMMCOR_URL=https://gitlab.com/qchem/gammcor.git
if [ ! -f gammcor-integrals/meson.build ]; then
    [ -e .git ] || die "Git submodule gammcor-integrals is missing.

This copy of GAMMCOR is not a git repository (downloaded as an archive?),
so the submodule cannot be fetched. Clone the repository with its submodules:

  git clone --recursive $GAMMCOR_URL"

    RUNLABEL="Fetching the gammcor-integrals submodule"
    if [ $LIVE = 1 ]; then spin_stop; render; spin_start; fi
    # GIT_TERMINAL_PROMPT=0: fail instead of waiting, invisibly, for a password
    GIT_TERMINAL_PROMPT=0 run git submodule update --init --recursive || true
    [ -f gammcor-integrals/meson.build ] || die "Git submodule gammcor-integrals is missing and could not be fetched
(git output is in install.log).

The repository was probably cloned without --recursive. Fetch the submodule
inside this directory:

  git submodule update --init --recursive

or clone again, with submodules:

  git clone --recursive $GAMMCOR_URL"
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

# Extra link flags for a conda environment (ifx and MKL from conda-forge):
#  - conda-forge MKL (libmkl_intel_thread.so, 2026.1) calls __atomic_load and
#    __atomic_compare_exchange from libatomic, which ifx -qmkl does not link;
#  - the MKL, OpenMP and Intel runtime libraries are in <conda prefix>/lib,
#    which is not on the loader path. conda puts -rpath in $LDFLAGS, but Meson
#    ignores $LDFLAGS when the profile sets fortran_link_args.
# Without them Meson's sanity check fails ("ifx cannot compile programs",
# "executables ... are not runnable").
# Meson links with --as-needed and places -latomic before the MKL libraries
# that -qmkl appends, so the linker would drop libatomic as unused. Linking it
# with --no-as-needed keeps it; one -Wl argument also stops Meson moving it.
# The conda prefix is the one holding MKL: $MKLROOT if set (the conda recipe
# sets MKLROOT=$PREFIX, the host env), else the active env. A oneAPI MKLROOT
# has no libatomic.so next to MKL, so this adds nothing there.
CONDA_LINK_ARGS=()
CONDA_LIB="${MKLROOT:-$CONDA_PREFIX}/lib"
if [ -f "$CONDA_LIB/libmkl_intel_thread.so" ] && [ -f "$CONDA_LIB/libatomic.so" ]; then
    CONDA_LINK_ARGS=(-Wl,-rpath,"$CONDA_LIB" -Wl,--push-state,--no-as-needed,-latomic,--pop-state)
fi

# Print a native file that repeats fortran_link_args of <profile> with <args> appended.
#   link_args_native_file <profile> <args>...
link_args_native_file() {
    python3 - "$@" <<'EOF'
import ast, configparser, sys
ini = configparser.ConfigParser(interpolation=None)
ini.read(sys.argv[1])
args = ast.literal_eval(ini['built-in options']['fortran_link_args'])
print('[built-in options]')
print(f"fortran_link_args = {args + sys.argv[2:]!r}")
EOF
}

# Configure a Meson build directory, or update its options if it exists.
# A directory left behind by a failed setup has meson-private/ but no
# coredata.dat; meson setup starts it over.
#   meson_setup <build dir> <source dir> [options...]
meson_setup() {
    local dir=$1 src=$2; shift 2
    local profile="$src/.meson/profiles/ifx.ini"
    if [ -f "$dir/meson-private/coredata.dat" ]; then
        meson configure "$dir" "$@"
    else
        local native=(--native-file "$profile")
        if [ ${#CONDA_LINK_ARGS[@]} -gt 0 ]; then
            mkdir -p "$dir"
            link_args_native_file "$profile" "${CONDA_LINK_ARGS[@]}" > "$dir/conda.ini"
            native+=(--native-file "$PWD/$dir/conda.ini")
        fi
        meson setup "$dir" "$src" "${native[@]}" "$@"
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
run "$PYTHON" -m pip install "${PIP_ARGS[@]}" || die "pip install failed (details in install.log).
If pip says 'externally-managed-environment', use a venv, a conda env,
or run: $PYTHON -m pip install --user ${PIP_ARGS[*]}"

ok "gammcor_tools, $(elapsed "$T_STEP")"

# ---------------------------------------------------------------------------
# 8. gammcor  ->  build/gammcor, bin/gammcor
# ---------------------------------------------------------------------------
begin $S_GAMMCOR "Building GAMMCOR"

mkdir -p OBJ   # object directory of the Makefile build, for those who use make

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
printf '%sNo electrons were harmed during compilation.%s\n' "$LAV" "$RESET"
printf '\n  gammcor  %s%s%s\n' "$CYAN" "$PWD/bin/gammcor" "$RESET"
if [ -n "$PREFIX" ]; then printf '           %s%s%s\n' "$CYAN" "$PREFIX/bin/gammcor" "$RESET"; fi
printf '  log      %s%s%s\n\n' "$DIM" "$LOG" "$RESET"
