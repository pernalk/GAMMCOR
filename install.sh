#!/bin/bash
# Builds GAMMCOR (dSRS branch) from source: TREXIO (if not found), XCFun,
# gammcor. Steps that are already done with the same toolchain are skipped, so
# running it again is cheap.
#
#   ./install.sh [--prefix DIR] [--cpu TARGET] [--jobs N] [--trexio DIR]
#                [--command NAME] [--clean] [--verbose]
#
#   --prefix DIR    also install gammcor into DIR/bin
#   --cpu TARGET    auto (default) | native | avx512 | avx512-amd | avx2 | dispatch
#   --jobs N        parallel build jobs (default: nproc)
#   --trexio DIR    use the TREXIO installed in DIR (DIR/include/trexio_f.f90,
#                   DIR/lib/libtrexio.so); default: $TREXIO_ROOT, see below
#   --command NAME  name of the command made at the end (default gammcor-dsrs,
#                   none: make no command), see below
#   --clean         rebuild everything
#   --verbose       show compiler output on screen (it always goes to install.log)
#
# Toolchain: ifx, icx, MKL and Python come from one place, never from a mix.
#   conda     the active conda env has ifx and MKL ($CONDA_PREFIX/bin/ifx and
#             $CONDA_PREFIX/lib/libmkl_rt.so), e.g. an env created from
#             environment.yml. Everything is taken from the env: $MKLROOT is
#             ignored, and directories outside the env are dropped from
#             LD_LIBRARY_PATH, LIBRARY_PATH and CPATH (for this script only).
#             ifx, icx and python3 must resolve to $CONDA_PREFIX/bin: do not
#             load compiler or MKL modules after `conda activate`. gammcor is
#             linked with an RPATH to the env, so it keeps using the env's MKL
#             even when such modules are loaded later.
#   modules   anything else (no conda env, or one without ifx, such as base):
#             ifx and icx from PATH, MKL from $MKLROOT (HPC modules or oneAPI
#             setvars.sh).
#
# TREXIO >= 2.5.0 (the first with transition RDMs), the first of these that
# has it:
#   1. --trexio DIR, or $TREXIO_ROOT (an older TREXIO there is an error)
#   2. conda mode: the conda env (conda install -c conda-forge trexio)
#      modules mode: pkg-config (e.g. installed into ~/.local as in the manual,
#      with PKG_CONFIG_PATH set)
#      An older TREXIO here is skipped, with a note in install.log.
#   3. trexio/install, built by an earlier run of this script
#   4. otherwise it is built from source into trexio/install: the release
#      tarball trexio/trexio-<version>.tar.gz is downloaded from GitHub, or
#      put there by hand on a machine without internet. Always with HDF5
#      (TREXIO files are HDF5 files): conda mode: the env's HDF5; modules
#      mode: the system's (pkg-config hdf5 or h5cc). When there is none,
#      HDF5 is built first into hdf5/ by hdf5/build_hdf5.sh.
#
# install.log has every command that was run, quoted so that it can be copied
# and run by hand; a failed one is also shown on screen.
#
# The toolchain (mode, ifx and its version, icx, MKL, python3, the GCC used by
# icx) is stored in
# .toolchain after a successful run. When it differs on the next run, the
# changed fields are shown and everything is rebuilt, as with --clean, so
# objects and .mod files from different compilers are never mixed.
#
# Without --prefix, the end of the script makes `gammcor-dsrs` (--command) a
# command: a symlink to bin/gammcor in $CONDA_PREFIX/bin (conda mode: there
# while the env is active) or in ~/.local/bin (modules mode). It is not called
# `gammcor`, so that it does not replace the command of another GAMMCOR
# branch installed into the same env. A rebuild needs no new link. A file of
# that name there that is not a symlink is left alone.
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
TREXIO_ROOT="${TREXIO_ROOT:-}"
COMMAND=gammcor-dsrs

while [ $# -gt 0 ]; do
    case "$1" in
        --prefix) PREFIX="$(realpath -m "$2")"; shift 2 ;;
        --cpu)    CPU="$2";  shift 2 ;;
        --jobs)   JOBS="$2"; shift 2 ;;
        --trexio) TREXIO_ROOT="$(realpath -m "$2")"; shift 2 ;;
        --command) COMMAND="$2"; shift 2 ;;
        --clean)  CLEAN=1;   shift ;;
        -v|--verbose) VERBOSE=1; shift ;;
        -h|--help) awk 'NR > 1 { if (!/^#/) exit; print }' "$0"; exit 0 ;;
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
STEP_VARS=()    # names of the index variables, for insert_step
add_step() {   # add_step <index variable> <name>
    printf -v "$1" %d ${#NAMES[@]}
    NAMES+=("$2"); STATE+=(pending); INFO+=(""); EXTRA+=("")
    STEP_VARS+=("$1")
}

# Add a step at <position> once the list is on screen; the steps from there
# on move down by one.
insert_step() {   # insert_step <index variable> <name> <position>
    local pos=$3 v
    for v in "${STEP_VARS[@]}"; do
        if [ "${!v}" -ge "$pos" ]; then printf -v "$v" %d $((${!v} + 1)); fi
    done
    NAMES=("${NAMES[@]:0:pos}" "$2" "${NAMES[@]:pos}")
    STATE=("${STATE[@]:0:pos}" pending "${STATE[@]:pos}")
    INFO=("${INFO[@]:0:pos}" "" "${INFO[@]:pos}")
    EXTRA=("${EXTRA[@]:0:pos}" "" "${EXTRA[@]:pos}")
    printf -v "$1" %d "$pos"
    STEP_VARS+=("$1")
}

CUR=-1          # running step
RUNLABEL=""     # its label, e.g. "Building XCFun"
LAST_COMMAND="" # command of run() that failed, see run()
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
    LAST_COMMAND=""
    LOG_MARK=$(wc -l < "$LOG")
    T_STEP=$(now)
    STEP_SECONDS=$SECONDS
    printf '\n=== [%s] Step %d of %d: %s\n' "$(date '+%F %T')" \
        $((CUR + 1)) ${#NAMES[@]} "$RUNLABEL" >> "$LOG"
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
    printf '=== [%s] %s: ok %s\n' "$(date '+%F %T')" "${NAMES[$CUR]}" "$1" >> "$LOG"
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

# The command of run() that failed, to copy and run by hand (nothing if none)
show_failed_command() {
    if [ -n "$LAST_COMMAND" ]; then
        printf '\n%sFailed command:%s\n  %s\n' "$DIM" "$RESET" "$LAST_COMMAND" >&2
    fi
}

die() {
    trap - ERR
    exec 1>&3 2>&4
    fail
    printf '\n%s%s%s\n' "$RED" "$*" "$RESET" >&2
    printf '%s\n' "$*" >> "$LOG"
    show_failed_command
    printf '\nFull log: %s\n' "$LOG" >&2
    exit 1
}

on_error() {
    local status=$?
    trap - ERR
    exec 1>&3 2>&4
    printf '\n[%s] FAILED with exit status %d%s\n' "$(date '+%F %T')" $status \
        "${LAST_COMMAND:+: $LAST_COMMAND}" >> "$LOG"
    fail
    if [ $VERBOSE = 0 ]; then
        printf '\n%sLast 30 lines of install.log:%s\n' "$DIM" "$RESET" >&2
        tail -n 30 "$LOG" | sed 's/^/  /' >&2
    fi
    show_failed_command
    printf '\nFull log: %s\n' "$LOG" >&2
}
trap on_error ERR
trap 'spin_stop; if [ $LIVE = 1 ]; then printf "\e[?25h"; fi' EXIT

# Run a build command: output to install.log (and to the screen with --verbose).
# The command goes to install.log first, quoted, so that it can be copied from
# there. A failure fires the ERR trap inside run() (unless run is followed by
# ||), with LAST_COMMAND still naming the command; a success clears it.
run() {
    local arg quoted=()
    for arg in "$@"; do quoted+=("$(printf '%q' "$arg")"); done
    LAST_COMMAND="${quoted[*]}"
    printf '\n[%s] $ %s\n' "$(date '+%F %T')" "$LAST_COMMAND" >> "$LOG"
    if [ $VERBOSE = 1 ]; then
        "$@" 2>&1 | tee -a "$LOG"
    else
        "$@" >> "$LOG" 2>&1
    fi
    local status=$?
    if [ $status = 0 ]; then LAST_COMMAND=""; fi
    return $status
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
# Toolchain helpers
# ---------------------------------------------------------------------------
NL=$'\n'
STAMP=.toolchain                                   # toolchain of the last successful run
FP_KEYS=(mode ifx ifx_version icx mkl python gcc)  # its fields, in this order

# Keep only the entries of a :-separated path variable that are inside
# $CONDA_PREFIX; the others go into DROPPED.
keep_conda_paths() {   # keep_conda_paths <variable>
    local e entries kept=()
    DROPPED=()
    IFS=: read -ra entries <<< "${!1:-}"
    for e in "${entries[@]}"; do
        case "$e" in
            "") ;;
            "$CONDA_PREFIX" | "$CONDA_PREFIX"/*) kept+=("$e") ;;
            *) DROPPED+=("$e") ;;
        esac
    done
    if [ ${#kept[@]} -gt 0 ]; then
        export "$1=$(IFS=:; echo "${kept[*]}")"
    else
        unset "$1"
    fi
}

# MKL version from <include dir>/mkl_version.h, e.g. 2025.0.1 (empty if unknown)
mkl_version() {   # mkl_version <include dir>
    [ -f "$1/mkl_version.h" ] || return 0
    awk '$1 == "#define" { v[$2] = $3 }
         END { if ("__INTEL_MKL__" in v)
                   printf "%s.%s.%s", v["__INTEL_MKL__"], v["__INTEL_MKL_MINOR__"], v["__INTEL_MKL_UPDATE__"] }' \
        "$1/mkl_version.h" || true
}

# The version at the start of the first word of `<program> --version` that
# starts with one, e.g. 2025.0.4, or 1.13.2 from 1.13.2.git.kitware.jobserver-1
# (ninja from PyPI)
version_of() {   # version_of <program>
    "$1" --version 2>&1 | awk '{
        for (i = 1; i <= NF; i++)
            if (match($i, /^[0-9]+\.[0-9]+(\.[0-9]+)*/)) { print substr($i, 1, RLENGTH); exit } }' || true
}

# Are there artifacts of an earlier build?
built_before() {
    [ -e build ] || [ -e trexio/install ] || [ -e hdf5/lib ] \
        || compgen -G 'xcfun/lib/libxcfun.*' > /dev/null
}

# The current toolchain, in the format of the stamp file
fingerprint() {
    local key
    for key in "${FP_KEYS[@]}"; do printf '%s=%s\n' "$key" "${FP[$key]}"; done
}

# TREXIO (see the header of this file)
TREXIO_VERSION=2.6.1
TREXIO_TARBALL=trexio/trexio-$TREXIO_VERSION.tar.gz
TREXIO_URL=https://github.com/TREX-CoE/trexio/releases/download/v$TREXIO_VERSION/trexio-$TREXIO_VERSION.tar.gz
TREXIO_LOCAL="$PWD/trexio/install"
TREXIO_MIN=2.5.0

# Does <prefix> have a TREXIO (Fortran interface and shared library)?
has_trexio() {   # has_trexio <prefix>
    [ -n "$1" ] && [ -f "$1/include/trexio_f.f90" ] \
        && { [ -e "$1/lib/libtrexio.so" ] || [ -e "$1/lib64/libtrexio.so" ]; }
}

# Is it new enough? TREXIO_MIN added the transition RDMs read by
# sapt_interface.F90; checked by function, not by version number.
trexio_new_enough() {   # trexio_new_enough <prefix>
    grep -qi 'trexio_read_rdm_2e_transition_size' "$1/include/trexio_f.f90"
}

# A TREXIO that gammcor can use; an older one found <where> goes into
# TREXIO_SKIPPED.
usable_trexio() {   # usable_trexio <prefix> <where>
    has_trexio "$1" || return 1
    trexio_new_enough "$1" && return 0
    TREXIO_SKIPPED+="${TREXIO_SKIPPED:+$NL}skipped TREXIO $(trexio_version "$1") < $TREXIO_MIN ($2): $1"
    return 1
}

# TREXIO version from <prefix>/include/trexio.h, e.g. 2.6.1 (empty if unknown)
trexio_version() {   # trexio_version <prefix>
    sed -n 's/^#define TREXIO_PACKAGE_VERSION "\(.*\)"/\1/p' "$1/include/trexio.h" 2> /dev/null || true
}

# Download the TREXIO release tarball with curl or wget
download_trexio() {
    if command -v curl > /dev/null; then
        run curl -fL --retry 2 -o "$TREXIO_TARBALL.part" "$TREXIO_URL"
    elif command -v wget > /dev/null; then
        run wget -O "$TREXIO_TARBALL.part" "$TREXIO_URL"
    else
        return 1
    fi
}

# ---------------------------------------------------------------------------
# Steps
# ---------------------------------------------------------------------------
add_step S_ENV     "Environment"
add_step S_CPU     "CPU target"
if [ $CLEAN = 1 ]; then add_step S_CLEAN "Clean"; fi
add_step S_HDF5    "HDF5"
add_step S_TREXIO  "TREXIO"
add_step S_XCFUN   "XCFun"
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

# Toolchain mode, decided here once (see the header of this file)
if [ -n "$CONDA_PREFIX" ] && [ -x "$CONDA_PREFIX/bin/ifx" ] \
        && [ -f "$CONDA_PREFIX/lib/libmkl_rt.so" ]; then
    MODE=conda
    MODE_INFO="conda ($CONDA_PREFIX)"
else
    MODE=modules
    MODE_INFO="modules"
fi

PATHS_DROPPED=""
if [ $MODE = conda ]; then
    # -qmkl links the MKL in $MKLROOT when it is set (an MKL module sets it);
    # without it, -qmkl finds the MKL next to ifx, i.e. the env's.
    unset MKLROOT
    MKL="$CONDA_PREFIX/lib/libmkl_rt.so"
    MKL_VERSION=$(mkl_version "$CONDA_PREFIX/include")

    # A module loaded after `conda activate` puts its bin/ before the env's.
    SHADOWED=""
    for tool in ifx icx python3; do
        path=$(command -v $tool || true)
        if [ -z "$path" ] || [ "$(realpath "$(dirname "$path")")" != "$(realpath "$CONDA_PREFIX/bin")" ]; then
            SHADOWED+="$NL  $(printf '%-8s %s' "$tool" "${path:-not found}")"
        fi
    done
    [ -z "$SHADOWED" ] || die "The active conda env is $CONDA_PREFIX,
but these tools do not come from its bin/:
$SHADOWED

A module loaded after \`conda activate\` is shadowing the env, and install.sh
does not mix a conda env with modules. Unload that module (\`module list\`
shows what is loaded):

  module unload <module>

or start again from a clean shell:

  module purge
  conda activate ${CONDA_DEFAULT_ENV:-$CONDA_PREFIX}
  ./install.sh

A tool that is 'not found' is missing from the env: update it from
environment.yml (conda env update -f environment.yml)."

    # Library and header directories of modules must not reach the compilers,
    # the linker, or the test programs that Meson and configure run. This
    # changes the environment of install.sh only, not of your shell.
    for var in LD_LIBRARY_PATH LIBRARY_PATH CPATH; do
        keep_conda_paths $var
        if [ ${#DROPPED[@]} -gt 0 ]; then
            printf 'dropped from %s: %s\n' "$var" "${DROPPED[*]}" >> "$LOG"
            PATHS_DROPPED+="${PATHS_DROPPED:+, }$var (${#DROPPED[@]})"
        fi
    done
else
    command -v ifx > /dev/null || die "ifx not found.

GAMMCOR is compiled with the Intel oneAPI compilers. Load Intel oneAPI into
your environment.

On a workstation this may be:

  source /opt/intel/oneapi/setvars.sh

On an HPC cluster, load the appropriate Intel compiler module
(e.g. \`module avail intel\`).

Or build with conda instead: create the env from environment.yml, activate
it (with no Intel modules loaded) and run ./install.sh again."

    command -v icx > /dev/null || die "icx not found (Intel oneAPI C/C++ compiler), see above for ifx."

    [ -n "$MKLROOT" ] || die "MKL not found: \$MKLROOT is not set.

Load the MKL module that goes with the Intel compiler (e.g. \`module avail mkl\`),
or source /opt/intel/oneapi/setvars.sh."
    MKL="$MKLROOT"
    MKL_VERSION=$(mkl_version "$MKLROOT/include")

    # A conda env with ifx but without MKL: its ifx would be used with the MKL
    # of a module. An MKLROOT that is itself a conda prefix (the conda recipe
    # sets MKLROOT=$PREFIX, the host env) is not a module, so that is fine.
    if [ -n "$CONDA_PREFIX" ] && [ ! -d "$MKLROOT/conda-meta" ]; then
        for tool in ifx icx; do
            case "$(command -v $tool)" in
                "$CONDA_PREFIX"/*) die "$tool comes from the conda env $CONDA_PREFIX,
but that env has no MKL ($CONDA_PREFIX/lib/libmkl_rt.so is missing),
so MKL would come from \$MKLROOT=$MKLROOT.
install.sh does not mix a conda env with modules. Either install MKL into the env:

  conda install -c conda-forge mkl-devel

or deactivate the env and use the modules only:

  conda deactivate" ;;
            esac
        done
    fi
fi

for tool in python3 meson ninja make tar; do
    command -v $tool > /dev/null || die "$tool not found."
done
PYTHON=$(command -v python3)   # conda mode: checked above, in $CONDA_PREFIX/bin

# The GCC whose C++ headers and libstdc++ icx uses for XCFun (passed to
# xcfun/Makefile.common as GCC_INSTALL_DIR). The Makefile's default is the one
# of g++ on PATH, which in conda mode is the system's: conda names its g++
# e.g. x86_64-conda-linux-gnu-g++. Empty (no g++): icx chooses.
if [ $MODE = conda ]; then
    GXX=$(compgen -G "$CONDA_PREFIX/bin/*-g++" | sed -n 1p || true)
else
    GXX=$(command -v g++ || true)
fi
GCC_INSTALL_DIR=""
if [ -n "$GXX" ]; then
    GCC_INSTALL_DIR=$(dirname "$("$GXX" -print-libgcc-file-name)")
fi

# ifx by its absolute path: it goes into the Meson build directories (see
# toolchain_native_file), so a later `meson compile -C build` does not take
# another ifx from PATH.
IFX=$(command -v ifx)
ICX=$(command -v icx)

# Compare the toolchain with the one of the last successful run
declare -A FP=(
    [mode]=$MODE
    [ifx]=$(realpath "$IFX")
    [ifx_version]=$(ifx --version 2>&1 | sed -n 1p)
    [icx]=$(realpath "$ICX")
    [mkl]="$MKL${MKL_VERSION:+ $MKL_VERSION}"
    [python]=$(realpath "$PYTHON")
    [gcc]=$GCC_INSTALL_DIR
)
CLEAN_REASON=""
STAMP_CHANGES=""
if [ -f "$STAMP" ]; then
    declare -A OLD=()
    while IFS='=' read -r key value; do
        case "$key" in "" | "#"*) continue ;; esac
        OLD[$key]=$value
    done < "$STAMP"
    for key in "${FP_KEYS[@]}"; do
        if [ "${OLD[$key]-}" != "${FP[$key]}" ]; then
            STAMP_CHANGES+="${STAMP_CHANGES:+$NL}$key: ${OLD[$key]:-(none)} → ${FP[$key]}"
        fi
    done
    if [ -n "$STAMP_CHANGES" ]; then CLEAN_REASON="toolchain changed"; fi
elif built_before; then
    # built before the stamp existed, or by a run that did not finish
    CLEAN_REASON="no $STAMP stamp"
fi
fingerprint >> "$LOG"

# The environment the build tools see (in conda mode after the clean-up above)
{
    echo
    printf '%-16s %s\n' "PATH:" "$PATH"
    for var in LD_LIBRARY_PATH LIBRARY_PATH CPATH C_INCLUDE_PATH CPLUS_INCLUDE_PATH \
            PYTHONPATH PYTHONHOME MKLROOT CC CXX FC CFLAGS FFLAGS LDFLAGS; do
        printf '%-16s %s\n' "$var:" "${!var:-<unset>}"
    done
    if [ $MODE = modules ]; then
        printf '%-16s %s\n' "modules:" "${LOADEDMODULES:-<not recorded by the module system>}"
    fi
    echo
} >> "$LOG"

# env_table <width>: lays out the "name: version path" rows read from stdin.
# When a row is wider than <width>, every path goes on a line of its own, so
# that a long path does not push the tool name off the screen:
#   ifx:     2026.1.1
#      └─ /long/path/to/envs/gammcor-dev/bin/ifx
env_table() {
    local w=$1 name ver path i wide=0 rows=()
    while read -r name ver path; do
        rows+=("$name" "$ver" "$path")
        if [ $((20 + ${#path})) -gt "$w" ]; then wide=1; fi
    done
    for ((i = 0; i < ${#rows[@]}; i += 3)); do
        if [ $wide = 0 ]; then
            printf '%-8s %-10s %s\n' "${rows[i]}" "${rows[i+1]}" "${rows[i+2]}"
        else
            printf '%-8s %s\n   └─ %s\n' "${rows[i]}" "${rows[i+1]}" "${rows[i+2]}"
        fi
    done
}

# name: version, path
ENV_ROWS=$(
    v() { local s; s=$(version_of "$1"); printf '%s' "${s:-?}"; }
    printf '%-8s %-10s %s\n' \
        "ifx:" "$(v "$IFX")" "$IFX" \
        "icx:" "$(v "$ICX")" "$ICX" \
        "MKL:" "${MKL_VERSION:-?}" "$MKL"
    if [ -n "$GXX" ]; then printf '%-8s %-10s %s\n' "g++:" "$(v "$GXX")" "$GXX"; fi
    printf '%-8s %-10s %s\n' "python:" "$(v "$PYTHON")" "$PYTHON"
    printf '%-8s %-10s %s\n' \
        "meson:" "$(v meson)" "$(command -v meson)" \
        "ninja:" "$(v ninja)" "$(command -v ninja)"
)
ENV_NOTES=$(
    if [ -n "$PATHS_DROPPED" ]; then printf 'non-conda entries dropped from %s\n' "$PATHS_DROPPED"; fi
)
printf '%s\n' "$ENV_ROWS" ${ENV_NOTES:+"$ENV_NOTES"} >> "$LOG"
ENV_SUMMARY=$(env_table $((COLS - 6)) <<< "$ENV_ROWS"; printf '%s' "$ENV_NOTES")

# --trexio DIR: checked here, not first after HDF5 is built
if [ -n "$TREXIO_ROOT" ]; then
    has_trexio "$TREXIO_ROOT" || die "No TREXIO in $TREXIO_ROOT (from --trexio or \$TREXIO_ROOT):
$TREXIO_ROOT/include/trexio_f.f90 or $TREXIO_ROOT/lib/libtrexio.so is missing."
    trexio_new_enough "$TREXIO_ROOT" || die "The TREXIO in $TREXIO_ROOT (from --trexio or \$TREXIO_ROOT)
is version $(trexio_version "$TREXIO_ROOT"), but GAMMCOR needs TREXIO >= $TREXIO_MIN (transition RDMs).
Without --trexio, install.sh builds TREXIO $TREXIO_VERSION from source."
fi

ok "" "$ENV_SUMMARY"

# Artifacts of another toolchain are removed before anything is built: the
# skip logic below (TREXIO already built, meson configure of an existing build
# directory, incremental XCFun make) is only safe with the same toolchain.
if [ -n "$CLEAN_REASON" ] && [ $CLEAN = 0 ]; then
    CLEAN=1
    insert_step S_CLEAN "Clean" $((S_CPU + 1))
fi

# ---------------------------------------------------------------------------
# 2. CPU target
# ---------------------------------------------------------------------------
begin $S_CPU "Detecting the CPU"

if [ "$CPU" = auto ]; then
    CPU=$(python3 .meson/scripts/detect_cpu.py | sed -n 1p)
fi

if [ "$JOBS" = 1 ]; then CORES="1 core"; else CORES="$JOBS cores"; fi
ok "$CPU, compiling on $CORES in parallel"

# ---------------------------------------------------------------------------
# --clean: remove everything built by this script
# ---------------------------------------------------------------------------
if [ $CLEAN = 1 ]; then
    begin $S_CLEAN "Removing previous builds"
    if [ -n "$STAMP_CHANGES" ]; then
        printf 'toolchain changed since the last build:\n%s\n' "$STAMP_CHANGES" >> "$LOG"
    fi
    # the stamp goes too: it is written again only when this run succeeds.
    # A downloaded TREXIO tarball is kept.
    run rm -rf build trexio/build trexio/install hdf5/lib hdf5/include "$STAMP"
    run make -C xcfun clean
    ok "${CLEAN_REASON:+$CLEAN_REASON, }previous builds removed" "$STAMP_CHANGES"
fi

# Extra link flags for MKL from conda-forge (with ifx from conda-forge):
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
# The conda prefix is the one holding MKL:
#  - conda mode: the active env. --disable-new-dtags writes the rpath as
#    RPATH, not RUNPATH: the loader searches LD_LIBRARY_PATH before RUNPATH
#    but after RPATH, so gammcor keeps using the env's MKL and Intel runtime
#    when Intel/MKL modules are loaded later.
#  - modules mode: $MKLROOT, when it is a conda prefix (the conda recipe sets
#    MKLROOT=$PREFIX, the host env). A oneAPI MKLROOT has no libatomic.so next
#    to MKL, so this adds nothing there.
ATOMIC_ARGS=(-Wl,--push-state,--no-as-needed,-latomic,--pop-state)
CONDA_LINK_ARGS=()
if [ $MODE = conda ]; then
    CONDA_LINK_ARGS=(-Wl,--disable-new-dtags -Wl,-rpath,"$CONDA_PREFIX/lib")
    if [ -f "$CONDA_PREFIX/lib/libatomic.so" ]; then CONDA_LINK_ARGS+=("${ATOMIC_ARGS[@]}"); fi
elif [ -f "$MKLROOT/lib/libmkl_intel_thread.so" ] && [ -f "$MKLROOT/lib/libatomic.so" ]; then
    CONDA_LINK_ARGS=(-Wl,-rpath,"$MKLROOT/lib" "${ATOMIC_ARGS[@]}")
fi

# Print a native file, read after <profile>, that sets fortran to <ifx> (an
# absolute path, kept by the build directory) and, when <args> are given,
# repeats fortran_link_args of <profile> with <args> appended.
#   toolchain_native_file <profile> <ifx> [<args>...]
toolchain_native_file() {
    python3 - "$@" <<'EOF'
import ast, configparser, sys
profile, ifx, extra = sys.argv[1], sys.argv[2], sys.argv[3:]
print('[binaries]')
print(f"fortran = {ifx!r}")
if extra:
    ini = configparser.ConfigParser(interpolation=None)
    ini.read(profile)
    args = ast.literal_eval(ini['built-in options']['fortran_link_args'])
    print('[built-in options]')
    print(f"fortran_link_args = {args + extra!r}")
EOF
}

# Configure a Meson build directory, or update its options if it exists
# (it was built with the same toolchain: the Clean step removes it otherwise).
# A directory left behind by a failed setup has meson-private/ but no
# coredata.dat; meson setup starts it over.
#   meson_setup <build dir> <source dir> [options...]
meson_setup() {
    local dir=$1 src=$2; shift 2
    local profile="$src/.meson/profiles/ifx.ini"
    if [ -f "$dir/meson-private/coredata.dat" ]; then
        run meson configure "$dir" "$@"
    else
        mkdir -p "$dir"
        toolchain_native_file "$profile" "$IFX" "${CONDA_LINK_ARGS[@]}" \
            > "$dir/toolchain.ini" 2>> "$LOG"
        run meson setup "$dir" "$src" --native-file "$profile" \
            --native-file "$PWD/$dir/toolchain.ini" "$@"
    fi
}

# ---------------------------------------------------------------------------
# 3. HDF5  ->  $HDF5, for a TREXIO built from source (see the header of this
#    file): nothing when TREXIO is found, else the env's or the system's HDF5,
#    else hdf5/lib, hdf5/include
# ---------------------------------------------------------------------------
begin $S_HDF5 "Looking for TREXIO and HDF5"

TREXIO=""
TREXIO_SKIPPED=""
if [ -n "$TREXIO_ROOT" ]; then   # checked in step 1
    TREXIO=$TREXIO_ROOT
    TREXIO_FROM="--trexio"
elif [ $MODE = conda ] && usable_trexio "$CONDA_PREFIX" "conda env"; then
    TREXIO=$CONDA_PREFIX
    TREXIO_FROM="conda env"
elif [ $MODE = modules ] && PC_PREFIX=$(pkg-config --variable=prefix trexio 2> /dev/null) \
        && usable_trexio "$PC_PREFIX" "pkg-config"; then
    TREXIO=$PC_PREFIX
    TREXIO_FROM="pkg-config"
elif has_trexio "$TREXIO_LOCAL" \
        && ! ldd "$TREXIO_LOCAL/lib/libtrexio.so" 2>&1 | grep -q 'not found'; then
    # same toolchain, see the Clean step; its HDF5 must still be there
    TREXIO=$TREXIO_LOCAL
    TREXIO_FROM="already built"
fi

# HDF5: --with-hdf5 prefix for TREXIO's configure (empty: configure finds
# the system's HDF5 itself, with pkg-config hdf5 or h5cc, as here);
# HDF5_LIBDIR: where its shared libhdf5 is, for the RPATH of libtrexio.so
HDF5=""
HDF5_LIBDIR=""
if [ -n "$TREXIO" ]; then
    ok "not needed (TREXIO from $TREXIO_FROM)"
elif [ $MODE = conda ] && [ -f "$CONDA_PREFIX/include/hdf5.h" ]; then
    HDF5=$CONDA_PREFIX
    HDF5_LIBDIR=$CONDA_PREFIX/lib
    HDF5_FROM="conda env"
    ok "$HDF5_FROM" "$HDF5"
elif [ $MODE = modules ] && HDF5_PC=$(pkg-config --modversion 'hdf5 >= 1.8' 2> /dev/null); then
    HDF5_LIBDIR=$(pkg-config --variable=libdir hdf5)
    HDF5_FROM="pkg-config"
    ok "$HDF5_PC, $HDF5_FROM" "$HDF5_LIBDIR"
elif [ $MODE = modules ] && command -v h5cc > /dev/null; then
    # the -L of `h5cc -show`, e.g. icx -I<prefix>/include -L<prefix>/lib -lhdf5
    HDF5_LIBDIR=$(h5cc -show | tr ' ' '\n' | sed -n 's/^-L//p' | sed -n 1p)
    HDF5_FROM="h5cc"
    ok "$(h5cc -showconfig 2> /dev/null | sed -n 's/.*HDF5 Version: *//p'), $HDF5_FROM" "$(command -v h5cc)"
elif [ -f hdf5/lib/libhdf5_fortran.a ]; then   # same toolchain, see the Clean step
    HDF5=$PWD/hdf5
    HDF5_FROM="hdf5/"
    ok "already built"
else
    RUNLABEL="Building HDF5"
    if [ $LIVE = 1 ]; then spin_stop; render; spin_start; fi
    # build_hdf5.sh prefers the system zlib (/usr) when ZLIB_ROOT is unset
    if [ $MODE = conda ]; then export ZLIB_ROOT="$CONDA_PREFIX"; fi
    run ./hdf5/build_hdf5.sh
    HDF5=$PWD/hdf5
    HDF5_FROM="hdf5/"
    ok "$(elapsed "$T_STEP")$(warnings)"
fi

# ---------------------------------------------------------------------------
# 4. TREXIO  ->  $TREXIO, a prefix with include/trexio_f.f90 and lib/libtrexio.so
#    (where it comes from: see the header of this file)
# ---------------------------------------------------------------------------
if [ -n "$TREXIO" ]; then
    begin $S_TREXIO "Looking for TREXIO"
else
    begin $S_TREXIO "Building TREXIO $TREXIO_VERSION from source"

    if [ ! -f "$TREXIO_TARBALL" ]; then
        if download_trexio; then
            mv "$TREXIO_TARBALL.part" "$TREXIO_TARBALL"
        else
            rm -f "$TREXIO_TARBALL.part"
            if [ $MODE = conda ]; then
                HOWTO="install it into the conda env:

  conda install -c conda-forge trexio

or download"
            else
                HOWTO="download"
            fi
            die "TREXIO not found, and it could not be downloaded to build it from source.
Either $HOWTO

  $TREXIO_URL

into $PWD/trexio/ and run ./install.sh again, or pass an installed TREXIO:

  ./install.sh --trexio <prefix>"
        fi
    fi

    # Always with HDF5 (step 3). The static HDF5 of hdf5/ goes into
    # libtrexio.so, with the libraries it needs itself (zlib, szip), which
    # its libhdf5.la lists.
    TREXIO_ARGS=(--prefix="$TREXIO_LOCAL" --disable-static CC="$ICX" FC="$IFX")
    if [ -n "$HDF5" ]; then TREXIO_ARGS+=(--with-hdf5="$HDF5"); fi
    if [ "$HDF5" = "$PWD/hdf5" ]; then
        TREXIO_ARGS+=(LIBS="$(sed -n "s/^dependency_libs='\(.*\)'$/\1/p" hdf5/lib/libhdf5.la)")
    fi
    # A shared HDF5 outside the loader's default directories must be found
    # through libtrexio.so, when gammcor is linked and when it runs.
    if [ -n "$HDF5_LIBDIR" ]; then TREXIO_ARGS+=(LDFLAGS="-Wl,-rpath,$HDF5_LIBDIR"); fi

    run rm -rf trexio/build trexio/install
    run mkdir -p trexio/build
    run tar -xzf "$TREXIO_TARBALL" -C trexio/build --strip-components=1
    run env -C trexio/build ./configure "${TREXIO_ARGS[@]}"
    run make -C trexio/build -j "$JOBS"
    run make -C trexio/build install
    has_trexio "$TREXIO_LOCAL" || die "TREXIO was built, but $TREXIO_LOCAL has no
include/trexio_f.f90 or lib/libtrexio.so (details in install.log)."
    MISSING_LIBS=$(ldd "$TREXIO_LOCAL/lib/libtrexio.so" 2>&1 | grep 'not found' || true)
    [ -z "$MISSING_LIBS" ] || die "TREXIO was built, but libraries of $TREXIO_LOCAL/lib/libtrexio.so
are not found:
$MISSING_LIBS"
    TREXIO=$TREXIO_LOCAL
    TREXIO_FROM="built from source, HDF5 from $HDF5_FROM, $(elapsed "$T_STEP")$(warnings)"
fi
printf 'TREXIO: %s (%s)\n%s' "$TREXIO" "$TREXIO_FROM" "${TREXIO_SKIPPED:+$TREXIO_SKIPPED$NL}" >> "$LOG"

TV=$(trexio_version "$TREXIO")
ok "${TV:+$TV, }$TREXIO_FROM" "$TREXIO"

# ---------------------------------------------------------------------------
# 5. XCFun  ->  xcfun/lib/libxcfun.a
# ---------------------------------------------------------------------------
begin $S_XCFUN "Building XCFun"

run make -C xcfun intel GCC_INSTALL_DIR="$GCC_INSTALL_DIR"

ok "$(elapsed "$T_STEP")$(warnings)"

# ---------------------------------------------------------------------------
# 6. gammcor  ->  build/gammcor, bin/gammcor
# ---------------------------------------------------------------------------
begin $S_GAMMCOR "Building GAMMCOR"

MESON_ARGS=(-Dcpu="$CPU" -Dtrexio="$TREXIO")
if [ -n "$PREFIX" ]; then MESON_ARGS+=(-Dprefix="$PREFIX"); fi

meson_setup build . "${MESON_ARGS[@]}"
run meson compile -C build -j "$JOBS"

ok "$(elapsed "$T_STEP")$(warnings)"

# ---------------------------------------------------------------------------
# 7. Install (only with --prefix)
# ---------------------------------------------------------------------------
if [ -n "$PREFIX" ]; then
    begin $S_INSTALL "Installing into $PREFIX"
    run meson install -C build
    ok "$PREFIX/bin/gammcor"
fi

# ---------------------------------------------------------------------------
# 8. `gammcor-dsrs` command (--command; not with --prefix): a symlink to
#     bin/gammcor in a directory on PATH, see the header of this file. No step
#     of its own: it is shown in the summary at the end, and a failure here is
#     only a note.
# ---------------------------------------------------------------------------
LINK=""
LINK_NOTE=""
if [ -z "$PREFIX" ] && [ "$COMMAND" != none ]; then
    if [ $MODE = conda ]; then LINK_DIR="$CONDA_PREFIX/bin"; else LINK_DIR="$HOME/.local/bin"; fi
    OLD_TARGET=$(readlink "$LINK_DIR/$COMMAND" || true)
    if [ -e "$LINK_DIR/$COMMAND" ] && [ ! -L "$LINK_DIR/$COMMAND" ]; then
        LINK_NOTE="$COMMAND command not created: $LINK_DIR/$COMMAND exists and is not a symlink."
    # a link from another GAMMCOR directory is replaced: the last install wins
    elif ! { run mkdir -p "$LINK_DIR" && run ln -sfn "$PWD/bin/gammcor" "$LINK_DIR/$COMMAND"; }; then
        LINK_NOTE="$COMMAND command not created: linking $LINK_DIR/$COMMAND failed (see install.log)."
    else
        LINK="$LINK_DIR/$COMMAND"
        printf '=== [%s] %s command: %s\n' "$(date '+%F %T')" "$COMMAND" "$LINK" >> "$LOG"
        if [ -n "$OLD_TARGET" ] && [ "$OLD_TARGET" != "$PWD/bin/gammcor" ]; then
            LINK_NOTE="it pointed to $OLD_TARGET before"
        fi
        hash -r
        case ":$PATH:" in
            *":$LINK_DIR:"*)
                FOUND=$(command -v "$COMMAND" || true)
                if [ "$FOUND" != "$LINK" ]; then
                    LINK_NOTE+="${LINK_NOTE:+$NL}another $COMMAND comes first on PATH: $FOUND"
                fi ;;
            *)
                LINK_NOTE+="${LINK_NOTE:+$NL}$LINK_DIR is not on PATH. Add to ~/.bashrc:$NL"
                LINK_NOTE+="  export PATH=\"\$HOME/.local/bin:\$PATH\"" ;;
        esac
    fi
fi

# ---------------------------------------------------------------------------
# Done
# ---------------------------------------------------------------------------
# Only a run that got here vouches for the artifacts: record its toolchain.
{
    echo "# toolchain of the last successful ./install.sh, $(date)"
    fingerprint
} > "$STAMP"

FOOTER="${GREEN}Done${RESET} in ${CYAN}$(elapsed "$T_TOTAL")${RESET}"
if [ $LIVE = 1 ]; then
    render
else
    printf '\n%s%s%s\n%s\n' "$DIM" "$RULE" "$RESET" "$FOOTER"
fi
printf '%sNo electrons were harmed during compilation.%s\n' "$LAV" "$RESET"
printf '\n  executable  %s%s%s\n' "$CYAN" "$PWD/bin/gammcor" "$RESET"
if [ -n "$PREFIX" ]; then printf '  installed   %s%s%s\n' "$CYAN" "$PREFIX/bin/gammcor" "$RESET"; fi
if [ -n "$LINK" ]; then printf '  command     %s%s%s\n' "$CYAN" "$COMMAND" "$RESET"; fi
if [ -n "$LINK_NOTE" ]; then
    while IFS= read -r line; do printf '              %s%s%s\n' "$YELLOW" "$line" "$RESET"; done <<< "$LINK_NOTE"
fi
printf '  log         %s%s%s\n\n' "$DIM" "$LOG" "$RESET"
