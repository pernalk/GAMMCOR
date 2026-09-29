project(
  'gammcor',
  'fortran',
  version : '4.1',
  meson_version : '>=1.1',
  default_options : [
    'warning_level=0',   # warnings are set per profile (-warn nounused)
    'buildtype=plain',   # Meson adds NO flags of its own: everything the
                         # compiler sees comes from .meson/profiles/*.ini
    'backend=ninja',
  ],
)

if meson.backend() != 'ninja'
  error('This project only supports the Ninja backend.')
endif

fs   = import('fs')
fc   = meson.get_compiler('fortran')
root = meson.project_source_root()

# Profile metadata, set in the [properties] section of the native file.
profile  = meson.get_external_property('profile',  'release')
exe_name = meson.get_external_property('exe_name', 'gammcor')

python = import('python').find_installation('python3')

cpu = get_option('cpu')
cpu_info = cpu
if cpu == 'auto'
  detect_cpu = root / '.meson/scripts/detect_cpu.py'
  cpu_lines = run_command(python, detect_cpu, check : true).stdout().strip().split('\n')
  cpu = cpu_lines[0]
  cpu_info = cpu + ' (auto: ' + cpu_lines[1] + ')'
endif
cpu_flags = {
  'avx512'     : ['-xCORE-AVX512', '-qopt-zmm-usage=high'],
  'avx512-amd' : ['-march=x86-64-v4'],
  'avx2'       : ['-march=x86-64-v3'],
  'native'     : ['-xHost'],
  'dispatch'   : ['-march=x86-64-v3', '-axCORE-AVX512', '-qopt-zmm-usage=high'],
}[cpu]
add_project_arguments(cpu_flags, language : 'fortran')
add_project_link_arguments(cpu_flags, language : 'fortran')

# ---------------------------------------------------------------------------
# Libraries
#
# These are project structure, identical in every profile, so they live here
# and not in the .ini files. Compiler/optimisation/MKL flags are the opposite:
# they are profile-specific and live in .meson/profiles/*.ini.
# ---------------------------------------------------------------------------

# XCFun: static lib plus two pre-built .o holding the Fortran wrappers,
# built by `make -C xcfun intel` (install.sh does it).
if not fs.exists(root / 'xcfun/lib/libxcfun.a')
  error('XCFun is not built. Run:  make -C xcfun intel   (or ./install.sh)')
endif
xcfun_dep = declare_dependency(
  include_directories : include_directories('xcfun/fortran'),
  dependencies : fc.find_library('xcfun', dirs : root / 'xcfun/lib', required : true),
)
xcfun_objects = files(
  'xcfun/fortran/xcfun_module.o',
  'xcfun/fortran/xcfun_autogen.o',
)

# TREXIO: the C library, plus its Fortran interface trexio_f.f90 (module
# `trexio`), which is compiled here as one of the sources.
#
# The Makefile/CMake builds compile SOURCE/trexio_module.F90 instead, which
# only does `#include "trexio_f.f90"`. Meson finds the compile order by
# scanning each source for `module`/`use` and does not follow #include, so
# it would not know that `use trexio` needs trexio_module.F90 first. Taking
# trexio_f.f90 itself as a source gives the same trexio.mod and a correct
# order.
trexio_prefix = get_option('trexio')
if trexio_prefix == ''
  trexio_pc = dependency('trexio', method : 'pkg-config', required : false)
  if not trexio_pc.found()
    error('TREXIO not found by pkg-config. Pass -Dtrexio=<prefix>, or run ./install.sh')
  endif
  trexio_incdir = trexio_pc.get_variable(pkgconfig : 'includedir')
  trexio_libdir = trexio_pc.get_variable(pkgconfig : 'libdir')
else
  trexio_incdir = trexio_prefix / 'include'
  trexio_libdir = trexio_prefix / 'lib'
  if not fs.exists(trexio_libdir / 'libtrexio.so') and fs.exists(trexio_prefix / 'lib64/libtrexio.so')
    trexio_libdir = trexio_prefix / 'lib64'
  endif
endif
if not fs.exists(trexio_incdir / 'trexio_f.f90')
  error('TREXIO: ' + trexio_incdir / 'trexio_f.f90' + ' not found')
endif
# TREXIO >= 2.5.0: the transition RDMs read by sapt_interface.F90
if not fs.read(trexio_incdir / 'trexio_f.f90').to_lower().contains('trexio_read_rdm_2e_transition_size')
  error('TREXIO in ' + trexio_incdir + ' is too old: TREXIO >= 2.5.0 is needed (transition RDMs)')
endif
trexio_dep = declare_dependency(
  dependencies : fc.find_library('trexio', dirs : trexio_libdir, required : true),
  # the directory is not on the loader path (conda env, trexio/install, ...)
  link_args : ['-Wl,-rpath,' + trexio_libdir],
)
trexio_f90 = configure_file(
  input  : trexio_incdir / 'trexio_f.f90',
  output : 'trexio_f.f90',
  copy   : true,
)

# ---------------------------------------------------------------------------
# Sources: the list of the Makefile of this branch (SOURCE/trexio_module.F90
# replaced by trexio_f.f90, see above).
#   Meson deliberately has no glob: an explicit list is what makes an
#   incremental build reproducible. Order does NOT matter -- Meson scans
#   `module`/`use` statements itself and derives the compile order.
# ---------------------------------------------------------------------------
sources = files(
  # -- core / driver
  'SOURCE/mainp.f',
  'SOURCE/initia.f',
  'SOURCE/dmscf.f',
  'SOURCE/misc.f',
  'SOURCE/optocc.f',
  'SOURCE/optnorb.f',
  'SOURCE/eigenno.f',
  'SOURCE/matvec.f',
  'SOURCE/matrix.f',
  'SOURCE/hf.f',
  'SOURCE/cpdmft.f',
  'SOURCE/nonadia.f',
  'SOURCE/projector.f',
  'SOURCE/ekt.f',
  'SOURCE/newton.f',
  'SOURCE/optapsg.f',
  'SOURCE/exact2el.f',
  'SOURCE/erpa.f',
  'SOURCE/interpa.f',
  'SOURCE/srlrdynamic.f',

  # -- infrastructure
  'SOURCE/priunit.f90',
  'SOURCE/types.f90',
  'SOURCE/blocktypes.f90',
  'SOURCE/timing.f90',
  'SOURCE/diis.f90',
  'SOURCE/inputfill.f90',
  'SOURCE/systemdef.f90',
  'SOURCE/rdwrext.f90',
  'SOURCE/sorter.f90',
  'SOURCE/git_info.f90',
  'SOURCE/build_info.f90',
  'SOURCE/arithmetic.f90',
  'SOURCE/io.f90',
  'SOURCE/clock.f90',
  'SOURCE/string.f90',
  'SOURCE/sort.f90',
  'SOURCE/display.f90',

  # -- TREXIO
  'SOURCE/write_trexio.f90',

  # -- integrals / transformation / Cholesky
  'SOURCE/tran.f90',
  'SOURCE/tran_Chol.f90',
  'SOURCE/sorter_Cholesky.f90',
  'SOURCE/Cholesky.f90',
  'SOURCE/Cholesky_old.f90',

  # -- grids and DFT
  'SOURCE/gridmolpro.f90',
  'SOURCE/dftgrid.f',
  'SOURCE/lsd_sr.f',
  'SOURCE/dftfun_exerfpbe.f',
  'SOURCE/dftfun_exerf.f',
  'SOURCE/dftfun_ecerfpbe.f',
  'SOURCE/dftfun_ecerf.f',
  'SOURCE/dftacg_pw92c.f',
  'SOURCE/caspidft.f',
  'SOURCE/vv10.f',
  'SOURCE/xcfun.f90',

  # -- A/B matrices
  'SOURCE/abmats.f90',
  'SOURCE/abfofo.f90',
  'SOURCE/ab0fofo.f90',
  'SOURCE/abchol.f90',
  'SOURCE/abcipsi.f90',

  # -- adiabatic connection
  'SOURCE/acfd.f',
  'SOURCE/accas.f',
  'SOURCE/ac_iter.f',
  'SOURCE/ac_fofo.f90',
  'SOURCE/ac_exact_2el.f',

  # -- iterative solvers
  'SOURCE/iter_stats.f90',
  'SOURCE/C_integrator.f90',
  'SOURCE/iter_algorithm.f90',
  'SOURCE/iter_algorithm_diis.f90',
  'SOURCE/iter_algorithm_damping.f90',
  'SOURCE/lambda_calculator.f90',
  'SOURCE/lambda_calculator_diag.f90',
  'SOURCE/lambda_calculator_block.f90',
  'SOURCE/lambda_calculator_projector.f90',

  # -- SAPT
  'SOURCE/sapt_main.f90',
  'SOURCE/sapt_exch.f90',
  'SOURCE/sapt_pol.f90',
  'SOURCE/sapt_utils.f90',
  'SOURCE/sapt_response.f90',
  'SOURCE/sapt_interface.F90',
  'SOURCE/sapt_dRS.F90',
  'SOURCE/sapt_os.f90',
  'SOURCE/sapt_Chol_pol.f90',
  'SOURCE/sapt_Chol_exch.f90',

  # -- exchange
  'SOURCE/exmisc.f90',
  'SOURCE/exdpino.f90',
  'SOURCE/exi.f90',
  'SOURCE/exappr.f90',
  'SOURCE/srefex.f90',
)

# ---------------------------------------------------------------------------
# Target
# ---------------------------------------------------------------------------
gammcor_exe = executable(
  'gammcor',
  sources, trexio_f90,
  objects : xcfun_objects,
  # SOURCE/ is on the include path so `include 'commons.inc'` resolves
  # regardless of the compiler's include-resolution rules.
  include_directories : include_directories('SOURCE'),
  dependencies : [xcfun_dep, trexio_dep],
  # `meson install -C build` puts it into <prefix>/bin/gammcor.
  install : true,
)

# ---------------------------------------------------------------------------
# Post-build: copy the binary to bin/<exe_name>, report build time, print the
# success banner. `console : true` makes Ninja stream its output straight to
# the terminal, and runs it last.
# ---------------------------------------------------------------------------
post_build = files('.meson/scripts/post_build.py')

post_build_target = custom_target('post_build_tasks',
  input   : gammcor_exe,
  output  : 'post_build.done',
  command : [
    python, post_build,
    '--exe',       '@INPUT@',
    '--bin-dir',   root / 'bin',
    '--build-dir', meson.project_build_root(),
    '--profile',   profile,
    '--exe-name',  exe_name,
    '--compiler',  fc.cmd_array()[0] + ' ' + fc.version(),
    '--stamp',     '@OUTPUT@',
  ],
  console          : true,
  build_by_default : true,
)

summary({
  'profile'    : profile,
  'compiler'   : fc.get_id() + ' ' + fc.version(),
  'cpu'        : cpu_info,
  'trexio'     : trexio_libdir,
  'executable' : 'bin' / exe_name,
}, section : 'gammcor')
