project(
  'gammcor',
  'fortran',
  version : '4.1',
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

fc   = meson.get_compiler('fortran')
root = meson.project_source_root()

# Profile metadata, set in the [properties] section of the native file.
profile  = meson.get_external_property('profile',  'release')
exe_name = meson.get_external_property('exe_name', 'gammcor')

# ---------------------------------------------------------------------------
# Pre-built libraries shipped inside the tree.
#
# These are project structure, identical in every profile, so they live here
# and not in the .ini files. Compiler/optimisation/MKL flags are the opposite:
# they are profile-specific and live in .meson/profiles/*.ini.
# ---------------------------------------------------------------------------

# XCFun: static lib plus two pre-built .o holding the Fortran wrappers.
xcfun_dep = declare_dependency(
  include_directories : include_directories('xcfun/fortran'),
  dependencies : fc.find_library('xcfun', dirs : root / 'xcfun/lib', required : true),
)
xcfun_objects = files(
  'xcfun/fortran/xcfun_module.o',
  'xcfun/fortran/xcfun_autogen.o',
)

# gammcor-integrals: cholesky.a has a non-standard name, so link it by path.
cholesky_dep = declare_dependency(
  include_directories : include_directories('gammcor-integrals/include'),
  link_args : [root / 'gammcor-integrals/lib/cholesky.a'],
)

# HDF5: built into hdf5/ (see hdf5/README.md).
hdf5_dep = declare_dependency(
  include_directories : include_directories('hdf5/include'),
  link_args : ['-L' + root / 'hdf5/lib', '-Wl,-rpath,' + root / 'hdf5/lib',
               '-lhdf5_hl_fortran', '-lhdf5_fortran', '-lhdf5_hl', '-lhdf5',
               '-lz', '-lsz', '-ldl', '-lm'],
)

# ---------------------------------------------------------------------------
# Sources
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
  'SOURCE/units.f90',
  'SOURCE/print_utils.f90',
  'SOURCE/types.f90',
  'SOURCE/blocktypes.f90',
  'SOURCE/memory.f90',
  'SOURCE/timing.f90',
  'SOURCE/diis.f90',
  'SOURCE/lin.f90',
  'SOURCE/inputfill.f90',
  'SOURCE/systemdef.f90',
  'SOURCE/rdwrext.f90',
  'SOURCE/sorter.f90',
  'SOURCE/git_info.f90',
  'SOURCE/build_info.f90',
  'SOURCE/h5_reader.f90',

  # -- integrals / transformation / Cholesky
  'SOURCE/gammint.f90',
  'SOURCE/geom_input.f90',
  'SOURCE/tran.f90',
  'SOURCE/tran_Chol.f90',
  'SOURCE/CholeskyOTF.F90',
  'SOURCE/interface_pp.f90',

  # -- grids and DFT
  'SOURCE/gridmolpro.f90',
  'SOURCE/griddalton.f90',
  'SOURCE/gridintern.f90',
  'SOURCE/dftgrid.f',
  'SOURCE/lsd_sr.f',
  'SOURCE/dftfun_exerfpbe.f',
  'SOURCE/dftfun_exerf.f',
  'SOURCE/dftfun_ecerfpbe.f',
  'SOURCE/dftfun_ecerf.f',
  'SOURCE/dftacg_pw92c.f',
  'SOURCE/caspidft.f',
  'SOURCE/vv10.f',
  'SOURCE/muerf.f',
  'SOURCE/xcfun.f90',

  # -- A/B matrices
  'SOURCE/abmats.f90',
  'SOURCE/abfofo.f90',
  'SOURCE/ab0fofo.f90',
  'SOURCE/abchol.f90',
  'SOURCE/abuks.f90',

  # -- adiabatic connection
  'SOURCE/acfd.f',
  'SOURCE/accas.f',
  'SOURCE/ac_iter.f',
  'SOURCE/ac_fofo.f90',
  'SOURCE/ac_exact_2el.f',
  'SOURCE/ppac_types.f90',

  # -- pp path: shared helpers
  'SOURCE/quadratures.f90',
  'SOURCE/pp_utils.f90',

  # -- pp path: THC back end
  'SOURCE/thc_auto_p0.f90',
  'SOURCE/thc_auto_p2.f90',
  'SOURCE/thc_auto_p3.f90',
  'SOURCE/thc_energy.f90',

  # -- pp path: AC0 kernels, one per integral back end
  'SOURCE/ppac0.f90',
  'SOURCE/ppac0_thc.f90',

  # -- pp path: drivers
  'SOURCE/ppac_subs.f90',
  'SOURCE/ppac_dispatch.f90',
  'SOURCE/ppac_simple_subs.f90',
  'SOURCE/ppac_simple.f90',

  # -- ph path and job routing
  'SOURCE/phac_spinres.f90',
  'SOURCE/mp.f90',
  'SOURCE/dispatcher.f90',
  'SOURCE/polari.F90',

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
  'SOURCE/sapt_interface.f90',
  'SOURCE/sapt_Chol_pol.f90',
  'SOURCE/sapt_Chol_exch.f90',
  'SOURCE/sapt_visual.f90',
  'SOURCE/sapt_files.f90',
  'SOURCE/sapt_os.f90',

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
  sources,
  objects : xcfun_objects,
  # SOURCE/ is on the include path so `include 'commons.inc'` resolves
  # regardless of the compiler's include-resolution rules.
  include_directories : include_directories('SOURCE'),
  dependencies : [xcfun_dep, cholesky_dep, hdf5_dep],
)

# ---------------------------------------------------------------------------
# Post-build: copy the binary to bin/<exe_name>, report build time, print the
# success banner. `console : true` makes Ninja stream its output straight to
# the terminal, and runs it last.
# ---------------------------------------------------------------------------
python     = import('python').find_installation('python3')
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
  'executable' : 'bin' / exe_name,
}, section : 'gammcor')
