# Building gammcor with Meson

The Meson build is independent of the Makefile. Both build systems describe the
same sources and can be used interchangeably; neither is required by the other,
and they write to separate locations, so builds do not interfere.

## Prerequisites

Meson ≥ 0.56 and Ninja.

This must be done inside the Python virtual environment. First activate it:

```bash
source ~/.virtualenvs/my_env/bin/activate
```

Then install Meson:

```bash
pip install meson
```

or

```bash
pip3 install meson
```

## Layout

The build is described by two kinds of file:

- `meson.build` -- the build script: the source list, the pre-built libraries to
  link, the executable target and the tests.
- `.meson/profiles/*.ini` -- the profiles, one per build variant. Compiler,
  optimisation, OpenMP, coarray and MKL settings are defined per profile.

Nothing in `meson.build` is compiler- or optimisation-specific, so a new build
variant means a new `.ini` file and no change to the build script.

## Configuration

Do only once, unless you delete the `build` or `build-debug` folders for some
reason. All flags and settings are set inside `.meson/profiles/ifx.ini` and
`.meson/profiles/ifx-debug.ini`.

```bash
meson setup build       --native-file .meson/profiles/ifx.ini
meson setup build-debug --native-file .meson/profiles/ifx-debug.ini
```

The profile is bound to the build directory at setup time. Re-run only after
deleting a build directory or editing a profile.

## Compilation

```bash
meson compile -C build
meson compile -C build-debug
```

The build directories are independent and can be rebuilt in any order.

## Artifacts

The executables are installed to:

```
path-to-gammcor/bin/gammcor          # release
path-to-gammcor/bin/gammcor-debug    # debug
```

Each is also left in its build directory, as `build/gammcor` and
`build-debug/gammcor` respectively.

## Tests

Not wired up yet. Once a test suite lands in the tree, `test()` targets go back
into `meson.build` and `meson test -C build` picks them up.

## Cleaning

`meson compile` does not require a `make clean` first. Ninja records what
produced every object -- the exact command line and every file it depends on,
`.mod` files included -- and recompiles whatever no longer matches. Changed
sources, changed compiler flags and changed module dependencies are therefore
all handled on their own. `make clean` exists because a Makefile misses those
cases and leaves stale objects behind; here there is nothing stale to clear.

If for some reason you want to do it anyway, simply delete the `build` or
`build-debug` folder:

```bash
rm -rf build
```

This discards the configuration along with the objects, so the build directory
must be set up again with its `--native-file` before it can be compiled:

```bash
meson setup build --native-file .meson/profiles/ifx.ini
```

## Adding a new module

Add the file to the `sources` list in `meson.build`, under the comment block it
belongs to: 

```meson
sources = files(
  ...
  'SOURCE/my_new_module.f90',
)
```

Then compile as usual:

```bash
meson compile -C build
```

There is no glob and no dependency list to maintain: the explicit list is the
only thing you edit. Order does not matter -- Meson reads the `module` and `use`
statements itself and works out the compile order. Adding a file changes
`meson.build`, which Ninja notices, so `meson setup` does not need to be re-run.
