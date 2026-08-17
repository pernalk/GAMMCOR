#!/usr/bin/env python3
"""Post-build step for gammcor.

Runs as the last edge of every `meson compile`:

  1. copies the freshly linked binary to  bin/<exe_name>
  2. works out how long the build took, from Ninja's own log
  3. prints a success banner

Invoked by the `post_build_tasks` custom_target in meson.build -- not by hand.
"""

import argparse
import shutil
import sys
from pathlib import Path

WIDTH = 78


def build_seconds(build_dir: Path):
    """Wall-clock duration of the Ninja invocation that is running us.

    .ninja_log records one tab-separated line per completed edge:

        start_ms  end_ms  output_mtime_ns  output_path  command_hash

    start_ms/end_ms are milliseconds since *that particular* Ninja invocation
    started, and lines are appended in completion order -- so within a single
    invocation end_ms only ever increases. The log is not truncated between
    builds, so walking backwards from the end of the file and stopping at the
    first place where end_ms jumps back up isolates the current run.

    Returns seconds as a float, or None if the log is missing or unparsable
    (a stale-looking number would be worse than no number).
    """
    log = build_dir / '.ninja_log'
    try:
        lines = log.read_text(errors='ignore').splitlines()
    except OSError:
        return None

    spans = []
    for line in reversed(lines):
        parts = line.split('\t')
        if len(parts) < 5:
            continue                      # header line, or a partial write
        try:
            start, end = int(parts[0]), int(parts[1])
        except ValueError:
            continue
        if spans and end > spans[-1][1]:
            break                         # crossed into the previous build
        spans.append((start, end))

    if not spans:
        return None
    return (max(e for _, e in spans) - min(s for s, _ in spans)) / 1000.0


def human(seconds):
    if seconds is None:
        return 'n/a'
    minutes, secs = divmod(seconds, 60)
    if minutes:
        return f'{int(minutes)} min {secs:04.1f} s'
    return f'{secs:.1f} s'


def banner(rows):
    rule = '=' * WIDTH
    print()
    print(rule)
    print('  gammcor -- BUILD SUCCESSFUL')
    print(rule)
    for key, value in rows:
        print(f'  {key:<13}: {value}')
    print(rule)
    print()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--exe',       required=True, type=Path)
    p.add_argument('--bin-dir',   required=True, type=Path)
    p.add_argument('--build-dir', required=True, type=Path)
    p.add_argument('--profile',   required=True)
    p.add_argument('--exe-name',  required=True)
    p.add_argument('--compiler',  default='')
    p.add_argument('--stamp',     required=True, type=Path)
    args = p.parse_args()

    if not args.exe.is_file():
        print(f'post_build: expected executable at {args.exe}', file=sys.stderr)
        return 1

    args.bin_dir.mkdir(parents=True, exist_ok=True)
    dest = args.bin_dir / args.exe_name
    shutil.copy2(args.exe, dest)

    banner([
        ('profile',    args.profile),
        ('compiler',   args.compiler.strip()),
        ('executable', dest),
        ('size',       f'{dest.stat().st_size / 1024**2:.1f} MiB'),
        ('build time', human(build_seconds(args.build_dir))),
    ])

    args.stamp.write_text('done\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
