#!/usr/bin/env python3
"""Pack a legacy GammCor job directory (.bin/.txt) into schema-v1 HDF5.

Run inside an old job folder:

    legacy2h5.py                    # -> pyscf_data.h5
    legacy2h5.py OLDJOB -o out.h5   # explicit dir / name
    legacy2h5.py --dry-run          # report only, write nothing

Reproduces the layout written by Pyscf2Gammcor_mini.py.  Anything the legacy
run never produced (SYSTEM/*, INTS/OVERLAP, AO_LABELS, ...) is simply absent
from the output -- the Fortran side must probe with h5ltpath_valid_f.

Requires numpy + h5py only.  No PySCF.
"""

import argparse
import os
import re
import sys
import time
from glob import glob

import numpy as np
import h5py

F8 = np.dtype(np.float64)
I4 = np.dtype(np.int32)

# --- auxdata*.txt field names, in the order Pyscf2Gammcor_mini writes them ---
AUX_BASE = ['NBASIS', 'NI', 'NA', 'NV', 'ECAS', 'ENUC', 'NEL', 'NATORB', 'FROZEN']
AUX_INT = {'NBASIS', 'NI', 'NA', 'NV', 'NEL', 'NATORB', 'FROZEN'}
PER_STATE_AUX = ('ECAS', 'EROHF')

# --- legacy file(s) -> (hdf5 path, shape spec, dtype, attrs) ----------------
# The first element is the candidate legacy name, or a tuple of them tried in
# order -- older jobs spell the same quantity differently (C.bin / CAONO.bin,
# rdm1.bin / occ.bin), exactly as interface_pp.f90 used to probe for them.
# shape spec entries are dimension keys resolved at runtime; 'na' = ncas,
# 'nb' = nbasis, 'eri' = the NAddr3-packed triangular length.
GLOBAL_MAP = [
    ('HCore.bin',     'INTS/CORE_HAMILTONIAN', ('nb', 'nb'), F8, {}),
    (('C.bin', 'CAONO.bin'), 'POSTHF/ORB_COEFF', ('nb', 'nb'), F8,
     {'DIMS': ['mo', 'ao']}),
    (('rdm1.bin', 'occ.bin'), 'POSTHF/OCC',      ('nb',),      F8, {}),
    ('occ_rohf.bin',  'SCF/MO_OCC',            ('nb',),      F8, {}),
    ('mo_occ_int.bin', 'SCF/MO_OCC_INT',       ('nb',),      I4, {}),
    ('TWOEl.bin',     'MOINTS/ERI',            ('eri',),     F8,
     {'PACKING': 'NADDR3'}),
]

RDM2_ATTRS_G = {'INDEX_ORDER': 'GAMMCOR', 'SPIN_AVERAGED': 1}

# stem (without the _<tag> suffix) -> (dataset name, shape, attrs)
STATE_MAP = [
    ('rdm2_full_reordered', 'RDM2_FULL_REORDERED', ('na',) * 4,
     {'INDEX_ORDER': 'GAMMCOR'}),
    ('rdm2_full',           'RDM2',      ('na',) * 4, {'INDEX_ORDER': 'PYSCF'}),
    ('rdm2_aaaa',           'RDM2_AAAA', ('na',) * 4, RDM2_ATTRS_G),
    ('rdm2_bbbb',           'RDM2_BBBB', ('na',) * 4, RDM2_ATTRS_G),
    ('rdm2_abab',           'RDM2_ABAB', ('na',) * 4, RDM2_ATTRS_G),
    ('rdm2_baba',           'RDM2_BABA', ('na',) * 4, RDM2_ATTRS_G),
    ('rdm1_full',           'RDM1',      ('na', 'na'), {}),
    ('rdm1a',               'RDM1_A',    ('na', 'na'), {}),
    ('rdm1b',               'RDM1_B',    ('na', 'na'), {}),
    ('rdm1p',               'OCC_A',     ('nb',),      {}),
    ('rdm1m',               'OCC_B',     ('nb',),      {}),
]

TAG_RE = re.compile(r'^(?P<stem>.+?)_(?P<tag>\d+\.\d+)\.bin$')


class Report:
    def __init__(self):
        self.written, self.skipped, self.problems = [], [], []

    def ok(self, path, src, shape):
        self.written.append((path, src, shape))

    def miss(self, src, why):
        self.skipped.append((src, why))

    def bad(self, src, why):
        self.problems.append((src, why))

    def show(self):
        for p, s, sh in self.written:
            print('  wrote  %-42s <- %-34s %s' % (p, s, sh))
        for s, w in self.skipped:
            print('  skip   %-42s    %s' % (s, w))
        for s, w in self.problems:
            print('  PROBLEM %-41s    %s' % (s, w), file=sys.stderr)


# ---------------------------------------------------------------------------
#  discovery
# ---------------------------------------------------------------------------
def find_states(jobdir):
    """Return sorted state tags, and whether the job used per-state filenames."""
    tags = set()
    for f in os.listdir(jobdir):
        m = TAG_RE.match(f)
        if m and any(m.group('stem') == s for s, _, _, _ in STATE_MAP):
            tags.add(m.group('tag'))
    if tags:
        # numeric sort: 1.1, 2.1, ... 10.1  (not lexicographic)
        return sorted(tags, key=lambda t: tuple(int(x) for x in t.split('.'))), True
    # state-specific run: untagged filenames, single state
    for stem, _, _, _ in STATE_MAP:
        if os.path.exists(os.path.join(jobdir, stem + '.bin')):
            return ['1.1'], False
    return [], False


def read_auxdata(path):
    """auxdata*.txt -> dict.  Field names inferred from line count."""
    with open(path) as f:
        raw = [ln.strip() for ln in f if ln.strip()]
    names = list(AUX_BASE[:len(raw)])
    if not names or len(names) != len(raw):
        raise ValueError('unexpected line count %d' % len(raw))
    if 'rohf' in os.path.basename(path).lower():
        names[4] = 'EROHF'
    out = {}
    for name, val in zip(names, raw):
        out[name] = int(float(val)) if name in AUX_INT else float(val)
    return out


def infer_dims(jobdir, tags, tagged, aux, override):
    """nbasis / ncas from auxdata; fall back to file sizes when it is missing."""
    dims = {}
    if aux:
        if 'NBASIS' in aux:
            dims['nb'] = aux['NBASIS']
        if 'NA' in aux:
            dims['na'] = aux['NA']

    def suffix(stem):
        return '%s_%s.bin' % (stem, tags[0]) if tagged else '%s.bin' % stem

    if 'na' not in dims:
        for stem in ('rdm2_full_reordered', 'rdm2_full', 'rdm2_aaaa'):
            p = os.path.join(jobdir, suffix(stem))
            if os.path.exists(p) and tags:
                n = round((os.path.getsize(p) // 8) ** 0.25)
                if n ** 4 * 8 == os.path.getsize(p):
                    dims['na'] = n
                    print('note: ncas=%d inferred from %s' % (n, suffix(stem)))
                    break
    if 'nb' not in dims:
        p = os.path.join(jobdir, 'HCore.bin')
        if os.path.exists(p):
            n = round((os.path.getsize(p) // 8) ** 0.5)
            if n * n * 8 == os.path.getsize(p):
                dims['nb'] = n
                print('note: nbasis=%d inferred from HCore.bin' % n)

    dims.update({k: v for k, v in override.items() if v})
    if 'nb' in dims:
        n1 = dims['nb'] * (dims['nb'] + 1) // 2
        dims['eri'] = n1 * (n1 + 1) // 2
    return dims


# ---------------------------------------------------------------------------
#  writing
# ---------------------------------------------------------------------------
def put(f, path, data, attrs, rep, src, dry):
    if not dry:
        if path in f:
            del f[path]
        d = f.create_dataset(path, data=data)
        for k, v in (attrs or {}).items():
            d.attrs[k] = v
    rep.ok(path, src, data.shape)


def load_bin(jobdir, fname, shape_keys, dtype, dims, rep):
    """Read a legacy stream file, C-order reshape.  Returns None if unusable."""
    src = os.path.join(jobdir, fname)
    if not os.path.exists(src):
        return None
    missing = [k for k in shape_keys if k not in dims]
    if missing:
        rep.bad(fname, 'cannot size it: %s unknown' % ', '.join(missing))
        return None
    shape = tuple(dims[k] for k in shape_keys)
    want = int(np.prod(shape)) * dtype.itemsize
    have = os.path.getsize(src)
    if have != want:
        rep.bad(fname, 'size %d B, expected %d B for %s' % (have, want, shape))
        return None
    # .bin was produced by ascontiguousarray(...).tofile(): C-order, no transpose
    return np.fromfile(src, dtype=dtype).reshape(shape)


def convert(jobdir, outfile, override, dry):
    rep = Report()
    tags, tagged = find_states(jobdir)

    # --- auxdata: per-state where tagged, plus a global fallback -----------
    aux_by_tag, aux_any = {}, {}
    for cand in sorted(glob(os.path.join(jobdir, 'auxdata*.txt'))):
        base = os.path.basename(cand)
        try:
            a = read_auxdata(cand)
        except Exception as err:
            rep.bad(base, 'unparsable (%s)' % err)
            continue
        m = re.match(r'auxdata_(\d+\.\d+)\.txt$', base)
        aux_by_tag[m.group(1) if m else '*'] = a
        aux_any = aux_any or a
    if not aux_any:
        rep.miss('auxdata*.txt', 'absent; dimensions inferred from file sizes')

    dims = infer_dims(jobdir, tags, tagged, aux_any, override)
    if not dims.get('nb') and not dims.get('na'):
        print('error: no recognizable GammCor data in %s' % jobdir, file=sys.stderr)
        return 1

    method = ('ROHF' if os.path.exists(os.path.join(jobdir, 'auxdata_rohf.txt'))
              else 'SA-CASSCF' if len(tags) > 1 else 'CASSCF')
    print('%s: %d state(s) %s, nbasis=%s, ncas=%s, method=%s'
          % (jobdir, len(tags), tags or '-', dims.get('nb', '?'),
             dims.get('na', '?'), method))

    f = None if dry else h5py.File(outfile, 'w')
    try:
        if not dry:
            f.attrs['FORMAT'] = 'GAMMCOR-PYSCF'
            f.attrs['VERSION'] = 1
            f.attrs['GENERATOR'] = 'legacy2h5.py'
            f.attrs['DATE'] = time.strftime('%Y-%m-%d %H:%M:%S')
            f.attrs['METHOD'] = method
            f.attrs['PROVENANCE'] = 'converted from legacy binaries in %s' \
                                    % os.path.abspath(jobdir)

        # --- global arrays -------------------------------------------------
        natorb = aux_any.get('NATORB', 0)
        for names, path, keys, dtype, attrs in GLOBAL_MAP:
            names = (names,) if isinstance(names, str) else names
            fname = next((n for n in names
                          if os.path.exists(os.path.join(jobdir, n))), None)
            if fname is None:
                rep.miss(' / '.join(names), 'not present')
                continue
            a = load_bin(jobdir, fname, keys, dtype, dims, rep)
            if a is None:
                continue
            attrs = dict(attrs)
            if path == 'POSTHF/ORB_COEFF':
                attrs['TYPE'] = 'NATURAL' if natorb else 'CANONICAL'
            if path == 'MOINTS/ERI':
                attrs['NORB'] = int(dims['nb'])
            put(f, path, a, attrs, rep, fname, dry)

        # --- per-state arrays ----------------------------------------------
        energies = []
        for tag in tags:
            grp = 'POSTHF/STATES/%s' % tag
            a_tag = aux_by_tag.get(tag) or aux_by_tag.get('*') or {}
            for stem, dsname, keys, attrs in STATE_MAP:
                fname = '%s_%s.bin' % (stem, tag) if tagged else '%s.bin' % stem
                a = load_bin(jobdir, fname, keys, F8, dims, rep)
                if a is None:
                    continue
                put(f, '%s/%s' % (grp, dsname), a, attrs, rep, fname, dry)

            # ASCII rdm2.dat only as a fallback for a missing RDM2 binary
            dat = os.path.join(jobdir, 'rdm2.dat' if not tagged
                               else 'rdm2_%s.dat' % tag)
            has_rdm2 = any(p.endswith('%s/RDM2' % grp) for p, _, _ in rep.written)
            if os.path.exists(dat) and not has_rdm2 and 'na' in dims:
                cols = np.loadtxt(dat)
                n = dims['na']
                if cols.shape[0] == n ** 4:
                    a = cols[:, 4].reshape((n,) * 4)
                    put(f, '%s/RDM2' % grp, a, {'INDEX_ORDER': 'PYSCF'},
                        rep, os.path.basename(dat), dry)
                else:
                    rep.bad(os.path.basename(dat),
                            '%d rows, expected %d' % (cols.shape[0], n ** 4))

            # per-state scalars + group attributes
            ecas = a_tag.get('ECAS', a_tag.get('EROHF'))
            if ecas is not None:
                put(f, '%s/ECAS' % grp, np.array([ecas], F8), {}, rep,
                    'auxdata', dry)
                energies.append(ecas)
            if not dry:
                g = f.require_group(grp)
                g.attrs['LABEL'] = tag
                g.attrs['STATE_NUMBER'] = int(tag.split('.')[0])
                g.attrs['IRREP_ID'] = int(tag.split('.')[1])
                if ecas is not None:
                    g.attrs['ENERGY'] = float(ecas)

        # --- REF/ scalars + reconstructed ORB_SPACE -------------------------
        for name, val in (aux_by_tag.get('*') or aux_any or {}).items():
            if name in PER_STATE_AUX:
                continue
            arr = np.array([val], I4 if name in AUX_INT else F8)
            put(f, 'REF/%s' % name, arr, {}, rep, 'auxdata', dry)
        if {'nb'} <= set(dims) and 'NI' in aux_any and 'NA' in aux_any:
            idx = np.full(dims['nb'], 2, dtype=I4)
            idx[:aux_any['NI']] = 0
            idx[aux_any['NI']:aux_any['NI'] + aux_any['NA']] = 1
            put(f, 'REF/ORB_SPACE', idx, {}, rep, 'derived', dry)

        # --- state index ----------------------------------------------------
        if tags and not dry:
            f.create_dataset('POSTHF/NSTATES', data=np.array([len(tags)], I4))
            f.create_dataset('POSTHF/STATE_LABELS',
                             data=np.array(tags,
                                           dtype=h5py.string_dtype('utf-8')))
            if len(energies) == len(tags):
                f.create_dataset('POSTHF/STATE_ENERGY',
                                 data=np.asarray(energies, F8))
    finally:
        if f is not None:
            f.close()

    print('')
    rep.show()
    print('\n%d dataset(s) %s, %d absent, %d problem(s)'
          % (len(rep.written), 'would be written' if dry else 'written',
             len(rep.skipped), len(rep.problems)))
    if rep.problems:
        print('re-run with --nbasis/--nact if the sizes above look wrong.')
    return 1 if rep.problems else 0


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('jobdir', nargs='?', default='.')
    p.add_argument('-o', '--output', default='pyscf_data.h5')
    p.add_argument('--nbasis', type=int, help='override NBASIS')
    p.add_argument('--nact', type=int, help='override NA (ncas)')
    p.add_argument('--dry-run', action='store_true')
    a = p.parse_args()

    if not os.path.isdir(a.jobdir):
        sys.exit('not a directory: %s' % a.jobdir)
    out = a.output if os.path.isabs(a.output) \
        else os.path.join(a.jobdir, a.output)
    tmp = out + '.tmp'

    rc = convert(a.jobdir, tmp, {'nb': a.nbasis, 'na': a.nact}, a.dry_run)
    if not a.dry_run:
        if os.path.exists(tmp):
            os.replace(tmp, out)       # atomic: never leave a half-written .h5
            print('-> %s' % out)
    sys.exit(rc)


if __name__ == '__main__':
    main()
