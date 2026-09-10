#!/usr/bin/env python3

import argparse
from pathlib import Path

import h5py
import numpy as np


def state_key(name):
    try:
        return tuple(int(part) for part in name.split('.'))
    except ValueError:
        return (name,)


def scalar(h5, state, name):
    state_path = f'{state}/{name}'
    path = state_path if state_path in h5 else f'REF/{name}'
    return np.asarray(h5[path][()]).reshape(-1)[0].item()


def write_array(h5, path, output, filename, dtype=np.float64):
    if path not in h5:
        return None
    np.ascontiguousarray(h5[path][()], dtype=dtype).tofile(output / filename)
    return filename


def write_aux(h5, state, names, output, filename, trailing_newline):
    if not all(f'{state}/{name}' in h5 or f'REF/{name}' in h5 for name in names):
        return None
    content = '\n'.join(str(scalar(h5, state, name)) for name in names)
    if trailing_newline:
        content += '\n'
    (output / filename).write_text(content)
    return filename


def write_rdm2_dat(dataset, output, filename, rohf=False):
    if dataset.ndim != 4 or len(set(dataset.shape)) != 1:
        raise ValueError(f'{dataset.name} must be a rank-4 square array')
    n = dataset.shape[0]
    with (output / filename).open('w') as handle:
        for i in range(n):
            block = np.asarray(dataset[i])
            lines = []
            for j in range(n):
                for k in range(n):
                    for l in range(n):
                        value = block[j, k, l]
                        if rohf:
                            lines.append(
                                f'{i + 1:5d}{j + 1:5d}{k + 1:5d}{l + 1:5d}'
                                f'    {value:20.12f}\n'
                            )
                        else:
                            lines.append(
                                f'{i + 1:4d} {j + 1:4d} {k + 1:4d} {l + 1:4d}'
                                f' {value:19.12f}\n'
                            )
            handle.writelines(lines)
    return filename


def get_method(h5):
    method = h5.attrs.get('METHOD')
    if isinstance(method, bytes):
        method = method.decode()
    if method:
        return str(method).upper()
    if 'SCF/MO_OCC_INT' in h5:
        return 'ROHF'
    if 'POSTHF/STATE_WEIGHTS' in h5:
        return 'SA-CASSCF'
    return 'CASSCF'


def dump_casscf(h5, output, state_names, state_averaged):
    written = []
    common = (
        ('POSTHF/ORB_COEFF', 'C.bin'),
        ('INTS/CORE_HAMILTONIAN', 'HCore.bin'),
        ('MOINTS/ERI', 'TWOEl.bin'),
    )
    for path, filename in common:
        result = write_array(h5, path, output, filename)
        if result:
            written.append(result)

    if state_averaged:
        for tag in state_names:
            state = f'POSTHF/STATES/{tag}'
            arrays = (
                ('RDM2_AAAA', f'rdm2_aaaa_{tag}.bin'),
                ('RDM2_ABAB', f'rdm2_abab_{tag}.bin'),
                ('RDM2', f'rdm2_full_{tag}.bin'),
                ('RDM2_FULL_REORDERED', f'rdm2_full_reordered_{tag}.bin'),
                ('RDM1_A', f'rdm1a_{tag}.bin'),
                ('RDM1_B', f'rdm1b_{tag}.bin'),
                ('RDM1', f'rdm1_full_{tag}.bin'),
            )
            for dataset, filename in arrays:
                result = write_array(h5, f'{state}/{dataset}', output, filename)
                if result:
                    written.append(result)
            result = write_aux(
                h5, state,
                ('NBASIS', 'NI', 'NA', 'NV', 'ECAS', 'ENUC', 'NEL', 'NATORB'),
                output, f'auxdata_{tag}.txt', False,
            )
            if result:
                written.append(result)
        return written

    if not state_names:
        return written
    state = f'POSTHF/STATES/{state_names[0]}'
    arrays = (
        ('RDM2_AAAA', 'rdm2_aaaa.bin'),
        ('RDM2_BBBB', 'rdm2_bbbb.bin'),
        ('RDM2_ABAB', 'rdm2_abab.bin'),
        ('RDM2_BABA', 'rdm2_baba.bin'),
        ('RDM2', 'rdm2_full.bin'),
        ('RDM2_FULL_REORDERED', 'rdm2_full_reordered.bin'),
        ('RDM1_A', 'rdm1a.bin'),
        ('RDM1_B', 'rdm1b.bin'),
        ('OCC_A', 'rdm1p.bin'),
        ('OCC_B', 'rdm1m.bin'),
    )
    for dataset, filename in arrays:
        result = write_array(h5, f'{state}/{dataset}', output, filename)
        if result:
            written.append(result)
    result = write_array(h5, 'POSTHF/OCC', output, 'rdm1.bin')
    if result:
        written.append(result)
    result = write_aux(
        h5, state,
        ('NBASIS', 'NI', 'NA', 'NV', 'ECAS', 'ENUC', 'NEL', 'NATORB', 'FROZEN'),
        output, 'auxdata.txt', True,
    )
    if result:
        written.append(result)
    if f'{state}/RDM2_AAAA' in h5 and f'{state}/RDM2' in h5:
        written.append(write_rdm2_dat(h5[f'{state}/RDM2'], output, 'rdm2.dat'))
    return written


def dump_rohf(h5, output, state_names):
    written = []
    arrays = (
        ('SCF/MO_OCC_INT', 'mo_occ_int.bin', np.int32),
        ('SCF/MO_OCC', 'occ_rohf.bin', np.float64),
        ('INTS/CORE_HAMILTONIAN', 'HCore.bin', np.float64),
        ('POSTHF/ORB_COEFF', 'C.bin', np.float64),
    )
    for path, filename, dtype in arrays:
        result = write_array(h5, path, output, filename, dtype)
        if result:
            written.append(result)
    if not state_names:
        return written
    state = f'POSTHF/STATES/{state_names[0]}'
    result = write_aux(
        h5, state,
        ('NBASIS', 'NI', 'NA', 'NV', 'EROHF', 'ENUC', 'NEL'),
        output, 'auxdata_rohf.txt', True,
    )
    if result:
        written.append(result)
    if f'{state}/RDM2' in h5:
        written.append(
            write_rdm2_dat(h5[f'{state}/RDM2'], output, 'rdm2_rohf.dat', True)
        )
    return written


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('hdf5')
    parser.add_argument('-o', '--output-dir', default='.')
    args = parser.parse_args()
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)

    with h5py.File(args.hdf5, 'r') as h5:
        states = h5.get('POSTHF/STATES')
        state_names = sorted(states.keys(), key=state_key) if states else []
        method = get_method(h5)
        if method == 'ROHF':
            written = dump_rohf(h5, output, state_names)
        elif method in ('CASSCF', 'SA-CASSCF'):
            written = dump_casscf(h5, output, state_names, method == 'SA-CASSCF')
        else:
            raise ValueError(f'unsupported METHOD: {method}')

    for filename in written:
        print(filename)
    print(f'Wrote {len(written)} files to {output.resolve()}')


if __name__ == '__main__':
    main()
