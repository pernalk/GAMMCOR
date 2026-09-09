#!/usr/bin/env python3
"""Dump the orbital basis of a PySCF script in EMSL/GAMESS-US format.

Reads `mol.basis` out of the script and writes it to baza.dat.  The script is
parsed, never executed, so pointing this at a full calculation is harmless.
ECPs are not written.
"""

import argparse
import ast
import sys
from pathlib import Path

from pyscf import gto
from pyscf.data import elements
from pyscf.lib.exceptions import BasisNotFoundError

ANGULAR_LABELS = 'SPDFGHIKLMNOQRTUVWXYZ'


def write_gamess_us_basis(basis, filename='baza.dat'):
    """Write only the orbital basis in EMSL/GAMESS-US format."""
    formatted_basis = gto.format_basis(basis, sort_basis=False)

    with open(filename, 'w') as output:
        output.write('$DATA\n\n')

        for symbol, shells in formatted_basis.items():
            atomic_number = elements.charge(symbol)
            output.write(f'{elements.ATOMIC_NAMES[atomic_number].upper()}\n')

            for shell in shells:
                angular_momentum = shell[0]
                primitives = [row for row in shell[1:] if isinstance(row, (list, tuple))]
                if not primitives:
                    continue
                ncontractions = len(primitives[0]) - 1

                # GAMESS-US stores every generalized contraction as a
                # separate shell.  Zero coefficients can be left out.
                for contraction in range(ncontractions):
                    contracted_primitives = [
                        (primitive[0], primitive[contraction + 1])
                        for primitive in primitives
                        if primitive[contraction + 1] != 0.0
                    ]
                    if not contracted_primitives:
                        continue

                    label = ANGULAR_LABELS[angular_momentum]
                    output.write(f'{label}   {len(contracted_primitives)}\n')
                    for index, (exponent, coefficient) in enumerate(
                            contracted_primitives, start=1):
                        output.write(
                            f'{index:<4d}{exponent:22.12E}'
                            f'{coefficient:22.12E}\n'
                        )

            output.write('\n')

        output.write('$END\n')


def assigned_literal(tree, attribute):
    """Value last assigned to `<something>.attribute` or a bare `attribute`."""
    found = None
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign):
            continue
        for target in node.targets:
            name = getattr(target, 'attr', None) or getattr(target, 'id', None)
            if name == attribute:
                found = node.value
    if found is None:
        return None
    try:
        return ast.literal_eval(found)
    except ValueError:
        raise ValueError(f'`{attribute}` is not a plain literal, cannot read it '
                         f'without running the script')


def elements_of(atom, script_dir):
    """Element symbols of a `mol.atom` value, in order of first appearance.

    Handles both a path to an .xyz file and inline geometry text.
    """
    xyz = Path(atom)
    if not xyz.is_absolute():
        xyz = script_dir / xyz

    if xyz.is_file():
        lines = xyz.read_text().splitlines()
        try:
            natoms = int(lines[0].split()[0])
        except (IndexError, ValueError):
            raise ValueError(f'{xyz}: first line is not an atom count')
        lines = lines[2:2 + natoms]
    else:
        lines = str(atom).replace(';', '\n').splitlines()

    symbols = []
    for line in lines:
        fields = line.replace(',', ' ').split()
        if fields and fields[0] not in symbols:
            symbols.append(fields[0])
    return symbols


def basis_of_script(path):
    """The `mol.basis` of a pyscf script, as a dict ready for format_basis."""
    path = Path(path)
    tree = ast.parse(path.read_text(), filename=str(path))

    basis = assigned_literal(tree, 'basis')
    if basis is None:
        raise ValueError(f'{path}: no assignment to `mol.basis` found')
    if isinstance(basis, dict):
        return basis

    # A bare basis name applies to every element, so the geometry is needed.
    atom = assigned_literal(tree, 'atom')
    if atom is None:
        raise ValueError(f'{path}: `mol.basis` is the name {basis!r}, but there '
                         f'is no `mol.atom` to say which elements it applies to')
    symbols = elements_of(atom, path.parent)
    if not symbols:
        raise ValueError(f'{path}: could not read any element from `mol.atom`')
    return {symbol: basis for symbol in symbols}


def main():
    parser = argparse.ArgumentParser(
        prog='basis_dump',
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='example:\n  basis_dump CO_NiAg_dmrgscf.py\n')
    parser.add_argument('script', help='pyscf script to read `mol.basis` from')
    parser.add_argument('-o', '--output', default='baza.dat',
                        help='output file (default: baza.dat)')
    args = parser.parse_args()

    try:
        basis = basis_of_script(args.script)
    except (OSError, SyntaxError, ValueError) as err:
        print(f'error: {err}', file=sys.stderr)
        return 1

    try:
        write_gamess_us_basis(basis, args.output)
    except BasisNotFoundError as err:
        print(f'error: {err}', file=sys.stderr)
        return 1

    print(f'Orbital basis for {", ".join(basis)} written to '
          f'{args.output} in GAMESS-US format (ECP omitted).')
    return 0


if __name__ == '__main__':
    sys.exit(main())
