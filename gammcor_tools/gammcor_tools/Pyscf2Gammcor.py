import h5py
import numpy as np
import time
import os
from pyscf import gto, scf, mcscf, ao2mo, symm, tools, fci
from pyscf.tools import mo_mapping
from collections import Counter
from pyscf.tools import molden
import sys
from copy import deepcopy
from functools import reduce


# ---------------------------------------------------------------------------
#  Organized HDF5 dump (schema v1)
#
#  pyscf_data.h5
#    SYSTEM/   geometry, nuclear repulsion, molecule metadata
#    INTS/     AO one-electron integrals   (CORE_HAMILTONIAN <- HCore.bin)
#    MOINTS/   ERI, NAddr3-packed          (<- TWOEl.bin)
#    SCF/      HF reference orbitals, energies, occupations
#    REF/      orbital-space dimensions and flags   (<- auxdata*.txt)
#    POSTHF/   ORB_COEFF (<- C.bin), OCC (<- rdm1.bin),
#              STATES/<state>.<sym>/  RDM1*, RDM2*, OCC_A, OCC_B, ...
#
#  Arrays are written exactly as they were handed to .tofile(), so every
#  dataset is byte-identical to the .bin file it replaces and the Fortran
#  side reads the same numbers it reads today.
#
#  get_data_for_gammcor(..., dump_hdf5=True)  -> one HDF5 file
#  get_data_for_gammcor(..., dump_hdf5=False) -> the legacy loose files
# ---------------------------------------------------------------------------

H5FILE = "pyscf_data.h5"


class _Dump:
    """Single sink for everything the interface exports."""

    def __init__(self, dump_hdf5, filename=H5FILE):
        self.hdf5 = dump_hdf5
        self.filename = filename
        if self.hdf5 and os.path.exists(filename):
            os.remove(filename)

    # -- low level ---------------------------------------------------------
    @staticmethod
    def _scalar(value):
        if isinstance(value, (bool, int, np.integer)):
            return np.array([value], dtype=np.int32)
        return np.array([value], dtype=np.float64)

    def _set(self, path, data, attrs=None):
        with h5py.File(self.filename, 'a') as f:
            if path in f:
                del f[path]
            d = f.create_dataset(path, data=data)
            for k, v in (attrs or {}).items():
                d.attrs[k] = v

    def _strings(self, path, values):
        data = np.array([str(v) for v in values],
                        dtype=h5py.string_dtype(encoding='utf-8'))
        self._set(path, data)

    # -- payload that also exists as a legacy file -------------------------
    def array(self, path, arr, legacy=None, attrs=None):
        data = np.ascontiguousarray(arr, dtype=np.float64)
        if self.hdf5:
            self._set(path, data, attrs)
        elif legacy is not None:
            data.tofile(legacy)

    def int_array(self, path, arr, legacy=None, attrs=None):
        data = np.ascontiguousarray(arr, dtype=np.int32)
        if self.hdf5:
            self._set(path, data, attrs)
        elif legacy is not None:
            data.tofile(legacy)

    def aux(self, names, values, legacy, state_group=None, legacy_text=None):
        """auxdata*.txt -> one shape-(1,) dataset per field under REF/.

        ECAS/EROHF is per-state, so it goes into the state group when one is
        given; everything else is state-independent and lands in REF/."""
        if not self.hdf5:
            with open(legacy, 'w') as f:
                f.write(legacy_text if legacy_text is not None
                        else '\n'.join(map(str, values)))
            return
        for name, value in zip(names, values):
            target = 'REF'
            if name in ('ECAS', 'EROHF') and state_group is not None:
                target = state_group
            self._set('%s/%s' % (target, name), self._scalar(value))

    def text(self, legacy, build):
        """Write the ASCII-only rdm2.dat payload (legacy mode only).

        In HDF5 mode the 2-RDM is already stored in the file (RDM2 and
        RDM2_FULL_REORDERED) and GAMMCOR writes rdm2.dat itself from there,
        after the MO->NO transformation, so an ASCII copy dumped here would
        be redundant and - for non-natural orbitals - in the wrong basis."""
        if self.hdf5:
            return
        with open(legacy, 'w') as f:
            f.write(build())

    # -- metadata and extras: HDF5 only, no legacy counterpart -------------
    def provenance(self, method):
        if not self.hdf5:
            return
        with h5py.File(self.filename, 'a') as f:
            f.attrs['FORMAT'] = 'GAMMCOR-PYSCF'
            f.attrs['VERSION'] = 1
            f.attrs['GENERATOR'] = 'Pyscf2Gammcor.py'
            f.attrs['DATE'] = time.strftime('%Y-%m-%d %H:%M:%S')
            f.attrs['METHOD'] = method

    def system(self, mol):
        if not self.hdf5:
            return
        self._set('SYSTEM/GEOMETRY',
                  np.ascontiguousarray(mol.atom_coords(), dtype=np.float64))
        self._set('SYSTEM/ATOM_CHARGES',
                  np.asarray(mol.atom_charges(), dtype=np.int32))
        self._strings('SYSTEM/ATOM_SYMBOLS',
                      [mol.atom_symbol(a) for a in range(mol.natm)])
        self._set('SYSTEM/ENERGY_NUCLEAR', self._scalar(float(mol.energy_nuc())))
        with h5py.File(self.filename, 'a') as f:
            g = f.require_group('SYSTEM')
            g.attrs['CHARGE'] = int(mol.charge)
            g.attrs['SPIN2'] = int(mol.spin)
            g.attrs['NELECTRON'] = int(mol.nelectron)
            g.attrs['NATM'] = int(mol.natm)
            g.attrs['POINT_GROUP'] = str(mol.groupname)
            g.attrs['BASIS_NAME'] = str(mol.basis)
            g.attrs['UNIT'] = 'bohr'

    def ao_extras(self, mol):
        """OVERLAP/KINETIC/POTENTIAL: computed today and thrown away."""
        if not self.hdf5:
            return
        self._set('INTS/OVERLAP',
                  np.ascontiguousarray(mol.intor('int1e_ovlp'), dtype=np.float64))
        self._set('INTS/KINETIC',
                  np.ascontiguousarray(mol.intor('int1e_kin'), dtype=np.float64))
        self._set('INTS/POTENTIAL',
                  np.ascontiguousarray(mol.intor('int1e_nuc'), dtype=np.float64))
        self._strings('INTS/AO_LABELS', mol.ao_labels())
        with h5py.File(self.filename, 'a') as f:
            g = f.require_group('INTS')
            g.attrs['NBASIS'] = int(mol.nao_nr())
            g.attrs['ORDERING'] = 'PYSCF'

    def scf_reference(self, mol, myhf):
        if not self.hdf5:
            return
        self._set('SCF/MO_COEFF',
                  np.ascontiguousarray(np.asarray(myhf.mo_coeff).T, dtype=np.float64),
                  {'DIMS': ['mo', 'ao']})
        self._set('SCF/MO_ENERGY',
                  np.ascontiguousarray(myhf.mo_energy, dtype=np.float64))
        self._set('SCF/MO_OCC',
                  np.ascontiguousarray(myhf.mo_occ, dtype=np.float64))
        self._set('SCF/TOTAL_ENERGY', self._scalar(float(myhf.e_tot)))
        if mol.symmetry and mol.symm_orb is not None:
            try:
                self._strings('SCF/MO_IRREP', symm.label_orb_symm(
                    mol, mol.irrep_name, mol.symm_orb, myhf.mo_coeff))
            except Exception as err:
                print('note: SCF/MO_IRREP not written (%s)' % err)

    def orb_space(self, nbasis, NI, NA):
        """IndAux, written explicitly so Python and Fortran cannot disagree."""
        if not self.hdf5:
            return
        idx = np.full(int(nbasis), 2, dtype=np.int32)
        idx[:int(NI)] = 0
        idx[int(NI):int(NI) + int(NA)] = 1
        self._set('REF/ORB_SPACE', idx)

    def state_meta(self, group, label, state_number, energy, irrep, irrep_id,
                   weight=None, ms2_norm=None, nat_occ=None, no_transform=None,
                   active_irreps=None):
        if not self.hdf5:
            return
        with h5py.File(self.filename, 'a') as f:
            g = f.require_group(group)
            g.attrs['LABEL'] = str(label)
            g.attrs['STATE_NUMBER'] = int(state_number)
            g.attrs['ENERGY'] = float(energy)
            g.attrs['IRREP'] = str(irrep)
            g.attrs['IRREP_ID'] = int(irrep_id)
            if weight is not None:
                g.attrs['WEIGHT'] = float(weight)
            if ms2_norm is not None:
                g.attrs['MS2_NORM'] = float(ms2_norm)
        if nat_occ is not None:
            self._set('%s/NAT_OCC' % group,
                      np.ascontiguousarray(nat_occ, dtype=np.float64))
        if no_transform is not None:
            self._set('%s/NO_TRANSFORM' % group,
                      np.ascontiguousarray(no_transform, dtype=np.float64))
        if active_irreps is not None:
            self._strings('%s/ACTIVE_IRREPS' % group, active_irreps)

    def states_index(self, labels, energies, weights=None):
        if not self.hdf5:
            return
        self._set('POSTHF/NSTATES', np.array([len(labels)], dtype=np.int32))
        self._set('POSTHF/STATE_ENERGY', np.asarray(energies, dtype=np.float64))
        self._strings('POSTHF/STATE_LABELS', labels)
        if weights is not None:
            self._set('POSTHF/STATE_WEIGHTS', np.asarray(weights, dtype=np.float64))


def check_dm_aaaa_minus_bbbb_norm(dm2_aaaa, dm2_bbbb):

    diff = dm2_aaaa - dm2_bbbb
    norm = np.sum(diff ** 2)
    return norm

def reorder_rdm(dm2):

    return dm2.transpose(0, 2, 1, 3)


def is_dmrgci(solver):
    """Recognize DMRGCI, including PySCF's dynamic state-average wrapper."""
    return any(cls.__name__ == 'DMRGCI' for cls in type(solver).__mro__)


def calc_full_dm2(dm2s, thresh=1e-8):
    """
    Build and dump the spin-traced 2-RDM from spin blocks:
      dm2s = (dm2_aaaa, dm2_abab, dm2_bbbb)
    Returns all entries in raw PySCF ordering (i j k l value), as expected
    by the legacy AB_CAS_FOFO rdm2.dat reader.
    """
    dm2_aaaa, dm2_abab, dm2_bbbb = dm2s
    # build the β→α block by permuting indices
    dm2_ba = dm2_abab.transpose(2, 3, 0, 1)

    n = dm2_aaaa.shape[0]
    print('nnn', n)
    content = ""
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for l in range(n):
                    val = (
                        dm2_aaaa[i, j, k, l]
                        + dm2_abab[i, j, k, l]
                        + dm2_ba[i, j, k, l]
                        + dm2_bbbb[i, j, k, l]
                    )
                    a = dm2_aaaa[i, j, k, l]
                    b = dm2_abab[i, j, k, l]
                    c= dm2_ba[i, j, k, l]
                    d=dm2_bbbb[i, j, k, l]

                        # use 1-based indices for output
                    line = f"{i+1:4d} {j+1:4d} {k+1:4d} {l+1:4d} {val:19.12f}\n"
                    content += line
    return content


def calc_full_dm2_nospin(dm2):
    """Format a spin-traced 2-RDM in raw PySCF order for rdm2.dat."""
    n = dm2.shape[0]
    return ''.join(
        f"{i+1:4d} {j+1:4d} {k+1:4d} {l+1:4d} {dm2[i,j,k,l]:19.12f}\n"
        for i in range(n) for j in range(n) for k in range(n) for l in range(n)
    )

def get_irrep_labels(irrep, point_group):
    irrep_map = {
        "D2h": ["Ag", "B3u", "B2u", "B1g", "B1u", "B2g", "B3g", "Au"],
        "C2v": ["A1", "B1", "B2", "A2"],
        "C2h": ["Ag", "Au", "Bu", "Ag"],
        "D2": ["A", "B3", "B2", "B1"],
        "Cs": ["A'", 'A"'],
        "C2": ["A", "B"],
        "Ci": ["Ag", "Au"],
        "C1": ['A']
    }

    irreps = irrep_map[point_group]
    return irreps.index(irrep) + 1


def get_data_for_gammcor(mol, myhf, mycas, mymp = None, dump_eri=False, simple = True,
                         dump_hdf5 = True):
    """Process and export CASSCF calculation data for GAMMCOR.

    dump_hdf5 = True  -> one organized HDF5 file (H5FILE), see _Dump
    dump_hdf5 = False -> the legacy loose .bin/.txt files, exactly as before
    """

    storage_name = H5FILE
    lll = 180
    c = "-"
    script_name = sys.argv[0].split('/')[-1].replace('.py', '')

    state_av = hasattr(mycas, 'weights') and mycas.weights is not None
    dump = _Dump(dump_hdf5)
    dump.provenance('SA-CASSCF' if state_av else 'CASSCF')
    dump.system(mol)
    dump.ao_extras(mol)
    dump.scf_reference(mol, myhf)
    state_tags = []
    state_energies = []

    if hasattr(mycas, 'weights') and mycas.weights is not None:
        print(" State-Averaged CASSCF Analysis ".center(lll, f"{c}"))

        # Check for incompatible settings
        # if mycas.natorb:
        #     print("ERROR: Natural orbitals with SA-CASSCF not supported!")
        #     print("Please correct your input configuration.\n")
        #     sys.exit(0)

        print(f"State averaging weights: {mycas.weights}\n")

        # Get state symmetries
        if hasattr(mycas.fcisolver, 'fcisolvers'):
            solvers = mycas.fcisolver.fcisolvers
            state_symmetries = []
            for solver in solvers:
                state_symmetries.extend([solver.wfnsym] * solver.nroots)
        else:
            solver = mycas.fcisolver
            state_symmetries = [solver.wfnsym] * solver.nroots

        print("State-by-State Analysis:")
        print("-" * 70)
        print(f"{'State':^8} {'Symmetry':^12} {'Energy':^20} {'State #':^8} {'Sym #':^8}")
        print("-" * 70)

        is_dmrgscf = is_dmrgci(mycas.fcisolver)
        is_dmrgscf = False
        if is_dmrgscf:
            dm1_states, dm2_states = mycas.fcisolver.states_make_rdm12(
                mycas.ci, mycas.ncas, mycas.nelecas)
        else:
            dm1s_states, dm2s_states = mycas.fcisolver.states_make_rdm12s(
                mycas.ci, mycas.ncas, mycas.nelecas)

        # Process each state
        symmetry_state_numbers = {}
        for i, (symmetry_in, energy) in enumerate(zip(state_symmetries, mycas.e_states)):
            if not mol.symmetry:
                symmetry = 'A'
            else:
                symmetry = symmetry_in
            count = symmetry_state_numbers.get(symmetry, 0)
            symmetry_state_numbers[symmetry] = count + 1
            print('symmetry', symmetry)
            print('symmetry_state_numbers', symmetry_state_numbers)
            state_number = count + 1
            sym_idx = get_irrep_labels(symmetry, mol.groupname)
            state_tag = f'{state_number}.{sym_idx}'
            sgrp = f'POSTHF/STATES/{state_tag}'
            if mol.irrep_name == None:
                ir_sym = 'A'
            else:
                ir_sym = symmetry
            print(f"{i+1:^8} {symmetry:^12} {energy:^20.10f} {state_number:^8} {sym_idx:^8}")

            # Get density matrices
            is_dmrgscf = False
            if is_dmrgscf:
                dm1 = dm1_states[i]
                dm1a = dm1 / 2.0
                dm1b = dm1 / 2.0
                dm2_full_bin = dm2_states[i]
                norm = None
            else:
                dm1a = dm1s_states[0][i]
                dm1b = dm1s_states[1][i]
                dm1 = dm1a + dm1b
                dm2s = tuple(block[i] for block in dm2s_states)
                norm = check_dm_aaaa_minus_bbbb_norm(dm2s[0], dm2s[2])

            ncore, ncas = mycas.ncore, mycas.ncas
            print(f"\nDensity matrix (State {i+1}):")
            for j in range(ncas):
                for jj in range(ncas):
                    print(j, jj, dm1[j,jj])

            # Natural orbital analysis
            eval, evec = np.linalg.eig(dm1)
            idx = eval.argsort()[::-1]
            eval_sorted = eval[idx]
            evec_sorted = evec[:, idx]

            # Transform orbitals
            ncore, ncas = mycas.ncore, mycas.ncas
            nmo = mycas.mo_coeff.shape[1]
            
            full_transform = np.eye(nmo)
            active_slice = slice(ncore, ncore + ncas)
            full_transform[active_slice, active_slice] = evec_sorted

            CAONO = np.asarray(mycas.mo_coeff)
            transformed_CAONO = full_transform @ CAONO @ full_transform

            # Set occupations
            full_occ = np.zeros(nmo)
            full_occ[:ncore] = 2.0
            full_occ[active_slice] = eval_sorted

            # Save molden file
            if not simple:
                mol_nosym = mol.copy()
                mol_nosym.symmetry = False
                tools.molden.from_mo(mol_nosym, f"molden_{state_number}.{sym_idx}.inp", 
                                     transformed_CAONO, occ=full_occ)


                

            # Print orbital analysis
            if mol.symmetry and mol.symm_orb is not None:
                mo_irreps = symm.label_orb_symm(mol, mol.irrep_name, mol.symm_orb, mycas.mo_coeff)
                active_irreps = mo_irreps[ncore:ncore+ncas]
                irreps_sorted = active_irreps[idx]

                print(f"\nActive Space Natural Orbitals (State {i+1}):")
                print("-" * 50)
                print(f"{'MO#':^6} {'Occupancy':^15} {'Symmetry':^8}")
                print("-" * 50)
                ncore, ncas = mycas.ncore, mycas.ncas
                for j in range(ncas):
                    print(f"{j+ncore+1:^6} {eval_sorted[j]:^15.6f} {irreps_sorted[j]:^8}")
                print()

            # Check for MS≠0 states
            if norm is not None and norm > 1e-5:
                print("\nWARNING: Possible MS≠0 state detected!")
                print(f"Norm of (dm_aaaa - dm_bbbb): {norm:10.6f}\n")

            # Save density matrices
            is_dmrgscf = False
            if is_dmrgscf:
                dm2_full_reordered = reorder_rdm(dm2_full_bin)
            else:
                # Prepare 2-RDM components (Reorder: pyscf->gammcor)
                dm2_aaaa = reorder_rdm(dm2s[0])
                dm2_abab = reorder_rdm(dm2s[1])
                dm2_bbbb = reorder_rdm(dm2s[2])
                dm2_baba = dm2_abab.transpose(1, 0, 3, 2)

                # Calculate Spin Averages (Match state-specific logic) and dump
                dm2_spin_avg = (dm2_aaaa + dm2_bbbb) / 2.0
                dm2_ab_avg   = (dm2_abab + dm2_baba) / 2.0
                dump.array(f'{sgrp}/RDM2_AAAA', dm2_spin_avg, f'rdm2_aaaa_{state_tag}.bin',
                           attrs={'INDEX_ORDER': 'GAMMCOR', 'SPIN_AVERAGED': 1})
                dump.array(f'{sgrp}/RDM2_ABAB', dm2_ab_avg, f'rdm2_abab_{state_tag}.bin',
                           attrs={'INDEX_ORDER': 'GAMMCOR', 'SPIN_AVERAGED': 1})

                dm2_full_bin = 2.0 * dm2s[1] + dm2s[0] + dm2s[2]
                dm2_full_reordered = 2.0 * (dm2_spin_avg + dm2_ab_avg)

            dump.array(f'{sgrp}/RDM2', dm2_full_bin, f'rdm2_full_{state_tag}.bin',
                       attrs={'INDEX_ORDER': 'PYSCF'})
            dump.array(f'{sgrp}/RDM2_FULL_REORDERED', dm2_full_reordered,
                       f'rdm2_full_reordered_{state_tag}.bin',
                       attrs={'INDEX_ORDER': 'GAMMCOR'})

            # 6. Dump 1-RDM Binaries (Separate Spins + Total)
            dm1_full = dm1a + dm1b

            dump.array(f'{sgrp}/RDM1_A', dm1a, f'rdm1a_{state_tag}.bin')
            dump.array(f'{sgrp}/RDM1_B', dm1b, f'rdm1b_{state_tag}.bin')
            dump.array(f'{sgrp}/RDM1', dm1_full, f'rdm1_full_{state_tag}.bin')

            # Save auxiliary data
            nbasis = mol.nao_nr()
            NI = (mol.nelectron - sum(mycas.nelecas)) // 2
            auxdata = [nbasis, NI, ncas, nbasis-NI-ncas, energy, 
                      mol.energy_nuc(), mol.nelectron, int(mycas.natorb)]
            
            dump.aux(['NBASIS', 'NI', 'NA', 'NV', 'ECAS', 'ENUC', 'NEL', 'NATORB'],
                     auxdata, f'auxdata_{state_tag}.txt', state_group=sgrp)
            dump.orb_space(nbasis, NI, ncas)
            dump.state_meta(sgrp, state_tag, state_number, energy, ir_sym, sym_idx,
                            weight=mycas.weights[i], ms2_norm=norm,
                            nat_occ=eval_sorted, no_transform=evec_sorted,
                            active_irreps=(irreps_sorted if (mol.symmetry and
                                           mol.symm_orb is not None) else None))
            state_tags.append(state_tag)
            state_energies.append(energy)

        # Save additional data
        HCore = myhf.get_hcore()

        if not simple:
            tools.molden.from_scf(myhf, "moldenhf.inp")
            tools.molden.from_mcscf(mycas, "molden_mcscf.inp")
        dump.array('POSTHF/ORB_COEFF', CAONO.T, 'C.bin',
                   attrs={'TYPE': 'CANONICAL', 'DIMS': ['mo', 'ao']})
        dump.array('INTS/CORE_HAMILTONIAN', HCore.T, 'HCore.bin')
        dump.states_index(state_tags, state_energies, mycas.weights)


        if dump_eri:
            
            doMOtrans = False
            
            # This is version with MO transoformation done in pyscf
            start_time = time.time()

            
            eri_ao = mol.intor('int2e')

            # This is version with MO transoformation done in pyscf
            if doMOtrans:
                eri_mo = ao2mo.incore.full(eri_ao, CAONO)
                eri_mo = ao2mo.restore(1, eri_mo, CAONO.shape[1])
                n_orb = CAONO.shape[1]
            else:
                # this is version where integrals are not transformed - do it in gammcor
                t1 = time.time()
                eri_mo = ao2mo.restore(1, eri_ao, mol.nao) 
                n_orb = mol.nao
                t2 = time.time()
                print(f"\nTime spent on ERI transformation and processing: {t2 - t1:.2f} seconds")


            mid_time = time.time()
            n = n_orb
            nbasis = n_orb
            ninte1 = nbasis * (nbasis + 1) // 2
            ninte2 = ninte1 * (ninte1 + 1) // 2
            twono = np.zeros(ninte2)
                


            p, q, r, s = np.meshgrid(np.arange(n_orb), np.arange(n_orb), np.arange(n_orb), np.arange(n_orb), indexing='ij')


            def naddr3_vec(i1, i2, i3, i4):
                i1, i2, i3, i4 = i1+1, i2+1, i3+1, i4+1
                addr12 = (np.maximum(i1,i2) * (np.maximum(i1,i2)-1))//2 + np.minimum(i1,i2)
                addr34 = (np.maximum(i3,i4) * (np.maximum(i3,i4)-1))//2 + np.minimum(i3,i4)
                addr = (np.maximum(addr12,addr34) * (np.maximum(addr12,addr34)-1))//2 + np.minimum(addr12,addr34)
                return addr - 1  


            addrs = naddr3_vec(p, q, r, s)

            
            twono[addrs.ravel()] = eri_mo.ravel()                
            # Wymiary tablicy
            print(f"Wymiary (shape): {twono.shape}")

            # Rozmiar w bajtach
            print(f"Rozmiar w bajtach: {twono.nbytes}")

            # Typ danych i rozmiar pojedynczego elementu
            print(f"Typ danych: {twono.dtype}")
            print(f"Bajtów na element: {twono.itemsize}")

            # Rozmiar w MB
            print(f"Rozmiar w MB: {twono.nbytes / (1024*1024):.2f}")
            dump.array('MOINTS/ERI', twono, 'TWOEl.bin',
                       attrs={'PACKING': 'NADDR3', 'NORB': int(n_orb),
                              'BASIS': 'MO' if doMOtrans else 'AO',
                              'TRANSFORMED': int(doMOtrans)})

            
            end_time = time.time()
            print(f"\nTime spent on ERI transformation and processing: {end_time - start_time:.2f} seconds")
            print(f"\nTime spent on ERI transformation and processing: {end_time - mid_time:.2f} seconds")



    else:
        print(" State-specific CASSCF Analysis ".center(lll, f"{c}"))

        occ = mycas.mo_occ  
        occ = occ / 2.0

        CAONO = mycas.mo_coeff  # Natural orbital coefficients

        XOne = mol.intor('int1e_kin')+ mol.intor('int1e_nuc')
        HCore = myhf.get_hcore()

        np.set_printoptions(precision=6, suppress=True)
        print(HCore[:5, :5])

        state_tag = '1.1'
        sgrp = f'POSTHF/STATES/{state_tag}'
        dump.array('INTS/CORE_HAMILTONIAN', HCore.T, 'HCore.bin')

        is_dmrgscf = is_dmrgci(mycas.fcisolver)
        is_dmrgscf = False
        if is_dmrgscf:
            dm1, dm2_full_bin = mycas.fcisolver.make_rdm12(
                mycas.ci, mycas.ncas, mycas.nelecas)
            dm2_full_reordered = reorder_rdm(dm2_full_bin)
            dm1s = (dm1 / 2.0, dm1 / 2.0)
            norm = None
            dump.text('rdm2.dat', lambda: calc_full_dm2_nospin(dm2_full_bin))
        else:
            dm1s, dm2s = mycas.fcisolver.make_rdm12s(
                mycas.ci, mycas.ncas, mycas.nelecas)
            dm1 = dm1s[0] + dm1s[1]
            norm = check_dm_aaaa_minus_bbbb_norm(dm2s[0], dm2s[2])
            dm2_aaaa = reorder_rdm(dm2s[0])
            dm2_abab = reorder_rdm(dm2s[1])
            dm2_baba = dm2_abab.transpose(1,0,3,2)
            dm2_bbbb = reorder_rdm(dm2s[2])
            dump.text('rdm2.dat', lambda: calc_full_dm2(dm2s))
            dm2_spin_avg = (dm2_aaaa + dm2_bbbb) / 2.0
            dm2_ab_avg = (dm2_abab + dm2_baba) / 2.0
            _rdm2_attrs = {'INDEX_ORDER': 'GAMMCOR', 'SPIN_AVERAGED': 1}
            dump.array(f'{sgrp}/RDM2_AAAA', dm2_spin_avg, 'rdm2_aaaa.bin', attrs=_rdm2_attrs)
            dump.array(f'{sgrp}/RDM2_BBBB', dm2_spin_avg, 'rdm2_bbbb.bin', attrs=_rdm2_attrs)
            dump.array(f'{sgrp}/RDM2_ABAB', dm2_ab_avg, 'rdm2_abab.bin', attrs=_rdm2_attrs)
            dump.array(f'{sgrp}/RDM2_BABA', dm2_ab_avg, 'rdm2_baba.bin', attrs=_rdm2_attrs)
            dm2_full_bin = 2.0*dm2s[1] + dm2s[0] + dm2s[2]
            dm2_full_reordered = 2.0 * (dm2_spin_avg + dm2_ab_avg)
        dump.array(f'{sgrp}/RDM2', dm2_full_bin, 'rdm2_full.bin',
                   attrs={'INDEX_ORDER': 'PYSCF'})
        dump.array(f'{sgrp}/RDM2_FULL_REORDERED', dm2_full_reordered,
                   'rdm2_full_reordered.bin', attrs={'INDEX_ORDER': 'GAMMCOR'})

        nbasis = mol.nao_nr()
        print('nbasis', nbasis)
        NI = (mol.nelectron - sum(mycas.nelecas))//2
        NA = mycas.ncas
        NV = nbasis - NI-NA
        casscf_energy = mycas.e_tot
        Enuc = mol.energy_nuc()
        NEL = mol.nelectron

        dm1_spin_ava = (dm1s[0]+dm1s[1])/2.0
        dm1_spin_avb = (dm1s[0]+dm1s[1])/2.0
        dump.array(f'{sgrp}/RDM1_A', dm1_spin_ava, 'rdm1a.bin',
                   attrs={'SPIN_AVERAGED': 1})
        dump.array(f'{sgrp}/RDM1_B', dm1_spin_avb, 'rdm1b.bin',
                   attrs={'SPIN_AVERAGED': 1})
        dump.array(f'{sgrp}/RDM1', dm1s[0] + dm1s[1])

        rdm1p = np.zeros(nbasis)
        rdm1m = np.zeros(nbasis)
        rdm1p[:NI] = 1.0
        rdm1m[:NI] = 1.0
        rdm1p[NI:NI+NA] = np.diag(dm1s[0])
        rdm1m[NI:NI+NA] = np.diag(dm1s[1])
        rdm1p[NI+NA:] = 0.0
        rdm1m[NI+NA:] = 0.0
        dump.array(f'{sgrp}/OCC_A', rdm1p, 'rdm1p.bin')
        dump.array(f'{sgrp}/OCC_B', rdm1m, 'rdm1m.bin')
        
        if (mycas.natorb == True):
            natorb = 1
            print("Natural orbitals used. mycas.natorb:", mycas.natorb)
        else:
            natorb = 0
            print("Natural orbitals NOT used. mycas.natorb=", mycas.natorb)
            
        frozen = 0
        if mymp is not None:
            frozen = mymp.frozen
            
        dump.aux(['NBASIS', 'NI', 'NA', 'NV', 'ECAS', 'ENUC', 'NEL', 'NATORB', 'FROZEN'],
                 [nbasis, NI, NA, NV, casscf_energy, Enuc, NEL, natorb, frozen],
                 'auxdata.txt', state_group=sgrp,
                 legacy_text=f"{nbasis}\n{NI}\n{NA}\n{NV}\n{casscf_energy}\n{Enuc}\n{NEL}\n{natorb}\n{frozen}\n")
        dump.orb_space(nbasis, NI, NA)
        dump.state_meta(sgrp, state_tag, 1, casscf_energy, 'A', 1,
                        ms2_norm=norm,
                        nat_occ=np.diag(dm1s[0]) + np.diag(dm1s[1]))
        dump.states_index([state_tag], [casscf_energy])

        atom_coords = mol.atom_coords()  # Shape: (n_atoms, 3)
        n_atoms = mol.natm
        if not simple:
            tools.molden.from_scf(myhf, "moldenhf.inp")
            tools.molden.from_mcscf(mycas, "molden_mcscf.inp")
        dump.array('POSTHF/ORB_COEFF', CAONO.T, 'C.bin',
                   attrs={'TYPE': 'NATURAL' if natorb else 'CANONICAL',
                          'DIMS': ['mo', 'ao']})
        dump.array('POSTHF/OCC', occ, 'rdm1.bin')
        dump.array('INTS/CORE_HAMILTONIAN', HCore.T, 'HCore.bin')

        if mol.symmetry and mol.symm_orb is not None:
            mo_irreps = symm.label_orb_symm(mol, mol.irrep_name, mol.symm_orb, mycas.mo_coeff)
        
            ncore = mycas.ncore
            ncas = mycas.ncas

            active_irreps = mo_irreps[ncore:ncore+ncas]
            dump.state_meta(sgrp, state_tag, 1, casscf_energy, 'A', 1,
                            active_irreps=active_irreps)

            occ = mycas.mo_occ
            
            print("\n=== Active Space Orbitals ===\n")
            print("{:<6} {:<12} {:<6}".format("MO#", "Occupancy", "Irrep"))
            for j in range(ncas):
                print("{:<6} {:<18.9f} {:<6}".format(
                    j + ncore + 1,
                    occ[j+ncore],
                    active_irreps[j]))
        else:
            # Print without symmetry labels
            print("\n===  Active Space Orbitals ===\n")
            print("{:<6} {:<12}".format("MO#", "Occupancy"))
            ncas = mycas.ncas
            ncore = mycas.ncore

            for j in range(ncas):
                print("{:<6} {:<18.9f}".format(
                    j + ncore + 1,
                    occ[j+ncore]))



        if dump_eri:
            
            doMOtrans = True
            
            # This is version with MO transoformation done in pyscf
            start_time = time.time()

            
            eri_ao = mol.intor('int2e')

            # This is version with MO transoformation done in pyscf
            if doMOtrans:
                eri_mo = ao2mo.incore.full(eri_ao, CAONO)
                eri_mo = ao2mo.restore(1, eri_mo, CAONO.shape[1])
                n_orb = CAONO.shape[1]
            else:
                # this is version where integrals are not transformed - do it in gammcor
                t1 = time.time()
                eri_mo = ao2mo.restore(1, eri_ao, mol.nao) 
                n_orb = mol.nao
                t2 = time.time()
                print(f"\nTime spent on ERI transformation and processing: {t2 - t1:.2f} seconds")


            mid_time = time.time()
            n = n_orb
            nbasis = n_orb
            ninte1 = nbasis * (nbasis + 1) // 2
            ninte2 = ninte1 * (ninte1 + 1) // 2
            twono = np.zeros(ninte2)
                


            p, q, r, s = np.meshgrid(np.arange(n_orb), np.arange(n_orb), np.arange(n_orb), np.arange(n_orb), indexing='ij')


            def naddr3_vec(i1, i2, i3, i4):
                i1, i2, i3, i4 = i1+1, i2+1, i3+1, i4+1
                addr12 = (np.maximum(i1,i2) * (np.maximum(i1,i2)-1))//2 + np.minimum(i1,i2)
                addr34 = (np.maximum(i3,i4) * (np.maximum(i3,i4)-1))//2 + np.minimum(i3,i4)
                addr = (np.maximum(addr12,addr34) * (np.maximum(addr12,addr34)-1))//2 + np.minimum(addr12,addr34)
                return addr - 1  


            addrs = naddr3_vec(p, q, r, s)

            
            twono[addrs.ravel()] = eri_mo.ravel()                
            # Wymiary tablicy
            print(f"Wymiary (shape): {twono.shape}")

            # Rozmiar w bajtach
            print(f"Rozmiar w bajtach: {twono.nbytes}")

            # Typ danych i rozmiar pojedynczego elementu
            print(f"Typ danych: {twono.dtype}")
            print(f"Bajtów na element: {twono.itemsize}")

            # Rozmiar w MB
            print(f"Rozmiar w MB: {twono.nbytes / (1024*1024):.2f}")
            dump.array('MOINTS/ERI', twono, 'TWOEl.bin',
                       attrs={'PACKING': 'NADDR3', 'NORB': int(n_orb),
                              'BASIS': 'MO' if doMOtrans else 'AO',
                              'TRANSFORMED': int(doMOtrans)})

            
            end_time = time.time()
            print(f"\nTime spent on ERI transformation and processing: {end_time - start_time:.2f} seconds")
            print(f"\nTime spent on ERI transformation and processing: {end_time - mid_time:.2f} seconds")


def get_rohf_for_gammcor(mol, myhf, dump_hdf5 = True):
    """ROHF/MP2 export. Same quantities and the same conventions as
    get_rohf_for_gammcor in Pyscf2Gammcor_legacy.py; dump_hdf5 only changes where
    they land.

    dump_hdf5 = True  -> appended to H5FILE (SCF/, INTS/, REF/, POSTHF/)
    dump_hdf5 = False -> mo_occ_int.bin, occ_rohf.bin, auxdata_rohf.txt,
                         HCore.bin, C.bin, rdm2_rohf.dat
    """
    dump = _Dump(dump_hdf5)
    dump.provenance('ROHF')
    dump.system(mol)
    dump.ao_extras(mol)
    dump.scf_reference(mol, myhf)

    hf_occ = myhf.mo_occ

    occ_integers = np.zeros(len(hf_occ), dtype=np.int32)
    for i, occ in enumerate(hf_occ):
        if occ == 2.0:
            occ_integers[i] = 2
        elif occ == 1.0:
            occ_integers[i] = 1
        else:  # occ == 0.0
            occ_integers[i] = 0

    dump.int_array('SCF/MO_OCC_INT', occ_integers, 'mo_occ_int.bin')
    print(f"Zapisano tablice okupacji: {occ_integers}")

    nbasis = mol.nao_nr()
    print('nbasis', nbasis)
    NI = np.sum(hf_occ == 2.0)
    NA = np.sum(hf_occ == 1.0)
    NV = nbasis - NI - NA
    NOccup = NI + NA

    ROHF = myhf.e_tot
    Enuc = mol.energy_nuc()
    NEL = mol.nelectron

    print(f"NI (podwojnie obsadzone): {NI}")
    print(f"NA (pojedynczo obsadzone): {NA}")
    print(f"NV (wirtualne): {NV}")

    occ = myhf.mo_occ
    CAONO = myhf.mo_coeff
    HCore = mol.intor('int1e_kin') + mol.intor('int1e_nuc')

    state_tag = '1.1'
    sgrp = f'POSTHF/STATES/{state_tag}'

    dump.array('INTS/CORE_HAMILTONIAN', HCore.T, 'HCore.bin')
    dump.array('POSTHF/ORB_COEFF', CAONO.T, 'C.bin',
               attrs={'TYPE': 'CANONICAL', 'DIMS': ['mo', 'ao']})
    dump.array('SCF/MO_OCC', occ, 'occ_rohf.bin')

    dump.aux(['NBASIS', 'NI', 'NA', 'NV', 'EROHF', 'ENUC', 'NEL'],
             [nbasis, int(NI), int(NA), int(NV), ROHF, Enuc, NEL],
             'auxdata_rohf.txt', state_group=sgrp,
             legacy_text=f"{nbasis}\n{NI}\n{NA}\n{NV}\n{ROHF}\n{Enuc}\n{NEL}\n")
    dump.orb_space(nbasis, NI, NA)
    dump.state_meta(sgrp, state_tag, 1, ROHF, 'A', 1)
    dump.states_index([state_tag], [ROHF])

    # ===================================================================
    # Efektywne obsadzenie alfa to 0.5 dla wszystkich NA orbitali
    occ_alpha_eff = np.full(NA, 0.5)
    # Efektywne obsadzenie beta to 0.0 dla wszystkich NA orbitali
    occ_beta_eff = np.full(NA, 0.0)

    # ===================================================================
    # Krok 2: Zbuduj 2-RDM. Konwencja dokladnie jak w Pyscf2Gammcor_legacy.py:
    #         X[I,J,K,L] = 2 * Gamma_JL,IK
    # ===================================================================
    rdm2_rohf = np.zeros((NA, NA, NA, NA))
    for i in range(NA):
        for j in range(NA):
            for k in range(NA):
                for l in range(NA):
                    gamma_JL_IK = 0.0

                    # Czesc kulombowska
                    if j == i and l == k:
                        gamma_JL_IK += (occ_alpha_eff[j] + occ_beta_eff[j]) * \
                                       (occ_alpha_eff[l] + occ_beta_eff[l])

                    # Czesc wymienna
                    if j == k and l == i:
                        gamma_JL_IK -= (occ_alpha_eff[j] * occ_alpha_eff[l] +
                                        occ_beta_eff[j] * occ_beta_eff[l])

                    rdm2_rohf[i, j, k, l] = 2.0 * gamma_JL_IK

    if dump_hdf5:
        dump.array(f'{sgrp}/RDM2', rdm2_rohf,
                   attrs={'INDEX_ORDER': 'ROHF_EFF'})
    else:
        with open("rdm2_rohf.dat", 'w') as f:
            for I_fortran in range(1, NA + 1):
                for J_fortran in range(1, NA + 1):
                    for K_fortran in range(1, NA + 1):
                        for L_fortran in range(1, NA + 1):
                            X = rdm2_rohf[I_fortran-1, J_fortran-1,
                                          K_fortran-1, L_fortran-1]
                            f.write(f"{I_fortran:5d}{J_fortran:5d}{K_fortran:5d}"
                                    f"{L_fortran:5d}    {X:20.12f}\n")
