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

def get_one_indices(binary_string):
    return [i+1 for i, bit in enumerate(reversed(binary_string)) if bit == '1']


def analyze_ci_coeffes(mycas):

    ci_coeffs = mycas.ci
    print(type(ci_coeffs))

    if hasattr(mycas.fcisolver, 'fcisolvers'):
        solver = mycas.fcisolver
        nroots = solver.nroots
    else:
        nroots = 1

    for root in range(0, nroots):
        for state_index in range(0, 1):
            #state_index = 0
            if nroots > 1:
                ci_state=ci_coeffs[root][state_index]
            else:
                ci_state=ci_coeffs[state_index]

            N = 20
            sorted_indices = np.argsort(np.abs(ci_state))[::-1]
            top_configs = sorted_indices[:N]
            top_ci_values = ci_state[top_configs]
            number_of_active_orbitals = mycas.ncas

            print("Top configurations for state {}: {}".format(state_index, top_configs))
            print("Corresponding CI coefficients: {}".format(top_ci_values))
    
            for idx in top_configs:
                binary = bin(idx)[2:].zfill(number_of_active_orbitals)
                print(idx, ' ', binary)
                print(get_one_indices(binary))


def irrep_analyze(mol, myhf):
    mo_irreps = symm.label_orb_symm(mol, mol.irrep_name, mol.symm_orb, myhf.mo_coeff)
    total_mo_counts = Counter(mo_irreps)
    mo_occ = myhf.mo_occ
    occupied_indices = np.where(mo_occ > 0)[0]
    occupied_irreps = mo_irreps[occupied_indices]
    counts_occ = Counter(occupied_irreps)
    print(occupied_irreps)
    print(counts_occ)
    print(mo_irreps)
    for i in range(0, len(mo_irreps)):
        print(i+1, mo_irreps[i])


    print('mol_groupname', mol.groupname)
    if mol.groupname == 'D2h':
        order = ['Ag', 'B3u', 'B2u', 'B1g', 'B1u', 'B2g', 'B3g', 'Au']
    elif mol.groupname == 'C2v':
        order = ['A1', 'A2', 'B1', 'B2']
    elif mol.groupname == 'Cs':
        order = ["A'", 'A"']
    elif mol.groupname == 'C2h':
        order = ['Ag', 'Bg', 'Au', 'Bu']
    else:
        # Dla innych grup po prostu użyj wszystkich znalezionych irrepsów
        order = sorted(set(mo_irreps))
    
    print()
    print('Symmetry occupancy after HF')
    print()
    print('  '.join(f'{x:3}' for x in order))
    values = [str(counts_occ.get(key, 0)) for key in order]
    print('  '.join(f'{x:3}' for x in values))
    

    print('\nTotal number of orbitals:')
    total_values = [str(total_mo_counts.get(key, 0)) for key in order]
    print('  '.join(f'{x:3}' for x in total_values))
    
def irrep_analyze_mo(mol, mo, occ):
    mo_irreps = symm.label_orb_symm(mol, mol.irrep_name, mol.symm_orb, mo)
    total_mo_counts = Counter(mo_irreps)

    occupied_indices = np.where(occ > 0)[0]
    occupied_irreps = mo_irreps[occupied_indices]
    counts_occ = Counter(occupied_irreps)
    print(occupied_irreps)
    print(counts_occ)
    print(mo_irreps)
    for i in range(0, len(mo_irreps)):
        print(i+1, mo_irreps[i])


    print('mol_groupname', mol.groupname)
    if mol.groupname == 'D2h':
        order = ['Ag', 'B3u', 'B2u', 'B1g', 'B1u', 'B2g', 'B3g', 'Au']
    elif mol.groupname == 'C2v':
        order = ['A1', 'A2', 'B1', 'B2']
    elif mol.groupname == 'Cs':
        order = ["A'", 'A"']
    elif mol.groupname == 'C2h':
        order = ['Ag', 'Bg', 'Au', 'Bu']
    else:
        # Dla innych grup po prostu użyj wszystkich znalezionych irrepsów                                                                                                                                                                                                                 
        order = sorted(set(mo_irreps))

    print()
    print('Symmetry occupancy after mp2')
    print()
    print('  '.join(f'{x:3}' for x in order))
    values = [str(counts_occ.get(key, 0)) for key in order]
    print('  '.join(f'{x:3}' for x in values))


    print('\nTotal number of orbitals:')
    total_values = [str(total_mo_counts.get(key, 0)) for key in order]
    print('  '.join(f'{x:3}' for x in total_values))


def analysis_of_mo(mol, myhf):

    mo_irreps = symm.label_orb_symm(mol, mol.irrep_name, mol.symm_orb, myhf.mo_coeff)

    mo_coeff = myhf.mo_coeff
    num_mo = mo_coeff.shape[1]
    
    print("\n=== Molecular Orbital Analysis ===\n")
    print("{:<6} {:<6} {:<15} {:<12}".format("MO#", "Irrep", "Energy (Hartree)", "Occupancy"))

    for i in range(num_mo):
        irrep = mo_irreps[i]
        energy = myhf.mo_energy[i]
        occupancy = myhf.mo_occ[i]
        print("{:<6} {:<6} {:<15.6f} {:<12}".format(i+1, irrep, energy, occupancy))
        
        print("  Significant Atomic Orbital Contributions:")
        significant_aos = mo_coeff[:, i][abs(mo_coeff[:, i]) > 0.2]  # Threshold can be adjusted
        for ao_idx, coef in enumerate(mo_coeff[:, i]):
            if abs(coef) > 0.2:
                ao_label = mol.ao_labels()[ao_idx]  # e.g., 'C 2s', 'C 2px'
                print(f"    {ao_label}: {coef:.3f}")
        print()
    
    molden_filename = 'molden_hf.molden'
    molden.from_mo(mol, molden_filename, myhf.mo_coeff)
    print(f"Molecular orbitals have been exported to '{molden_filename}' for visualization.\n")
    

def analysis_of_mo_cas(mol, mycas):
    mo_irreps = symm.label_orb_symm(mol, mol.irrep_name, mol.symm_orb, mycas.mo_coeff)

    mo_coeff = mycas.mo_coeff
    num_mo = mo_coeff.shape[1]

    print("\n=== Molecular Orbital Analysis ===\n")
    print("{:<6} {:<6} {:<15} {:<12}".format("MO#", "Irrep", "Energy (Hartree)", "Occupancy"))

    for i in range(num_mo):
        irrep = mo_irreps[i]
        energy = mycas.mo_energy[i]
        occupancy = mycas.mo_occ[i]  # Ensure mo_occ is available
        print("{:<6} {:<6} {:<15.6f} {:<20.15f}".format(i+1, irrep, energy, occupancy))

        # Print significant AO contributions
        print("  Significant Atomic Orbital Contributions:")
        significant_aos = mo_coeff[:, i][abs(mo_coeff[:, i]) > 0.2]  # Adjust threshold if needed
        for ao_idx, coef in enumerate(mo_coeff[:, i]):
            if abs(coef) > 0.2:
                ao_label = mol.ao_labels()[ao_idx]  # e.g., 'C 2s', 'C 2px'
                print(f"    {ao_label}: {coef:.3f}")
        print()

    # Export MOs to Molden file for visualization
    molden_filename = 'molden_casscf.molden'
    molden.from_mo(mol, molden_filename, mycas.mo_coeff)
    print(f"Molecular orbitals have been exported to '{molden_filename}' for visualization.\n")
    
    


def check_dm_aaaa_minus_bbbb_norm(dm2_aaaa, dm2_bbbb):

    diff = dm2_aaaa - dm2_bbbb
    norm = np.sum(diff ** 2)
    return norm

def myocc(mf):
    mol = mf.mol
    orbsym = symm.label_orb_symm(mol, mol.irrep_id, mol.symm_orb, mf.mo_coeff)
    doccsym = np.array(orbsym)[mf.mo_occ==2]
    soccsym = np.array(orbsym)[mf.mo_occ==1]
    for ir,irname in zip(mol.irrep_id, mol.irrep_name):
        print('%s, double-occ = %d, single-occ = %d' %
              (irname, sum(doccsym==ir), sum(soccsym==ir)))

def reorder_rdm(dm2):

    return dm2.transpose(0, 2, 1, 3)


def calc_full_dm(dm2s):

    dm2_aaaa, dm2_abab = dm2s
    dim = dm2s[0].shape[0]

    content = ""

    for l in range(dim):
        for k in range(dim):
            for j in range(dim):
                for i in range(dim):

                    # aa = dm2s[0][i, j, k, l]
                    # ab = dm2s[1][i, j, k, l]

                    aa = dm2_aaaa[i, j, k, l]
                    ab = dm2_abab[i, j, k, l]

                    #dm2_reorer[j, i, l, k] = 2.0*(aa+ab)
                    value   = 2.0*(aa+ab)
                    if (abs(value)>1.e-8):
                        content += f"{j+1:>4d} {i+1:>4d} {l+1:>4d} {k+1:>4d} {value:>19.12f} \n"
                        #print(f"{j+1:>4d} {i+1:>4d} {l+1:>4d} {k+1:>4d} {value:>19.12f}")
    return content


def calc_full_dm_nospin(dm2):

    dim = dm2.shape[0]

    content = ""

    for l in range(dim):
        for k in range(dim):
            for j in range(dim):
                for i in range(dim):

                    aa = dm2[i, j, k, l]

                    value   = aa
                    if (abs(value)>1.e-8):
                        content += f"{j+1:>4d} {i+1:>4d} {l+1:>4d} {k+1:>4d} {value:>19.12f} \n"
    return content


def calc_full_dm_ms(dm2s):

    dm2_aaaa, dm2_abab, dm2_bbbb = dm2s
    dim = dm2s[0].shape[0]

    content = ""

    for l in range(dim):
        for k in range(dim):
            for j in range(dim):
                for i in range(dim):

                    # aa = dm2s[0][i, j, k, l]
                    # ab = dm2s[1][i, j, k, l]

                    aa = dm2_aaaa[i, j, k, l]
                    ab = dm2_abab[i, j, k, l]
                    bb = dm2_bbbb[i, j, k, l]

                    #dm2_reorer[j, i, l, k] = 2.0*(aa+ab)
                    value   = 2.0*(ab) + aa + bb
                    if (abs(value)>1.e-8):
                        content += f"{j+1:>4d} {i+1:>4d} {l+1:>4d} {k+1:>4d} {value:>19.12f} \n"
                        print(f"{j+1:>4d} {i+1:>4d} {l+1:>4d} {k+1:>4d} {value:>19.12f}")
    return content


def calc_full_dm2_3(dm2s, thresh=1e-8):
    """
    Build and dump the spin-traced 2-RDM from spin blocks:
      dm2s = (dm2_aaaa, dm2_abab, dm2_bbbb)
    Returns a string of nonzero entries in chemists' ordering (i j k l value).
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
                    print(i, j, k, l)
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

                    print(val, a, b, c, d)
                    if abs(val) > thresh:
                        # use 1-based indices for output
                        print(  f"{i+1:4d} {j+1:4d} {k+1:4d} {l+1:4d} {a:19.12f} {b:19.12f} {c:19.12f} {d:19.12f}")
                    line = f"{i+1:4d} {j+1:4d} {k+1:4d} {l+1:4d} {val:19.12f}\n"
                    content += line
                    #print(line, end="")
    return content
                        

import numpy as np
import time
import os
from pyscf import gto, scf, mcscf, ao2mo, symm, tools, fci
from pyscf.tools import mo_mapping
from collections import Counter
from pyscf.tools import molden
import sys
from copy import deepcopy
import basis_set_exchange as bse
from functools import reduce

def get_one_indices(binary_string):
    return [i+1 for i, bit in enumerate(reversed(binary_string)) if bit == '1']


def analyze_ci_coeffes(mycas):

    ci_coeffs = mycas.ci
    print(type(ci_coeffs))

    if hasattr(mycas.fcisolver, 'fcisolvers'):
        solver = mycas.fcisolver
        nroots = solver.nroots
    else:
        nroots = 1

    for root in range(0, nroots):
        for state_index in range(0, 1):
            #state_index = 0
            if nroots > 1:
                ci_state=ci_coeffs[root][state_index]
            else:
                ci_state=ci_coeffs[state_index]

            N = 20
            sorted_indices = np.argsort(np.abs(ci_state))[::-1]
            top_configs = sorted_indices[:N]
            top_ci_values = ci_state[top_configs]
            number_of_active_orbitals = mycas.ncas

            print("Top configurations for state {}: {}".format(state_index, top_configs))
            print("Corresponding CI coefficients: {}".format(top_ci_values))
    
            for idx in top_configs:
                binary = bin(idx)[2:].zfill(number_of_active_orbitals)
                print(idx, ' ', binary)
                print(get_one_indices(binary))


def irrep_analyze(mol, myhf):
    mo_irreps = symm.label_orb_symm(mol, mol.irrep_name, mol.symm_orb, myhf.mo_coeff)
    total_mo_counts = Counter(mo_irreps)
    mo_occ = myhf.mo_occ
    occupied_indices = np.where(mo_occ > 0)[0]
    occupied_irreps = mo_irreps[occupied_indices]
    counts_occ = Counter(occupied_irreps)
    print(occupied_irreps)
    print(counts_occ)
    print(mo_irreps)
    for i in range(0, len(mo_irreps)):
        print(i+1, mo_irreps[i])


    print('mol_groupname', mol.groupname)
    if mol.groupname == 'D2h':
        order = ['Ag', 'B3u', 'B2u', 'B1g', 'B1u', 'B2g', 'B3g', 'Au']
    elif mol.groupname == 'C2v':
        order = ['A1', 'A2', 'B1', 'B2']
    elif mol.groupname == 'Cs':
        order = ["A'", 'A"']
    elif mol.groupname == 'C2h':
        order = ['Ag', 'Bg', 'Au', 'Bu']
    else:
        # Dla innych grup po prostu użyj wszystkich znalezionych irrepsów
        order = sorted(set(mo_irreps))
    
    print()
    print('Symmetry occupancy after HF')
    print()
    print('  '.join(f'{x:3}' for x in order))
    values = [str(counts_occ.get(key, 0)) for key in order]
    print('  '.join(f'{x:3}' for x in values))
    

    print('\nTotal number of orbitals:')
    total_values = [str(total_mo_counts.get(key, 0)) for key in order]
    print('  '.join(f'{x:3}' for x in total_values))
    
def irrep_analyze_mo(mol, mo, occ):
    mo_irreps = symm.label_orb_symm(mol, mol.irrep_name, mol.symm_orb, mo)
    total_mo_counts = Counter(mo_irreps)

    occupied_indices = np.where(occ > 0)[0]
    occupied_irreps = mo_irreps[occupied_indices]
    counts_occ = Counter(occupied_irreps)
    print(occupied_irreps)
    print(counts_occ)
    print(mo_irreps)
    for i in range(0, len(mo_irreps)):
        print(i+1, mo_irreps[i])


    print('mol_groupname', mol.groupname)
    if mol.groupname == 'D2h':
        order = ['Ag', 'B3u', 'B2u', 'B1g', 'B1u', 'B2g', 'B3g', 'Au']
    elif mol.groupname == 'C2v':
        order = ['A1', 'A2', 'B1', 'B2']
    elif mol.groupname == 'Cs':
        order = ["A'", 'A"']
    elif mol.groupname == 'C2h':
        order = ['Ag', 'Bg', 'Au', 'Bu']
    else:
        # Dla innych grup po prostu użyj wszystkich znalezionych irrepsów                                                                                                                                                                                                                 
        order = sorted(set(mo_irreps))

    print()
    print('Symmetry occupancy after mp2')
    print()
    print('  '.join(f'{x:3}' for x in order))
    values = [str(counts_occ.get(key, 0)) for key in order]
    print('  '.join(f'{x:3}' for x in values))


    print('\nTotal number of orbitals:')
    total_values = [str(total_mo_counts.get(key, 0)) for key in order]
    print('  '.join(f'{x:3}' for x in total_values))


def analysis_of_mo(mol, myhf):

    mo_irreps = symm.label_orb_symm(mol, mol.irrep_name, mol.symm_orb, myhf.mo_coeff)

    mo_coeff = myhf.mo_coeff
    num_mo = mo_coeff.shape[1]
    
    print("\n=== Molecular Orbital Analysis ===\n")
    print("{:<6} {:<6} {:<15} {:<12}".format("MO#", "Irrep", "Energy (Hartree)", "Occupancy"))

    for i in range(num_mo):
        irrep = mo_irreps[i]
        energy = myhf.mo_energy[i]
        occupancy = myhf.mo_occ[i]
        print("{:<6} {:<6} {:<15.6f} {:<12}".format(i+1, irrep, energy, occupancy))
        
        print("  Significant Atomic Orbital Contributions:")
        significant_aos = mo_coeff[:, i][abs(mo_coeff[:, i]) > 0.2]  # Threshold can be adjusted
        for ao_idx, coef in enumerate(mo_coeff[:, i]):
            if abs(coef) > 0.2:
                ao_label = mol.ao_labels()[ao_idx]  # e.g., 'C 2s', 'C 2px'
                print(f"    {ao_label}: {coef:.3f}")
        print()
    
    molden_filename = 'molden_hf.molden'
    molden.from_mo(mol, molden_filename, myhf.mo_coeff)
    print(f"Molecular orbitals have been exported to '{molden_filename}' for visualization.\n")
    

def analysis_of_mo_cas(mol, mycas):
    mo_irreps = symm.label_orb_symm(mol, mol.irrep_name, mol.symm_orb, mycas.mo_coeff)

    mo_coeff = mycas.mo_coeff
    num_mo = mo_coeff.shape[1]

    print("\n=== Molecular Orbital Analysis ===\n")
    print("{:<6} {:<6} {:<15} {:<12}".format("MO#", "Irrep", "Energy (Hartree)", "Occupancy"))

    for i in range(num_mo):
        irrep = mo_irreps[i]
        energy = mycas.mo_energy[i]
        occupancy = mycas.mo_occ[i]  # Ensure mo_occ is available
        print("{:<6} {:<6} {:<15.6f} {:<20.15f}".format(i+1, irrep, energy, occupancy))

        # Print significant AO contributions
        print("  Significant Atomic Orbital Contributions:")
        significant_aos = mo_coeff[:, i][abs(mo_coeff[:, i]) > 0.2]  # Adjust threshold if needed
        for ao_idx, coef in enumerate(mo_coeff[:, i]):
            if abs(coef) > 0.2:
                ao_label = mol.ao_labels()[ao_idx]  # e.g., 'C 2s', 'C 2px'
                print(f"    {ao_label}: {coef:.3f}")
        print()

    # Export MOs to Molden file for visualization
    molden_filename = 'molden_casscf.molden'
    molden.from_mo(mol, molden_filename, mycas.mo_coeff)
    print(f"Molecular orbitals have been exported to '{molden_filename}' for visualization.\n")
    
    


def check_dm_aaaa_minus_bbbb_norm(dm2_aaaa, dm2_bbbb):

    diff = dm2_aaaa - dm2_bbbb
    norm = np.sum(diff ** 2)
    return norm

def myocc(mf):
    mol = mf.mol
    orbsym = symm.label_orb_symm(mol, mol.irrep_id, mol.symm_orb, mf.mo_coeff)
    doccsym = np.array(orbsym)[mf.mo_occ==2]
    soccsym = np.array(orbsym)[mf.mo_occ==1]
    for ir,irname in zip(mol.irrep_id, mol.irrep_name):
        print('%s, double-occ = %d, single-occ = %d' %
              (irname, sum(doccsym==ir), sum(soccsym==ir)))

def reorder_rdm(dm2):

    return dm2.transpose(0, 2, 1, 3)


def calc_full_dm_mp2(dm2):

    dim = dm2[0].shape[0]

    content = ""

    for l in range(dim):
        for k in range(dim):
            for j in range(dim):
                for i in range(dim):

                    
                    aa = dm2[i, j, k, l]

                    value   = aa
                    if (abs(value)>1.e-8):
                        content += f"{j+1:>4d} {i+1:>4d} {l+1:>4d} {k+1:>4d} {value:>19.12f} \n"
                        #print(f"{j+1:>4d} {i+1:>4d} {l+1:>4d} {k+1:>4d} {value:>19.12f}")

    return content

def trtr(rdm1, mo_coeff):
    nmo = mo_coeff.shape[1]
    rdm1 = reduce(np.dot, (mo_coeff.T, rdm1, mo_coeff))
    return rdm1


def get_irrep_labels(irrep, point_group):
    irrep_map = {
        "D2h": ["Ag", "B3u", "B2u", "B1g", "B1u", "B2g", "B3g", "Au"],
        "C2v": ["A1", "B1", "B2", "A2"],
        "C2h": ["Ag", "Au", "Bu", "Ag"],
        "D2": ["A", "B3", "B2", "B1"],
        "Cs": ["A'", 'A"'],
        "C2": ["A", "B"],
        "Ci": ["Ag", "Au"]
    }

    irreps = irrep_map[point_group]
    return irreps.index(irrep) + 1



def get_data_for_gammcor(mol, myhf, mycas, mymp = None, dump_eri=False, simple = True):

    lll = 180
    c = "-"
    """Process and export CASSCF calculation data for GAMMCOR."""
    script_name = sys.argv[0].split('/')[-1].replace('.py', '')

    if hasattr(mycas, 'weights') and mycas.weights is not None:
        print(" State-Averaged CASSCF Analysis ".center(lll, f"{c}"))

        # Check for incompatible settings
        if mycas.natorb:
            print("ERROR: Natural orbitals with SA-CASSCF not supported!")
            print("Please correct your input configuration.\n")
            sys.exit(0)

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

        tdm1_ij = []
        tdm2_ij = []
        norb = mycas.ncas
        nelec = mycas.nelecas 

        for i in range(0, solver.nroots):
            for j in range(0, solver.nroots):
                tdm1, tdm2 = solver.trans_rdm12(solver.ci[i], solver.ci[j], norb, nelec)
           #     print('type of tdm1', type(tdm1), tdm1.shape)
            #    print('type of tdm2', type(tdm2), tdm2.shape)
                np.set_printoptions(precision=4, suppress=True)  # Adjust precision and suppress scientific notation
             #   print("tdm1:")
              #  print(tdm1)
               # print("\ntdm2:")
                #print(tdm2)
                tdm1_ij.append(tdm1)
                tdm2_ij.append(tdm2)


        # Process each state
        symmetry_state_numbers = {}
        for i, (symmetry, energy) in enumerate(zip(state_symmetries, mycas.e_states)):
            count = symmetry_state_numbers.get(symmetry, 0)
            symmetry_state_numbers[symmetry] = count + 1
            state_number = count + 1
            sym_idx = get_irrep_labels(symmetry, mol.groupname)

            print(f"{i+1:^8} {symmetry:^12} {energy:^20.10f} {state_number:^8} {sym_idx:^8}")

            # Get density matrices
            dm1s, dm2s = mycas.fcisolver.states_make_rdm12s(mycas.ci, mycas.ncas, mycas.nelecas)
            dm1 = dm1s[0][i] + dm1s[1][i]
            dm2 = dm2s[0][i] + dm2s[1][i]

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
                for j in range(ncas):
                    print(f"{j+ncore+1:^6} {eval_sorted[j]:^15.6f} {irreps_sorted[j]:^8}")
                print()

            # Check for MS≠0 states
            norm = check_dm_aaaa_minus_bbbb_norm(dm2s[0][i], dm2s[2][i])
            if norm > 1e-5:
                print("\nWARNING: Possible MS≠0 state detected!")
                print(f"Norm of (dm_aaaa - dm_bbbb): {norm:10.6f}\n")

            # Save density matrices
            dm2_aaaa = reorder_rdm(dm2s[0][i])
            dm2_abab = reorder_rdm(dm2s[1][i])
#            dm2_full = calc_full_dm((dm2s[0][i], dm2s[1][i]))
            dm2_full = calc_full_dm_ms((dm2s[0][i], dm2s[1][i], dm2s[2][i]))

            # Write files
            dm2_aaaa.tofile(f'rdm2_aaaa_{state_number}.{sym_idx}.bin')
            dm2_abab.tofile(f'rdm2_abab_{state_number}.{sym_idx}.bin')

            
            dm1.tofile(f'rdm1_{state_number}.{sym_idx}.bin')
            
            with open(f'rdm2_{state_number}.{sym_idx}.dat', 'w') as f:
                f.write(dm2_full)

            dm2_full_bin = 2.0*(dm2s[0][i] + dm2s[1][i])
            dm2_full_bin.tofile(f'rdm2_{state_number}.{sym_idx}.bin')

            # Save auxiliary data
            nbasis = mol.nao_nr()
            NI = (mol.nelectron - sum(mycas.nelecas)) // 2
            auxdata = [nbasis, NI, ncas, nbasis-NI-ncas, energy, 
                      mol.energy_nuc(), mol.nelectron, int(mycas.natorb)]
            
            with open(f'auxdata_{state_number}.{sym_idx}.txt', 'w') as f:
                f.write('\n'.join(map(str, auxdata)))

        # Save additional data
        HCore = mol.intor('int1e_kin') + mol.intor('int1e_nuc')
        if not simple:
            tools.molden.from_scf(myhf, "moldenhf.inp")
        CAONO.T.astype(np.float64).tofile('C.bin')
        HCore.T.astype(np.float64).tofile('HCore.bin')

    else:
        print(" Single-State CASSCF Analysis ".center(lll, f"{c}"))

        occ = mycas.mo_occ  
        occ = occ / 2.0          
        #print(occ)
        CAONO = mycas.mo_coeff  # Natural orbital coefficients

        #XOne = mol.intor('int1e_kin')+ mol.intor('int1e_nuc')
        # Transform one-electron integrals to MO basis
        #print('xone',  XOne[0,0], XOne[1,1])
        #print('cao', CAONO[0,0], CAONO[0,1], CAONO[1,0])
        #HCore = np.dot(CAONO.T, np.dot(XOne, CAONO))
        #print('hcore', HCore[0,0], HCore[1,1])

        HCore = mol.intor('int1e_kin') + mol.intor('int1e_nuc')
        print("\nZ określoną precyzją:")
        np.set_printoptions(precision=6, suppress=True)
        print(HCore[:5, :5])
        HCore.T.astype(np.float64).tofile('HCore.bin')



        dm1s, dm2s = mycas.fcisolver.make_rdm12s(mycas.ci, mycas.ncas, mycas.nelecas)
        print(type(dm1s), 'type-1')
        print(dm1s[0].shape, 'shape-1')
        print("Macierz dm1s[0]:")
        np.set_printoptions(precision=6, suppress=True, linewidth=100)
        print(dm1s[0])

        stdm1, stdm2 = mycas.fcisolver.make_rdm12(mycas.ci, mycas.ncas, mycas.nelecas)
        print(type(stdm2), 'rr2-type')
        print(stdm2.shape, 'rr2-shape')

        rdm1, rdm2, rdm3= mycas.fcisolver.make_rdm123(mycas.ci, mycas.ncas, mycas.nelecas)
        print(rdm3.shape, 'rr')

        rdm1, rdm2, rdm3= mycas.fcisolver.make_rdm123s(mycas.ci, mycas.ncas, mycas.nelecas)
        print(type(rdm3), 'rr2')
        print(rdm3[0].shape, 'rr3')
        

        dm1 = dm1s[0]+dm1s[1]
        for i in range(0, len(dm1s[0])):
            print('pluszek', dm1s[0][i], dm1s[1][i])
        
        print('dm1dm1dm1dm1dm1')
        print(np.array_str(dm1, precision=2, suppress_small=True))


        norm = check_dm_aaaa_minus_bbbb_norm(dm2s[0], dm2s[2])
        print('norm of ||dm_aaaa-dm_bbbb|| is', norm)
        if (norm > 1e-5):
            print('**********************************')
            print()
            print('WARNING: state might have Ms !=0')
            print()
            print('**********************************')


        # print()
        
        # dmfull = dm1s[0]+dm1s[1]
        # print(np.array_str(dmfull, precision=2, suppress_small=True))
        # print()
        
        # tr1 = np.trace(dm1s[0])
        # print('tr1', tr1)
        # tr2 = np.trace(dm1s[1])
    
        # nmo = mycas.mo_coeff.shape[1]
        # casdm1, casdm2 = mycas.fcisolver.make_rdm12(mycas.ci, mycas.ncas, mycas.nelecas)
        # rdm1, rdm2 = mcscf.addons._make_rdm12_on_mo(casdm1, casdm2, mycas.ncore, mycas.ncas, nmo)
        
        # tr = np.trace(rdm1)
        # print(tr)
        # print(rdm1)
        
        # rdm2_full = reorder_rdm(rdm2)
    
        dm2_aaaa = reorder_rdm(dm2s[0])
        dm2_abab = reorder_rdm(dm2s[1])
        dm2_baba = dm2_abab.transpose(2,3,0,1)
        dm2_bbbb = reorder_rdm(dm2s[2])
        dm2_full_zzz = calc_full_dm((dm2s[0], dm2s[1]))
#        dm2_full = calc_full_dm_ms((dm2s[0], dm2s[1], dm2s[2]))

        dm2_full = calc_full_dm2_3(dm2s)
        #dm2_full = calc_full_dm(dm2s)

        st_rdm2 = reorder_rdm(stdm2)
        st_rdm2 = calc_full_dm_nospin(st_rdm2)

        if not simple:
            with open('st_rdm2.dat', 'w') as f:
                f.write(st_rdm2)
            print()
            #print('im here', dm2_full)
            print()
        with open('rdm2.dat', 'w') as f:
            f.write(dm2_full)

        if not simple:
            with open('rdm2_zzz.dat', 'w') as f:
                f.write(dm2_full_zzz)

        # dm2_aaaa = np.asfortranarray(dm2_aaaa)
        # dm2_abab = np.asfortranarray(dm2_abab)
    
        dm2_aaaa.tofile('rdm2_aaaa.bin')
        dm2_abab.tofile('rdm2_abab.bin')
        dm2_bbbb.tofile('rdm2_bbbb.bin')

        
        dm2_full_bin = 2.0*(dm2s[1])+dm2s[0]+dm2s[2]
        dm2_full_bin.tofile('rdm2_full.bin')

        # rdm2_full.tofile('rdm2_full.bin')

        nbasis = mol.nao_nr()
        print('nbasis', nbasis)
        NI = (mol.nelectron - sum(mycas.nelecas))//2
        NA = mycas.ncas
        NV = nbasis - NI-NA
        casscf_energy = mycas.e_tot
        Enuc = mol.energy_nuc()
        NEL = mol.nelectron


        rdm1p = np.zeros(NI + NA)
        rdm1m = np.zeros(NI + NA)
        rdm1p[:NI] = 1.0
        rdm1m[:NI] = 1.0
        rdm1p[NI:NI+NA] = np.diag(dm1s[0])
        rdm1m[NI:NI+NA] = np.diag(dm1s[1])
        if not simple:
            rdm1p.astype(np.float64).tofile('rdm1p.bin')
            rdm1m.astype(np.float64).tofile('rdm1m.bin')
        
        if (mycas.natorb == True):
            natorb = 1
            print("Natural orbitals used. mycas.natorb:", mycas.natorb)
        else:
            natorb = 0
            print("Natural orbitals NOT used. mycas.natorb=", mycas.natorb)

        frozen = 0
        if mymp is not None:
            frozen = mymp.frozen

        with open('auxdata.txt', 'w') as f:
            f.write(f"{nbasis}\n{NI}\n{NA}\n{NV}\n{casscf_energy}\n{Enuc}\n{NEL}\n{natorb}\n{frozen}\n")

        # Extract atom coordinates
        atom_coords = mol.atom_coords()  # Shape: (n_atoms, 3)
        n_atoms = mol.natm
        if not simple:
            tools.molden.from_scf(myhf, "moldenhf.inp")
            tools.molden.from_mcscf(mycas, "molden_mcscf.inp")
        CAONO.T.astype(np.float64).tofile('C.bin')
        occ.astype(np.float64).tofile('rdm1.bin')
        HCore.T.astype(np.float64).tofile('HCore.bin')

        NIA = NI + NA
        orbital_energies = myhf.mo_energy
        occupied_energies = orbital_energies[:NIA]
        virtual_energies = orbital_energies[NIA:]
        if not simple:
            print(f"Zajęte (core+active): {len(occupied_energies)} orbitali")
            print(f"Wirtualne: {len(virtual_energies)} orbitali")
            occupied_energies.astype(np.float64).tofile('eorbi.bin')
            virtual_energies.astype(np.float64).tofile('eorba.bin')

        if mol.symmetry and mol.symm_orb is not None:
            mo_irreps = symm.label_orb_symm(mol, mol.irrep_name, mol.symm_orb, mycas.mo_coeff)
        
            ncore = mycas.ncore
            ncas = mycas.ncas

            active_irreps = mo_irreps[ncore:ncore+ncas]

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
            idx = naddr3_vec(0, 0, 0, 13)
            print(f"twono({idx}) = {twono[idx]}")                                
            print(f"Liczba elementów: {twono.size}")

            # Wymiary tablicy
            print(f"Wymiary (shape): {twono.shape}")

            # Rozmiar w bajtach
            print(f"Rozmiar w bajtach: {twono.nbytes}")

            # Typ danych i rozmiar pojedynczego elementu
            print(f"Typ danych: {twono.dtype}")
            print(f"Bajtów na element: {twono.itemsize}")

            # Rozmiar w MB
            print(f"Rozmiar w MB: {twono.nbytes / (1024*1024):.2f}")
            twono.astype(np.float64).tofile('TWOEl.bin')

            
            end_time = time.time()
            print(f"\nTime spent on ERI transformation and processing: {end_time - start_time:.2f} seconds")
            print(f"\nTime spent on ERI transformation and processing: {end_time - mid_time:.2f} seconds")


import numpy as np

def get_rohf_for_gammcor(mol, myhf):
    hf_occ = myhf.mo_occ
    
    occ_integers = np.zeros(len(hf_occ), dtype=np.int32)
    for i, occ in enumerate(hf_occ):
        if occ == 2.0:
            occ_integers[i] = 2
        elif occ == 1.0:
            occ_integers[i] = 1
        else:  # occ == 0.0
            occ_integers[i] = 0
    
    occ_integers.astype(np.int32).tofile('mo_occ_int.bin')
    print(f"Zapisano tablicę okupacji do mo_occ_int.bin: {occ_integers}")
    
    nbasis = mol.nao_nr()
    print('nbasis', nbasis)
    NI = np.sum(hf_occ == 2.0)
    NA = np.sum(hf_occ == 1.0)
    NV = nbasis - NI - NA
    NOccup = NI + NA
    
    ROHF = myhf.e_tot
    Enuc = mol.energy_nuc()
    NEL = mol.nelectron
    
    print(f"NI (podwójnie obsadzone): {NI}")
    print(f"NA (pojedynczo obsadzone): {NA}")
    print(f"NV (wirtualne): {NV}")

    occ = myhf.mo_occ
    CAONO = myhf.mo_coeff  
    HCore = mol.intor('int1e_kin') + mol.intor('int1e_nuc')
    HCore.T.astype(np.float64).tofile('HCore.bin')
    CAONO.T.astype(np.float64).tofile('C.bin')
    occ.astype(np.float64).tofile('occ_rohf.bin')
    

    with open('auxdata_rohf.txt', 'w') as f:
        f.write(f"{nbasis}\n{NI}\n{NA}\n{NV}\n{ROHF}\n{Enuc}\n{NEL}\n")

         # ===================================================================
    # Efektywne obsadzenie alfa to 0.5 dla wszystkich NA orbitali
    occ_alpha_eff = np.full(NA, 0.5)
    # Efektywne obsadzenie beta to 0.0 dla wszystkich NA orbitali
    occ_beta_eff = np.full(NA, 0.0)

    # ===================================================================
    # Krok 2: Zbuduj 2-RDM i zapisz do pliku z przenumerowanymi indeksami
    # ===================================================================
    filename="rdm2_rohf.dat"
    with open(filename, 'w') as f:
        # Pętle po indeksach w przestrzeni aktywnej, od 1 do NA
        for I_fortran in range(1, NA + 1):
            for J_fortran in range(1, NA + 1):
                for K_fortran in range(1, NA + 1):
                    for L_fortran in range(1, NA + 1):
                        
                        # Indeksy Pythona wewnątrz podprzestrzeni (od 0 do NA-1)
                        i, j, k, l = I_fortran-1, J_fortran-1, K_fortran-1, L_fortran-1
                        
                        # Obliczamy element Γ_JL,IK używając EFEKTYWNYCH obsadzeń
                        gamma_JL_IK = 0.0
                        
                        # Część kulombowska
                        if j == i and l == k:
                            gamma_JL_IK += (occ_alpha_eff[j] + occ_beta_eff[j]) * (occ_alpha_eff[l] + occ_beta_eff[l])
                            
                        # Część wymienna
                        if j == k and l == i:
                            gamma_JL_IK -= (occ_alpha_eff[j] * occ_alpha_eff[l] + occ_beta_eff[j] * occ_beta_eff[l])
                            
                        # Zgodnie z komentarzem w kodzie Fortran, X = 2 * Γ_JL,IK
                        X = 2.0 * gamma_JL_IK
                        
                       # if abs(X) > 1e-9:
                            # Zapisujemy indeksy z pętli (od 1 do NA)
                        f.write(f"{I_fortran:5d}{J_fortran:5d}{K_fortran:5d}{L_fortran:5d}    {X:20.12f}\n")
            


            
# def get_mp2_for_gammcor(mol, myhf, mymp):

# #    natorb_threshold = 1.e-8
    
#     aaa = True
#     if aaa == True:
        
#         dm1 = mymp.make_rdm1()
#         dm2 = mymp.make_rdm2()
#         # nbasis = mol.n_ao_nr()
#         noons, natorbs = np.linalg.eigh(dm1)
#         #       noons = noons[::-1]  # sort descending                                                                                                             
#         #       print("Natural orbital occupation numbers:", noons)

#         # occ = np.zeros(nbasis))
#         # for i in range(0, nbasis):
#         #     if n

# #        occ = mymp.mo_occ  
# #        occ = occ / 2.0          

#         CAONO = mymp.mo_coeff  # Natural orbital coefficients

#         # HCore = mol.intor('int1e_kin') + mol.intor('int1e_nuc')
#         # HCore.T.astype(np.float64).tofile('HCore.bin')

#  #       dm1 = mymp.make_rdm1()
#  #       dm2 = mymp.make_rdm2()

#         dm2_full = reorder_rdm(dm2)
#         dm2_full = calc_full_dm_mp2(dm2_full)
        
#         with open('rdm2_mp2.dat', 'w') as f:
#             f.write(dm2_full)
            
        # dm2_full_bin = 2.0*(dm2s[0]+dm2s[1])
        # dm2_full_bin.tofile('rdm2_full.bin')

        # rdm2_full.tofile('rdm2_full.bin')

#         nbasis = mol.nao_nr()
#         print('nbasis', nbasis)
#         NI = (mol.nelectron )//2
#         NA = 0
#         NV = nbasis - NI-NA
#         casscf_energy = mymp.e_tot
#         Enuc = mol.energy_nuc()
#         NEL = mol.nelectron
#         #      if (mycas.natorb == True):
#         natorb = 1
# #        print("Natural orbitals used. mycas.natorb:")
#    #     else:
#     #        natorb = 0
#      #       print("Natural orbitals NOT used. mycas.natorb=", mycas.natorb)        

#         with open('auxdata.txt', 'w') as f:
#             f.write(f"{nbasis}\n{NI}\n{NA}\n{NV}\n{casscf_energy}\n{Enuc}\n{NEL}\n{natorb}\n")

#         # Extract atom coordinates
#         atom_coords = mol.atom_coords()  # Shape: (n_atoms, 3)
#         n_atoms = mol.natm
#         tools.molden.from_scf(myhf, "moldenhf.inp")
# #        toaols.molden.from_mcscf(mycas, "molden_mcscf.inp")
#         CAONO.T.astype(np.float64).tofile('C.bin')
#         occ.astype(np.float64).tofile('rdm1.bin')
#         HCore.T.astype(np.float64).tofile('HCore.bin')

#         if mol.symmetry and mol.symm_orb is not None:
#             mo_irreps = symm.label_orb_symm(mol, mol.irrep_name, mol.symm_orb, mycas.mo_coeff)
        
#             ncore = mycas.ncore
#             ncas = mycas.ncas

#             active_irreps = mo_irreps[ncore:ncore+ncas]

#             occ = mycas.mo_occ
            
#             print("\n=== Active Space Orbitals ===\n")
#             print("{:<6} {:<12} {:<6}".format("MO#", "Occupancy", "Irrep"))
#             for j in range(ncas):
#                 print("{:<6} {:<18.9f} {:<6}".format(
#                     j + ncore + 1,
#                     occ[j+ncore],
#                     active_irreps[j]))
#         else:
#             # Print without symmetry labels
#             print("\n===  Active Space Orbitals ===\n")
#             print("{:<6} {:<12}".format("MO#", "Occupancy"))
#             for j in range(ncas):
#                 print("{:<6} {:<18.9f}".format(
#                     j + ncore + 1,
#                     occ[j+ncore]))





def analyze_configurations(mol, myhf, mycas):
    """
    Analyze and print dominant electronic configurations for each state in a CASSCF calculation.
    """

    # Determine if the calculation is State-Averaged or State-Specific
    is_sa = hasattr(mycas.fcisolver, 'nroots') and mycas.fcisolver.nroots > 1
    
    if is_sa:
        n_states = mycas.fcisolver.nroots
        print(f"\nState-Averaged CASSCF: Number of States = {n_states}\n")
    else:
        n_states = 1
        print(f"\nState-Specific CASSCF: Single State\n")
    
    # Get MO coefficients and energies
    mo_coeff = mycas.mo_coeff  # Shape: (num_AOs, num_MOs)
    mo_energies = mycas.mo_energy
    mo_irreps = symm.label_orb_symm(mol, mol.irrep_name, mol.symm_orb, mo_coeff)
    
    # Access CI coefficients
    if is_sa:
        ci_vectors = mycas.fcisolver.ci  # Shape: (n_states, n_configs)
    else:
        ci_vectors = [mycas.fcisolver.ci]  # List with one element
    
    for state in range(n_states):
        print(f"--- State {state + 1} ---")
        if is_sa:
            wfnsym = mycas.fcisolver.wfnsym[state]
        else:
            wfnsym = mycas.fcisolver.wfnsym
        print(f"Symmetry: {wfnsym}")
        print(f"Energy: {mycas.e_tot:.6f} Hartree\n")  # For SA, e_tot is averaged
        
        # Get CI coefficients for the state
        ci = ci_vectors[state]
        
        # Identify the dominant configurations (e.g., top 3)
        top_n = 3
        top_indices = np.argsort(np.abs(ci))[-top_n:][::-1]  # Indices of top_n configurations
        print(f"Top {top_n} Dominant Configurations:")
        
        for idx in top_indices:
            coef = ci[idx]

            # Attempt to retrieve the configuration bitstring
            if hasattr(mycas.fcisolver, 'get_config'):
                config = mycas.fcisolver.get_config(idx)
                print(f"    Configuration: {config}")
            else:
                print("    Configuration extraction not available.")
        print("\n")


