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
        occupancy = mycas.mo_occ[i]  
        print("{:<6} {:<6} {:<15.6f} {:<20.15f}".format(i+1, irrep, energy, occupancy))


        print("  Significant Atomic Orbital Contributions:")
        significant_aos = mo_coeff[:, i][abs(mo_coeff[:, i]) > 0.2] 
        for ao_idx, coef in enumerate(mo_coeff[:, i]):
            if abs(coef) > 0.2:
                ao_label = mol.ao_labels()[ao_idx] 
                print(f"    {ao_label}: {coef:.3f}")
        print()


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



def get_data_for_gammcor(mol, myhf, mycas, dump_eri=False):
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
            dm2_full = calc_full_dm((dm2s[0][i], dm2s[1][i]))

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
        HCore.T.astype(np.float64).tofile('HCore.bin')



        dm1s, dm2s = mycas.fcisolver.make_rdm12s(mycas.ci, mycas.ncas, mycas.nelecas)


        dm1 = dm1s[0]+dm1s[1]
        
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
        dm2_full = calc_full_dm((dm2s[0], dm2s[1]))
        #dm2_full = calc_full_dm(dm2s)
        
        with open('rdm2.dat', 'w') as f:
            f.write(dm2_full)
            
        # dm2_aaaa = np.asfortranarray(dm2_aaaa)
        # dm2_abab = np.asfortranarray(dm2_abab)
    
        dm2_aaaa.tofile('rdm2_aaaa.bin')
        dm2_abab.tofile('rdm2_abab.bin')

        dm2_full_bin = 2.0*(dm2s[0]+dm2s[1])
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
        if (mycas.natorb == True):
            natorb = 1
            print("Natural orbitals used. mycas.natorb:", mycas.natorb)
        else:
            natorb = 0
            print("Natural orbitals NOT used. mycas.natorb=", mycas.natorb)        

        with open('auxdata.txt', 'w') as f:
            f.write(f"{nbasis}\n{NI}\n{NA}\n{NV}\n{casscf_energy}\n{Enuc}\n{NEL}\n{natorb}\n")

        # Extract atom coordinates
        atom_coords = mol.atom_coords()  # Shape: (n_atoms, 3)
        n_atoms = mol.natm
        tools.molden.from_scf(myhf, "moldenhf.inp")
        tools.molden.from_mcscf(mycas, "molden_mcscf.inp")
        CAONO.T.astype(np.float64).tofile('C.bin')
        occ.astype(np.float64).tofile('rdm1.bin')
        HCore.T.astype(np.float64).tofile('HCore.bin')

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
                                
                            
            twono.astype(np.float64).tofile('TWOEl.bin')

            
            end_time = time.time()
            print(f"\nTime spent on ERI transformation and processing: {end_time - start_time:.2f} seconds")
            print(f"\nTime spent on ERI transformation and processing: {end_time - mid_time:.2f} seconds")



def analyze_configurations(mol, myhf, mycas):
    """
    Analyze and print dominant electronic configurations for each state in a CASSCF calculation.
    """


    is_sa = hasattr(mycas.fcisolver, 'nroots') and mycas.fcisolver.nroots > 1
    
    if is_sa:
        n_states = mycas.fcisolver.nroots
        print(f"\nState-Averaged CASSCF: Number of States = {n_states}\n")
    else:
        n_states = 1
        print(f"\nState-Specific CASSCF: Single State\n")
    
    mo_coeff = mycas.mo_coeff  # Shape: (num_AOs, num_MOs)
    mo_energies = mycas.mo_energy
    mo_irreps = symm.label_orb_symm(mol, mol.irrep_name, mol.symm_orb, mo_coeff)
    
    if is_sa:
        ci_vectors = mycas.fcisolver.ci  # Shape: (n_states, n_configs)
    else:
        ci_vectors = [mycas.fcisolver.ci]  
    
    for state in range(n_states):
        print(f"--- State {state + 1} ---")
        if is_sa:
            wfnsym = mycas.fcisolver.wfnsym[state]
        else:
            wfnsym = mycas.fcisolver.wfnsym
        print(f"Symmetry: {wfnsym}")
        print(f"Energy: {mycas.e_tot:.6f} Hartree\n")  
        
        ci = ci_vectors[state]
        
        top_n = 3
        top_indices = np.argsort(np.abs(ci))[-top_n:][::-1]
        print(f"Top {top_n} Dominant Configurations:")
        
        for idx in top_indices:
            coef = ci[idx]


            if hasattr(mycas.fcisolver, 'get_config'):
                config = mycas.fcisolver.get_config(idx)
                print(f"    Configuration: {config}")
            else:
                print("    Configuration extraction not available.")
        print("\n")


