#!/bin/bash

source ~/Programs/qp2/quantum_package.rc

# geom in bohr (-a)
qp create_ezfio -b 6-31g "bebe.xyz" -a -o "cas_A"
qp set_file "cas_A"

qp set electrons elec_alpha_num 2
qp set electrons elec_beta_num  2

qp set nuclei nucl_label  "['Be','x']"
qp set nuclei nucl_charge "[4.0, 0.0]"

qp set davidson_keywords distributed_davidson False

qp run scf | tee ${EZFIO_FILE}.scf.A.out
qp set_frozen_core 
sleep 1

qp set determinants n_states 4

qp run cis | tee ${EZFIO_FILE}.cis.A.out
sleep 1
qp run save_natorb | tee ${EZFIO_FILE}.natorb.A.out

#CIS within the active space
qp set_mo_class -c "[1]" -a "[2-5]" -d "[6-18]"
qp run cis | tee ${EZFIO_FILE}.cis_2states_active_spaceA.out
qp set determinants read_wf True

qp set_mo_class -c "[1]" -a "[2-5]" -v "[6-18]"

qp set casscf_cipsi small_active_space True
qp set mol_properties calc_tr_dipole_moment True
qp run casscf | tee ${EZFIO_FILE}.casscf_2states.A.out

# Export TREXIO using gammcor_plugin
qp set gammcor_plugin trexio_file \"A.hdf5\"
qp set gammcor_plugin cholesky_rdm False
qp set gammcor_plugin cholesky_tolerance 1.e-5
qp set gammcor_plugin export_tr_rdm True
qp run export_gammcor >> "export_A.out"
qp run gammcor_plugin >> "export_A.out"

# monomer B

# geom in bohr (-b)
qp create_ezfio -b 6-31g  "bebe.xyz" -a -o "cas_B"
qp set_file "cas_B"

qp set electrons elec_alpha_num 2
qp set electrons elec_beta_num  2

qp set nuclei nucl_label  "['x', 'Be']"
qp set nuclei nucl_charge "[0.0, 4.0]"

qp set davidson_keywords distributed_davidson False

qp run scf | tee ${EZFIO_FILE}.scf.B.out
qp set_frozen_core 
sleep 1

qp set determinants n_states 4

qp run cis | tee ${EZFIO_FILE}.cis.B.out
sleep 1
qp run save_natorb | tee ${EZFIO_FILE}.natorb.B.out

qp set_mo_class -c "[1]" -a "[2-5]" -d "[6-18]"
qp run cis | tee ${EZFIO_FILE}.cis_2states_active_space.out
qp set determinants read_wf True

qp set_mo_class -c "[1]" -a "[2-5]" -v "[6-18]"

qp set casscf_cipsi small_active_space True 
qp set mol_properties calc_tr_dipole_moment true
qp run casscf | tee ${EZFIO_FILE}.casscf.B.out

# Export TREXIO using gammcor_plugin
qp set gammcor_plugin trexio_file \"B.hdf5\"
qp set gammcor_plugin cholesky_rdm False
qp set gammcor_plugin cholesky_tolerance 1.e-5
qp set gammcor_plugin export_tr_rdm True
qp run export_gammcor >> "export_B.out"
qp run gammcor_plugin >> "export_B.out"

