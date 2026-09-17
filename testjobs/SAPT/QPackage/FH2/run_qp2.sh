#!/bin/bash

QP_ROOT=$HOME/qp2
source ${QP_ROOT}/quantum_package.rc 

geomfile=$1

## monomer A ###
qp create_ezfio -b aug-cc-pvtz -m 2 "$geomfile" -a -o "cas_A"
qp set_file "cas_A"

qp set electrons elec_alpha_num 5
qp set electrons elec_beta_num  4

qp set nuclei nucl_label  "['F', 'x', 'x']"
qp set nuclei nucl_charge "[9.0, 0.0, 0.0]"

qp set determinants s2_eig true
qp set determinants expected_s2 0.75

qp run scf | tee ${EZFIO_FILE}.scf.A.out
qp set_frozen_core
sleep 1

#qp set determinants n_states 3
#
#qp run cis | tee ${EZFIO_FILE}.cis.A.out
#sleep 1
#qp run save_natorb | tee ${EZFIO_FILE}.natorb.A.out
#
##CIS within the active space
#qp set_mo_class -c "[1-2]" -a "[3-5]" -d "[6-92]"
#qp run cis | tee ${EZFIO_FILE}.cis_3_states_active_spaceA.out
#qp set determinants read_wf True
#
#qp set_mo_class -c "[1-2]" -a "[3-5]" -v "[6-92]"
#qp set casscf_cipsi small_active_space True 
#qp set mol_properties calc_tr_dipole_moment True
#qp run casscf | tee ${EZFIO_FILE}.casscf.A.out

# Export TREXIO using gammcor_plugin
qp set gammcor_plugin trexio_file \"A.hdf5\"
qp set gammcor_plugin cholesky_rdm False
qp set gammcor_plugin cholesky_tolerance 1.e-5
qp set gammcor_plugin export_tr_rdm True
qp run export_gammcor >> "export_A.out"
qp run gammcor_plugin >> "export_A.out"


# monomer B
qp create_ezfio -b aug-cc-pvtz "$geomfile" -m 2 -a -o "cas_B"
qp set_file "cas_B"

qp set electrons elec_alpha_num 1
qp set electrons elec_beta_num  1

qp set nuclei nucl_label  "['x', 'H', 'H']"
qp set nuclei nucl_charge "[0.0, 1.0, 1.0]"

qp set davidson_keywords distributed_davidson False

qp run scf | tee ${EZFIO_FILE}.scf.B.out
sleep 1

qp set_mo_class -a "[1]" -v "[2-92]"

# Export TREXIO using gammcor_plugin
qp set gammcor_plugin trexio_file \"B.hdf5\"
qp set gammcor_plugin cholesky_rdm False
qp set gammcor_plugin cholesky_tolerance 1.e-5
qp set gammcor_plugin export_tr_rdm False
qp run export_gammcor >> "export_B.out"
qp run gammcor_plugin >> "export_B.out"

