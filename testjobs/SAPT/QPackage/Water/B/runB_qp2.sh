#!/bin/bash

QP_ROOT=$HOME/qp2
source ${QP_ROOT}/quantum_package.rc 

geomfile=$1

## monomer A ###
# -a : bohr
qp create_ezfio -b aug-cc-pvdz "$geomfile" -a -o "cas_B"
qp set_file "cas_B"

qp set electrons elec_alpha_num 5
qp set electrons elec_beta_num  5

qp set davidson_keywords distributed_davidson False

qp set nuclei nucl_label  "['x','x','x', 'O', 'H', 'H']"
qp set nuclei nucl_charge "[0.0, 0.0, 0.0, 8.0, 1.0, 1.0]"


qp run scf | tee ${EZFIO_FILE}.scf.B.out
qp set_frozen_core
sleep 1

qp set determinants n_states 2

qp run cis | tee ${EZFIO_FILE}.cis.B.out
sleep 1
qp run save_natorb | tee ${EZFIO_FILE}.natorb.B.out

##CIS within the active space
#qp set_mo_class -c "[1]" -a "[2-9]" -d "[10-82]"
#qp run cis | tee ${EZFIO_FILE}.cis_2_states_active_spaceA.out
#qp set determinants read_wf True

qp set determinants read_wf True
qp set_mo_class -c "[1]" -a "[2-9]" -v "[10-82]"

qp set casscf_cipsi small_active_space True 
qp set mol_properties calc_tr_dipole_moment True
qp run casscf | tee ${EZFIO_FILE}.casscf.B.out

# Export TREXIO using gammcor_plugin
qp set gammcor_plugin trexio_file \"B.hdf5\"
qp set gammcor_plugin cholesky_rdm False
qp set gammcor_plugin cholesky_tolerance 1.e-5
qp set gammcor_plugin export_tr_rdm True
qp run export_gammcor >> "export_B.out"
qp run gammcor_plugin >> "export_B.out"

