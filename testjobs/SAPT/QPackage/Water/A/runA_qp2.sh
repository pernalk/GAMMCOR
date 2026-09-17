#!/bin/bash

QP_ROOT=$HOME/qp2
source ${QP_ROOT}/quantum_package.rc 

geomfile=$1

## monomer A ###
# -a : bohr
qp create_ezfio -b cc-pvdz "$geomfile" -a -o "cas_A"
qp set_file "cas_A"

qp set electrons elec_alpha_num 5
qp set electrons elec_beta_num  5

qp set davidson_keywords distributed_davidson False

qp set nuclei nucl_label  "['O', 'H', 'H', 'x', 'x', 'x']"
qp set nuclei nucl_charge "[8.0, 1.0, 1.0, 0.0, 0.0, 0.0]"

qp run scf | tee ${EZFIO_FILE}.scf.A.out
qp set_frozen_core
sleep 1

qp set determinants n_states 2

qp run cis | tee ${EZFIO_FILE}.cis.A.out
sleep 1
qp run save_natorb | tee ${EZFIO_FILE}.natorb.A.out

qp set determinants read_wf True
qp set_mo_class -c "[1]" -a "[2-9]" -v "[10-82]"

qp set casscf_cipsi small_active_space True 
qp set mol_properties calc_tr_dipole_moment True
qp run casscf | tee ${EZFIO_FILE}.casscf.A.out

# Export TREXIO using gammcor_plugin
qp set gammcor_plugin trexio_file \"A.hdf5\"
qp set gammcor_plugin cholesky_rdm False
qp set gammcor_plugin cholesky_tolerance 1.e-5
qp set gammcor_plugin export_tr_rdm True
qp run export_gammcor >> "export_A.out"
qp run gammcor_plugin >> "export_A.out"

