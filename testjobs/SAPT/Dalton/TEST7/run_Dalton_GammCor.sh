#!/bin/bash

# ************ How to use this script ********
#  bash <path_to_dalton_executable> <path_to_gammcor_executable> 

DALTON_EXEC=$1
GAMMCOR_EXEC=$2

# ***********  Aim of the scipt is twofold ***
# 1) RUN DALTON  : CASSCF
# 2) RUN GAMMCOR : SAPT

export OMP_NUM_THREADS=1

mkdir MONOMER_A
mkdir MONOMER_B

### CALCULATION STARTS HERE
# ****************************************

cp template_A.mol  MONOMER_A/run.mol
cp template_A.dal  MONOMER_A/run.dal
cd MONOMER_A
$DALTON_EXEC -noarch -mb 100 -get "AOTWOINT AOONEINT SIRIFC SIRIUS.RST rdm2.dat " run  > skas

cd ..
cp template_B.mol  MONOMER_B/run.mol
cp template_B.dal  MONOMER_B/run.dal

cd MONOMER_B
$DALTON_EXEC -noarch -mb 100 -get "AOONEINT AMFI_SYMINFO.TXT SIRIFC SIRIUS.RST rdm2.dat " run  > skas

cd ..

#MONOMER A
mv MONOMER_A/run.AOONEINT    AOONEINT_A
mv MONOMER_A/run.AOTWOINT    AOTWOINT_A
mv MONOMER_A/run.SIRIFC      SIRIFC_A
mv MONOMER_A/run.SIRIUS.RST  SIRIUS_A.RST
mv MONOMER_A/run.rdm2.dat    rdm2_A.dat

#MONOMER_B
mv MONOMER_B/run.AOONEINT       AOONEINT_B
mv MONOMER_B/run.SIRIFC         SIRIFC_B
mv MONOMER_B/run.SIRIUS.RST     SIRIUS_B.RST
mv MONOMER_B/run.rdm2.dat       rdm2_B.dat 
mv MONOMER_B/run.AMFI_SYMINFO*  SYMINFO_B

# run GammCor
$GAMMCOR_EXEC > gammcor.out

