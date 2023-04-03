#!/bin/bash

function dir_clean {

 rm AO* SIRI* DALT* occ* enuc.dat rdm2* 
 rm -r MONOMER_*

}

function occ_read_sym {
INAct=$(grep "Number of closed shell electrons" $1 | awk '{print $7}')
Act=$(grep "@    Active orbitals" $1 | awk '{print $4}')
echo '      '$INAct' ' $Act > occupations.dat

nsym=$(grep "Total number of sym" $1 | awk '{print $6}')
nsymm1=$(awk -v ii="$nsym" 'BEGIN {print ii - 1} ')

j=1
for i in $(seq 1 $nsym) ; do

    let j+=1
    start=$(grep -n "Symmetry $i" $1 | cut -d':' -f1)
    if [[ $i -lt $nsym ]] ; then
       stmp=$(grep  -n "Symmetry $j" $1 | cut -d':' -f1)
       stop=$(awk -v ii="$stmp" 'BEGIN {print ii - 2} ')
    else
       stmp=$(grep  -n "Printout of final " $1 | cut -d':' -f1)
       stop=$(awk -v ii="$stmp" 'BEGIN {print ii - 3} ')
    fi

    str=$(sed -n "$start p" $1)

    if [[ "$str" == *"Total"* ]] ; then
#      echo $str
#      echo $start,$stop

      awk -v x=$start -v y=$stop 'NR==x+2,NR==y' $1 >> tmp.dat

    fi

done

readarray occ < tmp.dat
echo '   '${occ[*]} >> occupations.dat
rm -f tmp.dat

inact=$(grep "Inactive orbitals" $1 | cut -d "|" -f2)
act=$(grep "Active orbitals" $1 | cut -d "|" -f2)
echo '   '$act >> occupations.dat
echo '   '$inact >> occupations.dat

}

# ****************************************
# FIRST, RUN DALTON CALCULATIONS
# SAPT GOES SECOND (prdmt.exe)

dir_clean

DALTON_EXEC="/home/michalhapka/dalton_mod/build/dalton"
PRDMFT_EXEC="/home/michalhapka/pr-dmft-open/prdmft.exe"

mkdir MONOMER_A
mkdir MONOMER_B

### CALCULATION STARTS HERE
cp template_A.mol  MONOMER_A/run.mol
cp template_A.dal  MONOMER_A/run.dal
cd MONOMER_A
$DALTON_EXEC -mb 500 -get "AOTWOINT AOONEINT AMFI_SYMINFO.TXT SIRIFC SIRIUS.RST DALTON.MOPUN rdm2.dat " run  > skas
#occ_read_sym run.out

cd ..
cp template_B.mol  MONOMER_B/run.mol
cp template_B.dal  MONOMER_B/run.dal

cd MONOMER_B
$DALTON_EXEC -mb 500 -get "AOTWOINT AOONEINT AMFI_SYMINFO.TXT SIRIFC SIRIUS.RST DALTON.MOPUN rdm2.dat " run  > skas
#occ_read_sym run.out

cd ..

##MONOMER A
mv MONOMER_A/run.AOONEINT AOONEINT_A
mv MONOMER_A/run.AOTWOINT AOTWOINT_A
mv MONOMER_A/run.SIRIFC SIRIFC_A
mv MONOMER_A/run.SIRIUS.RST SIRIUS_A.RST
mv MONOMER_A/run.DALTON.MOPUN DALTON_A.MOPUN
mv MONOMER_A/run.rdm2.dat rdm2_A.dat
mv MONOMER_A/occupations.dat occupations_A.dat 

rm MONOMER_A/run.tar* 
#
#MONOMER_B
mv MONOMER_B/run.AOONEINT AOONEINT_B
mv MONOMER_B/run.SIRIFC SIRIFC_B
mv MONOMER_B/run.SIRIUS.RST SIRIUS_B.RST
mv MONOMER_B/run.DALTON.MOPUN DALTON_B.MOPUN
mv MONOMER_B/run.rdm2.dat rdm2_B.dat 
mv MONOMER_B/occupations.dat occupations_B.dat 
mv MONOMER_B/run.AMFI_SYMINFO*  SYMINFO_B

rm MONOMER_B/run.AO* 
rm MONOMER_B/run.tar* 

# run SAPT
#$PRDMFT_EXEC > sapt_dalton.out 

# remove integrals, SIRI and RDMs
#rm AO* SIRI* SYM* rdm* DALTON*
