#!/bin/bash
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 024:00:00
#SBATCH --mem=40GB
#SBATCH --nodelist=cn06

GAMMCOR_EXEC="/home/michalhapka/pr-dmft/build_dSRS/gammcor"

export OMP_NUM_THREADS=8

WRKPATH=$(pwd)

mkdir -p /tmp/michal/$$
cd /tmp/michal/$$

cp $WRKPATH/A0.hdf5   .

# check theta=60
for r in 60 ; do

 echo $r
 mkdir -p $r
 mkdir -p $WRKPATH/$r

 xyz_file=$WRKPATH"/geom/r6.2_th"$r".xyz"

    cd $r
    cp $WRKPATH/run_qp2.sh .
    cp $WRKPATH/input.inp  .

    bash run_qp2.sh $xyz_file
    cp cas_A.casscf.A.out $WRKPATH/$r/
    cp cas_B.scf.B.out $WRKPATH/$r/
    cp export_A.out  $WRKPATH/$r/
    cp export_B.out  $WRKPATH/$r/

    cp ../A0.hdf5 .
    $GAMMCOR_EXEC > $WRKPATH/$r/$r"_GAMMCOR.out"
    cd ../
#    rm -r $r
done

rm -r /tmp/michal/$$

