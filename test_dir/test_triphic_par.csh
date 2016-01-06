#!/bin/tcsh

# assumes cmake has been run with -DCMAKE_BUILD_TYPE=DEBUG - this dictates
# the path of the exe, nothing more than that.
# uses triphic's built-in SMARTS and points definitions
# uses openmpi

mpirun -np 4 ../exe_DEBUG/triphic \
    -Q 2Q61_lig_1.mol2 \
    -D bit_of_chembl_20.oeb.gz \
    -O triphic_test_par_out1.sdf \
    --protein 2Q61_prot.pdb

mpirun -np 4 ../exe_DEBUG/triphic \
    -Q 2Q61_lig_1.mol2 \
    -D chembl_20_first_100_ions_and_tauts.smi \
    -O triphic_test_par_out2.sdf \
    --protein 2Q61_prot.pdb \
    --use-omega --use-flipper
