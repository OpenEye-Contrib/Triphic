#!/bin/tcsh

# assumes cmake has been run with -DCMAKE_BUILD_TYPE=DEBUG - this dictates
# the path of the exe, nothing more than that.
# uses triphic's built-in SMARTS and points definitions
# uses openmpi

mpirun -np 4 ../exe_DEBUG/triphic \
    -Q 2Q61_lig_1.mol2 \
    -D bit_of_chembl_20.oeb.gz \
    -O r.sdf \
    --protein 2Q61_prot.pdb

