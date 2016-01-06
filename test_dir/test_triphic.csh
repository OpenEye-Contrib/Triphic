#!/bin/tcsh

# assumes cmake has been run with -DCMAKE_BUILD_TYPE=DEBUG
# uses triphic's built-in SMARTS and points definitions

../exe_DEBUG/triphic \
    -Q 2Q61_lig_1.mol2 \
    -D bit_of_chembl_20.oeb.gz \
    -O triphic_test_out.sdf \
    --protein 2Q61_prot.pdb
