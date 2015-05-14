#!/bin/tcsh

# assumes cmake has been run with -DCMAKE_BUILD_TYPE=DEBUG
../exe_DEBUG/plurality \
    -Q second_test.pphore \
    -D bit_of_chembl_20.oeb.gz \
    -O r.sdf \
    -S test.smt \
    -P TEST.points
