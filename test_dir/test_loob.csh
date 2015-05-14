#!/bin/tcsh

# assumes cmake has been run with -DCMAKE_BUILD_TYPE=DEBUG
../exe_DEBUG/loob \
    -S test.smt \
    -P TEST.points \
    -D bit_of_chembl_20.oeb.gz \
    --ascii-fps-file r.bits
