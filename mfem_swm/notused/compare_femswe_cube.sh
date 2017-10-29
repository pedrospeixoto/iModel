#! /bin/bash

test -n "$1" && LVL=$1

NX=$((3*(2**(LVL-1))))
CELLS=$((6*NX*NX))
CELLS=`printf "%010d" $CELLS`

IDENTIFIER="cube_$CELLS"

make clean
make femswe

echo "simulation of "$IDENTIFIER""

perf stat -e cycles,instructions,cache-references,cache-misses,branch-instructions,branch-misses ./femswe swenml/swenml_"$IDENTIFIER".in
mv run000001_restart_00000020.dat run000001_restart_00000020_"$IDENTIFIER".dat

echo "simulation of "$IDENTIFIER" SFC"

perf stat -e cycles,instructions,cache-references,cache-misses,branch-instructions,branch-misses ./femswe swenml/swenml_"$IDENTIFIER"_sfc.in
mv run000001_restart_00000020.dat run000001_restart_00000020_"$IDENTIFIER"_sfc.dat

echo "Difference?"
diff -q run000001_restart_00000020_"$IDENTIFIER"_sfc.dat run000001_restart_00000020_"$IDENTIFIER".dat

