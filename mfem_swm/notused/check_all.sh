#! /bin/bash

MIN=1
MAX=5

test -n "$1" && MIN=$1
test -n "$2" && MAX=$2

CMIN=$((MIN))
CMAX=$((MAX))

MIN=$((MIN+1))
MAX=$((MAX+1))

for i in `seq $CMIN $CMAX`; do
	NX=$((3*(2**(i-1))))
	CELLS=$((6*NX*NX))
	CELLS=`printf "%010d" $CELLS`

	IDENTIFIER="cube_$CELLS"

	echo "checking "$IDENTIFIER""

	diff restart_files/run000001_restart_00000020_"$IDENTIFIER".dat ../../trunk/original_version/restart_files/run000001_restart_00000020_"$IDENTIFIER".dat

	echo "checking "$IDENTIFIER" SFC"
	diff restart_files/run000001_restart_00000020_"$IDENTIFIER"_sfc.dat ../../trunk/original_version/restart_files/run000001_restart_00000020_"$IDENTIFIER".dat
done

for i in `seq $MIN $MAX`; do
	CELLS=$((5*2**(2*i-1)+2))
	CELLS=`printf "%010d" $CELLS`

	IDENTIFIER="hex_$CELLS"

	echo "checking "$IDENTIFIER""

	diff restart_files/run000001_restart_00000020_"$IDENTIFIER".dat ../../trunk/original_version/restart_files/run000001_restart_00000020_"$IDENTIFIER".dat

	echo "checking "$IDENTIFIER" SFC"

	diff restart_files/run000001_restart_00000020_"$IDENTIFIER"_sfc.dat ../../trunk/original_version/restart_files/run000001_restart_00000020_"$IDENTIFIER".dat
done
