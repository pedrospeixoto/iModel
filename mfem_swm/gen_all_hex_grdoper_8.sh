#! /bin/bash

MIN=1
MAX=9


for i in `seq $MIN $MAX`; do
    NX=$((3*(2**(i-1))))
    CELLS=$((6*NX*NX))
    CELLS=`printf "%010d" $CELLS`

    IDENTIFIER="cube_$CELLS"
    GRDFILE="grd/gridmap_cube_$CELLS.dat"

    echo 
    echo "---------------------------------------------"
    echo "           Running grid level" $i
    echo "---------------------------------------------"
    echo
    echo $GRDFILE
    ./buildop_fem cube 0 $GRDFILE
    echo
    echo
done