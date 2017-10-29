#! /bin/bash

MIN=1
MAX=9

for i in `seq $MIN $MAX`; do
    echo 
    echo "---------------------------------------------"
    echo "           Running grid level" $i
    echo "---------------------------------------------"
    echo
    ./gengrid_cube $i
    echo
    echo
done