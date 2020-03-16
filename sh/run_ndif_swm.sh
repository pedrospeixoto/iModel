#!/bin/bash

source inc_vars.sh 

# Script for convergence analysis on transport model for different times steps
#------------------------------------------

echo "            Basic Parameters     "
echo "---------------------------------"
echo " "
cat par/swm.par

for i in 1 2 3 4 5 6; do
    dif=$(((10) ** ($i)))
    #dif=$i
    echo $dif
    awk -v dif="$dif" '{ if ( NR == 25 ) { print dif;} else {print $0;} }'  par/swm.par > par/swm2.par
    cp par/swm2.par par/swm.par
    ./imodel
done

