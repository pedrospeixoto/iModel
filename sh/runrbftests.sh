#!/bin/bash

source inc_vars.sh 

# Script for multiple vector reconstruction tests
#------------------------------------------

rbf_h="0.02 0.04 0.06 0.08 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2 3 4 5 6 8"
  
for h in $rbf_h;
do
#Recon method
    #echo $h
    awk '{ if ( NR == 7 ) { print "rbfetr";} else {print $0;} }'  par/simul.par > par/simul2.par
    awk -v par="$h 0" '{ if ( NR == 19 ) { print par;} else {print $0;} }'  par/simul2.par > par/simul.par
    #./runngrids9.sh
    ./imodel
done




