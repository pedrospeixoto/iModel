#!/bin/bash

source sh/inc_vars.sh 

# Script for multiple shallow water equation tests

# HR95 grid
awk '{ if ( NR == 17 ) { print "HR95JT_007.xyz";} else {print $0;} }'  par/mesh.par > par/mesh2.par
cp par/mesh2.par par/mesh.par

#new method
cp par/swmnew.par par/swm.par
awk '{ if ( NR == 27 ) { print "none 0";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel


#trisk
cp par/swmtrsk.par par/swm.par
awk '{ if ( NR == 27 ) { print "none 0";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel


# SCVT grid
awk '{ if ( NR == 17 ) { print "icos_pol_scvt_h1_7.xyz";} else {print $0;} }'  par/mesh.par > par/mesh2.par
cp par/mesh2.par par/mesh.par

#new method
cp par/swmnew.par par/swm.par
awk '{ if ( NR == 27 ) { print "none 0";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel


#trisk
cp par/swmtrsk.par par/swm.par
awk '{ if ( NR == 27 ) { print "none 0";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

