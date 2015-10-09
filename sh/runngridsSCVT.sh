#!/bin/bash

source inc_vars.sh 

# Script for convergence analysis on HR95 grids
#------------------------------------------

echo "            Basic Parameters     "
echo "---------------------------------"
echo " "
cat par/mesh.par

awk '{ if ( NR == 5 ) { print "icos";} else {print $0;} }'  par/mesh.par > par/mesh2.par
awk '{ if ( NR == 9 ) { print "scvt";} else {print $0;} }'  par/mesh2.par > par/mesh.par
#cp par/mesh2.par par/mesh.par

awk '{ if ( NR == 17 ) { print "icos_pol_scvt_h0_0.gmt";} else {print $0;} }'  par/mesh.par > par/mesh2.par
cp par/mesh2.par par/mesh.par
./imodel

awk '{ if ( NR == 17 ) { print "icos_pol_scvt_h1_1.gmt";} else {print $0;} }'  par/mesh.par > par/mesh2.par
cp par/mesh2.par par/mesh.par
./imodel

awk '{ if ( NR == 17 ) { print "icos_pol_scvt_h1_2.gmt";} else {print $0;} }'  par/mesh.par > par/mesh2.par
cp par/mesh2.par par/mesh.par
./imodel

awk '{ if ( NR == 17 ) { print "icos_pol_scvt_h1_3.gmt";} else {print $0;} }'  par/mesh.par > par/mesh2.par
cp par/mesh2.par par/mesh.par
./imodel

awk '{ if ( NR == 17 ) { print "icos_pol_scvt_h1_4.gmt";} else {print $0;} }'  par/mesh.par > par/mesh2.par
cp par/mesh2.par par/mesh.par
./imodel

awk '{ if ( NR == 17 ) { print "icos_pol_scvt_h1_5.gmt";} else {print $0;} }'  par/mesh.par > par/mesh2.par
cp par/mesh2.par par/mesh.par
./imodel

awk '{ if ( NR == 17 ) { print "icos_pol_scvt_h1_6.gmt";} else {print $0;} }'  par/mesh.par > par/mesh2.par
cp par/mesh2.par par/mesh.par
./imodel

awk '{ if ( NR == 17 ) { print "icos_pol_scvt_h1_7.gmt";} else {print $0;} }'  par/mesh.par > par/mesh2.par
cp par/mesh2.par par/mesh.par
./imodel

#awk '{ if ( NR == 17 ) { print "icos_pol_scvt_h1_8.gmt";} else {print $0;} }'  par/mesh.par > par/mesh2.par
#cp par/mesh2.par par/mesh.par
#./imodel

#awk '{ if ( NR == 17 ) { print "icos_pol_scvt_h1_9.gmt";} else {print $0;} }'  par/mesh.par > par/mesh2.par
#cp par/mesh2.par par/mesh.par
#./imodel
