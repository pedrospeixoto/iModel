#!/bin/bash

source inc_vars.sh 

# Script for convergence analysis on Heikes SCVT grids
#------------------------------------------

echo "            Basic Parameters     "
echo "---------------------------------"
echo " "
cat par/mesh.par

awk '{ if ( NR == 5 ) { print "read";} else {print $0;} }'  par/mesh.par > par/mesh2.par
awk '{ if ( NR == 9 ) { print "nopt";} else {print $0;} }'  par/mesh2.par > par/mesh.par
#cp par/mesh2.par par/mesh.par

awk '{ if ( NR == 17 ) { print "point_centroidal_000000000642.xyz";} else {print $0;} }'  par/mesh.par > par/mesh2.par
cp par/mesh2.par par/mesh.par
./imodel

awk '{ if ( NR == 17 ) { print "point_centroidal_000000002562.xyz";} else {print $0;} }'  par/mesh.par > par/mesh2.par
cp par/mesh2.par par/mesh.par
./imodel

awk '{ if ( NR == 17 ) { print "point_centroidal_000000010242.xyz";} else {print $0;} }'  par/mesh.par > par/mesh2.par
cp par/mesh2.par par/mesh.par
./imodel

awk '{ if ( NR == 17 ) { print "point_centroidal_000000040962.xyz";} else {print $0;} }'  par/mesh.par > par/mesh2.par
cp par/mesh2.par par/mesh.par
./imodel

awk '{ if ( NR == 17 ) { print "point_centroidal_000000163842.xyz";} else {print $0;} }'  par/mesh.par > par/mesh2.par
cp par/mesh2.par par/mesh.par
./imodel

awk '{ if ( NR == 17 ) { print "point_centroidal_000000655362.xyz";} else {print $0;} }'  par/mesh.par > par/mesh2.par
cp par/mesh2.par par/mesh.par
./imodel

#awk '{ if ( NR == 17 ) { print "point_centroidal_000000000642.xyz";} else {print $0;} }'  par/mesh.par > par/mesh2.par
#cp par/mesh2.par par/mesh.par
#./imodel
