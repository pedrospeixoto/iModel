#!/bin/bash

# Script for multiple deformational tests
#------------------------------------------

source inc_vars.sh 

echo "            Basic Parameters     "
echo "---------------------------------"
echo " "
cat par/mesh.par

echo "-----------------NOPT - GRID-------------------------------------------"
#Grid
awk '{ if ( NR == 9 ) { print "nopt";} else {print $0;} }'  par/mesh.par > par/mesh2.par
#Test case
awk '{ if ( NR == 15 ) { print 3;} else {print $0;} }'  par/mesh2.par > par/mesh.par

cat par/mesh.par
echo

#./imodel
# runtests
./runntimes.sh

#Test case
awk '{ if ( NR == 15 ) { print 4;} else {print $0;} }'  par/mesh2.par > par/mesh.par
./runntimes.sh

echo "-----------------SCVT - GRID-------------------------------------------"
#Grid
awk '{ if ( NR == 9 ) { print "scvt";} else {print $0;} }'  par/mesh.par > par/mesh2.par
#Test case
awk '{ if ( NR == 15 ) { print 3;} else {print $0;} }'  par/mesh2.par > par/mesh.par

cat par/mesh.par

#./driver
# runtests
./runntimes.sh

#Test case
awk '{ if ( NR == 15 ) { print 4;} else {print $0;} }'  par/mesh2.par > par/mesh.par
./runntimes.sh

echo "-----------------SPRG - GRID-------------------------------------------"
#Grid
awk '{ if ( NR == 9 ) { print "sprg";} else {print $0;} }'  par/mesh.par > par/mesh2.par
#Test case
awk '{ if ( NR == 15 ) { print 3;} else {print $0;} }'  par/mesh2.par > par/mesh.par

cat par/mesh.par

#./driver
# runtests
./runntimes.sh

#Test case
awk '{ if ( NR == 15 ) { print 4;} else {print $0;} }'  par/mesh2.par > par/mesh.par
./runntimes.sh

echo "-----------------SALT - GRID-------------------------------------------"
#Grid
awk '{ if ( NR == 9 ) { print "salt";} else {print $0;} }'  par/mesh.par > par/mesh2.par
#Test case
awk '{ if ( NR == 15 ) { print 3;} else {print $0;} }'  par/mesh2.par > par/mesh.par

cat par/mesh.par

#./driver
# runtests
./runntimes.sh

#Test case
awk '{ if ( NR == 15 ) { print 4;} else {print $0;} }'  par/mesh2.par > par/mesh.par
./runntimes.sh
