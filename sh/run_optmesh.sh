#!/bin/bash

source inc_vars.sh 

# Script for multiple grid tests
#------------------------------------------

echo "            Basic Parameters     "
echo "---------------------------------"
echo " "
cat par/mesh.par

echo "-----------------NOPT - GRID-------------------------------------------"
#Grid
awk '{ if ( NR == 9 ) { print "nopt";} else {print $0;} }'  par/mesh.par > par/mesh2.par
#Test case mesh quality
awk '{ if ( NR == 15 ) { print 3;} else {print $0;} }'  par/mesh2.par > par/mesh.par

cat par/mesh.par
echo

# runtests
./runntimes.sh

#Test case divs
awk '{ if ( NR == 15 ) { print 4;} else {print $0;} }'  par/mesh2.par > par/mesh.par
./runntimes.sh

echo "-----------------SPRG - GRID-------------------------------------------"
#Grid
awk '{ if ( NR == 9 ) { print "sprg";} else {print $0;} }'  par/mesh.par > par/mesh2.par
#Hierarchy flag to symmetry
awk '{ if ( NR == 19 ) { print "3";} else {print $0;} }'  par/mesh2.par > par/mesh3.par
#Test case mesh quality
awk '{ if ( NR == 15 ) { print 3;} else {print $0;} }'  par/mesh3.par > par/mesh.par

cat par/mesh.par

# runtests
./runntimes.sh

#Test case divs
awk '{ if ( NR == 15 ) { print 4;} else {print $0;} }'  par/mesh3.par > par/mesh.par
./runntimes.sh

#Return Hierarchy flag to normal
awk '{ if ( NR == 19 ) { print "1";} else {print $0;} }'  par/mesh3.par > par/mesh.par


echo "-----------------SCVT - GRID-------------------------------------------"
#Grid
awk '{ if ( NR == 9 ) { print "scvt";} else {print $0;} }'  par/mesh.par > par/mesh2.par
#Test case mesh quality
awk '{ if ( NR == 15 ) { print 3;} else {print $0;} }'  par/mesh2.par > par/mesh.par

cat par/mesh.par

# runtests
./runntimes.sh

#Test case div
awk '{ if ( NR == 15 ) { print 4;} else {print $0;} }'  par/mesh2.par > par/mesh.par
./runntimes.sh


echo "-----------------HR95 grid - GRID-------------------------------------------"
#Grid
awk '{ if ( NR == 5 ) { print "read";} else {print $0;} }'  par/mesh.par > par/mesh2.par
awk '{ if ( NR == 9 ) { print "nopt";} else {print $0;} }'  par/mesh2.par > par/mesh3.par
#Test case mesh quality
awk '{ if ( NR == 15 ) { print 3;} else {print $0;} }'  par/mesh3.par > par/mesh.par

cat par/mesh.par

# runtests
./runntimesHR95.sh

#Test case divs
awk '{ if ( NR == 15 ) { print 4;} else {print $0;} }'  par/mesh3.par > par/mesh.par
./runntimes.sh
