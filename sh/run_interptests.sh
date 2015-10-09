#!/bin/bash

# Script for multiple interpolation tests
#------------------------------------------

source inc_vars.sh 

echo "            Basic Parameters     "
echo "---------------------------------"
echo " "
cat par/simul.par

#Interp method
awk '{ if ( NR == 3 ) { print "none";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 5 ) { print "HA";} else {print $0;} }'  par/simul2.par > par/simul.par

./imodel
./imodel
./imodel
./imodel
./imodel

#Interp method
awk '{ if ( NR == 3 ) { print "neartrv";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 5 ) { print "HA";} else {print $0;} }'  par/simul2.par > par/simul.par

./imodel
./imodel
./imodel
./imodel
./imodel

#Interp method
awk '{ if ( NR == 3 ) { print "lintrv";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 5 ) { print "HA";} else {print $0;} }'  par/simul2.par > par/simul.par

./imodel
./imodel
./imodel
./imodel
./imodel

#Interp method
awk '{ if ( NR == 3 ) { print "hermtrv";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 5 ) { print "HA";} else {print $0;} }'  par/simul2.par > par/simul.par

./imodel
./imodel
./imodel
./imodel
./imodel

#Interp method
awk '{ if ( NR == 3 ) { print "lsqhx";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 5 ) { print "HA";} else {print $0;} }'  par/simul2.par > par/simul.par

./imodel
./imodel
./imodel
./imodel
./imodel

#Interp method
awk '{ if ( NR == 3 ) { print "rbftr";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 5 ) { print "HA";} else {print $0;} }'  par/simul2.par > par/simul.par

./imodel
./imodel
./imodel
./imodel
./imodel

#Interp method
awk '{ if ( NR == 3 ) { print "rbfetr";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 5 ) { print "HA";} else {print $0;} }'  par/simul2.par > par/simul.par

./imodel
./imodel
./imodel
./imodel
./imodel

#Interp method
awk '{ if ( NR == 3 ) { print "rbfhx";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 5 ) { print "HA";} else {print $0;} }'  par/simul2.par > par/simul.par

./imodel
./imodel
./imodel
./imodel
./imodel

#Interp method
awk '{ if ( NR == 3 ) { print "natlap";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 5 ) { print "HA";} else {print $0;} }'  par/simul2.par > par/simul.par

./imodel
./imodel
./imodel
./imodel
./imodel

#Interp method
awk '{ if ( NR == 3 ) { print "natsib";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 5 ) { print "HA";} else {print $0;} }'  par/simul2.par > par/simul.par

./imodel
./imodel
./imodel
./imodel
./imodel

#Interp method
awk '{ if ( NR == 3 ) { print "natfar";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 5 ) { print "HA";} else {print $0;} }'  par/simul2.par > par/simul.par

./imodel
./imodel
./imodel
./imodel
./imodel


#Interp method
awk '{ if ( NR == 3 ) { print "lmshep";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 5 ) { print "HA";} else {print $0;} }'  par/simul2.par > par/simul.par

./imodel
./imodel
./imodel
./imodel
./imodel

#Interp method
awk '{ if ( NR == 3 ) { print "qdshep";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 5 ) { print "HA";} else {print $0;} }'  par/simul2.par > par/simul.par

./imodel
./imodel
./imodel
./imodel
./imodel


#Interp method
awk '{ if ( NR == 3 ) { print "none";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 5 ) { print "TA";} else {print $0;} }'  par/simul2.par > par/simul.par

./imodel
./imodel
./imodel
./imodel
./imodel

#Interp method
awk '{ if ( NR == 3 ) { print "neartrc";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 5 ) { print "TA";} else {print $0;} }'  par/simul2.par > par/simul.par

./imodel
./imodel
./imodel
./imodel
./imodel


#Interp method
awk '{ if ( NR == 3 ) { print "wach";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 5 ) { print "TA";} else {print $0;} }'  par/simul2.par > par/simul.par

./imodel
./imodel
./imodel
./imodel
./imodel


#Interp method
awk '{ if ( NR == 3 ) { print "none";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 5 ) { print "HC";} else {print $0;} }'  par/simul2.par > par/simul.par

./imodel
./imodel
./imodel
./imodel
./imodel


#Interp method
awk '{ if ( NR == 3 ) { print "wach";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 5 ) { print "HC";} else {print $0;} }'  par/simul2.par > par/simul.par

./imodel
./imodel
./imodel
./imodel
./imodel


#Interp method
awk '{ if ( NR == 3 ) { print "lmshep";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 5 ) { print "HC";} else {print $0;} }'  par/simul2.par > par/simul.par

./imodel
./imodel
./imodel
./imodel
./imodel

#Interp method
awk '{ if ( NR == 3 ) { print "p1nc";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 5 ) { print "HC";} else {print $0;} }'  par/simul2.par > par/simul.par

./imodel
./imodel
./imodel
./imodel
./imodel
