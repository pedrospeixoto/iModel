#!/bin/bash

source inc_vars.sh 

# Script for convergence analysis on transport model for different times steps
#------------------------------------------

echo "            Basic Parameters     "
echo "---------------------------------"
echo " "
cat par/swm.par

awk '{ if ( NR == 7 ) { print "0 1 0";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

awk '{ if ( NR == 7 ) { print "0 2 0";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

awk '{ if ( NR == 7 ) { print "0 4 0";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

awk '{ if ( NR == 7 ) { print "0 8 0";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

awk '{ if ( NR == 7 ) { print "0 16 0";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

awk '{ if ( NR == 7 ) { print "0 32 0";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

awk '{ if ( NR == 7 ) { print "0 64 0";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

awk '{ if ( NR == 7 ) { print "0 128 0";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

awk '{ if ( NR == 7 ) { print "0 256 0";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

awk '{ if ( NR == 7 ) { print "0 512 0";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

awk '{ if ( NR == 7 ) { print "0 1024 0";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

awk '{ if ( NR == 7 ) { print "0 2048 0";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel



