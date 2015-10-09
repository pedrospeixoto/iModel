#!/bin/bash

source inc_vars.sh 

# Script for convergence analysis on transport model for different times steps
#------------------------------------------

echo "            Basic Parameters     "
echo "---------------------------------"
echo " "
cat par/trans.par

awk '{ if ( NR == 9 ) { print "8 0";} else {print $0;} }'  par/trans.par > par/trans2.par
cp par/trans2.par par/trans.par
./runngrids.sh

awk '{ if ( NR == 9 ) { print "16 0";} else {print $0;} }'  par/trans.par > par/trans2.par
cp par/trans2.par par/trans.par
./runngrids.sh

awk '{ if ( NR == 9 ) { print "32 0";} else {print $0;} }'  par/trans.par > par/trans2.par
cp par/trans2.par par/trans.par
./runngrids.sh

awk '{ if ( NR == 9 ) { print "64 0";} else {print $0;} }'  par/trans.par > par/trans2.par
cp par/trans2.par par/trans.par
./runngrids.sh

awk '{ if ( NR == 9 ) { print "128 0";} else {print $0;} }'  par/trans.par > par/trans2.par
cp par/trans2.par par/trans.par
./runngrids.sh

awk '{ if ( NR == 9 ) { print "256 0";} else {print $0;} }'  par/trans.par > par/trans2.par
cp par/trans2.par par/trans.par
./runngrids.sh

awk '{ if ( NR == 9 ) { print "512 0";} else {print $0;} }'  par/trans.par > par/trans2.par
cp par/trans2.par par/trans.par
./runngrids.sh

awk '{ if ( NR == 9 ) { print "1024 0";} else {print $0;} }'  par/trans.par > par/trans2.par
cp par/trans2.par par/trans.par
./runngrids.sh

awk '{ if ( NR == 9 ) { print "2048 0";} else {print $0;} }'  par/trans.par > par/trans2.par
cp par/trans2.par par/trans.par
./runngrids.sh


