#!/bin/bash

source inc_vars.sh 

# Script for multiple vector reconstruction tests
#------------------------------------------

echo "            Basic Parameters     "
echo "---------------------------------"
echo " "

#Recon method
awk '{ if ( NR == 19 ) { print "perhx";} else {print $0;} }'  par/trans.par > par/trans2.par
awk '{ if ( NR == 23 ) { print "0";} else {print $0;} }'  par/trans2.par > par/trans.par
./imodel
#./runngrids9.sh

#Recon method
awk '{ if ( NR == 19 ) { print "perhx";} else {print $0;} }'  par/trans.par > par/trans2.par
awk '{ if ( NR == 23 ) { print "lsqhxe";} else {print $0;} }'  par/trans2.par > par/trans.par
./imodel

#Recon method
awk '{ if ( NR == 19 ) { print "rbfhx";} else {print $0;} }'  par/trans.par > par/trans2.par
awk '{ if ( NR == 23 ) { print "0";} else {print $0;} }'  par/trans2.par > par/trans.par
./imodel

#Recon method
awk '{ if ( NR == 19 ) { print "lsqhxe";} else {print $0;} }'  par/trans.par > par/trans2.par
awk '{ if ( NR == 23 ) { print "0";} else {print $0;} }'  par/trans2.par > par/trans.par
./imodel

#Recon method
awk '{ if ( NR == 19 ) { print "lsqtrc";} else {print $0;} }'  par/trans.par > par/trans2.par
awk '{ if ( NR == 23 ) { print "0";} else {print $0;} }'  par/trans2.par > par/trans.par
./imodel
