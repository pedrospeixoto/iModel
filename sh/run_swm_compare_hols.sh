#!/bin/bash

source sh/inc_vars.sh 

# Script for multiple shallow water equation tests

#awk '{ if ( NR == 17 ) { print "HR95JT_007.xyz";} else {print $0;} }'  par/mesh.par > par/mesh2.par
#cp par/mesh2.par par/mesh.par

#Test for different holinsgsworth parameters
cp par/swm.par par/swm_original.par
awk '{ if ( NR == 27 ) { print "0.001";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

cp par/swm.par par/swm2.par
awk '{ if ( NR == 27 ) { print "0.01";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

cp par/swm.par par/swm2.par
awk '{ if ( NR == 27 ) { print "0.1";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

cp par/swm.par par/swm_original.par
awk '{ if ( NR == 27 ) { print "0.2";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

cp par/swm.par par/swm2.par
awk '{ if ( NR == 27 ) { print "0.3";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

cp par/swm.par par/swm2.par
awk '{ if ( NR == 27 ) { print "0.4";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

cp par/swm.par par/swm2.par
awk '{ if ( NR == 27 ) { print "0.5";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

cp par/swm.par par/swm2.par
awk '{ if ( NR == 27 ) { print "1.0";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

cp par/swm.par par/swm2.par
awk '{ if ( NR == 27 ) { print "2.0";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

cp par/swm.par par/swm2.par
awk '{ if ( NR == 27 ) { print "3.0";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

cp par/swm.par par/swm2.par
awk '{ if ( NR == 27 ) { print "4.0";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

cp par/swm.par par/swm2.par
awk '{ if ( NR == 27 ) { print "5.0";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

cp par/swm.par par/swm2.par
awk '{ if ( NR == 27 ) { print "10.0";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

