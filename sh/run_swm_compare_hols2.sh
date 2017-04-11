#!/bin/bash

source inc_vars.sh 

# Script for multiple shallow water equation tests

#awk '{ if ( NR == 17 ) { print "HR95JT_007.xyz";} else {print $0;} }'  par/mesh.par > par/mesh2.par
#cp par/mesh2.par par/mesh.par

#Test for different holinsgsworth parameters
cp par/swm.par par/swm_original2.par
awk '{ if ( NR == 23 ) { print "1";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

cp par/swm.par par/swm2.par
awk '{ if ( NR == 23 ) { print "2";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

cp par/swm.par par/swm2.par
awk '{ if ( NR == 23 ) { print "3";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

cp par/swm.par par/swm_original.par
awk '{ if ( NR == 23 ) { print "5";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

cp par/swm.par par/swm2.par
awk '{ if ( NR == 23 ) { print "10";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

cp par/swm.par par/swm2.par
awk '{ if ( NR == 23 ) { print "20";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

cp par/swm.par par/swm2.par
awk '{ if ( NR == 23 ) { print "30";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

cp par/swm.par par/swm2.par
awk '{ if ( NR == 23 ) { print "40";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

cp par/swm.par par/swm2.par
awk '{ if ( NR == 23 ) { print "50";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

cp par/swm.par par/swm2.par
awk '{ if ( NR == 23 ) { print "50";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

cp par/swm.par par/swm2.par
awk '{ if ( NR == 23 ) { print "100";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

cp par/swm.par par/swm2.par
awk '{ if ( NR == 23 ) { print "1000";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

cp par/swm.par par/swm2.par
awk '{ if ( NR == 23 ) { print "10000";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

