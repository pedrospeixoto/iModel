#!/bin/bash

source inc_vars.sh 

# Script for multiple shallow water equation tests

#awk '{ if ( NR == 17 ) { print "HR95JT_007.xyz";} else {print $0;} }'  par/mesh.par > par/mesh2.par
#cp par/mesh2.par par/mesh.par

#Test for different holinsgsworth parameters
cp par/swm.par par/swm_original.par
awk '{ if ( NR == 11 ) { print "gass 0.0";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

cp par/swm.par par/swm_original.par
awk '{ if ( NR == 11 ) { print "gass 0.1";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

cp par/swm.par par/swm_original.par
awk '{ if ( NR == 11 ) { print "gass 0.2";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

cp par/swm.par par/swm_original.par
awk '{ if ( NR == 11 ) { print "gass 0.3";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

cp par/swm.par par/swm_original.par
awk '{ if ( NR == 11 ) { print "gass 0.4";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

cp par/swm.par par/swm_original.par
awk '{ if ( NR == 11 ) { print "gass 0.5";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

cp par/swm.par par/swm_original.par
awk '{ if ( NR == 11 ) { print "gass 0.6";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel


cp par/swm.par par/swm_original.par
awk '{ if ( NR == 11 ) { print "gass 0.625";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

cp par/swm.par par/swm_original.par
awk '{ if ( NR == 11 ) { print "gass 0.65";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

cp par/swm.par par/swm_original.par
awk '{ if ( NR == 11 ) { print "gass 0.7";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

cp par/swm.par par/swm_original.par
awk '{ if ( NR == 11 ) { print "gass 0.75";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

cp par/swm.par par/swm_original.par
awk '{ if ( NR == 11 ) { print "gass 0.8";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

cp par/swm.par par/swm_original.par
awk '{ if ( NR == 11 ) { print "gass 0.85";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

cp par/swm.par par/swm_original.par
awk '{ if ( NR == 11 ) { print "gass 0.9";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

cp par/swm.par par/swm_original.par
awk '{ if ( NR == 11 ) { print "gass 0.95";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel

cp par/swm.par par/swm_original.par
awk '{ if ( NR == 11 ) { print "gass 1.0";} else {print $0;} }'  par/swm.par > par/swm2.par
cp par/swm2.par par/swm.par
./imodel
