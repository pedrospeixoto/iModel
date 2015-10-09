#!/bin/bash

source inc_vars.sh 

# Script for multiple vector reconstruction tests
#------------------------------------------

echo "            Basic Parameters     "
echo "---------------------------------"
echo " "
#cat par/simul.par

#RosbyWave
awk '{ if ( NR == 21 ) { print "8";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 21 ) { print "8";} else {print $0;} }'  par/simul2.par > par/simul.par

#Recon method
awk '{ if ( NR == 7 ) { print "rbfhx";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 19 ) { print "1 0";} else {print $0;} }'  par/simul2.par > par/simul.par
#./imodel
#./runngrids9.sh

#Recon method
awk '{ if ( NR == 7 ) { print "rbfhxpc";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 19 ) { print "1 0";} else {print $0;} }'  par/simul2.par > par/simul.par
#./imodel
./runngrids9.sh

#Recon method
awk '{ if ( NR == 7 ) { print "rbfhx";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 19 ) { print "0.001 1";} else {print $0;} }'  par/simul2.par > par/simul.par
#./imodel
./runngrids9.sh

#Recon method
awk '{ if ( NR == 7 ) { print "rbfhxpc";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 19 ) { print "0.001 1";} else {print $0;} }'  par/simul2.par > par/simul.par
#./imodel
./runngrids9.sh

#Recon method
awk '{ if ( NR == 7 ) { print "rbfhx";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 19 ) { print "0.005 1";} else {print $0;} }'  par/simul2.par > par/simul.par
#./imodel
./runngrids9.sh

#Recon method
awk '{ if ( NR == 7 ) { print "rbfhxpc";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 19 ) { print "0.005 1";} else {print $0;} }'  par/simul2.par > par/simul.par
#./imodel
./runngrids9.sh

#Recon method
awk '{ if ( NR == 7 ) { print "rbfhx";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 19 ) { print "0.011 1";} else {print $0;} }'  par/simul2.par > par/simul.par
#./imodel
./runngrids9.sh

#Recon method
awk '{ if ( NR == 7 ) { print "rbfhxpc";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 19 ) { print "0.011 1";} else {print $0;} }'  par/simul2.par > par/simul.par
#./imodel
./runngrids9.sh

#Rotation
awk '{ if ( NR == 21 ) { print "5";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 21 ) { print "5";} else {print $0;} }'  par/simul2.par > par/simul.par

#Recon method
awk '{ if ( NR == 7 ) { print "rbfhx";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 19 ) { print "1 0";} else {print $0;} }'  par/simul2.par > par/simul.par
#./imodel
./runngrids9.sh

#Recon method
awk '{ if ( NR == 7 ) { print "rbfhxpc";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 19 ) { print "1 0";} else {print $0;} }'  par/simul2.par > par/simul.par
#./imodel
./runngrids9.sh

#Recon method
awk '{ if ( NR == 7 ) { print "rbfhx";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 19 ) { print "0.001 1";} else {print $0;} }'  par/simul2.par > par/simul.par
#./imodel
./runngrids9.sh

#Recon method
awk '{ if ( NR == 7 ) { print "rbfhxpc";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 19 ) { print "0.001 1";} else {print $0;} }'  par/simul2.par > par/simul.par
#./imodel
./runngrids9.sh

#Recon method
awk '{ if ( NR == 7 ) { print "rbfhx";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 19 ) { print "0.005 1";} else {print $0;} }'  par/simul2.par > par/simul.par
#./imodel
./runngrids9.sh

#Recon method
awk '{ if ( NR == 7 ) { print "rbfhxpc";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 19 ) { print "0.005 1";} else {print $0;} }'  par/simul2.par > par/simul.par
#./imodel
./runngrids9.sh

#Recon method
awk '{ if ( NR == 7 ) { print "rbfhx";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 19 ) { print "0.011 1";} else {print $0;} }'  par/simul2.par > par/simul.par
#./imodel
./runngrids9.sh

#Recon method
awk '{ if ( NR == 7 ) { print "rbfhxpc";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 19 ) { print "0.011 1";} else {print $0;} }'  par/simul2.par > par/simul.par
#./imodel
./runngrids9.sh



exit 0

#Recon method
awk '{ if ( NR == 3 ) { print "nonehx";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 5 ) { print "HC";} else {print $0;} }'  par/simul2.par > par/simul.par

./runngrids.sh

#Recon method
awk '{ if ( NR == 3 ) { print "nonetr";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 5 ) { print "HC";} else {print $0;} }'  par/simul2.par > par/simul.par

./runngrids.sh

#Recon method
awk '{ if ( NR == 3 ) { print "+1perhx+2lintrv+";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 5 ) { print "HC";} else {print $0;} }'  par/simul2.par > par/simul.par

./runngrids.sh

#Recon method
awk '{ if ( NR == 3 ) { print "+1lsqhxe+2lintrv+";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 5 ) { print "HC";} else {print $0;} }'  par/simul2.par > par/simul.par

./runngrids.sh

#Recon method
awk '{ if ( NR == 3 ) { print "+1perhx+2lintrv+3lsqhxe+";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 5 ) { print "HC";} else {print $0;} }'  par/simul2.par > par/simul.par

./runngrids.sh

#Recon method
awk '{ if ( NR == 3 ) { print "+1wht+2wach+";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 5 ) { print "HC";} else {print $0;} }'  par/simul2.par > par/simul.par

./runngrids.sh

#Recon method
awk '{ if ( NR == 3 ) { print "+1lsqtrc+2wach+";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 5 ) { print "HC";} else {print $0;} }'  par/simul2.par > par/simul.par

./runngrids.sh

#Recon method
awk '{ if ( NR == 3 ) { print "+1rbftr+2wach+";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 5 ) { print "HC";} else {print $0;} }'  par/simul2.par > par/simul.par

./runngrids.sh

#Recon method
awk '{ if ( NR == 3 ) { print "+1rbfetr+2wach+";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 5 ) { print "HC";} else {print $0;} }'  par/simul2.par > par/simul.par

./runngrids.sh

#Recon method
awk '{ if ( NR == 3 ) { print "+1pertr+2wach+";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 5 ) { print "HC";} else {print $0;} }'  par/simul2.par > par/simul.par

./runngrids.sh

#Recon method
awk '{ if ( NR == 3 ) { print "+1rbfhx+2lintrv+";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 5 ) { print "HC";} else {print $0;} }'  par/simul2.par > par/simul.par

./runngrids.sh

#Recon method
awk '{ if ( NR == 3 ) { print "+1kls+2lintrv+";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 5 ) { print "HC";} else {print $0;} }'  par/simul2.par > par/simul.par

./runngrids.sh
