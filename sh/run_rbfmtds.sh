#!/bin/bash

source inc_vars.sh 

# Script for multiple vector reconstruction tests
#------------------------------------------

echo "            Basic Parameters     "
echo "---------------------------------"
echo " "
#cat par/simul.par

#Recon method
awk '{ if ( NR == 5 ) { print "wach";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 7 ) { print "rbftr";} else {print $0;} }'  par/simul2.par > par/simul.par
sh/runrbftests.sh
#./runngrids9.sh

#Recon method
awk '{ if ( NR == 5 ) { print "wach";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 7 ) { print "rbfetr";} else {print $0;} }'  par/simul2.par > par/simul.par
sh/runrbftests.sh
#./runngrids9.sh


#Recon method
awk '{ if ( NR == 5 ) { print "lintrv";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 7 ) { print "rbfhx";} else {print $0;} }'  par/simul2.par > par/simul.par
awk '{ if ( NR == 19 ) { print "1 0";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ print $0; }'  par/simul2.par > par/simul.par
sh/runrbftests.sh
#./runngrids9.sh

#Recon method
awk '{ if ( NR == 5 ) { print "lintrv";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 7 ) { print "rbfhxpc";} else {print $0;} }'  par/simul2.par > par/simul.par
sh/runrbftests.sh
#./runngrids9.sh

#Recon method
awk '{ if ( NR == 5 ) { print "wach";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 7 ) { print "rbfetrp";} else {print $0;} }'  par/simul2.par > par/simul.par
sh/runrbftests.sh
#./runngrids9.sh

