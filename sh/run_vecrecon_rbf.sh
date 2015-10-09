#!/bin/bash

source inc_vars.sh 

# Script for multiple vector reconstruction tests RBF methods
#cat par/simul.par

#RosbyWave
awk '{ if ( NR == 21 ) { print "8";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 21 ) { print "8";} else {print $0;} }'  par/simul2.par > par/simul.par


#RBF Recon methods - Variable shape par 0.0755/h -> results in eps=2 for glevel 4
############################

awk '{ if ( NR == 5 ) { print "lintrv";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 7 ) { print "rbfhx";} else {print $0;} }'  par/simul2.par > par/simul.par
awk '{ if ( NR == 19 ) { print "0.0755 1";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ print $0; }'  par/simul2.par > par/simul.par

sh/runngrids9.sh

#Recon method
awk '{ if ( NR == 5 ) { print "lintrv";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 7 ) { print "rbfhxpc";} else {print $0;} }'  par/simul2.par > par/simul.par

sh/runngrids9.sh

#Recon method
awk '{ if ( NR == 5 ) { print "wach";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 7 ) { print "rbfetr";} else {print $0;} }'  par/simul2.par > par/simul.par

sh/runngrids9.sh

#Recon method
awk '{ if ( NR == 5 ) { print "wach";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 7 ) { print "rbfetrp";} else {print $0;} }'  par/simul2.par > par/simul.par

sh/runngrids9.sh

exit 0

#Recon method
awk '{ if ( NR == 5 ) { print "wach";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 7 ) { print "rbftr";} else {print $0;} }'  par/simul2.par > par/simul.par

sh/runngrids9.sh


#RBF Recon methods - Fixed shape par
############################

awk '{ if ( NR == 5 ) { print "lintrv";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 7 ) { print "rbfhx";} else {print $0;} }'  par/simul2.par > par/simul.par
awk '{ if ( NR == 19 ) { print "2 0";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ print $0; }'  par/simul2.par > par/simul.par

sh/runngrids9.sh

#Recon method
awk '{ if ( NR == 5 ) { print "lintrv";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 7 ) { print "rbfhxpc";} else {print $0;} }'  par/simul2.par > par/simul.par

sh/runngrids9.sh

#Recon method
awk '{ if ( NR == 5 ) { print "wach";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 7 ) { print "rbftr";} else {print $0;} }'  par/simul2.par > par/simul.par

sh/runngrids9.sh

#Recon method
awk '{ if ( NR == 5 ) { print "wach";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 7 ) { print "rbfetr";} else {print $0;} }'  par/simul2.par > par/simul.par

sh/runngrids9.sh

#Recon method
awk '{ if ( NR == 5 ) { print "wach";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 7 ) { print "rbfetrp";} else {print $0;} }'  par/simul2.par > par/simul.par

sh/runngrids9.sh

