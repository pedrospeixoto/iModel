#!/bin/bash

# Script for multiple vector reconstruction tests
#cat par/simul.par

source inc_vars.sh 

#Recon method - nonehx
awk '{ if ( NR == 5 ) { print "nonehx";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 7 ) { print "nonehx";} else {print $0;} }'  par/simul2.par > par/simul.par
#./imodel
sh/runngrids9.sh

#Recon method - nonetr
awk '{ if ( NR == 5 ) { print "nonetr";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 7 ) { print "nonetr";} else {print $0;} }'  par/simul2.par > par/simul.par
#./imodel
sh/runngrids9.sh


#Recon method - perot tr
awk '{ if ( NR == 5 ) { print "pertr";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 7 ) { print "pertr";} else {print $0;} }'  par/simul2.par > par/simul.par
#./imodel
sh/runngrids9.sh

#Recon method - perot hx
awk '{ if ( NR == 5 ) { print "perhx";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 7 ) { print "perhx";} else {print $0;} }'  par/simul2.par > par/simul.par
#./imodel
sh/runngrids9.sh

#Recon method - ebb / kls
awk '{ if ( NR == 5 ) { print "kls";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 7 ) { print "kls";} else {print $0;} }'  par/simul2.par > par/simul.par
#./imodel
sh/runngrids9.sh

#Recon method - whitney
awk '{ if ( NR == 5 ) { print "wht";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 7 ) { print "wht";} else {print $0;} }'  par/simul2.par > par/simul.par
#./imodel
sh/runngrids9.sh

#Recon method - lsqhxe
awk '{ if ( NR == 5 ) { print "lsqhxe";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 7 ) { print "lsqhxe";} else {print $0;} }'  par/simul2.par > par/simul.par
#./imodel
sh/runngrids9.sh

#Recon method - lsqtrc
awk '{ if ( NR == 5 ) { print "lsqtrc";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 7 ) { print "lsqtrc";} else {print $0;} }'  par/simul2.par > par/simul.par
#./imodel
sh/runngrids9.sh

#Recon method - hybrid
#############################
awk '{ if ( NR == 5 ) { print "lintrv";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 7 ) { print "perhx";} else {print $0;} }'  par/simul2.par > par/simul.par
awk '{ if ( NR == 9 ) { print "lsqhxe";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 17 ) { print "0.01 val";} else {print $0;} }'  par/simul2.par > par/simul.par
#./imodel
sh/runngrids9.sh

awk '{ if ( NR == 9 ) { print "0";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ print $0; }'  par/simul2.par > par/simul.par


#RBF Recon methods
############################3

awk '{ if ( NR == 5 ) { print "rbfhx";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 7 ) { print "rbfhx";} else {print $0;} }'  par/simul2.par > par/simul.par
awk '{ if ( NR == 19 ) { print "2 0";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ print $0; }'  par/simul2.par > par/simul.par
#./imodel
sh/runngrids9.sh

#Recon method
awk '{ if ( NR == 5 ) { print "rbfhxpc";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 7 ) { print "rbfhxpc";} else {print $0;} }'  par/simul2.par > par/simul.par
#./imodel
sh/runngrids9.sh

#Recon method
awk '{ if ( NR == 5 ) { print "rbftr";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 7 ) { print "rbftr";} else {print $0;} }'  par/simul2.par > par/simul.par
#./imodel
sh/runngrids9.sh

#Recon method
awk '{ if ( NR == 5 ) { print "rbfetr";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 7 ) { print "rbfetr";} else {print $0;} }'  par/simul2.par > par/simul.par
#./imodel
sh/runngrids9.sh

#Recon method
awk '{ if ( NR == 5 ) { print "rbfetrp";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 7 ) { print "rbfetrp";} else {print $0;} }'  par/simul2.par > par/simul.par
#./imodel
sh/runngrids9.sh

#Do not do plots for other runs
awk '{ if ( NR == 23 ) { print "0";} else {print $0;} }'  par/simul.par > par/simul2.par
awk '{ if ( NR == 23 ) { print "0";} else {print $0;} }'  par/simul2.par > par/simul.par

