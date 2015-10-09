#!/bin/bash

source inc_vars.sh 

# Script for multiple shallow water equation tests

#Test Case 2
awk '{ if ( NR == 3 ) { print "2";} else {print $0;} }'  par/swm.par > par/swm2.par
awk '{ if ( NR == 3 ) { print "2";} else {print $0;} }'  par/swm2.par > par/swm.par

#HTC TRSK Method
awk '{ if ( NR == 9 ) { print "HTC";} else {print $0;} }'  par/swm.par > par/swm2.par
awk '{ if ( NR == 11 ) { print "trsk";} else {print $0;} }'  par/swm2.par > par/swm.par
awk '{ if ( NR == 13 ) { print "trsk";} else {print $0;} }'  par/swm.par > par/swm2.par
awk '{ if ( NR == 15 ) { print "trsk";} else {print $0;} }'  par/swm2.par > par/swm.par
awk '{ if ( NR == 17 ) { print "trsk";} else {print $0;} }'  par/swm.par > par/swm2.par
awk '{ if ( NR == 17 ) { print "trsk";} else {print $0;} }'  par/swm2.par > par/swm.par

sh/runngrids9.sh
sh/runngridsHR95HK_9.sh
sh/runngridsSCVT_9.sh


#HC TRSK Method
awk '{ if ( NR == 9 ) { print "HC";} else {print $0;} }'  par/swm.par > par/swm2.par
awk '{ if ( NR == 15 ) { print "bary";} else {print $0;} }'  par/swm2.par > par/swm.par

sh/runngrids9.sh
sh/runngridsHR95HK_9.sh
sh/runngridsSCVT_9.sh

#HC TRSK Method
awk '{ if ( NR == 9 ) { print "HC";} else {print $0;} }'  par/swm.par > par/swm2.par
awk '{ if ( NR == 15 ) { print "trsk";} else {print $0;} }'  par/swm2.par > par/swm.par

sh/runngrids9.sh
sh/runngridsHR95HK_9.sh
sh/runngridsSCVT_9.sh


#HTC TRSK Method with perot
awk '{ if ( NR == 9 ) { print "HTC";} else {print $0;} }'  par/swm.par > par/swm2.par
awk '{ if ( NR == 11 ) { print "perhx";} else {print $0;} }'  par/swm2.par > par/swm.par

sh/runngrids9.sh
sh/runngridsHR95HK_9.sh
sh/runngridsSCVT_9.sh

#HC TRSK Method with perot
awk '{ if ( NR == 9 ) { print "HC";} else {print $0;} }'  par/swm.par > par/swm2.par
awk '{ if ( NR == 11 ) { print "trsk";} else {print $0;} }'  par/swm2.par > par/swm.par
awk '{ if ( NR == 13 ) { print "trsk";} else {print $0;} }'  par/swm.par > par/swm2.par
awk '{ if ( NR == 11 ) { print "perhx";} else {print $0;} }'  par/swm2.par > par/swm.par
sh/runngrids9.sh
sh/runngridsHR95HK_9.sh
sh/runngridsSCVT_9.sh

#HC PEROT Method
awk '{ if ( NR == 9 ) { print "HC";} else {print $0;} }'  par/swm.par > par/swm2.par
awk '{ if ( NR == 11 ) { print "perhx";} else {print $0;} }'  par/swm2.par > par/swm.par
awk '{ if ( NR == 13 ) { print "pered";} else {print $0;} }'  par/swm.par > par/swm2.par
awk '{ if ( NR == 15 ) { print "bary";} else {print $0;} }'  par/swm2.par > par/swm.par

sh/runngrids9.sh
sh/runngridsHR95HK_9.sh
sh/runngridsSCVT_9.sh

#HC DTRED Method
awk '{ if ( NR == 13 ) { print "dtred";} else {print $0;} }'  par/swm.par > par/swm2.par
awk '{ if ( NR == 15 ) { print "bary";} else {print $0;} }'  par/swm2.par > par/swm.par

sh/runngrids9.sh
sh/runngridsHR95HK_9.sh
sh/runngridsSCVT_9.sh

#HTC DTRED Method
awk '{ if ( NR == 9 ) { print "HTC";} else {print $0;} }'  par/swm.par > par/swm2.par
awk '{ if ( NR == 11 ) { print "trsk";} else {print $0;} }'  par/swm2.par > par/swm.par
awk '{ if ( NR == 13 ) { print "dtred";} else {print $0;} }'  par/swm.par > par/swm2.par
awk '{ if ( NR == 15 ) { print "trsk";} else {print $0;} }'  par/swm2.par > par/swm.par
awk '{ if ( NR == 17 ) { print "trsk";} else {print $0;} }'  par/swm.par > par/swm2.par
awk '{ if ( NR == 17 ) { print "trsk";} else {print $0;} }'  par/swm2.par > par/swm.par

sh/runngrids9.sh
sh/runngridsHR95HK_9.sh
sh/runngridsSCVT_9.sh

#HC DTRED Method
awk '{ if ( NR == 9 ) { print "HC";} else {print $0;} }'  par/swm.par > par/swm2.par
awk '{ if ( NR == 11 ) { print "trsk";} else {print $0;} }'  par/swm2.par > par/swm.par
awk '{ if ( NR == 13 ) { print "dtred";} else {print $0;} }'  par/swm.par > par/swm2.par
awk '{ if ( NR == 15 ) { print "trsk";} else {print $0;} }'  par/swm2.par > par/swm.par
awk '{ if ( NR == 17 ) { print "trsk";} else {print $0;} }'  par/swm.par > par/swm2.par
awk '{ if ( NR == 17 ) { print "trsk";} else {print $0;} }'  par/swm2.par > par/swm.par

sh/runngrids9.sh
sh/runngridsHR95HK_9.sh
sh/runngridsSCVT_9.sh

#HTC PERED Method
awk '{ if ( NR == 9 ) { print "HTC";} else {print $0;} }'  par/swm.par > par/swm2.par
awk '{ if ( NR == 11 ) { print "trsk";} else {print $0;} }'  par/swm2.par > par/swm.par
awk '{ if ( NR == 13 ) { print "pered";} else {print $0;} }'  par/swm.par > par/swm2.par
awk '{ if ( NR == 15 ) { print "trsk";} else {print $0;} }'  par/swm2.par > par/swm.par
awk '{ if ( NR == 17 ) { print "trsk";} else {print $0;} }'  par/swm.par > par/swm2.par
awk '{ if ( NR == 17 ) { print "trsk";} else {print $0;} }'  par/swm2.par > par/swm.par

sh/runngrids9.sh
sh/runngridsHR95HK_9.sh
sh/runngridsSCVT_9.sh

#HTC PERED Method
awk '{ if ( NR == 9 ) { print "HC";} else {print $0;} }'  par/swm.par > par/swm2.par
awk '{ if ( NR == 11 ) { print "trsk";} else {print $0;} }'  par/swm2.par > par/swm.par
awk '{ if ( NR == 13 ) { print "pered";} else {print $0;} }'  par/swm.par > par/swm2.par
awk '{ if ( NR == 15 ) { print "trsk";} else {print $0;} }'  par/swm2.par > par/swm.par
awk '{ if ( NR == 17 ) { print "trsk";} else {print $0;} }'  par/swm.par > par/swm2.par
awk '{ if ( NR == 17 ) { print "trsk";} else {print $0;} }'  par/swm2.par > par/swm.par

sh/runngrids9.sh
sh/runngridsHR95HK_9.sh
sh/runngridsSCVT_9.sh
