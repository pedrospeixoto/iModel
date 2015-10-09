#!/bin/bash

# Script for multiple deformational tests
#------------------------------------------
source inc_vars.sh 

echo "            Basic Parameters     "
echo "---------------------------------"
echo " "
cat deform.par

#Test case
awk '{ if ( NR == 3 ) { print 5;} else {print $0;} }'  deform.par > deform2.par
#Initial Condition Field
awk '{ if ( NR == 5 ) { print 4;} else {print $0;} }'  deform2.par > deform.par

cat deform.par

#./imodel
# runtests
./runntimes.sh

#Change trajectory
awk '{ if ( NR == 7 ) { print 0;} else {print $0;} }'  deform.par > deform2.par
cp deform2.par deform.par

cat deform.par
#./imodel
# runtests
./runntimes.sh


#Test case
awk '{ if ( NR == 3 ) { print 5;} else {print $0;} }'  deform.par > deform2.par
#Initial Condition Field
awk '{ if ( NR == 5 ) { print 4;} else {print $0;} }'  deform2.par > deform.par
#Change trajectory
awk '{ if ( NR == 7 ) { print 1;} else {print $0;} }'  deform.par > deform2.par
cp deform2.par deform.par

cat deform.par

#./imodel
# runtests
./runntimes.sh

#Change trajectory
awk '{ if ( NR == 7 ) { print 0;} else {print $0;} }'  deform.par > deform2.par
cp deform2.par deform.par

cat deform.par
#./imodel
# runtests
./runntimes.sh

