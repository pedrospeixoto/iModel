#!/bin/bash

source inc_vars.sh 

# Script for convergence analysis on transport model for different methods
#-------------------------------------------------------------------------


awk '{ if ( NR == 15 ) { print "lsqtrc";} else {print $0;} }'  par/trans.par > par/trans2.par
cp par/trans2.par par/trans.par
./runngrids.sh
#./runntimesgrids.sh

awk '{ if ( NR == 15 ) { print "lsqhxe";} else {print $0;} }'  par/trans.par > par/trans2.par
cp par/trans2.par par/trans.par
./runngrids.sh
#./runntimesgrids.sh

awk '{ if ( NR == 15 ) { print "perhx";} else {print $0;} }'  par/trans.par > par/trans2.par
cp par/trans2.par par/trans.par
./runngrids.sh
#./runntimesgrids.sh

awk '{ if ( NR == 15 ) { print "+1perhx+2lintrv+3lsqhxe+";} else {print $0;} }'  par/trans.par > par/trans2.par
cp par/trans2.par par/trans.par
./runngrids.sh
#./runntimesgrids.sh
