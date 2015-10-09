#!/bin/bash

# Script for building several meshs
#-------------------------------------------------------------------------

source sh/inc_vars.sh 


awk '{ if ( NR == 9 ) { print "nopt";} else {print $0;} }'  par/mesh.par > par/mesh2.par
cp par/mesh2.par par/mesh.par
./runngrids.sh

awk '{ if ( NR == 9 ) { print "scvt";} else {print $0;} }'  par/mesh.par > par/mesh2.par
cp par/mesh2.par par/mesh.par
./runngrids.sh
