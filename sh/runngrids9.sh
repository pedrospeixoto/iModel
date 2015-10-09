#!/bin/bash

source inc_vars.sh 

awk '{ if ( NR == 5 ) { print "icos";} else {print $0;} }'  par/mesh.par > par/mesh2.par
awk '{ if ( NR == 9 ) { print "nopt";} else {print $0;} }'  par/mesh2.par > par/mesh.par


./imodel par/mesh.par 12 # > log.txt
./imodel par/mesh.par 42 # > log.txt
./imodel par/mesh.par 162 # > log.txt
./imodel par/mesh.par 642 # > log.txt
./imodel par/mesh.par 2562 # > log.txt
./imodel par/mesh.par 10242 # > log.txt
./imodel par/mesh.par 40962 # > log.txt
./imodel par/mesh.par 160000 # > log.txt
./imodel par/mesh.par 600000 # > log.txt
./imodel par/mesh.par 2000000 # > log.txt
