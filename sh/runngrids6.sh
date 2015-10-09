#!/bin/bash

source inc_vars.sh 

./imodel par/mesh.par 12 # > log.txt
./imodel par/mesh.par 42 # > log.txt
./imodel par/mesh.par 162 # > log.txt
./imodel par/mesh.par 642 # > log.txt
./imodel par/mesh.par 2562 # > log.txt
./imodel par/mesh.par 10242 # > log.txt
./imodel par/mesh.par 40962 # > log.txt
#./imodel par/mesh.par 160000 # > log.txt
#./imodel par/mesh.par 600000 # > log.txt
