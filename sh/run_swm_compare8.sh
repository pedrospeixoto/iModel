#!/bin/bash

source inc_vars.sh 

# Script for multiple shallow water equation tests

cp par/swmnew.par par/swm.par
sh/runngridsHR95_8.sh
sh/runngridsSCVT_8.sh

cp par/swmtrsk.par par/swm.par
sh/runngridsHR95_8.sh
sh/runngridsSCVT_8.sh
