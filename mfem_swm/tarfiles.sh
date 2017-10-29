#!/bin/bash

# Srcipt to tar instalation files 

date=` date +%F `
version=` date +%y.%m.%d `
echo "Today: " $date

sourcefiles="*.f90 \
*.sh"

parfiles="nml/* \
*.in "

gmtscripts="gmt/*.sh \
gmt/*.cpt \
gmt/mapcenter.dat \
grd/*.gmt "

scripts="notused/*.sh "

docs="doc/* "

others="Makefile \
README.* \
matlab/* "

files="$sourcefiles $parfiles $scripts $gmtscripts $refs $docs $others"

output="femswe$version.tar.bz2"

tar cjfv $output $files

echo "File " $output " ready!"
echo

echo "-------------------------------------------"
