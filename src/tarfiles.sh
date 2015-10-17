#!/bin/bash

# Srcipt to tar instalation files of imodel

date=` date +%F `
version=` date +%y.%m.%d `
echo "Today: " $date

sourcefiles="src/*.f90 \
src/*.sh"

parfiles="par/*.par "

gmtscripts="gmt/*.sh \
gmt/*.cpt \
gmt/mapcenter.dat "

scripts="sh/*.sh "

refs="ref/*.f90 \
ref/*.pdf \
ref/Makefile "

docs="doc/*.tex \
doc/*.bib \
doc/*.pdf "

others="Makefile \
README.* \
matlab/*.m "

files="$sourcefiles $parfiles $scripts $gmtscripts $refs $docs $others"

output="imodel$version.tar.bz2"

tar cjfv $output $files

echo "File " $output " ready!"
echo

echo "-------------------------------------------"
