#!/bin/bash

# Creates a tarball

date=` date +%F `
version=` date +%y.%m.%d `
echo "Today: " $date

sourcefiles="../src/*"

reffiles="../ref/*"

parfiles="../par/* "

griddir="../grid/* "

altfiles="../altitude/* "

scripts="../sh/* "

pyscripts="../pyscripts/*.py "

others="../Makefile\
 ../README.* ../*.sh"

files="$sourcefiles $parfiles $scripts $pyscripts $altfiles $others $griddir"
#files="$sourcefiles $parfiles $scripts $pyscripts $others"

output="imodel.tar.bz2"

tar cjfv $output $files

echo "File " $output " ready!"
echo
echo "-------------------------------------------"
