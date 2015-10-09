#!/bin/bash

#Script to create a link for the main executable on the main directory
version=` date +%y.%m.%d `
binname="bin/imodel$version"
echo $binname
cp bin/imodel $binname
rm -rf imodel
if [[ -f imodel ]] ; then
    ln -n -s -f $binname imodel
    echo "A link for the main executable (named 'imodel') was updated " $binname  
else
    if [[ -f bin/imodel ]] ; then
	ln -s $binname imodel
	echo "A link for the main executable (named 'imodel') was created " $binname 
    else
	echo "Could not create link, bin/imodel does not exist" $binname 
    fi
fi
echo
