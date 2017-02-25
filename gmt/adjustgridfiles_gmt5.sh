#!/bin/bash
echo '============================'
echo 'GMT cannot have spaces in comment and annotation lines'
echo ' All spaces in front of # and > will be removed in all'
echo '     files in grid folder'
echo 'P. Peixoto - Feb 2017       '
echo '============================'
echo

cd ..
if [ -d grd ]; then
    cd grd
elif [ -d grid ]; then
    cd grid
else
    echo "No grid directory found"
    exit 0
fi

sed -i s/" #"/"#"/g *.gmt
sed -i s/" >"/">"/g *.gmt

echo "Done"
