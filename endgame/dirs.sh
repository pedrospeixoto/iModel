#!/bin/bash
#Script to set directory structure 

#echo
#echo "Setting directory structure..."
#echo
#if [ -d bin ] ; then
#    echo "Binary files will be placed inside          bin/" 
#else
#    mkdir bin
#    echo "Directory 'bin' created" 
#    echo "Binary files will be placed inside          bin/" 
#fi

if [ -d dump ] ; then
    echo "Output data files will be put inside        dump/"
else
    mkdir dump
    echo "Directory 'dump' created" 
    echo "Output data files will be put inside        dump/"
fi


if [ -d bin ] ; then
    echo "Binary files will be put inside        bin/"
else
    mkdir bin
    echo "Directory 'bin' created" 
    echo "Binary files will be put inside        bin/"
fi


#echo "End of directory structuring."
echo

