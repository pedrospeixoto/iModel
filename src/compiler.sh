#!/bin/bash
#------------------------------------------------------------
# Check Operating System
# Check compiler to be used
# Define Fortran flags
#
# Exports variables, should be run using double dots: . ./compiler.sh
#
# Author: P. Peixoto Jul 2012
#------------------------------------------------------------

#System check
echo
OS=`uname`
echo "Operating System:" $OS
echo

#Compiler Check
echo "Defining compiler..."
echo

F90=`command -v ifort`
if [ -f $F90 ] ; then
    echo "IFORT exists! " $F90
    F90="ifort"
else
    F90="gfortran"
fi
echo "Using compiler:" $F90
export F90=$F90

#Debug flag
#Set variable for debugging (1) or not (0)
DEBUG="0"

#Fortran flags depending on system and compiler
#LINUX
if test $OS = "Linux"  ; then
    echo "Using Linux Flags ..."
    #IFORT FLAGS
    if test $F90 = "ifort" ; then
	if test $DEBUG = "1" ; then
            #DEBUG FLAGS FOR ECLIPSE
            FFLAG=" -traceback -debug extended"
	else
	    #Optimized flags for ifort
	    FFLAG=" -openmp -O3 -ipo -parallel -xHOST"
	fi
    else
    #GFORTRAN FLAGS
	if test $DEBUG = "1" ; then
            #DEBUG FLAGS FOR ECLIPSE
            FFLAG=" -O0 -g"
	else
	    #Optimized flags for gfortran
	    FFLAG="-O3 -fopenmp"
	fi
    fi
else
#WINDOWS
    echo "Using Windows/Cygwin Flags ..."
    #IFORT FLAGS
    if test $F90 = "ifort" ; then
	if test $DEBUG = "1" ; then
            #DEBUG FLAGS FOR ECLIPSE
            FFLAG=" -debug extended " #Needs testing
	else
	    #Optimized flags for ifort
	    FFLAG="-Qopenmp -O3 -Qipo -Qparallel -QxHOST"
            #-Qguide4 -fast -Qparallel  -Qpar-report1 
            #DO NOT USE -fast : it disables prec-div and ruins accuracy
	fi
    else
    #GFORTRAN FLAGS
	if test $DEBUG = "1" ; then
            #DEBUG FLAGS FOR ECLIPSE
            FFLAG=" -O0 -g"
	else
	    #Optimized flags for gfortran
	    FFLAG="-O3 -fopenmp"
	fi
    fi
fi
export FFLAG="$FFLAG"
echo "  flags: " $FFLAG

#-------------------------------------------------------------------