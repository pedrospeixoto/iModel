#!/bin/bash

awk '{ if ( NR < 100 || NR%10==0 ) { print $0 ;} }' $1 > $2 
