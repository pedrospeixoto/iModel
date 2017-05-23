#!/bin/bash

echo '============================'
echo 'Convert pdf/eps to grayscale'
echo 'P. Peixoto - May 2017       '
echo '============================'
echo

#Close all open gv windows (to avoid a bunch of open stuff)
killall gv


#Arguments check
if [ $# -eq 0  ] ; then
    echo "Please enter the file that needs to be converted in either .pdf or .eps"
    exit 0
fi  

#Set naming variables
original=$1
filename=$(basename "$original")
extension="${filename##*.}"
filename="${filename%.*}"
in=$filename"."$extension
echo "IN: " $in
echo "basename: " $filename
echo "extension: " $extension
echo

outbase=$filename"_grey"
out=$filename"_grey".$extension
echo "OUT: " $out
echo "basename: " $outbase
cp $in $out

if [ "$extension" == "eps" ]; then
	echo "Converting to pdf..."
	epstopdf $out
fi

gs -dNOPAUSE -dBATCH -dSAFER -sDEVICE=pdfwrite -dProcessColorModel=/DeviceGray -dColorConversionStrategy=/Gray -dPDFUseOldCMS=false  -dNOCACHE -o $outbase"tmp.pdf" -f $outbase".pdf"
mv $outbase"tmp.pdf" $outbase".pdf"
if [ "$extension" == "eps" ]; then
	pdf2ps -eps $outbase".pdf" $outbase".eps"
fi


#if [ "$extension" == "pdf" ]; then
#    echo "Converting PDF file ..."
#    gs -dNOPAUSE -dBATCH -dSAFER -sDEVICE=pdfwrite -dProcessColorModel=/DeviceGray -dColorConversionStrategy=/Gray -dPDFUseOldCMS=false  -dNOCACHE -o $out -f $filename.pdf
#fi
#if [ "$extension" == "eps" ]; then
#    echo "Converting EPS file ..."
#    gs -dNOPAUSE -dBATCH -dSAFER -sDEVICE=pdfwrite -dProcessColorModel=/DeviceGray -dColorConversionStrategy=/Gray -dPDFUseOldCMS=false  -dNOCACHE -o $out -f $filename.eps
#fi
#gmt psconvert -Te -A -P $baseplot"_grey".ps  # EPS
#pdf2ps -eps $baseplot"_grey.ps" $baseplot"_grey.eps"

#pdftops  $baseplot"_grey.pdf" $baseplot"_grey.eps"


echo
echo "-----------------------------------------------------------"
