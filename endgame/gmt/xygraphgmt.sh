#!/bin/bash

echo '============================'
echo 'XY graph script for GMT'
echo 'P. Peixoto - Jul 2012       '
echo '============================'
echo

killall gv

#Directory for graphs
graphdir=../graphs

echo 

#KIND OF PLOT
#echo
#echo "Select kind of plot:"
#read kplot
#echo

#Arguments check
if [ $# -eq 0  ] ; then
    echo "Please enter arguments:"
    echo " Arguments (files) should be passed in the following order:"
    echo " Data file (.txt) in ASCII, output file (graph.ps)"
    exit 0
else
    data=$1
    plot=$2
fi  



#Set naming variables

tmp=$data
data=`basename $tmp .txt`
dirdata=`dirname $tmp`
echo "Dir and file to be graphed:" $dirdata $data

#Output graph file
if [ $plot ] ; then
    echo "Output filename:" $graphdir/$plot 
else
    plot=$graphdir/
    plot=$plot$data
    baseplot=$plot
    plot=$plot.ps
fi

echo "Graph to be build:" $plot

#Pre-calculus

dscale=0.00001
gmtconvert -F0,5 -V $data.txt > ermax.tmp
gmtconvert -F0,6 -V $data.txt > er2.tmp

reg=`minmax -I$dscale ermax.tmp`
minmaxdata=`minmax $data.txt -C`
minglevel=`echo $minmaxdata | awk '{printf "%d", $1 }' `
maxglevel=`echo $minmaxdata | awk '{printf "%d", $2 }' `
echo "Min/Max Glevel:" $minglevel $maxglevel
let "maxglevel = $maxglevel + 1"
echo "Min/Max Glevel:" $minglevel $maxglevel

minerror=`echo $minmaxdata | awk '{printf "%1.8e", $13 }' `
maxerror=`echo $minmaxdata | awk '{printf "%1.8e", $12 }' `
echo "Min/Max Error:" $minerror $maxerror
#minerror=`echo "scale=8; minerror-minerror/10." | bc -l`
#echo "Adjusted Min/Max Error:" $minerror $maxerror

reg=-R$minglevel/$maxglevel/$minerror/$maxerror

#MAP REGION
width=25c  
height=25c 
map=$reg" -JX10/10l"
# scalepos = xpos/ypos/length/width[h] 
# scalepos=23.2c/5.5c/12c/0.5c

echo Map : $map

#GMT set up
. gmt_shell_functions.sh
gmtset PAPER_MEDIA Custom_${height}x${width}
gmtset ANNOT_FONT_SIZE +12p


#Basemap - edit -Bxgx to get lines where wanted
psbasemap  $map -Ba1f1:"glevel":/a1f1:"Error":WS -K -V \
    --GRID_PEN_PRIMARY=faint,gray,-  --D_FORMAT="%0.1e" > $plot



#Error Max 
psxy $map -Wthin,black -O -K -V ermax.tmp >> $plot
psxy $map -W5,black -Sc0.2c -Gblack -O -K -V ermax.tmp >> $plot

#Error 2
psxy $map -Wthin,black -O -K -V er2.tmp >> $plot
psxy $map -W5,black -Ss0.2c -Gblack -O -K -V er2.tmp >> $plot


#Legend Elements
cat > legend.gmt << END
S 0.1c kcircleline 0.2c black 0.25p 0.5 Error Max
>
END

pslegend legend.gmt $map -O -Dx10c/10c/7c/1.7c/TC >> $plot



echo
echo "Plotted file:" $plot

echo "and pdf and eps files (for pdflatex and latex)."
ps2raster -Te -P -A $plot 
ps2pdf $plot $baseplot.pdf
rm $plot
gv $baseplot.eps &

echo
echo "-----------------------------------------------------------"