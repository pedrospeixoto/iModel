#!/bin/bash

#-----------------------------------------
# Script to animate scalar field on icosaedral mesh
#
# Pedro da Silva Peixoto
# Feb 2013
#
#-----------------------------------------

echo 
echo '   Animate Scalar Transport  '
echo '-----------------------------'

killall gv
. gmt_shell_functions.sh

#Directory for graphs
graphdir=../graphs

#First argument - what is to be animated
name=$1 
if [[ $1 ]] ; then
    echo "Plotting test for: " $name
else
    echo "Please use as first argument the kind of test from a file until _t"
    echo "Ex: '../data/transp_case3_...phi_t"
    exit 0
fi

#Second argument should be the mesh type
mesh=$2 # icoXXX
if [[ $2 ]] ; then
    echo "Mesh considered : " $mesh
else
    echo "Please use as second argument the kind of mesh that fills the rest of the filename"
    echo " without the extension"
    echo "Ex: 'icos_eqs_nopt_5' "
    exit 0
fi

#Set file names and directories
datadir=`dirname $name`
echo "Directory for data:" $datadir

#Get basename
base=`basename $name`
echo "Basename: " $base

#Make a list of the files to be plotted
ls $name*$mesh.dat -t -r > files.txt  

#first record is the newest file (lastest timestep)
timesteps=`grep -c "$name" files.txt` 
one=1
let timesteps=$timesteps-$one
echo "Timesteps:" $timesteps
if [[ $timesteps == 0 ]] ; then
    echo "No files to be plotted "
    exit 0
fi

#get file whit the last timestep
maxtimefile=`awk 'END{print}' files.txt `
echo "Max time file:" $maxtimefile

#remove the end and the directory
tmp=`basename $maxtimefile _$mesh.dat`
echo $tmp


#Get max time
maxtime=`echo $tmp | sed s:"$base"::g`
echo "Max time:" $maxtime

#maxtime=512
let "half= $maxtime/2 " 
let "dtime= $maxtime/$timesteps " #6 #120 #6 #120

if [ $dtime ] ; then
    echo "Dtime=" $dtime
else
    echo "Error in timestep. Dtime:" $dtime
    exit 0
fi

kone=0
echo "Subtract one on field values?"
echo " 0) NO"
echo " 1) YES"
#read kone
echo    

kabs=0
echo "Use absolute value?"
echo " 0) NO"
echo " 1) YES"
#read kabs
echo    

echo "Select color pallet:"
echo " 1) White-DarkBlue"
echo " 2) Blue-White-Red"
echo " 3) Multicolor jet"
echo " 4) White-Light Blue"
echo " 5) Seismology - Multicolor"
read kcolor
echo    

echo "Use fixed scale:"
echo " 0) NO"
echo " 1) YES"
read fixedscale
echo    

#   Create grid files
#-------------------------------------------------------

#Set Counters
time=0
minval=100000
maxval=-100000
#Creat grid files
while [ $time -le $maxtime ]; do

    filename=$base$time"_"$mesh

    echo Creating grid file: $datadir/$filename.grd 

    xyz2grd $datadir/$filename.dat -I0.25 -Rg -bi -F -G$datadir/$filename.grd 

    if [ $kone -eq 1 ] ; then
       #Subtract one from values
	grdmath $datadir/$filename.grd -1 ADD = $datadir/$filename.grd
    fi

    if [ $kabs -eq 1 ] ; then
       #Set absolute value 
	grdmath $datadir/$filename.grd ABS = $datadir/$filename.grd
    fi

    minvaltmp=`grdinfo $datadir/$filename.grd -C | awk '{printf "%12.8e", $6 }' `
    maxvaltmp=`grdinfo $datadir/$filename.grd -C | awk '{printf "%12.8e", $7 }' `
    echo Time: $time Min: $minvaltmp Max: $maxvaltmp
    minval=`perl -e "if ( $minval > $minvaltmp ) { print $minvaltmp; } else { print $minval; }"`
    maxval=`perl -e "if ( $maxval < $maxvaltmp ) { print $maxvaltmp; } else { print $maxval; }"`
    let time=$time+$dtime

done

echo "Min:" $minval "Max:" $maxval


#Map parameters
#---------------------------

#Map center
#read lon lat < mapcenter.dat
lat=0
lon=0
#echo Map Center: $lat $lon

#Dots per inch on tiff file
dpi=600 #175

#Page size
width=29c  
height=16c 

#GMT settings
. gmt_shell_functions.sh
gmtset PAPER_MEDIA Custom_${height}x${width}
gmtset ANNOT_FONT_SIZE +12p

#Region - Global
reg="-R-180/180/-90/90"
#reg="-Rg"

#Mercator
map="-R-180/180/-80/80 -JM16"

#Linear
map=$reg" -Jx0.06d"

#Map
echo "Map used:" $map
#echo

#      Set scale and color pallet
# --------------------------------------------

# scalepos = xpos/ypos/length/width[h] 
#scalepos=18c/6c/12c/0.5c
scalepos=23c/5.5c/12c/0.5c

minval2=`echo $minval | sed 's/e/\\*10\\^/' | sed 's/+//'`
maxval2=`echo $maxval | sed 's/e/\\*10\\^/' | sed 's/+//'`

scale=`echo "scale=10; ( $maxval2 - $minval2 )/100" | bc -l`
scaleopt="-T$minval/$maxval/$scale"
color=colorsanim.cpt

echo "Scale ticks:" $scale
echo "Scale min/max/d :" $scaleopt
echo


if [ $kcolor -eq 1 ] ; then  #White-Darkblue
    makecpt -Cdarkblue -I -T$minval/$maxval/$scale -Z \
    	--COLOR_BACKGROUND=blue  --COLOR_FOREGROUND=white > $color
    colscale=darkblue
fi

if [ $kcolor -eq 2 ] ; then #Blue - Red
    postest=`echo "scale=12; (-0.0000000001 <= $minval2)" | bc -l `
    if [ $postest -eq 0 ] ; then
	maxgtabsmin=0
	maxgtabsmin=`echo "scale=12; (-1 * $minval2 <= $minval2)" | bc -l `
	echo max gt abs min ? $maxgtabsmin
	if [ $maxgtabsmin -eq 0 ] ; then #abs negative value is larger
	    a=$minval
	    b=`echo "scale=8; ( -1* $minval2 )" | bc -l `
	else
	    a=`echo "scale=8; ( -1* $maxval2 )" | bc -l `
	    b=$maxval
	fi
	makecpt -Cbluered -T$a/$b/$scale -Z -V \
    	    --COLOR_BACKGROUND=blue  --COLOR_FOREGROUND=red > $color
    else
        makecpt -Cpolar -T$minval/$maxval/$scale -Z \
    	    --COLOR_BACKGROUND=blue  --COLOR_FOREGROUND=red > $color
    fi
    colscale=bluered
fi


if [ $kcolor -eq 3 ] ; then #Multicolor jet
    makecpt -Cjet -T$minval/$maxval/$scale -Z  > $color
    colscale=jet
fi

if [ $kcolor -eq 4 ] ; then #Light Blue -White
    makecpt -Cblue -T$minval/$maxval/$scale -Z  > $color
    colscale=blue
fi


if [ $kcolor -eq 5 ] ; then #Multicolor jet
    #makecpt -Cseis -I -T$minval/$maxval/$scale -Z   > $scalar.cpt
    makecpt -Cmulticol $scaleopt -Z -V  > $color
    colscale=multicol
    	#--COLOR_BACKGROUND=0/0/127  --COLOR_FOREGROUND=127/0/0
        #grd2cpt $scalar.grd -Cjet  -Z -V  > $scalar.cpt
fi    

scale=`echo "scale=8; ( $maxval2 - $minval2 )/10" | bc -l`

#  Create Images
#----------------------------------------------------

time=0
count=1000
while [ $time -le $maxtime ]; do
    filename=$base$time"_"$mesh
    echo Creating image for: $filename

    #Set basemap
    #psbasemap  $map -B40g40/20 -K -U"$filename" > $graphdir/$filename.ps
    psbasemap  $map -B40g40/20 -K  > $graphdir/$filename.ps

    #Variable color scale
    if [ $fixedscale -eq 0 ] ; then
	if [ $kcolor -eq 2 ] ; then #Blue - Red
	    grd2cpt $datadir/$filename.grd -C$colscale  -Z -T=  > $color
	else
	    grd2cpt $datadir/$filename.grd -C$colscale  -Z  > $color
	fi
    fi

    #Draw filled contourplot
    grdimage $datadir/$filename.grd -C$color $map -O -K >> $graphdir/$filename.ps

    #Set tickmarks distance for scale
    psscale -C$color -D$scalepos -B$scale -O -K  \
	--D_FORMAT="%0.2e">> $graphdir/$filename.ps
    #psscale -C$color -D$scalepos -B$scale -O -K  >> $graphdir/$filename.ps

    #Draw world map
    pscoast  $map -W0.25p  -O  >> $graphdir/$filename.ps

    #RIP to TIFF at specified dpi
    ps2raster -P -E$dpi -Tt $graphdir/$filename.ps -Ggs

    #Convert to jpeg
    convert -quality 100 -resize 16% $graphdir/$filename.tif $graphdir/$base"_"$mesh"_"$count.jpg
    
    #if [ "$time" -eq "0" -o  "$time" -eq "1" -o  "$time" -eq "2" -o  "$time" -eq "3" -o \
    #"$time" -eq "$half" -o "$time" -eq "$maxtime" ] ; then
    if [ "$time" -eq "0" -o "$time" -eq "$half" -o "$time" -eq "$maxtime" ] ; then
	echo "Keeping file (.eps and .pdf)" 
	ps2raster -Te -P -A $graphdir/$filename.ps
	ps2raster -Tf -P -A $graphdir/$filename.ps
    fi
    
    rm  -r $graphdir/$filename.ps
    rm  -r $datadir/$filename.grd
    rm  -r $graphdir/$filename.tif

    let time=$time+$dtime
    let count=$count+$one

done

echo "Creating high quality animated gif : " $graphdir/$base"_"$mesh.gif
echo "       please be patient ..."
#  Create animated GIF file 
convert -delay 20 -quality 100 $graphdir/$base"_"$mesh*.jpg $graphdir/$base"_"$mesh.gif

# echo "Creating high quality animated AVI : " $graphdir/$base"_"$mesh"_quality".avi
# echo "       please be patient ..."
# mencoder mf://$graphdir/*.jpg -mf w=800:h=600:fps=5:type=jpg -ovc copy -oac copy -quiet -o $graphdir/$base"_"$mesh"_quality".avi

# echo "Creating small size movie (avi) : " $graphdir/$base"_"$mesh.avi
# echo "       please be patient ..."
# avconv -loglevel quiet -r 5 -i $graphdir/$base"_"$mesh"_"1%03d.jpg -y -an $graphdir/$base"_"$mesh.avi

#rm -r $graphdir/$base"_"$mesh"_"*.jpg





