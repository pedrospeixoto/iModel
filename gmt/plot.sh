#!/bin/bash

echo '============================'
echo 'Plotting script for GMT'
echo 'P. Peixoto - Jun 2017       '
echo '============================'
echo

#Close all open gv windows (to avoid a bunch of open stuff)
killall gv

#Directory for graphs
#graphdir=graphs
graphdir=../graphs

gmt5='gmt'
gmt4='GMT'
if type "$gmt5" > /dev/null; then
    gmt=$gmt5
    echo "Using GMT 5"
elif type "$gmt4" > /dev/null; then
    gmt=$gmt4
    echo "Using GMT 4"
elif type "psxy" > /dev/null; then
    gmt=""
    echo "Legacy gmt4  version"
else
    echo "GMT not installed"
    exit 0
fi

echo 

#KIND OF PLOT
echo
echo "Select kind of plot:"
echo " 1) Mesh plot only"
echo " 2) Scalar field only"
echo " 3) Vector Field Only"
echo " 4) Mesh and scalar field"
echo " 5) Mesh and vector field"
echo " 6) Scalar and vector field"
echo " 7) Mesh, scalar and vector fields"
read kplot
echo

#Arguments check
if [ $# -eq 0  ] ; then
    echo "Please enter arguments:"
    echo " Arguments (files) should be passed in the following order:"
    echo " Mesh (nodes file .gmt), scalar (.dat), vector (.dat), output filename (.ps) "
    echo " Grid files and vector fields are ASCII. Scalar fields must be binary files"
    exit 0
fi  

#MAP PROJECTION
echo "Select map projection:"
echo " 1) Lambert - AZIMUTHAL (-JA) - Global"
echo " 2) Linear (-JX) - Global"  #Mercator - CYLINDRICAL (-JM)"
echo " 3) Mercator - CYLINDRICAL (-JM) - Regional only"
echo " 4) Polar Stereographic (-Js) - Regional only"
echo " 5) Linear (-JX) - North Hem"  #Mercator - CYLINDRICAL (-JM)"
read kmap
echo


#MAP CENTER
#lon=0 #-179.5 #1.23296 #31.7174
#lat=0 #-80 #-89.5 #-61.30099 #-0.0374
lon=0 #17.5
lat=0
#lon=30.125
#lat=-64.875
#lon=90.5 #1.23296 #31.7174
#lat=-0.5 #-89.5 #-61.30099 #-0.0374

echo $lon $lat > mapcenter.dat

#MAP REGION
if [ $kmap -eq 1 ] ;then #Azimutal projection - global
    if [ $kplot -eq 1 ] ;then #square box
	width=15c  
	height=15c 
    else #if colorbar scale needed, set space 
	width=22c  
	height=15c 
    fi
    # scale colorbar position
    # scalepos = xpos/ypos/length/width[h] 
    scalepos=14.5c/6.5c/12c/0.5c
    #set global region "g"
    reg="-Rg"
    #reg="-R-1/-89/3/-84r"
    map=$reg" -JA$lon/$lat/14"
    echo
    echo "Map Center (lon lat):" $lon $lat
    echo 
fi

if [ $kmap -eq 2 ] ;then #Linear projection - global
    if [ $kplot -eq 1 ] ;then #Only mesh plot
	width=16c  
	height=16c 
    else #Space for scale needed
	width=29c  
	height=16c 
    fi
    # scalepos = xpos/ypos/length/width[h] 
    scalepos=22.3c/5.5c/11c/0.5c
    #Global
    reg="-R-180/180/-90/90"
    #reg="-R-180/180/0/90"
    #map=$reg" -JX16/8"
    map=$reg" -Jx0.06d"    
fi

if [ $kmap -eq 3 ] ;then  #Local Mercator projection
    if [ $kplot -eq 1 ] ;then #Only mesh plot
	width=15c  
	height=15c 
    else      #Space for scale needed
	width=19c  
	#width=15c  
	height=15c 
    fi
    # scalepos = xpos/ypos/length/width[h] 
    scalepos=12.9c/5.5c/11c/0.5c

    #Range of map
    dlat=20
    dlon=20
    latmin=`echo "scale=12; ( $lat - $dlat)" | bc -l `
    latmax=`echo "scale=12; ( $lat + $dlat)" | bc -l `
    lonmin=`echo "scale=12; ( $lon - $dlon)" | bc -l `
    lonmax=`echo "scale=12; ( $lon + $dlon)" | bc -l `
    #Local - Prefer region with 30 degree square
    #reg="-R-65/65/-65/65"
    reg="-R$lonmin/$lonmax/$latmin/$latmax"
    map=$reg" -JM11"
    echo "Region:" $reg
    #echo "   Manualy define region in plot.sh file!!!."
fi

if [ $kmap -eq 4 ] ;then #Polar Stereog projection - local
    if [ $kplot -eq 1 ] ;then
	width=18c  
	height=18c 
    else
	width=26c  
	height=18c 
    fi
    pole=-90 #90
    latpol=-60 #70.5
    size=7
    # scalepos = xpos/ypos/length/width[h] 
    scalepos=15.5c/6.5c/12c/0.5c
    reg="-R-180/180/$pole/$latpol"
    map=$reg" -Js0/$pole/$size/$latpol"
    echo
    echo "Pole used (lat):" $pole
    echo 
fi

if [ $kmap -eq 5 ] ;then #Linear projection - global
    if [ $kplot -eq 1 ] ;then #Only mesh plot
	width=16c  
	height=8c 
    else #Space for scale needed
	width=29c  
	height=8c 
    fi
    # scalepos = xpos/ypos/length/width[h] 
    scalepos=22.3c/2.5c/6c/0.5c
    #Global
    #reg="-R0/360/-90/90"
    reg="-R-180/180/0/90"
    #map=$reg" -JX16/8"
    map=$reg" -Jx0.06d"    
fi


echo Map : $map

#GMT settings
gmt_shell_functions.sh
if [ "$gmt" == "$gmt5" ]; then
   $gmt gmtset PS_MEDIA Custom_${height}x${width}
   $gmt gmtset FONT_ANNOT_PRIMARY +18p
else
   $gmt gmtset PAPER_MEDIA Custom_${height}x${width}
   $gmt gmtset ANNOT_FONT_SIZE +18p
fi

#Map composition
mesh=""
case $kplot in
    1) 	mesh=$1 ; plot=$2 ;;
    2) 	scalar=$1 ; plot=$2 ;;
    3) 	vec=$1 ; plot=$2 ;;
    4) 	mesh=$1 ; scalar=$2 ; plot=$3 ;;
    5) 	mesh=$1 ; vec=$2 ; plot=$3 ;;
    6)	scalar=$1 ; vec=$2 ; plot=$3 ;;
    7) 	mesh=$1 ; scalar=$2 ; vec=$3	; plot=$4 ;;
esac

#Set naming variables

#Output graph file
if [ $plot ] ; then
    echo "Output filename:" $plot    
else
    plot=$graphdir/
fi

if [ $scalar ] ; then
    tmp=$scalar
    scalar=`basename $tmp .dat`
    dirscalar=`dirname $tmp`
    plot=$plot$scalar
fi


if [ $mesh ] ; then

    tmp=$mesh  
    mesh=`basename $tmp _nodes.gmt`
    dirmesh=`dirname $tmp `
    if [ $scalar ] ; then
	plot=$plot
	#plot=$plot"_"$mesh
    else
	plot=$plot$mesh
    fi
fi

if [ $vec ] ; then
    tmp=$vec
    vec=`basename $tmp .dat`
    dirvec=`dirname $tmp `
    plot=$plot$vec
    #echo $plot
fi

baseplot=$plot
plot=$plot.ps

echo $mesh

#Set mesh properties
if [ $mesh ] ; then
    if [ $kmap -ne 2 -a $kmap -ne 5 ] ; then
	mesh=$dirmesh/$mesh
	echo 'Mesh to be ploted:' $mesh
	echo "Select kind of mesh to be ploted:"
	echo " 1) - "
	echo " 2) Simple triangular/primal geodesic grid - no labels"
	echo " 3) Geodesic triangular/primal grid with labels"
	echo " 4) Geodesic triangular/primal grid with labels and edge vectors"
	echo " 5) Simple Geodesic voronoi/dual grid"
	echo " 6) Simple Geodesic voronoi/dual grid with triangulation/primal"
	echo " 7) Geodesic voronoi/dual grid with labels, edge vectors"
	read kmesh
        #Grid files
	nodes=$mesh'_nodes.gmt'
	ed=$mesh'_ed.gmt'
	edhx=$mesh'_edhx.gmt'
	edc=$mesh'_edc.gmt'
	ednr=$mesh'_ednr.gmt'
	edtg=$mesh'_edtg.gmt'
	edhxnr=$mesh'_edhxnr.gmt'
	edhxtg=$mesh'_edhxtg.gmt'
	trcc=$mesh'_trcc.gmt'
	#kmesh=1
    else
	echo "Global mesh plot not possible with Linear projection, choose another projection..."
	exit 0
    fi
fi


#      Start ploting 
#---------------------------------

#Set basemap
if [ "$gmt" == "$gmt5" ]; then
    grid_pen="--MAP_GRID_PEN_PRIMARY"
    base25="-BWSne -Bx40 -By40"
else
    grid_pen="--GRID_PEN_PRIMARY"
    base25="-BWSne40/40"
fi

if [ $kmap -eq 1 ] ; then   #Azimutal projection
    #Set basemap without gridlines
    if [ $kplot -eq 1 ] ; then	#Only mesh
	$gmt psbasemap  $map -B0 -K -V  -Xc -Yc > $plot
    else  #Scale will be needed
	$gmt psbasemap  $map -B0 -K -V -Yc > $plot
    fi
fi

if [ $kmap -eq 2 -o $kmap -eq 5 ] ; then #Linear projection
    #Basemap - edit -Bxgx to get lines where wanted

    #    $gmt psbasemap  $map -BWSne40/40 -Yc -K  \
    $gmt psbasemap  $map $base25 -Yc -K  \
	$grid_pen=faint,gray,-  > $plot
fi

if [ $kmap -eq 3 ] ; then  #Mercator projection, local
    $gmt psbasemap  $map -B5/5 -Yc -K  \
	$grid_pen=faint,gray,-  > $plot
fi

if [ $kmap -eq 4 ] ; then  #Polar plot
    if [ $kplot -eq 1 ] ; then	#Only mesh
	#                                 -B dlon/dlat 
	$gmt psbasemap  $map -K -V  -Xc -Yc -B20g20 > $plot
    else  #Scale will be needed
	$gmt psbasemap  $map -K -V -Yc -B20g20 > $plot
    fi
fi

#Scalar field plot
if [ $scalar ] ; then

    echo Basename for scalar: $scalar
    scalar=$dirscalar/$scalar
    
    # For endgame model, data file will end with grid resolution
    #   with 'x' separating the dimensions - search for this x and get resolution
    strindex() { 
	x="${1%%$2*}"
	[[ $x = $1 ]] && echo -1 || echo ${#x}
    }
    scal=${scalar: -5}
    pos=$( strindex "$scal" x )   # prints 4
    #echo $pos
    #echo $scal
    if [ $pos -eq -1 ] ; then	 # imodel and fem model - use 0.25
	dy=0.25      #1440x720 grid - std output of imodel and femswm
	#dy=1.40625  #256x128 grid
	#dy=0.17578125 #2048x1024 grid
    else
	len=${#scal}
	grid=${scal: $pos : $len}
	len=${#grid}
	grid=${grid: 1 : $len}
	#echo $grid
	dy=`echo "scale=12; (180.0/$grid)" | bc -l `
    fi
    echo "Grid spacing used: " $dy

    if [ "$gmt" == "$gmt5" ]; then
	xyz2grd_pars="-bif -r"
    else
	xyz2grd_pars="-bis -F"
    fi
    
    # Convert binary xyz values to grid structure -bi 
    #      (-bis means single precision, def is double)
    # Set dy above
    $gmt xyz2grd $scalar.dat -I$dy  $reg $xyz2grd_pars -G$scalar.grd -V

    #Check if field has negative values
    minval=`$gmt grdinfo $scalar.grd -C | awk '{printf "%12.8e", $6 }' `
    minval2=`echo $minval | sed 's/e/\\*10\\^/' | sed 's/+//'`
    #echo $minval $minval2
    maxval=`$gmt grdinfo $scalar.grd -C | awk '{printf "%12.8e", $7 }' `
    maxval2=`echo $maxval | sed 's/e/\\*10\\^/' | sed 's/+//'`

    #Ask if absolute value wanted
    postest=`echo "scale=12; (-0.0000000001 <= $minval2)" | bc -l `
    kabs=0
    #   if [ $postest -eq 0 ] ; then
    #	    echo "Field with negative values. Use absolute value?:"
    #	    echo " 0) NO"
    #	    echo " 1) YES"
    #	    read kabs
    #	    echo    
    #    fi

    if [ $kabs -eq 1 ] ; then
	#Set absolute value  - could set log scale
	$gmt grdmath $scalar.grd ABS = $scalar.grd
    fi

    #    #Ask if want to subtract one from values (density values)
    #    nearonetest1=`echo "scale=12; (0.9 <= $minval2)" | bc -l `
    #    nearonetest2=`echo "scale=12; (1.1 >= $maxval2)" | bc -l `
    kone=0
    #    if [ $nearonetest1 -eq 1 -a $nearonetest2 -eq 1 ] ; then
    #	    echo "Field near value one. Substract one?:"
    #	    echo " 0) NO"
    #	    echo " 1) YES"
    #	    read kone
    #	    echo    
    #    fi

    if [ $kone -eq 1 ] ; then
	#Subtract one from values
	$gmt grdmath $scalar.grd -1 ADD = $scalar.grd       
	postest=0
    fi

    #$gmt grdmath $scalar.grd $minxval SUB = $scalar.grd       

    maxgtabsmin=0
    maxgtabsmin=`echo "scale=12; (-1 * $minval2 <= $minval2)" | bc -l `
    echo max gt abs min ? $maxgtabsmin

    #$gmt grdmath $scalar.grd ABS 1 ADD LOG = $scalar.grd
    #$gmt grdmath $scalar.grd ABS EXP -1 ADD = $scalar.grd
    #$gmt grdmath $scalar.grd -1 ADD ABS +1 ADD LOG = $scalar.grd

    kg=0
    #    echo "Multiply/divide by gravity ?:"
    #    echo " 0) NO"
    #    echo " 1) Multiply (h->phi or x/phi->xi/h) "
    #    echo " 2) Divide   (phi->h or xi/h->xi/phi"
    #    read kg
    #    echo    
    
    if [ $kg -eq 1 ] ; then
	$gmt grdmath $scalar.grd 9.80616 MUL = $scalar.grd
    fi
    if [ $kg -eq 2 ] ; then
	$gmt grdmath $scalar.grd 9.80616 DIV = $scalar.grd
    fi

    #Scale marks
    minval=`$gmt grdinfo $scalar.grd -C | awk '{printf "%1.8e", $6 }' `
    maxval=`$gmt grdinfo $scalar.grd -C | awk '{printf "%1.8e", $7 }' `
    #minval=0 #-5.48408926e-02
    #maxval=1.00001  #1.15045559e+00
    #minval=-8.21844757e-01
    #maxval=7.33441532e-01

    kscale=0
    echo "Fixed scale?:"
    echo " 0) NO - Variable"
    echo " 1) YES - Values will be asked"
    read kscale
    echo    
    
    if [ $kscale -eq 1 ] ; then
	echo " Enter min scale value (actual: " $minval " )"	    
	read minval
	echo " Enter max scale value (actual: " $maxval " )"    
	read maxval
    fi

    echo "Min" $minval "Max" $maxval

    minval2=`echo $minval | sed 's/e/\\*10\\^/' | sed 's/+//'`
    maxval2=`echo $maxval | sed 's/e/\\*10\\^/' | sed 's/+//'`
    #maxval2=`echo "scale=10;  ($maxval2 + $maxval2)" | bc -l `
    #maxval2=`echo "scale=10;  ($maxval2 + 0.0000001)" | bc -l `
    if [ $kmap -eq 5 ] ;then #Linear projection - regional
	scale=`echo "scale=12; ( $maxval2 - $minval2 )/7" | bc -l `
	scale=`echo "scale=16; ( $maxval2 - $minval2 )/7" | bc -l `
    else
	scale=`echo "scale=12; ( $maxval2 - $minval2 )/10" | bc -l `
	scale=`echo "scale=16; ( $maxval2 - $minval2 )/10" | bc -l `
    fi

    echo "Min" $minval2 "Max" $maxval2
    echo "Scale:" $scale
    #echo "Scale OPT:" $scaleopt
    
    #Tiks for the scale
    
    #Create color pallet
    #------------------------------------
    
    kcolor=1

    echo "Select color pallet:"
    echo " 1) White-DarkBlue"
    echo " 2) Blue-White-Red"
    echo " 3) Multicolor jet"
    echo " 4) White-Light Blue"
    echo " 5) Seismology - Multicolor"
    echo " 6) Color - Seis + 0-white to 1-black"
    echo " 7) Color - Inverted Seis + 1-white to 0-black"
    echo " 8) Grayscale "
    echo " 9) Blue-White-Red smoother"
    read kcolor
    echo    

    if [ $kcolor -eq 1 ] ; then  #White-Darkblue
        $gmt makecpt -Cdarkblue -I -T$minval/$maxval/$scale -Z \
    	    --COLOR_BACKGROUND=blue  --COLOR_FOREGROUND=white > $scalar.cpt
    fi

    if [ $kcolor -eq 2 ] ; then #Blue - Red
	if [ $postest -eq 0 ] ; then
	    #grd2cpt $scalar.grd -Cpolar  -Z -T= \
	    $gmt grd2cpt $scalar.grd -Cbluered  -Z -T=  > $scalar.cpt
	else            
            $gmt makecpt -Cpolar -T$minval/$maxval/$scale -Z \
    		--COLOR_BACKGROUND=blue  --COLOR_FOREGROUND=red > $scalar.cpt
	fi
    fi

    if [ $kcolor -eq 9 ] ; then #Blue - Red
	if [ $postest -eq 0 ] ; then
	    $gmt grd2cpt $scalar.grd -Cbluered_zero  -Z -T=  > $scalar.cpt
	else            
            $gmt makecpt -Cpolar -T$minval/$maxval/$scale -Z \
    		--COLOR_BACKGROUND=blue  --COLOR_FOREGROUND=red > $scalar.cpt
	fi
    fi

    if [ $kcolor -eq 3 ] ; then #Multicolor jet
	$gmt makecpt -Cjet -T$minval/$maxval/$scale -Z  \
    	    --COLOR_BACKGROUND=0/0/127  --COLOR_FOREGROUND=127/0/0 > $scalar.cpt
    fi
    
    if [ $kcolor -eq 4 ] ; then #Light Blue -White
	$gmt makecpt -Cblue -T$minval/$maxval/$scale -Z  > $scalar.cpt
    fi

    if [ $kcolor -eq 5 ] ; then #Multicolor jet
	$gmt makecpt -Cmulticol  -T$minval/$maxval/$scale -Z   > $scalar.cpt
    fi    

    if [ $kcolor -eq 6 ] ; then #Multicolor
	#makecpt -Crainbow  -T$minval/$maxval/$scale -Z   > $scalar.cpt
	#makecpt -Chaxby -I -T$minval/$maxval/$scale -Z   > $scalar.cpt
	$gmt makecpt -Cseism_nonlin -I -T$minval/$maxval/$scale -Z   > $scalar.cpt
	#makecpt -Cwysiwyg -I -T$minval/$maxval/$scale -Z   > $scalar.cpt
    	#--COLOR_BACKGROUND=0/0/127  --COLOR_FOREGROUND=127/0/0
    fi    


    if [ $kcolor -eq 7 ] ; then #Multicolor
	#makecpt -Cjet  -T$minval/$maxval/$scale -Z   > $scalar.cpt
	#makecpt -Chaxby -I -T$minval/$maxval/$scale -Z   > $scalar.cpt
	#makecpt -Cseism -I -T$minval/$maxval/$scale -Z   > $scalar.cpt
	$gmt makecpt -Cseism_nonlin -T$minval/$maxval/$scale -Z   > $scalar.cpt
	#makecpt -Cwysiwyg -I -T$minval/$maxval/$scale -Z   > $scalar.cpt
    	#--COLOR_BACKGROUND=0/0/127  --COLOR_FOREGROUND=127/0/0
    fi
    
    if [ $kcolor -eq 8 ] ; then #Greyscale
	#makecpt -Cjet  -T$minval/$maxval/$scale -Z   > $scalar.cpt
	#makecpt -Chaxby -I -T$minval/$maxval/$scale -Z   > $scalar.cpt
	#makecpt -Cseism -I -T$minval/$maxval/$scale -Z   > $scalar.cpt
	$gmt grd2cpt $scalar.grd -Cgray  -Z   > $scalar.cpt
	#$gmt makecpt -Cseism_nonlin -T$minval/$maxval/$scale -Z   > $scalar.cpt
	#makecpt -Cwysiwyg -I -T$minval/$maxval/$scale -Z   > $scalar.cpt
    	#--COLOR_BACKGROUND=0/0/127  --COLOR_FOREGROUND=127/0/0
    fi    


    #grd2cpt $scalar.grd -Crainbow  -Z -V  > $scalar.cpt
    #grd2cpt $scalar.grd -Chot  -Z -V  -I > $scalar.cpt

    #-------------------------------------------------

    #Plot scalar field map
    $gmt grdimage $scalar.grd -C$scalar.cpt  $map  -O -K -V  >> $plot

    if [ "$gmt" == "$gmt5" ]; then
	format_scale='--FORMAT_FLOAT_MAP=%.1e --FONT_ANNOT_PRIMARY=20'
    else
	format_scale='--D_FORMAT=%0.1e --ANOT_FONT_SIZE=20 '
    fi
    
    #Set tickmarks distance for scale
    $gmt psscale -C$scalar.cpt -D$scalepos -B$scale -O -K -V  \
	$format_scale >> $plot   

    #Print coast
    #$gmt pscoast  $map -Wfaint,150 -K -O >> $plot

    rm -rf $scalar.grd
    rm -rf $scalar.cpt
fi

# PLOTING MESH
if [ $mesh ] ; then
    #Print coast
    #$gmt pscoast  $map -Wfaint,gray -K -O >> $plot

    if [ "$gmt" == "$gmt5" ]; then
	multiflag=""
    else
	multiflag="-m"
    fi
    
    if [ $kmesh -ne 5 ] ; then
	#Print nodes
	$gmt psxy $nodes  $map -Sc0.005c  -W0.5,gray -Ggray -O  -K -V  >> $plot
	#$gmt psxy $nodes  $map -Sc0.5c  -W0.5/blue -Ggray -O  -K -V  >> $plot

	#Print triangle edges
	#GMT5 cannot cope with comments in the same as gmt 4, so remove them
	#SED -i.bak '/>/d' $ed
	#sed -i.bak '/#/d' $ed
	#echo $map
	if [ $scalar ] ; then
	    $gmt psxy $ed $map  -Wthin,gray $multiflag  -O -K -V  >> $plot
	else
	    $gmt psxy $ed $map -Wthin,blue $multiflag  -O -K -V  >> $plot
	fi
    fi

    #Labels
    if [ $kmesh -eq 3 -o $kmesh -eq 4 -o $kmesh -eq 7 ] ; then
	echo "Labeling mesh ..."
	#Save node label file
	#  lon, lat, label
	if [ "$gmt" == "$gmt5" ]; then
	    awk '{print $1, $2, NR}' $nodes > 'nodes_ll.txt'
	else
	    awk '{print $1, $2, 10, 0, 0, "LT", NR}' $nodes > 'nodes_ll.txt'
	fi
	
	#Plot nodes's labels
	if [ "$gmt" == "$gmt5" ]; then
	    font="-F+f8p,Helvetica,black"
	else
	    font=""
	fi
	$gmt pstext 'nodes_ll.txt' $font  $map -O -K >> $plot
	
	#Triangle circumcenters
	$gmt psxy $trcc  $map -St0.05c  -W0.05,red -Gred  -O -K >> $plot
	
	#Save triangle labels
	if [ "$gmt" == "$gmt5" ]; then
	    awk '{print $1, $2, NR}' $trcc > 'triangle_cc.txt'
	else
	    awk '{print $1, $2, 5, 0, 0, "LT", NR}' $trcc > 'triangle_cc.txt'
	fi
	
	#Plot triangle labels
	if [ "$gmt" == "$gmt5" ]; then
	    font="-F+f6p,Helvetica,black"
	else
	    font=""
	fi
	$gmt pstext 'triangle_cc.txt' $font $map -O  -K >> $plot	
	rm 'nodes_ll.txt'  'triangle_cc.txt'
    fi
    
    #Normal and tangent edge vectors
    if [ $kmesh -eq 4  ] ; then
	echo "Ploting mesh normal and tangent vectors ..."

        # Plot normal and tangent vectors
	$gmt psxy $ed $map   -W1,blue   -O -K  >> $plot

	$gmt psxy $edc  $map -Sc0.0001c  -W0.01,blue -Ggreen  -O  -K >> $plot

        #  lon, lat, label 
	awk '{print $1, $2, NR-1}' $edc > 'edc_ll.txt'

	$gmt pstext 'edc_ll.txt' -F+f6p,Helvetica,black $map -O -K  >> $plot

	#vecstyle="0.01/0.05/0.05" #depreciated
	vecstyle="0.01i+bc+ea+g"
	# GMT 4
	#vecstyle="0.01/0.05/0.05"
	$gmt psxy $ednr $map -W0.5,green -Ggreen -SV$vecstyle -O -K  >> $plot
	$gmt psxy $edtg $map -W0.5,green -Ggreen -SV$vecstyle -O -K  >> $plot
	
        #Clean workspace
	rm  'edc_ll.txt'	    
    fi
    #Voronoi mesh, edges only 
    if [ $kmesh -ge 5  ] ; then
	if [ $scalar ] ; then
	    #$gmt psxy $edhx  $map   -Wthin,gray   -O -K   >> $plot
	    $gmt psxy $edhx  $map  $multiflag  -W0.5,black   -O -K   >> $plot
	else
	    #psxy $edhx -m  $map   -Wthin/black   -O -K   >> $plot
	    $gmt psxy $edhx  $map   $multiflag -W0.1,black   -O -K   >> $plot
	fi

    fi
    # Voronoi mesh with labels and normal/tangent edge vectors
    if [ $kmesh -eq 7  ] ; then
	# arrowwidth/headlength/headwidth  
	#vecstyle="0.008/0.03/0.03"
	#vecstyle="0.05/0.06/0.08" #depreciated in gmt5
	vecstyle="0.01i+bc+ea+g"
	
        # Set vectorcolor
	veccolor=0/100/0

	# -Wx is pen, x is brush size, -G is fill 
	$gmt psxy $edhxnr $map -W0.5,$veccolor -G$veccolor -SV$vecstyle -O -K  >> $plot
	$gmt psxy $edhxtg $map -W0.5,$veccolor -G$veccolor -SV$vecstyle -O -K  >> $plot

        #  lon, lat, font size, ? , ? , position, label 
	awk '{print $1, $2, NR-1}' $edc > 'edc_ll.txt'

	$gmt pstext 'edc_ll.txt' -F+f6p,Helvetica,black $map -O -K  >> $plot
    fi
fi

#Vector field plot
if [ $vec ] ; then
    echo "NOT converted to GMT5 and new GMT4"
    exit 0
    
    #pscoast  $map -Wfaint/150/150/150 -K -O >> $plot

    echo Basename for vectors: $vec
    vec=$dirvec/$vec

    #Use this AWK command to scale vectors
    # lon, lat, angle from north, length
    #awk '{ if ( NR % 18 == 0 ) { print $1, $2, $3, $4/200 ;} }' $vec.dat > 'vec.tmp'
    #awk '{ if ( NR % 1 == 0 ) { print $1, $2, $3, $4/200 ;} }' $vec.dat > 'vec.tmp'

    #Map for lon-lat vector grids
    #awk '{ if ( $1 % 3.125 == 0 && $2 % 3.125 == 0 ) { print $1, $2, $3, $4/200 ;} }' $vec.dat > 'vec.tmp'
    #awk '{ if ( NR < 1630 ) { print $1, $2, $3, $4/100 ;} }' $vec.dat > 'vec.tmp'
    #awk '{ { print $1 , $2, $3, $4*10000000 ;} }' $vec.dat > 'vec.tmp'
    awk '{ { print $1, $2, $3, $4/200 ;} }' $vec.dat > 'vec.tmp'

    #Draw vectors
    # arrowwidth/headlength/headwidth  
    vecstyle="0.01/0.03/0.05"

    # Set vector color and size
    veccolor=0/100/0

    # -W is pen, -G is fill
    $gmt psxy 'vec.tmp' $map -W1/$veccolor -G$veccolor -SV$vecstyle -O -K -V\
	--PAPER_MEDIA=Custom_${height}x${width} >> $plot

    #rm -rf 'vec.tmp'
fi

# Point to close the graph file
#echo "30 -65" > point1.dat
#echo "30 -64" > point2.dat
#echo "31 -65" > point3.dat
#echo "31 -64" > point4.dat
echo "-180 -90" > 'point.dat'
#echo "-23.8744234108428        34.7927455551923" > 'point.dat'
$gmt psxy  point.dat  $map -Sx0.0005c  -W0.0005 -Gred -O  >> $plot
#psxy  point.dat  $map -Sx0.1c  -W0.5/black -Gred -O -K >> $plot
#psxy point.dat  $map -Sx0.01c  -W0.05/blue -Gblue -O  >> $plot
#psxy point1.dat  $map -Sx0.1c  -W0.5/red -Gred -O  -V -K >> $plot
#psxy point2.dat  $map -Sx0.1c  -W0.5/red -Gred -O -K -V >> $plot
#psxy point3.dat  $map -Sx0.1c  -W0.5/red -Gred -O -K -V >> $plot
#psxy point4.dat  $map -Sx0.1c  -W0.5/red -Gred -O -K  -V >> $plot
#psxy mapcenter.dat  $map -Sx0.1c  -W0.5/green -Ggreen -O >> $plot
echo
echo "Plotted file:" $plot

newtitle=`basename $plot .ps`
sed -i s/"^%%Title:.*"/"%%Title:$newtitle"/ $plot

# Convert file to EPS and PDF and generate grayscale figures for papers
#-------------------------------------------------------------------------

if [ "$gmt" == "$gmt5" ]; then
    $gmt psconvert -Te -A -P $plot  # EPS
    $gmt psconvert -Tf -A -P $plot  # PDF
else
    ps2eps -R + -f $plot
    epstopdf $baseplot".eps"
fi
echo "and pdf and eps files (for pdflatex and latex)."

#Grayscale figures
cp $plot $baseplot"_gray".ps
gs -dNOPAUSE -dBATCH -dSAFER -sDEVICE=pdfwrite -dProcessColorModel=/DeviceGray -dColorConversionStrategy=/Gray -dPDFUseOldCMS=false  -dNOCACHE -o $baseplot"_gray.pdf" -f $baseplot.pdf
gs -dNOPAUSE -dBATCH -dSAFER -sDEVICE=pdfwrite -dProcessColorModel=/DeviceGray -dColorConversionStrategy=/Gray -dPDFUseOldCMS=false  -dNOCACHE -o $baseplot"_gray.ps" -f $baseplot.ps
pdf2ps -eps $baseplot"_gray.ps" $baseplot"_gray.eps"


#rm $plot
#rm $baseplot"_gray.ps"
gv $baseplot.ps &
#okular $baseplot.pdf &

echo
echo "-----------------------------------------------------------"
