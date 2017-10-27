#!/bin/bash

echo '============================'
echo 'Plotting script for GMT'
echo 'P. Peixoto - Aug 2011       '
echo '============================'
echo

killall gv


#Directory for graphs
graphdir=../graphs


echo 

#Arguments check
if [ $# -eq 0  ] ; then
    echo "Please enter arguments: 2 scalar fields"
    exit 0
fi  

#MAP PROJECTION
echo "Select map projection:"
echo " 1) Lambert - AZIMUTHAL (-JA) - Global"
echo " 2) Linear (-JX) - Global"  #Mercator - CYLINDRICAL (-JM)"
echo " 3) Mercator - CYLINDRICAL (-JM) - Regional only"
echo " 4) Polar Stereographic (-Js) - Regional only"
read kmap
echo


#MAP CENTER
#lon=0 #-179.5 #1.23296 #31.7174
#lat=0 #-80 #-89.5 #-61.30099 #-0.0374
lon=0 #17.5
lat=65
#lon=30.125
#lat=-64.875
#lon=90.5 #1.23296 #31.7174
#lat=-0.5 #-89.5 #-61.30099 #-0.0374

echo $lon $lat > mapcenter.dat

#MAP REGION
if [ $kmap -eq 1 ] ;then #Azimutal projection - global
    if [ $kplot -eq 1 ] ;then
	width=15c  
	height=15c 
    else
	width=22c  
	height=15c 
    fi
    # scalepos = xpos/ypos/length/width[h] 
    scalepos=14.5c/6.5c/12c/0.5c
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
        #map=$reg" -JX16/8"
    map=$reg" -Jx0.06d"    
fi

if [ $kmap -eq 3 ] ;then  #Local Mercator projection
    if [ $kplot -eq 1 ] ;then #Only mesh plot
	width=15c  
	height=15c 
    else      #Space for scale needed
	width=19c  
	height=15c 
    fi
     # scalepos = xpos/ypos/length/width[h] 
    scalepos=12.5c/5.5c/11c/0.5c

    #Range of map
    dlat=10
    dlon=50
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
    pole=-90
    latpol=-89.5
    size=7
    # scalepos = xpos/ypos/length/width[h] 
    scalepos=15.5c/6.5c/12c/0.5c
    reg="-R-180/180/$pole/$latpol"
    map=$reg" -Js0/$pole/$size/$latpol"
    echo
    echo "Pole used (lat):" $pole
    echo 
fi

echo Map : $map

#GMT settings
. gmt_shell_functions.sh
gmtset PAPER_MEDIA Custom_${height}x${width}
gmtset ANNOT_FONT_SIZE +18p

#Map composition
scalar1=$1 
scalar2=$2

#Set naming variables


tmp1=$scalar1
tmp2=$scalar2
scalar1=`basename $tmp1 .dat`
scalar2=`basename $tmp2 .dat`
dirscalar=`dirname $tmp1`
plot=$plot$scalar1$scalar2
scalar=$scalar1$scalar2

baseplot=$plot
plot=$plot.ps

echo $mesh


#      Start ploting 
#---------------------------------

#Set basemap
if [ $kmap -eq 1 ] ; then   #Azimutal projection
   #Set basemap without gridlines
    if [ $kplot -eq 1 ] ; then	#Only mesh
	psbasemap  $map -B0 -K -V  -Xc -Yc > $plot
    else  #Scale will be needed
	psbasemap  $map -B0 -K -V -Yc > $plot
    fi
fi

if [ $kmap -eq 2 ] ; then #Linear projection
    #Basemap - edit -Bxgx to get lines where wanted

    psbasemap  $map -BWSne40/40 -Yc -K  \
	--GRID_PEN_PRIMARY=faint,gray,-  > $plot
fi

if [ $kmap -eq 3 ] ; then  #Mercator projection, local
    psbasemap  $map -B5/5 -Yc -K  \
	--GRID_PEN_PRIMARY=faint,gray,-  > $plot
fi

if [ $kmap -eq 4 ] ; then  #Polar plot
    if [ $kplot -eq 1 ] ; then	#Only mesh
	psbasemap  $map -K -V  -Xc -Yc -B20g20/0.1g0.1 > $plot
    else  #Scale will be needed
	psbasemap  $map -K -V -Yc -B20g20 > $plot
    fi
fi


echo Basename for scalar1: $scalar1
scalar1=$dirscalar/$scalar1
echo Basename for scalar2: $scalar2
scalar2=$dirscalar/$scalar2
    
    # Convert binary xyz values to grid structure -bi
    #xyz2grd $scalar.dat -I0.125  $reg -bi -F -G$scalar.grd -V
    xyz2grd $scalar1.dat -I0.25  $reg -bi -F -G$scalar1.grd -V
    xyz2grd $scalar2.dat -I0.25  $reg -bi -F -G$scalar2.grd -V
    #xyz2grd $scalar.dat -I1  $reg -bi -F -G$scalar.grd -V    

    grdmath $scalar1.grd $scalar2.grd SUB = $scalar.grd
 
    #Check if field has negative values
    minval=`grdinfo $scalar.grd -C | awk '{printf "%12.8e", $6 }' `
    minval2=`echo $minval | sed 's/e/\\*10\\^/' | sed 's/+//'`
    #echo $minval $minval2
    maxval=`grdinfo $scalar.grd -C | awk '{printf "%12.8e", $7 }' `
    maxval2=`echo $maxval | sed 's/e/\\*10\\^/' | sed 's/+//'`

    #Ask if absolute value wanted
    postest=`echo "scale=12; (-0.0000000001 <= $minval2)" | bc -l `
    kabs=0
    if [ $postest -eq 0 ] ; then
	    echo "Field with negative values. Use absolute value?:"
	    echo " 0) NO"
	    echo " 1) YES"
	    read kabs
	    echo    
    fi

    if [ $kabs -eq 1 ] ; then
       #Set absolute value and log scale
	grdmath $scalar.grd ABS = $scalar.grd
    fi

 
 
    maxgtabsmin=0
    maxgtabsmin=`echo "scale=12; (-1 * $minval2 <= $minval2)" | bc -l `
    echo max gt abs min ? $maxgtabsmin

    #grdmath $scalar.grd ABS 1 ADD LOG = $scalar.grd
    #grdmath $scalar.grd ABS EXP -1 ADD = $scalar.grd
    #grdmath $scalar.grd -1 ADD ABS +1 ADD LOG = $scalar.grd

    #Scale marks
    minval=`grdinfo $scalar.grd -C | awk '{printf "%1.8e", $6 }' `
    maxval=`grdinfo $scalar.grd -C | awk '{printf "%1.8e", $7 }' `
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
    scale=`echo "scale=12; ( $maxval2 - $minval2 )/10" | bc -l `
    scale=`echo "scale=16; ( $maxval2 - $minval2 )/10" | bc -l `
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
    echo " 7) Color "

    read kcolor
    echo    

    if [ $kcolor -eq 1 ] ; then  #White-Darkblue
        makecpt -Cdarkblue -I -T$minval/$maxval/$scale -Z \
    	    --COLOR_BACKGROUND=blue  --COLOR_FOREGROUND=white > $scalar.cpt
    fi

    if [ $kcolor -eq 2 ] ; then #Blue - Red
	if [ $postest -eq 0 ] ; then
	    #grd2cpt $scalar.grd -Cpolar  -Z -T= \
	    grd2cpt $scalar.grd -Cbluered  -Z -T=  > $scalar.cpt
	else            
            makecpt -Cpolar -T$minval/$maxval/$scale -Z \
    		--COLOR_BACKGROUND=blue  --COLOR_FOREGROUND=red > $scalar.cpt
	fi
    fi
    
    if [ $kcolor -eq 3 ] ; then #Multicolor jet
	makecpt -Cjet -T$minval/$maxval/$scale -Z  \
    		--COLOR_BACKGROUND=0/0/127  --COLOR_FOREGROUND=127/0/0 > $scalar.cpt
        #grd2cpt $scalar.grd -Cjet  -Z -V  > $scalar.cpt
    fi
    
    if [ $kcolor -eq 4 ] ; then #Light Blue -White
	makecpt -Cblue -T$minval/$maxval/$scale -Z  > $scalar.cpt
    fi

    if [ $kcolor -eq 5 ] ; then #Multicolor jet
	#makecpt -Cseis -I -T$minval/$maxval/$scale -Z   > $scalar.cpt
	makecpt -Cmulticol  -T$minval/$maxval/$scale -Z   > $scalar.cpt
    	#--COLOR_BACKGROUND=0/0/127  --COLOR_FOREGROUND=127/0/0
    fi    

    if [ $kcolor -eq 6 ] ; then #Multicolor
	#makecpt -Crainbow  -T$minval/$maxval/$scale -Z   > $scalar.cpt
	#makecpt -Chaxby -I -T$minval/$maxval/$scale -Z   > $scalar.cpt
	makecpt -Cseism_nonlin -I -T$minval/$maxval/$scale -Z   > $scalar.cpt
	#makecpt -Cwysiwyg -I -T$minval/$maxval/$scale -Z   > $scalar.cpt
    	#--COLOR_BACKGROUND=0/0/127  --COLOR_FOREGROUND=127/0/0
    fi    

    if [ $kcolor -eq 7 ] ; then #Multicolor
	#makecpt -Cjet  -T$minval/$maxval/$scale -Z   > $scalar.cpt
	#makecpt -Chaxby -I -T$minval/$maxval/$scale -Z   > $scalar.cpt
	makecpt -Cseism -I -T$minval/$maxval/$scale -Z   > $scalar.cpt
	#makecpt -Cwysiwyg -I -T$minval/$maxval/$scale -Z   > $scalar.cpt
    	#--COLOR_BACKGROUND=0/0/127  --COLOR_FOREGROUND=127/0/0
    fi    

    #grd2cpt $scalar.grd -Crainbow  -Z -V  > $scalar.cpt
    #grd2cpt $scalar.grd -Chot  -Z -V  -I > $scalar.cpt

    #-------------------------------------------------

    #Plot scalar field map
    grdimage $scalar.grd -C$scalar.cpt  $map  -O -K -V  >> $plot

    #Set tickmarks distance for scale
    psscale -C$scalar.cpt -D$scalepos -B$scale -O -K -V  \
    --D_FORMAT="%0.1e" --ANOT_FONT_SIZE="20">> $plot
    #psscale -C$scalar.cpt -D$scalepos -O -K -V  >> $plot

    #Print coast
    #pscoast  $map -Wfaint/150/150/150 -K -O  >> $plot

    #rm -rf $scalar.grd




#echo "30 -65" > point1.dat
#echo "30 -64" > point2.dat
#echo "31 -65" > point3.dat
#echo "31 -64" > point4.dat
echo "-180 -90" > 'point.dat'
#echo "-23.8744234108428        34.7927455551923" > 'point.dat'
psxy  point.dat  $map -Sx0.0005c  -W0.0005/black -Gred -O  >> $plot
#psxy  point.dat  $map -Sx0.1c  -W0.5/black -Gred -O -K >> $plot

#Plot mapcenter - Just to end the ps file

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

echo "and pdf and eps files (for pdflatex and latex)."
ps2raster -Te -A -P $plot 
ps2raster -Tf -A -P $plot 
#ps2pdf14 -dEmbedAllFonts=true -dPDFSETTINGS=/prepress $plot $baseplot-tmp.pdf
#pdfcrop $baseplot-tmp.pdf $baseplot.pdf
rm $plot
gv $baseplot.eps &
#okular $baseplot.pdf &

echo
echo "-----------------------------------------------------------"