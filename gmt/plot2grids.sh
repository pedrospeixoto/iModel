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
    echo "Please enter arguments:"
    echo " Arguments (files) should be passed in the following order:"
    echo " 2 Mesh (nodes file .gmt)"
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
lon=0
lat=0
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
    dlat=20
    dlon=40
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
mesh=$1
mesh2=$2

#Set naming variables

#Output graph file
if [ $plot ] ; then
    echo "Output filename:" $plot    
else
    plot=$graphdir/
fi

if [ $mesh ] ; then
    tmp=$mesh  
    mesh=`basename $tmp _nodes.gmt`
    tmp2=$mesh2
    mesh2=`basename $tmp2 _nodes.gmt`
    dirmesh=`dirname $tmp `
	plot=$plot$mesh$mesh2
fi

baseplot=$plot
plot=$plot.ps

echo $mesh
echo $mesh2

#Set mesh properties
if [ $mesh ] ; then
    if [ $kmap -ne 2 ] ; then
	mesh=$dirmesh/$mesh
	echo 'Mesh to be ploted:' $mesh
	mesh2=$dirmesh/$mesh2
	echo 'Second mesh to be ploted:' $mesh2
	echo
	echo "Select kind of mesh to be ploted:"
	echo " 1) - "
	echo " 2) Simple triangular geodesic grid - no labels"
	echo " 3) Geodesic triangular grid with labels"
	echo " 4) Geodesic triangular grid with labels and edge vectors"
	echo " 5) Simple Geodesic voronoi grid"
	echo " 6) Simple Geodesic voronoi grid with triangulation"
	echo " 7) Geodesic voronoi grid with labels, edge vectors"
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

	nodes2=$mesh2'_nodes.gmt'
	ed2=$mesh2'_ed.gmt'
	edhx2=$mesh2'_edhx.gmt'
	edc2=$mesh2'_edc.gmt'
	ednr2=$mesh2'_ednr.gmt'
	edtg2=$mesh2'_edtg.gmt'
	edhxnr2=$mesh2'_edhxnr.gmt'
	edhxtg2=$mesh2'_edhxtg.gmt'
	trcc2=$mesh2'_trcc.gmt'
	#kmesh=1
    else
	echo "Global mesh plot not possible with Linear projection, choose another projection..."
	exit 0
    fi
fi


#      Start ploting 
#---------------------------------

#Set basemap
if [ $kmap -eq 1 ] ; then   #Azimutal projection
   #Set basemap without gridlines
	psbasemap  $map -B0 -K -V  -Xc -Yc > $plot
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
	psbasemap  $map -K -V  -Xc -Yc -B20g20/0.1g0.1 > $plot
fi

# PLOTING MESH
if [ $mesh ] ; then
	#Print coast
    #pscoast  $map -Wfaint/150/150/150 -K -O >> $plot

    if [ $kmesh -ne 5 ] ; then
	#Print nodes
	psxy $nodes  $map -Sc0.005c  -W0.5/blue -Gred -O  -K -V  >> $plot
	psxy $nodes2  $map -Sc0.005c  -W0.5/red -Gblue -O  -K -V  >> $plot

	#Print triangle edges
	psxy $ed -m  $map   -Wthin/blue   -O -K -V  >> $plot
	psxy $ed2 -m  $map   -Wthin/red   -O -K -V  >> $plot
    fi

	#Labels
    if [ $kmesh -eq 3 -o $kmesh -eq 4 -o $kmesh -eq 7 ] ; then
	echo "Labeling mesh ..."
	    #Save node label file
	    #  lon, lat, font size, ? , ? , position, label 
	awk '{print $1, $2, 6, 0, 0, "LT", NR}' $nodes > 'nodes_ll.txt'

	    #Plot nodes's labels
	 pstext 'nodes_ll.txt' $map -O -K >> $plot
	
	    #Triangle circumcenters
	 psxy $trcc  $map -St0.05c  -W0.05/red -Gred  -O -K >> $plot
	
	    #Save triangle labels
	awk '{print $1, $2, 5, 0, 0, "LT", NR}' $trcc > 'triangle_cc.txt'
	
	    #Plot triangle labels
	 pstext 'triangle_cc.txt' $map -O  -K >> $plot
	
	rm 'nodes_ll.txt'  'triangle_cc.txt'
    fi
    
	#Normal and tangent edge vectors
    if [ $kmesh -eq 4  ] ; then
	echo "Ploting mesh normal and tangent vectors ..."

             # Plot normal and tangent vectors
	psxy $ed -m  $map   -W1/blue   -O -K  >> $plot
	psxy $ed2 -m  $map   -W1/red   -O -K  >> $plot

	psxy $edc  $map -Sc0.1c  -W0.1/blue -Ggreen  -O  -K >> $plot
	psxy $edc  $map -Sc0.1c  -W0.1/red -Ggreen  -O  -K >> $plot

            #  lon, lat, font size, ? , ? , position, label 
	awk '{print $1, $2, 5, 0, 0, "LT", NR-1}' $edc > 'edc_ll.txt'

	pstext 'edc_ll.txt' $map -O -K  >> $plot

	vecstyle="0.01/0.05/0.05"
	psxy $ednr $map -W1/green -Ggreen -SV$vecstyle -O -K  >> $plot
	
	psxy $edtg $map -W1/green -Ggreen -SV$vecstyle -O -K  >> $plot
	
              #Clean workspace
	rm  'edc_ll.txt'	    
    fi
	#Voronoi mesh, edges only 
    if [ $kmesh -ge 5  ] ; then

	    #psxy $edhx -m  $map   -Wthin/black   -O -K   >> $plot
	psxy $edhx -m  $map   -W2/black   -O -K   >> $plot
	psxy $edhx2 -m  $map   -W2/blue   -O -K   >> $plot
    fi
	# Voronoi mesh with labels and normal/tangent edge vectors
    if [ $kmesh -eq 7  ] ; then
	    # arrowwidth/headlength/headwidth  
	#vecstyle="0.008/0.03/0.03"
	vecstyle="0.05/0.06/0.08"

            # Set vectorcolor
	veccolor=0/100/0

	    # -Wx is pen, x is brush size, -G is fill 
	psxy $edhxnr $map -W0.5/$veccolor -G$veccolor -SV$vecstyle -O -K  >> $plot
	psxy $edhxtg $map -W0.5/$veccolor -G$veccolor -SV$vecstyle -O -K  >> $plot

         #  lon, lat, font size, ? , ? , position, label 
	awk '{print $1, $2, 5, 0, 0, "LT", NR-1}' $edc > 'edc_ll.txt'

	pstext 'edc_ll.txt' $map -O -K  >> $plot
    fi
fi


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