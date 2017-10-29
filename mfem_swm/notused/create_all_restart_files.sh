#! /bin/bash

make
make sfc_optimize

MIN=1
MAX=5

test -n "$1" && MIN=$1
test -n "$2" && MAX=$2

CMIN=$((MIN))
CMAX=$((MAX))

MIN=$((MIN+1))
MAX=$((MAX+1))

#
# BASIC GRIDS
#
for i in `seq $CMIN $CMAX`; do
	echo "Cubed sphere level $i"
	./gengrid_cube $i > /dev/null
done

for i in `seq $MIN $MAX`; do
	echo "Hexagonal grid level $i"
	./gengrid_hex $i > /dev/null
done


#
# SFC OPTIMIZATIONS
#
for i in `seq $CMIN $CMAX`; do
	#
	# CUBED SPHERE
	#
	NX=$((3*(2**(i-1))))
	CELLS=$((6*NX*NX))
	CELLS=`printf "%010d" $CELLS`

	IDENTIFIER="cube_$CELLS"
	OUTPUTFILE="gridmap_""$IDENTIFIER""_sfc.dat"

	echo "SFC optimization for cubed sphere level $i"
	./sfc_optimize cube $CELLS 1 2 2 $OUTPUTFILE > /dev/null
done

for i in `seq $MIN $MAX`; do
	#
	# HEX GRIDS
	#
	CELLS=$((5*2**(2*i-1)+2))
	CELLS=`printf "%010d" $CELLS`

	IDENTIFIER="hex_$CELLS"
	OUTPUTFILE="gridmap_""$IDENTIFIER""_sfc.dat"

	echo "SFC optimization for hex sphere level $i"
	./sfc_optimize hex $CELLS 1 2 2 $OUTPUTFILE > /dev/null
done



#
# FEM stencils and values
#
for i in `seq $CMIN $CMAX`; do
	NX=$((3*(2**(i-1))))
	CELLS=$((6*NX*NX))
	CELLS=`printf "%010d" $CELLS`

	IDENTIFIER="cube_$CELLS"

	echo "Build operators for $IDENTIFIER"
	INPUTFILE="gridmap_""$IDENTIFIER"".dat"
	OUTPUTFILE="gridopermap_""$IDENTIFIER"".dat"

	./buildop_fem cube $CELLS $INPUTFILE $OUTPUTFILE > /dev/null


	echo "Build operators for $IDENTIFIER SFC"
	INPUTFILE="gridmap_""$IDENTIFIER""_sfc.dat"
	OUTPUTFILE="gridopermap_""$IDENTIFIER""_sfc.dat"

	./buildop_fem cube $CELLS $INPUTFILE $OUTPUTFILE > /dev/null
done
for i in `seq $MIN $MAX`; do
	CELLS=$((5*2**(2*i-1)+2))
	CELLS=`printf "%010d" $CELLS`

	IDENTIFIER="hex_$CELLS"

	echo "Build operators for $IDENTIFIER"
	INPUTFILE="gridmap_""$IDENTIFIER"".dat"
	OUTPUTFILE="gridopermap_""$IDENTIFIER"".dat"

	./buildop_fem hex $CELLS $INPUTFILE $OUTPUTFILE > /dev/null


	echo "Build operators for $IDENTIFIER SFC"
	INPUTFILE="gridmap_""$IDENTIFIER""_sfc.dat"
	OUTPUTFILE="gridopermap_""$IDENTIFIER""_sfc.dat"

	./buildop_fem hex $CELLS $INPUTFILE $OUTPUTFILE > /dev/null
done

for i in `seq $CMIN $CMAX`; do
	NX=$((3*(2**(i-1))))
	CELLS=$((6*NX*NX))
	CELLS=`printf "%010d" $CELLS`

	IDENTIFIER="cube_$CELLS"

	echo "simulation of "$IDENTIFIER""

	./femswe swenml/swenml_"$IDENTIFIER".in > /dev/null
	mv run000001_restart_00000020.dat restart_files/run000001_restart_00000020_"$IDENTIFIER".dat

	echo "simulation of "$IDENTIFIER" SFC"

	./femswe swenml/swenml_"$IDENTIFIER"_sfc.in > /dev/null
	mv run000001_restart_00000020.dat restart_files/run000001_restart_00000020_"$IDENTIFIER"_sfc.dat
done

for i in `seq $MIN $MAX`; do
	CELLS=$((5*2**(2*i-1)+2))
	CELLS=`printf "%010d" $CELLS`

	IDENTIFIER="hex_$CELLS"

	echo "simulation of "$IDENTIFIER""

	./femswe swenml/swenml_"$IDENTIFIER".in > /dev/null
	mv run000001_restart_00000020.dat restart_files/run000001_restart_00000020_"$IDENTIFIER".dat 

	echo "simulation of "$IDENTIFIER" SFC"

	./femswe swenml/swenml_"$IDENTIFIER"_sfc.in > /dev/null
	mv run000001_restart_00000020.dat restart_files/run000001_restart_00000020_"$IDENTIFIER"_sfc.dat 
done
