#! /bin/bash

for i in 1 2 3 4 5 6 7 8 9 10; do
	
	#
	# CUBED SPHERE
	#
	NX=$((3*(2**(i-1))))
	CELLS=$((6*NX*NX))
	CELLS=`printf "%010d" $CELLS`

	IDENTIFIER="cube_$CELLS"
	OUTPUT_FILE="swenml_""$IDENTIFIER"".in"

	IDENTIFIER_SFC="cube_$CELLS""_sfc"
	OUTPUT_FILE_SFC="swenml_""$IDENTIFIER""_sfc.in"

	TIMESTEPSIZE=$(((1800*(2**5))/(2**$i)))
	TIMESTEPS=20

	echo "Cubed Sphere Level $i:"
	echo " +    Identifier: $IDENTIFIER"
	echo " + Timestep size: $TIMESTEPSIZE"
	echo " +     Timesteps: $TIMESTEPS"
	echo

	cp template_swenml.in $OUTPUT_FILE

	sed -i "s/%IDENTIFIER%/$IDENTIFIER/"  $OUTPUT_FILE
	sed -i "s/%TIMESTEPSIZE%/$TIMESTEPSIZE/"  $OUTPUT_FILE
	sed -i "s/%TIMESTEPS%/$TIMESTEPS/"  $OUTPUT_FILE


	cp template_swenml.in $OUTPUT_FILE_SFC

	sed -i "s/%IDENTIFIER%/$IDENTIFIER_SFC/"  $OUTPUT_FILE_SFC
	sed -i "s/%TIMESTEPSIZE%/$TIMESTEPSIZE/"  $OUTPUT_FILE_SFC
	sed -i "s/%TIMESTEPS%/$TIMESTEPS/"  $OUTPUT_FILE_SFC
done


for i in 1 2 3 4 5 6 7 8 9 10; do
	#
	# HEXAGONAL GRID
	#
	CELLS=$((5*2**(2*i-1)+2))
	CELLS=`printf "%010d" $CELLS`

	IDENTIFIER="hex_$CELLS"
	OUTPUT_FILE="swenml_""$IDENTIFIER"".in"

	IDENTIFIER_SFC="hex_$CELLS""_sfc"
	OUTPUT_FILE_SFC="swenml_""$IDENTIFIER""_sfc.in"

	TIMESTEPSIZE=$(((1800*(2**5))/(2**$i)))
	TIMESTEPS=20

	echo "   Hexagonal Level $i:"
	echo " +    Identifier: $IDENTIFIER"
	echo " + Timestep size: $TIMESTEPSIZE"
	echo " +     Timesteps: $TIMESTEPS"
	echo


	cp template_swenml.in $OUTPUT_FILE

	sed -i "s/%IDENTIFIER%/$IDENTIFIER/"  $OUTPUT_FILE
	sed -i "s/%TIMESTEPSIZE%/$TIMESTEPSIZE/"  $OUTPUT_FILE
	sed -i "s/%TIMESTEPS%/$TIMESTEPS/"  $OUTPUT_FILE


	cp template_swenml.in $OUTPUT_FILE_SFC

	sed -i "s/%IDENTIFIER%/$IDENTIFIER_SFC/"  $OUTPUT_FILE_SFC
	sed -i "s/%TIMESTEPSIZE%/$TIMESTEPSIZE/"  $OUTPUT_FILE_SFC
	sed -i "s/%TIMESTEPS%/$TIMESTEPS/"  $OUTPUT_FILE_SFC
done


