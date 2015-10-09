#--------------------------------------
# Makefile 
#
# P. Peixoto - Jul 2012
#----------------------------------

F90 := ifort
#F90 := gfortran
#FFLAG := -O0 -traceback -debug extended -check noarg_temp_created -warn
FFLAG := -mcmodel large -shared-intel -shared-libgcc -O3 -xHOST -heap-arrays
#FFLAG := -openmp -O3 -parallel -xHOST -mcmodel large -shared-intel -shared-libgcc
 
all: endgame

run: endgame
	./endgame

endgame: eg_sw_ref_psp.f90 
	$(F90) $(FFLAG) eg_sw_ref_psp.f90 -o endgame

.PHONY: clean

clean:
	rm -f endgame
	rm -f *.mod

archive:
	tar cvf eg_sw_ref.tar.gz *.f90 Makefile 


