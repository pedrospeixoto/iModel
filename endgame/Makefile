#--------------------------------------
# Makefile 
#
# P. Peixoto - Jul 2012
#----------------------------------

#F90 := ifort
F90 := gfortran
#FFLAG := -O0 -traceback -debug extended -check noarg_temp_created -warn
#FFLAG := -mcmodel large -shared-intel -shared-libgcc -O3 -xHOST -heap-arrays
#FFLAG := -openmp -O3 -parallel -xHOST -mcmodel large -shared-intel -shared-libgcc

all: config endgame

run: endgame
	./endgame

endgame: src/eg_sw_ref_psp.f90 
	$(F90) $(FFLAG) src/eg_sw_ref_psp.f90 -o endgame -Ibin
	mv *.mod bin/
	
#Configure Enviroment (directories)
config:
	chmod +x *.sh
	. ./dirs.sh

.PHONY: clean

clean:
	rm -f endgame
	rm -f *.mod

archive:
	tar cvf eg_sw_ref.tar.gz src/*.f90 Makefile man/EG_refsoln_doc.pdf 


