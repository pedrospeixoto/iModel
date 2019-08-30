#--------------------------------------
# Makefile for IMODEL/IGRID routines
#  Geodesical mesh and atm modeling 
#     for the sphere
#
# P. Peixoto - Jul 2012
#----------------------------------

#Set 1 if you want to debug
# run as 'make DEBUG=1' to debug
DEBUG := 0
#DEBUG := 1

#Check for operating system 
OS := $(shell uname -s)

#Check if ifort exists
F90 := $(shell which ifort)
ifdef F90
	F90 = ifort
ifeq ($(OS), Linux)
ifeq ($(DEBUG), 1)
     #FFLAG := -O0 -traceback -debug extended -check all -warn -ftrapuv -check noarg_temp_created
     FFLAG := -qopenmp -traceback -debug extended -check noarg_temp_created -ftrapuv -check all
     #-fpp
else
     FFLAG := -O3 -qopenmp -parallel -heap-arrays
     #-fpp
     #FFLAG := -openmp -O0 -ipo
     #FFLAG := -O3 -xHOST
     #FFLAG := -openmp -O3 -ipo -parallel -xHOST -traceback -debug extended -check all -warn -ftrapuv -check noarg_temp_created 
endif
else #Windows/Cygwin
ifeq ($(DEBUG), 1)
     FFLAG := /traceback /check
else
     FFLAG := -Qopenmp -O3 -Qipo -Qparallel -QxHOST
endif
endif
endif

#If gfortran set from flag (make F90=gfortran)
ifeq ($(F90), gfortran)
	F90 = gfortran
	F90 := $(shell which gfortran)
ifeq ($(DEBUG), 1)
     FFLAG := -O0 -g -fbounds-check -Wall -Wno-unused -Wno-conversion
     #-cpp
else
     FFLAG := -O3 -fopenmp -march=native 
     #-cpp
     #FFLAG := -O3 -flto -fopenmp -march=native
endif
endif

#If ifort not available, use gfortran
#   same flags for Linux and Windows/Cygwin
ifndef F90
	F90 = gfortran
	F90 := $(shell which gfortran)
ifeq ($(DEBUG), 1)
     FFLAG := -O0 -g -fbounds-check -Wall -Wno-unused -Wno-conversion
else
     FFLAG := -O3 -fopenmp
endif
endif

#Include path - for modules created
IMOD=-Ibin 

#Objects
OBJ=bin/constants.obj \
bin/stripackm.obj \
bin/datastruct.obj \
bin/refinter.obj \
bin/smeshpack.obj \
bin/interpack.obj \
bin/diffoperpack.obj \
bin/simulpack.obj \
bin/swm_data.obj \
bin/swm_operators.obj \
bin/swm.obj \
bin/datastructmult.obj \
bin/transport.obj \
bin/poisson.obj \
bin/multigrid.obj \
bin/eispack.obj \
bin/highorder.obj \

#Compile and build all
all: header config bin/imodel ending

#Make all and run executable
run: all
	./imodel

#Print heading
header:
	@echo --------------------------------------------
	@echo Compiling and building the software   
	@echo --------------------------------------------
	@echo 
	@echo For DEBUG use DEBUG=1
	@echo For COMPILER use F90=gfortran or F90=ifort
	@echo 

#Configure Enviroment (directories)
config:
	chmod +x src/*.sh
	. src/dirs.sh


#Constants
bin/constants.obj: src/constants.f90
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv constants.mod bin/.

#Data struct
bin/datastruct.obj: src/datastruct.f90
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv datastruct.mod bin/.

#Data struct for multigrid
bin/datastructmult.obj: src/datastructmult.f90
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv datastructmult.mod bin/.

#Interpolation for local refinement pack #################################################################################
bin/refinter.obj: src/refinter.f90
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv refinter.mod bin/.

#Spherical triangulation pack
bin/stripackm.obj: src/stripackm.f90
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv stripackm.mod bin/.

#Interpolation pack
bin/interpack.obj: src/interpack.f90 
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv interpack.mod bin/.

#Numerical differential operators pack
bin/diffoperpack.obj: src/diffoperpack.f90 
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv diffoperpack.mod bin/.

#Spherical mesh pack
bin/smeshpack.obj: src/smeshpack.f90 
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv smeshpack.mod bin/.

#Simulation pack
bin/simulpack.obj: src/simulpack.f90 
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv simulpack.mod bin/.

#Multigrid pack
bin/poisson.obj: src/poisson.f90 
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv poisson.mod bin/.

#Multigrid pack
bin/multigrid.obj: src/multigrid.f90 
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv multigrid.mod bin/.

#Transport tests pack
bin/transport.obj: src/transport.f90 
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv transport.mod bin/.

#EISPACK - Eigenvalues calculation
bin/eispack.obj: src/eispack.f90 
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv eispack.mod bin/.

#Shallow water model
bin/swm_data.obj: src/swm_data.f90 
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv swm_data.mod bin/.

bin/swm_operators.obj: src/swm_operators.f90 
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv swm_operators.mod bin/.


bin/swm.obj: src/swm.f90 
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv swm.mod bin/.
	
#High-Order transport 	
bin/highorder.obj: src/highorder.f90 
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv highorder.mod bin/.


#Main executable
bin/imodel: src/imodel.f90 $(OBJ)
	$(F90) $(FFLAG)  src/imodel.f90 $(OBJ) -o $@ $(IMOD)

#Creates a link for executable and prints ending
ending: 
	chmod +x src/link.sh
	src/link.sh
	@echo End of compilation
	@echo
	@echo "Set parameter files (pars / *.par )" 
	@echo "   and then run 'imodel'"
	@echo "------------------------------------------------------------------"

#Clean targets
clean: 
	rm -rf bin/*.obj bin/*.o bin/*.mod	
	rm -rf bin/imodel*
	rm -rf *~

cleandata: clean
	rm -rf data/
	rm -rf graphs/
	rm -rf bin/
	rm imodel

cleangrids: clean
	rm -rf graphs/
	rm -rf grids/

cleanall: clean cleandata cleangrids

# Create a tar.bz2 file with all important files
archive: 
	chmod +x src/tarfiles.sh
	./src/tarfiles.sh

#Backup all important files
backup: 
	chmod +x src/backup.sh
	./src/backup.sh
