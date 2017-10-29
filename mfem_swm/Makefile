#TODO: this is just straight-forward makefile.
# Replace e.g. with autotools or scons

FC:=gfortran
FFLAGS=-O3 #-g3

CC:=gcc

all:	release

release:	FFLAGS=-O3 #-g3
release:	executables

openmp:	FFLAGS=-O3 -g3 -fopenmp
openmp:	executables

debug:		FFLAGS=-O0 -g3 -fcheck=all
debug:		executables

debug_openmp:		FFLAGS=-O0 -g3 -fcheck=all -fopenmp
debug_openmp:		executables


intel:		FFLAGS=-O3 -fopenmp -parallel
intel:		FC=ifort
intel:		executables

intel_debug:	FFLAGS=-O0 -g3 -check all
intel_debug:	FC=ifort
intel_debug:	executables

intel_serial:		FFLAGS=-O3 -g3
intel_serial:		FC=ifort
intel_serial:		executables


executables:	config	gengrid_hex gengrid_cube buildop_fem femswe

buildop_fem:	buildop_fem.f90
	$(FC) $(FFLAGS) grid.f90 buildop_fem.f90 -o buildop_fem

gengrid_hex:	gengrid_hex.f90
	$(FC) $(FFLAGS) grid.f90 gengrid_hex.f90 -o gengrid_hex


gengrid_cube:	gengrid_cube.f90
	$(FC) $(FFLAGS) grid.f90 gengrid_cube.f90 -o gengrid_cube

femswe:	femswe.f90
	$(FC) $(FFLAGS) grid.f90 femswe.f90 -o femswe

adjacency_matrices: adjacency_matrices.f90
	$(FC) $(FFLAGS) adjacency_matrices.f90 -o adjacency_matrices

sfc_optimize: sfc_optimize.f90
	$(CC) -c -o hilbert2d.o computeHilbert2d.c
	$(CC) -c -o hilbert3d.o computeHilbert3d.c
	$(FC) $(FFLAGS) -c grid.f90
	$(FC) $(FFLAGS) -c -o sfc_optimize.o sfc_optimize.f90
	$(FC) $(FFLAGS) hilbert2d.o hilbert3d.o grid.o sfc_optimize.o -o sfc_optimize

#Configure Enviroment (directories)
config:
	chmod +x *.sh
	. ./dirs.sh

run:
	./gengrid_hex
	./gengrid_cube
	./buildop_fem
	./femswe


clean:
	rm -f buildop_fem
	rm -f femswe
	rm -f gengrid_hex
	rm -f gengrid_cube
	rm -f adjacency_matrices
	rm -f sfc_optimize
	rm -f fort.43
	rm *.mod
