# iModel
Icosahedral grid tools for geophysical fluid dynamics modelling

Pedro Peixoto Oct 2015
pedrosp@ime.usp.br

-------------------------------------------------------

iModel is a pack of tools to work with icosahedral and Voronoi geodesic grid based geophysical fluid models. Fully Fortran written, with outputs using GMT/Matlab. It contains:

 -Grid generator and grid tools, including grid optimization
- Interpolation and vector reconstruction pack
- Multigrid solver
- Transport model (Semi-Lagrangian)
- Shallow water model

Please read doc/manual.pdf for further information.

- Runs on Linux (tested on Debian and Ubuntu) 
- Runs on Windows-Cygwin 
- Fully written in Fortran 90 

1) Use the Makefile to compile (just type 'make')

2) Run using "./imodel". These will call the necessary routines 
    for grid generation or modelling. Edit imodel.f90 
    for your purpuses. 

3) Mesh parameters must be set in par/mesh.par or other par/*.par files
   Other parameters must be set in par/ directory

3') Choose the simulation to be run in mesh.par (1, 2...)

4) It is recommented to have installed Intel Fortran Compiler (but gfortran may be used)

5) It is recommended to have Generic Mapping Tool installed as visualization tool

6) Problems? Send me an e-mail

