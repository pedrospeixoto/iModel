# iModel
Icosahedral grid tools for geophysical fluid dynamics modelling

Pedro Peixoto Fev 2019
(pedrosp@ime.usp.br)

With contributions from MSc/PhD students.

-------------------------------------------------------

iModel is a pack of tools to work with icosahedral and Voronoi geodesic grid based geophysical fluid models. Fully Fortran written, with outputs using GMT/Matlab/Python. It contains:

- Grid generator and grid tools, including grid optimization
- Locally refined grids, topography based (MSc L. Santos)
- Interpolation and vector reconstruction pack
- Multigrid solver (MSc M. Silva)
- Transport model (Semi-Lagrangian)
- High order finite volume transport model (PhD Bramatti)
- Finite Volume Shallow Water Model 
- Finite Volume Moist Shallow Water Model (MSc L. Santos)

Additionally includes (Developed by John Thuburn):
- A spherical shallow water model version of ENDGame model, which uses finite differences on lat-long grid (endgame/)
- A spherical shallow water model of a mixed finite elements schemes, uses cubed shere or hexagonal grid (mfem_swm/)

Additionally includes (Developed by Pedro Peixoto)
- A planar shallow water model writen in Matlab using finite differences regular C-grid energy enstrophy conserving schemes (fdc_een_swm)

Please read doc/manual.pdf for further information and doc/ in submodel folders.

iModel:
--------
 
- Runs on Linux (tested on Debian and Ubuntu) 
- Fortran 90 (tested with ifort and gfortran)

1) Use the Makefile to compile (just type 'make')

2) Run using "./imodel". These will call the necessary routines 
    for grid generation or modelling. Edit imodel.f90 
    for your purpuses. 

3) Mesh parameters must be set in par/mesh.par or other par/*.par files
   Other parameters must be set in par/ directory

4) Choose the simulation to be run in mesh.par (1, 2...)

5) Output is written in data/
 
6) 
  a) Use GMT visualization tool scripts from gmt/plot.sh to plot output (Generic Mapping Tool need to be installed separately) - depreceated 
  b) Use Python scripts from pyscript folder for visualization.

7) Problems? Send an e-mail to ppeixoto@usp.br

----------------------------------------------------------------------------

Main references for iModel (see www.ime.usp.br/~pedrosp for exact reference):

- Santos, Peixoto, 2021: Topography based local spherical Voronoi grid refinement on classical and moist shallow-water finite volume models (Geoscientific Model Development)

- Yonggang Yu; Ning Wang; Jacques Middlecoff; Pedro Peixoto; Mark Govett, 2020: Comparing Numerical Accuracy of Icosahedral A-grid and C-grid Schemes in Solving the Shallow-Water Model (Monthly Weather Review) 

- Peixoto, Thuburn and Bell, 2018: Numerical instabilities of spherical shallow water models considering small equivalent depths ( Quarterly Journal of the Royal Meteorological Society)

- Bell, Peixoto, Thuburn, 2017: Numerical instabilities of vector invariant momentum equations on rectangular C-grids (Quart. J. Roy. Meteorol. Soc.) 

- Peixoto, 2016: Accuracy analysis of mimetic finite volume operators on geodesic grids and a consistent alternative (Journal of Computational Physics)
 
- Peixoto, PS and Barros, SRM, 2014 : On vector field reconstructions for semi-Lagrangian transport methods on geodesic staggered grids (Journal of Computational Physics) 

- Peixoto, PS and Barros, SRM, 2013 : Analysis of grid imprinting on geodesic spherical icosahedral grids (Journal of Computational Physics)


 

