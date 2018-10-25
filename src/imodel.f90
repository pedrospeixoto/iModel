program imodel
  !--------------------------------------------------
  !  Icosahedral Model Main Program - Driver
  !
  !   Pedro Peixoto Sept 2015
  !--------------------------------------------------

  !Global constants
  use constants, only: &
       simulcase

  !Data structures
  use datastruct, only: &
       grid_structure

  !Main mesh build and structuring module
  use smeshpack, only: &
       getparameters, &
       meshbuild, &
       printending, &
       printheader

  !Test routines
  use simulpack, only: &
       divergence_test, &
       laplacian_test, &
       rotational_test, &
       meshquality, &
       meshquality_tiledareas, &
       test_geo2reg, &
       test_edgeconnections, &
       scalar_interpolation_test, &
       test_trsearch, &
       vector_interpolation_test, &
       vector_reconstruction_test, &
       vec_tg_reconstruction_test

  !Transport testing module
  use transport, only: &
       passive_advection, &
       transptests

  !Data structure for multigrid routines
  use datastructmult

  !Multigrid solver routines
  use multigrid, only: &
    mesh_generation, &
    laplacian, &
    relaxac, &
    relation_meshes, &
    interpol_linear, &
    multigrid_solver


  !Spurconvergence analysis of poisson equaiton
  use poisson, only: &
    poisson_main


  use swm, only: &
    swm_tests, &
    swm_horiz_loc_trunc_er, &
    swm_normal_mode_analysis

  !Variable declaration
  implicit none

  type(grid_structure) :: mesh
  type(multimesh) :: meshvet

  !Print a header on screen
  call printheader()

  !Create/Load mesh

   !Read mesh user defined parameters and simulation case
   call getparameters(mesh)

  if(simulcase /= 11)then
	  !Call mesh generation and structuring algorithm
    call meshbuild(mesh)
  end if


  !Do a simulation/test with the mesh loaded
  select case(simulcase)
  case(1) !Test geodesic to regular grid conversion tool
     call test_geo2reg(mesh)
  case(2) !Test Search methods
     call test_trsearch(mesh)
     call test_edgeconnections(mesh)
  case(3) !Mesh distortion tests
     call meshquality(mesh)
     call meshquality_tiledareas(mesh)
  case(4) !Divergence Tests
     call divergence_test(mesh)
  case(5) !Divergence Tests
     call laplacian_test(mesh)
  case(6)   !Test scalar interpolations
     call scalar_interpolation_test(mesh)
  case(7) !Test vector interpolation
     call vector_interpolation_test(mesh)
  case(8) !Test vector reconstruction
     call vector_reconstruction_test(mesh)
  case(9) !Passive advection simulation
     call passive_advection(mesh)
  case(10) !Transport deformational flow simulation
     call transptests(mesh)
  case(11) !Multigrid tests

    !Generate multigrid meshes
    call mesh_generation(meshvet)

    !Test laplacian
    call laplacian(meshvet)

     !call getparameters_mult()

    !Test relaxation solver
    call relaxac(meshvet)

    !Set mesh relations
    call relation_meshes(meshvet)

    !Test interpolation
    call interpol_linear(meshvet)

    !Test multigrid
    call multigrid_solver(meshvet)

  case(12)!Tg reconstruction test
    call vec_tg_reconstruction_test(mesh)

  case(13) !Curl/rotational discretization tests
    call rotational_test(mesh)

  case(14) !Horizontal Distc Shallow water model disgnostics
    call swm_horiz_loc_trunc_er(mesh)

  case(15) !Shallow water model test cases
    call swm_tests(mesh)

  case(16) !Superconvergence tests of poisson equation
    call poisson_main(mesh)

  case default
     print*, "Please select a proper simulation case ...:", simulcase
  end select
  !Print finishing line
  call printending()

end program imodel

