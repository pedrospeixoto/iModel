module datastructmult
  !====================================================================
  ! Data structures for multigrid method
  ! Marline I. Silva - Sept 2014
  !=====================================================================

 !Use global constants and kinds
  use constants, only: i4, r8

  use datastruct, only: &
    grid_structure

 
  type cell_relation
  !----------------------------------------------------
  ! Cellwise node indexes relating 2 different grid levels
  !----------------------------------------------------

     !two indices coarse mesh cell for localization index cell
     !fine mesh  
     integer(i4),dimension(1:2)::cel

     !index cell fine mesh between two cells 
     integer(i4)::celm
     
     !index cell coarse mesh and fine mesh
     integer(i4)::cmg
     
     !index cell fine mesh contained in the cell coarse mesh
     integer(i4)::cmf

  end type cell_relation

  type grid_relation
  !--------------------------------------------------------
  !  Full mesh node indexes relation between 2 grid levels
  !--------------------------------------------------------
     !list of nodes fine mesh not contained in coarse mesh
     type(cell_relation), allocatable::indmf(:)
     
     !list of nodes cell fine mesh
     type(cell_relation), allocatable::cf(:)

     !list of nodes cell coarse mesh
     type(cell_relation), allocatable::cng(:)

     !nodes cell coarse mesh
     integer(i4), allocatable::cg(:)
     
  
  end type grid_relation



  type multimesh
    !-------------------------------------------
    ! Hierarchy of grids with several levels
    !-------------------------------------------
     !Lista de malhas 
     type(grid_structure), allocatable::mesh(:)

     !Fine mesh
     type(grid_relation),allocatable::mf(:)

     !Coarser mesh
     type(grid_relation),allocatable::mg(:)

  end type multimesh


  type multisol2
  !----------------------------------------
  ! Function vectors for solution and residual
  !------------------------------------------
     !function value
     real(r8),allocatable::valf(:)

     !Residue value
     real(r8),allocatable::resf(:)

  end type multisol2

  type multisol
    !------------------------------------
    ! Solution for each grid level
    !-----------------------------------
     !value at each level
     type(multisol2), allocatable::solm(:)
     
   end type multisol


  type interpolable_meshes
     ! Values array, ordered in the same sequence as the
     !   valeus function

     real (r8), allocatable  :: f(:)
  end type interpolable_meshes

end module datastructmult
