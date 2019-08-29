module highorder
  !=============================================================================
  !  HIGH-ORDER module
  ! 
  !	Pack for several simulations on transport on deformational flow simulation
  !  on the sphere using Voronoi and Donalds grids
  ! 
  ! Jeferson Brambatti Granjeiro (brambatti@usp.br)
  !	Jun 2019
  !=============================================================================

  !Use global constants and kinds
  use constants, only: &
       erad, &
       i4, &
       pardir, &
       pi, &
       pi2, &
       r8, &
       unitspharea

  !Use main grid data structures
  use datastruct, only: &
       grid_structure, &
       scalar_field
       
  !Use routines from the spherical mesh pack
  use smeshpack
  
  implicit none

  !Global variables

  !Flags
  integer (i4):: testcase !Test case
  integer (i4):: initialfield !Initial condition

  !Time variables
  real (r8):: w   !Angular speed
  real (r8):: dt  !Time step
  real (r8):: T   !Period

  !Number of time steps and plots
  integer (i4):: ntime
  integer (i4):: nplots

  !Logical for plots or not
  logical:: plots

  !Plotsteps - will plot at every plotsteps timesteps
  integer (i4):: plotsteps

  !Method order 
  integer (i4):: order
 
  !Name for files and kind of staggering
  character (len=128)::  transpname
  character (len=8)::  stag
  character (len=8)::  radius 

  !Kind of interpolation for scalar and vector fields
  character (len=64):: ksinterpol
  character (len=64):: kvinterpol
  
  
  
  !======================================================================================
  ! ESTRUTURA DO METODO HIGH-ORDER

  !------------------------------------------------------------------
  ! node structure
  !------------------------------------------------------------------


  type,private  :: ngbr_structure
    ! Numero de vizinhos de cada no (primeiros vizinhos ou segundos vizinhos)
    integer(i4)       :: numberngbr
    
    ! Listagem dos primeiros vizinhos de um determinado no
    integer(i4),allocatable     :: lvv(:)
    
    ! Listagem das distancias dos primeiros vizinhos de um determinado no
    real(r8),allocatable        :: lvd(:)
  end type ngbr_structure  
  
  
  type,private  :: node_structure


  ! Informacoes de cada grupo de vizinhos
  type (ngbr_structure),allocatable  :: ngbr(:)
    
  end type node_structure

  

  
  
  
  
  
  
  
 
  !Criando a estrutura no
  type (node_structure),allocatable  :: node(:)  

  !======================================================================================  

contains   

  !======================================================================================
  !    HIGH-ORDER TESTS
  !======================================================================================

  subroutine highordertests(mesh)
  !-----------------------------------------
  !  Main transport tests routine
  !-----------------------------------------
 
  implicit none  
  integer(i4) :: i
  integer(i4) :: nodes
  integer(i4) :: ngbr
  integer(i4) :: nlines
  integer(i4) :: ncolumns

  !Numero de vizinhos e vizinhos de cada no
  integer(i4),allocatable   :: nbsv(:,:)   
  
  !Mesh
  type(grid_structure) :: mesh

  !Read parameters
  call gettransppars(mesh)    
    
  nodes = size(mesh%v(:)%p(1))

  !Alocando nos
  allocate(node(0:nodes))
  
  !Alocacao basica de memoria
  do i = 1,nodes
    allocate(node(i)%ngbr(1:order-1))
    ngbr = mesh%v(i)%nnb
       
    allocate(node(i)%ngbr(1)%lvv(1:ngbr+1))
    allocate(node(i)%ngbr(1)%lvd(1:ngbr+1))
  
    node(i)%ngbr(1)%numberngbr = ngbr
 
    !Armazena os primeiros vizinhos e suas respectivas distancias
    node(i)%ngbr(1)%lvv(1:ngbr+1) = (/i,mesh%v(i)%nb(1:ngbr)/)
    node(i)%ngbr(1)%lvd(1:ngbr+1) = (/0.0D0,mesh%v(i)%nbdg(1:ngbr)/)

    end do  

  
  print*, nodes , order
  read(*,*)
  
  if (order > 2) then
  nlines=maxval(mesh%v(:)%nnb)
  ncolumns=maxval(mesh%v(:)%nnb)+1  
  allocate(nbsv(13,nodes))
  
  nbsv = 0
  call find_neighbors(nbsv,nlines,ncolumns,nodes)
  
  do i=1,nodes
  print*,i
  end do 
  
  end if 

  end subroutine highordertests


  subroutine gettransppars(mesh)
    !---------------------------------------------------
    ! gettransppars
    !    Reads transport test parameters from file named "highorder.par"
    !    Saves parameters on global variables
    !--------------------------------------------------
    !Grid structure
    type(grid_structure), intent(in) :: mesh

    !Filename with parameters
    character (len=256):: filename

    !File unit
    integer (i4):: fileunit

    !Buffer for file reading
    character (len=300):: buffer

    !Flag to set num of time steps adjusted by
    !  grid level
    integer (i4):: adjustntime

    !Temp char
    character (len=64):: atmp

    !Couters
    integer(i4)::i
    integer(i4)::j

    !Standard parameters file
    filename=trim(pardir)//"highorder.par"
    print*,"Transport parameters (file): ", trim(filename)
    print*
    call getunit(fileunit)

    !A parameters file must exist
    open(fileunit,file=filename,status='old')

    read(fileunit,*)  buffer
    read(fileunit,*)  buffer
    read(fileunit,*)  order
    read(fileunit,*)  buffer     
    read(fileunit,*)  radius  
    read(fileunit,*)  buffer        
    read(fileunit,*)  testcase
    read(fileunit,*)  buffer
    read(fileunit,*)  initialfield
    read(fileunit,*)  buffer
    read(fileunit,*)  ntime, adjustntime
    read(fileunit,*)  buffer
    read(fileunit,*)  stag
   
    close(fileunit)

    !Number of iterations for a whole cycle of T
    !number of time steps
    if(adjustntime == 1)then
       if(mesh%glevel < 5 )then
          ntime= ntime / (2**(5-mesh%glevel))
       else
          ntime= ntime * (2**(mesh%glevel-5))
       end if
    end if

    !Set a standart name for files
    write(atmp,'(i8)') int(order)
    transpname="order"//trim(adjustl(trim(atmp)))    
    write(atmp,'(i8)') int(testcase)
    transpname=trim(adjustl(trim(transpname)))//"_v"//trim(adjustl(trim(atmp)))
    write(atmp,'(i8)') int(initialfield)
    transpname=trim(adjustl(trim(transpname)))//"_in"//trim(adjustl(trim(atmp)))    

    return
  end subroutine gettransppars
  
    subroutine find_neighbors(nbsv,nlines,ncolumns,nodes)  
    integer(i4),intent(in)    :: nlines
    integer(i4),intent(in)    :: ncolumns
    integer(i4),intent(in)    :: nodes
    integer(i4),dimension(13,nodes),intent(out)     :: nbsv
      
    integer     :: i
    integer     :: k
    integer     :: j
    integer     :: l
    integer     :: m
    integer     :: dimv
    integer     :: lend
    logical     :: logic
    integer(i4),allocatable   :: nbv(:,:)
    integer(i4),allocatable   :: vecv(:)
      
    allocate(nbv(nlines,ncolumns))
    allocate(vecv(size(nbv)))

    do i=1,nodes
      dimv=node(i)%ngbr(1)%numberngbr+1
      if (dimv-1 > nlines) then
        write(*,"(5x,'Problema com dimensao na matriz NB',/)")
        stop
      end if
        
    nbv = 0
    do k=2,dimv
      do j=1,node(node(i)%ngbr(1)%lvv(k))%ngbr(1)%numberngbr+1
        nbv(k-1,j) = node(node(i)%ngbr(1)%lvv(k))%ngbr(1)%lvv(j)
      end do
    end do

    do j=1,k-2
       l=1
       lend = node(node(i)%ngbr(1)%lvv(j+1))%ngbr(1)%numberngbr+1
       do while (l <= lend)
        logic = .false.
        m = 1
        do
          if (nbv(j,l) == node(i)%ngbr(1)%lvv(m)) then
            nbv(j,l) = 0
            logic = .true.
          end if
          m = m+1
            if (m >= 1+size(node(i)%ngbr(1)%lvv)) logic = .true.
              if (logic) exit
         end do
         l = l+1
       end do
    end do
        
    vecv = 0
    do k = 1,nlines
      do l = 1,ncolumns
        vecv(l+(k-1)*ncolumns) = nbv(k,l)
      end do
    end do
        
    !Eliminando os vertices repetidos
    do k = 1,size(nbv)-1
      l = k
      logic = .false.
      do 
        l = l+1
        if (vecv(k) /= 0) then
          if (vecv(k) == vecv(l)) then
            vecv(l) = 0
            logic = .true.
          end if
        else
          logic = .true.
        end if
        if (l >= size(nbv)) logic = .true.
        if (logic) exit
      end do
    end do

    !Contando o numero de vizinhos de cada no
    l = 0
    do k = 1,size(vecv)
      if (vecv(k) /= 0) then
        l = l+1
        if (l >= 13) then
          write(*,"(5x,'Mais vizinhos do que o esperado. Esperando no maximo',I8,/)")13
          stop
        end if
        nbsv(l+1,i) = vecv(k)
      end if
    end do
        
    !Salvando o numero de segundos vizinhos de cada no
    nbsv(1,i) = l
  end do

    deallocate(nbv,vecv)
    return
  end subroutine find_neighbors

  

end module highorder
