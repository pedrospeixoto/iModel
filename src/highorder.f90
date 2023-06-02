module highorder
  !=============================================================================
  ! HIGH-ORDER module
  ! 
  !	Pack for several simulations on transport on deformational flow simulation
  ! on the sphere using Voronoi and Donalds grids
  ! 
  ! Jeferson Brambatti Granjeiro (jbram@usp.br)
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
       unitspharea, &
       day2sec

  !Use main grid data structures
  use datastruct, only: &
       grid_structure, &
       scalar_field, &
       vector
  !Use routines from the spherical inter pack
  use smeshpack
  use interpack

  !Routines from shallow water operators
  use swm_operators, only: &
    scalar_hx2ed, &
    scalar_elem_product, &
    div_hx, &
    zero_vector
  implicit none

  !Global variables

  !Flags
  integer (i4):: testcase !Test case
  integer (i4):: initialfield !Initial condition

  !Time variables
  real (r8):: w   !Angular speed
  real (r8):: dt  !Time step
  real (r8):: T   !Period
  real (r8):: H   !Height 
  real (r8):: radius

  character(len=64) :: atime

  !Number of time steps and plots
  integer (i4):: ntime
  integer (i4):: nplots

  !Logical for plots or not
  logical:: plots

  !Logical for moist shallow water model
  logical:: moistswm

  !Logical for monotonic filter
  logical:: monotonicfilter

  ! High order advection scheme variable
  character (len=6):: advmtd

  ! Time integrator 
  character (len=32):: time_integrator


  !Plotsteps - will plot at every plotsteps timesteps
  integer (i4):: plotsteps

  !Method order 
  integer (i4):: order
  
  !Name for files and kind of staggering
  character (len=128)::  transpname
  character (len=8)::  stag
  character (len=8)::  charradius 
  character (len=8)::  method
  character (len=8)::  controlvolume

  !Kind of interpolation for scalar and vector fields
  character (len=64):: ksinterpol
  character (len=64):: kvinterpol

  !======================================================================================
  ! ESTRUTURA DO METODO HIGH-ORDER
  !======================================================================================

  ! Condition initial 
  type (scalar_field):: phi
  type (scalar_field):: phi_new


  !--------------------------------------------------------------------------------------
  ! Aux variables for RK3 when monotonic filter is employed (following Wang et al 2009 notation)
  type(scalar_field):: phi_star
  type(scalar_field):: phi_RKS0
  type(scalar_field):: phi_RKS2
  type(scalar_field):: phi_tilda
  type(scalar_field):: phi_tilda_min
  type(scalar_field):: phi_tilda_max
  type(scalar_field):: phi_min
  type(scalar_field):: phieq
  type(scalar_field):: phi_max 
  type(scalar_field):: div_uphi
  type(scalar_field):: phi_ed
  type(scalar_field):: uphi

  ! RK3 variables
  real(r8), dimension(:), allocatable:: phif0
  real(r8), dimension(:), allocatable:: phif1
  real(r8), dimension(:), allocatable:: phif2
  real(r8), dimension(:), allocatable:: phif3

  ! Fluxes
  real(r8), dimension(:,:), allocatable:: F_star
  real(r8), dimension(:,:), allocatable:: F_cor
  real(r8), dimension(:,:), allocatable:: F_step2
  !--------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------
  ! node structure
  !--------------------------------------------------------------------------------------

  type :: coords_structure
     real(r8),dimension(1:2):: xy
     real(r8),dimension(1:3):: xyz
  end type coords_structure

  type :: coords2_structure
     real(r8),dimension(:,:),allocatable:: xyz2
  end type coords2_structure

  type :: edges_structure
     real(r8),dimension(:,:),allocatable:: xyz2
     real(r8),allocatable:: PG(:)
  end type edges_structure

  type :: edgesg_structure
     real(r8),allocatable:: MRGG(:,:)
     real(r8),allocatable:: MPGG(:,:)
     real(r8),allocatable:: VBGG(:)
     real(r8),allocatable:: coefg(:)
  end type edgesg_structure


  type :: gauss_structure
     real(r8),allocatable:: lpg(:,:)
     real(r8),allocatable:: lwg(:)
     real(r8),allocatable:: lvn(:,:)
     type (vector), allocatable:: velocity_quadrature(:) ! Normal velocity at quadrature points - used in the moist swm
     integer(i4),allocatable:: upwind_donald(:)
  end type gauss_structure

  type :: ngbr_structure
     !Numero de vizinhos de cada no (primeiros vizinhos ou segundos vizinhos)
     integer(i4)       :: numberngbr

     !Listagem dos primeiros vizinhos de um determinado no
     integer(i4),allocatable     :: lvv(:)

     !Listagem das distancias dos primeiros vizinhos de um determinado no
     real(r8),allocatable        :: lvd(:)
  end type ngbr_structure

  type :: flux_structure
     real(r8),allocatable    :: flux
  end type flux_structure

  type :: node_structure

     !Informacoes de cada grupo de vizinhos
     type (ngbr_structure),allocatable  :: ngbr(:)

     type (coords_structure),allocatable:: stencilpln(:)

     type (coords_structure),allocatable:: bar(:)

     type (coords_structure),allocatable:: mid(:)

     type (coords_structure),allocatable:: circ(:)

     type (coords2_structure),allocatable:: Pbar(:)

     type (coords2_structure),allocatable:: Pmid(:)

     type (coords2_structure),allocatable:: Pcirc(:)

     type (flux_structure),allocatable:: S(:)

     type (edges_structure),allocatable:: edge(:)

     type (edgesg_structure),allocatable:: edgeg(:)

     type (gauss_structure),allocatable:: G(:)

     integer(i4),allocatable:: stencil(:)

     !Matriz de Reconstrucao e Pseudoinversa FV-OLG
     real(r8),allocatable:: MRO(:,:)
     real(r8),allocatable:: MPO(:,:)

     !Vetor B FV-OLG
     real(r8),allocatable:: VBO(:)

     !Matriz de Reconstrucao e Pseudoinversa FV-GAS
     real(r8),allocatable:: MRG(:,:)
     real(r8),allocatable:: MPG(:,:)

     !Vetor B FV-GAS
     real(r8),allocatable:: VBG(:)

     !Coeficientes do polinomio
     real(r8),allocatable:: coef(:)

     !Momentos
     real(r8),allocatable:: moment(:)

     !Termos Geometricos
     real(r8),allocatable:: geometric(:,:)

     !Flux at each edge
     real(r8),allocatable:: edge_flux(:)

     !Condition initial
     real(r8):: phi

     !Area volume control 
     real(r8):: area

     !Tamanho do passo de tempo
     real(r8):: dt

     !Numero de passos de tempo
     integer(i4):: ntime

     !Solucao phi
     real(r8):: phi_new
     real(r8):: phi_new2
     real(r8):: phi_exa
     real(r8):: phi_old

     !Erro
     real(r8):: erro

     !Upwind FV-GAS
     integer(i4), allocatable:: upwind_voronoi(:)

  end type node_structure

  !Criando a estrutura no
  type (node_structure),allocatable  :: node(:)  

  ! Used in vector reconstruction- store shared edges indexes
  integer(i4), allocatable:: sh_edges_indexes(:,:)
  integer(i4), allocatable:: edges_indexes(:,:,:)

  !===============================================================================================  

contains   

  subroutine gettransppars(mesh)
    !----------------------------------------------------------------------------------------------
    ! gettransppars
    !    Reads transport test parameters from file named "highorder.par"
    !    Saves parameters on global variables
    !-----------------------------------------------------------------------------------------------
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

    !Monotonic limiter var
    character (len=8):: mono

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
    read(fileunit,*)  controlvolume
    read(fileunit,*)  buffer
    read(fileunit,*)  advmtd
    read(fileunit,*)  buffer
    read(fileunit,*)  time_integrator
    read(fileunit,*)  buffer  
    read(fileunit,*)  mono     
    read(fileunit,*)  buffer     
    read(fileunit,*)  charradius  
    read(fileunit,*)  buffer        
    read(fileunit,*)  testcase
    read(fileunit,*)  buffer
    read(fileunit,*)  initialfield
    read(fileunit,*)  buffer
    read(fileunit,*)  ntime, adjustntime
    read(fileunit,*)  buffer
    read(fileunit,*)  stag

    if(charradius == "U")then
       radius=1.0D0
       T=5.0D0
       H=1.0D0
    else
       radius=6371220.0D0
       T=1036800.0D0
       H=1000.0D0
    end if

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

    !Mixing ration
    phi%pos=0
    phi%n=mesh%nv
    allocate(phi%f(1:phi%n)) 
    allocate(phi_new%f(1:phi%n)) 
    phi_new=phi

    !Set a standart name for files
    write(atmp,'(i8)') int(order)
    transpname="order"//trim(adjustl(trim(atmp)))    
    write(atmp,'(i8)') int(testcase)
    transpname=trim(adjustl(trim(transpname)))//"_v"//trim(adjustl(trim(atmp)))
    write(atmp,'(i8)') int(initialfield)
    transpname=trim(adjustl(trim(transpname)))//"_in"//trim(adjustl(trim(atmp)))    

    ! Advection scheme order
    select case(advmtd)
    case('upw1')
        order=2 ! just to avoid allocation problems
        method='O'

    case('og2')
        order=2
        method='O'

    case('og3') 
        order=3
        method='O'

    case('og4')
        order=4
        method='O'

    case('sg3')
        order=3
        method='G'

    case default
        print*, "Invalid advection method. Please select a proper order."
        stop
    end select

    transpname=trim(transpname)//"_advmethod_"//trim(advmtd)
    if (time_integrator == 'rk3' .or. time_integrator == 'rk4') then
        transpname=trim(transpname)//"_"//trim(time_integrator)
    else
        print*, "Invalid time integrator. Please select a proper method."
        stop
    endif

    ! Monotonic limiter
    select case(mono)
    case('0')
       monotonicfilter = .False.
    case('1')
       monotonicfilter = .True.
    case default
        print*, "Invalid monotonic limiter. Please select a proper one."
        stop
    end select
    transpname=trim(transpname)//"_mono"//trim(mono)

   return
  end subroutine gettransppars


  !======================================================================================
  !    HIGH-ORDER TESTS
  !======================================================================================

  subroutine highordertests(mesh)
    !-----------------------------------------
    !  Main transport tests routine
    !-----------------------------------------

    implicit none  
    integer(i4) :: i
    integer(i4) :: j
    integer(i4) :: k
    integer(i4) :: nodes
    integer(i4) :: ngbr
    integer(i4) :: nlines
    integer(i4) :: ncolumns
    real(r8):: time
    real(r8):: mass, mass0, mass_var

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

    if (order > 2) then
       nlines=maxval(mesh%v(:)%nnb)
       ncolumns=maxval(mesh%v(:)%nnb)+1  
       allocate(nbsv(13,nodes))
       nbsv = 0
       call find_neighbors(nbsv,nlines,ncolumns,nodes)
       do i=1,nodes
          allocate(node(i)%ngbr(2)%lvv(1:nbsv(1,i)+1))
          allocate(node(i)%ngbr(2)%lvd(1:nbsv(1,i)+1))
          node(i)%ngbr(2)%numberngbr = nbsv(1,i)
          node(i)%ngbr(2)%lvv = (/i,nbsv(2:,i)/)
       end do
       deallocate(nbsv)
    end if

    call allocation(nodes,mesh)

    call stencil(nodes,mesh)

    if(controlvolume=="D") then
       call coordsbarmid(nodes,mesh)  

       call edges_donalds(nodes)     

       call area_donald(nodes,mesh)

       call gaussedges(nodes,mesh)  

       call upwind_donald(nodes,mesh)

    else

       call coordscirc(nodes,mesh)

       call edges_voronoi(nodes)

       call gaussedges(nodes,mesh)

       call upwind_voronoi(nodes,mesh)

    endif

    call time_step(nodes,mesh,0.5785D0) 

    !Method FV-GAS
    if(method=='G')then  
       print*, 'Method FV-GAS'

       call condition_initial_gas(nodes,mesh) 

       call matrix_gas(nodes,mesh)

       call vector_gas(nodes,mesh)

       call reconstruction_gas(nodes,mesh)

!       call interpolation(nodes,mesh) 

       !call flux_edges_gas(nodes,mesh,0,0.0D0) 

       call  init_quadrature_edges(mesh)
       ! GASSMANN NO PLANO
       !call matrix_g(nodes,mesh) 
       !call vector_g(nodes,mesh)
       !call reconstruction_g(nodes,mesh)
       !call interpolation_g(mesh)
       !call flux_edges_g(nodes,mesh,0,0.0D0)  

       do i=1,nodes
          node(i)%phi_exa=node(i)%phi_new2
       enddo

       ! Initial mass
       mass0 = sum(node(1:nodes)%phi_new2*mesh%hx(1:nodes)%areag)

       !Avanço temporal 
       time=0.0D0 
       do j=0,node(0)%ntime
          !Atualizando a solucao
          do i=1,nodes
             node(i)%phi_old=node(i)%phi_new2
          enddo

          if (j==node(0)%ntime)then
             call erro(nodes,mesh,10.0D0)
          endif
      
          !Plotando Inicial e Final
          if(j==0.or.j==node(0)%ntime/2 .or.j==node(0)%ntime)then  
             do i=1,nodes
                phi%f(i) = node(i)%phi_new2
             end do
             write(atime,'(i8)') nint(time)
             phi%name=trim(adjustl(trim(transpname))//"_phi_t_"//trim(adjustl(trim(atime))))
             call plot_scalarfield(phi,mesh)
          end if

          if (time_integrator=='rk4') then
             do k=1,4
                call vector_gas(nodes,mesh)
                call reconstruction_gas(nodes,mesh)
                call flux_gas(nodes,mesh,k-1,time)
                call rungekutta4(nodes,mesh,k) 
             enddo
          elseif (time_integrator=='rk3') then
            call ode_rk3_adv(nodes, mesh, time)
            !call ode_rk3_advection (mesh, phi_new, phi, time, node(0)%dt, radius)
            !node(1:nodes)%phi_new2=phi_new%f
          endif

          ! Mass update
          mass = sum(node(1:nodes)%phi_new2*mesh%hx(1:nodes)%areag)
          mass_var = abs((mass0-mass)/mass0)
          time=time+node(0)%dt
          print*, "Time:",  time
          print '(a33, 3e16.8)','min, max, mass variation = ',&
          minval(node(1:nodes)%phi_new2), maxval(node(1:nodes)%phi_new2), mass_var
          print*,''
       enddo

       !Method FV-OLG
    else
       print*, 'Method FV-OLG'
       call matrix_olg(nodes,mesh)

       call  init_quadrature_edges(mesh)
       !call reconstruction_olg(nodes,mesh)  

       !call interpolation(nodes,mesh) 

       !call flux_edges_gas(nodes,mesh,0,0.0D0) 

       !call flux_edges_olg(nodes,mesh,0,0.0D0) 

       !call divergente_olg(nodes,mesh,0,0.0D0)

       do i=1,nodes
          node(i)%phi_exa=node(i)%phi_new2
       enddo

       ! Initial mass
       if(controlvolume=='V')then
         mass0 = sum(node(1:nodes)%phi_new2*mesh%hx(1:nodes)%areag)
       else
         mass0 = sum(node(1:nodes)%phi_new2*node(1:nodes)%area)
       end if

       !Avanço temporal 
       time=0.0D0

       do j=0,node(0)%ntime
          !Atualizando a solucao
          do i=1,nodes
             node(i)%phi_old=node(i)%phi_new2
          enddo
 
          if (j==node(0)%ntime)then
             call erro(nodes,mesh,10.0D0)
          endif

          !Plotando Inicial e Final
          if(j==0.or.j==node(0)%ntime/2 .or.j==node(0)%ntime)then  
             do i=1,nodes
                phi%f(i) = node(i)%phi_new2
             end do
             write(atime,'(i8)') nint(time)
             phi%name=trim(adjustl(trim(transpname))//"_phi_t_"//trim(adjustl(trim(atime))))
             call plot_scalarfield(phi,mesh)
          end if

          if (time_integrator=='rk4') then
            do k=1,4
               call vector_olg2(nodes)
               call reconstruction_olg(nodes,mesh) 
               call flux_olg(nodes,mesh,k-1,time)
               call rungekutta4(nodes,mesh,k) 
            enddo
          elseif (time_integrator=='rk3') then
            call ode_rk3_adv(nodes, mesh, time)
            !call ode_rk3_advection (mesh, phi_new, phi, time, node(0)%dt, radius)
            !node(1:nodes)%phi_new2=phi_new%f
          endif

          time=time+node(0)%dt

          ! Mass update
          if(controlvolume=='V')then
            mass = sum(node(1:nodes)%phi_new2*mesh%hx(1:nodes)%areag)
          else
            mass = sum(node(1:nodes)%phi_new2*node(1:nodes)%area)
          end if

          mass_var = abs((mass0-mass)/mass0)

          print*, "Time:",  time
          print '(a33, 3e16.8)','min, max, mass variation = ',&
          minval(node(1:nodes)%phi_new2), maxval(node(1:nodes)%phi_new2), mass_var
          print*,''
          enddo
    end if

    deallocate(node)
    return
  end subroutine highordertests



  !----------------------------------------------------------------------------------------
  !  Initial conditions and vector fields
  !----------------------------------------------------------------------------------------
  function f(lon, lat, time)
    !-----------------------------------
    !  F - initial conditions for scalar fields
    ! 
    !   Lat in [-pi/2,pi/2], Lon in [-pi,pi]
    ! 
    !  P. Peixoto - Feb2013
    !---------------------------------------------
    real (r8), intent(in) :: lon
    real (r8), intent(in) :: lat
    real (r8), intent(in), optional :: time

    real (r8):: f

    !General parameters
    real (r8):: h1
    real (r8):: h2

    !Cosbell parameters
    real (r8):: b
    real (r8):: c
    real (r8):: r
    real (r8):: r1
    real (r8):: r2

    !Gaussian parameters
    real (r8), dimension(3) :: p
    real (r8), dimension(3) :: p1
    real (r8), dimension(3) :: p2
    real (r8):: b0

    !Center of bell, slot, hill
    real (r8):: lat1
    real (r8):: lon1
    real (r8):: lat2
    real (r8):: lon2

    !Copy of point position
    real (r8):: lonp
    real (r8):: latp

    !Exponential/Gaussian parameters
    real (r8):: L

    !Keep a local copy of point position
    if(present(time))then !Exact scalar for rotation
       if(testcase==5.or. testcase==7)then
          lonp=lon-time*2.*pi/T          
       end if
       if(lonp<-pi)then
          lonp=lonp+pi2
       elseif(lonp>pi)then
          lonp=lonp-pi2
       end if
    else
       lonp=lon
    end if
    latp=lat

    !Center of bell, hill or cilinder
    ! It depends on the test case
    select case(testcase)
    case(1) !Nondivergent case-1
       lon1=0._r8
       lat1=pi/(3._r8)
       lon2=0._r8
       lat2=-pi/(3._r8)
    case(3) !Divergent case 3
       lon1=-pi/(4._r8)
       lat1=0._r8
       lon2=pi/(4._r8)
       lat2=0._r8
    case(2, 4) !Nondivergent case-2 and case 4 (rotation
       lon1=-pi/(6._r8)
       lat1=0._r8
       lon2=pi/(6._r8)
       lat2=0._r8
    case (5, 7) !Passive Adv - Rotation
       lon1=0.
       lat1=0.
    case default
       lon1=0.
       lat1=0.
       lon2=0.
       lat2=0.
    end select

    !Initial scalar fields
    select case(initialfield)
    case(0) !Constant
       f=1._r8
    case(1) !Cosine bell
       b=0.1_r8
       c=0.9_r8
       r=1._r8/2._r8
       r1= arcdistll(lonp, latp, lon1, lat1)
       if(r1<r)then
          h1=(1._r8/2._r8)*(1._r8+dcos(pi*r1/r))
          f=b+c*h1
          return
       end if
       r2= arcdistll(lonp, latp, lon2, lat2)
       if(r2<r)then
          h2=(1._r8/2._r8)*(1._r8+dcos(pi*r2/r))
          f=b+c*h2
          return
       end if
       f=b
    case(2) !Gaussian
       b0=5
       call sph2cart(lonp, latp, p(1), p(2), p(3))
       call sph2cart(lon1, lat1, p1(1), p1(2), p1(3))
       call sph2cart(lon2, lat2, p2(1), p2(2), p2(3))
       h1=dexp(-b0*norm(p-p1)**2)
       h2=dexp(-b0*norm(p-p2)**2)
       f=h1+h2
    case(3) !Non-smooth- Slotted cylinder
       b=0._r8
       c=1._r8 !0.9
       r=1._r8/2._r8
       f=b

       r1= arcdistll(lonp, latp, lon1, lat1)
       if(r1<=r.and.abs(lonp-lon1)>=r/(6._r8))then
          f=c
          return
       end if
       r2= arcdistll(lonp, latp, lon2, lat2)
       if(r2<=r.and.abs(lon-lon2)>=r/(6._r8))then
          f=c
          return
       end if
       if(r1<=r.and.abs(lonp-lon1) < r/(6._r8) .and. &
            (latp-lat1)<-(5._r8)*r/(12._r8))then
          f=c
          return
       end if
       if(r2<=r.and.abs(lonp-lon2)<r/(6._r8).and. &
            (latp-lat2)>(5._r8)*r/(12._r8))then
          f=c
          return
       end if
    case(4) !Passive Adv - Rotation
       r = (latp-lat1)**2+(lonp-lon1)**2
       L=5*pi/180._r8
       f=exp(-r/L)
    case(5) !Global positive function
       f=dsin(latp)
       !f=(dcos(latp))**3*(dsin(lonp))**2
    case(6) !One Gaussian  
       b0=5
       call sph2cart(lonp, latp, p(1), p(2), p(3))
       call sph2cart(lon1, lat1, p1(1), p1(2), p1(3))
       f=dexp(-b0*norm(p-p1)**2)
    case(7) !Gaussian
       b0=5
       call sph2cart(lonp, latp, p(1), p(2), p(3))
       call sph2cart(lon1, lat1, p1(1), p1(2), p1(3))
       f=dexp(-b0*norm(p-p1)**2)       
    case default
       f=0.

    end select

    return
  end function f

  function velocity(p, time)
    !-----------------------------------
    !  V - velocity in Cartesian coordinates
    !
    !   p in Cartesian coords
    !
    !  P. Peixoto - Feb2013
    !---------------------------------------------
    real (r8), intent(in) :: p(1:3)
    real (r8):: velocity(1:3)

    !Period, actual time
    real (r8), intent(in)::  time

    !Auxiliar variables
    real (r8):: lon
    real (r8):: lat
    real (r8):: utmp
    real (r8):: vtmp

    call cart2sph(p(1), p(2), p(3), lon, lat)
    utmp=u(lon, lat, time)
    vtmp=v(lon, lat, time)
    call convert_vec_sph2cart(utmp, vtmp, p, velocity)
    return
  end function velocity

  function u(lon, lat, time)
    !---------------------------------------------------------------------------------------------
    !  U - velocity in West-East direction
    !   for a given testcase
    !
    !   Lat in [-pi/2,pi/2], Lon in [-pi,pi]
    !
    !  P. Peixoto - Feb2012
    !---------------------------------------------------------------------------------------------
    real (r8), intent(in) :: lon
    real (r8), intent(in) :: lat
    real (r8):: u

    real (r8):: omega_m
    real (r8):: omega_0
    real (r8):: a
    integer (i4):: wave_m

    !Period, actual time
    real (r8):: time

    !Auxiliar variables
    real (r8):: k
    real (r8):: lonp

    u=0
    ! Velocity for each testcase
    select case(testcase)
    case(1) !Nondivergent case-1
       k=2.4
       u=k*(dsin((lon+pi)/2.)**2)*(dsin(2.*lat))*(dcos(pi*time/T))
    case(2) !Nondivergent case-2
       k=2.
       u=k*(dsin((lon+pi))**2)*(dsin(2.*lat))*(dcos(pi*time/T))
    case(3) !Divergent case-3
       k=1.
       u=-k*(dsin((lon+pi)/2._r8)**2)*(dsin(2.*lat))*(dcos(lat)**2)*(dcos(pi*time/T))
    case(4) !Nondivergent case-4 - case 2 + solid body rotation
       k=2.
       lonp=lon-2*pi*time/T
       u=k*(dsin((lonp+pi))**2)*(dsin(2.*lat))*(dcos(pi*time/T))+2.*pi*dcos(lat)/T
    case (5, 7) !Passive Adv - Rotation
       u=2.*radius*pi*dcos(lat)/T
       
    case(6) !Rossby-Harwitz
       a=erad
       omega_m=7.848e-6_r8
       omega_0=7.848e-6_r8
       wave_m=4_i4
       u=a*omega_0*dcos(lat)+ &
            a*omega_m*(dcos(lat))**(wave_m-1)* &
            (wave_m*dsin(lat)**2 - dcos(lat)**2) * dcos(wave_m*lon)
       u=u/50._r8
    case default
       print*, "Error in transport tests: Unknown vector field ", testcase
       stop
    end select

    return
  end function u

  function v(lon, lat, time)
    !----------------------------------------------------------------------------------------------
    !  V - velocity in South-North direction
    !   for a given testcase
    !
    !   Lat in [-pi/2,pi/2], Lon in [-pi,pi]
    !
    !  P. Peixoto - Feb2012
    !----------------------------------------------------------------------------------------------
    real (r8), intent(in) :: lon
    real (r8), intent(in) :: lat
    real (r8):: v
    real (r8):: omega_m
    real (r8):: omega_0
    real (r8):: a
    integer (i4):: wave_m

    !Period, actual time
    real (r8)::  time

    !Auxiliar variables
    real (r8):: k
    real (r8):: lonp

    ! Velocity for each testcase
    v=0
    select case(testcase)
    case(1) !Nondivergent case-1
       k=2.4
       v=k*(dsin(lon+pi))*(dcos(lat))*(dcos(pi*time/T))/2.
    case(2) !Nondivergent case-2
       k=2
       v=k*(dsin(2*(lon+pi)))*(dcos(lat))*(dcos(pi*time/T))
    case(3) !Divergent case-3
       k=1
       v=(k/2._r8)*(dsin((lon+pi)))*(dcos(lat)**3)*(dcos(pi*time/T))
    case(4) !Nondivergent case-2
       k=2
       lonp=lon-2*pi*time/T
       v=k*(dsin(2*(lonp+pi)))*(dcos(lat))*(dcos(pi*time/T))
    case(5,7) !Rotation
       v=0.
    case(6) !Rossby Harwitz
       a=erad
       omega_m=7.848e-6_r8
       omega_0=7.848e-6_r8
       wave_m=4_i4
       v = - a*omega_m*wave_m*((dcos(lat))**(wave_m-1))* &
            dsin(lat)*dsin(wave_m*lon)
       v=v/50._r8
    case default
       print*, "Error in transport tests: Unknown vector field ", testcase
       stop
    end select

    return
  end function v

  subroutine find_neighbors(nbsv,nlines,ncolumns,nodes)  
    !----------------------------------------------------------------------------------------------
    !    Procurando os primeiros e segundos vizinhos 
    !----------------------------------------------------------------------------------------------
    implicit none 
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

  subroutine allocation(nodes,mesh)
    !----------------------------------------------------------------------------------------------
    !    Allocation 
    !----------------------------------------------------------------------------------------------
    implicit none
    integer(i4),intent(in)         :: nodes
    type(grid_structure),intent(in):: mesh
    integer(i4):: i
    integer(i4):: j
    integer(i4):: k
    integer(i4):: l
    integer(i4):: VD
    integer(i4):: ngbr1
    integer(i4):: ngbr2
    !Error flag
    integer(i4):: ist


    do i=1,nodes
       ngbr1=node(i)%ngbr(1)%numberngbr
       allocate(node(i)%S(0:3))
       do k=0,3
          allocate(node(i)%S(k)%flux)
       end do

       if(controlvolume=="D")then
          VD=2
          allocate(node(i)%bar(1:ngbr1))
          allocate(node(i)%mid(1:ngbr1))            
       else
          VD=1
          allocate(node(i)%circ(1:ngbr1))      
       endif

       allocate(node(i)%edge(VD*ngbr1))
       allocate(node(i)%edgeg(ngbr1+1))

       allocate(node(i)%edge_flux(VD*ngbr1))

       do j=1,VD*ngbr1
          allocate(node(i)%edge(j)%xyz2(2,3))
          allocate(node(i)%edge(j)%pg(2))
       enddo


       if(order==2)then
          allocate(node(i)%stencil(0:ngbr1))
          allocate(node(i)%stencilpln(ngbr1))
          allocate(node(i)%MRO(ngbr1+1,3))
          allocate(node(i)%MPO(1:2,ngbr1))
          allocate(node(i)%VBO(1:ngbr1+1))
          allocate(node(i)%moment(1:3))
          allocate(node(i)%geometric(ngbr1,3))
          allocate(node(i)%coef(1:6))

          allocate(node(i)%G(1))
          allocate(node(i)%G(1)%lwg(2*ngbr1))
          allocate(node(i)%G(1)%lpg(2*ngbr1,3))
          allocate(node(i)%G(1)%velocity_quadrature(2*ngbr1))
          allocate(node(i)%G(1)%lvn(2*ngbr1,3))
          allocate(node(i)%G(1)%upwind_donald(2*ngbr1))

          allocate(node(i)%upwind_voronoi(1:ngbr1))

       end if

       if(order>2)then
          ngbr2=node(i)%ngbr(2)%numberngbr
          allocate(node(i)%stencil(0:ngbr1+ngbr2))
          allocate(node(i)%stencilpln(ngbr1+ngbr2))
          allocate(node(i)%G(2)) ! Original era G(3)
          allocate(node(i)%G(1)%lwg(2*ngbr1))
          allocate(node(i)%G(1)%lpg(2*ngbr1,3))
          allocate(node(i)%G(1)%lvn(2*ngbr1,3))
          allocate(node(i)%G(1)%velocity_quadrature(2*ngbr1))
          allocate(node(i)%G(2)%velocity_quadrature(2*ngbr1))
          allocate(node(i)%G(1)%upwind_donald(2*ngbr1))
          allocate(node(i)%G(2)%lwg(2*ngbr1))
          allocate(node(i)%G(2)%lpg(2*ngbr1,3))
          allocate(node(i)%G(2)%lvn(2*ngbr1,3))
          allocate(node(i)%G(2)%upwind_donald(4*ngbr1))
       end if

       if(order == 3)then
          allocate(node(i)%MRO(ngbr1+ngbr2+1,6))
          allocate(node(i)%MPO(1:5,ngbr1+ngbr2))
          allocate(node(i)%VBO(1:ngbr1+ngbr2+1))
          allocate(node(i)%moment(1:6))
          allocate(node(i)%geometric(ngbr1+ngbr2,6))
          allocate(node(i)%coef(1:6))

          allocate(node(i)%MRG(ngbr1,5))
          allocate(node(i)%MPG(5,ngbr1))
          allocate(node(i)%VBG(1:ngbr1))

          allocate(node(i)%upwind_voronoi(1:ngbr1))


          allocate(node(i)%edgeg(1)%MRGG(mesh%v(i)%nnb,5))
          allocate(node(i)%edgeg(1)%MPGG(5,mesh%v(i)%nnb))
          allocate(node(i)%edgeg(1)%VBGG(1:mesh%v(i)%nnb))
          allocate(node(i)%edgeg(1)%coefg(1:6))         
          do j=2,ngbr1+1
             allocate(node(i)%edgeg(j)%MRGG(mesh%v(mesh%v(i)%nb(j-1))%nnb,5))
             allocate(node(i)%edgeg(j)%MPGG(5,mesh%v(mesh%v(i)%nb(j-1))%nnb))
             allocate(node(i)%edgeg(j)%VBGG(1:mesh%v(mesh%v(i)%nb(j-1))%nnb))
             allocate(node(i)%edgeg(j)%coefg(1:6))
          enddo


       end if

       if(order == 4)then
          allocate(node(i)%MRO(ngbr1+ngbr2+1,10))
          allocate(node(i)%MPO(1:9,ngbr1+ngbr2))
          allocate(node(i)%VBO(1:ngbr1+ngbr2+1))
          allocate(node(i)%moment(1:10))
          allocate(node(i)%geometric(ngbr1+ngbr2,10))
          allocate(node(i)%coef(1:10))

          allocate(node(i)%upwind_voronoi(1:ngbr1))

       end if
     end do

     ! RK vars
     allocate(phif0(1:mesh%nv), stat=ist)
     allocate(phif1(1:mesh%nv), stat=ist)
     allocate(phif2(1:mesh%nv), stat=ist)
     allocate(phif3(1:mesh%nv), stat=ist)
     div_uphi%n = mesh%nv
     uphi%n = mesh%ne
     phi_ed%n = mesh%ne
     allocate(div_uphi%f(1:div_uphi%n), stat=ist)
     allocate(uphi%f(1:uphi%n), stat=ist)
     allocate(phi_ed%f(1:phi_ed%n), stat=ist)
     phieq%n=mesh%nv
     phieq%pos=0
     allocate(phieq%f(1:phieq%n), stat=ist)

     if(time_integrator=='rk3' .and.monotonicfilter)then
          phi_star%n = mesh%nv
          phi_RKS0%n = mesh%nv
          phi_RKS2%n = mesh%nv
          phi_tilda%n = mesh%nv
          phi_tilda_min%n = mesh%nv
          phi_tilda_max%n = mesh%nv
          phi_min%n = mesh%nv
          phi_max%n = mesh%nv
  
          allocate(phi_star%f(1:phi_star%n), stat=ist)
          allocate(phi_RKS0%f(1:phi_RKS0%n), stat=ist)
          allocate(phi_RKS2%f(1:phi_RKS2%n), stat=ist)
          allocate(phi_tilda%f(1:phi_tilda%n), stat=ist)
          allocate(phi_tilda_min%f(1:phi_tilda_min%n), stat=ist)
          allocate(phi_tilda_max%f(1:phi_tilda_max%n), stat=ist)
          allocate(phi_min%f(1:phi_min%n), stat=ist)
          allocate(phi_max%f(1:phi_max%n), stat=ist)

          if(controlvolume=='V')then
            allocate(F_star(1:mesh%nv, 1:maxval(mesh%v(:)%nnb)), stat=ist)
            allocate(F_step2(1:mesh%nv, 1:maxval(mesh%v(:)%nnb)), stat=ist)
            allocate(F_cor(1:mesh%nv, 1:maxval(mesh%v(:)%nnb)), stat=ist)

          else if(controlvolume=='D')then
            allocate(F_star(1:mesh%nv, 1:2*maxval(mesh%v(:)%nnb)), stat=ist)
            allocate(F_step2(1:mesh%nv, 1:2*maxval(mesh%v(:)%nnb)), stat=ist)
            allocate(F_cor(1:mesh%nv, 1:2*maxval(mesh%v(:)%nnb)), stat=ist)

          end if
      end if
 
    return
  end subroutine allocation

  subroutine stencil(nodes,mesh)
    !----------------------------------------------------------------------------------------------
    !   Determinando o estencil para os metodos de 2, 3 e 4 ordens 
    !----------------------------------------------------------------------------------------------  
    implicit none
    integer(i4),intent(in):: nodes
    type(grid_structure),intent(in):: mesh
    integer(i4):: i
    integer(i4):: j
    integer(i4):: l
    integer(i4):: m
    integer(i4):: jend
    real(r8):: cx
    real(r8):: cy
    real(r8):: sx
    real(r8):: sy
    real(r8):: x
    real(r8):: y
    real(r8):: z
    real(r8):: xp
    real(r8):: yp
    real(r8):: zp

    if(order==2)then
       do i=1,nodes
          x=mesh%v(i)%p(1)
          y=mesh%v(i)%p(2)
          z=mesh%v(i)%p(3) 
          call constr(x,y,z,cx,sx,cy,sy)
          jend=node(i)%ngbr(1)%numberngbr
          node(i)%stencil(0)=node(i)%ngbr(1)%lvv(1)
          do j=1,jend
             node(i)%stencil(j)=node(i)%ngbr(1)%lvv(j+1)
             x=mesh%v(node(i)%ngbr(1)%lvv(j+1))%p(1)
             y=mesh%v(node(i)%ngbr(1)%lvv(j+1))%p(2)
             z=mesh%v(node(i)%ngbr(1)%lvv(j+1))%p(3)
             call aplyr(x,y,z,cx,sx,cy,sy,xp,yp,zp)
             node(i)%stencilpln(j)%xy=(/xp,yp/)
          end do
       end do
    else
       do i=1,nodes
          x=mesh%v(i)%p(1)
          y=mesh%v(i)%p(2)
          z=mesh%v(i)%p(3) 
          call constr(x,y,z,cx,sx,cy,sy)
          node(i)%stencil(0) = node(i)%ngbr(1)%lvv(1)
          l=0
          do m=1,2
             jend=node(i)%ngbr(m)%numberngbr
             do j=1,jend 
                l=l+1
                node(i)%stencil(l) = node(i)%ngbr(m)%lvv(j+1)
                x=mesh%v(node(i)%ngbr(m)%lvv(j+1))%p(1)
                y=mesh%v(node(i)%ngbr(m)%lvv(j+1))%p(2)
                z=mesh%v(node(i)%ngbr(m)%lvv(j+1))%p(3)
                call aplyr(x,y,z,cx,sx,cy,sy,xp,yp,zp)
                node(i)%stencilpln(l)%xy=(/xp,yp/)
             end do
          end do
       end do
    end if
    return
  end subroutine stencil

  subroutine upwind_voronoi(nodes,mesh)  
    !----------------------------------------------------------------------------------------------
    !    Determinando os nos vizinhos para cada aresta 
    !----------------------------------------------------------------------------------------------
    implicit none
    integer(i4),intent(in):: nodes
    type(grid_structure),intent(in):: mesh
    integer(i4):: i
    integer(i4):: j
    integer(i4):: k
    integer(i4):: n
    integer(i4):: s
    integer(i4):: diml
    integer(i4):: jend
    real(r8):: aux_dist
    real(r8),allocatable:: dist(:)
    real(r8),allocatable:: pv(:)
    real(r8),allocatable:: pg(:)

    allocate(pv(3),pg(3),dist(6))
    dist=1.0D+10
    do i=1,nodes
       jend=node(i)%ngbr(1)%numberngbr
       diml=nint((order)/2.0D0)
       !Percorrendo todas as faces do volume de controle i
       do n=1,jend
          do s=1,1 
             !Determinando as coordenadas dos pontos de gauss da esfera para o plano 
             pg(1)=node(i)%G(s)%lpg(n,1)
             pg(2)=node(i)%G(s)%lpg(n,2)
             pg(3)=node(i)%G(s)%lpg(n,3)
             do k=1,node(i)%ngbr(1)%numberngbr
                pv=mesh%v(node(i)%ngbr(1)%lvv(k+1))%p
                dist(k)=arclen(pg,pv)
             enddo
             aux_dist=minval(dist)
             k=0
             do
                k=k+1
                if(k>6)exit
                if(aux_dist==dist(k))exit
             enddo
             dist=1.0D+10
             node(i)%upwind_voronoi(n) = node(i)%ngbr(1)%lvv(k+1)
          enddo
       enddo
    enddo
    deallocate(pv,pg,dist)
    return
  end subroutine upwind_voronoi

  subroutine upwind_donald(nodes,mesh)  
    !----------------------------------------------------------------------------------------------
    !    Determinando os nos vizinhos para cada aresta 
    !----------------------------------------------------------------------------------------------
    implicit none
    integer(i4),intent(in):: nodes
    type(grid_structure),intent(in):: mesh
    integer(i4):: i
    integer(i4):: j
    integer(i4):: k
    integer(i4):: n
    integer(i4):: s
    integer(i4):: diml
    integer(i4):: jend
    real(r8):: aux_dist
    real(r8),allocatable:: dist(:)
    real(r8),allocatable:: pv(:)
    real(r8),allocatable:: pg(:)

    allocate(pv(3),pg(3),dist(6))
    dist=1.0D+10
    do i=1,nodes
       jend=node(i)%ngbr(1)%numberngbr
       diml=nint((order)/2.0D0)
       !Percorrendo todas as faces do volume de controle i
       do n=1,2*jend
          do s=1,diml   
             !Determinando as coordenadas dos pontos de gauss da esfera para o plano 
             pg(1)=node(i)%G(s)%lpg(n,1)
             pg(2)=node(i)%G(s)%lpg(n,2)
             pg(3)=node(i)%G(s)%lpg(n,3)
             do k=1,node(i)%ngbr(1)%numberngbr
                pv=mesh%v(node(i)%ngbr(1)%lvv(k+1))%p
                dist(k)=arclen(pg,pv)
             enddo
             aux_dist=minval(dist)
             k=0
             do
                k=k+1
                if(k>6)exit
                if(aux_dist==dist(k))exit
             enddo
             dist=1.0D+10
             node(i)%G(s)%upwind_donald(n) = node(i)%ngbr(1)%lvv(k+1)
          enddo
       enddo
    enddo
    deallocate(pv,pg,dist)
    return
  end subroutine upwind_donald


  subroutine pseudoinversa(ll,cc,AA,PI)
    !----------------------------------------------------------------------------------------------
    !    Calculando a Pseudoinversa 
    !----------------------------------------------------------------------------------------------
    implicit none
    integer(i4):: ii,jj
    integer(i4):: ll,cc
    integer(i4):: m,n,i,j,k
    integer(i4):: LDA, LDU, LDVT
    integer(i4):: LWMAX
    integer(i4):: INFO, LWORK
    real(r8),intent(in)  :: AA(ll,cc)
    real(r8),intent(out) :: PI(cc,ll)
    real(r8),allocatable :: U(:,:),UT(:,:)
    real(r8),allocatable :: VT(:,:),V(:,:)
    real(r8),allocatable :: S(:),IS(:,:),SS(:,:)
    real(r8),allocatable :: WORK(:), A(:,:)
    real(r8),allocatable :: UU(:,:), VVTT(:,:)
    EXTERNAL         DGESVD
    INTRINSIC        INT, MIN
    M=ll
    N=cc 
    LDA = M
    LDU = M
    LDVT = N
    LWMAX = 1000
    allocate(UU(M,N),VVTT(N,N),SS(N,N), A(ll,cc))
    allocate(WORK(LWMAX),U(LDU,M),UT(N,M),V(N,N),VT(LDVT,N),S(N),IS(N,N))
    !Query the optimal workspace.
    LWORK = -1
    CALL DGESVD( 'All', 'All', M, N, AA, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
    LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
    !Compute SVD.
    CALL DGESVD( 'All', 'All', M, N, AA, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
    !Check for convergence.
    IF( INFO.GT.0 ) THEN
       WRITE(*,*)'The algorithm computing SVD failed to converge.'
       STOP
    END IF
    ! Construindo a matriz UU  
    do i=1,m
       do j=1,n
          UU(i,j) = U(i,j)
       end do
    end do
    ! Construindo a matriz S 
    do i = 1,n
       if (dabs(S(i)) .ge. 1.0D-15) then
          SS(i,i) = S(i)
       else 
          SS(i,i) = 0.0D0
       end if
    end do
    ! Construindo a matriz transposta de UU 
    UT = TRANSPOSE(UU)
    ! Construindo a matriz transposta de VT 
    VVTT = TRANSPOSE(VT)
    ! Construindo a matriz inversa de S   
    IS = 0.0D0    
    do i = 1,n
       if (dabs(S(i)) .ge. 1.0D-15) IS(i,i) = 1.0D0/S(i)
    end do
    ! Construindo a matriz pseudo-inversa de AA
    PI = (matmul((matmul(VVTT,IS)),UT))
    deallocate(v,u,ut,vt,s,is,work)
    return
  end subroutine pseudoinversa

  subroutine matrix_gas(nodes,mesh)  
    !----------------------------------------------------------------------------------------------
    !    Construcao do sistema A*X = B e determinando a Pseudoinversa
    !----------------------------------------------------------------------------------------------
    implicit none
    integer(i4),intent(in)    :: nodes
    type(grid_structure),intent(in)      :: mesh
    integer(i4):: i
    integer(i4):: j
    real(r8):: x
    real(r8):: y  
    !real(r8):: peso
    real(r8),allocatable:: p(:)

    do i=1,nodes
       do j=1,node(i)%ngbr(1)%numberngbr
          x=node(i)%stencilpln(j)%xy(1)
          y=node(i)%stencilpln(j)%xy(2)
          node(i)%MRG(j,1)=x
          node(i)%MRG(j,2)=y
          node(i)%MRG(j,3)=x*x
          node(i)%MRG(j,4)=x*y
          node(i)%MRG(j,5)=y*y
       enddo
       call pseudoinversa(node(i)%ngbr(1)%numberngbr,5,node(i)%MRG,node(i)%MPG)
    enddo
    return
  end subroutine matrix_gas

  subroutine vector_gas(nodes,mesh)  
    !----------------------------------------------------------------------------------------------
    !    Construcao do sistema A*X = B e determinando a Pseudoinversa
    !----------------------------------------------------------------------------------------------
    implicit none
    integer(i4),intent(in):: nodes
    type(grid_structure),intent(in):: mesh
    integer(i4):: i
    integer(i4):: j
    real(r8):: peso

    do i=1,nodes
       do j=1,node(i)%ngbr(1)%numberngbr
!          node(i)%VBG(j) = f(mesh%v(node(i)%stencil(j))%lon,mesh%v(node(i)%stencil(j))%lat) &
 !         - f(mesh%v(i)%lon,mesh%v(i)%lat)
           node(i)%VBG(j) = node(node(i)%stencil(j))%phi_new2 - node(i)%phi_new2
       end do
    end do
  end subroutine vector_gas

  subroutine reconstruction_gas(nodes,mesh)  
    !----------------------------------------------------------------------------------------------
    !    Reconstrucao da solucao
    !----------------------------------------------------------------------------------------------
    implicit none
    integer(i4),intent(in):: nodes
    type(grid_structure),intent(in):: mesh
    integer(i4):: i
    integer(i4):: j
    integer(i4):: m
    integer(i4):: n
    integer(i4):: l
    integer(i4):: c
    real(r8):: aux

    do i=1,nodes
       l=ubound(node(i)%MRG,1)
       c=ubound(node(i)%MRG,2)
!       node(i)%coef(1)=node(i)%phi_new
       node(i)%coef(1)=node(i)%phi_new2
       do m=1,c
          aux=0.0D0
          do n=1,l
             aux=aux+node(i)%MPG(m,n)*node(i)%VBG(n)
          end do
          node(i)%coef(m+1)=aux
       end do
    end do
    return 
  end subroutine reconstruction_gas

  subroutine matrix_olg(nodes,mesh)  
    !----------------------------------------------------------------------------------------------
    !    Matrix Method FV-OLG
    !----------------------------------------------------------------------------------------------
    implicit none
    integer,intent(in):: nodes
    type(grid_structure),intent(in):: mesh
    integer(i4):: i 
    integer(i4):: j 
    integer(i4):: l
    integer(i4):: c
    integer(i4):: m
    integer(i4):: n
    real(r8),allocatable:: MRG(:,:)

    if(order==2)then
       do i=1,nodes
          !Calculando os momentos para o node i 
          call moment(i,mesh)
          !Calculando os termos geometricos para o node i 
          call geometric(i,mesh)
          !Calculando phi para o node i 
          call condition_initial_olg(i,mesh)
          !Determinado o vetor B FV-OLG
          call vector_olg(i) 
          !Determinado a matriz de reconstrucao MRO FV-OLG 
          l=ubound(node(i)%geometric,1)+1
          allocate(MRG(l-1,2))
          node(i)%MRO(1,1)=node(i)%moment(1)/node(i)%moment(1)
          node(i)%MRO(1,2:)=node(i)%moment(2:)/node(i)%moment(1)
          node(i)%MRO(2:,1)=node(i)%geometric(1:,1)
          do j=2,l
             node(i)%MRO(j,2:) = node(i)%geometric(j-1,2:)*node(i)%geometric(j-1,1)
          end do
          l = ubound(node(i)%MRO,1)
          c = ubound(node(i)%MRO,2)
          !Determinando a Pseudoinversa MPI FV-OLG
          do m=2,l
             do n=2,c
                MRG(m-1,n-1)=node(i)%MRO(m,n) - node(i)%MRO(m,1)*node(i)%MRO(1,n)
             end do
          end do
          call pseudoinversa(l-1,c-1,MRG,node(i)%MPO)
          deallocate(MRG)
       enddo
    elseif(order==3)then
       do i=1,nodes
          !Calculando os momentos para o node i 
          call moment(i,mesh)
          !Calculando os termos geometricos para o node i 
          call geometric(i,mesh)
          !Calculando phi para o node i 
          call condition_initial_olg(i,mesh)
          !Determinado o vetor B FV-OLG
          call vector_olg(i) 
          !Determinado a matriz de reconstrucao MRO FV-OLG 
          l=ubound(node(i)%geometric,1)+1
          allocate(MRG(l-1,5))
          node(i)%MRO(1,1)=node(i)%moment(1)/node(i)%moment(1)
          node(i)%MRO(1,2:)=node(i)%moment(2:)/node(i)%moment(1)
          node(i)%MRO(2:,1)=node(i)%geometric(1:,1)
          do j=2,l
             node(i)%MRO(j,2:) = node(i)%geometric(j-1,2:)*node(i)%geometric(j-1,1)
          end do
          l = ubound(node(i)%MRO,1)
          c = ubound(node(i)%MRO,2)
          !Determinando a Pseudoinversa MPI FV-OLG
          do m=2,l
             do n=2,c
                MRG(m-1,n-1)=node(i)%MRO(m,n) - node(i)%MRO(m,1)*node(i)%MRO(1,n)
             end do
          end do
          call pseudoinversa(l-1,c-1,MRG,node(i)%MPO)
          deallocate(MRG)
       enddo
    else if (order == 4) then  
       do i=1,nodes
          !Calculando os momentos para o node i 
          call moment(i,mesh)
          !Calculando os termos geometricos para o node i 
          call geometric(i,mesh)
          !Calculando phi para o node i 
          call condition_initial_olg(i,mesh)         
          !Determinado o vetor B FV-OLG
          call vector_olg(i)          
          !Determinado a matriz de reconstrucao MRO FV-OLG 
          l=ubound(node(i)%geometric,1)+1
          allocate(MRG(l-1,9))
          node(i)%MRO(1,1)=node(i)%moment(1)/node(i)%moment(1)
          node(i)%MRO(1,2:)=node(i)%moment(2:)/node(i)%moment(1)
          node(i)%MRO(2:,1)=node(i)%geometric(1:,1)
          do j=2,l
             node(i)%MRO(j,2:) = node(i)%geometric(j-1,2:)*node(i)%geometric(j-1,1)
          end do
          l = ubound(node(i)%MRO,1)
          c = ubound(node(i)%MRO,2)
          !Determinando a Pseudoinversa MPI FV-OLG
          do m=2,l
             do n=2,c
                MRG(m-1,n-1)=node(i)%MRO(m,n) - node(i)%MRO(m,1)*node(i)%MRO(1,n)
             end do
          end do
          call pseudoinversa(l-1,c-1,MRG,node(i)%MPO)
          deallocate(MRG)
       enddo
    end if
    return
  end subroutine matrix_olg

  subroutine vector_olg(no)  
    !----------------------------------------------------------------------------------------------
    !    Construcao do vetor B FV-OLG
    !----------------------------------------------------------------------------------------------
    implicit none
    integer(i4),intent(in):: no
    !type(grid_structure),intent(in):: mesh
    integer(i4):: i
    integer(i4):: j
    integer(i4):: k
    integer(i4):: l
    real(r8):: m
    i=no
    if(order==2)then 
       node(i)%VBO(1)=node(i)%phi_new
       k=ubound(node(i)%geometric,1)+1
       do j=1,k-1
          l=node(i)%stencil(j)
          node(i)%VBO(j+1) = node(l)%phi_new*node(i)%geometric(j,1)
          m=node(i)%geometric(j,1)
          node(i)%VBO(j+1)=node(i)%VBO(j+1)-m*node(i)%VBO(1)
       end do
    else
       node(i)%VBO(1)=node(i)%phi_new
       k=ubound(node(i)%geometric,1)+1
       do j=1,k-1
          l=node(i)%stencil(j)
          node(i)%VBO(j+1) = node(l)%phi_new*node(i)%geometric(j,1)
          m=node(i)%geometric(j,1)
          node(i)%VBO(j+1)=node(i)%VBO(j+1)-m*node(i)%VBO(1)
       end do
    endif
    return
  end subroutine vector_olg

  subroutine vector_olg2(nodes)  
    !----------------------------------------------------------------------------------------------
    !    Construcao do vetor B FV-OLG
    !----------------------------------------------------------------------------------------------
    implicit none
    integer(i4),intent(in):: nodes
    !type(grid_structure),intent(in):: mesh
    integer(i4):: i
    integer(i4):: j
    integer(i4):: k
    integer(i4):: l,entrei
    real(r8):: m

    if(order==2)then 
       do i=1,nodes
          node(i)%VBO(1)=node(i)%phi_new2
          k=ubound(node(i)%geometric,1)+1
          do j=1,k-1
             l=node(i)%stencil(j)
             node(i)%VBO(j+1) = node(l)%phi_new2*node(i)%geometric(j,1)
             m=node(i)%geometric(j,1)
             node(i)%VBO(j+1)=node(i)%VBO(j+1)-m*node(i)%VBO(1)
          end do
       enddo
    else
       do i=1,nodes
          node(i)%VBO(1)=node(i)%phi_new2
          k=ubound(node(i)%geometric,1)+1
          do j=1,k-1
             l=node(i)%stencil(j)
             node(i)%VBO(j+1) = node(l)%phi_new2*node(i)%geometric(j,1)
             m=node(i)%geometric(j,1)
             node(i)%VBO(j+1)=node(i)%VBO(j+1)-m*node(i)%VBO(1)
          end do
       enddo
    endif
    return
  end subroutine vector_olg2

  subroutine reconstruction_olg(nodes,mesh)  
    !----------------------------------------------------------------------------------------------
    !    Reconstrucao da solucao
    !----------------------------------------------------------------------------------------------
    implicit none
    integer(i4),intent(in):: nodes
    type(grid_structure),intent(in):: mesh
    integer(i4):: i
    integer(i4):: m
    integer(i4):: n
    integer(i4):: l
    integer(i4):: c
    real(r8):: soma

    if(order==2)then
       do i=1,nodes
          l=ubound(node(i)%MRO,1)
          c=ubound(node(i)%MRO,2)
          do m=1,c-1
             soma=0.0D0
             do n=1,l-1
                soma=soma+node(i)%MPO(m,n)*node(i)%VBO(n+1)
             enddo
             node(i)%coef(m+1)=soma
          enddo
          soma=0.0D0
          do m=2,3
             soma=soma+node(i)%coef(m)*node(i)%MRO(1,m)
          enddo
          node(i)%coef(1)=node(i)%VBO(1)-soma 
       enddo
    elseif(order==3)then
       do i = 1,nodes
          !Calculando os coeficients (Pseudoinversa x B) 
          l=ubound(node(i)%MRO,1)
          c=ubound(node(i)%MRO,2)
          do m=1,c-1
             soma=0.0D0
             do n=1,l-1
                soma=soma+node(i)%MPO(m,n)*node(i)%VBO(n+1)
             end do
             node(i)%coef(m+1) = soma
          end do
          ! A restricao da media é utilizada para determinar o coeficiente c1  
          soma = 0.0D0
          do m = 2,6
             soma = soma + node(i)%coef(m)*node(i)%MRO(1,m)
          end do
          node(i)%coef(1) = node(i)%VBO(1) - soma 
       end do
    elseif(order==4)then  
       do i=1,nodes
          !Calculando os coeficients (Pseudoinversa x B) 
          l=ubound(node(i)%MRO,1)
          c=ubound(node(i)%MRO,2)
          do m=1,c-1
             soma=0.0D0
             do n=1,l-1
                soma=soma + node(i)%MPO(m,n) * node(i)%VBO(n+1)
             end do
             node(i)%coef(m+1) = soma
          end do
          !A restricao da media é utilizada para determinar o coeficiente c1  
          soma=0.0D0
          do m=2,10
             soma = soma + node(i)%coef(m)*node(i)%MRO(1,m)
          end do
          node(i)%coef(1) = node(i)%VBO(1) - soma 
       end do
    endif
    return 
  end subroutine reconstruction_olg


  subroutine coordsbarmid(nodes,mesh)  
    !----------------------------------------------------------------------------------------------
    !    Determinando as coordenadas do volume de controle de Donalds
    !----------------------------------------------------------------------------------------------
    implicit none
    integer,intent(in):: nodes
    type(grid_structure),intent(in):: mesh
    integer(i4):: i
    integer(i4):: j
    integer(i4):: k
    integer(i4):: l
    integer(i4):: jend
    real(r8)  :: cx       
    real(r8)  :: cy
    real(r8)  :: sx
    real(r8)  :: sy
    real(r8)  :: x
    real(r8)  :: y
    real(r8)  :: z
    real(r8)  :: xp
    real(r8)  :: yp
    real(r8)  :: zp

    do i=1,nodes
       x=mesh%v(i)%p(1)
       y=mesh%v(i)%p(2)
       z=mesh%v(i)%p(3)
       !Calculando as constantes da matriz de rotacao
       call constr(x,y,z,cx,sx,cy,sy)
       call aplyr(x,y,z,cx,sx,cy,sy,xp,yp,zp)
       if(order==2)jend=node(i)%ngbr(1)%numberngbr
       if(order>2)jend=node(i)%ngbr(1)%numberngbr+node(i)%ngbr(2)%numberngbr
       allocate(node(i)%Pbar(1:jend+1))
       allocate(node(i)%Pmid(1:jend+1))
       !Percorrrendo os primeiros vizinhos dos primeiros vizinhos node i, incluindo o node i 
       do j=1,jend+1
          k=node(i)%stencil(j-1) 
          jend=node(k)%ngbr(1)%numberngbr
          allocate(node(i)%Pbar(j)%xyz2(1:jend,2))
          allocate(node(i)%Pmid(j)%xyz2(1:jend,2))
          do l=1,jend
             node(k)%bar(l)%xyz=trbarycenter(mesh%v(k)%tr(l),mesh)
             call aplyr(node(k)%bar(l)%xyz(1),node(k)%bar(l)%xyz(2),node(k)%bar(l)%xyz(3),cx,sx,cy,sy,xp,yp,zp)
             node(i)%Pbar(j)%xyz2(l,:) = (/xp, yp/)
             node(k)%mid(l)%xyz=(mesh%v(k)%p + mesh%v(node(k)%stencil(l))%p) &
                  /norm((mesh%v(k)%p + mesh%v(node(k)%stencil(l))%p)) 
             call aplyr(node(k)%mid(l)%xyz(1),node(k)%mid(l)%xyz(2),node(k)%mid(l)%xyz(3),cx,sx,cy,sy,xp,yp,zp)      
             node(i)%Pmid(j)%xyz2(l,:) = (/xp, yp/)
          enddo
       enddo
    enddo
    return 
  end subroutine coordsbarmid

  subroutine coordscirc(nodes,mesh)  
    !----------------------------------------------------------------------------------------------
    !    Determinando as coordenadas do volume de controle de Donalds
    !----------------------------------------------------------------------------------------------
    implicit none
    integer,intent(in):: nodes
    type(grid_structure),intent(in):: mesh
    integer(i4):: i
    integer(i4):: j
    integer(i4):: k
    integer(i4):: l
    integer(i4):: jend
    real(r8)  :: cx       
    real(r8)  :: cy
    real(r8)  :: sx
    real(r8)  :: sy
    real(r8)  :: x
    real(r8)  :: y
    real(r8)  :: z
    real(r8)  :: xp
    real(r8)  :: yp
    real(r8)  :: zp

    do i=1,nodes
       x=mesh%v(i)%p(1)
       y=mesh%v(i)%p(2)
       z=mesh%v(i)%p(3)
       !Calculando as constantes da matriz de rotacao
       call constr(x,y,z,cx,sx,cy,sy)
       call aplyr(x,y,z,cx,sx,cy,sy,xp,yp,zp)
       if(order==2)jend=node(i)%ngbr(1)%numberngbr
       if(order>2)jend=node(i)%ngbr(1)%numberngbr+node(i)%ngbr(2)%numberngbr
       allocate(node(i)%Pcirc(1:jend+1))
       !Percorrrendo os primeiros vizinhos dos primeiros vizinhos node i, incluindo o node i 
       do j=1,jend+1
          k=node(i)%stencil(j-1) 
          jend=node(k)%ngbr(1)%numberngbr
          allocate(node(i)%Pcirc(j)%xyz2(1:jend,2))
          do l=1,jend
             node(k)%circ(l)%xyz=mesh%tr(mesh%v(k)%tr(l))%c%p
             call aplyr(node(k)%circ(l)%xyz(1),node(k)%circ(l)%xyz(2),node(k)%circ(l)%xyz(3),cx,sx,cy,sy,xp,yp,zp)
             node(i)%Pcirc(j)%xyz2(l,:) = (/xp, yp/)
          enddo
       enddo
    enddo
    return 
  end subroutine coordscirc

  subroutine area_donald(nodes,mesh)  
    !----------------------------------------------------------------------------------------------
    !    Calculando a area do volume de controle de Donalds
    !----------------------------------------------------------------------------------------------
    implicit none
    integer,intent(in):: nodes
    type(grid_structure),intent(in):: mesh
    integer(i4):: i
    integer(i4):: j
    integer(i4):: k
    integer(i4):: l
    integer(i4):: m
    integer(i4):: kk
    integer(i4):: jend
    real(r8):: atri
    real(r8),allocatable:: p1(:)
    real(r8),allocatable:: p2(:)      
    real(r8),allocatable:: p3(:)

    allocate(p1(3),p2(3),p3(3))
    do i=1,nodes
       p1=mesh%v(i)%p
       jend=node(i)%ngbr(1)%numberngbr
       l=1
       m=0
       atri=0.0D0
       do j=1,jend
          p2 = node(i)%bar(j)%xyz
          k=l-1
          do kk = 1,2
             m = m+1
             k = k+1
             if ((j == jend) .and. (kk == 2)) then
                p3 = node(i)%mid(1)%xyz
             else
                p3 = node(i)%mid(k)%xyz
             end if
             atri=atri+sphtriarea(p1,p2,p3)
          end do
          l = k
       end do
       node(i)%area=atri
    end do
    deallocate(p1,p2,p3)
    return 
  end subroutine area_donald
  
 
  subroutine gaussrootsweights(mm,x,w)  
    !----------------------------------------------------------------------------------------------
    !    Determinando as arestas do volume de controle de Donalds
    !----------------------------------------------------------------------------------------------
    implicit none
    integer,intent(in)  :: mm
    integer             :: i,j
    real(r8)            :: a,p1,p2,p3,xe
    real(r8)            :: pp,xd,xl,xm,z,z1
    real(r8),intent(out):: x(mm)
    real(r8),intent(out):: w(mm)
    real (r8)           :: eps
    eps=3.0D-14
    xe = -1.0d0
    xd =  1.0d0
    xl =  (xd-xe)*0.50d0
    xm =  (xd+xe)*0.50d0

    do i = 1,mm
       z  = dcos((dacos(-1.0D0))*(dble(i)-0.25d0)/(dble(mm)+0.25d0))
       a  = 2.0d0
       do while (a .gt. eps) 
          p1 = 1.0d0
          p2 = 0.0d0
          do j = 1,mm
             p3 = p2
             p2 = p1
             p1 = ((2.0d0*dble(j)-1.0d0)*z*p2-(dble(j)-1.0d0)*p3)/(dble(j))
          end do
          pp = (dble(mm)*(z*p1-p2))/(z*z-1.0d0)
          z1 = z
          z = z1-(p1/pp)
          a = dabs(z-z1)
       end do
       x(i)      = xm+xl*z
       x(mm+1-i) = xm-xl*z
       w(i)      = (2.0d0*xl)/((1.0d0-z*z)*pp*pp)
       w(mm+1-i) = w(i)
    enddo
    return
  end subroutine gaussrootsweights

  subroutine edges_donalds(nodes)  
    !----------------------------------------------------------------------------------------------
    !    Determinando as arestas do volume de controle de Donalds
    !----------------------------------------------------------------------------------------------
    implicit none
    integer,intent(in)    :: nodes
    integer(i4):: i
    integer(i4):: j
    integer(i4):: k
    integer(i4):: l
    integer(i4):: m
    integer(i4):: n
    integer(i4):: kk
    integer(i4):: jend
    integer(i4):: diml
    real(r8)  :: xi
    real(r8)  :: xf
    real(r8)  :: yi
    real(r8)  :: yf
    real(r8)  :: zi
    real(r8)  :: zf
    real(r8)  :: xb
    real(r8)  :: yb
    real(r8)  :: zb
    real(r8)  :: xc
    real(r8)  :: yc
    real(r8)  :: zc
    real(r8),allocatable  :: x(:) 
    real(r8),allocatable  :: w(:)   
    real(r8),allocatable  :: p1(:)
    real(r8),allocatable  :: p2(:) 

    allocate (p1(3),p2(3))
    diml = nint((order)/2.0D0)
    allocate(x(diml),w(diml))
    call gaussrootsweights(diml, x, w)

    do i=1,nodes
       jend = node(i)%ngbr(1)%numberngbr
       l=1
       m=0
       n=0
       do j=1,jend
          xb=node(i)%bar(j)%xyz(1)
          yb=node(i)%bar(j)%xyz(2)
          zb=node(i)%bar(j)%xyz(3)
          k=l-1
          do kk=1,2
             m=m+1
             k=k+1
             if((j == jend) .and. (kk == 2))then
                xc=node(i)%mid(1)%xyz(1)
                yc=node(i)%mid(1)%xyz(2)
                zc=node(i)%mid(1)%xyz(3)
             else
                xc=node(i)%mid(k)%xyz(1)
                yc=node(i)%mid(k)%xyz(2)
                zc=node(i)%mid(k)%xyz(3)
             end if
             n=n+1
             !Para n par
             if(mod(n,2)==0)then 
                xi=xb
                yi=yb
                zi=zb
                xf=xc
                yf=yc
                zf=zc              
                !Para n impar
             else
                xi=xc
                yi=yc
                zi=zc
                xf=xb
                yf=yb
                zf=zb
             end if
             p1=(/xi,yi,zi/)
             p2=(/xf,yf,zf/)
             node(i)%edge(n)%xyz2(1,:)=p1
             node(i)%edge(n)%xyz2(2,:)=p2
          end do
          l=k
       end do
    end do
    deallocate(x,w,p1,p2)
    return
  end subroutine edges_donalds


  subroutine edges_voronoi(nodes)  
    !----------------------------------------------------------------------------------------------
    !    Determinando as arestas do volume de controle de Voronoi
    !----------------------------------------------------------------------------------------------
    implicit none
    integer,intent(in)    :: nodes
    integer(i4):: i
    integer(i4):: j
    integer(i4):: k
    integer(i4):: l
    integer(i4):: m
    integer(i4):: n
    integer(i4):: jend
    real(r8)  :: xi
    real(r8)  :: xf
    real(r8)  :: yi
    real(r8)  :: yf
    real(r8)  :: zi
    real(r8)  :: zf
    real(r8),allocatable  :: p1(:)
    real(r8),allocatable  :: p2(:) 

    allocate (p1(3),p2(3))
    do i=1,nodes
       jend = node(i)%ngbr(1)%numberngbr
       do j=1,jend
          if(j==jend) then
             xi=node(i)%circ(j)%xyz(1)
             yi=node(i)%circ(j)%xyz(2)
             zi=node(i)%circ(j)%xyz(3)

             xf=node(i)%circ(1)%xyz(1)
             yf=node(i)%circ(1)%xyz(2)
             zf=node(i)%circ(1)%xyz(3)      
          else 
             xi=node(i)%circ(j)%xyz(1)
             yi=node(i)%circ(j)%xyz(2)
             zi=node(i)%circ(j)%xyz(3)

             xf=node(i)%circ(j+1)%xyz(1)
             yf=node(i)%circ(j+1)%xyz(2)
             zf=node(i)%circ(j+1)%xyz(3)          
          endif
          p1=(/xi,yi,zi/)
          p2=(/xf,yf,zf/)
          node(i)%edge(j)%xyz2(1,:)=p1
          node(i)%edge(j)%xyz2(2,:)=p2
       end do
    end do
    deallocate(p1,p2)
    return
  end subroutine edges_voronoi

  subroutine gaussedges(nodes,mesh)  
    !----------------------------------------------------------------------------------------------
    !    Determinando as raizes e pesos de gauss para cada aresta do volume de controle de Donalds
    !----------------------------------------------------------------------------------------------
    implicit none
    integer,intent(in):: nodes
    type(grid_structure),intent(in)      :: mesh       
    integer(i4):: i
    integer(i4):: VD
    integer(i4):: n
    integer(i4):: s
    integer(i4):: jend
    integer(i4):: diml
    real(r8):: alfa
    real(r8):: theta
    real(r8),allocatable:: x(:) 
    real(r8),allocatable:: w(:)   
    real(r8),allocatable:: p1(:)
    real(r8),allocatable:: p2(:) 
    real(r8),allocatable:: p3(:) 
    real(r8),allocatable:: paux(:)

    allocate (p1(3),p2(3),p3(3),paux(3))
    diml=nint((order)/2.0D0) 
    if(controlvolume=="D")then
       VD=2
    else
       VD=1
    endif
    if(method=='G') diml=1
    allocate(x(diml),w(diml))
    call gaussrootsweights(diml, x, w)
    do i=1,nodes
       jend=node(i)%ngbr(1)%numberngbr
       !Percorrendo todas as faces do volume de controle i
       do n=1,VD*jend
          p1=node(i)%edge(n)%xyz2(1,:) 
          p2=node(i)%edge(n)%xyz2(2,:) 
          do s=1,diml
             p3=(p2-p1*dot_product(p1,p2))/norm(p2-p1*dot_product(p1,p2))
             theta=dacos(dot_product(p1,p2)/(norm(p1)*norm(p2)))        
             alfa=theta*((x(s)+1)/2.0D0)  
             !Determinando as coordenadas dos pontos de gauss de cada face     
             node(i)%G(s)%lpg(n,1:3)= p1*dcos(alfa)+p3*dsin(alfa)/norm(p1*dcos(alfa)+p3*dsin(alfa))
             !Determinando do vetor unitario normal as faces 
             node(i)%G(s)%lvn(n,1:3)=-cross_product(p1,p2)/norm(cross_product(p1,p2))                   
             !Determinando o peso de gauss nas faces 
             node(i)%G(s)%lwg(n)=arclen(p1,p2)*(w(s)/2.0D0)
          end do
       end do
    end do

    deallocate(x,w,p1,p2,p3,paux)
    return
  end subroutine gaussedges


  subroutine condition_initial_gas(nodes,mesh)  
    !----------------------------------------------------------------------------------------------
    !    Determinando a condicao inical Dunavant  FV-OLG
    !----------------------------------------------------------------------------------------------
    implicit none
    integer(i4),intent(in):: nodes
    type(grid_structure),intent(in):: mesh
    integer(i4):: i

    do i=1,nodes
!       node(i)%phi_new = f(mesh%v(i)%lon,mesh%v(i)%lat)
       node(i)%phi_new2 = f(mesh%v(i)%lon,mesh%v(i)%lat)
    enddo

    return 
  end subroutine condition_initial_gas

  subroutine condition_initial_olg(no,mesh)  
    !----------------------------------------------------------------------------------------------
    !    Determinando a condicao inical Dunavant  FV-OLG
    !----------------------------------------------------------------------------------------------
    implicit none
    integer(i4),intent(in):: no
    type(grid_structure),intent(in):: mesh
    integer(i4):: i
    integer   :: ii
    integer   :: k
    integer   :: l
    integer   :: ll
    integer   :: m
    integer   :: mm
    integer   :: kk
    integer   :: jend
    integer   :: jj
    integer   :: VD
    real(r8)  :: cx
    real(r8)  :: cy
    real(r8)  :: sx
    real(r8)  :: sy
    real(r8)  :: x
    real(r8)  :: y
    real(r8)  :: z
    real(r8)  :: xp
    real(r8)  :: yp
    real(r8)  :: zp 
    real(r8)  :: xa
    real(r8)  :: ya
    real(r8)  :: za          
    real(r8)  :: xb
    real(r8)  :: yb
    real(r8)  :: xc
    real(r8)  :: yc
    real(r8)  :: xx
    real(r8)  :: yy
    real(r8)  :: zz
    real(r8)  :: soma
    real(r8)  :: atrinova, areatotal, somatotal
    real(r8)  :: tmpx
    real(r8)  :: tmpy
    real(r8)  :: w1,w2,w3
    real(r8)  :: num1,num2,num3,num4,num5
    real(r8)  :: xp1
    real(r8)  :: yp1
    real(r8)  :: zp1
    real(r8)  :: lon
    real(r8)  :: lat
    real(r8),dimension(7) :: c1 
    real(r8),dimension(7) :: c2 
    real(r8),dimension(7) :: c3 
    real(r8),dimension(7) :: wg

    w1 = 9.0D0 / 40.0D0 
    w2 = ( 155.0D0 - dsqrt(15.0D0) ) / 1200.0D0 
    w3 = ( 155.0D0 + dsqrt(15.0D0) ) / 1200.0D0 

    num1 = 1.0D0 / 3.0D0 
    num2 = ( 6.0D0 - dsqrt(15.0D0) ) / 21.0D0 
    num3 = ( 9.0D0 + 2.0D0*dsqrt(15.0D0) ) / 21.0D0 
    num4 = ( 6.0D0 + dsqrt(15.0D0) ) / 21.0D0 
    num5 = ( 9.0D0 - 2.0D0*dsqrt(15.0D0) ) / 21.0D0 
    !ATENCAO: peso do j-ésimo ponto de Gauss = wg[j]
    wg ( 1 ) = w1 
    wg ( 2 ) = w2 ; wg ( 3 ) = w2 ; wg ( 4 ) = w2 
    wg ( 5 ) = w3 ; wg ( 6 ) = w3 ; wg ( 7 ) = w3 
    !ATENCAO: j-ésimo ponto de Gauss = c1[j]*x1 + c2[j]*x2 + c3[j]*x3
    c1 ( 1 ) = num1 ; c2 ( 1 ) = num1 ; c3 ( 1 ) = num1  
    c1 ( 2 ) = num2 ; c2 ( 2 ) = num2 ; c3 ( 2 ) = num3 
    c1 ( 3 ) = num2 ; c2 ( 3 ) = num3 ; c3 ( 3 ) = num2 
    c1 ( 4 ) = num3 ; c2 ( 4 ) = num2 ; c3 ( 4 ) = num2 
    c1 ( 5 ) = num4 ; c2 ( 5 ) = num4 ; c3 ( 5 ) = num5 
    c1 ( 6 ) = num4 ; c2 ( 6 ) = num5 ; c3 ( 6 ) = num4 
    c1 ( 7 ) = num5 ; c2 ( 7 ) = num4 ; c3 ( 7 ) = num4 

    i = no

    x = mesh%v(i)%p(1)
    y = mesh%v(i)%p(2)
    z = mesh%v(i)%p(3)
    !Calculando as constantes da matriz de rotacao
    call constr(x,y,z,cx,sx,cy,sy)
    if(order==2)jend=node(i)%ngbr(1)%numberngbr
    if(order>2)jend=node(i)%ngbr(1)%numberngbr+node(i)%ngbr(2)%numberngbr
    do jj = 1,jend+1
       !Percorrrendo dos vizinhos do node i (incluindo o node i) e calculando a media da solucao exata
       ii = node(i)%stencil(jj-1) 
       jend = node(ii)%ngbr(1)%numberngbr
       xa = mesh%v(ii)%p(1)
       ya = mesh%v(ii)%p(2)
       za = mesh%v(ii)%p(3)
       call aplyr(xa,ya,za,cx,sx,cy,sy,xp,yp,zp)         
       xa = xp
       ya = yp
       za = zp
       l = 1
       m = 0
       node(ii)%phi_new=0.0D0
       areatotal = 0.0D0
       somatotal = 0.0D0
       tmpx = 0.0D0
       tmpy = 0.0D0
       if(controlvolume=="D")then
          VD=2
       else
          VD=1
       endif
       do ll = 1,VD*jend
          m = m+1
          call aplyr(node(ii)%edge(m)%xyz2(1,1),node(ii)%edge(m)%xyz2(1,2),node(ii)%edge(m)%xyz2(1,3),cx,sx,cy,sy,xp,yp,zp)    
          xc=xp
          yc=yp
          call aplyr(node(ii)%edge(m)%xyz2(2,1),node(ii)%edge(m)%xyz2(2,2),node(ii)%edge(m)%xyz2(2,3),cx,sx,cy,sy,xp,yp,zp)    
          xb=xp
          yb=yp                 
          soma = 0.0D0
          atrinova = 0.5D0*dabs(((xb-xc)*(ya-yc))-((xa-xc)*(yb-yc)))
          do mm=1,7
             xx = c1 (mm) * xa + c2 ( mm ) * xb + c3 ( mm ) * xc 
             yy = c1 (mm) * ya + c2 ( mm ) * yb + c3 ( mm ) * yc 
             !Determinando as coordenadas (x,y,z) na esfera   
             zz = dsqrt (1 - xx*xx - yy*yy)
             xp1 =  cy*xx+sy*zz
             yp1 = -sy*sx*xx + cx*yy + cy*sx*zz
             zp1 = -sy*cx*xx - sx*yy + cy*cx*zz
             !Convertendo coordenadas cartesianas em geograficas
             call cart2sph (xp1,yp1,zp1, lon, lat )
             soma = soma + atrinova * wg ( mm ) * f(lon,lat)
          end do
          node(ii)%phi_new=node(ii)%phi_new + soma
       end do
!       node(ii)%phi_new = f(mesh%v(ii)%lon,mesh%v(ii)%lat)
        node(ii)%phi_new=node(ii)%phi_new/node(ii)%moment(1)
       if (i==ii) then
          node(i)%phi_new2=node(ii)%phi_new
       end if
    end do
    return
  end subroutine condition_initial_olg


  subroutine moment(no,mesh)  
    !----------------------------------------------------------------------------------------------
    !     Momentos da Matrix Method FV-OLG
    !----------------------------------------------------------------------------------------------
    implicit none
    integer,intent(in):: no
    type(grid_structure),intent(in):: mesh 
    integer   :: i
    integer   :: l
    integer   :: k
    integer   :: m
    integer   :: n
    integer   :: ii
    integer   :: jj
    integer   :: kk
    integer   :: qq
    integer   :: ll
    integer   :: mm
    integer   :: nn
    integer   :: lll
    integer   :: VD
    integer   :: jend
    integer   :: col
    integer   :: diml
    real(r8)  :: xg
    real(r8)  :: yg
    real(r8)  :: xi
    real(r8)  :: xf
    real(r8)  :: yi
    real(r8)  :: yf
    real(r8)  :: aux1
    real(r8)  :: aux2
    real(r8)  :: aux3
    real(r8)  :: cx
    real(r8)  :: cy
    real(r8)  :: sx
    real(r8)  :: sy
    real(r8)  :: xx
    real(r8)  :: yy
    real(r8)  :: zz
    real(r8)  :: xp
    real(r8)  :: yp
    real(r8)  :: zp 
    real(r8)  :: xa
    real(r8)  :: ya
    real(r8)  :: za          
    real(r8)  :: xb
    real(r8)  :: yb
    real(r8)  :: xc
    real(r8)  :: yc
    real(r8),allocatable:: x(:) 
    real(r8),allocatable:: w(:)

    diml = nint((order+1)/2.0D0)
    allocate(x(diml),w(diml))
    call gaussrootsweights(diml, x, w)
    !Calculando o momento do node i 
    i = no
    xx = mesh%v(i)%p(1)
    yy = mesh%v(i)%p(2)
    zz = mesh%v(i)%p(3)
    call constr(xx,yy,zz,cx,sx,cy,sy)
    if(order == 2)jend = node(i)%ngbr(1)%numberngbr
    if(order >  2)jend = node(i)%ngbr(1)%numberngbr + node(i)%ngbr(2)%numberngbr
    do jj = 1,jend+1
       ii = node(i)%stencil(jj-1) 
       jend = node(ii)%ngbr(1)%numberngbr
       xa = mesh%v(ii)%p(1)
       ya = mesh%v(ii)%p(2)
       za = mesh%v(ii)%p(3)
       call aplyr(xa,ya,za,cx,sx,cy,sy,xp,yp,zp)         
       xa = xp
       ya = yp
       za = zp
       l = 1
       m = 0
       n = 0
       node(ii)%moment = 0.0D0
       if(controlvolume=="D")then
          VD=2
       else
          VD=1
       endif
       do lll = 1,VD*jend
          m = m+1
          call aplyr(node(ii)%edge(m)%xyz2(1,1),node(ii)%edge(m)%xyz2(1,2),node(ii)%edge(m)%xyz2(1,3),cx,sx,cy,sy,xp,yp,zp)    
          xi=xp
          yi=yp
          call aplyr(node(ii)%edge(m)%xyz2(2,1),node(ii)%edge(m)%xyz2(2,2),node(ii)%edge(m)%xyz2(2,3),cx,sx,cy,sy,xp,yp,zp)    
          xf=xp
          yf=yp      
          aux1 = 0.0D0 
          aux2 = 0.0D0
          aux3 = 0.0D0
          do qq = 1,diml
             xg = xi+0.50D0*(xf-xi)*(x(qq)+1.0D0)
             yg = yi+0.50D0*(yf-yi)*(x(qq)+1.0D0)
             aux1 = w(qq)*0.50D0*(yf-yi)
             aux2 = xg - xa
             aux3 = yg - ya
             col = 0
             do ll = 0,order-1
                do mm = 0,ll
                   nn = ll-mm
                   col = col+1
                   node(ii)%moment(col) = node(ii)%moment(col)+aux1*aux2**(nn+1)/(dble(nn)+ 1.0D0)*aux3**(mm)
                end do
             end do
          end do
       end do
    end do
    deallocate(x,w)
    return 
  end subroutine moment

  recursive function Factorial(n)  result(Fact)
    implicit none
    integer :: Fact
    integer,intent(in) :: n
    if (n == 0) then
       Fact = 1
    else
       Fact = n * Factorial(n-1)
    end if
  end function Factorial

  subroutine geometric(no,mesh)  
    !----------------------------------------------------------------------------------------------
    !    Termos Geometricos da Matrix Method FV-OLG
    !----------------------------------------------------------------------------------------------
    implicit none
    integer,intent(in):: no
    type(grid_structure),intent(in):: mesh 
    integer   :: i
    integer   :: j
    integer   :: l
    integer   :: m
    integer   :: r
    integer   :: s
    integer   :: nn
    integer   :: col
    integer   :: jend
    integer   :: coluna
    real(r8)  :: cx
    real(r8)  :: cy
    real(r8)  :: sx
    real(r8)  :: sy
    real(r8)  :: xx
    real(r8)  :: yy
    real(r8)  :: zz
    real(r8)  :: xp
    real(r8)  :: yp
    real(r8)  :: zp 
    real(r8)    :: bin1
    real(r8)    :: bin2
    real(r8)    :: p1
    real(r8)    :: p2
    real(r8)    :: soma
    real(r8),allocatable  :: xy(:,:)

    allocate(xy(0:order-1,0:order-1))
    !Determinando dos Termos Geometricos para Segunda, Terceira e Quarta Ordem   
    i = no
    xx = mesh%v(i)%p(1)
    yy = mesh%v(i)%p(2)
    zz = mesh%v(i)%p(3)
    call constr(xx,yy,zz,cx,sx,cy,sy)
    call aplyr(xx,yy,zz,cx,sx,cy,sy,xp,yp,zp)         
    xx = xp
    yy = yp
    if(order == 2) jend = node(i)%ngbr(1)%numberngbr
    if(order > 2)  jend = node(i)%ngbr(1)%numberngbr+node(i)%ngbr(2)%numberngbr
    do j = 1,jend
       col = 0
       do l = 0,order-1
          do m = 0,l
             nn = l-m
             col = col+1
             xy(nn,m) = (node(node(i)%stencil(j))%moment(col))/(node(node(i)%stencil(j))%moment(1))
          end do
       end do
       coluna = 1
       do l = 1,order-1
          do m = 0,l
             nn = l-m
             soma  = 0.0D0
             do r = 0,m
                do s = 0,nn
                   bin1 = dble(Factorial(m))/(dble(Factorial(m-r))*dble(Factorial(r)))
                   bin2 = dble(Factorial(nn))/(dble(Factorial(nn-s))*dble(Factorial(s)))
                   p1 = (node(i)%stencilpln(j)%xy(1) - xx)**s
                   p2 = (node(i)%stencilpln(j)%xy(2) - yy)**r
                   soma = soma+bin1*bin2*p1*p2*xy(nn-s,m-r)
                end do
             end do
             coluna = coluna+1
             node(i)%geometric(j,coluna) = soma
          end do
       end do
       node(i)%geometric(j,1)=1.0D0/dsqrt(node(i)%stencilpln(j)%xy(1)**2 + node(i)%stencilpln(j)%xy(2)**2)
    end do
    deallocate(xy)
    return 
  end subroutine geometric



  !======================================================================================
  !    ROTINES OF TESTS
  !======================================================================================
  subroutine interpolation(nodes,mesh)  
    implicit none
    type(grid_structure),intent(in)      :: mesh
    integer(i4),intent(in):: nodes
    integer(i4) :: ii
    integer(i4) :: jj
    integer (i4):: nlon !Total number of longitudes
    integer (i4):: nlat !Total number of latitudes
    !Auxiliar variables
    real (r8):: pp(1:3)  !General point in R3
    real (r8):: dlat    !Latitude step
    real (r8):: dlon    !Longitude step
    real (r8):: tlat    !Latitude
    real (r8):: tlon    !Longitude
    real(r8)  :: cx
    real(r8)  :: cy
    real(r8)  :: sx
    real(r8)  :: sy
    real(r8)  :: x
    real(r8)  :: y
    real(r8)  :: z
    real(r8)  :: xp
    real(r8)  :: yp
    real(r8)  :: zp
    integer(i4):: pol_ref
    integer(i4):: contador
    real(r8)  :: lon, lat
    real(r8)  :: erro 
    real(r8)  :: erro_Linf, erro_L2
    character (len=256):: filename
    real(r8),allocatable    :: p(:)
    real(r8),allocatable    :: fexact(:,:)
    real(r8),allocatable    :: fest(:,:)
    real(r8),allocatable    :: ferror(:,:)

    erro_Linf = 0.0D0
    erro_L2 = 0.0D0
    erro = 0.0D0
    contador = 0
    !Lat-lon grid sizes
    nlat=720 !180 !720
    nlon=1440 !360 !1440
    dlat=180._r8/dble(nlat)!real(nlat, r8)
    dlon=360._r8/dble(nlon)!real(nlon, r8)
    allocate (p(3))
    allocate(fexact(1:nlon,1:nlat))
    allocate(fest(1:nlon,1:nlat))
    allocate(ferror(1:nlon,1:nlat))

    !Pixel registration mode (for GMT graph)
    tlat=-90._r8+dlat/2._r8
    do jj=1,nlat
       tlon=-180._r8+dlon/2._r8
       do ii=1,nlon

          call sph2cart(tlon*deg2rad,tlat*deg2rad, pp(1), pp(2), pp(3))
          call cart2sph(pp(1), pp(2), pp(3),lon, lat)
          p = (/pp(1), pp(2), pp(3)/)

          !print*, lon,lat,f(lon,lat) 
          !read(*,*)           
          ! Determinando o polinomio interpolador mais proximo do ponto p 
          pol_ref = getnearnode(p, mesh)

          ! Determinando o valor da funcao f e o valor do polinomio interpolador    
          x = mesh%v(pol_ref)%p(1)
          y = mesh%v(pol_ref)%p(2)
          z = mesh%v(pol_ref)%p(3)
          call constr(x,y,z,cx,sx,cy,sy)
          call aplyr(pp(1), pp(2), pp(3),cx,sx,cy,sy,xp,yp,zp) 

          if(order==2)then         
             fest(ii,jj) = (node(pol_ref)%coef(1) + node(pol_ref)%coef(2)*xp + node(pol_ref)%coef(3)*yp)
          elseif(order==3)then
             fest(ii,jj) = node(pol_ref)%coef(1) + node(pol_ref)%coef(2)*xp + node(pol_ref)%coef(3)*yp + &
                  node(pol_ref)%coef(4)*xp*xp + node(pol_ref)%coef(5)*xp*yp + node(pol_ref)%coef(6)*yp*yp
          elseif(order==4)then  
             fest(ii,jj) = node(pol_ref)%coef(1) + node(pol_ref)%coef(2)*xp + node(pol_ref)%coef(3)*yp + &
                  node(pol_ref)%coef(4)*xp*xp + node(pol_ref)%coef(5)*xp*yp + node(pol_ref)%coef(6)*yp*yp + &                     
                  node(pol_ref)%coef(7)*xp*xp*xp + node(pol_ref)%coef(8)*xp*xp*yp + &
                  node(pol_ref)%coef(9)*xp*yp*yp + node(pol_ref)%coef(10)*yp*yp*yp
          end if

          fexact(ii,jj)=f(lon,lat) 
          ferror(ii,jj)=fest(ii,jj) - fexact(ii,jj)
          erro = ferror(ii,jj)
          erro_L2 = erro_L2 + erro*erro 
          contador = contador + 1
          erro = ferror(ii,jj)! f(lon,lat) - (node(pol_ref)%coef(1) + node(pol_ref)%coef(2)*xp + node(pol_ref)%coef(3)*yp)
          if (erro_Linf < abs (erro))  erro_Linf =  abs ( erro )        
          tlon=tlon+dlon
       end do
       tlat=tlat+dlat
    end do

    !    ! Plotando a funcao exata 
    !    filename = "cossine_exata_ordem4"
    !    call plot_scalarfield_sphgrid(nlon, nlat, fexact,filename)
    ! 
    !    ! Plotando o erro 
    !    filename = "cossine_erro_ordem4"
    !    call plot_scalarfield_sphgrid(nlon, nlat, ferror,filename)
    ! 
    !   ! Plotando a funcao estimada 
    !    filename = "cossine_estimada_ordem4"
    !    call plot_scalarfield_sphgrid(nlon, nlat, fest,filename)
    print*,erro_Linf,'erro_Linf'
    print*, dsqrt(erro_L2/contador), 'erro_L2'
    return
  end subroutine interpolation


  subroutine flux_edges_olg(nodes,mesh,z,time)  
    !----------------------------------------------------------------------------------------------
    !    Calculando o fluxo nas arestas
    !----------------------------------------------------------------------------------------------
    implicit none 
    integer,intent(in)    :: nodes
    type(grid_structure),intent(in)      :: mesh
    integer,intent(in),optional   :: z
    real(r8),intent(in),optional  :: time
    integer(i4)   :: i
    integer(i4)   :: j
    integer(i4)   :: k
    integer(i4)   :: l
    integer(i4)   :: m
    integer(i4)   :: n
    integer(i4)   :: w
    integer(i4)   :: s
    integer(i4)   :: kk
    integer(i4)   :: cc
    integer(i4)   :: diml
    integer(i4)   :: contador
    integer(i4)   :: jend
    real(r8)  :: dot
    real(r8)  :: cx
    real(r8)  :: cy
    real(r8)  :: sx
    real(r8)  :: sy
    real(r8)  :: xx
    real(r8)  :: yy
    real(r8)  :: zz
    real(r8)  :: xp
    real(r8)  :: yp
    real(r8)  :: zp
    real(r8)  :: xb
    real(r8)  :: yb
    real(r8)  :: zb
    real(r8)  :: xc
    real(r8)  :: yc
    real(r8)  :: zc
    real(r8)  :: xi
    real(r8)  :: yi
    real(r8)  :: zi
    real(r8)  :: xf
    real(r8)  :: yf
    real(r8)  :: zf
    real(r8)  :: FEPS, temp, aaxx, erro_Linf, erro, integral, lon, lat
    real(r8),allocatable  :: p1(:)      
    real(r8),allocatable  :: p2(:)      
    real(r8),allocatable  :: p3(:)      
    real(r8),allocatable  :: vn(:)      

    allocate (p1(3),p2(3),p3(3),vn(3))
    FEPS = 1.0D-8
    dot = 0.0D0
    temp = 0.0D0
    aaxx = 0.0D0
    erro_Linf = 0.0D0
    contador = 0
    do i = 1,nodes
       jend = node(i)%ngbr(1)%numberngbr
       l = 1
       m = 0
       n = 0
       diml = nint((order)/2.0D0)
       if(controlvolume=="D")then
          cc=2
       else
          cc=1
       endif
       if(method=='G') diml=1
       do n = 1,cc*jend     
          node(i)%S(z)%flux = 0.0D0
          ! Percorrendo todas as faces do volume de controle i
          p1 = node(i)%edge(n)%xyz2(1,:) 
          p2 = node(i)%edge(n)%xyz2(2,:) 
          do s=1,diml
             contador = contador + 1
             ! Calculando o produto interno entre o vetor velocidade e o vetor normal unitario a edge 
             dot = dot_product (velocity(node(i)%G(s)%lpg(n,1:3), time),node(i)%G(s)%lvn(n,1:3))
             ! Determinando os coeficientes da matriz de rotação para o node i 
             xx = mesh%v(i)%p(1)
             yy = mesh%v(i)%p(2)
             zz = mesh%v(i)%p(3)
             call constr(xx,yy,zz,cx,sx,cy,sy)
             ! Determinando as coordenadas dos pontos de gauss da esfera para o plano 
             p3(1) = node(i)%G(s)%lpg(n,1)
             p3(2) = node(i)%G(s)%lpg(n,2)
             p3(3) = node(i)%G(s)%lpg(n,3)
             call aplyr(p3(1), p3(2), p3(3),cx,sx,cy,sy,xp,yp,zp) 
             call cart2sph (xp,yp,zp,lon, lat)
             integral = f(lon,lat)*node(i)%G(s)%lwg(n)*dot
             ! Calculando o fluxo na esfera (Utilizando o polinomio interpolador no plano para estimar o valor na superficie esferica)
             if (order == 2) then           
                node(i)%S(z)%flux = (node(i)%coef(1) + node(i)%coef(2)*xp + & 
                     node(i)%coef(3)*yp)*node(i)%G(s)%lwg(n)*dot
             end if
             if (order == 3) then    
                node(i)%S(z)%flux =  (node(i)%coef(1) + node(i)%coef(2)*xp + node(i)%coef(3)*yp + &
                     node(i)%coef(4)*xp*xp + node(i)%coef(5)*xp*yp + node(i)%coef(6)*yp*yp)*node(i)%G(s)%lwg(n)*dot       
             end if
             if (order == 4) then     
                node(i)%S(z)%flux = (node(i)%coef(1) + node(i)%coef(2)*xp + node(i)%coef(3)*yp + &
                     node(i)%coef(4)*xp*xp + node(i)%coef(5)*xp*yp + node(i)%coef(6)*yp*yp + &                     
                     node(i)%coef(7)*xp*xp*xp + node(i)%coef(8)*xp*xp*yp + &
                     node(i)%coef(9)*xp*yp*yp + node(i)%coef(10)*yp*yp*yp)*node(i)%G(s)%lwg(n)*dot   
             end if
             !   vn = node(i)%G(s)%lvn(n,1:3)
          end do
          !call numerical_integration(p1,p2,vn,integral)
          erro = dabs(node(i)%S(z)%flux - integral)  
          aaxx = erro * erro
          if (erro_Linf < abs (erro))  erro_Linf =  abs ( erro )
          temp = temp +  aaxx
       end do
    end do
    print*, erro_Linf, 'erro_Linf'
    print*, dsqrt (temp/contador) , 'erro_L2'
    deallocate(p1,p2,p3,vn)
    return  
  end subroutine flux_edges_olg

  subroutine flux_edges_gas(nodes,mesh,z,time)  
    !----------------------------------------------------------------------------------------------
    !    Calculando o fluxo nas arestas
    !----------------------------------------------------------------------------------------------
    implicit none 
    integer,intent(in)    :: nodes
    type(grid_structure),intent(in)      :: mesh
    integer,intent(in),optional   :: z
    real(r8),intent(in),optional  :: time
    integer(i4)   :: i,ii
    integer(i4)   :: j,jj
    integer(i4)   :: k
    integer(i4)   :: l
    integer(i4)   :: m
    integer(i4)   :: n
    integer(i4)   :: w
    integer(i4)   :: s
    integer(i4)   :: kk
    integer(i4)   :: cc
    integer(i4)   :: diml
    integer(i4)   :: contador
    integer(i4)   :: jend
    real(r8)  :: dot
    real(r8)  :: cx
    real(r8)  :: cy
    real(r8)  :: sx
    real(r8)  :: sy
    real(r8)  :: xx
    real(r8)  :: yy
    real(r8)  :: zz
    real(r8)  :: xp
    real(r8)  :: yp
    real(r8)  :: zp
    real(r8)  :: xb
    real(r8)  :: yb
    real(r8)  :: zb
    real(r8)  :: xc
    real(r8)  :: yc
    real(r8)  :: zc
    real(r8)  :: xi
    real(r8)  :: yi
    real(r8)  :: zi
    real(r8)  :: xf
    real(r8)  :: yf
    real(r8)  :: zf, theta, alfa, div_est, sol_numerica, sol_exata, derivada_i, derivada_j, cosseno, seno, aux
    real(r8)  :: FEPS, temp, aaxx, erro_Linf, erro, integral,phi_i, phi_j,dist,sinal, aux1, aux2, aux3, lon, lat
    real(r8)  :: lat1, lat2, flux_numerico, flux_exato, lon1
    real(r8),allocatable  :: p1(:)      
    real(r8),allocatable  :: p2(:)      
    real(r8),allocatable  :: p3(:)      
    real(r8),allocatable  :: vn(:)      
    real(r8),allocatable  :: pm(:)      


    allocate (p1(3),p2(3),p3(3),vn(3),pm(3))
    FEPS = 1.0D-8
    dot = 0.0D0
    temp = 0.0D0
    flux_numerico = 0.0D0
    flux_exato = 0.0D0
    aaxx = 0.0D0
    erro_Linf = 0.0D0
    contador = 0  
    !stop
    do i = 132,132
         node(i)%S(z)%flux = 0.0D0
         jend = node(i)%ngbr(1)%numberngbr
         p3 =mesh%v(i)%p 
         call cart2sph (p3(1), p3(2), p3(3),lon,lat)
         lon1=lon
!         print*, lon, lat
!           if(dabs(lat)>=0.0D0 .and. dabs(lat)<=1.0D-15) then
!           print*, i,lon,lat
!           read(*,*)

        if(order >  2)jend = node(i)%ngbr(1)%numberngbr + node(i)%ngbr(2)%numberngbr
        do jj = 1,jend
           ii = node(i)%stencil(jj)          
           p3 =mesh%v(ii)%p 
           call cart2sph (p3(1), p3(2), p3(3),lon,lat)
           print*,ii,lat
           if(dabs(lat)>=0.0D0 .and. dabs(lat)<=1.0D-14) then
           print*, i,ii,lon,lat
           read(*,*)
 
           ! calcular aqui o fluxo....
           
          !Determinando os valores para o node i 
           xx = mesh%v(i)%p(1)
           yy = mesh%v(i)%p(2)
           zz = mesh%v(i)%p(3)
           call constr(xx,yy,zz,cx,sx,cy,sy)

           xx = mesh%v(ii)%p(1)
           yy = mesh%v(ii)%p(2)
           zz = mesh%v(ii)%p(3)
           call aplyr(xx,yy,zz,cx,sx,cy,sy,xp,yp,zp) 
           cosseno=(xp/dsqrt(xp**2+yp**2))
           seno=(yp/dsqrt(xp**2+yp**2))
           phi_i = node(i)%coef(1)
           derivada_i = 2.0D0*node(i)%coef(4)*(cosseno**2) + & 
           2.0D0*node(i)%coef(5)*cosseno*seno + 2.0D0*node(i)%coef(6)*(seno**2)
 
          !Determinando os valores dos vizinhos do node i 
           xx = mesh%v(ii)%p(1)
           yy = mesh%v(ii)%p(2)
           zz = mesh%v(ii)%p(3)
           call constr(xx,yy,zz,cx,sx,cy,sy)
           xx = mesh%v(i)%p(1)
           yy = mesh%v(i)%p(2)
           zz = mesh%v(i)%p(3)
           call aplyr(xx,yy,zz,cx,sx,cy,sy,xp,yp,zp) 
           cosseno=(xp/dsqrt(xp**2+yp**2))
           seno=(yp/dsqrt(xp**2+yp**2))
           phi_j = node(ii)%coef(1)
           derivada_j = 2.0D0*node(ii)%coef(4)*(cosseno**2) + &
           2.0D0*node(ii)%coef(5)*cosseno*seno + 2.0D0*node(ii)%coef(6)*(seno**2)           

           !Distancia entre o node i e seu respectivo vizinho   
           dist=arclen(mesh%v(i)%p,mesh%v(ii)%p)      

                      
           aux1 = (1.0D0/2.0D0)*(phi_i + phi_j)
           aux2 = (1.0D0/12.0D0)*((dist)**2)*(derivada_j + derivada_i)
!           aux3 = sinal*(1.0D0/48.0D0)*((dist)**2)*(derivada_j - derivada_i)
           
           sol_numerica = aux1 - aux2 !+ aux3
!           print*, sol_numerica , 'solucao numerica'

           pm=(mesh%v(i)%p + mesh%v(ii)%p)/2.0D0
           pm=pm/norm(pm)
           
           
           print*, arclen(mesh%v(i)%p,pm) , 'DISTANCIA'
           
           call cart2sph (pm(1),pm(2),pm(3),lon, lat)
           sol_exata = f(mesh%v(i)%lon,mesh%v(i)%lat)
           
           print*, sol_numerica,'sol_numerica antes'
           print*, sol_exata,'sol_exata antes'
                      
           if(lon>lon1) then 
           sol_numerica = -sol_numerica
           sol_exata = -sol_exata
           end if 

           print*, sol_numerica,'sol_numerica depois'
           print*, sol_exata,'sol_exata depois'
                   
           flux_numerico = flux_numerico +  sol_numerica  
           flux_exato = flux_exato + sol_exata
         end if     
         enddo


 enddo    

print*, flux_numerico,'fluxo numerico'
print*, flux_exato, 'fluxo exato' 
print*, dabs(flux_numerico-flux_exato), 'ERRO'
!read(*,*)


    deallocate(p1,p2,p3,vn,pm)
    return  
  end subroutine flux_edges_gas

  subroutine numerical_integration(p1,p2,vn,integral)  
    !----------------------------------------------------------------------------------------------
    !    Calculando a integral numerica utilizando regra de simpson composta 
    !----------------------------------------------------------------------------------------------
    implicit none
    real(r8),dimension (3),intent(in) :: p1
    real(r8),dimension (3),intent(in) :: p2
    real(r8),dimension (3),intent(in) :: vn
    real(r8),intent(out)              :: integral
    integer   :: k
    integer   :: n
    real(r8)  :: aa
    real(r8)  :: bb
    real(r8)  :: hh
    real(r8)  :: alfa,aux1
    real(r8)  :: somapar
    real(r8)  :: somaimpar
    real(r8)  :: lon, lat
    real(r8)  :: fa,fb
    real(r8)  :: time,dot
    real(r8),allocatable  :: a(:) 
    real(r8),allocatable  :: b(:) 
    real(r8),allocatable  :: p3(:) 
    real(r8),allocatable  :: auxiliar(:) 
    real(r8),allocatable  :: pontoauxiliar1(:) 
    real(r8),allocatable  :: pontoauxiliar2(:) 
    allocate (a(3),b(3),auxiliar(3),pontoauxiliar1(3),pontoauxiliar2(3),p3(3))
    a = p1
    b = p2
    aa = 0.0D0
    bb = 1.0D0
    time = 0.0D0
    !n quantidade de intervalos (numero par)
    n = 2000
    !calculando o angulo em radianos entre p1 e p2
    alfa = dacos(dot_product(p1,p2)/(norm(p1)*norm(p2)))
    !Calculando o tamanho de h  
    hh = (alfa)/dble(2*n)
    p3 = (p2 - p1*dot_product(p1,p2)) / norm(p2 - p1*dot_product(p1,p2))
    !Calculando os termos pares     
    somapar = 0.0D0
    do k=1, n-1
       aux1 = (2.0d0*dble(k))*hh
       auxiliar  = p1*dcos(aux1) + p3*dsin(aux1)
       dot = dot_product (velocity(auxiliar, time),vn)       
       call cart2sph (auxiliar(1), auxiliar(2),auxiliar(3), lon, lat )
       somapar = somapar + 2.0D0*f(lon,lat)!*dot  
    end do
    !Calculando os termos impares     
    somaimpar = 0.0D0
    do k=1, n
       aux1 = (2.0d0*dble(k)-1)*hh
       auxiliar  = p1*dcos(aux1) + p3*dsin(aux1)      
       dot = dot_product (velocity(auxiliar, time),vn) 
       call cart2sph (auxiliar(1), auxiliar(2),auxiliar(3), lon, lat )
       somaimpar = somaimpar + 4.0D0*f(lon,lat)!*dot   
    end do
    call cart2sph (a(1), a(2),a(3), lon, lat )
    dot = dot_product (velocity(a, time),vn)   
    fa = f(lon,lat)!*dot
    call cart2sph (b(1), b(2),b(3), lon, lat )
    dot = dot_product (velocity(b, time),vn)   
    fb = f(lon,lat)!*dot
    integral = (fa + fb + somapar + somaimpar)*hh/3.0D0
    deallocate(a,b,auxiliar,pontoauxiliar1,pontoauxiliar2,p3)
    return 
  end subroutine numerical_integration


  subroutine divergente_olg(nodes,mesh,z,time)
    implicit none
    integer,intent(in)    :: nodes
    type(grid_structure),intent(in)      :: mesh
    integer,intent(in),optional   :: z
    real(r8),intent(in),optional  :: time
    integer(i4)   :: i
    integer(i4)   :: j
    integer(i4)   :: k
    integer(i4)   :: l
    integer(i4)   :: m
    integer(i4)   :: n
    integer(i4)   :: w
    integer(i4)   :: s
    integer(i4)   :: kk
    integer(i4)   :: cc
    integer(i4)   :: diml
    integer(i4)   :: contador
    integer(i4)   :: jend
    real(r8)  :: dot
    real(r8)  :: cx
    real(r8)  :: cy
    real(r8)  :: sx
    real(r8)  :: sy
    real(r8)  :: xx
    real(r8)  :: yy
    real(r8)  :: zz
    real(r8)  :: xp
    real(r8)  :: yp
    real(r8)  :: zp
    real(r8)  :: xb
    real(r8)  :: yb
    real(r8)  :: zb
    real(r8)  :: xc
    real(r8)  :: yc
    real(r8)  :: zc
    real(r8)  :: xi
    real(r8)  :: yi
    real(r8)  :: zi
    real(r8)  :: xf
    real(r8)  :: yf
    real(r8)  :: zf
    real(r8)  :: FEPS, temp, aaxx, erro_Linf, erro, integral, div_est
    real(r8),allocatable  :: p1(:)      
    real(r8),allocatable  :: p2(:)      
    real(r8),allocatable  :: p3(:)      
    real(r8),allocatable  :: vn(:)      

    allocate (p1(3),p2(3),p3(3),vn(3))
    FEPS = 1.0D-8
    dot = 0.0D0
    temp = 0.0D0
    aaxx = 0.0D0
    erro_Linf = 0.0D0
    contador = 0
    do i = 1,nodes
       node(i)%S(z)%flux = 0.0D0
       jend = node(i)%ngbr(1)%numberngbr
       l = 1
       m = 0
       n = 0
       diml = nint((order)/2.0D0)
       if(controlvolume=="D")then
          cc=2
       else
          cc=1
       endif
       if(method=='G') diml=1
       do n = 1,cc*jend     

          ! Percorrendo todas as faces do volume de controle i
          p1 = node(i)%edge(n)%xyz2(1,:) 
          p2 = node(i)%edge(n)%xyz2(2,:) 
          do s=1,diml
             contador = contador + 1
             ! Calculando o produto interno entre o vetor velocidade e o vetor normal unitario a edge 
             dot = dot_product (velocity(node(i)%G(s)%lpg(n,1:3), time),node(i)%G(s)%lvn(n,1:3))
             ! Determinando os coeficientes da matriz de rotação para o node i 
             xx = mesh%v(i)%p(1)
             yy = mesh%v(i)%p(2)
             zz = mesh%v(i)%p(3)
             call constr(xx,yy,zz,cx,sx,cy,sy)
             ! Determinando as coordenadas dos pontos de gauss da esfera para o plano 
             p3(1) = node(i)%G(s)%lpg(n,1)
             p3(2) = node(i)%G(s)%lpg(n,2)
             p3(3) = node(i)%G(s)%lpg(n,3)
             call aplyr(p3(1), p3(2), p3(3),cx,sx,cy,sy,xp,yp,zp) 
             ! Calculando o fluxo na esfera (Utilizando o polinomio interpolador no plano para estimar o valor na superficie esferica)
             if (order == 2) then           
                node(i)%S(z)%flux = node(i)%S(z)%flux + (node(i)%coef(1) + node(i)%coef(2)*xp + & 
                     node(i)%coef(3)*yp)*node(i)%G(s)%lwg(n)*dot
             end if
             if (order == 3) then    
                node(i)%S(z)%flux = node(i)%S(z)%flux + (node(i)%coef(1) + node(i)%coef(2)*xp + node(i)%coef(3)*yp + &
                     node(i)%coef(4)*xp*xp + node(i)%coef(5)*xp*yp + node(i)%coef(6)*yp*yp)*node(i)%G(s)%lwg(n)*dot       
             end if
             if (order == 4) then     
                node(i)%S(z)%flux = node(i)%S(z)%flux  + (node(i)%coef(1) + node(i)%coef(2)*xp + node(i)%coef(3)*yp + &
                     node(i)%coef(4)*xp*xp + node(i)%coef(5)*xp*yp + node(i)%coef(6)*yp*yp + &                     
                     node(i)%coef(7)*xp*xp*xp + node(i)%coef(8)*xp*xp*yp + &
                     node(i)%coef(9)*xp*yp*yp + node(i)%coef(10)*yp*yp*yp)*node(i)%G(s)%lwg(n)*dot   
             end if
          end do
       end do
    end do

    do i = 1,nodes
       div_est = node(i)%S(z)%flux/node(i)%area
       erro = dabs(div_est)
       aaxx = erro * erro
       if (erro_Linf < abs (erro))  erro_Linf =  abs ( erro )
       temp = temp +  aaxx
    end do
    print*, erro_Linf, 'norma 1'
    print*, dsqrt (temp/nodes) , 'norma 2'
    deallocate(p1,p2,p3,vn)
    return
  end subroutine divergente_olg



  subroutine time_step(nodes,mesh,cfl)  
    !----------------------------------------------------------------------------------------------
    !    Calculando o passo de tempo
    !----------------------------------------------------------------------------------------------
    implicit none

    integer,intent(in):: nodes
    type(grid_structure),intent(in):: mesh
    real(r8),intent(in),optional:: cfl

    integer(i4):: i
    integer(i4):: aux
    real(r8):: time
    real(r8):: dx
    real(r8):: FEPS
    real(r8):: norma
    real(r8),allocatable:: time_local(:)
    real(r8),allocatable:: velocidade_total(:)
    real(r8),allocatable:: vel(:)
    FEPS = 1.0D-8
    time = 0.0D0  
    allocate (time_local(1:nodes))
    allocate (velocidade_total(1:nodes))
    time_local = 0.0D0
    allocate (vel(3))

    do i = 1,nodes
       vel = velocity(mesh%v(i)%p, time)
       norma = norm(vel)       
       if ( norma < FEPS ) norma = FEPS  

       if(controlvolume=="D") then    
          dx = sqrt (node(i)%area/node(i)%ngbr(1)%numberngbr) 
       else 
          dx = sqrt (mesh%hx(i)%areag/node(i)%ngbr(1)%numberngbr) 
          
       endif
       time_local(i)= cfl*dx/norma 
    end do
   
    node(0)%dt  =  minval (time_local)!/node(0)%div
    aux = int(T/node(0)%dt)

    if (mod(aux,5)==0) then
       node(0)%dt = minval (time_local)
    else  
       do while (mod(aux,5)/=0) 
          aux = aux + 1 
       end do
       node(0)%dt = T/aux
    end if
    node(0)%ntime = int(aux)
    node(0)%dt = T/node(0)%ntime
    deallocate (time_local,vel)     
    return
  end subroutine time_step

  subroutine flux_olg(nodes,mesh,z,time,uedges)  
    !----------------------------------------------------------------------------------------------
    !    Calculando o fluxo nas arestas
    !----------------------------------------------------------------------------------------------
    implicit none
    integer(i4),intent(in):: nodes
    type(grid_structure),intent(inout):: mesh
    integer,intent(in),optional   :: z
    real(r8),intent(in),optional  :: time
    type(scalar_field), intent(inout), optional :: uedges ! Only for moist shallow water model
    integer(i4):: i
    integer(i4):: j
    integer(i4):: k
    integer(i4):: w
    integer(i4):: l
    integer(i4):: n
    integer(i4):: e
    integer(i4):: i1
    integer(i4):: i2
    integer(i4):: s
    integer(i4):: cc
    integer(i4):: diml
    integer(i4):: jend
    real(r8):: dot
    real(r8):: xx
    real(r8):: yy
    real(r8):: zz
    real(r8):: xp
    real(r8):: yp
    real(r8):: zp
    real(r8):: cx
    real(r8):: cy
    real(r8):: sx
    real(r8):: sy
    real(r8):: pol, pol_filtered
    real(r8):: minval_i , maxval_i
    real(r8):: minval_i1, maxval_i1
    real(r8):: minval_i2, maxval_i2
    real(r8):: FEPS
    real(r8),allocatable:: p1(:),p2(:),p3(:)
    allocate(p1(3),p2(3),p3(3))
    FEPS = 1.0D-8


    do i=1,nodes
       node(i)%S(z)%flux=0.0D0
       jend=node(i)%ngbr(1)%numberngbr
       diml=nint((order)/2.0D0)
       !Percorrendo todas as faces do volume de controle i
       if(controlvolume=="D")then
          cc=2
       else
          cc=1
       endif

       ! Flux at edges
       node(i)%edge_flux(:) = 0._r8

       do n=1,cc*jend
          p1=node(i)%edge(n)%xyz2(1,:) 
          p2=node(i)%edge(n)%xyz2(2,:) 
          do s=1,diml
             !Calculando o produto interno entre o vetor velocidade e o vetor normal unitario a edge 
             if (.not. moistswm)then
               dot=dot_product(velocity(node(i)%G(s)%lpg(n,1:3), time), node(i)%G(s)%lvn(n,1:3))
             else
               dot=dot_product(node(i)%G(s)%velocity_quadrature(n)%v, node(i)%G(s)%lvn(n,1:3))
             end if
             !Determinando os coeficientes da matriz de rotação para o node i 
             xx=mesh%v(i)%p(1)
             yy=mesh%v(i)%p(2)
             zz=mesh%v(i)%p(3)
             call constr(xx,yy,zz,cx,sx,cy,sy)
             !Determinando as coordenadas dos pontos de gauss da esfera para o plano 
             p3(1)=node(i)%G(s)%lpg(n,1)
             p3(2)=node(i)%G(s)%lpg(n,2)
             p3(3)=node(i)%G(s)%lpg(n,3)
             call aplyr(p3(1),p3(2),p3(3),cx,sx,cy,sy,xp,yp,zp)
             if(dot<FEPS.and.dot>-1.0*FEPS)then
              node(i)%S(z)%flux=node(i)%S(z)%flux+0.0D0
              elseif(dot>0)then
                if(order==2)then
                   ! Polynomial evaluated at quadrature point
                   pol = node(i)%coef(1) + node(i)%coef(2)*xp + & 
                         node(i)%coef(3)*yp

                else if(order==3)then
                   ! Polynomial evaluated at quadrature point
                   pol = node(i)%coef(1) + node(i)%coef(2)*xp + node(i)%coef(3)*yp + &
                        node(i)%coef(4)*xp*xp + node(i)%coef(5)*xp*yp + node(i)%coef(6)*yp*yp

                else if(order==4) then     
                   ! Polynomial evaluated at quadrature point
                   pol = node(i)%coef(1) + node(i)%coef(2)*xp + node(i)%coef(3)*yp + &
                        node(i)%coef(4)*xp*xp + node(i)%coef(5)*xp*yp + node(i)%coef(6)*yp*yp + &                     
                        node(i)%coef(7)*xp*xp*xp + node(i)%coef(8)*xp*xp*yp + &
                        node(i)%coef(9)*xp*yp*yp + node(i)%coef(10)*yp*yp*yp
                end if

                ! Correct the value if monotonic limiter is applied
                if (monotonicfilter .and. time_integrator=='rk4')then
                  !if(moistswm)then ! To ensure mass conservation
                    e = edges_indexes(i,n,1)
                    i1 = mesh%edhx(e)%sh(1)
                    i2 = mesh%edhx(e)%sh(2)
                        
                    call monotonic_limiter(pol, i1, mesh, minval_i1, maxval_i1)
                    call monotonic_limiter(pol, i2, mesh, minval_i2, maxval_i2)
                        
                    minval_i = min(minval_i1, minval_i2)
                    maxval_i = max(maxval_i1, maxval_i2)
                        
                  !else
                  !  call monotonic_limiter(pol, i, mesh, minval_i, maxval_i)
                  !end if
                     
                  pol_filtered = max(pol,minval_i)
                  pol_filtered = min(pol_filtered,maxval_i)
                        
                  pol = pol_filtered
                end if

                ! Flux update
                node(i)%S(z)%flux = node(i)%S(z)%flux + pol*node(i)%G(s)%lwg(n)*dot

                ! Store flux at edge
                node(i)%edge_flux(n) = node(i)%edge_flux(n) + pol*node(i)%G(s)%lwg(n)*dot

             else
   
                if(controlvolume=="D")then
                  w=node(i)%G(s)%upwind_donald(n)
                else
                  w=node(i)%upwind_voronoi(n)
                endif

                !Determinando o valor do fluxo - UPWIND    
                xx=mesh%v(w)%p(1)
                yy=mesh%v(w)%p(2)
                zz=mesh%v(w)%p(3)
                call constr(xx,yy,zz,cx,sx,cy,sy)
                call aplyr(p3(1),p3(2),p3(3),cx,sx,cy,sy,xp,yp,zp) 

                if(order==2)then
                   ! Polynomial evaluated at quadrature point
                   pol = node(w)%coef(1) + node(w)%coef(2)*xp + & 
                         node(w)%coef(3)*yp

                else if(order==3)then
                   ! Polynomial evaluated at quadrature point
                   pol = node(w)%coef(1) + node(w)%coef(2)*xp + node(w)%coef(3)*yp + &
                        node(w)%coef(4)*xp*xp + node(w)%coef(5)*xp*yp + node(w)%coef(6)*yp*yp

                else if(order==4)then
                   ! Polynomial evaluated at quadrature point
                   pol = node(w)%coef(1) + node(w)%coef(2)*xp + node(w)%coef(3)*yp + &
                        node(w)%coef(4)*xp*xp + node(w)%coef(5)*xp*yp + node(w)%coef(6)*yp*yp + &                     
                        node(w)%coef(7)*xp*xp*xp + node(w)%coef(8)*xp*xp*yp + node(w)%coef(9)*xp*yp*yp + &
                        node(w)%coef(10)*yp*yp*yp
                end if

                ! Correct the value if monotonic limiter is applied
                if (monotonicfilter .and. time_integrator=='rk4')then
                  !if(moistswm)then ! To ensure mass conservation
                    e = edges_indexes(i,n,1)
                    i1 = mesh%edhx(e)%sh(1)
                    i2 = mesh%edhx(e)%sh(2)
                        
                    call monotonic_limiter(pol, i1, mesh, minval_i1, maxval_i1)
                    call monotonic_limiter(pol, i2, mesh, minval_i2, maxval_i2)
                        
                    minval_i = min(minval_i1, minval_i2)
                    maxval_i = max(maxval_i1, maxval_i2)
                        
                  !else
                  !    call monotonic_limiter(pol, i, mesh, minval_i, maxval_i)
                  !end if
                     
                  pol_filtered = max(pol,minval_i)
                  pol_filtered = min(pol_filtered,maxval_i)
 
                  pol = pol_filtered
                end if

                ! Flux update
                node(i)%S(z)%flux = node(i)%S(z)%flux  + pol*node(i)%G(s)%lwg(n)*dot 

                ! Store flux at edge
                node(i)%edge_flux(n) = node(i)%edge_flux(n) + pol*node(i)%G(s)%lwg(n)*dot


             endif
          enddo
       enddo
    enddo
    deallocate(p1,p2,p3)
    return  
  end subroutine flux_olg



    subroutine flux_gas(nodes,mesh,z,time, uedges)  
    !----------------------------------------------------------------------------------------------
    !    Calculando o fluxo nas arestas
    !----------------------------------------------------------------------------------------------

    implicit none 
    integer,intent(in)    :: nodes
    type(grid_structure),intent(inout)      :: mesh
    integer,intent(in),optional   :: z
    real(r8),intent(in),optional  :: time
    type(scalar_field), intent(inout), optional :: uedges ! Only for moist shallow water model
    integer(i4)   :: i
    integer(i4)   :: j
    integer(i4)   :: k
    integer(i4)   :: l
    integer(i4)   :: m
    integer(i4)   :: n
    integer(i4)   :: w
    integer(i4)   :: s
    integer(i4)   :: kk
    integer(i4)   :: cc
    integer(i4)   :: diml
    integer(i4)   :: contador
    integer(i4)   :: jend
    real(r8)  :: dot
    real(r8)  :: cx
    real(r8)  :: cy
    real(r8)  :: sx
    real(r8)  :: sy
    real(r8)  :: xx
    real(r8)  :: yy
    real(r8)  :: zz
    real(r8)  :: xp
    real(r8)  :: yp
    real(r8)  :: zp
    real(r8)  :: xb
    real(r8)  :: yb
    real(r8)  :: zb
    real(r8)  :: xc
    real(r8)  :: yc
    real(r8)  :: zc
    real(r8)  :: xi
    real(r8)  :: yi
    real(r8)  :: zi
    real(r8)  :: xf
    real(r8)  :: yf
    real(r8)  :: zf, theta, alfa, div_est, sol_numerica, sol_exata, derivada_i, derivada_j, cosseno, seno, aux
    real(r8)  :: FEPS, temp, aaxx, erro_Linf, erro, integral,phi_i, phi_j,dist,sinal, aux1, aux2, aux3, lon, lat
    real(r8)  :: lat1, lat2
    real(r8),allocatable  :: p1(:)      
    real(r8),allocatable  :: p2(:)      
    real(r8),allocatable  :: p3(:)      
    real(r8),allocatable  :: vn(:)      
    real(r8) :: dir_i(1:3), dir_j(1:3), rot_dir_i(1:3), rot_dir_j(1:3)

    allocate (p1(3),p2(3),p3(3),vn(3))
    FEPS = 1.0D-8
    dot = 0.0D0
    temp = 0.0D0
    aaxx = 0.0D0
    erro_Linf = 0.0D0
    contador = 0  

    do i = 1,nodes
         node(i)%S(z)%flux = 0.0D0
         jend = node(i)%ngbr(1)%numberngbr
         node(i)%edge_flux(:) = 0._r8
         do n = 1,jend     
           contador = contador + 1
           w=node(i)%upwind_voronoi(n)
           p3 =(mesh%v(i)%p+mesh%v(w)%p)/2.0D0 
           p3 = p3/norm(p3)       

           if (.not. moistswm)then
              dot = dot_product (velocity(p3, time),node(i)%G(1)%lvn(n,1:3))
           else
              dot = dot_product (node(i)%G(1)%velocity_quadrature(n)%v, node(i)%G(1)%lvn(n,1:3))
           endif

           if(dot>=0.0D0)then
              sinal=+1.0D0
           else
              sinal=-1.0D0
           endif
 
           !Determinando os valores para o node i 
           xx = mesh%v(i)%p(1)
           yy = mesh%v(i)%p(2)
           zz = mesh%v(i)%p(3)
           call constr(xx,yy,zz,cx,sx,cy,sy)

           xx = mesh%v(w)%p(1)
           yy = mesh%v(w)%p(2)
           zz = mesh%v(w)%p(3)
           call aplyr(xx,yy,zz,cx,sx,cy,sy,xp,yp,zp) 

 
           phi_i = node(i)%coef(1)

           ! Compute the second derivative in normal direction
           dir_i = node(i)%G(1)%lvn(n,1:3)
           dir_i = proj_vec_sphere(dir_i, mesh%v(i)%p)
           call aplyr(dir_i(1),dir_i(2),dir_i(3),cx,sx,cy,sy,rot_dir_i(1),rot_dir_i(2),rot_dir_i(3))
           rot_dir_i = rot_dir_i/norm(rot_dir_i)
           derivada_i = 2.d0*node(i)%coef(4)*rot_dir_i(1)*rot_dir_i(1) + &
                        2.d0*node(i)%coef(5)*rot_dir_i(1)*rot_dir_i(2) + &
                        2.d0*node(i)%coef(6)*rot_dir_i(2)*rot_dir_i(2)

           !Determinando os valores dos vizinhos do node i 
           xx = mesh%v(w)%p(1)
           yy = mesh%v(w)%p(2)
           zz = mesh%v(w)%p(3)
           call constr(xx,yy,zz,cx,sx,cy,sy)

           xx = mesh%v(i)%p(1)
           yy = mesh%v(i)%p(2)
           zz = mesh%v(i)%p(3)
           call aplyr(xx,yy,zz,cx,sx,cy,sy,xp,yp,zp) 

           phi_j = node(w)%coef(1)

           ! Compute second derivate in the normal direction to the edge
           dir_j = -node(i)%G(1)%lvn(n,1:3)
           dir_j = proj_vec_sphere(dir_j, mesh%v(w)%p)
           call aplyr(dir_j(1),dir_j(2),dir_j(3),cx,sx,cy,sy,rot_dir_j(1),rot_dir_j(2),rot_dir_j(3))
           rot_dir_j = rot_dir_j/norm(rot_dir_j)
           derivada_j = 2.d0*node(w)%coef(4)*rot_dir_j(1)*rot_dir_j(1) + &
                        2.d0*node(w)%coef(5)*rot_dir_j(1)*rot_dir_j(2) + &
                        2.d0*node(w)%coef(6)*rot_dir_j(2)*rot_dir_j(2)

           !Distancia entre o node i e seu respectivo vizinho   
           dist=arclen(mesh%v(i)%p,mesh%v(w)%p)           

           aux1 = (1.0D0/2.0D0)*(phi_i + phi_j)
           aux2 = (1.0D0/12.0D0)*((dist)**2)*(derivada_j + derivada_i)
           aux3 = sinal*(1.0D0/12.0D0)*((dist)**2)*(derivada_j - derivada_i)
           !aux3 = 0.d0
           !print*,aux2
           
           node(i)%S(z)%flux = node(i)%S(z)%flux  + (aux1-aux2+aux3)*node(i)%G(1)%lwg(n)*dot

           ! Store flux at edge
           node(i)%edge_flux(n) = node(i)%edge_flux(n) +  (aux1-aux2+aux3)*node(i)%G(1)%lwg(n)*dot

      end do
    end do
  deallocate(p1,p2,p3,vn)   
   return  
   end subroutine flux_gas   

   subroutine rungekutta3(nodes,mesh,k)  
      !----------------------------------------------------------------------------------------------
      !    AVançando no tempo
      !----------------------------------------------------------------------------------------------
      implicit none
      integer,intent(in):: nodes
      type(grid_structure),intent(in):: mesh
      integer,intent(in),optional:: k
      integer(i4):: i
      real(r8):: area
     
      select case (k)

      case  (1)
        do i=1,nodes
          if(controlvolume=='D') then
            area=node(i)%area
          else
            area=mesh%hx(i)%areag 
          endif
          node(i)%phi_new2 = node(i)%phi_old - 1.0D0/3.0D0 * (node(0)%dt/area)*node(i)%S(0)%flux 
        end do

      case  (2) 
        do i=1,nodes
          if(controlvolume=='D') then
            area=node(i)%area
          else
            area=mesh%hx(i)%areag
          endif
          node(i)%phi_new2 = node(i)%phi_old - 0.5D0 * (node(0)%dt/area)*node(i)%S(1)%flux
        end do

      case  (3) 
        do i=1,nodes
          if(controlvolume=='D') then
            area=node(i)%area
          else
            area=mesh%hx(i)%areag
          endif
          node(i)%phi_new2 = node(i)%phi_old  - (node(0)%dt/area) * node(i)%S(2)%flux     
        end do

      end select
      return
    end subroutine rungekutta3   

  subroutine ode_rk3_adv(nodes, mesh, time)
    integer(i4) :: k
    real(r8) :: time
    integer,intent(in):: nodes
    type(grid_structure),intent(inout):: mesh

    if(advmtd=='og2' .or. advmtd=='og3' .or. advmtd=='og4')then
      if(.not. monotonicfilter)then
        do k=1, 3
          call vector_olg2(nodes)
          call reconstruction_olg(nodes,mesh) 
          call flux_olg(nodes,mesh,k-1,time) 
          call rungekutta3(nodes,mesh,k) 
        enddo

      else
        do k=1, 2
          call vector_olg2(nodes)
          call reconstruction_olg(nodes,mesh) 
          call flux_olg(nodes,mesh,k-1,time) 
          call rungekutta3(nodes,mesh,k) 
        enddo

        ! Applies monotonic filter
        phi_RKS0%f(:) = node(1:nodes)%phi_old ! phi at time t
        phi_RKS2%f(:) = node(1:nodes)%phi_new2 ! phi at time t+dt/2 (second RK stage)
        call monotonicfilter_rk3(mesh, phi_RKS0, phi_RKS2, node(0)%dt, radius, time)
        node(1:nodes)%phi_new2 = phi_RKS2%f(1:nodes)

      end if

    else! Gassman scheme
      if(.not. monotonicfilter)then
        do k=1, 3
          call vector_gas(nodes,mesh)
          call reconstruction_gas(nodes,mesh)
          call flux_gas(nodes,mesh,k-1,time)
          call rungekutta3(nodes,mesh,k) 
        enddo

      else
        do k=1, 2
          call vector_gas(nodes,mesh)
          call reconstruction_gas(nodes,mesh)
          call flux_gas(nodes,mesh,k-1,time)
          call rungekutta3(nodes,mesh,k) 
        enddo

        ! Applies monotonic filter
        phi_RKS0%f(:) = node(1:nodes)%phi_old ! phi at time t
        phi_RKS2%f(:) = node(1:nodes)%phi_new2 ! phi at time t+dt/2 (second RK stage)
        call monotonicfilter_rk3(mesh, phi_RKS0, phi_RKS2, node(0)%dt, radius, time)
        node(1:nodes)%phi_new2 = phi_RKS2%f(1:nodes)
      end if

   end if

  end subroutine ode_rk3_adv


    subroutine ode_rk3_advection (mesh, phi_new, phi, time, dt, radius, u, u_new)
      !----------------------------------------------------------------------------------
      !! ode_rk3 takes one Runge-Kutta step for a vector ODE.
      !    t - time that will be calculated (t0+dt)
      !    h - scalar_field for thickness at current time
      !    u - scalar_field for velocities at current time
      !    dt - time step
      !    h_new and u_new - fields at t+dt
      !----------------------------------------------------------------------------------

      !Grid
      type(grid_structure), intent(inout) :: mesh

      !Scalar field (defined on voronoi centers)
      type(scalar_field), intent(inout):: phi  !General

      !Velocities (defined on edges - only normal component)
      type(scalar_field), optional, intent(inout):: u, u_new  !General

      !Time, time-step and sphere radius 
      real(r8):: time, dt, radius

      !Updated fields
      ! Scalar field (defined on voronoi centers)
      type(scalar_field), intent(inout):: phi_new  !General

      call zero_vector(phif0)
      call zero_vector(phif1)
      call zero_vector(phif2)
      call zero_vector(phif3)

      if(present(u))then
        call tendency_advection(phi, phif0, mesh, time, radius, u)
      else
        call tendency_advection(phi, phif0, mesh, time, radius)
      end if

      !First RK step
      phi_new%f(1:phi%n)   = phi%f(1:phi%n)    + dt * phif0(1:phi%n) / 3.0_r8

      if(present(u))then
        call tendency_advection(phi_new, phif1, mesh, time, radius, u)
      else
        call tendency_advection(phi_new, phif1, mesh, time, radius)
      end if

      !Second RK step
      phi_new%f(1:phi%n)   = phi%f(1:phi%n)    + dt * phif1(1:phi%n) / 2.0_r8

      if(present(u))then
        call tendency_advection(phi_new, phif2, mesh, time, radius, u)
      else
        call tendency_advection(phi_new, phif2, mesh, time, radius)
      end if

      ! Third  RK step
      ! Last RK step applies a different approach if the monotonic limiter is active 
      if(.not. monotonicfilter) then
        phi_new%f(1:phi%n)   = phi%f(1:phi%n)    + dt * phif2(1:phi%n)

      else if(monotonicfilter) then
        if(present(u))then
          call  monotonicfilter_rk3(mesh, phi, phi_new, dt, radius, time, u, u_new)
        else
          call  monotonicfilter_rk3(mesh, phi, phi_new, dt, radius, time)
        end if
      end if

      return
    end subroutine ode_rk3_advection

  subroutine tendency_advection(phi, phif, mesh, time, radius, u)
    !--------------------------------------
    !Calculates the Right Hand Side (spatial discret./tendency)
    !   of advection equation
    !-------------------------------------------------------------

    !Grid
    type(grid_structure), intent(inout) :: mesh

    !Tracer (defined on voronoi centers)
    type(scalar_field), intent(in):: phi  !General

    !Velocities (defined on edges - only normal component)
    type(scalar_field), optional, intent(inout):: u  !General

    !Time
    real(r8), intent(in):: time, radius

    !Right hand side of advction equation
    real(r8), intent(inout)::phif(:)

    !Initialize RHS
    call zero_vector(phif)

    !===============================================================
    !Calculate tracer phi tendency
    !===============================================================
    !Calculate divergence / tracer eq RHS
    if(present(u)) then
      call divhx(phi, div_uphi, mesh, radius, time, u)
    else
      call divhx(phi, div_uphi, mesh, radius, time)
    end if
    phif = -div_uphi%f
  return
  end subroutine tendency_advection

  subroutine rungekutta4(nodes,mesh,k)  
    !----------------------------------------------------------------------------------------------
    !    AVançando no tempo
    !----------------------------------------------------------------------------------------------
    implicit none
    integer,intent(in):: nodes
    type(grid_structure),intent(in):: mesh
    integer,intent(in),optional:: k
    integer(i4):: i
    real(r8):: area
   
    select case (k)
    case  (1)
    do i=1,nodes
      if(controlvolume=='D') then
        area=node(i)%area
      else
        area=mesh%hx(i)%areag 
      endif
        node(i)%phi_new2 = node(i)%phi_old - 0.5D0 * (node(0)%dt/area)*node(i)%S(0)%flux 
    end do

    case  (2) 
    do i=1,nodes
      if(controlvolume=='D') then
        area=node(i)%area
      else
        area=mesh%hx(i)%areag
      endif
        node(i)%phi_new2 = node(i)%phi_old - 0.5D0 * (node(0)%dt/area)*node(i)%S(1)%flux
    end do

    case  (3) 
    do i=1,nodes
      if(controlvolume=='D') then
        area=node(i)%area
      else
        area=mesh%hx(i)%areag
      endif
        node(i)%phi_new2 = node(i)%phi_old  - (node(0)%dt/area) * node(i)%S(2)%flux     
    end do

    case  (4) 
    do i=1,nodes
      if(controlvolume=='D') then
        area=node(i)%area
      else
        area=mesh%hx(i)%areag
      endif
        node(i)%phi_new2 = node(i)%phi_old  - (1.0D0/6.0D0)* (node(0)%dt/area) * (node(i)%S(0)%flux + &
                             2.0D0 * node(i)%S(1)%flux + 2.0D0 * node(i)%S(2)%flux + node(i)%S(3)%flux)  
    end do

    end select
    return
  end subroutine rungekutta4

  subroutine erro(nodes,mesh,time)  
    !----------------------------------------------------------------------------------------------
    !    AVançando no tempo
    !----------------------------------------------------------------------------------------------
    implicit none
    integer(i4),intent(in):: nodes
    type(grid_structure),intent(in):: mesh
    real(r8),intent(in),optional  :: time
    integer(i4):: i
    real(r8):: erro_Linf 
    real(r8):: erro_L1 
    real(r8):: erro_L2
    real(r8):: razao_L2
    real(r8):: razao_Linf
    real(r8):: area

    erro_Linf = 0.0D0
    erro_L1 = 0.0D0
    erro_L2 = 0.0D0
    razao_L2 = 0.0D0
    razao_Linf = 0.0D0

    do i=1,nodes
     if(method=='D') then
     area=node(i)%area
     else
     area=mesh%hx(i)%areag
     endif 
       !node(i)%erro=node(i)%phi_new2
       node(i)%erro=node(i)%phi_exa-node(i)%phi_new2
       phi%f(i) = node(i)%erro    
       if(erro_Linf<dabs(node(i)%erro)) erro_Linf= dabs(node(i)%erro) 
       if (razao_Linf < dabs (node(i)%phi_exa)) razao_Linf =  dabs (node(i)%phi_exa) 
       erro_L1=erro_L1 + dabs (node(i)%erro) * area
       erro_L2 = erro_L2 + (node(i)%erro) * (node(i)%erro) * area
       razao_L2 = razao_L2 + area*(node(i)%phi_exa)*(node(i)%phi_exa)
       write(atime,'(i8)') nint(time)
       phi%name=trim(adjustl(trim(transpname))//"_phi_t_"//trim(adjustl(trim(atime))))
    end do
    erro_L1 = erro_L1/razao_L2
    erro_L2 = dsqrt (erro_L2)
    razao_L2 = dsqrt (razao_L2)
    erro_L2 =  erro_L2/razao_L2
    erro_Linf = erro_Linf/razao_Linf
    call plot_scalarfield(phi,mesh)
    print*,'Erro_Linf',erro_Linf,'Erro L2',erro_L2
    return
  end subroutine erro

    !============================================================================================================
    ! TESTE GASSMANN NO PLANO 
    subroutine matrix_g(nodes,mesh)  
    !----------------------------------------------------------------------------------------------
    !    Construcao do sistema A*X = B e determinando a Pseudoinversa
    !----------------------------------------------------------------------------------------------
    implicit none
    integer(i4),intent(in)    :: nodes
    type(grid_structure),intent(in)      :: mesh
    integer(i4):: i
    integer(i4):: j
    integer(i4):: k
    integer(i4):: l
    integer(i4):: ii
    integer(i4):: ll
    integer(i4):: cc
    real(r8):: x
    real(r8):: y  
    real(r8):: xx
    real(r8):: yy
    real(r8):: zz
    real(r8):: cx
    real(r8):: sx
    real(r8):: cy
    real(r8):: sy
    real(r8):: xp
    real(r8):: yp
    real(r8):: zp
    real(r8),allocatable:: p(:)
    !real(r8),allocatable:: MRG(:,:)
    allocate(p(3))

    do i=1,nodes
       !Determinando os coeficientes da matriz de rotação para o node i 
       xx = mesh%v(i)%p(1)
       yy = mesh%v(i)%p(2)
       zz = mesh%v(i)%p(3)
       call constr(xx,yy,zz,cx,sx,cy,sy)
       !Percorrendo os primeiros vizinhos j-esimo do node i
       do j=1,node(i)%ngbr(1)%numberngbr+1
          k=node(i)%stencil(j-1)
          !Percorrendo os primeiros vizinhos l-esimo dos primeiros vizinhos j
          do l=1,node(k)%ngbr(1)%numberngbr
             p(1)=mesh%v(node(k)%stencil(l))%p(1)
             p(2)=mesh%v(node(k)%stencil(l))%p(2)
             p(3)=mesh%v(node(k)%stencil(l))%p(3)
             call aplyr(p(1),p(2),p(3),cx,sx,cy,sy,xp,yp,zp)          
             node(i)%edgeg(j)%MRGG(l,1)=xp 
             node(i)%edgeg(j)%MRGG(l,2)=yp
             node(i)%edgeg(j)%MRGG(l,3)=xp*xp
             node(i)%edgeg(j)%MRGG(l,4)=xp*yp
             node(i)%edgeg(j)%MRGG(l,5)=yp*yp
          enddo
          call pseudoinversa(node(node(i)%stencil(j-1))%ngbr(1)%numberngbr,5,node(i)%edgeg(j)%MRGG,node(i)%edgeg(j)%MPGG)
       end do
    end do
    return
    deallocate(p)
  end subroutine matrix_g



  subroutine vector_g(nodes,mesh)  
    !----------------------------------------------------------------------------------------------
    !    Construcao do sistema A*X = B e determinando a Pseudoinversa
    !----------------------------------------------------------------------------------------------
    implicit none
    integer(i4),intent(in):: nodes
    type(grid_structure),intent(in):: mesh
    integer(i4):: i
    integer(i4):: j
    integer(i4):: k
    integer(i4):: l

    do i=1,nodes
       do j=1,node(i)%ngbr(1)%numberngbr+1
          k=node(i)%stencil(j-1)
          do l=1,node(k)%ngbr(1)%numberngbr
             node(i)%edgeg(j)%VBGG(l) = f(mesh%v(node(k)%stencil(l))%lon,mesh%v(node(k)%stencil(l))%lat) &
                  - f(mesh%v(i)%lon,mesh%v(i)%lat)
          enddo
       end do
    end do
  end subroutine vector_g



  subroutine reconstruction_g(nodes,mesh)  
    !----------------------------------------------------------------------------------------------
    !    Reconstrucao da solucao
    !----------------------------------------------------------------------------------------------
    implicit none
    integer(i4),intent(in):: nodes
    type(grid_structure),intent(in):: mesh
    integer(i4):: i
    integer(i4):: j
    integer(i4):: m
    integer(i4):: n
    integer(i4):: l,ll
    integer(i4):: c
    integer(i4):: k
    integer(i4):: o
    real(r8):: aux

    do i=1,nodes
       do j=1,node(i)%ngbr(1)%numberngbr+1
          k=node(i)%stencil(j-1)
          node(i)%edgeg(j)%coefg(1)=node(k)%phi_new
          do ll=1,node(k)%ngbr(1)%numberngbr
             l=ubound(node(i)%edgeg(j)%MRGG,1)
             c=ubound(node(i)%edgeg(j)%MRGG,2)
             do m=1,c
                aux=0.0D0
                do n=1,l
                   aux=aux+node(i)%edgeg(j)%MPGG(m,n)*node(i)%edgeg(j)%VBGG(n)
                end do
                node(i)%edgeg(j)%coefg(m+1)=aux
             enddo
          enddo
       end do
    end do
    return 
  end subroutine reconstruction_g


  subroutine interpolation_g(mesh)  
    !----------------------------------------------------------------------------------------------
    !    Test of interpolation
    !----------------------------------------------------------------------------------------------
    implicit none
    type(grid_structure),intent(in):: mesh
    integer(i4):: i
    integer(i4):: j
    integer(i4):: nlon
    integer(i4):: nlat
    integer(i4):: pol_ref
    real(r8)  :: cx
    real(r8)  :: cy
    real(r8)  :: sx
    real(r8)  :: sy
    real(r8)  :: x
    real(r8)  :: y
    real(r8)  :: z
    real(r8)  :: xp
    real(r8)  :: yp
    real(r8)  :: zp
    real(r8):: lon
    real(r8):: lat
    real(r8):: dlon
    real(r8):: dlat
    real(r8):: tlon
    real(r8):: tlat
    real(r8):: erro
    real(r8):: erro_Linf
    real(r8):: erro_L2
    real(r8)::p(1:3)
    real(r8)::pp(1:3)
    real(r8),allocatable:: fexact(:,:)
    real(r8),allocatable:: fest(:,:)
    real(r8),allocatable:: ferror(:,:)

    !Lat-lon grid sizes
    nlat=720 !180 !720
    nlon=1440 !360 !1440
    dlat=180._r8/dble(nlat)!real(nlat, r8)
    dlon=360._r8/dble(nlon)!real(nlon, r8)

    allocate(fexact(1:nlon,1:nlat))
    allocate(fest(1:nlon,1:nlat))
    allocate(ferror(1:nlon,1:nlat))

    erro=0.0D0
    erro_Linf=0.0D0
    erro_L2=0.0D0

    !Pixel registration mode (for GMT graph)
    tlat=-90._r8+dlat/2._r8
    do j=1,nlat
       tlon=-180._r8+dlon/2._r8
       do i=1,nlon
          call sph2cart(tlon*deg2rad,tlat*deg2rad, pp(1), pp(2), pp(3))
          call cart2sph(pp(1), pp(2), pp(3),lon, lat)
          p = (/pp(1), pp(2), pp(3)/)

          ! Determinando o polinomio interpolador mais proximo do ponto p 
          pol_ref = getnearnode(p, mesh)

          ! Determinando o valor da funcao f e o valor do polinomio interpolador    
          x = mesh%v(pol_ref)%p(1)
          y = mesh%v(pol_ref)%p(2)
          z = mesh%v(pol_ref)%p(3)
          call constr(x,y,z,cx,sx,cy,sy)
          call aplyr(pp(1), pp(2), pp(3),cx,sx,cy,sy,xp,yp,zp) 

          !          fest(i,j) = node(pol_ref)%coef(1) + node(pol_ref)%coef(2)*xp + node(pol_ref)%coef(3)*yp + &
          !          node(pol_ref)%coef(4)*xp*xp + node(pol_ref)%coef(5)*xp*yp + node(pol_ref)%coef(6)*yp*yp

          fest(i,j) = node(pol_ref)%edgeg(1)%coefg(1) +  node(pol_ref)%edgeg(1)%coefg(2)*xp + &
               node(pol_ref)%edgeg(1)%coefg(3)*yp +  node(pol_ref)%edgeg(1)%coefg(4)*xp*xp + &
               node(pol_ref)%edgeg(1)%coefg(5)*xp*yp + node(pol_ref)%edgeg(1)%coefg(6)*yp*yp


          fexact(i,j)=f(lon,lat) 
          ferror(i,j)=fest(i,j) - fexact(i,j)

          erro = fest(i,j) - fexact(i,j)
          ! print*, fest(i,j)
          ! print*, fexact(i,j)
          ! print*, erro
          ! read(*,*)

          if (erro_Linf < abs (erro))  erro_Linf =  abs ( erro )  

          tlon=tlon+dlon
       end do
       tlat=tlat+dlat
    end do
    print*,erro_Linf,'erro_Linf'
    return
  end subroutine interpolation_g



  subroutine flux_edges_g(nodes,mesh,z,time)  
    !----------------------------------------------------------------------------------------------
    !    Calculando o fluxo nas arestas
    !----------------------------------------------------------------------------------------------
    implicit none 
    integer,intent(in)    :: nodes
    type(grid_structure),intent(in)      :: mesh
    integer,intent(in),optional   :: z
    real(r8),intent(in),optional  :: time
    integer(i4)   :: i
    integer(i4)   :: j
    integer(i4)   :: k
    integer(i4)   :: l
    integer(i4)   :: m
    integer(i4)   :: n
    integer(i4)   :: w
    integer(i4)   :: s
    integer(i4)   :: kk
    integer(i4)   :: cc
    integer(i4)   :: diml
    integer(i4)   :: contador
    integer(i4)   :: jend
    real(r8)  :: dot
    real(r8)  :: cx
    real(r8)  :: cy
    real(r8)  :: sx
    real(r8)  :: sy
    real(r8)  :: xx
    real(r8)  :: yy
    real(r8)  :: zz
    real(r8)  :: xp
    real(r8)  :: yp
    real(r8)  :: zp
    real(r8)  :: xb
    real(r8)  :: yb
    real(r8)  :: zb
    real(r8)  :: xc
    real(r8)  :: yc
    real(r8)  :: zc
    real(r8)  :: xi
    real(r8)  :: yi
    real(r8)  :: zi
    real(r8)  :: xf
    real(r8)  :: yf
    real(r8)  :: zf, theta, alfa,  div_est, delphi_i, delphi_j
    real(r8)  :: FEPS, temp, aaxx, erro_Linf, erro, integral,phi_i, phi_j,dist,sinal, aux1, aux2, aux3, lon, lat
    real(r8),allocatable  :: p1(:)      
    real(r8),allocatable  :: p2(:)      
    real(r8),allocatable  :: p3(:)      
    real(r8),allocatable  :: vn(:)      

    allocate (p1(3),p2(3),p3(3),vn(3))
    FEPS = 1.0D-8
    dot = 0.0D0
    temp = 0.0D0
    aaxx = 0.0D0
    erro_Linf = 0.0D0
    contador = 0
    do i = 1,nodes
       !       node(i)%S(z)%flux = 0.0D0
       jend = node(i)%ngbr(1)%numberngbr
       diml = 1
       cc=1
       m=1
       do n = 1,cc*jend     
          node(i)%S(z)%flux = 0.0D0
          ! Percorrendo todas as faces do volume de controle i
          p1 = node(i)%edge(n)%xyz2(1,:) 
          p2 = node(i)%edge(n)%xyz2(2,:) 
          !Determinando as coordenadas dos pontos de gauss de cada face     
          p3 =(p1+p2)/2.0D0 
          p3 = p3/norm(p3)
          !Determinando do vetor unitario normal as faces 
          vn=-cross_product(p1,p2)/norm(cross_product(p1,p2))   

          do s=1,diml
             contador = contador + 1
             ! Calculando o produto interno entre o vetor velocidade e o vetor normal unitario a edge 
             dot = dot_product (velocity(p3, time),vn)
             if(dot>=0.0D0)then
                sinal=+1.0D0
             else
                sinal=-1.0D0
             endif
             ! Determinando os coeficientes da matriz de rotação para o node i 
             xx = mesh%v(i)%p(1)
             yy = mesh%v(i)%p(2)
             zz = mesh%v(i)%p(3)
             call constr(xx,yy,zz,cx,sx,cy,sy)
             call aplyr(xx, yy, zz,cx,sx,cy,sy,xp,yp,zp)  

             ! Calculando o fluxo na esfera (Utilizando o polinomio interpolador no plano para estimar o valor na superficie esferica)
             if (order == 3) then   
                phi_i = node(i)%edgeg(1)%coefg(1) +  node(i)%edgeg(1)%coefg(2)*xp + &
                     node(i)%edgeg(1)%coefg(3)*yp +  node(i)%edgeg(1)%coefg(4)*xp*xp + &
                     node(i)%edgeg(1)%coefg(5)*xp*yp + node(i)%edgeg(1)%coefg(6)*yp*yp
                delphi_i=2.0D0*node(i)%edgeg(1)%coefg(3)


                !Determinando os coeficientes da matriz de rotação para vizinho do node i  
                w=node(i)%upwind_voronoi(n)
                xx = mesh%v(w)%p(1)
                yy = mesh%v(w)%p(2)
                zz = mesh%v(w)%p(3)
                !call constr(xx,yy,zz,cx,sx,cy,sy)
                !call aplyr(pg(1), pg(2), pg(3),cx,sx,cy,sy,xp,yp,zp)
                m=m+1
                call aplyr(xx, yy, zz,cx,sx,cy,sy,xp,yp,zp)          
                phi_j=node(i)%edgeg(m)%coefg(1) +  node(i)%edgeg(m)%coefg(2)*xp + &
                     node(i)%edgeg(m)%coefg(3)*yp +  node(i)%edgeg(m)%coefg(4)*xp*xp + &
                     node(i)%edgeg(m)%coefg(5)*xp*yp + node(i)%edgeg(m)%coefg(6)*yp*yp
                delphi_j=2.0D0*node(i)%edgeg(m)%coefg(3) 

                !print*,node(i)%edgeg(m)%coefg(1)
                !print*,phi_j
                !read(*,*)

                dist=arclen(mesh%v(i)%p,mesh%v(w)%p)

                node(i)%S(z)%flux = ((1.0D0/2.0D0)*(phi_i+phi_j) + (1.0D0/12.0D0)*(dist**2)*(delphi_i+delphi_j) + &
                sinal*(dist**2)*(0.25D0/12.0D0)*(delphi_j-delphi_i))*dist

             end if
          end do
          call numerical_integration(mesh%v(i)%p,mesh%v(w)%p,vn,integral)
          erro = dabs(node(i)%S(z)%flux - integral)  
          aaxx = erro * erro
          if (erro_Linf < abs (erro))  erro_Linf =  abs ( erro )
          temp = temp +  aaxx
       end do
    end do
    print*, erro_Linf, 'erro_Linf'
    print*, dsqrt (temp/contador) , 'erro_L2'

    !    do i = 1,nodes
    !       div_est = node(i)%S(z)%flux/mesh%hx(i)%areag
    !       erro = dabs(div_est)
    !       aaxx = erro * erro
    !       if (erro_Linf < abs (erro))  erro_Linf =  abs ( erro )
    !       temp = temp +  aaxx
    !    end do
    !    print*, erro_Linf, 'norma 1'
    !    print*, dsqrt (temp/nodes) , 'norma 2'    

    deallocate(p1,p2,p3,vn)
    return  
  end subroutine flux_edges_g


  subroutine init_quadrature_edges(mesh)  
    !----------------------------------------------------------------------------------------------
    ! Initialize some index variables regarding quadrature edges
    !----------------------------------------------------------------------------------------------
    implicit none

    type(grid_structure),intent(inout):: mesh
    integer(i4):: i
    integer(i4):: j
    integer(i4):: j1, j2
    integer(i4):: e
    integer(i4):: k, q, nquad, maxnnb
    real(r8):: p(1:3), lon, lat

    ! Number of quadrature points
    nquad = nint((order)/2.0D0)

    if (method=='G') then
      nquad = 1
    end if

    ! Maximum of neighbors in a Voronoi cell
    maxnnb = maxval(mesh%v(:)%nnb)
    ! Allocation
    allocate(sh_edges_indexes(1:mesh%ne, 1:4))
    allocate(edges_indexes(1:mesh%nv, 1:maxnnb, 1:nquad))
    sh_edges_indexes = 0

    do i = 1, mesh%nv
      do j = 1, mesh%v(i)%nnb
        do q = 1, nquad
          ! Get the quadrature points
          p = node(i)%G(q)%lpg(j,1:3)

          ! Get the nearest edge index
          e =  getnearhxedge(p, mesh)

          ! Store edge index
          edges_indexes(i,j,q) = e

          ! Save the index j
          if(sh_edges_indexes(e,1) == 0)then
            sh_edges_indexes(e,1) = j
          else
            sh_edges_indexes(e,2) = j
          end if
        end do
      end do
    end do
  end subroutine init_quadrature_edges

  subroutine reconstruct_velocity_quadrature(mesh, uedges)  
    !----------------------------------------------------------------------------------------------
    ! Reconstructs the normal velocity at quadrature points
    !----------------------------------------------------------------------------------------------
    implicit none

    type(grid_structure),intent(inout):: mesh
    type(scalar_field), intent(inout) :: uedges
    real(r8):: urecon(1:3)
    real(r8):: p(1:3), p1(1:3), p2(1:3)
    real(r8):: u0, utmp, vtmp, lat, lon, aux, aux1, aux2,error, uexact(1:3)
    integer(i4):: i, j, k, e, q
    integer(i4):: nquad
    integer(i4):: i1, i2, j1, j2, q1, q2
    character (len=8) :: reconadvmtd = "lsqhxe"

    ! Number of quadrature points
    nquad = nint((order)/2.0D0)

    !Alocate space if necessary
    if(.not.allocated(uedges%pol))then
      allocate(uedges%pol(1:mesh%nv))
    end if

    ! Olliver-Gooch method
    if (method=='O')then
      !$omp parallel do &
      !$omp default(none) &
      !$omp shared(mesh, uedges, node, sh_edges_indexes, nquad) &
      !$omp private(i1, i2, j1, j2, p1, p2, q1, q2, q,lon,lat,utmp,vtmp,u0,uexact) &
      !$omp private(urecon) &
      !$omp schedule(static)
      do e = 1, mesh%ne ! Edges loop
        do q = 1, nquad ! Quadrature points loop
          q1 = q
          if (q==1) then
            q2 = nquad
          else
            q2 = 1
          end if

          !Get neighbor Voronoi cells
          i1 = mesh%edhx(e)%sh(1)
          i2 = mesh%edhx(e)%sh(2)
          
          ! Get edge e indexes
          j1 = sh_edges_indexes(e,1)
          j2 = sh_edges_indexes(e,2)
          
          ! Quadrature points
          p1 = node(i1)%G(q1)%lpg(j1,1:3)
          p2 = node(i2)%G(q2)%lpg(j2,1:3)
         
          !print*,p1
          ! Reconstruct the velocity field at quadrature points
          urecon = vecrecon_lsq_ed(p1, uedges, mesh, e)
          !call cart2sph(p1(1), p1(2), p1(3), lon, lat)
          !u0 = 2._r8*pi*erad/(12._r8*day2sec)
          !utmp = u0*cos(lat)
          !vtmp = 0._r8
          !call convert_vec_sph2cart(utmp, vtmp, p1, uexact)
          !urecon=uexact 

          ! Store the velocity
          node(i1)%G(q1)%velocity_quadrature(j1)%v = urecon
          node(i2)%G(q2)%velocity_quadrature(j2)%v = urecon

          !p = p1
          !call cart2sph(p1(1), p1(2), p1(3), lon, lat)
          !u0 = 2._r8*pi*erad/(12._r8*day2sec)
          !utmp = u0*cos(lat)
          !vtmp = 0._r8
          !call convert_vec_sph2cart(utmp, vtmp, p1, uexact)
          !error = max(error, norm2(uexact-urecon)/norm2(uexact))
      end do    
    end do
    !$omp end parallel do

    else
      ! Gassman method
      !$omp parallel do &
      !$omp default(none) &
      !$omp shared(mesh, uedges, node, sh_edges_indexes, nquad, reconadvmtd) &
      !$omp private(j, k, p) &
      !$omp private(urecon) &
      !$omp schedule(static)
      do i = 1,mesh%nv
        do j = 1,node(i)%ngbr(1)%numberngbr
          k=node(i)%upwind_voronoi(j)

          !Reconstruction point
          p =(mesh%v(i)%p+mesh%v(k)%p)/2.0D0 
          p = p/norm(p)

          ! Reconstruct the vector field at given point
          urecon = vector_reconstruct (p, uedges, mesh, reconadvmtd)

          ! Store the velocity field
          node(i)%G(1)%velocity_quadrature(j)%v = urecon 
        end do
      end do
      !$omp end parallel do
    endif

  end subroutine reconstruct_velocity_quadrature

  subroutine monotonic_limiter(Pol,no,mesh,min_output, max_output)
      implicit none
      
      type(grid_structure),intent(in)      :: mesh
      integer (i4), intent(in) :: no
      real (r8), intent (in)   :: Pol
      real (r8), intent (out)  :: min_output, max_output
      real (r8), allocatable   :: limitador (:)
      integer   :: i,l,j,jend
      real(r8)  :: maximo,minimo
     
      i = no
       
      if (order == 2) then
        jend = node(i)%ngbr(1)%numberngbr
        allocate (limitador (1:jend+1))
      else 
        jend = node(i)%ngbr(1)%numberngbr+node(i)%ngbr(2)%numberngbr
        allocate (limitador (1:jend+1))
      end if       
  
      limitador(1) = node(i)%coef(1)    
      do l=1,jend 
        j = node(i)%stencil(l)
        limitador(l+1) = node(j)%coef(1)
      end do 
      
      max_output = maxval(limitador)
      min_output = minval(limitador)

      deallocate (limitador)
      return
  end subroutine monotonic_limiter
    
  !-----------------------------------------------------------------------------------
  ! Implementation of the flux correction method from the paper Wang et al 2009 -
  ! "Evaluation of Scalar Advection Schemes in the Advanced Research WRF
  ! Model Using Large-Eddy Simulations of Aerosol–Cloud Interactions"
  ! This routine is applied in the last stage of the RK3 scheme
  !-----------------------------------------------------------------------------------
  subroutine monotonicfilter_rk3(mesh, phi_step0, phi_step2,  dt, radius, time, u_step0, u_step2, hSphi)
      type(grid_structure), intent(inout) :: mesh
      type(scalar_field), intent(inout):: phi_step0 ! scalar field at time t
      type(scalar_field), intent(inout):: phi_step2 ! Scalar field from second step in RK3 (time t+dt/2)
      type(scalar_field), optional, intent(inout):: u_step0 ! velocity at time t
      type(scalar_field), optional, intent(inout):: u_step2 ! velocity at time t+dt/
      type(scalar_field), optional, intent(inout):: hSphi ! Source - optional
      real(r8), intent(in) :: dt, radius, time ! time-step and sphere radius 
      integer(i4) :: i, j, k, jj

      !-----------------------------------------------------------------------------------
      ! Flux for phi_star using 1st order upwind scheme at time t
      if (present(hSphi))then ! check if source was given
        !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
        !$OMP SHARED(phi_star, phi_step0, dt, hSphi)
        phi_star%f = phi_step0%f + dt*hSphi%f ! equation 3b from Wang et al 2009
        !$OMP END PARALLEL WORKSHARE
      else
        !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
        !$OMP SHARED(phi_star, phi_step0)
        phi_star%f = phi_step0%f ! equation 3b from Wang et al 2009
        !$OMP END PARALLEL WORKSHARE
      end if

      ! Compute upwind flux
      if(present(u_step0))then
        call flux_hx_upw1(phi_star, F_star, mesh, time, u_step0)
      else
        call flux_hx_upw1(phi_star, F_star, mesh, time)
      end if

      if(controlvolume=='V')then
        !-----------------------------------------------------------------------------------
        ! Flux for phi from the previous RK step using highorder scheme at time t+dt/2
        if (advmtd=='sg3' .or. advmtd=='og2' .or. advmtd=='og3' .or. advmtd=='og4')then
          !Calculate divergence and edges flux
          if (present(u_step2))then
            call divhx(phi_step2, div_uphi, mesh, radius, time, u_step2)
          else
            call divhx(phi_step2, div_uphi, mesh, radius, time)
          end if

          ! Store the flux
          !$omp parallel do &
          !$omp default(none) &
          !$omp shared(mesh, node, F_step2) &
          !$omp schedule(static)
          do i = 1, mesh%nv
            F_step2(i,1) = node(i)%edge_flux(mesh%v(i)%nnb)
            F_step2(i,2:mesh%v(i)%nnb) = node(i)%edge_flux(1:mesh%v(i)%nnb-1)
          end do
          !$omp end parallel do

        else if (advmtd=='upw1')then
          if (present(u_step2))then
            call flux_hx_upw1(phi_step2,F_step2, mesh, time, u_step2)
          else
            print*,'monotonicfilter_rk3 error in 1st order upwind: need u_step2'
            stop
          end if

        else if (advmtd=='trsk')then
          if (present(u_step2))then
            !Interpolate vapour to edges and calculate flux at edges
            call scalar_hx2ed(phi_step2, phi_ed, mesh)      !hQv: cell->edge
            call scalar_elem_product(u_step2, phi_ed, uphi) !Flux uhQv at edges
            call flux_hx(uphi, F_step2, mesh)
          else
            print*,'monotonicfilter_rk3 error in trsk: u_step2 missing'
            stop
          end if

        else
          print*, 'ERROR in monotonicfilter_rk3: invalid advmth:  ', advmtd
          stop
        end if
        !-----------------------------------------------------------------------------------


        !-----------------------------------------------------------------------------------
        ! Flux for phi corrected (equation 4 from wang et al 2009)
        !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
        !$OMP SHARED(F_cor, F_step2, F_star)      
        F_cor(:,:) = F_step2(:,:) - F_star(:,:) 
        !$OMP END PARALLEL WORKSHARE
        !-----------------------------------------------------------------------------------

        ! Equation 5 from Wang et al 2009 - 1st order upwind solution
        !$omp parallel do &
        !$omp default(none) &
        !$omp shared(mesh, phi_tilda, phi_star, F_star, dt, radius) &
        !$omp schedule(static)
        do i = 1, mesh%nv
          phi_tilda%f(i) = phi_star%f(i) - dt*sum(F_star(i,:))/mesh%hx(i)%areag/radius
        end do
        !$omp end parallel do

        ! Compute negative and positive flux updates
        !$omp parallel do &
        !$omp default(none) &
        !$omp shared(mesh, phi_tilda, phi_tilda_min, phi_tilda_max, F_cor, dt, radius) &
        !$omp schedule(static)
        do i = 1, mesh%nv
          phi_tilda_min%f(i) = phi_tilda%f(i)
          phi_tilda_max%f(i) = phi_tilda%f(i)
          do j = 1, mesh%v(i)%nnb
            if(F_cor(i,j)>0._r8)then
              ! Equation 6a from Wang et al 2009
              phi_tilda_min%f(i) = phi_tilda_min%f(i) - dt*F_cor(i,j)/mesh%hx(i)%areag/radius

            else
              ! Equation 6b from Wang et al 2009
              phi_tilda_max%f(i) = phi_tilda_max%f(i) - dt*F_cor(i,j)/mesh%hx(i)%areag/radius
            end if
          end do
        end do
        !$omp end parallel do

        ! Min/max in neighborhood at time t
        !$omp parallel do &
        !$omp default(none) &
        !$omp shared(mesh, phi_step0, phi_min, phi_max) &
        !$omp private(j, k) &
        !$omp schedule(static)
        do i = 1, mesh%nv
          phi_min%f(i) = phi_step0%f(i)
          phi_max%f(i) = phi_step0%f(i)
          do j = 1, mesh%v(i)%nnb
            k = mesh%v(i)%nb(j)
            if(phi_step0%f(k) < phi_min%f(i))then
              phi_min%f(i) = phi_step0%f(k)
            else if(phi_step0%f(k) > phi_max%f(i))then
              phi_max%f(i) = phi_step0%f(k)
            end if
         end do
        end do
        !$omp end parallel do

        ! Flux correction
        !$omp parallel do &
        !$omp default(none) &
        !$omp shared(mesh, F_cor, phi_tilda, phi_tilda_min, phi_tilda_max, phi_min, phi_max) &
        !$omp private(j, k) &
        !$omp schedule(static)
        do i = 1, mesh%nv
          do j = 1, mesh%v(i)%nnb
            k = mesh%v(i)%nb(j)
            if(F_cor(i,j)>0._r8) then
              if(abs(phi_tilda%f(i)-phi_tilda_min%f(i))>eps2 .and.  abs(phi_tilda%f(k)-phi_tilda_max%f(k))>eps2) then
              F_cor(i,j) = min(1._r8,  (phi_tilda%f(i)-phi_min%f(i))/(phi_tilda%f(i)-phi_tilda_min%f(i)),&
                                       (phi_tilda%f(k)-phi_max%f(k))/(phi_tilda%f(k)-phi_tilda_max%f(k)))*F_cor(i,j)
              end if
            else
              if(abs(phi_tilda%f(i)-phi_tilda_max%f(i))>eps2 .and.  abs(phi_tilda%f(k)-phi_tilda_min%f(k))>eps2) then
                F_cor(i,j) =  min(1._r8,  (phi_tilda%f(i)-phi_max%f(i))/(phi_tilda%f(i)-phi_tilda_max%f(i)),&
                                          (phi_tilda%f(k)-phi_min%f(k))/(phi_tilda%f(k)-phi_tilda_min%f(k)))*F_cor(i,j)
              end if
            end if
          end do
        end do
        !$omp end parallel do

        ! Final solution
        !$omp parallel do &
        !$omp default(none) &
        !$omp shared(mesh, phi_step2, phi_tilda, F_cor, radius, dt) &
        !$omp schedule(static)
        do i = 1, mesh%nv
          phi_step2%f(i) = phi_tilda%f(i) - dt*(sum(F_cor(i,:)))/mesh%hx(i)%areag/radius
        end do
        !$omp end parallel do

      !======================================================================================
      else ! Donald Diagram
        if (advmtd=='sg3' .or. advmtd=='og2' .or. advmtd=='og3' .or. advmtd=='og4')then
          call divhx(phi_step2, div_uphi, mesh, radius, time)
        end if

        ! Store the flux
        !$omp parallel do &
        !$omp default(none) &
        !$omp shared(mesh, node, F_step2) &
        !$omp schedule(static)
        do i = 1, mesh%nv
          F_step2(i,1:2) = node(i)%edge_flux(2*mesh%v(i)%nnb-1:2*mesh%v(i)%nnb)
          F_step2(i,3:2*mesh%v(i)%nnb) = node(i)%edge_flux(1:2*mesh%v(i)%nnb-2)
        end do
        !$omp end parallel do

        !-----------------------------------------------------------------------------------
        ! Flux for phi corrected (equation 4 from wang et al 2009)
        !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
        !$OMP SHARED(F_cor, F_step2, F_star)      
        F_cor(:,:) = F_step2(:,:) - F_star(:,:) 
        !$OMP END PARALLEL WORKSHARE
        !-----------------------------------------------------------------------------------

        ! Equation 5 from Wang et al 2009 - 1st order upwind solution
        !$omp parallel do &
        !$omp default(none) &
        !$omp shared(mesh, node, phi_tilda, phi_star, F_star, dt, radius) &
        !$omp schedule(static)
        do i = 1, mesh%nv
          phi_tilda%f(i) = phi_star%f(i) - dt*sum(F_star(i,:))/node(i)%area/radius
        end do
        !$omp end parallel do

        ! Compute negative and positive flux updates
        !$omp parallel do &
        !$omp default(none) &
        !$omp shared(mesh, node, phi_tilda, phi_tilda_min, phi_tilda_max, F_cor, dt, radius) &
        !$omp schedule(static)
        do i = 1, mesh%nv
          phi_tilda_min%f(i) = phi_tilda%f(i)
          phi_tilda_max%f(i) = phi_tilda%f(i)
          do j = 1, 2*mesh%v(i)%nnb
            if(F_cor(i,j)>0._r8)then
              ! Equation 6a from Wang et al 2009
              phi_tilda_min%f(i) = phi_tilda_min%f(i) - dt*F_cor(i,j)/node(i)%area/radius

            else
              ! Equation 6b from Wang et al 2009
              phi_tilda_max%f(i) = phi_tilda_max%f(i) - dt*F_cor(i,j)/node(i)%area/radius
            end if
          end do
        end do
        !$omp end parallel do

        ! Min/max in neighborhood at time t
        !$omp parallel do &
        !$omp default(none) &
        !$omp shared(mesh, phi_step0, phi_min, phi_max) &
        !$omp private(j, k) &
        !$omp schedule(static)
        do i = 1, mesh%nv
          phi_min%f(i) = phi_step0%f(i)
          phi_max%f(i) = phi_step0%f(i)
          do j = 1, mesh%v(i)%nnb
            k = mesh%v(i)%nb(j)
            if(phi_step0%f(k) < phi_min%f(i))then
              phi_min%f(i) = phi_step0%f(k)
            else if(phi_step0%f(k) > phi_max%f(i))then
              phi_max%f(i) = phi_step0%f(k)
            end if
         end do
        end do
        !$omp end parallel do

        ! Flux correction
        !$omp parallel do &
        !$omp default(none) &
        !$omp shared(mesh, F_cor, phi_tilda, phi_tilda_min, phi_tilda_max, phi_min, phi_max) &
        !$omp private(j, jj, k) &
        !$omp schedule(static)
        do i = 1, mesh%nv
          jj = 1
          do j = 1, mesh%v(i)%nnb
            k = mesh%v(i)%nb(j)
            if(F_cor(i,j)>0._r8) then
              if(abs(phi_tilda%f(i)-phi_tilda_min%f(i))>eps2 .and.  abs(phi_tilda%f(k)-phi_tilda_max%f(k))>eps2) then
                F_cor(i,jj) = min(1._r8,  (phi_tilda%f(i)-phi_min%f(i))/(phi_tilda%f(i)-phi_tilda_min%f(i)),&
                                       (phi_tilda%f(k)-phi_max%f(k))/(phi_tilda%f(k)-phi_tilda_max%f(k)))*F_cor(i,jj)
                F_cor(i,jj+1) = min(1._r8,  (phi_tilda%f(i)-phi_min%f(i))/(phi_tilda%f(i)-phi_tilda_min%f(i)),&
                                       (phi_tilda%f(k)-phi_max%f(k))/(phi_tilda%f(k)-phi_tilda_max%f(k)))*F_cor(i,jj+1)
 
              end if
            else
              if(abs(phi_tilda%f(i)-phi_tilda_max%f(i))>eps2 .and.  abs(phi_tilda%f(k)-phi_tilda_min%f(k))>eps2) then
                F_cor(i,jj) =  min(1._r8,  (phi_tilda%f(i)-phi_max%f(i))/(phi_tilda%f(i)-phi_tilda_max%f(i)),&
                                          (phi_tilda%f(k)-phi_min%f(k))/(phi_tilda%f(k)-phi_tilda_min%f(k)))*F_cor(i,jj)
                F_cor(i,jj+1) =  min(1._r8,  (phi_tilda%f(i)-phi_max%f(i))/(phi_tilda%f(i)-phi_tilda_max%f(i)),&
                                          (phi_tilda%f(k)-phi_min%f(k))/(phi_tilda%f(k)-phi_tilda_min%f(k)))*F_cor(i,jj+1)
              end if
            end if
            jj = jj+2
          end do
        end do
        !$omp end parallel do

        ! Final solution
        !$omp parallel do &
        !$omp default(none) &
        !$omp shared(mesh, node, phi_step2, phi_tilda, F_cor, F_star,F_step2, radius, dt) &
        !$omp schedule(static)
        do i = 1, mesh%nv
          phi_step2%f(i) = phi_tilda%f(i) - dt*(sum(F_cor(i,:)))/node(i)%area/radius
        end do
        !$omp end parallel do
      end if
  end subroutine monotonicfilter_rk3

  subroutine divhx(q, div, mesh, radius, time, u)
    !---------------------------------------------------------------
    !Calculate divergence at voronoi cells (hexagons)
    !   based on edge normal velocities
    !---------------------------------------------------------------
    type(grid_structure), intent(inout) :: mesh
    type(scalar_field), optional, intent(inout):: u ! velocity at cell edges
    type(scalar_field), intent(in):: q ! scalar at cell center
    type(scalar_field), intent(inout):: div !divergence - must be already allocated
    real(r8), intent(in) :: radius, time
    real(r8) :: area
    ! High order variables
    integer(i4) :: i
    
    if (advmtd=='sg3' .or. advmtd=='og2' .or. advmtd=='og3'.or. advmtd=='og4') then
      !$omp parallel do &
      !$omp default(none) &
      !$omp shared(mesh, q, node) &
      !$omp schedule(static)
      do i = 1, mesh%nv
        node(i)%phi_new2 = q%f(i)
        node(i)%phi_new  = q%f(i)
      end do
      !$omp end parallel do

      if (advmtd == 'og2' .or. advmtd=='og3'.or. advmtd=='og4')then ! Ollivier-Gooch method
        call vector_olg2(mesh%nv)
        call reconstruction_olg(mesh%nv, mesh) 
        if (present(u))then
          call flux_olg(mesh%nv, mesh, 0, time, u)
        else
          call flux_olg(mesh%nv, mesh, 0, time)
        end if

      else ! Gassman method
        call vector_gas(mesh%nv, mesh)
        call reconstruction_gas(mesh%nv, mesh) 
        if (present(u))then
          call flux_gas(mesh%nv, mesh, 0, time, u)
        else
          call flux_gas(mesh%nv, mesh, 0, time)
        end if
      end if

      !$omp parallel do &
      !$omp default(none) &
      !$omp shared(mesh, div, controlvolume, node, radius) &
      !$omp private(area) &
      !$omp schedule(static)
      do i = 1, mesh%nv
        if (controlvolume=='V')then
          area = mesh%hx(i)%areag 
        else
          area = node(i)%area
        end if

        div%f(i) = node(i)%S(0)%flux/area
        div%f(i) = div%f(i)/radius
      end do
      !$omp end parallel do

    else if(advmtd=='upw1') then
        call div_hx_upw1(div, q, u, mesh, radius)

    else if(advmtd=='trsk') then
        !Interpolate vapour to edges and calculate flux at edges
        call scalar_hx2ed(q, phi_ed, mesh)      !hQv: cell->edge
        call scalar_elem_product(u, phi_ed, uphi) !Flux uhQv at edges

        !Calculate divergence / vapour eq RHS
        call div_hx(uphi, div_uphi, mesh)
    end if

    return
  end subroutine divhx


  subroutine flux_hx_upw1(q, flux, mesh, time, u)
    !---------------------------------------------------------------
    !Calculate 1st order upwind flux  at voronoi edges of hexagons
    !---------------------------------------------------------------
    type(grid_structure), intent(in) :: mesh
    type(scalar_field), optional, intent(in):: u ! velocity at cell edges
    type(scalar_field), intent(in):: q ! scalar field at cell centers
    real(r8), intent(in) :: time
    real(r8), allocatable, intent(inout) :: flux(:,:)
    integer(i4):: i, j, k, l, ed, cc, jj
    real(r8):: signcor, u_edge

    if(present(u))then
      !$omp parallel do &
      !$omp default(none) &
      !$omp shared(mesh, u, q, flux) &
      !$omp private(j, l, k, signcor) &
      !$omp schedule(static)
      do i=1,mesh%nv
        !Divergence of uh on unit sphere using div_cell_Cgrig
        !For all edges forming the hexagon
        flux(i,:)=0._r8
        do j=1, mesh%v(i)%nnb
          !Get edge index
          l=mesh%v(i)%ed(j)

          !Get neighbor cell index
          k=mesh%v(i)%nb(j)

          !Get edge outer normal related to the hexagon
          signcor=real(mesh%hx(i)%ttgout(j), r8)

          if(signcor*u%f(l)>0._r8)then
            flux(i,j) = signcor*u%f(l)*mesh%edhx(l)%leng*q%f(i)
          else
            flux(i,j) = signcor*u%f(l)*mesh%edhx(l)%leng*q%f(k)
          end if
        end do
      end do

    else
      if(controlvolume=="V")then
        !$omp parallel do &
        !$omp default(none) &
        !$omp shared(mesh, q, time, flux) &
        !$omp private(j, l, k, u_edge, signcor) &
        !$omp schedule(static)
        do i=1,mesh%nv
          flux(i,:)=0._r8
          do j=1, mesh%v(i)%nnb
            !Get edge index
            l=mesh%v(i)%ed(j)

            !Get neighbor cell index
            k=mesh%v(i)%nb(j)

            !Get edge outer normal related to the hexagon
            signcor=real(mesh%hx(i)%ttgout(j), r8)

            u_edge = dot_product(velocity(mesh%ed(l)%c%p, time), mesh%ed(l)%tg)

            if(signcor*u_edge>0._r8)then
              flux(i,j) = signcor*u_edge*mesh%edhx(l)%leng*q%f(i)
            else
              flux(i,j) = signcor*u_edge*mesh%edhx(l)%leng*q%f(k)
            end if
          end do
        end do
        !$omp end parallel do

      else 
        !$omp parallel do &
        !$omp default(none) &
        !$omp shared(mesh, q, time, flux) &
        !$omp private(j, l, k, jj, u_edge, signcor) &
        !$omp schedule(static)
        do i=1,mesh%nv
          flux(i,:)=0._r8
          jj=1
          do j=1, mesh%v(i)%nnb
            !Get edge index
            l=mesh%v(i)%ed(j)

            !Get neighbor cell index
            k=mesh%v(i)%nb(j)

            !Get edge outer normal related to the hexagon
            signcor=real(mesh%hx(i)%ttgout(j), r8)
            u_edge = dot_product(velocity(mesh%ed(l)%c%p, time), mesh%ed(l)%tg)

            if(signcor*u_edge>0._r8)then
              flux(i,jj) = signcor*u_edge*mesh%edhx(l)%leng*q%f(i)*0.5_r8
            else
              flux(i,jj) = signcor*u_edge*mesh%edhx(l)%leng*q%f(k)*0.5_r8
            end if
            flux(i,jj+1) = flux(i,jj)
            jj = jj+2
          end do
        end do
        !$omp end parallel do
      end if
   end if


   return

  end subroutine flux_hx_upw1

  subroutine div_hx_upw1(div, q, u, mesh, radius)
    !---------------------------------------------------------------
    !Calculate divergence at voronoi cells (hexagons)
    !   based on edge normal velocities
    !---------------------------------------------------------------
    type(grid_structure), intent(in) :: mesh
    type(scalar_field), intent(in):: u ! velocity at cell edges
    type(scalar_field), intent(in):: q ! scalar field at cell centers
    type(scalar_field), intent(inout):: div !divergence - must be already allocated
    real(r8), intent(in) :: radius
    integer(i4):: i, j, k, l, ed
    real(r8):: signcor

    !$omp parallel do &
    !$omp default(none) &
    !$omp shared(mesh, q, u, div, radius) &
    !$omp private(j, l, k, signcor) &
    !$omp schedule(static)
    do i=1,mesh%nv
      !Divergence of uh on unit sphere using div_cell_Cgrig
      !For all edges forming the hexagon
      div%f(i)=0._r8
      do j=1, mesh%v(i)%nnb
        !Get edge index
        l=mesh%v(i)%ed(j)

        !Get neighbor cell index
        k=mesh%v(i)%nb(j)

        !Get edge outer normal related to the hexagon
        signcor=real(mesh%hx(i)%ttgout(j), r8)

        !Calculate numerical integration
        if(signcor*u%f(l)>0._r8)then
          div%f(i)=div%f(i)+signcor*u%f(l)*mesh%edhx(l)%leng*q%f(i)
        else
          div%f(i)=div%f(i)+signcor*u%f(l)*mesh%edhx(l)%leng*q%f(k)
        end if
      end do
      div%f(i)=div%f(i)/mesh%hx(i)%areag/radius
    end do
    !$omp end parallel do

    return

  end subroutine div_hx_upw1

  subroutine flux_hx(u, flux, mesh)
    !---------------------------------------------------------------
    !Calculate divergence at voronoi cells (hexagons)
    !   based on edge normal velocities
    !---------------------------------------------------------------
    type(grid_structure), intent(in) :: mesh
    type(scalar_field), intent(in):: u ! velocity at cell edges
    real(r8), allocatable, intent(inout) :: flux(:,:)
    integer(i4):: i, j, k, l, ed
    real(r8):: signcor

    !$omp parallel do &
    !$omp default(none) &
    !$omp shared(mesh, u, flux) &
    !$omp private(j, l, signcor) &
    !$omp schedule(static)
    do i=1,mesh%nv
      !Divergence of uh on unit sphere using div_cell_Cgrig
      !For all edges forming the hexagon
      flux(i,:)=0._r8
      do j=1, mesh%v(i)%nnb
        !Get edge index
        l=mesh%v(i)%ed(j)
        !Get edge outer normal related to the hexagon
        signcor=real(mesh%hx(i)%ttgout(j), r8)
        flux(i,j) = signcor*u%f(l)*mesh%edhx(l)%leng
      end do
    end do
    !$omp end parallel do

    return

  end subroutine flux_hx
end module highorder
