module transport
  !=============================================================================
  !  TRANSPORT module
  !
  !	Pack for several simulations on transport on deformational flow simulation
  !  on the sphere using Voronoi grids
  !
  !	Pedro da Silva Peixoto (pedrosp@ime.usp.br)
  !	Jan 2013
  !=============================================================================

  !Use global constants and kinds
  use constants, only: &
       datadir, &
       eps, &
       erad, &
       i4, &
       pardir, &
       pi, &
       pi2, &
       pio2, &
       r8, &
       rad2deg, &
       deg2rad, &
       unitspharea

  !Use main grid data structures
  use datastruct, only: &
       grid_structure, &
       scalar_field, &
       vectorinterpol_methods, &
       vector_field_cart

  !Use routines from the spherical mesh pack
  use smeshpack, only: &
       alignind, &
       alignindlimit, &
       arcdistll, &
       convert_vec_sph2cart, &
       error_norm_2, &
       error_norm_max, &
       error_norm_max_rel, &
       getunit, &
       norm, &
       sph2cart, &
       cart2sph, &
       proj_vec_sphere

  !Use interpolation routines
  use interpack, only: &
       getmultirecon, &
       gradcalc, &
       plot_cart_vectorfield, &
       plot_scalarfield, &
       precalcgradpol, &
       scalar_interpol, &
       vector_interpol, &
       vector_field_cart, &
       vrec_remap, &
       vector_reconstruct

  !Use differential operator routines
  use diffoperpack, only: &
       div_mesh_fullvec

  implicit none

  !Global variables

  !Flags
  integer (i4):: testcase !Test case
  integer (i4):: initialfield !Initial condition

  !Time variables
  real (r8):: w   !Angular speed
  real (r8):: dt  !Time step
  real (r8):: T   !Period

  !Trajectory Calculation Method
  ! 0 = Taylor Series or exact
  ! 1 =Ritchie
  integer (i4):: traj

  !Number of time steps and plots
  integer (i4):: ntime
  integer (i4):: nplots

  !Logical for plots or not
  logical:: plots

  !Discretizations to centroids
  logical:: discretcentroid

  !Plotsteps - will plot at every plotsteps timesteps
  integer (i4):: plotsteps

  !Name for files and kind of staggering
  character (len=128)::  transpname
  character (len=8)::  stag

  !Kind of interpolation for scalar and vector fields
  character (len=64):: ksinterpol
  character (len=64):: kvinterpol

  !Methods to be used for reconstruction and interpolation
  type(vectorinterpol_methods) :: recon_mtd

  !Inteprolation/Reconsgtruction method used for
  !     divergence/laplacian discretizations
  type(vectorinterpol_methods) :: discinterp_mtd

  !RBF shape parameter
  real (r8):: rbf_par

  !Alignment index cut off value
  real (r8) :: alignlimit

  !Percentage of aligned cells
  real (r8) :: alignpercent

  !RBF parameter
  real(r8) :: rbf_h !RBF shape parameter

  !Monotonicity flag
  integer (i4):: monot

  !Conservative flag
  integer (i4):: conserv

  !Paramter for variable rotation
  real(r8):: delta

contains 

  !======================================================================================
  !    TRANSPORT TESTS
  !======================================================================================

  subroutine gettransppars(mesh)
    !---------------------------------------------------
    ! gettransppars
    !    Reads transport test parameters from file named "trans.par"
    !    Saves parameters on global variables
    !
    ! Pedro Peixoto Jan 2012
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

    !Flag for reconstruction to cell centroids
    integer:: ireconcentroid

    !Flag for discretization to cell centroids
    integer:: idiscretcentroid

    !Alignement index type (value (val) or percentage (per))
    character (len=3):: aligntype

    !Alignement index threshold
    real (r8):: alignpar

    !Couters
    integer(i4)::i
    integer(i4)::j

    !Standard definition of the deformal tests
    testcase=1      !Non-div case 1
    initialfield=1  !1-Cosine bell ; 2-Gaussian ; 3-Slotted Cylinder
    plots=.true.   !Do plots or not
    T=5._r8      !Period
    w=pi/T     !Cicles per hour
    ntime= 15 * (2**(mesh%glevel-3)) !60  !Number of iterations for a whole cycle of T
    stag="HA"
    nplots=20
    traj=1 !Trajectory calculation - Ritchie
    delta=0. !Variable rotation speed

    !Standard parameters file
    filename=trim(pardir)//"trans.par"
    print*,"Transport parameters (file): ", trim(filename)
    print*
    call getunit(fileunit)

    !A parameters file must exist
    open(fileunit,file=filename,status='old')

    read(fileunit,*)  buffer
    read(fileunit,*)  buffer
    read(fileunit,*)  testcase
    read(fileunit,*)  buffer
    read(fileunit,*)  initialfield
    read(fileunit,*)  buffer
    read(fileunit,*)  traj
    read(fileunit,*)  buffer
    read(fileunit,*)  ntime, adjustntime
    read(fileunit,*)  buffer
    read(fileunit,*)  stag
    read(fileunit,*)  buffer
    read(fileunit,*)  monot
    read(fileunit,*)  buffer
    read(fileunit,*)  conserv
    read(fileunit,*)  buffer
    read(fileunit,*)  ksinterpol
    read(fileunit,*)  buffer
    read(fileunit,*)  recon_mtd%recon
    read(fileunit,*)  buffer
    read(fileunit,*)  recon_mtd%interp
    read(fileunit,*)  buffer
    read(fileunit,*)  recon_mtd%horecon
    read(fileunit,*)  buffer
    read(fileunit,*)  ireconcentroid
    read(fileunit,*)  buffer
    read(fileunit,*)  discinterp_mtd%interp
    read(fileunit,*)  buffer
    read(fileunit,*)  idiscretcentroid
    read(fileunit,*)  buffer
    read(fileunit,*)  alignpar, aligntype
    read(fileunit,*)  buffer
    read(fileunit,*)  rbf_par
    read(fileunit,*)  buffer
    read(fileunit,*)  nplots

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

    !Set number of times to plot
    if(nplots==0) then
       plots=.false.
    else
       plotsteps=ntime/nplots
    end if
    if(plotsteps==0)then
       !ntime too small - plot every timestep
       plotsteps=1
    end if

    !Time step
    dt=T/real(ntime, r8)

    print*, "Test Case (1-NonD1, 2-NonD2, 3-Div, 4-C2+SBRot)     : ", testcase
    print*, "Ini Field (1-cosbell, 2-twogauss, 3-cyl, 4-onegauss): ", initialfield
    print*, "Trajectory Calculation (0-Taylor 1-Ritchie)         : ", traj
    print*, "Staggering type                             : ", trim(stag)
    print*, "Plots? (nplots)                             : ", plots, nplots
    print*, "Total period                                : ", T
    print*, "Number of timesteps                         : ", ntime
    print*, "dt                                          : ",  dt
    print*, "Monotonicity?                               : ", monot
    print*, "Mass fixer?                                 : ", conserv

    !No higher order reconstruction method
    !    ("0" or "none" given)
    if(len(trim(recon_mtd%horecon))<2 .or. &
         trim(recon_mtd%horecon) == "none" )then
       recon_mtd%horecon=""
       recon_mtd%hyb=.false.
    else
       recon_mtd%hyb=.true.
    end if

    !Set vector interpolation name
    recon_mtd%name=trim(recon_mtd%recon)//&
         trim(recon_mtd%horecon)//trim(recon_mtd%interp)

    !Set reconstruction position (centroid?)
    if(ireconcentroid==1)then
       recon_mtd%massc=.true.
       recon_mtd%name=trim(recon_mtd%name)//"_b"
    else
       recon_mtd%massc=.false.
    end if
    recon_mtd%rbf_par=rbf_par

    !Adjust kind of interpolation to be used
    !call getmultirecon(kvinterpol,  recon_mtd)
    if(len(trim(recon_mtd%recon))>1)then
       print*, "Vector Reconstruction method      : ", trim(recon_mtd%recon)
    end if
    if(len(trim(recon_mtd%interp))>1)then
       print*, "Vector interpolation method       : ", trim(recon_mtd%interp)
    end if
    if(len(trim(recon_mtd%horecon))>1)then
       print*, "Vector higher order method        : ", trim(recon_mtd%horecon)
    end if
    if(len(trim(recon_mtd%name))>0)then
       print*, "Total vector interpolation method : ", trim(recon_mtd%name)
    end if
    if(recon_mtd%massc)then
       print*, "Reconstructing to mass centers"
    else
       print*, "Reconstructing to circumcenters or nodes"
    end if
    print*

    discinterp_mtd%hyb=.false.
    !Set interpolation/ reconstruction method for discrete div/lap
    if(len(trim(discinterp_mtd%interp))<=2 .or. &
         trim(discinterp_mtd%interp)=="none")then
       discinterp_mtd%recon=""
       discinterp_mtd%interp=""
       discinterp_mtd%horecon=""
    else
       if(trim(stag)=="HC".or.trim(stag)=="TC")then
          discinterp_mtd%recon=discinterp_mtd%interp
          discinterp_mtd%hyb=.true.
       end if
    end if
    !Set interpolation name for discretizations
    discinterp_mtd%name=trim(discinterp_mtd%interp)
    discinterp_mtd%rbf_par=rbf_par

    !Set discretization position (centroid?)
    if(idiscretcentroid==1)then
       discretcentroid=.true.
       if(trim(discinterp_mtd%name)=="")then
          discinterp_mtd%name=trim(discinterp_mtd%name)//"b"
       else
          discinterp_mtd%name=trim(discinterp_mtd%name)//"_b"
       end if
    else
       discretcentroid=.false.
    end if
    discinterp_mtd%massc=discretcentroid
    print*, "Interpolation for div mtd  : ", discinterp_mtd%interp
    if(discretcentroid)then
       print*,"Discretizing to mass centers"
    else
       print*,"Discretizing to nodes or circumcenters"
    end if

    !Set aligment index limit
    if(alignpar>=0 .and. recon_mtd%hyb)then
       if(trim(aligntype)=="per")then
          !Find index value relative to the percentage given
          alignlimit=alignindlimit(alignpar, mesh)
          alignpercent=alignpar
       else
          alignlimit=alignpar
          j=0
          do i=1, mesh%nv
             if(mesh%hx(i)%align<alignlimit+eps)then
                j=j+1
             end if
          end do
          alignpercent=real(j,r8)/real(mesh%nv,r8)
       end if
    else
       alignlimit=0
       alignpercent=1
    end if
    recon_mtd%alignlimit=alignlimit
    recon_mtd%alignpercent=alignpercent
    discinterp_mtd%alignlimit=alignlimit
    discinterp_mtd%alignpercent=alignpercent

    print*, "Alignment index threshold  : ", alignlimit
    print*, "Percentage of aligned cells: ", alignpercent

    !Set a standart name for files
    write(atmp,'(i8)') int(testcase)
    transpname="trans_v"//trim(adjustl(trim(atmp)))
    write(atmp,'(i8)') int(initialfield)
    transpname=trim(adjustl(trim(transpname)))//"_in"//trim(adjustl(trim(atmp)))
    write(atmp,'(i8)') int(traj)
    transpname=trim(adjustl(trim(transpname)))//"_tj"//trim(adjustl(trim(atmp)))
    write(atmp,'(i8)') int(ntime)
    transpname=trim(adjustl(trim(transpname)))//"_nt"//trim(adjustl(trim(atmp)))
    transpname=trim(adjustl(trim(transpname)))//"_"//trim(adjustl(trim(ksinterpol)))
    write(atmp,'(i8)') int(monot)
    transpname=trim(adjustl(trim(transpname)))//"_mt"//trim(adjustl(trim(atmp)))
    transpname=trim(adjustl(trim(transpname)))//"_"//trim(adjustl(trim(recon_mtd%name)))
    write(atmp,'(i8)') int(conserv)
    transpname=trim(adjustl(trim(transpname)))//"_mfx"//trim(adjustl(trim(atmp)))

    print*, "Transport Test Name for Plots: ", trim(transpname)
    print*

    return
  end subroutine gettransppars


  subroutine transptests(mesh)
    !-----------------------------------------
    !  Main transport tests routine
    !-----------------------------------------

    !Mesh
    type(grid_structure) :: mesh

    !Scalar fields for semi-lag
    type(scalar_field):: phi
    type(scalar_field):: phi_initial
    type(scalar_field):: phi_next
    type(scalar_field):: rho
    type(scalar_field):: rho_initial
    type(scalar_field):: rho_next
    type(scalar_field):: phi_rho
    type(scalar_field):: rho_div
    type(scalar_field):: div_vel
    type(scalar_field):: rho_dt2

    !General error field (multiuse) and aux field
    type(scalar_field):: errorfield
    type(scalar_field):: auxfield

    !Vector fields with velocities in cartesian coord
    type(vector_field_cart):: vel_ed  !edge
    type(vector_field_cart):: vel     !general
    type(vector_field_cart):: velerror     !general

    !Normal components of vector field, for edges
    type(scalar_field):: vn

    !File name for output
    character (len=256):: filename

    !Total Mass
    real(r8):: tmass_phi
    real(r8):: tmass_rho
    real(r8):: tmass

    !Total mass adjusted for conservation
    real(r8):: tmassc_rho
    real(r8):: tmassc !Total mass
    real(r8):: massc !Mass correction
    integer(i4):: massunit !Unit for output

    !Time vars
    real (r8):: time

    !Error variables
    real (r8):: error
    real (r8):: error2
    real (r8):: error2tmp
    real (r8):: errormax

    !Variables to follow trajectory of point
    !integer (i4) :: trajnode
    !real (r8) :: trajlat, trajlon, trajlonold, trajlatold

    !Semi-lag points
    real (r8):: dlat
    real (r8):: dlon
    real (r8):: pd(1:3)
    real (r8):: pdexact(1:3)

    !File units
    integer (i4):: errorsunit

    !Character vars
    character (len=256):: buffer
    character (len=64):: atime

    !Indexes
    integer (i4):: i
    integer (i4)::k
    logical::  ifile

    !Time counting variables
    real(r8):: elapsed_time
    real(r8):: start
    real(r8):: finish

    !Variaveis OpenMP
    !integer :: n_threads, id_thread, chunk,  OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM

    print*
    print*, "Transport Test...please be patient..."
    print*

    !Read parameters
    call gettransppars(mesh)

    !Check if exact scalar interpolation wanted
    if(trim(ksinterpol)=="exact" .and. testcase/=5 .and. testcase/=7)then
       print*, "Error on Transport Tests: I can only calculate the exact scalar field"
       print*, "   on the rotation test cases (5, 7)"
       stop
    end if

    if(testcase==7)then
       delta=pi2/(3.*T)
    else
       delta=0.
    end if

    !Set initial conditions and allocate space

    !Mixing ration
    phi%pos=0
    phi%n=mesh%nv
    allocate(phi%f(1:phi%n))

    !Density
    rho%pos=0
    rho%n=mesh%nv
    allocate(rho%f(1:phi%n))

    !Set scalar initial conditions
    do i=1,mesh%nv
       phi%f(i)=f(mesh%v(i)%lon,mesh%v(i)%lat)
       rho%f(i)=1._r8
    end do

    !Calculate gradients or polynomials of scalar field (for interpolation)
    call precalcgradpol(phi,mesh, ksinterpol)
    call precalcgradpol(rho,mesh, ksinterpol)

    !Setup initial velocity field for HC grids
    if(trim(stag)=="HC")then
       !Full velocity at hx edge midpoints
       vel_ed%pos=3
       vel_ed%n=mesh%ne
       allocate(vel_ed%p(1:vel_ed%n))

       !Normal components of velocities
       vn%pos=3
       vn%n=mesh%ne
       allocate(vn%f(1:vn%n))
    end if

    !Zero time counter
    k=0
    time=0.0_r8
    print*, "Initial time:", time

    !Full vectors on triangle vertices
    vel%pos=0  !Maybe chaged depending of reconstruction method
    call calc_vel(vel, time, mesh)
    velerror=vel

    !File for trajectory
    !filename=trim(datadir)//trim(transpname)//"_traj.gmt"
    !call getunit(trajunit)
    !open(trajunit,file=filename, status='replace')

    !Follow the trajectory of the node
    !call sph2cart(0._r8, pi/3._r8, p(1), p(2), p(3))
    !trajnode=getnearnode(p, mesh)
    !trajlon=0._r8
    !trajlat=pi/3._r8
    !Write trajectory points
    !write(trajunit, '(2f16.8)') trajlon*rad2deg, trajlat*rad2deg
    !write(*, *) k, time, trajlon*rad2deg, trajlat*rad2deg

    !Initialize other variables

    !Semi-lag variables
    phi_next=phi
    phi_initial=phi

    rho_next=rho
    rho_initial=rho
    rho_div=rho
    rho_dt2=rho

    phi_rho=phi

    errorfield=phi
    auxfield=phi

    errormax=0._r8
    error2=0._r8

    !Tracer density and mass calculations
    tmass=0.
    tmass_phi=0.
    tmass_rho=0.
    do i=1, mesh%nv
       phi_rho%f(i)=phi%f(i)*rho%f(i)
       tmass=tmass+phi_rho%f(i)*mesh%hx(i)%areag
       tmass_phi=tmass_phi+phi%f(i)*mesh%hx(i)%areag
       tmass_rho=tmass_rho+rho%f(i)*mesh%hx(i)%areag
    end do
    tmass=tmass/unitspharea
    tmass_phi=tmass_phi/unitspharea
    tmass_rho=tmass_rho/unitspharea

    !File for mass calculations
    filename=trim(datadir)//trim(transpname)//trim(mesh%name)//"_mass.txt"
    call getunit(massunit)
    open(massunit,file=filename, status='replace')
    write(massunit,'(a80)') "Time InitialmassRHO ActualMassRHO "//&
         "InitialmassPHIRHO ActualMassPHIRHO"
    write(massunit,'(f12.8,4f24.12)') time, &
         tmass_rho, tmass_rho, tmass, tmass

    !Pre calculate gradients or polinomials for interpolation
    call precalcgradpol(phi_rho, mesh, ksinterpol)

    !Plot initial conditions (if nplots>0)
    call plotfields()  !coment for DEBUG

    !Start counting time
    call cpu_time(start)

    !Start time iteration
    do k=1, ntime
       time=real(k, r8)*dt

       !Calculate velocities and divergence on grid
       !      at mid time step
       if(trim(stag)=="HA")then
          !Calculate velocities at nodes
          call calc_vel(vel, time-dt/2._r8, mesh)

          !Calculate divergence at mid time step
          ! for cell centers (nodes)
          div_vel=div_mesh_fullvec(vel, mesh)

       elseif(trim(stag)=="HC")then

          !Calculate full velocity at hx edge midpoints
          call calc_vel(vel_ed, time-dt/2._r8, mesh)

          !Calculate normal components of velocities
          !  The tangent velocities are thrown away, to simulate HC grid
          do i=1, vn%n
             !Calculate normal components of velocities in scalar variable
             vn%f(i)=dot_product(vel_ed%p(i)%v, mesh%edhx(i)%nr)

             !Calculate normal components of velocities in vector variable
             vel_ed%p(i)%v=vn%f(i)*mesh%edhx(i)%nr
          end do

          !Reconstruct velocity for either node or cell vertex
          !  Depending on the method this might change the position of vel%pos
          if(trim(recon_mtd%interp)/="exact")then
             call vrec_remap (vn, vel, mesh, recon_mtd)
             !Calculate errors on velocity
             !error=0._r8
             do i=1, mesh%nv
                !vtmp=vector_reconstruct (mesh%v(i)%p, vn, mesh, recon_mtd%recon)
                !error=max(0., norm(vtmp-velocity(mesh%v(i)%p, time)))
                velerror%p(i)%v=vel%p(i)%v-velocity(mesh%v(i)%p, time)
                velerror%p(i)%v=proj_vec_sphere(velerror%p(i)%v, mesh%v(i)%p)
                !print*, "Recon error:", i, error
             end do
          else !Exact vector calculation
             !Calculate exact velocities on nodes
             call calc_vel(vel, time-dt/2._r8, mesh)
          end if

          !Calculate divergence at mid time step
          ! for cell centers (nodes)
          div_vel=div_mesh_fullvec(vel_ed, mesh)

       end if

       !If exact wanted - force to be zero
       !  only work with non divergent field - obviously
       if(trim(recon_mtd%interp)=="exact")then
          div_vel%f(1:div_vel%n)=0._r8
          if(testcase==3 .and. k==1)then
             print*, "Cannot do exact calculation of divergent for a test case"
             print*, "   with divergence. Please hard code it.", testcase
             print*, "For now assuming zero!"
          end if
       end if
       !Calculate density * divergence
       rho_div%f = dt * rho_dt2%f * div_vel%f / 2._r8

       !Calculate field for departure point interpolation
       auxfield%f = rho%f - rho_div%f

       !Pre calculate polinomial for interpolation
       call precalcgradpol(auxfield, mesh, ksinterpol)

       error2tmp=0.

       !Calculate the departure point for each node in parallel

       !OPENMP PARALLEL REGION
       !$OMP PARALLEL  SHARED(phi, rho, phi_next, rho_next, &
       !$OMP       vel, errormax, error2tmp, mesh) &
       !$OMP DEFAULT(firstprivate)

       !Debug print for openmp
       !n_threads=omp_get_num_threads()
       !id_thread=omp_get_thread_num()
       !if(k==1)then
       !  print*, "Thread ",id_thread+1, " de ", n_threads, " comecou..."
       !end if

       !OPENMP PARALLEL DO REGION
       !$OMP DO SCHEDULE(auto) REDUCTION(max:errormax) &
       !$OMP       REDUCTION(+:error2tmp)
       do i=1, mesh%nv
          if(traj==1)then
             !Calculate departure point - semi-lag 2 time level
             pd=semilag2t_departpoint(mesh%v(i)%p, mesh, vel, dt, time)

             !Calculate semi-analytic or exact trajectory - just to get the error
             call departure_point(mesh%v(i)%lon, mesh%v(i)%lat, dlon, dlat, time, dt)
             call sph2cart(dlon, dlat, pdexact(1), pdexact(2), pdexact(3))

             error=norm(pd-pdexact)
             errorfield%f(i)=error
             errormax=max(errormax, error )
             error=error*error
             error2tmp=error2tmp+error

          else
             !Calculate semi-analytic or exact trajectory
             call departure_point(mesh%v(i)%lon, mesh%v(i)%lat, &
                  dlon, dlat, time, dt)
             call sph2cart(dlon, dlat, pd(1), pd(2), pd(3))
          end if


          if(trim(ksinterpol)=="exact")then
             !Calculate exact scalar value for rotation tests
             call cart2sph(pd(1), pd(2), pd(3), dlon, dlat)
             phi_next%f(i)=f(dlon, dlat, time-dt)
             !phi_next%f(i)=f(mesh%v(i)%lon, mesh%v(i)%lat, time)
             !phi_next%f(i)=f(dlon, dlat, time)
             !rho_next%f(i)=rho%f(i)
             rho_next%f(i)=1._r8-2.*rho_div%f(i)
          else !Interpolate scalar field
             phi_next%f(i)=scalar_interpol(pd, phi, mesh,  ksinterpol, monot=monot)
             rho_next%f(i)=scalar_interpol(pd, auxfield, mesh, ksinterpol, monot=monot)
             rho_next%f(i)=rho_next%f(i)-rho_div%f(i)
          end if


       end do
       !$OMP END DO
       !print*, "Thread ",id_thread+1, " de ", n_threads, " acabou!"
       !$OMP END PARALLEL

       !Trajectory error calculation
       error2tmp=dsqrt(error2tmp/real(mesh%nv, r8))
       error2=max(error2, error2tmp)

       !Calculate gradients for the next time step
       call precalcgradpol(phi_next, mesh, ksinterpol)

       !Calculate graidnets for next time step
       call precalcgradpol(rho_next, mesh, ksinterpol)

       !Uptade phi
       phi=phi_next

       !Update density
       rho=rho_next

       !Calculate total mass for unitary tracer
       tmassc_rho=0._r8
       do i=1, mesh%nv
          tmassc_rho=tmassc_rho+rho%f(i)*mesh%hx(i)%areag
       end do
       tmassc_rho=tmassc_rho/unitspharea

       !Conserve rho
       if(conserv==1)then
          massc=tmass_rho/tmassc_rho
          rho%f=massc*rho%f
       end if

       !OPENMP PARALLEL DO
       !$OMP PARALLEL  DO SHARED(phi, rho, phi_next, rho_next, &
       !$OMP       rho_dt2, phi_rho, mesh) &
       !$OMP       DEFAULT(firstprivate) &
       !$OMP       SCHEDULE(auto)
       do i=1, mesh%nv
          !Update the time extrapolation of density
          rho_dt2%f(i)=3.*rho_next%f(i)/2.-rho%f(i)/2.

          !Calculate the tracer density
          phi_rho%f(i)=phi%f(i)*rho%f(i)
       end do
       !$OMP END PARALLEL DO

       !Calculate tracer mass
       tmassc=0._r8
       do i=1, mesh%nv
          tmassc=tmassc+phi_rho%f(i)*mesh%hx(i)%areag
       end do
       tmassc=tmassc/unitspharea

       if(conserv==1)then
          massc=tmass/tmassc
          phi%f=massc*phi%f
       end if

       write(massunit,'(f12.8,6f24.12)') time, tmass_rho, tmassc_rho, tmass, tmassc

       !Calculate trajectory points using high order
       !trajlonold=trajlon
       !trajlatold=trajlat
       !Write trajectory points
       !call arrival_point(trajlonold, trajlatold, trajlon, trajlat, time-dt, dt)
       !write(trajunit, '(2f16.8)') trajlon*rad2deg, trajlat*rad2deg
       !write(*, *) k, time, trajlon*rad2deg, trajlat*rad2deg

       !Printable file at each plotsteps dtime
       print "(a8,i8,a12,i6,a12,i6,a8,f12.6)", "Mesh:",&
            mesh%nv, " Timesteps:", ntime, " k:", k, " Time:", time
       if(k==1 .or. k==ntime  .or. mod(k,plotsteps)==0 )then
          call plotfields()
       end if

       !Debug - do 1 time step only
       !if(k==ntime)then
       !exit !uncoment for DEBUG
       !end if
    end do
    !close(trajunit)
    close(massunit)

    call cpu_time(finish)
    elapsed_time= finish-start
    print*
    print '(a42, 2f12.6)', &
         " Transport test time (seg, seg/8): ", elapsed_time, elapsed_time/8
    print*

    !Plot errors
    call plot_errors()

    !File for errors
    filename=trim(datadir)//"transport_errors.txt"
    call getunit(errorsunit)

    inquire(file=filename, exist=ifile)
    if(ifile)then
       open(errorsunit,file=filename, status='old', position='append')
    else
       open(errorsunit,file=filename, status='replace')
       buffer="           n      mvdist case inif traj  nt       dt"//&
            "                scalar_interpol        vector_interpol      field"//&
            "     SLerror(inf)        SLerror(L2)            time      nplots  "//&
            "     Transpname"
       write(errorsunit, '(a)') buffer
    end if

    !Write trajectory errors
    call write_errors("TRAJ  ")
    print*, "Final trajectory error:", errormax, error2

    !Save final step error
    !PHI
    !error2=error_norm_2_rel(phi%f, phi_initial%f)
    error2=error_norm_2(phi%f, phi_initial%f, phi%n)
    errormax=error_norm_max_rel(phi%f, phi_initial%f, phi%n)
    errormax=error_norm_max(phi%f, phi_initial%f, phi%n)
    !field="PHI"
    call write_errors("PHI   ")  !coment for DEBUG
    print*, "PHI - Error (max, L2) :", errormax, error2

    !RHO
    !error2=error_norm_2_rel(rho%f, rho_initial%f)
    error2=error_norm_2(rho%f, rho_initial%f, rho%n)
    !errorsup=error_norm_max_rel(rho%f, rho_initial%f, rho%n)
    errormax=error_norm_max(rho%f, rho_initial%f, rho%n)
    !field="RHO"
    call write_errors("RHO   ") !coment for DEBUG
    print*, "RHO - Error (max, L2) :", errormax, error2

    !PHI*RHO
    !error2=error_norm_2_rel(phi_rho%f, phi_initial%f)
    error2=error_norm_2(phi_rho%f, phi_initial%f, phi%n)
    !errorsup=error_norm_max_rel(phi_rho%f, phi_initial%f, phi_rho%n)
    errormax=error_norm_max(phi_rho%f, phi_initial%f, phi_rho%n)
    !field="PHIRHO"
    call write_errors("PHIRHO") !coment for DEBUG
    print*, "PHI*RHO - Error (max, L2) :", errormax, error2

    if(plots)then
       print*
       print*, "Use ./gmt/anim.sh to see the fluid transporting"
       print*, "Use ./gmt/plot.sh to see the vector field, "//&
            " initial condition and errors"
    end if

    print*

    return

  contains

    subroutine plotfields()
      !Plot scalar fields to files
      if(plots)then

         !Time index
         write(atime,'(i8)') int(k)

         !Scalar field plots
         phi%name=trim(transpname)//"_phi_t"//trim(adjustl(trim(atime)))
         call plot_scalarfield(phi, mesh)

         rho%name=trim(transpname)//"_rho_t"//trim(adjustl(trim(atime)))
         call plot_scalarfield(rho, mesh)

         phi_rho%name=trim(transpname)//"_phi_rho_t"//trim(adjustl(trim(atime)))
         call plot_scalarfield(phi_rho, mesh)

         errorfield%name=trim(transpname)//"_traj_t"//trim(adjustl(trim(atime)))
         call plot_scalarfield(errorfield, mesh)

         !Velocity field
         !print*, "debug in plotfields: ", vel%pos
         vel%name=trim(transpname)//"_vel_t"//trim(adjustl(trim(atime)))
         call plot_cart_vectorfield(vel, mesh)

         !Error for velocity field
         !print*, "debug in plotfields: ", vel%pos
         velerror%name=trim(transpname)//"_veler_t"//trim(adjustl(trim(atime)))
         call plot_cart_vectorfield(velerror, mesh)
      end if

    end subroutine plotfields

    subroutine plot_errors()
      if(plots)then
         write(atime,'(i8)') int(ntime)

         !Calculate error at the end

         !Semi-lag
         errorfield%f=phi%f-phi_initial%f
         errorfield%name=trim(transpname)//"_phi_error_t"//trim(adjustl(trim(atime)))
         call plot_scalarfield(errorfield, mesh)

         errorfield%f=rho%f-rho_initial%f
         errorfield%name=trim(transpname)//"_rho_error_t"//trim(adjustl(trim(atime)))
         call plot_scalarfield(errorfield, mesh)

         errorfield%f=phi_rho%f-phi_initial%f
         errorfield%name=trim(transpname)//"_phi_rho_error_t"//trim(adjustl(trim(atime)))
         call plot_scalarfield(errorfield, mesh)

      end if
    end subroutine plot_errors

    subroutine write_errors(field)

      !Auxiliar string for field name
      character (len=6):: field

      !Write errors to file
      write(errorsunit, '(i12, 1f12.4, 3i4, i8, 1e16.8, 2a24, a10, 2e18.8, f18.6, i6, a60)') &
           mesh%nv,  mesh%meanvdist*rad2deg, testcase, &
           initialfield, traj, ntime, dt, trim(ksinterpol), &
           trim(recon_mtd%name), trim(field), errormax, error2, elapsed_time, nplots, &
           transpname

    end subroutine write_errors

  end subroutine transptests


  !---------------------------------
  !Calculates velocity field on nodes at a time step
  !  vel must have the position already set
  !---------------------------------
  subroutine calc_vel(vel, time, mesh)

    !Vector field with velocities at the  time
    type(vector_field_cart):: vel

    !Time
    real (r8):: time

    !Aux
    integer (i4):: j

    !Mesh
    type(grid_structure), intent(in) :: mesh

    !Auxiliar vars
    real (r8):: utmp
    real (r8):: vtmp

    select case(vel%pos)
    case(0) !Vectors on nodes
       if(.not. allocated(vel%p))then
          vel%n=mesh%nv
          allocate(vel%p(1:vel%n))
       end if
       do j=1, mesh%nv
          utmp=u(mesh%v(j)%lon,mesh%v(j)%lat, time)
          vtmp=v(mesh%v(j)%lon,mesh%v(j)%lat, time)
          call convert_vec_sph2cart(utmp, vtmp, mesh%v(j)%p, vel%p(j)%v)
       end do
    case(1) !Vectors on triangle circumcenters
       if(.not. allocated(vel%p))then
          vel%n=mesh%nt
          allocate(vel%p(1:vel%n))
       end if
       do j=1, mesh%nt
          utmp=u(mesh%tr(j)%c%lon,mesh%tr(j)%c%lat, time)
          vtmp=v(mesh%tr(j)%c%lon,mesh%tr(j)%c%lat, time)
          call convert_vec_sph2cart(utmp, vtmp, mesh%tr(j)%c%p, vel%p(j)%v)
       end do
    case(3) !Vectors on HX edges
       if(.not. allocated(vel%p))then
          vel%n=mesh%ne
          allocate(vel%p(1:vel%n))
       end if
       do j=1, mesh%ne
          utmp=u(mesh%edhx(j)%c%lon,mesh%edhx(j)%c%lat, time)
          vtmp=v(mesh%edhx(j)%c%lon,mesh%edhx(j)%c%lat, time)
          call convert_vec_sph2cart(utmp, vtmp, mesh%edhx(j)%c%p, vel%p(j)%v)
       end do
    case default
       print*, "calc_vel error: Please set a correct position for the velocity : ", vel%pos
       stop
    end select
    return

  end subroutine calc_vel

  !-------------------------------------------------------------------------------
  ! Semi-analytic trajectory calculations
  !-------------------------------------------------------------------------------

  !High order departure point calculation
  subroutine departure_point(alon, alat, dlon, dlat, time, dt)
    !Arrival point, time and dt
    real (r8), intent(in) :: alat
    real (r8), intent(in) :: alon
    real (r8), intent(in) :: time
    real (r8), intent(in) :: dt
    !Departure point
    real (r8), intent(out) :: dlat
    real (r8), intent(out) :: dlon
    !Intermidiate points
    real (r8):: lat
    real (r8):: lon
    !Intermediate time step
    real (r8):: tau
    !Index
    integer (i4):: k
    integer (i4):: m


    if(testcase==5 .or. testcase==7)then
       !In this case the exact trajectory is known
       m=1
       tau=dt
    else
       !In this case use the 3rd order traj approx
       !  with m substeps
       !m=10
       m=100
       tau=dt/real(m,r8)
    end if
    lon=alon
    lat=alat
    do k=1, m
       call departure_point_semiantraj(lon, lat, dlon, dlat, time-(k-1)*tau, tau)
       lon=dlon
       lat=dlat
    end do

    return

  end subroutine departure_point

  !Substeps of the high order method
  subroutine departure_point_semiantraj(lon, lat, dlon, dlat, time, dt)
    real (r8), intent(in) :: lat
    real (r8), intent(in) :: lon
    real (r8), intent(in) :: time
    real (r8), intent(in) :: dt
    real (r8), intent(out) :: dlat
    real (r8), intent(out) :: dlon
    real (r8):: u_tilde
    real (r8):: v_tilde
    real (r8):: sinlon
    real (r8):: coslon
    real (r8):: sinlon2
    real (r8):: coslon2
    real (r8):: coslat
    real (r8):: sinlat
    real (r8):: sinlonp
    real (r8):: coslonp
    real (r8):: k
    real (r8):: lonp

    sinlon=dsin(lon+pi)
    coslon=dcos(lon+pi)
    sinlat=dsin(lat)
    coslat=dcos(lat)

    if(abs(coslat)>eps)then
       u_tilde=u(lon, lat, time)/dcos(lat)
    else
       u_tilde=0
    end if
    v_tilde=v(lon, lat, time)

    select case(testcase)
    case(1)
       k=2.4
       sinlon2=dsin((lon+pi)/2)
       coslon2=dcos((lon+pi)/2)
       dlon=lon-(dt)*u_tilde-(dt**2)*k*sinlon2*(sinlon2*sinlat*dsin(w*time)*w &
            -u_tilde*sinlat*dcos(w*time)*coslon2 &
            -v_tilde*sinlon2*dcos(lat)*dcos(w*time))
       dlat=lat-dt*v_tilde-(k*(dt**2)/4.)*(sinlon*coslat*dsin(w*time)*w &
            -u_tilde*coslon*coslat*dcos(w*time)+v_tilde*sinlon*sinlat*dcos(w*time))
       !print*, lon*180/pi, lat*180/pi, dlon*180/pi, dlat*180/pi
    case(2)
       k=2.
       dlon=lon-(dt)*u_tilde-(dt**2)*k*sinlon*(sinlon*sinlat*dsin(w*time) &
            -2.*u_tilde*sinlat*dcos(w*time)*coslon &
            -v_tilde*sinlon*coslat*dcos(w*time))
       dlat=lat-dt*v_tilde-(k*(dt**2))*(sinlon*coslat*dsin(w*time)*coslon*w &
            -2.*u_tilde*(coslon**2)*coslat*dcos(w*time) &
            +u_tilde*coslat*dcos(w*time) &
            +v_tilde*sinlon*sinlat*coslon*dcos(w*time))
       !print*, lon*180/pi, lat*180/pi, dlon*180/pi, dlat*180/pi
    case(3) !Divergent case 3
       k=1.
       sinlon2=dsin((lon+pi)/2)
       coslon2=dcos((lon+pi)/2)
       dlon=lon-(dt)*u_tilde-(dt**2)*k*sinlon2*dcos(lat)*(sinlon2*dsin(w*time)*w*dcos(lat)*sinlat &
            -u_tilde*coslat*dcos(w*time)*coslon2*dsin(lat) &
            -3.*v_tilde*sinlon2*(coslat**2)*dcos(w*time)+2.*v_tilde*sinlon2*dcos(w*time))
       dlat=lat-dt*v_tilde-(k*(dt**2)/4.)*(dcos(lat)**2)*(sinlon*coslat*dsin(w*time)*w &
            -u_tilde*coslon*coslat*dcos(w*time) &
            +3.*v_tilde*sinlon*sinlat*dcos(w*time))
    case(4)
       k=2.
       lonp=lon-2._r8*pi*time/T
       sinlonp=dsin(lonp+pi)
       coslonp=dcos(lonp+pi)

       !u_tildep=u_tilde/dcos(lat)
       dlon=lon-(dt)*u_tilde-(dt**2)*k*sinlonp*(sinlonp*sinlat*dsin(w*time) &
            -2._r8*u_tilde*sinlat*dcos(w*time)*coslonp &
            -v_tilde*sinlonp*coslat*dcos(w*time))
       dlat=lat-dt*v_tilde-(k*(dt**2))*(sinlonp*coslat*dsin(w*time)*coslonp*w &
            -2._r8*u_tilde*(coslonp**2)*coslat*dcos(w*time) &
            +u_tilde*coslat*dcos(w*time) &
            +v_tilde*sinlonp*sinlat*coslonp*dcos(w*time))
       !print*, lon*180/pi, lat*180/pi, dlon*180/pi, dlat*180/pi
    case(5, 7) !Passive Advection - Rotation
       !dlon=lon-dt*(2._r8)*pi/T
       dlon=lon-(rotangle(time)-rotangle(time-dt)) !dt*(2._r8)*pi/T
       dlat=lat
    end select

    !print*, lon*rad2deg, dlon*rad2deg
    !Correct lat-lon frame
    if(dlon<=-pi)then
       dlon=dlon+pi2
    elseif(dlon>=pi)then
       dlon=dlon-pi2
    end if
    !This is an artifitial correction
    if(dlat<-pio2)then
       dlat=-pio2
    elseif(dlat>pio2)then
       dlat=pio2
    end if

    return
  end subroutine departure_point_semiantraj

  !High order arrival point calculation

  subroutine arrival_point(alon, alat, dlon, dlat, time, dt)
    !Arrival point, time and dt
    real (r8), intent(in) :: alat
    real (r8), intent(in) :: alon
    real (r8), intent(in) :: time
    real (r8), intent(in) :: dt
    !Departure point
    real (r8), intent(out) :: dlat
    real (r8), intent(out) :: dlon
    !Intermidiate points
    real (r8):: lat
    real (r8):: lon
    !Intermediate time step
    real (r8):: tau
    !Index
    integer (i4):: k
    integer (i4):: m

    m=3
    tau=dt/real(m,r8)
    lon=alon
    lat=alat
    do k=1, m
       call arrival_point_semiantraj(lon, lat, dlon, dlat, time+(k-1)*tau, tau)
       lon=dlon
       lat=dlat
    end do

    return

  end subroutine arrival_point

  subroutine arrival_point_semiantraj(lon, lat, dlon, dlat, time, dt)
    !Departure point
    real (r8), intent(in) :: lat
    real (r8), intent(in) :: lon
    real (r8), intent(in) :: time
    real (r8), intent(in) :: dt
    !Destination/arrival point
    real (r8), intent(out) :: dlat
    real (r8), intent(out) :: dlon
    real (r8):: u_tilde
    real (r8):: v_tilde
    real (r8):: sinlon
    real (r8):: coslon
    real (r8):: sinlon2
    real (r8):: coslon2
    real (r8):: coslat
    real (r8):: sinlat
    real (r8):: coslonp
    real (r8):: sinlonp
    real (r8):: k
    real (r8):: lonp

    sinlon=dsin(lon+pi)
    coslon=dcos(lon+pi)
    sinlat=dsin(lat)
    coslat=dcos(lat)

    if(abs(coslat)>eps)then
       u_tilde=u(lon, lat, time)/dcos(lat)
    else
       u_tilde=0
    end if
    v_tilde=v(lon, lat, time)

    select case(testcase)
    case(1)
       k=2.4
       sinlon2=dsin((lon+pi)/2)
       coslon2=dcos((lon+pi)/2)
       dlon=lon+(dt)*u_tilde-(dt**2)*k*sinlon2*(sinlon2*sinlat*dsin(w*time)*w &
            -u_tilde*sinlat*dcos(w*time)*coslon2 &
            -v_tilde*sinlon2*dcos(lat)*dcos(w*time))
       dlat=lat+dt*v_tilde-(k*(dt**2)/4)*(sinlon*coslat*dsin(w*time)*w &
            -u_tilde*coslon*coslat*dcos(w*time)+v_tilde*sinlon*sinlat*dcos(w*time))
       !print*, lon*180/pi, lat*180/pi, dlon*180/pi, dlat*180/pi
    case(2)
       k=2.
       dlon=lon+(dt)*u_tilde-(dt**2)*k*sinlon*(sinlon*sinlat*dsin(w*time) &
            -2*u_tilde*sinlat*dcos(w*time)*coslon &
            -v_tilde*sinlon*coslat*dcos(w*time))
       dlat=lat+dt*v_tilde-(k*(dt**2))*(sinlon*coslat*dsin(w*time)*coslon*w &
            -2*u_tilde*(coslon**2)*coslat*dcos(w*time) &
            +u_tilde*coslat*dcos(w*time) &
            +v_tilde*sinlon*sinlat*coslon*dcos(w*time))
       !print*, lon*180/pi, lat*180/pi, dlon*180/pi, dlat*180/pi
    case(3) !Divergent case 3
       k=1.
       sinlon2=dsin((lon+pi)/2)
       coslon2=dcos((lon+pi)/2)
       dlon=lon+(dt)*u_tilde-(dt**2)*k*sinlon2*dcos(lat)*(sinlon2*dsin(w*time)*w*dcos(lat)*sinlat &
            -u_tilde*coslat*dcos(w*time)*coslon2*dsin(lat) &
            -3.*v_tilde*sinlon2*(coslat**2)*dcos(w*time)+2.*v_tilde*sinlon2*dcos(w*time))
       dlat=lat+dt*v_tilde-(k*(dt**2)/4.)*(dcos(lat)**2)*(sinlon*coslat*dsin(w*time)*w &
            -u_tilde*coslon*coslat*dcos(w*time) &
            +3.*v_tilde*sinlon*sinlat*dcos(w*time))
    case(4)
       k=2.
       lonp=lon-2*pi*time/T
       sinlonp=dsin(lonp+pi)
       coslonp=dcos(lonp+pi)
       dlon=lon+(dt)*u_tilde-(dt**2)*k*sinlonp*(sinlonp*sinlat*dsin(w*time) &
            -2*u_tilde*sinlat*dcos(w*time)*coslonp &
            -v_tilde*sinlonp*coslat*dcos(w*time))
       dlat=lat+dt*v_tilde-(k*(dt**2))*(sinlonp*coslat*dsin(w*time)*coslonp*w &
            -2*u_tilde*(coslonp**2)*coslat*dcos(w*time) &
            +u_tilde*coslat*dcos(w*time) &
            +v_tilde*sinlonp*sinlat*coslonp*dcos(w*time))
       !print*, lon*180/pi, lat*180/pi, dlon*180/pi, dlat*180/pi
    end select
    return
  end subroutine arrival_point_semiantraj

  !========================================================================
  !    Semi-Lagrangian Routines
  !========================================================================

  function semilag2t_departpoint(pa, mesh, velm, dt, time)
    !-----------------------------------------------
    !  SEMILAG_DEPARTPOINT
    !   Calculates the departure point of a semilagrangean
    !   discretization of time as in Ritche 1987 and Garcia 2001
    !   with 2 time step levels
    !-----------------------------------------------

    !Point of arrival
    real (r8), intent(in) :: pa(1:3)

    !Vector field with velocities at the mean time step (t+dt/2)
    type(vector_field_cart), intent(in) :: velm

    !Mesh in which the variable belongs
    type(grid_structure), intent(in) :: mesh

    !Time step
    real (r8), intent(in) :: dt
    real (r8), intent(in) :: time

    !Medium point velocity
    real (r8):: vm(1:3)

    !Medium point
    real (r8):: pm(1:3)
    real (r8):: pmold(1:3)

    !Departure point
    real (r8):: semilag2t_departpoint(1:3)

    !Aux
    real(r8):: pmmove

    !Index
    !integer, optional:: node
    integer:: i
    integer:: n

    !Tmp vars
    !real (r8):: lon, lat

    n=2
    !n=1
    !n=slagiter !read from trans.par
    pm=pa
    i=0
    !vm=vector_interpol(pm, velm, mesh)
    !if(node==17)then
    !     call cart2sph(pm(1), pm(2), pm(3), lon, lat)
    !     print*, node, i, lon*rad2deg, lat*rad2deg
    !end if

    do i=1, n

       !Interpolate the velocity at medium point
       if(trim(recon_mtd%interp)=="exact")then
          vm=velocity(pm, time-dt/2._r8)
       else
          !if(i==1)then
          !vm=velocity(pm, time)
          !else
          vm=vector_interpol(pm, velm, mesh, recon_mtd%interp)

          !end if
          !vm=velocity(pm, time)
          !vm=velm%p(node)%v
          !vmtmp=velocity(pm, time)
          !print*, norm(vm-vmtmp)
       end if

       !Save point
       pmold=pm

       !Update point
       pm=pa-(dt/2._r8)*vm

       !Project to sphere
       pm=pm/norm(pm)

       !if(node==17)then
       !
       !  call cart2sph(pm(1), pm(2), pm(3), lon, lat)
       !  print*, node, i, lon*rad2deg, lat*rad2deg
       !  print*, "    Error Vel:", norm(vm-velocity(pmold, time))
       !end if

       !Check convergence
       pmmove=norm(pm-pmold)
       !print*, i, pmmove
       !if(pmmove<eps/100000._r8)then
       !print*, "Semi-Lag Traj. Stopping here. No movement."
       !   exit
       !end if
    end do

    semilag2t_departpoint=2._r8*(dot_product(pm, pa))*pm-pa

    return
  end function semilag2t_departpoint

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
          !lonp=lon-time*2.*pi/T
          lonp=lon-rotangle(time)
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
       b=0.1_r8
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
       f=2+dcos(3*latp)**3*dsin(5*lonp)**3
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
    !-----------------------------------
    !  U - velocity in West-East direction
    !   for a given testcase
    !
    !   Lat in [-pi/2,pi/2], Lon in [-pi,pi]
    !
    !  P. Peixoto - Feb2012
    !---------------------------------------------
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
       !u=2.*pi*dcos(lat)/T
       u=rotspeed(time)*dcos(lat)

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
    !-----------------------------------
    !  V - velocity in South-North direction
    !   for a given testcase
    !
    !   Lat in [-pi/2,pi/2], Lon in [-pi,pi]
    !
    !  P. Peixoto - Feb2012
    !---------------------------------------------
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

  function rotangle(time)
    ! ------------------------------------------
    ! ROT ANGLE
    ! Function for variable rotation velocity
    ! Returns the angle of rotation for a given time
    !----------------------------------------------

    !Period, actual time
    real (r8)::  time

    !Function value
    real (r8)::  rotangle
    real (r8):: tp

    tp = pi2*time/T
    rotangle= tp + delta*dsin(tp)*time
    !rotangle = tp !+ delta*(1.+time+50.*time**2+time**3+100000.*time**4) !dsin(tp)*time

    return
  end function rotangle

  function rotspeed(time)
    !--------------------------------------------------
    !  ROTSPEED
    ! Derivative of function for variable rotation velocity (rotangle)
    !----------------------------------------------

    !Period, actual time
    real (r8)::  time
    !Funcation value
    real (r8)::  rotspeed
    real (r8):: tp

    tp= pi2*time/T

    rotspeed= pi2/T + delta*(tp*dcos(tp)+dsin(tp))
    !rotspeed= pi2/T + delta*(-(4.*pi2/T)*dsin(tp)*dcos(tp)**3)
    ! rotspeed= pi2/T !+ delta*(1.+100.*time+3.*time**2+400000.*time**3)

    return
  end function rotspeed

  !===========================================================================
  !    PASSIVE ADVECTION
  !===========================================================================
  subroutine passive_advection(mesh)
    !-------------------------------------------------
    ! Main passive advection routine
    !------------------------------------------------

    !Mesh
    type(grid_structure), intent(in) :: mesh

    !Function/initial condition on the sphere
    type(scalar_field):: phi
    type(scalar_field):: phi_initial
    type(scalar_field):: phi_next
    type(scalar_field):: phi_error

    !Arrival and departure point (temporary variables)
    real (r8):: alat
    real (r8):: alon
    real (r8):: dlat
    real (r8):: dlon
    real (r8):: pd(1:3)

    !Reference lat, lon for rotated coordinates
    real (r8):: latref
    real (r8):: lonref

    !Angular speed, time step, max time, actual time
    real (r8):: w
    real (r8):: dt
    real (r8):: T
    real (r8):: time

    !Auxiliar variables
    !real (r8) :: error
    integer (i4):: i
    integer (i4)::k
    integer (i4):: ntime
    character (len=60)::  atime

    print*
    print*, "Passive Advection Test...please be patient..."
    print*

    !    alat=89*mesh%deg2rad
    !    alon=-110*mesh%deg2rad
    latref=-45*deg2rad !Rotation through pole
    !latref=0.*deg2rad !Simples rotation
    lonref=0.*deg2rad
    !T=24._r8*20._r8     !Hours
    T=64._r8     !Hours
    w=pi2/T      !cicles per hour
    dt=2._r8            !Hours
    ntime= int ( (1._r8)*T/dt , i4 )
    print*, "Number of time steps:", ntime

    !Set initial conditions
    time=0
    phi%pos=0
    phi%n=mesh%nv
    allocate(phi%f(1:phi%n))
    do i=1,mesh%nv
       phi%f(i)=g(mesh%v(i)%lon,mesh%v(i)%lat)
    end do
    call gradcalc(phi,mesh)

    !Plot initial condition
    write(atime,'(i8)') int(time)
    phi%name="passadv_t"//trim(adjustl(trim(atime)))
    print*, "Calculated time:", time, " hours "
    call plot_scalarfield(phi, mesh)

    phi_next=phi
    phi_initial=phi
    phi_error=phi
    do k=1,ntime
       time=time+dt
       do i=1,mesh%nv
          alat=mesh%v(i)%lat
          alon=mesh%v(i)%lon

          !Assume known trajectory
          call departure_point(alon, alat, dlon, dlat, lonref, latref, w, dt)
          call sph2cart(dlon, dlat, pd(1), pd(2), pd(3))
          !Interpol at departure point
          phi_next%f(i)=scalar_interpol(pd, phi, mesh, "hermtrv ")
       end do
       call gradcalc(phi_next, mesh)
       phi=phi_next
       !Printable file at each day
       !if(mod(k,ntime/40)==0.or.k==ntime)then
       if(mod(k,1)==0.or.k==ntime)then
          write(atime,'(i8)') int(time)
          phi%name="passadv_t"//trim(adjustl(trim(atime)))
          print*, "Calculated time:", time, " hours "
          call plot_scalarfield(phi, mesh)
       end if
       !Print error after 1 and 2 full rotations (20 and 40 days)
       if(k==ntime.or.k==ntime)then
          phi_error%f=phi%f-phi_initial%f
          call gradcalc(phi_error, mesh)
          write(atime,'(i8)') int(time)
          phi_error%name="passadv_error"//trim(adjustl(trim(atime)))
          call plot_scalarfield(phi_error, mesh)
       end if
    end do
    print*
    print*, "Use ./anim.sh passadv 'mesh' to see the fluid rotating"
    print*, "Use ./plot.sh to see the errors"
    print*
    return

  contains

    !---------------------------------
    !Initial conditions for passive advection
    !---------------------------------
    function g(lon,lat)
      real (r8), intent(in) :: lat !Radians
      real (r8), intent(in) :: lon !Radians
      real (r8):: g
      real (r8):: L
      real (r8):: r


      r = (lat-0)**2+(lon-0)**2
      L=5*pi/180._r8

      g=100*exp(-r/L)
      !f=dcos(lat) !5-sqrt(sigma)+2*(sin(lat)+sin(lon))

      return

    end function g

    !Departure point calculation
    subroutine departure_point(alon, alat, dlon, dlat, lonref, latref, w, dt)

      real (r8), intent(in)  :: alat    !Arrival position
      real (r8), intent(in)  :: alon    !Arrival position
      real (r8), intent(out) :: dlat    !Departure position
      real (r8), intent(out) :: dlon    !Departure position
      real (r8), intent(in)  :: w       !Angular speed
      real (r8), intent(in)  :: dt      !Time step
      real (r8), intent(in)  :: latref     !Reference for rotation
      real (r8), intent(in)  :: lonref     !Reference for rotation

      real (r8):: alatrot      !Rotated arrival position
      real (r8):: alonrot      !Rotated arrival position
      real (r8):: dlatrot      !Rotated departure position
      real (r8):: dlonrot      !Rotated departure position

      call convert_latlon2rotated(alon, alat, alonrot, alatrot, lonref, latref)
      dlatrot=alatrot
      dlonrot=alonrot-dt*w
      call convert_rotated2latlon(dlonrot, dlatrot, dlon, dlat, lonref, latref)

      return
    end subroutine departure_point


  end subroutine passive_advection

  subroutine convert_latlon2rotated(lon, lat, lonrot, latrot, lonref, latref)
    !-------------------------------------------------------------------
    ! CONVERT_LATLON2ROTATED
    !
    !   Converts latitude and longitude coordinates to rotated
    !   coordinates. All values must be in radians.
    !   Lat in [-pi/2,pi/2], Lon in [-pi,pi]
    !
    !   Based on Giraldo (1999) -Trajectory Calc. for Spherical Geodesic Grids
    !
    !   Pedro Peixoto Jan 2011
    !---------------------------------------------------------------------
    real (r8), intent(in)  :: lat     !Original latitude longitude
    real (r8), intent(in)  :: lon     !Original latitude longitude
    real (r8), intent(in)  :: latref  !Reference latitude longitude
    real (r8), intent(in)  :: lonref  !Reference latitude longitude
    real (r8), intent(out) :: latrot  !Rotated latitude longitude
    real (r8), intent(out) :: lonrot  !Rotated latitude longitude

    real (r8):: sinlat
    real (r8):: sinlatref
    real (r8):: coslat
    real (r8):: coslatref
    real (r8):: coslon_ref
    real (r8):: sinlon_ref

    sinlat=dsin(lat)
    sinlatref=dsin(latref)
    coslat=dcos(lat)
    coslatref=dcos(latref)
    coslon_ref=dcos(lon-lonref)
    sinlon_ref=dsin(lon-lonref)

    latrot=dasin(sinlat*coslatref-coslat*sinlatref*coslon_ref)
    lonrot=datan2(coslat*sinlon_ref, coslatref*coslon_ref*coslat+sinlatref*sinlat)
    if(lonrot>pi)then
       lonrot=lonrot-2*pi
    elseif(lonrot<-pi)then
       lonrot=lonrot+2*pi
    end if
    return
  end subroutine convert_latlon2rotated

  subroutine convert_rotated2latlon(lonrot, latrot, lon, lat, lonref, latref)
    !-------------------------------------------------------------------
    ! CONVERT_ROTATED2LATLON
    !
    !   Converts rotated latitude and longitude coordinates to original
    !   coordinates. All values must be in radians.
    !   Lat in [-pi/2,pi/2], Lon in [-pi,pi]
    !
    !   Based on Giraldo (1999) -Trajectory Calc. for Spherical Geodesic Grids
    !
    !   Pedro Peixoto Jan 2011
    !---------------------------------------------------------------------
    real (r8), intent(in) :: latrot !Rotated latitude longitude
    real (r8), intent(in) :: lonrot !Rotated latitude longitude
    real (r8), intent(in) :: latref !Reference latitude longitude
    real (r8), intent(in) :: lonref !Reference latitude longitude
    real (r8), intent(out) :: lat   !Original latitude longitude
    real (r8), intent(out) :: lon   !Original latitude longitude

    real (r8):: sinlatrot
    real (r8):: sinlatref
    real (r8):: coslatrot
    real (r8):: coslatref
    real (r8)::  tanlatrot
    real (r8):: sinlonrot
    real (r8):: coslonrot

    sinlatrot=dsin(latrot)
    sinlatref=dsin(latref)
    coslatrot=dcos(latrot)
    coslatref=dcos(latref)
    tanlatrot=dtan(latrot)
    sinlonrot=dsin(lonrot)
    coslonrot=dcos(lonrot)

    lat=dasin(coslatref*sinlatrot+coslatrot*sinlatref*coslonrot)
    lon=datan2(sinlonrot*coslatrot, (coslatrot*coslonrot*coslatref-sinlatref*sinlatrot) )+lonref

    if(lon>pi)then
       lon=lon-2*pi
    elseif(lon<-pi)then
       lon=lon+2*pi
    end if
    return
  end subroutine convert_rotated2latlon


end module transport
