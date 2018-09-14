module poisson
  !====================================================================
  ! Main module for superconveregence analysis of Poisson equation
  ! Pedro Peixoto Oct 2018
  ! partially based on multigrid solver of
  ! Marline I. Silva - Sept 2014
  !=====================================================================

  !Global constants
  use constants, only: &
    datadir, &
    deg2rad, &
    eps, &
    erad, &
    i4, &
    pardir, &
    pi, &
    r8, &
    rad2deg


  !Data structures
  use datastruct, only:  &
    grid_structure, &
    scalar_field, &
    vector, &
    vector_field_cart

  !Spherical mesh routines
  use smeshpack, only: &
    alignind, &
    alignindlimit, &
    arcdistll, &
    arcintersec, &
    arclen, &
    cart2sph, &
    choleskydecomp, &
    choleskysolve, &
    convert_vec_sph2cart, &
    convijll, &
    convllij, &
    distortionhx, &
    error_norm_2, &
    error_norm_max, &
    gcircarcintersec, &
    getnearnode, &
    getnearnodes, &
    gettr, &
    getunit, &
    hxtr_intersec_areas, &
    modint, &
    norm, &
    ortogonalarc, &
    proj_vec_sphere, &
    sph2cart, &
    sphpolarea, &
    sqtriintersec, &
    vorbarycenter



  !Interpolation pack
  use interpack, only: &
    plot_scalarfield, &
    plot_cart_vectorfield, &
    plot_grad_vectorfield, &
    plot_scalarfield_sphgrid, &
    plot_vectorfield_sphgrid

  implicit none	

  !Function for test
  integer (i4):: testfunc
  integer (i4):: testcase

  !Logical for other tests
  logical:: discretcentroid !Discretizations to centroids
  logical:: testgrad !Test gradient
  logical:: plots !Perform plots

  !Name for interpolation files
  character (len=128):: simulname
  character (len=128):: testcasename

  !SOR solver parameters
  real (r8):: w         !Relaxation parameter
  real (r8)::nrm_res    !Norm residue
  real (r8)::TOL        !Tolerance
  integer (i4):: numit  !max num of iterations


contains

  subroutine getsimulpars(mesh)
    !---------------------------------------------------
    ! GETSIMULPARS
    !    Reads simulation parameters from file named "poisson.par"
    !    Saves parameters on global variables
    !
    ! Pedro Peixoto Oct 2018
    !--------------------------------------------------

    !Mesh
    type(grid_structure) :: mesh

    !Buffer for strings
    character (len=300):: buffer

    !Flag for plots
    integer:: iplots

    !Flag for discretization to cell centroids
    integer:: idiscretcentroid

    !Parameters file
    character (len=256):: filename

    !Unit for input file
    integer (i4):: fileunit
    !integer :: error

    !Couters
    integer(i4)::i
    integer(i4)::j
    character (len=16):: atmp

    !Standard definition of the interpolation tests
    testfunc=1      !Non-div case 1
    plots=.true.   !Do plots or not

    w=1                   !Relaxation parameter
    tol=0.00000001        !Tolerance
    nrm_res=100           !Norm residue
    numit=20              !Maximum number of iterations

    !Standard interpolation parameters file
    filename=trim(pardir)//"poisson.par"

    print*,"Simulation parameters (file): ", trim(filename)
    print*
    call getunit(fileunit)

    !A parameters file already must exist
    open(fileunit,file=filename,status='old')

    read(fileunit,*)  buffer
    read(fileunit,*)  buffer
    read(fileunit,*)  testcase
    read(fileunit,*)  buffer
    read(fileunit,*)  testfunc
    read(fileunit,*)  buffer
    read(fileunit,*)  idiscretcentroid
    read(fileunit,*)  buffer
    read(fileunit,*)  iplots
    read(fileunit,*)  buffer
    read(fileunit,*)  w, numit, tol

    close(fileunit)

    !Set ploting parameter
    if(iplots==1)then
      plots=.true.
    else
      plots=.false.
    end if


    !Set discretization position (centroid?)
    if(idiscretcentroid==1)then
      discretcentroid=.true.
    else
      discretcentroid=.false.
    end if



    simulname="poisson"

    !Names depending on test cases
    select case(testcase)
      case(1) !Truncation errors
        simulname=trim(adjustl(trim(simulname)))//"_laptrunc"
        testcasename="laptrunc"
      case(2) !Poisson error
        simulname=trim(adjustl(trim(simulname)))//"_poiser"
        testcasename="poiser"
      case default
        simulname="_"
    end select

    !Test function (vector or scalar field)
    write(atmp,'(i8)') int(testfunc)
    simulname=trim(adjustl(trim(simulname)))//"_f"//trim(adjustl(trim(atmp)))

    print*, "Test case             : ", testcase, testcasename
    print*, "Function used         : ", testfunc
    print*, "Plots? (nplots)       : ", plots
    if(trim(testcasename)=="poiser")then
      print*, "Solver parameters (w,n,tol)   :", w, numit, tol
    end if

    print*

    print*, "Simulation Name for Plots: ", trim(simulname)
    print*

    return
  end subroutine getsimulpars

  subroutine poisson_main(mesh)
    !---------------------------------------------------
    ! poisson_main
    !    Main routine for test cases related to poisson superconv. analysis
    !
    ! Pedro Peixoto Oct 2018
    !--------------------------------------------------
    !Mesh
    type(grid_structure) :: mesh

    call getsimulpars(mesh)

    !Names depending on test cases
    select case(testcase)
      case(1) !Truncation errors
        call laplacian_truncation(mesh)
      case(2) !Poisson errors
        call poisson_error(mesh)
      case default
        print*, "ERROR on poisson_main: Don't know this testcase:", testcase
        stop
    end select

  end subroutine poisson_main


  !-----------------------------------------------------------------------------------------
  !    FUNCTIONS FOR TESTS
  !-----------------------------------------------------------------------------------------

  function f(p)
    !-----------------------------------
    !  F - initial conditions for scalar fields
    !
    !   p is a point in cartesian coords
    !
    !  P. Peixoto - Feb2012
    !---------------------------------------------
    real (r8), intent(in) :: p(1:3)
    real (r8):: lon
    real (r8):: lat

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
    real (r8), dimension(1:3) :: p1
    real (r8), dimension(1:3) :: p2
    real (r8):: b0

    !Center of bell, slot, hill
    real (r8):: lat1
    real (r8):: lon1
    real (r8):: lat2
    real (r8):: lon2

    !Wavelength of test function
    real (r8):: m
    real (r8):: n

    m=1
    n=1

    !convert to spherical coords
    call cart2sph ( p(1), p(2), p(3), lon, lat )

    if(testfunc>0 .and. testfunc <4)then
      !Deformational test functions
      lon1=-pi/(4._r8)
      lat1=0._r8
      lon2=pi/(4._r8)
      lat2=0._r8
    end if

    !scalar fields
    select case(testfunc)
      case(1) !Cosine bell
        b=0.1_r8
        c=1._r8 !0.9
        r=1._r8/2._r8
        r1= arcdistll(lon, lat, lon1, lat1)
        if(r1<r)then
          h1=(1._r8/2._r8)*(1._r8+dcos(pi*r1/r))
          f=b+c*h1
          return
        end if
        r2= arcdistll(lon, lat, lon2, lat2)
        if(r2<r)then
          h2=(1._r8/2._r8)*(1._r8+dcos(pi*r2/r))
          f=b+c*h2
          return
        end if
        f=b

      case(2) !Gaussian
        b0=5
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

        r1= arcdistll(lon, lat, lon1, lat1)
        if(r1<=r.and.abs(lon-lon1)>=r/(6._r8))then
          f=c
          return
        end if
        r2= arcdistll(lon, lat, lon2, lat2)
        if(r2<=r.and.abs(lon-lon2)>=r/(6._r8))then
          f=c
          return
        end if
        if(r1<=r.and.abs(lon-lon1) < r/(6._r8) .and. &
          (lat-lat1)<-(5._r8)*r/(12._r8))then
          f=c
          return
        end if
        if(r2<=r.and.abs(lon-lon2)<r/(6._r8).and. &
          (lat-lat2)>(5._r8)*r/(12._r8))then
          f=c
          return
        end if

      case(4) !Smooth Trigonometric
        f=((dcos(lat))**3)*((dsin(lon))**2)

      case(5) !Exponential
        f= (exp(p(1))+2*exp(p(2)+p(3)))/10

      case(6) !Linear in R3
        f=(1.0+2*p(1)+3*p(2)+4*p(3))/6.0

      case(7) !Constant 1
        f=1

      case(8) !Trigonom oscilation (beta of HR95)
        f=dcos(m * lon) * dcos(n * lat) ** 4
         !((dcos(lat))**4)*((dsin(lon))**7)

      case default
        print*, "F error: unknown scalar field (testfunc)", testfunc
        stop
    end select

    return
  end function f


  function df(p)
    !-------------------------------------------
    !  Grad of function "f"
    ! used to test interpolations
    ! in cartesian coords
    ! Returns the gradient in cartesian coords
    !------------------------------------------
    real (r8):: df(1:3)
    real (r8):: p(1:3)
    real (r8):: vlat
    real (r8):: vlon
    real (r8):: lat
    real (r8):: lon


    !Wavelength of test function
    real (r8):: m
    real (r8):: n

    m=1
    n=1

    call cart2sph (p(1), p(2), p(3), lon, lat)
    df=0.

    select case(testfunc)
      case(1) !Cosine bell

      case(2) !Gaussian

      case(3) !Non-smooth- Slotted cylinder
        df=0
      case(4) !Smooth Trigonometric
        !f=((dcos(lat))**3)*((dsin(lon))**2)
        vlon=2*dsin(lon)*dcos(lon)*((dcos(lat))**2)
        vlat=-3*((dcos(lat))**2)*((dsin(lon))**2)*dsin(lat)
        call convert_vec_sph2cart(vlon, vlat, p, df)
      case(5) !Exponential
        !f= (exp(p(1))+2*exp(p(2)+p(3)))/10
        df(1)=exp(p(1))/10
        df(2)=exp(p(2)+p(3))/5
        df(3)=exp(p(2)+p(3))/5
      case(6) !Linear in R3
        !f=(1.0+2*p(1)+3*p(2)+4*p(3))/6.0
        df(1)=2./6.
        df(2)=3./6.
        df(3)=4./6.
      case(7) !Constant
        df=0
      case(8) !Oscilating trig (beta of HR95)
        vlon = -dsin(m*lon)*m*dcos(n*lat)**4 / dcos(lat)
        vlat = -4._r8*dcos(m*lon)*dcos(n*lat)**3 * dsin(n*lat)* n
        !vlon=7*(dsin(lon)**6)*dcos(lon)*((dcos(lat))**3)
        !vlat=-4*((dcos(lat))**3)*((dsin(lon))**7)*dsin(lat)
        call convert_vec_sph2cart(vlon, vlat, p, df)
      case default
        print*, "DF error: unknown gradient for field (testfunc)", testfunc
        stop
    end select

    !Project vector to the sphere (if not already)
    df=proj_vec_sphere(df, p)

    return
  end function df

  function lap_exact(p)
    !-----------------------------------
    !  Laplacian of scalar field given in f()
    !
    !   p is a point in cartesian coords
    !---------------------------------------------
    real (r8), intent(in) :: p(1:3)
    real (r8):: lon
    real (r8):: lat

    real (r8):: lap_exact

    !Wavelength of test function
    real (r8):: m
    real (r8):: n

    !convert to spherical coords
    call cart2sph ( p(1), p(2), p(3), lon, lat )

    !Laplacian
    select case(testfunc)
      case(4) !Smooth Trigonometric
        lap_exact=dcos(lat) * (-5.0_r8 * dcos(lon) ** 2 + 7.0_r8 &
          - 12.0_r8 * dcos(lat)** 2 * dsin(lon) ** 2)
      case(7) !Constant 1
        lap_exact=0
      case (8) !Oscilat trig (beta of HR95)
        m=1
        n=1
        lap_exact = -dcos(m*lon)*dcos(n*lat)**2 * (m**2*dcos(n * lat) **2 &
          - 4._r8 * dsin(lat) * dcos(n*lat)*dsin(n*lat)*n*dcos(lat) &
          - 12._r8 * dcos(lat)** 2 * n** 2 + 16._r8*dcos(lat)**2 *dcos(n* lat)**2* n** 2) &
          / dcos(lat) ** 2
         !lap_exact= (dcos(lat)**2) * (dsin(lon) **5 )* &
         !    (33._r8 * dcos(lon) ** 2 + 9.0_r8 &
         !    - 20.0_r8 * dcos(lat)** 2 * dsin(lon) ** 2)
      case default
        print*, "Laplacian (lap_exact) error: unknown scalar field (testfunc) :", testfunc
        stop
    end select

    return
  end function lap_exact

  function func(p)
    real (r8)::p(1:3)
    real(r8)::func
    real(r8)::lat
    real(r8)::lon
    real(r8)::m
    real(r8)::n

    m=1
    n=1

    !TESTE 1
    !func=(1+2*p(1)+3*p(2)+4*p(3))/6._r8

    !TESTE 2
    !func = (exp(p(1))+ 2*(exp(p(2)+p(3))))/10

    call cart2sph (p(1), p(2), p(3), lon, lat)
    !TESTE 3
    !func=((dcos(lat))**3)*((dsin(lon))**2)

    !TESTE 4
    func=dcos(m*lon)*dcos(n*lat)**4

    !TESTE 5
    !func=dsin(lat)

    !TESTE 6
    !func=dsin(3*lon)*((dcos(3*lat)))**4

    !TESTE 7
    !func=-(1/2)*((dcos(lat))**3)*dcos(3*lon)

    !TESTE 8
    !func=dcos(lat)*dcos(lon)

    return
  end function func

  !----------------------------------------------------------------
  !     Exact Laplacian of the preceding functions
  !----------------------------------------------------------------

  function lap_ext(p)
    real (r8)::p(1:3)
    real (r8):: lat
    real (r8)::lon
    real (r8)::m
    real (r8)::n
    real(r8)::lap_ext

    m=1
    n=1

    !Converte para coordenadas esféricas
    call cart2sph (p(1), p(2), p(3), lon, lat)

    !TESTE 1 
    !lap_ext= (((3*dsin(lon)+2*dcos(lon))*(1-dcos(lat)**2+dsin(lat)**2))/6._r8*dcos(lat)) &
    !- 4._r8*dsin(lat)/3._r8

    !TESTE 2
    !lap_ext=(1.0_r8/10.0_r8*dcos(lat))*(2.0_r8*exp(dcos(lat)*dsin(lon)+dsin(lat))*(dcos(lat) &
    !*((dcos(lon))**2)-dsin(lon))+exp(dcos(lat)*dcos(lon))*(dcos(lat)*((dsin(lon))**2) &
    !-dcos(lon))-2.0_r8*dsin(lat)*exp(dcos(lat)*dsin(lon)+dsin(lat))*(dcos(lat)-dsin(lon) &
    !*dsin(lat))+dcos(lon)*((dsin(lat))**2)*exp(dcos(lat)*dcos(lon))) +0.1*(exp(dcos(lat) &
    !*dcos(lon))*(((dcos(lon))**2)*((dsin(lat))**2)-dcos(lon)*dcos(lat))+2.0_r8*exp(dcos(lat) &
    !*dsin(lon)+dsin(lat))*((dcos(lat)-dsin(lon)*(dsin(lat)))**2)-(dsin(lon)*dcos(lat) &
    !+dsin(lat)))

    !TESTE 3 (-lap u + u=f)
    !lap_ext=-dcos(lat)*(-5.0_r8*dcos(lon)**2+7.0_r8-12.0_r8*dcos(lat)**2*dsin(lon)**2)  &
    !    +((dcos(lat))**3)*((dsin(lon))**2)

    !-------------------------------------------------------------
    !TESTE 4
    ! lap_ext = -dcos(m*lon)*dcos(n*lat)**2 * (m**2*dcos(n * lat) **2 &
    !    - 4._r8 * dsin(lat) * dcos(n*lat)*dsin(n*lat)*n*dcos(lat) &
    !  -12._r8 * dcos(lat)**2*n**2+16._r8*dcos(lat)**2 *dcos(n* lat)**2* n**2)&
    !/dcos(lat)**2


    !TESTE 4 (-lap u +u =f)
    lap_ext = -(-dcos(m*lon)*dcos(n*lat)**2 * (m**2*dcos(n * lat) **2 &
      - 4._r8 * dsin(lat) * dcos(n*lat)*dsin(n*lat)*n*dcos(lat) &
      -12._r8 * dcos(lat)**2*n**2+16._r8*dcos(lat)**2 *dcos(n* lat)**2* n**2)&
      /dcos(lat)**2)   +  dcos(m*lon)*dcos(n*lat)**4

    !------------------------------------------------------------
    !TESTE 5
    !lap_ext=-2*dsin(lat)

    !TESTE 5 (-lap u+u=f)
    !lap_ext=2*dsin(lat)+dsin(lat)
    !----------------------------------------------------

    !TESTE 6
    !lap_ext=(-9*dsin(3*lon)*((dcos(3*lat))**4)/dcos(lat)**2) &
    ! -6*dsin(3*lon)*(dcos(lat))**2*(1-2*dcos(2*lat))**2*(2*dcos(2*lat) &
    !-2*dcos(4*lat)+13*dcos(6*lat)-7)   


    !TESTE 7
    !lap_ext=0!6*(dcos(lat)**3)*(dcos(3*lon))


    !TESTE 8
    !lap_ext=0.0

    return
  end function lap_ext




  !======================================================================
  !    LAPLACIAN TESTS
  !======================================================================
  subroutine laplacian_truncation(mesh)

    !Mesh
    type(grid_structure) :: mesh


    !Scalar fields
    type(scalar_field):: func
    type(scalar_field):: lap_ex
    type(scalar_field):: lap
    type(scalar_field):: error
    type(scalar_field):: ergrad
    type(scalar_field):: eddif

    !Counters
    integer:: i
    integer:: j
    integer:: k

    !Errors
    logical:: ifile
    integer:: iunit
    real (r8)::  error2
    real (r8):: errormax
    real (r8):: errorgrad2
    real (r8):: errorgradsup
    real (r8):: maxeddif
    character (len=256):: filename

    !Aux
    real(r8):: p(1:3)
    real(r8):: gradex
    real(r8):: gradest

    print*
    print*,"Laplacian Truncation Analysis "
    print*


    !Scalars on hexagon centers (nodes) - Exact
    lap_ex%pos=0
    lap_ex%n=mesh%nv
    lap_ex%name="lap_ex"
    allocate(lap_ex%f(1:lap_ex%n))
    do i=1, mesh%nv
      !Laplacian on node
      p=mesh%v(i)%p
      !lap_ex%f(i)=lap_exact(p)

      !Laplacian on barycenter
      !p=mesh%hx(i)%b%p
      lap_ex%f(i)=lap_exact(p)
       !print "(i8, 4f16.8)", i, lap_ex%f(i), p
    end do

    !Edge displacement
    eddif%n=mesh%nv
    eddif%name="edgedif"
    eddif%pos=0
    allocate(eddif%f(1:mesh%nv))
    do i=1, mesh%nv
      eddif%f(i)=0._r8
      do j=1,mesh%v(i)%nnb
        k=mesh%v(i)%ed(j)
        eddif%f(i)=max(eddif%f(i), (arclen(mesh%ed(k)%c%p, mesh%edhx(k)%c%p)/ &
          mesh%edhx(k)%leng))
      end do
       !eddif%f(i)=eddif%f(i)/mesh%v(i)%nnb
    end do

    !Edge intersection maximum difference
    maxeddif=maxval(abs(eddif%f(1:eddif%n)))

    !Gradient error
    ergrad%n=mesh%nv
    ergrad%name="ergrad"
    ergrad%pos=0
    allocate(ergrad%f(1:mesh%nv))
    do i=1, mesh%nv
      ergrad%f(i)=0._r8
    end do

    !Laplacian error
    error%n=mesh%nv
    error%name="erlap"
    error%pos=0
    allocate(error%f(1:mesh%nv))

    !Test function - Scalars on hexagon centers (nodes)
    func%pos=0
    func%n=mesh%nv
    func%name="func"
    allocate(func%f(1:func%n))
    do i=1, mesh%nv
      !Function on node
      p=mesh%v(i)%p
      func%f(i)=f(p)
    end do

    !Calculate Numerical Laplacian
    lap%pos=0
    lap%n=mesh%nv
    lap%name="lap_est"
    allocate(lap%f(1:lap%n))
    do i=1, mesh%nv
      lap%f(i)=0._r8
      ergrad%f(i)=0._r8
      do j=1, mesh%v(i)%nnb
        !Edge index
        k=mesh%v(i)%ed(j)
        !Hexagonal edge midpoint
        p=mesh%edhx(k)%c%p
        !ExactGrad=ExactGradVector*NormalVectorEdge*CorrectionOutHx
        gradex=dot_product(df(p),mesh%edhx(k)%nr)*mesh%hx(i)%nr(j)
        !Estimated Gradient Normal component
        gradest=(func%f(mesh%v(i)%nb(j))-func%f(i))/mesh%ed(k)%leng
        !arclen(p, mesh%v(i)%p) + &
        !arclen(p, mesh%v(mesh%v(i)%nb(j))%p)

        !Maximum gradient error for cell
        ergrad%f(i)=max(abs(gradest-gradex), ergrad%f(i))
        !Updade Laplacian
        lap%f(i)=lap%f(i)+gradest*mesh%edhx(k)%leng
      end do
      lap%f(i)=lap%f(i)/mesh%hx(i)%areag
      !Error in laplacian
      error%f(i)=lap_ex%f(i)-lap%f(i)
       !print"(i4, 3f16.8)", i, lap%f(i), lap_ex%f(i), error%f(i)
    end do

    !Global Errors

    error2=error_norm_2(lap%f, lap_ex%f, error%n)
    errormax=error_norm_max(lap%f, lap_ex%f, error%n)
    print*, "Error Lap (max, L2): ", errormax, error2

    !Global Errors for gradients
    errorgrad2=dsqrt(dot_product(ergrad%f,ergrad%f)/ergrad%n)
    errorgradsup=maxval(abs(ergrad%f(1:ergrad%n)))
    print*, "Error GRAD (max, L2): ", errorgradsup, errorgrad2
    print*

    ! Save error estimates
    !-------------------------------------------------------

    !For hexagonal methods
    filename=trim(datadir)//"lap_errors_"//trim(mesh%kind)// &
      "_"//trim(mesh%pos)//"_"//trim(mesh%optm)//".txt"
    call getunit(iunit)

    inquire(file=filename, exist=ifile)
    if(ifile)then
      open(iunit,file=filename, status='old', position='append')
    else
      open(iunit,file=filename, status='replace')
      write(iunit, '(a)') &
        " n    distance testfunc errormax error2 errorgradsup errorgrad2 maxeddif"
    end if

    write(iunit, '(i8, f18.8, i8, 5f22.12)') mesh%nv, mesh%meanvdist*rad2deg, &
      testfunc, errormax, error2, errorgradsup, errorgrad2, maxeddif

    close(iunit)

    ! Plot fields
    !-------------------------
    if(plots) then
      print*, "Plotting variables ... "

      lap_ex%name=trim(simulname)//"_exact"
      call plot_scalarfield(lap_ex, mesh)

      lap%name=trim(simulname)//"_est"
      call plot_scalarfield(lap, mesh)

      error%name=trim(simulname)//"_error"
      call plot_scalarfield(error, mesh)

      eddif%name=trim(simulname)//"_eddif"
      call plot_scalarfield(eddif, mesh)

      ergrad%name=trim(simulname)//"_ergrad"
      call plot_scalarfield(ergrad, mesh)

    end if

    return

  end subroutine laplacian_truncation

  !======================================================================
  !    POISSON TESTS
  !======================================================================
  subroutine poisson_error(mesh)

    !Mesh
    type(grid_structure) :: mesh


    !Scalar fields
    type(scalar_field):: func
    type(scalar_field):: lap_ex
    type(scalar_field):: lap
    type(scalar_field):: error
    type(scalar_field):: ergrad
    type(scalar_field):: eddif

    !Counters
    integer:: i
    integer:: j
    integer:: k

    !Errors
    logical:: ifile
    integer:: iunit
    real (r8)::  error2
    real (r8):: errormax
    real (r8):: errorgrad2
    real (r8):: errorgradsup
    real (r8):: maxeddif
    character (len=256):: filename

    !Aux
    real(r8):: p(1:3)
    real(r8):: gradex
    real(r8):: gradest

    print*
    print*,"Poisson Error Analysis "
    print*


    !Scalars on hexagon centers (nodes) - Exact
    lap_ex%pos=0
    lap_ex%n=mesh%nv
    lap_ex%name="lap_ex"
    allocate(lap_ex%f(1:lap_ex%n))
    do i=1, mesh%nv
      !Laplacian on node
      p=mesh%v(i)%p
      !lap_ex%f(i)=lap_exact(p)

      !Laplacian on barycenter
      !p=mesh%hx(i)%b%p
      lap_ex%f(i)=lap_exact(p)
       !print "(i8, 4f16.8)", i, lap_ex%f(i), p
    end do

    !Poisson solution error
    error%n=mesh%nv
    error%name="erpois"
    error%pos=0
    allocate(error%f(1:mesh%nv))

    !Test function - Scalars on hexagon centers (nodes)
    func%pos=0
    func%n=mesh%nv
    func%name="func"
    allocate(func%f(1:func%n))
    do i=1, mesh%nv
      !Function on node
      p=mesh%v(i)%p
      func%f(i)=f(p)
    end do

    !Calculate Numerical Laplacian
    lap%pos=0
    lap%n=mesh%nv
    lap%name="lap_est"
    allocate(lap%f(1:lap%n))
    do i=1, mesh%nv
      lap%f(i)=0._r8
      do j=1, mesh%v(i)%nnb
        !Edge index
        k=mesh%v(i)%ed(j)
        !Hexagonal edge midpoint
        p=mesh%edhx(k)%c%p


        !arclen(p, mesh%v(i)%p) + &
        !arclen(p, mesh%v(mesh%v(i)%nb(j))%p)

        !Updade Laplacian
        !lap%f(i)=lap%f(i)+gradest*mesh%edhx(k)%leng
      end do
      lap%f(i)=lap%f(i)/mesh%hx(i)%areag
      !Error in laplacian
      error%f(i)=lap_ex%f(i)-lap%f(i)
       !print"(i4, 3f16.8)", i, lap%f(i), lap_ex%f(i), error%f(i)
    end do

    !Global Errors
    error2=error_norm_2(lap%f, lap_ex%f, error%n)
    errormax=error_norm_max(lap%f, lap_ex%f, error%n)
    print*, "Error Lap (max, L2): ", errormax, error2


    ! Save error estimates
    !-------------------------------------------------------

    !For hexagonal methods
    filename=trim(datadir)//"lap_errors_"//trim(mesh%kind)// &
      "_"//trim(mesh%pos)//"_"//trim(mesh%optm)//".txt"
    call getunit(iunit)

    inquire(file=filename, exist=ifile)
    if(ifile)then
      open(iunit,file=filename, status='old', position='append')
    else
      open(iunit,file=filename, status='replace')
      write(iunit, '(a)') &
        " n    distance  testfunc errormax error2"
    end if

    write(iunit, '(i8, f18.8, i8, 5f22.12)') mesh%nv, mesh%meanvdist*rad2deg, &
      testfunc, errormax, error2

    close(iunit)

    ! Plot fields
    !-------------------------
    if(plots) then
      print*, "Plotting variables ... "

      lap_ex%name=trim(simulname)//"_exact"
      call plot_scalarfield(lap_ex, mesh)

      lap%name=trim(simulname)//"_est"
      call plot_scalarfield(lap, mesh)

      error%name=trim(simulname)//"_error"
      call plot_scalarfield(error, mesh)

    end if

    return

  end subroutine poisson_error

  subroutine sor(mesh, g, fap)
    !------------------------------------------------------------
    !RELAXAC_MULT_CELL
    ! Relaxation Scheme - Weighted Gauss-Seidel
    !
    !  return: "fap" that approximates u solution of
    !         for equation -Lap(u)=g
    !------------------------------------------------------------
    type(grid_structure)::mesh  !Mesh
    type(scalar_field)::fap     !Approximate value of the function on each grid point
    type(scalar_field)::g       !Independent term of equation

    !Counters
    integer::k
    integer::i
    integer::j

    real (r8)::ledhx     !Edge length
    real (r8)::lchx       !Distance between two neighboring centers
    real (r8)::laptmp
    real (r8)::laptmp1


    !For each iteration
    do k=1,numit

      !For each grid point
      do i=1,mesh%nv

        laptmp=0_r8
        laptmp1=0_r8
        ledhx=0_r8
        Lchx=0_r8

        !For each grid neighboring of i
        do j=1,mesh%v(i)%nnb

          !Edge length
          ledhx=mesh%edhx(mesh%v(i)%ed(j))%leng

          !Distance between two neighboring
          lchx=mesh%ed(mesh%v(i)%ed(j))%leng

          laptmp1=laptmp1+(ledhx/Lchx)
          laptmp=laptmp+(ledhx/Lchx)*(fap%f(mesh%v(i)%nb(j)))
        enddo

        laptmp=laptmp/mesh%hx(i)%areag
        laptmp1=laptmp1/mesh%hx(i)%areag

        !para o problema lap f =g
        !fap%f(i)=(1-w)*fap%f(i)+w*((laptmp-g%f(i))/(laptmp1))

        !para o problema -lap f+f =g
        fap%f(i)=(1-w)*(fap%f(i))+w*((laptmp+g%f(i))/(laptmp1+1))

      enddo

    enddo

    return
  end subroutine sor


  subroutine relaxac(mesh)
    !--------------------------------------------------------
    !RELAXAC test
    !
    ! Test for relaxation of the problem using the 
    !  method of Jacobi.
    ! For this test we use the subroutine "relaxac_mult_cell"
    !---------------------------------------------------------

    type(grid_structure)::mesh               !Mesh
    type(scalar_field)::fexact      !Exact value of the function on each grid point
    type(scalar_field)::fap         !Approximate value of the function on each grid point
    type(scalar_field)::fap_aux         !Approximate value of the function on each grid point
    type(scalar_field)::resid       !Residue
    type(scalar_field)::faperro     !Error between the exact value and the approximate value of the function
    type(scalar_field)::g           !Independent term of equation
    type(scalar_field)::lap_ap      !Approximate value of the laplacian on each grid point
    !type(interpolable_meshes)::fap_aux       !Auxiliary vector

    !Couter
    integer::i

    !Norm errors on each grid level
    real(r8)::nrm2
    real(r8)::nrmmax
    real(r8)::nrmres


    !A parameters file must exist
    character(len=60)::filename
    character (len=300):: buffer
    integer::fileunit
    integer::iunit


    !Approximate value of the function on each grid point
    allocate(fap%f(1:mesh%nv))
    allocate(fap_aux%f(1:mesh%nv))

    !Exact value of the function on each grid point
    allocate(fexact%f(1:mesh%nv))

    !Residue
    allocate(resid%f(1:mesh%nv))

    !Independent term of equation
    allocate(g%f(1:mesh%nv))


    !Approximate value of the laplacian on each grid point
    allocate(lap_ap%f(1:mesh%nv))

    !Error between the exact value and the approximate value of the function
    allocate(faperro%f(1:mesh%nv))


    !For each grid point
    do i=1,mesh%nv
      fap%f(i)=0!inic(mesh%v(i)%p)                        !initial condition
      !initial condition
      fexact%f(i)=func(mesh%v(i)%p)    !Exact value
      g%f(i)=lap_ext(mesh%v(i)%p)      !Independent term
    enddo



    !Calculate approximate value of the function on each grid point
    call relaxac_mult_cell_GS(mesh,g,fap,numit,w)

    !Calculation of the residue after each iteration
    call lap_cell(mesh,lap_ap,fap)
    resid%f=g%f-lap_ap%f

    !Norm residue
    nrm_res=error_nrm_max(resid%f,mesh%nv)


    !Norm errors on each grid level

    nrm2=error_norm_2(fexact%f,fap%f,mesh%nv)
    nrmmax=error_norm_max(fexact%f,fap%f,mesh%nv)
    nrmres=nrm_res


    !Deallocate
    deallocate(fap%f,fexact%f,faperro%f,g%f,resid%f,lap_ap%f)


    !Write values on file
    filename=trim(datadir)//"NRM_ERRO_RELAXAC.txt"
    call getunit ( iunit )
    open(iunit, file=filename, status='replace')
    write(iunit,*) "nrm2, nrmmax"

    write(iunit,"(3f32.16)") nrm2, nrmmax

    close (iunit)



    return
  end subroutine relaxac


  !==========================================================================
  ! Subroutines for relaxation, interpolation and transfer residue
  !==========================================================================

  subroutine lap_cell(mesh,lapc,fap)
    !------------------------------------------------------------
    !LAP_CELL
    !calculation of the Laplacian each point of mesh.
    !------------------------------------------------------------

    type(grid_structure)::mesh            !mesh
    !value function
    type(scalar_field)::fap
    !Approximate value of the laplacian on each grid point
    type(scalar_field)::lapc

    !counter
    integer::i
    integer::j
    integer::k

    !number points of mesh
    integer::n

    real (r8):: ledhx    !Edge length
    real (r8)::Lchx      !Distance between two neighboring
    real (r8)::laptmp


    !number points of the mesh
    n=mesh%nv

    !For each grid point
    do  i=1,n

      laptmp=0_r8
      ledhx=0_r8
      Lchx=0_r8

      !For each grid neighboring of i
      do j=1, mesh%v(i)%nnb

        !Edge index
        k=mesh%v(i)%ed(j)

        !Edge length
        ledhx=mesh%edhx(k)%leng

        !Distance between two neighboring
        Lchx=mesh%ed(k)%leng

        laptmp=laptmp+(ledhx/Lchx)*(fap%f(mesh%v(i)%nb(j))-fap%f(i))

      enddo


      !para o problema -lap f+f =g
      lapc%f(i)=-(laptmp/mesh%hx(i)%areag)+fap%f(i)

       !para o problema lap f =g
       !lapc%f(i)=laptmp

    enddo
    return 
  end subroutine lap_cell


  subroutine relaxac_mult_cell(mesh,g,fap,fap_aux,numit, w)
    !------------------------------------------------------------
    !RELAXAC_MULT_CELL
    ! Relaxation Scheme - Weighted Jacobi
    !------------------------------------------------------------
    type(grid_structure)::mesh             !Mesh
    type(scalar_field)::fap_aux     !Auxiliary vector
    !Approximate value of the function on each grid point
    type(scalar_field)::fap
    !Independent term of equation
    type(scalar_field)::g

    !parameter of relaxation
    real(r8)::w

    !Conter
    integer::k
    integer::i
    integer::j

    !Number of Iterations
    integer::numit

    !number of grid points
    integer::n

    real (r8):: ledhx     !Edge length
    real (r8)::Lchx       !Distance between two neighboring
    real (r8)::laptmp
    real (r8)::laptmp1

    !number points of the mesh
    n=mesh%nv

    !For each iteration
    do k=1,numit

      !For each grid point
      do i=1,n

        laptmp=0_r8
        laptmp1=0_r8
        ledhx=0_r8
        Lchx=0_r8

        !For each grid neighboring of i
        do j=1,mesh%v(i)%nnb

          !Edge length
          ledhx=mesh%edhx(mesh%v(i)%ed(j))%leng

          !Distance between two neighboring
          !arclen(mesh%v(mesh%v(i)%nb(j))%p,mesh%v(i)%p)
          Lchx=mesh%ed(mesh%v(i)%ed(j))%leng

          laptmp1=laptmp1+(ledhx/Lchx)
          laptmp=laptmp+(ledhx/Lchx)*(fap%f(mesh%v(i)%nb(j)))

        enddo

        laptmp=laptmp/mesh%hx(i)%areag
        laptmp1=laptmp1/mesh%hx(i)%areag

        !problem -lap f+f =g
        fap_aux%f(i)=(1-w)*(fap%f(i))+w*((laptmp+g%f(i))/(laptmp1+1))


         !problem lap f =g
         !fap_aux%f(i)=(1-w)*fap%f(i)+w*((laptmp-g%f(i))/(laptmp1))

      enddo

      fap%f=fap_aux%f

    enddo

    return 
  end subroutine relaxac_mult_cell

  subroutine relaxac_mult_cell_GS(mesh,g,fap,numit, w)
    !------------------------------------------------------------
    !RELAXAC_MULT_CELL
    ! Relaxation Scheme - Weighted Gauss-Seidel
    !------------------------------------------------------------
    type(grid_structure)::mesh  !Mesh
    type(scalar_field)::fap     !Approximate value of the function on each grid point
    type(scalar_field)::g       !Independent term of equation

    !parameter of relaxation
    real(r8)::w

    !Conter
    integer::k
    integer::i
    integer::j

    !Number of Iterations
    integer::numit

    !number of grid points
    integer::n

    real (r8):: ledhx     !Edge length
    real (r8)::Lchx       !Distance between two neighboring
    real (r8)::laptmp
    real (r8)::laptmp1

    !number points of the mesh
    n=mesh%nv

    !For each iteration
    do k=1,numit

      !For each grid point
      do i=1,n

        laptmp=0_r8
        laptmp1=0_r8
        ledhx=0_r8
        Lchx=0_r8

        !For each grid neighboring of i
        do j=1,mesh%v(i)%nnb

          !Edge length
          ledhx=mesh%edhx(mesh%v(i)%ed(j))%leng

          !Distance between two neighboring
          Lchx=mesh%ed(mesh%v(i)%ed(j))%leng
          laptmp1=laptmp1+(ledhx/Lchx)
          laptmp=laptmp+(ledhx/Lchx)*(fap%f(mesh%v(i)%nb(j)))
        enddo

        laptmp=laptmp/mesh%hx(i)%areag
        laptmp1=laptmp1/mesh%hx(i)%areag

        !para o problema lap f =g
        !fap%f(i)=(1-w)*fap%f(i)+w*((laptmp-g%f(i))/(laptmp1))


        !para o problema -lap f+f =g
        fap%f(i)=(1-w)*(fap%f(i))+w*((laptmp+g%f(i))/(laptmp1+1))


      enddo

    enddo

    return 
  end subroutine relaxac_mult_cell_GS


  function error_nrm_max(f, n)
    !-------------------------------------------
    !Calcula o maximo valor absoluto
    !  do vetor f de tamanho n
    !-------------------------------------------
    integer (i4):: n
    real (r8), dimension(1:n), intent(in) :: f
    real (r8):: error_nrm_max

    error_nrm_max=maxval(abs(f(1:n)))

    return
  end function error_nrm_max


  function error_norm_g(f,n)
    !-------------------------------------------
    !Calcula a raiz quadrada da soma
    !  dos elementos do vetor f
    !-------------------------------------------
    integer (i4):: n
    real (r8), dimension(1:n), intent(in) :: f
    real (r8):: error_norm_g
    real (r8):: sum_sq
    integer::i

    sum_sq=0.0

    do i=1,n 
      sum_sq=sum_sq+(f(i)**2)
    enddo

    error_norm_g=dsqrt(sum_sq)

    return
  end function error_norm_g


end module poisson

