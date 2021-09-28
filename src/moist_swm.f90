module moist_swm
  !=============================================================================
  !  Moist Shallow Water Model
  ! 
  !  Pack for simulations on the moist shallow water model developed
  !  by Zerroukat and Allen (2015) on the sphere using Voronoi grids
  ! 
  !  Luan da Fonseca Santos(luansantos@ime.usp.br)
  !  Oct 2019
  ! 
  !=============================================================================
    !Use global constants and kinds
    use constants, only: &
    i4, &
    r8, &
    erad, &
    rearth, &
    eps, &
    Omega, &
    datadir, &
    refdir, &
    rad2deg, &
    deg2rad, &
    pi, &
    pi2, &
    pio2, &
    piby2, &
    grav, &
    gravity, &
    gravi, &
    gravo2, &
    sec2day, &
    day2sec, &
    rotatn, &
    nlat_alt, &
    nlon_alt

    !Use main grid data structures
    use datastruct, only: &
    grid_structure, &
    scalar_field

    !Use routines from the spherical mesh pack
    use smeshpack, only: &
    convert_vec_sph2cart,&
    arcdistll, &
    sph2cart, &
    cart2sph, &
    convert_vec_cart2sph

    !Use routines from the interpolation pack
    use interpack, only: &
    plot_scalarfield, &
    plot_cart_vectorfield

    !Use routines from the SWM
    use swm_data      !Everything
    use swm_operators !Everything
    use swm, only: &
    initialize_gridprop, &
    sumf_areas, &
    calc_energies, &
    tendency, &
    error_calc

    !use refinter, only: &
    !andean_mountain_data, &
    !smooth_andean_mountain, &
    !earth_elevation, &
    !interpol_densf

    !=============================================================================
    !Moist shallow water model variables
    !=============================================================================

    !Temperature state variables (defined on voronoi centers)
    type(scalar_field):: theta       !temperature   - diagnostic variable
    type(scalar_field):: htheta      !h*temperature - prognostic variable
    type(scalar_field):: htheta_old  !h*temperature - prognostic variable 
    type(scalar_field):: h2theta     !h^2*temperature
    type(scalar_field):: hStheta     !Source for temperature equation
    type(scalar_field):: theta_exact !Only for test cases 2 and 3
    type(scalar_field):: theta_error !Only for test cases 2 and 3

    !Water (defined on voronoi centers)
    type(scalar_field):: water
    type(scalar_field):: hwater

    !Vapour state variables (defined on voronoi centers)
    type(scalar_field):: Qv          !vapour   - diagnostic variable
    type(scalar_field):: hQv         !h*vapour - prognostic variable
    type(scalar_field):: hQv_old     !h*vapour - prognostic variable
    type(scalar_field):: hSQv        !Source for vapour equation
    type(scalar_field):: delta_Qv    !Used in source computation 
    type(scalar_field):: Qv_exact    !Only for test cases 2 and 3
    type(scalar_field):: Qv_error    !Only for test cases 2 and 3

    !Cloud state variables  (defined on voronoi centers)
    type(scalar_field):: Qc          !cloud   - diagnostic variable
    type(scalar_field):: hQc         !h*cloud - prognostic variable
    type(scalar_field):: hQc_old     !h*cloud - prognostic variable
    type(scalar_field):: hSQc        !Source for cloud equation
    type(scalar_field):: delta_Qc    !used in source computation 
    type(scalar_field):: Qc_exact    !Only for test cases 2 and 3
    type(scalar_field):: Qc_error    !Only for test cases 2 and 3

    !Rain state variables  (defined on voronoi centers)
    type(scalar_field):: Qr          !rain   - diagnostic variable
    type(scalar_field):: hQr         !h*rain - prognostic variable
    type(scalar_field):: hQr_old     !h*rain - prognostic variable
    type(scalar_field):: hSQr        !Source for rain equation
    type(scalar_field):: delta_Qr    !used in source computation 
    type(scalar_field):: Qr_exact    !Only for test cases 2 and 3
    type(scalar_field):: Qr_error    !Only for test cases 2 and 3

    !Velocity source (defined on edges)
    type(scalar_field):: Su          !Source for momentum equation
    type(scalar_field):: ueast
    type(scalar_field):: unorth

    !Scalar fields from hx to ed (defined on edges)
    type(scalar_field):: theta_ed    !temperature
    type(scalar_field):: htheta_ed   !h*temperature
    type(scalar_field):: hQv_ed      !h*vapor
    type(scalar_field):: hQc_ed      !h*cloud
    type(scalar_field):: hQr_ed      !h*rain

    !velocity x scalar field (defined on edges - only normal component)
    type(scalar_field):: uhtheta     !velocity*h*temperature
    type(scalar_field):: uhQv        !velocity*h*vapour
    type(scalar_field):: uhQc        !velocity*h*cloud
    type(scalar_field):: uhQr        !velocity*h*rain

    !Divergences (defined on voronoi centers)
    type(scalar_field):: div_uhtheta !div of velocity*h*temperature
    type(scalar_field):: div_uhQv    !div of velocity*h*vapour
    type(scalar_field):: div_uhQc    !div of velocity*h*cloud
    type(scalar_field):: div_uhQr    !div of velocity*h*rain

    !Exact divergences (defined on voronoi centers) 
    !Only for test cases 2 and 3
    type(scalar_field):: div_uhtheta_exact !div of velocity*h*temperature
    type(scalar_field):: div_uhQv_exact    !div of velocity*h*temperature
    type(scalar_field):: div_uhQc_exact    !div of velocity*h*cloud
    type(scalar_field):: div_uhQr_exact    !div of velocity*h*rain

    !Gradients (defined on edges)
    type(scalar_field):: gradPI         !gradient of h^2*temperature
    type(scalar_field):: gradPI_oh      !gradient of h^2*temperature over h
    type(scalar_field):: grad_b         !gradient of topography
    type(scalar_field):: theta_grad_b   !temperature*gradient of topography

    !Exact gradients (defined on edges)
    !Only for test cases 2 and 3
    type(scalar_field):: gradPI_exact       !gradient of h^2*temperature
    type(scalar_field):: gradPI_oh_exact    !gradient of h^2*temperature over h
    type(scalar_field):: grad_b_exact       !gradient of topography
    type(scalar_field):: theta_grad_b_exact !temperature*gradient of topography

    !RHS
    type(scalar_field):: tempeq     !Temperature
    type(scalar_field):: vapoureq   !Vapour
    type(scalar_field):: cloudeq    !Cloud
    type(scalar_field):: raineq     !Rain

    !Runge-Kutta variables
    real(r8), dimension(:), allocatable:: tempf0
    real(r8), dimension(:), allocatable:: tempf1
    real(r8), dimension(:), allocatable:: tempf2
    real(r8), dimension(:), allocatable:: tempf3

    real(r8), dimension(:), allocatable:: vapourf0
    real(r8), dimension(:), allocatable:: vapourf1
    real(r8), dimension(:), allocatable:: vapourf2
    real(r8), dimension(:), allocatable:: vapourf3

    real(r8), dimension(:), allocatable:: cloudf0
    real(r8), dimension(:), allocatable:: cloudf1
    real(r8), dimension(:), allocatable:: cloudf2
    real(r8), dimension(:), allocatable:: cloudf3

    real(r8), dimension(:), allocatable:: rainf0
    real(r8), dimension(:), allocatable:: rainf1
    real(r8), dimension(:), allocatable:: rainf2
    real(r8), dimension(:), allocatable:: rainf3

    !Parameters of Zerroukat and Allen 2015 JCP paper
    real(r8)::q0
    real(r8)::Lscale
    real(r8)::gamma_r
    real(r8)::gamma_v
    real(r8)::q_precip
    real(r8)::Twater, iniwater


!=============================================================================
contains

subroutine allocate_global_moistswm_vars()
    !------------------------------------
    !Allocate fields for the moist SWM
    !------------------------------------

    !Error flag
    integer(i4):: ist

    !Allocate the SWM variables
    call allocate_globalswmvars()

    !Allocate the additional variables of the moist SWM
    tempeq%n=mesh%nv
    tempeq%name="tempeq"
    tempeq%pos=0
    allocate(tempeq%f(1:tempeq%n), stat=ist)

    vapoureq%n=mesh%nv
    vapoureq%name="vapoureq"
    vapoureq%pos=0
    allocate(vapoureq%f(1:vapoureq%n), stat=ist)

    cloudeq%n=mesh%nv
    cloudeq%name="cloudeq"
    cloudeq%pos=0
    allocate(cloudeq%f(1:cloudeq%n), stat=ist)

    raineq%n=mesh%nv
    raineq%name="raineq"
    raineq%pos=0
    allocate(raineq%f(1:raineq%n), stat=ist)

    !Water
    water%n=mesh%nv
    water%name="water"
    water%pos=0
    allocate(water%f(1:water%n), stat=ist)
    allocate(hwater%f(1:water%n), stat=ist)

    !Temperature
    theta%n=mesh%nv
    theta%name="temperature"
    theta%pos=0
    allocate(theta%f(1:theta%n), stat=ist)
    allocate(htheta_old%f(1:theta%n), stat=ist)
    allocate(htheta%f(1:theta%n), stat=ist)
    allocate(h2theta%f(1:theta%n), stat=ist)
    allocate(hStheta%f(1:theta%n), stat=ist)
    allocate(theta_exact%f(1:theta%n), stat=ist)
    allocate(theta_error%f(1:theta%n), stat=ist)

    !Vapour
    Qv%n=mesh%nv
    Qv%name="vapour"
    Qv%pos=0
    allocate(Qv%f(1:Qv%n), stat=ist)
    allocate(hQv_old%f(1:Qv%n), stat=ist)
    allocate(hQv%f(1:Qv%n), stat=ist)
    allocate(hSQv%f(1:Qv%n), stat=ist)
    allocate(delta_Qv%f(1:Qv%n), stat=ist)
    allocate(Qv_exact%f(1:Qv%n), stat=ist)
    allocate(Qv_error%f(1:Qv%n), stat=ist)

    !Cloud
    Qc%n=mesh%nv
    Qc%name="cloud"
    Qc%pos=0
    allocate(Qc%f(1:Qc%n), stat=ist)
    allocate(hQc%f(1:Qc%n), stat=ist)
    allocate(hQc_old%f(1:Qc%n), stat=ist)
    allocate(hSQc%f(1:Qc%n), stat=ist)
    allocate(delta_Qc%f(1:Qc%n), stat=ist)
    allocate(Qc_exact%f(1:Qc%n), stat=ist)
    allocate(Qc_error%f(1:Qc%n), stat=ist)

    !Rain
    Qr%n=mesh%nv
    Qr%name="rain"
    Qr%pos=0
    allocate(Qr%f(1:Qr%n), stat=ist)
    allocate(hQr%f(1:Qr%n), stat=ist)
    allocate(hQr_old%f(1:Qr%n), stat=ist)
    allocate(hSQr%f(1:Qr%n), stat=ist)
    allocate(delta_Qr%f(1:Qr%n), stat=ist)
    allocate(Qr_exact%f(1:Qr%n), stat=ist)
    allocate(Qr_error%f(1:Qr%n), stat=ist)

    !Divergence terms
    div_uhQv%n=mesh%nv
    div_uhQv%pos=0
    div_uhQv%name="divuqv"
    allocate(div_uhQv%f(1:div_uhQv%n), stat=ist)
    allocate(div_uhQc%f(1:div_uhQv%n), stat=ist)
    allocate(div_uhQr%f(1:div_uhQv%n), stat=ist)
    allocate(div_uhtheta%f(1:div_uhQv%n), stat=ist)

    allocate(div_uhQv_exact%f(1:div_uhQv%n), stat=ist)
    allocate(div_uhQc_exact%f(1:div_uhQv%n), stat=ist)
    allocate(div_uhQr_exact%f(1:div_uhQv%n), stat=ist)
    allocate(div_uhtheta_exact%f(1:div_uhQv%n), stat=ist)

    !Scalar fields on edges
    hQv_ed%n=mesh%ne
    hQv_ed%name="qv_ed"
    hQv_ed%pos=edpos
    allocate(hQv_ed%f(1:hQv_ed%n), stat=ist)
    allocate(hQc_ed%f(1:hQv_ed%n), stat=ist)
    allocate(hQr_ed%f(1:hQv_ed%n), stat=ist)
    allocate(htheta_ed%f(1:hQv_ed%n), stat=ist)
    allocate(theta_ed%f(1:hQv_ed%n), stat=ist)

    !Gradient terms
    gradPI%n=mesh%ne
    gradPI%name="gradPI"
    gradPI%pos=edpos
    allocate(gradPI%f(1:hQv_ed%n), stat=ist)
    allocate(gradPI_oh%f(1:hQv_ed%n), stat=ist)
    allocate(grad_b%f(1:hQv_ed%n), stat=ist)
    allocate(theta_grad_b%f(1:hQv_ed%n), stat=ist)

    allocate(gradPI_exact%f(1:hQv_ed%n), stat=ist)
    allocate(gradPI_oh_exact%f(1:hQv_ed%n), stat=ist)
    allocate(grad_b_exact%f(1:hQv_ed%n), stat=ist)
    allocate(theta_grad_b_exact%f(1:hQv_ed%n), stat=ist)

    !velocity x scalar field (defined on edges - only normal component)
    uhQv%n=mesh%ne
    uhQv%name="uQv"
    uhQv%pos=edpos
    allocate(uhQv%f(1:uhQv%n), stat=ist)
    allocate(uhQc%f(1:uhQv%n), stat=ist)
    allocate(uhQr%f(1:uhQv%n), stat=ist)
    allocate(uhtheta%f(1:uhQv%n), stat=ist)

    !Runge-Kutta variables 
    !Temperature
    allocate(tempf0(1:theta%n))
    allocate(tempf1(1:theta%n))
    allocate(tempf2(1:theta%n))
    allocate(tempf3(1:theta%n))
    
    !Vapour
    allocate(vapourf0(1:Qv%n))
    allocate(vapourf1(1:Qv%n))
    allocate(vapourf2(1:Qv%n))
    allocate(vapourf3(1:Qv%n))

    !Cloud
    allocate(cloudf0(1:Qc%n))
    allocate(cloudf1(1:Qc%n))
    allocate(cloudf2(1:Qc%n))
    allocate(cloudf3(1:Qc%n))

    !Rain
    allocate(rainf0(1:Qr%n))
    allocate(rainf1(1:Qr%n))
    allocate(rainf2(1:Qr%n))
    allocate(rainf3(1:Qr%n))

    !Velocities (defined on voronoi centers)
    ueast%n=mesh%nv
    ueast%name="ueast"
    ueast%pos=0
    allocate(ueast%f(1:ueast%n), stat=ist)
    allocate(unorth%f(1:ueast%n), stat=ist)

    !Source for momentum equation
    Su%n=mesh%ne
    Su%name="su"
    Su%pos=edpos
    allocate(Su%f(1:u%n), stat=ist)

    IF(ist>0) STOP 'Error in allocate_globalmoistswmvars'

    call initialize_global_moist_swm_vars()

end subroutine allocate_global_moistswm_vars


subroutine initialize_global_moist_swm_vars()
    !---------------------------------------------------
    !Initialize fields with zero in paralel
    !---------------------------------------------------

    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(theta, theta_exact, theta_error, htheta, htheta_old, h2theta) &
    !$OMP SHARED(water, hwater) &
    !$OMP SHARED(Qv, Qv_exact, Qv_error, hQv, hQv_old) &
    !$OMP SHARED(Qc, Qc_exact, Qc_error, hQc ,hQc_old) &
    !$OMP SHARED(Qr, Qr_exact, Qr_error, hQr, hQr_old) &
    !$OMP SHARED(hStheta, hSQv, hSQr, hSQc) &
    !$OMP SHARED(delta_Qv, delta_Qr, delta_Qc) &
    !$OMP SHARED(div_uhQv,div_uhQc,div_uhQr,div_uhtheta) &
    !$OMP SHARED(gradPI,gradPI_oh,grad_b,theta_grad_b) &
    !$OMP SHARED(div_uhQv_exact,div_uhQc_exact,div_uhQr_exact,div_uhtheta_exact) &
    !$OMP SHARED(gradPI_exact,gradPI_oh_exact,grad_b_exact,theta_grad_b_exact) &
    !$OMP SHARED(uhQv,uhQc,uhQr,uhtheta) &
    !$OMP SHARED(hqv_ed,hqc_ed,hqr_ed,htheta_ed,theta_ed) &
    !$OMP SHARED(tempeq, tempf0, tempf1, tempf2, tempf3) &
    !$OMP SHARED(vapoureq, vapourf0, vapourf1, vapourf2, vapourf3) &
    !$OMP SHARED(cloudeq, cloudf0, cloudf1, cloudf2, cloudf3) &
    !$OMP SHARED(raineq, rainf0, rainf1, rainf2, rainf3) &
    !$OMP SHARED(ueast, unorth,Su) 

    tempeq%f(1:theta%n)=0._r8
    vapoureq%f(1:Qv%n)=0._r8
    cloudeq%f(1:Qc%n)=0._r8
    raineq%f(1:Qr%n)=0._r8

    !Water
    water%f = 0._r8
    hwater = water

    !Temperature
    theta%f=0._r8
    theta_exact=theta
    theta_error=theta
    htheta=theta
    htheta_old=theta
    h2theta=theta
    hStheta=theta

    !Vapour
    Qv%f=0._r8
    Qv_exact=Qv
    Qv_error=Qv
    hQv=Qv
    hQv_old=Qv
    hSQv=Qv
    delta_Qv=Qv

    !Cloud
    Qc%f=0._r8
    Qc_exact=Qc
    Qc_error=Qc
    hQc=Qc
    hQc_old=Qc
    hSQc=Qc
    delta_Qc=Qc

    !Rain
    Qr%f=0._r8
    Qr_exact=Qr
    Qr_error=Qr
    hQr=Qr
    hQr_old=Qr
    hSQr=Qr
    delta_Qr=Qr

    !Runge kutta variables
    !Temperature
    tempf0(1:theta%n)=0._r8
    tempf1(1:theta%n)=0._r8
    tempf2(1:theta%n)=0._r8
    tempf3(1:theta%n)=0._r8

    !Vapour
    vapourf0(1:qv%n)=0._r8
    vapourf1(1:qv%n)=0._r8
    vapourf2(1:qv%n)=0._r8
    vapourf3(1:qv%n)=0._r8

    !Cloud
    cloudf0(1:qc%n)=0._r8
    cloudf1(1:qc%n)=0._r8
    cloudf2(1:qc%n)=0._r8
    cloudf3(1:qc%n)=0._r8

    !Rain
    rainf0(1:qr%n)=0._r8
    rainf1(1:qr%n)=0._r8
    rainf2(1:qr%n)=0._r8
    rainf3(1:qr%n)=0._r8

    !Divergences
    div_uhQv%f=0._r8
    div_uhQc=div_uhQv
    div_uhQr=div_uhQv
    div_uhtheta=div_uhQv

    div_uhQv_exact=div_uhQv
    div_uhQc_exact=div_uhQv
    div_uhQr_exact=div_uhQv
    div_uhtheta_exact=div_uhQv

    !Fields at edges
    hQv_ed%f=0._r8
    hQc_ed=hQv_ed
    hQr_ed=hQv_ed
    htheta_ed=hQv_ed
    theta_ed=hQv_ed

    !Gradients
    gradPI%f=0._r8
    gradPI_oh=gradPI
    grad_b=gradPI
    theta_grad_b=gradPI   
    
    gradPI_exact=gradPI
    gradPI_oh_exact=gradPI
    grad_b_exact=gradPI
    theta_grad_b_exact=gradPI   
    
    !velocity x scalar field
    uhQv%f=0._r8
    uhQc=uhQv
    uhQr=uhQv
    uhtheta=uhQv

    !Velocities
    ueast%f=0._r8
    unorth%f=0._r8
    Su%f=0._r8

    !$OMP END PARALLEL WORKSHARE
    
  end subroutine initialize_global_moist_swm_vars

  function qsat(T,H,B,q_0)
    real(r8)::qsat
    real(r8), intent(in)::T,H,B,q_0
    qsat = q_0*exp(20._r8*T)/((H+B)*grav)
  end function qsat

  function dqsat_dtheta(T,H,B,q_0)
    real(r8)::dqsat_dtheta
    real(r8), intent(in)::T,H,B,q_0
    dqsat_dtheta = 20._r8*qsat(T,H,B,q_0)
  end function dqsat_dtheta

  function F_quad(f1,f2,f3,lat)
    real(r8)::F_quad
    real(r8), intent(in)::f1,f2,f3,lat
    F_quad = lat*(lat-pio2)*f1 -2._r8*(lat+pio2)*(lat-pio2)*f2 + lat*(lat+pio2)*f3
    F_quad = F_quad*2._r8/(pi*pi)
  end function F_quad

  function bottom(lat,lon,latc,lonc,rlat,rlon)
    real(r8)::bottom
    real(r8), intent(in)::lat,lon,latc,lonc,rlat,rlon
    bottom = 1._r8 - min(1._r8,sqrt( ((lon-lonc)/rlon)**2 + ((lat-latc)/rlat)**2))
  end function bottom

  subroutine initialize_moist_swm_fields()
    !---------------------------------------------------
    ! initialize_fields
    ! see dataswm.f90 for field that will be initialized
    !---------------------------------------------------
  
    !Auxiliar variables
    real(r8)::lat
    real(r8)::lon
    real(r8)::lat0
    real(r8)::lon0
        real(r8)::lon1
                real(r8)::r1
    real(r8)::lat_rotated
    real(r8)::lon_rotated
    real(r8)::latd
    real(r8)::lond
    real(r8)::utmp
    real(r8)::vtmp
    real(r8)::rmax
    real(r8)::h0
    real(r8)::h_ct
    real(r8)::mu1
    real(r8)::mu2
    real(r8)::mu3
    real(r8)::mu4
    real(r8)::theta_sp
    real(r8)::theta_eq
    real(r8)::theta_np
    real(r8)::theta_min
    real(r8)::theta_max
    real(r8)::vectmp(1:3)
    real(r8)::postmp(1:3)
    real(r8)::u0
    real(r8)::w
    real(r8)::sigma
    real(r8)::phi0
    real(r8)::phi_ct
    real(r8)::temp0
    real(r8)::xsi
    real(r8)::Rmat(1:3,1:3) !Rotation matrix
        real(r8)::Rmat2(1:3,1:3) !Rotation matrix
    real(r8)::RmatT(1:3,1:3) !Transpose/inverse of Rotation matrix 
    real(r8)::h1
    real(r8)::sinn
    real(r8)::coss

    !Orography
    real(r8), allocatable :: alt_table(:,:)
    integer(i4):: iunit

    !Galewiski test case
    integer(i4):: j
    integer(i4):: jy
    integer(i4):: nygg
    real(r8):: u00
    real(r8):: alpha
    real(r8):: clat
    real(r8), allocatable :: hgg(:)
    real(r8):: l1
    real(r8):: l2
    real(r8):: lat1
    real(r8):: lat2
    real(r8):: en
    real(r8):: umen
    real(r8):: dygg
    real(r8):: beta
    real(r8):: totvol
    real(r8):: totarea
    real(r8):: den
    real(r8):: cc1
    real(r8):: cc2
    real(r8):: cc3
    real(r8):: cc4
    real(r8):: u1
    real(r8):: u2
    real(r8):: e1
    real(r8):: e2
    real(r8):: hpert
    real(r8):: long
      
    !Indexes
    integer(i4):: i !voronoi cell index
    integer(i4):: k !trinagle index
    integer(i4):: l !edge index

    Lscale=10._r8
    gamma_r=0.001_r8
    q_precip=0.0001_r8
    latd = 0._r8*deg2rad
    lond = 0._r8*deg2rad

    select case(testcase)
    !======================================================================================
    ! Steady state - from Zerroukat and Allen JCP 2015
    !======================================================================================
    case(2)
      !Field at cell's center
      u0    = 20._r8
      phi0  = 3._r8*10**4 
      w     = omega*erad*u0 + u0*u0*0.5_r8
      sigma = w/10._r8
      temp0 = phi0*phi0/300._r8
      xsi = 0.00_r8

      do i=1, mesh%nv
        lon = mesh%v(i)%lon
        lat = mesh%v(i)%lat
        sinn = sin(lat)
        coss = cos(lat)

        h%f(i) = (phi0 -(w+sigma)*sinn**2)*gravi
        bt%f(i) = 0._r8 

        theta%f(i) = temp0 + sigma*(coss**2)*((w+sigma)*(coss**2) + (phi0-w-sigma)*2._r8 )
        theta%f(i) = theta%f(i)/( phi0**2 + ((w+sigma)**2)*sinn**4 -2._r8*phi0*(w+sigma)*sinn**2)
        Qv%f(i) = (1._r8-xsi)*qsat(theta%f(i),h%f(i),bt%f(i),1._r8)

        !Fluxes
        hQv%f(i) = h%f(i)*Qv%f(i)
        hQc%f(i) = h%f(i)*Qc%f(i)
        hQr%f(i) = h%f(i)*Qr%f(i)
        htheta%f(i) = h%f(i)*theta%f(i)
      end do

      !print*,minval(theta%f),maxval(theta%f)
      !stop
      q0 = 0.02_r8/maxval(qv%f)
      !print*,maxval(qv%f)
      Qv%f = q0*Qv%f
      hQv%f = q0*hQv%f

      !Velocity
      if(useStagHTC)then
        do l=1, mesh%ne
          lat = mesh%ed(l)%c%lat
          lon = mesh%ed(l)%c%lon
          utmp = u0*cos(lat)
          vtmp = 0._r8
          call convert_vec_sph2cart(utmp, vtmp, mesh%ed(l)%c%p, vectmp)
          v_ed%p(l)%v=vectmp
          u%f(l)=dot_product(vectmp,mesh%ed(l)%tg)
        end do

      elseif(useStagHC)then
        do l=1, mesh%ne
          !vectmp=mesh%edhx(l)%c%p
          !vectmp = matmul(RmatT,vectmp)
          !call cart2sph(vectmp(1),vectmp(2),vectmp(3),lon,lat)
          !utmp = u0*dcos(lat)
          !vtmp = 0._r8
          !call convert_vec_sph2cart(utmp, vtmp, mesh%edhx(l)%c%p, vectmp)
          !v_ed%p(l)%v=vectmp
          !u%f(l)=dot_product(vectmp,mesh%edhx(l)%nr)
       end do
      end if

      h_exact = h
      u_exact = u
      theta_exact = theta
      Qv_exact = Qv
      Qc_exact = Qc
      Qr_exact = Qr

    !======================================================================================
    ! Flow over a mountain - from Zerroukat and Allen JCP 2015
    !======================================================================================
    case(4)

    !Parameters
    u0     = 20._r8
    phi0   = 5960._r8*grav
    h0     = 2000._r8
    phi_ct = (erad*omega*u0+(u0**2)/2._r8)
    w      = omega*erad*u0 + u0*u0*0.5_r8
    sigma  = 0._r8!w/10._r8
    xsi    = 0.001_r8
    mu1    = 0.05_r8
    mu2    = 0.98_r8
    theta_sp = -40._r8/300._r8 
    theta_eq =  30._r8/300._r8 
    theta_np = -20._r8/300._r8 
     
    !lon0  = -pi*0.5_r8 -0.52777d0*pi
    !lon1 = lon0 + 2.d0*pi
    !lat0  = -pi/6._r8  
    !lat1  = -pi/6._r8   
    rmax   = pi/9._r8
    lon0  = pi*0.5_r8 + 30.d0*deg2rad
    lat0  = pi/6._r8
          
    !Variables at Voronoi centers
    do i=1, mesh%nv
        vectmp = mesh%v(i)%p
        lon = mesh%v(i)%lon
        lat = mesh%v(i)%lat
        sinn = sin(lat)
        coss = cos(lat)
   
        h%f(i) = (phi0 -w*sinn**2)*gravi

        r  = dsqrt((lon-lon0)**2+(lat-lat0)**2)

        if(r<rmax)then
          bt%f(i)=2000._r8*(1._r8-r/rmax)
        else
          bt%f(i)=0.
        endif

        ! Correct h to allow orography
        h%f(i)=h%f(i)-bt%f(i)
    
        !Temperature 
        !lon = lon+1.d0*pi + 0.52777d0*pi !-95graus
        lon = lon+pi+pi + 30.d0*deg2rad !+ 0.52777d0*pi !-95graus
        theta%f(i)=F_quad(theta_sp,(1._r8-mu1)*theta_eq,theta_np,lat) + mu1*theta_eq*dcos(lat)*dsin(lon)
      
        !Vapour
        Qv%f(i) = mu2*qsat(theta%f(i),h%f(i),bt%f(i),1._r8)
        
        !Fluxes
        hQv%f(i) = h%f(i)*Qv%f(i)
        hQc%f(i) = h%f(i)*Qc%f(i)
        hQr%f(i) = h%f(i)*Qr%f(i)
        htheta%f(i) = h%f(i)*theta%f(i)
      end do

      q0 = 0.02_r8/maxval(qv%f)
      Qv%f = q0*Qv%f
      hQv%f = q0*hQv%f

      !print*,maxval(bt%f),minval(bt%f)
      !stop
      alpha = 0.d0
      if(useStagHTC)then
        do l=1, mesh%ne          
          lat = mesh%ed(l)%c%lat
          lon = mesh%ed(l)%c%lon
          utmp=u0*(cos(lat)*cos(alpha) + cos(lon)*sin(lat)*sin(alpha))
          vtmp=-u0*(sin(lon)*sin(alpha))
          call convert_vec_sph2cart(utmp, vtmp, mesh%ed(l)%c%p, vectmp)
          v_ed%p(l)%v=vectmp
          u%f(l)=dot_product(vectmp,mesh%ed(l)%tg)
        end do

      elseif(useStagHC)then
        do l=1, mesh%ne
          utmp=u0*dcos(mesh%edhx(l)%c%lat)
          vtmp=0._r8
          call convert_vec_sph2cart(utmp, vtmp, mesh%edhx(l)%c%p, vectmp)
          v_ed%p(l)%v=vectmp
          u%f(l)=dot_product(vectmp,mesh%edhx(l)%nr)
        end do
      end if

    !======================================================================================
    ! Galewski - Jet in Southern Hemisphere
    !======================================================================================
    case(7,8)
      !Parameters
      !mu1    = 0.05_r8
      mu1    = 0.00002_r8
      mu2    = 0.98_r8
      theta_sp = -40._r8/300._r8 
      theta_eq =  30._r8/300._r8 
      theta_np = -20._r8/300._r8 
      if(testcase==7)then
        u00 = 80.0
        lat0 = pi/7.0
        lat1 = pi/2.0 - lat0
      else
        if(testcase==8)then ! Jet in Southern Hemisphere
           u00 = 80.0
           lat0 = -5.d0*deg2rad
           lat1 = -45.d0*deg2rad
        end if
      end if
	
      en = exp(-4/(lat1 - lat0)**2)
      umen = u00/en
      totvol = 0.0D0
      totarea = 0.0D0      

      ! Integrate to tabulate h and psi as functions of geographical
      ! latitude
      nygg = 4*FLOOR(SQRT(REAL(mesh%nv)))
      allocate(hgg(nygg+1)) !, psigg(nygg+1))
      dygg = pi/nygg
      hgg(1) = 0.0D0
      !psigg(1) = 0.0D0
      do j = 2, nygg
        l1 = (j-2)*dygg - piby2
        den = (l1 - lat0)*(l1 - lat1)
        if (den .lt. 0.0D0) then
          u1 = umen*exp(1.0D0/den)
        else
          u1 = 0.0D0
        endif
        l2 = (j-1)*dygg - piby2
        den = (l2 - lat0)*(l2 - lat1)
        if (den .lt. 0.0D0) then
          u2 = umen*exp(1.0D0/den)
        else
          u2 = 0.0D0
        endif
        !print*, l2*rad2deg, " u: ", (u1+u2)*0.5_r8
        !psigg(j) = psigg(j-1) - 0.5d0*(u1 + u2)*dygg
        u1 = u1*(2.0d0*rotatn*SIN(l1) + TAN(l1)*u1/rearth)
        u2 = u2*(2.0d0*rotatn*SIN(l2) + TAN(l2)*u2/rearth)
        hgg(j) = hgg(j-1) - rearth*0.5d0*(u1 + u2)*dygg

        totarea = totarea + COS(l2)*dygg
        totvol = totvol + hgg(j)*COS(l2)*dygg

      enddo
      !psigg(nygg+1) = psigg(nygg)
      hgg(nygg+1) = hgg(nygg)
      totvol = totvol/(totarea*gravity)
      hgg = hgg + (1.0D4 - totvol)*gravity !potential phi2
      hgg=hgg*gravi !Height

      ! Now assign h as a function of geographical latitude
      ! using interpolation from tabulated values
      totvol = 0.00
      totarea = 0.0D0
      do i = 1, mesh%nv
        ! l1 = flat(if0,ngrids) + piby2
        !CALL centroid(if0,long,lat)
        lat=mesh%v(i)%lat !modif psp
        l1 = lat + piby2
        jy = floor(l1/dygg) + 1
        beta = (l1 - (jy - 1)*dygg)/dygg
        if (jy == 1 .or. jy == nygg) then
          ! Linear interpolation
          cc2 = 1.0D0 - beta
          cc3 = beta
          !phi2(if0) = (cc2*hgg(jy) + cc3*hgg(jy+1))*farea(if0,ngrids)
          h%f(i)=(cc2*hgg(jy) + cc3*hgg(jy+1)) !modif psp
        else
          ! Cubic interpolation
          cc1 = -beta*(beta - 1.0D0)*(beta - 2.0D0)/6.0D0
          cc2 = 0.5D0*(beta + 1.0D0)*(beta - 1.0D0)*(beta - 2.0D0)
          cc3 = -0.5D0*(beta + 1.0D0)*beta*(beta - 2.0D0)
          cc4 = (beta + 1.0D0)*beta*(beta - 1.0D0)/6.0D0
          !phi2(if0) = (cc1*hgg(jy-1) + cc2*hgg(jy) + cc3*hgg(jy+1) + cc4*hgg(jy+2))*farea(if0,ngrids)
          h%f(i) = (cc1*hgg(jy-1) + cc2*hgg(jy) + cc3*hgg(jy+1) + cc4*hgg(jy+2)) !*farea(if0,ngrids)
        endif
        totarea = totarea + mesh%hx(i)%areag !farea(if0,ngrids)

        totvol = totvol + h%f(i) ! phi2(if0)
         !print*, totarea, totvol
      enddo
      deallocate(hgg)!, psigg)


      !Set velocity field
      do l=1,mesh%ne
        utmp=0._r8
        vtmp=0._r8
        if(useStagHC)then
          lat = mesh%edhx(l)%c%lat
          den = (lat - lat0)*(lat - lat1)
          if (den .lt. 0.0D0) then
            utmp = umen*exp(1.0D0/den)
          end if
          call convert_vec_sph2cart(utmp, vtmp, mesh%edhx(l)%c%p, vectmp)
          v_ed%p(l)%v=vectmp
          u%f(l)=dot_product(vectmp,mesh%edhx(l)%nr)
        elseif(useStagHTC)then
          lat = mesh%ed(l)%c%lat
          den = (lat - lat0)*(lat - lat1)
          if (den .lt. 0.0D0) then
            utmp = umen*exp(1.0D0/den)
          end if
          call convert_vec_sph2cart(utmp, vtmp, mesh%ed(l)%c%p, vectmp)
          v_ed%p(l)%v=vectmp
          u%f(l)=dot_product(vectmp,mesh%ed(l)%tg)
        end if
      end do

      !Add perturbation 
      ! Geopotential perturbation
      hpert = 120.0D0
      alpha = 1.0D0/3.0D0
      beta = 1.0D0/15.0D0

      if(testcase == 8)then
        lat2 = -25.d0*deg2rad
      else
        lat2 = 0.5D0*piby2
      end if

      do i = 1, mesh%nv
        lat=mesh%v(i)%lat
        long=mesh%v(i)%lon
        l2 = lat
        !l1 = long - 0.7d0*pi
        l1 = long - 0.93d0*pi

        clat = COS(l2)
        if(testcase ==8)then
           l1 = l1 + 120.d0*deg2rad
        end if
        e1 = EXP(-(l1/alpha)**2)
        e2 = EXP(-((lat2 - l2)/beta)**2)
        h%f(i) = h%f(i)+hpert*clat*e1*e2
        !h%f(i) = 1.d0+hpert*clat*e1*e2
        !print*, h%f(i) 
      enddo
      !stop

      !Variables at Voronoi centers
      do i=1, mesh%nv
        lon    = mesh%v(i)%lon
        lat    = mesh%v(i)%lat
        !Temperature
        !lon = lon+pi
        !theta%f(i)= theta_sp*lat*(lat-pio2) - (1._r8-mu1)*theta_eq*(lat+pio2)*(lat-pio2) + theta_np*lat*(lat+pio2)
        theta%f(i)= theta_sp*lat*(lat+pio2) - (1._r8-mu1)*theta_eq*(lat+pio2)*(lat-pio2) + theta_np*lat*(lat-pio2)

        !theta%f(i)= theta_sp*lat*(lat-pio2) - (1._r8-mu1*theta_eq)*(lat+pio2)*(lat-pio2) + theta_np*lat*(lat+pio2)
        !theta%f(i)=F_quad(theta_sp,(1._r8-mu1)*theta_eq,theta_np,lat)! + mu1*theta_eq*dcos(lat)*dsin(lon)
        !theta%f(i)=F_quad(theta_sp,theta_eq,theta_np,lat)! + mu1*theta_eq*dcos(lat)*dsin(lon)

        !Vapour
        Qv%f(i) = mu2*qsat(theta%f(i),h%f(i),bt%f(i),1._r8)
        !Qv%f(i) = qsat(theta%f(i),h%f(i),bt%f(i),1._r8)

        !Fluxes
        hQv%f(i) = h%f(i)*Qv%f(i)
        hQc%f(i) = h%f(i)*Qc%f(i)
        hQr%f(i) = h%f(i)*Qr%f(i)
        htheta%f(i) = h%f(i)*theta%f(i)
      end do

      q0 = 0.02_r8/maxval(qv%f)
      Qv%f = q0*Qv%f
      hQv%f = q0*hQv%f

      !print*,q0,maxval(qv%f),minval(qv%f)
      !stop 

    case default
      print*, "SWM_initialize_fields error - please select a proper test case:", testcase
      stop
    end select

   !Check for CFL constraints
    maxvel=maxval(u%f(1:mesh%ne))
    cfl=abs(maxvel*dt/(mesh%minvdist*erad))

    print*, "CFL:", cfl
      if(cfl>2)then
        print*, "CFL too large, problems may occur"
    end if
  end subroutine initialize_moist_swm_fields


  subroutine tendency_moist_swm(h, u, htheta, hQv, hQc, hQr, masseq, momeq, tempeq, vapoureq, cloudeq, raineq)
    !--------------------------------------
    !Calculates the Right Hand Side (spatial discret./tendency)
    !   of mass, velocity and moist variables equations
    !-------------------------------------------------------------

    !Fluid thickness (defined on voronoi centers)
    type(scalar_field), intent(in):: h  !General

    !Velocities (defined on edges - only normal component)
    type(scalar_field), intent(in):: u  !General

    !Temperature (defined on voronoi centers)
    type(scalar_field), intent(in):: htheta  !General

    !Vapour (defined on voronoi centers)
    type(scalar_field), intent(in):: hQv  !General

    !Cloud (defined on voronoi centers)
    type(scalar_field), intent(in):: hQc  !General

    !Rain (defined on voronoi centers)
    type(scalar_field), intent(in):: hQr  !General

    !Time
    !real(r8), intent(in):: dtime

    !Right hand side of mass equation (number of cell equations)
    real(r8), intent(inout)::masseq(:)

    !Right hand side of momentum equation (number of edge equations)
    real(r8), intent(inout)::momeq(:)

    !Right hand side of temperature (number of cell equations)
    real(r8), intent(inout)::tempeq(:)

    !Right hand side of vapour (number of cell equations)
    real(r8), intent(inout)::vapoureq(:)
    
    !Right hand side of cloud (number of cell equations)
    real(r8), intent(inout)::cloudeq(:)
    
    !Right hand side of rain (number of cell equations)
    real(r8), intent(inout)::raineq(:)  


    !Compute the SWM tendency
    call tendency(h, u, masseq, momeq)

    !Initialize RHS of temperature and moist variables equations (in paralel)
    call zero_vector(tempeq)
    call zero_vector(vapoureq)
    call zero_vector(cloudeq)
    call zero_vector(raineq)

    !===============================================================
    !Calculate temperature tendency
    !===============================================================

    !Interpolate temperature to edges and calculate flux at edges
    call scalar_hx2ed(htheta, htheta_ed, mesh)      !htheta: cell->edge
    call scalar_elem_product(u, htheta_ed, uhtheta) !Flux uhtheta at edges

    !Calculate divergence / temp eq RHS
    call div_hx(uhtheta, div_uhtheta, mesh)  
    !Temp. eq. RHS
    tempeq = -div_uhtheta%f
    
    !===============================================================
    !Calculate vapour tendency
    !===============================================================

    !Interpolate vapour to edges and calculate flux at edges
    call scalar_hx2ed(hQv, hQv_ed, mesh)      !hQv: cell->edge
    call scalar_elem_product(u, hQv_ed, uhQv) !Flux uhQv at edges

    !Calculate divergence / vapour eq RHS
    call div_hx(uhQv, div_uhQv, mesh)

    !Vapour eq. RHS
    vapoureq = -div_uhQv%f

    !===============================================================
    !Calculate cloud tendency
    !===============================================================

    !Interpolate vapour to edges and calculate flux at edges
    call scalar_hx2ed(hQc, hQc_ed, mesh)      !hQc: cell->edge
    call scalar_elem_product(u, hQc_ed, uhQc) !Flux uhQc at edges

    !Calculate divergence / cloud eq RHS  
    call div_hx(uhQc, div_uhQc, mesh)

    !Cloud eq. RHS
    cloudeq = -div_uhQc%f

    !===============================================================
    !Calculate rain tendency
    !===============================================================

    !Interpolate vapour to edges and calculate flux at edges
    call scalar_hx2ed(hQr, hQr_ed, mesh)      !hQr: cell->edge
    call scalar_elem_product(u, hQr_ed, uhQr) !Flux uhQr at edges

    !Calculate divergence / rain eq RHS
    call div_hx(uhQr, div_uhQr, mesh)

    !Rain eq. RHS
    raineq = -div_uhQr%f

    !===============================================================
    !Compute and add the source
    !===============================================================
    call source(h, u, htheta, hQv, hQc, hQr, dt)

    momeq    = momeq    + Su%f 
    tempeq   = tempeq   + hStheta%f
    vapoureq = vapoureq + hSQv%f
    cloudeq  = cloudeq  + hSQc%f
    raineq   = raineq   + hSQr%f  
    return
  end subroutine tendency_moist_swm

  subroutine source(h, u, htheta, hQv, hQc, hQr, dtime)
    !Fluid thickness (defined on voronoi centers)
    type(scalar_field), intent(in):: h  !General

    !Velocities (defined on edges - only normal component)
    type(scalar_field), intent(in):: u  !General

    !Temperature (defined on voronoi centers)
    type(scalar_field), intent(in):: htheta  !General

    !Vapour (defined on voronoi centers)
    type(scalar_field), intent(in):: hQv  !General

    !Cloud (defined on voronoi centers)
    type(scalar_field), intent(in):: hQc  !General

    !Rain (defined on voronoi centers)
    type(scalar_field), intent(in):: hQr  !General

    !time step
    real(r8):: dtime
    integer(i4):: i

    !===============================================================
    !Calculate momentum source that comes from physics
    !===============================================================

    !Compute h^2*theta and its gradient
    call scalar_elem_product(h, htheta, h2theta)
    call grad_ed(h2theta, gradPI, mesh)
    call scalar_elem_divide(gradPI, h_ed, gradPI_oh)

    !Compute the topography gradient and multiply it by theta 
    call grad_ed(bt, grad_b, mesh)
    call scalar_hx2ed(theta, theta_ed, mesh)      !htheta: cell->edge
    call scalar_elem_product(theta_ed, grad_b, theta_grad_b)

    !Source for momentum equation
    Su%f = gravo2*gradPI_oh%f + grav*theta_grad_b%f

    !===============================================================
    !Calculate temperature, vapour, cloud and rain sources 
    !===============================================================

    call scalar_elem_divide(htheta, h, theta)
    call scalar_elem_divide(hQv, h, Qv)
    call scalar_elem_divide(hQc, h, Qc)
    call scalar_elem_divide(hQr, h, Qr)

    !$omp parallel do &
    !$omp default(none) &
    !$omp shared(mesh) &
    !$omp shared(delta_Qv, delta_Qc, delta_Qr) &
    !$OMP shared(hStheta, hSqv, hSqc, hSqr) &
    !$omp shared(Qv, Qc, Qr) &
    !$omp shared(theta, h, bt) &
    !$omp shared(dtime, gamma_v, gamma_r, q0, q_precip, Lscale) &
    !$omp schedule(static)

    do i=1, mesh%nv
      gamma_v = 1._r8+Lscale*dqsat_dtheta(theta%f(i),h%f(i), bt%f(i), q0)
      gamma_v = 1._r8/gamma_v  

      !Vapour
      delta_Qv%f(i) = gamma_v*(qv%f(i)-qsat(theta%f(i),h%f(i), bt%f(i), q0))
      delta_Qv%f(i) = max(0._r8,delta_Qv%f(i))/dtime

      !Cloud
      delta_Qc%f(i) = max(0._r8, -gamma_v*(qv%f(i)-qsat(theta%f(i),h%f(i), bt%f(i), q0)))
      delta_Qc%f(i) = min(qc%f(i),delta_Qc%f(i))/dtime

      !Rain
      delta_Qr%f(i) = max(0._r8,gamma_r*(Qc%f(i)-q_precip))/dtime

      !Computes the sources
      hSqv%f(i) = delta_Qc%f(i) - delta_Qv%f(i)
      hSqc%f(i) = delta_Qv%f(i) - delta_Qc%f(i) - delta_Qr%f(i)
      hSqr%f(i) = delta_Qr%f(i)
      hStheta%f(i) = Lscale*(delta_Qv%f(i)-delta_Qc%f(i))

      hSQv%f(i) = h%f(i)*hSQv%f(i)
      hSQc%f(i) = h%f(i)*hSQc%f(i)
      hSQr%f(i) = h%f(i)*hSQr%f(i)
      hStheta%f(i) =  h%f(i)*hStheta%f(i)
    end do

    !$omp end parallel do
  end subroutine source

  subroutine moist_swm_tests(meshtmp)
    !-----------------------------------------
    !  Main test routine tests routine
    !-----------------------------------------
    !Grid structure (incomming)
    type(grid_structure) :: meshtmp

    !Indexes
    integer(i4):: i !For node values
    integer(i4):: k !For triangles
    integer(i4):: l !For edges

    !Time in seconds
    real(r8)::time
    
    !Check for blow ups
    integer(i4)::blowup=0

    !Error variables - for tc2
    real(r8):: rel_error_h
    real(r8):: rel_error_u
    real(r8):: rel_error_Qv
    real(r8):: rel_error_Qc
    real(r8):: rel_error_Qr
    real(r8):: rel_error_theta

    !Total mass of cloud, rain and vapour
    real(r8):: Tcloud
    real(r8):: Train
    real(r8):: Tvapour

    !Save global variable mesh
    mesh=meshtmp

    !Get test case parameters
    call swm_phys_pars(usetime=.true.)

    !Allocate variables
    call allocate_global_moistswm_vars()

    !Pre calculate grid properties
    call initialize_gridprop()
   
    !Initialize fields
    call initialize_moist_swm_fields()

    !Compute the diffusion coeficient
    if(K2_max>0.d0)then    
        select case(diffus)
            case("const") !constant diffusion
                dif_coef_hx%f = K2_max
            case("align") !align based coefficient
                call alignment_coef(K2_max, dif_coef_hx, mesh)
            case("diam")  !diameter based coefficient
                call diameter_coef(K2_max, dif_coef_hx, mesh)
            case default
                print*,'ERROR: invalid diffusion coefficient function'
                stop
        end select
        
        !Interpolate to edges and triangles
        call scalar_hx2ed(dif_coef_hx, dif_coef_ed, mesh) !cell->edge
        call scalar_hx2tr(dif_coef_hx, dif_coef_tr, mesh) !cell->tr             
    end if
    
    !Compute the hyperdiffusion coef
    if(K4_max>0.d0)then
        select case(hyperdiffus)
            case("const") !constant diffusion
                hypdif_coef_hx%f = K4_max
            case("align") !align based coefficient
                call alignment_coef(K4_max, hypdif_coef_hx, mesh)
            case("diam")  !diameter based coefficient
                call diameter_coef(K4_max,hypdif_coef_hx, mesh)
            case default
                print*,'ERROR: invalid hyperdiffusion coefficient function'
                stop
        end select
        
        !Interpolate the hyperdiffusion coefficient from hx to ed
        call scalar_hx2ed(hypdif_coef_hx, hypdif_coef_ed, mesh) !cell->edge
        call scalar_hx2tr(hypdif_coef_hx, hypdif_coef_tr, mesh) !cell->tr
        
        dif_coef_hx%f = dsqrt(hypdif_coef_hx%f)
        dif_coef_ed%f = dsqrt(hypdif_coef_ed%f)
        dif_coef_tr%f = dsqrt(hypdif_coef_tr%f)
        !stop        
    end if
    
    !Calculate derived initial fields
    call tendency_moist_swm(h, u, htheta, hQv, hQc, hQr, masseq%f, momeq%f, tempeq%f, vapoureq%f, cloudeq%f, raineq%f)

    !Plot initial fields
    call plotfields_mswm(0, 0._r8)
    !stop 

    !Calculate total mass
    hwater%f = hQv%f + hQc%f + hQr%f
    inimass=sumf_areas(h)
    iniwater=sumf_areas(hwater)

    !Calculate energies
    call calc_energies(Penergy0, Kenergy0, Tenergy0, Availenergy0)

    u_old=u
    h_old=h
    htheta_old=htheta
    hQv_old=hQv
    hQc_old=hQc
    hQr_old=hQr
    
    !Time loop
    do k=1, ntime
      !Calculate u and h for time:
      time=real(k, r8)*dt
      call ode_rk4_moist_swm(time, h_old, u_old, htheta_old, hQv_old, hQc_old, hQr_old, &
                         h, u, htheta, hQv, hQc, hQr, dt)

      !Apply the monotonic filter for tracers
      call monotonic_filter(hQv)             
      call monotonic_filter(hQr)
      call monotonic_filter(hQc)

      !call mass_fixer(hQv,hQr,hQc, iniwater)
      !call scalar_elem_divide(htheta, h, theta)
      call scalar_elem_divide(hQv, h, Qv)
      call scalar_elem_divide(hQc, h, Qc)
      call scalar_elem_divide(hQr, h, Qr)

      !compute the mass of each tracer
      Train=sumf_areas(qr)
      Tcloud=sumf_areas(qc)
      Tvapour=sumf_areas(qv)
      
      if(testcase==2 .or. testcase==3)then
        h_error%f = h_exact%f - h%f
        u_error%f = u_exact%f - u%f
        Qv_error%f = Qv_exact%f - Qv%f
        Qc_error%f = Qc_exact%f - Qc%f
        Qr_error%f = Qr_exact%f - Qr%f
        theta_error%f = theta_exact%f - theta%f
        rel_error_h = maxval(abs(h_error%f))/maxval(abs(h_exact%f))
        rel_error_u = maxval(abs(u_error%f))/maxval(abs(u_exact%f))
        rel_error_Qv = maxval(abs(Qv_error%f))/maxval(abs(Qv_exact%f))
        rel_error_Qc = maxval(abs(Qc_error%f))
        rel_error_Qr = maxval(abs(Qr_error%f))
        rel_error_theta = maxval(abs(theta_error%f))/maxval(abs(theta_exact%f))
        !print*, k, ntime
        print*, "Time (dys) :",   k*dt*sec2day, " of ", ntime*dt*sec2day
        print*, "Step = ", k, " of ", ntime
        print '(a33, 3e16.8)','linf errors of (h, u, theta) = ',rel_error_h,rel_error_u,rel_error_theta
        print '(a33, 3e16.8)','linf errors of  (qv, qc, qr) = ',rel_error_Qv,rel_error_Qc,rel_error_Qr
        !print '(a22, 2e16.8)',' height = ',minval(h%f),maxval(h%f)
        !print '(a22, 2e16.8)',' velocity = ',minval(u%f),maxval(u%f)
        !print '(a22, 2e16.8)',' temperature = ',minval(theta%f),maxval(theta%f)
        !print '(a22, 2e16.8)',' vapour = ',minval(qv%f),maxval(qv%f)
        !print '(a22, 2e16.8)',' cloud = ',minval(qc%f),maxval(qc%f)
        !print '(a22, 2e16.8)',' rain = ',minval(qr%f),maxval(qr%f)
        call write_swmp_error_file(time)
      else 
        print*, "Time (dys) :",   k*dt*sec2day, " of ", ntime*dt*sec2day
        print*, "Step = ", k, " of ", ntime
        print*,'                          min               max               mass'
        print '(a22, 2e16.8)',' height = ',minval(h%f),maxval(h%f)
        print '(a22, 2e16.8)',' velocity = ',minval(u%f),maxval(u%f)
        print '(a22, 2e16.8)',' temperature = ',minval(theta%f),maxval(theta%f)
        print '(a22, 3e16.8)',' vapour = ',minval(qv%f),maxval(qv%f), Tvapour
        print '(a22, 3e16.8)',' cloud = ',minval(qc%f),maxval(qc%f), Tcloud
        print '(a22, 3e16.8)',' rain = ',minval(qr%f),maxval(qr%f), Train
      end if

      !print*,'CFL = ',cfl
      !Plot fields
      call plotfields_mswm(k, time)

      !Write errors in file
      call write_evol_file_cswm(time, iniwater, Twater, inimass, Penergy0, Kenergy0, Tenergy0, Availenergy0,&
                                   tmass, Penergy, Kenergy, Tenergy, Availenergy)

      call write_water_evol_file(time, Train, Tcloud, Tvapour)
      
      !Calculate total mass
      Tmass=sumf_areas(h)

      !Calculate total water
      hwater%f = hQr%f + hQv%f + hQc%f
      
      !call scalar_elem_divide(hwater, h, water)
      Twater = sumf_areas(hwater)

      !Calculate erngies
      call calc_energies(Penergy, Kenergy, Tenergy, Availenergy)
      print '(a33, 2e16.8)','Change in mass of h*(total water):', (Twater-iniwater)/iniwater
      print*,''
      !update fields
      u_old=u
      h_old=h
      htheta_old=htheta
      hQv_old=hQv
      hQc_old=hQc
      hQr_old=hQr
    end do
  end  subroutine moist_swm_tests


  subroutine write_evol_file_cswm(time, iniwater, Twater, inimass, Penergy0, Kenergy0, Tenergy0, Availenergy0,&
    tmass, Penergy, Kenergy, Tenergy, Availenergy)
  !----------------------------------------------------------
  !  write info to specific file for this specific model set up
  !    at defined time steps
  !----------------------------------------------------------
  !File name for output
  character (len=256):: filename

  !File units
  integer (i4):: iunit
  logical::  iopen

  !inputs
  real(r8), intent(in):: time
  real(r8), intent(in):: iniwater, Twater
  real(r8), intent(in):: inimass,  Tmass
  real(r8), intent(in):: Tenergy0, Tenergy
  real(r8), intent(in):: Kenergy0, Kenergy
  real(r8), intent(in):: Penergy0, Penergy
  real(r8), intent(in):: Availenergy0, Availenergy

  integer(i4) :: k,l
  integer(i4) :: errorsunit
  logical :: ifile

  !File for errors
  filename=trim(datadir)// trim(swmname)//"_"//trim(mesh%name)//"_evolution.txt"
  call getunit(errorsunit)
  inquire(file=filename, exist=ifile)
  if(ifile)then
    if(k>0)then
      open(errorsunit,file=filename, status='replace')
    else
      open(errorsunit,file=filename, status='old', position='append')
    end if
  else
    open(errorsunit,file=filename, status='replace')
  end if
  
  write(errorsunit, *) time*sec2day, Tenergy, Kenergy, Penergy, Availenergy, Tmass, Twater, &
  (Tenergy-Tenergy0)/Tenergy0, &
  (Kenergy-Kenergy0)/Kenergy0, &
  (Penergy-Penergy0)/Penergy0, &
  (Availenergy-Availenergy0)/Availenergy0, &
  (Tmass-inimass)/inimass, &
  (Twater-iniwater)/iniwater
  close(errorsunit)
end subroutine write_evol_file_cswm

subroutine write_water_evol_file(time, train, tcloud, tvapour)
!----------------------------------------------------------
!  write info about water evolution in the moist shallow water model
!  to specific file for this specific model set up   at defined time steps
!----------------------------------------------------------
!File name for output
character (len=256):: filename

!File units
integer (i4):: iunit
logical::  iopen

!inputs
real(r8), intent(in):: time
real(r8), intent(in):: train, tcloud, tvapour

integer(i4) :: k,l
integer(i4) :: errorsunit
logical :: ifile

!File for errors
filename=trim(datadir)// trim(swmname)//"_"//trim(mesh%name)//"_evolution_water.txt"
call getunit(errorsunit)
inquire(file=filename, exist=ifile)
if(ifile)then
  if(k>0)then
    open(errorsunit,file=filename, status='replace')
  else
    open(errorsunit,file=filename, status='old', position='append')
  end if
else
  open(errorsunit,file=filename, status='replace')
end if

write(errorsunit, *) time*sec2day, train, tcloud, tvapour
 close(errorsunit)
end subroutine write_water_evol_file

  
subroutine plotfields_mswm(k, time)
    !-------------------------------------------
    !  Plot fields
    !  k- index for time couting
    !  time - current time step
    !-------------------------------------------

    integer(i4)::k !Time index
    integer(i4)::i
    real(r8)::time
    character (len=60)::  atime

    write(atime,*) nint(time)

    if(.not.plots)return

    if( (k==ntime  .or. mod(k,plotsteps)==0 ) )then
        !Scalar field plots
        h%name=trim(swmname)//"_h_t"//trim(adjustl(trim(atime)))
        call plot_scalarfield(h, mesh)

        do i=1, mesh%nv
          call convert_vec_cart2sph(mesh%v(i)%p, v_hx%p(i)%v, ueast%f(i), unorth%f(i))
        end do

        !ueast%name=trim(swmname)//"_ueast_t"//trim(adjustl(trim(atime)))
        !call plot_scalarfield(ueast, mesh)

        !unorth%name=trim(swmname)//"_unorth_t"//trim(adjustl(trim(atime)))
        !call plot_scalarfield(unorth, mesh)

        theta%name=trim(swmname)//"_theta_t"//trim(adjustl(trim(atime)))
        call plot_scalarfield(theta, mesh)

        qv%name=trim(swmname)//"_qv_t"//trim(adjustl(trim(atime)))
        call plot_scalarfield(qv, mesh)

        qc%name=trim(swmname)//"_qc_t"//trim(adjustl(trim(atime)))
        call plot_scalarfield(qc, mesh)

        qr%name=trim(swmname)//"_qr_t"//trim(adjustl(trim(atime)))
        call plot_scalarfield(qr, mesh)

        if(maxval(bt%f(1:bt%n)) > eps)then
          hbt%name=trim(swmname)//"_hbt_t"//trim(adjustl(trim(atime)))
          call plot_scalarfield(hbt, mesh)

          bt%name=trim(swmname)//"_bt_t"//trim(adjustl(trim(atime)))
          call plot_scalarfield(bt, mesh)
        end if
     
        if(testcase==2 .or. testcase==3)then
          h_error%name=trim(swmname)//"_h_error_t"//trim(adjustl(trim(atime)))
          call plot_scalarfield(h_error, mesh)

          u_error%name=trim(swmname)//"_u_error_t"//trim(adjustl(trim(atime)))
          call plot_scalarfield(u_error, mesh)

          theta_error%name=trim(swmname)//"_theta_error_t"//trim(adjustl(trim(atime)))
          call plot_scalarfield(theta_error, mesh)

          qv_error%name=trim(swmname)//"_qv_error_t"//trim(adjustl(trim(atime)))
          call plot_scalarfield(qv_error, mesh)
        end if
        !eta%name=trim(swmname)//"_eta_t"//trim(adjustl(trim(atime)))
        !call plot_scalarfield(eta, mesh)

        !zeta%name=trim(swmname)//"_zeta_t"//trim(adjustl(trim(atime)))
        !call plot_scalarfield(zeta, mesh)

        !ke_hx%name=trim(swmname)//"_Kenergy_t"//trim(adjustl(trim(atime)))
        !call plot_scalarfield(ke_hx, mesh)

        q_tr%name=trim(swmname)//"_pv_t"//trim(adjustl(trim(atime)))
        call plot_scalarfield(q_tr, mesh)

        q_ed%name=trim(swmname)//"_pv_ed_t"//trim(adjustl(trim(atime)))
        call plot_scalarfield(q_ed, mesh)

        !divuh%name=trim(swmname)//"_divuh_t"//trim(adjustl(trim(atime)))
        !call plot_scalarfield(divuh, mesh)
    end if
  end subroutine plotfields_mswm

  subroutine ode_rk4_moist_swm ( t, h, u, htheta, hqv, hqc, hqr, h_new, u_new, htheta_new,&
                                 hqv_new, hqc_new, hqr_new, dt)
    !----------------------------------------------------------------------------------
    !! ode_rk4 takes one Runge-Kutta step for a vector ODE.
    !    t - time that will be calculated (t0+dt)
    !    h - scalar_field for thickness at current time
    !    u - scalar_field for velocities at current time
    !    dt - time step
    !    h_new and u_new - fields at t+dt
    !----------------------------------------------------------------------------------

    !Fluid thickness (defined on voronoi centers)
    type(scalar_field), intent(in):: h  !General

    !Velocities (defined on edges - only normal component)
    type(scalar_field), intent(in):: u  !General

    !Temperature (defined on voronoi centers) 
    type(scalar_field), intent(in):: htheta  !General

    !Vapour (defined on voronoi centers)
    type(scalar_field), intent(in):: hQv  !General

    !Cloud (defined on voronoi centers)
    type(scalar_field), intent(in):: hQc !General

    !Rain (defined on voronoi centers)
    type(scalar_field), intent(in):: hQr !General

    !Time and time step
    real(r8):: t
    real(r8):: dt

    !Updated fields
    !Fluid thickness (defined on voronoi centers)
    type(scalar_field):: h_new  !General

    !Velocities (defined on edges - only normal component)
    type(scalar_field):: u_new  !General

    !Temperature (defined on voronoi centers)
    type(scalar_field):: htheta_new  !General

    !Vapour (defined on voronoi centers)
    type(scalar_field):: hqv_new  !General

    !Cloud (defined on voronoi centers)
    type(scalar_field):: hqc_new  !General

    !Rain (defined on voronoi centers)
    type(scalar_field):: hqr_new  !General

    !Times
    real(r8):: t0
    real(r8):: t1
    real(r8):: t2
    real(r8):: t3

    u_new      = u
    h_new      = h
    htheta_new = htheta
    hQv_new    = hQv
    hQc_new    = hQc
    hQr_new    = hQr

    masseq%f   = massf0
    momeq%f    = momf0
    tempeq%f   = tempf0
    vapoureq%f = vapourf0
    cloudeq%f  = cloudf0
    raineq%f   = rainf0

    !Initial f (f0)
    t0=t-dt
    call tendency_moist_swm(h, u, htheta, hQv, hQc, hQr, massf0, momf0, tempf0, vapourf0, cloudf0, rainf0)

    !First RK step
    t1 = t0 + dt/2._r8

    u_new%f(1:u%n)          = u%f(1:u%n)           + dt * momf0(1:u%n) / 2.0_r8
    h_new%f(1:h%n)          = h%f(1:h%n)           + dt * massf0(1:h%n) / 2.0_r8
    htheta_new%f(1:theta%n) = htheta%f(1:theta%n)  + dt * tempf0(1:theta%n) / 2.0_r8
    hQv_new%f(1:qv%n)       = hQv%f(1:qv%n)        + dt * vapourf0(1:qv%n) / 2.0_r8
    hQc_new%f(1:qc%n)       = hQc%f(1:qc%n)        + dt * cloudf0(1:qc%n) / 2.0_r8
    hQr_new%f(1:qr%n)       = hQr%f(1:qr%n)        + dt * rainf0(1:qr%n) / 2.0_r8

    call tendency_moist_swm(h_new, u_new, htheta_new, hQv_new, hQc_new, hQr_new, &
    massf1, momf1, tempf1, vapourf1, cloudf1, rainf1)

    !Second RK step
    t2 = t0 + dt/2._r8

    u_new%f(1:u%n)          = u%f(1:u%n)           + dt * momf1(1:u%n) / 2.0_r8
    h_new%f(1:h%n)          = h%f(1:h%n)           + dt * massf1(1:h%n) / 2.0_r8
    htheta_new%f(1:theta%n) = htheta%f(1:theta%n)  + dt * tempf1(1:theta%n) / 2.0_r8
    hQv_new%f(1:qv%n)       = hQv%f(1:qv%n)        + dt * vapourf1(1:qv%n) / 2.0_r8
    hQc_new%f(1:qc%n)       = hQc%f(1:qc%n)        + dt * cloudf1(1:qc%n) / 2.0_r8
    hQr_new%f(1:qr%n)       = hQr%f(1:qr%n)        + dt * rainf1(1:qr%n) / 2.0_r8

    call tendency_moist_swm(h_new, u_new, htheta_new, hQv_new, hQc_new, hQr_new, &
    massf2, momf2, tempf2, vapourf2, cloudf2, rainf2)


    !Third  RK step
    t3 = t0 + dt
    u_new%f(1:u%n)          = u%f(1:u%n)           + dt * momf2(1:u%n)
    h_new%f(1:h%n)          = h%f(1:h%n)           + dt * massf2(1:h%n) 
    htheta_new%f(1:theta%n) = htheta%f(1:theta%n)  + dt * tempf2(1:theta%n)
    hQv_new%f(1:qv%n)       = hQv%f(1:qv%n)        + dt * vapourf2(1:qv%n) 
    hQc_new%f(1:qc%n)       = hQc%f(1:qc%n)        + dt * cloudf2(1:qc%n)
    hQr_new%f(1:qr%n)       = hQr%f(1:qr%n)        + dt * rainf2(1:qr%n)

    call tendency_moist_swm(h_new, u_new, htheta_new, hQv_new, hQc_new, hQr_new, &
    massf3, momf3, tempf3, vapourf3, cloudf3, rainf3)

    !
    ! Combine them to estimate the solution at time t+dt
    !
    u_new%f(1:u%n) = u%f(1:u%n) + dt * (momf0(1:u%n)+2._r8*momf1(1:u%n) &
    +2._r8*momf2(1:u%n)+momf3(1:u%n))/6._r8

    h_new%f(1:h%n) = h%f(1:h%n) + dt * (massf0(1:h%n)+2._r8*massf1(1:h%n) &
    +2._r8*massf2(1:h%n)+massf3(1:h%n))/6._r8

    htheta_new%f(1:theta%n) = htheta%f(1:theta%n) + dt * (tempf0(1:theta%n)+2._r8*tempf1(1:theta%n) &
    +2._r8*tempf2(1:theta%n)+tempf3(1:theta%n))/6._r8

    hQv_new%f(1:Qv%n) = hQv%f(1:qv%n) + dt * (vapourf0(1:qv%n)+2._r8*vapourf1(1:qv%n) &
    +2._r8*vapourf2(1:qv%n)+vapourf3(1:qv%n))/6._r8

    hQc_new%f(1:Qc%n) = hQc%f(1:qc%n) + dt * (cloudf0(1:qc%n)+2._r8*cloudf1(1:qc%n) &
    +2._r8*cloudf2(1:qc%n)+cloudf3(1:qc%n))/6._r8

    hQr_new%f(1:Qr%n) = hQr%f(1:qr%n) + dt * (rainf0(1:qr%n)+2._r8*rainf1(1:qr%n) &
    +2._r8*rainf2(1:qr%n)+rainf3(1:qr%n))/6._r8
    return
  end subroutine ode_rk4_moist_swm


  subroutine swm_phys_pars(usetime)
    !---------------------------------------------------
    ! swmpars
    !    Reads swm test parameters from file named "swm.par"
    !    Saves parameters on global variables
    !
    !--------------------------------------------------

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

    !Loggical if time variables need or not to be set up
    logical::usetime

    !Temporary Integer fo nopv flag
    integer::inoPV

    !Couters
    !integer(i4)::i
    !integer(i4)::j

    !Standard definition of the deformal tests
    testcase=2      !Williamson 1992
    plots=.true.   !Do plots or not
    nplots=20      !number of plots to output
    ntime= 15 * (2**(mesh%glevel-3)) !60  !Number time steps
    adjustntime=0   !Adjust time steps to grid level 0 or 1


    !Standard parameters file
    filename=trim(pardir)//"moist_swm.par"
    print*,"Shallow Water Model parameters (file): ", trim(filename)
    print*
    call getunit(fileunit)

    !A parameters file must exist
    open(fileunit,file=filename,status='old')

    read(fileunit,*)  buffer
    read(fileunit,*)  buffer
    read(fileunit,*)  testcase, test_lterror
    read(fileunit,*)  buffer
    read(fileunit,*)  period, maxtime
    read(fileunit,*)  buffer
    read(fileunit,*)  dt, ntime, adjustntime
    read(fileunit,*)  buffer
    read(fileunit,*)  stag
    read(fileunit,*)  buffer
    read(fileunit,*)  mtdwrapper
    read(fileunit,*)  buffer
    read(fileunit,*)  reconmtd, gasscoef
    read(fileunit,*)  buffer
    read(fileunit,*)  coriolis_reconmtd
    read(fileunit,*)  buffer
    read(fileunit,*)  sinterpol
    read(fileunit,*)  buffer
    read(fileunit,*)  gradmtd
    read(fileunit,*)  buffer
    read(fileunit,*)  areamtd
    read(fileunit,*)  buffer
    read(fileunit,*)  nplots, nprints, iploterrors
    read(fileunit,*)  buffer
    read(fileunit,*)  K2_max, diffus
    read(fileunit,*)  buffer
    read(fileunit,*)  K4_max, hyperdiffus
    read(fileunit,*)  buffer
    read(fileunit,*)  hollgw
    read(fileunit,*)  buffer
    read(fileunit,*)  trskindmax
    read(fileunit,*)  buffer
    read(fileunit,*)  pv_stab, pvspar
    read(fileunit,*)  buffer
    read(fileunit,*)  inoPV
    close(fileunit)

    !Flag to test linearized version of equation
    ! used for normal modes calculation
    fsphere=0

    !Check staggering
    useStagHTC=trim(stag)=="HTC"
    useStagHC=trim(stag)=="HC"
    if((.not.useStagHC) .and. (.not.useStagHTC))then
      print*, "Unknown staggering", stag
      stop
    end if

    if(.not.(trim(mtdwrapper)=="trsk10".or. trim(mtdwrapper)=="pxt16".or.trim(mtdwrapper)=="gass18"))then
      if(trim(mtdwrapper)/="none")then
        print*, "No valid method wrapper selected (trsk10, pxt16, gass18)", mtdwrapper
        mtdwrapper="none"
      endif
    endif

    !Check scalar interpol
    if(trim(mtdwrapper)=="trsk10")then
      useSinterpolTrisk=.true.
      sinterpol="trsk"
    elseif(trim(mtdwrapper)=="pxt16")then
      useSinterpolBary=.true.
      sinterpol="bary"
    elseif(trim(mtdwrapper)=="gass18")then
      useSinterpolGass=.true.
      sinterpol="gass"
    else
      useSinterpolTrisk=trim(sinterpol)=="trsk"
      useSinterpolBary=trim(sinterpol)=="bary"
      useSinterpolGass=trim(sinterpol)=="gass"
    endif
    if((.not.useSinterpolTrisk).and.(.not.useSinterpolBary).and.(.not.useSinterpolGass))then
      print*, "Unknown interpolation", sinterpol, mtdwrapper
      stop
    end if

    !Check  Vector reconstruction method
    if(trim(mtdwrapper)=="trsk10")then
      useCoriolisMtdTrisk=.true.
      coriolis_reconmtd="trsk"
    elseif(trim(mtdwrapper)=="pxt16")then
      useCoriolisMtdPered=.true.
      coriolis_reconmtd="pered"
    elseif(trim(mtdwrapper)=="gass18")then
      useCoriolisMtdGass=.true.
      coriolis_reconmtd="gass"
    else
      useCoriolisMtdPered=trim(coriolis_reconmtd)=="pered"
      useCoriolisMtdDtred=trim(coriolis_reconmtd)=="dtred"
      useCoriolisMtdTrisk=trim(coriolis_reconmtd)=="trsk"
      useCoriolisMtdHyb=trim(coriolis_reconmtd)=="hyb"
      useCoriolisMtdGass=trim(coriolis_reconmtd)=="gass"
      useCoriolisMtdExact=trim(coriolis_reconmtd)=="exact"
    endif
    if((.not. useCoriolisMtdPered).and.(.not. useCoriolisMtdTrisk) &
      .and.(.not. useCoriolisMtdDtred).and.(.not. useCoriolisMtdHyb) &
      .and.(.not.useCoriolisMtdGass).and.(.not. useCoriolisMtdExact))then
      print*, "Unknown Coriolis vector reconstruction method", coriolis_reconmtd, mtdwrapper
      stop
    end if

    noPV=inoPV>0
    if(noPV)then
      if(useCoriolisMtdDtred.or.useCoriolisMtdPered)then
        print*, "Cannot use this Coriolis vector reconstruction method with level model (noPV)", coriolis_reconmtd, noPV
        stop
      endif
    endif

    !Check  Vector reconstruction method
    if(trim(mtdwrapper)=="trsk10")then
      useReconmtdTrisk=.true.
      reconmtd="trsk"
    elseif(trim(mtdwrapper)=="pxt16")then
      useReconmtdPerhx=.true.
      reconmtd="perhx"
    elseif(trim(mtdwrapper)=="gass18")then
      useReconmtdTrisk=.true.
      useReconmtdGass=.false.
      reconmtd="trsk"
    else
      useReconmtdPerhx=trim(reconmtd)=="perhx"
      useReconmtdTrisk=trim(reconmtd)=="trsk"
      useReconmtdGass=trim(reconmtd)=="gass"
      useReconmtdMelv=trim(reconmtd)=="melv"
      if(trim(reconmtd)=="dubos")then
        useReconmtdGass=.true.
        gasscoef=0.0
      end if
    endif
    if((.not. useReconmtdPerhx).and.(.not. useReconmtdTrisk) &
      .and.(.not. useReconmtdGass).and.(.not. useReconmtdMelv))then
      print*, "Unknown vector reconstruction method", reconmtd
      stop
    end if

    !Check gradmetd
    if (trim(mtdwrapper)=="none") then

      useGradmtdTrisk=trim(gradmtd)=="trsk"
      useGradmtdBary=trim(gradmtd)=="bary"
    elseif(trim(mtdwrapper)=="trsk10" .or. trim(mtdwrapper)=="pxt16" .or. trim(mtdwrapper)=="gass18")then
      useGradmtdTrisk=.true.
    endif
    if((.not.useGradmtdTrisk).and.(.not.useGradmtdBary))then
      print*, "Unknown gradient discretization method", gradmtd
      stop
    end if

    !Check areamtd
    useGeoAreas=trim(areamtd(1:3))=="geo"
    useTiledAreas=trim(areamtd(1:4))=="tile"
    usePlannarAreas=.false. !trim(areamtd(1:4))=="plan" !not implmented
    if((.not.useGeoAreas).and.(.not.useTiledAreas).and.(.not.usePlannarAreas))then
      print*, "Area not well specified, using geodesic areas", areamtd
      useGeoAreas=.true.
    end if

    !Pontential Vorticity correction
    useAPVM=trim(pv_stab)=="apvm"
    useCLUST=trim(pv_stab)=="clust"
    useOrigPV=trim(pv_stab)=="none"
    if((.not.useCLUST).and.(.not.useAPVM).and.(.not.useOrigPV))then
      print*, "Not using any correction for PV", pv_stab
    end if
    if(useCLUST)then
      if(pvspar<=0)then
        print*, "Please provide a parameter (0-1) for CLUST in swm.par", pvspar
        stop
      end if
    end if

    !Number of time steps
    if(period <= 0)then
      stop "swmpars error: period should be positive and must be given in days"
    endif
    !period and maxtime should be inputed in days, so we have to convert to seconds
    period=period*day2sec
    maxtime=maxtime*day2sec

    if(dt>0)then !use the given time step
      if(adjustntime == 1)then
        dt = dt * (2**(5_i4-mesh%glevel))
      end if
      ntime=ceiling(maxtime/dt,i4)
       !print*, dt, ntime
    elseif(dt==0 .and. ntime>0)then !use number of time steps
      if(adjustntime == 1)then
        ntime = ntime * (2**(mesh%glevel-5_i4))
      end if
      dt=maxtime/real(ntime, r8)
    elseif(ntime<=0 .and. dt<=0)then !use dt=50s
      dt=50 !seconds
      if(adjustntime == 1)then
        dt = dt * (2**(5_i4-mesh%glevel))
      end if
      ntime=ceiling(maxtime/dt,i4)
    endif

    !Set number of times to plot
    if(nplots<=0) then
      plots=.false.
    else
      plotsteps=ntime/nplots
    end if
    if(plotsteps<=0)then
      !ntime too small or nplots too large - plot every timestep
      plotsteps=1
    end if

    if(nprints<=0)then
      nprints=1000000
    end if

    if(iploterrors<=0)then
      ploterrors=.false.
    else
      ploterrors=.true.
    end if

    !Set a standart name for files
    write(atmp,'(i8)') int(testcase)
    swmname="moist_swm_tc"//trim(adjustl(trim(atmp)))

    !Print information
    print*, "Test Case     : ", testcase
    if(usetime)then
      print*, "Integration period (dys): ", period*sec2day
      print*, "Stopping time (dys)     : ", maxtime*sec2day
      print*, "dt  (sec)               : ",  dt
      print*, "Number of timesteps     : ", ntime
      write(atmp,'(i8)') nint(dt)
      swmname=trim(swmname)//"_dt"//trim(adjustl(trim(atmp)))
    else
      ntime=1
    end if
    print*, "Staggering              : ", stag
    print*, "Method wrapper          : ", mtdwrapper
    print*, "Scalar interpolation    : ", sinterpol
    print*, "Vector recon method     : ", reconmtd
    print*, "Coriolis recon method   : ", coriolis_reconmtd
    print*, "Gradient method         : ", gradmtd
    print*, "Area method             : ", areamtd
    print*, "Hollingsworth depth     : ", hollgw
    print*, "PV stabilization mtd    : ", pv_stab
    if(noPV)then
      print*, "Using level model (no PV)  "
    else
      print*, "Using layer model (with PV)  "
    end if
    print*
    swmname=trim(swmname)//"_"//trim(adjustl(trim(stag)))

    ! If wrapper set, shorten name
    if(trim(mtdwrapper)/="none")then
      swmname=trim(swmname)//"_"//trim(adjustl(trim(mtdwrapper)))
    else
      swmname=trim(swmname)//"_vrec"//trim(adjustl(trim(reconmtd)))
      if(useReconmtdGass)then
        write(atmp,'(f4.2)') real(gasscoef)
        swmname=trim(swmname)//trim(adjustl(trim(atmp)))
      end if
      swmname=trim(swmname)//"_crec"//trim(adjustl(trim(coriolis_reconmtd)))
      swmname=trim(swmname)//"_sint"//trim(adjustl(trim(sinterpol)))
      swmname=trim(swmname)//"_grd"//trim(adjustl(trim(gradmtd)))
    end if

    swmname=trim(swmname)//"_area"//trim(adjustl(trim(areamtd)))

    if((.not.useOrigPV))then
      write(atmp,'(f5.2)') real(pvspar)
      swmname=trim(swmname)//"_pvs"//trim(adjustl(trim(pv_stab)))
      swmname=trim(swmname)//trim(adjustl(trim(atmp)))
    end if

    if(noPV)then
      swmname=trim(swmname)//"_nopv"
    endif

    if(hollgw>0 .and. (testcase==32 .or. testcase==33 .or. testcase==34 .or. testcase==35) )then
      write(atmp,'(f6.2)') real(hollgw)
      swmname=trim(swmname)//"_hol"//trim(adjustl(trim(atmp)))
       !print*, atmp
    end if


    if( K2_max>0 ) then
      write(atmp,'(f10.3)') real(dlog(K2_max)/dlog(10._r8))
      swmname=trim(swmname)//"_"//trim(diffus)//"_diffusion_10to"//trim(adjustl(trim(atmp)))
    end if

    if( K4_max>0 ) then
      write(atmp,'(f10.3)') real(dlog(K4_max)/dlog(10._r8)) 
      swmname=trim(swmname)//"_"//trim(hyperdiffus)//"_hyperdiffusion_10to"//trim(adjustl(trim(atmp)))     
    end if 
    
    !RefSolRead=testcase==5.or. testcase==51.or.testcase==6.or.testcase==21.or.testcase==23.or.testcase==60
    !RefSolAnal= testcase==1.or.testcase==2.or. testcase==22.or. testcase==24 &
    !  .or. testcase==32.or. testcase==33 .or. testcase==34 .or. testcase==35 .or. &
    !  testcase==40 .or. testcase==41.or. testcase==42 .or.testcase==56 .or. testcase==57

    print*, "SWM Name for Plots: ", trim(swmname)
    print*

    return
  end subroutine swm_phys_pars

  subroutine write_swmp_error_file(time)
    !----------------------------------------------------------
    ! Write the errors in a file
    !----------------------------------------------------------
    !Scalar fields
    !type(scalar_field), intent(in):: h, h_exact         !Height field
    !type(scalar_field), intent(in):: u, u_exact          !Velocity - normal component
    type(scalar_field) :: grad_gh,grad_gh_error,grad_gh_exact
    type(scalar_field) :: grad_K,grad_K_error,grad_K_exact
    
    real(r8), intent(in) :: time

    !Error and norm
    real(r8):: errormax_h, errormaxrel_h, error2_h
    real(r8):: errormax_u, errormaxrel_u, error2_u
    real(r8):: errormax_theta, errormaxrel_theta, error2_theta
    real(r8):: errormax_qv, errormaxrel_qv, error2_qv
    real(r8):: errormax_qc, errormaxrel_qc, error2_qc
    real(r8):: errormax_qr, errormaxrel_qr, error2_qr

    integer(i4) :: k,l
    integer(i4) :: errorsunit
    logical :: ifile

    character (len=256):: filename
    
    !Computes the errors
    call error_calc(h, h_exact, h_error, errormaxrel_h, error2_h, errormax_h)
    call error_calc(u, u_exact, u_error, errormaxrel_u, error2_u, errormax_u)
    call error_calc(theta, theta_exact, theta_error, errormaxrel_theta, error2_theta, errormax_theta)
    call error_calc(qv, qv_exact, qv_error, errormaxrel_qv, error2_qv, errormax_qv)
    call error_calc(qc, qc_exact, qc_error, errormaxrel_qc, error2_qc, errormax_qc)
    call error_calc(qr, qr_exact, qr_error, errormaxrel_qr, error2_qr, errormax_qr)
    
    !File for errors
    filename=trim(datadir)//trim(swmname)//"_errors_h_"//trim(mesh%name)//".txt"
    call getunit(errorsunit)
    inquire(file=filename, exist=ifile)
    if(ifile)then
      if(k>0)then
        open(errorsunit,file=filename, status='replace')
      else
        open(errorsunit,file=filename, status='old', position='append')
      end if
    else
      open(errorsunit,file=filename, status='replace')
    end if
    write(errorsunit, *) time*sec2day, errormaxrel_h, error2_h, errormax_h
    close(errorsunit)

    !File for errors
    filename=trim(datadir)//trim(swmname)//"_errors_u_"//trim(mesh%name)//".txt"
    call getunit(errorsunit)
    inquire(file=filename, exist=ifile)
    if(ifile)then
      if(k>0)then
        open(errorsunit,file=filename, status='replace')
      else
        open(errorsunit,file=filename, status='old', position='append')
      end if
    else
      open(errorsunit,file=filename, status='replace')
    end if
    write(errorsunit, *) time*sec2day, errormaxrel_u, error2_u, errormax_u
    close(errorsunit)

        !File for errors
    filename=trim(datadir)//trim(swmname)//"_errors_theta_"//trim(mesh%name)//".txt"
    call getunit(errorsunit)
    inquire(file=filename, exist=ifile)
    if(ifile)then
      if(k>0)then
        open(errorsunit,file=filename, status='replace')
      else
        open(errorsunit,file=filename, status='old', position='append')
      end if
    else
      open(errorsunit,file=filename, status='replace')
    end if
    write(errorsunit, *) time*sec2day, errormaxrel_theta, error2_theta, errormax_theta
    close(errorsunit)


   !File for errors
    filename=trim(datadir)//trim(swmname)//"_errors_qv_"//trim(mesh%name)//".txt"
    call getunit(errorsunit)
    inquire(file=filename, exist=ifile)
    if(ifile)then
      if(k>0)then
        open(errorsunit,file=filename, status='replace')
      else
        open(errorsunit,file=filename, status='old', position='append')
      end if
    else
      open(errorsunit,file=filename, status='replace')
    end if
    write(errorsunit, *) time*sec2day, errormaxrel_qv, error2_qv, errormax_qv
    close(errorsunit)

    !File for errors
    filename=trim(datadir)//trim(swmname)//"_errors_qc_"//trim(mesh%name)//".txt"
    call getunit(errorsunit)
    inquire(file=filename, exist=ifile)
    if(ifile)then
      if(k>0)then
        open(errorsunit,file=filename, status='replace')
      else
        open(errorsunit,file=filename, status='old', position='append')
      end if
    else
      open(errorsunit,file=filename, status='replace')
    end if
    write(errorsunit, *) time*sec2day, errormaxrel_qc, error2_qc, errormax_qc
    close(errorsunit)

    !File for errors
    filename=trim(datadir)//trim(swmname)//"_errors_qr_"//trim(mesh%name)//".txt"
    call getunit(errorsunit)
    inquire(file=filename, exist=ifile)
    if(ifile)then
      if(k>0)then
        open(errorsunit,file=filename, status='replace')
      else
        open(errorsunit,file=filename, status='old', position='append')
      end if
    else
      open(errorsunit,file=filename, status='replace')
    end if
    write(errorsunit, *) time*sec2day, errormaxrel_qr, error2_qr, errormax_qr
    close(errorsunit)

  end subroutine write_swmp_error_file

  
  !================================================
  ! Mononotic filter for tracer phi
  !================================================
  subroutine monotonic_filter(phi) 
    type(scalar_field), intent(inout):: phi
    type(scalar_field) :: phi_mass
    real(kind=8):: mass_water_initial
    integer(i4) :: i,j,k, nnb, iter
    real(r8):: mass_hQc, mass_hQv, mass_hQr, modified_mass_hQr, sumareas, summass,summass2
    real(r8):: eps, eps2

    phi_mass= phi
    !stop
    !$omp parallel do &
    !$omp default(none) &
    !$omp shared(mesh) &
    !$OMP shared(phi,phi_mass)&
    !$omp schedule(static)

    !Compute the mass for each cell
    do i = 1, mesh%nv
      phi_mass%f(i) = phi%f(i)*mesh%hx(i)%areag
    end do
    !$omp end parallel do
    summass = sum(phi_mass%f)

    iter = 0
    !eps = maxval(phi_mass%f)*0.00000001d0
    eps2 = maxval(phi_mass%f)*0.00001d0
    eps = max(eps2,0.000000000001)
    do while (minval(phi_mass%f) < -eps) 
    !!$omp parallel do &
    !!$omp default(none) &
    !!$omp shared(mesh) &
    !!$OMP shared(phi_mass) &
    !!$omp private(nnb,k) &
    !!$omp schedule(static)
      do i = 1, mesh%nv
        if(phi_mass%f(i) < 0.d0)then
          nnb = size(mesh%v(i)%nb(:))
          do j = 1, nnb
            k = mesh%v(i)%nb(j)
            phi_mass%f(k) = phi_mass%f(k) + phi_mass%f(i)/nnb
          end do
          phi_mass%f(i) = 0.d0
        end if
      end do
    !!$omp end parallel do
    iter = iter+1
    end do

    !Compute the new tracer value  for each cell
    !$omp parallel do &
    !$omp default(none) &
    !$omp shared(mesh) &
    !$OMP shared(phi_mass)&
    !$omp schedule(static)

    do i = 1, mesh%nv
      if(phi_mass%f(i) < 0.d0)then
        phi_mass%f(i) = 0.d0
      end if
    end do
    !$omp end parallel do

    !Compute the new tracer value  for each cell
    !$omp parallel do &
    !$omp default(none) &
    !$omp shared(mesh) &
    !$OMP shared(phi,phi_mass)&
    !$omp schedule(static)
    do i = 1, mesh%nv
      phi%f(i) = phi_mass%f(i)/mesh%hx(i)%areag
    end do
    !$omp end parallel do

    !print*, 'rain', summass, sum(phi_mass%f(:)),summass-sum(phi_mass%f(:))
    !print*,iter
  end subroutine

end module moist_swm
