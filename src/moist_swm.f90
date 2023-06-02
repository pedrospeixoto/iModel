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
    convert_vec_cart2sph, &
    norm, &
    getnearhxedge

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

    !Use routine from highorder module
    use highorder, only: &
    node, &
    controlvolume, &
    method, &
    advmtd, &
    time_integrator, &
    monotonicfilter, &
    moistswm, &
    order, &
    monotonicfilter_rk3, & 
    find_neighbors, &
    allocation, &
    stencil, &
    coordscirc, &
    edges_voronoi, &
    gaussedges, &
    upwind_voronoi, &
    matrix_olg, &
    matrix_gas, &
    condition_initial_gas, &
    gaussrootsweights, &
    reconstruct_velocity_quadrature, &
    edges_indexes, &
    sh_edges_indexes, &
    init_quadrature_edges, &
    divhx, &
    phi_star, phi_min, phi_max, &
    phi_tilda, phi_tilda_min, phi_tilda_max, div_uphi, &
    phi_ed, uphi, &
    F_star, F_cor, F_step2,&
    tendency_advection

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
    type(scalar_field):: qv          !vapour   - diagnostic variable
    type(scalar_field):: hqv         !h*vapour - prognostic variable
    type(scalar_field):: hqv_old     !h*vapour - prognostic variable
    type(scalar_field):: hSqv        !Source for vapour equation
    type(scalar_field):: delta_qv    !Used in source computation 
    type(scalar_field):: qv_exact    !Only for test case 2
    type(scalar_field):: qv_error    !Only for test case 2

    !Cloud state variables  (defined on voronoi centers)
    type(scalar_field):: qc          !cloud   - diagnostic variable
    type(scalar_field):: hqc         !h*cloud - prognostic variable
    type(scalar_field):: hqc_old     !h*cloud - prognostic variable
    type(scalar_field):: hSqc        !Source for cloud equation
    type(scalar_field):: delta_qc    !used in source computation 
    type(scalar_field):: qc_exact    !Only for test case 2
    type(scalar_field):: qc_error    !Only for test case 2

    !Rain state variables  (defined on voronoi centers)
    type(scalar_field):: qr          !rain   - diagnostic variable
    type(scalar_field):: hqr         !h*rain - prognostic variable
    type(scalar_field):: hqr_old     !h*rain - prognostic variable
    type(scalar_field):: hSqr        !Source for rain equation
    type(scalar_field):: delta_qr    !used in source computation 
    type(scalar_field):: qr_exact    !Only for test case 2
    type(scalar_field):: qr_error    !Only for test case 2

    !Tracer state variables  (defined on voronoi centers)
    type(scalar_field):: tracer          !tracer   - diagnostic variable
    type(scalar_field):: htracer         !h*tracer - prognostic variable
    type(scalar_field):: htracer_old     !h*tracer - prognostic variable
    type(scalar_field):: tracer_exact    !Only for test case 2
    type(scalar_field):: tracer_error    !Only for test case 2

    !Source for momentum equation
    type(scalar_field):: Su

    !Scalar fields from hx to ed (defined on edges)
    type(scalar_field):: theta_ed    !temperature

    !Gradients (defined on edges)
    type(scalar_field):: gradPI         !gradient of h^2*temperature
    type(scalar_field):: gradPI_oh      !gradient of h^2*temperature over h
    type(scalar_field):: grad_b         !gradient of topography
    type(scalar_field):: theta_grad_b   !temperature*gradient of topography

    !RHS
    type(scalar_field):: tempeq     !Temperature
    type(scalar_field):: vapoureq   !Vapour
    type(scalar_field):: cloudeq    !Cloud
    type(scalar_field):: raineq     !Rain
    type(scalar_field):: tracereq   !tracer

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

    real(r8), dimension(:), allocatable:: tracerf0
    real(r8), dimension(:), allocatable:: tracerf1
    real(r8), dimension(:), allocatable:: tracerf2
    real(r8), dimension(:), allocatable:: tracerf3

    !Parameters of Zerroukat and Allen 2015 JCP paper
    real(r8)::q0
    real(r8)::Lscale
    real(r8)::gamma_r
    real(r8)::gamma_v
    real(r8)::q_precip
    real(r8)::Twater, iniwater

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

    ! Tracer
    tracereq%n=mesh%nv
    tracereq%name="tracer"
    tracereq%pos=0
    allocate(tracereq%f(1:tracer%n), stat=ist)

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
    qv%n=mesh%nv
    qv%name="vapour"
    qv%pos=0
    allocate(qv%f(1:qv%n), stat=ist)
    allocate(hqv_old%f(1:qv%n), stat=ist)
    allocate(hqv%f(1:qv%n), stat=ist)
    allocate(hSqv%f(1:qv%n), stat=ist)
    allocate(delta_qv%f(1:qv%n), stat=ist)
    allocate(qv_exact%f(1:qv%n), stat=ist)
    allocate(qv_error%f(1:qv%n), stat=ist)

    !Cloud
    qc%n=mesh%nv
    qc%name="cloud"
    qc%pos=0
    allocate(qc%f(1:qc%n), stat=ist)
    allocate(hqc%f(1:qc%n), stat=ist)
    allocate(hqc_old%f(1:qc%n), stat=ist)
    allocate(hSqc%f(1:qc%n), stat=ist)
    allocate(delta_qc%f(1:qc%n), stat=ist)
    allocate(qc_exact%f(1:qc%n), stat=ist)
    allocate(qc_error%f(1:qc%n), stat=ist)

    !Rain
    qr%n=mesh%nv
    qr%name="rain"
    qr%pos=0
    allocate(qr%f(1:qr%n), stat=ist)
    allocate(hqr%f(1:qr%n), stat=ist)
    allocate(hqr_old%f(1:qr%n), stat=ist)
    allocate(hSqr%f(1:qr%n), stat=ist)
    allocate(delta_qr%f(1:qr%n), stat=ist)
    allocate(qr_exact%f(1:qr%n), stat=ist)
    allocate(qr_error%f(1:qr%n), stat=ist)

    !Tracer
    tracer%n=mesh%nv
    tracer%name="tracer"
    tracer%pos=0
    allocate(tracer%f(1:tracer%n), stat=ist)
    allocate(htracer%f(1:tracer%n), stat=ist)
    allocate(htracer_old%f(1:tracer%n), stat=ist)
    allocate(tracer_exact%f(1:tracer%n), stat=ist)
    allocate(tracer_error%f(1:tracer%n), stat=ist)

    !Scalar fields on edges
    theta_ed%n=mesh%ne
    theta_ed%pos=edpos
    allocate(theta_ed%f(1:mesh%ne), stat=ist)

    !Gradient terms
    gradPI%n=mesh%ne
    gradPI%name="gradPI"
    gradPI%pos=edpos
    allocate(gradPI%f(1:gradPI%n), stat=ist)
    allocate(gradPI_oh%f(1:gradPI%n), stat=ist)
    allocate(grad_b%f(1:gradPI%n), stat=ist)
    allocate(theta_grad_b%f(1:gradPI%n), stat=ist)

    !Runge-Kutta variables 
    !Temperature
    allocate(tempf0(1:theta%n))
    allocate(tempf1(1:theta%n))
    allocate(tempf2(1:theta%n))
    allocate(tempf3(1:theta%n))
    
    !Vapour
    allocate(vapourf0(1:qv%n))
    allocate(vapourf1(1:qv%n))
    allocate(vapourf2(1:qv%n))
    allocate(vapourf3(1:qv%n))

    !Cloud
    allocate(cloudf0(1:qc%n))
    allocate(cloudf1(1:qc%n))
    allocate(cloudf2(1:qc%n))
    allocate(cloudf3(1:qc%n))

    !Rain
    allocate(rainf0(1:qr%n))
    allocate(rainf1(1:qr%n))
    allocate(rainf2(1:qr%n))
    allocate(rainf3(1:qr%n))

    !Tracer
    allocate(tracerf0(1:tracer%n))
    allocate(tracerf1(1:tracer%n))
    allocate(tracerf2(1:tracer%n))
    allocate(tracerf3(1:tracer%n))

    !Source for momentum equation
    Su%n=mesh%ne
    Su%name="su"
    Su%pos=edpos
    allocate(Su%f(1:u%n), stat=ist)

    IF(ist>0) STOP 'Error in allocate_globalmoistswmvars'

    call initialize_global_moist_swm_vars()
    if (advmtd=='sg3' .or. advmtd=='og2' .or. advmtd=='og3' .or. advmtd=='og4' .or. advmtd=='upw1' .or. advmtd=='trsk')then
      call highorder_adv_vars()
      if (advmtd .ne. 'upw1')then
        call init_quadrature_edges(mesh)
      end if
    end if
  end subroutine allocate_global_moistswm_vars


!--------------------------------------------------------------------------
! High order advection variables allocation and initialization
!--------------------------------------------------------------------------
subroutine highorder_adv_vars()
    integer(i4) :: i
    integer(i4) :: ngbr
    integer(i4) :: nlines
    integer(i4) :: ncolumns

    !Numero de vizinhos e vizinhos de cada no
    integer(i4),allocatable   :: nbsv(:,:)   

    controlvolume = "V"
    moistswm = .True.

    if (advmtd=='sg3' .or. advmtd=='og2' .or. advmtd=='og3' .or. advmtd=='og4' &
        .or. advmtd=='trsk' .or. advmtd=='upw1') then
        !Alocando nos
        allocate(node(0:mesh%nv))

        !Alocacao basica de memoria
        do i = 1,mesh%nv
          allocate(node(i)%ngbr(1:order-1))
          ngbr = mesh%v(i)%nnb
          allocate(node(i)%ngbr(1)%lvv(1:ngbr+1))
          allocate(node(i)%ngbr(1)%lvd(1:ngbr+1))
          node(i)%ngbr(1)%numberngbr = ngbr

          !Armazena os primeiros vizinhos e suas respectivas distancias
          node(i)%ngbr(1)%lvv(1:ngbr+1) = (/i,mesh%v(i)%nb(1:ngbr)/)
          node(i)%ngbr(1)%lvd(1:ngbr+1) = (/0.0D0,mesh%v(i)%nbdg(1:ngbr)/)
        end do

        if (advmtd=='sg3' .or. advmtd=='og3' .or. advmtd=='og4') then
           nlines=maxval(mesh%v(:)%nnb)
           ncolumns=maxval(mesh%v(:)%nnb)+1  
           allocate(nbsv(13,mesh%nv))
           nbsv = 0
           call find_neighbors(nbsv,nlines,ncolumns,mesh%nv)
           do i=1,mesh%nv
             allocate(node(i)%ngbr(2)%lvv(1:nbsv(1,i)+1))
             allocate(node(i)%ngbr(2)%lvd(1:nbsv(1,i)+1))
             node(i)%ngbr(2)%numberngbr = nbsv(1,i)
             node(i)%ngbr(2)%lvv = (/i,nbsv(2:,i)/)
           end do
           deallocate(nbsv)
       end if
       call allocation(mesh%nv, mesh)

       call stencil(mesh%nv,mesh)

       ! Voronoi's cell - not Donald's cell!
       call coordscirc(mesh%nv,mesh)

       call edges_voronoi(mesh%nv)

       call gaussedges(mesh%nv,mesh)

       call upwind_voronoi(mesh%nv,mesh)

       if (advmtd=='og2' .or. advmtd=='og3' .or. advmtd=='og4')then
          call matrix_olg(mesh%nv,mesh)
       else if(advmtd=='sg3') then!Gassman
          call condition_initial_gas(mesh%nv,mesh) 
          call matrix_gas(mesh%nv,mesh)
       endif

       do i=1,mesh%nv
          node(i)%phi_exa=node(i)%phi_new2
       enddo
    end if

    div_uphi%n = mesh%nv
    uphi%n = mesh%ne
    phi_ed%n = mesh%ne
    phi_ed%pos = 6
    allocate(div_uphi%f(1:div_uphi%n), stat=ist)
    allocate(uphi%f(1:uphi%n), stat=ist)
    allocate(phi_ed%f(1:phi_ed%n), stat=ist)

    if(time_integrator=='rk3' .and. monotonicfilter) then
      phi_star%n = mesh%nv
      phi_tilda%n = mesh%nv
      phi_tilda_min%n = mesh%nv
      phi_tilda_max%n = mesh%nv
      phi_min%n = mesh%nv
      phi_max%n = mesh%nv
      allocate(phi_star%f(1:phi_star%n), stat=ist)
      allocate(phi_tilda%f(1:phi_tilda%n), stat=ist)
      allocate(phi_tilda_min%f(1:phi_tilda_min%n), stat=ist)
      allocate(phi_tilda_max%f(1:phi_tilda_max%n), stat=ist)
      allocate(phi_min%f(1:phi_min%n), stat=ist)
      allocate(phi_max%f(1:phi_max%n), stat=ist)
      allocate(F_star(1:mesh%nv, 1:maxval(mesh%v(:)%nnb)), stat=ist)
      allocate(F_step2(1:mesh%nv, 1:maxval(mesh%v(:)%nnb)), stat=ist)
      allocate(F_cor(1:mesh%nv, 1:maxval(mesh%v(:)%nnb)), stat=ist)
    end if
 
end subroutine highorder_adv_vars

subroutine initialize_global_moist_swm_vars()
    !---------------------------------------------------
    !Initialize fields with zero in paralel
    !---------------------------------------------------

    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(theta, theta_exact, theta_error, htheta, htheta_old, h2theta) &
    !$OMP SHARED(water, hwater) &
    !$OMP SHARED(qv, qv_exact, qv_error, hqv, hqv_old) &
    !$OMP SHARED(qc, qc_exact, qc_error, hqc ,hqc_old) &
    !$OMP SHARED(qr, qr_exact, qr_error, hqr, hqr_old) &
    !$OMP SHARED(hStheta, hSqv, hSqr, hSqc) &
    !$OMP SHARED(delta_qv, delta_qr, delta_qc) &
    !$OMP SHARED(gradPI,gradPI_oh,grad_b,theta_grad_b) &
    !$OMP SHARED(theta_ed) &
    !$OMP SHARED(tempeq, tempf0, tempf1, tempf2, tempf3) &
    !$OMP SHARED(vapoureq, vapourf0, vapourf1, vapourf2, vapourf3) &
    !$OMP SHARED(cloudeq, cloudf0, cloudf1, cloudf2, cloudf3) &
    !$OMP SHARED(raineq, rainf0, rainf1, rainf2, rainf3) &
    !$OMP SHARED(tracer, tracer_exact, tracer_error, htracer, htracer_old) &
    !$OMP SHARED(tracereq, tracerf0, tracerf1, tracerf2, tracerf3) &
    !$OMP SHARED(Su)

    tempeq%f(1:theta%n)=0._r8
    vapoureq%f(1:qv%n)=0._r8
    cloudeq%f(1:qc%n)=0._r8
    raineq%f(1:qr%n)=0._r8
    tracereq%f(:)=0._r8

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
    qv%f=0._r8
    qv_exact=qv
    qv_error=qv
    hqv=qv
    hqv_old=qv
    hSqv=qv
    delta_qv=qv

    !Cloud
    qc%f=0._r8
    qc_exact=qc
    qc_error=qc
    hqc=qc
    hqc_old=qc
    hSqc=qc
    delta_qc=qc

    !Rain
    qr%f=0._r8
    qr_exact=qr
    qr_error=qr
    hqr=qr
    hqr_old=qr
    hSqr=qr
    delta_qr=qr

    !Tracer
    tracer%f=0._r8
    tracer_exact=tracer
    tracer_error=tracer
    htracer=tracer
    htracer_old=tracer

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

    !tracer
    tracerf0(1:qr%n)=0._r8
    tracerf1(1:qr%n)=0._r8
    tracerf2(1:qr%n)=0._r8
    tracerf3(1:qr%n)=0._r8

    !Fields at edges
    theta_ed%f=0._r8

    !Gradients
    gradPI%f=0._r8
    gradPI_oh=gradPI
    grad_b=gradPI
    theta_grad_b=gradPI   
    
    ! Source for momentum equation
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
        qv%f(i) = (1._r8-xsi)*qsat(theta%f(i),h%f(i),bt%f(i),1._r8)

        !Fluxes
        hqv%f(i) = h%f(i)*qv%f(i)
        hqc%f(i) = h%f(i)*qc%f(i)
        hqr%f(i) = h%f(i)*qr%f(i)
        htheta%f(i) = h%f(i)*theta%f(i)
      end do

      q0 = 0.02_r8/maxval(qv%f)
      qv%f = q0*qv%f
      hqv%f = q0*hqv%f

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
          lat = mesh%edhx(l)%c%lat
          lon = mesh%edhx(l)%c%lon
          utmp=u0*dcos(lat)
          vtmp=0._r8
          call convert_vec_sph2cart(utmp, vtmp, mesh%edhx(l)%c%p, vectmp)
          v_ed%p(l)%v=vectmp
          u%f(l)=dot_product(vectmp,mesh%edhx(l)%nr)
        end do
      end if

      h_exact = h
      u_exact = u
      theta_exact = theta
      qv_exact = qv
      qc_exact = qc
      qr_exact = qr

    !======================================================================================
    ! Flow over a mountain - from Zerroukat and Allen JCP 2015
    !======================================================================================
    case(3)

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
        qv%f(i) = mu2*qsat(theta%f(i),h%f(i),bt%f(i),1._r8)
        
        !Fluxes
        hqv%f(i) = h%f(i)*qv%f(i)
        hqc%f(i) = h%f(i)*qc%f(i)
        hqr%f(i) = h%f(i)*qr%f(i)
        htheta%f(i) = h%f(i)*theta%f(i)
      end do

      q0 = 0.02_r8/maxval(qv%f)
      qv%f = q0*qv%f
      hqv%f = q0*hqv%f

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
          lat = mesh%edhx(l)%c%lat
          lon = mesh%edhx(l)%c%lon
          utmp=u0*dcos(lat)
          vtmp=0._r8
          call convert_vec_sph2cart(utmp, vtmp, mesh%edhx(l)%c%p, vectmp)
          v_ed%p(l)%v=vectmp
          u%f(l)=dot_product(vectmp,mesh%edhx(l)%nr)
        end do
      end if

    !======================================================================================
    ! Galewski - Jet in Southern Hemisphere
    !======================================================================================
    case(4,5)
      !Parameters
      !mu1    = 0.05_r8
      mu1    = 0.00002_r8
      mu2    = 0.98_r8
      theta_sp = -40._r8/300._r8 
      theta_eq =  30._r8/300._r8 
      theta_np = -20._r8/300._r8 
      if(testcase==4)then
        u00 = 80.0
        lat0 = pi/7.0
        lat1 = pi/2.0 - lat0
      else
        if(testcase==5)then ! Jet in Southern Hemisphere
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
      nygg = 4*FLOOR(SqrT(REAL(mesh%nv)))
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

      if(testcase == 5)then
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
        if(testcase ==4)then
           l1 = l1 + 120.d0*deg2rad
        end if
        e1 = EXP(-(l1/alpha)**2)
        e2 = EXP(-((lat2 - l2)/beta)**2)
        h%f(i) = h%f(i)+hpert*clat*e1*e2
        !h%f(i) = 1.d0+hpert*clat*e1*e2
        !print*, h%f(i) 
      enddo

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
        qv%f(i) = mu2*qsat(theta%f(i),h%f(i),bt%f(i),1._r8)
        !qv%f(i) = qsat(theta%f(i),h%f(i),bt%f(i),1._r8)

        !Fluxes
        hqv%f(i) = h%f(i)*qv%f(i)
        hqc%f(i) = h%f(i)*qc%f(i)
        hqr%f(i) = h%f(i)*qr%f(i)
        htheta%f(i) = h%f(i)*theta%f(i)
      end do

      q0 = 0.02_r8/maxval(qv%f)
      qv%f = q0*qv%f
      hqv%f = q0*hqv%f

    case default
      print*, "initialize_moist_swm_fields error - please select a proper test case:", testcase
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


  subroutine tendency_moist_swm(h, u, htheta, hqv, hqc, hqr, htracer, &
      masseq, momeq, tempeq, vapoureq, cloudeq, raineq, tracereq, time)
    !--------------------------------------
    !Calculates the Right Hand Side (spatial discret./tendency)
    !   of mass, velocity and moist variables equations
    !-------------------------------------------------------------

    !Fluid thickness (defined on voronoi centers)
    type(scalar_field), intent(in):: h  !General

    !Velocities (defined on edges - only normal component)
    type(scalar_field), intent(inout):: u  !General

    !Temperature (defined on voronoi centers)
    type(scalar_field), intent(in):: htheta  !General

    !Vapour (defined on voronoi centers)
    type(scalar_field), intent(in):: hqv  !General

    !Cloud (defined on voronoi centers)
    type(scalar_field), intent(in):: hqc  !General

    !Rain (defined on voronoi centers)
    type(scalar_field), intent(in):: hqr  !General

    !Tracer (defined on voronoi centers)
    type(scalar_field), intent(in):: htracer  !General

    !Time
    real(r8), intent(in):: time

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

    !Right hand side of tracer (number of cell equations)
    real(r8), intent(inout)::tracereq(:) 

    !Compute the SWM tendency
    call tendency(h, u, masseq, momeq)

    !Initialize RHS of temperature and moist variables equations (in paralel)
    call zero_vector(tempeq)
    call zero_vector(vapoureq)
    call zero_vector(cloudeq)
    call zero_vector(raineq)
    call zero_vector(tracereq)

    !===============================================================
    ! Reconstructs to normal velocity at quadrature points
    !===============================================================
    if (advmtd=='og2' .or. advmtd=='og3' .or. advmtd=='og4' .or. advmtd=='sg3')then
       call reconstruct_velocity_quadrature(mesh, u)
    end if

    !===============================================================
    !Calculate temperature tendency
    !===============================================================
    !Calculate divergence / temp eq RHS
    !call tendency_advection(htheta, tempeq, mesh, time, erad, u)

    !===============================================================
    !Calculate vapour tendency
    !===============================================================
    !Calculate divergence / vapour eq RHS
    !call tendency_advection(hqv, vapoureq, mesh, time, erad, u)

    !===============================================================
    !Calculate cloud tendency
    !===============================================================
    !Calculate divergence / cloud eq RHS
    !call tendency_advection(hqc, cloudeq, mesh, time, erad, u)

    !===============================================================
    !Calculate rain tendency
    !===============================================================
    !Calculate divergence / rain eq RHS
    !call tendency_advection(hqr, raineq, mesh, time, erad, u)

    !===============================================================
    !Calculate tracer tendency
    !===============================================================
    !Calculate divergence / tracer eq RHS
    call tendency_advection(htracer, tracereq, mesh, time, erad, u)

    !===============================================================
    !Compute and add the source
    !===============================================================
    if (time_integrator == 'rk4') then ! RK3 handles with the source in a different way
      call source(h, u, htheta, hqv, hqc, hqr, dt)
      momeq    = momeq    + Su%f 
      tempeq   = tempeq   + hStheta%f
      vapoureq = vapoureq + hSqv%f
      cloudeq  = cloudeq  + hSqc%f
      raineq   = raineq   + hSqr%f  
    end if 
    return

  end subroutine tendency_moist_swm

  subroutine source(h, u, htheta, hqv, hqc, hqr, dtime)
    !Fluid thickness (defined on voronoi centers)
    type(scalar_field), intent(in):: h  !General

    !Velocities (defined on edges - only normal component)
    type(scalar_field), intent(in):: u  !General

    !Temperature (defined on voronoi centers)
    type(scalar_field), intent(in):: htheta  !General

    !Vapour (defined on voronoi centers)
    type(scalar_field), intent(in):: hqv  !General

    !Cloud (defined on voronoi centers)
    type(scalar_field), intent(in):: hqc  !General

    !Rain (defined on voronoi centers)
    type(scalar_field), intent(in):: hqr  !General

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
    call scalar_elem_divide(hqv, h, qv)
    call scalar_elem_divide(hqc, h, qc)
    call scalar_elem_divide(hqr, h, qr)

    !$omp parallel do &
    !$omp default(none) &
    !$omp shared(mesh) &
    !$omp shared(delta_qv, delta_qc, delta_qr) &
    !$OMP shared(hStheta, hSqv, hSqc, hSqr) &
    !$omp shared(qv, qc, qr) &
    !$omp shared(theta, h, bt) &
    !$omp shared(dtime, gamma_v, gamma_r, q0, q_precip, Lscale) &
    !$omp schedule(static)

    do i=1, mesh%nv
      gamma_v = 1._r8+Lscale*dqsat_dtheta(theta%f(i),h%f(i), bt%f(i), q0)
      gamma_v = 1._r8/gamma_v  

      !Vapour
      delta_qv%f(i) = gamma_v*(qv%f(i)-qsat(theta%f(i),h%f(i), bt%f(i), q0))
      delta_qv%f(i) = max(0._r8,delta_qv%f(i))/dtime

      !Cloud
      delta_qc%f(i) = max(0._r8, -gamma_v*(qv%f(i)-qsat(theta%f(i),h%f(i), bt%f(i), q0)))
      delta_qc%f(i) = min(qc%f(i),delta_qc%f(i))/dtime

      !Rain
      delta_qr%f(i) = max(0._r8,gamma_r*(qc%f(i)-q_precip))/dtime

      !Computes the sources
      hSqv%f(i) = delta_qc%f(i) - delta_qv%f(i)
      hSqc%f(i) = delta_qv%f(i) - delta_qc%f(i) - delta_qr%f(i)
      hSqr%f(i) = delta_qr%f(i)
      hStheta%f(i) = Lscale*(delta_qv%f(i)-delta_qc%f(i))

      hSqv%f(i) = h%f(i)*hSqv%f(i)
      hSqc%f(i) = h%f(i)*hSqc%f(i)
      hSqr%f(i) = h%f(i)*hSqr%f(i)
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
    real(r8):: rel_error_qv
    real(r8):: rel_error_qc
    real(r8):: rel_error_qr
    real(r8):: rel_error_theta

    !Total mass of cloud, rain and vapour
    real(r8):: Tcloud
    real(r8):: Train
    real(r8):: Tvapour
    real(r8):: Ttracer
    real(r8):: Ttracer0
    real(r8):: Ttracer_change

    real(r8):: lon, lat, p(1:3), pc(1:3)
    real(r8):: b, c, f, r1, r2, r, lonp, latp, lon1, lat1, lon2, lat2, lon3, lat3


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
    end if

    !Calculate derived initial fields
    call tendency_moist_swm(h, u, htheta, hqv, hqc, hqr, htracer, masseq%f, momeq%f, tempeq%f, vapoureq%f, &
    cloudeq%f, raineq%f, tracereq%f, 0._r8)

    !Plot initial fields
    call plotfields_mswm(0, 0._r8)

    !Calculate total mass
    hwater%f = hqv%f + hqc%f + hqr%f
    inimass=sumf_areas(h)
    iniwater=sumf_areas(hwater)

    !Calculate energies
    call calc_energies(Penergy0, Kenergy0, Tenergy0, Availenergy0)

    u_old=u
    h_old=h
    htheta_old=htheta
    hqv_old=hqv
    hqc_old=hqc
    hqr_old=hqr
    htracer_old=htracer
 
    Ttracer0 = sumf_areas(htracer)

    !Time loop
    do k=1, ntime
        !Calculate u and h for time:
        time=real(k, r8)*dt

        if (time_integrator=='rk4') then
          call ode_rk4_moist_swm(time, h_old, u_old, htheta_old, hqv_old, hqc_old, hqr_old, htracer_old,&
                           h, u, htheta, hqv, hqc, hqr, htracer, dt)
        else if(time_integrator=='rk3') then
          call ode_rk3_moist_swm(time, h_old, u_old, htheta_old, hqv_old, hqc_old, hqr_old, htracer_old,&
                           h, u, htheta, hqv, hqc, hqr, htracer, dt)
        end if

        !Apply the monotonic filter for tracers-only for RK4
        if(time_integrator=='rk4' .and. monotonicfilter)then
          call monotonic_filter(hQv)
          call monotonic_filter(hQr)
          call monotonic_filter(hQc)
        end if

        call scalar_elem_divide(hqv, h, qv)
        call scalar_elem_divide(hqc, h, qc)
        call scalar_elem_divide(hqr, h, qr)

        !compute the mass of each tracer
        Train=sumf_areas(qr)
        Tcloud=sumf_areas(qc)
        Tvapour=sumf_areas(qv)

        if(testcase==2)then
          h_error%f = h_exact%f - h%f
          u_error%f = u_exact%f - u%f
          qv_error%f = qv_exact%f - qv%f
          qc_error%f = qc_exact%f - qc%f
          qr_error%f = qr_exact%f - qr%f
          theta_error%f = theta_exact%f - theta%f
          rel_error_h = maxval(abs(h_error%f))/maxval(abs(h_exact%f))
          rel_error_u = maxval(abs(u_error%f))/maxval(abs(u_exact%f))
          rel_error_qv = maxval(abs(qv_error%f))/maxval(abs(qv_exact%f))
          rel_error_qc = maxval(abs(qc_error%f))
          rel_error_qr = maxval(abs(qr_error%f))
          rel_error_theta = maxval(abs(theta_error%f))/maxval(abs(theta_exact%f))
          !print*, k, ntime
          print*, "Time (dys) :",   k*dt*sec2day, " of ", ntime*dt*sec2day
          print*, "Step = ", k, " of ", ntime
          print '(a33, 3e16.8)','linf errors of (h, u, theta) = ',rel_error_h,rel_error_u,rel_error_theta
          print '(a33, 3e16.8)','linf errors of  (qv, qc, qr) = ',rel_error_qv,rel_error_qc,rel_error_qr
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

        if (rel_error_h>1000.d0 .or. maxval(abs(u%f))>100000._r8) stop
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
        hwater%f = hqr%f + hqv%f + hqc%f
        
        !call scalar_elem_divide(hwater, h, water)
        Twater = sumf_areas(hwater)

        !Calculate erngies
        call calc_energies(Penergy, Kenergy, Tenergy, Availenergy)

        Ttracer = sumf_areas(htracer)
        print '(a33, 2e16.8)','Min/max:', minval(tracer%f),  maxval(tracer%f)
        print '(a33, 2e16.8)','Change in mass of tracer:', (Ttracer-Ttracer0)/Ttracer0
        print '(a33, 2e16.8)','Change in mass of h*(total water):', (Twater-iniwater)/iniwater
        print*,''
        !stop

        !update fields
        u_old=u
        h_old=h
        htheta_old=htheta
        hqv_old=hqv
        hqc_old=hqc
        hqr_old=hqr
        htracer_old=htracer
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

          theta%name=trim(swmname)//"_theta_t"//trim(adjustl(trim(atime)))
          call plot_scalarfield(theta, mesh)

          qv%name=trim(swmname)//"_qv_t"//trim(adjustl(trim(atime)))
          call plot_scalarfield(qv, mesh)

          qc%name=trim(swmname)//"_qc_t"//trim(adjustl(trim(atime)))
          call plot_scalarfield(qc, mesh)

          qr%name=trim(swmname)//"_qr_t"//trim(adjustl(trim(atime)))
          call plot_scalarfield(qr, mesh)

          tracer%name=trim(swmname)//"_tracer_t"//trim(adjustl(trim(atime)))
          call plot_scalarfield(tracer, mesh)

          htracer%name=trim(swmname)//"_htracer_t"//trim(adjustl(trim(atime)))
          call plot_scalarfield(htracer, mesh)

          tracer_error%f = tracer_exact%f - tracer%f
          tracer_error%name=trim(swmname)//"_tracer_error_t"//trim(adjustl(trim(atime)))
          call plot_scalarfield(tracer_error, mesh)

          if(maxval(bt%f(1:bt%n)) > eps)then
            hbt%name=trim(swmname)//"_hbt_t"//trim(adjustl(trim(atime)))
            call plot_scalarfield(hbt, mesh)

            bt%name=trim(swmname)//"_bt_t"//trim(adjustl(trim(atime)))
            call plot_scalarfield(bt, mesh)
          end if
       
          if(testcase==2)then
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

          !q_tr%name=trim(swmname)//"_pv_t"//trim(adjustl(trim(atime)))
          !call plot_scalarfield(q_tr, mesh)

          !q_ed%name=trim(swmname)//"_pv_ed_t"//trim(adjustl(trim(atime)))
          !call plot_scalarfield(q_ed, mesh)

          !divuh%name=trim(swmname)//"_divuh_t"//trim(adjustl(trim(atime)))
          !call plot_scalarfield(divuh, mesh)
      end if
    end subroutine plotfields_mswm



    subroutine ode_rk4_moist_swm ( t, h, u, htheta, hqv, hqc, hqr, htracer, h_new, u_new, htheta_new,&
                                   hqv_new, hqc_new, hqr_new, htracer_new, dt)
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
      type(scalar_field), intent(inout):: u  !General

      !Temperature (defined on voronoi centers) 
      type(scalar_field), intent(in):: htheta  !General

      !Vapour (defined on voronoi centers)
      type(scalar_field), intent(in):: hqv  !General

      !Cloud (defined on voronoi centers)
      type(scalar_field), intent(in):: hqc !General

      !Rain (defined on voronoi centers)
      type(scalar_field), intent(in):: hqr !General

      !Tracer (defined on voronoi centers)
      type(scalar_field), intent(in):: htracer !General

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

      !Tracer (defined on voronoi centers)
      type(scalar_field):: htracer_new  !General
      
      !Times
      real(r8):: t0
      real(r8):: t1
      real(r8):: t2
      real(r8):: t3

      u_new      = u
      h_new      = h
      htheta_new = htheta
      hqv_new    = hqv
      hqc_new    = hqc
      hqr_new    = hqr
      htracer_new= htracer

      masseq%f   = massf0
      momeq%f    = momf0
      tempeq%f   = tempf0
      vapoureq%f = vapourf0
      cloudeq%f  = cloudf0
      raineq%f   = rainf0
      tracereq%f   = tracerf0

      !Initial f (f0)
      t0=t-dt
      call tendency_moist_swm(h, u, htheta, hqv, hqc, hqr, htracer, massf0, momf0, tempf0, vapourf0, cloudf0, rainf0, tracerf0, t0)

      !First RK step
      t1 = t0 + dt/2._r8

      !u_new%f(1:u%n)          = u%f(1:u%n)           + dt * momf0(1:u%n) / 2.0_r8
      !h_new%f(1:h%n)          = h%f(1:h%n)           + dt * massf0(1:h%n) / 2.0_r8
      !htheta_new%f(1:theta%n) = htheta%f(1:theta%n)  + dt * tempf0(1:theta%n) / 2.0_r8
      !hqv_new%f(1:qv%n)       = hqv%f(1:qv%n)        + dt * vapourf0(1:qv%n) / 2.0_r8
      !hqc_new%f(1:qc%n)       = hqc%f(1:qc%n)        + dt * cloudf0(1:qc%n) / 2.0_r8
      !hqr_new%f(1:qr%n)       = hqr%f(1:qr%n)        + dt * rainf0(1:qr%n) / 2.0_r8
      htracer_new%f(1:qr%n)   = htracer%f(1:qr%n)    + dt * tracerf0(1:qr%n) / 2.0_r8

      call tendency_moist_swm(h_new, u_new, htheta_new, hqv_new, hqc_new, hqr_new, htracer_new, &
      massf1, momf1, tempf1, vapourf1, cloudf1, rainf1, tracerf1, t1)

      !Second RK step
      t2 = t0 + dt/2._r8

      !u_new%f(1:u%n)          = u%f(1:u%n)           + dt * momf1(1:u%n) / 2.0_r8
      !h_new%f(1:h%n)          = h%f(1:h%n)           + dt * massf1(1:h%n) / 2.0_r8
      !htheta_new%f(1:theta%n) = htheta%f(1:theta%n)  + dt * tempf1(1:theta%n) / 2.0_r8
      !hqv_new%f(1:qv%n)       = hqv%f(1:qv%n)        + dt * vapourf1(1:qv%n) / 2.0_r8
      !hqc_new%f(1:qc%n)       = hqc%f(1:qc%n)        + dt * cloudf1(1:qc%n) / 2.0_r8
      !hqr_new%f(1:qr%n)       = hqr%f(1:qr%n)        + dt * rainf1(1:qr%n) / 2.0_r8
      htracer_new%f(1:qr%n)   = htracer%f(1:qr%n)    + dt * tracerf1(1:qr%n) / 2.0_r8

      call tendency_moist_swm(h_new, u_new, htheta_new, hqv_new, hqc_new, hqr_new, htracer_new,&
      massf2, momf2, tempf2, vapourf2, cloudf2, rainf2, tracerf2, t2)


      !Third  RK step
      t3 = t0 + dt
      !u_new%f(1:u%n)          = u%f(1:u%n)           + dt * momf2(1:u%n)
      !h_new%f(1:h%n)          = h%f(1:h%n)           + dt * massf2(1:h%n) 
      !htheta_new%f(1:theta%n) = htheta%f(1:theta%n)  + dt * tempf2(1:theta%n)
      !hqv_new%f(1:qv%n)       = hqv%f(1:qv%n)        + dt * vapourf2(1:qv%n) 
      !hqc_new%f(1:qc%n)       = hqc%f(1:qc%n)        + dt * cloudf2(1:qc%n)
      !hqr_new%f(1:qr%n)       = hqr%f(1:qr%n)        + dt * rainf2(1:qr%n)
      htracer_new%f(1:qr%n)   = htracer%f(1:qr%n)    + dt * tracerf2(1:qr%n)

      call tendency_moist_swm(h_new, u_new, htheta_new, hqv_new, hqc_new, hqr_new, htracer_new, &
      massf3, momf3, tempf3, vapourf3, cloudf3, rainf3, tracerf3, t3)

      !
      ! Combine them to estimate the solution at time t+dt
      !
      !u_new%f(1:u%n) = u%f(1:u%n) + dt * (momf0(1:u%n)+2._r8*momf1(1:u%n) &
      !+2._r8*momf2(1:u%n)+momf3(1:u%n))/6._r8

      !h_new%f(1:h%n) = h%f(1:h%n) + dt * (massf0(1:h%n)+2._r8*massf1(1:h%n) &
      !+2._r8*massf2(1:h%n)+massf3(1:h%n))/6._r8

      !htheta_new%f(1:theta%n) = htheta%f(1:theta%n) + dt * (tempf0(1:theta%n)+2._r8*tempf1(1:theta%n) &
      !+2._r8*tempf2(1:theta%n)+tempf3(1:theta%n))/6._r8

      !hqv_new%f(1:qv%n) = hqv%f(1:qv%n) + dt * (vapourf0(1:qv%n)+2._r8*vapourf1(1:qv%n) &
      !+2._r8*vapourf2(1:qv%n)+vapourf3(1:qv%n))/6._r8

      !hqc_new%f(1:qc%n) = hqc%f(1:qc%n) + dt * (cloudf0(1:qc%n)+2._r8*cloudf1(1:qc%n) &
      !+2._r8*cloudf2(1:qc%n)+cloudf3(1:qc%n))/6._r8

      !hqr_new%f(1:qr%n) = hqr%f(1:qr%n) + dt * (rainf0(1:qr%n)+2._r8*rainf1(1:qr%n) &
      !+2._r8*rainf2(1:qr%n)+rainf3(1:qr%n))/6._r8

      htracer_new%f(1:qr%n) = htracer%f(1:qr%n) + dt * (tracerf0(1:qr%n)+2._r8*tracerf1(1:qr%n) &
      +2._r8*tracerf2(1:qr%n)+tracerf3(1:qr%n))/6._r8

      return
    end subroutine ode_rk4_moist_swm


    subroutine ode_rk3_moist_swm ( t, h, u, htheta, hqv, hqc, hqr, htracer, h_new, u_new, htheta_new,&
                                   hqv_new, hqc_new, hqr_new, htracer_new, dt)
      !----------------------------------------------------------------------------------
      !! ode_rk3 takes one Runge-Kutta step for a vector ODE.
      !    t - time that will be calculated (t0+dt)
      !    h - scalar_field for thickness at current time
      !    u - scalar_field for velocities at current time
      !    dt - time step
      !    h_new and u_new - fields at t+dt
      !----------------------------------------------------------------------------------

      !Fluid thickness (defined on voronoi centers)
      type(scalar_field), intent(inout):: h  !General

      !Velocities (defined on edges - only normal component)
      type(scalar_field), intent(inout):: u  !General

      !Temperature (defined on voronoi centers) 
      type(scalar_field), intent(inout):: htheta  !General

      !Vapour (defined on voronoi centers)
      type(scalar_field), intent(inout):: hqv  !General

      !Cloud (defined on voronoi centers)
      type(scalar_field), intent(inout):: hqc !General

      !Rain (defined on voronoi centers)
      type(scalar_field), intent(inout):: hqr !General

      !Tracer (defined on voronoi centers)
      type(scalar_field), intent(inout):: htracer !General

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

      !Tracer (defined on voronoi centers)
      type(scalar_field):: htracer_new  !General
      
      !Times
      real(r8):: t0
      real(r8):: t1
      real(r8):: t2
      real(r8):: t3
      integer(i4) :: i, j, k

      u_new      = u
      h_new      = h
      htheta_new = htheta
      hqv_new    = hqv
      hqc_new    = hqc
      hqr_new    = hqr
      htracer_new= htracer

      masseq%f   = massf0
      momeq%f    = momf0
      tempeq%f   = tempf0
      vapoureq%f = vapourf0
      cloudeq%f  = cloudf0
      raineq%f   = rainf0
      tracereq%f   = tracerf0

      ! Compute source
      call source(h, u, htheta, hqv, hqc, hqr, dt)
 
      call tendency_moist_swm(h, u, htheta, hqv, hqc, hqr, htracer, massf0, momf0, tempf0, vapourf0, cloudf0, rainf0, tracerf0, t)

      !First RK step
      u_new%f(1:u%n)          = u%f(1:u%n)           + dt * (momf0(1:u%n) + Su%f(1:u%n)) / 3.0_r8
      h_new%f(1:h%n)          = h%f(1:h%n)           + dt *  massf0(1:h%n) / 3.0_r8
      htheta_new%f(1:theta%n) = htheta%f(1:theta%n)  + dt * (tempf0(1:theta%n) + hStheta%f(1:theta%n)) / 3.0_r8
      hqv_new%f(1:qv%n)       = hqv%f(1:qv%n)        + dt * (vapourf0(1:qv%n)  + hSqv%f(1:qv%n)      ) / 3.0_r8
      hqc_new%f(1:qc%n)       = hqc%f(1:qc%n)        + dt * (cloudf0(1:qc%n)   + hSqc%f(1:qc%n)      ) / 3.0_r8
      hqr_new%f(1:qr%n)       = hqr%f(1:qr%n)        + dt * (rainf0(1:qr%n)    + hSqr%f(1:qr%n)      ) / 3.0_r8
      htracer_new%f(1:qr%n)   = htracer%f(1:qr%n)    + dt * tracerf0(1:qr%n) / 3.0_r8

      call tendency_moist_swm(h_new, u_new, htheta_new, hqv_new, hqc_new, hqr_new, htracer_new, &
      massf1, momf1, tempf1, vapourf1, cloudf1, rainf1, tracerf1, t+dt/3._r8)

      !Second RK step
      u_new%f(1:u%n)          = u%f(1:u%n)           + dt * (momf1(1:u%n) + Su%f(1:u%n)) / 2.0_r8
      h_new%f(1:h%n)          = h%f(1:h%n)           + dt * massf1(1:h%n) / 2.0_r8
      htheta_new%f(1:theta%n) = htheta%f(1:theta%n)  + dt * (tempf1(1:theta%n) + hStheta%f(1:theta%n)) / 2.0_r8
      hqv_new%f(1:qv%n)       = hqv%f(1:qv%n)        + dt * (vapourf1(1:qv%n)  + hSqv%f(1:qv%n)      ) / 2.0_r8
      hqc_new%f(1:qc%n)       = hqc%f(1:qc%n)        + dt * (cloudf1(1:qc%n)   + hSqc%f(1:qc%n)      ) / 2.0_r8
      hqr_new%f(1:qr%n)       = hqr%f(1:qr%n)        + dt * (rainf1(1:qr%n)    + hSqr%f(1:qr%n)      ) / 2.0_r8
      htracer_new%f(1:qr%n)   = htracer%f(1:qr%n)    + dt * tracerf1(1:qr%n) / 2.0_r8

      call tendency_moist_swm(h_new, u_new, htheta_new, hqv_new, hqc_new, hqr_new, htracer_new,&
      massf2, momf2, tempf2, vapourf2, cloudf2, rainf2, tracerf2, t+dt/2._r8)

      ! Third  RK step
      ! Last RK step applies a different approach if the monotonic limiter is active 
      if(.not. monotonicfilter) then
        u_new%f(1:u%n)          = u%f(1:u%n)           + dt * (momf2(1:u%n) + Su%f(1:u%n))
        h_new%f(1:h%n)          = h%f(1:h%n)           + dt * massf2(1:h%n)
        htheta_new%f(1:theta%n) = htheta%f(1:theta%n)  + dt * (tempf2(1:theta%n) + hStheta%f(1:theta%n))
        hqv_new%f(1:qv%n)       = hqv%f(1:qv%n)        + dt * (vapourf2(1:qv%n)  + hSqv%f(1:qv%n))
        hqc_new%f(1:qc%n)       = hqc%f(1:qc%n)        + dt * (cloudf2(1:qc%n)   + hSqc%f(1:qc%n))
        hqr_new%f(1:qr%n)       = hqr%f(1:qr%n)        + dt * (rainf2(1:qr%n)    + hSqr%f(1:qr%n))
        htracer_new%f(1:qr%n)   = htracer%f(1:qr%n)    + dt * tracerf2(1:qr%n)

      else if(monotonicfilter) then
        ! Update tracers
        call  monotonicfilter_rk3(mesh, htheta , htheta_new , dt, erad, t, u, u_new, hStheta)
        call  monotonicfilter_rk3(mesh, hqv    , hqv_new    , dt, erad, t, u, u_new, hSqv)
        call  monotonicfilter_rk3(mesh, hqr    , hqr_new    , dt, erad, t, u, u_new, hSqr)
        call  monotonicfilter_rk3(mesh, hqc    , hqc_new    , dt, erad, t, u, u_new, hSqc)
        call  monotonicfilter_rk3(mesh, htracer, htracer_new, dt, erad, t, u, u_new)
 
        !Update momentum and mass conservation equations
        u_new%f(1:u%n)          = u%f(1:u%n)           + dt * (momf2(1:u%n) + Su%f(1:u%n))
        h_new%f(1:h%n)          = h%f(1:h%n)           + dt * massf2(1:h%n)
      end if
      return
    end subroutine ode_rk3_moist_swm

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

    !Monotonization flag
    character (len=6)::mono
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
    read(fileunit,*)  advmtd
    read(fileunit,*)  buffer
    read(fileunit,*)  time_integrator
    read(fileunit,*)  buffer
    read(fileunit,*)  mono
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
    print*, "Advection method        : ", advmtd
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
   
    ! Advection scheme order
    select case(advmtd)
    case('upw1')
        order=2 ! just to avoid allocation problems in highorder_adv_vars
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

    case('trsk')
        order=2 ! just to avoid allocation problems in highorder_adv_vars
        method='O'


    case default
        print*, "Invalid advection method. Please select a proper order."
        stop
    end select

    swmname=trim(swmname)//"_advmethod_"//trim(advmtd)

    if (time_integrator == 'rk3' .or. time_integrator == 'rk4') then
        swmname=trim(swmname)//"_"//trim(time_integrator)
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
    swmname=trim(swmname)//"_mono"//trim(mono)

    print*, "SWM Name for Plots: ", trim(swmname)

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
    real(r8):: mass_hqc, mass_hqv, mass_hqr, modified_mass_hqr, sumareas, summass,summass2
    real(r8):: eps, eps2

    phi_mass= phi
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
