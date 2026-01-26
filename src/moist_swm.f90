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
    nlon_alt, &
    eps2

    !Use main grid data structures
    use datastruct, only: &
    grid_structure, &
    scalar_field, &
    gaussian_quadrature

    !Use routines from the spherical mesh pack
    use smeshpack, only: &
    convert_vec_sph2cart,&
    arcdistll, &
    sph2cart, &
    cart2sph, &
    convert_vec_cart2sph, &
    norm, &
    getnearhxedge, &
    modint

    !Use routines from the interpolation pack
    use interpack, only: &
    plot_scalarfield, &
    plot_cart_vectorfield, &
    aplyr, &
    aplyrt, &
    constr, &
    vecrecon_lsqeval, &
    vecrecon_lsqfitpol
 

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
    node, order, controlvolume, &
    find_neighbors, &
    allocation, &
    stencil, &
    coordscirc, &
    edges_voronoi, &
    gaussedges, &
    upwind_voronoi, &
    matrix_og, &
    init_quadrature_edges, &
    pseudoinversa

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

    !Tracer state variables  (defined on voronoi centers)
    type(scalar_field):: tracer          !tracer   - diagnostic variable
    type(scalar_field):: htracer         !h*tracer - prognostic variable
    type(scalar_field):: htracer_old     !h*tracer - prognostic variable
    type(scalar_field):: tracer_exact    !Only for test case 2
    type(scalar_field):: tracer_error    !Only for test case 2

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
    type(scalar_field):: htracer_ed  !h*tracer

    !velocity x scalar field (defined on edges - only normal component)
    type(scalar_field):: uhtheta     !velocity*h*temperature
    type(scalar_field):: uhQv        !velocity*h*vapour
    type(scalar_field):: uhQc        !velocity*h*cloud
    type(scalar_field):: uhQr        !velocity*h*rain
    type(scalar_field):: uhtracer    !velocity*h*tracer

    !Divergences (defined on voronoi centers)
    type(scalar_field):: div_uhtheta  !div of velocity*h*temperature
    type(scalar_field):: div_uhQv     !div of velocity*h*vapour
    type(scalar_field):: div_uhQc     !div of velocity*h*cloud
    type(scalar_field):: div_uhQr     !div of velocity*h*rain
    type(scalar_field):: div_uhtracer !div of velocity*h*tracer

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
    type(scalar_field):: tracereq   !tracer

    ! Variables for transport of a generic tracer
    type(scalar_field):: phi      ! tracer
    type(scalar_field):: uphi      ! tracer flux at edges
    type(scalar_field):: phi_ed    ! tracer at edges
    type(scalar_field):: div_uphi ! divergence

    ! variables for flux corrected transport scheme (flux_limiter)
    type(scalar_field):: phi_lo_upw   ! tracer computed using 1st order upwind
    type(scalar_field):: phi_min   !
    type(scalar_field):: phi_max   ! 
    type(scalar_field):: phi_plus  ! 
    type(scalar_field):: phi_minus ! 
    type(scalar_field):: uphi_lo_upw  ! tracer flux (1st order upwind) at edges
    type(scalar_field):: uphi_cor  ! corrected flux
    type(scalar_field):: r_limiter ! limiting factor

    ! Quadrature - used by OG schemes
    type(gaussian_quadrature) :: gauss_quad

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

    ! Monotonic filter
    integer(i4):: mono_filter

    ! Advection method for tracers
    character (len=6):: advmtd

    ! Time integrator 
    character (len=32):: time_integrator
    
    ! wind reconstruction
    character (len=5):: urecon_mtd = "ed1" ! edge centered linear LSQ reconstruction
    !character (len=5):: urecon_mtd = "ed2" ! edge centered quadratic LSQ reconstruction
    !character (len=5):: urecon_mtd = "hx" ! hexagon centered LSQ reconstruction


!=============================================================================
contains

subroutine allocate_global_moistswm_vars()
    !------------------------------------
    !Allocate fields for the moist SWM
    !------------------------------------

    !Error flag
    integer(i4):: ist
    integer(i4):: e

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

    !Tracer
    tracer%n=mesh%nv
    tracer%name="tracer"
    tracer%pos=0
    allocate(tracer%f(1:tracer%n), stat=ist)
    allocate(htracer%f(1:tracer%n), stat=ist)
    allocate(htracer_old%f(1:tracer%n), stat=ist)
    allocate(tracer_exact%f(1:tracer%n), stat=ist)
    allocate(tracer_error%f(1:tracer%n), stat=ist)

    !Divergence terms
    div_uhQv%n=mesh%nv
    div_uhQv%pos=0
    div_uhQv%name="divuqv"
    allocate(div_uhQv%f(1:div_uhQv%n), stat=ist)
    allocate(div_uhQc%f(1:div_uhQv%n), stat=ist)
    allocate(div_uhQr%f(1:div_uhQv%n), stat=ist)
    allocate(div_uhtheta%f(1:div_uhQv%n), stat=ist)
    allocate(div_uhtracer%f(1:div_uhQv%n), stat=ist)

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
    allocate(htracer_ed%f(1:hQv_ed%n), stat=ist)
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
    allocate(uhtracer%f(1:uhQv%n), stat=ist)

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

    !Tracer
    allocate(tracerf0(1:tracer%n))
    allocate(tracerf1(1:tracer%n))
    allocate(tracerf2(1:tracer%n))
    allocate(tracerf3(1:tracer%n))

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

    !----------------------------------------------------------
    ! Variables for transport of a tracer
    phi%n = mesh%nv
    div_uphi%n = mesh%nv
    uphi%n = mesh%ne
    phi_ed%n = mesh%ne

    phi%pos = 0
    div_uphi%pos = 0
    uphi%pos = 6
    phi_ed%pos = 6

    allocate(phi%f(1:phi%n), stat=ist)
    allocate(phi_ed%f(1:phi_ed%n), stat=ist)
    allocate(div_uphi%f(1:div_uphi%n), stat=ist)
    allocate(uphi%f(1:uphi%n), stat=ist)
    !----------------------------------------------------------


    !----------------------------------------------------------
    ! First order upwind scheme variables
    ! Needed by flux_limiter routine
    if (mono_filter==2) then
       phi_lo_upw%n = mesh%nv
       phi_min%n = mesh%nv
       phi_max%n = mesh%nv
       phi_plus%n = mesh%nv
       phi_minus%n = mesh%nv
       uphi_lo_upw%n = mesh%ne
       uphi_cor%n = mesh%ne
       r_limiter%n = mesh%ne

       phi_lo_upw%pos = 0
       phi_min%pos = 0
       phi_max%pos = 0
       phi_plus%pos = 0
       phi_minus%pos = 0
       uphi_lo_upw%pos = 6
       uphi_cor%pos = 6
       r_limiter%n = mesh%ne

       allocate(phi_lo_upw%f(1:phi_lo_upw%n), stat=ist)
       allocate(phi_min%f(1:phi_min%n), stat=ist)
       allocate(phi_max%f(1:phi_max%n), stat=ist)
       allocate(phi_plus %f(1:phi_plus%n), stat=ist)
       allocate(phi_minus%f(1:phi_minus%n), stat=ist)
       allocate(uphi_lo_upw%f(1:uphi_lo_upw%n), stat=ist)
       allocate(uphi_cor%f(1:uphi_cor%n), stat=ist)
       allocate(r_limiter%f(1:r_limiter%n), stat=ist)
    endif
    !----------------------------------------------------------

    IF(ist>0) STOP 'Error in allocate_globalmoistswmvars'

    call initialize_global_moist_swm_vars()


    !----------------------------------------------------------------------------------------------
    call map_edge_to_cell_local_index(mesh)  

    ! Least squares parameters for sg3/sg4
    if (advmtd=='sg3' .or. advmtd=='sg4') then
       ! 2nd Order polinomial least square fit on node
       allocate(phi%pol(1:mesh%nv))
       do i = 1, mesh%nv
          allocate(phi%pol(i)%c(1:6))
          allocate(phi%pol(i)%tgplane_nr(1:mesh%v(i)%nnb))
          allocate(phi%pol(i)%tgplane_coords(1:mesh%v(i)%nnb))
          allocate(phi%pol(i)%n2d(1:mesh%v(i)%nnb))
          allocate(phi%pol(i)%rhs(1:mesh%v(i)%nnb))
          allocate(phi%pol(i)%lsq_matrix(1:mesh%v(i)%nnb,5))
          allocate(phi%pol(i)%lsq_matrix_pinv(5,1:mesh%v(i)%nnb))
       end do
       call init_lsq_rotation_sg(phi, mesh)

    else if (advmtd=='og2' .or. advmtd=='og3' .or. advmtd=='og4') then
        call highorder_adv_vars()
        allocate(phi%pol(1:mesh%nv))

        if (advmtd=='og2') then !1st order poly + 1 quadrature point
           do i = 1, mesh%nv
              allocate(phi%pol(i)%c(1:3))
           end do
           nquad = 1

        else if(advmtd=='og3') then !2nd order poly + 2 quadrature points
           do i = 1, mesh%nv
              allocate(phi%pol(i)%c(1:6))
           end do
           nquad = 2

        else if(advmtd=='og4') then !3rd order poly + 2 quadrature points
           do i = 1, mesh%nv
              allocate(phi%pol(i)%c(1:10))
           end do
           nquad = 2

        end if

        ! Gaussian quadrature allocation
        ne = mesh%ne
        gauss_quad%n = nquad
        allocate(gauss_quad%edge(1:ne))
        do e = 1, mesh%ne
           allocate(gauss_quad%edge(e)%node(1:nquad))
           allocate(gauss_quad%edge(e)%tgp_p1(1:nquad))
           allocate(gauss_quad%edge(e)%tgp_p2(1:nquad))
           allocate(gauss_quad%edge(e)%nr(1:nquad))
           allocate(gauss_quad%edge(e)%q(1:nquad))
           allocate(gauss_quad%edge(e)%u(1:nquad))
           allocate(gauss_quad%edge(e)%w(1:nquad))
        enddo 
        call gaussedges2(mesh, gauss_quad%n)
        call init_lsq_rotation_og(phi, u, mesh)
    endif


end subroutine allocate_global_moistswm_vars

!--------------------------------------------------------------------------
! High order advection variables allocation and initialization
!--------------------------------------------------------------------------
subroutine highorder_adv_vars()
    integer(i4) :: i
    integer(i4) :: ngbr
    integer(i4) :: nlines
    integer(i4) :: ncolumns
    integer(i4) :: nodes

    !Numero de vizinhos e vizinhos de cada no
    integer(i4),allocatable   :: nbsv(:,:)   

    controlvolume = "V"
    !moistswm = .True.
    !print*, nodes

    if (advmtd=='og2' .or. advmtd=='og3' .or. advmtd=='og4') then
        if (advmtd=='og2') then
            order=2
        else if (advmtd=='og3') then
            order=3
        else
            order=4
        endif

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
           allocate(nbsv(15,nodes))
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
 
        call coordscirc(nodes,mesh)

        call edges_voronoi(nodes)

        call gaussedges(nodes)

        call upwind_voronoi(nodes,mesh)

        call matrix_og(nodes,mesh)

        call init_quadrature_edges(mesh)

       !if (advmtd=='og2' .or. advmtd=='og3' .or. advmtd=='og4')then
       !   call init_lsq_rotation_og(phi, u, mesh)
       !endif
 
    end if
end subroutine highorder_adv_vars

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
    !$OMP SHARED(div_uhQv,div_uhQc,div_uhQr,div_uhtheta,div_uhtracer) &
    !$OMP SHARED(gradPI,gradPI_oh,grad_b,theta_grad_b) &
    !$OMP SHARED(div_uhQv_exact,div_uhQc_exact,div_uhQr_exact,div_uhtheta_exact) &
    !$OMP SHARED(gradPI_exact,gradPI_oh_exact,grad_b_exact,theta_grad_b_exact) &
    !$OMP SHARED(uhQv,uhQc,uhQr,uhtheta,uhtracer) &
    !$OMP SHARED(hqv_ed,hqc_ed,hqr_ed,htheta_ed,htracer_ed,theta_ed) &
    !$OMP SHARED(tempeq, tempf0, tempf1, tempf2, tempf3) &
    !$OMP SHARED(vapoureq, vapourf0, vapourf1, vapourf2, vapourf3) &
    !$OMP SHARED(cloudeq, cloudf0, cloudf1, cloudf2, cloudf3) &
    !$OMP SHARED(raineq, rainf0, rainf1, rainf2, rainf3) &
    !$OMP SHARED(tracer, tracer_exact, tracer_error, htracer, htracer_old) &
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

    !Divergences
    div_uhQv%f=0._r8
    div_uhQc=div_uhQv
    div_uhQr=div_uhQv
    div_uhtheta=div_uhQv
    div_uhtracer=div_uhQv

    div_uhQv_exact=div_uhQv
    div_uhQc_exact=div_uhQv
    div_uhQr_exact=div_uhQv
    div_uhtheta_exact=div_uhQv

    !Fields at edges
    hQv_ed%f=0._r8
    hQc_ed=hQv_ed
    hQr_ed=hQv_ed
    htheta_ed=hQv_ed
    htracer_ed=hQv_ed
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
    uhtracer=uhQv

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
    real(r8)::p(1:3), pc(1:3)
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
    ! Pure advection of a Gaussian hill using a zonal wind
    !======================================================================================
    case(1)
        ! Fields at centers
        do i=1, mesh%nv
            lon = mesh%v(i)%lon
            lat = mesh%v(i)%lat
            h%f(i) = 1._r8
            bt%f(i) = 0._r8 
            theta%f(i) = 0._r8
            qv%f(i) = 0._r8
            hqv%f(i) = h%f(i)*qv%f(i)
            hqc%f(i) = h%f(i)*qc%f(i)
            hqr%f(i) = h%f(i)*qr%f(i)
            htheta%f(i) = h%f(i)*theta%f(i)

            ! tracer
            call sph2cart(lon, lat, p(1), p(2), p(3))
            call sph2cart(0._r8, 0._r8, pc(1), pc(2), pc(3))
            tracer%f(i) = f(lon, lat)
            htracer%f(i) = tracer%f(i)*h%f(i)
            tracer_exact%f(i) = tracer%f(i)
        end do

        !Velocity
        u0 = 2._r8*pi*erad/(12._r8*day2sec)
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


    !======================================================================================
    ! Steady state - from Zerroukat and Allen JCP 2015
    !======================================================================================
    case(2)
      !Field at cell's center
      !u0    = 20._r8
      u0 = 2._r8*pi*erad/(12._r8*day2sec)
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

        ! tracer
        call sph2cart(lon, lat, p(1), p(2), p(3))
        call sph2cart(0._r8, 0._r8, pc(1), pc(2), pc(3))
        tracer%f(i) = 1_r8*dexp(-5._r8*norm(p-pc)**2)
        tracer%f(i) = f(lon, lat)
        htracer%f(i) = tracer%f(i)*h%f(i)
        tracer_exact%f(i) = tracer%f(i)
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
          lat = mesh%edhx(l)%c%lat
          lon = mesh%edhx(l)%c%lon
          utmp=u0*cos(lat)
          vtmp=0._r8
          call convert_vec_sph2cart(utmp, vtmp, mesh%edhx(l)%c%p, vectmp)
          v_ed%p(l)%v=vectmp
          u%f(l)=dot_product(vectmp,mesh%edhx(l)%nr)
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
    !stop
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
    type(scalar_field), intent(inout):: htheta  !General

    !Vapour (defined on voronoi centers)
    type(scalar_field), intent(inout):: hQv  !General

    !Cloud (defined on voronoi centers)
    type(scalar_field), intent(inout):: hQc  !General

    !Rain (defined on voronoi centers)
    type(scalar_field), intent(inout):: hQr  !General

    !Tracer (defined on voronoi centers)
    type(scalar_field), intent(inout):: htracer  !General

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

    real(r8) :: t0

    t0 = time
    !Compute the SWM tendency
    if(testcase>1) then
       call tendency(h, u, masseq, momeq)
    endif
    !===============================================================
    ! Reconstructs to normal velocity at quadrature points
    !===============================================================
    if (advmtd=='og2' .or. advmtd=='og3' .or. advmtd=='og4') then
        call reconstruct_velocity_quadrature(mesh, u)
    end if

    if(testcase>1) then
       !Initialize RHS of temperature and moist variables equations (in paralel)
       call zero_vector(tempeq)
       call zero_vector(vapoureq)
       call zero_vector(cloudeq)
       call zero_vector(raineq)

       if(testcase==2) call zero_vector(tracereq)
       !===============================================================
       !Calculate temperature tendency
       !===============================================================
       call tendency_advection(htheta, tempeq, u, uhtheta, mesh)

       !===============================================================
       !Calculate vapour tendency
       !===============================================================
       call tendency_advection(hqv, vapoureq, u, uhQv, mesh)

       !===============================================================
       !Calculate cloud tendency
       !===============================================================
       call tendency_advection(hqc, cloudeq, u, uhQc, mesh)

       !===============================================================
       !Calculate rain tendency
       !===============================================================
       call tendency_advection(hqr, raineq, u, uhQr, mesh)

       !===============================================================
       !Calculate tracer tendency for TC2
       !===============================================================
       if (testcase==2) call tendency_advection(htracer, tracereq, u, uhtracer, mesh)

       !===============================================================
       !Compute and add the source
       !===============================================================

       call source(h, htheta, hQv, hQc, hQr, dt)

       momeq    = momeq    + Su%f 
       tempeq   = tempeq   + hStheta%f
       vapoureq = vapoureq + hSQv%f
       cloudeq  = cloudeq  + hSQc%f
       raineq   = raineq   + hSQr%f  


    else if(testcase==1)then
        !===============================================================
        !Calculate tracer tendency
        !===============================================================
        call tendency_advection(htracer, tracereq, u, uhtracer, mesh)
    end if


    return
  end subroutine tendency_moist_swm

  subroutine source(h, htheta, hQv, hQc, hQr, dtime)
    !Fluid thickness (defined on voronoi centers)
    type(scalar_field), intent(in):: h  !General

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

  subroutine source_u(h, htheta)
    !Fluid thickness (defined on voronoi centers)
    type(scalar_field), intent(in):: h  !General

    !Temperature (defined on voronoi centers)
    type(scalar_field), intent(in):: htheta  !General


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

    call scalar_elem_divide(htheta, h, theta)
  end subroutine source_u


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
    real(r8):: Ttracer
    real(r8):: Ttracer0

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

    call tendency_moist_swm(h, u, htheta, hQv, hQc, hQr, htracer, &
    masseq%f, momeq%f, tempeq%f, vapoureq%f, cloudeq%f, raineq%f, tracereq%f, 0._r8)

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


      !Apply the monotonic filter for tracers (negative mass distribution)
      if(mono_filter==1) then
        if(testcase>1)then 
          call monotonic_filter(hQv)             
          call monotonic_filter(hQr)
          call monotonic_filter(hQc)
          if(testcase==2) call monotonic_filter(htracer)
        else
          call monotonic_filter(htracer)
        endif
      endif

      !call mass_fixer(hQv,hQr,hQc, iniwater)
      !call scalar_elem_divide(htheta, h, theta)
      call scalar_elem_divide(hQv, h, Qv)
      call scalar_elem_divide(hQc, h, Qc)
      call scalar_elem_divide(hQr, h, Qr)

      !compute the mass of each tracer
      Train=sumf_areas(qr)
      Tcloud=sumf_areas(qc)
      Tvapour=sumf_areas(qv)

      if(testcase==2) then
        Ttracer = sumf_areas(htracer)
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
        print '(a33, 2e16.8)','tracer Min/max = ', minval(htracer%f/h%f),  maxval(htracer%f/h%f)
        print '(a33, 2e16.8)','Change in mass of tracer:', (Ttracer-Ttracer0)/Ttracer0 
        !print '(a22, 2e16.8)',' height = ',minval(h%f),maxval(h%f)
        !print '(a22, 2e16.8)',' velocity = ',minval(u%f),maxval(u%f)
        !print '(a22, 2e16.8)',' temperature = ',minval(theta%f),maxval(theta%f)
        !print '(a22, 2e16.8)',' vapour = ',minval(qv%f),maxval(qv%f)
        !print '(a22, 2e16.8)',' cloud = ',minval(qc%f),maxval(qc%f)
        !print '(a22, 2e16.8)',' rain = ',minval(qr%f),maxval(qr%f)
        call write_swmp_error_file(time, k)
      else 
        print*, "Time (dys) :",   k*dt*sec2day, " of ", ntime*dt*sec2day
        print*, "Step = ", k, " of ", ntime
        print*,'                          min               max               mass'
        if (testcase==1)then
            Ttracer = sumf_areas(htracer)
            print '(a33, 2e16.8)','tracer Min/max = ', minval(htracer%f/h%f),  maxval(htracer%f/h%f)
            print '(a33, 2e16.8)','Change in mass of tracer:', (Ttracer-Ttracer0)/Ttracer0
     
        else
            print '(a22, 2e16.8)',' height = ',minval(h%f),maxval(h%f)
            print '(a22, 2e16.8)',' velocity = ',minval(u%f),maxval(u%f)
            print '(a22, 2e16.8)',' temperature = ',minval(theta%f),maxval(theta%f)
            print '(a22, 3e16.8)',' vapour = ',minval(qv%f),maxval(qv%f), Tvapour
            print '(a22, 3e16.8)',' cloud = ',minval(qc%f),maxval(qc%f), Tcloud
            print '(a22, 3e16.8)',' rain = ',minval(qr%f),maxval(qr%f), Train
        end if

      end if

      !print*,'CFL = ',cfl
      !Plot fields
      call plotfields_mswm(k, time)

      !Write errors in file
      call write_evol_file_mswm(k, time, iniwater, Twater, inimass, Penergy0, Kenergy0, Tenergy0, Availenergy0,&
                                   tmass, Penergy, Kenergy, Tenergy, Availenergy)

      call write_water_evol_file(k, time, Train, Tcloud, Tvapour)
      
      !Calculate total mass
      Tmass=sumf_areas(h)

      !Calculate total water
      hwater%f = hQr%f + hQv%f + hQc%f
      
      !call scalar_elem_divide(hwater, h, water)
      Twater = sumf_areas(hwater)

      !Calculate erngies
      call calc_energies(Penergy, Kenergy, Tenergy, Availenergy)
      if(testcase>=2) then
         print '(a33, 2e16.8)','Change in mass of h*(total water):', (Twater-iniwater)/iniwater
         print*,''
      end if
      print*, 'adv=', advmtd, mesh%nv

      !update fields
      u_old=u
      h_old=h
      htheta_old=htheta
      hQv_old=hQv
      hQc_old=hQc
      hQr_old=hQr
      htracer_old=htracer
    end do
    if(testcase<2) then
       print '(a33, 2e16.8)','Final error:', maxval(abs(htracer%f-tracer_exact%f))
       print*,''
    end if
 
  end  subroutine moist_swm_tests


  subroutine write_evol_file_mswm(k, time, iniwater, Twater, inimass, Penergy0, Kenergy0, Tenergy0, Availenergy0,&
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
  integer(i4), intent(in) :: k
  integer(i4) :: l
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
end subroutine write_evol_file_mswm

subroutine write_water_evol_file(k,time, train, tcloud, tvapour)
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
integer(i4), intent(in) :: k

integer(i4) :: l
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
          if(testcase>1)then
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

              if(maxval(bt%f(1:bt%n)) > eps)then
                hbt%name=trim(swmname)//"_hbt_t"//trim(adjustl(trim(atime)))
                call plot_scalarfield(hbt, mesh)

                bt%name=trim(swmname)//"_bt_t"//trim(adjustl(trim(atime)))
                call plot_scalarfield(bt, mesh)
              end if

              if(testcase==2)then
                  tracer%f = htracer%f/h%f
                  tracer%name=trim(swmname)//"_tracer_t"//trim(adjustl(trim(atime)))
                  call plot_scalarfield(tracer, mesh)

                  tracer_error%f = tracer_exact%f - tracer%f
                  tracer_error%name=trim(swmname)//"_tracer_error_t"//trim(adjustl(trim(atime)))
                  call plot_scalarfield(tracer_error, mesh)
              endif
          else
            tracer%f = htracer%f/h%f
            tracer%name=trim(swmname)//"_tracer_t"//trim(adjustl(trim(atime)))
            call plot_scalarfield(tracer, mesh)

            tracer_error%f = tracer_exact%f - tracer%f
            tracer_error%name=trim(swmname)//"_tracer_error_t"//trim(adjustl(trim(atime)))
            call plot_scalarfield(tracer_error, mesh)
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
    type(scalar_field), intent(inout):: htheta  !General

    !Vapour (defined on voronoi centers)
    type(scalar_field), intent(inout):: hQv  !General

    !Cloud (defined on voronoi centers)
    type(scalar_field), intent(inout):: hQc !General

    !Rain (defined on voronoi centers)
    type(scalar_field), intent(inout):: hQr !General

    !Tracer (defined on voronoi centers)
    type(scalar_field), intent(inout):: htracer !General

    !Time and time step
    real(r8):: t
    real(r8):: dt

    !Updated fields
    !Fluid thickness (defined on voronoi centers)
    type(scalar_field):: h_new  !General

    !Velocities (defined on edges - only normal component)
    type(scalar_field), intent(inout):: u_new  !General

    !Temperature (defined on voronoi centers)
    type(scalar_field), intent(inout):: htheta_new  !General

    !Vapour (defined on voronoi centers)
    type(scalar_field), intent(inout):: hqv_new  !General

    !Cloud (defined on voronoi centers)
    type(scalar_field), intent(inout):: hqc_new  !General

    !Rain (defined on voronoi centers)
    type(scalar_field), intent(inout):: hqr_new  !General

    !Tracer (defined on voronoi centers)
    type(scalar_field), intent(inout):: htracer_new  !General
 
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
    tracereq%f = tracerf0

    !Initial f (f0)
    t0=t-dt
    call tendency_moist_swm(h, u, htheta, hqv, hqc, hqr, htracer, massf0, momf0, tempf0, vapourf0, cloudf0, rainf0, tracerf0, t0)

    !First RK step
    t1 = t0 + dt/2._r8
 
    if(testcase>1)then
        u_new%f(1:u%n)          = u%f(1:u%n)           + dt * momf0(1:u%n) / 2.0_r8
        h_new%f(1:h%n)          = h%f(1:h%n)           + dt * massf0(1:h%n) / 2.0_r8
        htheta_new%f(1:theta%n) = htheta%f(1:theta%n)  + dt * tempf0(1:theta%n) / 2.0_r8
        hqv_new%f(1:qv%n)       = hqv%f(1:qv%n)        + dt * vapourf0(1:qv%n) / 2.0_r8
        hqc_new%f(1:qc%n)       = hqc%f(1:qc%n)        + dt * cloudf0(1:qc%n) / 2.0_r8
        hqr_new%f(1:qr%n)       = hqr%f(1:qr%n)        + dt * rainf0(1:qr%n) / 2.0_r8
        if(testcase==2) htracer_new%f(1:qr%n)   = htracer%f(1:qr%n)    + dt * tracerf0(1:qr%n) / 2.0_r8

    else if(testcase==1)then
        htracer_new%f(1:qr%n)   = htracer%f(1:qr%n)    + dt * tracerf0(1:qr%n) / 2.0_r8
    end if
    call tendency_moist_swm(h_new, u_new, htheta_new, hqv_new, hqc_new, hqr_new, htracer_new, &
    massf1, momf1, tempf1, vapourf1, cloudf1, rainf1, tracerf1, t1)

    !Second RK step
    t2 = t0 + dt/2._r8

    if(testcase>1)then
        u_new%f(1:u%n)          = u%f(1:u%n)           + dt * momf1(1:u%n) / 2.0_r8
        h_new%f(1:h%n)          = h%f(1:h%n)           + dt * massf1(1:h%n) / 2.0_r8
        htheta_new%f(1:theta%n) = htheta%f(1:theta%n)  + dt * tempf1(1:theta%n) / 2.0_r8
        hqv_new%f(1:qv%n)       = hqv%f(1:qv%n)        + dt * vapourf1(1:qv%n) / 2.0_r8
        hqc_new%f(1:qc%n)       = hqc%f(1:qc%n)        + dt * cloudf1(1:qc%n) / 2.0_r8
        hqr_new%f(1:qr%n)       = hqr%f(1:qr%n)        + dt * rainf1(1:qr%n) / 2.0_r8
        if(testcase==2) htracer_new%f(1:qr%n)   = htracer%f(1:qr%n)    + dt * tracerf1(1:qr%n) / 2.0_r8

    else if(testcase==1)then
        htracer_new%f(1:qr%n)   = htracer%f(1:qr%n)    + dt * tracerf1(1:qr%n) / 2.0_r8
    end if

    call tendency_moist_swm(h_new, u_new, htheta_new, hqv_new, hqc_new, hqr_new, htracer_new,&
    massf2, momf2, tempf2, vapourf2, cloudf2, rainf2, tracerf2, t2)


    !Third  RK step
    t3 = t0 + dt

    if(testcase>1)then
      u_new%f(1:u%n)          = u%f(1:u%n)           + dt * momf2(1:u%n)
      h_new%f(1:h%n)          = h%f(1:h%n)           + dt * massf2(1:h%n) 
      htheta_new%f(1:theta%n) = htheta%f(1:theta%n)  + dt * tempf2(1:theta%n)
      hqv_new%f(1:qv%n)       = hqv%f(1:qv%n)        + dt * vapourf2(1:qv%n) 
      hqc_new%f(1:qc%n)       = hqc%f(1:qc%n)        + dt * cloudf2(1:qc%n)
      hqr_new%f(1:qr%n)       = hqr%f(1:qr%n)        + dt * rainf2(1:qr%n)
      if(testcase==2) htracer_new%f(1:qr%n)   = htracer%f(1:qr%n)    + dt * tracerf2(1:qr%n)

    else if(testcase==1)then
      htracer_new%f(1:qr%n)   = htracer%f(1:qr%n)    + dt * tracerf2(1:qr%n)
    end if

    call tendency_moist_swm(h_new, u_new, htheta_new, hqv_new, hqc_new, hqr_new, htracer_new, &
    massf3, momf3, tempf3, vapourf3, cloudf3, rainf3, tracerf3, t3)


    !
    ! Combine them to estimate the solution at time t+dt
    !
    if(mono_filter<=1) then
        if(testcase>1)then
            u_new%f(1:u%n) = u%f(1:u%n) + dt * (momf0(1:u%n)+2._r8*momf1(1:u%n) &
            +2._r8*momf2(1:u%n)+momf3(1:u%n))/6._r8

            h_new%f(1:h%n) = h%f(1:h%n) + dt * (massf0(1:h%n)+2._r8*massf1(1:h%n) &
            +2._r8*massf2(1:h%n)+massf3(1:h%n))/6._r8

            htheta_new%f(1:theta%n) = htheta%f(1:theta%n) + dt * (tempf0(1:theta%n)+2._r8*tempf1(1:theta%n) &
            +2._r8*tempf2(1:theta%n)+tempf3(1:theta%n))/6._r8

            hqv_new%f(1:qv%n) = hqv%f(1:qv%n) + dt * (vapourf0(1:qv%n)+2._r8*vapourf1(1:qv%n) &
            +2._r8*vapourf2(1:qv%n)+vapourf3(1:qv%n))/6._r8

            hqc_new%f(1:qc%n) = hqc%f(1:qc%n) + dt * (cloudf0(1:qc%n)+2._r8*cloudf1(1:qc%n) &
            +2._r8*cloudf2(1:qc%n)+cloudf3(1:qc%n))/6._r8

            hqr_new%f(1:qr%n) = hqr%f(1:qr%n) + dt * (rainf0(1:qr%n)+2._r8*rainf1(1:qr%n) &
            +2._r8*rainf2(1:qr%n)+rainf3(1:qr%n))/6._r8

            if(testcase==2) then
                htracer_new%f(1:qr%n) = htracer%f(1:qr%n) + dt * (tracerf0(1:qr%n)+2._r8*tracerf1(1:qr%n) &
                +2._r8*tracerf2(1:qr%n)+tracerf3(1:qr%n))/6._r8
            endif

        else if (testcase==1) then
          htracer_new%f(1:qr%n) = htracer%f(1:qr%n) + dt * (tracerf0(1:qr%n)+2._r8*tracerf1(1:qr%n) &
          +2._r8*tracerf2(1:qr%n)+tracerf3(1:qr%n))/6._r8
        end if
    else
        print*, 'ERROR on ode_rk4_moist_swm: mono_filter=2 is not implemented for RK4, only for RK3.'
        stop
    end if
 
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
      type(scalar_field), intent(inout):: htheta_new  !General

      !Vapour (defined on voronoi centers)
      type(scalar_field), intent(inout):: hqv_new  !General

      !Cloud (defined on voronoi centers)
      type(scalar_field), intent(inout):: hqc_new  !General

      !Rain (defined on voronoi centers)
      type(scalar_field), intent(inout):: hqr_new  !General

      !Tracer (defined on voronoi centers)
      type(scalar_field), intent(inout):: htracer_new  !General
      
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

      ! Compute source
      if(testcase>1) call source(h, htheta, hqv, hqc, hqr, dt)

      !-------------------------------------------------------------------------
      call tendency_moist_swm(h, u, htheta, hqv, hqc, hqr, htracer, massf0, momf0, tempf0, vapourf0, cloudf0, rainf0, tracerf0, t)

      !First RK step
      if(testcase>1)then
          u_new%f(1:u%n)          = u%f(1:u%n)           + dt * momf0(1:u%n) / 3.0_r8
          h_new%f(1:h%n)          = h%f(1:h%n)           + dt * massf0(1:h%n) / 3.0_r8
          htheta_new%f(1:theta%n) = htheta%f(1:theta%n)  + dt * tempf0(1:theta%n) / 3.0_r8
          hqv_new%f(1:qv%n)       = hqv%f(1:qv%n)        + dt * vapourf0(1:qv%n)  / 3.0_r8
          hqc_new%f(1:qc%n)       = hqc%f(1:qc%n)        + dt * cloudf0(1:qc%n)   / 3.0_r8
          hqr_new%f(1:qr%n)       = hqr%f(1:qr%n)        + dt * rainf0(1:qr%n)    / 3.0_r8
 
          if (testcase==2) htracer_new%f(1:htracer%n) = htracer%f(1:htracer%n)+dt*tracerf0(1:htracer%n)/3.0_r8

      else if (testcase==1)then
          htracer_new%f(1:htracer%n) = htracer%f(1:htracer%n)+dt*tracerf0(1:htracer%n)/3.0_r8
      end if

      !-------------------------------------------------------------------------
      ! update u source
      !if(testcase>1) call source_u(h_new, u_new, htheta_new)

      !-------------------------------------------------------------------------
      call tendency_moist_swm(h_new, u_new, htheta_new, hqv_new, hqc_new, hqr_new, htracer_new, &
      massf1, momf1, tempf1, vapourf1, cloudf1, rainf1, tracerf1, t+dt/3._r8)

      !Second RK step
      if(testcase>1)then
          u_new%f(1:u%n)          = u%f(1:u%n)           + dt * momf1(1:u%n) / 2.0_r8
          h_new%f(1:h%n)          = h%f(1:h%n)           + dt * massf1(1:h%n) / 2.0_r8
          htheta_new%f(1:theta%n) = htheta%f(1:theta%n)  + dt * tempf1(1:theta%n) / 2.0_r8
          hqv_new%f(1:qv%n)       = hqv%f(1:qv%n)        + dt * vapourf1(1:qv%n)  / 2.0_r8
          hqc_new%f(1:qc%n)       = hqc%f(1:qc%n)        + dt * cloudf1(1:qc%n)   / 2.0_r8
          hqr_new%f(1:qr%n)       = hqr%f(1:qr%n)        + dt * rainf1(1:qr%n)    / 2.0_r8
 
          if (testcase==2) htracer_new%f(1:htracer%n) = htracer%f(1:htracer%n)+dt*tracerf1(1:htracer%n)/2.0_r8

      else if (testcase==1)then
          htracer_new%f(1:htracer%n) = htracer%f(1:htracer%n)+dt*tracerf1(1:htracer%n)/2.0_r8
      end if

      !-------------------------------------------------------------------------
      ! update u source
      !if(testcase>1)  call source_u(h_new, u_new, htheta_new)

      !-------------------------------------------------------------------------
      call tendency_moist_swm(h_new, u_new, htheta_new, hqv_new, hqc_new, hqr_new, htracer_new,&
      massf2, momf2, tempf2, vapourf2, cloudf2, rainf2, tracerf2, t+dt/2._r8)
 
      !-------------------------------------------------------------------------
      ! Third  RK step
      ! Last RK step applies a different approach if the monotonic limiter is active 
      if(mono_filter<=1) then
        if(testcase>1)then
            u_new%f(1:u%n)          = u%f(1:u%n)           + dt * momf2(1:u%n)
            h_new%f(1:h%n)          = h%f(1:h%n)           + dt * massf2(1:h%n)
            htheta_new%f(1:theta%n) = htheta%f(1:theta%n)  + dt * tempf2(1:theta%n)
            hqv_new%f(1:qv%n)       = hqv%f(1:qv%n)        + dt * vapourf2(1:qv%n)
            hqc_new%f(1:qc%n)       = hqc%f(1:qc%n)        + dt * cloudf2(1:qc%n)
            hqr_new%f(1:qr%n)       = hqr%f(1:qr%n)        + dt * rainf2(1:qr%n)


            if(testcase==2) htracer_new%f(1:htracer%n) = htracer%f(1:htracer%n)+dt*tracerf2(1:htracer%n)

        else if (testcase==1)then
          htracer_new%f(1:htracer%n) = htracer%f(1:htracer%n)+dt*tracerf2(1:htracer%n)
        end if

      else 
        ! Update tracers
        if(testcase>1)then
            call flux_limiter(mesh, htheta , htheta_new , htheta_new, u, uhtheta, dt, erad, hStheta)
            call flux_limiter(mesh, hqv    , hqv_new    , hqv_new   , u, uhQv   , dt, erad, hSqv   )
            call flux_limiter(mesh, hqr    , hqr_new    , hqr_new   , u, uhQr   , dt, erad, hSqr   )
            call flux_limiter(mesh, hqc    , hqc_new    , hqc_new   , u, uhQc   , dt, erad, hSqc   )
            if (testcase==2) call flux_limiter(mesh, htracer, htracer_new, htracer_new, u, uhtracer, dt, erad)

            !Update momentum and mass conservation equations
            u_new%f(1:u%n)          = u%f(1:u%n)           + dt * momf2(1:u%n)
            h_new%f(1:h%n)          = h%f(1:h%n)           + dt * massf2(1:h%n)
            !theta_new%f(1:theta%n) = htheta%f(1:theta%n)  + dt * (tempf2(1:theta%n) + hStheta%f(1:theta%n))
        else if (testcase==1)then
            call  flux_limiter(mesh, htracer, htracer_new, htracer_new, u, uhtracer, dt, erad)
        end if

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

    !Monotonic filter flag
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

    ! Advection scheme
    select case(advmtd)
    case('upw','trsk','og2','og3','og4','sg2','sg3','sg4')
        swmname=trim(swmname)//"_advmethod_"//trim(advmtd)
    case default
        print*, "Invalid advection method. Please select a proper order."
        stop
    end select

    if (time_integrator == 'rk3' .or. time_integrator == 'rk4') then
        swmname=trim(swmname)//"_"//trim(time_integrator)
    else
        print*, "Invalid time integrator. Please select a proper method."
        stop
    endif

    ! Monotonic limiter
    select case(mono)
    case('0')
       mono_filter=0
    case('1')
       mono_filter=1
    case('2')
       mono_filter=2
    case default
        print*, "Invalid monotonic filter. Please select a proper one."
        stop
    end select
    swmname=trim(swmname)//"_mono"//trim(mono)
    print*, "SWM Name for Plots: ", trim(swmname)
    print*

    return
  end subroutine swm_phys_pars

  subroutine write_swmp_error_file(time, k)
    !----------------------------------------------------------
    ! Write the errors in a file
    !----------------------------------------------------------
    !Scalar fields
    !type(scalar_field), intent(in):: h, h_exact         !Height field
    !type(scalar_field), intent(in):: u, u_exact          !Velocity - normal component
    type(scalar_field) :: grad_gh,grad_gh_error,grad_gh_exact
    type(scalar_field) :: grad_K,grad_K_error,grad_K_exact
    
    real(r8), intent(in) :: time
    integer(i4), intent(in) :: k


    !Error and norm
    real(r8):: errormax_h, errormaxrel_h, error2_h
    real(r8):: errormax_u, errormaxrel_u, error2_u
    real(r8):: errormax_theta, errormaxrel_theta, error2_theta
    real(r8):: errormax_qv, errormaxrel_qv, error2_qv
    real(r8):: errormax_qc, errormaxrel_qc, error2_qc
    real(r8):: errormax_qr, errormaxrel_qr, error2_qr

    integer(i4) :: l
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

  subroutine tendency_advection(q, rhs_q, u, uq, mesh)
    !--------------------------------------
    !Calculates the Right Hand Side (spatial discret./tendency)
    !   of advection equation
    !-------------------------------------------------------------

    !Grid
    type(grid_structure), intent(inout) :: mesh

    !Tracer (defined on voronoi centers)
    type(scalar_field), intent(inout):: q  !General

    !Velocities (defined on edges - only normal component)
    type(scalar_field), intent(inout):: u  !General

    ! Normal flux (defined on edges - only normal component)
    type(scalar_field), intent(inout):: uq  !General

    !Right hand side of advction equation
    real(r8), intent(inout)::rhs_q(:)

    !===============================================================
    !Calculate tracer phi tendency
    !===============================================================
    !Calculate divergence / tracer eq RHS
    phi%f = q%f
    call divhx(phi, u, div_uphi, uq, phi_ed, mesh)
    rhs_q = -div_uphi%f

  return
  end subroutine tendency_advection


  subroutine divhx(q, u, div, uq, q_ed, mesh)
    !---------------------------------------------------------------
    !Calculate divergence at voronoi cells (hexagons)
    !   based on edge normal velocities
    !---------------------------------------------------------------
    type(grid_structure), intent(inout) :: mesh
    type(scalar_field), intent(inout):: u ! velocity at cell edges
    type(scalar_field), intent(inout):: q ! scalar at cell center
    type(scalar_field), intent(inout):: div   ! divergence - must be already allocated
    type(scalar_field), intent(inout):: uq    ! flux at cell edges
    type(scalar_field), intent(inout):: q_ed  ! q at cell edges
    integer(i4) :: i, j, e
    real(r8):: signcor
    real(r8):: area

    if (advmtd=='og2' .or. advmtd=='og3' .or. advmtd=='og4') then ! Ollivier-Gooch method
        call flux_og(q, uq, mesh) 

    else if (advmtd=='trsk' .or. advmtd=='sg2' .or. advmtd=='sg3' .or. advmtd=='sg4') then ! Skamarock-Gassman method
        call flux_sg(q, u, q_ed, uq, mesh) 
    
    else if(advmtd=='upw') then
        call flux_lo_upw(q, u, uq, mesh)

    else
        print*, 'ERROR on tendency_advection. Invalid advection method: ', advmtd
        stop
    endif

    !$omp parallel do &
    !$omp default(none) &
    !$omp shared(mesh, uq, div) &
    !$omp private(j, e, signcor) &
    !$omp schedule(static)
    !Divergence of uq
    do i = 1, mesh%nv
       !For all edges forming the hexagon
       div%f(i) = 0._r8
       do j = 1, mesh%v(i)%nnb
          !Get edge index
          e = mesh%v(i)%ed(j)

          !Get edge outer normal related to the hexagon
          signcor = real(mesh%hx(i)%ttgout(j), r8)

          !Calculate numerical integration
          div%f(i) = div%f(i) + signcor*uq%f(e)
       end do
       div%f(i)=div%f(i)/mesh%hx(i)%areag/erad
    end do
    !$omp end parallel do

    return
  end subroutine divhx

  subroutine flux_sg(q, u, q_ed, uq, mesh)  
    !----------------------------------------------------------------------------------------------
    ! Compute the Skamarock and Gassman flux
    !----------------------------------------------------------------------------------------------
    implicit none 
    type(grid_structure), intent(inout) :: mesh
    type(scalar_field), intent(inout) :: q    ! u at Voronoi centers
    type(scalar_field), intent(inout) :: q_ed ! q at edges
    type(scalar_field), intent(inout) :: u    ! u at edges
    type(scalar_field), intent(inout) :: uq   ! uq at edges (flux)
    integer(i4) :: i, j, k, l, n, e
    integer(i4) :: i1, i2
    integer(i4) :: j1, j2
    real(r8) :: signcor, sign_u
    real(r8) :: normal_2n_derivative_i, normal_2n_derivative_k
    real(r8) :: dist, aux1, aux2, aux3
    real(r8) :: beta
    real(r8) :: p1(1:3), p2(1:3), nr(1:3), dp

    ! Interpolate q to edges and calculate flux at edges
    call scalar_hx2ed(q, q_ed, mesh)      !cell->edge
    beta = 0.25_r8
    !beta = 1._r8
   
    ! Polynomial reconstruction of q
    if(advmtd=='sg3' .or. advmtd=='sg4') then
      ! right hand side for least squares fitting
      call rhs_sg(q, mesh)

      ! least squares fitting
      call least_squares_fitting_sg(q, mesh) 

      ! compute 2nd order derivatives using the least squares fitting
      call normal_2nd_derivative(q, mesh)
    endif

    !$OMP PARALLEL DO &
    !$OMP DEFAULT(NONE) & 
    !$OMP SHARED(mesh, uq, q, q_ed, u, advmtd, beta) &
    !$OMP PRIVATE(i1, i2, j1, j2, e) &
    !$OMP PRIVATE(dist, sign_u, signcor, normal_2n_derivative_i, normal_2n_derivative_k) &
    !$OMP PRIVATE(aux1, aux2, aux3) &
    !$OMP SCHEDULE(static)
    do e = 1, mesh%ne
      ! Get neighboring cells (i1 and i2) that share edge e
      i1 = mesh%edhx(e)%sh(1)
      i2 = mesh%edhx(e)%sh(2)

      ! Edge index on each cell
      j1 = mesh%edhx(e)%sh_local_index(1)
      j2 = mesh%edhx(e)%sh_local_index(2)

      ! Distance between the vertice i1 and its neighbor
      dist = mesh%ed(e)%leng
      
      aux1 = q_ed%f(e)

      !nr = mesh%ed(e)%tg
      !p1 = mesh%v(i1)%p
      !p2 = mesh%v(i2)%p
      !dp = dot_product(nr,p2-p1)
      !signcor = dsign(1._r8, dp)
      signcor = real(mesh%hx(i1)%ttgout(j1), r8) ! is this correct?
      !print*, signcor, dsign(1._r8, dp)
 

      ! sign of u normal component
      if (advmtd=='sg3' .or. advmtd=='sg4') then
         sign_u = sign(1._r8, signcor*u%f(e))
 
         !------------------------------------------------------------------------
         !node i1 rotation values and derivatives
         normal_2n_derivative_i = q%pol(i1)%n2d(j1)
         normal_2n_derivative_k = q%pol(i2)%n2d(j2)

         !------------------------------------------------------------------------
         ! Compute the flux
         if(advmtd=='sg3') then
            aux2 = -(1._r8/12._r8)*((dist)**2)*(normal_2n_derivative_k + normal_2n_derivative_i)
            aux3 = beta*sign_u*(1._r8/12._r8)*((dist)**2)*(normal_2n_derivative_k - normal_2n_derivative_i)
         else if (advmtd=='sg4')  then
            aux3 = 0._r8
            aux2 = -(1._r8/12._r8)*((dist)**2)*(normal_2n_derivative_k + normal_2n_derivative_i)
         endif
      end if

      if (advmtd=='sg3' .or. advmtd=='sg4') then
         uq%f(e) = (aux1+aux2+aux3)*mesh%edhx(e)%leng*u%f(e)
      else !sg2
         uq%f(e) = aux1*mesh%edhx(e)%leng*u%f(e)
      end if

      ! Store flux at edge
      !node(i1)%edge_flux(j1) = signcor*uq%f(e)
      !node(i2)%edge_flux(j2) = -node(i1)%edge_flux(j1)
    end do
    !$OMP END PARALLEL DO

    return  
  end subroutine flux_sg


  subroutine normal_2nd_derivative(q, mesh)  
    !----------------------------------------------------------------------------------------------
    ! Compute 2nd order derivatives at Voronoi centers on each normal direction
    ! using the lsq polynomial
    !----------------------------------------------------------------------------------------------
    implicit none 
    type(grid_structure),intent(inout) :: mesh
    type(scalar_field), intent(inout) :: q    ! q at Voronoi centers
    integer(i4) :: i, j
    real(r8) :: nr(1:2)

    !$OMP PARALLEL DO &
    !$OMP DEFAULT(NONE) & 
    !$OMP SHARED(mesh, q) &
    !$OMP PRIVATE(i, j, nr) &
    !$OMP SCHEDULE(static)
    do i = 1, mesh%nv
       do j = 1, mesh%v(i)%nnb
          !------------------------------------------------------------------------
          !node i rotation values and evaluate 2nd order normal derivatives
          nr = q%pol(i)%tgplane_nr(j)%v(1:2)
          q%pol(i)%n2d(j)  = 2._r8*(q%pol(i)%c(4)*nr(1)*nr(1) + &
                                    q%pol(i)%c(5)*nr(1)*nr(2) + &
                                    q%pol(i)%c(6)*nr(2)*nr(2))
          !------------------------------------------------------------------------
       end do
    end do
    !$OMP END PARALLEL DO
    return  
  end subroutine normal_2nd_derivative  

  subroutine rhs_sg(q, mesh)  
    !--------------------------------------------------------------------------------
    ! Compute the right-hand side (RHS) of the least squares system
    ! for the Skamarock and Gassmann scheme
    !----------------------------------------------------------------------------------------------
    implicit none
    type(grid_structure), intent(in) :: mesh
    type(scalar_field)  , intent(inout) :: q
    integer(i4):: i, j

    !-------------------------------------------------------------------------
    !$OMP PARALLEL DO &
    !$OMP DEFAULT(NONE) & 
    !$OMP SHARED(q, mesh) &
    !$OMP PRIVATE(i, j) &
    !$OMP SCHEDULE(static)
    do i = 1, mesh%nv
       do j = 1, mesh%v(i)%nnb
          q%pol(i)%rhs(j) = q%f(mesh%v(i)%nb(j)) - q%f(i)
       end do
    end do
    !$OMP END PARALLEL DO
  end subroutine rhs_sg



  subroutine least_squares_fitting_sg(q, mesh)  
    !----------------------------------------------------------------------------------------------
    ! Compute the least squares polynomial fitting for the Skamarock and Gassman scheme
    ! 2nd order polynomial
    !----------------------------------------------------------------------------------------------
    implicit none
    type(grid_structure), intent(in) :: mesh
    type(scalar_field)  , intent(inout) :: q
    integer(i4):: i, j
    integer(i4):: m, n
    integer(i4):: l, c
    real(r8):: aux

    !$OMP PARALLEL DO &
    !$OMP DEFAULT(NONE) & 
    !$OMP SHARED(mesh, q) &
    !$OMP PRIVATE(i, l, c, m, aux) &
    !$OMP SCHEDULE(static)
    do i=1,mesh%nv
       l=ubound(q%pol(i)%lsq_matrix,1)
       c=ubound(q%pol(i)%lsq_matrix,2)
       q%pol(i)%c(1) = q%f(i)
       do m=1,c
          aux=0.0D0
          do n=1,l
             aux = aux + q%pol(i)%lsq_matrix_pinv(m,n)*q%pol(i)%rhs(n)
          end do
          q%pol(i)%c(m+1)=aux
       end do
    end do
    !$OMP END PARALLEL DO
    return 
  end subroutine least_squares_fitting_sg

  subroutine flux_og(q, uq, mesh)  
     !----------------------------------------------------------------------------------------------
     ! Compute the Ollivier-Gooch flux
     !----------------------------------------------------------------------------------------------
     implicit none 
     type(grid_structure), intent(inout) :: mesh
     type(scalar_field), intent(inout) :: q    ! q at Voronoi centers
     type(scalar_field), intent(inout) :: uq   ! uq at edges (flux)
     integer(i4):: e, k, i1, i2, iupw
     integer(i4):: nquad
     real(r8) :: pol, dot_prod
     real(r8) :: x, y, z
     real(r8) :: xp, yp, zp
     real(r8) :: cx, cy, sx, sy
     real(r8) :: p(1:3), nr(1:3)
 
     nquad = gauss_quad%n
     !node(1:mesh%nv)%phi_new2 = q%f

     call rhs_og(q, mesh)

     call reconstruction_og(mesh) 

     !$OMP PARALLEL DO &
     !$OMP DEFAULT(NONE) & 
     !$OMP SHARED(node, mesh, q, uq, nquad, gauss_quad, order) &
     !$OMP PRIVATE(i1, i2, iupw) &
     !$OMP PRIVATE(pol, dot_prod, nr) &
     !$OMP PRIVATE(cx, sx, cy, sy, xp, yp) &
     !$OMP SCHEDULE(static)
     do e = 1, mesh%ne
        ! Flux at edges
        uq%f(e) = 0._r8

        !Get neighbor Voronoi cells
        i1 = mesh%edhx(e)%sh(1)
        i2 = mesh%edhx(e)%sh(2)
 
        ! Get normal vector
        nr = mesh%ed(e)%tg

        do k = 1, nquad
           ! Compute the dot product of the normal vector with the velocity vector
           !p  = gauss_quad%edge(e)%node(k)%p 
           dot_prod = dot_product(gauss_quad%edge(e)%u(k)%v, nr)
           !dot_prod = gauss_quad%edge(e)%u(k)%v(1)

           if (dot_prod>=0._r8) then
               iupw = i1
               xp = gauss_quad%edge(e)%tgp_p1(k)%p(1)
               yp = gauss_quad%edge(e)%tgp_p1(k)%p(2)
           else
               iupw = i2
               xp = gauss_quad%edge(e)%tgp_p2(k)%p(1)
               yp = gauss_quad%edge(e)%tgp_p2(k)%p(2)
           end if

           ! Get rotation coefs
           cx = q%pol(iupw)%cx
           sx = q%pol(iupw)%sx
           cy = q%pol(iupw)%cy
           sy = q%pol(iupw)%sy

           ! Polynomial evaluated at the quadrature point
           if (order==2) then
              pol = node(iupw)%coef(1) + node(iupw)%coef(2)*xp + & 
                    node(iupw)%coef(3)*yp

           else if(order==3) then
              pol = node(iupw)%coef(1) + node(iupw)%coef(2)*xp + node(iupw)%coef(3)*yp + &
                    node(iupw)%coef(4)*xp*xp + node(iupw)%coef(5)*xp*yp + node(iupw)%coef(6)*yp*yp

           else if(order==4) then     
              pol = node(iupw)%coef(1) + node(iupw)%coef(2)*xp + node(iupw)%coef(3)*yp + &
                    node(iupw)%coef(4)*xp*xp + node(iupw)%coef(5)*xp*yp + node(iupw)%coef(6)*yp*yp + &                     
                    node(iupw)%coef(7)*xp*xp*xp + node(iupw)%coef(8)*xp*xp*yp + &
                    node(iupw)%coef(9)*xp*yp*yp + node(iupw)%coef(10)*yp*yp*yp
           end if

           ! Flux update
           uq%f(e) = uq%f(e) + pol*gauss_quad%edge(e)%w(k)*dot_prod

        enddo
     enddo
     !$OMP END PARALLEL DO
     return  
  end subroutine flux_og


  subroutine rhs_og(q, mesh)
    !----------------------------------------------------------------------------------------------
    ! Compute the RHS for the OG reconstruction
    !----------------------------------------------------------------------------------------------
    implicit none
    type(scalar_field), intent(inout) :: q  ! q at Voronoi centers
    type(grid_structure), intent(in) :: mesh
    integer(i4) :: i, j, k, l
    real(r8) :: m
 
    !$OMP PARALLEL DO &
    !$OMP DEFAULT(NONE) & 
    !$OMP SHARED(node, mesh, q) & 
    !$OMP PRIVATE(i, j, k, l, m) &
    !$OMP SCHEDULE(static)
    do i = 1, mesh%nv
       node(i)%VBO(1) = q%f(i)
       k = ubound(node(i)%geometric,1)+1
       do j = 1,k-1
          l = node(i)%stencil(j)
          node(i)%VBO(j+1) = q%f(l)*node(i)%geometric(j,1)
          m = node(i)%geometric(j,1)
          node(i)%VBO(j+1) = node(i)%VBO(j+1) - m*node(i)%VBO(1)
       end do
    enddo
    !$OMP END PARALLEL DO
    return
  end subroutine rhs_og

  subroutine reconstruction_og(mesh)  
    !----------------------------------------------------------------------------------------------
    ! Reconstruction for OG scheme
    !----------------------------------------------------------------------------------------------
    implicit none
    type(grid_structure),intent(in):: mesh
    integer(i4):: i, m, n, l ,c, max_m
    real(r8):: soma

    select case(order)
        case(2); max_m = 3
        case(3); max_m = 6
        case(4); max_m = 10
    end select

    !$OMP PARALLEL DO &
    !$OMP DEFAULT(NONE) & 
    !$OMP SHARED(mesh, node, max_m) & 
    !$OMP PRIVATE(i, l, c, m, soma) &
    !$OMP SCHEDULE(static)
    do i = 1, mesh%nv
       l = ubound(node(i)%MRO,1)
       c = ubound(node(i)%MRO,2)
       do m = 1, c-1
          soma = 0.d0
          do n = 1, l-1
             soma = soma + node(i)%MPO(m,n)*node(i)%VBO(n+1)
          enddo
          node(i)%coef(m+1) = soma
       enddo
       soma = 0.d0
       do m = 2, max_m
          soma = soma + node(i)%coef(m)*node(i)%MRO(1,m)
       enddo
       node(i)%coef(1) = node(i)%VBO(1)-soma 
    enddo
    !$OMP END PARALLEL DO
    return 
  end subroutine reconstruction_og


subroutine map_edge_to_cell_local_index(mesh)  
    !----------------------------------------------------------------------------------------------
    ! Initialize the local edge indices (sh_local_index) within each neighboring cell.
    ! For every edge 'e', the routine determines the position of that edge
    ! in the neighbor list of the two adjacent cells (i and k).
    !----------------------------------------------------------------------------------------------
    implicit none 
    type(grid_structure),intent(inout) :: mesh
    integer(i4) :: i, j, k, e, l, n

    ! Calculate sh_local_index
    do e = 1, mesh%ne
      ! Get neighboring cells (i and k) that share edge e
      i = mesh%edhx(e)%sh(1)
      k = mesh%edhx(e)%sh(2)

      j = -10
      do n = 1, mesh%v(i)%nnb
         if( e == mesh%v(i)%ed(n) ) j=n
      end do
      if(j==-10) then
         print*, 'Error: couldnt find edge ',e,' in cell ', i
         stop
      endif
      mesh%edhx(e)%sh_local_index(1) = j

      l = -10
      do n = 1, mesh%v(k)%nnb
         if( e == mesh%v(k)%ed(n) ) l=n
      end do
      if(l==-10) then
         print*, 'Error: couldnt find edge ',e,' in cell ', k
         stop
      endif
      mesh%edhx(e)%sh_local_index(2) = l
    end do

   return  
   end subroutine map_edge_to_cell_local_index




  subroutine init_lsq_rotation_sg(q, mesh)  
    !----------------------------------------------------------------------------------------------
    ! Initialize some Least Squares rotation parametes
    ! normal directions (needed by 2nd order derivatives) and 
    !----------------------------------------------------------------------------------------------
    implicit none 
    type(grid_structure),intent(inout) :: mesh
    type(scalar_field), intent(inout) :: q
    integer(i4) :: i, j, k, e, l, n
    real(r8) :: cx, cy, sx, sy
    real(r8) :: x, y, z
    real(r8) :: xp, yp, zp
    real(r8) :: p(1:3), normal_vector(1:3)
    real(r8) :: dir_i(1:3), rot_dir_i(1:3)
    real(r8) :: signcor

    ! Calculate sh_local_index
    do e = 1, mesh%ne
      ! Get neighboring cells (i and k) that share edge e
      i = mesh%edhx(e)%sh(1)
      k = mesh%edhx(e)%sh(2)

      j = -10
      do n = 1, mesh%v(i)%nnb
         if( e == mesh%v(i)%ed(n) ) j=n
      end do
      if(j==-10) then
         print*, 'Error: couldnt find edge ',e,' in cell ', i
         stop
      endif
      mesh%edhx(e)%sh_local_index(1) = j

      l = -10
      do n = 1, mesh%v(k)%nnb
         if( e == mesh%v(k)%ed(n) ) l=n
      end do
      if(l==-10) then
         print*, 'Error: couldnt find edge ',e,' in cell ', k
         stop
      endif
      mesh%edhx(e)%sh_local_index(2) = l
    end do

    !$OMP PARALLEL DO &
    !$OMP DEFAULT(NONE) & 
    !$OMP SHARED(mesh, q) &
    !$OMP PRIVATE(i, j, e, k, p) &
    !$OMP PRIVATE(x, y, z, cx, sx, cy, sy, xp, yp, zp) &
    !$OMP PRIVATE(rot_dir_i, dir_i, signcor, normal_vector) &
    !$OMP SCHEDULE(static)
    ! Rotate Voronoi center points
    do i = 1, mesh%nv
       !------------------------------------------------------------------------
       !node i rotation values
       x = mesh%v(i)%p(1)
       y = mesh%v(i)%p(2)
       z = mesh%v(i)%p(3)

       call constr(x, y, z, q%pol(i)%cx, q%pol(i)%sx, q%pol(i)%cy, q%pol(i)%sy)
       cx = q%pol(i)%cx
       sx = q%pol(i)%sx
       cy = q%pol(i)%cy
       sy = q%pol(i)%sy

       do j = 1, mesh%v(i)%nnb
           ! get edge
           k = mesh%v(i)%nb(j)
           e = mesh%v(i)%ed(j)
           p = mesh%ed(e)%c%p

           ! signcor for normal vector
           signcor = real(mesh%hx(i)%ttgout(j), r8)

           ! Normal vector - pointing outwards
           normal_vector = signcor*mesh%ed(e)%tg

           ! rotate node k
           x = mesh%v(k)%p(1)
           y = mesh%v(k)%p(2)
           z = mesh%v(k)%p(3)
           call aplyr(x,y,z,cx,sx,cy,sy,xp,yp,zp)
           q%pol(i)%tgplane_coords(j)%p(1) = xp
           q%pol(i)%tgplane_coords(j)%p(2) = yp

           ! Compute normal direction
           dir_i = normal_vector
           dir_i = proj_vec_sphere(dir_i, mesh%v(i)%p)
           call aplyr(dir_i(1),dir_i(2),dir_i(3),cx,sx,cy,sy,rot_dir_i(1),rot_dir_i(2),rot_dir_i(3))
           rot_dir_i = rot_dir_i/norm(rot_dir_i)
           q%pol(i)%tgplane_nr(j)%v(1) = rot_dir_i(1)
           q%pol(i)%tgplane_nr(j)%v(2) = rot_dir_i(2)

           ! Compute least squares matrix
           q%pol(i)%lsq_matrix(j,1) = xp
           q%pol(i)%lsq_matrix(j,2) = yp
           q%pol(i)%lsq_matrix(j,3) = xp*xp
           q%pol(i)%lsq_matrix(j,4) = xp*yp
           q%pol(i)%lsq_matrix(j,5) = yp*yp

       end do
       call pseudoinversa(mesh%v(i)%nnb, 5, q%pol(i)%lsq_matrix, q%pol(i)%lsq_matrix_pinv)
    end do
    !$OMP END PARALLEL DO
   return  
   end subroutine init_lsq_rotation_sg


  subroutine flux_lo_upw(q, u, uq, mesh)  
    !----------------------------------------------------------------------------------------------
    ! Compute the upwind flux
    !----------------------------------------------------------------------------------------------
    implicit none 
    type(grid_structure), intent(inout) :: mesh
    type(scalar_field), intent(inout) :: q    ! u at Voronoi centers
    type(scalar_field), intent(inout) :: u    ! u at edges
    type(scalar_field), intent(inout) :: uq   ! uq at edges (flux)
    integer(i4) :: i1, i2
    integer(i4) :: j1, j2
    integer(i4) :: e, k
    integer(i4) :: i, j
    real(r8) :: signcor
    real(r8) :: p1(1:3), p2(1:3), nr(1:3), dp

    !$OMP PARALLEL DO &
    !$OMP DEFAULT(NONE) & 
    !$OMP SHARED(mesh, uq, q, u) &
    !$OMP PRIVATE(i1, i2) &
    !$OMP PRIVATE(p1, p2, nr, dp) &
    !$OMP PRIVATE(signcor) &
    !$OMP SCHEDULE(static)
    do e = 1, mesh%ne
      ! Get neighboring cells (i1 and i2) that share edge e
      i1 = mesh%edhx(e)%sh(1)
      i2 = mesh%edhx(e)%sh(2)

      nr = mesh%ed(e)%tg
      p1 = mesh%v(i1)%p
      p2 = mesh%v(i2)%p
      dp = dot_product(nr,p2-p1)
      signcor = dsign(1._r8, dp)
      !signcor = real(mesh%hx(i1)%ttgout(j1), r8) ! is this correct?

      ! Sign correction 
      if(signcor*u%f(e)>0._r8)then
        uq%f(e) = u%f(e)*mesh%edhx(e)%leng*q%f(i1)
      else
        uq%f(e) = u%f(e)*mesh%edhx(e)%leng*q%f(i2)
      end if
    end do
    !$OMP END PARALLEL DO

    return  
  end subroutine flux_lo_upw

  !-----------------------------------------------------------------------------------
  ! Implementation of the flux correction method from the paper Wang et al 2009 -
  ! "Evaluation of Scalar Advection Schemes in the Advanced Research WRF
  ! Model Using Large-Eddy Simulations of Aerosol–Cloud Interactions"
  ! This routine is applied in the last stage of the RK3 scheme
  !-----------------------------------------------------------------------------------
  subroutine flux_limiter(mesh, phi_init_step, phi_prev_step, phi_new_step, u_init_step, uphi_ho, dt, radius, hSphi)
      type(grid_structure), intent(inout) :: mesh
      type(scalar_field), intent(inout):: phi_init_step ! scalar field at time t
      type(scalar_field), intent(inout):: phi_prev_step ! Scalar field from second step in RK3 (time t+dt/2)
      type(scalar_field), intent(inout):: phi_new_step  ! scalar field at time t+dt
      type(scalar_field), intent(inout):: u_init_step ! velocity at time t
      type(scalar_field), intent(inout):: uphi_ho     ! high order flux 
      type(scalar_field), optional, intent(inout):: hSphi ! Source - optional
      real(r8), intent(in) :: dt, radius ! time-step and sphere radius 
      real(r8) :: flux, div, signcor, r1, r2
      integer(i4) :: i, j, k, e
      integer(i4) :: i1, i2, j1, j2, iupw

      !-----------------------------------------------------------------------------------
      ! Flux for phi_star using 1st order upwind scheme at time t
      if (present(hSphi))then ! check if source was given
        !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
        !$OMP SHARED(phi_lo_upw, phi_init_step, phi_prev_step, phi, dt, hSphi)
        phi_lo_upw%f = phi_init_step%f + dt*hSphi%f ! equation 3b from Wang et al 2009
        phi%f = phi_prev_step%f
        !$OMP END PARALLEL WORKSHARE
      else
        !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
        !$OMP SHARED(phi_lo_upw, phi_init_step, phi_prev_step, phi)
        phi_lo_upw%f = phi_init_step%f ! equation 3b from Wang et al 2009
        phi%f = phi_prev_step%f
        !$OMP END PARALLEL WORKSHARE
      end if
      !-----------------------------------------------------------------------------------


      !-----------------------------------------------------------------------------------
      ! Compute 1st order upwind flux
      call flux_lo_upw(phi_lo_upw, u_init_step, uphi_lo_upw, mesh)  

      !-----------------------------------------------------------------------------------
      ! Equation 5 from Wang et al 2009 - 1st order upwind solution
      !$omp parallel do &
      !$omp default(none) &
      !$omp shared(mesh, phi_lo_upw, uphi_lo_upw, radius, dt) &
      !$omp private (flux, signcor, div, e, j) &
      !$omp schedule(static)
      do i = 1, mesh%nv
         flux = 0._r8
         do j = 1, mesh%v(i)%nnb 
            e = mesh%v(i)%ed(j)
            signcor = real(mesh%hx(i)%ttgout(j), r8)
            flux = flux + signcor*uphi_lo_upw%f(e)
         end do
         div = flux/mesh%hx(i)%areag/radius
         phi_lo_upw%f(i) = phi_lo_upw%f(i) - dt*div
      end do
      !$omp end parallel do
      !-----------------------------------------------------------------------------------


      !-----------------------------------------------------------------------------------
      ! Corrected flux for phi (equation 4 from wang et al 2009)
      ! High order flux (uphi_ho) minus 1st order upwind flux (uphi_lo_upw)
      !$omp parallel do &
      !$omp default(none) &
      !$omp shared(mesh, uphi_ho, uphi_lo_upw, uphi_cor) 
      do e = 1, mesh%ne
         uphi_cor%f(e) = uphi_ho%f(e) - uphi_lo_upw%f(e) 
      enddo
      !$omp end parallel do
      !-----------------------------------------------------------------------------------


      !-----------------------------------------------------------------------------------
      ! Compute negative and positive flux updates
      !$omp parallel do &
      !$omp default(none) &
      !$omp shared(mesh, phi_lo_upw, phi_plus, phi_minus, uphi_cor, dt, radius) &
      !$omp private (flux, e, signcor) &
      !$omp schedule(static)
      do i = 1, mesh%nv
         phi_minus%f(i) = phi_lo_upw%f(i)
         phi_plus %f(i) = phi_lo_upw%f(i)
         do j = 1, mesh%v(i)%nnb
            e = mesh%v(i)%ed(j)
            signcor = real(mesh%hx(i)%ttgout(j), r8)
            if (signcor*uphi_cor%f(e)>0._r8)then
               ! Equation 6a from Wang et al 2009
               phi_plus %f(i) = phi_plus %f(i) - dt*signcor*uphi_cor%f(e)/mesh%hx(i)%areag/radius
            else
               ! Equation 6b from Wang et al 2009
               phi_minus%f(i) = phi_minus%f(i) - dt*signcor*uphi_cor%f(e)/mesh%hx(i)%areag/radius
            end if
        end do
      end do
      !$omp end parallel do
      !-----------------------------------------------------------------------------------

      !-----------------------------------------------------------------------------------
      ! Min/max in neighborhood at time t
      !$omp parallel do &
      !$omp default(none) &
      !$omp shared(mesh, phi_init_step, phi_min, phi_max) &
      !$omp private(j, k) &
      !$omp schedule(static)
      do i = 1, mesh%nv
         phi_min%f(i) = phi_init_step%f(i)
         phi_max%f(i) = phi_init_step%f(i)
         do j = 1, mesh%v(i)%nnb
            k = mesh%v(i)%nb(j)
            if (phi_init_step%f(k) < phi_min%f(i))then
               phi_min%f(i) = phi_init_step%f(k)
            else if(phi_init_step%f(k) > phi_max%f(i))then
               phi_max%f(i) = phi_init_step%f(k)
            end if
         end do
      end do
      !$omp end parallel do
      !-----------------------------------------------------------------------------------

      !-----------------------------------------------------------------------------------
      ! Compute the limiting factor
      ! 
      !$omp parallel do &
      !$omp default(none) &
      !$omp shared(mesh, uphi_cor, u_init_step, r_limiter) &
      !$omp shared(phi_min, phi_max, phi_plus, phi_minus, phi_lo_upw) &
      !$omp private(i1, i2, j1, j2, signcor, iupw, r1, r2)
      do e = 1, mesh%ne
         ! Get neighboring cells (i1 and i2) that share edge e
         i1 = mesh%edhx(e)%sh(1)
         i2 = mesh%edhx(e)%sh(2)

         ! Edge index on each cell
         j1 = mesh%edhx(e)%sh_local_index(1)
         j2 = mesh%edhx(e)%sh_local_index(2)

         signcor = real(mesh%hx(i1)%ttgout(j1), r8) ! is this correct?
         if(signcor<=0.0001) then
            print*, 'ERROR in flux_limiter, signcor is non positive!'
            stop
         endif

         if (u_init_step%f(e)>0._r8) then !normal vector always points from i1 towards i2
            iupw = i1
         else
            iupw = i2
         endif

         r1 = 1._r8
         r2 = 1._r8
  
         if (uphi_cor%f(e)>0._r8)then
           if (abs(phi_lo_upw%f(i1)-phi_plus%f(i1))>eps2*abs(phi_lo_upw%f(i1))) then
              r1 = (phi_lo_upw%f(i1)-phi_min%f(i1))/(phi_lo_upw%f(i1)-phi_plus%f(i1))
           endif

           if (abs(phi_lo_upw%f(i2)-phi_minus%f(i2))>eps2*abs(phi_lo_upw%f(i2))) then
              r2 = (phi_lo_upw%f(i2)-phi_max%f(i2))/(phi_lo_upw%f(i2)-phi_minus%f(i2))
           endif

         else
           if (abs(phi_lo_upw%f(i1)-phi_minus%f(i1))>eps2*abs(phi_lo_upw%f(i1))) then
              r1 = (phi_lo_upw%f(i1)-phi_max%f(i1))/(phi_lo_upw%f(i1)-phi_minus%f(i1))
           endif

           if (abs(phi_lo_upw%f(i2)-phi_plus%f(i2))>eps2*abs(phi_lo_upw%f(i2))) then
              r2 = (phi_lo_upw%f(i2)-phi_min%f(i2))/(phi_lo_upw%f(i2)-phi_plus%f(i2))
           endif

         endif
         r_limiter%f(e) = min(1._r8, r1, r2)

      enddo
      !$omp end parallel do
      !-----------------------------------------------------------------------------------


      !-----------------------------------------------------------------------------------

      !-----------------------------------------------------------------------------------
      ! Corrected divergence
      ! Final solution
      !$omp parallel do &
      !$omp default(none) &
      !$omp shared(mesh, phi_prev_step, phi_lo_upw, phi_init_step, uphi_cor, uphi_ho, radius, dt) &
      !$omp shared(phi_new_step, uphi_lo_upw, r_limiter) &
      !$omp private (flux, signcor, div, e) &
      !$omp schedule(static)
      do i = 1, mesh%nv
         flux = 0._r8
         do j = 1, mesh%v(i)%nnb 
            e = mesh%v(i)%ed(j)
            signcor = real(mesh%hx(i)%ttgout(j), r8)
            flux = flux + signcor*r_limiter%f(e)*uphi_cor%f(e)
            !flux = flux + signcor*(uphi_lo_upw%f(e) + r_limiter%f(e)*uphi_cor%f(e))
            !print*, uphi_cor%f(e), uphi_ho%f(e)-uphi_lo_upw%f(e)
            !flux = flux + signcor*(uphi_ho%f(e))
            !flux = flux + signcor*(uphi_lo_upw%f(e))
            !flux = flux + signcor*(uphi_cor%f(e))
         end do
         flux = flux/mesh%hx(i)%areag/radius
         phi_new_step%f(i) = phi_lo_upw%f(i) - dt*flux
         !phi_new_step%f(i) = phi_init_step%f(i) - dt*flux
      end do
      !$omp end parallel do
      !-----------------------------------------------------------------------------------

  end subroutine flux_limiter


subroutine init_vecrecon_lsqfitpol_ed(e, u, mesh)
!-----------------------------------------------------------------------
! Vector reconstruction at an edge using least-squares polynomial fitting
!-----------------------------------------------------------------------
  integer (i4), intent(in)    :: e
  type(grid_structure), intent(in) :: mesh
  type(scalar_field), intent(inout) :: u

  !--------------------------------------------------
  ! Local variables
  !--------------------------------------------------
  integer (i4) :: i, j, n, l, m, k, ed
  integer (i4) :: v1, v2
  integer (i4) :: n_edges
  integer (i4), allocatable :: ed_list(:)

  real (r8) :: cx, sx, cy, sy
  real (r8) :: x, y, z
  real (r8) :: xp, yp, zp

  real (r8) :: nx, ny, nz
  real (r8) :: nr(3)

  integer (i4) :: n_coefs
  real (r8), allocatable :: a(:,:)

  ! --- geometric weights
  real (r8) :: pref(3)
  real (r8) :: d, dsum, dmax
  real (r8) :: av, rmax, rinv, wt

  real(r8) :: sn, st, Ln, Lt, sn2sum, st2sum, dloc
  ! Linear or quadratic least squares
  if (urecon_mtd=='ed1') then
     n_coefs = 6
  else if (urecon_mtd=='ed2') then
     n_coefs = 12
  else
    print*, 'ERROR on init_vecrecon_lsqfitpol_ed: invalid urecon_mtd, ', urecon_mtd
    stop
  endif

  !Alocate space if necessary
  if (.not. allocated(u%pol)) then
     allocate(u%pol(1:mesh%ne))
  end if

  !==================================================
  ! Build edge stencil (two cells sharing edge e)
  !==================================================

  ! Primal edge endpoints
  i = mesh%edhx(e)%sh(1)
  j = mesh%edhx(e)%sh(2)


  if (urecon_mtd=='ed2') then ! add more points for quadratic LSQ
    !--------------------------------------------------
    ! Dual edge endpoints (Voronoi vertices)
    ! Correspond to circumcenters of the two
    ! Delaunay triangles adjacent to edge e
    !--------------------------------------------------
    v1 = mesh%edhx(e)%v(1)
    v2 = mesh%edhx(e)%v(2)

    !--------------------------------------------------
    ! For each adjacent triangle, find the vertex
    ! opposite to primal edge (i,j)
    !--------------------------------------------------
    l = -1
    do n = 1, 3
       k = mesh%tr(v1)%v(n)
       if (k /= i .and. k /= j) then
          l = k
          exit
       end if
    end do

    m = -1
    do n = 1, 3
       k = mesh%tr(v2)%v(n)
       if (k /= i .and. k /= j) then
          m = k
          exit
       end if
    end do

    if (l<0 .or. m<0) then
      print*, 'ERROR on init_vecrecon_lsqfitpol_ed: couldnt find opposite vertex'
      stop
    endif
  endif

  !--------------------------------------------------
  ! Count stencil edges:
  ! - edge e
  ! - edges incident to vertices i and j
  ! - edges incident to opposite vertices l and m
  !   excluding those already counted via i or j
  !--------------------------------------------------

  n_edges = 1

  do n = 1, mesh%v(i)%nnb
     if (mesh%v(i)%ed(n) /= e) n_edges = n_edges + 1
  end do

  do n = 1, mesh%v(j)%nnb
     if (mesh%v(j)%ed(n) /= e) n_edges = n_edges + 1
  end do

  if (urecon_mtd=='ed2') then ! add more points for quadratic LSQ
    do n = 1, mesh%v(l)%nnb
       if ( .not. any(mesh%v(i)%ed(1:mesh%v(i)%nnb) == mesh%v(l)%ed(n)) .and. &
            .not. any(mesh%v(j)%ed(1:mesh%v(j)%nnb) == mesh%v(l)%ed(n)) ) then
          n_edges = n_edges + 1
       end if
    end do

    do n = 1, mesh%v(m)%nnb
       if ( .not. any(mesh%v(i)%ed(1:mesh%v(i)%nnb) == mesh%v(m)%ed(n)) .and. &
            .not. any(mesh%v(j)%ed(1:mesh%v(j)%nnb) == mesh%v(m)%ed(n)) ) then
          n_edges = n_edges + 1
       end if
    end do
  endif


  allocate(ed_list(n_edges))
  ed_list(1) = e

  k = 1

  !----------------------------------
  ! Edges incident to vertex i
  !----------------------------------
  do n = 1, mesh%v(i)%nnb
     ed = mesh%v(i)%ed(n)
     if (ed /= e) then
        k = k + 1
        ed_list(k) = ed
     end if
  end do


  !----------------------------------
  ! Edges incident to vertex j
  !----------------------------------
  do n = 1, mesh%v(j)%nnb
     ed = mesh%v(j)%ed(n)
     if (ed /= e) then
        k = k + 1
        ed_list(k) = ed
     end if
  end do

  if (urecon_mtd=='ed2') then ! add more points for quadratic LSQ
    !----------------------------------
    ! Edges incident to vertex l
    !----------------------------------
    do n = 1, mesh%v(l)%nnb
       ed = mesh%v(l)%ed(n)
       if (.not. any(ed_list(1:k) == ed)) then
          k = k + 1
          ed_list(k) = ed
       end if
    end do

    !----------------------------------
    ! Edges incident to vertex m
    !----------------------------------
    do n = 1, mesh%v(m)%nnb
       ed = mesh%v(m)%ed(n)
       if (.not. any(ed_list(1:k) == ed)) then
          k = k + 1
          ed_list(k) = ed
       end if
    end do
  end if
     


  !==================================================
  ! Allocate LSQ structures
  !==================================================


  !Check if coeficients already exist,
  ! if so, they will be overwritten
  if (.not. allocated(u%pol(e)%c)) then
     allocate(u%pol(e)%c(1:n_coefs))
  end if

  if (.not. allocated(u%pol(e)%eds)) then
     allocate(u%pol(e)%eds(1:n_edges))
  endif

  ! store edges used by stencil
  do k = 1, n_edges
     u%pol(e)%eds(k) = ed_list(k) 
  end do

  if (.not. allocated(u%pol(e)%lsq_matrix_pinv)) then
     allocate(u%pol(e)%lsq_matrix_pinv(1:n_coefs,1:n_edges))
  endif

 
  allocate(a(n_edges, n_coefs))

  !==================================================
  ! Reference point = center of edge e
  !==================================================

  pref = mesh%ed(e)%c%p
  !pref = mesh%edhx(e)%c%p

  !--------------------------------------------------
  ! Compute angular distances
  !--------------------------------------------------
  if (.not. allocated(u%pol(e)%wt))then
     allocate(u%pol(e)%wt(1:n_edges))
  endif

  if (.not. allocated(u%pol(e)%rhs))then
     allocate(u%pol(e)%rhs(1:n_edges))
  endif

  !==================================================
  ! Rotation to tangent plane at edge e
  !==================================================
  call constr(pref(1), pref(2), pref(3), cx, sx, cy, sy)

  !==================================================
  ! Compute anisotropic length scales in tangent plane
  !==================================================

  sn2sum = 0.0_r8
  st2sum = 0.0_r8

  do l = 1, n_edges

     ed = ed_list(l)

     ! Position in local tangent plane
     call aplyr(mesh%ed(ed)%c%p(1), mesh%ed(ed)%c%p(2), &
                mesh%ed(ed)%c%p(3), cx, sx, cy, sy, xp, yp, zp)

     ! Edge normal (Voronoi normal) in tangent plane
     nr = mesh%ed(ed)%tg
     call aplyr(nr(1), nr(2), nr(3), cx, sx, cy, sy, nx, ny, nz)

     dloc = dsqrt(nx*nx + ny*ny)
     nx = nx / dloc
     ny = ny / dloc

     ! Normal and tangential coordinates
     sn =  xp*nx + yp*ny
     st = -xp*ny + yp*nx

     if (ed /= e) then
        sn2sum = sn2sum + sn*sn
        st2sum = st2sum + st*st
     end if

  end do

  Ln = dsqrt( sn2sum / real(max(1,n_edges-1), r8) )
  Lt = dsqrt( st2sum / real(max(1,n_edges-1), r8) )

  ! Avoid zero scales
  Ln = max(Ln, 1.0e-12_r8)
  Lt = max(Lt, 1.0e-12_r8)

  !==================================================
  ! Assemble weighted LSQ system
  !==================================================

  do l = 1, n_edges

     ed = ed_list(l)

     !-----------------------------
     ! Position in tangent plane
     !-----------------------------
     call aplyr(mesh%ed(ed)%c%p(1), mesh%ed(ed)%c%p(2), &
                mesh%ed(ed)%c%p(3), cx, sx, cy, sy, xp, yp, zp)

     !-----------------------------
     ! Edge normal in tangent plane
     !-----------------------------
     nr = mesh%ed(ed)%tg
     call aplyr(nr(1), nr(2), nr(3), cx, sx, cy, sy, nx, ny, nz)

     dloc = dsqrt(nx*nx + ny*ny)
     nx = nx / dloc
     ny = ny / dloc

     !-----------------------------
     ! Normal / tangential coords
     !-----------------------------
     sn =  xp*nx + yp*ny
     st = -xp*ny + yp*nx

     !-----------------------------
     ! Anisotropic geometric weight
     !-----------------------------
     wt = 1.0_r8 / (1.0_r8 + (sn/Ln)**2 + (st/Lt)**2)
     u%pol(e)%wt(l) = wt

     !-----------------------------
     ! LSQ matrix
     !-----------------------------
     if (urecon_mtd == 'ed1') then
        ! Linear reconstruction (6 DOF)

        a(l,1) = wt * nx
        a(l,2) = wt * ny

        a(l,3) = wt * (xp * nx)
        a(l,4) = wt * (xp * ny)

        a(l,5) = wt * (yp * nx)
        a(l,6) = wt * (yp * ny)

     else if (urecon_mtd == 'ed2') then
        ! Quadratic reconstruction (12 DOF)

        a(l, 1) = wt * nx
        a(l, 2) = wt * ny

        a(l, 3) = wt * (xp * nx)
        a(l, 4) = wt * (xp * ny)

        a(l, 5) = wt * (yp * nx)
        a(l, 6) = wt * (yp * ny)

        a(l, 7) = wt * (xp*xp * nx)
        a(l, 8) = wt * (xp*xp * ny)

        a(l, 9) = wt * (xp*yp * nx)
        a(l,10) = wt * (xp*yp * ny)

        a(l,11) = wt * (yp*yp * nx)
        a(l,12) = wt * (yp*yp * ny)

     end if

  end do

  !==================================================
  ! Solve LSQ
  !==================================================

  call pseudoinversa(n_edges, n_coefs, a, u%pol(e)%lsq_matrix_pinv)

  !Save the rotation parameters
  u%pol(e)%cx=cx
  u%pol(e)%cy=cy
  u%pol(e)%sx=sx
  u%pol(e)%sy=sy

 
  !==================================================
  ! Cleanup
  !==================================================

  deallocate(ed_list, a)

end subroutine init_vecrecon_lsqfitpol_ed


  subroutine init_vecrecon_lsqfitpol_hxe(kc, var, mesh)
    !----------------------------------------------------------
    ! POLINOMIAL 1st ORDER LEAST SQUARE vector reconstruction
    !
    ! On input:
    !       kc : central triangle/node index for reconstruction
    !            calculation
    !       stencil:
    !           - hxe : 12 point centered on node
    !       mesh : with at least 12 nodes
    !       var  : with var%f filled with normal components at edges
    ! On output:
    !      var%pol%c is filled with the coeficients
    !      var%pol%cx/cy/sx/sy is filled with the rotaion parameters
    !
    !----------------------------------------------------------------

    !Node or triangle index serving as center of reconstruction
    integer (i4), intent(in) :: kc

    !Mesh structure
    type(grid_structure), intent(in) :: mesh

    !Variable used in reconstruction
    type(scalar_field), intent(inout) :: var

    !Globally fixed minimum and maximal lengths
    integer, parameter :: lmx=30

    !List of neighbour nodes
    integer (i4):: eds(1:lmx)

    !Floating point tolerance for ill conditioning detection
    real (r8), parameter :: dtol=0.0001

    ! Distances
    real (r8):: d
    real (r8):: dsum
    real (r8):: rmax
    real (r8):: dmax

    !Normal vectors
    real (r8):: nx
    real (r8):: ny
    real (r8):: nz
    real (r8):: nr(1:3)

    !Indexes
    integer (i4):: i
    integer (i4):: j
    integer (i4):: k
    integer (i4):: l
    integer (i4):: n
    integer (i4):: ed
    !integer :: ier

    ! Weighting Average
    real (r8):: av

    ! Inverse of a radius of influence R which
    ! enters into the weights, R = 1+RF unless all ele-
    ! ments of NPTS are used in the fit (LNP =
    ! LMAX+1), in which case R is the distance
    ! function associated with some point more
    ! distant from node than NPTS(LMAX)
    real(r8):: rinv

    ! Components of a plane rotation about the X-
    ! axis and Y-axis. Then define a rotation
    ! of 'node' to the north pole
    real (r8):: cx
    real (r8):: sx
    real (r8):: cy
    real (r8):: sy

    !Givens rotation parameters
    real (r8):: c
    real (r8):: s

    !Coordinates of NP in the rotated coordinate
    !system unless ZP < 0, in which case
    !(XP,YP,0) lies on the equator
    real (r8):: xp
    real (r8):: yp
    real (r8):: zp

    ! Weight for the equation coresponding to a point
    ! WT = (R-D)/(R*D) = 1/D - RIN, where D = 1-ZP
    real (r8):: wt

    ! Polinomial coeficients
    real (r8):: coef(1:6)

    !Augmented Matrix for the least square problem
    real (r8), allocatable :: a(:,:)

    !Reference point
    real (r8):: pref(1:3)

    !Auxiliar variables
    real (r8):: dmin

    !Test for the values position
    if(var%pos/=6)then
       print*, "init_vecrecon_lsqfitpol_hxe ERROR : Variable not on hexagon edge midpoints"
       print*, "pos:", var%pos
       stop
    end if

    !Save the nearest edges to be used
    dsum = 0._r8
    dmax=-10._r8
    dmin=100000.0_r8

    !Check if index given is ok
    if(kc>mesh%nv .or. kc <1)then
       print*, "init_vecrecon_lsqfitpol_hxe error: Index given is not correct", kc, mesh%nv
       stop
    end if

    !Alocate space if necessary
    if(.not.allocated(var%pol))then
       allocate(var%pol(1:mesh%nv))
    end if


    !Check if coeficients already exist,
    ! if so, they will be overwritten
    if(.not.allocated(var%pol(kc)%c))then
       allocate(var%pol(kc)%c(1:6))
       !else !nothing to be done
       !   return
    end if

    pref=mesh%v(kc)%p

    !For each neighbour node
    do i = 1, mesh%v(kc)%nnb
       !Save index of edge
       eds(i)=mesh%v(kc)%ed(i)
       !Calculate distance (-cos(angle))
       d=-dot_product(mesh%v(kc)%p, mesh%ed(eds(i))%c%p)
       dmax=max(d, dmax)
       dsum = dsum + 1. - d ** 2
    end do

    !Now add the surounding edges depending on the stencil
    call addextraeqs()

    !Number of edges used
    n=j

    if(.not.allocated(var%pol(kc)%eds))then
       allocate(var%pol(kc)%eds(1:n))
    endif

    ! store edges used by stencil
    do i = 1, n
      var%pol(kc)%eds(i) = eds(i) 
    end do

    if(.not.allocated(var%pol(kc)%lsq_matrix_pinv))then
       allocate(var%pol(kc)%lsq_matrix_pinv(1:6,1:n))
    endif
    allocate(a(1:n,1:6))

    if(.not.allocated(var%pol(kc)%wt))then
       allocate(var%pol(kc)%wt(1:n))
    endif

    if(.not.allocated(var%pol(kc)%rhs))then
       allocate(var%pol(kc)%rhs(1:n))
    endif

    !Average of distances
    av = dsqrt (dsum / real(n, r8))
    var%pol(kc)%av = av

    !The radius of action of the reconstruction is arbitrarily
    !  increased by 25 percent. This is done because the R
    !  distance must be larger then
    ! the nodes distance, otherwise weights would be null
    rmax = dmax + .25*abs(dmax)
    rinv = 1. / (1. + rmax)
    !print*, node, n
    !print*, "EDS:", eds
    !print*, "AV DMAX RMAX RINV", av, dmax, rmax, rinv

    ! Construct the rotation coeficients
    call constr(pref(1), pref(2), pref(3), cx, sx, cy, sy)

    ! Compute wt
    do  i = 1, n
       ed = var%pol(kc)%eds(i)

       !Apply the rotation for this point
       call aplyr (mesh%ed(ed)%c%p(1), mesh%ed(ed)%c%p(2), &
                   mesh%ed(ed)%c%p(3), cx, sx, cy, sy, xp, yp, zp)

       !Define the weight of the edge
       var%pol(kc)%wt(i) = 1. / (1. - zp) - rinv
       wt = var%pol(kc)%wt(i)

       !rhs
       var%pol(kc)%rhs(i) = (var%f(ed))*wt
    end do

    ! Compute LSQ matrix
    do i=1, n
       ed = var%pol(kc)%eds(i)

       !Apply the rotation for this point
       call aplyr (mesh%ed(ed)%c%p(1), mesh%ed(ed)%c%p(2), &
                   mesh%ed(ed)%c%p(3), cx, sx, cy, sy, xp, yp, zp)

       !Define the weight of the edge
       wt = var%pol(kc)%wt(i) 

       !Get Normal vector component at this edge and rotate it
       nr=mesh%ed(ed)%tg
       call aplyr (mesh%ed(ed)%tg(1), mesh%ed(ed)%tg(2), &
                   mesh%ed(ed)%tg(3), cx, sx, cy, sy, nx, ny, nz)

       !Ignore z normal component and correct norm of 2d vector
       nx=nx/dsqrt(1-nz**2)
       ny=ny/dsqrt(1-nz**2)

       ! Store line
       a(i,1:6)=(/ (xp*nx)*wt/av, (xp*ny)*wt/av, (yp*nx)*wt/av, &
                   (yp*ny)*wt/av, (nx)*wt, (ny)*wt /)
    end do

    !Save the rotation parameters
    var%pol(kc)%cx=cx
    var%pol(kc)%cy=cy
    var%pol(kc)%sx=sx
    var%pol(kc)%sy=sy

    !Compute the pseudoinverse
    call pseudoinversa(n, 6, a, var%pol(kc)%lsq_matrix_pinv)

    deallocate(a)
    return
    contains

    subroutine addextraeqs()
      !---------------------------------------------------
      ! Add 6 extra edges to the least squares
      !---------------------------------------------------
      j=n
      !For each neighbour node
      do i=1, mesh%v(kc)%nnb
         !Edge counter
         j=i+mesh%v(kc)%nnb
         !Get the triangle
         k=mesh%v(kc)%tr(i)
         !Get the edge index that is missing from the list
         eds(j)=0
         do l=1,3
            if(mesh%tr(k)%ed(l) /= eds(i) .and. &
                 mesh%tr(k)%ed(l) /= eds(modint(i+1, mesh%v(kc)%nnb)) ) then
               !Add edge to list
               eds(j)=mesh%tr(k)%ed(l)
               exit
            end if
         end do
         if(eds(j)==0)then
            print*, "vecrecon_lsfit ERROR: could not add surrounding edge"
            stop
         end if
         !Calculate distance (-cos(angle))
         d=-dot_product(mesh%v(kc)%p, mesh%ed(eds(j))%c%p)
         dmax=max(d, dmax)
         dsum = dsum + 1. - d ** 2
      end do

      !set n
      n=j

      !Average of distances
      av = dsqrt (dsum / real(n, r8))

      !The radius
      rmax = dmax + .25*abs(dmax)
      rinv = 1. / (1. + rmax)

    end subroutine addextraeqs

  end subroutine init_vecrecon_lsqfitpol_hxe



  subroutine vecrecon_lsqfitpol_hxe(kc, var)
    !----------------------------------------------------------
    ! POLINOMIAL 1st ORDER LEAST SQUARE vector reconstruction
    !
    ! On input:
    !       kc : central triangle/node index for reconstruction
    !            calculation
    !       stencil:
    !           - hxe : 12 point centered on node
    !       var  : with var%f filled with normal components at edges
    ! On output:
    !      var%pol%c is filled with the coeficients
    !      var%pol%cx/cy/sx/sy is filled with the rotaion parameters
    !
    !----------------------------------------------------------------

    !Node or triangle index serving as center of reconstruction
    integer (i4), intent(in) :: kc

    !Variable used in reconstruction
    type(scalar_field), intent(inout) :: var

    !Indexes
    integer (i4):: i
    integer (i4):: n
    integer (i4):: ed
    !integer :: ier

    ! Polinomial coeficients
    real (r8):: coef(1:6)

    !Test for the values position
    if(var%pos/=6)then
       print*, "vecrecon_lsfitpol ERROR : Variable not on hexagon edge midpoints"
       print*, "pos:", var%pos
       stop
    end if

    ! Compute RHS
    n = size(var%pol(kc)%eds(:))
    do  i = 1, n
       ed = var%pol(kc)%eds(i)

       !rhs
       var%pol(kc)%rhs(i) = (var%f(ed))*var%pol(kc)%wt(i) 
    end do

    coef = matmul(var%pol(kc)%lsq_matrix_pinv,var%pol(kc)%rhs)
    coef(1:4)=coef(1:4)/var%pol(kc)%av 

    !Save the coeficients
    var%pol(kc)%c=coef
  end subroutine vecrecon_lsqfitpol_hxe


  !=====================================================================================
  !    POLINOMIAL Least Square Vector Reconstruction
  !=====================================================================================

  function vecrecon_lsq_hxe_ed (p, var, mesh, e)
    !----------------------------------------------------------
    !
    ! This function assumes that the point p lives in Voronoi edge "e"
    !
    ! Evaluate the reconstructed vector using a least square fit
    ! using the two surrounding Voronoi cells
    !
    ! It uses 12 hexagon midpoint normal vector components for 
    ! each one of the two surrounding Voronoi cells and average them
    ! at the end
    !----------------------------------------------------------------
    !Approximation point
    real (r8), intent(in) :: p(1:3)

    !Variable
    type(scalar_field), intent(inout) :: var

    !Stencil used
    character(len=4) :: stencil="hxe"

    !Mesh
    type(grid_structure), intent(in) :: mesh

    !Edge where p lives
    integer (i4), intent(in):: e
    
    ! Voronoi cells sharing edge e
    integer (i4) :: k1, k2

    !Returning value of approximation
    real (r8):: vecrecon_lsq_hxe_ed(1:3), vecrecon_lsq_ed_1(1:3), vecrecon_lsq_ed_2(1:3)

    if(var%pos/=6)then
       print*, "vecrecon_lsq_hxe_ed warning: vector field given in incorrect position", var%pos
       stop
    end if

    ! Get surrounding Voronoi cells
    k1 = mesh%edhx(e)%sh(1)
    k2 = mesh%edhx(e)%sh(2)

    ! Least reconstruction using Voronoi cell k1
    call vecrecon_lsqfitpol_hxe (k1, var)
    vecrecon_lsq_ed_1 = vecrecon_lsqeval(k1, p, var)

    ! Least reconstruction using Voronoi cell k2
    call vecrecon_lsqfitpol_hxe (k2, var)
    vecrecon_lsq_ed_2 = vecrecon_lsqeval(k2, p, var)

    ! Force vectors to be tangent to the sphere
    vecrecon_lsq_ed_1 = proj_vec_sphere(vecrecon_lsq_ed_1, p)
    vecrecon_lsq_ed_2 = proj_vec_sphere(vecrecon_lsq_ed_2, p)

    ! Average the vectors
    vecrecon_lsq_hxe_ed = (vecrecon_lsq_ed_1+vecrecon_lsq_ed_2)*0.5_r8
    return

  end function vecrecon_lsq_hxe_ed

function vecrecon_lsq_ed(p, u, e) result(vecrecon)
!-----------------------------------------------------------------------
! Reconstruct vector field at point p using LSQ polynomial coefficients
! precomputed at edge e (tangent-plane formulation).
!-----------------------------------------------------------------------

  !-----------------------------
  ! Input
  !-----------------------------
  real (r8), intent(in)            :: p(1:3)     ! Target point on sphere
  integer (i4), intent(in)         :: e           ! Central edge
  type(scalar_field), intent(inout):: u           ! Edge-based scalar field

  !-----------------------------
  ! Output
  !-----------------------------
  real (r8) :: vecrecon(1:3)       ! Reconstructed vector at p

  !-----------------------------
  ! Local variables
  !-----------------------------
  real (r8) :: urecon(1:3)
  real (r8) :: cx, sx, cy, sy      ! Rotation parameters
  real (r8) :: xp, yp, zp          ! Rotated coordinates
  real (r8) :: nx, ny, nz
  real (r8) :: nr(3)

  integer (i4) :: n_edges, l, ed

  !==================================================
  ! Assemble RHS using stored stencil and weights
  !==================================================

  n_edges = size(u%pol(e)%eds)

  do l = 1, n_edges
     ed = u%pol(e)%eds(l)
     u%pol(e)%rhs(l) = u%f(ed) * u%pol(e)%wt(l)
  end do

  !==================================================
  ! Solve LSQ system (coefficients already projected)
  !==================================================

  u%pol(e)%c = matmul(u%pol(e)%lsq_matrix_pinv, u%pol(e)%rhs)

  !==================================================
  ! Load rotation parameters (tangent plane at edge e)
  !==================================================

  cx = u%pol(e)%cx
  cy = u%pol(e)%cy
  sx = u%pol(e)%sx
  sy = u%pol(e)%sy

  !==================================================
  ! Evaluate reconstructed field at target point p
  !==================================================

  ! Rotate target point to tangent plane
  call aplyr(p(1), p(2), p(3), cx, sx, cy, sy, xp, yp, zp)

  if (urecon_mtd=='ed1') then
    ! Polynomial evaluation (linear)
    nr(1) = u%pol(e)%c(1) + u%pol(e)%c(3)*xp + u%pol(e)%c(5)*yp
    nr(2) = u%pol(e)%c(2) + u%pol(e)%c(4)*xp + u%pol(e)%c(6)*yp
    nr(3) = 0.0_r8

  else if (urecon_mtd=='ed2') then
    !----------------------------------
    ! Quadratic polynomial evaluation
    !----------------------------------
    nr(1) = u%pol(e)%c(1)  &
          + u%pol(e)%c(3)  * xp &
          + u%pol(e)%c(5)  * yp &
          + u%pol(e)%c(7)  * xp*xp &
          + u%pol(e)%c(9)  * xp*yp &
          + u%pol(e)%c(11) * yp*yp

    nr(2) = u%pol(e)%c(2)  &
          + u%pol(e)%c(4)  * xp &
          + u%pol(e)%c(6)  * yp &
          + u%pol(e)%c(8)  * xp*xp &
          + u%pol(e)%c(10) * xp*yp &
          + u%pol(e)%c(12) * yp*yp

    nr(3) = 0.0_r8
  endif

  ! Rotate back to physical space and project onto sphere
  call aplyrt(nr(1), nr(2), cx, sx, cy, sy, urecon)
  vecrecon = proj_vec_sphere(urecon, p)

end function vecrecon_lsq_ed


  subroutine reconstruct_velocity_quadrature(mesh, u)  
    !----------------------------------------------------------------------------------------------
    ! Reconstructs the normal velocity at quadrature points
    !----------------------------------------------------------------------------------------------
    implicit none

    type(grid_structure),intent(inout):: mesh
    type(scalar_field), intent(inout) :: u ! normal u at edges center
    real(r8):: urecon(1:3), dmax, emax
    real(r8):: p(1:3)
    real(r8):: u0, utmp, vtmp, lat, lon, aux, aux1, aux2,error, uexact(1:3), soma
    integer(i4):: i, j, k, e, q, jj, jend, s
    integer(i4):: nquad
    integer(i4):: i1, i2

    ! Number of quadrature points
    nquad = gauss_quad%n
    error = 0.d0

    if(useStagHTC)then

    !$omp parallel do &
    !$omp default(none) &
    !$omp shared(mesh, u, nquad, gauss_quad, urecon_mtd) &
    !$omp private(i1, i2, p) &
    !$omp private(urecon) &
    !$omp schedule(static)
    do e = 1, mesh%ne ! Edges loop
       do q = 1, nquad ! Quadrature points loop
          !Get neighbor Voronoi cells
          i1 = mesh%edhx(e)%sh(1)
          i2 = mesh%edhx(e)%sh(2)

          ! Quadrature points
          p = gauss_quad%edge(e)%node(q)%p 
     
          ! Reconstruct the velocity field at quadrature points
          if (urecon_mtd == "hx") then
             urecon = vecrecon_lsq_hxe_ed(p, u, mesh, e)

          else if (urecon_mtd == "ed1" .or. urecon_mtd == 'ed2') then
             urecon = vecrecon_lsq_ed (p, u, e)
          endif


          !call cart2sph(p(1), p(2), p(3), lon, lat)
          !u0 = 2._r8*pi*erad/(12._r8*day2sec)
          !utmp = u0*dcos(lat)
          !vtmp = 0._r8
          !call convert_vec_sph2cart(utmp, vtmp, p, uexact)
          !!urecon=uexact 
          !error = max(error, maxval(abs(uexact-urecon))/maxval(abs(uexact)))

          ! Store the velocity
          gauss_quad%edge(e)%u(q)%v = urecon
          !gauss_quad%edge(e)%u(q)%v(1) = u%f(e)
       end do    
    end do
    !$omp end parallel do

    else if (advmtd=='og2') then

    !$omp parallel do &
    !$omp default(none) &
    !$omp shared(mesh, u, nquad, gauss_quad) &
    !$omp private(i1, i2, p) &
    !$omp private(urecon) &
    !$omp schedule(static)
    do e = 1, mesh%ne ! Edges loop
       do q = 1, nquad ! Quadrature points loop
          !Get neighbor Voronoi cells
          !i1 = mesh%edhx(e)%sh(1)
          !i2 = mesh%edhx(e)%sh(2)
      
          ! Quadrature points
          !p = gauss_quad%edge(e)%node(q)%p 
     
          ! Reconstruct the velocity field at quadrature points
          !urecon = vecrecon_lsq_ed(p, u, mesh, e)

          !call cart2sph(p(1), p(2), p(3), lon, lat)
          !u0 = 2._r8*pi*erad/(12._r8*day2sec)
          !utmp = u0*dcos(lat)
          !vtmp = 0._r8
          !call convert_vec_sph2cart(utmp, vtmp, p, uexact)
          !urecon=uexact 
          !error = max(error, maxval(abs(uexact-urecon))/maxval(abs(uexact)))

          ! Store the velocity
          gauss_quad%edge(e)%u(q)%v = u%f(e)
       end do    
    end do
    !$omp end parallel do

    endif
    !print*, error
    !stop

  end subroutine reconstruct_velocity_quadrature

   subroutine init_lsq_rotation_og(q, u, mesh)  
    !----------------------------------------------------------------------------------------------
    ! Initialize some Least Squares rotation parametes (rotation matrices,
    ! rotated quadrature points, ...) needed by OG schemes. 
    !----------------------------------------------------------------------------------------------
    implicit none 
    type(grid_structure),intent(inout) :: mesh
    type(scalar_field), intent(inout) :: q
    type(scalar_field), intent(inout) :: u ! normal velocity
    integer(i4) :: i, j, k, e, l, n, n1, n2
    integer(i4) :: i1, i2, j1, j2
    real(r8)  :: cx, cy, sx, sy
    real(r8)  :: x, y, z
    real(r8)  :: xp, yp, zp
    real(r8) :: p(1:3), normal_vector(1:3)
    real(r8) :: dir_i(1:3), rot_dir_i(1:3)
    real(r8) :: p1(1:3), p2(1:3), nr(1:3)
    real(r8) :: signcor

    !$OMP PARALLEL DO &
    !$OMP DEFAULT(NONE) & 
    !$OMP SHARED(mesh, q) &
    !$OMP PRIVATE(i, j, e, k, p) &
    !$OMP PRIVATE(x, y, z, cx, sx, cy, sy, xp, yp, zp) &
    !$OMP PRIVATE(rot_dir_i, dir_i, signcor, normal_vector) &
    !$OMP SCHEDULE(static)
    ! Rotate Voronoi center points
    do i = 1, mesh%nv
       !------------------------------------------------------------------------
       !node i rotation values
       x = mesh%v(i)%p(1)
       y = mesh%v(i)%p(2)
       z = mesh%v(i)%p(3)

       call constr(x, y, z, q%pol(i)%cx, q%pol(i)%sx, q%pol(i)%cy, q%pol(i)%sy)
    end do
    !$OMP END PARALLEL DO

    do e = 1, mesh%ne ! Edges loop
       do n = 1, gauss_quad%n ! Quadrature points loop
          !Get surronding Voronoi cells that share edge e
          i1 = mesh%edhx(e)%sh(1)
          i2 = mesh%edhx(e)%sh(2)
      
          ! Quadrature points
          p1 = gauss_quad%edge(e)%node(n)%p
          nr = gauss_quad%edge(e)%nr(n)%v

          !----------------------------------------------------------
          ! apply rotations of cell i and j at the quadrature points
          !----------------------------------------------------------
          x = p1(1)
          y = p1(2)
          z = p1(3)

          !----------------------------------------------------------
          ! cell i
          cx = q%pol(i1)%cx
          sx = q%pol(i1)%sx
          cy = q%pol(i1)%cy
          sy = q%pol(i1)%sy

          call aplyr(x,y,z,cx,sx,cy,sy,xp,yp,zp)

          gauss_quad%edge(e)%tgp_p1(n)%p(1) = xp
          gauss_quad%edge(e)%tgp_p1(n)%p(2) = yp

          !----------------------------------------------------------
          ! cell j
          cx = q%pol(i2)%cx
          sx = q%pol(i2)%sx
          cy = q%pol(i2)%cy
          sy = q%pol(i2)%sy

          call aplyr(x,y,z,cx,sx,cy,sy,xp,yp,zp)

          gauss_quad%edge(e)%tgp_p2(n)%p(1) = xp
          gauss_quad%edge(e)%tgp_p2(n)%p(2) = yp

       end do
    end do 
 
    ! init velocity field reconstruction parameters
    if (useStagHTC) then
       if (urecon_mtd=="hx") then
         do i = 1, mesh%nv
            call init_vecrecon_lsqfitpol_hxe (i, u, mesh)
         end do

       else if (urecon_mtd=="ed1" .or. urecon_mtd=="ed2" ) then
         do e = 1, mesh%ne
            call init_vecrecon_lsqfitpol_ed(e, u, mesh)
         end do

       endif
    endif
   return  
   end subroutine init_lsq_rotation_og


  subroutine gauss_legendre_quadrature(x, w, nquad)  
    !----------------------------------------------------------------------------------------------
    ! Computes the nodes (roots) and weights for Gauss–Legendre quadrature.
    !----------------------------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: nquad
    integer :: i, j
    real(r8), intent(out):: x(1:nquad)
    real(r8), intent(out):: w(1:nquad)
    real(r8) :: a, p1, p2, p3, xe
    real(r8) :: pp, xd, xl, xm, z, z1
    real(r8) :: eps
 
    eps = 3.0D-14
    xe  = -1.0d0
    xd  =  1.0d0
    xl  =  (xd-xe)*0.50d0
    xm  =  (xd+xe)*0.50d0

    do i = 1, nquad
       z  = dcos((dacos(-1.0D0))*(dble(i)-0.25d0)/(dble(nquad)+0.25d0))
       a  = 2.0d0
       do while (a .gt. eps) 
          p1 = 1.0d0
          p2 = 0.0d0
          do j = 1, nquad
             p3 = p2
             p2 = p1
             p1 = ((2.0d0*dble(j)-1.0d0)*z*p2-(dble(j)-1.0d0)*p3)/(dble(j))
          end do
          pp = (dble(nquad)*(z*p1-p2))/(z*z-1.0d0)
          z1 = z
          z = z1-(p1/pp)
          a = dabs(z-z1)
       end do
       x(i)      = xm+xl*z
       x(nquad+1-i) = xm-xl*z
       w(i)      = (2.0d0*xl)/((1.0d0-z*z)*pp*pp)
       w(nquad+1-i) = w(i)
    enddo
    return
  end subroutine gauss_legendre_quadrature

  subroutine gaussedges2(mesh, nquad)  
    !----------------------------------------------------------------------------------------------
    ! Computes Gauss-Legendre quadrature points and weights for each edge of a Voronoi control volume.
    !----------------------------------------------------------------------------------------------
    implicit none
    integer, intent(in):: nquad
    type(grid_structure),intent(in) :: mesh       
    integer(i4):: e, n
    integer(i4):: i1, i2
    real(r8), allocatable:: x(:) 
    real(r8), allocatable:: w(:)   
    real(r8) :: p1(1:3), p2(1:3), p3(1:3)
    real(r8) :: nr(1:3), tg(1:3)
    real(r8):: alfa, theta

    ! Allocate arrays for quadrature points and weights
    allocate(x(1:nquad), w(1:nquad))
    call gauss_legendre_quadrature(x, w, nquad)

    ! Loop over all edges
    do e = 1, mesh%ne
       ! Get vertex indices of the edge
       i1 = mesh%ed(e)%sh(1)
       i2 = mesh%ed(e)%sh(2)

       ! Get circumcenter positions of neighboring triangles (vertices of Voronoi polygon)
       p1 = mesh%tr(i1)%c%p
       p2 = mesh%tr(i2)%c%p

       ! Original tangent/normal vector associated with the edge
       nr = mesh%ed(e)%tg

       ! Loop over Gauss points along the edge
       do n = 1, nquad
          ! Compute tangent vector along the arc (project p2 onto tangent of p1)
          p3 = (p2 - p1*dot_product(p1, p2))/norm(p2 - p1*dot_product(p1,p2))
        
          ! Angle between p1 and p2 (arc length on unit sphere)
          theta = dacos(dot_product(p1, p2)/(norm(p1)*norm(p2)))        
          alfa = theta*((x(n)+1) / 2.d0)  

          ! Compute position of Gauss point on the unit sphere along the edge
          gauss_quad%edge(e)%node(n)%p = p1*dcos(alfa) + p3*dsin(alfa)/norm(p1*dcos(alfa) + p3*dsin(alfa))
          !print*, norm2(gauss_quad%edge(e)%node(n)%p-mesh%ed(e)%c%p)
          !print*, norm2(gauss_quad%edge(e)%node(n)%p-mesh%edhx(e)%c%p)

          ! Assign edge normal vector (fixed for the edge)
          gauss_quad%edge(e)%nr(n)%v = nr! -cross_product(p1, p2)/norm(cross_product(p1, p2))                   

          ! Assign quadrature weight scaled by arc length of the edge
          gauss_quad%edge(e)%w(n) = arclen(p1, p2)*(w(n)/2.d0)
       enddo
       !print*
    enddo
    !stop
    deallocate(x, w)
    return
  end subroutine gaussedges2



  !----------------------------------------------------------------------------------------
  !  Initial condition for tracer
  !----------------------------------------------------------------------------------------
  function f(lon, lat)
    real (r8), intent(in) :: lon
    real (r8), intent(in) :: lat
    real (r8), dimension(3) :: p
    real (r8), dimension(3) :: pc
    real (r8):: latc, lonc
    real (r8):: b, c, r, r1, f0
    real (r8):: f
    integer (i4):: tracer_ic

    tracer_ic = 1

    select case(tracer_ic)
    case(0)
       b=0._r8
       c=0.9_r8
       f0 = 1000._r8/2._r8
       r=1._r8/3._r8
       latc = 0._r8
       lonc = 0._r8
       r1= arcdistll(lon, lat, lonc, latc)
       if(r1<r)then
          h1=f0*(1._r8+dcos(pi*r1/r))
          f=h1
          return
       end if
       f=b

    case(1)
       latc = 0._r8
       lonc = 0._r8
       call sph2cart(lon, lat, p(1), p(2), p(3))
       call sph2cart(lonc, latc, pc(1), pc(2), pc(3))
       f = 1.0_r8*dexp(-5._r8*norm(p-pc)**2)
    end select

  end function f

end module moist_swm
