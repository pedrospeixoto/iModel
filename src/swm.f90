module swm
  !=============================================================================
  !  Shallow Water Model
  !
  ! Pack for several simulations on a shallow water model
  !  on the sphere using Voronoi grids
  !
  ! Pedro da Silva Peixoto (pedrosp@ime.usp.br)
  ! Oct 2018
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
    sec2day, &
    day2sec, &
    rotatn

  !Global variables and operators for shallow water model
  use swm_data !Everything
  use swm_operators !Everything

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
    arclen, &
    bar_coord_tr, &
    calc_tiled_areas, &
    convert_vec_sph2cart, &
    cross_product, &
    error_norm_2, &
    error_norm_max, &
    error_norm_max_rel, &
    getunit, &
    getedindexonhx, &
    insidetrmesh, &
    norm, &
    sph2cart, &
    trhx_intersec_areas, &
    cart2sph, &
    proj_vec_sphere, &
    sphtriarea

  use interpack, only: &
    calc_trisk_weights, &
    wachspress_coords_hxv, &
    bar_coord_tr, &
    trsk_order_index, &
    plot_scalarfield, &
    plot_cart_vectorfield

  !Use differential operator routines
  use diffoperpack, only: &
    div_cell_Cgrid, &
    grad_edge_Cgrid


  !use eispack, only: &
  !  rg

  implicit none


contains 

  !======================================================================================
  !    SWM TESTS
  !======================================================================================

  subroutine swm_tests(meshtmp)
    !-----------------------------------------
    !  Main test routine tests routine
    !-----------------------------------------
    !Grid structure (incomming)
    type(grid_structure) :: meshtmp

    !Right hand side of mass equation (number of cell equations)
    !real(r8)::masseq(1:meshtmp%nv)
    !real(r8)::masseq(1:meshtmp%nv)

    !Right hand side of momentum equation (number of edge equations)
    !real(r8)::momeq(1:meshtmp%ne)
    !real(r8)::momeq(1:meshtmp%ne)

    !Errors
    real(r8):: errormaxrel_h=0.
    real(r8):: errormaxrel_u=0.
    real(r8):: errormax_h=0.
    real(r8):: errormax_u=0.
    real(r8):: error2_h=0.
    real(r8):: error2_u=0.
    real(r8):: errormax=0.
    real(r8):: errormaxrel=0.
    real(r8):: error2=0.
    real(r8):: RMSdiv=0.
    real(r8):: maxdiv
    real(r8):: max_gradke
    real(r8):: gradke_tmp
    real(r8):: signcor

    !Indexes
    integer(i4):: i !For node values
    integer(i4):: k !For triangles
    integer(i4):: l !For edges

    !Time in seconds
    real(r8)::time

    !logical :: nonlin_pert=.true.
    real(r8):: nonlin_alpha=1.0
    real(r8):: u00=1.0
    real(r8), dimension(:), allocatable:: h_force, u_force

    !Check for blow ups
    integer(i4)::blowup=0

    !blowup=0

    !Save global variable mesh
    mesh=meshtmp

    !Get test case parameters
    call swmpars(usetime=.true.)

    !Allocate variables
    call allocate_globalswmvars()

    !Read reference data - Spectral model
    !  This is just to get reference times
    call read_refdata_endgame(0._r8)

    !Pre calculate grid properties
    call initialize_gridprop()

    !Initialize fields
    call initialize_fields()

    !Calculate derived initial fields
    !call eqs_spatial_disc(0._r8, dt, h, u, masseq%f, momeq%f)
    call tendency(h, u, masseq%f, momeq%f)
    !h=h_old
    !h=u_old

    !Plot initial fields
    call plotfields(k=0, time=0._r8)

    u_old=u
    h_old=h
    !Calculate total mass
    inimass=sumf_areas(h)

    !Calculate energies
    call calc_energies(Penergy0, Kenergy0, Tenergy0, Availenergy0)

    call write_evol_error(0, 0._r8, inimass, errormax_h, error2_h, errormax_u, error2_u,  errormax, error2, &
      Penergy0, Kenergy0, Tenergy0, Availenergy0, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8)

    call write_evol(0, 0._r8, inimass, Penergy0, Kenergy0, Tenergy0, Availenergy0, 0._r8)

    allocate(u_force(1:u%n))
    allocate(h_force(1:h%n))

    !Time loop
    do k=1, ntime
      !Calculate u and h for time:
      time=real(k, r8)*dt
      ! print*
      !print*, "Time: ", time
      !do l=1, mesh%ne
      ! print*, l, grad_ghbK%f(l), grad_ghbK%f(l)-grad_ghbK_exact%f(l), uhq_perp%f(l), uhq_perp%f(l)-uhq_perp_exact%f(l), momeq(l)
      !end do

      !u%f(1:u%n) = u_old%f(1:u%n) + dt * momeq(1:u%n)
      !h%f(1:h%n) = h_old%f(1:h%n) + dt * masseq(1:h%n)
      call ode_rk4 (time, h_old, u_old, h, u, dt)

      !Linear analysis for Hollingsworth problem
      if(testcase==34 .or. testcase==35)then
        !Set forcing
        if(k==1)then
          u_force=u_exact%f-u%f
          h_force=h_exact%f-h%f

          !Add perturbation to h
          i=18
          h%f(i)=h%f(i)+h%f(i)/100.0
          h%f(mesh%v(i)%nb(1))=h%f(mesh%v(i)%nb(1))+h%f(mesh%v(i)%nb(1))/1000.0
          h%f(mesh%v(i)%nb(2))=h%f(mesh%v(i)%nb(2))+h%f(mesh%v(i)%nb(2))/1000.0
          h%f(mesh%v(i)%nb(3))=h%f(mesh%v(i)%nb(3))+h%f(mesh%v(i)%nb(3))/1000.0
           !call plotfields(0, 1.0_r8)
           !print*, mesh%v(i)%lat
        end if
        u%f=u%f+u_force
        h%f=h%f+h_force

        call error_calc(h, h_exact, h_error, errormaxrel_h, error2_h, errormax_h)
        call error_calc(u, u_exact, u_error, errormaxrel_u, error2_u, errormax_u)

        !alpha*max|u_pert|=u00
        !alpha exp(lambda t) ==> how does alpha behave? alpha=nonlin_alpha
        !nonlin_alpha=0.0001/errormax_u
        nonlin_alpha=0.000001/error2_u

        !Rebuild the total h and u with
        h%f=h_exact%f+nonlin_alpha*(h%f-h_exact%f)
        u%f=u_exact%f+nonlin_alpha*(u%f-u_exact%f)
      end if

      if(RefSolRead)then
        call read_refdata_endgame(time)
         !Sets useRefSol
      else if(RefSolAnal)then
        useRefSol=.true.
      else
        useRefSol=.false.
      end if

      !print*, time, useRefSol
      if( useRefSol  )then
        !Calculate total mass
        Tmass=sumf_areas(h)

        !Calculate erngies
        call calc_energies(Penergy, Kenergy, Tenergy, Availenergy)

        !Calculate Errors and output them
        call error_calc(h, h_exact, h_error, errormaxrel_h, error2_h, errormax_h)
        call error_calc(u, u_exact, u_error, errormaxrel_u, error2_u, errormax_u)
        call error_calc(q_tr, q_tr_exact, q_tr_error, errormaxrel, error2, errormax)
        RMSdiv=dsqrt(sumfsq_areas(divu))
        maxdiv=maxval(abs(divu%f))
        max_gradke=0.0_r8
        do l=1, mesh%ne
          !print*, l
          !Calculate gradient of Ke on edges
          if(useStagHC)then
            signcor=dsign( 1._r8, dot_product(mesh%edhx(l)%nr, &
              mesh%v(mesh%edhx(l)%sh(2))%p-mesh%v(mesh%edhx(l)%sh(1))%p ))
            gradke_tmp=signcor*(ke_hx%f(mesh%edhx(l)%sh(2))-ke_hx%f(mesh%edhx(l)%sh(1)))
          elseif(useStagHTC)then
            ! Obs: the tangent of a triangle edge (which is used as normal
            !   of the voronoi cell edges) is always defined such
            !   that it point from its vertex 1 to its vertex 2
            !   Therefore, the gradient needs no correction term (n_ei)
            gradke_tmp=(ke_hx%f(mesh%ed(l)%v(2))-ke_hx%f(mesh%ed(l)%v(1)))
          end if
          max_gradke=max(max_gradke, abs(gradke_tmp/mesh%ed(l)%leng/erad))
        end do

        call write_evol_error(k, time, (Tmass-inimass)/inimass, errormaxrel_h, error2_h, &
          errormaxrel_u, error2_u,  errormaxrel, error2, &
          (Penergy-Penergy0)/Penergy0, (Kenergy-Kenergy0)/Kenergy0, (Tenergy-Tenergy0)/Tenergy0, &
          (Availenergy-Availenergy0)/Availenergy0, RMSdiv, maxdiv, max_gradke, nonlin_alpha)

        if((errormaxrel_h > 20.0 .or. isnan(errormaxrel_h)).and. k > 3 )then

          blowup=blowup+1
          print*, "System might be unstable, large errors:", errormaxrel_h
          !Plot fields
          call plotfields(ntime, time)
          if(blowup >= 5 .or. isnan(errormaxrel_h) .or. errormaxrel_h > 100000000.0 )then
            print*, "Stopping due to large errors", errormaxrel_h
            exit
          end if
        end if

        if(k==ntime)then
          call write_errors(time, errormaxrel_h, error2_h, errormaxrel_u, error2_u, (Tmass-inimass)/inimass, &
            (Penergy-Penergy0)/Penergy0, (Kenergy-Kenergy0)/Kenergy0, (Tenergy-Tenergy0)/Tenergy0, nonlin_alpha)
        end if

       !No reference or analytic solution
      elseif( k<=2 .or. k==ntime  .or. mod(k,nprints)==0 )then
        !Calculate total mass
        Tmass=sumf_areas(h)

        !Calculate erngies
        call calc_energies(Penergy, Kenergy, Tenergy, Availenergy)
        RMSdiv=dsqrt(sumfsq_areas(divu))

        call write_evol(k, time, (Tmass-inimass)/inimass, &
          (Penergy-Penergy0)/Penergy0, (Kenergy-Kenergy0)/Kenergy0, (Tenergy-Tenergy0)/Tenergy0, &
          (Availenergy-Availenergy0)/Availenergy0, RMSdiv)

      end if

      !Plot fields
      call plotfields(k, time)

      !update fields
      u_old=u
      h_old=h

    end do

    !Print user help for ploting
    if(plots)then
      print*
      print*, "Use ./gmt/anim.sh to see the fluid animation"
      print*, "Use ./gmt/plot.sh to see the scalar/vector field, "//&
        " initial condition and errors"
    end if

  end subroutine swm_tests

  subroutine swm_horiz_loc_trunc_er(meshtmp)
    !-----------------------------------------
    !  Local Truncation Error Analysis for SWM
    !-----------------------------------------
    !Grid structure (incomming)
    type(grid_structure), intent(in) :: meshtmp

    !Errors
    real(r8):: errormaxrel
    real(r8):: error2rel
    real(r8):: errormax

    !Mimetic properties
    real(r8):: prodrule
    real(r8):: hdivu
    real(r8):: ugradh
    real(r8):: coriolis_engy
    real(r8):: crossrule
    real(r8):: crossruletmp
    real(r8):: sumdivuh
    real(r8):: signcor

    !Indexes
    integer(i4):: i
    !integer(i4):: j
    integer(i4):: k
    integer(i4):: l
    integer(i4):: ed
    !integer(i4):: n
    !integer ( i4 ) ierr

    !File name for output
    character (len=256):: filename
    character (len=100):: fmt
    !character (len=100):: atmp

    !File units
    integer (i4):: errorsunit
    !integer (i4):: iunit
    logical::  ifile

    !Save global variable mesh
    mesh=meshtmp

    !Set flag for local truncation error test

    !Get test case parameters
    call swmpars(usetime=.false.)
    test_lterror=1

    !Allocate variables
    call allocate_globalswmvars()

    !Pre calculate grid properties
    call initialize_gridprop()

    !Initialize fields
    call initialize_fields()

    !Calculate initial derived fields
    !call eqs_spatial_disc(0._r8, dt,  h, u, masseq%f, momeq%f)
    call tendency( h, u, masseq%f, momeq%f)

    !File for errors
    filename=trim(datadir)//"swm_loc_trunc_errors.txt"
    call getunit(errorsunit)

    inquire(file=filename, exist=ifile)
    if(ifile)then
      open(errorsunit,file=filename, status='old', position='append')
    else
      open(errorsunit,file=filename, status='replace')
      write(errorsunit, '(a140)') "     Operator       Mesh       MaxErrorRel    "//&
        "    RMSError    MaxError     Methods        Grid"
    end if
    fmt="(a16, i12,  3e18.8, a80)"

    !LTE analysis
    print*, "     Operator       Mesh       MaxErrorRel        RMSError         MaxError"
    !---------------------------------------------------
    ! H interpolation to edges - should be second order on HTC grid
    !---------------------------------------------------
    call error_calc(h_ed, h_ed_exact, h_ed_error, errormaxrel, error2rel, errormax)
    print '(a16, i12,  3e18.8)', " h_ed ", mesh%nv,  errormaxrel, error2rel, errormax
    write(errorsunit, trim(fmt)) " h_ed ", mesh%nv,  errormaxrel, error2rel, errormax, &
      trim(swmname)//" "//trim(mesh%name)

    !---------------------------------------------------
    ! Div of uh
    !---------------------------------------------------
    call error_calc(divuh, divuh_exact, divuh_error, errormaxrel, error2rel, errormax)
    print '(a16, i12,  3e18.8)', " divuh ", mesh%nv,  errormaxrel, error2rel, errormax
    write(errorsunit,trim(fmt)) " divuh ", mesh%nv,  errormaxrel, error2rel, errormax, &
      trim(swmname)//" "//trim(mesh%name)

    !---------------------------------------------------
    ! Kinectic Energy calculation
    !---------------------------------------------------
    call error_calc(ke_hx, ke_hx_exact, ke_hx_error, errormaxrel, error2rel, errormax)
    print '(a16, i12,  3e18.8)', " kenergy ", mesh%nv,  errormaxrel, error2rel, errormax
    write(errorsunit,trim(fmt)) " kenergy ", mesh%nv,  errormaxrel, error2rel, errormax, &
      trim(swmname)//" "//trim(mesh%name)

    !---------------------------------------------------
    ! Kinectic Energy calculation on Triangles
    !---------------------------------------------------
    if(useReconmtdGass)then
      call error_calc(ke_tr, ke_tr_exact, ke_tr_error, errormaxrel, error2rel, errormax)
      print '(a16, i12,  3e18.8)', " kenergy_tr", mesh%nv,  errormaxrel, error2rel, errormax
      write(errorsunit,trim(fmt)) " kenergy_tr", mesh%nv,  errormaxrel, error2rel, errormax, &
        trim(swmname)//" "//trim(mesh%name)
       !print*, Kin_energy_tr%f(1:10)
       !print*, Kin_energy_tr_exact%f(1:10)
    end if

    !---------------------------------------------------
    ! Absolute vorticity
    !---------------------------------------------------
    call error_calc(eta, eta_exact, eta_error, errormaxrel, error2rel, errormax)
    print '(a16, i12,  3e18.8)', " abs_vort ", mesh%nv,  errormaxrel, error2rel, errormax
    write(errorsunit,trim(fmt)) " abs_vort ", mesh%nv,  errormaxrel, error2rel, errormax, &
      trim(swmname)//" "//trim(mesh%name)

    !---------------------------------------------------
    ! Thickness interpolation to tr center
    !---------------------------------------------------
    call error_calc(h_tr, h_tr_exact, h_tr_error, errormaxrel, error2rel, errormax)
    print '(a16, i12,  3e18.8)', " h_tr ", mesh%nv,  errormaxrel, error2rel, errormax
    write(errorsunit,trim(fmt)) " h_tr ", mesh%nv,  errormaxrel, error2rel, errormax, &
      trim(swmname)//" "//trim(mesh%name)

    !---------------------------------------------------
    ! Potential vorticity at tr center
    !---------------------------------------------------
    call error_calc(q_tr, q_tr_exact, q_tr_error, errormaxrel, error2rel, errormax)
    print '(a16, i12,  3e18.8)', " q_tr ", mesh%nv,  errormaxrel, error2rel, errormax
    write(errorsunit,trim(fmt)) " q_tr ", mesh%nv,  errormaxrel, error2rel, errormax, &
      trim(swmname)//" "//trim(mesh%name)

    !---------------------------------------------------
    ! Potential vorticity at cell center
    !---------------------------------------------------
    call error_calc(q_hx, q_hx_exact, q_hx_error, errormaxrel, error2rel, errormax)
    print '(a16, i12,  3e18.8)', " q_hx ", mesh%nv,  errormaxrel, error2rel, errormax
    write(errorsunit,trim(fmt)) " q_hx ", mesh%nv,  errormaxrel, error2rel, errormax, &
      trim(swmname)//" "//trim(mesh%name)

    !---------------------------------------------------
    ! Potential vorticity at edges
    !---------------------------------------------------
    call error_calc(q_ed, q_ed_exact, q_ed_error, errormaxrel, error2rel, errormax)
    print '(a16, i12,  3e18.8)', " q_ed ", mesh%nv,  errormaxrel, error2rel, errormax
    write(errorsunit,trim(fmt)) " q_ed ", mesh%nv,  errormaxrel, error2rel, errormax, &
      trim(swmname)//" "//trim(mesh%name)

    !---------------------------------------------------
    ! Gradient of ghbK - Kenergy calculated
    !---------------------------------------------------
    call error_calc(grad_ghbK, grad_ghbK_exact, grad_ghbK_error, errormaxrel, error2rel, errormax)
    print '(a16, i12,  3e18.8)', " grad_ghbK ", mesh%nv,  errormaxrel, error2rel, errormax
    write(errorsunit,trim(fmt)) " grad_ghbK ", mesh%nv,  errormaxrel, error2rel, errormax, &
      trim(swmname)//" "//trim(mesh%name)

    !---------------------------------------------------
    ! Gradient of ghbK - Kenergy given exactly
    !---------------------------------------------------
    ghbK%f(1:ghbK%n)=grav*(h_exact%f(1:ghbK%n)+bt%f(1:ghbK%n))+ke_hx_exact%f(1:ghbK%n)
    do l=1, mesh%ne
      grad_ghbK%f(l)=grad_edge_Cgrid(l, ghbk, u%pos, gradmtd, mesh)/erad
       !   print*, l, grad_ghbK%f(l)
    end do
    call error_calc(grad_ghbK, grad_ghbK_exact, grad_ghbK_error, errormaxrel, error2rel, errormax)
    print '(a16, i12,  3e18.8)', " grad_ghbK_exact ", mesh%nv,  errormaxrel, error2rel, errormax
    write(errorsunit,trim(fmt)) " grad_ghbK_exact ", mesh%nv,  errormaxrel, error2rel, errormax, &
      trim(swmname)//" "//trim(mesh%name)

    !---------------------------------------------------
    ! Vector reconstruction to triangle centers - vhq_tr
    !---------------------------------------------------
    call error_calc_vec(vhq_tr, vhq_tr_exact, vhq_tr_error, errormaxrel, error2rel, errormax)
    print '(a16, i12,  3e18.8)', " vhq_tr ", mesh%nv,  errormaxrel, error2rel, errormax
    write(errorsunit,trim(fmt)) " vhq_tr ", mesh%nv,  errormaxrel, error2rel, errormax, &
      trim(swmname)//" "//trim(mesh%name)

    !---------------------------------------------------
    ! PV gradient - q_grad_tr
    !---------------------------------------------------
    call error_calc_vec(q_grad_tr, q_grad_tr_exact, q_grad_tr_error, errormaxrel, error2rel, errormax)
    print '(a16, i12,  3e18.8)', " q_grad_tr ", mesh%nv,  errormaxrel, error2rel, errormax
    write(errorsunit,trim(fmt)) " q_grad_tr ", mesh%nv,  errormaxrel, error2rel, errormax, &
      trim(swmname)//" "//trim(mesh%name)

    !---------------------------------------------------
    ! Vector reconstruction to cell centers - vh_hx
    !---------------------------------------------------
    call error_calc_vec(vh_hx, vh_hx_exact, vh_hx_error, errormaxrel, error2rel, errormax)
    print '(a16, i12,  3e18.8)', " vh_hx ", mesh%nv,  errormaxrel, error2rel, errormax
    write(errorsunit,trim(fmt)) " vh_hx ", mesh%nv,  errormaxrel, error2rel, errormax, &
      trim(swmname)//" "//trim(mesh%name)

    !---------------------------------------------------
    ! Vector reconstruction to cell centers multiplied by interpolated PV - vhq_hx
    !---------------------------------------------------
    call error_calc_vec(vhq_hx, vhq_hx_exact, vhq_hx_error, errormaxrel, error2rel, errormax)
    print '(a16, i12,  3e18.8)', " vhq_hx ", mesh%nv,  errormaxrel, error2rel, errormax
    write(errorsunit,trim(fmt)) " vhq_hx ", mesh%nv,  errormaxrel, error2rel, errormax, &
      trim(swmname)//" "//trim(mesh%name)

    !---------------------------------------------------
    ! Perp operator (uhq)
    !---------------------------------------------------
    call error_calc(uhq_perp, uhq_perp_exact, uhq_perp_error, errormaxrel, error2rel, errormax)
    print '(a16, i12,  3e18.8)', " uhq_perp ", mesh%nv,  errormaxrel, error2rel, errormax
    write(errorsunit,trim(fmt)) " uhq_perp ", mesh%nv,  errormaxrel, error2rel, errormax, &
      trim(swmname)//" "//trim(mesh%name)


    !---------------------------------------------------
    ! Laplacian on edges
    !---------------------------------------------------
    call error_calc(lapu, lapu_exact, lapu_error, errormaxrel, error2rel, errormax)
    print '(a16, i12,  3e18.8)', " lapu ", mesh%nv,  errormaxrel, error2rel, errormax
    write(errorsunit,trim(fmt)) " lapu ", mesh%nv,  errormaxrel, error2rel, errormax, &
      trim(swmname)//" "//trim(mesh%name)


    !---------------------------------------------------
    ! Mass equation (overall LTE)
    !---------------------------------------------------
    call error_calc(masseq, masseq_exact, masseq_error, errormaxrel, error2rel, errormax)
    print '(a16, i12,  3e18.8)', " masseq ", mesh%nv,  errormaxrel, error2rel, errormax
    write(errorsunit,trim(fmt)) " masseq ", mesh%nv,  errormaxrel, error2rel, errormax, &
      trim(swmname)//" "//trim(mesh%name)

    !---------------------------------------------------
    ! Momentum equation (overall LTE)
    !---------------------------------------------------
    call error_calc(momeq, momeq_exact, momeq_error, errormaxrel, error2rel, errormax)
    print '(a16, i12,  3e18.8)', " momeq ", mesh%nv,  errormaxrel, error2rel, errormax
    write(errorsunit,trim(fmt)) " momeq ", mesh%nv,  errormaxrel, error2rel, errormax, &
      trim(swmname)//" "//trim(mesh%name)

    !---------------------------------------------------
    ! Mimetic property verification
    !---------------------------------------------------
    print*
    print*, "Verify mimetic properties (normalized)"
    print*, "-------------------------"
    print*
    print*, "1) Mass Conservation"
    print*, " To be analysed on the evolving system"
    print*

    print*, "2) Product rule - normalized fields"
    ! Product rule verification
    !Sum all h*div(u)
    prodrule=0.
    do i=1, mesh%nv
      hdivu=h%f(i)*div_cell_Cgrid(i,u,mesh) /maxh/maxvel !/erad
      divuh_error%f(i)=div_cell_Cgrid(i,u,mesh)/erad
      !print*, i, h%f(i),  div_cell_Cgrid(i,u,mesh), hdivu
      !print*,i, hdivu
      prodrule=prodrule+mesh%hx(i)%areag*hdivu !*erad**2
       !print*, i, hdivu, mesh%hx(i)%areag, erad**2
    end do
    print*, "Sum of h*Div(u): ", prodrule
    !Now sum all u.Grad(h)
    do l=1, mesh%ne
      ugradh=u%f(l)*grad_edge_Cgrid(l, h, u%pos, gradmtd, mesh)/maxh/maxvel !/erad
      prodrule=prodrule+ugradh*mesh%ed(l)%leng*mesh%edhx(l)%leng !*erad**2
    end do
    !prodrule=prodrule*erad**2
    !This should be equal to the sum of div(h*U), which should be zero
    ! by divergence theorem
    print*, "Sum of h*Div(u)+u*Grad(h): ", prodrule
    sumdivuh=0.
    do i=1, mesh%nv
      sumdivuh=sumdivuh+mesh%hx(i)%areag*divuh%f(i)/maxh/maxvel !*erad**2
    end do
    print*, "Sum of Div(hu): ", sumdivuh
    print*, "    Difference:", dabs(sumdivuh-prodrule)
    print*


    print*, "3) Cross rule - normalized fields"
    ! Cross rule verification
    ! It must be valid for any scalar field
    !  we will use the gradient of the h field
    !  that was already calculated with the scalar field
    !  assumed exact
    !Calculate Gradient of h
    do l=1, mesh%ne
      grad_h%f(l)=grad_edge_Cgrid(l, h_ed, u%pos, gradmtd, mesh) !/erad
    end do

    !Calculate cross product of gradient
    crossrule=0.
    do k=1, mesh%nt
      crossruletmp=0.
      !loop over triangle edges
      do l=1, 3
        ed=mesh%tr(k)%ed(l)
        if(useStagHTC)then
          signcor=mesh%tr(k)%tg(l)
        elseif(useStagHC)then
          signcor=dsign( 1._r8, real(mesh%tr(k)%tg(l), r8)* &
            dot_product(mesh%edhx(ed)%nr, mesh%ed(ed)%tg))
           !crossruletmp=crossruletmp+grad_h%f(ed)*mesh%ed(ed)%leng*signcor
        end if
        crossruletmp=crossruletmp+grad_h%f(ed)*mesh%ed(ed)%leng*signcor
         !print*,k, l, grad_h%f(ed)*mesh%ed(ed)%leng*signcor, crossruletmp
      end do
      crossruletmp=crossruletmp/mesh%tr(k)%areag / maxh /maxh !/erad/erad
      crossrule=crossrule+mesh%tr(k)%areag*crossruletmp !*erad**2
       !print*,k, crossruletmp, crossrule
       !crossrule=crossrule+mesh%tr(k)%areag*crossruletmp*erad**2
    end do
    print*, "Sum of k.Rot(Grad(h)):", crossrule
    print*

    call calc_energies(Penergy0, Kenergy0, Tenergy0, Availenergy0)
    print*, "4) Coriolis spurious energy budget - normalized fields "
    coriolis_engy=0.
    do l=1, mesh%ne
      !coriolis_engy=coriolis_engy+ &
      !    uh%f(l)*uhq_perp%f(l)*mesh%edhx(l)%leng*mesh%ed(l)%leng! *0.5_r8 *erad/(maxh)/(maxvel**3)
      coriolis_engy=coriolis_engy+ (uh%f(l)*uhq_perp%f(l))* & !/(maxh)/(maxvel**3))* &
        mesh%edhx(l)%leng*mesh%ed(l)%leng *0.5_r8 !/(maxh)/(maxvel**3)
       ! 2.*sphtriarea(mesh%tr(mesh%ed(l)%sh(1))%c%p, &
       ! mesh%tr(mesh%ed(l)%sh(2))%c%p, mesh%v(mesh%ed(l)%v(1))%p)

       !print*, uh%f(l)*uhq_perp%f(l)*mesh%edhx(l)%leng*mesh%ed(l)%leng *0.5_r8 !*erad/(maxh)/(maxvel**3)
       !print*, uh%f(l)*uhq_perp%f(l)*mesh%edhx(l)%leng*mesh%ed(l)%leng*0.5_r8*erad/(maxh)/(maxvel**3)
    end do
    print*, "Sum of Coriolis energy (Ae*q*h^2*u*u_perp):", coriolis_engy !, &
    print*, "Sum of Coriolis energy / Tenergy:", coriolis_engy/Tenergy0 !, &
    !coriolis_engy*0.5_r8 *erad**2/(maxh)/(maxvel**3)
    !print*, erad, (maxh), (maxvel), (maxvel**3), 0.5_r8 *erad/(maxh)/(maxvel**3)

    !if(testcase==11 .or. testcase==12)then
    !print*
    !print*, "5) Normal modes "
    !call swm_normal_mode_analysis(mesh)
    !print*, "Uncomment call to see results"
    !print*
    !end if

    !Plot fields
    call plotfields(0, 0._r8)

    return
  end subroutine swm_horiz_loc_trunc_er



  subroutine swm_normal_mode_analysis(mesh)
    !-----------------------------------------
    !  Main test routine tests routine
    !-----------------------------------------
    !Grid structure (incomming)
    type(grid_structure), intent(in) :: mesh

    !Right hand side of mass equation (number of cell equations)
    real(r8), allocatable::masseq(:)
    real(r8), allocatable:: masseq0(:)
    !real(r8), allocatable:: masseq(:), masseq0(:)

    !Right hand side of momentum equation (number of edge equations)
    real(r8), allocatable::momeq(:)
    real(r8), allocatable:: momeq0(:)
    !real(r8), allocatable:: masseq(:), masseq0(:)
    !type(scalar_field):: momeq

    !Linearized model matrix
    real(r8), allocatable:: M(:,:)
    !real(r8), allocatable:: Z(:,:)

    type(scalar_field):: h_0  !Initial h
    type(scalar_field):: u_0  !Initial u

    !Eigenvalues
    !real(r8), allocatable:: wr(:) !real
    !real(r8), allocatable:: wi(:) !imag
    !integer (i4) matz


    !Perturbation
    real(r8)::pert

    !Indexes
    integer(i4):: i
    integer(i4):: j
    !integer(i4):: k
    !integer(i4):: l
    !integer(i4):: ed
    integer(i4):: n

    !integer ( kind = 4 ) ierr

    !File name for output
    character (len=256):: filename
    character (len=100):: fmt
    character (len=100):: atmp

    !File units
    integer (i4):: iunit

    !Save global variable mesh
    !mesh=meshtmp

    !System dimension
    n=mesh%nv+mesh%ne
    allocate(masseq(1:mesh%nv))
    allocate(masseq0(1:mesh%nv))
    allocate(momeq(1:mesh%ne))
    allocate(momeq0(1:mesh%ne))
    allocate(M(1:n, 1:n))
    !allocate(Z(1:n, 1:n))
    !allocate(wr(1:n))
    !allocate(wi(1:n))

    allocate(h_0%f(1:h%n))
    allocate(u_0%f(1:u%n))

    !Get test case parameters
    !call swmpars()

    !Allocate variables
    !call allocate_globalswmvars()

    !Pre calculate grid properties
    !call initialize_gridprop()

    !Initialize fields
    !testcase=11
    !call initialize_fields()

    h_0=h
    u_0=u
    print*, "f-sphere:", fsphere, "testcase:", testcase

    !Fix perturbation
    pert = mesh%meanhxarea !*erad

    if(test_lterror==1)then
      test_lterror=0
    end if

    !Calculate initial nonlinear tendency
    !call eqs_spatial_disc(0._r8, dt,  h_0, u_0, masseq0, momeq0)
    call tendency(h_0, u_0, masseq0, momeq0)

    print*
    print*, "Generating Matrices" ! - mass field"

    ! Perturb each input variable in turn
    ! First the mass field
    write(*,'(a9, $)') "Mass:"
    do i=1, mesh%nv
      h=h_0
      u=u_0
      h%f(i)=h_0%f(i)+pert
      !call eqs_spatial_disc(0._r8, 0._r8,  h, u, masseq, momeq)
      call tendency(h, u, masseq, momeq)
      ! And save as column of system matrix
      !$OMP PARALLEL WORKSHARE DEFAULT(shared)
      M(1:mesh%nv,i) = (masseq(1:mesh%nv) - masseq0(1:mesh%nv))/pert
      M(mesh%nv+1:n,i) = (momeq(1:mesh%ne) - momeq0(1:mesh%ne))/pert
      !$OMP END PARALLEL WORKSHARE
      if(mod(i,mesh%nv/10)==0)then
        write(*,'(a1, $)') "*" !print*, bar(1:12)
      end if
    end do
    print*


    !print*, "Generating Matrix - velocity field"
    write(*,'(a9, $)') "Velocity:" !print*, bar(1:12)
    !Now the velocity
    do i=1, mesh%ne
      !print*, "u", i
      h=h_0
      u=u_0
      u%f(i)=u_0%f(i)+pert
      !print*, "h"
      !print*, h%f
      !print*, "u"
      !print*, u%f
      !read*, j
      !call eqs_spatial_disc(0._r8, 0._r8,  h, u, masseq, momeq)
      call tendency( h, u, masseq, momeq)
      !print*, masseq
      !print*, momeq
      ! And save as column of system matrix
      M(1:mesh%nv,mesh%nv+i) = (masseq(1:mesh%nv) - masseq0(1:mesh%nv))/pert
      M(mesh%nv+1:n,mesh%nv+i) = (momeq(1:mesh%ne) - momeq0(1:mesh%ne))/pert
      !print*, "mass"
      !print*, masseq
      !print*, "mom"
      !print*, momeq
      !print*
      if(mod(i,mesh%ne/10)==0)then
        write(*,'(a1, $)') "*" !print*, bar(1:12)
      end if
    end do
    print*

    !print*, "Writing Matrix..."
    !Output matrix
    if(n<100)then
      filename=trim(datadir)//trim(swmname)//"_linearmatrix_"//trim(mesh%name)//".txt"
      call getunit(iunit)
      open(iunit,file=filename, status='replace')

      write(atmp,'(i10)') n
      fmt='('//trim(adjustl(atmp))//'f24.8)'
      !print*, fmt
      do i=1, n
        print *, M(i,1:n)
        write(iunit, fmt) M(i,1:n)
         !read*,j
      end do

    else

      !Output matrix - Sparse print
      filename=trim(datadir)//trim(swmname)//"_linearmatrixsparse_"//trim(mesh%name)//".txt"
      call getunit(iunit)
      open(iunit,file=filename, status='replace')
      write(*,'(a9, $)') "Writing:" !print*, bar(1:12)
      !print*, fmt
      do i=1, n
        !print *, M(i,1:n)
        do j=1, n
          if(dabs(M(i,j))>eps/100._r8)then
            write(iunit, *) i, j, M(i,j)
          end if
        end do
        if(mod(i,n/10)==0)then
          write(*,'(a1, $)') "*" !print*, bar(1:12)
        end if
         !read*,j
      end do
    end if
    print*
    print*, "Matrix printed as file: ", trim(filename)



    print*, "Use MATLAB to calculate the eigenvalues or uncomment the EISPACK-rg call (slow)"
    print*

    !   matz=0
    !   call rg ( n, M(1:n, 1:n), wr(1:n), wi(1:n), matz, z(1:n,1:n), ierr )
    !
    !    !Output eigenvalues
    !    filename=trim(datadir)//trim(swmname)//"_linearmatrixeigen_"//trim(mesh%name)//".txt"
    !    call getunit(iunit)
    !    open(iunit,file=filename, status='replace')
    !    write(iunit, *) "Real    Imag"
    !    do i=1, n
    !       write(iunit, *) wr(i), wi(i)
    !    end do
    !
    !    print*, "Eigenvalues"
    !    j=0 !counter of nonzero eigenvalues
    !    k=0 !counter of non zero real parts
    !    l=0 !counter of non zero imag parts
    !
    !    do i=1,n
    !       if(abs(wr(i))>eps.or.abs(wi(i))>eps)then
    !          j=j+1
    !          if(abs(wr(i))>eps)then
    !             k=k+1
    !          end if
    !          if(abs(wi(i))>eps)then
    !             l=l+1
    !          end if
    !          !print*,i, j, wr(i), wi(i)
    !       end if
    !    end do
    !
    !    print*
    !    print*, "Nonzero     :", j, " of ", n
    !    print*, "Nonzero real:", k, " of ", n
    !    print*, "Nonzero imag:", l, " of ", n
    !
    !    !print*, wr(1:n)
    !    !print*
    !    !print*, wi(1:n)

    test_lterror=1
    return
  end subroutine swm_normal_mode_analysis


  subroutine initialize_gridprop()
    !-----------------------------------------------------
    !  Pre calculate grid dependent only value
    !-----------------------------------------------------

    integer(i4)::i, j, k, l
    real(r8)::ratiotrsk

    !Areas of intersection between primal and dual grids
    ! these are stored in mesh
    call trhx_intersec_areas(mesh)

    !Calculate trisk weights - required if any "trsk" flag set
    call calc_trisk_weights(mesh)

    !Calculate tiled areas
    if(useTiledAreas)then
      call calc_tiled_areas(mesh)
    end if

    !Calculate Wachspress coordinates - necessary only to have pv in hexagons
    !  not in use any more
    if(useSinterpolBary .and. test_lterror==1)then
      do i=1, mesh%nv
        wachc_tr2v(i)=wachspress_coords_hxv(mesh%v(i)%p, mesh, i)
      end do
    end if

    !Pre calculate position of midpoints of hexagonal edges relative to triangles
    if(useSinterpolBary.and.useStagHC)then
      do l=1, mesh%ne
        if(insidetrmesh(mesh%edhx(l)%c%p, mesh%ed(l)%sh(1), mesh))then
          mesh%edhx(l)%c%kt=mesh%ed(l)%sh(1)
        else !if(insidetrmesh(mesh%edhx(l)%c%p, mesh%ed(l)%sh(2), mesh))then
          mesh%edhx(l)%c%kt=mesh%ed(l)%sh(2)
        end if
        allocate(mesh%edhx(l)%c%b(1:3))
        mesh%edhx(l)%c%b=bar_coord_tr(mesh%edhx(l)%c%p, mesh%edhx(l)%c%kt, mesh)
      end do
    end if

    !Pre calculate barycentric coordinates of triangle circumcenters
    if(useSinterpolBary)then
      do k=1, mesh%nt
        allocate(mesh%tr(k)%c%b(1:3))
        mesh%tr(k)%c%b=bar_coord_tr(mesh%tr(k)%c%p, k, mesh)
      end do
    end if

    !Hybrid coriolis recon scheme - possibly to be unstable
    if(useCoriolisMtdHyb)then
      do l=1, mesh%ne
        trskind%f(l)=trsk_order_index(l, u%pos, mesh)
        isTrskindLow(l)=.false.
        if(trskind%f(l)<trskindmax)then
          isTrskindLow(l)=.true.
          ratiotrsk=ratiotrsk+1._r8/mesh%ne
        end if
      end do
      print*, "Percent of TRISK usage in hybrid scheme: ", ratiotrsk*100., "%"
    end if

  end subroutine initialize_gridprop

  subroutine initialize_fields()
    !---------------------------------------------------
    ! initialize_fields
    ! see dataswm.f90 for field that will be initialized
    !---------------------------------

    !Auxiliar variables
    real(r8)::u0
    real(r8):: h0
    real(r8):: utmp
    real(r8):: vtmp
    real(r8):: atmp
    real(r8):: btmp
    real(r8):: ctmp
    !real(r8):: maxvel
    !real(r8):: maxh
    !real(r8):: htmp
    !real(r8):: lon
    !real(r8):: lat
    real(r8):: h_ct
    real(r8):: eta_ct
    real(r8):: grad_ct
    real(r8):: lon0
    real(r8):: lat0
    real(r8):: lon
    real(r8):: lat
    real(r8):: p0(1:3)
    real(r8):: p(1:3)
    real(r8):: r
    real(r8):: rmax
    real(r8):: rmaxc
    real(r8):: vectmp(1:3)
    real(r8):: K_RH !parameters of RH wave
    real(r8):: w !parameters of RH wave
    real(r8):: ctmp1 !constant tmp parameter
    real(r8):: Rmat(1:3,1:3) !Rotation matrix
    real(r8):: RmatT(1:3,1:3) !Transpose/inverse of Rotation matrix
    real(r8):: F, C

    !Indexes
    integer(i4):: i !voronoi cell index
    integer(i4):: k !trinagle index
    integer(i4):: l !edge index
    integer(i4):: m !RH wave number
    integer(i4):: n !Wave number for test case 40

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

    !File name for input
    character (len=256):: filename
    character (len=256):: tctmp, atime

    !File units
    integer (i4):: iunit
    logical::  ifile


    maxvel=1.
    maxh=1.


    select case(testcase)
      case(0) !Validate time discret.
        !      (solid body rot with trigonom h field)
        h%f=0._r8
        do i=1, mesh%nv
          h%f(i)=dcos(mesh%v(i)%lon) &
            *dcos(mesh%v(i)%lat)**2
           !h_exact%f(i)=dcos(pi+mesh%v(i)%lon) &
           !     *dcos(mesh%v(i)%lat)**2
        end do
        !h_0=h
        h_exact=h
        u%f=0

      case(1) !Solid body rotation

        u0=pi2*erad/(12._r8*day2sec)
        h0=1000._r8
        lon0=0._r8 !3._r8*pi*0.5_r8
        lat0=0._r8
        rmax=erad/3._r8
        call sph2cart(lon0, lat0, p0(1), p0(2), p0(3))
        h%f=0._r8

        do i=1, mesh%nv
          r=erad*arclen(mesh%v(i)%p, p0)
          !h%f(i)=0._r8
          if(r<rmax)then
            h%f(i)=(h0*0.5_r8)*(1._r8+dcos(pi*r/rmax))
          endif
          !print*, h%f(i)
          utmp=u0*dcos(mesh%v(i)%lat)
          vtmp=0._r8
          call convert_vec_sph2cart(utmp, vtmp, mesh%v(i)%p, vectmp)
        end do
        !h_0=h
        h_exact=h

        h_ed%f=0._r8
        if(useStagHTC)then
          do l=1, mesh%ne
            utmp=u0*dcos(mesh%ed(l)%c%lat)
            vtmp=0._r8
            call convert_vec_sph2cart(utmp, vtmp, mesh%ed(l)%c%p, vectmp)
            v_ed%p(l)%v=vectmp
            u%f(l)=dot_product(vectmp,mesh%ed(l)%tg)
            r=erad*arclen(mesh%ed(l)%c%p, p0)
            !h%f(i)=0._r8
            if(r<rmax)then
              h_ed%f(i)=(h0*0.5_r8)*(1._r8+dcos(pi*r/rmax))
            endif
          end do
        elseif(useStagHC)then
          do l=1, mesh%ne
            utmp=u0*dcos(mesh%edhx(l)%c%lat)
            vtmp=0._r8
            call convert_vec_sph2cart(utmp, vtmp, mesh%edhx(l)%c%p, vectmp)
            v_ed%p(l)%v=vectmp
            u%f(l)=dot_product(vectmp,mesh%edhx(l)%nr)
            r=erad*arclen(mesh%edhx(l)%c%p, p0)
            !h%f(i)=0._r8
            if(r<rmax)then
              h_ed%f(i)=(h0*0.5_r8)*(1._r8+dcos(pi*r/rmax))
            endif
          end do
        end if
        !u_0=u
        u_exact=u
        v_ed_exact=v_ed
        maxvel=u0

      case(2) !Global Steady State Zonal Geo Flow
        u0=pi2*erad/(12._r8*day2sec)
        h0=2.94e4_r8*gravi
        h_ct=(erad*omega*u0+(u0**2)/2._r8)*gravi
        eta_ct=(2.*u0/erad+2.*omega)
        grad_ct=u0*eta_ct
        do i=1, mesh%nv
          !print*, i
          h%f(i)=h0-h_ct*dsin(mesh%v(i)%lat)**2
          bt%f(i)=0._r8
          divuh%f(i)=0._r8
          !print*, h%f(i)
          !utmp=u0*dcos(mesh%v(i)%lat)
          !vtmp=0._r8
          !call convert_vec_sph2cart(utmp, vtmp, mesh%v(i)%p, vectmp)
          !Calculate exact Kinectic energy
          if(test_lterror==1)then
            ke_hx_exact%f(i)=(u0*dcos(mesh%v(i)%lat))**2/2._r8
            utmp=u0*dcos(mesh%v(i)%lat)
            vtmp=0._r8
            call convert_vec_sph2cart(utmp, vtmp, mesh%v(i)%p, vectmp)
            vh_hx_exact%p(i)%v=vectmp*h%f(i)
            q_hx_exact%f(i)=eta_ct*dsin(mesh%v(i)%lat)/h%f(i)
            vhq_hx_exact%p(i)%v=vh_hx_exact%p(i)%v*q_hx_exact%f(i)
             !print*, i, q_hx_exact%f(i)
          end if
        end do
        !h_0=h
        h_exact=h
        maxh=h0
        if(test_lterror==1)then
          divuh_exact=divuh
        end if

        if(useStagHTC)then
          do l=1, mesh%ne
            lon=mesh%ed(l)%c%lon
            lat=mesh%ed(l)%c%lat
            utmp=u0*dcos(lat)
            vtmp=0._r8
            call convert_vec_sph2cart(utmp, vtmp, mesh%ed(l)%c%p, vectmp)
            v_ed%p(l)%v=vectmp
            u%f(l)=dot_product(vectmp,mesh%ed(l)%tg)
            if(test_lterror==1)then
              h_ed_exact%f(l)=h0-h_ct*dsin(lat)**2
              call convert_vec_sph2cart(0._r8, &
                -grad_ct*dsin(lat)*dcos(lat), mesh%ed(l)%c%p, vectmp)

              grad_ghbK_exact%f(l)=dot_product(vectmp,mesh%ed(l)%tg)
              !print*, eta_ct*dsin(mesh%ed(l)%c%lat), h_ed_exact%f(l)
              q_ed_exact%f(l)=eta_ct*dsin(lat)/h_ed_exact%f(l)
              uhq_perp_exact%f(l)=-eta_ct*dsin(lat)*dot_product(v_ed%p(l)%v,mesh%ed(l)%nr)
              !Laplacian
              utmp=u0*(1.-2.*dcos(lat)**2)/(erad**2*dcos(lat))
              vtmp=0._r8
              call convert_vec_sph2cart(utmp, vtmp, mesh%ed(l)%c%p, vectmp)
              lapu_exact%f(l)=dot_product(vectmp,mesh%ed(l)%tg)
            end if
          end do
        elseif(useStagHC)then
          do l=1, mesh%ne
            lon=mesh%edhx(l)%c%lon
            lat=mesh%edhx(l)%c%lat
            utmp=u0*dcos(lat)
            vtmp=0._r8
            call convert_vec_sph2cart(utmp, vtmp, mesh%edhx(l)%c%p, vectmp)
            v_ed%p(l)%v=vectmp
            u%f(l)=dot_product(vectmp,mesh%edhx(l)%nr)
            if(test_lterror==1)then
              h_ed_exact%f(l)=h0-h_ct*dsin(lat)**2
              call convert_vec_sph2cart(0._r8, &
                -grad_ct*dsin(lat)*dcos(lat), mesh%edhx(l)%c%p, vectmp)
              grad_ghbK_exact%f(l)=dot_product(vectmp,mesh%edhx(l)%nr)
              q_ed_exact%f(l)=eta_ct*dsin(lat)/h_ed_exact%f(l)
              uhq_perp_exact%f(l)=+eta_ct*dsin(lat)*dot_product(v_ed%p(l)%v,mesh%edhx(l)%tg)
              !Laplacian
              utmp=u0*(1.-2.*dcos(lat)**2)/(erad**2*dcos(lat))
              vtmp=0._r8
              call convert_vec_sph2cart(utmp, vtmp, mesh%edhx(l)%c%p, vectmp)
              lapu_exact%f(l)=dot_product(vectmp,mesh%edhx(l)%nr)
            end if
          end do
        end if
        !u_0=u
        u_exact=u
        v_ed_exact=v_ed
        maxvel=u0

        do k=1, mesh%nt
          !Exact PV
          q_tr_exact%f(k)=eta_ct*dsin(mesh%tr(k)%c%lat)/(h0-h_ct*dsin(mesh%tr(k)%c%lat)**2)
          if(test_lterror==1)then
            !Kin energy
            ke_tr_exact%f(k)=(u0*dcos(mesh%tr(k)%c%lat))**2/2._r8
            !Grad of PV
            vtmp=eta_ct*dcos(mesh%tr(k)%c%lat)*(h0+h_ct*dsin(mesh%tr(k)%c%lat)**2) &
              /(erad*(h0-h_ct*dsin(mesh%tr(k)%c%lat)**2)**2)
            utmp=0.
            call convert_vec_sph2cart(utmp, vtmp, mesh%tr(k)%c%p,  q_grad_tr_exact%p(k)%v)
             !print*, k,  vtmp, q_grad_tr_exact%p(k)%v, eta_ct, (h0+h_ct*dsin(mesh%tr(k)%c%lat)**2), &
             !((h0-h_ct*dsin(mesh%tr(k)%c%lat)**2)**2)
          end if
        end do

        if(test_lterror==1)then
          !Loop over triangles
          do k=1, mesh%nt
            !Absolute vorticity
            eta_exact%f(k)=eta_ct*dsin(mesh%tr(k)%c%lat)
            h_tr_exact%f(k)=h0-h_ct*dsin(mesh%tr(k)%c%lat)**2
            q_tr_exact%f(k)=eta_exact%f(k)/h_tr_exact%f(k)
            !print*, k, q_tr_exact%f(k)
            utmp=u0*dcos(mesh%tr(k)%c%lat)
            vtmp=0._r8
            call convert_vec_sph2cart(utmp, vtmp, mesh%tr(k)%c%p, vectmp)
            vhq_tr_exact%p(k)%v=vectmp*eta_exact%f(k)

          end do
          masseq_exact%f=0._r8
          momeq_exact%f=0._r8
        end if

      case(5, 51) !Flow over mountain
        u0=20._r8
        h0=5960._r8
        h_ct=(erad*omega*u0+(u0**2)/2._r8)*gravi

        lon0=-pi*0.5_r8
        lat0=pi/6._r8
        rmax=pi/9._r8

        if(testcase==5)then
          do i=1, mesh%nv
            !print*, i
            lon=mesh%v(i)%lon
            lat=mesh%v(i)%lat
            r=dsqrt((lon-lon0)**2+(lat-lat0)**2)
            h%f(i)=h0-h_ct*dsin(mesh%v(i)%lat)**2


            !print*, r, rmax
            if(r<rmax)then
              bt%f(i)=2000._r8*(1-r/rmax)
               !print*, i,bt%f(i), h%f(i)
            else
              bt%f(i)=0.
            endif
            ! Correct h to allow orography
            h%f(i)=h%f(i)-bt%f(i)
          end do

        elseif(testcase==51)then
          !u0=pi2*erad/(12._r8*day2sec)
          !h0=2.94e4_r8*gravi
          h_ct=(erad*omega*u0+(u0**2)/2._r8)*gravi

          do i=1, mesh%nv

            bt%f(i)=h0-h_ct*dsin(mesh%v(i)%lat)**2
            !print*, i, h0, h_ct, h%f(i)
            call sph2cart (lon0, lat0, p0(1), p0(2), p0(3))

            !Add mountain
            r=norm(mesh%v(i)%p-p0)
            rmaxc=2*(1-cos(rmax))
            bt%f(i)=bt%f(i)+20._r8*exp(-0.4*r**2/rmaxc**2)
            !bt%f(i)=2000._r8*exp(-0.4*r**2/rmaxc**2)

            !  h thin layer to allow orography
            h%f(i)=hollgw-20._r8*exp(-0.4*r**2/rmaxc**2)
          end do
        endif

        !h_0=h
        maxh=h0

        if(useStagHTC)then
          do l=1, mesh%ne
            utmp=u0*dcos(mesh%ed(l)%c%lat)
            vtmp=0._r8
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
        !u_0=u
        u_exact=u
        maxvel=u0

      case(6)

        K_RH=7.848e-6_r8
        w=7.848e-6_r8
        m=4_i4
        !wave_m=8_i4
        h0=8.0e3_r8
        ctmp1=50.001334560000004

        !Velocity field

        if(useStagHTC)then
          do l=1, mesh%ne
            lat=mesh%ed(l)%c%lat
            lon=mesh%ed(l)%c%lon
            utmp= erad*w*dcos(lat)+ &
              erad*K_RH*(dcos(lat))**(m-1)* &
              (m*dsin(lat)**2 - dcos(lat)**2) * dcos(m*lon)
            vtmp = - erad*K_RH*m*((dcos(lat))**(m-1))* &
              dsin(lat)*dsin(m*lon)
            call convert_vec_sph2cart(utmp, vtmp, mesh%ed(l)%c%p, vectmp)
            v_ed%p(l)%v=vectmp
            u%f(l)=dot_product(vectmp,mesh%ed(l)%tg)
            if(test_lterror==1)then
              !Thickness
              h_ed_exact%f(l)=h_tc6(lon, lat)
              !Potential vorticity

              q_ed_exact%f(l)=1.5695581066106648e-7*Sec(lat)*(-800.0213529600001*Cos(lat)**3*Cos(4*lon)*Sin(lat) + &
                Sin(lat)*(ctmp1*Cos(lat) + ctmp1*Cos(lat)**3*Cos(4*lon)*(-Cos(lat)**2 + 4*Sin(lat)**2)) - &
                Cos(lat)*(-ctmp1*Sin(lat) + 500.01334560000004*Cos(lat)**4*Cos(4*lon)*Sin(lat) - &
                150.00400368*Cos(lat)**2*Cos(4*lon)*Sin(lat)*(-Cos(lat)**2 + 4*Sin(lat)**2)))+2.*Omega*dsin(lat)
              q_ed_exact%f(l)=q_ed_exact%f(l)/h_ed_exact%f(l)

              !Calculate perp operator - the normal to triangle edge is the tangent to the voronoi edge
              ! that is defined to be 90deg counter-clockwise rotated from to the tangent of the triangle
              ! So no sign correction needed (see datastruct.f90)
              uhq_perp_exact%f(l)=q_ed_exact%f(l)*h_ed_exact%f(l)*dot_product(v_ed%p(l)%v,mesh%ed(l)%nr)

              !Grad of gh+K
              utmp=7.847790533053324e-8*Sec(lat)*(320017.0825959745*Cos(lat)**6*Cos(4*lon)*Sin(lat)**2*Sin(4*lon) - &
                400.01067648000003*Cos(lat)**3*(-Cos(lat)**2 + 4*Sin(lat)**2)* &
                (ctmp1*Cos(lat) + ctmp1*Cos(lat)**3*Cos(4*lon)*(-Cos(lat)**2 + 4*Sin(lat)**2))*Sin(4*lon)) + &
                6.371219999999999e6*Sec(lat)*(-1.690312704e-10*Cos(lat)**4*(26. - 25*Cos(lat)**2)*Sin(4*lon) - &
                1.23182208e-10*Cos(lat)**8*(-6 + 5*Cos(lat)**2)*Sin(8.*lon))
              vtmp=Cos(lat)*Sin(lat)*(-0.008077015579484158 + Cos(lat)**8* &
                (-0.004905130920335999 - 0.005101336157149441*Cos(4*lon)**2 - 0.004905130920335999*Cos(8.*lon)) + &
                Cos(lat)**2*Cos(4*lon)*(-0.028000320675545084 - 0.006278567578030081*Sin(lat)**2) + &
                Cos(lat)**6*(-0.020405344628597756 + 0.004708925683522559*Cos(8.*lon) + &
                0.025114270312120324*Cos(4*lon)**2*Sin(lat)**2 + &
                0.006278567578030081*Sin(4*lon)**2) + Cos(lat)**4* &
                (0.018835702734090236 + 0.045878824528197124*Cos(4*lon) - &
                0.018835702734090243*Cos(4*lon)**2*Sin(lat)**4 - &
                0.018835702734090243*Sin(lat)**2*Sin(4*lon)**2))
              call convert_vec_sph2cart(utmp, vtmp, mesh%ed(l)%c%p, vectmp)
              grad_ghbK_exact%f(l)=dot_product(vectmp,mesh%ed(l)%tg)
            end if
          end do
        elseif(useStagHC)then
          do l=1, mesh%ne
            lat=mesh%edhx(l)%c%lat
            lon=mesh%edhx(l)%c%lon

            utmp= erad*w*dcos(lat)+ &
              erad*K_RH*(dcos(lat))**(m-1)* &
              (m*dsin(lat)**2 - dcos(lat)**2) * dcos(m*lon)
            vtmp = - erad*K_RH*m*((dcos(lat))**(m-1))* &
              dsin(lat)*dsin(m*lon)
            call convert_vec_sph2cart(utmp, vtmp, mesh%edhx(l)%c%p, vectmp)
            v_ed%p(l)%v=vectmp
            u%f(l)=dot_product(vectmp,mesh%edhx(l)%nr)
            if(test_lterror==1)then

              !Thickness
              h_ed_exact%f(l)=h_tc6(lon, lat)
              !Potential vorticity
              q_ed_exact%f(l)=1.5695581066106648e-7*Sec(lat)*(-800.0213529600001*Cos(lat)**3*Cos(4*lon)*Sin(lat) + &
                Sin(lat)*(ctmp1*Cos(lat) + ctmp1*Cos(lat)**3*Cos(4*lon)*(-Cos(lat)**2 + 4*Sin(lat)**2)) - &
                Cos(lat)*(-ctmp1*Sin(lat) + 500.01334560000004*Cos(lat)**4*Cos(4*lon)*Sin(lat) - &
                150.00400368*Cos(lat)**2*Cos(4*lon)*Sin(lat)*(-Cos(lat)**2 + 4*Sin(lat)**2)))+2.*Omega*dsin(lat)
              q_ed_exact%f(l)=q_ed_exact%f(l)/h_ed_exact%f(l)

              !call convert_vec_sph2cart(-vtmp, utmp, mesh%edhx(l)%c%p, vectmp)
              !Calculate perp operator - the minus sign is the tangent vector is to make the perp always counter-clockwise
              !   with respect to the normal (it is defined in datastruct.f90 as tang 90deg clockwise from normal)
              uhq_perp_exact%f(l)=q_ed_exact%f(l)*h_ed_exact%f(l)*dot_product(v_ed%p(l)%v,-mesh%edhx(l)%tg)


              !Grad of gh+K
              utmp=7.847790533053324e-8*Sec(lat)*(320017.0825959745*Cos(lat)**6*Cos(4*lon)*Sin(lat)**2*Sin(4*lon) - &
                400.01067648000003*Cos(lat)**3*(-Cos(lat)**2 + 4*Sin(lat)**2)* &
                (ctmp1*Cos(lat) + ctmp1*Cos(lat)**3*Cos(4*lon)*(-Cos(lat)**2 + 4*Sin(lat)**2))*Sin(4*lon)) + &
                6.371219999999999e6*Sec(lat)*(-1.690312704e-10*Cos(lat)**4*(26. - 25*Cos(lat)**2)*Sin(4*lon) - &
                1.23182208e-10*Cos(lat)**8*(-6 + 5*Cos(lat)**2)*Sin(8.*lon))
              vtmp=Cos(lat)*Sin(lat)*(-0.008077015579484158 + Cos(lat)**8* &
                (-0.004905130920335999 - 0.005101336157149441*Cos(4*lon)**2 - 0.004905130920335999*Cos(8.*lon)) + &
                Cos(lat)**2*Cos(4*lon)*(-0.028000320675545084 - 0.006278567578030081*Sin(lat)**2) + &
                Cos(lat)**6*(-0.020405344628597756 + 0.004708925683522559*Cos(8.*lon) + &
                0.025114270312120324*Cos(4*lon)**2*Sin(lat)**2 + &
                0.006278567578030081*Sin(4*lon)**2) + Cos(lat)**4* &
                (0.018835702734090236 + 0.045878824528197124*Cos(4*lon) - &
                0.018835702734090243*Cos(4*lon)**2*Sin(lat)**4 - &
                0.018835702734090243*Sin(lat)**2*Sin(4*lon)**2))
              call convert_vec_sph2cart(utmp, vtmp, mesh%edhx(l)%c%p, vectmp)
              grad_ghbK_exact%f(l)=dot_product(vectmp,mesh%edhx(l)%nr)
            end if
          end do
        end if

        !Height field
        do i=1, mesh%nv
          lon=mesh%v(i)%lon
          lat=mesh%v(i)%lat

          !Thickness h
          h%f(i)=h_tc6(lon, lat)

          if(test_lterror==1)then
            !Exact initial divergence of h*u
            if(dabs(lat)>pio2-eps.or. dabs(lon)<eps )then
              divuh_exact%f(i)=0!avoid division by zero
            elseif(dabs(lat)<eps/10000)then
              divuh_exact%f(i)=1.5695581066106648e-7*(0. + 4.1394841903864507e12*(ctmp1 - ctmp1*Cos(4*lon))* &
                (-1.690312704e-10*Sin(4*lon) + 1.23182208e-10*Sin(8.*lon)))
            else
              divuh_exact%f(i)=0.15673434284549875*Cos(lat)**14*(-0.9109225617267884*(1/cos(lat))**10*Sin(4*lon) - &
                56.0567730293408*(1/sin(2*lat))**6*Sin(lat)**8*Sin(8*lon) + &
                0.1531934829004216*(1/cos(lat))**6*Sin(8.*lon) + &
                0.6383061787517568*Sin(4*lon)*Tan(lat)**2 + 0.6383061787517568*Cos(8.*lon)*Sin(4*lon)*Tan(lat)**2 + &
                (2.655353703607308 - 0.6127739316016865*Cos(8.*lon))*(1/cos(lat))**2*Sin(4*lon)*Tan(lat)**2 + &
                (1/cos(lat))**8*Sin(4*lon)*(0.8758870785834504 + 1.*Tan(lat)**2) + &
                (1/cos(lat))**4*(-0.12766123575035135*Sin(8.*lon) - 2.4510957264067463*Sin(4*lon)*Tan(lat)**2) + &
                Cos(4*lon)*(-0.8758870785834504*(1/cos(lat))**4*Sin(4*lon) + &
                0.9109225617267884*(1/cos(lat))**6*Sin(4*lon) + &
                Sin(8.*lon)*(0.12766123575035135 - 0.5106449430014054*Tan(lat)**2) + &
                (1/cos(lat))**2*Sin(8.*lon)*(-0.1531934829004216 + 0.6127739316016864*Tan(lat)**2)))
            end if

            ke_hx_exact%f(i)=((ctmp1*Cos(lat) + ctmp1*Cos(lat)**3*Cos(4*lon)*(-Cos(lat)**2 + 4*Sin(lat)**2))**2 + &
              40002.13532449681*Cos(lat)**6*Sin(lat)**2*Sin(4*lon)**2)/2.
            !Velocities
            utmp= erad*w*dcos(lat)+ &
              erad*K_RH*(dcos(lat))**(m-1)* &
              (m*dsin(lat)**2 - dcos(lat)**2) * dcos(m*lon)
            vtmp = - erad*K_RH*m*((dcos(lat))**(m-1))* &
              dsin(lat)*dsin(m*lon)
            call convert_vec_sph2cart(utmp, vtmp, mesh%v(i)%p, vectmp)
            vh_hx_exact%p(i)%v=vectmp*h%f(i)

            !Potential vorticity
            q_hx_exact%f(i)= 1.5695581066106648e-7*Sec(lat)*(-800.0213529600001*Cos(lat)**3*Cos(4*lon)*Sin(lat) + &
              Sin(lat)*(ctmp1*Cos(lat) + ctmp1*Cos(lat)**3*Cos(4*lon)*(-Cos(lat)**2 + 4*Sin(lat)**2)) - &
              Cos(lat)*(-ctmp1*Sin(lat) + 500.01334560000004*Cos(lat)**4*Cos(4*lon)*Sin(lat) - &
              150.00400368*Cos(lat)**2*Cos(4*lon)*Sin(lat)*(-Cos(lat)**2 + 4*Sin(lat)**2)))+2.*Omega*dsin(lat)
            q_hx_exact%f(i)=q_hx_exact%f(i)/h%f(i)

            vhq_hx_exact%p(i)%v=vh_hx_exact%p(i)%v*q_hx_exact%f(i)
             !print*, i, q_hx_exact%f(i)
          end if

        end do

        if(test_lterror==1)then
          !Loop over triangles
          do k=1, mesh%nt
            lon=mesh%tr(k)%c%lon
            lat=mesh%tr(k)%c%lat
            !Absolute vorticity
            eta_exact%f(k)= 1.5695581066106648e-7*Sec(lat)*(-800.0213529600001*Cos(lat)**3*Cos(4*lon)*Sin(lat) + &
              Sin(lat)*(ctmp1*Cos(lat) + ctmp1*Cos(lat)**3*Cos(4*lon)*(-Cos(lat)**2 + 4*Sin(lat)**2)) - &
              Cos(lat)*(-ctmp1*Sin(lat) + 500.01334560000004*Cos(lat)**4*Cos(4*lon)*Sin(lat) - &
              150.00400368*Cos(lat)**2*Cos(4*lon)*Sin(lat)*(-Cos(lat)**2 + 4*Sin(lat)**2)))
            eta_exact%f(k)=eta_exact%f(k)+2.*Omega*dsin(lat)
             !h_tr_exact%f(k)=h0-h_ct*dsin(mesh%tr(k)%c%lat)**2
             !q_tr_exact%f(k)=eta_exact%f(k)/h_tr_exact%f(k)
             !print*, k, q_tr_exact%f(k)
             !utmp=u0*dcos(mesh%tr(k)%c%lat)
             !vtmp=0._r8
             !call convert_vec_sph2cart(utmp, vtmp, mesh%tr(k)%c%p, vectmp)
             !vhq_tr_exact%p(k)%v=vectmp*eta_exact%f(k)
             !
          end do
          do i=1, mesh%nv
            masseq_exact%f(i)=-divuh_exact%f(i)
          end do
          do l=1, mesh%ne
            momeq_exact%f(l)=-uhq_perp_exact%f(l)-grad_ghbK_exact%f(l)
          end do
        end if

        bt%f=0.
        !h_0=h
        h_exact=h
        u_0=u
        maxh=maxval(h%f(1:mesh%nv))
        maxvel=maxval(u%f(1:mesh%ne))


      case(11, 12) !Linearized equations  - for normal mode calculations
        h%f=1e5_r8*gravi
        !u_0%f=0.
        u%f=0.
        !h=h_0
        bt%f=0.
        if(testcase==11)then
          fsphere=1
          fcte=1.4584e-4_r8
        elseif(testcase==12)then
          fsphere=0
        end if

        !print*, 1e5_r8*gravi,10e5_r8, gravi
        !print*
        !print*, h_0%f
        u_exact=u
        h_exact=h
        masseq_exact%f=0
        momeq_exact%f=0


      case(21, 22, 23) !  ! Galewsky et al test case - From J. Thuburns code

        if(testcase==23)then
          u00 = 200.0
          lat0 = pi/7.0
          lat1 = pi/2.0 - lat0
        else
          u00 = 80.0
          lat0 = pi/7.0
          lat1 = pi/2.0 - lat0
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

        !Add perturbation only in test case 21 or 23
        if(testcase==21 .or. testcase==23)then
          ! Geopotential perturbation
          !if(testcase==21)then
          hpert = 120.0D0
          alpha = 1.0D0/3.0D0
          beta = 1.0D0/15.0D0
          !elseif(testcase==23)then
          ! hpert = 120.0D0
          ! alpha = 1.0D0/2.0D0
          ! beta = 1.0D0/2.0D0
          !end if
          lat2 = 0.5D0*piby2
          do i = 1, mesh%nv
            ! l2 = flat(if0,ngrids)
            ! l1 = flong(if0,ngrids)
            ! CALL centroid(if0,long,lat)
            lat=mesh%v(i)%lat
            long=mesh%v(i)%lon
            l2 = lat
            l1 = long
            clat = COS(l2)
            !IF (l1 > pi) l1 = l1 - 2.0d0*pi
            e1 = EXP(-(l1/alpha)**2)
            e2 = EXP(-((lat2 - l2)/beta)**2)
            h%f(i) = h%f(i)+hpert*clat*e1*e2
          enddo
        end if
        !h_0=h
        !u_0=u
        h_exact=h
        u_exact=u

      case(24) !GLW test case with initial condition read for day 4 form endgame

        print*, "Reading initial data for day 5 from file "
        !Read delayed GLW test case (on day 5)
        tctmp="GLW"
        write(atime,'(I10.10)') 432000 !345600 !5 days
        !Read thickness
        filename=trim(refdir)//trim(tctmp)//"_"// &
          trim(mesh%name)//"_"//trim(atime)//"h.dat"
        print*, trim(filename)
        inquire(file=filename, exist=ifile)
        if(ifile)then
          call getunit(iunit)
          open(iunit,file=filename, status='old',FORM='UNFORMATTED')
          bt%f=0
          do j=1,mesh%nv
            read(iunit) h_exact%f(j)
            h%f(j)=h_exact%f(j)-bt%f(j)
             !print*, j, h_exact%f(j)
          end do
        else
          print*, "Cannot find reference file:", trim(filename)

        end if

        !Read velocities depending on the grid choice
        if(useStagHTC)then
          call getunit(iunit)

          filename=trim(refdir)//trim(tctmp)//"_"// &
            trim(mesh%name)//"_"//trim(atime)//"uv_ed.dat"
          open(iunit,file=filename, status='old',FORM='UNFORMATTED')
          print*, trim(filename)
          do l=1, mesh%ne

            read(iunit) utmp, vtmp

            call convert_vec_sph2cart(utmp, vtmp, mesh%ed(l)%c%p, vectmp)
            u%f(l)=dot_product(vectmp,mesh%ed(l)%tg)
          end do
        elseif(useStagHC)then
          call getunit(iunit)
          filename=trim(refdir)//trim(tctmp)//"_"// &
            trim(mesh%name)//"_"//trim(atime)//"uv_edhx.dat"
          open(iunit,file=filename, status='old',FORM='UNFORMATTED')
          print*, trim(filename)
          do l=1, mesh%ne
            read(iunit) utmp, vtmp
            call convert_vec_sph2cart(utmp, vtmp, mesh%edhx(l)%c%p, vectmp)
            u%f(l)=dot_product(vectmp,mesh%edhx(l)%nr)
          end do
        end if


        print*, "Reading reference solution for day 6"
        !Read delayed GLW test case (on day 6)
        tctmp="GLW"
        write(atime,'(I10.10)') 518400 !6 days

        !Read thickness
        filename=trim(refdir)//trim(tctmp)//"_"// &
          trim(mesh%name)//"_"//trim(atime)//"h.dat"
        print*, trim(filename)
        inquire(file=filename, exist=ifile)
        if(ifile)then
          call getunit(iunit)
          open(iunit,file=filename, status='old',FORM='UNFORMATTED')
          bt%f=0
          do j=1,mesh%nv
            read(iunit) h_exact%f(j)
            h_exact%f(j)=h_exact%f(j)-bt%f(j)
             !print*, j, h_exact%f(j)
          end do
        else
          print*, "Cannot find reference file:", trim(filename)

        end if

        !Read velocities depending on the grid choice
        if(useStagHTC)then
          call getunit(iunit)

          filename=trim(refdir)//trim(tctmp)//"_"// &
            trim(mesh%name)//"_"//trim(atime)//"uv_ed.dat"
          open(iunit,file=filename, status='old',FORM='UNFORMATTED')
          print*, trim(filename)
          do l=1, mesh%ne

            read(iunit) utmp, vtmp

            call convert_vec_sph2cart(utmp, vtmp, mesh%ed(l)%c%p, vectmp)
            u_exact%f(l)=dot_product(vectmp,mesh%ed(l)%tg)
          end do
        elseif(useStagHC)then
          call getunit(iunit)
          filename=trim(refdir)//trim(tctmp)//"_"// &
            trim(mesh%name)//"_"//trim(atime)//"uv_edhx.dat"
          open(iunit,file=filename, status='old',FORM='UNFORMATTED')
          print*, trim(filename)
          do l=1, mesh%ne
            read(iunit) utmp, vtmp
            call convert_vec_sph2cart(utmp, vtmp, mesh%edhx(l)%c%p, vectmp)
            u_exact%f(l)=dot_product(vectmp,mesh%edhx(l)%nr)
          end do
        end if

      case(32, 33, 34, 35) !Global Steady State Zonal Geo Flow Hollingsworth instability

        u0=pi2*erad/(12._r8*day2sec)
        h0=2.94e4_r8*gravi


        if(hollgw <eps)then
          print*, "Please add a positive layer thickness in swm.par for the Hollingsworth test case"
          stop
        end if


        ! For tc33, h is contant, u is rotation and the f=0
        if(testcase==33.or.testcase==35)then
          fcte=0.0_r8
          fsphere=2
          !u0=200
          h_ct=((u0**2)/2._r8)*gravi
          eta_ct=(2.*u0/erad)
        elseif(testcase==34)then
       	                                !F-sphere
       	                                !fcte=1.4584e-4_r8
          !fsphere=2
          !u0=200
          !h_ct=((u0**2)/2._r8)*gravi
          !eta_ct=(2.*u0/erad)

          !Normal run
          fsphere=0
          h_ct=(erad*omega*u0+(u0**2)/2._r8)*gravi
          eta_ct=(2.*u0/erad+2.*omega)
        else
          h_ct=(erad*omega*u0+(u0**2)/2._r8)*gravi
          eta_ct=(2.*u0/erad+2.*omega)
        end if

        grad_ct=u0*eta_ct
        print*, "Test case parameters"

        print*, "u0", u0
        print*, "h_ct", h_ct


        do i=1, mesh%nv
          !print*, i
          h%f(i)=hollgw
          !if(testcase==32)then !Only put bottom topograpy in tc32 - geostrophic balanced flow
          bt%f(i)=h0-h_ct*dsin(mesh%v(i)%lat)**2
          if(fsphere==2)then
            bt%f(i)=bt%f(i)-(erad/grav)*fcte*dsin(mesh%v(i)%lat)
          end if
          !elseif(testcase==33)then
          !   bt%f(i)=h0-((u0**2)/2._r8)*gravi*dsin(mesh%v(i)%lat)**2
          !end if
          divuh%f(i)=0._r8
          !print*, h%f(i)
          !utmp=u0*dcos(mesh%v(i)%lat)
          !vtmp=0._r8
          !call convert_vec_sph2cart(utmp, vtmp, mesh%v(i)%p, vectmp)
          !Calculate exact Kinectic energy
          if(test_lterror==1)then
            ke_hx_exact%f(i)=(u0*dcos(mesh%v(i)%lat))**2/2._r8
            utmp=u0*dcos(mesh%v(i)%lat)
            vtmp=0._r8
            call convert_vec_sph2cart(utmp, vtmp, mesh%v(i)%p, vectmp)
            vh_hx_exact%p(i)%v=vectmp*h%f(i)
            q_hx_exact%f(i)=eta_ct*dsin(mesh%v(i)%lat)/h%f(i)
            vhq_hx_exact%p(i)%v=vh_hx_exact%p(i)%v*q_hx_exact%f(i)
             !print*, i, q_hx_exact%f(i)
          end if
        end do
        !h_0=h
        h_exact=h
        maxh=h0
        if(test_lterror==1)then
          divuh_exact=divuh
        end if

        if(useStagHTC)then
          do l=1, mesh%ne
            lon=mesh%ed(l)%c%lon
            lat=mesh%ed(l)%c%lat
            utmp=u0*dcos(lat)
            vtmp=0._r8
            call convert_vec_sph2cart(utmp, vtmp, mesh%ed(l)%c%p, vectmp)
            v_ed%p(l)%v=vectmp
            u%f(l)=dot_product(vectmp,mesh%ed(l)%tg)
            if(test_lterror==1)then
              h_ed_exact%f(l)=hollgw
              call convert_vec_sph2cart(0._r8, &
                -grad_ct*dsin(lat)*dcos(lat), mesh%ed(l)%c%p, vectmp)

              grad_ghbK_exact%f(l)=dot_product(vectmp,mesh%ed(l)%tg)
              !print*, eta_ct*dsin(mesh%ed(l)%c%lat), h_ed_exact%f(l)
              q_ed_exact%f(l)=eta_ct*dsin(lat)/h_ed_exact%f(l)
              uhq_perp_exact%f(l)=-eta_ct*dsin(lat)*dot_product(v_ed%p(l)%v,mesh%ed(l)%nr)
              !Laplacian
              utmp=u0*(1.-2.*dcos(lat)**2)/(erad**2*dcos(lat))
              vtmp=0._r8
              call convert_vec_sph2cart(utmp, vtmp, mesh%ed(l)%c%p, vectmp)
              lapu_exact%f(l)=dot_product(vectmp,mesh%ed(l)%tg)
            end if
          end do
        elseif(useStagHC)then
          do l=1, mesh%ne
            lon=mesh%edhx(l)%c%lon
            lat=mesh%edhx(l)%c%lat
            utmp=u0*dcos(lat)
            vtmp=0._r8
            call convert_vec_sph2cart(utmp, vtmp, mesh%edhx(l)%c%p, vectmp)
            v_ed%p(l)%v=vectmp
            u%f(l)=dot_product(vectmp,mesh%edhx(l)%nr)
            if(test_lterror==1)then
              h_ed_exact%f(l)=hollgw
              call convert_vec_sph2cart(0._r8, &
                -grad_ct*dsin(lat)*dcos(lat), mesh%edhx(l)%c%p, vectmp)
              grad_ghbK_exact%f(l)=dot_product(vectmp,mesh%edhx(l)%nr)
              q_ed_exact%f(l)=eta_ct*dsin(lat)/h_ed_exact%f(l)
              uhq_perp_exact%f(l)=+eta_ct*dsin(lat)*dot_product(v_ed%p(l)%v,mesh%edhx(l)%tg)
              !Laplacian
              utmp=u0*(1.-2.*dcos(lat)**2)/(erad**2*dcos(lat))
              vtmp=0._r8
              call convert_vec_sph2cart(utmp, vtmp, mesh%edhx(l)%c%p, vectmp)
              lapu_exact%f(l)=dot_product(vectmp,mesh%edhx(l)%nr)
            end if
          end do
        end if
        !u_0=u
        u_exact=u
        v_ed_exact=v_ed
        maxvel=u0

        if(test_lterror==1)then
          !Loop over triangles
          do k=1, mesh%nt
            !Kin energy
            ke_tr_exact%f(k)=(u0*dcos(mesh%tr(k)%c%lat))**2/2._r8
            !Absolute vorticity
            eta_exact%f(k)=eta_ct*dsin(mesh%tr(k)%c%lat)
            h_tr_exact%f(k)=hollgw
            q_tr_exact%f(k)=eta_exact%f(k)/h_tr_exact%f(k)
            !print*, k, q_tr_exact%f(k)
            utmp=u0*dcos(mesh%tr(k)%c%lat)
            vtmp=0._r8
            call convert_vec_sph2cart(utmp, vtmp, mesh%tr(k)%c%p, vectmp)
            vhq_tr_exact%p(k)%v=vectmp*eta_exact%f(k)

          end do
          masseq_exact%f=0._r8
          momeq_exact%f=0._r8
        end if

      case(40) !Steady state localised test

        h0=10000*gravi !2.94e4_r8*gravi
        k=120 ! good for glevel 5
        n=k*2+2
        u0=sqrt(n*h0*grav) !pi2*erad/(12._r8*day2sec)
        fsphere=2
        fcte=0_r8
        print*, "h0:", h0, "k:", k
        !print*, "hiiii"
        !Set height field
        do i=1,mesh%nv
          !print*, i
          h%f(i)=h0*(2-dsin(mesh%v(i)%lat)**n)
          bt%f(i)=0.
          divuh%f(i)=0._r8

          !Calculate exact Kinectic energy
          if(test_lterror==1)then
            ke_hx_exact%f(i)=(u0*dsin(mesh%v(i)%lat)**k*dcos(mesh%v(i)%lat))**2/2._r8
            utmp=u0*dsin(mesh%v(i)%lat)**k*dcos(mesh%v(i)%lat)
            vtmp=0._r8
            call convert_vec_sph2cart(utmp, vtmp, mesh%v(i)%p, vectmp)
            vh_hx_exact%p(i)%v=vectmp*h%f(i)
            !q_hx_exact%f(i)=eta_ct*dsin(mesh%v(i)%lat)/h%f(i)
            vhq_hx_exact%p(i)%v=vh_hx_exact%p(i)%v*q_hx_exact%f(i)
             !print*, i, q_hx_exact%f(i)
          end if
        end do

        !Set velocity field
        do l=1,mesh%ne
          utmp=0._r8
          vtmp=0._r8
          !print*, l
          if(useStagHC)then
            lat = mesh%edhx(l)%c%lat
            lon = mesh%edhx(l)%c%lon
            utmp=u0*dsin(lat)**k*dcos(lat)
            vtmp=0.
            ! print*, lat, utmp
            !end if
            call convert_vec_sph2cart(utmp, vtmp, mesh%edhx(l)%c%p, vectmp)
            v_ed%p(l)%v=vectmp
            u%f(l)=dot_product(vectmp,mesh%edhx(l)%nr)
          elseif(useStagHTC)then
            lat = mesh%ed(l)%c%lat
            lon = mesh%ed(l)%c%lon
            utmp=u0*dsin(lat)**k*dcos(lat)
            vtmp=0.
            call convert_vec_sph2cart(utmp, vtmp, mesh%ed(l)%c%p, vectmp)
            v_ed%p(l)%v=vectmp
            u%f(l)=dot_product(vectmp,mesh%ed(l)%tg)
          end if
        end do

        h_exact=h
        u_exact=u

      case(41) !Rotated Steady state localised test

        h0=100000*gravi !2.94e4_r8*gravi
        k=120 ! good for glevel 5
        n=k*2+2
        u0=sqrt(n*h0*grav) !pi2*erad/(12._r8*day2sec)
        fsphere=2
        fcte=0_r8

        ! Set point of peak for blob (north hemisphere)
        lon0=1*deg2rad
        lat0=3*deg2rad
        call sph2cart(lon0, lat0, p0(1), p0(2), p0(3))
        print*, "Peak center (lon, lat):"
        print*, lon0*rad2deg, lat0*rad2deg
        print*
        ! Calculate rotation matrix
        Rmat(1:3,1:3)=0._r8
        Rmat(1,1)=dsin(lat0)*dcos(lon0)
        Rmat(1,2)=dsin(lat0)*dsin(lon0)
        Rmat(1,3)=-dcos(lat0)
        Rmat(2,1)=-dsin(lon0)
        Rmat(2,2)=dcos(lon0)
        Rmat(3,1)=dcos(lat0)*dcos(lon0)
        Rmat(3,2)=dcos(lat0)*dsin(lon0)
        Rmat(3,3)=dsin(lat0)

        !Inverse rotation
        RmatT=transpose(Rmat)

        !print*, p0
        !print*, Rmat
        !print*
        !  lon=0*deg2rad
        !  lat=0*deg2rad
        !  call sph2cart(lon, lat, p0(1), p0(2), p0(3))
        !print*, matmul(Rmat, p0)

        !Set height field
        do i=1,mesh%nv
          !print*, i
          p=matmul(Rmat, mesh%v(i)%p)
          call cart2sph(p(1), p(2), p(3), lon, lat)
          h%f(i)=h0*(2-dsin(lat)**n)
          bt%f(i)=0.
          divuh%f(i)=0._r8

          !Calculate exact Kinectic energy
          if(test_lterror==1)then
            ke_hx_exact%f(i)=(u0*dsin(lat)**k*dcos(lat))**2/2._r8
            utmp=u0*dsin(lat)**k*dcos(lat)
            vtmp=0._r8
            call convert_vec_sph2cart(utmp, vtmp, p, vectmp)
            vh_hx_exact%p(i)%v=matmul(RmatT, vectmp)*h%f(i)
             !q_hx_exact%f(i)=eta_ct*dsin(mesh%v(i)%lat)/h%f(i)
             !vhq_hx_exact%p(i)%v=vh_hx_exact%p(i)%v*q_hx_exact%f(i)
             !print*, i, q_hx_exact%f(i)
          end if
        end do

        !Set velocity field
        do l=1,mesh%ne
          utmp=0._r8
          vtmp=0._r8
          !print*, l
          if(useStagHC)then
            p=matmul(Rmat, mesh%edhx(l)%c%p)
            call cart2sph(p(1), p(2), p(3), lon, lat)
            utmp=u0*dsin(lat)**k*dcos(lat)
            vtmp=0.
            ! print*, lat, utmp
            !end if
            call convert_vec_sph2cart(utmp, vtmp, p, vectmp)
            v_ed%p(l)%v=matmul(RmatT, vectmp)
            u%f(l)=dot_product(v_ed%p(l)%v,mesh%edhx(l)%nr)
          elseif(useStagHTC)then
            p=matmul(Rmat, mesh%ed(l)%c%p)
            call cart2sph(p(1), p(2), p(3), lon, lat)
            utmp=u0*dsin(lat)**k*dcos(lat)
            vtmp=0.
            call convert_vec_sph2cart(utmp, vtmp, p, vectmp)
            v_ed%p(l)%v=matmul(RmatT, vectmp)
            u%f(l)=dot_product(v_ed%p(l)%v,mesh%ed(l)%tg)
          end if
        end do

        h_exact=h
        u_exact=u

      case(42) !Rotated Steady state localised test on f-sphere

        h0=100000*gravi !2.94e4_r8*gravi
        k=160 ! good for glevel 5
        n=k*2+2
        u0=sqrt(n*h0*grav) !pi2*erad/(12._r8*day2sec)
        fsphere=2
        fcte=1.4584e-4_r8
        print*, "f: ", fcte
        ! Set point of peak for blob (north hemisphere)
        lon0=1*deg2rad
        lat0=3*deg2rad
        call sph2cart(lon0, lat0, p0(1), p0(2), p0(3))
        print*, "Peak center (lon, lat):"
        print*, lon0*rad2deg, lat0*rad2deg
        print*
        ! Calculate rotation matrix
        Rmat(1:3,1:3)=0._r8
        Rmat(1,1)=dsin(lat0)*dcos(lon0)
        Rmat(1,2)=dsin(lat0)*dsin(lon0)
        Rmat(1,3)=-dcos(lat0)
        Rmat(2,1)=-dsin(lon0)
        Rmat(2,2)=dcos(lon0)
        Rmat(3,1)=dcos(lat0)*dcos(lon0)
        Rmat(3,2)=dcos(lat0)*dsin(lon0)
        Rmat(3,3)=dsin(lat0)

        !Inverse rotation
        RmatT=transpose(Rmat)

        !print*, p0
        !print*, Rmat
        !print*
        !  lon=0*deg2rad
        !  lat=0*deg2rad
        !  call sph2cart(lon, lat, p0(1), p0(2), p0(3))
        !print*, matmul(Rmat, p0)

        !Set height field
        do i=1,mesh%nv
          !print*, i
          p=matmul(Rmat, mesh%v(i)%p)
          call cart2sph(p(1), p(2), p(3), lon, lat)
          if(lat>0)then
            h%f(i)=h0*(2-dsin(lat)**n)
          else
            h%f(i)=2.*h0
          end if

          bt%f(i)=0.
          divuh%f(i)=0._r8

           !          !Calculate exact Kinectic energy
           !          if(test_lterror==1)then
           !             Kin_energy_exact%f(i)=(u0*dsin(lat)**k*dcos(lat))**2/2._r8
           !             utmp=u0*dsin(lat)**k*dcos(lat)
           !             vtmp=0._r8
           !             call convert_vec_sph2cart(utmp, vtmp, p, vectmp)
           !             vh_hx_exact%p(i)%v=matmul(RmatT, vectmp)*h%f(i)
           !             !q_hx_exact%f(i)=eta_ct*dsin(mesh%v(i)%lat)/h%f(i)
           !             !vhq_hx_exact%p(i)%v=vh_hx_exact%p(i)%v*q_hx_exact%f(i)
           !             !print*, i, q_hx_exact%f(i)
           !          end if
        end do

        !Set velocity field
        do l=1,mesh%ne
          utmp=0._r8
          vtmp=0._r8
          !print*, l
          if(useStagHC)then
            p=matmul(Rmat, mesh%edhx(l)%c%p)
            call cart2sph(p(1), p(2), p(3), lon, lat)
            if(lat>eps)then
              F=erad*fcte*dcos(lat)/dsin(lat)
            else
              F=0
            end if
            C=-grav*h0*n*dsin(lat)**(n-2)*dcos(lat)**2
            if(lat>0)then
              utmp=(-F+dsqrt(F**2-4*C))/2
            else
              utmp=0
            end if
            vtmp=0.
            ! print*, lat, utmp
            !end if
            call convert_vec_sph2cart(utmp, vtmp, p, vectmp)
            v_ed%p(l)%v=matmul(RmatT, vectmp)
            u%f(l)=dot_product(v_ed%p(l)%v,mesh%edhx(l)%nr)
          elseif(useStagHTC)then
            p=matmul(Rmat, mesh%ed(l)%c%p)
            call cart2sph(p(1), p(2), p(3), lon, lat)
            if(abs(lat)>eps)then
              F=erad*fcte*dcos(lat)/dsin(lat)
            else
              F=0
            end if
            C=-grav*h0*n*dsin(lat)**(n-2)*dcos(lat)**2
            if(lat>0)then
              utmp=(-F+dsqrt(F**2-4*C))/2
            else
              utmp=0
            end if
            vtmp=0.
            call convert_vec_sph2cart(utmp, vtmp, p, vectmp)
            v_ed%p(l)%v=matmul(RmatT, vectmp)
            u%f(l)=dot_product(v_ed%p(l)%v,mesh%ed(l)%tg)
          end if
        end do

        h_exact=h
        u_exact=u

      case(43) !Play toy

        u00 = 20.0
        lat0 = 0._r8 !pi/7.0
        !lat1 = 0._r8 !pi/7.0 ! - lat0

        en = 1. !exp(-4/(lat1 - lat0)**2)
        umen = -u00/en

        !Set velocity field
        do l=1,mesh%ne
          utmp=0._r8
          vtmp=0._r8
          if(useStagHC)then
            lat = mesh%edhx(l)%c%lat
            den = (lat - lat0)**2 !(lat - lat1)
            !if (den .lt. 0.0D0) then
            utmp = umen*exp(-den*50._r8)
            ! print*, lat, utmp
            !end if
            call convert_vec_sph2cart(utmp, vtmp, mesh%edhx(l)%c%p, vectmp)
            v_ed%p(l)%v=vectmp
            u%f(l)=dot_product(vectmp,mesh%edhx(l)%nr)
          elseif(useStagHTC)then
            lat = mesh%ed(l)%c%lat
            den = (lat - lat0)**2 !(lat - lat1)
            !if (den .lt. 0.0D0) then
            utmp = umen*exp(-den*50._r8)
            ! print*, lat, utmp
            !end if
            call convert_vec_sph2cart(utmp, vtmp, mesh%ed(l)%c%p, vectmp)
            v_ed%p(l)%v=vectmp
            u%f(l)=dot_product(vectmp,mesh%ed(l)%tg)
          end if
        end do

        h%f=100.
        h_exact=h
        u_exact=u

      case default
        print*, "SWM_initialize_fields error - please select a proper test case:", testcase
        stop
    end select

    !Check for CFL constraints
    cfl=abs(maxvel*dt/(mesh%minvdist*erad))
    if(testcase/=11)then
      print*, "CFL:", cfl
      if(cfl>2)then
        print*, "CFL too large, problems may occur"
      end if
    end if
  contains

    function h_tc6(lon, lat)
      real(r8):: h_tc6
      real(r8)::lon
      real(r8):: lat

      if(dabs(lat)>pio2-eps)then
        atmp=0. !avoid division by zero
      else
        atmp=0.5*w*(2.*omega+w)*dcos(lat)**2+&
          0.25*K_RH**2*dcos(lat)**(m*2)* &
          ((m+1)*dcos(lat)**2+(2.*m**2-m-2.)-2.*m**2/(dcos(lat)**2))
      end if

      btmp=2.*((omega+w)*K_RH/((m+1)*(m+2)))* &
        dcos(lat)**(m)*((m**2+2.*m+2.)-(m+1)**2*dcos(lat)**2)

      ctmp=0.25*(K_RH)**2*dcos(lat)**(2*m)*((m+1)*dcos(lat)**2-(m+2))

      h_tc6=h0+gravi*(erad**2)*(atmp + btmp*dcos(m*lon) + ctmp*dcos(2.*m*lon))

      return
    end function h_tc6
  end subroutine initialize_fields

  subroutine tendency(h, u, masseq, momeq)
    !--------------------------------------
    !Calculates the Right Hand Side (spatial discret./tendency)
    !   of mass and velocity equations
    !-------------------------------------------------------------

    !Fluid thickness (defined on voronoi centers)
    type(scalar_field), intent(in):: h  !General

    !Velocities (defined on edges - only normal component)
    type(scalar_field), intent(in):: u  !General

    !Time
    !real(r8), intent(in):: t, dt

    !Right hand side of mass equation (number of cell equations)
    real(r8), intent(inout)::masseq(:)

    !Right hand side of momentum equation (number of edge equations)
    real(r8), intent(inout)::momeq(:)

    !call eqs_spatial_disc(0.0_r8, 0.0_r8, h, u, masseq, momeq)
    !return

    !Initialize RHS (in paralel)
    call zero_vector(momeq)
    call zero_vector(masseq)

    !---------------------------------------------------------------
    !Interpolate thickness to edges and calculate flux at edges
    !---------------------------------------------------------------
    call scalar_hx2ed(h, h_ed, mesh)      !h: cell->edge
    call scalar_elem_product(u, h_ed, uh) !Flux uh at edges

    !---------------------------------------------------------------
    !Calculate vorticity/PV related fields
    !---------------------------------------------------------------

    !Absolute vorticity
    call vort_tr(u, eta, mesh) !Vorticity at triangles

    !Depth at triangles/Rhombi
    if(useSinterpolGass)then !Gass
      call scalar_edtr2trcc(h_ed, h_tr, mesh) !h: ed->tr
      call scalar_tr2ed(h_tr, h_rhb, mesh)    !h: tr->rhombi
    else !TRSK or bary/PXT
      call scalar_hx2tr(h, h_tr, mesh) !h: cell->tr
    end if

    !Interpolate vort from tr to edges (or rhombi in case of Gass)
    call scalar_tr2ed(eta, eta_ed, mesh) ! vort:tr->ed/rhombi

    !Calculate the potential vorticity at triangles
    if(testcase<2)then
      q_tr%f=0.
    else
      call scalar_elem_divide(eta, h_tr, q_tr) !PV on triangles
    end if

    !Calculate vort/PV on edges/rhombi
    if(useSinterpolGass)then
      call scalar_elem_divide(eta_ed, h_rhb, q_ed) !PV on rhombi
    else
      call scalar_tr2ed(q_tr, q_ed, mesh) !PV: tr->ed
    end if

    ! Calculate Grad of PV, if needed
    if((.not.useOrigPV).or. test_lterror==1)then
      call grad_tr2ed(q_tr, q_grad_ed, mesh)
    end if

    !Reconstruct grad of PV to tr center Perot style, if needed
    if(useAPVM.or.useCLUST)then
      !TODO: check if this is correct (/erad) and matches 2016 paper
      ! Fig 18 and 19 of paper do not seem correct)
      call vector_edtr2tr_perot(q_grad_ed, q_grad_tr, mesh)
    end if


    !--------------------------------------------------------
    ! Vector reconstructions
    !--------------------------------------------------------

    !Reconstruct vhq to triangle circumcenter (use Perot's method)
    !  Only used in dtred Coriolis method
    if(useCoriolisMtdDtred .or. test_lterror==1)then
      call vector_edhx2tr_perot(uh, vhq_tr, mesh) !vhq is actualy only vh here
      call vector_elem_product(vhq_tr, q_tr, vhq_tr) !now q is multiplied into vhq
    endif

    ! Velocity reconstruction use Perot's method
    call vector_edhx2hx_perot(u, v_hx, mesh)

    !---------------------------------------------------------------
    ! Calculate divergence / mas eq RHS
    !---------------------------------------------------------------
    call div_hx(uh, divuh, mesh)
    !Mass eq. RHS
    masseq = -divuh%f !-divuh_exact%f !

    !---------------------------------------------------------------
    ! Kinetic energy evaluations
    !---------------------------------------------------------------

    !Calculate K energy at triangles (for Gassman schemes)
    if(useReconmtdGass)then
      call kinetic_energy_tr(u, ke_tr, mesh)
    end if

    !Kinetic energy calculation
    call kinetic_energy_hx(u, v_hx, ke_tr, uh, h, ke_hx, mesh)

    !---------------------------------------------------------------
    ! Grad of bernoulli potential
    !---------------------------------------------------------------

    !Add bottom topography to depth
    hbt%f=h%f+bt%f

    !Calculate g(h+b)+K at cell centers (Bernoulli potential)
    ghbK%f=grav*hbt%f+ke_hx%f

    !Grad g(h+b)+K
    call grad_ed(ghbk, grad_ghbK, mesh)

    !Update momentum eq tendency
    momeq= -grad_ghbK%f !-grad_ghbK_exact%f(l)

    !---------------------------------------------------------------
    ! Calculate PV at cells
    !---------------------------------------------------------------
    if(test_lterror==1)then
      !Calculate pv on cells - just for analysis
      call vector_elem_product(v_hx, h, vh_hx)
      call scalar_tr2hx(q_tr, q_hx, mesh)
      call vector_elem_product(vh_hx, q_hx, vhq_hx)
    end if

    !--------------------------------------------------------------
    !   APVM PV - Filtering of PV noise in triangles
    !---------------------------------------------------------------
    if(useAPVM)then
      call apvm(q_grad_tr, v_hx, pvspar, dt, q_ed, mesh)
    end if

    !---------------------------------------------------------------
    !Calculate coriolis term
    !---------------------------------------------------------------
    call coriolis_ed(u, uh, eta_ed, q_ed, vhq_tr, uhq_perp, mesh)

    !Update mom eq tendency
    momeq=momeq - uhq_perp%f  !-uhq_perp_exact%f(l) !

    !Calculate divergence of velocity - used in difusion and for diagnostics
    call div_hx(u, divu, mesh)

    !------------------------------------------
    !Diffusion
    !------------------------------------------
    if(difus>0)then !Needs verification!!!
      call laplacian_hx(u, eta, lapu, mesh)
      momeq=momeq + difus*lapu%f
    end if

    if(testcase<=1)then
      call zero_vector(momeq)
    end if

    return

  end subroutine tendency

  subroutine ode_rk4 ( t, h, u, h_new, u_new, dt)
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

    !Time and time step
    real(r8):: t
    real(r8):: dt

    !Updated fields
    !Fluid thickness (defined on voronoi centers)
    type(scalar_field):: h_new  !General

    !Velocities (defined on edges - only normal component)
    type(scalar_field):: u_new  !General

    !Times
    real(r8):: t0
    real(r8):: t1
    real(r8):: t2
    real(r8):: t3

    u_new=u
    h_new=h
    masseq%f=massf0
    momeq%f=momf0

    !Initial f (f0)
    t0=t-dt
    !call eqs_spatial_disc(t0, dt, h, u, massf0, momf0)
    call tendency(h, u, massf0, momf0)
    !print*, "Mass, mom: ", maxval(massf0), maxval(momf0)
    !print*, "h, u: ", maxval(h_new%f), maxval(u_new%f)

    !First RK step
    t1 = t0 + dt/2._r8
    u_new%f(1:u%n) = u%f(1:u%n) + dt * momf0(1:u%n) / 2.0_r8
    h_new%f(1:h%n) = h%f(1:h%n) + dt * massf0(1:h%n) / 2.0_r8
    !call eqs_spatial_disc(t1, dt, h_new, u_new, massf1, momf1)
    call tendency(h_new, u_new, massf1, momf1)
    !print*, "Mass, mom: ", maxval(massf1), maxval(momf1)
    !print*, "h, u: ", maxval(h_new%f), maxval(u_new%f)

    !Second RK step
    t2 = t0 + dt/2._r8
    u_new%f(1:u%n) = u%f(1:u%n) + dt * momf1(1:u%n) / 2.0_r8
    h_new%f(1:h%n) = h%f(1:h%n) + dt * massf1(1:h%n) / 2.0_r8
    !call eqs_spatial_disc(t2, dt, h_new, u_new, massf2, momf2)
    call tendency(h_new, u_new, massf2, momf2)
    !print*, "Mass, mom: ", maxval(massf2), maxval(momf2)
    !print*, "h, u: ", maxval(h_new%f), maxval(u_new%f)

    !Third  RK step
    t3 = t0 + dt
    u_new%f(1:u%n) = u%f(1:u%n) + dt * momf2(1:u%n)
    h_new%f(1:h%n) = h%f(1:h%n) + dt * massf2(1:h%n)
    !call eqs_spatial_disc(t3, dt, h_new, u_new, massf3, momf3)
    call tendency(h_new, u_new, massf3, momf3)
    !print*, "Mass, mom: ", maxval(massf3), maxval(momf3)
    !print*, "h, u: ", maxval(h_new%f), maxval(u_new%f)

    !
    ! Combine them to estimate the solution at time t+dt
    !
    u_new%f(1:u%n) = u%f(1:u%n) + dt * (momf0(1:u%n)+2._r8*momf1(1:u%n) &
      +2._r8*momf2(1:u%n)+momf3(1:u%n))/6._r8
    h_new%f(1:h%n) = h%f(1:h%n) + dt * (massf0(1:h%n)+2._r8*massf1(1:h%n) &
      +2._r8*massf2(1:h%n)+massf3(1:h%n))/6._r8

    !print*, "h, u: ", maxval(h_new%f), maxval(u_new%f)

    return
  end subroutine ode_rk4


  subroutine calc_energies(Penergy, Kenergy, Tenergy, Availenergy)
    !Calculate Kinetic, Potential and Total energies

    integer(i4):: i
    integer(i4):: l
    real(r8):: Tenergy
    real(r8):: Kenergy
    real(r8):: Kenergy_new
    real(r8):: Penergy
    real(r8):: Availenergy
    real(r8):: UnavailPenergy
    real(r8):: mean_height
    real(r8):: mean_height_bt

    !Calculate energies
    Penergy=0.
    do i=1, mesh%nv
      Penergy=Penergy+mesh%hx(i)%areag*grav*h%f(i)*(h%f(i)*0.5+bt%f(i))
    end do
    !print*, "Potential:",Penergy

    UnavailPenergy=0.
    !Calculate mean height
    mean_height=sumf_areas(hbt)
    !print*, mean_height, maxval(h%f(1:mesh%nv)), mean_height_bt, maxval(bt%f(1:mesh%nv))
    do i=1, mesh%nv
      mean_height_bt=mean_height-bt%f(i)

      !Unavailable energy
      UnavailPenergy=UnavailPenergy+mesh%hx(i)%areag*grav*mean_height_bt*(mean_height_bt*0.5+bt%f(i))
    end do
    !print*, "Unavailable:", UnavailPenergy

    !print*, "Potential:",Penergy

    UnavailPenergy=0.
    !Calculate mean height
    mean_height=sumf_areas(hbt)
    !print*, mean_height, maxval(h%f(1:mesh%nv)), mean_height_bt, maxval(bt%f(1:mesh%nv))
    do i=1, mesh%nv
      mean_height_bt=mean_height-bt%f(i)

      !Unavailable energy
      UnavailPenergy=UnavailPenergy+mesh%hx(i)%areag*grav*mean_height_bt*(mean_height_bt*0.5+bt%f(i))
    end do
    !print*, "Unavailable:", UnavailPenergy

    Kenergy=0.
    if(useReconmtdPerhx)then
      do i=1, mesh%nv
        Kenergy=Kenergy+ke_hx%f(i)*mesh%hx(i)%areag
      end do
    else       
      do l=1, mesh%ne
        Kenergy=Kenergy+0.5*mesh%ed(l)%leng*mesh%edhx(l)%leng*h_ed%f(l)*u%f(l)**2
      end do
    end if

    Tenergy=Penergy+Kenergy
    Availenergy=Tenergy-UnavailPenergy
    !print*, "Available:", Availenergy

    Availenergy=Tenergy-UnavailPenergy
    !print*, "Available:", Availenergy

    return
  end subroutine calc_energies

  subroutine plotfields(k, time)
    !-------------------------------------------
    !  Plot fields
    !  k- index for time couting
    !  time - current time step
    !-------------------------------------------

    integer(i4)::k !Time index
    real(r8)::time
    character (len=60)::  atime

    write(atime,*) nint(time)

    if(.not.plots)return

    !Plot Fields
    !print*, k, ploterrors, useRefSol, (k==ntime  .or. mod(k,plotsteps)==0 ).or. (ploterrors .and. useRefSol)
    !if( (k==ntime  .or. mod(k,plotsteps)==0 ).or. (ploterrors .and. useRefSol))then
    if( (k==ntime  .or. mod(k,plotsteps)==0 ).or. (ploterrors .and. useRefSol .and. mod(k,plotsteps)==0))then

      if(test_lterror/=1)then
        !Scalar field plots
        h%name=trim(swmname)//"_h_t"//trim(adjustl(trim(atime)))
        call plot_scalarfield(h, mesh)

        if(maxval(bt%f(1:bt%n)) > eps)then
          hbt%name=trim(swmname)//"_hbt_t"//trim(adjustl(trim(atime)))
          call plot_scalarfield(hbt, mesh)
        end if

        eta%name=trim(swmname)//"_eta_t"//trim(adjustl(trim(atime)))
        call plot_scalarfield(eta, mesh)

        ke_hx%name=trim(swmname)//"_Kenergy_t"//trim(adjustl(trim(atime)))
        call plot_scalarfield(ke_hx, mesh)

        q_tr%name=trim(swmname)//"_pv_t"//trim(adjustl(trim(atime)))
        call plot_scalarfield(q_tr, mesh)

        q_ed%name=trim(swmname)//"_pv_ed_t"//trim(adjustl(trim(atime)))
        call plot_scalarfield(q_ed, mesh)

        divuh%name=trim(swmname)//"_divuh_t"//trim(adjustl(trim(atime)))
        call plot_scalarfield(divuh, mesh)

        if(useRefSol)then !Calculate only u and h errors
          h_error%f=h%f-h_exact%f
          h_error%name=trim(swmname)//"_h_error_t"//trim(adjustl(trim(atime)))
          call plot_scalarfield(h_error, mesh)

          u_error%f=u%f-u_exact%f
          u_error%name=trim(swmname)//"_u_error_t"//trim(adjustl(trim(atime)))
          call plot_scalarfield(u_error, mesh)

          q_tr_error%f=q_tr%f-q_tr_exact%f
          q_tr_error%name=trim(swmname)//"_q_tr_error_t"//trim(adjustl(trim(atime)))
          call plot_scalarfield(q_tr_error, mesh)
        end if

      elseif(test_lterror==1)then !Calculate all Errors

        !Scalar field plots
        h%name=trim(swmname)//"_h_t"//trim(adjustl(trim(atime)))
        call plot_scalarfield(h, mesh)

        if(maxval(bt%f(1:bt%n)) > eps)then
          hbt%name=trim(swmname)//"_hbt_t"//trim(adjustl(trim(atime)))
          call plot_scalarfield(hbt, mesh)
        end if

        eta%name=trim(swmname)//"_eta_t"//trim(adjustl(trim(atime)))
        call plot_scalarfield(eta, mesh)

        ke_hx%name=trim(swmname)//"_Kenergy_t"//trim(adjustl(trim(atime)))
        call plot_scalarfield(ke_hx, mesh)


        v_hx%name=trim(swmname)//"_v_hx_t"//trim(adjustl(trim(atime)))
        call plot_cart_vectorfield(v_hx, mesh)

        v_ed%name=trim(swmname)//"_v_ed_t"//trim(adjustl(trim(atime)))
        call plot_cart_vectorfield(v_ed, mesh)

        q_grad_tr%name=trim(swmname)//"_q_grad_tr_t"//trim(adjustl(trim(atime)))
        call plot_cart_vectorfield(q_grad_tr, mesh)

        q_tr%name=trim(swmname)//"_pv_t"//trim(adjustl(trim(atime)))
        call plot_scalarfield(q_tr, mesh)

        divuh%name=trim(swmname)//"_divuh_t"//trim(adjustl(trim(atime)))
        call plot_scalarfield(divuh, mesh)

        uhq_perp%name=trim(swmname)//"_uhq_perp_t"//trim(adjustl(trim(atime)))
        call plot_scalarfield(uhq_perp, mesh)

        ghbK%name=trim(swmname)//"_ghbK_t"//trim(adjustl(trim(atime)))
        call plot_scalarfield(ghbK, mesh)

        momeq%name=trim(swmname)//"_momeq_t"//trim(adjustl(trim(atime)))
        call plot_scalarfield(momeq, mesh)

        masseq%name=trim(swmname)//"_masseq_t"//trim(adjustl(trim(atime)))
        call plot_scalarfield(masseq, mesh)

        h_error%f=h%f-h_exact%f
        h_error%name=trim(swmname)//"_h_error_t"//trim(adjustl(trim(atime)))
        call plot_scalarfield(h_error, mesh)

        u_error%f=u%f-u_exact%f
        u_error%name=trim(swmname)//"_u_error_t"//trim(adjustl(trim(atime)))
        call plot_scalarfield(u_error, mesh)

        ke_hx_error%f=ke_hx%f-ke_hx_exact%f
        ke_hx_error%name=trim(swmname)//"_Kenergy_error_t"//trim(adjustl(trim(atime)))
        call plot_scalarfield(ke_hx, mesh)

        eta_error%f=eta%f-eta_exact%f
        eta_error%name=trim(swmname)//"_eta_error_t"//trim(adjustl(trim(atime)))
        call plot_scalarfield(eta_error, mesh)

        q_tr_error%f=q_tr%f-q_tr_exact%f
        q_tr_error%name=trim(swmname)//"_pv_error_t"//trim(adjustl(trim(atime)))
        call plot_scalarfield(q_tr_error, mesh)

        q_grad_tr_exact%name=trim(swmname)//"_q_grad_tr_exact_t"//trim(adjustl(trim(atime)))
        call plot_cart_vectorfield(q_grad_tr_exact, mesh)


        divuh_error%f=divuh%f-divuh_exact%f
        divuh_error%name=trim(swmname)//"_divuh_error_t"//trim(adjustl(trim(atime)))
        call plot_scalarfield(divuh_error, mesh)

        uhq_perp_error%f=uhq_perp%f-uhq_perp_exact%f
        uhq_perp_error%name=trim(swmname)//"_uhq_perp_error_t"//trim(adjustl(trim(atime)))
        call plot_scalarfield(uhq_perp_error, mesh)

        grad_ghbK_error%f=grad_ghbK%f-grad_ghbK_exact%f
        grad_ghbK_error%name=trim(swmname)//"_grad_ghbK_t"//trim(adjustl(trim(atime)))
        call plot_scalarfield(grad_ghbK_error, mesh)

      end if

    end if
  end subroutine plotfields

  subroutine plot_errors_uh()
    !-------------------------------------------
    !  Plot error fields for u and h
    !-------------------------------------------

    if(.not.plots)return

    !Scalar field plots
    h_error%name=trim(swmname)//"_h_error"
    call plot_scalarfield(h_error, mesh)

    u_error%name=trim(swmname)//"_u_error"
    call plot_scalarfield(u_error, mesh)

  end subroutine plot_errors_uh

  subroutine write_errors(time, errormax_h, error2_h, errormax_u, error2_u, error_mass, &
    Penergy, Kenergy, Tenergy, tmp)
    !----------------------------------------------------------
    !  write errors to fixed general file for all model parameters
    !----------------------------------------------------------
    !File name for output
    character (len=256):: filename

    !File units
    integer (i4):: iunit
    logical::  ifile

    !Errors
    real(r8):: time
    real(r8):: errormax_h
    real(r8):: errormax_u
    real(r8):: error2_h
    real(r8):: error2_u
    real(r8):: error_mass
    real(r8):: Tenergy
    real(r8):: Kenergy
    real(r8):: Penergy
    real(r8):: tmp

    !Char buffer
    character (len=512):: buffer

    !File for errors
    filename=trim(datadir)//"swm_errors.txt"
    call getunit(iunit)

    buffer="        n        mvdist  tcase       ntime"//&
      "     dt(s)                 cfl                time(dys)  "//&
      "        errormax_h       error2_h               errormax_u     "//&
      "   error2_u               mass              Penergy      "//&
      "     Kenergy             Tenergy       SWMparameters"

    inquire(file=filename, exist=ifile)
    if(ifile)then
      open(iunit,file=filename, status='old', position='append')
    else
      open(iunit,file=filename, status='replace')
      write(iunit, '(a)') buffer
    end if

    !Write errors to file
    write(iunit, '(i12, 1f12.4, i4, i10, 12e20.10 , a120)') &
      mesh%nv,  mesh%meanvdist*rad2deg, testcase, &
      ntime, dt, cfl, time*sec2day, &
      errormax_h, error2_h, errormax_u, error2_u, error_mass, &
      Penergy, Kenergy, Tenergy, tmp, &
      trim(adjustl(trim(swmname)))//"_"//trim(adjustl(trim(mesh%name)))

    print*
    print*, "Final Errors"
    buffer="         n    mvdist  tcase   nt     dt(s) "//&
      "       cfl     time(dys) errormax_h      error2_h     "//&
      "   errormax_u      error2_u        mass             "//&
      " Penergy           Kenergy     Tenergy     tmp"

    write(*,'(a)') trim(buffer)
    write(*, '(i12, 1f8.3, i4, i10, 3f10.4, 9e16.6)') &
      mesh%nv,  mesh%meanvdist*rad2deg, testcase, &
      ntime, dt, cfl, time*sec2day, &
      errormax_h, error2_h, errormax_u, error2_u, error_mass, &
      Penergy, Kenergy, Tenergy, tmp
    print*
  end subroutine write_errors

  subroutine write_evol_error(k, time, tmass, errormax_h, error2_h, &
    errormax_u, error2_u, &
    errormax_pv, error2_pv,&
    Penergy, Kenergy, Tenergy, Availenergy, &
    RMSdiv, maxdiv, max_gradke, tmp)
    !----------------------------------------------------------
    !  write errors to specific file for this specific model set up
    !    at defined time steps
    !----------------------------------------------------------
    !File name for output
    character (len=256):: filename

    !File units
    integer (i4):: iunit
    logical::  iopen

    !Errors
    real(r8):: tmass
    real(r8):: errormax_h
    real(r8):: errormax_u
    real(r8):: error2_h
    real(r8):: error2_u
    real(r8):: errormax_pv
    real(r8):: error2_pv
    !real(r8):: error_mass
    real(r8):: Tenergy
    real(r8):: Kenergy
    real(r8):: Penergy
    real(r8):: Availenergy
    real(r8):: RMSdiv, maxdiv, max_gradke, tmp
    integer(i4):: k
    !Time
    real(r8):: time

    !Char buffer
    character (len=512):: buffer

    !File for errors
    filename=trim(datadir)//trim(swmname)//"_evolution_"//trim(mesh%name)//".txt"
    buffer="        n        mvdist  tcase       k       nt    "//&
      "     dt(s)                 cfl                time(dys)  "//&
      "        errormax_h       error2_h               errormax_u     "//&
      "   error2_u               errormax_pv       error2_pv            "//&
      " mass              Penergy      "//&
      "     Kenergy             Tenergy          Availenergy      "//&
      "     RMSdiv              maxdiv    "//&
      "     max_gradke          tmp        SWMparameters"

    inquire(file=filename, opened=iopen)
    if(iopen)then
      inquire(file=filename, number=iunit)
    else
      call getunit(iunit)
      open(iunit,file=filename, status='replace')
      !buffer="      time     mass       Penergy     Kenergy          Tenergy"
      write(iunit, '(a)') trim(buffer)
    end if



    if(k==0 .or. mod(k,nprints*40)==0 )then
      buffer="     n    mvdist  tcase     k     nt   dt(s) "//&
        "       cfl   time(dys) errormax_h   error2_h  "//&
        "   errormax_u  error2_u    mass       "//&
        " Penergy           Kenergy     Tenergy    Availenergy  errormax_pv   error2_pv  RMSdiv  maxdiv  max_gradke    tmp"
      write(*, '(a)') trim(buffer)
    end if



    !if( k==1 .or. k==ntime  .or. mod(k,nprints)==0 )then
    !Write errors to file
    write(iunit, '(i12, 1f12.4, i4, 2i10, 18e20.10 , a100)') &
      mesh%nv,  mesh%meanvdist*rad2deg, testcase, &
      k, ntime, dt, cfl, time*sec2day, &
      errormax_h, error2_h, errormax_u, error2_u, errormax_pv, error2_pv, tmass, &
      Penergy, Kenergy, Tenergy , Availenergy, RMSdiv, maxdiv, max_gradke, tmp, &
      trim(adjustl(trim(swmname)))//"_"//trim(adjustl(trim(mesh%name)))
    write(*, '(i8, 1f8.3, i4, 2i8, 3f8.3, 14e12.4, 1e16.8)') &
      mesh%nv,  mesh%meanvdist*rad2deg, testcase, &
      k, ntime, dt, cfl, time*sec2day, &
      errormax_h, error2_h, errormax_u, error2_u, tmass, &
      Penergy, Kenergy, Tenergy, Availenergy, errormax_pv, error2_pv, RMSdiv, maxdiv, max_gradke, tmp
    !if(k==0)then
    !Buffer for write_evol output - neat output this way
    !buffer="        n    mvdist  tcase       k       nt    dt(s)  "//&
    ! "    cfl     time(dys) "//&
    ! "         mass "//&
    ! "        Penergy          Kenergy       Tenergy "
    !write(*, '(a)') trim(buffer)
    !end if
    !end if
    !print*

  end subroutine write_evol_error

  subroutine write_evol(k, time, tmass, Penergy, Kenergy, Tenergy, Availenergy, RMSdiv)
    !----------------------------------------------------------
    !  write info to specific file for this specific model set up
    !    at defined time steps
    !----------------------------------------------------------
    !File name for output
    character (len=256):: filename

    !File units
    integer (i4):: iunit
    logical::  iopen

    !Errors
    real(r8):: tmass
    real(r8):: Tenergy
    real(r8):: Kenergy
    real(r8):: Penergy
    real(r8):: Availenergy
    real(r8):: RMSdiv
    integer(i4):: k
    !Time
    real(r8):: time

    !Char buffer
    character (len=256):: buffer

    !File for errors
    filename=trim(datadir)//trim(swmname)//"_mass_energy_"//trim(mesh%name)//".txt"
    buffer="        n        mvdist  tcase       k       nt        dt(s)       "//&
      "          cfl                 time(dys)             mass               Penergy    "//&
      "          Kenergy            Tenergy            Availenergy           RMSdiv"

    inquire(file=filename, opened=iopen)
    if(iopen)then
      inquire(file=filename, number=iunit)
    else
      call getunit(iunit)
      open(iunit,file=filename, status='replace')
      !buffer="      time     mass       Penergy     Kenergy          Tenergy"
      write(iunit, '(a)') trim(buffer)
    end if

    if(k==0 .or. mod(k,nprints*40)==0 )then
      buffer="        n    mvdist  tcase       k       nt    dt(s)  "//&
        "    time(dys) "//&
        "     mass "//&
        "        Penergy          Kenergy       Tenergy        Availenergy       RMSdiv"
      write(*, '(a)') trim(buffer)
    end if

    !if( k==1 .or. k==ntime  .or. mod(k,nprints)==0 )then
    !Write errors to file
    write(iunit, '(i12, 1f12.4, i4, 2i10, 9e28.16)') & !, a60, a60)') &
      mesh%nv,  mesh%meanvdist*rad2deg, testcase, &
      k, ntime, dt, cfl, time*sec2day, tmass, &
      Penergy, Kenergy, Tenergy, Availenergy, RMSdiv !, &
    !trim(adjustl(trim(swmname))), trim(adjustl(trim(mesh%name)))
    write(*, '(i12, 1f8.3, i4, 2i10, 2f10.4, 6e16.6)') &
      mesh%nv,  mesh%meanvdist*rad2deg, testcase, &
      k, ntime, dt, time*sec2day, tmass, &
      Penergy, Kenergy, Tenergy, Availenergy, RMSdiv
    !end if
    !print*

  end subroutine write_evol

  subroutine read_refdata()
    !----------------------------------------------------------
    !  Read reference model data (Spectral SWM)
    !  Reads h, u, v for certains timesteps
    !----------------------------------------------------------
    !File name for output
    character (len=256):: filename

    !File units
    integer (i4):: iunit
    logical::  ifile

    !Time


    !Indexes
    integer(i4):: i
    integer(i4):: j
    integer(i4):: k

    print*, "Reading reference data..."

    filename=trim(refdir)//"graphinfo"
    call getunit(iunit)

    inquire(file=filename, exist=ifile)
    if(ifile)then
      open(iunit,file=filename, status='old')
    else
      print*, "ERROR: Could not find graphinfo in ref/ folder"
      stop
    end if

    read(iunit, *) ntimeref, nlonref, nlatref
    print*, ntimeref, nlonref, nlatref
    allocate(latlevels(1:nlatref))

    read(iunit, *) latlevels(1:nlatref)
    print*, latlevels

    allocate(href(1:nlonref, 1:nlatref, 1:ntimeref))
    allocate(uref(1:nlonref, 1:nlatref, 1:ntimeref))
    allocate(vref(1:nlonref, 1:nlatref, 1:ntimeref))

    filename=trim(refdir)//"graphout"
    call getunit(iunit)

    inquire(file=filename, exist=ifile)
    if(ifile)then
      !open(iunit,file=filename, status='old')
      open(iunit,file=filename,access='sequential',form='unformatted')
    else
      print*, "ERROR: Could not find graphout in ref/ folder"
      stop
    end if

    do k =1,ntimeref
      read(iunit) ((href(i,j, k),i=1,nlonref),j=1,nlatref)
      read(iunit) ((uref(i,j, k),i=1,nlonref),j=1,nlatref)
      read(iunit) ((vref(i,j, k),i=1,nlonref),j=1,nlatref)
    end do

    !print*, href(320, 180, 1), uref(nlonref/2,nlatref/2, 1), vref(nlonref/2,nlatref/2, 1)
    !stop

    !print*
    !print*, "Total mass:", tmass

  end subroutine read_refdata

  subroutine read_refdata_endgame(time)
    !----------------------------------------------------------
    !  Read reference model data from END GAME
    !  Reads h, u, v for certains timesteps
    !----------------------------------------------------------
    !File name for input
    character (len=256):: filename
    character (len=256):: atmp, atime

    !File units
    integer (i4):: iunit
    logical::  ifile

    !Time
    real(r8)::time

    !Indexes
    integer(i4):: i
    integer(i4):: j
    integer(i4):: k
    integer(i4):: l


    !Auxiliar velocity veriables
    real(r8):: utmp
    real(r8):: vtmp
    real(r8):: vectmp(1:3)

    useRefSol=.false.

    !Read reference times (only once)
    if(.not. allocated(reftimes))then

      write(atmp, '(i2.2)') testcase
      if(testcase==21)then
        filename=trim(refdir)//"GLW_"//"reftimes.dat"
      else
        filename=trim(refdir)//"TC"//trim(atmp)//"_"//"reftimes.dat"
      end if
      inquire(file=filename, exist=ifile)

      if(.not.ifile)then
        print*, "There is no reference file to read"
        return
      end if
      print*, "Reading reference data on times given from file:", &
        trim(filename)
      call getunit(iunit)
      open(iunit, FILE=trim(filename))
      read(iunit, *) nreftime
      allocate(reftimes(1:nreftime))
      do i=1, nreftime
        read(iunit, *) reftimes(i)
      end do
       !print*, "Times:", reftimes(1:nreftime)
    end if
    !write(*,*) time, nint(time), nint(reftimes(1:nreftime))
    !Check if it is necessary to read reference data
    do i=1,nreftime
      !if(nint(time) == nint(reftimes(i)).or. .true.)then
      !print*, nint(time), nint(reftimes(i))
      if(nint(time) == nint(reftimes(i)))then

        print*, "Reading reference data from file"
        write(atmp, '(i2.2)') testcase
        write(atime,'(I10.10)') NINT(reftimes(i))
        !Read thickness
        if(testcase==21)then
          filename=trim(refdir)//"GLW_"// &
            trim(mesh%name)//"_"//trim(atime)//"h.dat"
        else
          filename=trim(refdir)//"TC"//trim(atmp)//"_"// &
            trim(mesh%name)//"_"//trim(atime)//"h.dat"
        end if
        print*, trim(filename)
        inquire(file=filename, exist=ifile)
        if(ifile)then
          call getunit(iunit)
          open(iunit,file=filename, status='old',FORM='UNFORMATTED')
          do j=1,mesh%nv
            read(iunit) h_exact%f(j)
            h_exact%f(j)=h_exact%f(j)-bt%f(j)
             !print*, j, h_exact%f(j)
          end do
        else
          print*, "Cannot find reference file:", trim(filename)
          cycle
        end if
        useRefSol=.true.
        !Read velocities depending on the grid choice
        if(useStagHTC)then
          call getunit(iunit)
          if(testcase==21)then
            filename=trim(refdir)//"GLW_"// &
              trim(mesh%name)//"_"//trim(atime)//"uv_ed.dat"
          else
            filename=trim(refdir)//"TC"//trim(atmp)//"_"// &
              trim(mesh%name)//"_"//trim(atime)//"uv_ed.dat"
          end if
          open(iunit,file=filename, status='old',FORM='UNFORMATTED')
          print*, trim(filename)
          do l=1, mesh%ne
            read(iunit) utmp, vtmp
            call convert_vec_sph2cart(utmp, vtmp, mesh%ed(l)%c%p, vectmp)
            u_exact%f(l)=dot_product(vectmp,mesh%ed(l)%tg)
          end do
        elseif(useStagHC)then
          call getunit(iunit)
          if(testcase==21)then
            filename=trim(refdir)//"GLW_"// &
              trim(mesh%name)//"_"//trim(atime)//"uv_edhx.dat"
          else
            filename=trim(refdir)//"TC"//trim(atmp)//"_"// &
              trim(mesh%name)//"_"//trim(atime)//"uv_edhx.dat"
          end if
          open(iunit,file=filename, status='old',FORM='UNFORMATTED')
          print*, trim(filename)
          do l=1, mesh%ne
            read(iunit) utmp, vtmp
            call convert_vec_sph2cart(utmp, vtmp, mesh%edhx(l)%c%p, vectmp)
            u_exact%f(l)=dot_product(vectmp,mesh%edhx(l)%nr)
          end do
        end if
        if(plots)then
          write(atime,*) NINT(reftimes(i))
          h_exact%name=trim(swmname)//"_h_exact_t"//trim(adjustl(trim(atime)))
          call plot_scalarfield(h_exact, mesh)
          u_exact%name=trim(swmname)//"_u_exact_t"//trim(adjustl(trim(atime)))
          call plot_scalarfield(u_exact, mesh)
        end if
      end if

    end do

    !print*
    !print*, "Total mass:", tmass

  end subroutine read_refdata_endgame

  function refdata2icos(time)
    !-------------------------------------------
    !  Convert reference solution form Gaussian grid to icosahedral grid
    !   for a given time (givens in days)
    !-------------------------------------------
    real(r8)::time
    real(r8):: timeref
    integer(i4):: ktime
    integer(i4):: i
    real(r8):: refdata2icos
    !The total running time of the reference solution is 14 days
    !The reference data will have ntimeref timesteps equaly spaced
    !  within 0 and 14 days
    !
    timeref=14_i4/ntimeref !data for every ktimeref days

    ktime=nint(time/timeref)+1

    if(ktime>ntimeref)then
      print*, "Warning: Cannot evaluate reference solution at this time"
      stop
    end if

    !Get h field
    do i=1, mesh%nv
       !lat=mesh%v(i)%lat
    end do
    refdata2icos=0.

  end function refdata2icos

  function sumf_areas(sf)
    !---------------------------------------
    !Calculates Sum f(i)A_i / Sum A_i
    ! where A_i is the area of the cell
    ! it will depend on the position of the scalar field f
    !-------------------------------------------------------

    !Scalar field
    type(scalar_field), intent(in):: sf  !General

    !Sum
    real(r8)::sumf_areas
    real(r8):: sumareas

    !Indexes
    integer(i4)::i

    sumf_areas=0.
    sumareas=0.

    select case(sf%pos)
      case(0) !Scalars on Voronoi cell nodes
        do i=1, sf%n
          sumareas=sumareas+mesh%hx(i)%areag
          sumf_areas=sumf_areas+mesh%hx(i)%areag*sf%f(i)
        end do

      case(1) !Scalars on triangles
        do i=1, sf%n
          sumareas=sumareas+mesh%tr(i)%areag
          sumf_areas=sumf_areas+mesh%tr(i)%areag*sf%f(i)
        end do

      case(3, 6) !Scalars on Voronoi cell edges
        do i=1, sf%n
          sumareas=sumareas+mesh%ed(i)%leng*mesh%edhx(i)%leng*0.5_r8
          sumf_areas=sumf_areas+mesh%hx(i)%areag*sf%f(i)
        end do

      case default
        print*, "Don't know how to calculate sumf_areas on this position:", sf%pos
        stop
    end select

    sumf_areas=sumf_areas/sumareas

  end function sumf_areas

  function sumfsq_areas(sf)
    !---------------------------------------
    !Calculates Sum f(i)**2A_i / Sum A_i
    ! where A_i is the area of the cell
    ! it will depend on the position of the scalar field f
    !-------------------------------------------------------

    !Scalar field
    type(scalar_field), intent(in):: sf  !General

    !Sum
    real(r8):: sumfsq_areas
    real(r8):: sumareas
    real(r8):: area

    !Indexes
    integer(i4)::i

    sumfsq_areas=0.
    sumareas=0.

    select case(sf%pos)
      case(0) !Scalars on Voronoi cell nodes
        do i=1, sf%n
          sumareas=sumareas+mesh%hx(i)%areag
          sumfsq_areas=sumfsq_areas+mesh%hx(i)%areag*sf%f(i)**2
        end do

      case(1) !Scalars on triangles
        do i=1, sf%n
          sumareas=sumareas+mesh%tr(i)%areag
          sumfsq_areas=sumfsq_areas+mesh%tr(i)%areag*sf%f(i)**2
        end do

      case(3, 6) !Scalars on Voronoi cell edges
        do i=1, sf%n
          area=mesh%ed(i)%leng*mesh%edhx(i)%leng*0.5_r8
          sumareas=sumareas+area
          sumfsq_areas=sumfsq_areas+area*sf%f(i)**2
        end do

      case default
        print*, "Don't know how to calculate sumfsq_areas on this position:", sf%pos
        stop
    end select
    !print*, sumareas, unitspharea
    sumfsq_areas=sumfsq_areas/sumareas

  end function sumfsq_areas

  subroutine error_calc(f, fref, fdif, errormaxrel, error2rel, errormax)
    !----------------------------------------------------------
    !Calculate error norms of the field "f" relative to the
    ! reference "fref" using norm "norm"
    ! error2 : L2 norm
    ! errormax : Linfinity norm
    ! (Williamson 1992)
    !
    ! Returns:
    !  - errormax and error2
    !  - fdif (pointwise difference between f and fref)
    !----------------------------------------------------------

    !Scalar fields
    type(scalar_field), intent(in):: f     !Estimated field
    type(scalar_field), intent(in):: fref  !Reference values
    type(scalar_field), intent(inout):: fdif  !Error

    !Error and norm
    real(r8)::errormaxrel
    real(r8):: error2rel
    real(r8):: errormax
    !integer(i4)::norm !norm used

    !Auxiliar vars
    real(r8):: refsum
    real(r8):: refmax

    fdif=f
    fdif%f=f%f-fref%f

    !L2 error
    refsum=sumfsq_areas(fref)
    if(refsum>0)then
      error2rel=dsqrt(sumfsq_areas(fdif))/dsqrt(refsum)
    else
      error2rel=dsqrt(sumfsq_areas(fdif))
    end if

    !Maximum error
    refmax=maxval(abs(fref%f(1:fref%n)))
    if(refmax>0)then
      !print*, "using rel:", trim(f%name), refmax, maxval(fdif%f), maxval(fref%f), maxval(f%f)
      errormax=maxval(abs(fdif%f(1:fdif%n)))
      errormaxrel=errormax/refmax
    else
      !print*, "using abs:", trim(f%name), refmax
      errormax=maxval(abs(fdif%f(1:fdif%n)))
      errormaxrel=errormax
    end if

  end subroutine error_calc

  subroutine error_calc_vec(v, vref, vdif, errormaxrel, error2rel, errormax)
    !----------------------------------------------------------
    !Calculate error norms of the vector field "v" relative to the
    ! reference "vref" using norm "norm"
    ! error2rel : L2 norm
    ! errormaxrel : Linfinity norm
    ! errormax : Maximum difference
    ! (Williamson 1992)
    !
    ! Returns:
    !  - errormax, errormaxrel and error2rel
    !  - vdif (pointwise difference between v and vref)
    !----------------------------------------------------------

    !Vector fields
    type(vector_field_cart), intent(in):: v     !Estimated field
    type(vector_field_cart), intent(in):: vref  !Reference values
    type(vector_field_cart), intent(inout):: vdif  !Error

    !Scalar fields (norms of vectors
    type(scalar_field):: fref  !Reference values
    type(scalar_field):: fdif  !Error

    !Error and norm
    real(r8)::errormax !errors
    real(r8):: error2rel !errors
    real(r8):: errormaxrel !errors
    !integer(i4)::norm !norm used

    !Auxiliar vars
    real(r8):: refsum
    real(r8):: refmax

    integer (i4)::i

    allocate(fref%f(1:v%n))
    allocate(fdif%f(1:v%n))
    fref%n=v%n
    fdif%n=v%n
    fref%pos=v%pos
    fdif%pos=v%pos

    vdif=v
    do i=1, v%n
      vdif%p(i)%v=v%p(i)%v-vref%p(i)%v
      fdif%f(i)=norm(vdif%p(i)%v)
      fref%f(i)=norm(vref%p(i)%v)
    end do

    !L2 error
    refsum=sumfsq_areas(fref)
    if(refsum>0)then
      error2rel=dsqrt(sumfsq_areas(fdif))/dsqrt(refsum)
    else
      error2rel=dsqrt(sumfsq_areas(fdif))
    end if

    !Maximum error
    refmax=maxval(abs(fref%f(1:fref%n)))
    if(refmax>0)then
      errormax=maxval(abs(fdif%f(1:fdif%n)))
      errormaxrel=errormax/refmax
    else
      errormaxrel=errormax
    end if

  end subroutine error_calc_vec

  function sec(x)
    ! Secant calculation
    real(r8):: x
    real(r8)::sec
    if(dabs(dabs(x)-pi/2)<eps/1000.)then
      sec=0.
    else
      sec=1./dcos(x)
    end if
    return
  end function sec


end module swm
!-------------------------------------------------


