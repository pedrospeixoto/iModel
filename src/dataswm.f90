module dataswm
  !=============================================================================
  !  Global data for shallow water model
  !
  ! Pedro da Silva Peixoto (pedrosp@ime.usp.br)
  ! Apr 2015
  !=============================================================================

  !Use global constants and kinds
  use constants

  !Data structures for geodesic grids
  use datastruct, only: &
       general_coords, &
       grid_structure, &
       scalar_field, &
       vector_field_cart

  use smeshpack, only: &
       getunit

  !-----------------------------------
  !Global variables
  !-----------------------------------

  !Grid structure
  type(grid_structure) :: mesh

  !Flags
  integer (i4):: testcase !Test case
  integer (i4):: initialfield !Initial condition

  !Local Truncation Error tests
  integer (i4):: test_lterror !0=false, 1=true

  !Linearized SW equation tests constant f
  integer (i4):: fsphere !0=false, 1=true and linearized eq, 2=true but non linear eq
  real(r8):: fcte !Value of f for f-sphere=1

  !Time variables
  real (r8):: dt  !Time step
  real (r8):: period   !Period
  real (r8):: maxtime   !stopping time
  integer (i4):: ntime !number of time steps
  real(r8):: cfl

  !Difusion coeficient
  real(r8):: difus

  !Hollingsworth coeficient
  real(r8):: hollgw

  ! Potential vorticity correction/stabilization flag
  character (len=6) :: pv_stab
  logical :: useAPVM
  logical :: useCLUST
  logical :: useOrigPV
  real(r8) :: pvspar

  !Use reference solution
  logical::useRefSol
  logical::RefSolRead, RefSolAnal

  !Coriolis Vector reconstruction method
  ! trsk = TRiSK scheme
  ! dtred  = Use dual triangles
  ! pered  = Modified Perot
  ! hyb  = Hybrid Perod and Trisk
  character (len=6)::  coriolis_reconmtd
  logical :: useCoriolisMtdPered !Modified Perot
  logical :: useCoriolisMtdTrisk !Trisk
  logical :: useCoriolisMtdDtred !Dtred
  logical :: useCoriolisMtdHyb   !Hybrid
  logical :: useCoriolisMtdExact !Exact - if possible

  !General Vector reconstruction method
  character (len=6)::  reconmtd
  logical :: useReconmtdPerhx !Perot
  logical :: useReconmtdTrisk !Trisk
  logical :: useReconmtdGass !Gassman's version to avoid Hol. instabil
  logical :: useReconmtdMelv !Melvin's version to avoid Hol. instabil
  real(r8):: gasscoef

  !Staggering type
  ! HC - velocities on Voronoi edge midpoints (just normals)- edpos=3
  ! HTC - velocities on the intersection of Voronoi-Triangle edges
  !            (just normals - relative to Voronoi cells) - edpos=6
  character (len=6)::  stag
  logical :: useStagHTC
  logical :: useStagHC

  !Edge positioning index
  integer(i4):: edpos

  !Scalar interpolations
  ! trsk - use trisk interpolations (simple area/length averages)
  ! bary - use 2nd order barycentric interpolation
  character (len=6)::  sinterpol
  logical :: useSinterpolTrisk !Trisk interpolation methods (area averaging)
  logical :: useSinterpolBary  !Linear barycentric interpolation

  !Scalar interpolations
  ! trsk - use trisk interpolations (simple area/length averages)
  ! bary - use 2nd order barycentric interpolation
  character (len=6)::  gradmtd
  logical :: useGradmtdTrisk !Simple diference
  logical :: useGradmtdBary  !Linear barycentric gradient

  !Plotsteps - will plot at every plotsteps timesteps
  integer (i4):: plotsteps
  logical:: plots   !Logical for plots or not
  integer (i4):: nplots !Number of time to plot
  integer (i4):: nprints !Number of time to plot
  integer (i4):: iploterrors !Logical for plots of errors 0 or 1!
  logical:: ploterrors   !Logical for plots of errors

  !Name for files and kind of staggering
  character (len=128)::  swmname

  !Total and initial mass
  real(r8):: tmass, inimass

  !Total energy
  real(r8):: Tenergy, Kenergy, Penergy, Availenergy
  real(r8):: Tenergy0, Kenergy0, Penergy0, Availenergy0

  !Reference Maximum values - used for normalization
  real(r8):: maxvel !Velocities
  real(r8):: maxh   !Fluid thickness


  !-------------------------------------
  !Fields given at grid points
  !-------------------------------------

  !Tendencies
  type(scalar_field):: momeq !Momentum
  type(scalar_field):: masseq !Mass
  type(scalar_field):: masseq_exact
  type(scalar_field):: masseq_error
  type(scalar_field):: momeq_exact
  type(scalar_field):: momeq_error

  !Velocities (defined on edges - only normal component)
  type(scalar_field):: u  !General
  type(scalar_field):: u_old  !Temporary
  type(scalar_field):: u_0 !Initial
  type(scalar_field):: u_error !Error
  type(scalar_field):: u_exact  !Exact
  type(scalar_field):: uh  !Velocity X thickness
  type(scalar_field):: uh_exact  !Velocity X thickness
  type(scalar_field):: uh_error  !Velocity X thickness

  type(scalar_field):: uhq_perp  !Perpendicular component (tangent)
  type(scalar_field):: uhq_perp_error  !Perpendicular component (tangent)
  type(scalar_field):: uhq_perp_exact  !Perpendicular component (tangent)

  type(vector_field_cart):: v_hx  !Full vector velocity at cell centers
  type(vector_field_cart):: v_ed  !Full vector velocity at edges
  type(vector_field_cart):: v_ed_exact  !Full exact vector velocity at edge
  type(vector_field_cart):: vhq_tr  !Velocity X thickness x pv (tr)
  type(vector_field_cart):: vhq_tr_exact
  type(vector_field_cart):: vhq_tr_error
  type(vector_field_cart):: vhq_hx  !Velocity X thickness x pv (hx)
  type(vector_field_cart):: vhq_hx_exact
  type(vector_field_cart):: vhq_hx_error
  type(vector_field_cart):: vh_hx  !Velocity X thickness  (hx)
  type(vector_field_cart):: vh_hx_exact
  type(vector_field_cart):: vh_hx_error

  !Fluid thickness (defined on voronoi centers)
  type(scalar_field):: h  !General
  type(scalar_field):: h_old  !Temporary
  type(scalar_field):: h_error  !Error
  type(scalar_field):: h_exact  !Exact
  type(scalar_field):: h_ed  !on edges
  type(scalar_field):: h_ed_exact  !Exact
  type(scalar_field):: h_ed_error  !Error
  type(scalar_field):: h_tr  !on tr cc
  type(scalar_field):: h_tr_exact  !Exact
  type(scalar_field):: h_tr_error  !Error

  !Vorticity
  type(scalar_field):: eta   !absolute vorticity - on triangles
  type(scalar_field):: eta_exact   !absolute vorticity - on triangles
  type(scalar_field):: eta_error   !absolute vorticity - on triangles
  type(scalar_field):: q_tr  !pv on triangles
  type(scalar_field):: q_tr_exact  !pv on triangles
  type(scalar_field):: q_tr_error  !pv on triangles
  type(scalar_field):: q_hx  !pv on hexs
  type(scalar_field):: q_hx_exact  !pv on hexs
  type(scalar_field):: q_hx_error  !pv on hexs
  type(scalar_field):: q_ed  !pv on edges
  type(scalar_field):: q_ed_exact  !pv on edges
  type(scalar_field):: q_ed_error  !pv on edges

  type(scalar_field):: q_grad_ed  !grad pv on edges
  type(vector_field_cart):: q_grad_tr  !graf pv on hexs
  type(vector_field_cart):: q_grad_tr_exact  !pv on hexs
  type(vector_field_cart):: q_grad_tr_error  !pv on hexs

  !Kinectic energy
  type(scalar_field):: Kin_energy  !on cell center
  type(scalar_field):: Kin_energy_exact  !on cell center
  type(scalar_field):: Kin_energy_error  !on cell center
  type(scalar_field):: Kin_energy_tr  !on tr center
  type(scalar_field):: Kin_energy_tr_exact  !on tr center
  type(scalar_field):: Kin_energy_tr_error  !on tr center
  !type(scalar_field):: Kin_energy_ed  !on edges

  !Auxiliar field
  type(scalar_field):: ghbK  !g x (h + b) + K field at Cell centers
  type(scalar_field):: grad_ghbK  !Gradient of ghbk - at edges
  type(scalar_field):: grad_ghbK_exact  !Gradient of ghbk
  type(scalar_field):: grad_ghbK_error  !Gradient of ghbk
  type(scalar_field):: grad_h  !Gradient of h - at edges

  !Total energy
  type(scalar_field):: Tot_energy  !on cell center
  type(scalar_field):: Tot_energy_ed  !on edges

  !Bottom topography (at cell nodes)
  type(scalar_field):: bt
  type(scalar_field):: hbt  !h+bt

  !Divergence at cell
  type(scalar_field):: divuh
  type(scalar_field):: divuh_exact
  type(scalar_field):: divuh_error
  type(scalar_field):: divu

  !Divergence at dual (with vorticity)
  type(scalar_field):: divueta
  type(scalar_field):: divueta_exact
  type(scalar_field):: divueta_error

  !Vector Laplacian (at edges)
  type(scalar_field):: lapu
  type(scalar_field):: lapu_exact
  type(scalar_field):: lapu_error

  !Trsk consistency index
  type(scalar_field):: trskind !Low values indicate godd for trisk
  logical, allocatable:: isTrskindLow(:) !isTrskindHigh=true if > trskindmax
  real(r8):: trskindmax  !Maximum - above that - use Pered

  !Wachspress coordinates vector for interpolation TR -> V
  type(general_coords), allocatable:: wachc_tr2v(:)

  !Reference solution from END GAME
  integer(i4):: nreftime
  real(r8), allocatable:: reftimes(:)

  !Reference solution from Spectral model
  ! (nlon, nlat, ntime)
  integer(i4):: nlonref, nlatref, ntimeref
  real(r8), allocatable:: href(:,:,:), uref(:,:,:), vref(:,:,:)

  !Latitude levels for Gaussian grid
  real(r8), allocatable:: latlevels(:)

  !Runge Kutta variables
  !Right hand side of mass equation (number of cell equations)
  real(r8), dimension(:), allocatable:: massf0
  real(r8), dimension(:), allocatable:: massf1
  real(r8), dimension(:), allocatable:: massf2
  real(r8), dimension(:), allocatable:: massf3

  !Right hand side of momentum equation (number of edge equations)
  real(r8), dimension(:), allocatable:: momf0
  real(r8), dimension(:), allocatable:: momf1
  real(r8), dimension(:), allocatable:: momf2
  real(r8), dimension(:), allocatable:: momf3

contains

  subroutine swmpars(usetime)
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
    filename=trim(pardir)//"swm.par"
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
    read(fileunit,*)  reconmtd, gasscoef
    read(fileunit,*)  buffer
    read(fileunit,*)  coriolis_reconmtd
    read(fileunit,*)  buffer
    read(fileunit,*)  sinterpol
    read(fileunit,*)  buffer
    read(fileunit,*)  gradmtd
    read(fileunit,*)  buffer
    read(fileunit,*)  nplots, nprints, iploterrors
    read(fileunit,*)  buffer
    read(fileunit,*)  difus
    read(fileunit,*)  buffer
    read(fileunit,*)  hollgw
    read(fileunit,*)  buffer
    read(fileunit,*)  trskindmax
    read(fileunit,*)  buffer
    read(fileunit,*)  pv_stab, pvspar
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

    !Check scalar interpol
    useSinterpolTrisk=trim(sinterpol)=="trsk"
    useSinterpolBary=trim(sinterpol)=="bary"
    if((.not.useSinterpolTrisk).and.(.not.useSinterpolBary))then
       print*, "Unknown interpolation", sinterpol
       stop
    end if

    !Check  Vector reconstruction method
    useCoriolisMtdPered=trim(coriolis_reconmtd)=="pered"
    useCoriolisMtdDtred=trim(coriolis_reconmtd)=="dtred"
    useCoriolisMtdTrisk=trim(coriolis_reconmtd)=="trsk"
    useCoriolisMtdHyb=trim(coriolis_reconmtd)=="hyb"
    useCoriolisMtdExact=trim(coriolis_reconmtd)=="exact"
    if((.not. useCoriolisMtdPered).and.(.not. useCoriolisMtdTrisk) &
         .and.(.not. useCoriolisMtdDtred).and.(.not. useCoriolisMtdHyb) &
         .and.(.not. useCoriolisMtdExact))then
       print*, "Unknown Coriolis vector reconstruction method", coriolis_reconmtd
       stop
    end if

    !Check  Vector reconstruction method
    useReconmtdPerhx=trim(reconmtd)=="perhx"
    useReconmtdTrisk=trim(reconmtd)=="trsk"
    useReconmtdGass=trim(reconmtd)=="gass"
    useReconmtdMelv=trim(reconmtd)=="melv"
    if((.not. useReconmtdPerhx).and.(.not. useReconmtdTrisk) &
         .and.(.not. useReconmtdGass).and.(.not. useReconmtdMelv))then
       print*, "Unknown vector reconstruction method", reconmtd
       stop
    end if

    !Check gradmetd
    useGradmtdTrisk=trim(gradmtd)=="trsk"
    useGradmtdBary=trim(gradmtd)=="bary"
    if((.not.useGradmtdTrisk).and.(.not.useGradmtdBary))then
       print*, "Unknown gradient discretization method", gradmtd
       stop
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
    swmname="swm_tc"//trim(adjustl(trim(atmp)))

    !Print information
    print*, "Test Case (Will1994)    : ", testcase
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
    print*, "Scalar interpolation    : ", sinterpol
    print*, "Vector recon method     : ", reconmtd
    print*, "Coriolis recon method   : ", coriolis_reconmtd
    print*, "Gradient method         : ", gradmtd
    print*, "Hollingsworth parameter : ", hollgw
    print*, "PV stabilization mtd    : ", pv_stab
    print*
    swmname=trim(swmname)//"_"//trim(adjustl(trim(stag)))
    swmname=trim(swmname)//"_vrec"//trim(adjustl(trim(reconmtd)))
    swmname=trim(swmname)//"_crec"//trim(adjustl(trim(coriolis_reconmtd)))
    swmname=trim(swmname)//"_sint"//trim(adjustl(trim(sinterpol)))
    swmname=trim(swmname)//"_grd"//trim(adjustl(trim(gradmtd)))
    if((.not.useOrigPV))then
       write(atmp,'(f5.2)') real(pvspar)
       swmname=trim(swmname)//"_pvs"//trim(adjustl(trim(pv_stab)))
       swmname=trim(swmname)//trim(adjustl(trim(atmp)))
    end if
    if(hollgw>0 .and. (testcase==32 .or. testcase==33 .or. testcase==34) )then
       write(atmp,'(f6.2)') real(hollgw)
       swmname=trim(swmname)//"_hol"//trim(adjustl(trim(atmp)))
       !print*, atmp
    end if

    RefSolRead=testcase==5.or. testcase==51.or.testcase==6.or.testcase==21.or.testcase==23
    RefSolAnal= testcase==1.or.testcase==2.or. testcase==22.or. testcase==24 &
         .or. testcase==32.or. testcase==33 .or. testcase==34 .or. &
         testcase==40 .or. testcase==41.or. testcase==42

    print*, "SWM Name for Plots: ", trim(swmname)
    print*

    return
  end subroutine swmpars

  subroutine allocate_globalswmvars()
    !--------------------------------
    !Allocate fields for SWM
    !------------------------------------

    !Error flag
    integer(i4):: ist

    if(useStagHTC)then
       !Instersection of tr edge and hx edge
       edpos=6
    elseif(useStagHC)then
       !Midpoint of Voronoi edges
       edpos=3
    else
       print*, "Staggered undefined", stag
       stop
    end if

    !Tendencies
    momeq%n=mesh%ne
    momeq%name="momeq"
    momeq%pos=edpos
    allocate(momeq%f(1:momeq%n), stat=ist)  !General
    if(test_lterror==1)then
       allocate(momeq_exact%f(1:momeq%n), stat=ist)
       allocate(momeq_error%f(1:momeq%n), stat=ist)
    end if

    masseq%n=mesh%nv
    masseq%name="masseq"
    masseq%pos=0
    allocate(masseq%f(1:masseq%n), stat=ist)  !General
    if(test_lterror==1)then
       allocate(masseq_exact%f(1:masseq%n), stat=ist)
       allocate(masseq_error%f(1:masseq%n), stat=ist)
    end if

    !Velocities (defined on edges - only normal component)
    u%n=mesh%ne
    u%name="u"
    u%pos=edpos
    allocate(u%f(1:u%n), stat=ist)  !General
    allocate(u_old%f(1:u%n), stat=ist)
    allocate(u_0%f(1:u%n), stat=ist)
    allocate(u_error%f(1:u%n), stat=ist)
    allocate(u_exact%f(1:u%n), stat=ist)
    allocate(uh%f(1:u%n), stat=ist)

    !Perp velocities (defined on edges - only normal component)
    uhq_perp%n=mesh%ne
    uhq_perp%name="uhq_perp"
    uhq_perp%pos=edpos
    allocate(uhq_perp%f(1:uhq_perp%n), stat=ist)
    if(test_lterror==1)then
       allocate(uhq_perp_exact%f(1:u%n), stat=ist)
       allocate(uhq_perp_error%f(1:u%n), stat=ist)
    end if

    v_hx%n=mesh%nv
    v_hx%name="v"
    v_hx%pos=0
    allocate(v_hx%p(1:v_hx%n), stat=ist)

    v_ed%n=mesh%ne
    v_ed%name="v_ed"
    v_ed%pos=edpos
    allocate(v_ed%p(1:v_ed%n), stat=ist)

    !Vector velocities
    if(test_lterror==1)then
       v_ed_exact%n=mesh%ne
       v_ed_exact%name="v_ed_exact"
       v_ed_exact%pos=edpos
       allocate(v_ed_exact%p(1:v_ed_exact%n), stat=ist)
    end if

    h%n=mesh%nv
    h%name="h"
    h%pos=0
    allocate(h%f(1:h%n), stat=ist)
    allocate(h_old%f(1:h%n), stat=ist)
    allocate(h_error%f(1:h%n), stat=ist)
    allocate(h_exact%f(1:h%n), stat=ist)

    h_ed%n=mesh%ne
    h_ed%name="h_ed"
    h_ed%pos=edpos
    allocate(h_ed%f(1:h_ed%n), stat=ist)
    if(test_lterror==1)then
       allocate(h_ed_exact%f(1:h_ed%n), stat=ist)
       allocate(h_ed_error%f(1:h_ed%n), stat=ist)
    end if

    h_tr%n=mesh%nt
    h_tr%name="h_tr"
    h_tr%pos=1
    allocate(h_tr%f(1:h_tr%n), stat=ist)
    if(test_lterror==1)then
       allocate(h_tr_exact%f(1:h_tr%n), stat=ist)
       allocate(h_tr_error%f(1:h_tr%n), stat=ist)
    end if

    vhq_tr%n=mesh%nt
    vhq_tr%name="vhq_tr"
    vhq_tr%pos=1
    allocate(vhq_tr%p(1:vhq_tr%n), stat=ist)

    if(test_lterror==1)then
       allocate(vhq_tr_exact%p(1:vhq_tr%n), stat=ist)
       allocate(vhq_tr_error%p(1:vhq_tr%n), stat=ist)
    end if

    vhq_hx%n=mesh%nv
    vhq_hx%name="vhq_hx"
    vhq_hx%pos=0
    allocate(vhq_hx%p(1:vhq_hx%n), stat=ist)

    if(test_lterror==1)then
       allocate(vhq_hx_exact%p(1:vhq_hx%n), stat=ist)
       allocate(vhq_hx_error%p(1:vhq_hx%n), stat=ist)
    end if

    if(test_lterror==1)then
       vh_hx%n=mesh%nv
       vh_hx%name="vh_hx"
       vh_hx%pos=0
       allocate(vh_hx%p(1:vh_hx%n), stat=ist)
       allocate(vh_hx_exact%p(1:vhq_hx%n), stat=ist)
       allocate(vh_hx_error%p(1:vhq_hx%n), stat=ist)
    end if

    !Absolute vorticity
    eta%n=mesh%nt
    eta%pos=1
    eta%name="eta"
    allocate(eta%f(1:eta%n), stat=ist)

    if(test_lterror==1)then
       allocate(eta_exact%f(1:eta%n), stat=ist)
       allocate(eta_error%f(1:eta%n), stat=ist)
    end if

    !Potential vorticity
    q_tr%n=mesh%nt
    q_tr%pos=1
    q_tr%name="q_tr"
    allocate(q_tr%f(1:q_tr%n), stat=ist)

    !if(test_lterror==1)then
    allocate(q_tr_exact%f(1:q_tr%n), stat=ist)
    allocate(q_tr_error%f(1:q_tr%n), stat=ist)
    !end if

    if(test_lterror==1)then
       q_hx%n=mesh%nv
       q_hx%pos=0
       q_hx%name="q_hx"
       allocate(q_hx%f(1:q_hx%n), stat=ist)
       allocate(q_hx_exact%f(1:q_hx%n), stat=ist)
       allocate(q_hx_error%f(1:q_hx%n), stat=ist)

    end if

    q_ed%n=mesh%ne
    q_ed%pos=edpos
    q_ed%name="q_ed"
    allocate(q_ed%f(1:q_ed%n), stat=ist)
    if(test_lterror==1)then
       allocate(q_ed_exact%f(1:q_ed%n), stat=ist)
       allocate(q_ed_error%f(1:q_ed%n), stat=ist)
    end if

    q_grad_ed%n=mesh%ne
    q_grad_ed%pos=edpos
    q_grad_ed%name="q_grad_ed"
    allocate(q_grad_ed%f(1:q_ed%n), stat=ist)

    q_grad_tr%n=mesh%nt
    q_grad_tr%pos=1
    q_grad_tr%name="q_grad_tr"
    allocate(q_grad_tr%p(1:q_grad_tr%n), stat=ist)

    if(test_lterror==1)then
       allocate(q_grad_tr_exact%p(1:q_grad_tr%n), stat=ist)
       allocate(q_grad_tr_error%p(1:q_grad_tr%n), stat=ist)
    end if

    !Kinectic energy
    Kin_energy%n=mesh%nv
    Kin_energy%pos=0
    Kin_energy%name="Kin_energy"
    allocate(Kin_energy%f(1:Kin_energy%n), stat=ist)

    if(test_lterror==1)then
       allocate(Kin_energy_exact%f(1:Kin_energy%n), stat=ist)
       allocate(Kin_energy_error%f(1:Kin_energy%n), stat=ist)
    end if

    Kin_energy_tr%n=mesh%nt
    Kin_energy_tr%pos=1
    Kin_energy_tr%name="Kin_energy_tr"
    allocate(Kin_energy_tr%f(1:Kin_energy_tr%n), stat=ist)
    if(test_lterror==1)then
       allocate(Kin_energy_tr_exact%f(1:Kin_energy_tr%n), stat=ist)
       allocate(Kin_energy_tr_error%f(1:Kin_energy_tr%n), stat=ist)
    end if

    !Auxiliar field
    ghbK%n=mesh%nv
    ghbK%pos=0
    ghbK%name="ghbK"
    allocate(ghbK%f(1:ghbK%n), stat=ist)

    grad_ghbK%n=mesh%ne
    grad_ghbK%pos=edpos
    grad_ghbK%name="grad_ghbK"
    allocate(grad_ghbK%f(1:grad_ghbK%n), stat=ist)

    if(test_lterror==1)then
       allocate(grad_ghbK_exact%f(1:grad_ghbK%n), stat=ist)
       allocate(grad_ghbK_error%f(1:grad_ghbK%n), stat=ist)

       grad_h%n=mesh%ne
       grad_h%pos=edpos
       grad_h%name="grad_ed"
       allocate(grad_h%f(1:grad_h%n), stat=ist)
    end if

    !Bottom topography
    bt%n=mesh%nv
    bt%pos=0
    bt%name="bt"
    allocate(bt%f(1:bt%n), stat=ist)

    !H+bt
    hbt%n=mesh%nv
    hbt%pos=0
    hbt%name="hbt"
    allocate(hbt%f(1:bt%n), stat=ist)

    !Divergence
    divuh%n=mesh%nv
    divuh%pos=0
    divuh%name="divuh"
    allocate(divuh%f(1:divuh%n), stat=ist)

    if(test_lterror==1)then
       allocate(divuh_exact%f(1:divuh%n), stat=ist)
       allocate(divuh_error%f(1:divuh%n), stat=ist)
    end if

    !Divergence
    divu%n=mesh%nv
    divu%pos=0
    divu%name="divu"
    allocate(divu%f(1:divu%n), stat=ist)

    if(test_lterror==1)then
       divueta%n=mesh%nt
       divueta%pos=1
       divueta%name="divueta"
       allocate(divueta%f(1:divueta%n), stat=ist)
       allocate(divueta_exact%f(1:divueta%n), stat=ist)
       allocate(divueta_error%f(1:divueta%n), stat=ist)
    end if

    lapu%n=mesh%ne
    lapu%pos=edpos
    lapu%name="lapu"
    allocate(lapu%f(1:lapu%n), stat=ist)

    if(test_lterror==1)then
       allocate(lapu_exact%f(1:lapu%n), stat=ist)
       allocate(lapu_error%f(1:lapu%n), stat=ist)

    end if

    if(useCoriolisMtdHyb)then
       trskind%n=mesh%ne
       trskind%pos=u%pos
       trskind%name="trskind"
       allocate(trskind%f(1:trskind%n), stat=ist)
       allocate(isTrskindLow(1:trskind%n), stat=ist)
    end if

    !Wachspress coordinates
    if(useSinterpolBary)then
       allocate(wachc_tr2v(1 :mesh%nv), stat=ist)
    end if

    !Runge kutta variables
    allocate(massf0(1:h%n))
    allocate(massf1(1:h%n))
    allocate(massf2(1:h%n))
    allocate(massf3(1:h%n))
    allocate(momf0(1:u%n))
    allocate(momf1(1:u%n))
    allocate(momf2(1:u%n))
    allocate(momf3(1:u%n))

    IF(ist>0) STOP 'Error in allocate_globalswmvars'

    call initialize_globalvars()

    return
  end subroutine allocate_globalswmvars

  subroutine initialize_globalvars()
    !---------------------------------------------------
    !Initialize fields with zero in paralel
    !---------------------------------------------------

    integer(i4)::i, k, l

    !First initialize vector variable -so that they can be used later in
    !  workshare enviroment

    !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
    !$OMP PRIVATE(i) SHARED(v_hx)
    do i=1, v_hx%n
       v_hx%p(i)%v=0._r8
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
    !$OMP PRIVATE(k) SHARED(vhq_tr, q_grad_tr)
    do k=1, vhq_tr%n
       vhq_tr%p(k)%v(1:3)=0.
       q_grad_tr%p(k)%v(1:3)=0.
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
    !$OMP PRIVATE(l) SHARED(v_ed)
    do l=1, v_ed%n
       v_ed%p(l)%v=0.
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(test_lterror, momeq, masseq, u, u_old, u_0) &
    !$OMP SHARED(u_error, u_exact, uh, uhq_perp) &
    !$OMP SHARED(v_hx, v_ed, vhq_hx) &
    !$OMP SHARED(h, h_old, h_0, h_error, h_exact, h_ed) &
    !$OMP SHARED(h_tr, eta, q_tr, q_ed) &
    !$OMP SHARED(Kin_energy, Kin_energy_tr, ghbK, grad_ghbK, hbt, bt) &
    !$OMP SHARED(divuh, q_tr_exact, q_tr_error, q_grad_ed, divu, lapu ) &
    !$OMP SHARED(massf0, massf1, massf2, massf3) &
    !$OMP SHARED(momf0, momf1, momf2, momf3)
    momeq%f(1:u%n)=0._r8
    masseq%f(1:h%n)=0._r8
    !Velocities
    u%f(1:u%n)=0._r8
    u_old=u
    u_0=u
    u_error=u
    u_exact=u
    uh=u
    uhq_perp%f=0._r8
    vhq_hx=v_hx

    !Mass-Thickness
    h%f=0._r8
    h_old=h
    !h_0=h
    h_error=h
    h_exact=h
    h_ed%f=0._r8
    h_tr%f=0._r8

    !Vorticity
    eta%f=0._r8
    q_tr%f=0._r8
    q_ed%f=0._r8

    !Energy
    Kin_energy%f=0._r8
    Kin_energy_tr%f=0._r8
    ghbK%f=0._r8
    grad_ghbK%f=0._r8

    !Bottom Topography
    bt%f=0._r8
    hbt%f=0._r8

    !Divergence
    divuh%f=0._r8
    divu%f=0._r8

    !Laplacian
    lapu%f=0._r8

    !PV analysis
    q_tr_exact=q_tr
    q_tr_error=q_tr

    q_grad_ed%f=0._r8

    !Runge kutta variables
    massf0(1:h%n)=0._r8
    massf1(1:h%n)=0._r8
    massf2(1:h%n)=0._r8
    massf3(1:h%n)=0._r8
    momf0(1:u%n)=0._r8
    momf1(1:u%n)=0._r8
    momf2(1:u%n)=0._r8
    momf3(1:u%n)=0._r8
    !$OMP END PARALLEL WORKSHARE


    if(test_lterror==1)then
       !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
       !$OMP SHARED(test_lterror, uhq_perp, uhq_perp_exact, uhq_perp_error) &
       !$OMP SHARED(v_ed, v_ed_exact, h_ed, h_ed_exact, h_ed_error) &
       !$OMP SHARED(h_tr, h_tr_exact, h_tr_error) &
       !$OMP SHARED(vhq_tr, vhq_tr_exact, vhq_tr_error) &
       !$OMP SHARED(vhq_hx, vhq_hx_exact, vhq_hx_error) &
       !$OMP SHARED(vh_hx, vh_hx_exact, vh_hx_error) &
       !$OMP SHARED(eta, eta_exact, eta_error) &
       !$OMP SHARED(q_tr, q_tr_exact, q_tr_error) &
       !$OMP SHARED(q_hx, q_hx_exact, q_hx_error) &
       !$OMP SHARED(q_ed, q_ed_exact, q_ed_error) &
       !$OMP SHARED( q_grad_tr, q_grad_tr_exact, q_grad_tr_error) &
       !$OMP SHARED(Kin_energy, Kin_energy_exact, Kin_energy_error) &
       !$OMP SHARED(kin_energy_tr_exact, kin_energy_tr, kin_energy_tr_error) &
       !$OMP SHARED(grad_ghbK, grad_ghbK_exact, grad_ghbK_error, grad_h) &
       !$OMP SHARED(divuh, divuh_exact, divuh_error) &
       !$OMP SHARED(divueta, divueta_exact, divueta_error) &
       !$OMP SHARED(lapu, lapu_exact, lapu_error) &
       !$OMP SHARED(momeq,  masseq, momeq_exact, masseq_exact, momeq_error, masseq_error)
       momeq_exact=momeq
       masseq_exact=masseq
       momeq_error=momeq
       masseq_error=masseq
       uhq_perp_exact=uhq_perp
       uhq_perp_error=uhq_perp
       h_ed_exact=h_ed
       h_ed_error=h_ed
       h_tr_exact=h_tr
       h_tr_error=h_tr
       v_ed_exact=v_ed
       vhq_tr_exact=vhq_tr
       vhq_tr_error=vhq_tr
       vhq_hx_exact=vhq_hx
       vhq_hx_error=vhq_hx
       vh_hx=vhq_hx
       vh_hx_exact=vh_hx
       vh_hx_error=vh_hx
       eta_exact=eta
       eta_error=eta
       q_tr_exact=q_tr
       q_tr_error=q_tr
       q_grad_tr_exact=q_grad_tr
       q_grad_tr_error=q_grad_tr
       q_hx%f=0._r8
       q_hx_exact=q_hx
       q_hx_error=q_hx
       q_ed_exact=q_ed
       q_ed_error=q_ed
       Kin_energy_exact=Kin_energy
       Kin_energy_error=Kin_energy
       Kin_energy_tr_exact=Kin_energy_tr
       Kin_energy_tr_error=Kin_energy_tr
       grad_ghbK_exact=grad_ghbK
       grad_ghbK_error=grad_ghbK
       grad_h=grad_ghbK
       divuh_exact=divuh
       divuh_error=divuh
       divueta=q_tr
       divueta_exact=divueta
       divueta_error=divueta
       lapu_exact=lapu
       lapu_error=lapu
       !$OMP END PARALLEL WORKSHARE

    end if

  end subroutine initialize_globalvars

end module dataswm
