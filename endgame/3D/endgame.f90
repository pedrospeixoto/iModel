! LOOSE ENDS
!
!
! 3. Proper output routines
!
! 11. Implement SLICE !!!
!
!
! =====================================================

MODULE runtype

! Basic information on the run
! Switches are read from namelist rundata;
! default values are set here.

! Five character run ID
CHARACTER*5 :: runid = '00001'

! True for restart, false to set initial data in code
LOGICAL :: lrestart = .false.

! Name of file from which to restart if lrestart = .true.
CHARACTER*25 :: yrestart = 'run00001restart0000000192'


END MODULE runtype

! ======================================================

MODULE switches

! Switches controlling options
! Switches are read from namelist switchdata;
! default values are set here.

! True for SLICE treatment of density, false for
! semi-Lagrangian treatment of density.
! Only the latter is currently implemented
LOGICAL :: slice = .FALSE.

! 1 for non-hydrostatic, 0 for quasi-hydrostatic
! Note: for a hydrostatic solution alpha_w should also be set to 1.0d0
INTEGER :: delta_v = 1

! True for monotone advection of theta, false for
! standard cubic Lagrange.
LOGICAL :: monotone = .TRUE.


END MODULE switches

! ======================================================

MODULE grid

! Information about the grid size and resolution

! Number of points in each direction
! INTEGER, PARAMETER :: p = 7, nx = 2**p, ny = 2**(p-1), nyp = ny + 1
INTEGER, PARAMETER :: p = 7, nx = 192, ny = 96, nyp = ny + 1
INTEGER, PARAMETER :: nz = 30, nzp = nz + 1

! Useful constants
REAL*8, PARAMETER :: pi = 3.141592653589793d0, twopi = 2.0d0*pi, &
                     piby2 = pi/2.0d0, dx = 2.0d0*pi/nx, dy=pi/ny

! Rotation of model coordinates relative to geographical coordinates
REAL*8, PARAMETER :: rotgrid = 0.0d0*pi

! One-dimensional arrays of coordinates and related quantities
REAL*8 :: xp(nx), yp(ny), xu(nx), yu(ny), xv(nx), yv(nyp), &
          cosp(ny), cosv(nyp), sinp(ny), sinv(nyp), area(ny), &
          areatot, &
	  etap(nz), etaw(nzp)

! Geographical longitude and latitude and
! components of 2 x Earth's rotation vector
REAL*8 :: geolonu(nx,ny),  geolatu(nx,ny), &
          geolonv(nx,nyp), geolatv(nx,nyp), &
          geolonp(nx,ny),  geolatp(nx,ny), &
          coriol1(nx,ny),  coriol2(nx,ny),  coriol3(nx,ny)

! Surface and 3D fields of distance from Earth's centre
REAL*8 :: rgeoid(nx,ny), rsurf(nx,ny), rp(nz,nx,ny), rw(nzp,nx,ny)

! Areas of grid cell faces (below, to the West, and to the South
! of each cell) and grid cell volume
REAL*8 :: areab(nzp,nx,ny), areaw(nz,nx,ny), areas(nz,nx,nyp), &
          volume(nz,nx,ny)

! Coefficients used to interpolate from w levels to p levels
! and from p levels to w levels
REAL*8 :: above1(nz), below1(nz), above2(nzp), below2(nzp)


END MODULE grid

! ======================================================

MODULE state

! Variables describing the current model state

USE grid

INTEGER, PARAMETER :: ntracers = 2

! Current time step
REAL*8 :: rho0(nz,nx,ny), theta0(nzp,nx,ny), &
          u0(nz,nx,ny), v0(nz,nx,nyp), w0(nzp,nx,ny), &
	  etadot0(nzp,nx,ny)
REAL*8, ALLOCATABLE :: tracers0(:,:,:,:)

! New time step
REAL*8 :: rho(nz,nx,ny), theta(nzp,nx,ny), &
          u(nz,nx,ny), v(nz,nx,nyp), w(nzp,nx,ny), &
	  etadot(nzp,nx,ny)
REAL*8, ALLOCATABLE :: tracers(:,:,:,:)

! Geopotential
REAL*8 :: phi(nz,nx,ny)

END MODULE state

! ======================================================

MODULE trajectories

! Variables describing Lagrangian trajectories

INTEGER, PARAMETER :: ntrajectories = 44310
INTEGER, PARAMETER :: nlabels = 2

! Coordinates of current locations
REAL*8, ALLOCATABLE :: xtraj(:), ytraj(:), etraj(:), labels(:,:)

END MODULE trajectories

! ======================================================

MODULE departure

! Departure points for u points, v points, w points and
! p points

USE grid

REAL*8 :: xdepu(nz,nx,ny),   ydepu(nz,nx,ny),   edepu(nz,nx,ny),  &
          xdepv(nz,nx,nyp),  ydepv(nz,nx,nyp),  edepv(nz,nx,nyp), &
          xdepw(nzp,nx,ny),  ydepw(nzp,nx,ny),  edepw(nzp,nx,ny), &
          xdepp(nz,nx,ny),   ydepp(nz,nx,ny),   edepp(nz,nx,ny)

END MODULE departure

! ======================================================

MODULE work

! Arrays used as workspace for RHS calculations

USE grid

REAL*8 :: rhsgridrho(nz,nx,ny), rhsgridtheta(nzp,nx,ny), &
          rhsgridu(nz,nx,ny), rhsgridv(nz,nx,nyp), rhsgridw(nzp,nx,ny)
REAL*8 :: rhsdeprho(nz,nx,ny), rhsdeptheta(nzp,nx,ny), &
          rhsdepu(nz,nx,ny), rhsdepv(nz,nx,nyp), rhsdepw(nzp,nx,ny)
REAL*8 :: rhsrho(nz,nx,ny), rhstheta(nzp,nx,ny), &
          rhsu(nz,nx,ny), rhsv(nz,nx,nyp), rhsw(nzp,nx,ny)
REAL*8 :: q16(nzp,nx,ny), rhs_helm(nz,nx,ny)
REAL*8 :: dv3d(nzp,nx,ny)
LOGICAL :: ldiag


END MODULE work

! ======================================================

MODULE increments

! Increments computed by Helmholtz solver and back substitution

USE grid

REAL*8 :: exner_inc(nz,nx,ny), rho_inc(nz,nx,ny), theta_inc(nzp,nx,ny), &
          u_inc(nz,nx,ny), v_inc(nz,nx,nyp), w_inc(nzp,nx,ny)

END MODULE increments

! ======================================================

MODULE timestep

! Data related to time stepping

! Implicitness parameters (0.5 is centred)
REAL*8, PARAMETER :: alpha_x     = 0.5d0, beta_x     = 1.0d0 - alpha_x
REAL*8, PARAMETER :: alpha_rho   = 0.501d0, beta_rho   = 1.0d0 - alpha_rho
REAL*8, PARAMETER :: alpha_u     = 0.501d0, beta_u     = 1.0d0 - alpha_u
REAL*8, PARAMETER :: alpha_w     = 0.501d0, beta_w     = 1.0d0 - alpha_w

! Values read from namelist timedata. Default values are set here.
REAL*8 :: dt = 3600.0d0          ! Time step
INTEGER :: nstop = 10            ! Length of run
INTEGER :: noutput = 10          ! Number of steps between standard output
INTEGER :: nrestart = 10         ! Number of steps between restart dumps

! Current time step
INTEGER :: istep

! *** Could make alpha's namelist variables ***
! *** Could save some flops by precomputing alpha*dt and beta*dt ***


END MODULE timestep

! ======================================================

MODULE constants

! Physical constants:

! Earth's radius
REAL*8, PARAMETER :: rearth = 6371229.0d0
! Spherical planet
REAL*8, PARAMETER :: requat = rearth
REAL*8, PARAMETER :: rpole  = rearth
! Earth-like spheroid
!REAL*8, PARAMETER :: requat = 6378100.0d0
!REAL*8, PARAMETER :: rpole  = 6356800.0d0
! Jupiter-like spheroid
!REAL*8, PARAMETER :: requat = 6584800.0d0
!REAL*8, PARAMETER :: rpole  = 6157600.0d0

! Earth's rotation rate
REAL*8, PARAMETER :: rotatn = 7.29212d-5

! Domain depth
REAL*8, PARAMETER :: domain = 3.0d4

! Gravitational acceleration at the Earth's surface
REAL*8, PARAMETER :: gravity = 9.80616d0

! Gas constant for dry air
REAL*8, PARAMETER :: gascon = 278.04d0

! Specific heat capacity at constant pressure
REAL*8, PARAMETER :: cp = 1004.64d0

! kappa = R / cp
REAL*8, PARAMETER :: kappa = gascon/cp

! (1 - kappa) / kappa
REAL*8, PARAMETER :: onemkbyk = (1.0d0 - kappa)/kappa

! kappa / (1 - kappa)
REAL*8, PARAMETER :: kbyonemk = 1.0d0/onemkbyk

! Reference pressure used in Exner
REAL*8, PARAMETER :: p00 = 1.0d5


END MODULE constants

! ======================================================

MODULE refstate

! Reference state for iterative solver, and coefficients
! used for discrete Helmholtz problem

USE grid

! Reference state
REAL*8 :: rho_ref(nz,nx,ny), theta_ref(nzp,nx,ny), exner_ref(nz,nx,ny)

! Helmholtz coefficients
REAL*8 :: rtw(nz,nx,ny), rts(nz,nx,nyp), &
          cabove(nz,nx,ny), cbelow(nz,nx,ny), &
	  ccell(nz,nx,ny), bb_ref(nzp,nx,ny), &
	  dexnerdr_ref(nzp,nx,ny), thetabar_ref(nz,nx,ny), &
	  dthetabardr_ref(nzp,nx,ny), drhobardr_ref(nz,nx,ny)

END MODULE refstate

! ======================================================

MODULE lagdiagnostics

! Fields used in calculation of Lagrangian diagnostics

USE grid

INTEGER, PARAMETER :: nlayers = 32, nlayersp = nlayers + 1
REAL*8 :: thetalayers(nlayersp)
REAL*8 :: mpt(nlayers), massbelow(nlayersp)
REAL*8 :: pv(nz,nx,nyp)

END MODULE lagdiagnostics

! ======================================================

MODULE util

! Useful function

CONTAINS

FUNCTION near(dist,domain)

! To correct for wrap-around when evaluating distances
! between points

IMPLICIT NONE
REAL*8, INTENT(IN) :: dist, domain
REAL*8 :: near
REAL*8 :: hdom

hdom = 0.5*domain
near = dist
IF (dist .ge.  hdom) near = near - domain
IF (dist .le. -hdom) near = near + domain

END FUNCTION near

END MODULE util

! ====================================================

MODULE channels

! Tidy list of all I/O channels in one place to avoid accidental
! overuse of any channel number

INTEGER, PARAMETER :: channml = 20          ! For reading namelists
INTEGER, PARAMETER :: chanout = 30          ! Standard output files for IDL
INTEGER, PARAMETER :: chanoutm = 40         ! Standard output files for Matlab
INTEGER, PARAMETER :: chanresin = 50        ! Input channel for restart run
INTEGER, PARAMETER :: chanresout = 60       ! Restart dumps
INTEGER, PARAMETER :: chandump1 = 70        ! Quick look dump file for IDL (1D fields)
INTEGER, PARAMETER :: chandump2 = 71        ! Quick look dump file for IDL (2D fields)
INTEGER, PARAMETER :: chandumpm = 80        ! Quick look dump file for Matlab
INTEGER, PARAMETER :: chann48 = 48          ! Channel for reading N48 data
INTEGER, parameter :: chanmass = 49         ! Channel for global mass diagnostic
INTEGER, PARAMETER :: chanmasspt = 85       ! Channel for mass per unit theta diagnostic
INTEGER, PARAMETER :: chandiag = 84         ! Channel for global diagnostics
INTEGER, PARAMETER :: chanscatterm = 86     ! Channel for lagrangian conservation scatter plots
INTEGER, PARAMETER :: chantraj = 87         ! Channel for trajectory data

END MODULE channels

! ====================================================

MODULE globdiag

! Save some global diagnostics, allowing their change over
! one step to be computed

REAL*8 :: totalm0 = 0.0d0, internaleng0 = 0.0d0, potentialeng0 = 0.0d0, &
          kineticeng0 = 0.0d0, totaleng0 = 0.0d0, totalent0 = 0.0d0,    &
          totalam0 = 0.0d0, imbalance = 0.0d0,                          &
          unavailie = 0.0d0, unavailpe = 0.0d0

END MODULE globdiag

! ====================================================

PROGRAM endgame

! John Thuburn 19/9/2010
! Independent implementation of ENDGame version 3.02
! 
! There are some differences from the Version 3.02 documentation
!
! 1. Only deep atmosphere spherical geometry is implemented.
!
! 2. Only uniform spacing in model coordinate longitude and latitude
!    is implemented, though the model coordinates may be rotated relative
!    to geographical coordinates.
!
! 3. Only `Version 2' is implemented, i.e. when iterated to
!    convergence the solution is independent of the reference
!    state.
!
! 4. An incremental formulation of the iterative solver is used.
!    This allows the Helmholtz problem to be solved to a less tight
!    tolerance. It also simplifies some of the RHS calculations.
!
! 5. The elliptic solver (conditional semi-coarsening multigrid
!    in the horizontal, line solve in the vertical) can use a fully
!    three-dimensional reference state.
!
! 6. The simpler (rather than `flux-form') form of the gradient
!    operator in bent coordinates is used. We can then ensure (by
!    appropriate discretization of curl) that curl(grad()) = 0.
!
! 7. A simpler calculation of departure point velocities is used.
!    E.g. for u departure points, v and w are averaged to u points
!    using simple averaging; all velocity components are then interpolated
!    from u points to departure points before applying the rotation matrix.
!    The calculation at v and w points is similar. This potentially slightly
!    less accurate near the poles than the scheme documented in v3.02.
!
! 8. The gravitational term is coded as grad(Phi). This ensures no spurious
!    vorticity source provided grad and curl are implemented mimetically.
!    It also allows us to test the effect of horizontal components of gravity
!    as an alternative to the spheroidal coordinate approach. It also allows
!    a gravitational tide to be forced.
!
! 9. Top and bottom values of Psi_w are obtained differently - under review.
!

IMPLICIT NONE

! ----------------------------------------------------

! Computations and actions that are performed once at the start of an integration
CALL preliminary
print *,'Done preliminary'

! Time integration
CALL integrate

! Final actions at the end of the run
CALL epilogue

! -----------------------------------------------------

END PROGRAM endgame

! =====================================================

SUBROUTINE preliminary

! Computations and actions that are performed once at the start of an integration

USE channels
USE runtype
use grid
use constants
use state

IMPLICIT NONE
integer :: k, j, i

! Open dump file for matlab plotting
OPEN(chandumpm,FILE='dump.m')

! Open dump files for IDL plotting
OPEN(chandump1,FILE='dump1d.dat')
OPEN(chandump2,FILE='dump2d.dat')

! Read namelist data
CALL readnml

IF (.NOT. lrestart) THEN

  ! Initial run

  ! Set up grid
  CALL setupgrid
  print *,'Grid set up'

  ! Geopotential
  CALL setphi
  print *,'Geopotential set up'

  ! Initial condition
  CALL initial
  print *,'Initial data set'

  ! Initialize departure points
  CALL inidep
  print *,'Departure points initialized'

  ! Perform departure point iterations to obtain a reasonable
  ! initial estimate of departure points
  CALL depart_all2
  CALL depart_all2
  print *,'Preliminary departure iterations done twice'

  ! Check that initial data is reasonably well balanced and perform
  ! balancing iterations if necessary
  CALL checkbal

  ! Perturbation for baroclinic wave test
  ! print *,'*** U perturbation for baroclinic wave test ***'
  ! CALL inibwpert
  print *,'*** Perturbation for UMJS baroclinic wave test ***'
  CALL iniumjspert

ELSE

  ! Restart run

  print *,'**** NEED TO INCLUDE TRACERS AND TRAJECTORIES IN RESTART ****'
!  stop
  
  ! Set up grid   *** Ideally grid info should be read from restart file ***
  CALL setupgrid
  print *,'Grid set up'
  
  ! Geopotential
  CALL setphi
  print *,'Geopotential set up'
  
  ! Read in restart information
  CALL restart
  
ENDIF

! Initialize Lagrangian diagnostics
CALL start_diagnostics


END SUBROUTINE preliminary

! =====================================================

SUBROUTINE integrate

! Master routine controlling time integration

USE timestep
USE state           ! *** not needed here once a proper output routine is written
USE refstate        ! *** ditto ***
USE trajectories    ! *** ditto ***
USE lagdiagnostics  ! *** ditto ***
use constants

IMPLICIT NONE
INTEGER :: ii, i, ix, iy, iz
REAL*8 :: maxw, pmin, lonmin, latmin, ubar(nz,ny), thetabar(nzp,ny)

! Output inital state...
CALL masspertheta
CALL unavailPE2
CALL globaldiags
CALL potvort
print *,' '
call bwout
call dumpm(u(:,10,:),'ulon10',nz,ny)
!call dumpm(v(:,10,:),'vlon10',nz,nyp)
!call dumpm(w(:,10,:),'wlon10',nzp,ny)
!call dumpm(rho(:,10,:),'rholon10',nz,ny)
!call dumpm(theta(:,64,:),'thetalon64',nzp,ny)
!call dumpm(theta(:,:,64),'thetalat64',nzp,nx)
!call dumpm(u(1,:,:),'ulev1',nx,ny)
!call dumpm(v(10,:,:),'vlev10',nx,nyp)
!call dumpm(w(10,:,:),'wlev10',nx,ny)
!call dumpm(w(10,:,:),'wlev10',nx,ny)
!call dumpm(theta(20,:,:),'thetalev20',nx,ny)
!call dumpm(u(20,:,:),'ulev20',nx,ny)
! call dumpm(v(15,:,:),'vlev15',nx,nyp)
! call dumpm(w(15,:,:),'wlev15',nx,ny)
!call dumpm(theta(8,:,:),'thetalev8',nx,ny)
!call dumptrajm(xtraj,ytraj,etraj,ntrajectories,'xy')
!call dumpm(tracers(8,:,:,1),'trcrlev8',nx,ny)
!call dumpm(pv(8,:,:),'pvlev8',nx,nyp)
!call dumpm(tracers(8,:,:,2),'trcr2lev8',nx,ny)
!call dumpm(u(36,:,:),'ulev36',nx,ny)
!call dumpm(u(40,:,:),'ulev40',nx,ny)
!call dumpm(v(30,:,:),'vlev30',nx,nyp)
!call dumpm(w(30,:,:),'wlev30',nx,ny)
!call dumpm(u(:,:,ny/2),'ueq',nz,nx)
! call dumpm(theta(30,:,:),'thetalev30',nx,ny)
! call dumpm(theta(30,:,:),'thetalev30',nx,ny)
!call dumpm(w(3,:,:),'wlev3',nx,ny)
!call dumpm(w(:,:,1),'wlat1',nzp,nx)
! call dumpm(w(4,:,:),'wlev4',nx,ny)
! call dumpm(w(6,:,:),'wlev6',nx,ny)
! call dumpm(w(8,:,:),'wlev8',nx,ny)
! call dumpm(w(10,:,:),'wlev10',nx,ny)
! call dumpm(w(14,:,:),'wlev14',nx,ny)
! call dumpm(w(18,:,:),'wlev18',nx,ny)
! call dumpm(w(22,:,:),'wlev22',nx,ny)
! call dumpm(w(26,:,:),'wlev26',nx,ny)
!call dumpm(w(68,:,:),'wlev68',nx,ny)

!CALL dumpscatterm
CALL dumptrajectories

! Loop over time steps
DO ii = 1, nstop

  istep = istep + 1
  CALL step
  
  maxw = MAXVAL(ABS(w))
  print *,'Done step ',istep
  print *,'max w = ',maxw
  CALL masspertheta
  CALL unavailPE2
  CALL globaldiags
  print *,' '
  write(44,*) istep, maxw
  write(46,*) istep, MINVAL(theta(1,:,:)), MAXVAL(theta(1,:,:))
  call locatepmin(pmin,lonmin,latmin)
  write(47,'(I5,3E14.5)') istep, pmin, lonmin, latmin
  call locatemax(w,maxw,ix,iy,iz,nx,ny,nzp)
  print *,'Max abs w = ',maxw,'  at coords ',iz, ix, iy

write(87,*) istep, maxval(abs(u(1,:,:))), maxval(abs(v(1,:,:))), maxval(abs(w(2,:,:)))


  ! Output ...
  IF (MOD(istep,noutput) == 0) THEN
    CALL potvort
    !call bwout
    !call dumpm(u(:,10,:),'ulon10',nz,ny)
    !call dumpm(v(:,10,:),'vlon10',nz,nyp)
    !call dumpm(w(:,10,:),'wlon10',nzp,ny)
    !call dumpm(rho(:,10,:),'rholon10',nz,ny)
    !call dumpm(theta(:,30,:),'thetalon30',nzp,ny)
    !call dumpm(theta(:,:,95),'thetalat95',nzp,nx)
    !call dumpm(w(:,:,95),'wlat95',nzp,nx)
    !call dumpm(theta(:,:,2),'thetalat2',nzp,nx)
    !call dumpm(w(:,:,2),'wlat2',nzp,nx)
    !call dumpm(v(:,128,:),'vlon128',nz,nyp)
    !call dumpm(theta(25,:,:),'thetalev25',nx,ny)
    !call dumpm(theta(50,:,:),'thetalev50',nx,ny)
    !call dumpm(theta(75,:,:),'thetalev75',nx,ny)
    !call dumpm(theta(25,:,:),'thetalev25',nx,ny)
    !call dumpm(u(1,:,:),'ulev1',nx,ny)
    !call dumpm(v(4,:,:),'vlev4',nx,nyp)
    !call dumpm(w(10,:,:),'wlev10',nx,ny)
    call dumpm(u(1,:,:),'ulev1',nx,ny)
    call dumpm(v(1,:,:),'vlev1',nx,nyp)
    !call dumpm(w(4,:,:),'wlev4',nx,ny)
    !call dumpm(u(18,:,:),'ulev18',nx,ny)
    !call dumpm(u(36,:,:),'ulev36',nx,ny)
    !call dumpm(u(40,:,:),'ulev40',nx,ny)
    !call dumpm(v(3,:,:),'vlev3',nx,nyp)
    ! call dumpm(w(14,:,:),'wlev14',nx,ny)
    !call dumpm(w(69,:,:),'wlev69',nx,ny)
    !call dumpm(w(:,:,ny),'wnp',nzp,nx)
    !call dumpm(w(:,:,50),'wlat50',nzp,nx)
    !call dumpm(theta(:,:,64),'thetalat64',nzp,nx)
    !call dumpm(u(:,:,ny/2),'ueq',nz,nx)
    call dumpm(u(:,:,71),'ulat71',nz,nx)
    call dumpm(v(:,:,71),'vlat71',nz,nx)
    !call dumpm(w(:,:,64),'wlat64',nzp,nx)
    !call dumpm(theta(8,:,:),'thetalev8',nx,ny)
    !call dumptrajm(xtraj,ytraj,etraj,ntrajectories,'xy')
    !call dumpm(tracers(8,:,:,1),'trcrlev8',nx,ny)
    !call dumpm(pv(8,:,:),'pvlev8',nx,nyp)
    !call dumpm(tracers(8,:,:,2),'trcr2lev8',nx,ny)
    !call dumpm(tracers(:,:,48,2),'trcr2lat48',nz,nx)
    !call dumpm(bb_ref(:,:,50),'bblat50',nzp,nx)
    ubar = 0.0d0
    thetabar = 0.0d0
    DO i = 1, nx
      ubar = ubar + u(:,i,:)
      thetabar = thetabar + theta(:,i,:)
    ENDDO
    ubar = ubar/nx
    thetabar = thetabar/nx
    !call dumpm(ubar,'ubar',nz,ny)
    !call dumpm(thetabar,'thetabar',nzp,ny)

!    CALL dumpscatterm
    CALL dumptrajectories

  ENDIF

  ! Restart dump
  IF (MOD(istep,nrestart) == 0) THEN
    CALL restartdump
  ENDIF

ENDDO


END SUBROUTINE integrate

! =====================================================

SUBROUTINE epilogue

! Final actions that are performed at the end of the integration

USE channels

IMPLICIT NONE

! Close off dump files
WRITE(chandump1,*) -1
CLOSE(chandump1)
WRITE(chandump2,*) -1, -1
CLOSE(chandump2)
CLOSE(chandumpm)


END SUBROUTINE epilogue

! =====================================================

SUBROUTINE readnml

! Read namelists

USE runtype
USE switches
USE timestep
USE channels

IMPLICIT NONE


NAMELIST /rundata/ runid, lrestart, yrestart
NAMELIST /switchdata/ slice, delta_v, monotone
NAMELIST /timedata/ dt, nstop, noutput, nrestart


OPEN(channml,FILE='endgamenml.in',DELIM='APOSTROPHE')
READ(channml,rundata)
READ(channml,switchdata)
READ(channml,timedata)
CLOSE(channml)


END SUBROUTINE readnml

! =====================================================

SUBROUTINE setupgrid

! Set up grid information

! Note longitude, latitude and area information
! is dimensionless, but height and radius are in SI units

USE grid
USE constants

IMPLICIT NONE

INTEGER :: i, im, j, jm, k, km, kp, hny
REAL*8 :: a, b, sinr, cosr, num, den, sinlat, rwa, rwb
REAL*8 :: u00, cs32ev, f1, f2, stretch, rthick, deta
REAL*8 :: lonc, latc, rc, hill, rr, pp, mu


! Special code for N48 restart
! PRINT *,'Reading N48 grid data'
! CALL readn48grid

! Horizontal grid information (model coordinates)
hny = ny/2
DO i = 1, nx
  xp(i) = (i-0.5d0)*dx
  xu(i) = (i-1)*dx
  xv(i) = (i-0.5d0)*dx
ENDDO
DO j = 1, ny
  yp(j) = (j-hny-0.5d0)*dy
  yu(j) = (j-hny-0.5d0)*dy
  yv(j) = (j-hny-1)*dy
  cosp(j) = COS(yp(j))
  cosv(j) = COS(yv(j))
  sinp(j) = SIN(yp(j))
  sinv(j) = SIN(yv(j))
ENDDO
yv(1) = -piby2
yv(nyp) = piby2
cosv(1) = 0.0
cosv(nyp) = 0.0
sinv(1) = 0.0
sinv(nyp) = 1.0

! Geographical coordinates on model grid
sinr = SIN(rotgrid)
cosr = COS(rotgrid)
DO i = 1, nx
  DO j = 1, ny
    ! u points
    sinlat = cosr*sinp(j) - sinr*cosp(j)*SIN(xu(i))
    geolatu(i,j) = ASIN(sinlat)
    num = cosr*cosp(j)*SIN(xu(i)) + sinr*sinp(j)
    den = cosp(j)*COS(xu(i))
    geolonu(i,j) = MODULO(ATAN2(num,den),twopi)
    ! p points
    sinlat = cosr*sinp(j) - sinr*cosp(j)*SIN(xp(i))
    geolatp(i,j) = ASIN(sinlat)
    num = cosr*cosp(j)*SIN(xp(i)) + sinr*sinp(j)
    den = cosp(j)*COS(xp(i))
    geolonp(i,j) = MODULO(ATAN2(num,den),twopi)
    ! Coriolis terms
    coriol1(i,j) = -2.0d0*rotatn*sinr*COS(xp(i))
    coriol2(i,j) =  2.0d0*rotatn*(cosr*cosp(j) + sinr*sinp(j)*SIN(xp(i)))
    coriol3(i,j) =  2.0d0*rotatn*(cosr*sinp(j) - sinr*cosp(j)*SIN(xp(i)))
  ENDDO
  DO j = 1,nyp
    ! v points
    sinlat = cosr*sinv(j) - sinr*cosv(j)*SIN(xv(i))
    geolatv(i,j) = ASIN(sinlat)
    num = cosr*cosv(j)*SIN(xv(i)) + sinr*sinv(j)
    den = cosv(j)*COS(xv(i))
    geolonv(i,j) = MODULO(ATAN2(num,den),twopi)
  ENDDO
ENDDO

! Grid cell areas
areatot = 0.0d0
DO j = 1, ny
  ! area(j) = dx*(sinv(j+1) - sinv(j))
  area(j) = dx*dy*cosp(j)
  areatot = areatot + area(j)
ENDDO
areatot = areatot*nx

!! Vertical coordinate
!etaw(1) = 0.0
!DO k = 1, nz
!  ! Uniform grid
!  etap(k) = domain*((k - 0.5d0)/nz)
!  etaw(k+1) = domain*((1.0d0*k)/nz)
!  ! Quadratically stretched grid
!  !etap(k) = domain*((k - 0.5d0)/nz)**2
!  !etaw(k+1) = domain*((1.0d0*k)/nz)**2
!ENDDO
!! Grid with stretching factor stretch^2
!Set ratio of finest to coarsest layer thickness
!rthick = 10.0d0
!stretch = SQRT(rthick**(1.0d0/(nz-1)))
!etaw(1) = 0.0d0
!! stretch = 1.03d0
!deta = (stretch - 1.0d0)*domain/(stretch**(2*nz) - 1.0d0)
!DO k = 1, nz
!  etap(k) = etaw(k) + deta
!  deta = deta*stretch
!  etaw(k+1) = etap(k) + deta
!  deta = deta*stretch
!ENDDO
!print *,'*** stretched grid ***  stretching factor ',stretch**2
! UMJS formula
mu = 15.0d0
etaw(1) = 0.0d0
DO k = 1, nz
  etap(k)   = domain*(SQRT(mu*((k-0.5d0)/nz)**2 + 1.0d0) - 1.0d0)/(SQRT(mu + 1.0d0) - 1.0d0)
  etaw(k+1) = domain*(SQRT(mu*((k-0.0d0)/nz)**2 + 1.0d0) - 1.0d0)/(SQRT(mu + 1.0d0) - 1.0d0)
ENDDO
print *,'*** UJMS vertical grid ***'
print *,'etap = ',etap
print *,'etaw = ',etaw
print *,'Ratio of top to bottom layer thicknesses = ',(etaw(nzp)-etaw(nz))/(etaw(2)-etaw(1))

! Shape of geoid
DO j = 1, ny
  DO i = 1, nx
    pp = geolatp(i,j)
    rgeoid(i,j) = 1.0/SQRT((COS(pp)/requat)**2 + (SIN(pp)/rpole)**2)
  ENDDO
ENDDO

! Surface height
!print *,'*** Surface height set for baroclinic instability test ***'
DO j = 1, ny
  DO i = 1, nx
    rsurf(i,j) = rgeoid(i,j)
    u00 = 35.0d0
    cs32ev = (COS((1.0d0 - 0.252d0)*piby2))**1.5d0
    f1 = 10.0d0/63.0d0 - 2.0d0*(sinp(j)**6)*(cosp(j)**2 + 1.0d0/3.0d0)
    f2 = 1.6d0*(cosp(j)**3)*(sinp(j)**2 + 2.0d0/3.0d0) - 0.25d0*pi
!    rsurf(i,j) = rsurf(i,j) + u00*cs32ev*(f1*u00*cs32ev + f2*rearth*rotatn)/gravity
  ENDDO
ENDDO
!print *,'*** Pointy mountain added ***'
!lonc = 1.5d0*pi
!latc = pi/6.0d0
!rc = pi/9.0d0
!hill = 2000.0d0
!DO j = 1, ny
!  DO i = 1, nx
!    rr = SQRT(min(rc**2,(geolonp(i,j) - lonc)**2 + (geolatp(i,j) - latc)**2))
!    rsurf(i,j) = rsurf(i,j) + hill*(1.0d0 - rr/rc)
!  ENDDO
!ENDDO
!print *,'*** COS**2 mountain added ***'
!lonc = 0.5d0*pi
!latc = -pi/6.0d0
!rc = pi/9.0d0
!hill = 3000.0d0
!DO j = 1, ny
!  DO i = 1, nx
!    rr = SQRT(min(rc**2,(geolonp(i,j) - lonc)**2 + (geolatp(i,j) - latc)**2))
!    rsurf(i,j) = rsurf(i,j) + hill*COS(piby2*rr/rc)**2
!  ENDDO
!ENDDO
!print *,'*** Gaussian mountain added ***'
!lonc = 0.5d0*pi
!latc = 0.0d0
!rc = 0.3d0
!hill = 1000.0d0
!DO j = 1, ny
!  DO i = 1, nx
!    rr = (geolonp(i,j) - lonc)**2 + (geolatp(i,j) - latc)**2
!    rsurf(i,j) = rsurf(i,j) + hill*EXP(-rr/(rc*rc))
!  ENDDO
!ENDDO

call dumpm(rsurf-rearth,'zsurf',nx,ny)

! Height in all grid cells
DO j = 1, ny
  DO i = 1, nx
    DO k = 1, nz
      a = etap(k)/domain
      b = 1.0d0 - a
      rp(k,i,j) = a*(domain + rgeoid(i,j)) + b*rsurf(i,j)
    ENDDO
  ENDDO
ENDDO

! Coefficients for interpolation from w levels to p levels
DO k = 1, nz
  kp = k + 1
  above1(k) = (etap(k) - etaw(k))/(etaw(kp) - etaw(k))
  below1(k) = 1.0d0 - above1(k)
ENDDO
! and from p levels to w levels
DO k = 2, nz
  km = k - 1
  above2(k) = (etaw(k) - etap(km))/(etap(k) - etap(km))
  below2(k) = 1.0d0 - above2(k)
ENDDO
! Boundary values should never be used
above2(1) = 0.0d0
below2(1) = 0.0d0
above2(nzp) = 0.0d0
below2(nzp) = 0.0d0

! r at w points.
DO j = 1, ny
  DO i = 1, nx
    rw(1,i,j) = rsurf(i,j)
    DO k = 2, nz
      km = k - 1
      rw(k,i,j) = below2(k)*rp(km,i,j) + above2(k)*rp(k,i,j)
    ENDDO
    rw(nzp,i,j) = rgeoid(i,j) + domain
  ENDDO
ENDDO

! Areas of all cell faces
! Faces below cells
DO j = 1, ny
  DO i = 1, nx
    DO k = 1, nzp
      areab(k,i,j) = rw(k,i,j)*rw(k,i,j)*area(j)
    ENDDO
  ENDDO
ENDDO
! Faces to West of cells
DO j = 1, ny
  DO i = 1, nx
    im = i - 1
    IF (i == 1) im = nx
    DO k = 1, nz
      kp = k + 1
      rwa = 0.5d0*(rw(kp,im,j) + rw(kp,i,j))
      rwb = 0.5d0*(rw(k, im,j) + rw(k, i,j))
      areaw(k,i,j) = (rwa*rwa - rwb*rwb)*0.5d0*dy
    ENDDO
  ENDDO
ENDDO
! Faces to South of cells
DO j = 2, ny
  jm = j - 1
  DO i = 1, nx
    DO k = 1, nz
      kp = k + 1
      rwa = 0.5d0*(rw(kp,i,jm) + rw(kp,i,j))
      rwb = 0.5d0*(rw(k, i,jm) + rw(k, i,j))
      areas(k,i,j) = (rwa*rwa - rwb*rwb)*0.5d0*cosv(j)*dx
    ENDDO
  ENDDO
ENDDO
areas(:,:,1) = 0.0d0
areas(:,:,nyp) = 0.0d0

! Volumes of all cells
DO j = 1, ny
  DO i = 1, nx
    DO k = 1, nz
      kp = k + 1
      rwa = rw(kp,i,j)
      rwb = rw(k, i,j)
      volume(k,i,j) = (rwa*rwa*rwa - rwb*rwb*rwb)*area(j)/3.0d0
    ENDDO
  ENDDO
ENDDO

! -------------------------------------------------


END SUBROUTINE setupgrid

! =====================================================

SUBROUTINE initial

! Set up initial conditions

USE  constants
USE timestep
USE state

IMPLICIT NONE


! Initial state
! CALL iniacoustic
! CALL inisbr
! CALL inibw
CALL iniumjs
! CALL inibubble
! CALL readn48data

! Compute eta coordinate vertical velocity
CALL findetadot(u0,v0,w0,etadot0)

! Allocate space for tracers and initialize them
CALL initracers

! Allocate space for trajectories and initialize them
CALL initrajectories

! First guess for next time level
rho = rho0
theta = theta0
u = u0
v = v0
w = w0
etadot = etadot0
IF (ntracers > 0) tracers = tracers0

! Initialize step count
istep = 0


END SUBROUTINE initial

! =======================================================

SUBROUTINE setphi

! Set up geopotential field

USE constants
USE state

IMPLICIT NONE
INTEGER :: k

! *** Should be consistent with finphicol ***

! Inverse square law for a spherical earth
! *** Note care is needed in defining the vertical coordinate
! to ensure that gravity acts only in the vertical ***
phi = gravity*rearth*(1.0d0 - rearth/rp)

! Constant gravitational acceleration
! phi = gravity*(rp - rearth)

! Spheroidal planet with constant vertical gradient
! DO k = 1, nz
!   phi(k,:,:) = gravity*(rp(k,:,:) - rgeoid)
! ENDDO


END SUBROUTINE setphi

! =====================================================

SUBROUTINE step

! Computations for one time step

USE timestep
USE refstate
USE state
USE work
USE increments
USE switches
USE trajectories, ONLY : ntrajectories

IMPLICIT NONE

INTEGER, PARAMETER :: nouter = 4, ninner = 1
INTEGER :: iinner, iouter, ng
REAL*8 :: usp(nz,nx), unp(nz,nx), wsp(nz), wnp(nz)

! Variables used for diagnostics
real*8 :: trurho(nz,nx,ny), trutheta(nzp,nx,ny), truu(nz,nx,ny), truv(nz,nx,nyp), truw(nzp,nx,ny)
real*8 :: errrho(nz,nx,ny), errtheta(nzp,nx,ny), erru(nz,nx,ny), errv(nz,nx,nyp), errw(nzp,nx,ny)
integer :: count
real*8 :: maxrhsw(20)


! Reset delta_v to unregularized value
dv3d = REAL(delta_v)

! Use old time level solution as reference profile for Helmholtz problem
rho_ref = rho0
theta_ref = theta0
CALL findexner(rho_ref,theta_ref,exner_ref)


! Compute RHS terms on the model grid involving old time level
CALL rhsgrid


! First guess for new time level solution is old time level solution

ldiag = .false.
count = 0

! Outer loop
DO iouter = 1, nouter

  print *,'Outer iteration ',iouter
  ! Outer loop update
  CALL up_outer

  ! Build coefficients for Helmholtz operator
  CALL build_helm

  ! Inner loop
  DO iinner = 1, ninner
  
    count = count + 1

    print *,'Inner iteration ',iinner  
    ! Inner loop update
    CALL up_inner

    ! Build RHS of Helmholtz problem from individual RHSs
    CALL build_rhs_helm
    call dumprms(rhsrho,'rhsrho',nz,nx,ny)
    call dumprms(rhstheta,'rhstheta',nzp,nx,ny)
    call dumprms(rhsu,'rhsu',nz,nx,ny)
    call dumprms(rhsv,'rhsv',nz,nx,nyp)
    call dumprms(rhsw,'rhsw',nzp,nx,ny)
    call dumprms(rhs_helm,'rhshelm',nz,nx,ny)
    maxrhsw(count) = MAXVAL(ABS(rhsw))
    
    ! Solve Helmholtz
    ng = p - 1
    CALL mgsolve(exner_inc,rhs_helm,ng)

    ! Backsubstitution to find all increments
    CALL backsub

    ! Increment solution
    rho = rho + rho_inc
    theta = theta + theta_inc
    u = u + u_inc
    v = v + v_inc
    w = w + w_inc
    
    ! Update polar values of v
    CALL polar(u,w,usp,v(1,1,1),wsp,unp,v(1,1,nyp),wnp)

    ! Update top and bottom boundary values of w
    CALL wtopnbottom

    ! Compute up-to-date etadot
    CALL findetadot(u,v,w,etadot)

  ENDDO

ENDDO

write(45,888) istep,maxrhsw(1:count)
888 format(i5,20e16.4)

! Jump convergence diagnostics
goto 999

! --------------------------------------------------------------------------

! Convergence diagnostics

trurho = rho
trutheta = theta
truu = u
truv = v
truw = w

!print *,'Go again with pertn...'
! theta(5,60,64) = theta(5,60,64) + 1.0d0
!call findetadot(u,v,w,etadot)

print *,'Go again with diagnostics this time'
rho = rho0
theta = theta0
u = u0
v = v0
w = w0
etadot = etadot0

count = 0

! Outer loop
DO iouter = 1, 20
  print *,'Outer iteration ',iouter
  ldiag = .true.

  ! Outer loop update
  CALL up_outer

  ! Inner loop
  DO iinner = 1, 1
    print *,'Inner iteration ',iinner
    count = count + 1

    if (ldiag) then
      !call dumpm((rho(:,46,:)-trurho(:,46,:))/rho(:,46,:),'rhoerr',nz,ny)
      !call dumpm(theta(:,46,:)-trutheta(:,46,:),'thetaerr',nzp,ny)
      !call dumpm(u(:,46,:)-truu(:,46,:),'uerr',nz,ny)
      !call dumpm(v(:,46,:)-truv(:,46,:),'verr',nz,nyp)
      !call dumpm(w(:,46,:)-truw(:,46,:),'werr',nzp,ny)
      !call dumpm(rho(:,:,64)-trurho(:,:,64),'rhoerr',nz,nx)
      !call dumpm(theta(:,:,64)-trutheta(:,:,64),'thetaerr',nzp,nx)
      !call dumpm(u(:,:,64)-truu(:,:,64),'uerr',nz,nx)
      !call dumpm(v(:,:,64)-truv(:,:,64),'verr',nz,nx)
      !call dumpm(w(:,:,64)-truw(:,:,64),'werr',nzp,nx)
      !call dumpm((rho(2,:,:)-trurho(2,:,:)),'rhoerr',nx,ny)
      !call dumpm(theta(2,:,:)-trutheta(2,:,:),'thetaerr',nx,ny)
      !call dumpm(u(2,:,:)-truu(2,:,:),'uerr',nx,ny)
      !call dumpm(v(2,:,:)-truv(2,:,:),'verr',nx,nyp)
      !call dumpm(w(2,:,:)-truw(2,:,:),'werr',nx,ny)
    endif

    ! Inner loop update (RHSs dumped in up_inner)
    CALL up_inner

    ! Build RHS of Helmholtz problem from individual RHSs
    CALL build_rhs_helm
    !call dumprms(rhsrho,'rhsrho',nz,nx,ny)
    !call dumprms(rhstheta,'rhstheta',nzp,nx,ny)
    !call dumprms(rhsu,'rhsu',nz,nx,ny)
    !call dumprms(rhsv,'rhsv',nz,nx,nyp)
    !call dumprms(rhsw,'rhsw',nzp,nx,ny)
    !call dumprms(rhs_helm,'rhshelm',nz,nx,ny)
    
    ! Solve Helmholtz
    ng = p - 1
    CALL mgsolve(exner_inc,rhs_helm,ng)

    ! Backsubstitution to find all increments
    CALL backsub

    ! Increment solution
    rho = rho + rho_inc
    theta = theta + theta_inc
    u = u + u_inc
    v = v + v_inc
    w = w + w_inc
    
    ! if (count .ge. 5) ldiag = .false.
    if (ldiag) then
      !call dumpm(rho_inc(:,64,:),'rhoinc',nz,ny)
      !call dumpm(theta_inc(:,64,:),'thetainc',nzp,ny)
      !call dumpm(u_inc(:,64,:),'uinc',nz,ny)
      !call dumpm(v_inc(:,64,:),'vinc',nz,nyp)
      !call dumpm(w_inc(:,64,:),'winc',nzp,ny)
      !call dumpm(rho_inc(:,:,95),'rhoinc',nz,nx)
      !call dumpm(theta_inc(:,:,95),'thetainc',nzp,nx)
      !call dumpm(u_inc(:,:,95),'uinc',nz,nx)
      !call dumpm(v_inc(:,:,95),'vinc',nz,nx)
      !call dumpm(w_inc(:,:,95),'winc',nzp,nx)
      !call dumpm(rho_inc(2,:,:),'rhoinc',nx,ny)
      !call dumpm(theta_inc(2,:,:),'thetainc',nx,ny)
      !call dumpm(u_inc(2,:,:),'uinc',nx,ny)
      !call dumpm(v_inc(2,:,:),'vinc',nx,nyp)
      !call dumpm(w_inc(2,:,:),'winc',nx,ny)
    endif

    ! Update polar values of v
    CALL polar(u,w,usp,v(1,1,1),wsp,unp,v(1,1,nyp),wnp)

    ! Update top and bottom boundary values of w
    CALL wtopnbottom

    ! Compute up-to-date etadot
    CALL findetadot(u,v,w,etadot)

  ENDDO

ENDDO

! --------------------------------------------------------------------------

999 continue

! --------------------------------------------------------------------------

! Advance tracers and trajectories
IF (ntracers > 0 .OR. ntrajectories > 0) THEN
  CALL advance_tracers
ENDIF

! --------------------------------------------------------------------------

! Add in forcing terms
! print *,'*** Held-Suarez forcing added ***'
! CALL heldsuarez
! CALL tlheldsuarez

! --------------------------------------------------------------------------

! Diagnose a measure of departure from hydrostatic balance
CALL diag_imbal

! --------------------------------------------------------------------------

! Copy fields ready for start of next step.
rho0 = rho
theta0 = theta
u0 = u
v0 = v
w0 = w
etadot0 = etadot
IF (ntracers > 0) tracers0 = tracers


END SUBROUTINE step

! =====================================================

SUBROUTINE advance_tracers

! Step forward tracers. Also step forward trajectories.
! This routine is called from the end of step, so both
! the old time level and new time level winds are available.

USE state
USE departure
USE trajectories

IMPLICIT NONE

INTEGER :: itracer, itraj
REAL*8 :: usp0(nz,nx), unp0(nz,nx), usp(nz,nx), unp(nz,nx), &
          esp0(nz), enp0(nz), esp(nz), enp(nz)
REAL*8, ALLOCATABLE :: uatp0(:,:,:), vatp0(:,:,:), eatp0(:,:,:), &
		       uatp(:,:,:),  vatp(:,:,:),  eatp(:,:,:)

! -------------------------------------------------------

! Compute polar values of u and v velocity components and etadot
CALL polar(u0,etadot0,usp0,v0(1,1,1),esp0,unp0,v0(1,1,nyp),enp0)
CALL polar(u ,etadot ,usp ,v(1,1,1) ,esp ,unp ,v(1,1,nyp) ,enp )

! -------------------------------------------------------

! Allocate space for p point velocities
ALLOCATE(uatp0(nz,nx,ny),vatp0(nz,nx,ny),eatp0(nz,nx,ny))
ALLOCATE(uatp(nz,nx,ny) ,vatp(nz,nx,ny) ,eatp(nz,nx,ny) )

! Obtain velocity components at p points
CALL vecatp(u0,v0,etadot0,uatp0,vatp0,eatp0)
CALL vecatp(u ,v ,etadot ,uatp ,vatp ,eatp )

! -------------------------------------------------------

IF (ntracers > 0) THEN

  ! Compute p depature points
  CALL departurep2(uatp0,vatp0,eatp0,uatp,vatp,eatp)

  ! Advect tracers
  ! (Note it is inefficient to do them one at a time because
  ! all the interpolation coefficients have to be recomputed;
  ! all at once would be much cheaper.)
  DO itracer = 1, ntracers
    CALL lagrangep(xdepp,ydepp,edepp,tracers0(1,1,1,itracer),tracers(1,1,1,itracer))
  ENDDO

  tracers0 = tracers

ENDIF

! -------------------------------------------------------

IF (ntrajectories > 0) THEN

  DO itraj = 1, ntrajectories
    CALL steptrajectory(xtraj(itraj),ytraj(itraj),etraj(itraj), &
                        uatp0,vatp0,eatp0,                      &
                        uatp ,vatp ,eatp )
  ENDDO

ENDIF

! -------------------------------------------------------

! Release space
DEALLOCATE(uatp0,vatp0,eatp0)
DEALLOCATE(uatp ,vatp ,eatp )


END SUBROUTINE advance_tracers

! =====================================================

SUBROUTINE rhsgrid

! Compute RHS terms in the grid at the old time level

USE switches
USE constants
USE timestep
USE state
USE work

IMPLICIT NONE

INTEGER :: i, ip, j, jp
REAL*8 :: temp1(nz,nx,ny), temp2(nz,nx,nyp), temp3(nzp,nx,ny)
integer :: nzm
real*8 :: a1, b1, a2, b2, uex, vex


! RHS of rho equation
IF (slice) THEN
  print *,'SLICE option is not implemented yet'
  STOP
ELSE
  ! Compute divergence
  CALL divergence(u0,v0,etadot0,temp1)
  ! Finalize computation of old time level terms in rho equation
  rhsgridrho = rho0*(1.0d0 - beta_rho*dt*temp1)
ENDIF

! RHS of theta equation
rhsgridtheta = theta0

! RHS of momentum equations
! Exner gradient
CALL findexner(rho0,theta0,temp1)
CALL grad(temp1,rhsgridu,rhsgridv,rhsgridw)

! Theta at u and v points
CALL thetaatuv(theta0,temp1,temp2)

! -Cp theta grad(exner)
rhsgridu = -cp*temp1*rhsgridu
rhsgridv = -cp*temp2*rhsgridv
rhsgridw = -cp*theta0*rhsgridw

! Geopotential gradient
CALL grad(phi,temp1,temp2,temp3)
rhsgridu = rhsgridu - temp1
rhsgridv = rhsgridv - temp2
rhsgridw = rhsgridw - temp3

! Coriolis terms
CALL coriolis(u0,v0,w0,rho0,temp1,temp2,temp3)
rhsgridu = rhsgridu + temp1
rhsgridv = rhsgridv + temp2
rhsgridw = rhsgridw + temp3

! u + beta*dt*psi
rhsgridu =      u0 + beta_u*dt*rhsgridu
rhsgridv =      v0 + beta_u*dt*rhsgridv
rhsgridw = dv3d*w0 + beta_w*dt*rhsgridw


! Need to set sensible values at top and bottom boundaries as they're needed
! for averaging to u and v points in up_outer
rhsgridw(1,:,:) = rhsgridw(2,:,:)
rhsgridw(nzp,:,:) = rhsgridw(nz,:,:)


!print *,'*** rhsgridw = - u^2 + v^2 / r  at bdies ***'
!! Boundary value of u^2 + v^2 / r using
!! linear extrapolation of u and v to the boundary.
!nzm = nz - 1
!a1 = (etaw(1) - etap(1))/(etap(2) - etap(1))
!b1 = 1.0d0 - a1
!a2 = (etaw(nzp) - etap(nzm))/(etap(nz) - etap(nzm))
!b2 = 1.0d0 - a2
!DO j = 1, ny
!  jp = j + 1
!  DO i = 1, nx
!    ip = i + 1
!    IF (i == nx) ip = 1
!    uex = a1*0.5d0*(u0(2,i,j) + u0(2,ip,j)) &
!        + b1*0.5d0*(u0(1,i,j) + u0(1,ip,j))
!    vex = a1*0.5d0*(v0(2,i,j) + v0(2,i,jp)) &
!        + b1*0.5d0*(v0(1,i,j) + v0(1,i,jp))
!    rhsgridw(1,i,j)   = -( uex**2 + vex**2 ) / rsurf(i,j)
!    uex = a2*0.5d0*(u0(nz ,i,j) + u0(nz ,ip,j)) &
!        + b2*0.5d0*(u0(nzm,i,j) + u0(nzm,ip,j))
!    vex = a2*0.5d0*(v0(nz ,i,j) + v0(nz ,i,jp)) &
!        + b2*0.5d0*(v0(nzm,i,j) + v0(nzm,i,jp))
!    rhsgridw(nzp,i,j) = -( uex**2 + vex**2 ) / (rearth + domain)
!  ENDDO
!ENDDO


END SUBROUTINE rhsgrid

! =====================================================

SUBROUTINE up_outer

! Outer loop update of RHS terms, mostly departure point
! terms

USE switches
USE state
USE departure
USE work
use refstate ! temporary for diagnostics

IMPLICIT NONE
INTEGER :: i, j, k
REAL*8, ALLOCATABLE :: temp1(:,:,:), temp2(:,:,:), temp3(:,:,:), &
                       rvatu(:,:,:), rwatu(:,:,:), &
		       ruatv(:,:,:), rwatv(:,:,:), &
		       ruatw(:,:,:), rvatw(:,:,:)
REAL*8 :: rhsusp(nz,nx), rhsunp(nz,nx), rhswsp(nz), rhswnp(nz)
REAL*8 :: sina, cosa, sind, cosd, sasd, sacd, casd, cacd, sdl, cdl, dlambda, &
          m11, m12, m13, m21, m22, m23, m31, m32, m33
	  

! Compute all departure points
CALL depart_all2

! RHS of rho equation
IF (slice) THEN
  print *,'SLICE option is not implemented yet'
  STOP
ELSE
  CALL lagrangep(xdepp,ydepp,edepp,rhsgridrho,rhsdeprho)
ENDIF


! Polar values of velocity equation RHS terms
CALL polar(rhsgridu,rhsgridw,rhsusp,rhsgridv(1,1,1),rhswsp,rhsunp,rhsgridv(1,1,nyp),rhswnp)

! RHS of u equation
! Allocate space for u point RHS
ALLOCATE(rvatu(nz,nx,ny),rwatu(nz,nx,ny),temp1(nz,nx,ny),temp2(nz,nx,ny),temp3(nz,nx,ny))
! Obtain RHS of v and w equations at u points
CALL vecatu(rhsgridv,rhsgridw,rvatu,rwatu)
! Interpolate the three vector components to departure points
CALL lagrangeu(xdepu,ydepu,edepu,rhsgridu,rvatu,rwatu,temp1,temp2,temp3)
! Rotate to arrival point coordinate system
DO  j = 1, ny
  ! Trig factors at arrival point
  sina = sinp(j)
  cosa = cosp(j)
  DO i = 1, nx
    DO k = 1, nz
      sind = SIN(ydepu(k,i,j))
      cosd = COS(ydepu(k,i,j))
      ! sasd = sina*sind
      ! sacd = sina*cosd
      ! casd = cosa*sind
      ! cacd = cosa*cosd
      dlambda = xu(i) - xdepu(k,i,j)
      sdl = SIN(dlambda)
      cdl = COS(dlambda)
      m11 = cdl
      m12 = sind*sdl
      m13 = -cosd*sdl
      ! m21 = -sina*sdl
      ! m22 = cacd + sasd*cdl
      ! m23 = casd - sacd*cdl
      ! m31 = cosa*sdl
      ! m32 = sacd - casd*cdl
      ! m33 = sasd + cacd*cdl
      rhsdepu(k,i,j) = temp1(k,i,j)*m11 + temp2(k,i,j)*m12 + temp3(k,i,j)*m13
    ENDDO
  ENDDO
ENDDO
! Deallocate space
DEALLOCATE(rvatu,rwatu,temp1,temp2,temp3)

! RHS of v equation
! Allocate space for v point RHS
ALLOCATE(ruatv(nz,nx,nyp),rwatv(nz,nx,nyp),temp1(nz,nx,nyp),temp2(nz,nx,nyp),temp3(nz,nx,nyp))
! Obtain RHS of u and w equations at v points
CALL vecatv(rhsgridu,rhsgridw,ruatv,rwatv)
! Polar values of ruatv and rwatv
DO i = 1, nx
  ruatv(:,i,1) = rhsusp(:,i)
  ruatv(:,i,nyp) = rhsunp(:,i)
  rwatv(:,i,1) = rhswsp(:)
  rwatv(:,i,nyp) = rhswnp(:)
ENDDO
! Interpolate the three vector components to departure points
CALL lagrangev(xdepv,ydepv,edepv,ruatv,rhsgridv,rwatv,temp1,temp2,temp3)
! Rotate to arrival point coordinate system
DO  j = 2, ny
  ! Trig factors at arrival point
  sina = sinv(j)
  cosa = cosv(j)
  DO i = 1, nx
    DO k = 1, nz
      sind = SIN(ydepv(k,i,j))
      cosd = COS(ydepv(k,i,j))
      sasd = sina*sind
      sacd = sina*cosd
      casd = cosa*sind
      cacd = cosa*cosd
      dlambda = xv(i) - xdepv(k,i,j)
      sdl = SIN(dlambda)
      cdl = COS(dlambda)
      ! m11 = cdl
      ! m12 = sind*sdl
      ! m13 = -cosd*sdl
      m21 = -sina*sdl
      m22 = cacd + sasd*cdl
      m23 = casd - sacd*cdl
      ! m31 = cosa*sdl
      ! m32 = sacd - casd*cdl
      ! m33 = sasd + cacd*cdl
      rhsdepv(k,i,j) = temp1(k,i,j)*m21 + temp2(k,i,j)*m22 + temp3(k,i,j)*m23
    ENDDO
  ENDDO
ENDDO
! Set polar values to zero
rhsdepv(:,:,1) = 0.0d0
rhsdepv(:,:,nyp) = 0.0d0
! Deallocate space
DEALLOCATE(ruatv,rwatv,temp1,temp2,temp3)

! RHS of w and theta equations
! Allocate space for w point RHS
ALLOCATE(ruatw(nzp,nx,ny),rvatw(nzp,nx,ny),temp1(nzp,nx,ny),temp2(nzp,nx,ny),temp3(nzp,nx,ny))
! Obtain RHS of u and v equations at w points
CALL vecatw(rhsgridu,rhsgridv,ruatw,rvatw)
! Interpolate the three vector components and theta to departure points
! CALL lagrangew(xdepw,ydepw,edepw,ruatw,rvatw,rhsgridw,rhsgridtheta,temp1,temp2,temp3,rhsdeptheta)
CALL lagrangewd(xdepw,ydepw,edepw,ruatw,rvatw,rhsgridw,rhsgridtheta, &
               temp1,temp2,temp3,rhsdeptheta,dthetabardr_ref)
! Rotate to arrival point coordinate system
DO  j = 1, ny
  ! Trig factors at arrival point
  sina = sinp(j)
  cosa = cosp(j)
  DO i = 1, nx
    DO k = 2, nz
      sind = SIN(ydepw(k,i,j))
      cosd = COS(ydepw(k,i,j))
      sasd = sina*sind
      sacd = sina*cosd
      casd = cosa*sind
      cacd = cosa*cosd
      dlambda = xp(i) - xdepw(k,i,j)
      sdl = SIN(dlambda)
      cdl = COS(dlambda)
      ! m11 = cdl
      ! m12 = sind*sdl
      ! m13 = -cosd*sdl
      ! m21 = -sina*sdl
      ! m22 = cacd + sasd*cdl
      ! m23 = casd - sacd*cdl
      m31 = cosa*sdl
      m32 = sacd - casd*cdl
      m33 = sasd + cacd*cdl
      rhsdepw(k,i,j) = temp1(k,i,j)*m31 + temp2(k,i,j)*m32 + temp3(k,i,j)*m33
    ENDDO
  ENDDO
ENDDO
! Set bottom and top values to zero
rhsdepw(1,:,:) = 0.0d0
rhsdepw(nzp,:,:) = 0.0d0
! Deallocate space
DEALLOCATE(ruatw,rvatw,temp1,temp2,temp3)

! LHS of theta equation must be updated in outer loop
! for optimal convergence
! LHS of theta equation is just latest guess for theta itself
rhstheta = rhsdeptheta - theta

! IF using slice then LHS of rho equation must be updated in
! outer loop for optimal convergence
! LHS of rho equation is just the latest guess for rho itself
IF (slice) THEN
  rhsrho = rhsdeprho - rho
ENDIF


END SUBROUTINE up_outer

! =====================================================

SUBROUTINE up_inner

! Inner loop update of RHS terms

USE switches
USE constants
USE timestep
USE state
USE work

use refstate


IMPLICIT NONE
REAL*8 :: lhsu(nz,nx,ny), lhsv(nz,nx,nyp), lhsw(nzp,nx,ny), &
          temp1(nz,nx,ny), temp2(nz,nx,nyp), temp3(nzp,nx,ny)


! Rho equation
IF (slice) THEN
  ! LHS and hence RHS is updated in up_outer
ELSE
  ! Divergence field
  CALL divergence(u,v,etadot,temp1)
  ! LHS terms
  temp1 = rho*(1.0d0 + alpha_rho*dt*temp1)
  rhsrho = rhsdeprho - temp1
ENDIF


! Theta equation
! LHS and hence RHS are updated in outer loop


! velocity equations
! Exner gradient
CALL findexner(rho,theta,temp1)
CALL grad(temp1,lhsu,lhsv,lhsw)

! Theta at u and v points
CALL thetaatuv(theta,temp1,temp2)

! -Cp theta grad(exner)
lhsu = -cp*temp1*lhsu
lhsv = -cp*temp2*lhsv
lhsw = -cp*theta*lhsw

! Geopotential gradient
CALL grad(phi,temp1,temp2,temp3)
lhsu = lhsu - temp1
lhsv = lhsv - temp2
lhsw = lhsw - temp3

! Coriolis terms
CALL coriolis(u,v,w,rho,temp1,temp2,temp3)
lhsu = lhsu + temp1
lhsv = lhsv + temp2
lhsw = lhsw + temp3

! u - alpha*dt*psi
lhsu =      u - alpha_u*dt*lhsu
lhsv =      v - alpha_u*dt*lhsv
lhsw = dv3d*w - alpha_w*dt*lhsw
! Tidy up top and bottom levels (shouldn't be necessary)
lhsw(1,:,:) = 0.0d0
lhsw(nzp,:,:) = 0.0d0
! And polar values
lhsv(:,:,1) = 0.0d0
lhsv(:,:,nyp) = 0.0d0

! Finalize RHS terms of velocity equations
rhsu = rhsdepu - lhsu
rhsv = rhsdepv - lhsv
rhsw = rhsdepw - lhsw

! Dump fields for convergence diagnostics if required
if (ldiag) then
  !call dumpm(rhsrho(:,64,:),'rhsrho',nz,ny)
  !call dumpm(rhstheta(:,64,:),'rhstheta',nzp,ny)
  !call dumpm(rhsu(:,64,:),'rhsu',nz,ny)
  !call dumpm(rhsv(:,64,:),'rhsv',nz,nyp)
  !call dumpm(rhsw(:,64,:),'rhsw',nzp,ny)
  !call dumpm(rhsrho(:,:,95),'rhsrho',nz,nx)
  !call dumpm(rhstheta(:,:,95),'rhstheta',nzp,nx)
  !call dumpm(rhsu(:,:,95),'rhsu',nz,nx)
  !call dumpm(rhsv(:,:,95),'rhsv',nz,nx)
  !call dumpm(rhsw(:,:,95),'rhsw',nzp,nx)
  !call dumpm(rhsrho(2,:,:),'rhsrho',nx,ny)
  !call dumpm(rhstheta(2,:,:),'rhstheta',nx,ny)
  !call dumpm(rhsu(2,:,:),'rhsu',nx,ny)
  !call dumpm(rhsv(2,:,:),'rhsv',nx,nyp)
  !call dumpm(rhsw(2,:,:),'rhsw',nx,ny)
endif


END SUBROUTINE up_inner

! =====================================================

SUBROUTINE build_rhs_helm

! Build RHS of Helmholtz problem from individual equation RHSs

USE constants
USE timestep
USE state
USE refstate
USE work

IMPLICIT NONE
INTEGER :: i, im, ip, j, jm, jp, k, kp
REAL*8 :: const, const2, rhobar
REAL*8 :: temp1(nz,nx,ny), temp2(nz,nx,nyp), temp3(nzp,nx,ny), div(nz,nx,ny)


! Useful constants
const = 1.0d0/(alpha_w*dt*cp)
const2 = alpha_rho*dt


! Leading order inclusion of Coriolis terms in Helmholtz problem
!  print *,'No Coriolis correction in Helmholtz'
!  print *,'Deep Coriolis correction in Helmholtz'
!  CALL coriolis(rhsu,rhsv,rhsw,rho,temp1,temp2,temp3)
!  rhsu = rhsu + alpha_u*dt*temp1
!  rhsv = rhsv + alpha_u*dt*temp2
!  rhsw = rhsw + alpha_w*dt*temp3
  print *,'Shallow Coriolis correction in Helmholtz'
  ! Borrow q16 array temporarily
  q16 = 0.0d0 ! Temporary use of array
  CALL coriolis(rhsu,rhsv,q16,rho,temp1,temp2,temp3)
  rhsu = rhsu + alpha_u*dt*temp1
  rhsv = rhsv + alpha_u*dt*temp2


! RHS of equation (16) divided by (alpha_w*dt*cp*theta_ref)
! Note this is needed for back substitution
q16 = (const*rhsw - dexnerdr_ref*rhstheta)/theta_ref


DO k = 1, nz
  kp = k + 1
  ! Average rhstheta to p points
  temp1(k,:,:) = above1(k)*rhstheta(kp,:,:) + below1(k)*rhstheta(k,:,:)
  ! and evaluate D_1 (R23) ( = C(Q16) )
  rhs_helm(k,:,:) = cabove(k,:,:)*q16(kp,:,:) - cbelow(k,:,:)*q16(k,:,:)
ENDDO


! Subtract off RHS of equation (19)
temp1 = rhsrho + rho_ref*temp1/thetabar_ref
rhs_helm = rhs_helm - temp1


! Divergence of rho_ref times RHS of u and v equations.
! First compute 'mass fluxes'
DO j = 1, ny
  DO i = 1, nx
    im = i - 1
    IF (i == 1) im = nx
    DO k = 1,nz
      rhobar = 0.5d0*(rho_ref(k,i,j) + rho_ref(k,im,j))
      temp1(k,i,j) = rhobar*areaw(k,i,j)*rhsu(k,i,j)
    ENDDO
  ENDDO
ENDDO
! Zero northward mass flux at poles
temp2(:,:,1) = 0.0d0
temp2(:,:,nyp) = 0.0d0
DO j = 2, ny
  jm = j - 1
  DO i = 1, nx
    DO k = 1,nz
      rhobar = 0.5d0*(rho_ref(k,i,j) + rho_ref(k,i,jm))
      temp2(k,i,j) = rhobar*areas(k,i,j)*rhsv(k,i,j)
    ENDDO
  ENDDO
ENDDO
! Now compute horizontal divergence
DO k = 1, nz
  DO j = 1, ny
    jp = j + 1
    DO i = 1, nx
      ip = i + 1
      IF (i == nx) ip = 1
      div(k,i,j) = ( (temp1(k,ip,j) - temp1(k,i,j))   &
                   + (temp2(k,i,jp) - temp2(k,i,j)) ) &
		   / volume(k,i,j)
    ENDDO
  ENDDO
ENDDO

rhs_helm = rhs_helm + const2*div


END SUBROUTINE build_rhs_helm

! =====================================================

SUBROUTINE backsub

! Perform back substitution to compute increments to all variables


USE switches
USE constants
USE timestep
USE refstate
USE work
USE increments

use state

IMPLICIT NONE

INTEGER :: k, kp
REAL*8 :: temp1(nz,nx,ny), temp2(nz,nx,nyp), temp2d(nx,ny), &
          gradpx(nz,nx,ny), gradpy(nz,nx,nyp), gradpz(nzp,nx,ny) 
REAL*8 :: const
REAL*8, ALLOCATABLE :: uatw(:,:,:), vatw(:,:,:), uincatw(:,:,:), vincatw(:,:,:)
REAL*8, ALLOCATABLE :: uatp(:,:,:), vatp(:,:,:), watp(:,:,:), &
                       uincatp(:,:,:), vincatp(:,:,:), wincatp(:,:,:)


! Compute gradient of exner increment
CALL grad(exner_inc,gradpx,gradpy,gradpz)

! Theta_ref at u and v points
CALL thetaatuv(theta_ref,temp1,temp2)
! Hence compute  - cp * theta_ref * grad(exner_inc)
gradpx = - cp*temp1*gradpx
gradpy = - cp*temp2*gradpy

! Increment to u and v
const = alpha_u*dt
u_inc = rhsu + const*gradpx
v_inc = rhsv + const*gradpy

! Increment to w
const = alpha_w*dt*cp
w_inc = const*theta_ref*bb_ref*(q16 - gradpz)
! Make sure it's zero at top and bottom boundaries
w_inc(1,:,:) = 0.0d0
w_inc(nzp,:,:) = 0.0d0

! Increment to theta
const = alpha_x*dt
theta_inc = rhstheta - const*w_inc*dthetabardr_ref

! Increment to rho using equation of state
print *,'rho_inc via equation of state'
DO k = 1, nz
  kp = k + 1
  temp2d = (above1(k)*theta_inc(kp,:,:) + below1(k)*theta_inc(k ,:,:))/thetabar_ref(k ,:,:)
  rho_inc(k,:,:) = rho_ref(k,:,:)*(onemkbyk*exner_inc(k,:,:)/exner_ref(k,:,:) - temp2d)
ENDDO


! Alternative increment to rho using mass equation
! (Might be problematic at high resolution if Helmholtz
! is not solved accurately)
!IF (slice) THEN
!  PRINT *,'Check that the form of massdiv is appropriate'
!  STOP
!ELSE
!  print *,'rho_inc via mass equation'
!  ! Approximate alpha * divergence of (rho_ref * u_inc)
!  CALL massdiv2(u_inc,v_inc,w_inc,temp1)
!  rho_inc = rhsrho - dt*temp1
!ENDIF


! `Accelerator' terms - approximate allowance for vertical change
! in departure points within inner loop calculation

! 1. Correct rhstheta (because otherwise rhstheta is only updated in
! the outer loop)
rhstheta = 0.0d0

! 2. Correct rhsrho or rhsdeprho to allow for vertical change in
! departure points
ALLOCATE(uincatp(nz,nx,ny), vincatp(nz,nx,ny), wincatp(nz,nx,ny))
IF (slice) THEN
  print *,'Need to correct rhsrho in backsub'
ELSE
  CALL vecatp(u_inc,v_inc,w_inc,uincatp,vincatp,wincatp)
  rhsdeprho = rhsdeprho - alpha_x*dt*wincatp*drhobardr_ref
ENDIF
DEALLOCATE(uincatp,vincatp,wincatp)


END SUBROUTINE backsub

! =====================================================

SUBROUTINE findetadot(u,v,w,etadot)

! Compute etadot from 3D wind field

USE grid

IMPLICIT NONE

REAL*8, INTENT(IN) :: u(nz,nx,ny), v(nz,nx,nyp), w(nzp,nx,ny)
REAL*8, INTENT(OUT) :: etadot(nzp,nx,ny)

INTEGER :: i, im, ip, i1, j, jm, jp, k, km, hnx
REAL*8 :: h1, h2, h3, d1xi3, d2xi3, d3eta
REAL*8 :: tempu(nz,nx,ny), tempv(nz,nx,nyp)


hnx = nx/2

! First compute the contributions to v . grad(r)
DO k = 1, nz
  DO i = 1, nx
    im = i - 1
    IF (i == 1) im = nx
    ! For polar points i1 is index of the point on the
    ! other side of the pole
    i1 = MODULO(i - hnx - 1, nx) + 1
    DO j = 1, ny
      h1 = 0.5d0*(rp(k,i,j) + rp(k,im,j))*cosp(j)
      d1xi3 = (rp(k,i,j) - rp(k,im,j))/dx
      tempu(k,i,j) = u(k,i,j)*d1xi3/h1
    ENDDO
    j = 1
    jm = j
    h2 = 0.5d0*(rp(k,i,j) + rp(k,i1,jm))
    d2xi3 = (rp(k,i,j) - rp(k,i1,jm))/dy
    tempv(k,i,j) = v(k,i,j)*d2xi3/h2
    DO j = 2, ny
      jm = j - 1
      h2 = 0.5d0*(rp(k,i,j) + rp(k,i,jm))
      d2xi3 = (rp(k,i,j) - rp(k,i,jm))/dy
      tempv(k,i,j) = v(k,i,j)*d2xi3/h2
    ENDDO
    j = nyp
    jm = j - 1
    h2 = 0.5d0*(rp(k,i,jm) + rp(k,i1,jm))
    d2xi3 = (rp(k,i1,jm) - rp(k,i,jm))/dy
    tempv(k,i,j) = v(k,i,j)*d2xi3/h2
  ENDDO
ENDDO

! Now average to w points and compute etadot
h3 = 1.0d0
DO k = 2, nz
  km = k - 1
  DO j = 1, ny
    jp = j + 1
    DO i = 1, nx
      ip = i + 1
      IF (i == nx) ip = 1
      d3eta = (etap(k) - etap(km))/(rp(k,i,j) - rp(km,i,j))
      etadot(k,i,j) = (w(k,i,j)/h3                             &
                    - 0.5d0*(                                  &
                   above2(k)*(tempu(k ,i,j) + tempu(k ,ip,j)   &
			    + tempv(k ,i,j) + tempv(k ,i,jp))  &
	         + below2(k)*(tempu(km,i,j) + tempu(km,ip,j)   &
			    + tempv(km,i,j) + tempv(km,i,jp))  &
                               ) )*d3eta
    ENDDO
  ENDDO
ENDDO

! Top and bottom boundaries
etadot(1,:,:) = 0.0d0
etadot(nzp,:,:) = 0.0d0


END SUBROUTINE findetadot

! ===============================================================

SUBROUTINE inidep

! Initialize departure points to equal arival points

USE grid
USE departure

IMPLICIT NONE

INTEGER :: i, j, k

DO j = 1, ny
  DO k = 1, nz
    xdepu(k,:,j) = xu(:)
    xdepv(k,:,j) = xv(:)
    xdepw(k,:,j) = xp(:)
    xdepp(k,:,j) = xp(:)
  ENDDO
  xdepw(nz,:,j) = xp(:)
ENDDO

DO i = 1, nx
  DO k = 1, nz
    ydepu(k,i,:) = yu(:)
    ydepv(k,i,:) = yv(:)
    ydepw(k,i,:) = yp(:)
    ydepp(k,i,:) = yp(:)
  ENDDO
  ydepw(nzp,i,:) = yp(:)
ENDDO

DO j = 1, ny
  DO i = 1, nx
    edepu(:,i,j) = etap(:)
    edepv(:,i,j) = etap(:)
    edepw(:,i,j) = etaw(:)
    edepp(:,i,j) = etap(:)
  ENDDO
ENDDO


END SUBROUTINE inidep

! ===============================================================

SUBROUTINE depart_all2

! Master subroutine controlling calculation of departure points
! This routine does a quasi-shallow rotation matrix treatment.
! To use the full rotation matrix treatment call depart_all instead
! (now deleted).

USE state

IMPLICIT NONE

INTEGER :: i
REAL*8 :: usp0(nz,nx), unp0(nz,nx), usp(nz,nx), unp(nz,nx), &
          esp0(nz), enp0(nz), esp(nz), enp(nz)
REAL*8, ALLOCATABLE :: vatu0(:,:,:), eatu0(:,:,:), vatu(:,:,:), eatu(:,:,:), &
                       uatv0(:,:,:), eatv0(:,:,:), uatv(:,:,:), eatv(:,:,:), &
		       uatw0(:,:,:), vatw0(:,:,:), uatw(:,:,:), vatw(:,:,:), &
		       uatp0(:,:,:), vatp0(:,:,:), eatp0(:,:,:), &
		       uatp(:,:,:),  vatp(:,:,:),  eatp(:,:,:)


! NOTE The amount of calculation could be reduced at the expense
! of storing several 3D arrays by saving some averaged velocity 
! values from the previous timestep

! -------------------------------------------------------

! Compute polar values of u and v velocity components and etadot
CALL polar(u0,etadot0,usp0,v0(1,1,1),esp0,unp0,v0(1,1,nyp),enp0)
CALL polar(u ,etadot ,usp ,v(1,1,1) ,esp ,unp ,v(1,1,nyp) ,enp )

! -------------------------------------------------------

! Allocate space for u point velocities
ALLOCATE(vatu0(nz,nx,ny),eatu0(nz,nx,ny),vatu(nz,nx,ny),eatu(nz,nx,ny))

! Obtain velocity components at u points
CALL vecatu(v0,etadot0,vatu0,eatu0)
CALL vecatu(v ,etadot ,vatu ,eatu )

! And hence compute u depature points
CALL departureu2(u0,vatu0,eatu0,u,vatu,eatu)

DEALLOCATE(vatu0,eatu0,vatu,eatu)

! -------------------------------------------------------

! Allocate space for v point velocities
ALLOCATE(uatv0(nz,nx,nyp),eatv0(nz,nx,nyp),uatv(nz,nx,nyp),eatv(nz,nx,nyp))

! Obtain velocity components at v points
CALL vecatv(u0,etadot0,uatv0,eatv0)
CALL vecatv(u ,etadot ,uatv ,eatv )
! Include polar values of u and etadot
DO i = 1, nx
  uatv0(:,i,1) = usp0(:,i)
  uatv0(:,i,nyp) = unp0(:,i)
  uatv(:,i,1) = usp(:,i)
  uatv(:,i,nyp) = unp(:,i)
  eatv0(:,i,1) = esp0(:)
  eatv0(:,i,nyp) = enp0(:)
  eatv(:,i,1) = esp(:)
  eatv(:,i,nyp) = enp(:)
ENDDO

! And hence compute v depature points
CALL departurev2(uatv0,v0,eatv0,uatv,v,eatv)

DEALLOCATE(uatv0,eatv0,uatv,eatv)

! -------------------------------------------------------

! Allocate space for w point velocities
ALLOCATE(uatw0(nzp,nx,ny),vatw0(nzp,nx,ny),uatw(nzp,nx,ny),vatw(nzp,nx,ny))

! Obtain velocity components at w points
CALL vecatw(u0,v0,uatw0,vatw0)
CALL vecatw(u ,v ,uatw ,vatw )

! And hence compute w depature points
CALL departurew2(uatw0,vatw0,etadot0,uatw,vatw,etadot)

DEALLOCATE(uatw0,vatw0,uatw,vatw)

! -------------------------------------------------------

! Allocate space for p point velocities
ALLOCATE(uatp0(nz,nx,ny),vatp0(nz,nx,ny),eatp0(nz,nx,ny))
ALLOCATE(uatp(nz,nx,ny) ,vatp(nz,nx,ny) ,eatp(nz,nx,ny) )

! Obtain velocity components at p points
CALL vecatp(u0,v0,etadot0,uatp0,vatp0,eatp0)
CALL vecatp(u ,v ,etadot ,uatp ,vatp ,eatp )

! And hence compute p depature points
CALL departurep2(uatp0,vatp0,eatp0,uatp,vatp,eatp)

DEALLOCATE(uatp0,vatp0,eatp0)
DEALLOCATE(uatp ,vatp ,eatp )


END SUBROUTINE depart_all2

! ===============================================================

SUBROUTINE polar(u,w,usp,vsp,wsp,unp,vnp,wnp)

! Determine polar values of u, v and w from u and w at nearest
! u-latitude. Also used for RHS of momentum equation.

USE grid

IMPLICIT NONE

REAL*8, INTENT(IN) :: u(nz,nx,ny), w(nzp,nx,ny)
REAL*8, INTENT(OUT) :: usp(nz,nx), vsp(nz,nx), wsp(nz), &
                       unp(nz,nx), vnp(nz,nx), wnp(nz)
REAL*8 :: tempwn(nzp), tempws(nzp)
INTEGER :: k, kp
REAL*8 :: a, b, lambda, vp


DO k = 1, nzp
  tempws(k) = SUM(w(k,:,1))/nx
  tempwn(k) = SUM(w(k,:,ny))/nx
ENDDO


DO k = 1, nz
  kp = k + 1

  ! South pole
  a = SUM(u(k,:,1)*SIN(xu(:)))
  b = SUM(u(k,:,1)*COS(xu(:)))
  lambda = ATAN2(b,-a)
  vp = (-a*COS(lambda) + b*SIN(lambda))*2.0/nx
  vsp(k,:) = vp*COS(xv(:) - lambda)
  usp(k,:) = -vp*SIN(xu(:) - lambda)
  wsp(k) = below1(k)*tempws(k) + above1(k)*tempws(kp)

  ! North pole
  a = SUM(u(k,:,ny)*SIN(xu(:)))
  b = SUM(u(k,:,ny)*COS(xu(:)))
  lambda = ATAN2(-b,a)
  vp = (a*COS(lambda) - b*SIN(lambda))*2.0/nx
  vnp(k,:) = vp*COS(xv(:) - lambda)
  unp(k,:) = vp*SIN(xu(:) - lambda)
  wnp(k) = below1(k)*tempwn(k) + above1(k)*tempwn(kp)

ENDDO


END SUBROUTINE polar

! ========================================================

SUBROUTINE vecatu(v,w,vatu,watu)

! To average v and w to u points

USE grid

IMPLICIT NONE

REAL*8, INTENT(IN) :: v(nz,nx,nyp), w(nzp,nx,ny)
REAL*8, INTENT(OUT) :: vatu(nz,nx,ny), watu(nz,nx,ny)
INTEGER :: k, kp, i, im, j, jp
REAL*8 :: a, b


DO k = 1, nz
  kp = k + 1
  a = above1(k)
  b = below1(k)
  DO j = 1, ny
    jp = j + 1
    DO i = 1, nx
      im = i - 1
      IF (i == 1) im = nx
      vatu(k,i,j) = 0.25d0*(v(k,im,j) + v(k,im,jp) + v(k,i,j) + v(k,i,jp))
      watu(k,i,j) = 0.5d0*(a*(w(kp,im,j) + w(kp,i,j)) &
                         + b*(w(k ,im,j) + w(k ,i,j)))
    ENDDO
  ENDDO
ENDDO


END SUBROUTINE vecatu

! ========================================================

SUBROUTINE vecatv(u,w,uatv,watv)

! To average u and w to v points
! (Note polar values of u must be computed by a
! separate call to polar()

USE grid

IMPLICIT NONE

REAL*8, INTENT(IN) :: u(nz,nx,ny), w(nzp,nx,ny)
REAL*8, INTENT(OUT) :: uatv(nz,nx,nyp), watv(nz,nx,nyp)
INTEGER :: k, kp, i, ip, j, jm
REAL*8 :: a, b


DO k = 1, nz
  kp = k + 1
  a = above1(k)
  b = below1(k)
  DO j = 2, ny
    jm = j - 1
    DO i = 1, nx
      ip = i + 1
      IF (i == nx) ip = 1
      uatv(k,i,j) = 0.25d0*(u(k,i,jm) + u(k,i,j) + u(k,ip,jm) + u(k,ip,j))
      watv(k,i,j) = 0.5d0*(a*(w(kp,i,jm) + w(kp,i,j)) &
                         + b*(w(k ,i,jm) + w(k ,i,j)))
    ENDDO
  ENDDO
ENDDO


END SUBROUTINE vecatv

! ========================================================

SUBROUTINE vecatw(u,v,uatw,vatw)

! To average u and v to w points

USE grid

IMPLICIT NONE

REAL*8, INTENT(IN) :: u(nz,nx,ny), v(nz,nx,nyp)
REAL*8, INTENT(OUT) :: uatw(nzp,nx,ny), vatw(nzp,nx,ny)
INTEGER :: k, km, i, ip, j, jp
REAL*8 :: a, b


DO k = 2, nz
  km = k - 1
  a = above2(k)
  b = below2(k)
  DO j = 1, ny
    jp = j + 1
    DO i = 1, nx
      ip = i + 1
      IF (i == nx) ip = 1
      uatw(k,i,j) = 0.5d0*(a*(u(k ,i,j) + u(k ,ip,j)) &
                         + b*(u(km,i,j) + u(km,ip,j)))
      vatw(k,i,j) = 0.5d0*(a*(v(k ,i,j) + v(k ,i,jp)) &
                         + b*(v(km,i,j) + v(km,i,jp)))
    ENDDO
  ENDDO
ENDDO

! Top and bottom boundaries: constant extrapolation
uatw(1,:,:) = uatw(2,:,:)
vatw(1,:,:) = vatw(2,:,:)
uatw(nzp,:,:) = uatw(nz,:,:)
vatw(nzp,:,:) = vatw(nz,:,:)


END SUBROUTINE vecatw

! ========================================================

SUBROUTINE vecatp(u,v,w,uatp,vatp,watp)

! To average u, v and w to p points

USE grid

IMPLICIT NONE

REAL*8, INTENT(IN) :: u(nz,nx,ny), v(nz,nx,nyp), w(nzp,nx,ny)
REAL*8, INTENT(OUT) :: uatp(nz,nx,ny), vatp(nz,nx,ny), watp(nz,nx,ny)
INTEGER :: k, kp, i, ip, j, jp
REAL*8 :: a, b


DO k = 1, nz
  kp = k + 1
  a = above1(k)
  b = below1(k)
  DO j = 1, ny
    jp = j + 1
    DO i = 1, nx
      ip = i + 1
      IF (i == nx) ip = 1
      uatp(k,i,j) = 0.5d0*(u(k,i,j) + u(k,ip,j))
      vatp(k,i,j) = 0.5d0*(v(k,i,j) + v(k,i,jp))
      watp(k,i,j) = a*w(kp,i,j) + b*w(k,i,j)
    ENDDO
  ENDDO
ENDDO


END SUBROUTINE vecatp

! ========================================================

SUBROUTINE thetaatuv(theta,thetau,thetav)

! To average theta to u and v points

USE grid

IMPLICIT NONE

REAL*8, INTENT(IN) ::  theta(nzp,nx,ny)
REAL*8, INTENT(OUT) :: thetau(nz,nx,ny), thetav(nz,nx,nyp)
REAL*8 :: thetap(nz,nx,ny)
INTEGER :: k, kp, i, im, j, jm
REAL*8 :: a, b

! Vertical interpolation
DO k = 1, nz
  kp = k + 1
  a = above1(k)
  b = below1(k)
  DO j = 1, ny
    DO i = 1, nx
      thetap(k,i,j) = a*theta(kp,i,j) + b*theta(k,i,j)
    ENDDO
  ENDDO
ENDDO

! Horizontal averaging
DO k = 1, nz
  DO i = 1, nx
    im = i - 1
    IF (i == 1) im = nx
    thetau(k,i,1) = 0.5d0*(thetap(k,im,1) + thetap(k,i,1))
    DO j = 2, ny
      jm = j - 1
      thetau(k,i,j) = 0.5d0*(thetap(k,im,j) + thetap(k,i,j))
      thetav(k,i,j) = 0.5d0*(thetap(k,i,jm) + thetap(k,i,j))
    ENDDO
  ENDDO
ENDDO
! Set polar values to zero
thetav(:,:,1) = 0.0d0
thetav(:,:,nyp) = 0.0d0


END SUBROUTINE thetaatuv

! ========================================================

SUBROUTINE findetacell(eta,etad,n,ix)

! Given eta values at cell top and bottom boundaries and a
! departure point eta value etad, find the index of the
! cell that etad lies in.
!
! If etad lies below the first cell then ix is set to 0
! If etad lies above the top cell then ix is set to n

IMPLICIT NONE

INTEGER, PARAMETER :: maxkount = 10
INTEGER, INTENT(IN) :: n
REAL*8, INTENT(IN) :: eta(n), etad
INTEGER, INTENT(OUT) :: ix
INTEGER :: ix1, ix2, ixnew, itkount
LOGICAL :: converged
REAL*8 :: eta1, eta2, rixnew


! First guess is based on assuming eta is linear in level number.
! A better first guess could be obtained by using more information about
! how eta is distributed
!e.g.    ix = FLOOR(n*SQRT(etad/(eta(n) - eta(1))))
ix = FLOOR(n*(etad/(eta(n) - eta(1))))
ix1 = MIN(MAX(ix,1),n-1)
ix2 = ix1 + 1

! Now iterate
converged = .FALSE.
itkount = 0
DO WHILE (.NOT. converged)
  itkount = itkount + 1
  eta1 = eta(ix1)
  eta2 = eta(ix2)
  rixnew = ((etad - eta1)*ix2 + (eta2 - etad)*ix1)/(eta2 - eta1)
  ixnew = FLOOR(rixnew)
  IF (ixnew == ix .OR. itkount .ge. maxkount) converged = .TRUE.
  ix = ixnew
  ix1 = MIN(MAX(ix,1),n-1)
  ix2 = ix1 + 1
ENDDO

ix = MIN(MAX(ix,0),n)


END SUBROUTINE findetacell

! ========================================================

SUBROUTINE departureu2(u0,vatu0,eatu0,u,vatu,eatu)

! Compute departure points for u points
! Assume a first guess is already provided

USE grid
USE constants
USE timestep
USE departure

IMPLICIT NONE

INTEGER, PARAMETER :: ndepit = 2
REAL*8, INTENT(INOUT) :: u0(nz,nx,ny), vatu0(nz,nx,ny), eatu0(nz,nx,ny), &
                         u(nz,nx,ny),  vatu(nz,nx,ny),  eatu(nz,nx,ny)
REAL*8 :: u0byr(nz,nx,ny), ubyr(nz,nx,ny), rtemp(nz,ny)
INTEGER :: i, im, j, k, id, idp, jd, jdp, kd, kdp, &
           id1, id1p, id2, id2p, hnx, idepit
REAL*8 :: sina, cosa, a1, a2, b1, b2, c1, c2, c1w, c2w, flip1, flip2, &
          b1flip, b2flip, rid, rjd, ud, vd, ed, x, y, z, &
	  sind, cosd, sasd, cacd, sdl, cdl, dlambda, &
	  urot, vrot, wrot, m11, m12, m21, m22, cosalphap1

! --------------------------------------------------------

! First divide horizontal velocities by r
DO i = 1, nx
  im = i - 1
  IF (i == 1) im = nx
  rtemp = 0.5d0*(rp(:,im,:) + rp(:,i,:))
  u0byr(:,i,:) = u0(:,i,:)/rtemp
  ubyr(:,i,:) = u(:,i,:)/rtemp
  vatu0(:,i,:) = vatu0(:,i,:)/rtemp
  vatu(:,i,:) = vatu(:,i,:)/rtemp
ENDDO


hnx = nx/2

! Loop over iterations
DO idepit = 1, ndepit

  DO  j = 1, ny

    ! Trig factors at arrival point
    sina = sinp(j)
    cosa = cosp(j)

    DO i = 1, nx
    
      im = i - 1
      IF (i == 1) im = nx

      DO k = 1, nz

        ! Determine departure cell indices based on latest estimate
        ! xdepu, ydepu, edepu
        rid = (xdepu(k,i,j) - xu(1))/dx
        id = FLOOR(rid)
        a2 = rid - id
        a1 = 1.0 - a2
        id = MODULO(id, nx) + 1
	idp = id + 1
	IF (id == nx) idp = 1
        rjd = (ydepu(k,i,j) - yu(1))/dy
        jd = FLOOR(rjd)
        b2 = rjd - jd
        b1 = 1.0 - b2
        jd = jd + 1
        jdp = jd + 1
	CALL findetacell(etap,edepu(k,i,j),nz,kd)
        IF (kd == 0) THEN
          ! Constant extrapolation below lowest u level
	  ! etadot linearly to zero
          kd = 1
          kdp = 2
          c1 = 1.0d0
          c2 = 0.0d0
	  c1w = (edepu(k,i,j) - etaw(1))/(etap(1) - etaw(1))
	  c2w = 0.0d0
        ELSEIF (kd == nz) THEN
          ! Constant extrapolation above highest u level
	  ! etadot linearly to zero
          kd = nz - 1
          kdp = nz
          c1 = 0.0d0
          c2 = 1.0d0
	  c1w = 0.0d0
	  c2w = (edepu(k,i,j) - etaw(nzp))/(etap(nz) - etaw(nzp))
        ELSE
          kdp = kd + 1
          c1 = (etap(kdp) - edepu(k,i,j))/(etap(kdp) - etap(kd))
          c2 = 1.0d0 - c1
	  c1w = c1
	  c2w = c2
	ENDIF

        ! Tricks to handle polar case
        ! Here we interpolate across the pole.
        ! (An alternative would be to use the polar values of u and v
        ! and interpolate between nearest u latitude and the pole.)
        id1 = id
        id1p = idp
        id2 = id
        id2p = idp
        flip1 = 1.0d0
        flip2 = 1.0d0
        IF (jd == 0) THEN ! South pole
          jd = 1
          id1  = MODULO(id1  + hnx - 1, nx) + 1
          id1p = MODULO(id1p + hnx - 1, nx) + 1
          flip1 = -1.0d0
        ELSEIF(jd == ny) THEN ! North pole
          jdp = ny
          id2  = MODULO(id2  + hnx - 1, nx) + 1
          id2p = MODULO(id2p + hnx - 1, nx) + 1
          flip2 = -1.0d0
        ENDIF
	 
        ! Linearly interpolate  velocity to estimated departure point
	b1flip = b1*flip1
	b2flip = b2*flip2
        ud = c1*((a1*u0byr(kd,id1,jd)            &
                + a2*u0byr(kd,id1p,jd))*b1flip   &
	       + (a1*u0byr(kd,id2,jdp)           &
	        + a2*u0byr(kd,id2p,jdp))*b2flip) &
	   + c2*((a1*u0byr(kdp,id1,jd)           &
                + a2*u0byr(kdp,id1p,jd))*b1flip  &
	       + (a1*u0byr(kdp,id2,jdp)          &
	        + a2*u0byr(kdp,id2p,jdp))*b2flip)
        vd = c1*((a1*vatu0(kd,id1,jd)            &
                + a2*vatu0(kd,id1p,jd))*b1flip   &
	       + (a1*vatu0(kd,id2,jdp)           &
	        + a2*vatu0(kd,id2p,jdp))*b2flip) &
	   + c2*((a1*vatu0(kdp,id1,jd)           &
                + a2*vatu0(kdp,id1p,jd))*b1flip  &
	       + (a1*vatu0(kdp,id2,jdp)          &
	        + a2*vatu0(kdp,id2p,jdp))*b2flip)
        ed = c1w*((a1*eatu0(kd,id1,jd)        &
                 + a2*eatu0(kd,id1p,jd))*b1   &
	        + (a1*eatu0(kd,id2,jdp)       &
	         + a2*eatu0(kd,id2p,jdp))*b2) &
	   + c2w*((a1*eatu0(kdp,id1,jd)       &
                 + a2*eatu0(kdp,id1p,jd))*b1  &
	        + (a1*eatu0(kdp,id2,jdp)      &
	         + a2*eatu0(kdp,id2p,jdp))*b2)
        
        ! Rotate to arrival point Cartesian system
        sind = SIN(ydepu(k,i,j))
        cosd = COS(ydepu(k,i,j))
        sasd = sina*sind
        cacd = cosa*cosd
        dlambda = xu(i) - xdepu(k,i,j)
        sdl = SIN(dlambda)
        cdl = COS(dlambda)
        cosalphap1 = 1.0d0 + sasd + cacd*cdl
        m11 = (cacd + (1.0d0 + sasd)*cdl)/cosalphap1
        m12 = (sina + sind)*sdl/cosalphap1
        m21 = -m12
        m22 = m11
        urot = ud*m11 + vd*m12
        vrot = ud*m21 + vd*m22

        ! Hence calculate better estimate of departure point
        ! in arrival point Cartesian system
        x =    - dt*(alpha_x*ubyr(k,i,j) + beta_x*urot)*0.5d0*cosalphap1
        y =    - dt*(alpha_x*vatu(k,i,j) + beta_x*vrot)*0.5d0*cosalphap1

        ! Project back to spherical coordinate system
        z = SQRT(1.0d0 - (x*x + y*y))
        sind = y*cosa + z*sina
        ydepu(k,i,j) = ASIN(sind)
        dlambda = ATAN2(x,z*cosa - y*sina)
        xdepu(k,i,j) = MODULO(xu(i) + dlambda, twopi)

        ! Calculation of etadep
        edepu(k,i,j) = etap(k) - dt*(alpha_x*eatu(k,i,j) + beta_x*ed)


      ENDDO

    ENDDO
  ENDDO

ENDDO


END SUBROUTINE departureu2

! ========================================================

SUBROUTINE departurev2(uatv0,v0,eatv0,uatv,v,eatv)

! Compute departure points for v points
! Assume a first guess is already provided

USE grid
USE constants
USE timestep
USE departure

IMPLICIT NONE

INTEGER, PARAMETER :: ndepit = 2
REAL*8, INTENT(INOUT) :: uatv0(nz,nx,nyp), v0(nz,nx,nyp), eatv0(nz,nx,nyp), &
                         uatv(nz,nx,nyp),  v(nz,nx,nyp),  eatv(nz,nx,nyp)
REAL*8 :: v0byr(nz,nx,nyp), vbyr(nz,nx,nyp), rtemp(nz,nx)
INTEGER :: i, j, jm, k, id, idp, jd, jdp, kd, kdp, &
           hnx, idepit
REAL*8 :: sina, cosa, a1, a2, b1, b2, c1, c2, c1w, c2w, &
          rid, rjd, ud, vd, ed, x, y, z, &
	  sind, cosd, sasd, cacd, sdl, cdl, dlambda, &
	  urot, vrot, m11, m12, m21, m22, cosalphap1

! --------------------------------------------------------

! First divide horizontal velocities by r
DO j = 2, ny
  jm = j - 1
  rtemp = 0.5d0*(rp(:,:,jm) + rp(:,:,j))
  uatv0(:,:,j) = uatv0(:,:,j)/rtemp
  uatv(:,:,j) = uatv(:,:,j)/rtemp
  v0byr(:,:,j) = v0(:,:,j)/rtemp
  vbyr(:,:,j) = v(:,:,j)/rtemp
ENDDO
DO k = 1, nz
  rtemp(k,:) = SUM(rp(k,:,1))/nx
ENDDO
uatv0(:,:,1) = uatv0(:,:,1)/rtemp
uatv(:,:,1) = uatv(:,:,1)/rtemp
v0byr(:,:,1) = v0(:,:,1)/rtemp
vbyr(:,:,1) = v(:,:,1)/rtemp
DO k = 1, nz
  rtemp(k,:) = SUM(rp(k,:,ny))/nx
ENDDO
uatv0(:,:,nyp) = uatv0(:,:,nyp)/rtemp
uatv(:,:,nyp) = uatv(:,:,nyp)/rtemp
v0byr(:,:,nyp) = v0(:,:,nyp)/rtemp
vbyr(:,:,nyp) = v(:,:,nyp)/rtemp

hnx = nx/2

! Loop over iterations
DO idepit = 1, ndepit
 
  DO  j = 2, ny

    jm = j - 1

    ! Trig factors at arrival point
    sina = sinv(j)
    cosa = cosv(j)

    DO i = 1, nx

      DO k = 1, nz

        ! Determine departure cell indices based on latest estimate
        ! xdepv, ydepv, edepv
        rid = (xdepv(k,i,j) - xv(1))/dx
        id = FLOOR(rid)
        a2 = rid - id
        a1 = 1.0 - a2
        id = MODULO(id, nx) + 1
	idp = id + 1
	IF (id == nx) idp = 1
        rjd = (ydepv(k,i,j) - yv(1))/dy
        jd = FLOOR(rjd)
        b2 = rjd - jd
        b1 = 1.0 - b2
        jd = jd + 1
        jdp = jd + 1
	CALL findetacell(etap,edepv(k,i,j),nz,kd)
        IF (kd == 0) THEN
          ! Constant extrapolation below lowest u level
	  ! etadot linearly to zero
          kd = 1
          kdp = 2
          c1 = 1.0d0
          c2 = 0.0d0
	  c1w = (edepv(k,i,j) - etaw(1))/(etap(1) - etaw(1))
	  c2w = 0.0d0
        ELSEIF (kd == nz) THEN
          ! Constant extrapolation above highest u level
	  ! etadot linearly to zero
          kd = nz - 1
          kdp = nz
          c1 = 0.0d0
          c2 = 1.0d0
	  c1w = 0.0d0
	  c2w = (edepv(k,i,j) - etaw(nzp))/(etap(nz) - etaw(nzp))
        ELSE
          kdp = kd + 1
          c1 = (etap(kdp) - edepv(k,i,j))/(etap(kdp) - etap(kd))
          c2 = 1.0d0 - c1
	  c1w = c1
	  c2w = c2
	ENDIF
	
        ! Linearly interpolate  velocity to estimated departure point
        ud = c1*((a1*uatv0(kd,id,jd)        &
                + a2*uatv0(kd,idp,jd))*b1   &
	       + (a1*uatv0(kd,id,jdp)       &
	        + a2*uatv0(kd,idp,jdp))*b2) &
	   + c2*((a1*uatv0(kdp,id,jd)       &
                + a2*uatv0(kdp,idp,jd))*b1  &
	       + (a1*uatv0(kdp,id,jdp)      &
	        + a2*uatv0(kdp,idp,jdp))*b2)
        vd = c1*((a1*v0byr(kd,id,jd)        &
                + a2*v0byr(kd,idp,jd))*b1   &
	       + (a1*v0byr(kd,id,jdp)       &
	        + a2*v0byr(kd,idp,jdp))*b2) &
	   + c2*((a1*v0byr(kdp,id,jd)       &
                + a2*v0byr(kdp,idp,jd))*b1  &
	       + (a1*v0byr(kdp,id,jdp)      &
	        + a2*v0byr(kdp,idp,jdp))*b2)
        ed = c1w*((a1*eatv0(kd,id,jd)        &
                 + a2*eatv0(kd,idp,jd))*b1   &
	        + (a1*eatv0(kd,id,jdp)       &
	         + a2*eatv0(kd,idp,jdp))*b2) &
	   + c2w*((a1*eatv0(kdp,id,jd)       &
                 + a2*eatv0(kdp,idp,jd))*b1  &
	        + (a1*eatv0(kdp,id,jdp)      &
	         + a2*eatv0(kdp,idp,jdp))*b2)
		 
        ! Rotate to arrival point Cartesian system
        sind = SIN(ydepv(k,i,j))
        cosd = COS(ydepv(k,i,j))
        sasd = sina*sind
        cacd = cosa*cosd
        dlambda = xv(i) - xdepv(k,i,j)
        sdl = SIN(dlambda)
        cdl = COS(dlambda)
        cosalphap1 = 1.0d0 + sasd + cacd*cdl
        m11 = (cacd + (1.0d0 + sasd)*cdl)/cosalphap1
        m12 = (sina + sind)*sdl/cosalphap1
        m21 = -m12
        m22 = m11
        urot = ud*m11 + vd*m12
        vrot = ud*m21 + vd*m22

        ! Hence calculate better estimate of departure point
        ! in arrival point Cartesian system
        x =    - dt*(alpha_x*uatv(k,i,j) + beta_x*urot)*0.5d0*cosalphap1
        y =    - dt*(alpha_x*vbyr(k,i,j) + beta_x*vrot)*0.5d0*cosalphap1

        ! Project back to spherical coordinate system
        z = SQRT(1.0d0 - (x*x + y*y))
        sind = y*cosa + z*sina
        ydepv(k,i,j) = ASIN(sind)
        dlambda = ATAN2(x,z*cosa - y*sina)
        xdepv(k,i,j) = MODULO(xv(i) + dlambda, twopi)

        ! Calculation of etadep
        edepv(k,i,j) = etap(k) - dt*(alpha_x*eatv(k,i,j) + beta_x*ed)

      ENDDO

    ENDDO
  ENDDO

ENDDO


END SUBROUTINE departurev2

! ========================================================

SUBROUTINE departurew2(uatw0,vatw0,etadot0,uatw,vatw,etadot)

! Compute departure points for w points
! Assume a first guess is already provided

USE grid
USE constants
USE timestep
USE departure

IMPLICIT NONE

INTEGER, PARAMETER :: ndepit = 2
REAL*8, INTENT(INOUT) :: uatw0(nzp,nx,ny), vatw0(nzp,nx,ny), etadot0(nzp,nx,ny), &
                         uatw(nzp,nx,ny),  vatw(nzp,nx,ny),  etadot(nzp,nx,ny)
INTEGER :: i, j, k, id, idp, jd, jdp, kd, kdp, &
           id1, id1p, id2, id2p, hnx, idepit
REAL*8 :: sina, cosa, a1, a2, b1, b2, c1, c2, flip1, flip2, &
          b1flip, b2flip, rid, rjd, ud, vd, ed, x, y, z, &
	  sind, cosd, sasd, cacd, sdl, cdl, dlambda, &
	  urot, vrot, wrot, m11, m12, m21, m22, cosalphap1

! --------------------------------------------------------

! First divide horizontal velocities by r
uatw0 = uatw0/rw
uatw = uatw/rw
vatw0 = vatw0/rw
vatw = vatw/rw


hnx = nx/2

! Loop over iterations
DO idepit = 1, ndepit

  DO  j = 1, ny

    ! Trig factors at arrival point
    sina = sinp(j)
    cosa = cosp(j)

    DO i = 1, nx

      DO k = 1, nzp

        ! Determine departure cell indices based on latest estimate
        ! xdepw, ydepw, edepw
        rid = (xdepw(k,i,j) - xp(1))/dx
        id = FLOOR(rid)
        a2 = rid - id
        a1 = 1.0 - a2
        id = MODULO(id, nx) + 1
	idp = id + 1
	IF (id == nx) idp = 1
        rjd = (ydepw(k,i,j) - yp(1))/dy
        jd = FLOOR(rjd)
        b2 = rjd - jd
        b1 = 1.0 - b2
        jd = jd + 1
        jdp = jd + 1
	IF (k == 1) THEN
	  ! Bottom boundary - departure point is on boundary
	  kd = 1
	  kdp = 2
	  c1 = 1.0d0
	  c2 = 0.0d0
	ELSEIF (k == nzp) THEN
	  ! Top boundary - departure point is on boundary
	  kd = nz
	  kdp = nzp
	  c1 = 0.0d0
	  c2 = 1.0d0
	ELSE
          CALL findetacell(etaw,edepw(k,i,j),nzp,kd)
          kd = MIN(MAX(kd,1),nz)
          kdp = kd + 1
          c1 = (etaw(kdp) - edepw(k,i,j))/(etaw(kdp) - etaw(kd))
          c2 = 1.0d0 - c1
        ENDIF
	
        ! Tricks to handle polar case
        ! Here we interpolate across the pole.
        ! (An alternative would be to use the polar values of u and v
        ! and interpolate between nearest u latitude and the pole.)
        id1 = id
        id1p = idp
        id2 = id
        id2p = idp
        flip1 = 1.0d0
        flip2 = 1.0d0
        IF (jd == 0) THEN ! South pole
          jd = 1
          id1  = MODULO(id1  + hnx - 1, nx) + 1
          id1p = MODULO(id1p + hnx - 1, nx) + 1
          flip1 = -1.0d0
        ELSEIF(jd == ny) THEN ! North pole
          jdp = ny
          id2  = MODULO(id2  + hnx - 1, nx) + 1
          id2p = MODULO(id2p + hnx - 1, nx) + 1
          flip2 = -1.0d0
        ENDIF

        ! Linearly interpolate  velocity to estimated departure point
	b1flip = b1*flip1
	b2flip = b2*flip2
        ud = c1*((a1*uatw0(kd,id1,jd)            &
                + a2*uatw0(kd,id1p,jd))*b1flip   &
	       + (a1*uatw0(kd,id2,jdp)           &
	        + a2*uatw0(kd,id2p,jdp))*b2flip) &
	   + c2*((a1*uatw0(kdp,id1,jd)           &
                + a2*uatw0(kdp,id1p,jd))*b1flip  &
	       + (a1*uatw0(kdp,id2,jdp)          &
	        + a2*uatw0(kdp,id2p,jdp))*b2flip)
        vd = c1*((a1*vatw0(kd,id1,jd)            &
                + a2*vatw0(kd,id1p,jd))*b1flip   &
	       + (a1*vatw0(kd,id2,jdp)           &
	        + a2*vatw0(kd,id2p,jdp))*b2flip) &
	   + c2*((a1*vatw0(kdp,id1,jd)           &
                + a2*vatw0(kdp,id1p,jd))*b1flip  &
	       + (a1*vatw0(kdp,id2,jdp)          &
	        + a2*vatw0(kdp,id2p,jdp))*b2flip)
        ed = c1*((a1*etadot0(kd,id1,jd)        &
                + a2*etadot0(kd,id1p,jd))*b1   &
	       + (a1*etadot0(kd,id2,jdp)       &
	        + a2*etadot0(kd,id2p,jdp))*b2) &
	   + c2*((a1*etadot0(kdp,id1,jd)       &
                + a2*etadot0(kdp,id1p,jd))*b1  &
	       + (a1*etadot0(kdp,id2,jdp)      &
	        + a2*etadot0(kdp,id2p,jdp))*b2)
        
        ! Rotate to arrival point Cartesian system
        sind = SIN(ydepw(k,i,j))
        cosd = COS(ydepw(k,i,j))
        sasd = sina*sind
        cacd = cosa*cosd
        dlambda = xp(i) - xdepw(k,i,j)
        sdl = SIN(dlambda)
        cdl = COS(dlambda)
        cosalphap1 = 1.0d0 + sasd + cacd*cdl
        m11 = (cacd + (1.0d0 + sasd)*cdl)/cosalphap1
        m12 = (sina + sind)*sdl/cosalphap1
        m21 = -m12
        m22 = m11
        urot = ud*m11 + vd*m12
        vrot = ud*m21 + vd*m22

        ! Hence calculate better estimate of departure point
        ! in arrival point Cartesian system
        x =    - dt*(alpha_x*uatw(k,i,j) + beta_x*urot)*0.5d0*cosalphap1
        y =    - dt*(alpha_x*vatw(k,i,j) + beta_x*vrot)*0.5d0*cosalphap1

        ! Project back to spherical coordinate system
        z = SQRT(1.0d0 - (x*x + y*y))
        sind = y*cosa + z*sina
        ydepw(k,i,j) = ASIN(sind)
        dlambda = ATAN2(x,z*cosa - y*sina)
        xdepw(k,i,j) = MODULO(xp(i) + dlambda, twopi)

        ! Calculation of etadep
        edepw(k,i,j) = etaw(k) - dt*(alpha_x*etadot(k,i,j) + beta_x*ed)

      ENDDO

    ENDDO
  ENDDO



ENDDO


END SUBROUTINE departurew2

! ========================================================

SUBROUTINE departurep2(uatp0,vatp0,eatp0,uatp,vatp,eatp)

! Compute departure points for p points
! Assume a first guess is already provided

USE grid
USE constants
USE timestep
USE departure

IMPLICIT NONE

INTEGER, PARAMETER :: ndepit = 2
REAL*8, INTENT(INOUT) :: uatp0(nz,nx,ny), vatp0(nz,nx,ny), eatp0(nz,nx,ny), &
                         uatp(nz,nx,ny),  vatp(nz,nx,ny),  eatp(nz,nx,ny)
INTEGER :: i, j, k, id, idp, jd, jdp, kd, kdp, &
           id1, id1p, id2, id2p, hnx, idepit
REAL*8 :: sina, cosa, a1, a2, b1, b2, c1, c2, c1w, c2w, flip1, flip2, &
          b1flip, b2flip, rid, rjd, ud, vd, ed, x, y, z, rd, ra, &
	  sind, cosd, sasd, cacd, sdl, cdl, dlambda, &
	  urot, vrot, wrot, m11, m12, m21, m22, cosalphap1

! --------------------------------------------------------

! First divide horizontal velocities by r
uatp0 = uatp0/rp
uatp = uatp/rp
vatp0 = vatp0/rp
vatp = vatp/rp


hnx = nx/2

! Loop over iterations
DO idepit = 1, ndepit

  DO  j = 1, ny

    ! Trig factors at arrival point
    sina = sinp(j)
    cosa = cosp(j)

    DO i = 1, nx
      
      DO k = 1, nz

        ! Determine departure cell indices based on latest estimate
        ! xdepp, ydepp, edepp
        rid = (xdepp(k,i,j) - xp(1))/dx
        id = FLOOR(rid)
        a2 = rid - id
        a1 = 1.0 - a2
        id = MODULO(id, nx) + 1
	idp = id + 1
	IF (id == nx) idp = 1
        rjd = (ydepp(k,i,j) - yp(1))/dy
        jd = FLOOR(rjd)
        b2 = rjd - jd
        b1 = 1.0 - b2
        jd = jd + 1
        jdp = jd + 1
	CALL findetacell(etap,edepp(k,i,j),nz,kd)
        IF (kd == 0) THEN
          ! Constant extrapolation below lowest u level
	  ! etadot linearly to zero
          kd = 1
          kdp = 2
          c1 = 1.0d0
          c2 = 0.0d0
	  c1w = (edepp(k,i,j) - etaw(1))/(etap(1) - etaw(1))
	  c2w = 0.0d0
        ELSEIF (kd == nz) THEN
          ! Constant extrapolation above highest u level
	  ! etadot linearly to zero
          kd = nz - 1
          kdp = nz
          c1 = 0.0d0
          c2 = 1.0d0
	  c1w = 0.0d0
	  c2w = (edepp(k,i,j) - etaw(nzp))/(etap(nz) - etaw(nzp))
        ELSE
          kdp = kd + 1
          c1 = (etap(kdp) - edepp(k,i,j))/(etap(kdp) - etap(kd))
          c2 = 1.0d0 - c1
	  c1w = c1
	  c2w = c2
	ENDIF
	
        ! Tricks to handle polar case
        ! Here we interpolate across the pole.
        ! (An alternative would be to use the polar values of u and v
        ! and interpolate between nearest u latitude and the pole.)
        id1 = id
        id1p = idp
        id2 = id
        id2p = idp
        flip1 = 1.0d0
        flip2 = 1.0d0
        IF (jd == 0) THEN ! South pole
          jd = 1
          id1  = MODULO(id1  + hnx - 1, nx) + 1
          id1p = MODULO(id1p + hnx - 1, nx) + 1
          flip1 = -1.0d0
        ELSEIF(jd == ny) THEN ! North pole
          jdp = ny
          id2  = MODULO(id2  + hnx - 1, nx) + 1
          id2p = MODULO(id2p + hnx - 1, nx) + 1
          flip2 = -1.0d0
        ENDIF

        ! Linearly interpolate  velocity to estimated departure point
	b1flip = b1*flip1
	b2flip = b2*flip2
        ud = c1*((a1*uatp0(kd,id1,jd)            &
                + a2*uatp0(kd,id1p,jd))*b1flip   &
	       + (a1*uatp0(kd,id2,jdp)           &
	        + a2*uatp0(kd,id2p,jdp))*b2flip) &
	   + c2*((a1*uatp0(kdp,id1,jd)           &
                + a2*uatp0(kdp,id1p,jd))*b1flip  &
	       + (a1*uatp0(kdp,id2,jdp)          &
	        + a2*uatp0(kdp,id2p,jdp))*b2flip)
        vd = c1*((a1*vatp0(kd,id1,jd)            &
                + a2*vatp0(kd,id1p,jd))*b1flip   &
	       + (a1*vatp0(kd,id2,jdp)           &
	        + a2*vatp0(kd,id2p,jdp))*b2flip) &
	   + c2*((a1*vatp0(kdp,id1,jd)           &
                + a2*vatp0(kdp,id1p,jd))*b1flip  &
	       + (a1*vatp0(kdp,id2,jdp)          &
	        + a2*vatp0(kdp,id2p,jdp))*b2flip)
        ed = c1w*((a1*eatp0(kd,id1,jd)        &
                 + a2*eatp0(kd,id1p,jd))*b1   &
	        + (a1*eatp0(kd,id2,jdp)       &
	         + a2*eatp0(kd,id2p,jdp))*b2) &
	   + c2w*((a1*eatp0(kdp,id1,jd)       &
                 + a2*eatp0(kdp,id1p,jd))*b1  &
	        + (a1*eatp0(kdp,id2,jdp)      &
	         + a2*eatp0(kdp,id2p,jdp))*b2)
        
        ! Rotate to arrival point Cartesian system
        sind = SIN(ydepp(k,i,j))
        cosd = COS(ydepp(k,i,j))
        sasd = sina*sind
        cacd = cosa*cosd
        dlambda = xp(i) - xdepp(k,i,j)
        sdl = SIN(dlambda)
        cdl = COS(dlambda)
        cosalphap1 = 1.0d0 + sasd + cacd*cdl
        m11 = (cacd + (1.0d0 + sasd)*cdl)/cosalphap1
        m12 = (sina + sind)*sdl/cosalphap1
        m21 = -m12
        m22 = m11
        urot = ud*m11 + vd*m12
        vrot = ud*m21 + vd*m22

        ! Hence calculate better estimate of departure point
        ! in arrival point Cartesian system
        x =    - dt*(alpha_x*uatp(k,i,j) + beta_x*urot)*0.5d0*cosalphap1
        y =    - dt*(alpha_x*vatp(k,i,j) + beta_x*vrot)*0.5d0*cosalphap1

        ! Project back to spherical coordinate system
        z = SQRT(1.0d0 - (x*x + y*y))
        sind = y*cosa + z*sina
        ydepp(k,i,j) = ASIN(sind)
        dlambda = ATAN2(x,z*cosa - y*sina)
        xdepp(k,i,j) = MODULO(xp(i) + dlambda, twopi)

        ! Calculation of etadep
        edepp(k,i,j) = etap(k) - dt*(alpha_x*eatp(k,i,j) + beta_x*ed)


      ENDDO

    ENDDO
  ENDDO

ENDDO


END SUBROUTINE departurep2

! ========================================================

SUBROUTINE steptrajectory(xt,yt,et,uatp0,vatp0,eatp0,uatp,vatp,eatp)

! Step forward a Lagrangian trajectory
! The algorithm is similar to that for computing semi-Lagrangian departure
! points except that here the departure point is known and the arrival point
! is to be found.

USE grid
USE constants
USE timestep

IMPLICIT NONE

INTEGER, PARAMETER :: niter = 3
REAL*8, INTENT(INOUT) :: xt, yt, et
REAL*8, INTENT(IN) :: uatp0(nz,nx,ny), vatp0(nz,nx,ny), eatp0(nz,nx,ny), &
                      uatp(nz,nx,ny),  vatp(nz,nx,ny),  eatp(nz,nx,ny)
INTEGER :: id, idp, jd, jdp, kd, kdp, &
           id1, id1p, id2, id2p, hnx, iter
REAL*8 :: sina, cosa, a1, a2, b1, b2, c1, c2, c1w, c2w, flip1, flip2, &
          b1flip, b2flip, rid, rjd, ud, vd, ed, ua, va, ea, &
          x, y, z, rd, ra, xtarr, ytarr, etarr, &
	  sind, cosd, sasd, cacd, sdl, cdl, dlambda, &
	  urot, vrot, wrot, m11, m12, m21, m22, cosalphap1

! --------------------------------------------------------

! Note horizontal velocities were divided by r in departurep2
! for the tracer update. Therefore this does not need to be
! done here.


hnx = nx/2

! First, since departure point is not generally on the grid,
! we need to interpolate the velocity to the departure point.

! Determine departure cell indices based on xt, yt, et
rid = (xt - xp(1))/dx
id = FLOOR(rid)
a2 = rid - id
a1 = 1.0 - a2
id = MODULO(id, nx) + 1
idp = id + 1
IF (id == nx) idp = 1
rjd = (yt - yp(1))/dy
jd = FLOOR(rjd)
b2 = rjd - jd
b1 = 1.0 - b2
jd = jd + 1
jdp = jd + 1
CALL findetacell(etap,et,nz,kd)
IF (kd == 0) THEN
  ! Constant extrapolation below lowest u level
  ! etadot linearly to zero
  kd = 1
  kdp = 2
  c1 = 1.0d0
  c2 = 0.0d0
  c1w = (et - etaw(1))/(etap(1) - etaw(1))
  c2w = 0.0d0
ELSEIF (kd == nz) THEN
  ! Constant extrapolation above highest u level
  ! etadot linearly to zero
  kd = nz - 1
  kdp = nz
  c1 = 0.0d0
  c2 = 1.0d0
  c1w = 0.0d0
  c2w = (et - etaw(nzp))/(etap(nz) - etaw(nzp))
ELSE
  kdp = kd + 1
  c1 = (etap(kdp) - et)/(etap(kdp) - etap(kd))
  c2 = 1.0d0 - c1
  c1w = c1
  c2w = c2
ENDIF
	
! Tricks to handle polar case
! Here we interpolate across the pole.
id1 = id
id1p = idp
id2 = id
id2p = idp
flip1 = 1.0d0
flip2 = 1.0d0
IF (jd == 0) THEN ! South pole
  jd = 1
  id1  = MODULO(id1  + hnx - 1, nx) + 1
  id1p = MODULO(id1p + hnx - 1, nx) + 1
  flip1 = -1.0d0
ELSEIF(jd == ny) THEN ! North pole
  jdp = ny
  id2  = MODULO(id2  + hnx - 1, nx) + 1
  id2p = MODULO(id2p + hnx - 1, nx) + 1
  flip2 = -1.0d0
ENDIF

! Linearly interpolate  velocity to departure point
b1flip = b1*flip1
b2flip = b2*flip2
ud = c1*((a1*uatp0(kd,id1,jd)            &
        + a2*uatp0(kd,id1p,jd))*b1flip   &
       + (a1*uatp0(kd,id2,jdp)           &
	+ a2*uatp0(kd,id2p,jdp))*b2flip) &
   + c2*((a1*uatp0(kdp,id1,jd)           &
        + a2*uatp0(kdp,id1p,jd))*b1flip  &
       + (a1*uatp0(kdp,id2,jdp)          &
        + a2*uatp0(kdp,id2p,jdp))*b2flip)
vd = c1*((a1*vatp0(kd,id1,jd)            &
        + a2*vatp0(kd,id1p,jd))*b1flip   &
       + (a1*vatp0(kd,id2,jdp)           &
        + a2*vatp0(kd,id2p,jdp))*b2flip) &
   + c2*((a1*vatp0(kdp,id1,jd)           &
        + a2*vatp0(kdp,id1p,jd))*b1flip  &
       + (a1*vatp0(kdp,id2,jdp)          &
        + a2*vatp0(kdp,id2p,jdp))*b2flip)
ed = c1w*((a1*eatp0(kd,id1,jd)        &
         + a2*eatp0(kd,id1p,jd))*b1   &
	+ (a1*eatp0(kd,id2,jdp)       &
         + a2*eatp0(kd,id2p,jdp))*b2) &
   + c2w*((a1*eatp0(kdp,id1,jd)       &
         + a2*eatp0(kdp,id1p,jd))*b1  &
        + (a1*eatp0(kdp,id2,jdp)      &
         + a2*eatp0(kdp,id2p,jdp))*b2)


! Set first guess for the arrival point to be the departure point
xtarr = xt
ytarr = yt
etarr = et

! Trig factors at departure point
sind = SIN(yt)
cosd = COS(yt)

! Loop over iterations
DO iter = 1, niter

  ! Determine arrival cell indices based on latest estimate
  ! xtarr, ytarr, etarr
  rid = (xtarr - xp(1))/dx
  id = FLOOR(rid)
  a2 = rid - id
  a1 = 1.0 - a2
  id = MODULO(id, nx) + 1
  idp = id + 1
  IF (id == nx) idp = 1
  rjd = (ytarr - yp(1))/dy
  jd = FLOOR(rjd)
  b2 = rjd - jd
  b1 = 1.0 - b2
  jd = jd + 1
  jdp = jd + 1
  CALL findetacell(etap,etarr,nz,kd)
  IF (kd == 0) THEN
    ! Constant extrapolation below lowest u level
    ! etadot linearly to zero
    kd = 1
    kdp = 2
    c1 = 1.0d0
    c2 = 0.0d0
    c1w = (etarr - etaw(1))/(etap(1) - etaw(1))
    c2w = 0.0d0
  ELSEIF (kd == nz) THEN
    ! Constant extrapolation above highest u level
    ! etadot linearly to zero
    kd = nz - 1
    kdp = nz
    c1 = 0.0d0
    c2 = 1.0d0
    c1w = 0.0d0
    c2w = (etarr - etaw(nzp))/(etap(nz) - etaw(nzp))
  ELSE
    kdp = kd + 1
    c1 = (etap(kdp) - etarr)/(etap(kdp) - etap(kd))
    c2 = 1.0d0 - c1
    c1w = c1
    c2w = c2
  ENDIF
	
  ! Tricks to handle polar case
  ! Here we interpolate across the pole.
  id1 = id
  id1p = idp
  id2 = id
  id2p = idp
  flip1 = 1.0d0
  flip2 = 1.0d0
  IF (jd == 0) THEN ! South pole
    jd = 1
    id1  = MODULO(id1  + hnx - 1, nx) + 1
    id1p = MODULO(id1p + hnx - 1, nx) + 1
    flip1 = -1.0d0
  ELSEIF(jd == ny) THEN ! North pole
    jdp = ny
    id2  = MODULO(id2  + hnx - 1, nx) + 1
    id2p = MODULO(id2p + hnx - 1, nx) + 1
    flip2 = -1.0d0
  ENDIF

  ! Linearly interpolate velocity to estimated arrival point
  b1flip = b1*flip1
  b2flip = b2*flip2
  ua = c1*((a1*uatp(kd,id1,jd)            &
          + a2*uatp(kd,id1p,jd))*b1flip   &
         + (a1*uatp(kd,id2,jdp)           &
          + a2*uatp(kd,id2p,jdp))*b2flip) &
     + c2*((a1*uatp(kdp,id1,jd)           &
          + a2*uatp(kdp,id1p,jd))*b1flip  &
         + (a1*uatp(kdp,id2,jdp)          &
          + a2*uatp(kdp,id2p,jdp))*b2flip)
  va = c1*((a1*vatp(kd,id1,jd)            &
          + a2*vatp(kd,id1p,jd))*b1flip   &
         + (a1*vatp(kd,id2,jdp)           &
          + a2*vatp(kd,id2p,jdp))*b2flip) &
     + c2*((a1*vatp(kdp,id1,jd)           &
          + a2*vatp(kdp,id1p,jd))*b1flip  &
         + (a1*vatp(kdp,id2,jdp)          &
          + a2*vatp(kdp,id2p,jdp))*b2flip)
  ea = c1w*((a1*eatp(kd,id1,jd)        &
           + a2*eatp(kd,id1p,jd))*b1   &
          + (a1*eatp(kd,id2,jdp)       &
           + a2*eatp(kd,id2p,jdp))*b2) &
     + c2w*((a1*eatp(kdp,id1,jd)       &
           + a2*eatp(kdp,id1p,jd))*b1  &
          + (a1*eatp(kdp,id2,jdp)      &
           + a2*eatp(kdp,id2p,jdp))*b2)
        
  ! Rotate to departure point Cartesian system
  sina = SIN(ytarr)
  cosa = COS(ytarr)
  sasd = sina*sind
  cacd = cosa*cosd
  dlambda = xt - xtarr
  sdl = SIN(dlambda)
  cdl = COS(dlambda)
  cosalphap1 = 1.0d0 + sasd + cacd*cdl
  m11 = (cacd + (1.0d0 + sasd)*cdl)/cosalphap1
  m12 = (sina + sind)*sdl/cosalphap1
  m21 = -m12
  m22 = m11
  urot = ua*m11 + va*m12
  vrot = ua*m21 + va*m22

  ! Hence calculate better estimate of arrival point
  ! in departure point Cartesian system
  x = dt*(alpha_x*urot + beta_x*ud)*0.5d0*cosalphap1
  y = dt*(alpha_x*vrot + beta_x*vd)*0.5d0*cosalphap1

  ! Project back to spherical coordinate system
  z = SQRT(1.0d0 - (x*x + y*y))
  sina = y*cosd + z*sind
  ytarr = ASIN(sina)
  dlambda = ATAN2(x,z*cosd - y*sind)
  xtarr = MODULO(xt + dlambda, twopi)

  ! Calculation of etarr
  etarr = et + dt*(alpha_x*ea + beta_x*ed)

ENDDO

xt = xtarr
yt = ytarr
et = etarr


END SUBROUTINE steptrajectory

! ========================================================

SUBROUTINE interptraj(xt,yt,et,qp,qt)

! Interpolate a scalar field, stored at p points, to the trajectory
! parcel location (xt,yt,et)

USE grid

IMPLICIT NONE

REAL*8, INTENT(IN) :: xt, yt, et
REAL*8, INTENT(IN) :: qp(nz,nx,ny)
REAL*8, INTENT(OUT) :: qt

INTEGER :: id, idp, jd, jdp, kd, kdp, &
           id1, id1p, id2, id2p, hnx
REAL*8 :: a1, a2, b1, b2, c1, c2, rid, rjd

! --------------------------------------------------------

hnx = nx/2

! Determine cell indices based on xt, yt, et
rid = (xt - xp(1))/dx
id = FLOOR(rid)
a2 = rid - id
a1 = 1.0 - a2
id = MODULO(id, nx) + 1
idp = id + 1
IF (id == nx) idp = 1
rjd = (yt - yp(1))/dy
jd = FLOOR(rjd)
b2 = rjd - jd
b1 = 1.0 - b2
jd = jd + 1
jdp = jd + 1
CALL findetacell(etap,et,nz,kd)
IF (kd == 0) THEN
  ! Constant extrapolation below lowest p level
  kd = 1
  kdp = 2
  c1 = 1.0d0
  c2 = 0.0d0
ELSEIF (kd == nz) THEN
  ! Constant extrapolation above highest p level
  kd = nz - 1
  kdp = nz
  c1 = 0.0d0
  c2 = 1.0d0
ELSE
  kdp = kd + 1
  c1 = (etap(kdp) - et)/(etap(kdp) - etap(kd))
  c2 = 1.0d0 - c1
ENDIF
	
! Tricks to handle polar case
! Here we interpolate across the pole.
id1 = id
id1p = idp
id2 = id
id2p = idp
IF (jd == 0) THEN ! South pole
  jd = 1
  id1  = MODULO(id1  + hnx - 1, nx) + 1
  id1p = MODULO(id1p + hnx - 1, nx) + 1
ELSEIF(jd == ny) THEN ! North pole
  jdp = ny
  id2  = MODULO(id2  + hnx - 1, nx) + 1
  id2p = MODULO(id2p + hnx - 1, nx) + 1
ENDIF

! Linearly interpolate qp to trajectory point
qt = c1*((a1*qp(kd,id1,jd)            &
        + a2*qp(kd,id1p,jd))*b1       &
       + (a1*qp(kd,id2,jdp)           &
	+ a2*qp(kd,id2p,jdp))*b2)     &
   + c2*((a1*qp(kdp,id1,jd)           &
        + a2*qp(kdp,id1p,jd))*b1      &
       + (a1*qp(kdp,id2,jdp)          &
        + a2*qp(kdp,id2p,jdp))*b2)


END SUBROUTINE interptraj

! ========================================================

SUBROUTINE labeltraj(ntraj,xtraj,ytraj,etraj,qp,qtraj)

! Compute a set of trajectory labels qtraj by interpolating the field
! qp from p points to the trajectory locations xtraj, ytraj, etraj

USE grid, ONLY : nz, nx, ny

IMPLICIT NONE

INTEGER, INTENT(IN) :: ntraj
REAL*8, INTENT(IN) :: xtraj(ntraj), ytraj(ntraj), etraj(ntraj), &
                      qp(nz,nx,ny)
REAL*8, INTENT(OUT) :: qtraj(ntraj)
INTEGER :: i


DO i = 1, ntraj
  CALL interptraj(xtraj(i),ytraj(i),etraj(i),qp,qtraj(i))
ENDDO


END SUBROUTINE labeltraj

! ========================================================

SUBROUTINE lagrangep(xd,yd,ed,q,qnew)

! 3D Cubic Lagrange interpolation for the p-points on the
! Spherical C-grid
!
! It is assumed, when interpolating across the pole, that q is a scalar
! such as density

USE grid
USE util

IMPLICIT NONE

REAL*8, INTENT(IN) :: xd(nz,nx,ny), yd(nz,nx,ny), ed(nz,nx,ny), &
		      q(nz,nx,ny)
REAL*8, INTENT(OUT) :: qnew(nz,nx,ny)
INTEGER :: i, j, id, idm, idp, idpp, jd, jdm, jdp, jdpp, &
           i1, i1m, i1p, i1pp, hnx, k, kd, kdm, kdp, kdpp
REAL*8 :: xdd, ydd, edd, d12, d13, d14, d23, d24, d34, &
          denx1, denx2, denx3, denx4, deny1, deny2, deny3, deny4, &
	  dene1, dene2, dene3, dene4, &
	  dd1, dd2, dd3, dd4, &
	  facx1, facx2, facx3, facx4, facy1, facy2, facy3, facy4, &
	  face1, face2, face3, face4, &
	  q1, q2, q3, q4, qqm, qq, qqp, qqpp, &
	  x(nx), y(ny), e(nz)


! Regular gridded values are at p points
x = xp
y = yp
e = etap

! Handy value for handling poles
hnx = nx/2

! Some horizontal interpolation factors can be pre-computed on uniform grid

! For longitude direction
d12 = -dx
d13 = -2*dx
d14 = -3*dx
d23 = -dx
d24 = -2*dx
d34 = -dx
denx1 =  d12*d13*d14
denx2 = -d12*d23*d24
denx3 =  d13*d23*d34
denx4 = -d14*d24*d34

! For latitude direction
d12 = -dy
d13 = -2*dy
d14 = -3*dy
d23 = -dy
d24 = -2*dy
d34 = -dy
deny1 =  d12*d13*d14
deny2 = -d12*d23*d24
deny3 =  d13*d23*d34
deny4 = -d14*d24*d34


DO j = 1, ny
  DO i = 1, nx
    DO k = 1, nz

      ! Find indices of departure point and stencil
      xdd = xd(k,i,j)
      ydd = yd(k,i,j)
      edd = ed(k,i,j)

      id =   MODULO(FLOOR((xdd-x(1))/dx),nx) + 1
      idm =  MODULO(id-2,nx) + 1
      idp =  MODULO(id,nx) + 1
      idpp = MODULO(idp,nx) + 1
      jd =   FLOOR((ydd-y(1))/dy) + 1
      jdm =  jd - 1
      jdp =  jd + 1
      jdpp = jdp + 1
      CALL findetacell(e,edd,nz,kd)
      kdm = kd - 1
      kdp = kd + 1
      kdpp = kdp + 1

      ! Factors for x-interpolation
      dd1 = near(xdd - x(idm), twopi)
      dd2 = near(xdd - x(id), twopi)
      dd3 = near(xdd - x(idp), twopi)
      dd4 = near(xdd - x(idpp), twopi)
      facx1 = dd2*dd3*dd4/denx1
      facx2 = dd1*dd3*dd4/denx2
      facx3 = dd1*dd2*dd4/denx3
      facx4 = dd1*dd2*dd3/denx4
    
      ! Factors for y-interpolation
      dd2 = ydd - (y(1) + jdm*dy)
      dd1 = dd2 + dy
      dd3 = dd2 - dy
      dd4 = dd2 - 2*dy
      facy1 = dd2*dd3*dd4/deny1
      facy2 = dd1*dd3*dd4/deny2
      facy3 = dd1*dd2*dd4/deny3
      facy4 = dd1*dd2*dd3/deny4
      
      ! Factors for eta-interpolation
!      IF (kd == 0) THEN
!        ! Below lowest eta point linear extrapolation
!        kdm = 1
!        kd = 1
!        kdp = 1
!        kdpp = 2
!        face1 = 0.0d0
!        face2 = 0.0d0
!        face3 = 1.0d0
!        face4 = 0.0d0
!      ELSEIF (kd == nz) THEN
!        ! Above highest eta point linear extrapolation
!        kdm = nz - 1
!        kd = nz
!        kdp = nz
!        kdpp = nz
!        face1 = 0.0d0
!        face2 = 1.0d0
!        face3 = 0.0d0
!        face4 = 0.0d0
      ! In or below lowest eta interval linear interpolation or extrapolation
      IF (kd .LE. 1) THEN
        kdm = 1
        kd = 1
        kdp = 2
        kdpp = 3
        dd3 = edd - e(kdp)
        d23 = e(kd) - e(kdp)
        face1 = 0.0d0
        face2 = dd3/d23
        face3 = 1.0d0 - face2
        face4 = 0.0d0
      ! In or above highest eta interval linear interpolation or extrapolation
      ELSEIF (kd .GE. nz-1) THEN
        kdm = nz-2
        kd = nz-1
        kdp = nz
        kdpp = nz
        dd3 = edd - e(kdp)
        d23 = e(kd) - e(kdp)
        face1 = 0.0d0
        face2 = dd3/d23
        face3 = 1.0d0 - face2
        face4 = 0.0d0
      ELSE
        ! Otherwise cubic interpolation
        dd1 = edd - e(kdm)
        dd2 = edd - e(kd)
        dd3 = edd - e(kdp)
        dd4 = edd - e(kdpp)
        d12 = e(kdm) - e(kd)
        d13 = e(kdm) - e(kdp)
        d14 = e(kdm) - e(kdpp)
        d23 = e(kd) - e(kdp)
        d24 = e(kd) - e(kdpp)
        d34 = e(kdp) - e(kdpp)
        dene1 = d12*d13*d14
        dene2 = -d12*d23*d24
        dene3 = d13*d23*d34
        dene4 = -d14*d24*d34
        face1 = dd2*dd3*dd4/dene1
        face2 = dd1*dd3*dd4/dene2
        face3 = dd1*dd2*dd4/dene3
        face4 = dd1*dd2*dd3/dene4
      ENDIF

      ! Interpolate at four rows
      ! First
      i1m = idm
      i1 = id
      i1p = idp
      i1pp = idpp
      IF (jd .le. 1) THEN
        jdm = 1 - jdm
        i1m  = MODULO(i1m - hnx - 1, nx) + 1
        i1   = MODULO(i1 - hnx - 1, nx) + 1
        i1p  = MODULO(i1p - hnx - 1, nx) + 1
        i1pp = MODULO(i1pp - hnx - 1, nx) + 1
      ENDIF
      q1 = face1*q(kdm,i1m ,jdm) + face2*q(kd,i1m ,jdm) + face3*q(kdp,i1m ,jdm) + face4*q(kdpp,i1m ,jdm)
      q2 = face1*q(kdm,i1  ,jdm) + face2*q(kd,i1  ,jdm) + face3*q(kdp,i1  ,jdm) + face4*q(kdpp,i1  ,jdm)
      q3 = face1*q(kdm,i1p ,jdm) + face2*q(kd,i1p ,jdm) + face3*q(kdp,i1p ,jdm) + face4*q(kdpp,i1p ,jdm)
      q4 = face1*q(kdm,i1pp,jdm) + face2*q(kd,i1pp,jdm) + face3*q(kdp,i1pp,jdm) + face4*q(kdpp,i1pp,jdm)
      qqm = q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4

      ! Second
      i1m = idm
      i1 = id
      i1p = idp
      i1pp = idpp
      IF (jd .eq. 0) THEN
        jd = 1
        i1m  = MODULO(i1m - hnx - 1, nx) + 1
        i1   = MODULO(i1 - hnx - 1, nx) + 1
        i1p  = MODULO(i1p - hnx - 1, nx) + 1
        i1pp = MODULO(i1pp - hnx - 1, nx) + 1
      ENDIF
      q1 = face1*q(kdm,i1m ,jd) + face2*q(kd,i1m ,jd) + face3*q(kdp,i1m ,jd) + face4*q(kdpp,i1m ,jd)
      q2 = face1*q(kdm,i1  ,jd) + face2*q(kd,i1  ,jd) + face3*q(kdp,i1  ,jd) + face4*q(kdpp,i1  ,jd)
      q3 = face1*q(kdm,i1p ,jd) + face2*q(kd,i1p ,jd) + face3*q(kdp,i1p ,jd) + face4*q(kdpp,i1p ,jd)
      q4 = face1*q(kdm,i1pp,jd) + face2*q(kd,i1pp,jd) + face3*q(kdp,i1pp,jd) + face4*q(kdpp,i1pp,jd)
      qq = q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4

      ! Third
      i1m = idm
      i1 = id
      i1p = idp
      i1pp = idpp
      IF (jd .eq. ny) THEN
        jdp = ny
        i1m  = MODULO(i1m - hnx - 1, nx) + 1
        i1   = MODULO(i1 - hnx - 1, nx) + 1
        i1p  = MODULO(i1p - hnx - 1, nx) + 1
        i1pp = MODULO(i1pp - hnx - 1, nx) + 1
      ENDIF
      q1 = face1*q(kdm,i1m ,jdp) + face2*q(kd,i1m ,jdp) + face3*q(kdp,i1m ,jdp) + face4*q(kdpp,i1m ,jdp)
      q2 = face1*q(kdm,i1  ,jdp) + face2*q(kd,i1  ,jdp) + face3*q(kdp,i1  ,jdp) + face4*q(kdpp,i1  ,jdp)
      q3 = face1*q(kdm,i1p ,jdp) + face2*q(kd,i1p ,jdp) + face3*q(kdp,i1p ,jdp) + face4*q(kdpp,i1p ,jdp)
      q4 = face1*q(kdm,i1pp,jdp) + face2*q(kd,i1pp,jdp) + face3*q(kdp,i1pp,jdp) + face4*q(kdpp,i1pp,jdp)
      qqp = q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4

      ! Fourth
      i1m = idm
      i1 = id
      i1p = idp
      i1pp = idpp
      IF (jd .ge. ny - 1) THEN
        jdpp = 2*ny + 1 - jdpp
        i1m  = MODULO(i1m - hnx - 1, nx) + 1
        i1   = MODULO(i1 - hnx - 1, nx) + 1
        i1p  = MODULO(i1p - hnx - 1, nx) + 1
        i1pp = MODULO(i1pp - hnx - 1, nx) + 1
      ENDIF
      q1 = face1*q(kdm,i1m ,jdpp) + face2*q(kd,i1m ,jdpp) + face3*q(kdp,i1m ,jdpp) + face4*q(kdpp,i1m ,jdpp)
      q2 = face1*q(kdm,i1  ,jdpp) + face2*q(kd,i1  ,jdpp) + face3*q(kdp,i1  ,jdpp) + face4*q(kdpp,i1  ,jdpp)
      q3 = face1*q(kdm,i1p ,jdpp) + face2*q(kd,i1p ,jdpp) + face3*q(kdp,i1p ,jdpp) + face4*q(kdpp,i1p ,jdpp)
      q4 = face1*q(kdm,i1pp,jdpp) + face2*q(kd,i1pp,jdpp) + face3*q(kdp,i1pp,jdpp) + face4*q(kdpp,i1pp,jdpp)
      qqpp = q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4

      ! Interpolate in y
      qnew(k,i,j) = qqm*facy1 + qq*facy2 + qqp*facy3 + qqpp*facy4

    ENDDO
  ENDDO
ENDDO


END SUBROUTINE lagrangep

! ========================================================

SUBROUTINE lagrangeu(xd,yd,ed,qu,qv,qw,qunew,qvnew,qwnew)

! 3D Cubic Lagrange interpolation for the three components of a vector field
! stored at the u-points on the spherical C-grid. Note the `flip' factors
! used when interpolating across the poles


USE grid
USE util

IMPLICIT NONE

REAL*8, INTENT(IN) :: xd(nz,nx,ny), yd(nz,nx,ny), ed(nz,nx,ny), &
		      qu(nz,nx,ny), qv(nz,nx,ny), qw(nz,nx,ny)
REAL*8, INTENT(OUT) :: qunew(nz,nx,ny), qvnew(nz,nx,ny), qwnew(nz,nx,ny)
INTEGER :: i, j, id, idm, idp, idpp, jd, jdm, jdp, jdpp, &
           i1, i1m, i1p, i1pp, hnx, k, kd, kdm, kdp, kdpp
REAL*8 :: xdd, ydd, edd, d12, d13, d14, d23, d24, d34, &
          denx1, denx2, denx3, denx4, deny1, deny2, deny3, deny4, &
	  dene1, dene2, dene3, dene4, &
	  dd1, dd2, dd3, dd4, &
	  facx1, facx2, facx3, facx4, facy1, facy2, facy3, facy4, &
	  face1, face2, face3, face4, &
	  q1, q2, q3, q4, qqum, qqu, qqup, qqupp, &
          qqvm, qqv, qqvp, qqvpp, qqwm, qqw, qqwp, qqwpp, &
	  flip, x(nx), y(ny), e(nz)

! Regular gridded values are at u points
x = xu
y = yu
e = etap

! Handy value for handling poles
hnx = nx/2

! Some horizontal interpolation factors can be pre-computed on uniform grid

! For longitude direction
d12 = -dx
d13 = -2*dx
d14 = -3*dx
d23 = -dx
d24 = -2*dx
d34 = -dx
denx1 =  d12*d13*d14
denx2 = -d12*d23*d24
denx3 =  d13*d23*d34
denx4 = -d14*d24*d34

! For latitude direction
d12 = -dy
d13 = -2*dy
d14 = -3*dy
d23 = -dy
d24 = -2*dy
d34 = -dy
deny1 =  d12*d13*d14
deny2 = -d12*d23*d24
deny3 =  d13*d23*d34
deny4 = -d14*d24*d34


DO j = 1, ny
  DO i = 1, nx
    DO k = 1, nz

      ! Find indices of departure point and stencil
      xdd = xd(k,i,j)
      ydd = yd(k,i,j)
      edd = ed(k,i,j)

      id =   MODULO(FLOOR((xdd-x(1))/dx),nx) + 1
      idm =  MODULO(id-2,nx) + 1
      idp =  MODULO(id,nx) + 1
      idpp = MODULO(idp,nx) + 1
      jd =   FLOOR((ydd-y(1))/dy) + 1
      jdm =  jd - 1
      jdp =  jd + 1
      jdpp = jdp + 1
      CALL findetacell(e,edd,nz,kd)
      kdm = kd - 1
      kdp = kd + 1
      kdpp = kdp + 1

      ! Factors for x-interpolation
      dd1 = near(xdd - x(idm), twopi)
      dd2 = near(xdd - x(id), twopi)
      dd3 = near(xdd - x(idp), twopi)
      dd4 = near(xdd - x(idpp), twopi)
      facx1 = dd2*dd3*dd4/denx1
      facx2 = dd1*dd3*dd4/denx2
      facx3 = dd1*dd2*dd4/denx3
      facx4 = dd1*dd2*dd3/denx4

      ! Factors for y-interpolation
      dd2 = ydd - (y(1) + jdm*dy)
      dd1 = dd2 + dy
      dd3 = dd2 - dy
      dd4 = dd2 - 2*dy
      facy1 = dd2*dd3*dd4/deny1
      facy2 = dd1*dd3*dd4/deny2
      facy3 = dd1*dd2*dd4/deny3
      facy4 = dd1*dd2*dd3/deny4

      ! Factors for eta-interpolation
!      IF (kd == 0) THEN
!        ! Below lowest eta point constant extrapolation
!        kdm = 1
!        kd = 1
!        kdp = 1
!        kdpp = 2
!        face1 = 0.0d0
!        face2 = 0.0d0
!        face3 = 1.0d0
!        face4 = 0.0d0
!      ELSEIF (kd == nz) THEN
!        ! Above highest eta point constant extrapolation
!        kdm = nz - 1
!        kd = nz
!        kdp = nz
!        kdpp = nz
!        face1 = 0.0d0
!        face2 = 1.0d0
!        face3 = 0.0d0
!        face4 = 0.0d0
      ! In lowest eta interval linear interpolation
      IF (kd .LE. 1) THEN
        kdm = 1
        kd = 1
        kdp = 2
        kdpp = 3
        dd3 = edd - e(kdp)
        d23 = e(kd) - e(kdp)
        face1 = 0.0d0
        face2 = dd3/d23
        face3 = 1.0d0 - face2
        face4 = 0.0d0
      ! In highest eta interval linear interpolation
      ELSEIF (kd .GE. nz-1) THEN
        kdm = nz-2
        kd = nz-1
        kdp = nz
        kdpp = nz
        dd3 = edd - e(kdp)
        d23 = e(kd) - e(kdp)
        face1 = 0.0d0
        face2 = dd3/d23
        face3 = 1.0d0 - face2
        face4 = 0.0d0
      ELSE
        ! Otherwise cubic interpolation
        dd1 = edd - e(kdm)
        dd2 = edd - e(kd)
        dd3 = edd - e(kdp)
        dd4 = edd - e(kdpp)
        d12 = e(kdm) - e(kd)
        d13 = e(kdm) - e(kdp)
        d14 = e(kdm) - e(kdpp)
        d23 = e(kd) - e(kdp)
        d24 = e(kd) - e(kdpp)
        d34 = e(kdp) - e(kdpp)
        dene1 = d12*d13*d14
        dene2 = -d12*d23*d24
        dene3 = d13*d23*d34
        dene4 = -d14*d24*d34
        face1 = dd2*dd3*dd4/dene1
        face2 = dd1*dd3*dd4/dene2
        face3 = dd1*dd2*dd4/dene3
        face4 = dd1*dd2*dd3/dene4
      ENDIF

      ! Interpolate at four rows
      ! First
      i1m = idm
      i1 = id
      i1p = idp
      i1pp = idpp
      flip = 1.0d0
      IF (jd .le. 1) THEN
        jdm = 1 - jdm
        i1m  = MODULO(i1m - hnx - 1, nx) + 1
        i1   = MODULO(i1 - hnx - 1, nx) + 1
        i1p  = MODULO(i1p - hnx - 1, nx) + 1
        i1pp = MODULO(i1pp - hnx - 1, nx) + 1
        flip = -1.0d0
      ENDIF
      q1 = face1*qu(kdm,i1m ,jdm) + face2*qu(kd,i1m ,jdm) + face3*qu(kdp,i1m ,jdm) + face4*qu(kdpp,i1m ,jdm)
      q2 = face1*qu(kdm,i1  ,jdm) + face2*qu(kd,i1  ,jdm) + face3*qu(kdp,i1  ,jdm) + face4*qu(kdpp,i1  ,jdm)
      q3 = face1*qu(kdm,i1p ,jdm) + face2*qu(kd,i1p ,jdm) + face3*qu(kdp,i1p ,jdm) + face4*qu(kdpp,i1p ,jdm)
      q4 = face1*qu(kdm,i1pp,jdm) + face2*qu(kd,i1pp,jdm) + face3*qu(kdp,i1pp,jdm) + face4*qu(kdpp,i1pp,jdm)
      qqum = (q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4)*flip
      q1 = face1*qv(kdm,i1m ,jdm) + face2*qv(kd,i1m ,jdm) + face3*qv(kdp,i1m ,jdm) + face4*qv(kdpp,i1m ,jdm)
      q2 = face1*qv(kdm,i1  ,jdm) + face2*qv(kd,i1  ,jdm) + face3*qv(kdp,i1  ,jdm) + face4*qv(kdpp,i1  ,jdm)
      q3 = face1*qv(kdm,i1p ,jdm) + face2*qv(kd,i1p ,jdm) + face3*qv(kdp,i1p ,jdm) + face4*qv(kdpp,i1p ,jdm)
      q4 = face1*qv(kdm,i1pp,jdm) + face2*qv(kd,i1pp,jdm) + face3*qv(kdp,i1pp,jdm) + face4*qv(kdpp,i1pp,jdm)
      qqvm = (q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4)*flip
      q1 = face1*qw(kdm,i1m ,jdm) + face2*qw(kd,i1m ,jdm) + face3*qw(kdp,i1m ,jdm) + face4*qw(kdpp,i1m ,jdm)
      q2 = face1*qw(kdm,i1  ,jdm) + face2*qw(kd,i1  ,jdm) + face3*qw(kdp,i1  ,jdm) + face4*qw(kdpp,i1  ,jdm)
      q3 = face1*qw(kdm,i1p ,jdm) + face2*qw(kd,i1p ,jdm) + face3*qw(kdp,i1p ,jdm) + face4*qw(kdpp,i1p ,jdm)
      q4 = face1*qw(kdm,i1pp,jdm) + face2*qw(kd,i1pp,jdm) + face3*qw(kdp,i1pp,jdm) + face4*qw(kdpp,i1pp,jdm)
      qqwm = q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4

      ! Second
      i1m = idm
      i1 = id
      i1p = idp
      i1pp = idpp
      flip = 1.0d0
      IF (jd .eq. 0) THEN
        jd = 1
        i1m  = MODULO(i1m - hnx - 1, nx) + 1
        i1   = MODULO(i1 - hnx - 1, nx) + 1
        i1p  = MODULO(i1p - hnx - 1, nx) + 1
        i1pp = MODULO(i1pp - hnx - 1, nx) + 1
        flip = -1.0d0
      ENDIF
      q1 = face1*qu(kdm,i1m ,jd) + face2*qu(kd,i1m ,jd) + face3*qu(kdp,i1m ,jd) + face4*qu(kdpp,i1m ,jd)
      q2 = face1*qu(kdm,i1  ,jd) + face2*qu(kd,i1  ,jd) + face3*qu(kdp,i1  ,jd) + face4*qu(kdpp,i1  ,jd)
      q3 = face1*qu(kdm,i1p ,jd) + face2*qu(kd,i1p ,jd) + face3*qu(kdp,i1p ,jd) + face4*qu(kdpp,i1p ,jd)
      q4 = face1*qu(kdm,i1pp,jd) + face2*qu(kd,i1pp,jd) + face3*qu(kdp,i1pp,jd) + face4*qu(kdpp,i1pp,jd)
      qqu = (q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4)*flip
      q1 = face1*qv(kdm,i1m ,jd) + face2*qv(kd,i1m ,jd) + face3*qv(kdp,i1m ,jd) + face4*qv(kdpp,i1m ,jd)
      q2 = face1*qv(kdm,i1  ,jd) + face2*qv(kd,i1  ,jd) + face3*qv(kdp,i1  ,jd) + face4*qv(kdpp,i1  ,jd)
      q3 = face1*qv(kdm,i1p ,jd) + face2*qv(kd,i1p ,jd) + face3*qv(kdp,i1p ,jd) + face4*qv(kdpp,i1p ,jd)
      q4 = face1*qv(kdm,i1pp,jd) + face2*qv(kd,i1pp,jd) + face3*qv(kdp,i1pp,jd) + face4*qv(kdpp,i1pp,jd)
      qqv = (q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4)*flip
      q1 = face1*qw(kdm,i1m ,jd) + face2*qw(kd,i1m ,jd) + face3*qw(kdp,i1m ,jd) + face4*qw(kdpp,i1m ,jd)
      q2 = face1*qw(kdm,i1  ,jd) + face2*qw(kd,i1  ,jd) + face3*qw(kdp,i1  ,jd) + face4*qw(kdpp,i1  ,jd)
      q3 = face1*qw(kdm,i1p ,jd) + face2*qw(kd,i1p ,jd) + face3*qw(kdp,i1p ,jd) + face4*qw(kdpp,i1p ,jd)
      q4 = face1*qw(kdm,i1pp,jd) + face2*qw(kd,i1pp,jd) + face3*qw(kdp,i1pp,jd) + face4*qw(kdpp,i1pp,jd)
      qqw = q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4

      ! Third
      i1m = idm
      i1 = id
      i1p = idp
      i1pp = idpp
      flip = 1.0d0
      IF (jd .eq. ny) THEN
        jdp = ny
        i1m  = MODULO(i1m - hnx - 1, nx) + 1
        i1   = MODULO(i1 - hnx - 1, nx) + 1
        i1p  = MODULO(i1p - hnx - 1, nx) + 1
        i1pp = MODULO(i1pp - hnx - 1, nx) + 1
        flip = -1.0d0
      ENDIF
      q1 = face1*qu(kdm,i1m ,jdp) + face2*qu(kd,i1m ,jdp) + face3*qu(kdp,i1m ,jdp) + face4*qu(kdpp,i1m ,jdp)
      q2 = face1*qu(kdm,i1  ,jdp) + face2*qu(kd,i1  ,jdp) + face3*qu(kdp,i1  ,jdp) + face4*qu(kdpp,i1  ,jdp)
      q3 = face1*qu(kdm,i1p ,jdp) + face2*qu(kd,i1p ,jdp) + face3*qu(kdp,i1p ,jdp) + face4*qu(kdpp,i1p ,jdp)
      q4 = face1*qu(kdm,i1pp,jdp) + face2*qu(kd,i1pp,jdp) + face3*qu(kdp,i1pp,jdp) + face4*qu(kdpp,i1pp,jdp)
      qqup = (q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4)*flip
      q1 = face1*qv(kdm,i1m ,jdp) + face2*qv(kd,i1m ,jdp) + face3*qv(kdp,i1m ,jdp) + face4*qv(kdpp,i1m ,jdp)
      q2 = face1*qv(kdm,i1  ,jdp) + face2*qv(kd,i1  ,jdp) + face3*qv(kdp,i1  ,jdp) + face4*qv(kdpp,i1  ,jdp)
      q3 = face1*qv(kdm,i1p ,jdp) + face2*qv(kd,i1p ,jdp) + face3*qv(kdp,i1p ,jdp) + face4*qv(kdpp,i1p ,jdp)
      q4 = face1*qv(kdm,i1pp,jdp) + face2*qv(kd,i1pp,jdp) + face3*qv(kdp,i1pp,jdp) + face4*qv(kdpp,i1pp,jdp)
      qqvp = (q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4)*flip
      q1 = face1*qw(kdm,i1m ,jdp) + face2*qw(kd,i1m ,jdp) + face3*qw(kdp,i1m ,jdp) + face4*qw(kdpp,i1m ,jdp)
      q2 = face1*qw(kdm,i1  ,jdp) + face2*qw(kd,i1  ,jdp) + face3*qw(kdp,i1  ,jdp) + face4*qw(kdpp,i1  ,jdp)
      q3 = face1*qw(kdm,i1p ,jdp) + face2*qw(kd,i1p ,jdp) + face3*qw(kdp,i1p ,jdp) + face4*qw(kdpp,i1p ,jdp)
      q4 = face1*qw(kdm,i1pp,jdp) + face2*qw(kd,i1pp,jdp) + face3*qw(kdp,i1pp,jdp) + face4*qw(kdpp,i1pp,jdp)
      qqwp = q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4

      ! Fourth
      i1m = idm
      i1 = id
      i1p = idp
      i1pp = idpp
      flip = 1.0d0
      IF (jd .ge. ny - 1) THEN
        jdpp = 2*ny + 1 - jdpp
        i1m  = MODULO(i1m - hnx - 1, nx) + 1
        i1   = MODULO(i1 - hnx - 1, nx) + 1
        i1p  = MODULO(i1p - hnx - 1, nx) + 1
        i1pp = MODULO(i1pp - hnx - 1, nx) + 1
        flip = -1.0d0
      ENDIF
      q1 = face1*qu(kdm,i1m ,jdpp) + face2*qu(kd,i1m ,jdpp) + face3*qu(kdp,i1m ,jdpp) + face4*qu(kdpp,i1m ,jdpp)
      q2 = face1*qu(kdm,i1  ,jdpp) + face2*qu(kd,i1  ,jdpp) + face3*qu(kdp,i1  ,jdpp) + face4*qu(kdpp,i1  ,jdpp)
      q3 = face1*qu(kdm,i1p ,jdpp) + face2*qu(kd,i1p ,jdpp) + face3*qu(kdp,i1p ,jdpp) + face4*qu(kdpp,i1p ,jdpp)
      q4 = face1*qu(kdm,i1pp,jdpp) + face2*qu(kd,i1pp,jdpp) + face3*qu(kdp,i1pp,jdpp) + face4*qu(kdpp,i1pp,jdpp)
      qqupp = (q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4)*flip
      q1 = face1*qv(kdm,i1m ,jdpp) + face2*qv(kd,i1m ,jdpp) + face3*qv(kdp,i1m ,jdpp) + face4*qv(kdpp,i1m ,jdpp)
      q2 = face1*qv(kdm,i1  ,jdpp) + face2*qv(kd,i1  ,jdpp) + face3*qv(kdp,i1  ,jdpp) + face4*qv(kdpp,i1  ,jdpp)
      q3 = face1*qv(kdm,i1p ,jdpp) + face2*qv(kd,i1p ,jdpp) + face3*qv(kdp,i1p ,jdpp) + face4*qv(kdpp,i1p ,jdpp)
      q4 = face1*qv(kdm,i1pp,jdpp) + face2*qv(kd,i1pp,jdpp) + face3*qv(kdp,i1pp,jdpp) + face4*qv(kdpp,i1pp,jdpp)
      qqvpp = (q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4)*flip
      q1 = face1*qw(kdm,i1m ,jdpp) + face2*qw(kd,i1m ,jdpp) + face3*qw(kdp,i1m ,jdpp) + face4*qw(kdpp,i1m ,jdpp)
      q2 = face1*qw(kdm,i1  ,jdpp) + face2*qw(kd,i1  ,jdpp) + face3*qw(kdp,i1  ,jdpp) + face4*qw(kdpp,i1  ,jdpp)
      q3 = face1*qw(kdm,i1p ,jdpp) + face2*qw(kd,i1p ,jdpp) + face3*qw(kdp,i1p ,jdpp) + face4*qw(kdpp,i1p ,jdpp)
      q4 = face1*qw(kdm,i1pp,jdpp) + face2*qw(kd,i1pp,jdpp) + face3*qw(kdp,i1pp,jdpp) + face4*qw(kdpp,i1pp,jdpp)
      qqwpp = q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4

      ! Interpolate in y
      qunew(k,i,j) = qqum*facy1 + qqu*facy2 + qqup*facy3 + qqupp*facy4
      qvnew(k,i,j) = qqvm*facy1 + qqv*facy2 + qqvp*facy3 + qqvpp*facy4
      qwnew(k,i,j) = qqwm*facy1 + qqw*facy2 + qqwp*facy3 + qqwpp*facy4

    ENDDO
  ENDDO
ENDDO


END SUBROUTINE lagrangeu

! ========================================================

SUBROUTINE lagrangev(xd,yd,ed,qu,qv,qw,qunew,qvnew,qwnew)

! 3D Cubic Lagrange interpolation for the three components of a vector field
! stored at the v-points on the spherical C-grid. Polar values of the components
! must be computed before calling this routine. Note the `flip' factors used
! when interpolating across the poles. Results are calculated only for non-polar rows.

USE grid
USE util

IMPLICIT NONE

REAL*8, INTENT(IN) :: xd(nz,nx,nyp), yd(nz,nx,nyp), ed(nz,nx,nyp), &
		      qu(nz,nx,nyp), qv(nz,nx,nyp), qw(nz,nx,nyp)
REAL*8, INTENT(OUT) :: qunew(nz,nx,nyp), qvnew(nz,nx,nyp), qwnew(nz,nx,nyp)
INTEGER :: i, j, id, idm, idp, idpp, jd, jdm, jdp, jdpp, &
           i1, i1m, i1p, i1pp, hnx, k, kd, kdm, kdp, kdpp
REAL*8 :: xdd, ydd, edd, d12, d13, d14, d23, d24, d34, &
          denx1, denx2, denx3, denx4, deny1, deny2, deny3, deny4, &
	  dene1, dene2, dene3, dene4, &
	  dd1, dd2, dd3, dd4, &
	  facx1, facx2, facx3, facx4, facy1, facy2, facy3, facy4, &
	  face1, face2, face3, face4, &
	  q1, q2, q3, q4, qqum, qqu, qqup, qqupp, &
          qqvm, qqv, qqvp, qqvpp, qqwm, qqw, qqwp, qqwpp, &
	  flip, x(nx), y(nyp), e(nz)

! Regular gridded values are at v points
x = xv
y = yv
e = etap

! Handy value for handling poles
hnx = nx/2

! Some horizontal interpolation factors can be pre-computed on uniform grid

! For longitude direction
d12 = -dx
d13 = -2*dx
d14 = -3*dx
d23 = -dx
d24 = -2*dx
d34 = -dx
denx1 =  d12*d13*d14
denx2 = -d12*d23*d24
denx3 =  d13*d23*d34
denx4 = -d14*d24*d34

! For latitude direction
d12 = -dy
d13 = -2*dy
d14 = -3*dy
d23 = -dy
d24 = -2*dy
d34 = -dy
deny1 =  d12*d13*d14
deny2 = -d12*d23*d24
deny3 =  d13*d23*d34
deny4 = -d14*d24*d34


DO j = 2, ny
  DO i = 1, nx
    DO k = 1, nz

      ! Find indices of departure point and stencil
      xdd = xd(k,i,j)
      ydd = yd(k,i,j)
      edd = ed(k,i,j)

      id =   MODULO(FLOOR((xdd-x(1))/dx),nx) + 1
      idm =  MODULO(id-2,nx) + 1
      idp =  MODULO(id,nx) + 1
      idpp = MODULO(idp,nx) + 1
      jd =   FLOOR((ydd-y(1))/dy) + 1
      ! Trap cases that might occur because of roundoff
      IF (jd .lt. 1) jd = 1
      IF (jd .gt. ny) jd = ny
      jdm =  jd - 1
      jdp =  jd + 1
      jdpp = jdp + 1
      CALL findetacell(e,edd,nz,kd)
      kdm = kd - 1
      kdp = kd + 1
      kdpp = kdp + 1
      
      ! Factors for x-interpolation
      dd1 = near(xdd - x(idm), twopi)
      dd2 = near(xdd - x(id), twopi)
      dd3 = near(xdd - x(idp), twopi)
      dd4 = near(xdd - x(idpp), twopi)
      facx1 = dd2*dd3*dd4/denx1
      facx2 = dd1*dd3*dd4/denx2
      facx3 = dd1*dd2*dd4/denx3
      facx4 = dd1*dd2*dd3/denx4

      ! Factors for y-interpolation
      dd2 = ydd - (y(1) + jdm*dy)
      dd1 = dd2 + dy
      dd3 = dd2 - dy
      dd4 = dd2 - 2*dy
      facy1 = dd2*dd3*dd4/deny1
      facy2 = dd1*dd3*dd4/deny2
      facy3 = dd1*dd2*dd4/deny3
      facy4 = dd1*dd2*dd3/deny4

      ! Factors for eta-interpolation
!      IF (kd == 0) THEN
!        ! Below lowest eta point constant extrapolation
!        kdm = 1
!        kd = 1
!        kdp = 1
!        kdpp = 2
!        face1 = 0.0d0
!        face2 = 0.0d0
!        face3 = 1.0d0
!        face4 = 0.0d0
!      ELSEIF (kd == nz) THEN
!        ! Above highest eta point constant extrapolation
!        kdm = nz - 1
!        kd = nz
!        kdp = nz
!        kdpp = nz
!        face1 = 0.0d0
!        face2 = 1.0d0
!        face3 = 0.0d0
!        face4 = 0.0d0
!      ! In lowest eta interval linear interpolation
      IF (kd .LE. 1) THEN
        kdm = 1
        kd = 1
        kdp = 2
        kdpp = 3
        dd3 = edd - e(kdp)
        d23 = e(kd) - e(kdp)
        face1 = 0.0d0
        face2 = dd3/d23
        face3 = 1.0d0 - face2
        face4 = 0.0d0
      ! In highest eta interval linear interpolation
      ELSEIF (kd .GE. nz-1) THEN
        kdm = nz-2
        kd = nz-1
        kdp = nz
        kdpp = nz
        dd3 = edd - e(kdp)
        d23 = e(kd) - e(kdp)
        face1 = 0.0d0
        face2 = dd3/d23
        face3 = 1.0d0 - face2
        face4 = 0.0d0
      ELSE
        ! Otherwise cubic interpolation
        dd1 = edd - e(kdm)
        dd2 = edd - e(kd)
        dd3 = edd - e(kdp)
        dd4 = edd - e(kdpp)
        d12 = e(kdm) - e(kd)
        d13 = e(kdm) - e(kdp)
        d14 = e(kdm) - e(kdpp)
        d23 = e(kd) - e(kdp)
        d24 = e(kd) - e(kdpp)
        d34 = e(kdp) - e(kdpp)
        dene1 = d12*d13*d14
        dene2 = -d12*d23*d24
        dene3 = d13*d23*d34
        dene4 = -d14*d24*d34
        face1 = dd2*dd3*dd4/dene1
        face2 = dd1*dd3*dd4/dene2
        face3 = dd1*dd2*dd4/dene3
        face4 = dd1*dd2*dd3/dene4
      ENDIF

      ! Interpolate at four rows
      ! First
      i1m = idm
      i1 = id
      i1p = idp
      i1pp = idpp
      flip = 1.0d0
      IF (jd == 1) THEN
        jdm = 2
        i1m  = MODULO(i1m - hnx - 1, nx) + 1
        i1   = MODULO(i1 - hnx - 1, nx) + 1
        i1p  = MODULO(i1p - hnx - 1, nx) + 1
        i1pp = MODULO(i1pp - hnx - 1, nx) + 1
        flip = -1.0d0
      ENDIF
      q1 = face1*qu(kdm,i1m ,jdm) + face2*qu(kd,i1m ,jdm) + face3*qu(kdp,i1m ,jdm) + face4*qu(kdpp,i1m ,jdm)
      q2 = face1*qu(kdm,i1  ,jdm) + face2*qu(kd,i1  ,jdm) + face3*qu(kdp,i1  ,jdm) + face4*qu(kdpp,i1  ,jdm)
      q3 = face1*qu(kdm,i1p ,jdm) + face2*qu(kd,i1p ,jdm) + face3*qu(kdp,i1p ,jdm) + face4*qu(kdpp,i1p ,jdm)
      q4 = face1*qu(kdm,i1pp,jdm) + face2*qu(kd,i1pp,jdm) + face3*qu(kdp,i1pp,jdm) + face4*qu(kdpp,i1pp,jdm)
      qqum = (q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4)*flip
      q1 = face1*qv(kdm,i1m ,jdm) + face2*qv(kd,i1m ,jdm) + face3*qv(kdp,i1m ,jdm) + face4*qv(kdpp,i1m ,jdm)
      q2 = face1*qv(kdm,i1  ,jdm) + face2*qv(kd,i1  ,jdm) + face3*qv(kdp,i1  ,jdm) + face4*qv(kdpp,i1  ,jdm)
      q3 = face1*qv(kdm,i1p ,jdm) + face2*qv(kd,i1p ,jdm) + face3*qv(kdp,i1p ,jdm) + face4*qv(kdpp,i1p ,jdm)
      q4 = face1*qv(kdm,i1pp,jdm) + face2*qv(kd,i1pp,jdm) + face3*qv(kdp,i1pp,jdm) + face4*qv(kdpp,i1pp,jdm)
      qqvm = (q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4)*flip
      q1 = face1*qw(kdm,i1m ,jdm) + face2*qw(kd,i1m ,jdm) + face3*qw(kdp,i1m ,jdm) + face4*qw(kdpp,i1m ,jdm)
      q2 = face1*qw(kdm,i1  ,jdm) + face2*qw(kd,i1  ,jdm) + face3*qw(kdp,i1  ,jdm) + face4*qw(kdpp,i1  ,jdm)
      q3 = face1*qw(kdm,i1p ,jdm) + face2*qw(kd,i1p ,jdm) + face3*qw(kdp,i1p ,jdm) + face4*qw(kdpp,i1p ,jdm)
      q4 = face1*qw(kdm,i1pp,jdm) + face2*qw(kd,i1pp,jdm) + face3*qw(kdp,i1pp,jdm) + face4*qw(kdpp,i1pp,jdm)
      qqwm = q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4

      ! Second
      q1 = face1*qu(kdm,idm ,jd) + face2*qu(kd,idm ,jd) + face3*qu(kdp,idm ,jd) + face4*qu(kdpp,idm ,jd)
      q2 = face1*qu(kdm,id  ,jd) + face2*qu(kd,id  ,jd) + face3*qu(kdp,id  ,jd) + face4*qu(kdpp,id  ,jd)
      q3 = face1*qu(kdm,idp ,jd) + face2*qu(kd,idp ,jd) + face3*qu(kdp,idp ,jd) + face4*qu(kdpp,idp ,jd)
      q4 = face1*qu(kdm,idpp,jd) + face2*qu(kd,idpp,jd) + face3*qu(kdp,idpp,jd) + face4*qu(kdpp,idpp,jd)
      qqu = (q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4)
      q1 = face1*qv(kdm,idm ,jd) + face2*qv(kd,idm ,jd) + face3*qv(kdp,idm ,jd) + face4*qv(kdpp,idm ,jd)
      q2 = face1*qv(kdm,id  ,jd) + face2*qv(kd,id  ,jd) + face3*qv(kdp,id  ,jd) + face4*qv(kdpp,id  ,jd)
      q3 = face1*qv(kdm,idp ,jd) + face2*qv(kd,idp ,jd) + face3*qv(kdp,idp ,jd) + face4*qv(kdpp,idp ,jd)
      q4 = face1*qv(kdm,idpp,jd) + face2*qv(kd,idpp,jd) + face3*qv(kdp,idpp,jd) + face4*qv(kdpp,idpp,jd)
      qqv = (q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4)
      q1 = face1*qw(kdm,idm ,jd) + face2*qw(kd,idm ,jd) + face3*qw(kdp,idm ,jd) + face4*qw(kdpp,idm ,jd)
      q2 = face1*qw(kdm,id  ,jd) + face2*qw(kd,id  ,jd) + face3*qw(kdp,id  ,jd) + face4*qw(kdpp,id  ,jd)
      q3 = face1*qw(kdm,idp ,jd) + face2*qw(kd,idp ,jd) + face3*qw(kdp,idp ,jd) + face4*qw(kdpp,idp ,jd)
      q4 = face1*qw(kdm,idpp,jd) + face2*qw(kd,idpp,jd) + face3*qw(kdp,idpp,jd) + face4*qw(kdpp,idpp,jd)
      qqw = q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4

      ! Third
      q1 = face1*qu(kdm,idm ,jdp) + face2*qu(kd,idm ,jdp) + face3*qu(kdp,idm ,jdp) + face4*qu(kdpp,idm ,jdp)
      q2 = face1*qu(kdm,id  ,jdp) + face2*qu(kd,id  ,jdp) + face3*qu(kdp,id  ,jdp) + face4*qu(kdpp,id  ,jdp)
      q3 = face1*qu(kdm,idp ,jdp) + face2*qu(kd,idp ,jdp) + face3*qu(kdp,idp ,jdp) + face4*qu(kdpp,idp ,jdp)
      q4 = face1*qu(kdm,idpp,jdp) + face2*qu(kd,idpp,jdp) + face3*qu(kdp,idpp,jdp) + face4*qu(kdpp,idpp,jdp)
      qqup = (q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4)
      q1 = face1*qv(kdm,idm ,jdp) + face2*qv(kd,idm ,jdp) + face3*qv(kdp,idm ,jdp) + face4*qv(kdpp,idm ,jdp)
      q2 = face1*qv(kdm,id  ,jdp) + face2*qv(kd,id  ,jdp) + face3*qv(kdp,id  ,jdp) + face4*qv(kdpp,id  ,jdp)
      q3 = face1*qv(kdm,idp ,jdp) + face2*qv(kd,idp ,jdp) + face3*qv(kdp,idp ,jdp) + face4*qv(kdpp,idp ,jdp)
      q4 = face1*qv(kdm,idpp,jdp) + face2*qv(kd,idpp,jdp) + face3*qv(kdp,idpp,jdp) + face4*qv(kdpp,idpp,jdp)
      qqvp = (q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4)
      q1 = face1*qw(kdm,idm ,jdp) + face2*qw(kd,idm ,jdp) + face3*qw(kdp,idm ,jdp) + face4*qw(kdpp,idm ,jdp)
      q2 = face1*qw(kdm,id  ,jdp) + face2*qw(kd,id  ,jdp) + face3*qw(kdp,id  ,jdp) + face4*qw(kdpp,id  ,jdp)
      q3 = face1*qw(kdm,idp ,jdp) + face2*qw(kd,idp ,jdp) + face3*qw(kdp,idp ,jdp) + face4*qw(kdpp,idp ,jdp)
      q4 = face1*qw(kdm,idpp,jdp) + face2*qw(kd,idpp,jdp) + face3*qw(kdp,idpp,jdp) + face4*qw(kdpp,idpp,jdp)
      qqwp = q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4

      ! Fourth
      i1m = idm
      i1 = id
      i1p = idp
      i1pp = idpp
      flip = 1.0d0
      IF (jd == ny) THEN
        jdpp = ny
        i1m  = MODULO(i1m - hnx - 1, nx) + 1
        i1   = MODULO(i1 - hnx - 1, nx) + 1
        i1p  = MODULO(i1p - hnx - 1, nx) + 1
        i1pp = MODULO(i1pp - hnx - 1, nx) + 1
        flip = -1.0d0
      ENDIF
      q1 = face1*qu(kdm,i1m ,jdpp) + face2*qu(kd,i1m ,jdpp) + face3*qu(kdp,i1m ,jdpp) + face4*qu(kdpp,i1m ,jdpp)
      q2 = face1*qu(kdm,i1  ,jdpp) + face2*qu(kd,i1  ,jdpp) + face3*qu(kdp,i1  ,jdpp) + face4*qu(kdpp,i1  ,jdpp)
      q3 = face1*qu(kdm,i1p ,jdpp) + face2*qu(kd,i1p ,jdpp) + face3*qu(kdp,i1p ,jdpp) + face4*qu(kdpp,i1p ,jdpp)
      q4 = face1*qu(kdm,i1pp,jdpp) + face2*qu(kd,i1pp,jdpp) + face3*qu(kdp,i1pp,jdpp) + face4*qu(kdpp,i1pp,jdpp)
      qqupp = (q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4)*flip
      q1 = face1*qv(kdm,i1m ,jdpp) + face2*qv(kd,i1m ,jdpp) + face3*qv(kdp,i1m ,jdpp) + face4*qv(kdpp,i1m ,jdpp)
      q2 = face1*qv(kdm,i1  ,jdpp) + face2*qv(kd,i1  ,jdpp) + face3*qv(kdp,i1  ,jdpp) + face4*qv(kdpp,i1  ,jdpp)
      q3 = face1*qv(kdm,i1p ,jdpp) + face2*qv(kd,i1p ,jdpp) + face3*qv(kdp,i1p ,jdpp) + face4*qv(kdpp,i1p ,jdpp)
      q4 = face1*qv(kdm,i1pp,jdpp) + face2*qv(kd,i1pp,jdpp) + face3*qv(kdp,i1pp,jdpp) + face4*qv(kdpp,i1pp,jdpp)
      qqvpp = (q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4)*flip
      q1 = face1*qw(kdm,i1m ,jdpp) + face2*qw(kd,i1m ,jdpp) + face3*qw(kdp,i1m ,jdpp) + face4*qw(kdpp,i1m ,jdpp)
      q2 = face1*qw(kdm,i1  ,jdpp) + face2*qw(kd,i1  ,jdpp) + face3*qw(kdp,i1  ,jdpp) + face4*qw(kdpp,i1  ,jdpp)
      q3 = face1*qw(kdm,i1p ,jdpp) + face2*qw(kd,i1p ,jdpp) + face3*qw(kdp,i1p ,jdpp) + face4*qw(kdpp,i1p ,jdpp)
      q4 = face1*qw(kdm,i1pp,jdpp) + face2*qw(kd,i1pp,jdpp) + face3*qw(kdp,i1pp,jdpp) + face4*qw(kdpp,i1pp,jdpp)
      qqwpp = q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4

      ! Interpolate in y
      qunew(k,i,j) = qqum*facy1 + qqu*facy2 + qqup*facy3 + qqupp*facy4
      qvnew(k,i,j) = qqvm*facy1 + qqv*facy2 + qqvp*facy3 + qqvpp*facy4
      qwnew(k,i,j) = qqwm*facy1 + qqw*facy2 + qqwp*facy3 + qqwpp*facy4

    ENDDO
  ENDDO
ENDDO


END SUBROUTINE lagrangev

! ========================================================

SUBROUTINE lagrangew(xd,yd,ed,qu,qv,qw,qt,qunew,qvnew,qwnew,qtnew)

! 3D Cubic Lagrange interpolation for the three components of a vector field
! plus a scalar stored at the w-points on the spherical C-grid. Note the `flip' factors
! used when interpolating across the poles

USE switches
USE grid
USE util

IMPLICIT NONE

REAL*8, INTENT(IN) :: xd(nzp,nx,ny), yd(nzp,nx,ny), ed(nzp,nx,ny), &
		      qu(nzp,nx,ny), qv(nzp,nx,ny), qw(nzp,nx,ny), qt(nzp,nx,ny)
REAL*8, INTENT(OUT) :: qunew(nzp,nx,ny), qvnew(nzp,nx,ny), qwnew(nzp,nx,ny), qtnew(nzp,nx,ny)
INTEGER :: i, j, id, idm, idp, idpp, jd, jdm, jdp, jdpp, &
           i1, i1m, i1p, i1pp, hnx, k, kd, kdm, kdp, kdpp
REAL*8 :: xdd, ydd, edd, d12, d13, d14, d23, d24, d34, &
          denx1, denx2, denx3, denx4, deny1, deny2, deny3, deny4, &
	  dene1, dene2, dene3, dene4, &
	  dd1, dd2, dd3, dd4, &
	  facx1, facx2, facx3, facx4, facy1, facy2, facy3, facy4, &
	  face1, face2, face3, face4, &
	  q1, q2, q3, q4, qqum, qqu, qqup, qqupp, &
          qqvm, qqv, qqvp, qqvpp, qqwm, qqw, qqwp, qqwpp, &
          qqtm, qqt, qqtp, qqtpp, flip, x(nx), y(ny), e(nzp), &
	  qtmin, qtmax


! Regular gridded values are at w points
x = xp
y = yp
e = etaw

! Handy value for handling poles
hnx = nx/2

! Some horizontal interpolation factors can be pre-computed on uniform grid

! For longitude direction
d12 = -dx
d13 = -2*dx
d14 = -3*dx
d23 = -dx
d24 = -2*dx
d34 = -dx
denx1 =  d12*d13*d14
denx2 = -d12*d23*d24
denx3 =  d13*d23*d34
denx4 = -d14*d24*d34

! For latitude direction
d12 = -dy
d13 = -2*dy
d14 = -3*dy
d23 = -dy
d24 = -2*dy
d34 = -dy
deny1 =  d12*d13*d14
deny2 = -d12*d23*d24
deny3 =  d13*d23*d34
deny4 = -d14*d24*d34


DO j = 1, ny
  DO i = 1, nx
    DO k = 1, nzp

      ! Find indices of departure point and stencil 
      xdd = xd(k,i,j)
      ydd = yd(k,i,j)
      edd = ed(k,i,j)

      id =   MODULO(FLOOR((xdd-x(1))/dx),nx) + 1
      idm =  MODULO(id-2,nx) + 1
      idp =  MODULO(id,nx) + 1
      idpp = MODULO(idp,nx) + 1
      jd =   FLOOR((ydd-y(1))/dy) + 1
      jdm =  jd - 1
      jdp =  jd + 1
      jdpp = jdp + 1
      CALL findetacell(e,edd,nzp,kd)
      kd = MIN(MAX(kd,1),nz)
      kdm = kd - 1
      kdp = kd + 1
      kdpp = kdp + 1

      ! Factors for x-interpolation
      dd1 = near(xdd - x(idm), twopi)
      dd2 = near(xdd - x(id), twopi)
      dd3 = near(xdd - x(idp), twopi)
      dd4 = near(xdd - x(idpp), twopi)
      facx1 = dd2*dd3*dd4/denx1
      facx2 = dd1*dd3*dd4/denx2
      facx3 = dd1*dd2*dd4/denx3
      facx4 = dd1*dd2*dd3/denx4

      ! Factors for y-interpolation
      dd2 = ydd - (y(1) + jdm*dy)
      dd1 = dd2 + dy
      dd3 = dd2 - dy
      dd4 = dd2 - 2*dy
      facy1 = dd2*dd3*dd4/deny1
      facy2 = dd1*dd3*dd4/deny2
      facy3 = dd1*dd2*dd4/deny3
      facy4 = dd1*dd2*dd3/deny4

      ! Factors for eta-interpolation
      IF (k == 1) THEN
        ! On the bottom boundary; just use boundary values
        kdm = 1
        kd = 1
        kdp = 1
        kdpp = 2
        face1 = 0.0d0
        face2 = 0.0d0
        face3 = 1.0d0
        face4 = 0.0d0
      ELSEIF (k == nzp) THEN
        ! On the top boundary; just use boundary values
        kdm = nz
        kd = nzp
        kdp = nzp
        kdpp = nzp
        face1 = 0.0d0
        face2 = 1.0d0
        face3 = 0.0d0
        face4 = 0.0d0
      ! In lowest eta interval linear interpolation
      ELSEIF (kd == 1) THEN
        kdm = 1
        kd = 1
        kdp = 2
        kdpp = 3
        dd3 = edd - e(kdp)
        d23 = e(kd) - e(kdp)
        face1 = 0.0d0
        face2 = dd3/d23
        face3 = 1.0d0 - face2
        face4 = 0.0d0
      ! In highest eta interval linear interpolation
      ELSEIF (kd == nz) THEN
        kdm = nz-1
        kd = nz
        kdp = nzp
        kdpp = nzp
        dd3 = edd - e(kdp)
        d23 = e(kd) - e(kdp)
        face1 = 0.0d0
        face2 = dd3/d23
        face3 = 1.0d0 - face2
        face4 = 0.0d0
      ELSE
        ! Otherwise cubic interpolation
        dd1 = edd - e(kdm)
        dd2 = edd - e(kd)
        dd3 = edd - e(kdp)
        dd4 = edd - e(kdpp)
        d12 = e(kdm) - e(kd)
        d13 = e(kdm) - e(kdp)
        d14 = e(kdm) - e(kdpp)
        d23 = e(kd) - e(kdp)
        d24 = e(kd) - e(kdpp)
        d34 = e(kdp) - e(kdpp)
        dene1 = d12*d13*d14
        dene2 = -d12*d23*d24
        dene3 = d13*d23*d34
        dene4 = -d14*d24*d34
        face1 = dd2*dd3*dd4/dene1
        face2 = dd1*dd3*dd4/dene2
        face3 = dd1*dd2*dd4/dene3
        face4 = dd1*dd2*dd3/dene4
      ENDIF

      ! Interpolate at four rows
      ! First
      i1m = idm
      i1 = id
      i1p = idp
      i1pp = idpp
      flip = 1.0d0
      IF (jd .le. 1) THEN
        jdm = 1 - jdm
        i1m  = MODULO(i1m - hnx - 1, nx) + 1
        i1   = MODULO(i1 - hnx - 1, nx) + 1
        i1p  = MODULO(i1p - hnx - 1, nx) + 1
        i1pp = MODULO(i1pp - hnx - 1, nx) + 1
        flip = -1.0d0
      ENDIF
      q1 = face1*qu(kdm,i1m ,jdm) + face2*qu(kd,i1m ,jdm) + face3*qu(kdp,i1m ,jdm) + face4*qu(kdpp,i1m ,jdm)
      q2 = face1*qu(kdm,i1  ,jdm) + face2*qu(kd,i1  ,jdm) + face3*qu(kdp,i1  ,jdm) + face4*qu(kdpp,i1  ,jdm)
      q3 = face1*qu(kdm,i1p ,jdm) + face2*qu(kd,i1p ,jdm) + face3*qu(kdp,i1p ,jdm) + face4*qu(kdpp,i1p ,jdm)
      q4 = face1*qu(kdm,i1pp,jdm) + face2*qu(kd,i1pp,jdm) + face3*qu(kdp,i1pp,jdm) + face4*qu(kdpp,i1pp,jdm)
      qqum = (q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4)*flip
      q1 = face1*qv(kdm,i1m ,jdm) + face2*qv(kd,i1m ,jdm) + face3*qv(kdp,i1m ,jdm) + face4*qv(kdpp,i1m ,jdm)
      q2 = face1*qv(kdm,i1  ,jdm) + face2*qv(kd,i1  ,jdm) + face3*qv(kdp,i1  ,jdm) + face4*qv(kdpp,i1  ,jdm)
      q3 = face1*qv(kdm,i1p ,jdm) + face2*qv(kd,i1p ,jdm) + face3*qv(kdp,i1p ,jdm) + face4*qv(kdpp,i1p ,jdm)
      q4 = face1*qv(kdm,i1pp,jdm) + face2*qv(kd,i1pp,jdm) + face3*qv(kdp,i1pp,jdm) + face4*qv(kdpp,i1pp,jdm)
      qqvm = (q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4)*flip
      q1 = face1*qw(kdm,i1m ,jdm) + face2*qw(kd,i1m ,jdm) + face3*qw(kdp,i1m ,jdm) + face4*qw(kdpp,i1m ,jdm)
      q2 = face1*qw(kdm,i1  ,jdm) + face2*qw(kd,i1  ,jdm) + face3*qw(kdp,i1  ,jdm) + face4*qw(kdpp,i1  ,jdm)
      q3 = face1*qw(kdm,i1p ,jdm) + face2*qw(kd,i1p ,jdm) + face3*qw(kdp,i1p ,jdm) + face4*qw(kdpp,i1p ,jdm)
      q4 = face1*qw(kdm,i1pp,jdm) + face2*qw(kd,i1pp,jdm) + face3*qw(kdp,i1pp,jdm) + face4*qw(kdpp,i1pp,jdm)
      qqwm = q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4
      q1 = face1*qt(kdm,i1m ,jdm) + face2*qt(kd,i1m ,jdm) + face3*qt(kdp,i1m ,jdm) + face4*qt(kdpp,i1m ,jdm)
      q2 = face1*qt(kdm,i1  ,jdm) + face2*qt(kd,i1  ,jdm) + face3*qt(kdp,i1  ,jdm) + face4*qt(kdpp,i1  ,jdm)
      q3 = face1*qt(kdm,i1p ,jdm) + face2*qt(kd,i1p ,jdm) + face3*qt(kdp,i1p ,jdm) + face4*qt(kdpp,i1p ,jdm)
      q4 = face1*qt(kdm,i1pp,jdm) + face2*qt(kd,i1pp,jdm) + face3*qt(kdp,i1pp,jdm) + face4*qt(kdpp,i1pp,jdm)
      IF (monotone) THEN
        qtmin = MIN(qt(kd,i1m ,jdm),qt(kdp,i1m ,jdm))
	qtmax = MAX(qt(kd,i1m ,jdm),qt(kdp,i1m ,jdm))
	q1 = MIN(MAX(qtmin,q1),qtmax)
	qtmin = MIN(qt(kd,i1  ,jdm),qt(kdp,i1  ,jdm))
	qtmax = MAX(qt(kd,i1  ,jdm),qt(kdp,i1  ,jdm))
	q2 = MIN(MAX(qtmin,q2),qtmax)
	qtmin = MIN(qt(kd,i1p ,jdm),qt(kdp,i1p ,jdm))
	qtmax = MAX(qt(kd,i1p ,jdm),qt(kdp,i1p ,jdm))
	q3 = MIN(MAX(qtmin,q3),qtmax)
	qtmin = MIN(qt(kd,i1pp,jdm),qt(kdp,i1pp,jdm))
	qtmax = MAX(qt(kd,i1pp,jdm),qt(kdp,i1pp,jdm))
	q4 = MIN(MAX(qtmin,q4),qtmax)
      ENDIF
      qqtm = q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4
      IF (monotone) THEN
        qtmin = MIN(q2,q3)
	qtmax = MAX(q2,q3)
	qqtm = MIN(MAX(qtmin,qqtm),qtmax)
      ENDIF
      
      ! Second
      i1m = idm
      i1 = id
      i1p = idp
      i1pp = idpp
      flip = 1.0d0
      IF (jd .eq. 0) THEN
        jd = 1
        i1m  = MODULO(i1m - hnx - 1, nx) + 1
        i1   = MODULO(i1 - hnx - 1, nx) + 1
        i1p  = MODULO(i1p - hnx - 1, nx) + 1
        i1pp = MODULO(i1pp - hnx - 1, nx) + 1
        flip = -1.0d0
      ENDIF
      q1 = face1*qu(kdm,i1m ,jd) + face2*qu(kd,i1m ,jd) + face3*qu(kdp,i1m ,jd) + face4*qu(kdpp,i1m ,jd)
      q2 = face1*qu(kdm,i1  ,jd) + face2*qu(kd,i1  ,jd) + face3*qu(kdp,i1  ,jd) + face4*qu(kdpp,i1  ,jd)
      q3 = face1*qu(kdm,i1p ,jd) + face2*qu(kd,i1p ,jd) + face3*qu(kdp,i1p ,jd) + face4*qu(kdpp,i1p ,jd)
      q4 = face1*qu(kdm,i1pp,jd) + face2*qu(kd,i1pp,jd) + face3*qu(kdp,i1pp,jd) + face4*qu(kdpp,i1pp,jd)
      qqu = (q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4)*flip
      q1 = face1*qv(kdm,i1m ,jd) + face2*qv(kd,i1m ,jd) + face3*qv(kdp,i1m ,jd) + face4*qv(kdpp,i1m ,jd)
      q2 = face1*qv(kdm,i1  ,jd) + face2*qv(kd,i1  ,jd) + face3*qv(kdp,i1  ,jd) + face4*qv(kdpp,i1  ,jd)
      q3 = face1*qv(kdm,i1p ,jd) + face2*qv(kd,i1p ,jd) + face3*qv(kdp,i1p ,jd) + face4*qv(kdpp,i1p ,jd)
      q4 = face1*qv(kdm,i1pp,jd) + face2*qv(kd,i1pp,jd) + face3*qv(kdp,i1pp,jd) + face4*qv(kdpp,i1pp,jd)
      qqv = (q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4)*flip
      q1 = face1*qw(kdm,i1m ,jd) + face2*qw(kd,i1m ,jd) + face3*qw(kdp,i1m ,jd) + face4*qw(kdpp,i1m ,jd)
      q2 = face1*qw(kdm,i1  ,jd) + face2*qw(kd,i1  ,jd) + face3*qw(kdp,i1  ,jd) + face4*qw(kdpp,i1  ,jd)
      q3 = face1*qw(kdm,i1p ,jd) + face2*qw(kd,i1p ,jd) + face3*qw(kdp,i1p ,jd) + face4*qw(kdpp,i1p ,jd)
      q4 = face1*qw(kdm,i1pp,jd) + face2*qw(kd,i1pp,jd) + face3*qw(kdp,i1pp,jd) + face4*qw(kdpp,i1pp,jd)
      qqw = q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4
      q1 = face1*qt(kdm,i1m ,jd) + face2*qt(kd,i1m ,jd) + face3*qt(kdp,i1m ,jd) + face4*qt(kdpp,i1m ,jd)
      q2 = face1*qt(kdm,i1  ,jd) + face2*qt(kd,i1  ,jd) + face3*qt(kdp,i1  ,jd) + face4*qt(kdpp,i1  ,jd)
      q3 = face1*qt(kdm,i1p ,jd) + face2*qt(kd,i1p ,jd) + face3*qt(kdp,i1p ,jd) + face4*qt(kdpp,i1p ,jd)
      q4 = face1*qt(kdm,i1pp,jd) + face2*qt(kd,i1pp,jd) + face3*qt(kdp,i1pp,jd) + face4*qt(kdpp,i1pp,jd)
      IF (monotone) THEN
        qtmin = MIN(qt(kd,i1m ,jd),qt(kdp,i1m ,jd))
	qtmax = MAX(qt(kd,i1m ,jd),qt(kdp,i1m ,jd))
	q1 = MIN(MAX(qtmin,q1),qtmax)
	qtmin = MIN(qt(kd,i1  ,jd),qt(kdp,i1  ,jd))
	qtmax = MAX(qt(kd,i1  ,jd),qt(kdp,i1  ,jd))
	q2 = MIN(MAX(qtmin,q2),qtmax)
	qtmin = MIN(qt(kd,i1p ,jd),qt(kdp,i1p ,jd))
	qtmax = MAX(qt(kd,i1p ,jd),qt(kdp,i1p ,jd))
	q3 = MIN(MAX(qtmin,q3),qtmax)
	qtmin = MIN(qt(kd,i1pp,jd),qt(kdp,i1pp,jd))
	qtmax = MAX(qt(kd,i1pp,jd),qt(kdp,i1pp,jd))
	q4 = MIN(MAX(qtmin,q4),qtmax)
      ENDIF
      qqt = q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4
      IF (monotone) THEN
        qtmin = MIN(q2,q3)
	qtmax = MAX(q2,q3)
	qqt = MIN(MAX(qtmin,qqt),qtmax)
      ENDIF
      
      ! Third
      i1m = idm
      i1 = id
      i1p = idp
      i1pp = idpp
      flip = 1.0d0
      IF (jd .eq. ny) THEN
        jdp = ny
        i1m  = MODULO(i1m - hnx - 1, nx) + 1
        i1   = MODULO(i1 - hnx - 1, nx) + 1
        i1p  = MODULO(i1p - hnx - 1, nx) + 1
        i1pp = MODULO(i1pp - hnx - 1, nx) + 1
        flip = -1.0d0
      ENDIF
      q1 = face1*qu(kdm,i1m ,jdp) + face2*qu(kd,i1m ,jdp) + face3*qu(kdp,i1m ,jdp) + face4*qu(kdpp,i1m ,jdp)
      q2 = face1*qu(kdm,i1  ,jdp) + face2*qu(kd,i1  ,jdp) + face3*qu(kdp,i1  ,jdp) + face4*qu(kdpp,i1  ,jdp)
      q3 = face1*qu(kdm,i1p ,jdp) + face2*qu(kd,i1p ,jdp) + face3*qu(kdp,i1p ,jdp) + face4*qu(kdpp,i1p ,jdp)
      q4 = face1*qu(kdm,i1pp,jdp) + face2*qu(kd,i1pp,jdp) + face3*qu(kdp,i1pp,jdp) + face4*qu(kdpp,i1pp,jdp)
      qqup = (q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4)*flip
      q1 = face1*qv(kdm,i1m ,jdp) + face2*qv(kd,i1m ,jdp) + face3*qv(kdp,i1m ,jdp) + face4*qv(kdpp,i1m ,jdp)
      q2 = face1*qv(kdm,i1  ,jdp) + face2*qv(kd,i1  ,jdp) + face3*qv(kdp,i1  ,jdp) + face4*qv(kdpp,i1  ,jdp)
      q3 = face1*qv(kdm,i1p ,jdp) + face2*qv(kd,i1p ,jdp) + face3*qv(kdp,i1p ,jdp) + face4*qv(kdpp,i1p ,jdp)
      q4 = face1*qv(kdm,i1pp,jdp) + face2*qv(kd,i1pp,jdp) + face3*qv(kdp,i1pp,jdp) + face4*qv(kdpp,i1pp,jdp)
      qqvp = (q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4)*flip
      q1 = face1*qw(kdm,i1m ,jdp) + face2*qw(kd,i1m ,jdp) + face3*qw(kdp,i1m ,jdp) + face4*qw(kdpp,i1m ,jdp)
      q2 = face1*qw(kdm,i1  ,jdp) + face2*qw(kd,i1  ,jdp) + face3*qw(kdp,i1  ,jdp) + face4*qw(kdpp,i1  ,jdp)
      q3 = face1*qw(kdm,i1p ,jdp) + face2*qw(kd,i1p ,jdp) + face3*qw(kdp,i1p ,jdp) + face4*qw(kdpp,i1p ,jdp)
      q4 = face1*qw(kdm,i1pp,jdp) + face2*qw(kd,i1pp,jdp) + face3*qw(kdp,i1pp,jdp) + face4*qw(kdpp,i1pp,jdp)
      qqwp = q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4
      q1 = face1*qt(kdm,i1m ,jdp) + face2*qt(kd,i1m ,jdp) + face3*qt(kdp,i1m ,jdp) + face4*qt(kdpp,i1m ,jdp)
      q2 = face1*qt(kdm,i1  ,jdp) + face2*qt(kd,i1  ,jdp) + face3*qt(kdp,i1  ,jdp) + face4*qt(kdpp,i1  ,jdp)
      q3 = face1*qt(kdm,i1p ,jdp) + face2*qt(kd,i1p ,jdp) + face3*qt(kdp,i1p ,jdp) + face4*qt(kdpp,i1p ,jdp)
      q4 = face1*qt(kdm,i1pp,jdp) + face2*qt(kd,i1pp,jdp) + face3*qt(kdp,i1pp,jdp) + face4*qt(kdpp,i1pp,jdp)
      IF (monotone) THEN
        qtmin = MIN(qt(kd,i1m ,jdp),qt(kdp,i1m ,jdp))
	qtmax = MAX(qt(kd,i1m ,jdp),qt(kdp,i1m ,jdp))
	q1 = MIN(MAX(qtmin,q1),qtmax)
	qtmin = MIN(qt(kd,i1  ,jdp),qt(kdp,i1  ,jdp))
	qtmax = MAX(qt(kd,i1  ,jdp),qt(kdp,i1  ,jdp))
	q2 = MIN(MAX(qtmin,q2),qtmax)
	qtmin = MIN(qt(kd,i1p ,jdp),qt(kdp,i1p ,jdp))
	qtmax = MAX(qt(kd,i1p ,jdp),qt(kdp,i1p ,jdp))
	q3 = MIN(MAX(qtmin,q3),qtmax)
	qtmin = MIN(qt(kd,i1pp,jdp),qt(kdp,i1pp,jdp))
	qtmax = MAX(qt(kd,i1pp,jdp),qt(kdp,i1pp,jdp))
	q4 = MIN(MAX(qtmin,q4),qtmax)
      ENDIF
      qqtp = q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4
      IF (monotone) THEN
        qtmin = MIN(q2,q3)
	qtmax = MAX(q2,q3)
	qqtp = MIN(MAX(qtmin,qqtp),qtmax)
      ENDIF
      
      ! Fourth
      i1m = idm
      i1 = id
      i1p = idp
      i1pp = idpp
      flip = 1.0d0
      IF (jd .ge. ny - 1) THEN
        jdpp = 2*ny + 1 - jdpp
        i1m  = MODULO(i1m - hnx - 1, nx) + 1
        i1   = MODULO(i1 - hnx - 1, nx) + 1
        i1p  = MODULO(i1p - hnx - 1, nx) + 1
        i1pp = MODULO(i1pp - hnx - 1, nx) + 1
        flip = -1.0d0
      ENDIF
      q1 = face1*qu(kdm,i1m ,jdpp) + face2*qu(kd,i1m ,jdpp) + face3*qu(kdp,i1m ,jdpp) + face4*qu(kdpp,i1m ,jdpp)
      q2 = face1*qu(kdm,i1  ,jdpp) + face2*qu(kd,i1  ,jdpp) + face3*qu(kdp,i1  ,jdpp) + face4*qu(kdpp,i1  ,jdpp)
      q3 = face1*qu(kdm,i1p ,jdpp) + face2*qu(kd,i1p ,jdpp) + face3*qu(kdp,i1p ,jdpp) + face4*qu(kdpp,i1p ,jdpp)
      q4 = face1*qu(kdm,i1pp,jdpp) + face2*qu(kd,i1pp,jdpp) + face3*qu(kdp,i1pp,jdpp) + face4*qu(kdpp,i1pp,jdpp)
      qqupp = (q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4)*flip
      q1 = face1*qv(kdm,i1m ,jdpp) + face2*qv(kd,i1m ,jdpp) + face3*qv(kdp,i1m ,jdpp) + face4*qv(kdpp,i1m ,jdpp)
      q2 = face1*qv(kdm,i1  ,jdpp) + face2*qv(kd,i1  ,jdpp) + face3*qv(kdp,i1  ,jdpp) + face4*qv(kdpp,i1  ,jdpp)
      q3 = face1*qv(kdm,i1p ,jdpp) + face2*qv(kd,i1p ,jdpp) + face3*qv(kdp,i1p ,jdpp) + face4*qv(kdpp,i1p ,jdpp)
      q4 = face1*qv(kdm,i1pp,jdpp) + face2*qv(kd,i1pp,jdpp) + face3*qv(kdp,i1pp,jdpp) + face4*qv(kdpp,i1pp,jdpp)
      qqvpp = (q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4)*flip
      q1 = face1*qw(kdm,i1m ,jdpp) + face2*qw(kd,i1m ,jdpp) + face3*qw(kdp,i1m ,jdpp) + face4*qw(kdpp,i1m ,jdpp)
      q2 = face1*qw(kdm,i1  ,jdpp) + face2*qw(kd,i1  ,jdpp) + face3*qw(kdp,i1  ,jdpp) + face4*qw(kdpp,i1  ,jdpp)
      q3 = face1*qw(kdm,i1p ,jdpp) + face2*qw(kd,i1p ,jdpp) + face3*qw(kdp,i1p ,jdpp) + face4*qw(kdpp,i1p ,jdpp)
      q4 = face1*qw(kdm,i1pp,jdpp) + face2*qw(kd,i1pp,jdpp) + face3*qw(kdp,i1pp,jdpp) + face4*qw(kdpp,i1pp,jdpp)
      qqwpp = q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4
      q1 = face1*qt(kdm,i1m ,jdpp) + face2*qt(kd,i1m ,jdpp) + face3*qt(kdp,i1m ,jdpp) + face4*qt(kdpp,i1m ,jdpp)
      q2 = face1*qt(kdm,i1  ,jdpp) + face2*qt(kd,i1  ,jdpp) + face3*qt(kdp,i1  ,jdpp) + face4*qt(kdpp,i1  ,jdpp)
      q3 = face1*qt(kdm,i1p ,jdpp) + face2*qt(kd,i1p ,jdpp) + face3*qt(kdp,i1p ,jdpp) + face4*qt(kdpp,i1p ,jdpp)
      q4 = face1*qt(kdm,i1pp,jdpp) + face2*qt(kd,i1pp,jdpp) + face3*qt(kdp,i1pp,jdpp) + face4*qt(kdpp,i1pp,jdpp)
      IF (monotone) THEN
        qtmin = MIN(qt(kd,i1m ,jdpp),qt(kdp,i1m ,jdpp))
	qtmax = MAX(qt(kd,i1m ,jdpp),qt(kdp,i1m ,jdpp))
	q1 = MIN(MAX(qtmin,q1),qtmax)
	qtmin = MIN(qt(kd,i1  ,jdpp),qt(kdp,i1  ,jdpp))
	qtmax = MAX(qt(kd,i1  ,jdpp),qt(kdp,i1  ,jdpp))
	q2 = MIN(MAX(qtmin,q2),qtmax)
	qtmin = MIN(qt(kd,i1p ,jdpp),qt(kdp,i1p ,jdpp))
	qtmax = MAX(qt(kd,i1p ,jdpp),qt(kdp,i1p ,jdpp))
	q3 = MIN(MAX(qtmin,q3),qtmax)
	qtmin = MIN(qt(kd,i1pp,jdpp),qt(kdp,i1pp,jdpp))
	qtmax = MAX(qt(kd,i1pp,jdpp),qt(kdp,i1pp,jdpp))
	q4 = MIN(MAX(qtmin,q4),qtmax)
      ENDIF
      qqtpp = q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4
      IF (monotone) THEN
        qtmin = MIN(q2,q3)
	qtmax = MAX(q2,q3)
	qqtpp = MIN(MAX(qtmin,qqtpp),qtmax)
      ENDIF

      ! Interpolate in y
      qunew(k,i,j) = qqum*facy1 + qqu*facy2 + qqup*facy3 + qqupp*facy4
      qvnew(k,i,j) = qqvm*facy1 + qqv*facy2 + qqvp*facy3 + qqvpp*facy4
      qwnew(k,i,j) = qqwm*facy1 + qqw*facy2 + qqwp*facy3 + qqwpp*facy4
      qtnew(k,i,j) = qqtm*facy1 + qqt*facy2 + qqtp*facy3 + qqtpp*facy4
      
      IF (monotone) THEN
!        ! Simple limiter for theta advection
!        qtmin = MIN(qt(kd,i1,jd ),qt(kd,i1p,jd ),qt(kdp,i1,jd ),qt(kdp,i1p,jd ), &
!                    qt(kd,i1,jdp),qt(kd,i1p,jdp),qt(kdp,i1,jdp),qt(kdp,i1p,jdp))
!        qtmax = MAX(qt(kd,i1,jd ),qt(kd,i1p,jd ),qt(kdp,i1,jd ),qt(kdp,i1p,jd ), &
!                    qt(kd,i1,jdp),qt(kd,i1p,jdp),qt(kdp,i1,jdp),qt(kdp,i1p,jdp))
!        qtnew(k,i,j) = MIN(MAX(qtmin,qtnew(k,i,j)),qtmax)
        qtmin = MIN(qqt,qqtp)
	qtmax = MAX(qqt,qqtp)
	qtnew(k,i,j) = MIN(MAX(qtmin,qtnew(k,i,j)),qtmax)
      ENDIF



    ENDDO
  ENDDO
ENDDO


END SUBROUTINE lagrangew

! ========================================================

SUBROUTINE lagrangewd(xd,yd,ed,qu,qv,qw,qt,qunew,qvnew,qwnew,qtnew,qdtnew)

! 3D Cubic Lagrange interpolation for the three components of a vector field
! plus a scalar stored at the w-points on the spherical C-grid. Note the `flip' factors
! used when interpolating across the poles
!
! This alternative routine also calculates the vertical derivative of the scalar.
!

USE switches
USE grid
USE util

IMPLICIT NONE

REAL*8, INTENT(IN) :: xd(nzp,nx,ny), yd(nzp,nx,ny), ed(nzp,nx,ny), &
		      qu(nzp,nx,ny), qv(nzp,nx,ny), qw(nzp,nx,ny), qt(nzp,nx,ny)
REAL*8, INTENT(OUT) :: qunew(nzp,nx,ny), qvnew(nzp,nx,ny), qwnew(nzp,nx,ny), &
                       qtnew(nzp,nx,ny), qdtnew(nzp,nx,ny)
INTEGER :: i, j, id, idm, idp, idpp, jd, jdm, jdp, jdpp, &
           i1, i1m, i1p, i1pp, hnx, k, kd, kdm, kdp, kdpp
REAL*8 :: xdd, ydd, edd, d12, d13, d14, d23, d24, d34, &
          denx1, denx2, denx3, denx4, deny1, deny2, deny3, deny4, &
	  dene1, dene2, dene3, dene4, &
	  dd1, dd2, dd3, dd4, &
	  facx1, facx2, facx3, facx4, facy1, facy2, facy3, facy4, &
	  face1, face2, face3, face4, &
	  faced1, faced2, faced3, faced4, &
	  q1, q2, q3, q4, qqum, qqu, qqup, qqupp, &
          qqvm, qqv, qqvp, qqvpp, qqwm, qqw, qqwp, qqwpp, &
          qqtm, qqt, qqtp, qqtpp, flip, x(nx), y(ny), e(nzp), &
	  qqdtm, qqdt, qqdtp, qqdtpp, &
	  qtmin, qtmax, ql1, ql2, ql3, ql4, qd1, qd2, qd3, qd4


if (.not. monotone) print *,'lagrangewd only works monotone !!!!!!!!!!!!!!!!'

! Regular gridded values are at w points
x = xp
y = yp
e = etaw

! Handy value for handling poles
hnx = nx/2

! Some horizontal interpolation factors can be pre-computed on uniform grid

! For longitude direction
d12 = -dx
d13 = -2*dx
d14 = -3*dx
d23 = -dx
d24 = -2*dx
d34 = -dx
denx1 =  d12*d13*d14
denx2 = -d12*d23*d24
denx3 =  d13*d23*d34
denx4 = -d14*d24*d34

! For latitude direction
d12 = -dy
d13 = -2*dy
d14 = -3*dy
d23 = -dy
d24 = -2*dy
d34 = -dy
deny1 =  d12*d13*d14
deny2 = -d12*d23*d24
deny3 =  d13*d23*d34
deny4 = -d14*d24*d34


DO j = 1, ny
  DO i = 1, nx
    DO k = 1, nzp

      ! Find indices of departure point and stencil 
      xdd = xd(k,i,j)
      ydd = yd(k,i,j)
      edd = ed(k,i,j)

      id =   MODULO(FLOOR((xdd-x(1))/dx),nx) + 1
      idm =  MODULO(id-2,nx) + 1
      idp =  MODULO(id,nx) + 1
      idpp = MODULO(idp,nx) + 1
      jd =   FLOOR((ydd-y(1))/dy) + 1
      jdm =  jd - 1
      jdp =  jd + 1
      jdpp = jdp + 1
      CALL findetacell(e,edd,nzp,kd)
      kd = MIN(MAX(kd,1),nz)
      kdm = kd - 1
      kdp = kd + 1
      kdpp = kdp + 1

      ! Factors for x-interpolation
      dd1 = near(xdd - x(idm), twopi)
      dd2 = near(xdd - x(id), twopi)
      dd3 = near(xdd - x(idp), twopi)
      dd4 = near(xdd - x(idpp), twopi)
      facx1 = dd2*dd3*dd4/denx1
      facx2 = dd1*dd3*dd4/denx2
      facx3 = dd1*dd2*dd4/denx3
      facx4 = dd1*dd2*dd3/denx4

      ! Factors for y-interpolation
      dd2 = ydd - (y(1) + jdm*dy)
      dd1 = dd2 + dy
      dd3 = dd2 - dy
      dd4 = dd2 - 2*dy
      facy1 = dd2*dd3*dd4/deny1
      facy2 = dd1*dd3*dd4/deny2
      facy3 = dd1*dd2*dd4/deny3
      facy4 = dd1*dd2*dd3/deny4

      ! Factors for eta-interpolation
      IF (k == 1) THEN
        ! On the bottom boundary; just use boundary values
	! Set vertical derivative to zero
        kdm = 1
        kd = 1
        kdp = 1
        kdpp = 2
        face1 = 0.0d0
        face2 = 0.0d0
        face3 = 1.0d0
        face4 = 0.0d0
	faced1 = 0.0d0
	faced2 = 0.0d0
	faced3 = 0.0d0
	faced4 = 0.0d0
      ELSEIF (k == nzp) THEN
        ! On the top boundary; just use boundary values
	! Set vertical derivative to zero
        kdm = nz
        kd = nzp
        kdp = nzp
        kdpp = nzp
        face1 = 0.0d0
        face2 = 1.0d0
        face3 = 0.0d0
        face4 = 0.0d0
	faced1 = 0.0d0
	faced2 = 0.0d0
	faced3 = 0.0d0
	faced4 = 0.0d0
      ! In lowest eta interval linear interpolation
      ELSEIF (kd == 1) THEN
        kdm = 1
        kd = 1
        kdp = 2
        kdpp = 3
        dd3 = edd - e(kdp)
        d23 = e(kd) - e(kdp)
        face1 = 0.0d0
        face2 = dd3/d23
        face3 = 1.0d0 - face2
        face4 = 0.0d0
	faced1 = 0.0d0
	faced2 = 1.0d0/d23
	faced3 = -faced2
	faced4 = 0.0d0
      ! In highest eta interval linear interpolation
      ELSEIF (kd == nz) THEN
        kdm = nz-1
        kd = nz
        kdp = nzp
        kdpp = nzp
        dd3 = edd - e(kdp)
        d23 = e(kd) - e(kdp)
        face1 = 0.0d0
        face2 = dd3/d23
        face3 = 1.0d0 - face2
        face4 = 0.0d0
	faced1 = 0.0d0
	faced2 = 1.0d0/d23
	faced3 = -faced2
	faced4 = 0.0d0
      ELSE
        ! Otherwise cubic interpolation
        dd1 = edd - e(kdm)
        dd2 = edd - e(kd)
        dd3 = edd - e(kdp)
        dd4 = edd - e(kdpp)
        d12 = e(kdm) - e(kd)
        d13 = e(kdm) - e(kdp)
        d14 = e(kdm) - e(kdpp)
        d23 = e(kd) - e(kdp)
        d24 = e(kd) - e(kdpp)
        d34 = e(kdp) - e(kdpp)
        dene1 = d12*d13*d14
        dene2 = -d12*d23*d24
        dene3 = d13*d23*d34
        dene4 = -d14*d24*d34
        face1 = dd2*dd3*dd4/dene1
        face2 = dd1*dd3*dd4/dene2
        face3 = dd1*dd2*dd4/dene3
        face4 = dd1*dd2*dd3/dene4
	faced1 = (dd2*dd3 + dd2*dd4 + dd3*dd4)/dene1
	faced2 = (dd1*dd3 + dd1*dd4 + dd3*dd4)/dene2
	faced3 = (dd1*dd2 + dd1*dd4 + dd2*dd4)/dene3
	faced4 = (dd1*dd2 + dd1*dd3 + dd2*dd3)/dene4
      ENDIF

      ! Interpolate at four rows
      ! First
      i1m = idm
      i1 = id
      i1p = idp
      i1pp = idpp
      flip = 1.0d0
      IF (jd .le. 1) THEN
        jdm = 1 - jdm
        i1m  = MODULO(i1m - hnx - 1, nx) + 1
        i1   = MODULO(i1 - hnx - 1, nx) + 1
        i1p  = MODULO(i1p - hnx - 1, nx) + 1
        i1pp = MODULO(i1pp - hnx - 1, nx) + 1
        flip = -1.0d0
      ENDIF
      q1 = face1*qu(kdm,i1m ,jdm) + face2*qu(kd,i1m ,jdm) + face3*qu(kdp,i1m ,jdm) + face4*qu(kdpp,i1m ,jdm)
      q2 = face1*qu(kdm,i1  ,jdm) + face2*qu(kd,i1  ,jdm) + face3*qu(kdp,i1  ,jdm) + face4*qu(kdpp,i1  ,jdm)
      q3 = face1*qu(kdm,i1p ,jdm) + face2*qu(kd,i1p ,jdm) + face3*qu(kdp,i1p ,jdm) + face4*qu(kdpp,i1p ,jdm)
      q4 = face1*qu(kdm,i1pp,jdm) + face2*qu(kd,i1pp,jdm) + face3*qu(kdp,i1pp,jdm) + face4*qu(kdpp,i1pp,jdm)
      qqum = (q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4)*flip
      q1 = face1*qv(kdm,i1m ,jdm) + face2*qv(kd,i1m ,jdm) + face3*qv(kdp,i1m ,jdm) + face4*qv(kdpp,i1m ,jdm)
      q2 = face1*qv(kdm,i1  ,jdm) + face2*qv(kd,i1  ,jdm) + face3*qv(kdp,i1  ,jdm) + face4*qv(kdpp,i1  ,jdm)
      q3 = face1*qv(kdm,i1p ,jdm) + face2*qv(kd,i1p ,jdm) + face3*qv(kdp,i1p ,jdm) + face4*qv(kdpp,i1p ,jdm)
      q4 = face1*qv(kdm,i1pp,jdm) + face2*qv(kd,i1pp,jdm) + face3*qv(kdp,i1pp,jdm) + face4*qv(kdpp,i1pp,jdm)
      qqvm = (q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4)*flip
      q1 = face1*qw(kdm,i1m ,jdm) + face2*qw(kd,i1m ,jdm) + face3*qw(kdp,i1m ,jdm) + face4*qw(kdpp,i1m ,jdm)
      q2 = face1*qw(kdm,i1  ,jdm) + face2*qw(kd,i1  ,jdm) + face3*qw(kdp,i1  ,jdm) + face4*qw(kdpp,i1  ,jdm)
      q3 = face1*qw(kdm,i1p ,jdm) + face2*qw(kd,i1p ,jdm) + face3*qw(kdp,i1p ,jdm) + face4*qw(kdpp,i1p ,jdm)
      q4 = face1*qw(kdm,i1pp,jdm) + face2*qw(kd,i1pp,jdm) + face3*qw(kdp,i1pp,jdm) + face4*qw(kdpp,i1pp,jdm)
      qqwm = q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4
      q1 = face1*qt(kdm,i1m ,jdm) + face2*qt(kd,i1m ,jdm) + face3*qt(kdp,i1m ,jdm) + face4*qt(kdpp,i1m ,jdm)
      q2 = face1*qt(kdm,i1  ,jdm) + face2*qt(kd,i1  ,jdm) + face3*qt(kdp,i1  ,jdm) + face4*qt(kdpp,i1  ,jdm)
      q3 = face1*qt(kdm,i1p ,jdm) + face2*qt(kd,i1p ,jdm) + face3*qt(kdp,i1p ,jdm) + face4*qt(kdpp,i1p ,jdm)
      q4 = face1*qt(kdm,i1pp,jdm) + face2*qt(kd,i1pp,jdm) + face3*qt(kdp,i1pp,jdm) + face4*qt(kdpp,i1pp,jdm)
      IF (monotone) THEN
        qtmin = MIN(qt(kd,i1m ,jdm),qt(kdp,i1m ,jdm))
	qtmax = MAX(qt(kd,i1m ,jdm),qt(kdp,i1m ,jdm))
	ql1 = MIN(MAX(qtmin,q1),qtmax)
	qtmin = MIN(qt(kd,i1  ,jdm),qt(kdp,i1  ,jdm))
	qtmax = MAX(qt(kd,i1  ,jdm),qt(kdp,i1  ,jdm))
	ql2 = MIN(MAX(qtmin,q2),qtmax)
	qtmin = MIN(qt(kd,i1p ,jdm),qt(kdp,i1p ,jdm))
	qtmax = MAX(qt(kd,i1p ,jdm),qt(kdp,i1p ,jdm))
	ql3 = MIN(MAX(qtmin,q3),qtmax)
	qtmin = MIN(qt(kd,i1pp,jdm),qt(kdp,i1pp,jdm))
	qtmax = MAX(qt(kd,i1pp,jdm),qt(kdp,i1pp,jdm))
	ql4 = MIN(MAX(qtmin,q4),qtmax)
      ELSE
        ql1 = q1
	ql2 = q2
	ql3 = q3
	ql4 = q4
      ENDIF
      qqtm = ql1*facx1 + ql2*facx2 + ql3*facx3 + ql4*facx4
      IF (monotone) THEN
        qtmin = MIN(ql2,ql3)
	qtmax = MAX(ql2,ql3)
	qqtm = MIN(MAX(qtmin,qqtm),qtmax)
      ENDIF
      qd1 = faced1*qt(kdm,i1m ,jdm) + faced2*qt(kd,i1m ,jdm) + faced3*qt(kdp,i1m ,jdm) + faced4*qt(kdpp,i1m ,jdm)
      qd2 = faced1*qt(kdm,i1  ,jdm) + faced2*qt(kd,i1  ,jdm) + faced3*qt(kdp,i1  ,jdm) + faced4*qt(kdpp,i1  ,jdm)
      qd3 = faced1*qt(kdm,i1p ,jdm) + faced2*qt(kd,i1p ,jdm) + faced3*qt(kdp,i1p ,jdm) + faced4*qt(kdpp,i1p ,jdm)
      qd4 = faced1*qt(kdm,i1pp,jdm) + faced2*qt(kd,i1pp,jdm) + faced3*qt(kdp,i1pp,jdm) + faced4*qt(kdpp,i1pp,jdm)
      IF (ql1 .NE. q1) qd1 = 0.0d0
      IF (ql2 .NE. q2) qd2 = 0.0d0
      IF (ql3 .NE. q3) qd3 = 0.0d0
      IF (ql4 .NE. q4) qd4 = 0.0d0
      qqdtm = qd1*facx1 + qd2*facx2 + qd3*facx3 + qd4*facx4
      
      ! Second
      i1m = idm
      i1 = id
      i1p = idp
      i1pp = idpp
      flip = 1.0d0
      IF (jd .eq. 0) THEN
        jd = 1
        i1m  = MODULO(i1m - hnx - 1, nx) + 1
        i1   = MODULO(i1 - hnx - 1, nx) + 1
        i1p  = MODULO(i1p - hnx - 1, nx) + 1
        i1pp = MODULO(i1pp - hnx - 1, nx) + 1
        flip = -1.0d0
      ENDIF
      q1 = face1*qu(kdm,i1m ,jd) + face2*qu(kd,i1m ,jd) + face3*qu(kdp,i1m ,jd) + face4*qu(kdpp,i1m ,jd)
      q2 = face1*qu(kdm,i1  ,jd) + face2*qu(kd,i1  ,jd) + face3*qu(kdp,i1  ,jd) + face4*qu(kdpp,i1  ,jd)
      q3 = face1*qu(kdm,i1p ,jd) + face2*qu(kd,i1p ,jd) + face3*qu(kdp,i1p ,jd) + face4*qu(kdpp,i1p ,jd)
      q4 = face1*qu(kdm,i1pp,jd) + face2*qu(kd,i1pp,jd) + face3*qu(kdp,i1pp,jd) + face4*qu(kdpp,i1pp,jd)
      qqu = (q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4)*flip
      q1 = face1*qv(kdm,i1m ,jd) + face2*qv(kd,i1m ,jd) + face3*qv(kdp,i1m ,jd) + face4*qv(kdpp,i1m ,jd)
      q2 = face1*qv(kdm,i1  ,jd) + face2*qv(kd,i1  ,jd) + face3*qv(kdp,i1  ,jd) + face4*qv(kdpp,i1  ,jd)
      q3 = face1*qv(kdm,i1p ,jd) + face2*qv(kd,i1p ,jd) + face3*qv(kdp,i1p ,jd) + face4*qv(kdpp,i1p ,jd)
      q4 = face1*qv(kdm,i1pp,jd) + face2*qv(kd,i1pp,jd) + face3*qv(kdp,i1pp,jd) + face4*qv(kdpp,i1pp,jd)
      qqv = (q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4)*flip
      q1 = face1*qw(kdm,i1m ,jd) + face2*qw(kd,i1m ,jd) + face3*qw(kdp,i1m ,jd) + face4*qw(kdpp,i1m ,jd)
      q2 = face1*qw(kdm,i1  ,jd) + face2*qw(kd,i1  ,jd) + face3*qw(kdp,i1  ,jd) + face4*qw(kdpp,i1  ,jd)
      q3 = face1*qw(kdm,i1p ,jd) + face2*qw(kd,i1p ,jd) + face3*qw(kdp,i1p ,jd) + face4*qw(kdpp,i1p ,jd)
      q4 = face1*qw(kdm,i1pp,jd) + face2*qw(kd,i1pp,jd) + face3*qw(kdp,i1pp,jd) + face4*qw(kdpp,i1pp,jd)
      qqw = q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4
      q1 = face1*qt(kdm,i1m ,jd) + face2*qt(kd,i1m ,jd) + face3*qt(kdp,i1m ,jd) + face4*qt(kdpp,i1m ,jd)
      q2 = face1*qt(kdm,i1  ,jd) + face2*qt(kd,i1  ,jd) + face3*qt(kdp,i1  ,jd) + face4*qt(kdpp,i1  ,jd)
      q3 = face1*qt(kdm,i1p ,jd) + face2*qt(kd,i1p ,jd) + face3*qt(kdp,i1p ,jd) + face4*qt(kdpp,i1p ,jd)
      q4 = face1*qt(kdm,i1pp,jd) + face2*qt(kd,i1pp,jd) + face3*qt(kdp,i1pp,jd) + face4*qt(kdpp,i1pp,jd)
      IF (monotone) THEN
        qtmin = MIN(qt(kd,i1m ,jd),qt(kdp,i1m ,jd))
	qtmax = MAX(qt(kd,i1m ,jd),qt(kdp,i1m ,jd))
	ql1 = MIN(MAX(qtmin,q1),qtmax)
	qtmin = MIN(qt(kd,i1  ,jd),qt(kdp,i1  ,jd))
	qtmax = MAX(qt(kd,i1  ,jd),qt(kdp,i1  ,jd))
	ql2 = MIN(MAX(qtmin,q2),qtmax)
	qtmin = MIN(qt(kd,i1p ,jd),qt(kdp,i1p ,jd))
	qtmax = MAX(qt(kd,i1p ,jd),qt(kdp,i1p ,jd))
	ql3 = MIN(MAX(qtmin,q3),qtmax)
	qtmin = MIN(qt(kd,i1pp,jd),qt(kdp,i1pp,jd))
	qtmax = MAX(qt(kd,i1pp,jd),qt(kdp,i1pp,jd))
	ql4 = MIN(MAX(qtmin,q4),qtmax)
      ELSE
        ql1 = q1
	ql2 = q2
	ql3 = q3
	ql4 = q4
      ENDIF
      qqt = ql1*facx1 + ql2*facx2 + ql3*facx3 + ql4*facx4
      IF (monotone) THEN
        qtmin = MIN(ql2,ql3)
	qtmax = MAX(ql2,ql3)
	qqt = MIN(MAX(qtmin,qqt),qtmax)
      ENDIF
      qd1 = faced1*qt(kdm,i1m ,jdm) + faced2*qt(kd,i1m ,jdm) + faced3*qt(kdp,i1m ,jdm) + faced4*qt(kdpp,i1m ,jdm)
      qd2 = faced1*qt(kdm,i1  ,jdm) + faced2*qt(kd,i1  ,jdm) + faced3*qt(kdp,i1  ,jdm) + faced4*qt(kdpp,i1  ,jdm)
      qd3 = faced1*qt(kdm,i1p ,jdm) + faced2*qt(kd,i1p ,jdm) + faced3*qt(kdp,i1p ,jdm) + faced4*qt(kdpp,i1p ,jdm)
      qd4 = faced1*qt(kdm,i1pp,jdm) + faced2*qt(kd,i1pp,jdm) + faced3*qt(kdp,i1pp,jdm) + faced4*qt(kdpp,i1pp,jdm)
      IF (ql1 .NE. q1) qd1 = 0.0d0
      IF (ql2 .NE. q2) qd2 = 0.0d0
      IF (ql3 .NE. q3) qd3 = 0.0d0
      IF (ql4 .NE. q4) qd4 = 0.0d0
      qqdt = qd1*facx1 + qd2*facx2 + qd3*facx3 + qd4*facx4

      ! Third
      i1m = idm
      i1 = id
      i1p = idp
      i1pp = idpp
      flip = 1.0d0
      IF (jd .eq. ny) THEN
        jdp = ny
        i1m  = MODULO(i1m - hnx - 1, nx) + 1
        i1   = MODULO(i1 - hnx - 1, nx) + 1
        i1p  = MODULO(i1p - hnx - 1, nx) + 1
        i1pp = MODULO(i1pp - hnx - 1, nx) + 1
        flip = -1.0d0
      ENDIF
      q1 = face1*qu(kdm,i1m ,jdp) + face2*qu(kd,i1m ,jdp) + face3*qu(kdp,i1m ,jdp) + face4*qu(kdpp,i1m ,jdp)
      q2 = face1*qu(kdm,i1  ,jdp) + face2*qu(kd,i1  ,jdp) + face3*qu(kdp,i1  ,jdp) + face4*qu(kdpp,i1  ,jdp)
      q3 = face1*qu(kdm,i1p ,jdp) + face2*qu(kd,i1p ,jdp) + face3*qu(kdp,i1p ,jdp) + face4*qu(kdpp,i1p ,jdp)
      q4 = face1*qu(kdm,i1pp,jdp) + face2*qu(kd,i1pp,jdp) + face3*qu(kdp,i1pp,jdp) + face4*qu(kdpp,i1pp,jdp)
      qqup = (q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4)*flip
      q1 = face1*qv(kdm,i1m ,jdp) + face2*qv(kd,i1m ,jdp) + face3*qv(kdp,i1m ,jdp) + face4*qv(kdpp,i1m ,jdp)
      q2 = face1*qv(kdm,i1  ,jdp) + face2*qv(kd,i1  ,jdp) + face3*qv(kdp,i1  ,jdp) + face4*qv(kdpp,i1  ,jdp)
      q3 = face1*qv(kdm,i1p ,jdp) + face2*qv(kd,i1p ,jdp) + face3*qv(kdp,i1p ,jdp) + face4*qv(kdpp,i1p ,jdp)
      q4 = face1*qv(kdm,i1pp,jdp) + face2*qv(kd,i1pp,jdp) + face3*qv(kdp,i1pp,jdp) + face4*qv(kdpp,i1pp,jdp)
      qqvp = (q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4)*flip
      q1 = face1*qw(kdm,i1m ,jdp) + face2*qw(kd,i1m ,jdp) + face3*qw(kdp,i1m ,jdp) + face4*qw(kdpp,i1m ,jdp)
      q2 = face1*qw(kdm,i1  ,jdp) + face2*qw(kd,i1  ,jdp) + face3*qw(kdp,i1  ,jdp) + face4*qw(kdpp,i1  ,jdp)
      q3 = face1*qw(kdm,i1p ,jdp) + face2*qw(kd,i1p ,jdp) + face3*qw(kdp,i1p ,jdp) + face4*qw(kdpp,i1p ,jdp)
      q4 = face1*qw(kdm,i1pp,jdp) + face2*qw(kd,i1pp,jdp) + face3*qw(kdp,i1pp,jdp) + face4*qw(kdpp,i1pp,jdp)
      qqwp = q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4
      q1 = face1*qt(kdm,i1m ,jdp) + face2*qt(kd,i1m ,jdp) + face3*qt(kdp,i1m ,jdp) + face4*qt(kdpp,i1m ,jdp)
      q2 = face1*qt(kdm,i1  ,jdp) + face2*qt(kd,i1  ,jdp) + face3*qt(kdp,i1  ,jdp) + face4*qt(kdpp,i1  ,jdp)
      q3 = face1*qt(kdm,i1p ,jdp) + face2*qt(kd,i1p ,jdp) + face3*qt(kdp,i1p ,jdp) + face4*qt(kdpp,i1p ,jdp)
      q4 = face1*qt(kdm,i1pp,jdp) + face2*qt(kd,i1pp,jdp) + face3*qt(kdp,i1pp,jdp) + face4*qt(kdpp,i1pp,jdp)
      IF (monotone) THEN
        qtmin = MIN(qt(kd,i1m ,jdp),qt(kdp,i1m ,jdp))
	qtmax = MAX(qt(kd,i1m ,jdp),qt(kdp,i1m ,jdp))
	ql1 = MIN(MAX(qtmin,q1),qtmax)
	qtmin = MIN(qt(kd,i1  ,jdp),qt(kdp,i1  ,jdp))
	qtmax = MAX(qt(kd,i1  ,jdp),qt(kdp,i1  ,jdp))
	ql2 = MIN(MAX(qtmin,q2),qtmax)
	qtmin = MIN(qt(kd,i1p ,jdp),qt(kdp,i1p ,jdp))
	qtmax = MAX(qt(kd,i1p ,jdp),qt(kdp,i1p ,jdp))
	ql3 = MIN(MAX(qtmin,q3),qtmax)
	qtmin = MIN(qt(kd,i1pp,jdp),qt(kdp,i1pp,jdp))
	qtmax = MAX(qt(kd,i1pp,jdp),qt(kdp,i1pp,jdp))
	ql4 = MIN(MAX(qtmin,q4),qtmax)
      ELSE
        ql1 = q1
	ql2 = q2
	ql3 = q3
	ql4 = q4
      ENDIF
      qqtp = ql1*facx1 + ql2*facx2 + ql3*facx3 + ql4*facx4
      IF (monotone) THEN
        qtmin = MIN(ql2,ql3)
	qtmax = MAX(ql2,ql3)
	qqtp = MIN(MAX(qtmin,qqtp),qtmax)
      ENDIF
      qd1 = faced1*qt(kdm,i1m ,jdm) + faced2*qt(kd,i1m ,jdm) + faced3*qt(kdp,i1m ,jdm) + faced4*qt(kdpp,i1m ,jdm)
      qd2 = faced1*qt(kdm,i1  ,jdm) + faced2*qt(kd,i1  ,jdm) + faced3*qt(kdp,i1  ,jdm) + faced4*qt(kdpp,i1  ,jdm)
      qd3 = faced1*qt(kdm,i1p ,jdm) + faced2*qt(kd,i1p ,jdm) + faced3*qt(kdp,i1p ,jdm) + faced4*qt(kdpp,i1p ,jdm)
      qd4 = faced1*qt(kdm,i1pp,jdm) + faced2*qt(kd,i1pp,jdm) + faced3*qt(kdp,i1pp,jdm) + faced4*qt(kdpp,i1pp,jdm)
      IF (ql1 .NE. q1) qd1 = 0.0d0
      IF (ql2 .NE. q2) qd2 = 0.0d0
      IF (ql3 .NE. q3) qd3 = 0.0d0
      IF (ql4 .NE. q4) qd4 = 0.0d0
      qqdtp = qd1*facx1 + qd2*facx2 + qd3*facx3 + qd4*facx4

      ! Fourth
      i1m = idm
      i1 = id
      i1p = idp
      i1pp = idpp
      flip = 1.0d0
      IF (jd .ge. ny - 1) THEN
        jdpp = 2*ny + 1 - jdpp
        i1m  = MODULO(i1m - hnx - 1, nx) + 1
        i1   = MODULO(i1 - hnx - 1, nx) + 1
        i1p  = MODULO(i1p - hnx - 1, nx) + 1
        i1pp = MODULO(i1pp - hnx - 1, nx) + 1
        flip = -1.0d0
      ENDIF
      q1 = face1*qu(kdm,i1m ,jdpp) + face2*qu(kd,i1m ,jdpp) + face3*qu(kdp,i1m ,jdpp) + face4*qu(kdpp,i1m ,jdpp)
      q2 = face1*qu(kdm,i1  ,jdpp) + face2*qu(kd,i1  ,jdpp) + face3*qu(kdp,i1  ,jdpp) + face4*qu(kdpp,i1  ,jdpp)
      q3 = face1*qu(kdm,i1p ,jdpp) + face2*qu(kd,i1p ,jdpp) + face3*qu(kdp,i1p ,jdpp) + face4*qu(kdpp,i1p ,jdpp)
      q4 = face1*qu(kdm,i1pp,jdpp) + face2*qu(kd,i1pp,jdpp) + face3*qu(kdp,i1pp,jdpp) + face4*qu(kdpp,i1pp,jdpp)
      qqupp = (q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4)*flip
      q1 = face1*qv(kdm,i1m ,jdpp) + face2*qv(kd,i1m ,jdpp) + face3*qv(kdp,i1m ,jdpp) + face4*qv(kdpp,i1m ,jdpp)
      q2 = face1*qv(kdm,i1  ,jdpp) + face2*qv(kd,i1  ,jdpp) + face3*qv(kdp,i1  ,jdpp) + face4*qv(kdpp,i1  ,jdpp)
      q3 = face1*qv(kdm,i1p ,jdpp) + face2*qv(kd,i1p ,jdpp) + face3*qv(kdp,i1p ,jdpp) + face4*qv(kdpp,i1p ,jdpp)
      q4 = face1*qv(kdm,i1pp,jdpp) + face2*qv(kd,i1pp,jdpp) + face3*qv(kdp,i1pp,jdpp) + face4*qv(kdpp,i1pp,jdpp)
      qqvpp = (q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4)*flip
      q1 = face1*qw(kdm,i1m ,jdpp) + face2*qw(kd,i1m ,jdpp) + face3*qw(kdp,i1m ,jdpp) + face4*qw(kdpp,i1m ,jdpp)
      q2 = face1*qw(kdm,i1  ,jdpp) + face2*qw(kd,i1  ,jdpp) + face3*qw(kdp,i1  ,jdpp) + face4*qw(kdpp,i1  ,jdpp)
      q3 = face1*qw(kdm,i1p ,jdpp) + face2*qw(kd,i1p ,jdpp) + face3*qw(kdp,i1p ,jdpp) + face4*qw(kdpp,i1p ,jdpp)
      q4 = face1*qw(kdm,i1pp,jdpp) + face2*qw(kd,i1pp,jdpp) + face3*qw(kdp,i1pp,jdpp) + face4*qw(kdpp,i1pp,jdpp)
      qqwpp = q1*facx1 + q2*facx2 + q3*facx3 + q4*facx4
      q1 = face1*qt(kdm,i1m ,jdpp) + face2*qt(kd,i1m ,jdpp) + face3*qt(kdp,i1m ,jdpp) + face4*qt(kdpp,i1m ,jdpp)
      q2 = face1*qt(kdm,i1  ,jdpp) + face2*qt(kd,i1  ,jdpp) + face3*qt(kdp,i1  ,jdpp) + face4*qt(kdpp,i1  ,jdpp)
      q3 = face1*qt(kdm,i1p ,jdpp) + face2*qt(kd,i1p ,jdpp) + face3*qt(kdp,i1p ,jdpp) + face4*qt(kdpp,i1p ,jdpp)
      q4 = face1*qt(kdm,i1pp,jdpp) + face2*qt(kd,i1pp,jdpp) + face3*qt(kdp,i1pp,jdpp) + face4*qt(kdpp,i1pp,jdpp)
      IF (monotone) THEN
        qtmin = MIN(qt(kd,i1m ,jdpp),qt(kdp,i1m ,jdpp))
	qtmax = MAX(qt(kd,i1m ,jdpp),qt(kdp,i1m ,jdpp))
	ql1 = MIN(MAX(qtmin,q1),qtmax)
	qtmin = MIN(qt(kd,i1  ,jdpp),qt(kdp,i1  ,jdpp))
	qtmax = MAX(qt(kd,i1  ,jdpp),qt(kdp,i1  ,jdpp))
	ql2 = MIN(MAX(qtmin,q2),qtmax)
	qtmin = MIN(qt(kd,i1p ,jdpp),qt(kdp,i1p ,jdpp))
	qtmax = MAX(qt(kd,i1p ,jdpp),qt(kdp,i1p ,jdpp))
	ql3 = MIN(MAX(qtmin,q3),qtmax)
	qtmin = MIN(qt(kd,i1pp,jdpp),qt(kdp,i1pp,jdpp))
	qtmax = MAX(qt(kd,i1pp,jdpp),qt(kdp,i1pp,jdpp))
	ql4 = MIN(MAX(qtmin,q4),qtmax)
      ELSE
        ql1 = q1
	ql2 = q2
	ql3 = q3
	ql4 = q4
      ENDIF
      qqtpp = ql1*facx1 + ql2*facx2 + ql3*facx3 + ql4*facx4
      IF (monotone) THEN
        qtmin = MIN(ql2,ql3)
	qtmax = MAX(ql2,ql3)
	qqtpp = MIN(MAX(qtmin,qqtpp),qtmax)
      ENDIF
      qd1 = faced1*qt(kdm,i1m ,jdm) + faced2*qt(kd,i1m ,jdm) + faced3*qt(kdp,i1m ,jdm) + faced4*qt(kdpp,i1m ,jdm)
      qd2 = faced1*qt(kdm,i1  ,jdm) + faced2*qt(kd,i1  ,jdm) + faced3*qt(kdp,i1  ,jdm) + faced4*qt(kdpp,i1  ,jdm)
      qd3 = faced1*qt(kdm,i1p ,jdm) + faced2*qt(kd,i1p ,jdm) + faced3*qt(kdp,i1p ,jdm) + faced4*qt(kdpp,i1p ,jdm)
      qd4 = faced1*qt(kdm,i1pp,jdm) + faced2*qt(kd,i1pp,jdm) + faced3*qt(kdp,i1pp,jdm) + faced4*qt(kdpp,i1pp,jdm)
      IF (ql1 .NE. q1) qd1 = 0.0d0
      IF (ql2 .NE. q2) qd2 = 0.0d0
      IF (ql3 .NE. q3) qd3 = 0.0d0
      IF (ql4 .NE. q4) qd4 = 0.0d0
      qqdtpp = qd1*facx1 + qd2*facx2 + qd3*facx3 + qd4*facx4

      ! Interpolate in y
      qunew(k,i,j) = qqum*facy1 + qqu*facy2 + qqup*facy3 + qqupp*facy4
      qvnew(k,i,j) = qqvm*facy1 + qqv*facy2 + qqvp*facy3 + qqvpp*facy4
      qwnew(k,i,j) = qqwm*facy1 + qqw*facy2 + qqwp*facy3 + qqwpp*facy4
      qtnew(k,i,j) = qqtm*facy1 + qqt*facy2 + qqtp*facy3 + qqtpp*facy4
      qdtnew(k,i,j) = qqdtm*facy1 + qqdt*facy2 + qqdtp*facy3 + qqdtpp*facy4
      
      IF (monotone) THEN
!        ! Simple limiter for theta advection
!        qtmin = MIN(qt(kd,i1,jd ),qt(kd,i1p,jd ),qt(kdp,i1,jd ),qt(kdp,i1p,jd ), &
!                    qt(kd,i1,jdp),qt(kd,i1p,jdp),qt(kdp,i1,jdp),qt(kdp,i1p,jdp))
!        qtmax = MAX(qt(kd,i1,jd ),qt(kd,i1p,jd ),qt(kdp,i1,jd ),qt(kdp,i1p,jd ), &
!                    qt(kd,i1,jdp),qt(kd,i1p,jdp),qt(kdp,i1,jdp),qt(kdp,i1p,jdp))
!        qtnew(k,i,j) = MIN(MAX(qtmin,qtnew(k,i,j)),qtmax)
        qtmin = MIN(qqt,qqtp)
	qtmax = MAX(qqt,qqtp)
	qtnew(k,i,j) = MIN(MAX(qtmin,qtnew(k,i,j)),qtmax)
      ENDIF

    ENDDO
  ENDDO
ENDDO


END SUBROUTINE lagrangewd

! ========================================================

SUBROUTINE grad(q,gradqx,gradqy,gradqe)

! Compute the three components of the gradient of a scalar field q.
! q is stored at p points; the components of the gradient are stored
! at u, v and w points. The y-component of the gradient is not computed
! at the poles. If required, it must be computed by calling polar after
! this routine. Similarly, the eta-component of the gradient is extrapolated
! to the bottom and top boundaries.

USE grid

IMPLICIT NONE
REAL*8, INTENT(IN) :: q(nz,nx,ny)
REAL*8, INTENT(OUT) :: gradqx(nz,nx,ny), gradqy(nz,nx,nyp), gradqe(nzp,nx,ny)
INTEGER :: i, im, j, jm, k, km, kp, nzm
REAL*8 :: temp(nz,nx,ny), h1, h2

nzm = nz - 1

! Vertical derivative of q
DO j = 1, ny
  DO i = 1, nx
    ! dq/dr at w points
    DO k = 2, nz
      km = k - 1
      gradqe(k,i,j) = (q(k,i,j) - q(km,i,j))/(rp(k,i,j) - rp(km,i,j))
    ENDDO
    ! Average to p points; linear extrapolation to level 1, linear to level nz
    temp(1,i,j) = 0.5d0*(3.0d0*gradqe(2,i,j) - gradqe(3,i,j))
    DO k = 2, nz-1
      kp = k + 1
      temp(k,i,j) = above1(k)*gradqe(kp,i,j) + below1(k)*gradqe(k,i,j)
    ENDDO
    temp(nz,i,j) = 0.5d0*(3.0d0*gradqe(nz,i,j) - gradqe(nz-1,i,j))
  ENDDO
ENDDO
! Constant extrapolation to bottom and top boundaries
gradqe(1,:,:) = gradqe(2,:,:)
gradqe(nzp,:,:) = gradqe(nz,:,:)

! Horizontal derivatives
DO k = 1, nz
  DO i = 1, nx
    im = i - 1
    IF (i == 1) im = nx
    j = 1
    h1 = 0.5d0*(rp(k,i,j) + rp(k,im,j))*cosp(j)
    gradqx(k,i,j) = ( (q(k,i,j) - q(k,im,j)) &
                    - (rp(k,i,j) - rp(k,im,j))*0.5d0*(temp(k,i,j) + temp(k,im,j)) ) &
                                              /(dx*h1)
    DO j = 2, ny
      jm = j - 1
      h1 = 0.5d0*(rp(k,i,j) + rp(k,im,j))*cosp(j)
      gradqx(k,i,j) = ( (q(k,i,j) - q(k,im,j)) &
                      - (rp(k,i,j) - rp(k,im,j))*0.5d0*(temp(k,i,j) + temp(k,im,j)) ) &
                                                /(dx*h1)
      h2 = 0.5d0*(rp(k,i,j) + rp(k,i,jm))
      gradqy(k,i,j) = ( (q(k,i,j) - q(k,i,jm)) &
                      - (rp(k,i,j) - rp(k,i,jm))*0.5d0*(temp(k,i,j) + temp(k,i,jm)) ) &
                                                /(dy*h2)
    ENDDO
  ENDDO
ENDDO
! Set northward derivative to zero at poles to avoid NaN in dummy computations
gradqy(:,:,1) = 0.0d0
gradqy(:,:,nyp) = 0.0d0

END SUBROUTINE grad

! ========================================================

SUBROUTINE gradpi(q,gradqx,gradqy,gradqe)

! Compute the three components of the gradient of a scalar field q.
! q is stored at p points; the components of the gradient are stored
! at u, v and w points. The y-component of the gradient is not computed
! at the poles. If required, it must be computed by calling polar after
! this routine. Similarly, the eta-component of the gradient is extrapolated
! to the bottom and top boundaries.
!
! Modified version of routine grad that uses extrapolation of d/dr (log q)
! to the top and bottom full levels, on the assumption that q is exner.


USE grid

IMPLICIT NONE
REAL*8, INTENT(IN) :: q(nz,nx,ny)
REAL*8, INTENT(OUT) :: gradqx(nz,nx,ny), gradqy(nz,nx,nyp), gradqe(nzp,nx,ny)
INTEGER :: i, im, j, jm, k, km, kp, nzm
REAL*8 :: temp(nz,nx,ny), h1, h2

nzm = nz - 1

! Vertical derivative of q
DO j = 1, ny
  DO i = 1, nx
    ! dq/dr at w points
    DO k = 2, nz
      km = k - 1
      gradqe(k,i,j) = (q(k,i,j) - q(km,i,j))/(rp(k,i,j) - rp(km,i,j))
    ENDDO
    ! Average to p points; constant extrapolation of d(log(q))/dr to level 1 and nz
    DO k = 2, nz-1
      kp = k + 1
      temp(k,i,j) = above1(k)*gradqe(kp,i,j) + below1(k)*gradqe(k,i,j)
    ENDDO
    temp(1,i,j) = temp(2,i,j)*q(1,i,j)/q(2,i,j)
    temp(nz,i,j) = temp(nzm,i,j)*q(nz,i,j)/q(nzm,i,j)
  ENDDO
ENDDO
! Constant extrapolation to bottom and top boundaries
gradqe(1,:,:) = gradqe(2,:,:)
gradqe(nzp,:,:) = gradqe(nz,:,:)

! Horizontal derivatives
DO k = 1, nz
  DO i = 1, nx
    im = i - 1
    IF (i == 1) im = nx
    j = 1
    h1 = 0.5d0*(rp(k,i,j) + rp(k,im,j))*cosp(j)
    gradqx(k,i,j) = ( (q(k,i,j) - q(k,im,j)) &
                    - (rp(k,i,j) - rp(k,im,j))*0.5d0*(temp(k,i,j) + temp(k,im,j)) ) &
                                              /(dx*h1)
    DO j = 2, ny
      jm = j - 1
      h1 = 0.5d0*(rp(k,i,j) + rp(k,im,j))*cosp(j)
      gradqx(k,i,j) = ( (q(k,i,j) - q(k,im,j)) &
                      - (rp(k,i,j) - rp(k,im,j))*0.5d0*(temp(k,i,j) + temp(k,im,j)) ) &
                                                /(dx*h1)
      h2 = 0.5d0*(rp(k,i,j) + rp(k,i,jm))
      gradqy(k,i,j) = ( (q(k,i,j) - q(k,i,jm)) &
                      - (rp(k,i,j) - rp(k,i,jm))*0.5d0*(temp(k,i,j) + temp(k,i,jm)) ) &
                                                /(dy*h2)
    ENDDO
  ENDDO
ENDDO
! Set northward derivative to zero at poles to avoid NaN in dummy computations
gradqy(:,:,1) = 0.0d0
gradqy(:,:,nyp) = 0.0d0

END SUBROUTINE gradpi

! ========================================================

SUBROUTINE wtopnbottom

! Compute the w field at the top and bottom boundaries

USE state

IMPLICIT NONE

INTEGER :: i, im, ip, i1, j, jm, jp, hnx
REAL*8 :: h1, h2, h3, d1xi3, d2xi3, &
          tempu1(nx,ny), tempv1(nx,nyp), tempu2(nx,ny), tempv2(nx,nyp)

hnx = nx/2

! w = - v . grad (r_surf).
! First compute the contributions to v . grad(r)
! at the surface underneath u and v points, with
! constant extrapolation of u and v to the surface.
DO i = 1, nx
  im = i - 1
  IF (i == 1) im = nx
  ! For polar points i1 is index of the point on the
  ! other side of the pole
  i1 = MODULO(i - hnx - 1, nx) + 1
  DO j = 1, ny
    h1 = 0.5d0*(rsurf(i,j) + rsurf(im,j))*cosp(j)
    d1xi3 = (rsurf(i,j) - rsurf(im,j))/dx
    tempu1(i,j) = u(1,i,j)*d1xi3/h1
    h1 = 0.5d0*(rw(nzp,i,j) + rw(nzp,im,j))*cosp(j)
    d1xi3 = (rw(nzp,i,j) - rw(nzp,im,j))/dx
    tempu2(i,j) = u(nz,i,j)*d1xi3/h1
  ENDDO
  j = 1
  jm = j
  h2 = 0.5d0*(rsurf(i,j) + rsurf(i1,jm))
  d2xi3 = (rsurf(i,j) - rsurf(i1,jm))/dy
  tempv1(i,j) = v(1,i,j)*d2xi3/h2
  h2 = 0.5d0*(rw(nzp,i,j) + rw(nzp,i1,jm))
  d2xi3 = (rw(nzp,i,j) - rw(nzp,i1,jm))/dy
  tempv2(i,j) = v(nz,i,j)*d2xi3/h2
  DO j = 2, ny
    jm = j - 1
    h2 = 0.5d0*(rsurf(i,j) + rsurf(i,jm))
    d2xi3 = (rsurf(i,j) - rsurf(i,jm))/dy
    tempv1(i,j) = v(1,i,j)*d2xi3/h2
    h2 = 0.5d0*(rw(nzp,i,j) + rw(nzp,i,jm))
    d2xi3 = (rw(nzp,i,j) - rw(nzp,i,jm))/dy
    tempv2(i,j) = v(nz,i,j)*d2xi3/h2
  ENDDO
  j = nyp
  jm = j - 1
  h2 = 0.5d0*(rsurf(i,jm) + rsurf(i1,jm))
  d2xi3 = (rsurf(i1,jm) - rsurf(i,jm))/dy
  tempv1(i,j) = v(1,i,j)*d2xi3/h2
  h2 = 0.5d0*(rw(nzp,i,jm) + rw(nzp,i1,jm))
  d2xi3 = (rw(nzp,i1,jm) - rw(nzp,i,jm))/dy
  tempv2(i,j) = v(nz,i,j)*d2xi3/h2
ENDDO

! Now average to w points and compute w
h3 = 1.0d0
DO j = 1, ny
  jp = j + 1
  DO i = 1, nx
    ip = i + 1
    IF (i == nx) ip = 1
    w(1  ,i,j) = h3*0.5d0*(tempu1(i,j) + tempu1(ip,j)   &
		         + tempv1(i,j) + tempv1(i,jp))
    w(nzp,i,j) = h3*0.5d0*(tempu2(i,j) + tempu2(ip,j)   &
		         + tempv2(i,j) + tempv2(i,jp))
  ENDDO
ENDDO


! Top boundary is flat: w = 0.
! w(nzp,:,:) = 0.0d0


END SUBROUTINE wtopnbottom

! ========================================================

SUBROUTINE coriolis(u,v,w,rho,cu,cv,cw)

! Evaluate Coriolis terms: -2 \Omega \times \mathbf{u}
!
! Note that no coriolis term is evaluated at the poles.
! Constant extrapolation is used for the bottom and top boundaries.

! *** Could save a few flops by removing factors of 2.0d0 in one
! place and factors of 0.5d0 in another ***

! *** Note a more natural and symmetrical form of the Coriolis terms is possible
! if we store the components of the Earth's rotation vector at the corresponding
! u, v and w points ***


USE constants
USE grid

IMPLICIT NONE

REAL*8, INTENT(IN) :: u(nz,nx,ny), v(nz,nx,nyp), w(nzp,nx,ny), rho(nz,nx,ny)
REAL*8, INTENT(out) :: cu(nz,nx,ny), cv(nz,nx,nyp), cw(nzp,nx,ny)
INTEGER :: i, im, ip, j, jm, jp, k, km, kp
REAL*8 :: mfu(nz,nx,ny), mfv(nz,nx,nyp), mfw(nzp,nx,ny), &
          tempu(nz,nx,ny), tempv(nz,nx,ny), tempw(nz,nx,ny)
REAL*8 :: h1, h2, h12, dr, rhobar, f1, f2, f3, &
          mfubar, mfvbar, mfwbar, r


! First compute mass fluxes
DO j = 1, ny
  DO i = 1, nx
    im = i - 1
    IF (i == 1) im = nx
    DO k = 1, nz
      rhobar = 0.5d0*(rho(k,im,j) + rho(k,i,j))
      mfu(k,i,j) = rhobar*u(k,i,j)*areaw(k,i,j)
    ENDDO
  ENDDO
ENDDO
! Zero northward mass flux at poles
mfv(:,:,1) = 0.0d0
mfv(:,:,nyp) = 0.0d0
DO j = 2, ny
  jm = j - 1
  DO i = 1, nx
    DO k = 1, nz
      rhobar = 0.5d0*(rho(k,i,jm) + rho(k,i,j))
      mfv(k,i,j) = rhobar*v(k,i,j)*areas(k,i,j)
    ENDDO
  ENDDO
ENDDO
DO j = 1, ny
  DO i = 1, nx
    DO k = 2, nz
      km = k - 1
      rhobar =  below2(k)*rho(km,i,j) + above2(k)*rho(k,i,j)
      mfw(k,i,j) = rhobar*w(k,i,j)*areab(k,i,j)
    ENDDO
    k = 1
    rhobar = rho(k,i,j)
    mfw(k,i,j) = rhobar*w(k,i,j)*areab(k,i,j)
    k = nzp
    mfw(k,i,j) = 0.0d0
  ENDDO
ENDDO

! Compute required terms at p points
DO j = 1, ny
  jp = j + 1
  DO i = 1, nx
    ip = i + 1
    IF (i == nx) ip = 1
    DO k = 1, nz
      kp = k + 1
      dr = rw(kp,i,j) - rw(k,i,j)
      mfubar = 0.5d0*(mfu(k,i,j) + mfu(k,ip,j))
      mfvbar = 0.5d0*(mfv(k,i,j) + mfv(k,i,jp))
      mfwbar = 0.5d0*(mfw(k,i,j) + mfw(kp,i,j))
      h1 = rp(k,i,j)*cosp(j)
      f1 = coriol1(i,j)/(rho(k,i,j)*h1*dx)
      h2 = rp(k,i,j)
      f2 = coriol2(i,j)/(rho(k,i,j)*h2*dy)
      f3 = coriol3(i,j)/(rho(k,i,j)*dr) ! h3 is included in dr
      tempu(k,i,j) = mfvbar*f3 - mfwbar*f2
      tempv(k,i,j) = mfwbar*f1 - mfubar*f3
      tempw(k,i,j) = mfubar*f2 - mfvbar*f1
    ENDDO
  ENDDO
ENDDO

! Average to u, v, and w points and complete calculation
DO j = 1, ny
  DO i = 1, nx
    im = i - 1
    IF (i == 1) im = nx
    DO k = 1, nz
      h1 = 0.5d0*(rp(k,i,j) + rp(k,im,j))*cosp(j)
      cu(k,i,j) = 0.5d0*(tempu(k,i,j) + tempu(k,im,j))/(h1*dx)
    ENDDO
  ENDDO
ENDDO
DO j = 2, ny
  jm = j - 1
  DO i = 1, nx
    DO k = 1, nz
      h2 = 0.5d0*(rp(k,i,j) + rp(k,i,jm))
      cv(k,i,j) = 0.5d0*(tempv(k,i,j) + tempv(k,i,jm))/(h2*dy)
    ENDDO
  ENDDO
ENDDO
DO j = 1, ny
  DO i = 1, nx
    DO k = 2, nz
      km = k - 1
      dr = rp(k,i,j) - rp(km,i,j)
      cw(k,i,j) = 0.5d0*(tempw(k,i,j) + tempw(km,i,j))/dr
    ENDDO
  ENDDO
ENDDO
! Constant extrapolation to bottom and top boundaries
cw(1  ,:,:) = cw(2 ,:,:)
cw(nzp,:,:) = cw(nz,:,:)

! Set polar values to zero to avoid problems with dummy computations
cv(:,:,1  ) = 0.0d0
cv(:,:,nyp) = 0.0d0


END SUBROUTINE coriolis

! ========================================================

SUBROUTINE findexner(rho,theta,exner)

! Compute the Exner presure from density rho and
! potential temperature theta

USE constants
USE grid

IMPLICIT NONE

REAL*8, INTENT(IN) :: rho(nz,nx,ny), theta(nzp,nx,ny)
REAL*8, INTENT(OUT) :: exner(nz,nx,ny)
INTEGER :: i, j, k, kp
REAL*8 :: rbyp, thetabar, a, b


rbyp = gascon/p00

DO k = 1, nz
  kp = k + 1
  a = above1(k)
  b = below1(k)
  exner(k,:,:) = (rbyp*rho(k,:,:)*(a*theta(kp,:,:) + b*theta(k,:,:)))**kbyonemk
ENDDO

END SUBROUTINE findexner

! ========================================================

SUBROUTINE divergence(u,v,etadot,div)

! Compute the 3D divergence of the velocity field expressed
! in terms of u, v and etadot

USE constants
USE grid

IMPLICIT NONE
REAL*8, INTENT(IN) :: u(nz,nx,ny), v(nz,nx,nyp), etadot(nzp,nx,ny)
REAL*8, INTENT(OUT) :: div(nz,nx,ny)
INTEGER :: i, ip, j, jp, k, km, kp
REAL*8 :: vfu(nz,nx,ny), vfv(nz,nx,nyp), vfw(nzp,nx,ny), drdeta


! First compute volume fluxes
vfu = u*areaw
vfv = v*areas
! Zero volume flux across bottom and top boundaries
vfw(1,:,:) = 0.0d0
vfw(nzp,:,:) = 0.0d0
DO j = 1, ny
  DO i = 1, nx
    DO k = 2, nz
      km = k - 1
      drdeta = (rp(k,i,j) - rp(km,i,j))/(etap(k) - etap(km))
      vfw(k,i,j) = drdeta*etadot(k,i,j)*areab(k,i,j)
    ENDDO
  ENDDO
ENDDO

! Now compute divergence
DO j = 1, ny
  jp = j + 1
  DO i = 1, nx
    ip = i + 1
    IF (i == nx) ip = 1
    DO k = 1, nz
      kp = k + 1
      div(k,i,j) = ( (vfu(k,ip,j) - vfu(k,i,j))  &
                   + (vfv(k,i,jp) - vfv(k,i,j))  &
		   + (vfw(kp,i,j) - vfw(k,i,j)) )  &
		   / volume(k,i,j)
    ENDDO
  ENDDO
ENDDO


END SUBROUTINE divergence

! ========================================================

SUBROUTINE massdiv(u,v,w,rho,div)

! Compute the 3D divergence of the mass flux, where velocity field 
! is expressed in terms of u, v and w, and w approx drdeta*etadot is assumed

USE constants
USE grid

IMPLICIT NONE
REAL*8, INTENT(IN) :: u(nz,nx,ny), v(nz,nx,nyp), w(nzp,nx,ny), rho(nz,nx,ny)
REAL*8, INTENT(OUT) :: div(nz,nx,ny)
INTEGER :: i, im, ip, j, jm, jp, k, km, kp
REAL*8 :: rhobar
REAL*8 :: mfu(nz,nx,ny), mfv(nz,nx,nyp), mfw(nzp,nx,ny)



print *,'*** Need alpha factors in routine massdiv ***'

! First compute mass fluxes
DO j = 1, ny
  DO i = 1, nx
    im = i - 1
    IF (i == 1) im = nx
    DO k = 1,nz
      rhobar = 0.5d0*(rho(k,i,j) + rho(k,im,j))
      mfu(k,i,j) = rhobar*areaw(k,i,j)*u(k,i,j)
    ENDDO
  ENDDO
ENDDO
! Zero northward mass flux at poles
mfv(:,:,1) = 0.0d0
mfv(:,:,nyp) = 0.0d0
DO j = 2, ny
  jm = j - 1
  DO i = 1, nx
    DO k = 1,nz
      rhobar = 0.5d0*(rho(k,i,j) + rho(k,i,jm))
      mfv(k,i,j) = rhobar*areas(k,i,j)*v(k,i,j)
    ENDDO
  ENDDO
ENDDO
! Zero volume flux across bottom and top boundaries
mfw(1,:,:) = 0.0d0
mfw(nzp,:,:) = 0.0d0
DO j = 1, ny
  DO i = 1, nx
    DO k = 2,nz
      km = k - 1
      rhobar = above2(k)*rho(k,i,j) + below2(k)*rho(km,i,j)
      mfw(k,i,j) = rhobar*areab(k,i,j)*w(k,i,j)
    ENDDO
  ENDDO
ENDDO

! Now compute divergence
DO j = 1, ny
  jp = j + 1
  DO i = 1, nx
    ip = i + 1
    IF (i == nx) ip = 1
    DO k = 1, nz
      kp = k + 1
      div(k,i,j) = ( (mfu(k,ip,j) - mfu(k,i,j))  &
                   + (mfv(k,i,jp) - mfv(k,i,j))  &
		   + (mfw(kp,i,j) - mfw(k,i,j)) )  &
		   / volume(k,i,j)
    ENDDO
  ENDDO
ENDDO


END SUBROUTINE massdiv

! ========================================================

SUBROUTINE massdiv2(u,v,w,div)

! Compute the 3D divergence of reference density times the velocity
! increment, where w approx drdeta*etadot is assumed.
!
! This version calculates 1/r^2 d/dr (r^2 rho_ref w) so as to mimic
! a finite difference rho_ref 1/r^2 d/dr (r^2 w) + w d/dr rho_ref

USE constants
USE grid
USE timestep
USE refstate

IMPLICIT NONE
REAL*8, INTENT(IN) :: u(nz,nx,ny), v(nz,nx,nyp), w(nzp,nx,ny)
REAL*8, INTENT(OUT) :: div(nz,nx,ny)
INTEGER :: i, im, ip, j, jm, jp, k, km, kp
REAL*8 :: rhobar, h123, h12
REAL*8 :: mfu(nz,nx,ny), mfv(nz,nx,nyp), wbar(nz,nx,ny)


! First compute mass fluxes
DO j = 1, ny
  DO i = 1, nx
    im = i - 1
    IF (i == 1) im = nx
    DO k = 1,nz
      rhobar = 0.5d0*(rho_ref(k,i,j) + rho_ref(k,im,j))
      mfu(k,i,j) = rhobar*areaw(k,i,j)*u(k,i,j)
    ENDDO
  ENDDO
ENDDO
! Zero northward mass flux at poles
mfv(:,:,1) = 0.0d0
mfv(:,:,nyp) = 0.0d0
DO j = 2, ny
  jm = j - 1
  DO i = 1, nx
    DO k = 1,nz
      rhobar = 0.5d0*(rho_ref(k,i,j) + rho_ref(k,i,jm))
      mfv(k,i,j) = rhobar*areas(k,i,j)*v(k,i,j)
    ENDDO
  ENDDO
ENDDO
DO j = 1, ny
  DO i = 1, nx
    DO k = 1, nz
      kp = k + 1
      wbar(k,i,j) = above1(k)*w(kp,i,j) + below1(k)*w(k,i,j)
    ENDDO
  ENDDO
ENDDO

! Now compute divergence
DO j = 1, ny
  jp = j + 1
  DO i = 1, nx
    ip = i + 1
    IF (i == nx) ip = 1
    DO k = 1, nz
      kp = k + 1
      div(k,i,j) = alpha_rho*(                             &
                               (mfu(k,ip,j) - mfu(k,i,j))  &
                             + (mfv(k,i,jp) - mfv(k,i,j))  &
                             + rho_ref(k,i,j)*             &
		               (w(kp,i,j)*areab(kp,i,k) - w(k,i,j)*areab(k,i,j)) &
                              )  / volume(k,i,j)           &
                 + alpha_x*wbar(k,i,j)*drhobardr_ref(k,i,j)
    ENDDO
  ENDDO
ENDDO


END SUBROUTINE massdiv2

! ========================================================

SUBROUTINE hydrostaticrho

! Construct rho field to be in hydrostatic balance with theta

USE grid
USE state
USE constants

IMPLICIT NONE

INTEGER :: i, j, k, kp
REAL*8 :: exnera, exnerb


DO i = 1, nx
  DO j = 1, ny
    exnera = (0.5d0*gascon*(theta(nz,i,j) + theta(nzp,i,j))*rho(nz,i,j)/p00)**kbyonemk
    DO k = nz - 1, 1, -1
      kp = k + 1
      exnerb = exnera + (phi(kp,i,j) - phi(k,i,j))/(cp*theta(kp,i,j))
      rho(k,i,j) = (p00/(0.5d0*gascon*(theta(k,i,j) + theta(kp,i,j))))*(exnerb**onemkbyk)
      exnera = exnerb
    ENDDO
  ENDDO
ENDDO


END SUBROUTINE hydrostaticrho

! ========================================================

SUBROUTINE build_helm

! Build coefficients used in Helmholtz problem

USE grid
USE constants
USE timestep
USE refstate
USE switches
USE work

IMPLICIT NONE

INTEGER :: i, j, k, im, jm, km, kp
REAL*8 :: c1, c2, c3, bbmin, dvmax, dvmin

REAL*8, ALLOCATABLE :: rhobar(:,:,:), temp2d(:,:)


! Allocate memory for temporary arrays
ALLOCATE(rhobar(nzp,nx,ny), temp2d(nx,ny))


! Useful constants
c1 = alpha_rho*alpha_u*dt*dt*cp
c2 = alpha_rho*alpha_w*dt*dt*cp
c3 = alpha_x*alpha_w*dt*dt*cp


! Coefficient of undifferentiated term
ccell = onemkbyk*rho_ref/exner_ref

! Theta in cells
DO k = 1, nz
  kp = k + 1
  thetabar_ref(k,:,:) = above1(k)*theta_ref(kp,:,:) + below1(k)*theta_ref(k,:,:)
ENDDO

! c1 * rho * theta * area / r at western cell boundaries
DO i = 1, nx
  im = i - 1
  IF (i == 1 ) im = nx
  DO j = 1, ny
    DO k = 1, nz
      rtw(k,i,j) = 0.5d0*c1*areaw(k,i,j) &
                           *(rho_ref(k,im,j) + rho_ref(k,i,j)) &
                           *(thetabar_ref(k,im,j) + thetabar_ref(k,i,j)) &
			   /(rp(k,im,j) + rp(k,i,j))
    ENDDO
  ENDDO
ENDDO


! c1 * rho * theta at southern cell boundaries
DO i = 1, nx
  DO j = 2, ny
    jm = j - 1
    DO k = 1, nz
      rts(k,i,j) = 0.5d0*c1*areas(k,i,j) &
                           *(rho_ref(k,i,jm) + rho_ref(k,i,j)) &
			   *(thetabar_ref(k,i,jm) + thetabar_ref(k,i,j)) &
                           /(rp(k,i,jm) + rp(k,i,j))
    ENDDO
  ENDDO
ENDDO
! Zero polar values (SP values are used in a dummy calculation in mgsolve)
rts(:,:,1) = 0.0d0
rts(:,:,nyp) = 0.0d0


! Vertical derivative of rho consistent with semi-Lagrangian advection
CALL refdbydr


! Vertical derivatives of exner at top and bottom of cells,
! bound dtheta/dr away from zero. *** No, just bound bb_ref ***
! Also rho at top and bottom of cells, and buoyancy term
! and drhobar/dr
! Values at bottom of domain
dexnerdr_ref(1,:,:) = 0.0d0
rhobar(1,:,:) = 0.0d0
bb_ref(1,:,:) = 0.0d0
DO k = 2, nz
  km = k - 1
  dexnerdr_ref(k,:,:) = (exner_ref(k,:,:) - exner_ref(km,:,:)) / (rp(k,:,:) - rp(km,:,:))
  rhobar(k,:,:) = above2(k)*rho_ref(k,:,:) + below2(k)*rho_ref(km,:,:)
  temp2d = c3*dthetabardr_ref(k,:,:)*dexnerdr_ref(k,:,:)
  dv3d(k,:,:) = MAX(dv3d(k,:,:),1.0d0 + temp2d)
  bb_ref(k,:,:) = dv3d(k,:,:) - temp2d
ENDDO

bbmin = MINVAL(bb_ref(2:nz,:,:))
dvmin = MINVAL(dv3d)
dvmax = MAXVAL(dv3d)
print *,'Minimum value of delta_v + alpha^2 dt^2 N^2 is ',bbmin
print *,'Maximum value of delta_v = ',dvmax,' Minimum value = ',dvmin
write(48,*) bbmin, dvmin, dvmax

bb_ref(2:nz,:,:) = 1.0d0/(bb_ref(2:nz,:,:))


! Values at top of domain
dexnerdr_ref(nzp,:,:) = 0.0d0
rhobar(nzp,:,:) = 0.0d0
bb_ref(nzp,:,:) = 0.0d0


! Coefficients for cells above and below cell (k,i,j)
! The boundary values set in the loop above ensure that
! cbelow(1,:,:) = 0 and cabove(nz,:,:) = 0, which ensures the
! correct boundary conditions
IF (slice) THEN
  print *,'Check discretization of C() for SLICE '
  STOP
!  DO k = 1, nz
!    kp = k + 1
!    temp2d = rho_ref(k,:,:)/thetabar_ref(k,:,:)
!    cabove(k,:,:) =  bb_ref(kp,:,:)*theta_ref(kp,:,:)*( &
!                     (c2*rrefw(kp)*rrefw(kp)*rhobar(kp,:,:)) &
!                     / (rrefp(k)*rrefp(k)*dzrefp(k)) &
!	             + &
!	             above1(k)*c3*temp2d*dthetabardr_ref(kp,:,:) )
!    cbelow(k,:,:) =  bb_ref(k ,:,:)*theta_ref(k,:,:)*( &
!                     (c2*rrefw(k )*rrefw(k )*rhobar(k ,:,:)) &
!                     / (rrefp(k)*rrefp(k)*dzrefp(k)) &
!	             - &
!		     below1(k)*c3*temp2d*dthetabardr_ref(k ,:,:) )
!  ENDDO
ELSE
  DO k = 1, nz
    kp = k + 1
    temp2d = rho_ref(k,:,:)/thetabar_ref(k,:,:)
    cabove(k,:,:) =  bb_ref(kp,:,:)*theta_ref(kp,:,:)*( &
                     (c2*areab(kp,:,:)*rho_ref(k,:,:)) &
                     / volume(k,:,:) &
	             + &
	             above1(k)*c3*(drhobardr_ref(k,:,:) + &
		                   temp2d*dthetabardr_ref(kp,:,:)) )
    cbelow(k,:,:) =  bb_ref(k ,:,:)*theta_ref(k,:,:)*( &
                     (c2*areab(k ,:,:)*rho_ref(k,:,:)) &
                     / volume(k,:,:) &
	             - &
		     below1(k)*c3*(drhobardr_ref(k,:,:) + &
		                   temp2d*dthetabardr_ref(k ,:,:)) )
  ENDDO
ENDIF



! Tidy temporary memory
DEALLOCATE(rhobar, temp2d)

END SUBROUTINE build_helm

! ===============================================================

SUBROUTINE refdbydr

! For consistency with semi-Lagrangian advection, the reference
! profiles of drho/dr and dtheta/dr should be evaluated using
! cubic Lagrange interpolation. This is done here.
! The vertical derivative is estimated as the average of the derivatives
! given by cubic fit above the level in question and cubic
! fit below the level in question. Linear interpolation is assumed near
! the bottom and top boundaries. This should be consistent with what
! is done in routines lagrangew and lagrangep; and changes there
! should be reflected here.

! NOTE
! dthetabardr_ref is now computed in computed in lagrangewd'

USE grid
USE refstate

IMPLICIT NONE

INTEGER :: kmm, km, k, kp, kpp, i, j
REAL*8 :: dd0, dd1, dd3, dd4, d01, d02, d03, d12, d13, d14, &
          d23, d24, d34, face0, face1, face2, face3, face4



! First, derivative of reference theta

!DO j = 1, ny
!
!DO i = 1, nx
!
!DO k = 2, nz
!
!  kmm = k - 2
!  km = k - 1
!  kp = k + 1
!  kpp = k + 2
!
!  IF (k == 2) THEN
!    ! Linear fit in interval below, cubic above
!    kmm = 1 ! to avoid accessing array index 0
!    dd0 = rw(k,i,j) - rw(kmm,i,j)
!    dd1 = rw(k,i,j) - rw(km,i,j)
!    dd3 = rw(k,i,j) - rw(kp,i,j)
!    dd4 = rw(k,i,j) - rw(kpp,i,j)
!    d01 = rw(kmm,i,j) - rw(km,i,j)
!    d02 = -dd0
!    d03 = rw(kmm,i,j) - rw(kp,i,j)
!    d12 = -dd1
!    d13 = rw(km,i,j) - rw(kp,i,j)
!    d14 = rw(km,i,j) - rw(kpp,i,j)
!    d23 = dd3
!    d24 = dd4
!    d34 = rw(kp,i,j) - rw(kpp,i,j)
!    face0 =  0.0d0
!    face1 =  0.5d0*(  1.0d0/d12 + dd3*dd4/(d12*d13*d14))
!    face2 =  0.5d0*( -1.0d0/d12 &
!                     -(dd1*dd3 + dd1*dd4 + dd3*dd4)/(d12*d23*d24) )
!    face3 =  0.5d0*dd1*dd4/(d13*d23*d34)
!    face4 = -0.5d0*dd1*dd3/(d14*d24*d34)
!  ELSEIF (k == nz) THEN
!    ! Linear fit in the interval above, cubic below
!    kpp = nzp ! to avoid accessing array index nz+2
!    dd0 = rw(k,i,j) - rw(kmm,i,j)
!    dd1 = rw(k,i,j) - rw(km,i,j)
!    dd3 = rw(k,i,j) - rw(kp,i,j)
!    dd4 = rw(k,i,j) - rw(kpp,i,j)
!    d01 = rw(kmm,i,j) - rw(km,i,j)
!    d02 = -dd0
!    d03 = rw(kmm,i,j) - rw(kp,i,j)
!    d12 = -dd1
!    d13 = rw(km,i,j) - rw(kp,i,j)
!    d14 = rw(km,i,j) - rw(kpp,i,j)
!    d23 = dd3
!    d24 = dd4
!    d34 = rw(kp,i,j) - rw(kpp,i,j)
!    face0 =  0.5d0*dd1*dd3/(d01*d02*d03)
!    face1 =  0.5d0*(-dd0*dd3/(d01*d12*d13))
!    face2 =  0.5d0*(  (dd0*dd1 + dd0*dd3 + dd1*dd3)/(d02*d12*d23) &
!                     +1.0d0/d23 )
!    face3 =  0.5d0*(-dd0*dd1/(d03*d13*d23) - 1.0d0/d23)
!    face4 =  0.0d0
!  ELSE
!    ! Cubic fit above and below
!    dd0 = rw(k,i,j) - rw(kmm,i,j)
!    dd1 = rw(k,i,j) - rw(km,i,j)
!    dd3 = rw(k,i,j) - rw(kp,i,j)
!    dd4 = rw(k,i,j) - rw(kpp,i,j)
!    d01 = rw(kmm,i,j) - rw(km,i,j)
!    d02 = -dd0
!    d03 = rw(kmm,i,j) - rw(kp,i,j)
!    d12 = -dd1
!    d13 = rw(km,i,j) - rw(kp,i,j)
!    d14 = rw(km,i,j) - rw(kpp,i,j)
!    d23 = dd3
!    d24 = dd4
!    d34 = rw(kp,i,j) - rw(kpp,i,j)
!    face0 =  0.5d0*dd1*dd3/(d01*d02*d03)
!    face1 =  0.5d0*(-dd0*dd3/(d01*d12*d13) + dd3*dd4/(d12*d13*d14))
!    face2 =  0.5d0*(  (dd0*dd1 + dd0*dd3 + dd1*dd3)/(d02*d12*d23) &
!                     -(dd1*dd3 + dd1*dd4 + dd3*dd4)/(d12*d23*d24) )
!    face3 =  0.5d0*(-dd0*dd1/(d03*d13*d23) + dd1*dd4/(d13*d23*d34))
!    face4 = -0.5d0*dd1*dd3/(d14*d24*d34)
!  ENDIF
!
!  dthetabardr_ref(k,i,j) = face0*theta_ref(kmm,i,j) &
!                         + face1*theta_ref(km ,i,j) &
!                         + face2*theta_ref(k  ,i,j) &
!                         + face3*theta_ref(kp ,i,j) &
!		         + face4*theta_ref(kpp,i,j)
!
!ENDDO

!ENDDO

!ENDDO

! Boundary values are not used
!dthetabardr_ref(1  ,:,:) = 0.0d0
!dthetabardr_ref(nzp,:,:) = 0.0d0



! Now, derivative of reference rho
DO j = 1, ny

DO i = 1, nx

DO k = 1, nz

  kmm = k - 2
  km = k - 1
  kp = k + 1
  kpp = k + 2

  IF (k == 1) THEN
    ! Linear fit using levels 1 and 2
    kmm = 1
    km = 1
    dd3 = rp(k,i,j) - rp(kp,i,j)
    face0 = 0.0d0
    face1 = 0.0d0
    face2 =  1.0d0/dd3
    face3 = -1.0d0/dd3
    face4 = 0.0d0
  ELSEIF (k == nz) THEN
    ! Linear fit using levels nz-1 and nz
    kpp = nz
    kp = nz
    dd1 = rp(k,i,j) - rp(km,i,j)
    face0 = 0.0d0
    face1 = -1.0d0/dd1
    face2 =  1.0d0/dd1
    face3 = 0.0d0
    face4 = 0.0d0
  ELSEIF (k == 2) THEN
    ! Linear fit in interval below, cubic above
    kmm = 1 ! to avoid accessing array index 0
    dd0 = rp(k,i,j) - rp(kmm,i,j)
    dd1 = rp(k,i,j) - rp(km,i,j)
    dd3 = rp(k,i,j) - rp(kp,i,j)
    dd4 = rp(k,i,j) - rp(kpp,i,j)
    d01 = rp(kmm,i,j) - rp(km,i,j)
    d02 = -dd0
    d03 = rp(kmm,i,j) - rp(kp,i,j)
    d12 = -dd1
    d13 = rp(km,i,j) - rp(kp,i,j)
    d14 = rp(km,i,j) - rp(kpp,i,j)
    d23 = dd3
    d24 = dd4
    d34 = rp(kp,i,j) - rp(kpp,i,j)
    face0 =  0.0d0
    face1 =  0.5d0*(  1.0d0/d12 + dd3*dd4/(d12*d13*d14))
    face2 =  0.5d0*( -1.0d0/d12 &
                     -(dd1*dd3 + dd1*dd4 + dd3*dd4)/(d12*d23*d24) )
    face3 =  0.5d0*dd1*dd4/(d13*d23*d34)
    face4 = -0.5d0*dd1*dd3/(d14*d24*d34)
  ELSEIF (k == nz-1) THEN
    ! Linear fit in the interval above, cubic below
    kpp = nz ! to avoid accessing array index nz+1
    dd0 = rp(k,i,j) - rp(kmm,i,j)
    dd1 = rp(k,i,j) - rp(km,i,j)
    dd3 = rp(k,i,j) - rp(kp,i,j)
    dd4 = rp(k,i,j) - rp(kpp,i,j)
    d01 = rp(kmm,i,j) - rp(km,i,j)
    d02 = -dd0
    d03 = rp(kmm,i,j) - rp(kp,i,j)
    d12 = -dd1
    d13 = rp(km,i,j) - rp(kp,i,j)
    d14 = rp(km,i,j) - rp(kpp,i,j)
    d23 = dd3
    d24 = dd4
    d34 = rp(kp,i,j) - rp(kpp,i,j)
    face0 =  0.5d0*dd1*dd3/(d01*d02*d03)
    face1 =  0.5d0*(-dd0*dd3/(d01*d12*d13))
    face2 =  0.5d0*(  (dd0*dd1 + dd0*dd3 + dd1*dd3)/(d02*d12*d23) &
                     +1.0d0/d23 )
    face3 =  0.5d0*(-dd0*dd1/(d03*d13*d23) - 1.0d0/d23)
    face4 =  0.0d0
  ELSE
    ! Cubic fit above and below
    dd0 = rp(k,i,j) - rp(kmm,i,j)
    dd1 = rp(k,i,j) - rp(km,i,j)
    dd3 = rp(k,i,j) - rp(kp,i,j)
    dd4 = rp(k,i,j) - rp(kpp,i,j)
    d01 = rp(kmm,i,j) - rp(km,i,j)
    d02 = -dd0
    d03 = rp(kmm,i,j) - rp(kp,i,j)
    d12 = -dd1
    d13 = rp(km,i,j) - rp(kp,i,j)
    d14 = rp(km,i,j) - rp(kpp,i,j)
    d23 = dd3
    d24 = dd4
    d34 = rp(kp,i,j) - rp(kpp,i,j)
    face0 =  0.5d0*dd1*dd3/(d01*d02*d03)
    face1 =  0.5d0*(-dd0*dd3/(d01*d12*d13) + dd3*dd4/(d12*d13*d14))
    face2 =  0.5d0*(  (dd0*dd1 + dd0*dd3 + dd1*dd3)/(d02*d12*d23) &
                     -(dd1*dd3 + dd1*dd4 + dd3*dd4)/(d12*d23*d24) )
    face3 =  0.5d0*(-dd0*dd1/(d03*d13*d23) + dd1*dd4/(d13*d23*d34))
    face4 = -0.5d0*dd1*dd3/(d14*d24*d34)
  ENDIF

  drhobardr_ref(k,i,j) = face0*rho_ref(kmm,i,j) &
                       + face1*rho_ref(km ,i,j) &
                       + face2*rho_ref(k  ,i,j) &
                       + face3*rho_ref(kp ,i,j) &
                       + face4*rho_ref(kpp,i,j)

ENDDO

ENDDO

ENDDO


END SUBROUTINE refdbydr

! ===============================================================

SUBROUTINE mgsolve(phi,rr,ng)

! Multigrid solver for elliptic equation
!
! const*Div (rho_ref theta_ref grad(phi) )
! + D_1( D_2( phi ) )
! - c0 phi
! = RHS
!
! using full multigrid algorithm with conditional semi-coarsening
! and vertical line solve.
!
! Coefficients are contained in module refstate.


USE grid
USE constants
USE refstate

IMPLICIT NONE

! Numbers of iterations on coarsest grid and other grids
INTEGER, PARAMETER :: niterc = 10, niter = 2, npass = 1


! Conditional coarsening parameter: coarsen in both directions only if
! coupling ratio > alphac
REAL*8, PARAMETER :: alphac = 1.0

INTEGER, INTENT(IN) :: ng
REAL*8, INTENT(IN) :: rr(nz,nx,ny)
REAL*8, INTENT(OUT) :: phi(nz,nx,ny)

INTEGER :: nnx(ng), nny(ng), igrid, igridm, jgrid, jgridm, iter, &
           nnxf, nnyf, nnxc, nnyc, nnxj, nnyj, ipass, sum, ptr(ng+1), &
	   i, j, k, ip, jm, jp, km, kp, ic, jc, ix, &
	   ixij, ixipj, ixijp, ixipjp
INTEGER :: imerge(ny,ng), ixmerge(ny,ng), ixsplit(ny,ng)
LOGICAL :: lsplit(ny,ng), lcns

REAL*8, ALLOCATABLE :: ff(:,:), rf(:,:), &
                       hhw(:,:), hhe(:,:), &
                       hhs(:,:), hhn(:,:), &
                       hha(:,:), hhb(:,:), hhc(:,:), hhv(:,:) 

REAL*8 :: temp1(nz,nx*ny), ddx, dra, drb, temp2(nz,nx,ny)

REAL*8 :: clatc(ny,ng), clats(ny,ng), clatn(ny,ng), &
          clens(ny,ng), clenn(ny,ng), clenew(ny,ng), &
	  cdists(ny,ng), cdistn(ny,ng), cdistew(ny,ng), &
	  carea(ny,ng), careainv(ny,ng), ratio


! ------------------------------------------------------------------------

! One pass should be enough. Warn user if npass is set to
! some other value for testing
IF (npass > 1) PRINT *,'mgsolve: npass = ',npass


! Construct information needed for each grid in the hierarchy
! First determine coarsening pattern

! Finest grid
nnx(ng) = nx
nny(ng) = ny
ddx = dx
DO j = 1, ny
  ! Latitude of cell centre
  clatc(j,ng) = yp(j)
  ! Latitude of northern edge
  clats(j,ng) = yv(j)
  ! Latitude of southern edge
  clatn(j,ng) = yv(j+1)
  ! Length of southern edge
  clens(j,ng) = cosv(j)*dx
  ! Length of northern edge
  clenn(j,ng) = cosv(j+1)*dx
  ! Length of lateral edges
  clenew(j,ng) = dy
  ! Distance to southern neighbour centre
  cdists(j,ng) = dy
  ! Distance to northern neighbour centre
  cdistn(j,ng) = dy
  ! Distance to EW neighbour centre
  cdistew(j,ng) = cosp(j)*dx
  ! Cell area
  carea(j,ng) = cosp(j)*dx*dy
  careainv(j,ng) = 1.0d0/carea(j,ng)
ENDDO

! Coarser grids
DO igrid = ng, 2, -1
  igridm = igrid - 1
  ddx = 2.0*ddx
  j = 1
  jc = 0
  DO WHILE (j <= nny(igrid))
    jp = j + 1
    lcns = .false.
    IF (j < nny(igrid)) THEN
      ratio = clenn(j,igrid)/cdistn(j,igrid)
      IF (ratio > alphac) THEN
        lcns = .true.
      ENDIF
    ENDIF
    IF (lcns) THEN
      ! Coarsen both directions
      jc = jc + 1
      jp = j + 1
      clatc(jc,igridm) = 0.5*(clats(j,igrid) + clatn(jp,igrid))
      clats(jc,igridm) = clats(j,igrid)
      clatn(jc,igridm) = clatn(jp,igrid)
      clens(jc,igridm) = 2.0*clens(j,igrid)
      clenn(jc,igridm) = 2.0*clenn(jp,igrid)
      clenew(jc,igridm) = clenew(j,igrid) + clenew(j+1,igrid)
      cdistew(jc,igridm) = cos(clatc(jc,igridm))*ddx
      carea(jc,igridm) = 2.0*(carea(j,igrid) + carea(j+1,igrid))
      careainv(jc,igridm) = 1.0d0/carea(jc,igridm)
      imerge(j,igrid) = 1
      imerge(jp,igrid) = 2
      ixmerge(j,igrid) = jc
      ixmerge(jp,igrid) = jc
      lsplit(jc,igridm) = .true.
      ixsplit(jc,igridm) = j
      j = j + 2
    ELSE
      ! Coarsen only in EW direction
      jc = jc + 1
      clatc(jc,igridm) = clatc(j,igrid)
      clats(jc,igridm) = clats(j,igrid)
      clatn(jc,igridm) = clatn(j,igrid)
      clens(jc,igridm) = 2.0*clens(j,igrid)
      clenn(jc,igridm) = 2.0*clenn(j,igrid)
      clenew(jc,igridm) = clenew(j,igrid)
      cdistew(jc,igridm) = 2.0*cdistew(j,igrid)
      carea(jc,igridm) = 2.0*carea(j,igrid)
      careainv(jc,igridm) = 1.0d0/carea(jc,igridm)
      imerge(j,igrid) = 0
      ixmerge(j,igrid) = jc
      lsplit(jc,igridm) = .false.
      ixsplit(jc,igridm) = j
      j = j + 1
    ENDIF
  ENDDO
  nnx(igridm) = nnx(igrid)/2
  nny(igridm) = jc
  ! Finally, calculate NS distance between cell centres
  DO jc = 1, nny(igridm)
    jm = jc - 1
    IF (jm > 0) THEN
      cdists(jc,igridm) = clatc(jc,igridm) - clatc(jm,igridm)
    ELSE
      ! Dummy value
      cdists(jc,igridm) = 1.0
    ENDIF
    jp = jc + 1
    IF (jc < nny(igridm)) THEN
      cdistn(jc,igridm) = clatc(jp,igridm) - clatc(jc,igridm)
    ELSE
      ! Dummy value
      cdistn(jc,igridm) = 1.0
    ENDIF
  ENDDO
ENDDO


! Compute pointers and allocate storage for 3D arrays that are needed at all resolutions
sum = 0
DO igrid = 1, ng
  ptr(igrid) = sum + 1
  sum = sum + nnx(igrid)*nny(igrid)
ENDDO
ptr(ng+1) = sum + 1
ALLOCATE(ff(nz,sum), rf(nz,sum))
ALLOCATE(hhw(nz,sum), hhe(nz,sum), hhs(nz,sum), hhn(nz,sum), &
         hha(nz,sum), hhb(nz,sum), hhc(nz,sum), hhv(nz,sum))


! Now compute coefficients needed at all resolutions
! First the finest grid
ix = ptr(ng)
DO j = 1, ny
  jp = j + 1
  DO i = 1, nx
    ip = i + 1
    IF (i == nx) ip = 1
    DO k = 1, nz
      km = k - 1
      kp = k + 1
      IF (k == 1) THEN
        drb = 1.0d0
      ELSE
        drb = rp(k,i,j) - rp(km,i,j)
      ENDIF
      IF (k == nz) THEN
        dra = 1.0d0
      ELSE
        dra = rp(kp,i,j) - rp(k,i,j)
      ENDIF
      hhv(k,ix) = volume(k,i,j)
      hhw(k,ix) = rtw(k,i,j)
      hhe(k,ix) = rtw(k,ip,j)
      hhs(k,ix) = rts(k,i,j)
      hhn(k,ix) = rts(k,i,jp)
      hha(k,ix) = carea(j,ng)*cabove(k,i,j)/dra
      hhb(k,ix) = carea(j,ng)*cbelow(k,i,j)/drb
      hhc(k,ix) = carea(j,ng)*ccell(k,i,j)
    ENDDO
    ix = ix + 1
  ENDDO
ENDDO

! Now successively coarser grids
DO igrid = ng, 2, -1
  igridm = igrid - 1
  ix = ptr(igridm)
  DO jc = 1, nny(igridm)
    j = ixsplit(jc,igridm)
    jp = j + 1
    IF (lsplit(jc,igridm)) THEN
      DO ic = 1, nnx(igridm)
        ip = ic + ic
	i = ip - 1
	ixij = ptr(igrid) - 1 + i + nnx(igrid)*(j-1)
	ixipj = ixij + 1
	ixijp = ixij + nnx(igrid)
	ixipjp = ixijp + 1
        DO k = 1, nz
	  hhv(k,ix) = hhv(k,ixij) + hhv(k,ixijp) + hhv(k,ixijp) + hhv(k,ixipjp)
          hhw(k,ix) = hhw(k,ixij) + hhw(k,ixijp)
          hhe(k,ix) = hhe(k,ixij) + hhe(k,ixijp)
	  hhs(k,ix) = hhs(k,ixij) + hhs(k,ixipj)
	  hhn(k,ix) = hhn(k,ixij) + hhn(k,ixipj)
	  hha(k,ix) = hha(k,ixij) + hha(k,ixipj) + hha(k,ixijp) + hha(k,ixipjp)
 	  hhb(k,ix) = hhb(k,ixij) + hhb(k,ixipj) + hhb(k,ixijp) + hhb(k,ixipjp)
	  hhc(k,ix) = hhc(k,ixij) + hhc(k,ixipj) + hhc(k,ixijp) + hhc(k,ixipjp)
        ENDDO
	ix = ix + 1
      ENDDO
    ELSE
      DO ic = 1, nnx(igridm)
        ip = ic + ic
	i = ip - 1
	ixij = ptr(igrid) - 1 + i + nnx(igrid)*(j-1)
	ixipj = ixij + 1
        DO k = 1, nz
	  hhv(k,ix) = hhv(k,ixij) + hhv(k,ixipj)
          hhw(k,ix) = hhw(k,ixij)
          hhe(k,ix) = hhe(k,ixij)
	  hhs(k,ix) = hhs(k,ixij) + hhs(k,ixipj)
	  hhn(k,ix) = hhn(k,ixij) + hhn(k,ixipj)
	  hha(k,ix) = hha(k,ixij) + hha(k,ixipj)
 	  hhb(k,ix) = hhb(k,ixij) + hhb(k,ixipj)
	  hhc(k,ix) = hhc(k,ixij) + hhc(k,ixipj)
        ENDDO
	ix = ix + 1
      ENDDO
    ENDIF
  ENDDO
ENDDO

! Finally include distance and area factors on all grids
DO igrid = ng, 1, -1
  igridm = igrid - 1
  ix = ptr(igrid)
  DO j = 1, nny(igrid)
      DO i = 1, nnx(igrid)
        DO k = 1, nz
          hhw(k,ix) = hhw(k,ix)/(hhv(k,ix)*cdistew(j,igrid))
          hhe(k,ix) = hhe(k,ix)/(hhv(k,ix)*cdistew(j,igrid))
	  hhs(k,ix) = hhs(k,ix)/(hhv(k,ix)*cdists(j,igrid))
	  hhn(k,ix) = hhn(k,ix)/(hhv(k,ix)*cdistn(j,igrid))
	  hha(k,ix) = hha(k,ix)*careainv(j,igrid)
 	  hhb(k,ix) = hhb(k,ix)*careainv(j,igrid)
	  hhc(k,ix) = hhc(k,ix)*careainv(j,igrid)
        ENDDO
	ix = ix + 1
      ENDDO
  ENDDO
ENDDO


! ----------------------------------------------------------------------------------


! Initialize solution to zero
phi = 0.0d0


DO ipass = 1, npass

! Initialize rhs as residual using latest estimate
IF (ipass == 1) THEN
  ! No need to do the calculation
  ix = 0
  DO j = 1, ny
    DO i = 1, nx
      rf(:,ptr(ng)+ix) = rr(:,i,j)
      ix = ix + 1
    ENDDO
  ENDDO
ELSE
  ix = ptr(ng)
  CALL residual(phi,rr,rf(1,ix), &
		hhw(1,ix),hhe(1,ix),hhs(1,ix),hhn(1,ix), &
		hha(1,ix),hhb(1,ix),hhc(1,ix),nz,nx,ny)
ENDIF

! Initialize solution to zero
ff = 0.0

! Inject right hand side to each grid in the hierarchy
DO igrid = ng, 2, -1
  igridm = igrid - 1
  nnxf = nnx(igrid)
  nnyf = nny(igrid)
  nnxc = nnx(igridm)
  nnyc = nny(igridm)
  CALL inject(rf(1,ptr(igrid)),rf(1,ptr(igridm)),lsplit(1,igridm),ixsplit(1,igridm), &
              carea(1,igrid),careainv(1,igridm),nx,ny,nz,nnxf,nnyf,nnxc,nnyc)
ENDDO

! Iterate to convergence on coarsest grid
nnxj = nnx(1)
nnyj = nny(1)
ff(1,ptr(1):ptr(2)-1) = 0.0d0
ix = ptr(1)
DO iter = 1, niterc
  CALL relax(ff(1,ix),rf(1,ix), &
             hhw(1,ix),hhe(1,ix),hhs(1,ix),hhn(1,ix), &
	     hha(1,ix),hhb(1,ix),hhc(1,ix),nz,nnxj,nnyj)	     
ENDDO


! Sequence of growing V-cycles
DO igrid = 2, ng

  igridm = igrid - 1
  nnxf = nnx(igrid)
  nnyf = nny(igrid)
  nnxc = nnx(igridm)
  nnyc = nny(igridm)

  ! Accurately prolong solution to grid igrid
  ! and execute one V-cycle starting from grid igrid

  !Accurately prolong
!  CALL prolong(ff(1,ptr(igridm)),ff(1,ptr(igrid)),imerge(1,igrid),ixmerge(1,igrid), &
!               clatc(1,igridm),clatc(1,igrid),nx,ny,nz,nnxf,nnyf,nnxc,nnyc)
  CALL prolong2(ff(1,ptr(igridm)),ff(1,ptr(igrid)),imerge(1,igrid),ixmerge(1,igrid), &
                clatc(1,igridm),clatc(1,igrid),nx,ny,nz,nnxf,nnyf,nnxc,nnyc)

  ! Descending part of V-cycle
  DO jgrid = igrid, 2, -1
    
    jgridm = jgrid - 1
    nnxf = nnx(jgrid)
    nnyf = nny(jgrid)
    nnxc = nnx(jgridm)
    nnyc = nny(jgridm)
    ix = ptr(jgrid)

    ! Relax on grid jgrid
    DO iter = 1, niter
      CALL relax(ff(1,ix),rf(1,ix), &
                 hhw(1,ix),hhe(1,ix),hhs(1,ix),hhn(1,ix), &
		 hha(1,ix),hhb(1,ix),hhc(1,ix),nz,nnxf,nnyf)
    ENDDO

    ! Calculate residual on jgrid
    temp1 = 0.0d0
    CALL residual(ff(1,ix),rf(1,ix),temp1, &
		  hhw(1,ix),hhe(1,ix),hhs(1,ix),hhn(1,ix), &
		  hha(1,ix),hhb(1,ix),hhc(1,ix),nz,nnxf,nnyf)
  
    ! Inject residual to jgrid-1
    CALL inject(temp1,rf(1,ptr(jgridm)),lsplit(1,jgridm),ixsplit(1,jgridm), &
                carea(1,jgrid),careainv(1,jgridm),nx,ny,nz,nnxf,nnyf,nnxc,nnyc)

    ! Set correction first guess to zero on grid jgrid-1
    ff(1:nz,ptr(jgridm):ptr(jgrid)-1) = 0.0d0

  ENDDO

  ! Relax to convergence on grid 1
  nnxj = nnx(1)
  nnyj = nny(1)
  ix = ptr(1)
  DO iter = 1, niterc
    CALL relax(ff(1,ix),rf(1,ix), &
               hhw(1,ix),hhe(1,ix),hhs(1,ix),hhn(1,ix), &
	       hha(1,ix),hhb(1,ix),hhc(1,ix),nz,nnxj,nnyj)
  ENDDO


  ! Ascending part of V-cycle
  DO jgrid = 2, igrid

    jgridm = jgrid - 1
    nnxf = nnx(jgrid)
    nnyf = nny(jgrid)
    nnxc = nnx(jgridm)
    nnyc = nny(jgridm)
    ix = ptr(jgrid)

    ! Prolong correction to grid jgrid
    CALL prolong(ff(1,ptr(jgridm)),temp1,imerge(1,jgrid),ixmerge(1,jgrid), &
                 clatc(1,jgridm),clatc(1,jgrid),nx,ny,nz,nnxf,nnyf,nnxc,nnyc)

    ! Add correction to solution on jgrid
    ff(1:nz,ptr(jgrid):ptr(jgrid+1)-1) = ff(1:nz,ptr(jgrid):ptr(jgrid+1)-1) &
                                       + temp1(1:nz,1:ptr(jgrid+1)-ptr(jgrid))

    ! Relax on grid jgrid
    DO iter = 1, niter
      CALL relax(ff(1,ix),rf(1,ix), &
                 hhw(1,ix),hhe(1,ix),hhs(1,ix),hhn(1,ix), &
		 hha(1,ix),hhb(1,ix),hhc(1,ix),nz,nnxf,nnyf)
    ENDDO

  ENDDO

ENDDO

! Add correction to phi
ix = 0
DO j = 1, ny
  DO i = 1, nx
    phi(:,i,j) = phi(:,i,j) + ff(:,ptr(ng)+ix)
    ix = ix + 1
  ENDDO
ENDDO



  ! Just for diagnostics...
  ix = ptr(ng)
  CALL residual(phi,rr,temp2, &
		hhw(1,ix),hhe(1,ix),hhs(1,ix),hhn(1,ix), &
		hha(1,ix),hhb(1,ix),hhc(1,ix),nz,nx,ny)
  call dumprms(temp2,'RMSresidual',nz,nx,ny)
  ! call dumpm(temp2(:,:,ny),'finalresnp',nz,nx)

ENDDO


DEALLOCATE(ff,rf)


END SUBROUTINE mgsolve

! ==========================================================

SUBROUTINE inject(ff,cf,lsplit,ixsplit,caf,cacinv,nx,ny,nz,nnxf,nnyf,nnxc,nnyc)

! Inject data from a fine grid to a coarser grid
! using full area weighting for phi

IMPLICIT NONE

INTEGER, INTENT(IN) :: nx, ny, nz, nnxf, nnyf, nnxc, nnyc, ixsplit(ny)
LOGICAL, INTENT(IN) :: lsplit(ny)
REAL*8, INTENT(IN) :: ff(nz,nnxf,nnyf),caf(ny),cacinv(ny)
REAL*8, INTENT(OUT) :: cf(nz,nnxc,nnyc)
INTEGER :: ic, i, ip, jc, j, jp, k


DO jc = 1, nnyc

  IF (lsplit(jc)) THEN
    ! Merge cells in both latitude and longitude
    j = ixsplit(jc)
    jp = j + 1
    DO ic = 1, nnxc
      ip = ic + ic
      i = ip - 1
      ! Weight by area
      DO k = 1, nz
        cf(k,ic,jc) = ((ff(k,i,j)  + ff(k,ip,j) )*caf(j)  &
	         +     (ff(k,i,jp) + ff(k,ip,jp))*caf(jp) &
                                                         )*cacinv(jc)
      ENDDO
    ENDDO
  ELSE
    ! Merge cells only in longitude
    j = ixsplit(jc)
    DO ic = 1, nnxc
      ip = ic + ic
      i = ip - 1
      ! Weight by area
      DO k = 1, nz
        cf(k,ic,jc) = ((ff(k,i,j) + ff(k,ip,j))*caf(j)  &
	                                              )*cacinv(jc)
      ENDDO
    ENDDO
  ENDIF

ENDDO

END SUBROUTINE inject

! ==========================================================

SUBROUTINE prolong(cf,ff,imerge,ixmerge,clatc,clatf,nx,ny,nz,nnxf,nnyf,nnxc,nnyc)

! Prolong phi field from a coarse grid to a fine grid:
! Cheap version using linear fitting

IMPLICIT NONE

INTEGER, INTENT(IN) :: nx, ny, nz, nnxf, nnyf, nnxc, nnyc, imerge(ny), ixmerge(ny)
REAL*8, INTENT(IN) :: cf(nz,nnxc,nnyc), clatc(ny), clatf(ny)
REAL*8, INTENT(OUT) :: ff(nz,nnxf,nnyf)
REAL*8 :: ws, wn
INTEGER :: i, im, ic, icm, j, jc, jcm, jcp, k


DO j = 1, nnyf

  jc = ixmerge(j)
  IF (imerge(j) == 0  .OR. j == 1 .OR. j == nnyf) THEN
    ! This row of cells was not merged in latitude
    ! or it is next to one pole
    DO ic = 1, nnxc
      icm = MODULO(ic-2,nnxc) + 1
      i = ic + ic - 1
      im = MODULO(i-2,nnxf) + 1
      DO k = 1, nz
        ff(k,im,j) = 0.75*cf(k,icm,jc) + 0.25*cf(k,ic,jc)
        ff(k,i,j)  = 0.25*cf(k,icm,jc) + 0.75*cf(k,ic,jc)
      ENDDO
    ENDDO
  ELSEIF (imerge(j) == 1) THEN
    ! This row is the southern row of a merged pair
    ! Interpolate using southern neighbour
    jcm = jc - 1
    ws = (clatc(jc) - clatf(j))/(clatc(jc) - clatc(jcm))
    wn = 1.0 - ws
    DO ic = 1, nnxc
      icm = MODULO(ic-2,nnxc) + 1
      i = ic + ic - 1
      im = MODULO(i-2,nnxf) + 1
      DO k = 1, nz
        ff(k,im,j) = ws*(0.75*cf(k,icm,jcm) + 0.25*cf(k,ic,jcm)) &
                   + wn*(0.75*cf(k,icm,jc)  + 0.25*cf(k,ic,jc))
        ff(k,i,j)  = ws*(0.25*cf(k,icm,jcm) + 0.75*cf(k,ic,jcm)) &
                   + wn*(0.25*cf(k,icm,jc)  + 0.75*cf(k,ic,jc))
      ENDDO
    ENDDO
  ELSE
    ! This row is the northern row of a merged pair
    ! Interpolate using northern neighbour
    jcp = jc + 1
    wn = (clatf(j) - clatc(jc))/(clatc(jcp) - clatc(jc))
    ws = 1.0 - wn
    DO ic = 1, nnxc
      icm = MODULO(ic-2,nnxc) + 1
      i = ic + ic - 1
      im = MODULO(i-2,nnxf) + 1
      DO k = 1,  nz
        ff(k,im,j) = ws*(0.75*cf(k,icm,jc)  + 0.25*cf(k,ic,jc)) &
                   + wn*(0.75*cf(k,icm,jcp) + 0.25*cf(k,ic,jcp))
        ff(k,i,j)  = ws*(0.25*cf(k,icm,jc)  + 0.75*cf(k,ic,jc)) &
                   + wn*(0.25*cf(k,icm,jcp) + 0.75*cf(k,ic,jcp))
      ENDDO
    ENDDO
  ENDIF

ENDDO


END SUBROUTINE prolong

! ==========================================================

SUBROUTINE prolong2(cf,ff,imerge,ixmerge,clatc,clatf,nx,ny,nz,nnxf,nnyf,nnxc,nnyc)

! Prolong phi field from a coarse grid to a fine grid:
! More accurate version using quadratic fitting

IMPLICIT NONE

INTEGER, INTENT(IN) :: nx, ny, nz, nnxf, nnyf, nnxc, nnyc, imerge(ny), ixmerge(ny)
REAL*8, INTENT(IN) :: cf(nz,nnxc,nnyc), clatc(ny), clatf(ny)
REAL*8, INTENT(OUT) :: ff(nz,nnxf,nnyf)
REAL*8 :: ws, wc, wn, w1, w2, w3
INTEGER :: i, ip, ic, icm, icp, j, jc, jcm, jcp, k


! Weights for EW interpolation
w1 = 5.0D0/32.0D0
w2 = 15.0D0/16.0D0
w3 = -3.0/32.0D0

DO j = 1, nnyf

  jc = ixmerge(j)
  IF (imerge(j) == 0) THEN
    ! This row of cells was not merged in latitude
    DO ic = 1, nnxc
      icm = MODULO(ic-2,nnxc) + 1
      icp = MODULO(ic,  nnxc) + 1
      i = ic + ic - 1
      ip = MODULO(i,nnxf) + 1
      DO k = 1, nz
        ff(k,i ,j) = w1*cf(k,icm,jc) + w2*cf(k,ic,jc) + w3*cf(k,icp,jc)
        ff(k,ip,j) = w3*cf(k,icm,jc) + w2*cf(k,ic,jc) + w1*cf(k,icp,jc)
      ENDDO
    ENDDO
  ELSE
    ! This row is one of a merged pair
    ! Interpolate using southern and northern neighbours
    jcm = jc - 1
    jcp = jc + 1
    if (jcm < 1 ) then
      print *,'Interp past SP '
      STOP
    endif
    ! Coefficients for NS interpolation
    ws = ((clatc(jcp) - clatf(j  ))*(clatc(jc ) - clatf(j  ))) / &
         ((clatc(jcp) - clatc(jcm))*(clatc(jc ) - clatc(jcm)))
    wc = ((clatc(jcp) - clatf(j  ))*(clatc(jcm) - clatf(j  ))) / &
         ((clatc(jcp) - clatc(jc ))*(clatc(jcm) - clatc(jc )))
    wn = ((clatc(jc ) - clatf(j  ))*(clatc(jcm) - clatf(j  ))) / &
         ((clatc(jc ) - clatc(jcp))*(clatc(jcm) - clatc(jcp)))
    DO ic = 1, nnxc
      icm = MODULO(ic-2,nnxc) + 1
      icp = MODULO(ic,  nnxc) + 1
      i = ic + ic - 1
      ip = MODULO(i,nnxf) + 1
      DO k = 1, nz
        ff(k,i ,j) = ws*(w1*cf(k,icm,jcm) + w2*cf(k,ic,jcm) + w3*cf(k,icp,jcm)) &
                   + wc*(w1*cf(k,icm,jc ) + w2*cf(k,ic,jc ) + w3*cf(k,icp,jc )) &
                   + wn*(w1*cf(k,icm,jcp) + w2*cf(k,ic,jcp) + w3*cf(k,icp,jcp))
        ff(k,ip,j) = ws*(w3*cf(k,icm,jcm) + w2*cf(k,ic,jcm) + w1*cf(k,icp,jcm)) &
                   + wc*(w3*cf(k,icm,jc ) + w2*cf(k,ic,jc ) + w1*cf(k,icp,jc )) &
                   + wn*(w3*cf(k,icm,jcp) + w2*cf(k,ic,jcp) + w1*cf(k,icp,jcp))
      ENDDO
    ENDDO
  ENDIF

ENDDO


END SUBROUTINE prolong2

! ==========================================================

SUBROUTINE relax(ff,rf,hhw,hhe,hhs,hhn,hha,hhb,hhc,nz,nnx,nny)

! Red-black relaxation in the horizontal; line solve in the vertical

IMPLICIT NONE

INTEGER, INTENT(IN) :: nz, nnx, nny
REAL*8, INTENT(IN) :: rf(nz,nnx,nny), &
                      hhw(nz,nnx,nny), hhe(nz,nnx,nny), &
		      hhs(nz,nnx,nny), hhn(nz,nnx,nny), &
		      hha(nz,nnx,nny), hhb(nz,nnx,nny), hhc(nz,nnx,nny)
REAL*8, INTENT(INOUT) :: ff(nz,nnx,nny)

INTEGER :: i, j, k, im, jm, ip, jp, iparity, istart

REAL*8 :: sub(nz), sup(nz), diag(nz), rhs(nz)


! Red-black relaxation of phi


DO iparity=1, 2

  ! First row
  j = 1
  jp = j + 1
  istart = MODULO(j + iparity,2) + 1
  DO i = istart, nnx, 2
    im = MODULO(i - 2,nnx) + 1
    ip = MODULO(i,nnx) + 1
    DO k = 1, nz
      diag(k) = -(  hhw(k,i,j) + hhe(k,i,j) + hhn(k,i,j) &
                  + hha(k,i,j) + hhb(k,i,j) + hhc(k,i,j) )
      sub(k) = hhb(k,i,j)
      sup(k) = hha(k,i,j)
      rhs(k) = rf(k,i,j) &
             - (hhw(k,i,j)*ff(k,im,j) + hhe(k,i,j)*ff(k,ip,j) &
	                              + hhn(k,i,j)*ff(k,i,jp))
    ENDDO
    CALL trisolveb(ff(:,i,j),sub,diag,sup,rhs,nz)
  ENDDO

  ! Middle rows
  DO j = 2, nny - 1
    jm = j - 1
    jp = j + 1
    istart = MODULO(j + iparity,2) + 1
    DO i = istart, nnx, 2
      im = MODULO(i - 2,nnx) + 1
      ip = MODULO(i,nnx) + 1
      DO k = 1, nz
        diag(k) = -( hhw(k,i,j) + hhe(k,i,j) + hhs(k,i,j) + hhn(k,i,j) &
	           + hha(k,i,j) + hhb(k,i,j) + hhc(k,i,j) )
	sub(k) = hhb(k,i,j)
	sup(k) = hha(k,i,j)
        rhs(k) = rf(k,i,j) &
               - (hhw(k,i,j)*ff(k,im,j) + hhe(k,i,j)*ff(k,ip,j) &
	        + hhs(k,i,j)*ff(k,i,jm) + hhn(k,i,j)*ff(k,i,jp))
      ENDDO
      CALL trisolveb(ff(:,i,j),sub,diag,sup,rhs,nz)
    ENDDO
  ENDDO
  
  ! Last row
  j = nny
  jm = j - 1
  istart = MODULO(j + iparity,2) + 1
  DO i = istart, nnx, 2
    im = MODULO(i - 2,nnx) + 1
    ip = MODULO(i,nnx) + 1
    DO k = 1, nz
        diag(k) = -( hhw(k,i,j) + hhe(k,i,j) + hhs(k,i,j) &
	           + hha(k,i,j) + hhb(k,i,j) + hhc(k,i,j) )
	sub(k) = hhb(k,i,j)
	sup(k) = hha(k,i,j)
        rhs(k) = rf(k,i,j) &
               - (hhw(k,i,j)*ff(k,im,j) + hhe(k,i,j)*ff(k,ip,j) &
	        + hhs(k,i,j)*ff(k,i,jm)                         )
    ENDDO
    CALL trisolveb(ff(:,i,j),sub,diag,sup,rhs,nz)
  ENDDO

ENDDO


END SUBROUTINE relax

! ==========================================================

SUBROUTINE residual(ff,rf,resf,hhw,hhe,hhs,hhn,hha,hhb,hhc,nz,nnx,nny)

! Calculate residual for elliptic problem

IMPLICIT NONE

INTEGER, INTENT(IN) :: nz, nnx, nny
REAL*8, INTENT(IN) :: ff(nz,nnx,nny), rf(nz,nnx,nny), &
                      hhw(nz,nnx,nny), hhe(nz,nnx,nny), &
		      hhs(nz,nnx,nny), hhn(nz,nnx,nny), &
		      hha(nz,nnx,nny), hhb(nz,nnx,nny), hhc(nz,nnx,nny)

REAL*8, INTENT(OUT) :: resf(nz,nnx,nny)

INTEGER :: i, j, im, jm, ip, jp, k
REAL*8 :: fc, fw, fe, fs, fn, fd, fu


! Residual in phi equation

! First row
j = 1
jp = MODULO(j,nny) + 1
DO i = 1, nnx
  im = MODULO(i - 2,nnx) + 1
  ip = MODULO(i,nnx) + 1
  DO k = 1, nz
    fc = ff(k,i,j)
    fw = ff(k,im,j)
    fe = ff(k,ip,j)
    fn = ff(k,i,jp)
    IF (k == 1) THEN
      fd = 0.0d0
    ELSE
      fd = ff(k-1,i,j)
    ENDIF
    IF (k == nz) THEN
      fu = 0.0d0
    ELSE
      fu = ff(k+1,i,j)
    ENDIF
      resf(k,i,j) = rf(k,i,j) &
                  - ( &
                          hhw(k,i,j)*(fw - fc) &
			+ hhe(k,i,j)*(fe - fc) &
			+ hhn(k,i,j)*(fn - fc) &
		        + hha(k,i,j)*(fu - fc) &
		        + hhb(k,i,j)*(fd - fc) &
		        - hhc(k,i,j)*fc )
  ENDDO
ENDDO

! Middle rows
DO j = 2, nny - 1
  jm = MODULO(j - 2,nny) + 1
  jp = MODULO(j,nny) + 1
  DO i = 1, nnx
    im = MODULO(i - 2,nnx) + 1
    ip = MODULO(i,nnx) + 1
    DO k = 1, nz
      fc = ff(k,i,j)
      fw = ff(k,im,j)
      fe = ff(k,ip,j)
      fs = ff(k,i,jm)
      fn = ff(k,i,jp)
      IF (k == 1) THEN
        fd = 0.0d0
      ELSE
        fd = ff(k-1,i,j)
      ENDIF
      IF (k == nz) THEN
        fu = 0.0d0
      ELSE
        fu = ff(k+1,i,j)
      ENDIF
      resf(k,i,j) = rf(k,i,j) &
                  - ( &
                          hhw(k,i,j)*(fw - fc) &
			+ hhe(k,i,j)*(fe - fc) &
			+ hhs(k,i,j)*(fs - fc) &
			+ hhn(k,i,j)*(fn - fc) &
		        + hha(k,i,j)*(fu - fc) &
		        + hhb(k,i,j)*(fd - fc) &
		        - hhc(k,i,j)*fc )
    ENDDO
  ENDDO
ENDDO

! Last row
j = nny
jm = MODULO(j - 2,nny) + 1
DO i = 1, nnx
  im = MODULO(i - 2,nnx) + 1
  ip = MODULO(i,nnx) + 1
  DO k = 1, nz
    fc = ff(k,i,j)
    fw = ff(k,im,j)
    fe = ff(k,ip,j)
    fs = ff(k,i,jm)
    IF (k == 1) THEN
      fd = 0.0d0
    ELSE
      fd = ff(k-1,i,j)
    ENDIF
    IF (k == nz) THEN
      fu = 0.0d0
    ELSE
      fu = ff(k+1,i,j)
    ENDIF
      resf(k,i,j) = rf(k,i,j) &
                  - ( &
                          hhw(k,i,j)*(fw - fc) &
			+ hhe(k,i,j)*(fe - fc) &
			+ hhs(k,i,j)*(fs - fc) &
		        + hha(k,i,j)*(fu - fc) &
		        + hhb(k,i,j)*(fd - fc) &
		        - hhc(k,i,j)*fc )
  ENDDO
ENDDO


END SUBROUTINE residual

! ==========================================================

SUBROUTINE trisolve(x,a,b,c,r,n)

! To solve the constant coefficient, periodic domain,
! tridiagonal linear system
! Ax = r
! where a is the value below the diagonal of A,
! b is the value on the diagonal of A,
! c is the value above the diagonal of A,
! and r is the known vector right hand side.

IMPLICIT NONE

INTEGER,INTENT(IN) :: n
REAL*8, INTENT(IN) :: a(n), b(n), c(n), r(n)
REAL*8, INTENT(OUT) :: x(n)
INTEGER :: j
REAL*8 :: q(n), s(n), rmx, p


rmx=r(n)

! Forward elimination sweep
q(1) = -c(1)/b(1)
x(1) = r(1)/b(1)
s(1) = -a(1)/b(1)
DO j = 2, n
  p = 1.0/(b(j)+a(j)*q(j-1))
  q(j) = -c(j)*p
  x(j) = (r(j)-a(j)*x(j-1))*p
  s(j) = -a(j)*s(j-1)*p
ENDDO

! Backward pass
q(n) = 0.0
s(n) = 1.0
DO j = n-1, 1, -1
  s(j) = s(j)+q(j)*s(j+1)
  q(j) = x(j)+q(j)*q(j+1)
ENDDO

! Final pass
x(n) = (rmx-c(n)*q(1)-a(n)*q(n-1))/(c(n)*s(1)+a(n)*s(n-1)+b(n))
DO j = 1, n-1
  x(j) = x(n)*s(j)+q(j)
ENDDO


END SUBROUTINE trisolve

! =========================================================

SUBROUTINE trisolveb(x,a,b,c,r,n)

! To solve the constant coefficient, bounded domain,
! tridiagonal linear system
! Ax = r
! where a is the value below the diagonal of A,
! b is the value on the diagonal of A,
! c is the value above the diagonal of A,
! and r is the known vector right hand side.
! a(1) and c(n) must be zero.

IMPLICIT NONE

INTEGER,INTENT(IN) :: n
REAL*8, INTENT(IN) :: a(n), b(n), c(n), r(n)
REAL*8, INTENT(OUT) :: x(n)
INTEGER :: j
REAL*8 :: q(n), p


! Forward elimination sweep
q(1) = -c(1)/b(1)
x(1) = r(1)/b(1)
DO j = 2, n
  p = 1.0/(b(j)+a(j)*q(j-1))
  q(j) = -c(j)*p
  x(j) = (r(j)-a(j)*x(j-1))*p
ENDDO

! Backward pass
DO j = n-1, 1, -1
  x(j) = x(j)+q(j)*x(j+1)
ENDDO


END SUBROUTINE trisolveb

! =========================================================

SUBROUTINE inisbr

! Set up initial conditions: solid body rotation

USE constants
USE state

IMPLICIT NONE

INTEGER :: i, im, j, jm, k, km, kp
REAL*8 :: t_iso, h_rho, a, b, exner(nz),dpidr(nz+1),thetabar(nz), &
          uoa, c1, u00, ptl, dr, cossq, sinr, cosr, u99, rtop, r99


t_iso = 270.0d0
h_rho = gascon*t_iso/gravity
u00 = 40.0d0
uoa = u00/rearth
c1 = uoa*(uoa + 2.0d0*rotatn)
DO i = 1, nx
  DO j = 1, ny
    cossq = COS(geolatp(i,j))**2
    ! Assign top level theta
    rtop = rw(nzp,i,j)
    ptl = 0.5d0*(rtop**2)*cossq*c1/(cp*t_iso)
    theta0(nzp,i,j) = t_iso*EXP(kappa*domain/h_rho)*EXP(-ptl)
    exner(nz) = EXP(-kappa*(rp(nz,i,j)-rgeoid(i,j))/h_rho)
    ptl = 0.5d0*rp(nz,i,j)*rp(nz,i,j)*cossq*c1/(cp*t_iso)
    exner(nz) = exner(nz)*EXP(ptl)
    rho0(nz,i,j) = (p00/(gascon*t_iso))*(exner(nz)**(onemkbyk + 1.0d0))
    ! Now integrate top down to get close to hydrostatic balance
    DO k = nz, 1, -1
      kp = k + 1
      km = k - 1
      a = above1(k)
      b = below1(k)
      thetabar(k) = t_iso/exner(k)
      theta0(k,i,j) = (thetabar(k) - a*theta0(kp,i,j))/b
      IF (k > 1) THEN
        a = above2(k)
        b = below2(k)
        dr = rp(k,i,j) - rp(km,i,j)
	ptl = c1*(a*rp(k,i,j) + b*rp(km,i,j))*cossq
        exner(km) = exner(k) + (-ptl*dr + phi(k,i,j) - phi(km,i,j))/(cp*theta0(k,i,j))
        rho0(km,i,j) = (p00/(gascon*t_iso))*(exner(km)**(onemkbyk + 1.0d0))
      ENDIF
    ENDDO
  ENDDO
ENDDO

sinr = SIN(rotgrid)
cosr = COS(rotgrid)
DO i = 1, nx
  DO j = 1, ny
    DO k = 1, nz
      r99 = 0.5*(rp(k,im,j) + rp(k,i,j))
      u99 = u00*(r99/rearth)
      u0(k,i,j) = (cosr*cosp(j) + sinr*sinp(j)*SIN(xu(i)))*u99
    ENDDO
  ENDDO
ENDDO
DO i = 1, nx
  DO j = 1, nyp
    DO k = 1, nz
      r99 = 0.5*(rp(k,i,jm) + rp(k,i,j))
      u99 = u00*(r99/rearth)
      v0(k,i,j) = sinr*COS(xv(i))*u99
    ENDDO
  ENDDO
ENDDO
DO i = 1, nx
  DO j = 1, ny
    DO k = 1, nzp
      w0(k,i,j) = 0.0d0
    ENDDO
  ENDDO
ENDDO


END SUBROUTINE inisbr

! =====================================================

SUBROUTINE inibw

! Set up initial conditions: Jablonowski and Williamson
! baroclinic instability test case
!
! Note we use `sigma' here in place of the `eta' used in
! J&W's paper to avoid confusion with the ENDGame vertical
! coordinate eta

USE  constants
USE state

IMPLICIT NONE

INTEGER :: i, j, k, kp, iter
REAL*8 :: a, b, exner(nz,ny),dpidr(nz+1),thetabar, &
          t00, gamma, deltatemp, sigmat, u00, sigma0, sigma, sigmav, &
	  tmean, latfac, oneby3, twoby3, eightby5, tenby63, slat, clat, &
	  piby4, one37by60, tenby3, f1, f2, phijw, phimean, temp, &
	  ssv, rcsv, zw, dra


! Handy constants
oneby3 = 1.0d0/3.0d0
twoby3 = 2.0d0/3.0d0
eightby5 = 8.0d0/5.0d0
tenby63 = 10.0d0/63.0d0
piby4 = 0.25d0*pi
one37by60 = 137.0d0/60.0d0
tenby3 = 10.0d0/3.0d0


! Surface temperature
t00 = 288.0d0
! Mean lapse rate
gamma = 0.005d0
! Stratospheric amendment to temperature
deltatemp = 4.8d5
! Tropopause value of sigma
sigmat = 0.2d0
! Maximum zonal wind
u00 = 35.0d0
! Displacement to sigma value
sigma0 = 0.252d0


DO i = 1, nx


! Surface theta
DO j = 1, ny
  sigma = 1.0d0
  slat = sinp(j)
  clat = cosp(j)
  f1 = tenby63 - 2.0d0*(slat**6)*(clat**2 + oneby3)
  f2 = eightby5*(clat**3)*(slat**2 + twoby3) - piby4
  DO iter = 1, 25
    sigmav = (sigma - sigma0)*piby2
    ssv = SIN(sigmav)
    rcsv = SQRT(COS(sigmav))
    tmean = t00*sigma**(gascon*gamma/gravity) &
          + deltatemp*MAX(0.0d0,(sigmat - sigma)**5)
    phimean = (1.0d0 - sigma**(gascon*gamma/gravity))*t00*gravity/gamma
    IF (sigma .lt. sigmat) THEN
      phimean = phimean  &
              - gascon*deltatemp*(   (LOG(sigma/sigmat) + one37by60)*sigmat**5 &
			           - 5.0d0*(sigmat**4)*sigma &
			           + 5.0d0*(sigmat**3)*(sigma**2) &
		                   - tenby3*(sigmat**2)*(sigma**3) &
		                   + 1.25d0*sigmat*(sigma**4) &
			           - 0.2d0*sigma**5 )
    ENDIF
    latfac = f1*u00*rcsv**3 + f2*rearth*rotatn
    phijw = phimean + u00*(rcsv**3)*latfac
    latfac = f1*2.0d0*u00*rcsv**3 + f2*rearth*rotatn
    temp = tmean + (0.75d0*sigma*pi*u00/gascon)*ssv*rcsv*latfac
    zw = rsurf(i,j) - rgeoid(i,j)
    sigma = sigma + (phijw - zw*gravity)/(gascon*temp/sigma)
    sigma = MAX(MIN(sigma,1.1D0),1.0d-8)
  ENDDO
  theta0(1,i,j) = temp*sigma**(-kappa)
ENDDO

DO k = 1, nz
  kp = k + 1
  
  DO j = 1, ny
  
    slat = sinp(j)
    clat = cosp(j)
    f1 = tenby63 - 2.0d0*(slat**6)*(clat**2 + oneby3)
    f2 = eightby5*(clat**3)*(slat**2 + twoby3) - piby4
    
    ! theta levels
    ! First determine sigma 
    sigma = 1.0d-6
    DO iter = 1, 25
      sigmav = (sigma - sigma0)*piby2
      ssv = SIN(sigmav)
      rcsv = SQRT(COS(sigmav))
      tmean = t00*sigma**(gascon*gamma/gravity) &
            + deltatemp*MAX(0.0d0,(sigmat - sigma)**5)
      phimean = (1.0d0 - sigma**(gascon*gamma/gravity))*t00*gravity/gamma
      IF (sigma .lt. sigmat) THEN
        phimean = phimean  &
                - gascon*deltatemp*(   (LOG(sigma/sigmat) + one37by60)*sigmat**5 &
			             - 5.0d0*(sigmat**4)*sigma &
			             + 5.0d0*(sigmat**3)*(sigma**2) &
		        	     - tenby3*(sigmat**2)*(sigma**3) &
		        	     + 1.25d0*sigmat*(sigma**4) &
			             - 0.2d0*sigma**5 )
      ENDIF
      latfac = f1*u00*rcsv**3 + f2*rearth*rotatn
      phijw = phimean + u00*(rcsv**3)*latfac
      latfac = f1*2.0d0*u00*rcsv**3 + f2*rearth*rotatn
      temp = tmean + (0.75d0*sigma*pi*u00/gascon)*ssv*rcsv*latfac
      zw = rw(kp,i,j) - rgeoid(i,j)
      sigma = sigma + (phijw - zw*gravity)/(gascon*temp/sigma)
      sigma = MAX(MIN(sigma,1.1D0),1.0d-8)
    ENDDO
    ! Now assign theta
    theta0(kp,i,j) = temp*sigma**(-kappa)
        
    ! rho levels
    ! First determine sigma 
    sigma = 1.0d-6
    DO iter = 1, 25
      sigmav = (sigma - sigma0)*piby2
      ssv = SIN(sigmav)
      rcsv = SQRT(COS(sigmav))
      tmean = t00*sigma**(gascon*gamma/gravity) &
            + deltatemp*MAX(0.0d0,(sigmat - sigma)**5)
      phimean = (1.0d0 - sigma**(gascon*gamma/gravity))*t00*gravity/gamma
      IF (sigma .lt. sigmat) THEN
        phimean = phimean  &
                - gascon*deltatemp*(   (LOG(sigma/sigmat) + one37by60)*sigmat**5 &
			             - 5.0d0*(sigmat**4)*sigma &
			             + 5.0d0*(sigmat**3)*(sigma**2) &
		        	     - tenby3*(sigmat**2)*(sigma**3) &
		        	     + 1.25d0*sigmat*(sigma**4) &
			             - 0.2d0*sigma**5 )
      ENDIF
      latfac = f1*u00*rcsv**3 + f2*rearth*rotatn
      phijw = phimean + u00*(rcsv**3)*latfac
      latfac = f1*2.0d0*u00*rcsv**3 + f2*rearth*rotatn
      temp = tmean + (0.75d0*sigma*pi*u00/gascon)*ssv*rcsv*latfac
      sigma = sigma + (phijw - phi(k,i,j))/(gascon*temp/sigma)
      sigma = MAX(MIN(sigma,1.1D0),1.0d-8)
    ENDDO
    ! Now assign rho and u
    !! rho0(k,:,j) = p00*sigma/(gascon*temp)
    u0(k,i,j) = u00*(rcsv**3)*(2.0d0*sinp(j)*cosp(j))**2
    ! Determine rho from hydrostatic balance
    IF (k == nz) THEN
      thetabar = above1(k)*theta0(kp,i,j) + below1(k)*theta0(k,i,j)
      exner(k,j) = sigma**kappa
      rho0(k,i,j) = (p00/(gascon*thetabar))*exner(k,j)**onemkbyk
    ENDIF
    
  ENDDO
ENDDO

! Now construct rho to satisfy hydrostatic balance
DO k = nz-1, 1, -1
  kp = k + 1
  DO j = 1, ny
    dra = rp(kp,i,j) - rp(k,i,j)
    exner(k,j) = exner(kp,j)                                                             &
               + ((phi(kp,i,j) - phi(k,i,j))                                             &
               -  ((above2(kp)*u0(kp,i,j) + below2(kp)*u0(k,i,j))**2)*dra/rearth  &
	       -   coriol2(i,j)*(above2(kp)*u0(kp,i,j) + below2(kp)*u0(k,i,j))*dra) &
	       / (cp*theta0(kp,i,j))
    thetabar = above1(k)*theta0(kp,i,j) + below1(k)*theta0(k,i,j)
    rho0(k,i,j) = (p00/(gascon*thetabar))*exner(k,j)**onemkbyk
  ENDDO
ENDDO

ENDDO

! v and w velocity components
v0 = 0.0d0
w0 = 0.0d0




END SUBROUTINE inibw

! =====================================================

SUBROUTINE inibwpert

! Initialize perturbation for baroclinic wave test case.
! This is in a separate subroutine from the initialization
! of the mean state so that balancing iterations can be done
! in between.

USE state

IMPLICIT NONE
INTEGER :: i, j
REAL*8 :: xc = pi/9.0d0, yc = 2.0d0*pi/9.0d0, r0 = 0.1d0, up = 1.0d0
REAL*8 :: cyc, syc, dlon, r, sy, cy



cyc = COS(yc)
syc = SIN(yc)

DO j = 1, ny
  DO  i = 1, nx
    cy = COS(geolatu(i,j))
    sy = SIN(geolatu(i,j))
    dlon = geolonu(i,j) - xc
    r = ACOS(syc*sy + cyc*cy*COS(dlon))
    u0(:,i,j) = u0(:,i,j) + up*EXP(-(r/r0)**2)
  ENDDO
ENDDO

u = u0


END SUBROUTINE inibwpert

! =====================================================

SUBROUTINE iniumjs

! Set up initial conditions: Ullrich et al. baroclinic
! instability test case
!
! Note an inverse square law gravity should be used for this
! initial state: check setphi and findphicol


USE constants
USE state

IMPLICIT NONE

INTEGER :: kk, i, j, k, im, ip, kp
REAL*8 :: b, psfc, te0, tp0, t00, gamma, aa, bb, cc, zz, rbya, zzbh2, &
          e1, e2, f2, tau1, tau2, i1, i2, kk2, tt, rca, rcafac, e3, pp, &
          rcc, exner(nz), uu, thetabar, dra, hscale, uatw, rr, aaa, bbb, &
          dex


! Basic state parameters
b = 2.0d0            ! Half width
kk = 3               ! Exponent for temperature
psfc = 1.0d5         ! Surface pressure
te0 = 310.0d0        ! Surface T at equator
tp0 = 240.0d0        ! Surface T at pole
t00 = 0.5d0*(te0 + tp0)
hscale = gascon*t00/gravity
gamma = 0.005d0      ! Lapse rate

! Derived constants
aa = 1.0d0/gamma
bb = (te0 - tp0)/((te0 + tp0)*tp0)
cc = 0.5d0*(kk  + 2)*(te0 - tp0)/(te0*tp0)
kk2 = kk/(kk + 2.0d0)


! First set theta
DO j = 1, ny
  DO i = 1, nx
    DO k = 1, nzp
      zz = rw(k,i,j) - rearth
      rbya = rw(k,i,j)/rearth
      zzbh2 = (zz/(b*hscale))**2
      e1 = EXP(zz*gamma/t00)
      e2 = EXP(-zzbh2)
      f2 = 1.0d0 - 2.0d0*zzbh2
      tau1 = aa*(gamma/t00)*e1 + bb*f2*e2
      tau2 = cc*f2*e2
      i1 = aa*(e1 - 1.0d0) + bb*zz*e2
      i2 = cc*zz*e2
      rca = rbya*COS(geolatp(i,j))
      rcafac = rca**kk - kk2*rca**(kk+2)
      ! Temperature
      tt = 1.0d0/((tau1 - tau2*rcafac)*rbya*rbya)
      e3 = (i2*rcafac - i1)*gravity/gascon
      pp = psfc*EXP(e3)
      ! theta
      theta0(k,i,j) = tt*(p00/pp)**kappa
    ENDDO
  ENDDO
ENDDO

! Now set zonal wind
! (This needs to be modified for a rotated grid)
DO j = 1, ny
  DO i = 1, nx
    im = i - 1
    IF (i == 1) im = nx
    DO k = 1, nzp
      rr = 0.5d0*(rp(k,im,j) + rp(k,i,j))
      zz = rr - rearth
      rbya = rr/rearth
      zzbh2 = (zz/(b*hscale))**2
      e1 = EXP(zz*gamma/t00)
      e2 = EXP(-zzbh2)
      f2 = 1.0d0 - 2.0d0*zzbh2
      tau1 = aa*(gamma/t00)*e1 + bb*f2*e2
      tau2 = cc*f2*e2
      i1 = aa*(e1 - 1.0d0) + bb*zz*e2
      i2 = cc*zz*e2
      rca = rbya*COS(geolatu(i,j))
      rcafac = rca**kk - kk2*rca**(kk+2)
      ! Temperature
      tt = 1.0d0/((tau1 - tau2*rcafac)*rbya*rbya)
      rcafac = rca**(kk-1) - rca**(kk+1)
      ! `Wind proxy'
      uu = kk*tt*i2*rcafac*gravity/rearth
      ! Zonal wind
      rcc = rr*COS(geolatu(i,j))
      u0(k,i,j) = SQRT((rotatn*rcc)**2 + rcc*uu) - rotatn*rcc
    ENDDO
  ENDDO
ENDDO

! Finally integrate from top to bottom to obtain exner and hence rho that
! give balance
! Coefficients to extrapolate exner to the surface
aaa = -etap(1)/(etap(2) - etap(1))
bbb = 1.0d0 - aaa
DO j = 1, ny
  DO i = 1, nx
    ip = i + 1
    IF (i == nx) ip = 1
    zz = rp(nz,i,j) - rearth
    rbya = rp(nz,i,j)/rearth
    zzbh2 = (zz/(b*hscale))**2
    e1 = EXP(zz*gamma/t00)
    e2 = EXP(-zzbh2)
    f2 = 1.0d0 - 2.0d0*zzbh2
    tau1 = aa*(gamma/t00)*e1 + bb*f2*e2
    tau2 = cc*f2*e2
    i1 = aa*(e1 - 1.0d0) + bb*zz*e2
    i2 = cc*zz*e2
    rca = rbya*COS(geolatp(i,j))
    rcafac = rca**kk - kk2*rca**(kk+2)
    ! Pressure
    e3 = (i2*rcafac - i1)*gravity/gascon
    pp = psfc*EXP(e3)
    exner(nz) = (pp/p00)**kappa
    DO k = nz-1, 1, -1
      kp = k + 1
      dra = rp(kp,i,j) - rp(k,i,j)
      uatw = 0.5d0*( above2(kp)*u0(kp,i ,j) + below2(kp)*u0(k,i ,j)     &
                   + above2(kp)*u0(kp,ip,j) + below2(kp)*u0(k,ip,j) )
      exner(k) = exner(kp)                        &
               + ( (phi(kp,i,j) - phi(k,i,j))     &
               -   (uatw**2)*dra/rw(kp,i,j)       &
	       -   coriol2(i,j)*uatw*dra      )   &
	       / (cp*theta0(kp,i,j))
    ENDDO
    ! Correction to exner to ensure surface pressure equals psfc
    dex = (psfc/p00)**kappa - (aaa*exner(2) + bbb*exner(1))
    exner = exner + dex
    ! Compute density
    DO k = 1, nz
      kp = k + 1
      thetabar = above1(k)*theta0(kp,i,j) + below1(k)*theta0(k,i,j)
      rho0(k,i,j) = (p00/(gascon*thetabar))*exner(k)**onemkbyk
    ENDDO
  ENDDO
ENDDO


! v and w velocity components
v0 = 0.0d0
w0 = 0.0d0


END SUBROUTINE iniumjs

! =====================================================

SUBROUTINE iniumjspert

! Initialize perturbation for Ullrich et al baroclinic wave test case.
! This is in a separate subroutine from the initialization
! of the mean state so that balancing iterations can be done
! in between.

USE state
USE constants

IMPLICIT NONE
INTEGER :: i, im, j, jm, k
REAL*8 :: d0, vp, zt, xc, yc
REAL*8 :: cyc, syc, dlon, r, sy, cy, cdx, d, taper, fac, sdx, zz, mag



! Set parameters
d0 = 1.0d0/6.0d0     ! Radius of perturbation
vp = 1.0d0           ! Wind amplitude of perturbation
zt = 1.5d4           ! Top of perturbation
xc = pi/9.0d0        ! Coordinates of centre of
yc = twopi/9.0d0     ! perturbation

mag = 16.0d0*vp/(3.0d0*SQRT(3.0d0))

cyc = COS(yc)
syc = SIN(yc)

! Set u
DO j = 1, ny
  DO  i = 1, nx
    im = i - 1
    IF (i == 1) im = nx
    cy = COS(geolatu(i,j))
    sy = SIN(geolatu(i,j))
    dlon = geolonu(i,j) - xc
    cdx = COS(dlon)
    d = ACOS(syc*sy + cyc*cy*cdx)
    IF (d > 0.0d0 .AND. d < d0) THEN
      DO k = 1, nz
        zz = 0.5d0*(rp(k,im,j) + rp(k,i,j)) - rearth
        IF (zz < zt) THEN
          taper = 1.0d0 - 3.0d0*(zz/zt)**2 + 2.0d0*(zz/zt)**3
          fac = (-syc*cy + cyc*sy*cdx)/SIN(d)
          u0(k,i,j) = u0(k,i,j) - mag*taper*(COS(piby2*d/d0)**3)*SIN(piby2*d/d0)*fac
        ENDIF
      ENDDO
    ENDIF
  ENDDO
ENDDO

! Set v
DO j = 2, ny
  jm = j - 1
  DO  i = 1, nx
    im = i - 1
    IF (i == 1) im = nx
    cy = COS(geolatv(i,j))
    sy = SIN(geolatv(i,j))
    dlon = geolonv(i,j) - xc
    cdx = COS(dlon)
    sdx = SIN(dlon)
    d = ACOS(syc*sy + cyc*cy*cdx)
    IF (d > 0.0d0 .AND. d < d0) THEN
      DO k = 1, nz
        zz = 0.5d0*(rp(k,i,jm) + rp(k,i,j)) - rearth
        IF (zz < zt) THEN
          taper = 1.0d0 - 3.0d0*(zz/zt)**2 + 2.0d0*(zz/zt)**3
          fac = cyc*sdx/SIN(d)
          v0(k,i,j) = v0(k,i,j) + mag*taper*(COS(piby2*d/d0)**3)*SIN(piby2*d/d0)*fac
        ENDIF
      ENDDO
    ENDIF
  ENDDO
ENDDO


u = u0
v = v0
print *,'max initial v = ',maxval(abs(v0))

END SUBROUTINE iniumjspert

! =====================================================

SUBROUTINE iniacoustic

! Set up initial conditions: internal acoustic wave on a resting
! hydrostatic basic state

USE  constants
USE state

IMPLICIT NONE

INTEGER :: i, j, k, km, kp
REAL*8 :: t_iso, h_rho, a, b, exner(nz),dpidr(nz+1),thetabar(nz)


! Hydrostatic balance
t_iso = 250.0d0
h_rho = gascon*t_iso/gravity
DO i = 1, nx
  DO j = 1, ny
    theta0(1,i,j) = t_iso
    exner(1) = EXP(-kappa*(rp(1,i,j)-rgeoid(i,j))/h_rho)
    rho0(1,i,j) = (p00/(gascon*t_iso))*(exner(1)**(onemkbyk + 1.0d0))
    DO k = 1, nz
      kp = k + 1
      a = above1(k)
      b = below1(k)
      thetabar(k) = t_iso/exner(k)
      theta0(kp,i,j) = (thetabar(k) - b*theta0(k,i,j))/a
      IF (k < nz) THEN
        exner(kp) = exner(k) - (phi(kp,i,j) - phi(k,i,j))/(cp*theta0(kp,i,j))
        rho0(kp,i,j) = (p00/(gascon*t_iso))*(exner(kp)**(onemkbyk + 1.0d0))
        dpidr(kp) = (exner(kp) - exner(k))/(rp(kp,i,j) - rp(k,i,j))
      ENDIF
    ENDDO
  ENDDO
ENDDO


DO i = 1, nx
  DO j = 1, ny
    DO k = 1, nz
      u0(k,i,j) = 0.0d0
    ENDDO
  ENDDO
ENDDO
DO i = 1, nx
  DO j = 1, nyp
    DO k = 1, nz
      v0(k,i,j) = 0.0d0
    ENDDO
  ENDDO
ENDDO
DO i = 1, nx
  DO j = 1, ny
    DO k = 1, nzp
      w0(k,i,j) = 0.0d0
!      w0(k,i,j) = 0.01d0*SIN(pi*etaw(k)/domain)*EXP(etaw(k)/(2.0*h_rho))
    ENDDO
  ENDDO
ENDDO


END SUBROUTINE iniacoustic

! =====================================================

SUBROUTINE inibubble

! Initial state for a warm bubble test case.
! Suggest using a 10km planetary radius

USE constants
USE state

IMPLICIT NONE

INTEGER :: i, j, k, km, kp
REAL*8 :: theta00, rbub, hbub, thetapert, x, y, z, x0, d, thetabar, &
          a, b, exner(nz)


! Isentropic basic state
theta00 = 300.0d0
theta0 = theta00
    
! With a warm bubble centred on the equator
! Bubble radius
rbub = 1000.0d0    
! Height of bubble centre
hbub = 2000.0d0
x0 = rearth + hbub
! Amplitude of potential temperature perturbation
thetapert = 2.0d0

DO j = 1, ny
  DO i = 1, nx
    DO k = 1, nzp
      x = rw(k,i,j)*cosp(j)*COS(xp(i))
      y = rw(k,i,j)*cosp(j)*SIN(xp(i))
      z = rw(k,i,j)*sinp(j)
      d = SQRT((x + x0)**2 + y**2 + z**2)
      d = MIN(d,rbub)
      theta0(k,i,j) = theta0(k,i,j) + thetapert*COS(piby2*d/rbub)**2
    ENDDO
    
    ! Now integrate top down to get close to hydrostatic balance
    exner(nz) = 1.0d0 - gravity*domain/(cp*theta00)
    a = above1(nz)
    b = below1(nz)
    thetabar = a*theta0(nzp,i,j) + b*theta0(nz,i,j)
    rho0(nz,i,j) = (p00/(gascon*thetabar))*exner(nz)**onemkbyk
    DO k = nz-1, 1, -1
      kp = k + 1
      exner(k) = exner(kp) + (phi(kp,i,j) - phi(k,i,j))/(cp*theta0(kp,i,j))
      a = above1(k)
      b = below1(k)
      thetabar = a*theta0(kp,i,j) + b*theta0(k,i,j)
      rho0(k,i,j) = (p00/(gascon*thetabar))*exner(k)**onemkbyk
    ENDDO
  ENDDO
ENDDO

u0 = 0.0d0
v0 = 0.0d0
w0 = 0.0d0


END SUBROUTINE inibubble

! =====================================================

SUBROUTINE initracers

! Allocate space for tracers and initialize

USE state

IMPLICIT NONE

INTEGER :: k

IF (ntracers > 0) THEN

  ALLOCATE(tracers0(nz,nx,ny,ntracers),tracers(nz,nx,ny,ntracers))

  ! First two tracers are used for theta and PV; set in start_diagnostics

  ! Other tracers...

ENDIF


END SUBROUTINE initracers

! =====================================================

SUBROUTINE initrajectories

! Allocate space for trajectories and initialize

USE trajectories
USE grid, ONLY : nz, pi, etap

IMPLICIT NONE

INTEGER, PARAMETER :: ntlat = 24, ntlong = 96
INTEGER :: i, j, k, itraj
REAL*8 :: dx, dy, y, x

IF (ntrajectories > 0) THEN

  ALLOCATE(xtraj(ntrajectories),ytraj(ntrajectories),etraj(ntrajectories))

  ! Quasi-uniform distribution of parcels in NH, on model p-levels
  dy = 0.5d0*pi/ntlat
  itraj = 0
  DO j = 1, ntlat
    y = (j - 0.5d0)*dy
    x = 0.0d0
    dx = 2.0d0*pi/(ntlong*COS(y))
    DO WHILE (x < 2.0d0*pi)
      DO k = 1, nz
        itraj = itraj + 1
        xtraj(itraj) = x
        ytraj(itraj) = y
        etraj(itraj) = etap(k)
      ENDDO
      x = x + dx
    ENDDO
  ENDDO

  IF (itraj .ne. ntrajectories) THEN
    PRINT *,'Number of trajectories initialized ',itraj, &
            ' different from ntrajectories ',ntrajectories
    STOP
  ENDIF

  IF (nlabels > 0) THEN
    ALLOCATE(labels(ntrajectories,nlabels))
  ENDIF

ENDIF


END SUBROUTINE initrajectories

! =====================================================

SUBROUTINE checkbal

! Check whether initial data are reasonably well balanced

USE constants
USE switches
USE state
USE work
USE refstate
USE timestep
use departure  ! temporary to check departures

IMPLICIT NONE

INTEGER :: igo, iter, i, jm, j, km, k, kp
LOGICAL :: done
REAL*8 :: maxu, maxv, maxw
REAL*8 :: r1(nzp,ny), r2(nz,nyp), q(nzp,ny), pp(nz,ny), r(nz,ny), &
          tt(nz), dpi, dra


! Define reference profiles
rho_ref = rho0
theta_ref = theta0
CALL findexner(rho_ref,theta_ref,exner_ref)


done = .false.
DO WHILE (.not. done)

  ! Compute RHS terms on the model grid
  CALL rhsgrid

  ! Compute departure points and departure point terms
  CALL up_outer
  
  ! Build Helmholtz coefficients
  ! First reset delta_v to unregularized value
  dv3d = REAL(delta_v)
  CALL build_helm

  ! Compute arrival point terms and hence right hand sides
  CALL up_inner

  maxu = MAXVAL(ABS(rhsu))
  maxv = MAXVAL(ABS(rhsv))
  maxw = MAXVAL(ABS(rhsw))  
  PRINT *,' '
  PRINT *,'Checking balance: '
  PRINT *,'Maximum imbalance'
  PRINT *,'in u  ',maxu
  PRINT *,'in v  ',maxv
  PRINT *,'in w  ',maxw
  PRINT *,' '
  PRINT *,'*********************************************************'
  PRINT *,'* Enter 0 to stop now.                                  *'
  PRINT *,'* Enter 1 to continue with the integration.             *'
  PRINT *,'* Enter 2 to initialize in quasi-hydrostatic balance    *'
  PRINT *,'* Enter 3 to carry out some balancing iterations        *'
  PRINT *,'*     (only suitable for zonally symmetric states       *'
  PRINT *,'*     with v=0 and w=0).                                *'
  PRINT *,'*********************************************************'
  PRINT *,' '
  READ(5,*) igo

  IF (igo == 0) THEN
    CALL epilogue
    STOP
  ENDIF

  IF (igo == 1) done = .true.
  
  IF (igo == 2) THEN
  
    ! Increment exner to get close to quasi-hydrostatic balance
  
    DO i = 1, nx
    
      r1 = rhsw(:,i,:)/(cp*theta0(:,i,:)*dt)   
      DO j = 1, ny
        dpi = 0.0d0
	pp(nz,:) = 0.0d0
        DO k = nz-1, 1, -1
          kp = k + 1
	  dra = rp(kp,i,j) - rp(k,i,j)
          dpi = dpi - r1(kp,j)*dra
	  pp(k,j) = dpi
        ENDDO
      ENDDO
       
      DO k = 1, nz
        kp = k + 1
	r(k,:) = rho0(k,i,:)*(onemkbyk*pp(k,:)/exner_ref(k,i,:))
      ENDDO
      rho(:,i,:) = rho(:,i,:) + r
       
    ENDDO
 
    rho0 = rho
 
  ENDIF

  IF (igo == 3) THEN

      ! Increment exner and theta to get closer to balance
      
      r1 = rhsw(:,1,:)/(cp*theta0(:,1,:)*dt)
      DO j = 2, ny
        jm = j - 1
	tt = 0.5d0*(thetabar_ref(:,1,jm) + thetabar_ref(:,1,j))
        r2(:,j) = rhsv(:,1,j)*rearth/(cp*tt*dt)
      ENDDO
      r2(:,1) = 0.0d0
      r2(:,nyp) = 0.0d0

      ! First chose Pi increment to satisfy v equation
      DO k = 1, nz
        pp(k,1) = 0.0d0
        DO j = 2, ny
          jm = j - 1
          pp(k,j) = pp(k,jm) + r2(k,j)*dy
        ENDDO
        pp(k,:) = pp(k,:) - SUM(pp(k,:))/ny
      ENDDO

      ! Latitudinally independent increment to Pi to correct
      ! w equation as far as possible
      dpi = 0.0d0
      DO k = nz-1, 1, -1
        kp = k + 1
	dpi = dpi - SUM(r1(kp,:)*(rp(kp,1,:) - rp(k,1,:)))/ny
	pp(k,:) = pp(k,:) + dpi
      ENDDO

      ! Now choose theta increment to satisfy w equation
      q(1,:) = 0.0d0
      DO k = 2, nz
        km = k - 1
        DO j = 1, ny
	  q(k,j) = r1(k,j) - (pp(k,j) - pp(km,j))/(rp(k,1,j) - rp(km,1,j))
	ENDDO
      ENDDO
      q(nzp,:) = 0.0d0
      q = q/MIN(dexnerdr_ref(:,1,:),-1.0d-12)
      q = q*theta0(:,1,:)
      
      DO k = 1, nz
        kp = k + 1
	r(k,:) = rho0(k,1,:)* &
	         (onemkbyk*pp(k,:)/exner_ref(k,1,:) &
                - (above1(k)*q(kp,:) + below1(k)*q(k,:))/thetabar_ref(k,1,:))
      ENDDO
       
      rho(:,1,:) = rho(:,1,:) + r
      theta(:,1,:) = theta(:,1,:) + q
      
      DO i = 2, nx
        rho(:,i,:) = rho(:,1,:)
	theta(:,i,:) = theta(:,1,:)
      ENDDO	
      
      ! Finally copy back into initial fields
      rho0 = rho
      theta0 = theta
      u0 = u
      v0 = v
      w0 = w

  ENDIF

ENDDO


END SUBROUTINE checkbal

! ============================================================

SUBROUTINE restartdump

! Write out all information needed to restart this integration

USE runtype
USE timestep
USE departure
USE state
USE channels

IMPLICIT NONE
CHARACTER*25 :: yfile

! Note: A more sophisticated version would write out all grid
! information too.  For now we trust the user not to mess this up.

! Construct file name
WRITE(yfile,100) runid, istep
100 FORMAT('run',A5,'restart',I10.10)

OPEN(chanresout,FILE=yfile,FORM='UNFORMATTED')
WRITE(chanresout) istep
WRITE(chanresout) xdepu
WRITE(chanresout) ydepu
WRITE(chanresout) edepu
WRITE(chanresout) xdepv
WRITE(chanresout) ydepv
WRITE(chanresout) edepv
WRITE(chanresout) xdepw
WRITE(chanresout) ydepw
WRITE(chanresout) edepw
WRITE(chanresout) xdepp
WRITE(chanresout) ydepp
WRITE(chanresout) edepp
WRITE(chanresout) rho
WRITE(chanresout) theta
WRITE(chanresout) u
WRITE(chanresout) v
WRITE(chanresout) w

CLOSE(chanresout)

END SUBROUTINE restartdump

! ============================================================

SUBROUTINE restart

! Read in restart data ready to continue an integration

USE runtype
USE timestep
USE departure
USE state
USE channels

IMPLICIT NONE


OPEN(chanresin,FILE=yrestart,FORM='UNFORMATTED')
READ(chanresin) istep
READ(chanresin) xdepu
READ(chanresin) ydepu
READ(chanresin) edepu
READ(chanresin) xdepv
READ(chanresin) ydepv
READ(chanresin) edepv
READ(chanresin) xdepw
READ(chanresin) ydepw
READ(chanresin) edepw
READ(chanresin) xdepp
READ(chanresin) ydepp
READ(chanresin) edepp
READ(chanresin) rho0
READ(chanresin) theta0
READ(chanresin) u0
READ(chanresin) v0
READ(chanresin) w0

! Compute eta coordinate vertical velocity
CALL findetadot(u0,v0,w0,etadot0)

rho = rho0
theta = theta0
u = u0
v = v0
w = w0
etadot = etadot0


END SUBROUTINE restart

! ============================================================

SUBROUTINE readn48grid

! Read grid data from N48 dump.

! xp and yp are read in here for comparison with the values
! subsequently set in setupgrid.

USE grid
USE constants
USE channels
IMPLICIT NONE

INTEGER :: k, i, j, l

! Note single precision!
REAL*4 :: griddata(3*nx*ny*nz)

! p points
OPEN(chann48,FILE="p_grid-000.dat",FORM="UNFORMATTED",CONVERT="BIG_ENDIAN")
READ(chann48) griddata
CLOSE(chann48)
l = 1
DO k = 1, nz
  DO j = 1, ny
    DO i = 1, nx
      xp(i) = griddata(l)
      yp(j) = griddata(l+1)
      rp(k,i,j) = griddata(l+2) + rearth
      l = l + 3
    ENDDO
  ENDDO
  etap(k) = MINVAL(rp(k,:,:)) - rearth
ENDDO

! w points
OPEN(chann48,FILE="w_grid-000.dat",FORM="UNFORMATTED",CONVERT="BIG_ENDIAN")
READ(chann48) griddata
CLOSE(chann48)
l = 1
DO k = 1, nzp
  DO j = 1, ny
    DO i = 1, nx
      xp(i) = griddata(l)
      yp(j) = griddata(l+1)
      rw(k,i,j) = griddata(l+2) + rearth
      l = l + 3
    ENDDO
  ENDDO
  IF (k == nzp) rw(k,:,:) = 80000.0D0 + rearth
  etaw(k) = MINVAL(rw(k,:,:)) - rearth
ENDDO


! Remove orography
DO k = 1, nz
  rp(k,:,:) = etap(k) + rearth
ENDDO
DO k = 1, nzp
  rw(k,:,:) = etaw(k) + rearth
ENDDO
PRINT *,'N48 orography removed'


! Set surface height
rsurf = rw(1,:,:)


END SUBROUTINE readn48grid

! ============================================================

SUBROUTINE readn48data

! Read fields from N48 dump

USE channels
USE state
IMPLICIT NONE

INTEGER :: i, j, k

REAL*8, ALLOCATABLE :: fld(:,:,:)


! Rho
ALLOCATE(fld(nx,ny,nz))
OPEN(chann48,FILE="rho-000000-000.dat",FORM="UNFORMATTED",CONVERT="BIG_ENDIAN")
READ(chann48) fld
CLOSE(chann48)
DO k = 1, nz
  DO j = 1, ny
    DO i = 1, nx
      rho0(k,i,j) = fld(i,j,k)
    ENDDO
  ENDDO
ENDDO
DEALLOCATE(fld)
!call dumpm(rho0(1,:,:),'rholev1',nx,ny)

! Theta
ALLOCATE(fld(nx,ny,nzp))
OPEN(chann48,FILE="theta-000000-000.dat",FORM="UNFORMATTED",CONVERT="BIG_ENDIAN")
READ(chann48) fld
CLOSE(chann48)
DO k = 1, nzp
  DO j = 1, ny
    DO i = 1, nx
      theta0(k,i,j) = fld(i,j,k)
    ENDDO
  ENDDO
ENDDO
DEALLOCATE(fld)

! u
ALLOCATE(fld(nx,ny,nz))
OPEN(chann48,FILE="Uvel-000000-000.dat",FORM="UNFORMATTED",CONVERT="BIG_ENDIAN")
READ(chann48) fld
CLOSE(chann48)
DO k = 1, nz
  DO j = 1, ny
    DO i = 1, nx
      u0(k,i,j) = fld(i,j,k)
    ENDDO
  ENDDO
ENDDO
DEALLOCATE(fld)

! v
ALLOCATE(fld(nx,nyp,nz))
OPEN(chann48,FILE="Vvel-000000-000.dat",FORM="UNFORMATTED",CONVERT="BIG_ENDIAN")
READ(chann48) fld
CLOSE(chann48)
DO k = 1, nz
  DO j = 1, nyp
    DO i = 1, nx
      v0(k,i,j) = fld(i,j,k)
    ENDDO
  ENDDO
ENDDO
DEALLOCATE(fld)

! w
ALLOCATE(fld(nx,ny,nzp))
OPEN(chann48,FILE="Wvel-000000-000.dat",FORM="UNFORMATTED",CONVERT="BIG_ENDIAN")
READ(chann48) fld
CLOSE(chann48)
DO k = 1, nzp
  DO j = 1, ny
    DO i = 1, nx
      w0(k,i,j) = fld(i,j,k)
    ENDDO
  ENDDO
ENDDO
DEALLOCATE(fld)


END SUBROUTINE readn48data

! ============================================================

SUBROUTINE bwout

! Output for baroclinic wave test case

USE constants
USE state

IMPLICIT NONE

INTEGER :: i, j, k, km, kp
REAL*8 :: exner(nz,nx,ny), psurf(nx,ny), t850(nx,ny), t2d(nz,ny)
REAL*8 :: a, b, exner0, exnera, exnerb, thetabar, meanps
LOGICAL :: done


! First compute exner
CALL findexner(rho,theta,exner)

! Extrapolate to surface
a = -etap(1)/(etap(2) - etap(1))
b = 1.0d0 - a
psurf = a*exner(2,:,:) + b*exner(1,:,:)


! Construct T at 850 hPa

exner0 = (85000.0d0/p00)**kappa

DO j = 1, ny
  DO i = 1, nx

    ! Find first theta level with exner < exner0
    k = 2
    km = k - 1
    exnerb = psurf(i,j)
    exnera = above2(k)*exner(k,i,j) + below2(k)*exner(km,i,j)
    done = .false.
    DO WHILE (.not. done)
      IF (exnera < exner0) THEN
        done = .true.
      ELSE
        km = k
        k = k + 1
        exnerb = exnera
        exnera = above2(k)*exner(k,i,j) + below2(k)*exner(km,i,j)
      ENDIF
    ENDDO
    a = (exner0 - exnerb)/(exnera - exnerb)
    b = 1.0d0 - a
    t850(i,j) = (a*theta(k,i,j) + b*theta(km,i,j))*exner0

  ENDDO
ENDDO

! Finally convert surface exner to pressure
! and subtract 1d5 to avoid loss of precision in output
psurf = p00*psurf**(1.0d0/kappa) - 1.0d5

! and output
call dumpm(psurf,'Psurf',nx,ny)
call dumpm(t850,'T850',nx,ny)


t2d = 0.0d0
DO i = 1, nx
  DO j = 1, ny
    DO k = 1, nz
      kp = k + 1
      thetabar = below1(k)*theta(k,i,j) + above1(k)*theta(kp,i,j)
      t2d(k,j) = t2d(k,j) + thetabar*exner(k,i,j)
    ENDDO
  ENDDO
ENDDO
t2d = t2d/nx

! call dumpm(t2d,'ZonalmeanT',nz,ny)

! Diagnose global mean surface pressure
meanps = 0.0d0
DO j = 1, ny
  DO i = 1, nx
    meanps = meanps + psurf(i,j)*area(j)
  ENDDO
ENDDO
meanps = meanps/areatot
print *,'Mean surface pressure = 1e5 + ',meanps




END SUBROUTINE bwout

! ============================================================

SUBROUTINE locatepmin(pmin,lonmin,latmin)

! Determine the coordinates of the minimum surface pressure

USE constants
USE state

IMPLICIT NONE

REAL *8, INTENT(OUT) :: pmin, lonmin, latmin
INTEGER :: i, j, ix, jx
REAL*8 :: exner(nz,nx,ny), psurf(nx,ny)
REAL*8 :: a, b


! First compute exner
CALL findexner(rho,theta,exner)

! Extrapolate to surface
a = -etap(1)/(etap(2) - etap(1))
b = 1.0d0 - a
psurf = a*exner(2,:,:) + b*exner(1,:,:)

! Convert surface exner to pressure
psurf = p00*psurf**(1.0d0/kappa)

! Find and dump location of minimum surface pressure
ix = 0
jx = 0
pmin = 1.0d6
DO j = 1, ny
  DO i = 1, nx
    IF (psurf(i,j) .lt. pmin) THEN
      ix = i
      jx = j
      pmin = psurf(i,j)
    ENDIF
  ENDDO
ENDDO
lonmin = geolonp(ix,jx)
latmin = geolatp(ix,jx)

END SUBROUTINE locatepmin

! ============================================================

SUBROUTINE locatemax(q,qmax,ix,iy,iz,nx,ny,nz)

! Determine the value and indices of the max abs of field q(nz,nx,ny)

IMPLICIT NONE

INTEGER, INTENT(IN) :: nx, ny, nz
REAL*8, INTENT(IN) :: q(nz,nx,ny)
INTEGER, INTENT(OUT) :: ix, iy, iz
REAL*8, INTENT(OUT) :: qmax
INTEGER :: i, j, k
REAL*8 :: aq


! Find and location of max abs
qmax = 0.0d0
ix = 0
iy = 0
iz = 0
DO j = 1, ny
  DO i = 1, nx
    DO k = 1, nz
      aq = ABS(q(k,i,j))
      IF (aq .gt. qmax) THEN
        ix = i
        iy = j
        iz = k
        qmax = aq
      ENDIF
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE locatemax

! ============================================================

SUBROUTINE potvort

! Compute potential vorticity

USE state
USE lagdiagnostics
!
IMPLICIT NONE
!
INTEGER :: j, k
REAL*8 :: gradthetar(nz,nx,nyp), gradthetal(nz,nx,nyp), gradthetap(nz,nx,nyp), &
          zetal(nz,nx,nyp), zetap(nz,nx,nyp), zetar(nz,nx,nyp), &
          drvbydr(nz,nx,nyp), drucbydr(nz,nx,nyp), dwbydr(nz,nx,nyp), &
          drdlam(nz,nx,nyp), drdphi(nz,nx,nyp), &
          rp_pv(nz,nx,nyp), rpvcosv(nz,nx,nyp), cortemp(nx,nyp), &
          pvmin, pvmax

REAL*8, ALLOCATABLE :: ratu(:,:,:), ratv(:,:,:), &
                       raboveu(:,:,:), rabovev(:,:,:)
REAL*8, ALLOCATABLE :: tempp(:,:,:), tempu(:,:,:), tempv(:,:,:), &
                       tempw(:,:,:), tempt(:,:,:), tempu2D(:,:)
!
! Allocate arrays
ALLOCATE(tempp(nz,nx,ny), tempu(nz,nx,ny), tempv(nz,nx,nyp), &
         tempw(nzp,nx,ny), tempt(nzp,nx,nyp), tempu2D(nx,nyp))

!
!
!  *****  PV_R TERM  *****
!
! Arrays needed only here
ALLOCATE(ratu(nz,nx,ny), ratv(nz,nx,nyp), &
         raboveu(nzp,nx,ny), rabovev(nzp,nx,nyp))

! Calculate gradthetar, and then average to PV points
CALL dwbydrw(nx,ny,rw,theta,tempp)                  ! dtheta/dr at p points
CALL averp2u(nz,ny,tempp,tempu)
CALL averp2v(nz,nx,tempu,gradthetar)                ! dtheta/dr at pv points

! Calculate and average zetar to PV points 
CALL averp2v(nz,nx,rp,ratv)                         ! r at v points
tempv = ratv*v                                      ! r*v at v points
CALL dpbydlam(nz,nyp,tempv,zetar)                   ! d(rv)/dlam at pv points
CALL averp2w(nx,nyp,tempv,tempt)                    ! r*v above v points
CALL averp2v(nzp,nx,rw,rabovev)                     ! r above v points
CALL dwbydrw(nx,nyp,rabovev,tempt,tempv)            ! d(rv)/dr at v points
CALL averp2u(nz,nyp,tempv,drvbydr)                  ! d(rv)/dr at pv points
CALL dpbydlam(nz,nyp,ratv,drdlam)                   ! dr/dlam at pv points
zetar = zetar - drdlam*drvbydr                      ! total v contrib to r*r*cos*zetar

CALL averp2u(nz,ny,rp,ratu)                         ! r at u points
DO j = 1, ny
   tempu(:,:,j) = ratu(:,:,j)*u(:,:,j)*cosp(j)      ! r*u*cosp at u points
END DO
CALL dpbydphi(nz,nx,tempu,tempv)                    ! d(rucosp)/dphi at pv points
zetar = zetar - tempv
CALL averp2w(nx,ny,tempu,tempw)                     ! r*u*cosp above u points
CALL averp2u(nzp,ny,rw,raboveu)                     ! r above u points
CALL dwbydrw(nx,ny,raboveu,tempw,tempu)             ! d(rucos)/dr at u points
CALL averp2v(nz,nx,tempu,drucbydr)                  ! d(rucos)/dr at pv points
CALL dpbydphi(nz,nx,ratu,drdphi)                    ! dr/dphi at pv points
zetar = zetar + drdphi*drucbydr                     ! total u and v contrib to r*r*cos*zetar

! Update zetar with metric (r*r*cosv)
CALL averp2v(nz,nx,ratu,rp_pv)                      ! r at pv points
DO j = 2, ny
   rpvcosv(:,:,j) = rp_pv(:,:,j)*cosv(j)            ! r*cos at pv points
   zetar(:,:,j) = zetar(:,:,j) &
                /( rp_pv(:,:,j)*rpvcosv(:,:,j) )    ! total zetar
END DO
rpvcosv(:,:,1) = rp_pv(:,:,1)*cosv(1)
rpvcosv(:,:,nyp) = rp_pv(:,:,j)*cosv(nyp)

! Update PV with vertical Coriolis component
CALL averp2u2D(ny,coriol3,tempu2D)
CALL averp2v2D(nx,tempu2D,cortemp)                  ! Vertical Coriolis at pv points
DO k = 1, nz
  pv(k,:,:) = gradthetar(k,:,:)*(zetar(k,:,:) + cortemp)
END DO

! Release arrays no longer needed
DEALLOCATE(ratu, ratv, raboveu, rabovev)

!
!
!  *****  PV_LAMBDA TERM  *****
!
! Calculate gradtheta lambda, and then average to PV points
CALL dpbydlam(nzp,ny,theta,tempw)                   ! dtheta/dlam above u points
CALL averw2p(nx,ny,tempw,tempu)                     ! dtheta/dlam at u points
CALL averp2v(nz,nx,tempu,gradthetal)                ! dtheta/dlam at pv points
gradthetal = gradthetal - drdlam*gradthetar         ! dtheta/dlam with bent term at pv points
DO j = 2, ny
   gradthetal(:,:,j) = gradthetal(:,:,j)/rpvcosv(:,:,j) ! zonal grad theta at pv points
END DO

! Calculate and average zeta lambda to PV points 
CALL dpbydphi(nzp,nx,w,tempt)                       ! dw/dphi above v points
CALL averw2p(nx,nyp,tempt,tempv)                    ! dw/dphi at v points
CALL averp2u(nz,nyp,tempv,zetal)                    ! dw/dphi at pv points
CALL dwbydrw(nx,ny,rw,w,tempp)                      ! dw/dr at p points
CALL averp2u(nz,ny,tempp,tempu)                     ! dw/dr at u points
CALL averp2v(nz,nx,tempu,dwbydr)                    ! dw/dr at pv points
zetal = zetal - drvbydr - drdphi*dwbydr             
zetal = zetal/rp_pv                                 ! total zetal

! Update PV with zonal Coriolis component
CALL averp2u2D(ny,coriol1,tempu2D)
CALL averp2v2D(nx,tempu2D,cortemp)                  ! zonal Coriolis at pv points
DO k = 1, nz
  pv(k,:,:) = pv(k,:,:) + gradthetal(k,:,:)*(zetal(k,:,:) + cortemp)
END DO

!
!
!  *****  PV_PHI TERM  *****
!
! Calculate gradtheta phi, and then average to PV points
CALL dpbydphi(nzp,nx,theta,tempt)                   ! dtheta/dphi above v points
CALL averw2p(nx,nyp,tempt,tempv)                    ! dtheta/dphi at v points
CALL averp2u(nz,nyp,tempv,gradthetap)               ! dtheta/dphi at pv points
gradthetap = gradthetap - drdphi*gradthetar         ! dtheta/dphi with bent term at pv points
gradthetap = gradthetap/rp_pv                       ! meridional grad theta at pv points

! Calculate and average zeta phi to PV points 
CALL dpbydlam(nzp,ny,w,tempw)                       ! dw/dlam above u points
CALL averw2p(nx,ny,tempw,tempu)                     ! dw/dlam at u points
CALL averp2v(nz,nx,tempu,zetap)                     ! dw/dlam at pv points
zetap =  drucbydr - zetap - drdlam*dwbydr
DO j = 2, ny
   zetap(:,:,j) = zetap(:,:,j)/rpvcosv(:,:,j)       ! total zetap
END DO

! Update PV with meridional Coriolis component
CALL averp2u2D(ny,coriol2,tempu2D)
CALL averp2v2D(nx,tempu2D,cortemp)                  ! meridional Coriolis at pv points
DO k = 1, nz
  pv(k,:,:) = pv(k,:,:) + gradthetap(k,:,:)*(zetap(k,:,:) + cortemp)
END DO

!
!
!  *****  RHO TERM  *****
!
! Calculate average of rho in PV points and put it in tempv array
CALL averp2u(nz,ny,rho,tempu)
CALL averp2v(nz,nx,tempu,tempv)

! Update PV, all except polar values
DO j = 2, ny
   pv(:,:,j) = pv(:,:,j)/tempv(:,:,j)
END DO

! Polar values of PV
DO k = 1, nz
  pv(k,:,1) = SUM(pv(k,:,2))/nx 
  pv(k,:,nyp) = SUM(pv(k,:,ny))/nx
ENDDO

!
! Maxima and minima
pvmax = MAXVAL(pv)
pvmin = MINVAL(pv)

!print *,'max zetar,x,y = ',maxval(zetar),maxval(zetal),maxval(zetap)
!print *,'max zetar,x,y = ',maxval(zetar(:,:,2:ny)),maxval(zetal(:,:,2:ny)),maxval(zetap(:,:,2:ny))
!print *,'max pv = ',maxval(pv)
!
! Dealocate arrays
DEALLOCATE(tempp, tempu, tempv, tempw, tempt, tempu2D)
!
RETURN
END SUBROUTINE potvort

! ==========================================================
!!!
SUBROUTINE averp2w(nxa,nya,q,averq_pw)
!
! IK: Compute the VERTICAL average of a scalar field Q.
! Q is stored at P (U,V) points; the values of the average 
! are stored at W (above U, above V) points . The values are 
! extrapolated to the bottom and top boundaries
!
USE grid, ONLY: nz, nzp, above2, below2, above1, below1
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nxa, nya
REAL*8, INTENT(IN) :: q(nz,nxa,nya)
REAL*8, INTENT(OUT) :: averq_pw(nzp,nxa,nya)
INTEGER :: k, km
!
! Vertical average of Q
! Q_bar at W (above U, above V) points
DO k = 2, nz
  km = k - 1
  averq_pw(k,:,:) = above2(k)*q(k,:,:) + below2(k)*q(km,:,:)
ENDDO
! ! ! ! Constant extrapolation to bottom and top boundaries
! ! ! averq_pw(1,:,:) = q(1,:,:)
! ! ! averq_pw(nzp,:,:) = q(nz,:,:)
! Linear extrapolation to bottom and top boundaries
averq_pw(1,:,:) = (q(1,:,:) - above1(1)*averq_pw(2,:,:))/below1(1)
averq_pw(nzp,:,:) = (q(nz,:,:) - below1(nz)*averq_pw(nz,:,:))/above1(nz)
!
RETURN
END SUBROUTINE averp2w
!!!
! ==========================================================
!!!
SUBROUTINE averw2p(nxa,nya,q,averq_wp)
!
! IK: Compute the VERTICAL average of a scalar field Q.
! Q is stored at W (above U, above V) points; the values of 
! the average are stored at P (U, V) points.
!
USE grid, ONLY: nz, nzp, above1, below1
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nxa, nya
REAL*8, INTENT(IN) :: q(nzp,nxa,nya)
REAL*8, INTENT(OUT) :: averq_wp(nz,nxa,nya)
INTEGER :: k, kp
!
! Vertical average of Q
! Q_bar at P (U, V) points
DO k = 1, nz
  kp = k + 1
  averq_wp(k,:,:) = above1(k)*q(kp,:,:) + below1(k)*q(k,:,:)
ENDDO
!
RETURN
END SUBROUTINE averw2p
!!!
! ==========================================================
!!!
SUBROUTINE averp2u(nza,nya,q,averq_pu)
!
! IK: Compute the HORIZONTAL average in LAMBDA direction of a 
! scalar field Q. Q is stored at P (W) points; the values of the 
! average are stored at U (above U) points.
!
USE grid, ONLY: nx
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nza, nya
REAL*8, INTENT(IN) :: q(nza,nx,nya)
REAL*8, INTENT(OUT) :: averq_pu(nza,nx,nya)
INTEGER :: i, im
!
! Horizontal average of Q
! Q_bar at U (above U) points
DO i = 1, nx
  im = i - 1
  IF (i == 1) im = nx
  averq_pu(:,i,:) = 0.5d0*(q(:,im,:) + q(:,i,:))
ENDDO
!
RETURN
END SUBROUTINE averp2u
!!!
! ==========================================================
!!!
SUBROUTINE averu2p(nza,nya,q,averq_up)
!
! IK: Compute the HORIZONTAL average in LAMBDA direction of a 
! scalar field Q. Q is stored at U (above U) points; the values of the 
! average are stored at P (W) points.
!
USE grid, ONLY: nx
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nza, nya
REAL*8, INTENT(IN) :: q(nza,nx,nya)
REAL*8, INTENT(OUT) :: averq_up(nza,nx,nya)
INTEGER :: i, ip
!
! Horizontal average of Q
! Q_bar at P (W) points
DO i = 1, nx
  ip = i + 1
  IF (i == nx) ip = 1
  averq_up(:,i,:) = 0.5d0*(q(:,i,:) + q(:,ip,:))
ENDDO
!
RETURN
END SUBROUTINE averu2p
!!!
! ==========================================================
!!!
SUBROUTINE averp2v(nza,nxa,q,averq_pv)
!
! IK: Compute the HORIZONTAL average in PHI direction of a 
! scalar field Q. Q is stored at P (W) points; the values of the 
! average are stored at V (above V) points.
!
USE grid, ONLY: ny, nyp, nx
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nza, nxa
REAL*8, INTENT(IN) :: q(nza,nxa,ny)
REAL*8, INTENT(OUT) :: averq_pv(nza,nxa,nyp)
INTEGER :: j, jm, k
!
! Horizontal average of Q
! Q_bar at V (above V) points
DO j = 2, ny
  jm = j - 1
  averq_pv(:,:,j) = 0.5d0*(q(:,:,jm) + q(:,:,j))
ENDDO
! Update polar values, can be modified outside the routine
DO k = 1, nza
  averq_pv(k,:,1) = SUM(averq_pv(k,:,2))/nx
  averq_pv(k,:,nyp) = SUM(averq_pv(k,:,ny))/nx
ENDDO
!
RETURN
END SUBROUTINE averp2v
!!!
! ==========================================================
!!!
SUBROUTINE averv2p(nza,nxa,q,averq_vp)
!
! IK: Compute the HORIZONTAL average in PHI direction of a 
! scalar field Q. Q is stored at V (above V) points; the values of the 
! average are stored at P (W) points.
!
USE grid, ONLY: ny, nyp
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nza, nxa
REAL*8, INTENT(IN) :: q(nza,nxa,nyp)
REAL*8, INTENT(OUT) :: averq_vp(nza,nxa,ny)
INTEGER :: j, jp
!
! Horizontal average of Q
! Q_bar at P (W) points
DO j = 1, ny
  jp = j + 1
  averq_vp(:,:,j) = 0.5d0*(q(:,:,j) + q(:,:,jp))
ENDDO
!
RETURN
END SUBROUTINE averv2p
!!!
! ==========================================================
!!!
SUBROUTINE dwbydrw(nxa,nya,rwq,q,dqbydrw)
!
! IK: Compute the VERTICAL derivative of a scalar field Q.
! Q is stored at W (above U, above V) points; the values of 
! the derivative are stored at P (U, V) points.
!
USE grid, ONLY: nz, nzp
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nxa, nya
REAL*8, INTENT(IN) :: rwq(nzp,nxa,nya), q(nzp,nxa,nya)
REAL*8, INTENT(OUT) :: dqbydrw(nz,nxa,nya)
REAL*8 :: dr
INTEGER :: i, j, k, kp
!
! Vertical derivative of Q
! Q_der at P (U, V) points
DO k = 1, nz
  kp = k + 1
  DO i = 1, nxa
    DO j = 1, nya
      dr = rwq(kp,i,j) - rwq(k,i,j)
      dqbydrw(k,i,j) = (q(kp,i,j) - q(k,i,j))/dr
    END DO
  END DO
END DO
!
RETURN
END SUBROUTINE dwbydrw
!!!
! ==========================================================
!!!
SUBROUTINE dpbydrp(nxa,nya,rpq,q,dqbydrp)
!
! IK: Compute the VERTICAL derivative of a scalar field Q.
! Q is stored at P (U,V) points; the values of the derivative 
! are stored at W (above U, above V) points . The values are 
! extrapolated to the bottom and top boundaries
!
USE grid, ONLY: nz, nzp, above1, below1
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nxa, nya
REAL*8, INTENT(IN) :: rpq(nz,nxa,nya), q(nz,nxa,nya)
REAL*8, INTENT(OUT) :: dqbydrp(nzp,nxa,nya)
REAL*8 :: dr
INTEGER :: i, j, k, km
!
! Vertical derivative of Q
! Q_der at W (above U, above V) points
DO k = 2, nz
  km = k - 1
  DO i = 1, nxa
    DO j = 1, nya
      dr = rpq(k,i,j) - rpq(km,i,j)
      dqbydrp(k,i,j) = (q(k,i,j) - q(km,i,j))/dr
    END DO
  END DO
END DO
! Constant extrapolation to bottom and top boundaries
dqbydrp(1,:,:) = dqbydrp(2,:,:)
dqbydrp(nzp,:,:) = dqbydrp(nz,:,:)
!
RETURN
END SUBROUTINE dpbydrp
!!!
! ==========================================================
!!!
SUBROUTINE dpbydlam(nza,nya,q,dqdlam)
!
! IK: Compute the HORIZONTAL derivative in LAMBDA direction of a 
! scalar field Q. Q is stored at P (W) points; the values of the 
! derivative are stored at U (above U) points.
!
USE grid, ONLY: nx, dx
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nza, nya
REAL*8, INTENT(IN) :: q(nza,nx,nya)
REAL*8, INTENT(OUT) :: dqdlam(nza,nx,nya)
INTEGER :: i, im
!
! Horizontal average of Q
! Q_der at U (above U) points
DO i = 1, nx
  im = i - 1
  IF (i == 1) im = nx
  dqdlam(:,i,:) = (q(:,i,:) - q(:,im,:))/dx
ENDDO
!
RETURN
END SUBROUTINE dpbydlam
!!!
! ==========================================================
!!!
SUBROUTINE dubydlam(nza,nya,q,dqdlam)
!
! IK: Compute the HORIZONTAL derivative in LAMBDA direction of a 
! scalar field Q. Q is stored at U (above U) points; the values of the 
! derivative are stored at P (W) points.
!
USE grid, ONLY: nx, dx
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nza, nya
REAL*8, INTENT(IN) :: q(nza,nx,nya)
REAL*8, INTENT(OUT) :: dqdlam(nza,nx,nya)
INTEGER :: i, ip
!
! Horizontal derivative of Q
! Q_der at P (W) points
DO i = 1, nx
  ip = i + 1
  IF (i == nx) ip = 1
  dqdlam(:,i,:) = (q(:,ip,:) - q(:,i,:))/dx
ENDDO
!
RETURN
END SUBROUTINE dubydlam
!!!
! ==========================================================
!!!
SUBROUTINE dpbydphi(nza,nxa,q,dqdphi)
!
! IK: Compute the HORIZONTAL derivative in PHI direction of a 
! scalar field Q. Q is stored at P (W) points; the values of the 
! derivative are stored at V (above V) points.
!
USE grid, ONLY: ny, nyp, dy
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nza, nxa
REAL*8, INTENT(IN) :: q(nza,nxa,ny)
REAL*8, INTENT(OUT) :: dqdphi(nza,nxa,nyp)
INTEGER :: j, jm
!
! Horizontal average of Q
! Q_der at V (above V) points
DO j = 2, ny
  jm = j - 1
  dqdphi(:,:,j) = (q(:,:,j) - q(:,:,jm))/dy
ENDDO
! Set polar values to zero, can be modified outside the routine
dqdphi(:,:,1) = 0.0d0
dqdphi(:,:,nyp) = 0.0d0
!
RETURN
END SUBROUTINE dpbydphi
!!!
! ==========================================================
!!!
SUBROUTINE dvbydphi(nza,nxa,q,dqdphi)
!
! IK: Compute the HORIZONTAL derivative in PHI direction of a 
! scalar field Q. Q is stored at V (above V) points; the values of the 
! derivative are stored at P (W) points.
!
USE grid, ONLY: ny, nyp, dy
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nza, nxa
REAL*8, INTENT(IN) :: q(nza,nxa,nyp)
REAL*8, INTENT(OUT) :: dqdphi(nza,nxa,ny)
INTEGER :: j, jp
!
! Horizontal average of Q
! Q_der at P (W) points
DO j = 1, ny
  jp = j + 1
  dqdphi(:,:,j) = (q(:,:,jp) - q(:,:,j))/dy
ENDDO
!
RETURN
END SUBROUTINE dvbydphi
!!!
! ==========================================================
!!!
SUBROUTINE averp2v2D(nxa,q,averq_pv)
!
! IK: Compute the HORIZONTAL average in PHI direction of a 
! 2D scalar field Q. Q is stored at P (W) points on one layer; 
! the values of the average are stored at V (above V) points,
! also on one layer.
!
USE grid, ONLY: ny, nyp, nx
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nxa
REAL*8, INTENT(IN) :: q(nxa,ny)
REAL*8, INTENT(OUT) :: averq_pv(nxa,nyp)
INTEGER :: j, jm
!
! Horizontal average of Q
! Q_bar at V (above V) points
DO j = 2, ny
  jm = j - 1
  averq_pv(:,j) = 0.5d0*(q(:,jm) + q(:,j))
ENDDO
! Update polar values, can be modified outside the routine
averq_pv(:,1) = SUM(averq_pv(:,2))/nx
averq_pv(:,nyp) = SUM(averq_pv(:,ny))/nx
!
RETURN
END SUBROUTINE averp2v2D
!!!
! ==========================================================
!!!
SUBROUTINE averp2u2D(nya,q,averq_pu)
!
! IK: Compute the HORIZONTAL average in LAMBDA direction of a 
! 2D scalar field Q. Q is stored at P (W) pointson one layer; 
! the values of the average are stored at U (above U) points,
! also on one layer.
!
USE grid, ONLY: nx
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nya
REAL*8, INTENT(IN) :: q(nx,nya)
REAL*8, INTENT(OUT) :: averq_pu(nx,nya)
INTEGER :: i, im
!
! Horizontal average of Q
! Q_bar at U (above U) points
DO i = 1, nx
  im = i - 1
  IF (i == 1) im = nx
  averq_pu(i,:) = 0.5d0*(q(im,:) + q(i,:))
ENDDO
!
RETURN
END SUBROUTINE averp2u2D
!!!
! ==========================================================

SUBROUTINE totalmass

! To diagnose  total atmospheric mass

USE channels
USE timestep
USE state

IMPLICIT NONE

INTEGER :: i, j, k
REAL*8 :: totalm

totalm = 0.0d0
DO i = 1, nx
  DO j = 1, ny
    DO k = 1, nz
      totalm = totalm + rho(k,i,j)*volume(k,i,j)
    ENDDO
  ENDDO
ENDDO
write(chanmass,'(I10,E15.8)') istep, totalm
print *,'Step ',istep,'  total mass = ',totalm


END SUBROUTINE totalmass

! ============================================================

SUBROUTINE globaldiags

! To diagnose total atmospheric mass, energy, entropy,
! and angular momentum. Energy is decomposed into internal,
! potential and kinetic
! The ke contributions are evaluated at u, v, w
! points.


USE channels
USE timestep
USE state
USE constants
USE switches
USE globdiag

IMPLICIT NONE

INTEGER :: i, j, k, ip, jp, kp, im, jm, km
REAL*8 :: cv, exner(nz,nx,ny), dm, ke, thetabar, urcosp, wr2c2, &
          totalm, internaleng, potentialeng, kineticeng, totaleng, &
          totalent, totalam, sinr, cosr, ra, rb, vv, totvol


cv = cp - gascon

CALL findexner(rho,theta,exner)

totvol = 0.0d0
totalm = 0.0d0
internaleng = 0.0d0
potentialeng = 0.0d0
kineticeng = 0.0d0
totalent = 0.0d0
totalam = 0.0d0

sinr = SIN(rotgrid)
cosr = COS(rotgrid)

DO j = 1, ny
  jp = j + 1
  DO i = 1, nx
    ip = i + 1
    IF (i == nx) ip = 1
    DO k = 1, nz
      kp = k + 1
      totvol = totvol + volume(k,i,j)
      dm = rho(k,i,j)*volume(k,i,j)
      totalm = totalm + dm
      ke = 0.5d0*( (0.5d0*(u(k,i,j) + u(k,ip,j)))**2               &
                 + (0.5d0*(v(k,i,j) + v(k,i,jp)))**2               &
                 + delta_v*(above1(k)*w(kp,i,j) + below1(k)*w(k,i,j))**2 )
      thetabar = above1(k)*theta(kp,i,j) + below1(k)*theta(k,i,j)
      internaleng = internaleng + cv*exner(k,i,j)*thetabar*dm
      potentialeng = potentialeng + phi(k,i,j)*dm
      kineticeng = kineticeng + ke*dm
      totalent = totalent + cp*log(thetabar)*dm
      urcosp = ( 0.5d0*(u(k,i,j) + u(k,ip,j))*(cosp(j)*cosr + sinp(j)*SIN(xp(i))*sinr)   &
               + 0.5d0*(v(k,i,j) + v(k,i,jp))*COS(xp(i))*sinr )                          &
             *rp(k,i,j)
      wr2c2 = rotatn*(cosp(j)*rp(k,i,j))**2
      totalam = totalam + (wr2c2 + urcosp)*dm
    ENDDO
  ENDDO
ENDDO

print *,'ke approx = ',kineticeng
kineticeng = 0.0d0

! Kinetic energy contributions from u, v and w points
DO j = 1, ny
  DO i = 1, nx
    im = i - 1
    IF (i == 1) im = nx
    DO k = 1, nz
      vv  = 0.5d0*(volume(k,im,j) + volume(k,i,j))
      kineticeng = kineticeng + (0.5d0*u(k,i,j)**2)*0.5d0*(rho(k,i,j) + rho(k,im,j))*vv
    ENDDO
  ENDDO
ENDDO
DO j = 2, ny
  jm = j - 1
  DO i = 1, nx
    DO k = 1, nz
      vv  = 0.5d0*(volume(k,i,jm) + volume(k,i,j))
      kineticeng = kineticeng + (0.5d0*v(k,i,j)**2)*0.5d0*(rho(k,i,j) + rho(k,i,jm))*vv
    ENDDO
  ENDDO
ENDDO
DO j = 1, ny
  DO i = 1, nx
    DO k = 2, nz
      km = k - 1
      ra = rp(k ,i,j)
      rb = rp(km,i,j)
      vv  = (ra*ra*ra - rb*rb*rb)*area(j)/3.0d0
      kineticeng = kineticeng + (0.5d0*w(k,i,j)**2)*(above2(k)*rho(k,i,j) + below2(k)*rho(km,i,j))*vv
    ENDDO
  ENDDO
ENDDO


totaleng = internaleng + potentialeng + kineticeng



WRITE(chandiag,'(I10,10E22.14)') istep, totalm, internaleng, potentialeng,      &
                                       kineticeng, totaleng, unavailie, unavailpe, &
                                       totalent, totalam, imbalance

PRINT *,'  Total volume           = ',totvol
PRINT *,'  Total mass             = ',totalm,       '      rate of change = ',(totalm - totalm0)/dt
PRINT *,'  Internal energy        = ',internaleng,  '      rate of change = ',1d-15*(internaleng - internaleng0)/dt
PRINT *,'  Potential energy       = ',potentialeng, '      rate of change = ',1d-15*(potentialeng - potentialeng0)/dt
PRINT *,'  Kinetic energy         = ',kineticeng,   '      rate of change = ',1d-15*(kineticeng - kineticeng0)/dt
PRINT *,'  Total energy           = ',totaleng,     '      rate of change = ',1d-15*(totaleng - totaleng0)/dt
PRINT *,'  Total entropy          = ',totalent,     '      rate of change = ',(totalent - totalent0)/dt
PRINT *,'  Total angular momentum = ',totalam,      '      rate of change = ',(totalam - totalam0)/dt
PRINT *,'  Hydrostatic imbalance  = ',imbalance

! Save values for next step
totalm0 = totalm
internaleng0 = internaleng
potentialeng0 = potentialeng
kineticeng0 = kineticeng
totaleng0 = totaleng
totalent0 = totalent
totalam0 = totalam


END SUBROUTINE globaldiags

! ============================================================

SUBROUTINE diag_imbal

! Compute and save a measure of the departure from hydrostatic balance

USE timestep
USE state
USE departure
USE work
USE refstate
USE globdiag

IMPLICIT NONE

INTEGER :: k, km, i, j
REAL*8 :: vol, ra, rb, vv, accel, vvx, aax
REAL*8, ALLOCATABLE :: temp1(:,:,:), temp2(:,:,:), &
                       ruatw(:,:,:), rvatw(:,:,:)


! Calculate rate of change of w along trajectory.
! Assume w0 is still the old value of w, and assume w departure points
! are still valid.

! Allocate space
ALLOCATE(ruatw(nzp,nx,ny),rvatw(nzp,nx,ny),temp1(nzp,nx,ny),temp2(nzp,nx,ny))
! To avoid the need for a separate interpolation routine, dummy calculations
! are carried out for u, v and theta fields
CALL vecatw(u0,v0,ruatw,rvatw)
! Interpolate to departure points
CALL lagrangewd(xdepw,ydepw,edepw,ruatw,rvatw,w0,theta0, &
               temp1,temp2,rhsdepw,rhsdeptheta,dthetabardr_ref)

! Deallocate space
DEALLOCATE(ruatw,rvatw,temp1,temp2)

! Integrate vertical acceleration term over the domain to obtain RMS value
vvx = 0.0d0
aax = 0.0d0
vol = 0.0d0
imbalance = 0.0d0
DO j = 1, ny
  DO i = 1, nx
    DO k = 2, nz
      km = k - 1
      ra = rp(k ,i,j)
      rb = rp(km,i,j)
      vv  = (ra*ra*ra - rb*rb*rb)*area(j)/3.0d0
      accel = (w(k,i,j) - rhsdepw(k,i,j))/dt
      vvx = max(vvx,vv)
      aax = max(aax,abs(accel))
      vol = vol + vv
      imbalance = imbalance + accel*accel*vv
    ENDDO
  ENDDO
ENDDO
imbalance = SQRT(imbalance/vol)


END SUBROUTINE diag_imbal

! ============================================================

SUBROUTINE heldsuarez

! Add in forcing/relaxation terms for Held-Suarez dynamical
! core test case

USE constants
USE timestep
USE state

IMPLICIT NONE

INTEGER :: i, im, j, jm, k, km
REAL*8 :: exner(nz,nx,ny), exsurf(nx,ny), psurf(nx,ny), sigma(nz,nx,ny)
REAL*8 :: a, b, kf, kv, ka, ks, kt, sigmab, deltat, deltatheta, sig, exn, thetaeq, &
          tequat, tstrat


! Define constants
sigmab = 0.7d0
kf = 1.0d0/86400.0d0
ka = 1.0d0/(40.0d0*86400.0d0)
ks = 1.0d0/(4.0d0*86400.0d0)
deltat = 60.0d0
deltatheta = 10.0d0
tequat = 315.0d0
tstrat = 200.0d0


! Compute exner
CALL findexner(rho,theta,exner)

! Extrapolate to surface
a = -etap(1)/(etap(2) - etap(1))
b = 1.0d0 - a
exsurf = a*exner(2,:,:) + b*exner(1,:,:)

! Convert surface exner to pressure/p00
psurf = exsurf**(1.0d0/kappa)

! Compute 3D sigma field
DO k = 1, nz
  sigma(k,:,:) = (exner(k,:,:)**(1.0d0/kappa))/psurf
ENDDO

! Relax u
DO j = 1, ny
  DO i = 1, nx
    im = i - 1
    IF (i == 1) im = nx
    DO k = 1, nz
      sig = 0.5d0*(sigma(k,im,j) + sigma(k,i,j))
      kv = kf*MAX(0.0d0,(sig - sigmab)/(1.0d0 - sigmab))
      u(k,i,j) = u(k,i,j)/(1.0 + kv*dt)
    ENDDO
  ENDDO
ENDDO
    
! Relax v
DO j = 2, ny
  jm = j - 1
  DO i = 1, nx
    DO k = 1, nz
      sig = 0.5d0*(sigma(k,i,jm) + sigma(k,i,j))
      kv = kf*MAX(0.0d0,(sig - sigmab)/(1.0d0 - sigmab))
      v(k,i,j) = v(k,i,j)/(1.0 + kv*dt)
    ENDDO
  ENDDO
ENDDO

! Relax theta
DO j = 1, ny
  DO i = 1, nx
    ! Bottom level
    k = 1
    sig = 1.0d0
    exn = exsurf(i,j)
    kt = ka + (ks - ka)*MAX(0.0d0,(sig - sigmab)/(1.0d0 - sigmab))*COS(geolatp(i,j))**4
    thetaeq = tequat - deltat*SIN(geolatp(i,j))**2 &
                      - deltatheta*(LOG(exn)/kappa)*COS(geolatp(i,j))**2
    thetaeq = MAX(tstrat/exn, thetaeq)
    theta(k,i,j) = (theta(k,i,j) + kt*dt*thetaeq)/(1.0 + kt*dt)
    ! Middle levels
    DO k = 2, nz
      km = k - 1
      sig = below2(k)*sigma(km,i,j) + above2(k)*sigma(k,i,j)
      exn = below2(k)*exner(km,i,j) + above2(k)*exner(k,i,j)
      kt = ka + (ks - ka)*MAX(0.0d0,(sig - sigmab)/(1.0d0 - sigmab))*COS(geolatp(i,j))**4
      thetaeq = tequat - deltat*SIN(geolatp(i,j))**2 &
                       - deltatheta*(LOG(exn)/kappa)*COS(geolatp(i,j))**2
      thetaeq = MAX(tstrat/exn, thetaeq)
      theta(k,i,j) = (theta(k,i,j) + kt*dt*thetaeq)/(1.0 + kt*dt)
    ENDDO
    ! Top level
    k = nzp
    a = (etaw(nzp) - etap(nz-1))/(etap(nz) - etap(nz-1))
    b = 1.0d0 - a
    sig = a*sigma(nz,i,j) + b*sigma(nz-1,i,j)
    exn = a*exner(nz,i,j) + b*exner(nz-1,i,j)
    kt = ka + (ks - ka)*MAX(0.0d0,(sig - sigmab)/(1.0d0 - sigmab))*COS(geolatp(i,j))**4
    thetaeq = tequat - deltat*SIN(geolatp(i,j))**2 &
                      - deltatheta*(LOG(exn)/kappa)*COS(geolatp(i,j))**2
    thetaeq = MAX(tstrat/exn, thetaeq)
    theta(k,i,j) = (theta(k,i,j) + kt*dt*thetaeq)/(1.0 + kt*dt)
  ENDDO
ENDDO


END SUBROUTINE heldsuarez

! ============================================================

SUBROUTINE tlheldsuarez

! Add in forcing/relaxation terms for tidally-locked Held-Suarez dynamical
! core test case

USE constants
USE timestep
USE state

IMPLICIT NONE

INTEGER :: i, im, j, jm, k, km
REAL*8 :: exner(nz,nx,ny), exsurf(nx,ny), psurf(nx,ny), sigma(nz,nx,ny)
REAL*8 :: a, b, kf, kv, ka, ks, kt, sigmab, deltat, deltatheta, sig, exn, thetaeq, &
          tequat, tstrat


! Define constants
sigmab = 0.7d0
kf = 1.0d0/86400.0d0
ka = 1.0d0/(40.0d0*86400.0d0)
ks = 1.0d0/(4.0d0*86400.0d0)
deltat = 60.0d0
deltatheta = 10.0d0
tequat = 315.0d0
tstrat = 200.0d0


! Compute exner
CALL findexner(rho,theta,exner)

! Extrapolate to surface
a = -etap(1)/(etap(2) - etap(1))
b = 1.0d0 - a
exsurf = a*exner(2,:,:) + b*exner(1,:,:)

! Convert surface exner to pressure/p00
psurf = exsurf**(1.0d0/kappa)

! Compute 3D sigma field
DO k = 1, nz
  sigma(k,:,:) = (exner(k,:,:)**(1.0d0/kappa))/psurf
ENDDO

! Relax u
DO j = 1, ny
  DO i = 1, nx
    im = i - 1
    IF (i == 1) im = nx
    DO k = 1, nz
      sig = 0.5d0*(sigma(k,im,j) + sigma(k,i,j))
      kv = kf*MAX(0.0d0,(sig - sigmab)/(1.0d0 - sigmab))
      u(k,i,j) = u(k,i,j)/(1.0 + kv*dt)
    ENDDO
  ENDDO
ENDDO
    
! Relax v
DO j = 2, ny
  jm = j - 1
  DO i = 1, nx
    DO k = 1, nz
      sig = 0.5d0*(sigma(k,i,jm) + sigma(k,i,j))
      kv = kf*MAX(0.0d0,(sig - sigmab)/(1.0d0 - sigmab))
      v(k,i,j) = v(k,i,j)/(1.0 + kv*dt)
    ENDDO
  ENDDO
ENDDO

! Relax theta
DO j = 1, ny
  DO i = 1, nx
    ! Bottom level
    k = 1
    sig = 1.0d0
    exn = exsurf(i,j)
    kt = ka + (ks - ka)*MAX(0.0d0,(sig - sigmab)/(1.0d0 - sigmab))*COS(geolatp(i,j))**4
    thetaeq = tequat + deltat*COS(geolonp(i,j) - pi)*COS(geolatp(i,j)) &
                      - deltatheta*(LOG(exn)/kappa)*COS(geolatp(i,j))**2
    thetaeq = MAX(tstrat/exn, thetaeq)
    theta(k,i,j) = (theta(k,i,j) + kt*dt*thetaeq)/(1.0 + kt*dt)
    ! Middle levels
    DO k = 2, nz
      km = k - 1
      sig = below2(k)*sigma(km,i,j) + above2(k)*sigma(k,i,j)
      exn = below2(k)*exner(km,i,j) + above2(k)*exner(k,i,j)
      kt = ka + (ks - ka)*MAX(0.0d0,(sig - sigmab)/(1.0d0 - sigmab))*COS(geolatp(i,j))**4
      thetaeq = tequat + deltat*COS(geolonp(i,j) - pi)*COS(geolatp(i,j)) &
                       - deltatheta*(LOG(exn)/kappa)*COS(geolatp(i,j))**2
      thetaeq = MAX(tstrat/exn, thetaeq)
      theta(k,i,j) = (theta(k,i,j) + kt*dt*thetaeq)/(1.0 + kt*dt)
    ENDDO
    ! Top level
    k = nzp
    a = (etaw(nzp) - etap(nz-1))/(etap(nz) - etap(nz-1))
    b = 1.0d0 - a
    sig = a*sigma(nz,i,j) + b*sigma(nz-1,i,j)
    exn = a*exner(nz,i,j) + b*exner(nz-1,i,j)
    kt = ka + (ks - ka)*MAX(0.0d0,(sig - sigmab)/(1.0d0 - sigmab))*COS(geolatp(i,j))**4
    thetaeq = tequat + deltat*COS(geolonp(i,j) - pi)*COS(geolatp(i,j)) &
                      - deltatheta*(LOG(exn)/kappa)*COS(geolatp(i,j))**2
    thetaeq = MAX(tstrat/exn, thetaeq)
    theta(k,i,j) = (theta(k,i,j) + kt*dt*thetaeq)/(1.0 + kt*dt)
  ENDDO
ENDDO


END SUBROUTINE tlheldsuarez

! ============================================================

SUBROUTINE masspertheta

! Compute and output mass per unit theta for a range of isentropic layers.
! Assume interfaces are given in ascending order in thetalayers.

USE lagdiagnostics
USE channels, ONLY : chanmasspt
USE timestep, ONLY : dt, istep

IMPLICIT NONE

INTEGER :: k


! Main calculation is done in masspertheta2
CALL masspertheta2(thetalayers,massbelow,nlayers)

! Convert to mass in theta layers
mpt(1:nlayers) = massbelow(2:nlayers+1) - massbelow(1:nlayers)

!print *,' '
!print *,'Mass below theta layers'
!print *, massbelow

!WRITE(chanmasspt,'(A,I10,A,F20.15,1X)') 'ISTEP = ', istep,' and T/hours = ',dt*istep/3600.0d0
!DO k = 1, nlayersp
!  WRITE(chanmasspt,'(I10,1x,10(ES25.15E3,1X))') k, thetalayers(k), massbelow(k)
!END DO
!WRITE(chanmasspt,'(A)') '   '

IF (istep == 0) then
  WRITE(chanmasspt,888) thetalayers
ENDIF
WRITE(chanmasspt,888) massbelow

888 FORMAT(33E15.7)


END SUBROUTINE masspertheta

! ============================================================

SUBROUTINE masspertheta2(thetalayers,massbelow,nlayers)

! Compute mass per unit theta for a range of isentropic layers.
! Assume interfaces are given in ascending order in thetalayers.
! If an interface intersects a grid box then count a fraction of its
! mass assuming mass is distributed linearly in theta

USE state

IMPLICIT NONE

INTEGER, INTENT(IN) :: nlayers
REAL*8, INTENT(IN) :: thetalayers(nlayers+1)
REAL*8, INTENT(OUT) :: massbelow(nlayers+1)
INTEGER :: i, j, k, kp, ilayer, nlayersp
REAL*8 :: tlev, tmin, tmax, mcell, mbelow


nlayersp = nlayers + 1

! First calculate mass below given theta levels
DO ilayer = 1, nlayersp
  tlev = thetalayers(ilayer)
  mbelow = 0.0d0
  DO j = 1, ny
    DO i = 1, nx
      DO k = 1, nz
        mcell = rho(k,i,j)*volume(k,i,j)
        kp = k + 1
	tmin = MIN(theta(k,i,j),theta(kp,i,j))
	tmax = MAX(theta(k,i,j),theta(kp,i,j))
	IF (tmax < tlev) THEN
	  mbelow = mbelow + mcell
	ELSEIF (tmin < tlev) THEN
	  mbelow = mbelow + mcell*(tlev - tmin)/(tmax - tmin)
	ELSE
	  ! No contribution from this cell
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  massbelow(ilayer) = mbelow
ENDDO

print *,'Range of theta = ',minval(theta),maxval(theta)


END SUBROUTINE masspertheta2

! ============================================================
!
SUBROUTINE unavailPE2
!
! Compute unavailable potential energy
! Orography is taken into account. However, we assume a spherical geopotential;
! for a spheroidal planet the calculation would need to be modified
!
USE grid
USE state, ONLY : theta
USE constants
USE globdiag
!
IMPLICIT NONE
!
INTEGER, PARAMETER :: nlevupe = 160, nlevupep = nlevupe+1
INTEGER :: k, iter, kmin, kmax, kmaxp
REAL*8 :: rr(nlevupep), rrbar(nlevupe), dr, vv(nlevupep), aa(nlevupep), &
          rho(nlevupe), thetabar, exner(nlevupe), a, b, phi(nlevupe), &
          imb(nlevupep), dibydr, aaa(nlevupep), bbb(nlevupep), ccc(nlevupep), &
          rrr(nlevupep), rinc(nlevupep), const, ddrr(nlevupe), ddrinc(nlevupe), &
          alpha, cv, upe, uie, utote, mass, entropy, dm, massbelow(nlevupep), &
          thetalayers(nlevupep), dtheta, thetamax, thetamin


cv = cp - gascon


! Set theta levels to span range of theta to ensure smooth variation of
! layer thickness
thetamin = minval(theta)
thetamax = maxval(theta)
! Use levels evenly spaced in theta
dtheta = (thetamax - thetamin)/(nlevupe-2)
DO k = 1, nlevupep
  thetalayers(k) = thetamin + (k-2)*dtheta
ENDDO
!print *,'thetalayers = ',thetalayers

! Determine mass below each theta level
CALL masspertheta2(thetalayers,massbelow,nlevupe)

! Sanity checks: must be zero mass below bottom theta level, and all mass
! below top theta level; the only way to be sure of this last one is if
! there's no mass in the top layer
IF (massbelow(1) > 0.0d0) THEN
  PRINT *,'Bottom diagnostic theta level should have a lower value'
  STOP
ENDIF
IF (massbelow(nlevupep) > massbelow(nlevupe)) THEN
  PRINT *,'Top diagnostic theta level should have a higher value'
  STOP
ENDIF

! First determine range of theta layers containing nonzero mass
DO k = 1, nlevupe
  IF (massbelow(k) == 0.0d0) kmin = k
  IF (massbelow(k+1) > massbelow(k)) kmax = k
ENDDO
kmaxp = kmax + 1

! First guess for locations of theta levels
rr(1:kmin) = MINVAL(rsurf)
rr(kmaxp:nlevupep) = rearth + domain
dr = (rr(nlevupep) - rr(1))/(kmaxp - kmin)
DO k = kmin, kmaxp
  rr(k) = rr(1) + (k-kmin)*dr
ENDDO

! Iterate to move levels towards a hydrostatically balanced state
DO iter = 1, 25

  ! Find level spacing
  ddrr(kmin:kmax) = rr(kmin+1:kmaxp) - rr(kmin:kmax)

  ! Interpolate to find locations of rho levels
  ! Uniform stretching at bottom and top, otherwise cubic interpolation
  alpha = SQRT(ddrr(kmin+1)/ddrr(kmin))
  b = 0.5d0 !alpha/(1.0d0 + alpha)
  a = 1.0d0 - b
  rrbar(kmin) = b*rr(kmin) + a*rr(kmin+1)
  DO k = kmin+1, kmax-1
    rrbar(k) = (9.0d0*(rr(k) + rr(k+1)) - (rr(k-1) + rr(k+2)))/16.0d0
  ENDDO
  alpha = SQRT(ddrr(kmax)/ddrr(kmax-1))
  b = 0.5d0 !alpha/(1.0d0 + alpha)
  a = 1.0d0 - b
  rrbar(kmax) = b*rr(kmax) + a*rr(kmaxp)

  ! Find volumes below each theta level and area of each theta level
  DO k = 1, nlevupep
    CALL volbelow(rr(k),vv(k),aa(k))
  ENDDO

  ! Compute rho and theta at rho levels, hence exner
  DO k = kmin, kmax
    rho(k) = (massbelow(k+1) - massbelow(k))/(vv(k+1) - vv(k))
    a = (rrbar(k) - rr(k))/(rr(k+1) - rr(k))
    b = 1.0d0 - a
    thetabar = a*thetalayers(k+1) + b*thetalayers(k)
    exner(k) = (gascon*rho(k)*thetabar/p00)**kbyonemk
  ENDDO

  ! Set up geopotentials
  CALL findphicol(kmaxp-kmin,rrbar(kmin:kmax),phi(kmin:kmax))

  ! Compute hydrostatic imbalance
  DO k = kmin+1, kmax
    dr = rrbar(k) - rrbar(k-1)
    imb(k) = (cp*thetalayers(k)*(exner(k) - exner(k-1)) + (phi(k) - phi(k-1)))/dr
  ENDDO
  ! print *,'Iteration ',iter,' Max imbalance ',MAXVAL(ABS(imb(kmin+1:kmax)))

  ! Set up a tridiagonal problem for increments to r to reduce imbalance
  aaa(kmin) = 0.0d0
  bbb(kmin) = 0.0d0
  ccc(kmin) = 1.0d0
  rrr(kmin) = 0.0d0
  DO k = kmin+1, kmax
    dr = rrbar(k) - rrbar(k-1)
    const = kbyonemk*cp*thetalayers(k)/dr
    aaa(k) = -const*aa(k+1)*exner(k  )/(vv(k+1) - vv(k  ))
    bbb(k) = -const*aa(k-1)*exner(k-1)/(vv(k  ) - vv(k-1))
    ccc(k) =  const*aa(k)*(exner(k  )/(vv(k+1) - vv(k  ))    &
                          +exner(k-1)/(vv(k  ) - vv(k-1)))
    rrr(k) = -imb(k)
  ENDDO
  aaa(kmaxp) = 0.0d0
  bbb(kmaxp) = 0.0d0
  ccc(kmaxp) = 1.0d0
  rrr(kmaxp) = 0.0d0
  CALL trisolveb(rinc(kmin:kmaxp),bbb(kmin:kmaxp),ccc(kmin:kmaxp),aaa(kmin:kmaxp),rrr(kmin:kmaxp),kmaxp+1-kmin)
  ddrinc(kmin:kmax) = rinc(kmin+1:kmaxp) - rinc(kmin:kmax)
  ! Limit increments to prevent spacing reducing to less than half current spacing
  alpha = MIN(1.0d0, -0.5d0/MINVAL(ddrinc(kmin:kmax)/ddrr(kmin:kmax)))
  rr(kmin:kmaxp) = rr(kmin:kmaxp) + alpha*rinc(kmin:kmaxp)

ENDDO

! Now compute the PE and IE of this balanced state;
! also mass and entropy as a cross check

! Interpolate to find locations of rho levels
! Uniform stretching at bottom and top, otherwise cubic interpolation
alpha = SQRT(ddrr(kmin+1)/ddrr(kmin))
b = 0.5d0 !alpha/(1.0d0 + alpha)
a = 1.0d0 - b
rrbar(kmin) = b*rr(kmin) + a*rr(kmin+1)
DO k = kmin+1, kmax-1
  rrbar(k) = (9.0d0*(rr(k) + rr(k+1)) - (rr(k-1) + rr(k+2)))/16.0d0
ENDDO
alpha = SQRT(ddrr(kmax)/ddrr(kmax-1))
b = 0.5d0 !alpha/(1.0d0 + alpha)
a = 1.0d0 - b
rrbar(kmax) = b*rr(kmax) + a*rr(kmaxp)

! Find volumes below each theta level and area of each theta level
DO k = 1, nlevupep
  CALL volbelow(rr(k),vv(k),aa(k))
ENDDO

! Set up geopotentials
CALL findphicol(kmaxp-kmin,rrbar(kmin:kmax),phi(kmin:kmax))

! Compute rho and theta at rho levels, hence exner, T, and PE
mass = 0.0d0
upe = 0.0d0
uie = 0.0d0
entropy = 0.0d0
DO k = kmin, kmax
  rho(k) = (massbelow(k+1) - massbelow(k))/(vv(k+1) - vv(k))
  a = (rrbar(k) - rr(k))/(rr(k+1) - rr(k))
  b = 1.0d0 - a
  thetabar = a*thetalayers(k+1) + b*thetalayers(k)
  exner(k) = (gascon*rho(k)*thetabar/p00)**kbyonemk
  dm = massbelow(k+1) - massbelow(k)
  mass = mass + dm
  upe = upe + phi(k)*dm
  uie = uie + cv*exner(k)*thetabar*dm
  entropy = entropy + cp*LOG(thetabar)*dm
ENDDO

print *,'unavailPE2:'
print *,'mass = ',mass,' volume = ',vv(nlevupep)
print *,'uie = ',uie
print *,'upe = ',upe
print *,'total ue = ',upe + uie
print *,'entropy = ',entropy

! Save for output along with other global diags
unavailie = uie
unavailpe = upe

!
RETURN
END SUBROUTINE unavailPE2
!
! ============================================================
!
SUBROUTINE volbelow(r,v,a)

! Compute the volume below radius r and the exposed area at radius r

USE grid, ONLY : rsurf, area, nx, ny

IMPLICIT NONE

REAL*8, INTENT(IN) :: r
REAL*8, INTENT(OUT) :: v, a

INTEGER :: i, j
REAL*8 :: rs

v = 0.0d0
a = 0.0d0

DO j = 1, ny
  DO i = 1, nx
    rs = rsurf(i,j)
    IF (r >= rsurf(i,j)) THEN
      a = a + r*r*area(j)
      v = v + area(j)*(r*r*r - rs*rs*rs)/3.0d0
    ENDIF
  ENDDO
ENDDO


END SUBROUTINE volbelow
!
! ============================================================
!
SUBROUTINE findphicol(n,rphi,pphi)
!
! IK/JT: calculates geopotential for a single column
! Phi is assumed to be a function only of r (spherical
! geopotential approximation).
! This routine should be consistent with routine setphi
!
USE constants
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: n
REAL*8, INTENT(IN)  :: rphi(n)
REAL*8, INTENT(OUT) :: pphi(n)
INTEGER :: k
!
! Calculate geopotential
!
! Spherical planet with inverse square law
pphi = gravity*rearth*(1.0d0 - rearth/rphi)
!
! Spheroidal planet with constant vertical gradient
! pphi = gravity*(rphi - rearth)
!
RETURN
END SUBROUTINE findphicol
!
! ==========================================================
!
SUBROUTINE start_diagnostics
!
! IK/JT: Subroutine to start lagrangian diagnostics
!
USE state
USE lagdiagnostics
USE trajectories
USE channels, ONLY : chanmasspt, chanscatterm

IMPLICIT NONE

INTEGER :: i, k, ilayer
REAL*8 :: thetamin = 210.0d0, thetamax = 810.0d0, dtheta
REAL*8 :: tempu(nz,nx,ny), tempp(nz,nx,ny)

!
! Initialize theta levels for cross-theta fluxes
! Use levels evenly spaced in theta
dtheta = (thetamax - thetamin)/nlayers
DO k = 1, nlayersp
  thetalayers(k) = thetamin + (k-1)*dtheta
ENDDO
!print *,'thetalayers = ',thetalayers

OPEN(chanmasspt,FILE = 'masspertheta.txt')


! Potential vorticity (to initialize tracers and labels)
CALL potvort

! Initialize tracers to agree with theta and with PV
IF (ntracers > 0) THEN
  CALL averw2p(nx,ny,theta0,tracers0(:,:,:,1))
  tracers(:,:,:,1) = tracers0(:,:,:,1)
ENDIF
IF (ntracers > 1) THEN
  CALL averv2p(nz,nx,pv,tempu)
  CALL averu2p(nz,ny,tempu,tracers0(:,:,:,2))
  tracers(:,:,:,2) = tracers0(:,:,:,2)
ENDIF


! Set trajectory labels
IF (ntrajectories > 0) THEN
  IF (nlabels > 0) THEN
    ! First label is theta
    CALL averw2p(nx,ny,theta0,tempp)
    CALL labeltraj(ntrajectories,xtraj,ytraj,etraj,tempp,labels(1,1))
  ENDIF
  IF (nlabels > 1) THEN
    ! Second label is PV
    CALL averv2p(nz,nx,pv,tempu)
    CALL averu2p(nz,ny,tempu,tempp)
    CALL labeltraj(ntrajectories,xtraj,ytraj,etraj,tempp,labels(1,2))
  ENDIF
ENDIF

OPEN(chanscatterm,FILE='scatter.m')


RETURN
END SUBROUTINE start_diagnostics
!
! ============================================================

SUBROUTINE dump(f,ytitle,nx,ny)

! Dump a 2D field for IDL plotting

USE channels

IMPLICIT NONE

INTEGER, INTENT(IN) :: nx, ny
REAL*8, INTENT(IN) :: f(nx,ny)
REAL :: fstar4(nx,ny)
CHARACTER*(*), INTENT(IN) :: ytitle


! Convert to single precision to reduce size of output
! and improve readability!
fstar4 = f
WRITE(chandump2,*) nx,ny
WRITE(chandump2,*) ytitle
WRITE(chandump2,*) fstar4


END SUBROUTINE dump

! ===============================================================

SUBROUTINE dump1d(z,f,ytitle,n)

! Dump a 1D field for IDL plotting

USE channels

IMPLICIT NONE

INTEGER, INTENT(IN) :: n
REAL*8, INTENT(IN) :: z(n), f(n)
REAL :: zstar4(n), fstar4(n)
CHARACTER*(*), INTENT(IN) :: ytitle

! Convert to single precision to reduce size of output
! and improve readability!
zstar4 = z
fstar4 = f

WRITE(chandump1,*) n
WRITE(chandump1,*) ytitle
WRITE(chandump1,*) zstar4
WRITE(chandump1,*) fstar4

END SUBROUTINE dump1d

! ===============================================================

SUBROUTINE dumpm(f,ytitle,nx,ny)

! Dump a 2D field for matlab plotting

USE channels

IMPLICIT NONE

INTEGER, INTENT(IN) :: nx, ny
REAL*8, INTENT(IN) :: f(nx,ny)
REAL :: fstar4(nx,ny)
CHARACTER*(*), INTENT(IN) :: ytitle

! Convert to single precision to reduce size of output
! and improve readability!
fstar4 = f

WRITE(chandumpm,*) 'nx = ',nx,';'
WRITE(chandumpm,*) 'ny = ',ny,';'
WRITE(chandumpm,*)  ytitle//' = [ ...'
WRITE(chandumpm,888) fstar4
WRITE(chandumpm,*) ' ];'
WRITE(chandumpm,*) '[ min('//ytitle//') max('//ytitle//') ]'
WRITE(chandumpm,*) 'z = reshape('//ytitle//',nx,ny);'
WRITE(chandumpm,*) 'contour(z'')'
WRITE(chandumpm,*) 'title('''//ytitle//''')'
WRITE(chandumpm,*) 'pause'


888 FORMAT(E16.4)

END SUBROUTINE dumpm

! ===============================================================

SUBROUTINE dumpma(f,ytitle,nx,ny)

! Dump a 2D field for matlab plotting as one of an
! array of plots

USE channels

IMPLICIT NONE

INTEGER, INTENT(IN) :: nx, ny
REAL*8, INTENT(IN) :: f(nx,ny)
REAL :: fstar4(nx,ny)
CHARACTER*(*), INTENT(IN) :: ytitle
CHARACTER*5 :: rowhead(5)
CHARACTER*3 :: colhead(3)
INTEGER :: plotcount = 0, colcount = 0, row, colm1, pindex, prow, pcol

! Convert to single precision to reduce size of output
! and improve readability!
fstar4 = f

prow = 5
pcol = 8
rowhead(1) = 'rho'
rowhead(2) = 'theta'
rowhead(3) = 'u'
rowhead(4) = 'v'
rowhead(5) = 'w'
colhead(1) = 'err'
colhead(2) = 'rhs'
colhead(3) = 'inc'


plotcount = MODULO(plotcount,prow*pcol) + 1
colm1 = (plotcount-1)/prow
row = plotcount - colm1*prow
IF (row == 1) colcount = MODULO(colcount,3) + 1
pindex = pcol*(row-1) + colm1 + 1
WRITE(chandumpm,*) 'nx = ',nx,';'
WRITE(chandumpm,*) 'ny = ',ny,';'
WRITE(chandumpm,*)  ytitle//' = [ ...'
WRITE(chandumpm,888) fstar4
WRITE(chandumpm,*) ' ];'
WRITE(chandumpm,*) 'ymin = num2str(min('//ytitle//'));'
WRITE(chandumpm,*) 'ymax = num2str(max('//ytitle//'));'
WRITE(chandumpm,*) '[ min('//ytitle//') max('//ytitle//') ]'
WRITE(chandumpm,*) 'z = reshape('//ytitle//',nx,ny);'
IF (plotcount == 1) WRITE(chandumpm,*) 'subplot(1,1,1)'
WRITE(chandumpm,*) 'subplot(',prow,',',pcol,',',pindex,')'
WRITE(chandumpm,*) 'contour(z,4)'
!WRITE(chandumpm,*) 'title(['''//ytitle//''',''  '',ymin,''  '',ymax],''Fontsize'',7)'
WRITE(chandumpm,*) 'title([ymin,''  '',ymax],''Fontsize'',7)'
IF (colm1 == 0) WRITE(chandumpm,*) 'ylabel(''',rowhead(row),''')'
IF (row == 1) WRITE(chandumpm,*) 'text(0.4*ny,1.35*nx,''',colhead(colcount),''')'
WRITE(chandumpm,*) 'pause'


888 FORMAT(E16.4)

END SUBROUTINE dumpma

! ===============================================================

SUBROUTINE dumpxm(f,clatc,ytitle,nnx,nny,nx,ny)

! Dump for matlab plotting - NS cross section

USE channels

IMPLICIT NONE

INTEGER, INTENT(IN) :: nnx, nny, nx, ny
REAL*8, INTENT(IN) :: f(nx,ny), clatc(ny)
REAL :: fstar4(nny), lat(nny)
CHARACTER*(*), INTENT(IN) :: ytitle

! Convert to single precision to reduce size of output
! and improve readability!
lat = clatc(1:nny)
fstar4 = f(nnx/4,1:nny)

WRITE(chandumpm,*) 'ny = ',ny,';'
WRITE(chandumpm,*) 'lat = [ ...'
WRITE(chandumpm,888) lat
WRITE(chandumpm,*) ' ];'
WRITE(chandumpm,*)  ytitle//' = [ ...'
WRITE(chandumpm,888) fstar4
WRITE(chandumpm,*) ' ];'
WRITE(chandumpm,*) '[ min('//ytitle//') max('//ytitle//') ]'
WRITE(chandumpm,*) 'plot(lat,'//ytitle//')'
! WRITE(chandumpm,*) 'axis([0,2*pi,-1,1])
WRITE(chandumpm,*) 'title('''//ytitle//''')'
WRITE(chandumpm,*) 'pause'

888 FORMAT(E16.4)

END SUBROUTINE dumpxm

! ===============================================================

SUBROUTINE dump1dm(z,f,ytitle,n)

! Dump a 1D field for matlab plotting

USE channels

IMPLICIT NONE

INTEGER, INTENT(IN) :: n
REAL*8, INTENT(IN) :: z(n), f(n)
REAL :: zstar4(n), fstar4(n)
CHARACTER*(*), INTENT(IN) :: ytitle

! Convert to single precision to reduce size of output
! and improve readability!
zstar4 = z
fstar4 = f

WRITE(chandumpm,*) 'n = ',n,';'
WRITE(chandumpm,*) 'z = [ ...'
WRITE(chandumpm,888) zstar4
WRITE(chandumpm,*) ' ];'
WRITE(chandumpm,*)  ytitle//' = [ ...'
WRITE(chandumpm,888) fstar4
WRITE(chandumpm,*) ' ];'
WRITE(chandumpm,*) '[ min('//ytitle//') max('//ytitle//') ]'
WRITE(chandumpm,*) 'plot('//ytitle//',z)'
WRITE(chandumpm,*) 'title('''//ytitle//''')'
WRITE(chandumpm,*) 'pause'

888 FORMAT(E16.4)

END SUBROUTINE dump1dm

! ===============================================================

SUBROUTINE dumptrajm(x,y,e,n,ysec)

! Dump trajectories for matlab plotting
! They will be plotted on top of the previous plot plotted

! Options for ysec are 'xy', 'zx', or 'zy', which should be chosen to
! match the contour plot just plotted

USE grid, ONLY : ny, dx, dy
USE channels, ONLY : chandumpm

IMPLICIT NONE

INTEGER, INTENT(IN) :: n
REAL*8, INTENT(IN) :: x(n), y(n), e(n)
REAL :: star4(n)
CHARACTER*(*), INTENT(IN) :: ysec

print *,'*** dumptrajm needs to be updated when plotting routines use'
print *,'proper coordinates ***'

WRITE(chandumpm,*) 'xtraj = [ ...'
IF (ysec == 'xy') THEN
  star4 = x/dx + 0.5d0
  WRITE(chandumpm,889) star4
ELSE
  star4 = e
  WRITE(chandumpm,889) star4
ENDIF
WRITE(chandumpm,*) ' ];'
WRITE(chandumpm,*) 'ytraj = [ ...'
IF (ysec == 'zx') THEN
  star4 = x/dx +0.5d0
  WRITE(chandumpm,889) star4
ELSE
  star4 = y/dy + 0.5d0*(ny + 1)
  WRITE(chandumpm,889) star4
ENDIF
WRITE(chandumpm,*) ' ];'
WRITE(chandumpm,*) 'hold on'
WRITE(chandumpm,*) 'plot(xtraj,ytraj,''+'')'
WRITE(chandumpm,*) 'hold off'
WRITE(chandumpm,*) 'pause'

889 FORMAT(E16.4)


END SUBROUTINE dumptrajm

! ===============================================================

SUBROUTINE dumpscatterm

! Create scatter plots of dynamical fields vs tracers vs trajectory labels

USE state
USE trajectories
USE lagdiagnostics, ONLY : pv
USE channels, ONLY : chanscatterm

IMPLICIT NONE

REAL*8 :: tempp(nz,nx,ny), tempu(nz,nx,ny)
REAL*8 :: dyndata(ntrajectories), trcdata(ntrajectories)
REAL*4 :: star4(ntrajectories)

! It might be desirable to select a subset of trajectories for producing
! these diagnostics...

IF (ntrajectories > 0) THEN

  ! First theta
  ! Trajectory first label, if it exists, is the initial theta
  ! Average dynamical theta to p points then sample at trajectory locations
  CALL averw2p(nx,ny,theta,tempp)
  CALL labeltraj(ntrajectories,xtraj,ytraj,etraj,tempp,dyndata)
  IF (ntracers > 0) THEN
    ! And sample theta-like tracer at trajectory locations
    CALL labeltraj(ntrajectories,xtraj,ytraj,etraj,tracers(1,1,1,1),trcdata)
  ENDIF

  IF (nlabels > 0) THEN
    ! label vs dynamical field
    WRITE(chanscatterm,*) 'field1 = [ ...'
    star4 = labels(:,1)
    WRITE(chanscatterm,887) star4
    WRITE(chanscatterm,*) ' ];'
    WRITE(chanscatterm,*) 'field2 = [ ...'
    star4 = dyndata
    WRITE(chanscatterm,887) star4
    WRITE(chanscatterm,*) ' ];'
    WRITE(chanscatterm,*) 'plot(field1,field2,''+'')'
    WRITE(chanscatterm,*) 'xlabel(''Initial \theta'')'
    WRITE(chanscatterm,*) 'ylabel(''\theta'')'
    WRITE(chanscatterm,*) 'pause'
  ENDIF
  IF (nlabels > 0 .AND. ntracers > 0) THEN
    ! label vs tracer
    WRITE(chanscatterm,*) 'field1 = [ ...'
    star4 = labels(:,1)
    WRITE(chanscatterm,887) star4
    WRITE(chanscatterm,*) ' ];'
    WRITE(chanscatterm,*) 'field2 = [ ...'
    star4 = trcdata
    WRITE(chanscatterm,887) star4
    WRITE(chanscatterm,*) ' ];'
    WRITE(chanscatterm,*) 'plot(field1,field2,''+'')'
    WRITE(chanscatterm,*) 'xlabel(''Initial \theta'')'
    WRITE(chanscatterm,*) 'ylabel(''\theta tracer'')'
    WRITE(chanscatterm,*) 'pause'
  ENDIF
  IF (ntracers > 0) THEN
    ! dynamical field vs tracer
    WRITE(chanscatterm,*) 'field1 = [ ...'
    star4 = dyndata
    WRITE(chanscatterm,887) star4
    WRITE(chanscatterm,*) ' ];'
    WRITE(chanscatterm,*) 'field2 = [ ...'
    star4 = trcdata
    WRITE(chanscatterm,887) star4
    WRITE(chanscatterm,*) ' ];'
    WRITE(chanscatterm,*) 'plot(field1,field2,''+'')'
    WRITE(chanscatterm,*) 'xlabel(''\theta'')'
    WRITE(chanscatterm,*) 'ylabel(''\theta tracer'')'
    WRITE(chanscatterm,*) 'pause'
  ENDIF

  ! Next PV
  ! Trajectory second label is initial PV
  ! Compute current pv, average to p points, then sample at
  ! trajectory locations
  CALL potvort
  CALL averv2p(nz,nx,pv,tempu)
  CALL averu2p(nz,ny,tempu,tempp)
  CALL labeltraj(ntrajectories,xtraj,ytraj,etraj,tempp,dyndata)
  IF (ntracers > 1) THEN
    ! And sample pv-like tracer at trajectory locations
    CALL labeltraj(ntrajectories,xtraj,ytraj,etraj,tracers(1,1,1,2),trcdata)
  ENDIF

  IF (nlabels > 1) THEN
    ! label vs dynamical field
    WRITE(chanscatterm,*) 'field1 = [ ...'
    star4 = labels(:,2)
    WRITE(chanscatterm,887) star4
    WRITE(chanscatterm,*) ' ];'
    WRITE(chanscatterm,*) 'field2 = [ ...'
    star4 = dyndata
    WRITE(chanscatterm,887) star4
    WRITE(chanscatterm,*) ' ];'
    WRITE(chanscatterm,*) 'plot(field1,field2,''+'')'
    WRITE(chanscatterm,*) 'xlabel(''Initial PV'')'
    WRITE(chanscatterm,*) 'ylabel(''PV'')'
    WRITE(chanscatterm,*) 'pause'
  ENDIF
  IF (nlabels > 1 .AND. ntracers > 1) THEN
    ! label vs tracer
    WRITE(chanscatterm,*) 'field1 = [ ...'
    star4 = labels(:,2)
    WRITE(chanscatterm,887) star4
    WRITE(chanscatterm,*) ' ];'
    WRITE(chanscatterm,*) 'field2 = [ ...'
    star4 = trcdata
    WRITE(chanscatterm,887) star4
    WRITE(chanscatterm,*) ' ];'
    WRITE(chanscatterm,*) 'plot(field1,field2,''+'')'
    WRITE(chanscatterm,*) 'xlabel(''Initial PV'')'
    WRITE(chanscatterm,*) 'ylabel(''PV tracer'')'
    WRITE(chanscatterm,*) 'pause'
  ENDIF
  IF (ntracers > 1) THEN
    ! dynamical field vs tracer
    WRITE(chanscatterm,*) 'field1 = [ ...'
    star4 = dyndata
    WRITE(chanscatterm,887) star4
    WRITE(chanscatterm,*) ' ];'
    WRITE(chanscatterm,*) 'field2 = [ ...'
    star4 = trcdata
    WRITE(chanscatterm,887) star4
    WRITE(chanscatterm,*) ' ];'
    WRITE(chanscatterm,*) 'plot(field1,field2,''+'')'
    WRITE(chanscatterm,*) 'xlabel(''PV'')'
    WRITE(chanscatterm,*) 'ylabel(''PV tracer'')'
    WRITE(chanscatterm,*) 'pause'
  ENDIF

ENDIF


887 FORMAT(E16.4)


END SUBROUTINE dumpscatterm

! ===============================================================

SUBROUTINE dumptrajectories

! Output label and tracer information interpolate to trajectories
! ready for subsequent processing and plotting

! nlabels should equal 2 and ntracers should equal 2, otherwise an error
! message is given

USE state
USE trajectories
USE timestep, ONLY : istep, dt
USE lagdiagnostics, ONLY : pv
USE channels, ONLY : chantraj

INTEGER :: itraj
REAL*8 :: tempp(nz,nx,ny), tempu(nz,nx,ny)
REAL*8 :: dyndata1(ntrajectories), trcdata1(ntrajectories)
REAL*8 :: dyndata2(ntrajectories), trcdata2(ntrajectories)
CHARACTER*19 :: yfilename

! -----


IF (ntrajectories > 0) THEN


  IF (ntracers < 2 .OR. nlabels < 2) THEN
    PRINT *,'To use routine dumptrajectories both ntracers and nlabels'
    PRINT *,'must be at least 2'
    STOP
  ENDIF


  WRITE(yfilename,'(''trajectories_'',i6.6)') istep
  OPEN(chantraj,FILE=yfilename)


  ! First theta
  ! Trajectory first label, if it exists, is the initial theta
  ! Average dynamical theta to p points then sample at trajectory locations
  CALL averw2p(nx,ny,theta,tempp)
  CALL labeltraj(ntrajectories,xtraj,ytraj,etraj,tempp,dyndata1)
  IF (ntracers > 0) THEN
    ! And sample theta-like tracer at trajectory locations
    CALL labeltraj(ntrajectories,xtraj,ytraj,etraj,tracers(1,1,1,1),trcdata1)
  ENDIF

  ! Next PV
  ! Trajectory second label is initial PV
  ! Compute current pv, average to p points, then sample at
  ! trajectory locations
  CALL potvort
  CALL averv2p(nz,nx,pv,tempu)
  CALL averu2p(nz,ny,tempu,tempp)
  CALL labeltraj(ntrajectories,xtraj,ytraj,etraj,tempp,dyndata2)
  IF (ntracers > 1) THEN
    ! And sample pv-like tracer at trajectory locations
    CALL labeltraj(ntrajectories,xtraj,ytraj,etraj,tracers(1,1,1,2),trcdata2)
  ENDIF

  DO itraj = 1, ntrajectories
    WRITE(chantraj,888) labels(itraj,1), dyndata1(itraj), trcdata1(itraj), &
                        labels(itraj,2), dyndata2(itraj), trcdata2(itraj)
  ENDDO


  CLOSE(chantraj)


ENDIF


888 FORMAT(6E16.6)



END SUBROUTINE dumptrajectories

! ===============================================================

SUBROUTINE dumprms(q,title,nz,nx,ny)

! Calculate and print the rms value of a field

INTEGER, INTENT(IN) :: nx, ny, nz
REAL*8, INTENT(IN) :: q(nz,nx,ny)
CHARACTER*(*), INTENT(IN) :: title
REAL*8 :: rms
real*8 :: qmax
integer :: kx, ix, jx

qmax = 0.0d0
rms = 0.0d0
DO i = 1, nx
  DO j = 1, ny
    DO k = 1, nz
      rms = rms + q(k,i,j)*q(k,i,j)
      if (abs(q(k,i,j)) .gt. qmax) then
        kx = k
	ix = i
	jx = j
	qmax = abs(q(k,i,j))
      endif
    ENDDO
  ENDDO
ENDDO
rms = sqrt(rms/(nz*nx*ny))

PRINT *,'RMS value of '//title//' = ',rms
! print *,'Max abs = ',qmax,'  at coords ',kx,ix,jx

END SUBROUTINE dumprms

! =========================================================

SUBROUTINE testdep

! Test departure pint calculation

USE grid
USE departure

IMPLICIT NONE
INTEGER :: i, j, k

CALL depart_all2

i = 48
j = 24
k = 2
write(22,99) xp(i), yp(j), etap(k), xdepp(k,i,j), ydepp(k,i,j), edepp(k,i,j)
DO k = 1, nz
  DO j = 1, ny
    DO i = 1, nx
      write(22,99) xp(i), yp(j), etap(k), xdepp(k,i,j), ydepp(k,i,j), edepp(k,i,j)
    ENDDO
  ENDDO
ENDDO

99 FORMAT(6F16.5)


END SUBROUTINE testdep

! =========================================================

SUBROUTINE testlagrange

! Test semi-Lagrangian advection based on cubic Lagrange interpolation

USE grid
USE departure

IMPLICIT NONE
INTEGER :: i, j, k, istep
REAL*8 :: q(nz,nx,ny), qnew(nz,nx,ny), temp(nx,ny), temp2(nz,ny)

! Compute departure points once (since wind does not change)
CALL depart_all2

! Initial condition for advected scalar
DO k = 1, nz
  DO j = 1, ny
    DO i = 1, nx
      q(k,i,j) = etap(k)*sinp(j)
    ENDDO
  ENDDO
ENDDO

temp(:,:) = q(5,1:nx,1:ny)
call dumpm(temp,'qlev5',nz,ny)
DO istep = 1, 100
  CALL lagrangep(xdepp,ydepp,edepp,q,qnew)
  q = qnew
  temp(:,:) = q(5,1:nx,1:ny)
  call dumpm(temp,'qlev5',nx,ny)
  print *,'Done step ',istep
ENDDO

END SUBROUTINE testlagrange

! =========================================================

SUBROUTINE testgrad

! Test calculation of gradient of a scalar

USE constants
USE grid

IMPLICIT NONE
INTEGER :: i, j, k
REAL*8 :: q(nz,nx,ny), dqdx(nz,nx,ny), dqdy(nz,nx,nyp), dqdz(nzp,nx,ny), &
          temp(nx,ny), temp2(nz,ny)


! Scalar field
DO k = 1, nz
  DO j = 1, ny
    DO i = 1, nx
      q(k,i,j) = -gravity*rearth*rearth/rp(k,i,j)
      ! q(k,i,j) = rearth*cosp(j)*SIN(xp(i))
    ENDDO
  ENDDO
ENDDO


CALL grad(q,dqdx,dqdy,dqdz)

temp(:,:) = dqdx(5,1:nx,1:ny)
call dumpm(temp,'dqdx',nx,ny)
temp(:,:) = dqdy(5,1:nx,1:ny)
call dumpm(temp,'dqdy',nx,ny)
temp(:,:) = dqdz(5,1:nx,1:ny)
call dumpm(temp,'dqdz',nx,ny)


END SUBROUTINE testgrad

! =========================================================

SUBROUTINE testdiv

! Test calculation of divergence of a velocity field

USE constants
USE state

IMPLICIT NONE
REAL*8 :: div(nz,nx,ny), temp(nx,ny), temp2(nz,ny)


! Velocity field is set up in initial

! Compute etadot
CALL findetadot(u,v,w,etadot)

! Now compute divergence
CALL divergence(u,v,etadot,div)

temp(:,:) = div(5,1:nx,1:ny)
call dumpm(temp,'div5',nx,ny)
temp(:,:) = div(1,1:nx,1:ny)
call dumpm(temp,'div1',nx,ny)
temp2(:,:) = div(1:nz,1,1:ny)
call dumpm(temp2,'divlong1',nz,ny)


END SUBROUTINE testdiv

! =========================================================

SUBROUTINE testcoriol

! Test calculation of gradient of a scalar

USE constants
USE state

IMPLICIT NONE
INTEGER :: i, j, k
REAL*8 :: cu(nz,nx,ny), cv(nz,nx,nyp), cw(nzp,nx,ny), &
          temp(nx,ny), temp2(nz,ny)



! Velocities and density set up in initial

! Compute Coriolis terms
CALL coriolis(u,v,w,rho,cu,cv,cw)

temp(:,:) = cu(1,1:nx,1:ny)
call dumpm(temp,'cu',nx,ny)
temp(:,:) = cv(1,1:nx,1:ny)
call dumpm(temp,'cv',nx,ny)
temp(:,:) = cw(2,1:nx,1:ny)
call dumpm(temp,'cw',nx,ny)


END SUBROUTINE testcoriol

! =========================================================

SUBROUTINE testmg

! To run the MG solver on a simple test case

USE constants
USE switches
USE work
USE increments

IMPLICIT NONE
INTEGER :: ng, i, j, k


! ----------------------------------------------------------

! Build the coefficients of the Helmholtz operator  
dv3d = REAL(delta_v)
CALL build_helm

! Set up a simple RHS field
rhs_helm = 0.0d0
DO i = 1, nx
  DO j = 1, ny
    DO k = 1, nz
      ! xp(i), yp(j), rp(k,i,j) are the longitude, latitude and radius coordinate of
      ! of the cell with indices i,j,k 
!      IF ((rp(k,i,j) - rearth) >= 1.6e4 .AND. yp(j) >= 1.0d0) THEN
      IF ((rp(k,i,j) - rearth) >= 1.6e4) THEN
        rhs_helm(k,i,j) = 1.0d0
      ENDIF
    ENDDO
  ENDDO
ENDDO


! Define number of grids to use
ng = p - 2

! Call solver
CALL mgsolve(exner_inc,rhs_helm,ng)


END SUBROUTINE testmg

! =========================================================

