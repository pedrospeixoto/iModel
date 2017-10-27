MODULE version

! Choice of numerical schemes

! 1. Semi-Lagrangian treatment of phi, following ENDGame
! 2. Semi-Lagrangian treatment of phi, in a form that
!    resembles the SLICE version as closely as possble
! 3. SLICE treatment of phi

INTEGER, PARAMETER :: ischeme = 1

! 0. C-grid correction to SLICE switched off
! 1. C-grid correction to SLICE included

INTEGER, PARAMETER :: cgridcorr = 0

! 0. No fix to SLICE departure areas
! 7. Direct fix to departure areas and area-coord remap;
!    use SLICE itself to iteratively estimate areas.
! (This cleaned up code has only these two alternatives.)

INTEGER, PARAMETER :: areafix = 7

END MODULE version

! ======================================================

MODULE grid

! Information about the grid size and resolution

! Set ng = p - 2 with the following
 INTEGER, PARAMETER :: p = 8, nx = 2**p, ny = 2**(p-1)
! Or set ng = p - 1 with the following
! INTEGER, PARAMETER :: p = 3, nx = 3*2**p, ny = 3*2**(p-1)
!INTEGER, PARAMETER :: p = 6, nx = 5*2**p, ny = 5*2**(p-1)

REAL*8, PARAMETER :: pi = 3.141592653589793D0, twopi = 2.0*pi, &
                     piby2 = pi/2.0, dx = 2.0*pi/nx, dy=pi/ny

REAL*8 :: rotgrid = 0.0d0*piby2

REAL*8 :: xp(nx), yp(ny), xu(nx), yu(ny), xv(nx), yv(ny+1), &
          cosp(ny), cosv(ny+1), sinp(ny), sinv(ny+1), area(ny), &
          dareadx, areatot, singeolatp(nx,ny), singeolatz(nx,ny+1), &
          geolatp(nx,ny), geolonp(nx,ny)



END MODULE grid

! ======================================================

MODULE constants

! Physical constants:

! Earth's radius
REAL*8 :: rearth

! Twice Earth's rotation rate
REAL*8 :: twoomega

! Gravity
REAL*8 :: gravity

! Numerical constants: Reference geopotential
REAL*8 :: phiref


END MODULE constants

! ======================================================

MODULE force

! Fields used for forced integrations, such as equilibrium
! phi fields for forced dissipative Galewsky

USE grid

REAL*8 :: phieqm(nx,ny), phieqmmean, &
          ueqm(nx,ny), veqm(nx,ny+1)
REAL*8 :: tauh = 100.0D0*86400.0D0, &
          tauv = 100.0D0*86400.0D0

! For physics-dynamics coupling
! lpdc          Logical switch to switch physics on or off
! taum          Timescale for relaxation of phimem to phi
! taul          Timescale for relaxation of latent energy
! phileq        Equilibrium value of latent energy
! philrecharge  Set phil to this value when convection completes
! thrconv       Threshold value of phi - phimem for convection
! tauconv       Timescale for convection
! betal         Coefficient for depletion of latent energy
LOGICAL :: lpdc = .false.
REAL*8, PARAMETER :: taum = 10.0D0*86400.0D0,    &
                     taul =  5.0D0*86400.0D0,    &
                     phileq = 4.0D0,             &
                     philrecharge = -1.0D0,      &
                     thrconv = 350.0D0,          &
                     betal = 0.1D0


END MODULE force

! ======================================================

MODULE thresholds

! Temporary module to compute diagnostics on phi - phimem

INTEGER :: tcount(10)
REAL*8 :: tvalues(10) = (/325.0D0, 350.0D0, 375.0D0, 400.0D0, &
                          425.0D0, 450.0D0, 475.0D0, 500.0D0, &
                          525.0D0, 550.0D0/), &
          tarea(10)


END MODULE thresholds

! ======================================================

MODULE state

! The state of the model's prognostic variables
! at the current and new time levels

USE grid

REAL*8 :: phi0(nx,ny), u0(nx,ny), v0(nx,ny+1), &
          phi(nx,ny), u(nx,ny), v(nx,ny+1), &
	  phis(nx,ny)
REAL*8 :: phimem0(nx,ny), philat0(nx,ny), &
          phimem(nx,ny), philat(nx,ny)


END MODULE state

! ======================================================

MODULE accumdiag

! Fields needed for accumulated diagnostics

USE grid

INTEGER :: naccum

REAL*8 :: accumphi(nx,ny), accumdiv(nx,ny), accumxi(nx, ny+1), &
          accumphi2(nx,ny), accumdiv2(nx,ny), accumxi2(nx,ny+1), &
          accumcon(nx,ny)

END MODULE accumdiag

! ======================================================

MODULE timeinfo

! Information related to run length, time step, etc

INTEGER :: nstop, idump, istep

REAL*8 :: tstop, dt, hdt, adt, bdt


END MODULE timeinfo

! ======================================================

MODULE util

! Useful function

CONTAINS

FUNCTION near(dx,domain)

! To correct for wrap-around when evaulating distances
! between points

IMPLICIT NONE
REAL*8, INTENT(IN) :: dx, domain
! REAL*8, INTENT(OUT) :: near
REAL*8 :: near
REAL*8 :: hdom

hdom = 0.5*domain
near = dx
IF (dx .ge. hdom) near = near - domain
IF (dx .le. -hdom) near = near + domain

END FUNCTION near


END MODULE util

! ====================================================

MODULE errdiag

! Save initial state for use in computing error diagnostics
! for steady state test cases like TC2

USE grid

INTEGER :: ccount = 0
REAL*8 :: phi_init(nx,ny), u_init(nx,ny), v_init(nx,ny+1)


END MODULE errdiag

! ====================================================

MODULE work

! Variables used in intermediate calculations

USE grid

! Departure point coordinates for u-points and v-points and phi-points
REAL*8 :: xdu(nx,ny), ydu(nx,ny), xdv(nx,ny+1), ydv(nx,ny+1), &
          xdp(nx,ny), ydp(nx,ny), aread(nx,ny)
	  
! Modified winds, divergence and departure points for SLICE
INTEGER :: jmods, jmodn
REAL*8 :: umod(nx,ny), vmod(nx,ny+1), u0mod(nx,ny), v0mod(nx,ny+1), &
          umodbar(nx,ny+1), vmodbar(nx,ny), u0modbar(nx,ny+1), v0modbar(nx,ny), &
          divmod(nx,ny), div0mod(nx,ny), &
          xdumod(nx,ny), ydumod(nx,ny), xdvmod(nx,ny+1), ydvmod(nx,ny+1)

! u at v-points and v at u-points
REAL*8 :: u0bar(nx,ny+1), v0bar(nx,ny), ubar(nx,ny+1), vbar(nx,ny)

! u and v at phi points
REAL*8 :: u0p(nx,ny), v0p(nx,ny), up(nx,ny), vp(nx,ny)

! Coriolis terms: fu at v-points and fv at u-points
REAL*8 :: fu0(nx,ny+1), fv0(nx,ny), fu(nx,ny+1), fv(nx,ny)

! Terms in momentum equations at u and v points
! and departure points at current time level
REAL*8 :: ru0(nx,ny), rv0(nx,ny+1), ru0bar(nx,ny+1), rv0bar(nx,ny), &
          rud(nx,ny), rvd(nx,ny+1)

! Terms in semi-Lagrangian form of phi equation at phi points
! and departure points at current time level
REAL*8 :: rphi0(nx,ny), rphid(nx,ny)

! Derivative of density, upon arrival at arrival cell,
! wrt departure point locations
REAL*8 :: dmdxu(nx,ny), dmdyu(nx,ny), dmdxv(nx,ny), dmdyv(nx,ny)

! Gradients of current-time-level terms in momentum equations
REAL*8 :: drudx(nx,ny), drudy(nx,ny), drvdx(nx,ny), drvdy(nx,ny)

! Gradient of orography
REAL*8 :: dphisdx(nx,ny), dphisdy(nx,ny+1)

! Coefficients for elliptic problem, and right hand side
REAL*8 :: rphi(nx,ny), rhs(nx,ny), ru(nx,ny), rv(nx,ny+1)

! Divergence at current time level
REAL*8 :: div0(nx,ny)

! Data used to diagnose convergence of iterations
REAL*8 :: oldu(nx,ny), oldv(nx,ny+1), oldphi(nx,ny), &
          deltau(nx,ny), deltav(nx,ny+1), deltaphi(nx,ny), &
	  olddiv(nx,ny), oldphid(nx,ny)
	  
! Modification to phi to represent convection
REAL*8 :: phiphys(nx,ny)


END MODULE work

! ======================================================

MODULE contest

! Variables used to diagnose convergence of iterations
INTEGER :: itcount
REAL*8:: rmsu(20), rmsv(20), rmsphi(20)

END MODULE contest

! ======================================================

MODULE alldata

! Collect all variables in one module

USE state
USE force
USE constants
USE timeinfo
USE work
USE version

END MODULE alldata

! =======================================================

PROGRAM sw

! John Thuburn 18/01/08
! Code to solve spherical
! shallow water equations using a fully
! implicit semi-Lagrangian scheme (including
! trajectories) with conservative (SLICE)
! treatment of phi.
!

IMPLICIT NONE

! ----------------------------------------------------

! PRINT *,'TRY INCREMENTAL NL SOLVER'
! print *,'Try interpolating without averaging first'

! Open dump file for matlab plotting
OPEN(34,FILE='swdump.m')

! Set up grid
CALL setupgrid
print *,'Grid set up'

! Output grid for use in tabulating reference solutions
! CALL dumpgrid
! STOP

! Set up physical constants
CALL setconst
print *,'Constants set up'

! Set up time step and run length
CALL timing
print *,'Timing info set up'

! Initialize accumulated diagnostics
CALL iniaccum
print *,'Accumulated diagnostics initialized'

! Set up initial conditions and orography
CALL initial
print *,'Initial data set up'

! Integrate
CALL integrate

! Output convergence test diagnostics
! CALL outcon

print *,'Finished'

! Close off dump files
WRITE(33,*) -1, -1
CLOSE(33)
CLOSE(34)

! -----------------------------------------------------

END PROGRAM sw

! =====================================================

SUBROUTINE setupgrid

USE grid

IMPLICIT NONE

INTEGER :: i, j, hny
REAL*8 :: sinr, cosr, sinlat, num, den

! Set up grid information

hny = ny/2

DO i = 1, nx
  xp(i) = (i-0.5)*dx
  xu(i) = (i-1)*dx
  xv(i) = (i-0.5)*dx
ENDDO

DO j = 1, ny
  yp(j) = (j-hny-0.5)*dy
  yu(j) = (j-hny-0.5)*dy
  yv(j) = (j-hny-1)*dy
  cosp(j) = COS(yp(j))
  cosv(j) = COS(yv(j))
  sinp(j) = SIN(yp(j))
  sinv(j) = SIN(yv(j))
ENDDO

yv(1) = -piby2
yv(ny+1) = piby2
cosv(1) = 0.0
cosv(ny+1) = 0.0
sinv(1) = -1.0
sinv(ny+1) = 1.0

dareadx = 0.0
areatot = 0.0
DO j = 1, ny
  ! area(j) = dx*(sinv(j+1) - sinv(j))
  area(j) = dx*dy*cosp(j)
  areatot = areatot + area(j)
  dareadx = dareadx + dy*cosp(j)
ENDDO
areatot = areatot*nx

! print *,'dareadx = ',dareadx,' areatot = ',areatot

! Sine of geographical latitude for rotated model grid,
! and geographical latitude and longitude
sinr = SIN(rotgrid)
cosr = COS(rotgrid)
! Phi points
DO i = 1, nx
  DO j = 1, ny
    sinlat = cosr*sinp(j) - sinr*cosp(j)*SIN(xp(i))
    singeolatp(i,j) = sinlat
    geolatp(i,j) = ASIN(sinlat)
    num = cosr*cosp(j)*SIN(xp(i)) + sinr*sinp(j)
    den = cosp(j)*COS(xp(i))
    geolonp(i,j) = MODULO(ATAN2(num,den),twopi)
  ENDDO
ENDDO
! Vorticity points
DO i = 1, nx
  DO j = 1, ny+1
    sinlat = cosr*sinv(j) - sinr*cosv(j)*SIN(xu(i))
    singeolatz(i,j) = sinlat
  ENDDO
ENDDO



END SUBROUTINE setupgrid

! =====================================================

SUBROUTINE setconst

USE constants

IMPLICIT NONE

! Set Earth's radius
rearth = 6371220.0d0

! Set rotation rate
twoomega = 1.4584d-4
! twoomega = 0.0

! Gravity
gravity = 9.80616D0

! Set reference geopotential for phi
! phiref = 2.3d4    ! mean for case 2
! phiref = 37979.96 ! min for case 5
! phiref = 58439.00 ! max for case 5
! phiref = 48209.48 ! minimax for case 5
! phiref = 78464.01 ! min for case 6
! phiref = 103502.53! max for case 6
! phiref = 90983.27 ! minimax for case 6
! phiref = 1.0d5    ! suitable for Galewsky
! phiref = 1.0d4    ! for Held-Suarez-like case
phiref = 1.0d+2    ! for Hollingswoth test

END SUBROUTINE setconst

! =====================================================

SUBROUTINE iniaccum

! Initialize accumulated diagnostics

USE accumdiag

IMPLICIT NONE


naccum = 0
accumphi = 0.0d0
accumdiv = 0.0d0
accumxi = 0.0d0
accumphi2 = 0.0d0
accumdiv2 = 0.0d0
accumxi2 = 0.0d0
accumcon = 0.0d0


END SUBROUTINE iniaccum

! =====================================================

SUBROUTINE initial

! Set up initial conditions

USE state
USE force
USE constants
USE timeinfo
use work ! for testing advection
USE errdiag

IMPLICIT NONE

INTEGER, PARAMETER :: nygg = 2*nx

INTEGER :: i, j, rrw, ip, jp, jy, ic
REAL*8 :: gwamp = 100.0d2, phimeangw = 1.0d5
REAL*8 :: u00, phi00, sinr, cosr
REAL*8 :: phis0, rr0, rr, latc, longc
REAL*8 :: wrw, krw, phi0rw, a, b, c
REAL*8 :: psigg(nygg+1), hgg(nygg+1), u1, u2, l1, l2, &
          lat0, lat1, umax, en, umen, dygg, psi1, psi2, &
          alpha, beta, totvol, totarea, hbar, den, num, hpert, &
          lat2, e1, e2, coslat, c1, c2, c3, c4


! ic = 1   Resting, constant phi
! ic = 2   Williamson et al. test case 2
! ic = 3   Williamson et al. test case 5
! ic = 5   Galewsky et al. barotropically unstable jet - specify stream fn
! ic = 6   SW Held-Suarez - specify stream function and state to relax to
! ic = 7   Rossby Haurwitz wave
! ic = 8   Like case 2 but with orography balancing u
ic = 8


! Departure point and advection tests

!DO j = 1, ny
!  DO i = 1, nx
!    if (abs(xu(i)-0.3) .lt. 0.3) Then
!      phi0(i,j) = cos(yu(j))*(cos(0.5*pi*(xu(i)-0.3)/0.3))**2
!    else
!      phi0(i,j) = 0.0
!    endif
!    phi0(i,j) = sin(yp(j))
!    u0(i,j) = 100.0*cosp(j)*(sin(5.0*yp(j))**2)
!    phi0(i,j) = 1000.0
!    u0(i,j) = (rearth*dx/dt)*cosp(j)*(-1)**(i+j)
!  ENDDO
!ENDDO

!DO j = 1, ny+1
!  DO i = 1, nx 
!    if (abs(xv(i)-0.3) .lt. 0.3) Then
!      rv0(i,j) = cos(yv(j))*(cos(0.5*pi*(xv(i)-0.3)/0.3))**2
!    else
!      rv0(i,j) = 0.0
!    endif
!    v0(i,j) = - 000.0*SIN(xv(i)) &
!              + 000.0*COS(xv(i))
!    v0(i,j) = cosv(j)*rearth*dy/(2.0*dt)
!     v0(i,j) = 0.0
!  ENDDO
!ENDDO


! Global non-rotating gravity wave

!DO j = 1, ny
!  DO i = 1, nx
!     phi0(i,j) = phimeangw + gwamp*COS(yp(j))*COS(xp(i))   ! Case 1
!      phi0(i,j) = phimeangw + gwamp*COS(yp(j))*SIN(xp(i))  ! Case 2
!      u0(i,j) = gwamp*COS(xu(i))/(SQRT(2.0*phimeangw))     ! Case 1
!     u0(i,j) = 0.0                                         ! Case 2
!  ENDDO
!ENDDO    

!DO j = 1, ny+1
!  DO i = 1, nx 
!    v0(i,j) = -gwamp*SIN(xv(i))*SIN(yv(j))/(SQRT(2.0*phimeangw)) ! Case 1
!     v0(i,j) = gwamp*COS(yv(j))/(SQRT(2.0*phimeangw))            ! Case 2
!  ENDDO
!ENDDO


! Zero orography - modify below for TC 5 and Held-Suarez
phis = 0.0d0


IF (ic == 1) THEN

  ! Resting state with uniform phi
  phi0 = 1.0d5
  u0 = 0.0d0
  v0 = 0.0d0

ELSEIF (ic == 2 .OR. ic == 3) THEN

    ! Balanced solid body rotation
    IF (ic == 2) THEN
    ! Test case 2
    u00 = 2.0d0*pi*rearth/(12.0d0*86400.0d0)
    phi00 = 2.94d4
  ELSEIF (ic == 3) THEN
    ! Test case 5
    u00 = 20.0d0
    phi00 = 5960.0d0*gravity
  ENDIF

  cosr = COS(rotgrid)
  sinr = SIN(rotgrid)

  DO j = 1, ny
    DO i = 1, nx
       phi0(i,j) = phi00 - 0.5*(rearth*twoomega*u00 + u00*u00)*singeolatp(i,j)**2
       u0(i,j) = u00*(cosp(j)*cosr + SIN(xu(i))*sinp(j)*sinr)
    ENDDO
  ENDDO    

  DO j = 1, ny+1
    DO i = 1, nx 
      v0(i,j) = u00*COS(xv(i))*sinr
    ENDDO
  ENDDO

  IF (ic == 3) THEN
    ! Include mountain
    phis0 = 2000.0d0*9.80616d0
    rr0 = pi/9.0d0
    latc = pi/6.0d0
    longc = 0.5d0*pi
    DO j = 1, ny
      DO i = 1, nx  
        ! Isolated mountain
        rr = SQRT(MIN(rr0*rr0,(geolonp(i,j)-longc)**2 + (geolatp(i,j)-latc)**2))
        phis(i,j) = phis0*(1.0 - rr/rr0)   
      ENDDO
    ENDDO
    ! Correct phi to allow for orography
    phi0 = phi0 - phis
  ENDIF

ELSEIF (ic == 8) THEN

  u00 = 2.0d0*pi*rearth/(12.0d0*86400.0d0)
  !phi00 = 1.0d-2
  phi00 = 1.0d+2

  cosr = COS(rotgrid)
  sinr = SIN(rotgrid)

  DO j = 1, ny
    DO i = 1, nx
       phis(i,j) = - 0.5*(rearth*twoomega*u00 + u00*u00)*singeolatp(i,j)**2
       phi0(i,j) = phi00
       u0(i,j) = u00*(cosp(j)*cosr + SIN(xu(i))*sinp(j)*sinr)
    ENDDO
  ENDDO    

  DO j = 1, ny+1
    DO i = 1, nx 
      v0(i,j) = u00*COS(xv(i))*sinr
    ENDDO
  ENDDO

ELSEIF (ic == 5) THEN

  ! Galewsky test initialized using stream function

  umax = 80.0D0
  lat0 = pi/7.0D0
  lat1 = pi/2.0D0 - lat0
  en = exp(-4.0D0/(lat1 - lat0)**2)
  umen = umax/en
  totvol = 0.0D0
  totarea = 0.0D0
  ! Integrate to tabulate h and psi as functions of geographical
  ! latitude
  dygg = pi/nygg
  hgg(1) = 0.0D0
  psigg(1) = 0.0D0
  DO j = 2, nygg
    l1 = (j-2)*dygg - piby2
    den = (l1 - lat0)*(l1 - lat1)
    IF (den .lt. 0.0D0) THEN
      u1 = umen*exp(1.0D0/den)
    ELSE
      u1 = 0.0D0
    ENDIF
    l2 = (j-1)*dygg - piby2
    den = (l2 - lat0)*(l2 - lat1)
    IF (den .lt. 0.0D0) THEN
      u2 = umen*exp(1.0D0/den)
    ELSE
      u2 = 0.0D0
    ENDIF
    psigg(j) = psigg(j-1) - 0.5*(u1 + u2)*dygg
    u1 = u1*(twoomega*SIN(l1) + TAN(l1)*u1/rearth)
    u2 = u2*(twoomega*SIN(l2) + TAN(l2)*u2/rearth)
    hgg(j) = hgg(j-1) - rearth*0.5*(u1 + u2)*dygg
    totarea = totarea + cos(l2)*dygg
    totvol = totvol + hgg(j)*cos(l2)*dygg
  ENDDO
  psigg(nygg+1) = psigg(nygg)
  hgg(nygg+1) = hgg(nygg)
  totvol = totvol/(totarea*gravity)
  hgg = hgg + (1.0D4 - totvol)*gravity

  ! Now assign h as a function of geographical latitude
  ! using interpolation from tabulated values
  totvol = 0.0D0
  totarea = 0.0D0
  DO j = 1, ny
    DO i = 1, nx
      l1 = ASIN(singeolatp(i,j)) + piby2
      jy = FLOOR(l1/dygg) + 1
      beta = (l1 - (jy - 1)*dygg)/dygg
      IF (jy == 1 .OR. jy == nygg) THEN
        ! Linear interpolation
        c2 = 1.0D0 - beta
        c3 = beta
        phi0(i,j) = c2*hgg(jy) + c3*hgg(jy+1)
      ELSE
        ! Cubic interpolation
        c1 = -beta*(beta - 1.0D0)*(beta - 2.0D0)/6.0D0
        c2 = 0.5D0*(beta + 1.0D0)*(beta - 1.0D0)*(beta - 2.0D0)
        c3 = -0.5D0*(beta + 1.0D0)*beta*(beta - 2.0D0)
        c4 = (beta + 1.0D0)*beta*(beta - 1.0D0)/6.0D0
        phi0(i,j) = c1*hgg(jy-1) + c2*hgg(jy) + c3*hgg(jy+1) + c4*hgg(jy+2)
      ENDIF
      totarea = totarea + area(j)
      totvol = totvol + area(j)*phi0(i,j)
    ENDDO
  ENDDO
  ! Now calculate velocity components by interpolating
  ! stream function to each vorticity point
  ! u field
  DO j = 1, ny
    jp = j + 1
    DO i = 1, nx
      l1 = ASIN(singeolatz(i,j)) + piby2
      jy = FLOOR(l1/dygg) + 1
      beta = (l1 - (jy - 1)*dygg)/dygg
      IF (jy == 1 .OR. jy == nygg) THEN
        ! Linear interpolation
        c2 = 1.0D0 - beta
        c3 = beta
        psi1 = c2*psigg(jy) + c3*psigg(jy+1)
      ELSE
        ! Cubic interpolation
        c1 = -beta*(beta - 1.0D0)*(beta - 2.0D0)/6.0D0
        c2 = 0.5D0*(beta + 1.0D0)*(beta - 1.0D0)*(beta - 2.0D0)
        c3 = -0.5D0*(beta + 1.0D0)*beta*(beta - 2.0D0)
        c4 = (beta + 1.0D0)*beta*(beta - 1.0D0)/6.0D0
        psi1 = c1*psigg(jy-1) + c2*psigg(jy) + c3*psigg(jy+1) + c4*psigg(jy+2)
      ENDIF
      l2 = ASIN(singeolatz(i,jp)) + piby2
      jy = FLOOR(l2/dygg) + 1
      beta = (l2 - (jy - 1)*dygg)/dygg
      IF (jy == 1 .OR. jy == nygg) THEN
        ! Linear interpolation
        c2 = 1.0D0 - beta
        c3 = beta
        psi2 = c2*psigg(jy) + c3*psigg(jy+1)
      ELSE
        ! Cubic interpolation
        c1 = -beta*(beta - 1.0D0)*(beta - 2.0D0)/6.0D0
        c2 = 0.5D0*(beta + 1.0D0)*(beta - 1.0D0)*(beta - 2.0D0)
        c3 = -0.5D0*(beta + 1.0D0)*beta*(beta - 2.0D0)
        c4 = (beta + 1.0D0)*beta*(beta - 1.0D0)/6.0D0
        psi2 = c1*psigg(jy-1) + c2*psigg(jy) + c3*psigg(jy+1) + c4*psigg(jy+2)
      ENDIF
      u0(i,j) = -(psi2 - psi1)/dy
    ENDDO
  ENDDO
  ! v field
  DO j = 2, ny
    DO i = 1,nx
      ip = i + 1
      IF (i .eq. nx) ip = 1
      l1 = ASIN(singeolatz(i,j)) + piby2
      jy = FLOOR(l1/dygg) + 1
      beta = (l1 - (jy - 1)*dygg)/dygg
      IF (jy == 1 .OR. jy == nygg) THEN
        ! Linear interpolation
        c2 = 1.0D0 - beta
        c3 = beta
        psi1 = c2*psigg(jy) + c3*psigg(jy+1)
      ELSE
        ! Cubic interpolation
        c1 = -beta*(beta - 1.0D0)*(beta - 2.0D0)/6.0D0
        c2 = 0.5D0*(beta + 1.0D0)*(beta - 1.0D0)*(beta - 2.0D0)
        c3 = -0.5D0*(beta + 1.0D0)*beta*(beta - 2.0D0)
        c4 = (beta + 1.0D0)*beta*(beta - 1.0D0)/6.0D0
        psi1 = c1*psigg(jy-1) + c2*psigg(jy) + c3*psigg(jy+1) + c4*psigg(jy+2)
      ENDIF
      l2 = ASIN(singeolatz(ip,j)) + piby2
      jy = FLOOR(l2/dygg) + 1
      beta = (l2 - (jy - 1)*dygg)/dygg
      IF (jy == 1 .OR. jy == nygg) THEN
        ! Linear interpolation
        c2 = 1.0D0 - beta
        c3 = beta
        psi2 = c2*psigg(jy) + c3*psigg(jy+1)
      ELSE
        ! Cubic interpolation
        c1 = -beta*(beta - 1.0D0)*(beta - 2.0D0)/6.0D0
        c2 = 0.5D0*(beta + 1.0D0)*(beta - 1.0D0)*(beta - 2.0D0)
        c3 = -0.5D0*(beta + 1.0D0)*beta*(beta - 2.0D0)
        c4 = (beta + 1.0D0)*beta*(beta - 1.0D0)/6.0D0
        psi2 = c1*psigg(jy-1) + c2*psigg(jy) + c3*psigg(jy+1) + c4*psigg(jy+2)
      ENDIF
      v0(i,j) = (psi2 - psi1)/(dx*cosv(j))
    ENDDO
  ENDDO
  v0(:,1) = 0.0
  v0(:,ny+1) = 0.0

  ! Geopotential perturbation
  alpha = 1.0D0/3.0D0
  beta = 1.0D0/15.0D0
  hpert = 120.0D0
  ! hpert = 0.0d0
  ! print *,'*** No Galewsky perturbation ***'
  lat2 = 0.5D0*piby2
  DO j = 1, ny
    DO i = 1, nx
      l2 = geolatp(i,j)
      coslat = COS(l2)
      l1 = geolonp(i,j)
!      l1 = geolonp(i,j) - 0.5d0*pi ! Shift pertn for ease of plotting
      IF (l1 > pi) l1 = l1 - twopi
      IF (l1 < -pi) l1 = l1 + twopi
      e1 = EXP(-(l1/alpha)**2)
      e2 = EXP(-((lat2 - l2)/beta)**2)
      phi0(i,j) = phi0(i,j) + gravity*hpert*coslat*e1*e2
    ENDDO
  ENDDO

ELSEIF (ic == 6) THEN

  ! Held-Suarez-like test initialized using stream function

  ! Equator to pole orography difference
!  phis0 = 1.4d4
!  DO j = 1, ny
!    DO i = 1, nx
!      phis(i,j) = phis0*cos(geolatp(i,j))**2
!    ENDDO
!  ENDDO

  umax = 120.0D0
  num = 12.0D0
  totvol = 0.0D0
  totarea = 0.0D0
  ! Integrate to tabulate h and psi as functions of geographical
  ! latitude
  dygg = pi/nygg
  hgg(1) = 0.0D0
  psigg(1) = 0.0D0
  DO j = 2, nygg
    l1 = (j-2)*dygg - piby2
    u1 = umax*(SIN(2.0D0*l1)**2)*((SIN(num*l1))**2 - 0.5D0)
    l2 = (j-1)*dygg - piby2
    u2 = umax*(SIN(2.0D0*l2)**2)*((SIN(num*l2))**2 - 0.5D0)
    psigg(j) = psigg(j-1) - 0.5*(u1 + u2)*dygg
    u1 = u1*(twoomega*SIN(l1) + TAN(l1)*u1/rearth)
    u2 = u2*(twoomega*SIN(l2) + TAN(l2)*u2/rearth)
    hgg(j) = hgg(j-1) - rearth*0.5*(u1 + u2)*dygg
    totarea = totarea + cos(l2)*dygg
    totvol = totvol + hgg(j)*cos(l2)*dygg
  ENDDO
  psigg(nygg+1) = psigg(nygg)
  hgg(nygg+1) = hgg(nygg)
!  totvol = totvol/(totarea*gravity)
!  hgg = hgg + (1.0D3 - totvol)*gravity

  ! Now assign h as a function of geographical latitude
  ! using interpolation from tabulated values
  DO j = 1, ny
    DO i = 1, nx
      l1 = ASIN(singeolatp(i,j)) + piby2
      jy = FLOOR(l1/dygg) + 1
      beta = (l1 - (jy - 1)*dygg)/dygg
      IF (jy == 1 .OR. jy == nygg) THEN
        ! Linear interpolation
        c2 = 1.0D0 - beta
        c3 = beta
        phi0(i,j) = c2*hgg(jy) + c3*hgg(jy+1)
      ELSE
        ! Cubic interpolation
        c1 = -beta*(beta - 1.0D0)*(beta - 2.0D0)/6.0D0
        c2 = 0.5D0*(beta + 1.0D0)*(beta - 1.0D0)*(beta - 2.0D0)
        c3 = -0.5D0*(beta + 1.0D0)*beta*(beta - 2.0D0)
        c4 = (beta + 1.0D0)*beta*(beta - 1.0D0)/6.0D0
        phi0(i,j) = c1*hgg(jy-1) + c2*hgg(jy) + c3*hgg(jy+1) + c4*hgg(jy+2)
      ENDIF
      ! Overwrite with smooth pole-equator gradient
      l1 = ASIN(singeolatp(i,j))
      phi0(i,j) = 8000.0d0*(1.0d0 - sin(l1)**2)
    ENDDO
  ENDDO
  ! Correct for orography
  phi0 = phi0 - phis
  ! Correct to obtain desired mean depth
  totvol = 0.0D0
  totarea = 0.0D0
  DO j = 1, ny
    DO i = 1, nx
      totarea = totarea + area(j)
      totvol = totvol + area(j)*phi0(i,j)
    ENDDO
  ENDDO
  totvol = totvol/(totarea*gravity)
  phi0 = phi0 + (1.0D3 - totvol)*gravity  
  ! Now calculate velocity components by interpolating
  ! stream function to each vorticity point
  ! u field
  DO j = 1, ny
    jp = j + 1
    DO i = 1, nx
      l1 = ASIN(singeolatz(i,j)) + piby2
      jy = FLOOR(l1/dygg) + 1
      beta = (l1 - (jy - 1)*dygg)/dygg
      IF (jy == 1 .OR. jy == nygg) THEN
        ! Linear interpolation
        c2 = 1.0D0 - beta
        c3 = beta
        psi1 = c2*psigg(jy) + c3*psigg(jy+1)
      ELSE
        ! Cubic interpolation
        c1 = -beta*(beta - 1.0D0)*(beta - 2.0D0)/6.0D0
        c2 = 0.5D0*(beta + 1.0D0)*(beta - 1.0D0)*(beta - 2.0D0)
        c3 = -0.5D0*(beta + 1.0D0)*beta*(beta - 2.0D0)
        c4 = (beta + 1.0D0)*beta*(beta - 1.0D0)/6.0D0
        psi1 = c1*psigg(jy-1) + c2*psigg(jy) + c3*psigg(jy+1) + c4*psigg(jy+2)
      ENDIF
      l2 = ASIN(singeolatz(i,jp)) + piby2
      jy = FLOOR(l2/dygg) + 1
      beta = (l2 - (jy - 1)*dygg)/dygg
      IF (jy == 1 .OR. jy == nygg) THEN
        ! Linear interpolation
        c2 = 1.0D0 - beta
        c3 = beta
        psi2 = c2*psigg(jy) + c3*psigg(jy+1)
      ELSE
        ! Cubic interpolation
        c1 = -beta*(beta - 1.0D0)*(beta - 2.0D0)/6.0D0
        c2 = 0.5D0*(beta + 1.0D0)*(beta - 1.0D0)*(beta - 2.0D0)
        c3 = -0.5D0*(beta + 1.0D0)*beta*(beta - 2.0D0)
        c4 = (beta + 1.0D0)*beta*(beta - 1.0D0)/6.0D0
        psi2 = c1*psigg(jy-1) + c2*psigg(jy) + c3*psigg(jy+1) + c4*psigg(jy+2)
      ENDIF
      u0(i,j) = -(psi2 - psi1)/dy
    ENDDO
  ENDDO
  ! v field
  DO j = 2, ny
    DO i = 1,nx
      ip = i + 1
      IF (i .eq. nx) ip = 1
      l1 = ASIN(singeolatz(i,j)) + piby2
      jy = FLOOR(l1/dygg) + 1
      beta = (l1 - (jy - 1)*dygg)/dygg
      IF (jy == 1 .OR. jy == nygg) THEN
        ! Linear interpolation
        c2 = 1.0D0 - beta
        c3 = beta
        psi1 = c2*psigg(jy) + c3*psigg(jy+1)
      ELSE
        ! Cubic interpolation
        c1 = -beta*(beta - 1.0D0)*(beta - 2.0D0)/6.0D0
        c2 = 0.5D0*(beta + 1.0D0)*(beta - 1.0D0)*(beta - 2.0D0)
        c3 = -0.5D0*(beta + 1.0D0)*beta*(beta - 2.0D0)
        c4 = (beta + 1.0D0)*beta*(beta - 1.0D0)/6.0D0
        psi1 = c1*psigg(jy-1) + c2*psigg(jy) + c3*psigg(jy+1) + c4*psigg(jy+2)
      ENDIF
      l2 = ASIN(singeolatz(ip,j)) + piby2
      jy = FLOOR(l2/dygg) + 1
      beta = (l2 - (jy - 1)*dygg)/dygg
      IF (jy == 1 .OR. jy == nygg) THEN
        ! Linear interpolation
        c2 = 1.0D0 - beta
        c3 = beta
        psi2 = c2*psigg(jy) + c3*psigg(jy+1)
      ELSE
        ! Cubic interpolation
        c1 = -beta*(beta - 1.0D0)*(beta - 2.0D0)/6.0D0
        c2 = 0.5D0*(beta + 1.0D0)*(beta - 1.0D0)*(beta - 2.0D0)
        c3 = -0.5D0*(beta + 1.0D0)*beta*(beta - 2.0D0)
        c4 = (beta + 1.0D0)*beta*(beta - 1.0D0)/6.0D0
        psi2 = c1*psigg(jy-1) + c2*psigg(jy) + c3*psigg(jy+1) + c4*psigg(jy+2)
      ENDIF
      v0(i,j) = (psi2 - psi1)/(dx*cosv(j))
    ENDDO
  ENDDO
  v0(:,1) = 0.0
  v0(:,ny+1) = 0.0

  ! Equilibrium phi for forced dissipative case
  phieqmmean = 1.0D3*gravity
  phieqm = 1.0d0*(phi0 - phieqmmean) + phieqmmean
  ueqm   = 1.0d0*u0
  veqm   = 1.0d0*v0

  ! Now reset initial state
  phi0 = 1.0D3*gravity
  u0 = 0.0D0
  v0 = 0.0D0

  ! And initialize physical fields
  phimem0 = phi0
  philat0 = phileq

  ! call dumpm(phieqm,'phieqm',nx,ny)

  ! Geopotential perturbation to kick off turbulence
  alpha = 1.0D0/3.0D0
  beta = 1.0D0/15.0D0
  hpert = 1.20D0
  ! hpert = 0.0d0
  ! print *,'*** No Held-Suarez perturbation ***'
  lat2 = 0.5D0*piby2
  DO j = 1, ny
    DO i = 1, nx
      l2 = geolatp(i,j)
      coslat = COS(l2)
      l1 = geolonp(i,j)
      IF (l1 > pi) l1 = l1 - twopi
      e1 = EXP(-(l1/alpha)**2)
      e2 = EXP(-((lat2 - l2)/beta)**2)
      phi0(i,j) = phi0(i,j) + gravity*hpert*coslat*e1*e2
    ENDDO
  ENDDO

ELSEIF (ic == 7) THEN

  ! Rossby Harwitz wave
  rrw = 4
  wrw = 1.0*7.848e-6
  ! wrw = twoomega/28.0  ! To obtain a stationary RH wave !
  krw = 1.0*7.848e-6
  phi0rw = 8.0e3*9.80616
  ! phi0rw = 8.0e5*9.80616 ! Large rossby radius ~ bve

  DO j = 1, ny
    a = 0.5*wrw*(twoomega + wrw)*cosp(j)*cosp(j) &
      + 0.25*(krw*krw*(cosp(j))**(2*rrw)) &
            *((rrw+1)*cosp(j)*cosp(j) + (2*rrw*rrw - rrw - 2) - &
	       2*rrw*rrw*(cosp(j))**(-2))
    b = (((twoomega + 2*wrw)*krw)/((rrw+1)*(rrw+2))) &
          *((cosp(j))**rrw) &
          *((rrw*rrw + 2*rrw + 2) - (rrw + 1)*(rrw + 1)*cosp(j)*cosp(j))
    c = 0.25*krw*krw*(cosp(j)**(2*rrw))*((rrw + 1)*cosp(j)*cosp(j) - (rrw+2))
    DO i = 1, nx
      phi0(i,j) = phi0rw &
                + rearth*rearth*( a + b*COS(rrw*xp(i)) + c*COS(2*rrw*xp(i)) )
      u0(i,j) = rearth*( wrw*cosp(j) &
                    +  krw*((cosp(j))**(rrw-1)) &
                    *(rrw*sinp(j)*sinp(j) - cosp(j)*cosp(j)) &
                       *COS(rrw*xu(i)) )
    ENDDO
  ENDDO    

  DO j = 1, ny+1
    DO i = 1, nx
      v0(i,j) = - rearth*krw*rrw*(cosv(j)**(rrw-1))*sinv(j)*SIN(rrw*xv(i))
    ENDDO
  ENDDO

ELSE

  PRINT *,'Initial condition ic = ',ic,' not coded.'
  STOP

ENDIF





! Save initial state for computing errors in steady test cases
phi_init = phi0
u_init = u0
v_init = v0


END SUBROUTINE initial

! =====================================================

SUBROUTINE timing

! Set up info on run length, time step, etc

USE timeinfo
USE version

IMPLICIT NONE

! Length of run
! tstop = 220*86400.0d0
! tstop = 12000.0
! tstop = 900.0*2740
tstop = 5*86400.0d0

! Size of time step
dt = 900.0d0
hdt = 0.5d0*dt

! Off-centring parameters
adt = 0.5d0*dt
bdt = dt - adt

IF ((ischeme == 3) .AND. (adt .NE. hdt)) THEN
  PRINT *,'*** Time scheme must be centred for SLICE ***'
  PRINT *,'*** Think carefully about coupling to departure '
  PRINT *,'point calculation before changing ***.'
  STOP
ENDIF

! Total number of steps to take
nstop = nint(tstop/dt)

! Number of steps between output dumps
idump = 480


END SUBROUTINE timing

! ======================================================

SUBROUTINE integrate

! Controlling subroutine for the time integration

USE timeinfo
USE alldata

IMPLICIT NONE

! -------------------------------------------------------

! Preliminary calculations
CALL prelim
print *,'Done preliminaries'

! Output initial data
CALL diagnostics(0)
CALL outstatem(0)
! CALL outthresh(0)
PRINT *,' '

! Loop over steps
DO istep = 1, nstop

  CALL step
  
  CALL diagnostics(istep)
  IF (istep*dt > 17280000.0D0) THEN
    CALL accumulate
  ENDIF

  IF (MOD(istep,idump) == 0) THEN
    CALL outstatem(istep)
    CALL outaccm(istep)
    ! CALL outthresh(istep)
  ENDIF

  ! Frequent dumps to make a movie!
!  if (istep.ge.19200 .or. MOD(istep,idump) == 0) then
!  if ((istep.ge.48000).and.(istep.le.51840).and.(mod(istep,1) == 0)) then
!    call outstatem(istep)
!  endif

  ! If required write out current phi field as reference solution
  ! for other models
  ! CALL writeref

  ! If required read in reference solution from high resolution
  ! run and compute and write errors
  ! CALL diffref

  IF (MOD(istep,1) == 0) THEN
    PRINT *,'Done step ',istep
    PRINT *,' '
  ENDIF
 
ENDDO
! call dumpm(phi+phis,'phitot',nx,ny)



END SUBROUTINE Integrate

! ======================================================

SUBROUTINE step

! Subroutine to take one time step

USE alldata
USE contest

IMPLICIT NONE

INTEGER :: nouter = 2, ninner = 2, iouter, iinner, ng, count

REAL*8 :: nu
REAL*8 :: div(nx,ny)

! Convergence diagnostics at final step
!if (istep == nstop) then
!  print *,'Doing convergence test'
!  nouter = 20
!  ninner = 20
!endif


! Constant for Helmholtz problem
nu = 1.0/(adt*adt*phiref)

! Terms in momentum equation at current time level
CALL momentum

! For semi-Lagrangian version:
! Terms in phi equation at current time level
IF (ischeme == 1 .OR. ischeme == 2 .OR. ischeme == 3) THEN
  CALL slphieqn
ENDIF

! Use current values of state variables as first guess
! for new value
phi = phi0 
u = u0
v = v0
!print *,'zero first guess'
!phi = 0.0
!u = 0.0
!v = 0.0

! C-grid average velocity components
CALL cgridave(nx,ny,u,v,ubar,vbar)
CALL uvatphi(nx,ny,u,v,up,vp)

count = 0

olddiv = 0.0
oldphid = 0.0

! Outer loop
DO iouter = 1, nouter

  itcount = iouter

  ! Save current values to examine convergence
  ! oldu = u
  ! oldv = v
  ! oldphi = phi

  ! Recalculate departure points
  CALL departure
  CALL departurep
  
  ! and modified departure points in polar regions
  CALL departuremod
  
  ! Recalculate departure point terms and RHS of
  ! Phi equation
  CALL up_outer
  
  ! Inner loop
  DO iinner = 1, ninner
  
    ! Save current values to examine convergence
    !oldu = u
    !oldv = v
    !oldphi = phi
  
    ! Update Coriolis terms and hence RHS of
    ! Helmholtz problem
    CALL up_inner
    
    ! Solve Helmholtz
    ! ng = p - 2
    ng = p - 1
    rhs = -nu*rhs
    CALL mgsolve(phi,rhs,nu,ng)
    
    ! Back substitute to find new u and v
    CALL backsub

    ! And compute corresponding c-grid average velocities
    CALL cgridave(nx,ny,u,v,ubar,vbar)
    CALL uvatphi(nx,ny,u,v,up,vp)
    
    ! Calculate modified velocities for polar regions
    CALL modifywind
    CALL cgridave(nx,ny,umod,vmod,umodbar,vmodbar)
    
    !deltau = u - oldu
    !deltav = v - oldv
    !deltaphi = phi - oldphi

    !count = count + 1
    !rmsu(count) = SQRT(SUM(deltau*deltau)/(nx*ny))
    !rmsv(count) = SQRT(SUM(deltav*deltav)/(nx*ny))
    !rmsphi(count) = SQRT(SUM(deltaphi*deltaphi)/(nx*ny))
    
    ! print *,'count = ',count
    ! print *,' diffs = ',rmsphi(count),rmsu(count),rmsv(count)
    
  ENDDO
   
  !deltau = u - oldu
  !deltav = v - oldv
  !deltaphi = phi - oldphi

  !if (istep == nstop) then
    ! call dump(deltaphi,'deltaphi',nx,ny)
    ! call dump(deltau,'deltau',nx,ny)
    ! call dump(deltav,'deltav',nx,ny+1)
  !endif

  !count = count + 1
  !rmsu(count) = SQRT(SUM(deltau*deltau)/(nx*ny))
  !rmsv(count) = SQRT(SUM(deltav*deltav)/(nx*ny))
  !rmsphi(count) = SQRT(SUM(deltaphi*deltaphi)/(nx*ny))
  
  !print *,'count = ',count
  !print *,' diffs = ',rmsphi(count),rmsu(count),rmsv(count)
  
ENDDO

! print *,'Could do a final phi update to ensure conservation'

IF (lpdc) THEN
  ! Advect physics tracers using latest available departure points
  CALL lagrangep(xdp,ydp,phimem0,phimem)
  CALL lagrangep(xdp,ydp,philat0,philat)
ENDIF

! Include forcing terms
! CALL forcing


phi0 = phi
u0 = u
v0 = v
u0bar = ubar
v0bar = vbar
u0p = up
v0p = vp
u0mod = umod
v0mod = vmod
u0modbar = umodbar
v0modbar = vmodbar
phimem0 = phimem
philat0 = philat


END SUBROUTINE step

! ======================================================

SUBROUTINE prelim

! Preliminary calculations before commencing time stepping.

USE alldata
USE contest

IMPLICIT NONE


! Gradient of orography
CALL gradphis

! All velocity components at all points
CALL cgridave(nx,ny,u0,v0,u0bar,v0bar)
CALL uvatphi(nx,ny,u0,v0,u0p,v0p)

! First guess for new time level fields
phi = phi0
u = u0
v = v0

! All velocity components at all points
CALL cgridave(nx,ny,u,v,ubar,vbar)
CALL uvatphi(nx,ny,u,v,up,vp)

! First guess for departure points at first step
! and departure areas
! (At subsequent steps the first guess is taken from
! the previous step).
CALL departurefg
CALL departurefgp

! Modified winds in polar cap region
CALL modifywind
u0mod = umod
v0mod = vmod
CALL cgridave(nx,ny,u0mod,v0mod,u0modbar,v0modbar)
CALL cgridave(nx,ny,umod,vmod,umodbar,vmodbar)

! First guess for modified departure points
CALL departuremodfg



END SUBROUTINE prelim

! ========================================================

SUBROUTINE departurefg

! First guess for departure point calculation.
! Use velocity at old time level.

! Also first guess for departure areas

USE alldata

IMPLICIT NONE
INTEGER :: i, j
REAL*8 :: sina, cosa, x, y, r, sind, dlambda

! u-point departure points
DO  j = 1, ny
  sina = sinp(j)
  cosa = cosp(j)
  DO i = 1, nx
    ! Displacement in local Cartesian system
    x = -u0(i,j)*dt
    y = -v0bar(i,j)*dt
    ! Project back to spherical coordinate system
    r = SQRT(x*x + y*y + rearth*rearth)
    sind = (y*cosa + rearth*sina)/r
    ydu(i,j) = ASIN(sind)
    dlambda = ATAN2(x,rearth*cosa - y*sina)
    xdu(i,j) = MODULO(xu(i) + dlambda, twopi)
  ENDDO
ENDDO

! v-point departure points
DO  j = 2, ny
  sina = sinv(j)
  cosa = cosv(j)
  DO i = 1, nx
    ! Displacement in local Cartesian system
    x = -u0bar(i,j)*dt
    y = -v0(i,j)*dt
    ! Project back to spherical coordinate system
    r = SQRT(x*x + y*y + rearth*rearth)
    sind = (y*cosa + rearth*sina)/r
    ydv(i,j) = ASIN(sind)
    dlambda = ATAN2(x,rearth*cosa - y*sina)
    xdv(i,j) = MODULO(xv(i) + dlambda, twopi)
  ENDDO
ENDDO

! Departure areas. First guess equals arrival areas
! DO j = 1, ny
!   aread(:,j) = area(j)
! ENDDO


END SUBROUTINE departurefg

! ========================================================

SUBROUTINE departuremodfg

! First guess for modified departure point calculation.
! Use velocity at old time level.

! Also first guess for departure areas

USE alldata

IMPLICIT NONE
INTEGER :: i, j
REAL*8 :: sina, cosa, x, y, r, sind, dlambda

! u-point departure points
DO  j = 1, ny
  sina = sinp(j)
  cosa = cosp(j)
  DO i = 1, nx
    ! Displacement in local Cartesian system
    x = -u0mod(i,j)*dt
    y = -v0modbar(i,j)*dt
    ! Project back to spherical coordinate system
    r = SQRT(x*x + y*y + rearth*rearth)
    sind = (y*cosa + rearth*sina)/r
    ydu(i,j) = ASIN(sind)
    dlambda = ATAN2(x,rearth*cosa - y*sina)
    xdumod(i,j) = MODULO(xu(i) + dlambda, twopi)
  ENDDO
ENDDO

! v-point departure points
DO  j = 2, ny
  sina = sinv(j)
  cosa = cosv(j)
  DO i = 1, nx
    ! Displacement in local Cartesian system
    x = -u0modbar(i,j)*dt
    y = -v0mod(i,j)*dt
    ! Project back to spherical coordinate system
    r = SQRT(x*x + y*y + rearth*rearth)
    sind = (y*cosa + rearth*sina)/r
    ydv(i,j) = ASIN(sind)
    dlambda = ATAN2(x,rearth*cosa - y*sina)
    xdvmod(i,j) = MODULO(xv(i) + dlambda, twopi)
  ENDDO
ENDDO

! Departure areas. First guess equals arrival areas
DO j = 1, ny
  aread(:,j) = area(j)
ENDDO


END SUBROUTINE departuremodfg

! ========================================================

SUBROUTINE cgridave(nx,ny,u,v,ubar,vbar)

! To average u to v points and v to u points
! on the C-grid

IMPLICIT NONE

INTEGER, INTENT(IN) :: nx, ny
REAL*8, INTENT(IN) :: u(nx,ny)
REAL*8 ,INTENT(INOUT) :: v(nx,ny+1)
REAL*8, INTENT(OUT) :: ubar(nx,ny+1), vbar(nx,ny)

INTEGER :: i, im, ip, j, jm, jp

! Deal with polar values first, as we need to reconstruct
! polar v values as well as u values
CALL polar(u,ubar(:,1),v(:,1),ubar(:,ny+1),v(:,ny+1))

! u at v points
DO j = 2, ny
  jm = j-1
  DO i = 1, nx
    ip = i+1
    IF (i == nx) ip = 1
    ubar(i,j) = 0.25*(u(i,jm)+u(i,j)+u(ip,jm)+u(ip,j))
  ENDDO
ENDDO

! v at u points
DO j = 1, ny
  jp = j+1
  DO i = 1, nx
    im = i-1
    IF (im == 0) im = nx
    vbar(i,j) = 0.25*(v(im,j)+v(im,jp)+v(i,j)+v(i,jp))
  ENDDO
ENDDO


END SUBROUTINE cgridave

! ========================================================

SUBROUTINE polar(u,usp,vsp,unp,vnp)

! Determine polar values of u and v from u at nearest
! u-latitude. Also used for RHS of momentum equation.

USE grid

IMPLICIT NONE

REAL*8, INTENT(IN) :: u(nx,ny)
REAL*8, INTENT(OUT) :: usp(nx), vsp(nx), unp(nx), vnp(nx)
REAL*8 :: a, b, lambda, vp

! South pole
a = SUM(u(:,1)*SIN(xu(:)))
b = SUM(u(:,1)*COS(xu(:)))
lambda = ATAN2(b,-a)
vp = (-a*COS(lambda) + b*SIN(lambda))*2.0/nx
vsp(:) = vp*COS(xv(:) - lambda)
usp(:) = -vp*SIN(xv(:) - lambda)

! North pole
a = SUM(u(:,ny)*SIN(xu(:)))
b = SUM(u(:,ny)*COS(xu(:)))
lambda = ATAN2(-b,a)
vp = (a*COS(lambda) - b*SIN(lambda))*2.0/nx
vnp(:) = vp*COS(xv(:) - lambda)
unp(:) = vp*SIN(xv(:) - lambda)


END SUBROUTINE polar

! ========================================================

SUBROUTINE coriolis(u,v,phi,fu,fv)

! To evaluate the Coriolis terms on the C-grid,
! taking account of energy conservation and improved
! Rossby mode dispersion

USE grid
USE constants

IMPLICIT NONE

! INTEGER, INTENT(IN) :: nx, ny
REAL*8, INTENT(IN) :: u(nx,ny), v(nx,ny), phi(nx,ny)
!REAL*8, INTENT(OUT) :: fu(nx,ny), fv(nx,ny)
REAL*8, INTENT(OUT) :: fu(nx,ny+1), fv(nx,ny)
REAL*8 :: tempv(nx,ny+1), tempu(nx,ny), tempp(nx,ny)

INTEGER :: i, im, ip, j, jm, jp

! ------

! phi v cos(lat) at v points
DO j = 2, ny
  jm = j - 1
  DO i = 1, nx
    tempv(i,j) = 0.5*(phi(i,jm) + phi(i,j))*cosv(j)*v(i,j)
  ENDDO
ENDDO
! zero at polar latitudes
tempv(:,1) = 0.0
tempv(:,ny+1) = 0.0

! Average to phi points and times f / phi
DO j = 1, ny
  jp = j + 1
  DO i = 1, nx
    tempp(i,j) = 0.5*(tempv(i,j) + tempv(i,jp))*twoomega*singeolatp(i,j)/phi(i,j)
  ENDDO
ENDDO

! Average to u points and divide by cos(lat) to get
! fv at u points
DO j = 1, ny
  DO i = 1, nx
    im = i-1
    IF (im == 0) im = nx
    fv(i,j) = 0.5*(tempp(im,j) + tempp(i,j))/cosp(j)
  ENDDO
ENDDO

! ------

! phi u at u points
DO j = 1, ny
  DO i = 1, nx
    im = i - 1
    IF (im == 0) im = nx
    tempu(i,j) = 0.5*(phi(im,j) + phi(i,j))*u(i,j)
  ENDDO
ENDDO

! Average to phi points and times f / phi
DO j = 1, ny
  DO i = 1, nx
    ip = i + 1
    IF (i == nx) ip = 1
    tempp(i,j) = 0.5*(tempu(i,j) + tempu(ip,j))*twoomega*singeolatp(i,j)/phi(i,j)
  ENDDO
ENDDO

! Average to v points to get
! fu at v points
DO j = 2, ny
  jm = j - 1
  DO i = 1, nx
    fu(i,j) = 0.5*(tempp(i,jm) + tempp(i,j))
  ENDDO
ENDDO
! zero at polar latitudes
fu(:,1) = 0.0
fu(:,ny+1) = 0.0

! ------

END SUBROUTINE coriolis

! ========================================================

SUBROUTINE departure

! Calculate departure points for a given estimate
! of the current and new u and v and their spatial
! averages


USE alldata

IMPLICIT NONE

INTEGER :: ndepit = 2, idepit, i, j, k, l, kp, lp, k1, k1p, k2, k2p, &
           hnx

REAL*8 :: a1, a2, b1, b2, ud, vd, sina, cosa, sind, cosd, rk, rl, &
          flip1, flip2, sinad, cosad, sdl, cdl, den, urot, vrot, &
	  m11, m12, m21, m22, x, y, r, dlambda


! ---------------------------------------------------------


! Handy quantity for polar interpolation
hnx = nx/2

! Loop over iterations
DO idepit = 1, ndepit
  
  ! u-point departure points
  DO  j = 1, ny
  
    sina = sinp(j)
    cosa = cosp(j)
    DO i = 1, nx
    
      ! Trig factors at estimated departure point
      sind = SIN(ydu(i,j))
      cosd = COS(ydu(i,j))
    
      ! Determine departure cell index based on latest estimate
      ! xdu, ydu
      rk = (xdu(i,j) - xu(1))/dx
      k = FLOOR(rk)
      a2 = rk - k
      a1 = 1.0 - a2
      k = MODULO(k, nx) + 1
      rl = (ydu(i,j) - yu(1))/dy
      l = FLOOR(rl)
      b2 = rl - l
      b1 = 1.0 - b2
      l = l + 1
      ! l = 0 or l = ny ??
      kp = k+1
      IF (k == nx) kp = 1
      lp = l+1
      
      ! Tricks to handle polar case
      ! Here we interpolate across the pole.
      ! (An alternative would be to use the polar values of u and v
      ! and interpolate between nearest u latitude and the pole.)
      k1 = k
      k1p = kp
      k2 = k
      k2p = kp
      flip1 = 1.0
      flip2 = 1.0
      IF (l == 0) THEN ! South pole
        l = 1
        k1 = MODULO(k1 + hnx - 1, nx) + 1
        k1p = MODULO(k1p + hnx - 1, nx ) + 1
        flip1 = -1.0
      ELSEIF(l == ny) THEN ! North pole
        lp = ny
        k2 = MODULO(k2 + hnx - 1, nx) + 1
        k2p = MODULO(k2p + hnx - 1, nx ) + 1
        flip2 = -1.0
      ENDIF
      
      ! Linearly interpolate  velocity to estimated departure point
      ud = (a1*b1*u0(k1,l) &
          + a2*b1*u0(k1p,l))*flip1 &
	 + (a1*b2*u0(k2,lp) &
	  + a2*b2*u0(k2p,lp))*flip2
      vd = (a1*b1*v0bar(k1,l) &
          + a2*b1*v0bar(k1p,l))*flip1 &
	 + (a1*b2*v0bar(k2,lp) &
	  + a2*b2*v0bar(k2p,lp))*flip2
	 
      ! Rotate to arrival point Cartesian system
      sinad = sina*sind
      cosad = cosa*cosd
      sdl = SIN(xu(i) - xdu(i,j))
      cdl = COS(xu(i) - xdu(i,j))
      den = 1.0 + sinad + cosad*cdl
      m11 = (cosad + (1.0 + sinad)*cdl) / den
      m12 = (sina + sind)*sdl / den
      m21 = -m12
      m22 = m11
      urot = ud*m11 + vd*m12
      vrot = ud*m21 + vd*m22

      ! Hence calculate better estimate of departure point
      ! in arrival point Cartesian system
      x = -hdt*(u(i,j) + urot)
      y = -hdt*(vbar(i,j) + vrot)

      ! Project back to spherical coordinate system
      r = SQRT(x*x + y*y + rearth*rearth)
      sind = (y*cosa + rearth*sina)/r
      ydu(i,j) = ASIN(sind)
      dlambda = ATAN2(x,rearth*cosa - y*sina)
      xdu(i,j) = MODULO(xu(i) + dlambda, twopi)

    ENDDO
  ENDDO


  ! v-point departure points
  DO  j = 2, ny
    sina = sinv(j)
    cosa = cosv(j)
    DO i = 1, nx
    
      ! Trig factors at estimated departure point
      sind = SIN(ydv(i,j))
      cosd = COS(ydv(i,j))
    
      ! Determine departure cell index based on latest estimate
      ! xdv, ydv
      rk = (xdv(i,j) - xv(1))/dx
      k = FLOOR(rk)
      a2 = rk - k
      a1 = 1.0 - a2
      k = MODULO(k, nx) + 1
      rl = (ydv(i,j) - yv(1))/dy
      l = FLOOR(rl)
      b2 = rl - l
      b1 = 1.0 - b2
      l = l + 1
      kp = k+1
      IF (k == nx) kp = 1
      lp = l+1
      
      ! Linearly interpolate  velocity to estimated departure point
      ud = a1*b1*u0bar(k,l) &
         + a2*b1*u0bar(kp,l) &
	 + a1*b2*u0bar(k,lp) &
	 + a2*b2*u0bar(kp,lp)
      vd = a1*b1*v0(k,l) &
         + a2*b1*v0(kp,l) &
	 + a1*b2*v0(k,lp) &
	 + a2*b2*v0(kp,lp)

      ! Rotate to arrival point Cartesian system
      sinad = sina*sind
      cosad = cosa*cosd
      sdl = SIN(xv(i) - xdv(i,j))
      cdl = COS(xv(i) - xdv(i,j))
      den = 1.0 + sinad + cosad*cdl
      m11 = (cosad + (1.0 + sinad)*cdl) / den
      m12 = (sina + sind)*sdl / den
      m21 = -m12
      m22 = m11
      urot = ud*m11 + vd*m12
      vrot = ud*m21 + vd*m22 

      ! Hence calculate better estimate of departure point
      ! in arrival point Cartesian system
      x = -hdt*(ubar(i,j) + urot)
      y = -hdt*(v(i,j) + vrot)

      ! Project back to spherical coordinate system
      r = SQRT(x*x + y*y + rearth*rearth)
      sind = (y*cosa + rearth*sina)/r
      ydv(i,j) = ASIN(sind)
      dlambda = ATAN2(x,rearth*cosa - y*sina)
      xdv(i,j) = MODULO(xv(i) + dlambda, twopi)
      
    ENDDO
  ENDDO

ENDDO


END SUBROUTINE departure

! ========================================================

SUBROUTINE departuremod

! Calculate departure points for a given estimate
! of the current and new u and v and their spatial
! averages


USE alldata

IMPLICIT NONE

INTEGER :: ndepit = 2, idepit, i, j, k, l, kp, lp, k1, k1p, k2, k2p, &
           hnx

REAL*8 :: a1, a2, b1, b2, ud, vd, sina, cosa, sind, cosd, rk, rl, &
          flip1, flip2, sinad, cosad, sdl, cdl, den, urot, vrot, &
	  m11, m12, m21, m22, x, y, r, dlambda


! ---------------------------------------------------------


! Handy quantity for polar interpolation
hnx = nx/2

! Loop over iterations
DO idepit = 1, ndepit
  
  ! u-point departure points
  DO  j = 1, ny
  
    sina = sinp(j)
    cosa = cosp(j)
    DO i = 1, nx
    
      ! Trig factors at estimated departure point
      sind = SIN(ydumod(i,j))
      cosd = COS(ydumod(i,j))
    
      ! Determine departure cell index based on latest estimate
      ! xdumod, ydumod
      rk = (xdumod(i,j) - xu(1))/dx
      k = FLOOR(rk)
      a2 = rk - k
      a1 = 1.0 - a2
      k = MODULO(k, nx) + 1
      rl = (ydumod(i,j) - yu(1))/dy
      l = FLOOR(rl)
      b2 = rl - l
      b1 = 1.0 - b2
      l = l + 1
      ! l = 0 or l = ny ??
      kp = k+1
      IF (k == nx) kp = 1
      lp = l+1
      
      ! Tricks to handle polar case
      ! Here we interpolate across the pole.
      ! (An alternative would be to use the polar values of u and v
      ! and interpolate between nearest u latitude and the pole.)
      k1 = k
      k1p = kp
      k2 = k
      k2p = kp
      flip1 = 1.0
      flip2 = 1.0
      IF (l == 0) THEN ! South pole
        l = 1
        k1 = MODULO(k1 + hnx - 1, nx) + 1
        k1p = MODULO(k1p + hnx - 1, nx ) + 1
        flip1 = -1.0
      ELSEIF(l == ny) THEN ! North pole
        lp = ny
        k2 = MODULO(k2 + hnx - 1, nx) + 1
        k2p = MODULO(k2p + hnx - 1, nx ) + 1
        flip2 = -1.0
      ENDIF
      
      ! Linearly interpolate  velocity to estimated departure point
      ud = (a1*b1*u0mod(k1,l) &
          + a2*b1*u0mod(k1p,l))*flip1 &
	 + (a1*b2*u0mod(k2,lp) &
	  + a2*b2*u0mod(k2p,lp))*flip2
      vd = (a1*b1*v0modbar(k1,l) &
          + a2*b1*v0modbar(k1p,l))*flip1 &
	 + (a1*b2*v0modbar(k2,lp) &
	  + a2*b2*v0modbar(k2p,lp))*flip2
	 
      ! Rotate to arrival point Cartesian system
      sinad = sina*sind
      cosad = cosa*cosd
      sdl = SIN(xu(i) - xdumod(i,j))
      cdl = COS(xu(i) - xdumod(i,j))
      den = 1.0 + sinad + cosad*cdl
      m11 = (cosad + (1.0 + sinad)*cdl) / den
      m12 = (sina + sind)*sdl / den
      m21 = -m12
      m22 = m11
      urot = ud*m11 + vd*m12
      vrot = ud*m21 + vd*m22

      ! Hence calculate better estimate of departure point
      ! in arrival point Cartesian system
      x = -hdt*(umod(i,j) + urot)
      y = -hdt*(vmodbar(i,j) + vrot)

      ! Project back to spherical coordinate system
      r = SQRT(x*x + y*y + rearth*rearth)
      sind = (y*cosa + rearth*sina)/r
      ydumod(i,j) = ASIN(sind)
      dlambda = ATAN2(x,rearth*cosa - y*sina)
      xdumod(i,j) = MODULO(xu(i) + dlambda, twopi)

    ENDDO
  ENDDO


  ! v-point departure points
  DO  j = 2, ny
    sina = sinv(j)
    cosa = cosv(j)
    DO i = 1, nx
    
      ! Trig factors at estimated departure point
      sind = SIN(ydvmod(i,j))
      cosd = COS(ydvmod(i,j))
    
      ! Determine departure cell index based on latest estimate
      ! xdv, ydv
      rk = (xdvmod(i,j) - xv(1))/dx
      k = FLOOR(rk)
      a2 = rk - k
      a1 = 1.0 - a2
      k = MODULO(k, nx) + 1
      rl = (ydvmod(i,j) - yv(1))/dy
      l = FLOOR(rl)
      b2 = rl - l
      b1 = 1.0 - b2
      l = l + 1
      kp = k+1
      IF (k == nx) kp = 1
      lp = l+1
      
      ! Linearly interpolate  velocity to estimated departure point
      ud = a1*b1*u0modbar(k,l) &
         + a2*b1*u0modbar(kp,l) &
	 + a1*b2*u0modbar(k,lp) &
	 + a2*b2*u0modbar(kp,lp)
      vd = a1*b1*v0mod(k,l) &
         + a2*b1*v0mod(kp,l) &
	 + a1*b2*v0mod(k,lp) &
	 + a2*b2*v0mod(kp,lp)

      ! Rotate to arrival point Cartesian system
      sinad = sina*sind
      cosad = cosa*cosd
      sdl = SIN(xv(i) - xdvmod(i,j))
      cdl = COS(xv(i) - xdvmod(i,j))
      den = 1.0 + sinad + cosad*cdl
      m11 = (cosad + (1.0 + sinad)*cdl) / den
      m12 = (sina + sind)*sdl / den
      m21 = -m12
      m22 = m11
      urot = ud*m11 + vd*m12
      vrot = ud*m21 + vd*m22 

      ! Hence calculate better estimate of departure point
      ! in arrival point Cartesian system
      x = -hdt*(umodbar(i,j) + urot)
      y = -hdt*(vmod(i,j) + vrot)

      ! Project back to spherical coordinate system
      r = SQRT(x*x + y*y + rearth*rearth)
      sind = (y*cosa + rearth*sina)/r
      ydvmod(i,j) = ASIN(sind)
      dlambda = ATAN2(x,rearth*cosa - y*sina)
      xdvmod(i,j) = MODULO(xv(i) + dlambda, twopi)
      
    ENDDO
  ENDDO

ENDDO

! Assign polar departure points to equal arrival points
xdvmod(:,1) = xv(:)
ydvmod(:,1) = -piby2
xdvmod(:,ny+1) = xv(:)
ydvmod(:,ny+1) = piby2


END SUBROUTINE departuremod

! ========================================================

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

SUBROUTINE slice1d(xg,xd,domain,q,qnew,qg,n)

! To use the SLICE-1D algorithm to conservatively
! advect a quantity q on a one-dimensional periodic
! domain

IMPLICIT NONE

INTEGER, INTENT(IN) :: n
REAL*8, INTENT(IN) :: xg(n), xd(n), domain, q(n)
REAL*8, INTENT(OUT) :: qnew(n), qg(n)
INTEGER :: i, idx(n), im, ip, k, j, length, jj
REAL*8 :: a(n), b(n), c(n), r(n), dx(n), rdx(n), xi(n), &
        part(n), a0(n), a1(n), a2(n), xx, sum


! Grid intervals and reciprocals
DO i = 1, n
  ip = MODULO(i,n)+1
  dx(i) = MODULO(xg(ip) - xg(i), domain)
  rdx(i) = 1.0/dx(i)
ENDDO

! Find indices to departure cells
! and fractions of cells
DO i = 1, n
  IF (xd(i) .ge. xg(n) .or. xd(i) .le. xg(1)) THEN
    k = n
  ELSE
    k = CEILING((xd(i) - xg(1))/dx(1))
    k = MAX(MIN(k,n),1)
    ! Safety check for irregular grid
    DO WHILE(xd(i) .lt. xg(k))
      k = k - 1
    ENDDO
    DO WHILE(xd(i) .gt. xg(k+1))
      k = k + 1
    ENDDO
  ENDIF
  idx(i) = k
  xi(i) = MODULO(xd(i) - xg(k),domain)*rdx(k)
ENDDO

! Set up coefficients for tridiagonal problem
! to determine parabolic spline fit
DO i = 1, n
  im = i-1
  IF (i == 1) im = n
  a(i) = rdx(im)
  b(i) = 2.0*(rdx(im) + rdx(i))
  c(i) = rdx(i)
  r(i) = 3.0*(q(im)*rdx(im) + q(i)*rdx(i))
ENDDO

! Solve tridiagonal problem
! to obtain cell edge values qg
CALL trisolve(qg,a,b,c,r,n)

! Hence find coefficients of parabolas
DO i = 1, n
  ip = i + 1
  IF (i == n) ip = 1
  a0(i) = qg(i)
  a1(i) = -2*qg(i) - qg(ip) + 3*q(i)    ! ZWS coeff / 2
  a2(i) = qg(i) + qg(ip) - 2*q(i)       ! ZWS coeff / 3
ENDDO

! Compute partial integrals for each departure point
! and grid point value of q at departure point
DO i = 1, n
  k = idx(i)
  xx = xi(i)
  part(i) = ((((a2(k)*xx + a1(k))*xx) + a0(k))*xx)*dx(k)
  qg(i) = (3.0*a2(k)*xx + 2.0*a1(k))*xx + a0(k)
ENDDO

! Finally compute integrals between departure points
! and update values of q
DO i = 1, n
  ip = i + 1
  IF (i == n) ip = 1
  sum = part(ip) - part(i)
  length = MODULO(idx(ip) - idx(i), n)
  DO j = 1,length
    jj = MODULO(idx(i) + j - 2,n) + 1
    sum = sum + q(jj)*dx(jj)
  ENDDO
  qnew(i) = sum*rdx(i)
ENDDO



END SUBROUTINE slice1d

! ==========================================================

SUBROUTINE slice1db(xg,xd,q,qnew,qg,n)

! To use the SLICE-1D algorithm to conservatively
! advect a quantity q on a one-dimensional bounded
! domain

IMPLICIT NONE

INTEGER, INTENT(IN) :: n
REAL*8, INTENT(IN) :: xg(n+1), xd(n+1), q(n)
REAL*8, INTENT(OUT) :: qnew(n), qg(n+1)
INTEGER :: i, idx(n+1), im, ip, k, j, length, jj
REAL*8 :: a(n+1), b(n+1), c(n+1), r(n+1), dx(n), rdx(n), xi(n+1), &
        part(n+1), a0(n), a1(n), a2(n), xx, sum


! Grid intervals and reciprocals
DO i = 1, n
  ip = i + 1
  dx(i) = xg(ip) - xg(i)
  rdx(i) = 1.0/dx(i)
ENDDO

! Find indices to departure cells
! and fractions of cells
! xd(1) should equal xg(1) and
! xd(n+1) should equal xg(n+1)

DO i = 1, n+1
  IF (xd(i) .ge. xg(n)) THEN
    k = n
  ELSEIF (xd(i) .le. xg(1)) THEN
    k = 1
  ELSE
    k = CEILING((xd(i) - xg(1))/dx(1))
    k = MAX(MIN(k,n),1)
    ! Safety check for irregular grid
    ! Ensure xd(i) lies between xg(k) and xg(k+1)
    DO WHILE(xd(i) .lt. xg(k))
      k = k - 1
    ENDDO
    DO WHILE(xd(i) .gt. xg(k+1))
      k = k + 1
    ENDDO
  ENDIF
  idx(i) = k
  xi(i) = (xd(i) - xg(k))*rdx(k)
ENDDO

! Set up coefficients for tridiagonal problem
! to determine parabolic spline fit
!
! For  N-S sweep on sphere bc's should be zero when using distance
! as the coordinate (1/8/08)
! But should remain zero gradient when using area as coordinate
! (3/8/08)
!
a(1) = 0.0
b(1) = 2.0
c(1) = 1.0
r(1) = 3.0*q(1)
! a(1) = 0.0
! b(1) = 1.0
! c(1) = 0.0
! r(1) = 0.0
DO i = 2, n
  im = i-1
  a(i) = rdx(im)
  b(i) = 2.0*(rdx(im) + rdx(i))
  c(i) = rdx(i)
  r(i) = 3.0*(q(im)*rdx(im) + q(i)*rdx(i))
ENDDO
a(n+1) = 1.0
b(n+1) = 2.0
c(n+1) = 0.0
r(n+1) = 3.0*q(n)
! a(n+1) = 0.0
! b(n+1) = 1.0
! c(n+1) = 0.0
! r(n+1) = 0.0
! Solve tridiagonal problem
! to obtain cell edge values qg
CALL trisolveb(qg,a,b,c,r,n+1)

! Hence find coefficients of parabolas
DO i = 1, n
  ip = i + 1
  a0(i) = qg(i)
  a1(i) = -2*qg(i) - qg(ip) + 3*q(i)    ! ZWS coeff / 2
  a2(i) = qg(i) + qg(ip) - 2*q(i)       ! ZWS coeff / 3
ENDDO

! Compute partial integrals for each departure point
! and grid point value of q at departure point
DO i = 1, n+1
  k = idx(i)
  xx = xi(i)
  part(i) = ((((a2(k)*xx + a1(k))*xx) + a0(k))*xx)*dx(k)
  qg(i) = (3.0*a2(k)*xx + 2.0*a1(k))*xx + a0(k)
ENDDO

! Finally compute integrals between departure points
! and update values of q
DO i = 1, n
  ip = i + 1
  sum = part(ip) - part(i)
  length = idx(ip) - idx(i)
  DO j = 0,length - 1
    jj = idx(i) + j
    sum = sum + q(jj)*dx(jj)
  ENDDO
  qnew(i) = sum*rdx(i)
ENDDO


END SUBROUTINE slice1db

! ==========================================================

SUBROUTINE slice2d(xdu,ydu,xdv,ydv,q,qnew)

! To advect a quantity q using the 2D SLICE algorithm

USE version
USE grid
USE util

IMPLICIT NONE

REAL*8, INTENT(IN) :: xdu(nx,ny), ydu(nx,ny), xdv(nx,ny+1), ydv(nx,ny+1), &
		      q(nx,ny)
REAL*8, INTENT(OUT) :: qnew(nx,ny)
INTEGER :: i, im, j, jm, jp
REAL*8 :: xiecv(nx,ny), a1, a2, yy, domain, hdomain, &
        qx(nx), qxnew(nx), qy(ny), qynew(ny), qgx(nx), qgy(ny+1), &
        s(ny+1), ds(ny), silcv(ny+1), x1, x2, s1, s2, deltax, deltay, &
	slen, hslen, dsi(ny), den, cl, &
	xdc(nx,ny+1), ydc(nx,ny+1), y1, y2, &
	xtemp(ny), xdu2(ny), dxdu(ny), dq, qcorr(nx,ny), &
	dmdxu(nx,ny)


! Useful constants
domain = twopi
hdomain = 0.5*domain


! Estimate departure points for cell corners
! Used for C-grid correction to SLICE
DO i = 1, nx
  im = i - 1
  IF (i == 1) im = nx
  DO j = 2, ny
    jm = j - 1
    x1 = xdv(im,j)
    x2 = xdv(i,j)
    y1 = ydu(i,jm)
    y2 = ydu(i,j)
    IF (x1 .gt. x2) x2 = x2 + domain
    xdc(i,j) = MODULO(0.5*(x1 + x2),domain)
    ydc(i,j) = MODULO(0.5*(y1 + y2),domain)
  ENDDO
ENDDO
! Polar corners
xdc(:,1) = xdu(:,1)
ydc(:,1) = -piby2
xdc(:,ny+1) = xdu(:,ny)
ydc(:,ny+1) = piby2


! Loop over rows
DO j = 1, ny

  ! Find x values of ends of Intermediate Eulerian Control Volumes
  yy = yu(j)
  DO i = 1, nx
    ! Find indices to u departure points to
    ! North and South of yu(j)
    jm = j
    jp = j
    deltay = ydu(i,j) - yy
    IF (deltay .gt. 0.0) THEN
      DO WHILE (deltay .gt. 0.0)
        jm = jm - 1
	IF (jm > 0 ) THEN
          deltay = ydu(i,jm) - yy
	ELSE
	  deltay = -1.0
	ENDIF
      ENDDO
      jp = jm + 1
    ELSE
      DO WHILE (deltay .le. 0.0)
        jp = jp + 1
	IF (jp < ny + 1) THEN
          deltay = ydu(i,jp) - yy
	ELSE
	  deltay = 1.0
	ENDIF
      ENDDO
      jm = jp - 1
    ENDIF
    ! Assume polar points fixed
    IF (jm == 0) THEN
      y1 = -piby2
      x1 = xu(i)
      y2 = ydu(i,jp)
      x2 = xdu(i,jp)
    ELSEIF (jp == ny+1) THEN
      y1 = ydu(i,jm)
      x1 = xdu(i,jm)
      y2 = piby2
      x2 = xu(i)
    ELSE
      y1 = ydu(i,jm)
      x1 = xdu(i,jm)
      y2 = ydu(i,jp)
      x2 = xdu(i,jp)
    ENDIF
    den = y2 - y1
    a1 = (y2 - yy)/den
    a2 = 1.0 - a1
    x1 = xu(i) + near(x1 - xu(i),domain)
    x2 = xu(i) + near(x2 - xu(i),domain)
    xiecv(i,j) = MODULO(a1*x1 + a2*x2,domain)
  ENDDO
  
  ! Remap to Intermediate Control Volumes
  qx = q(:,j)
  CALL slice1d(xu,xiecv(:,j),domain,qx,qxnew,qgx,nx)


  ! Mass in intermediate control volumes
  qnew(:,j) = qxnew*area(j)
  
  ! q at u-edges of intermediate control volumes
  dmdxu(:,j) = qgx

ENDDO


! Initialize correction to q
qcorr = 0.0


! Loop over columns
DO i = 1, nx
  
  ! Define coordinate s at top and bottom of Lagrangian
  ! Control Volumes (assume polar departure points are
  ! defined to be at the poles)
  s(1) = 0.0
  DO j = 2, ny+1
    jm = j - 1
    deltax = xdv(i,j) - xdv(i,jm)
    IF (deltax .gt. hdomain) THEN
      deltax = deltax - domain
    ELSEIF (deltax .lt. -hdomain) THEN
      deltax = deltax + domain
    ENDIF
    deltay = ydv(i,j) - ydv(i,jm)
    ! Approximation to spherical distances
    cl = COS(0.5*(ydv(i,j) + ydv(i,jm)))
    deltax=cl*deltax
    ds(jm) = SQRT(deltax*deltax + deltay*deltay)
    s(j) = s(jm) + ds(jm)
  ENDDO
  slen = s(ny+1)
  hslen = 0.5*slen
    
  ! Find s values of ends of
  ! Intermediate Lagrangian Control Volumes
  DO j = 2, ny
    yy = yv(j)
    ! Find indices to v departure points to
    ! North and South of yv(j)
    jm = j
    jp = j
    deltay = ydv(i,j) - yy
    IF (deltay .gt. 0.0) THEN
      DO WHILE (deltay .gt. 0.0)
        jm = jm - 1
	IF (jm > 0 ) THEN
          deltay = ydv(i,jm) - yy
	ELSE
	  deltay = -1.0
	ENDIF
      ENDDO
      jp = jm + 1
    ELSE
      DO WHILE (deltay .le. 0.0)
        jp = jp + 1
	IF (jp < ny + 1) THEN
          deltay = ydu(i,jp) - yy
	ELSE
	  deltay = 1.0
	ENDIF
      ENDDO
      jm = jp - 1
    ENDIF
    den = ydv(i,jp) - ydv(i,jm)
    a1 = (ydv(i,jp) - yy)/den
    a2 = 1.0 - a1
    s1 = s(jm)
    s2 = s(jp)
    silcv(j) = a1*s1 + a2*s2
  ENDDO
  silcv(1) = s(1)
  silcv(ny+1) = s(ny+1)
  DO j = 1, ny
    jp = j + 1
    dsi(j) = silcv(jp) - silcv(j)
  ENDDO
  
  ! Remap to Lagrangian Control Volumes
  qy(:) = qnew(i,:)/dsi(:)
  CALL slice1db(silcv,s,qy,qynew,qgy,ny)

  ! Mass in Lagrangian control volumes
  ! equals mass in Eulerian arrival control volumes.
  ! Divide by area to get density
  qnew(i,:) = qynew(:)*dsi(:)/area(:)

ENDDO


! C-grid correction to SLICE
! Effect should be small so use linear interpolation,
! and interpolate in y rather than s

DO i = 1, nx

  im = i - 1
  IF (i == 1) im = nx  

  ! Estimate q at cell edge as average of cell values
  ! either side
  DO j = 1, ny
    ! jp = j + 1
    ! deltay = MODULO(ydc(i,jp) - ydc(i,j),domain)
    ! dmdxu(i,j) = 0.5*(qnew(im,j) + qnew(i,j))*deltay/area(j)
    ! print *,'Check C-grid correction '
    dmdxu(i,j) = 0.5*(qnew(im,j) + qnew(i,j))/dx
  ENDDO
  
  ! Interpolate xd from edges of intermediate Eulerian CVs
  ! to edges of Lagrangian CVs
  ! being careful about wrap-around
  xtemp = xiecv(i,1:ny)
  DO j = 1, ny
    xtemp(j) = xu(i) + near(xtemp(j) - xu(i), domain)
  ENDDO
  CALL lin1db(yu,ydu(i,1:ny),xtemp,xdu2,ny)

  ! Compare estimate with the true xdu
  DO j = 1, ny
    dxdu(j) = near(xdu(i,j) - xdu2(j), domain)
  ENDDO

  ! Correction to masses in cells either side
  DO j = 1, ny
    dq = dxdu(j)*dmdxu(i,j)
    qcorr(i,j) = qcorr(i,j) - dq
    qcorr(im,j) = qcorr(im,j) + dq 
  ENDDO
  
ENDDO


IF (cgridcorr == 1) THEN
  ! Correct q
  qnew = qnew + qcorr
ENDIF


END SUBROUTINE slice2d

! ==========================================================

SUBROUTINE slice2da(xdu,ydu,xdv,ydv,aread,q,qnew)

! To advect a quantity q using the 2D SLICE algorithm
! In this version the departure cell areas are specified and
! integrated cell area is used as the coordinate for remapping
! in the N-S sweeps.

! No need for C-grid correction in this version (departure areas
! already feel the C-grid divergence, and any correction would
! compromise the departure cell areas).

USE version
USE grid
USE util

IMPLICIT NONE

! Note xdv and ydv are not actually used 
REAL*8, INTENT(IN) :: xdu(nx,ny), ydu(nx,ny), xdv(nx,ny+1), ydv(nx,ny+1), &
		      q(nx,ny), aread(nx,ny)
REAL*8, INTENT(OUT) :: qnew(nx,ny)
INTEGER :: i, ip, j, jm, jp
REAL*8 :: xiecv(nx,ny), a1, a2, yy, domain, hdomain, &
        qx(nx), qxnew(nx), qy(ny), qynew(ny), qgx(nx), qgy(ny+1), &
        s(ny+1), silcv(ny+1), x1, x2, deltay, &
	slen, den, &
	y1, y2, aicv(nx,ny), fac, astripd(nx), astripi(nx), &
        acorr(nx), xcorr(nx), dabar, dxbar


! Useful constants
domain = twopi
hdomain = 0.5*domain

! First find ends of Intermediate Control Volumes and corresponding
! Intermediate areas; hence find corresponding total area
! in each N-S strip for comparison with total arrival cell area
! in each N-S strip.
astripi = 0.0
astripd = 0.0
DO j = 1, ny

  yy = yu(j)
  DO i = 1, nx
    ! Find indices to u departure points to
    ! North and South of yu(j)
    jm = j
    jp = j
    deltay = ydu(i,j) - yy
    IF (deltay .gt. 0.0) THEN
      DO WHILE (deltay .gt. 0.0)
        jm = jm - 1
	IF (jm > 0 ) THEN
          deltay = ydu(i,jm) - yy
	ELSE
	  deltay = -1.0
	ENDIF
      ENDDO
      jp = jm + 1
    ELSE
      DO WHILE (deltay .le. 0.0)
        jp = jp + 1
	IF (jp < ny + 1) THEN
          deltay = ydu(i,jp) - yy
	ELSE
	  deltay = 1.0
	ENDIF
      ENDDO
      jm = jp - 1
    ENDIF
    ! Assume polar points fixed
    IF (jm == 0) THEN
      y1 = -piby2
      x1 = xu(i)
      y2 = ydu(i,jp)
      x2 = xdu(i,jp)
    ELSEIF (jp == ny+1) THEN
      y1 = ydu(i,jm)
      x1 = xdu(i,jm)
      y2 = piby2
      x2 = xu(i)
    ELSE
      y1 = ydu(i,jm)
      x1 = xdu(i,jm)
      y2 = ydu(i,jp)
      x2 = xdu(i,jp)
    ENDIF
    den = y2 - y1
    a1 = (y2 - yy)/den
    a2 = 1.0 - a1
    x1 = xu(i) + near(x1 - xu(i),domain)
    x2 = xu(i) + near(x2 - xu(i),domain)
    xiecv(i,j) = MODULO(a1*x1 + a2*x2,domain)
  ENDDO
  
  ! Note areas of intermediate control volumes and accumulate
  ! areas in strips
  DO i = 1, nx
    ip = i + 1
    IF (i == nx) ip = 1
    aicv(i,j) = near((xiecv(ip,j) - xiecv(i,j)),domain)*dy*cosp(j)
    astripi(i) = astripi(i) + aicv(i,j)
    astripd(i) = astripd(i) + aread(i,j)
  ENDDO

ENDDO


! Calculate correction to xiecv to make astripi coinicide
! as closely as possible with astripd. (They can only coincide
! exactly if the global sum of departure areas exactly agrees with
! area of sphere.)

! Required fractional in area
acorr = astripd - astripi

! Tweak xiecv
dabar = SUM(acorr)/nx
acorr = acorr - dabar
xcorr(1) = 0.0
DO i = 2, nx
  xcorr(i) = xcorr(i-1) + acorr(i-1)/dareadx
ENDDO
dxbar = SUM(xcorr)/nx
xcorr = xcorr - dxbar
DO j = 1, ny
  xiecv(:,j) = MODULO(xiecv(:,j) + xcorr,domain)
  aicv(:,j) = aicv(:,j) + acorr*cosp(j)*dy/dareadx
ENDDO

! Main loop over rows
DO j = 1, ny

  ! Remap to Intermediate Control Volumes
  qx = q(:,j)
  CALL slice1d(xu,xiecv(:,j),domain,qx,qxnew,qgx,nx)

  ! Mass in intermediate control volumes
  qnew(:,j) = qxnew*area(j)
  
ENDDO

! Convert mass to density in intermediate control volumes
qnew = qnew/aicv


! Loop over columns
DO i = 1, nx
  
  ! Define coordinate at top and bottom of Lagrangian control volumes
  ! and top and bottom of intermediate Lagrangian control volumes
  ! Assume polar departure points are at the poles
  s(1) = 0.0
  silcv(1) = 0.0
  DO j = 2, ny + 1
    jm = j - 1
    s(j) = s(jm) + aread(i,jm)
    silcv(j) = silcv(jm) + aicv(i,jm)
  ENDDO

  ! Check that total area in the strip is correct
  fac = silcv(ny+1)/s(ny+1)
  s = s*fac
  
  ! Remap to Lagrangian Control Volumes
  qy(:) = qnew(i,:)
  CALL slice1db(silcv,s,qy,qynew,qgy,ny)

  ! Mass in Lagrangian control volumes
  ! equals mass in Eulerian arrival control volumes.
  ! Divide by area to get density
  qnew(i,:) = qynew(:)*aicv(i,:)/area(:)

ENDDO




END SUBROUTINE slice2da

! ==========================================================


SUBROUTINE lin1db(xg,xd,q,qnew,n)

! To do 1d linear interpolation to
! advect a quantity q on a one-dimensional bounded
! domain

IMPLICIT NONE

INTEGER, INTENT(IN) :: n
REAL*8, INTENT(IN) :: xg(n), xd(n), q(n)
REAL*8, INTENT(OUT) :: qnew(n)
INTEGER :: i, idx(n), ip, k
REAL*8 :: dx(n - 1), rdx(n - 1), xx(n)



! Grid intervals and reciprocals
DO i = 1, n - 1
  ip = i + 1
  dx(i) = xg(ip) - xg(i)
  rdx(i) = 1.0/dx(i)
ENDDO

! Find indices to departure points
! and fractions of cells
DO i = 1, n
  IF (xd(i) .ge. xg(n)) THEN
    idx(i) = n - 1
    xx(i) = 1.0
  ELSEIF (xd(i) .le. xg(1)) THEN
    idx(i) = 1
    xx(i) = 0.0
  ELSE
    k = CEILING((xd(i) - xg(1))/dx(1))
    k = MAX(MIN(k,n-1),1)
    ! Safety check for irregular grid
    ! Ensure xd(i) lies between xg(k) and xg(k+1)
    DO WHILE(xd(i) .lt. xg(k))
      k = k - 1
    ENDDO
    DO WHILE(xd(i) .gt. xg(k+1))
      k = k + 1
    ENDDO
    idx(i) = k
    xx(i) = (xd(i) - xg(k))*rdx(k)
  ENDIF
ENDDO

! Interpolate
DO i = 1, n
  qnew(i) = (1.0 - xx(i))*q(idx(i)) + xx(i)*q(idx(i)+1)
ENDDO


END SUBROUTINE lin1db

! ==========================================================

SUBROUTINE momentum

! Calulate terms needed in the momentum equations
! at departure points

USE alldata

IMPLICIT NONE

INTEGER :: i, im, j, jm
REAL*8 :: phitot(nx,ny), &
          rdx, rdy, cdx

rdx = rearth*dx
rdy = rearth*dy

! Total geopotential
phitot = phis + phi0

! Calculate Coriolis terms
!CALL coriolis(nx,ny,u0,v0,phi0,fu,fv)
CALL coriolis(u0,v0,phi0,fu,fv)

! At u-points
DO j = 1, ny
  cdx = cosp(j)*rdx
  DO i = 1, nx
    im = i - 1
    IF (i == 1) im = nx
    ru0(i,j) = u0(i,j) &
             - bdt*( (phitot(i,j) - phitot(im,j))/cdx - fv(i,j) )
  ENDDO
ENDDO

! At v-points
DO j = 2, ny
  jm = j - 1
  DO i = 1, nx 
    rv0(i,j) = v0(i,j) &
             - bdt*( (phitot(i,j) - phitot(i,jm))/rdy + fu(i,j) )
  ENDDO
ENDDO

! Perform C-grid average and find polar values of rv0
CALL cgridave(nx,ny,ru0,rv0,ru0bar,rv0bar)



END SUBROUTINE momentum

! ==========================================================

SUBROUTINE diagnostics(istep)

! Output some basic diagnostics and error measures

USE state
USE work
USE errdiag
USE constants

IMPLICIT NONE
INTEGER, INTENT(IN) :: istep

INTEGER :: i, j, ip, jp
REAL*8 :: phierr(nx,ny), uerr(nx,ny), verr(nx,ny+1), l1phi, l2phi, linfphi, &
          l1u, l2u, linfu, l1v, l2v, linfv, aphi, au, av, aerr, da, mass, &
          div(nx,ny), rmsdiv, tota, uatp, vatp, ke, ape, angmom, totphitot, &
          phitot(nx,ny), wc2, dm, cosr, sinr, ucoslat, meanphitot




! Find polar values of v
CALL polar(u,ubar(:,1),v(:,1),ubar(:,ny+1),v(:,ny+1))

! Compute divergence
CALL divergence(u,v,div)

! Compute L1, L2 and Linf error norms for phi and velocity,
! Assuming steady state is the truth
phierr = phi - phi_init
uerr = u - u_init
verr = v - v_init
l1phi = 0.0d0
l2phi = 0.0d0
linfphi = 0.0d0
l1u = 0.0d0
l2u = 0.0d0
linfu = 0.0d0
l1v = 0.0d0
l2v = 0.0d0
linfv = 0.0d0
aphi = 0.0d0
au = 0.0d0
av = 0.0d0
DO j = 1, ny
  DO i = 1, nx
    da = area(j)
    aphi = aphi + da
    aerr = ABS(phierr(i,j))
    l1phi = l1phi + da*aerr
    l2phi = l2phi + da*aerr*aerr
    linfphi = MAX(linfphi,aerr)
    au = au + da
    aerr = ABS(uerr(i,j))
    l1u = l1u + da*aerr
    l2u = l2u + da*aerr*aerr
    linfu = MAX(linfu,aerr)
  ENDDO
ENDDO
l1phi = l1phi/aphi
l2phi = SQRT(l2phi/aphi)
DO j = 2, ny
  DO i = 1, nx
    da = dx*dy*cosv(j)
    av = av + da
    aerr = ABS(verr(i,j))
    l1v = l1v + da*aerr
    l2v = l2v + da*aerr*aerr
    linfv = MAX(linfv,aerr)
  ENDDO
ENDDO
l1u = (l1u + l1v)/(au + av)
l2u = SQRT((l2u + l2v)/(au + av))
!PRINT *,'Step ',istep
!PRINT *,'l1, l2, linf (phi) = ',l1phi, l2phi, linfphi
!PRINT *,'l1, l2, linf (u  ) = ',l1u  , l2u  , linfu
!PRINT *,' '


phitot = phis + phi0
cosr = COS(rotgrid)
sinr = SIN(rotgrid)
tota = 0.0d0
mass = 0.0d0
totphitot = 0.0d0
rmsdiv = 0.0d0
DO j = 1, ny
  da = area(j)*rearth*rearth
  DO i = 1, nx
    tota = tota + da
    mass = mass + phi(i,j)*da
    totphitot = totphitot + phitot(i,j)*da
    rmsdiv = rmsdiv + div(i,j)*div(i,j)*da
  ENDDO
ENDDO
rmsdiv = SQRT(rmsdiv/tota)
meanphitot = totphitot/tota
ke = 0.0d0
ape = 0.0d0
angmom = 0.0d0
DO j = 1, ny
  jp = j + 1
  da = area(j)*rearth*rearth
  DO i = 1, nx
    ip = i + 1
    IF (i == nx) ip = 1
    dm = phi(i,j)*da
    uatp = 0.5d0*(u(i,j) + u(ip,j))
    vatp = 0.5d0*(v(i,j) + v(i,jp))
    ke = ke + 0.5d0*(uatp**2 + vatp**2)*dm
    ape = ape + 0.5d0*da*(phitot(i,j) - meanphitot)**2
    ucoslat = uatp*(cosr*cosp(j) + sinr*SIN(xp(i)*sinp(j))) + vatp*sinr*COS(xp(i))
    wc2 = 0.5d0*twoomega*COS(geolatp(i,j))**2
    angmom = angmom + (wc2*rearth*rearth + ucoslat*rearth)*dm
  ENDDO
ENDDO


PRINT *,'Step ',istep
PRINT *,'Mass = ',mass,'  ke = ',ke,'  APE = ',ape,'  Angmom = ',angmom
WRITE(61,555) istep, mass, rmsdiv, ccount, ke, ape, angmom
555 FORMAT(I8.8,E20.8,E24.15,I12.8,3E24.15)


END SUBROUTINE diagnostics

! ==========================================================

SUBROUTINE accumulate

USE state
USE work
USE constants
USE force
USE accumdiag

IMPLICIT NONE
INTEGER :: i, im, j, jm
REAL*8 :: temp, xi(nx,ny+1), q(nx,ny+1), div(nx,ny), phiprime(nx,ny)

! Compute divergence
CALL divergence(u,v,div)

! Compute vorticity
! Find polar values of v
CALL polar(u,ubar(:,1),v(:,1),ubar(:,ny+1),v(:,ny+1))
DO j = 2, ny
  jm = j - 1
  DO i = 1, nx
    im = MODULO(i - 2,nx) + 1
    xi(i,j) = ((v(i,j) - v(im,j))/dx - (u(i,j)*cosp(j) - u(i,jm)*cosp(jm))/dy)/(cosv(j)*rearth)
  ENDDO
ENDDO
! South pole
temp = SUM(u(:,1))*cosp(1)*8/(dy*dy*nx*rearth)
xi(:,1) = -temp
! North pole
temp = SUM(u(:,ny))*cosp(ny)*8/(dy*dy*nx*rearth)
xi(:,ny+1) = temp

! Accumulate
naccum = naccum + 1
phiprime = phi - 1.0d3*gravity
accumphi = accumphi + phiprime
accumdiv = accumdiv + div
accumxi = accumxi + xi
accumphi2 = accumphi2 + phiprime*phiprime
accumdiv2 = accumdiv2 + div*div
accumxi2 = accumxi2 + xi*xi

IF (lpdc) THEN
  DO j = 1, ny
    DO i = 1, nx
      IF ((phi(i,j) - phimem(i,j) > thrconv).and.(philat(i,j) > 0.0D0)) &
          accumcon(i,j) = accumcon(i,j) + 1.0d0
    ENDDO
  ENDDO
ENDIF


END SUBROUTINE accumulate

! ==========================================================

SUBROUTINE outstate(istep)

! Output the model state

USE state
USE work
USE errdiag
USE constants

IMPLICIT NONE
INTEGER, INTENT(IN) :: istep
CHARACTER*6 :: ystep
CHARACTER*10 :: yname
CHARACTER*24 :: ytitle
INTEGER :: i, im, j, jm, ip, jp
REAL*8 :: xi(nx,ny+1), q(nx,ny+1), &
          mbar, dm, tote, da, source, div(nx,ny), ros(nx,ny), s1, s2, temp


! Find polar values of v
CALL polar(u,ubar(:,1),v(:,1),ubar(:,ny+1),v(:,ny+1))

DO j = 2, ny
  jm = j - 1
  DO i = 1, nx
   im = MODULO(i - 2,nx) + 1
   xi(i,j) = ((v(i,j) - v(im,j))/dx - (u(i,j)*cosp(j) - u(i,jm)*cosp(jm))/dy)/(cosv(j)*rearth)
   mbar = 0.25*(phi(i,j) + phi(im,j) + phi(i,jm) + phi(im,jm))
   q(i,j) = (twoomega*singeolatz(i,j) + xi(i,j)) / mbar
  ENDDO
ENDDO
! South pole
temp = SUM(u(:,1))*cosp(1)*8/(dy*dy*nx*rearth)
xi(:,1) = -temp
temp = SUM(phi(:,1))/nx
q(:,1) = (twoomega*singeolatz(1,1)+xi(:,1))/temp
! North pole
temp = SUM(u(:,ny))*cosp(ny)*8/(dy*dy*nx*rearth)
xi(:,ny+1) = temp
temp = SUM(phi(:,ny))/nx
q(:,ny+1) = (twoomega*singeolatz(1,ny+1) + xi(:,ny+1))/temp



WRITE(ystep,'(I6.6)') istep
yname = 'dump'//ystep
OPEN(23,FILE=yname)
WRITE(23,*) nx, ny
WRITE(23,*) xu, yu, xv, yv
ytitle = 'u    step '//ystep
WRITE(23,*) ytitle
WRITE(23,*) u
ytitle = 'v    step '//ystep
WRITE(23,*) ytitle
WRITE(23,*) v
ytitle = 'phi + phis step '//ystep
WRITE(23,*) ytitle
WRITE(23,*) phi + phis
ytitle = 'xi  step '//ystep
WRITE(23,*) ytitle
WRITE(23,*) xi
!ytitle = 'PV  step '//ystep
!WRITE(23,*) ytitle
!WRITE(23,*) q
!ytitle = 'Div step '//ystep
!WRITE(23,*) ytitle
!WRITE(23,*) div


CLOSE(23)


END SUBROUTINE outstate

! ==========================================================

SUBROUTINE outstatem(istep)

! Output the model state for matlab plotting

USE state
USE work
USE errdiag
USE constants
USE force

IMPLICIT NONE
INTEGER, INTENT(IN) :: istep
CHARACTER*8 :: ystep
CHARACTER*14 :: yname
CHARACTER*24 :: ytitle
INTEGER :: i, im, j, jm, ip, jp
REAL*8 :: xi(nx,ny+1), q(nx,ny+1), &
          mbar, dm, tote, da, source, div(nx,ny), ros(nx,ny), s1, s2, temp


! Write header information
WRITE(ystep,'(I8.8)') istep
yname = 'dump'//ystep//'.m'
OPEN(23,FILE=yname)
WRITE(23,*) 'nx = ',nx,';'
WRITE(23,*) 'ny = ',ny,';'
WRITE(23,*) 'xu = [...'
WRITE(23,888) xu
WRITE(23,*) '];'
WRITE(23,*) 'yu = [...'
WRITE(23,888) yu
WRITE(23,*) '];'
WRITE(23,*) 'xv = [...'
WRITE(23,888) xv
WRITE(23,*) '];'
WRITE(23,*) 'yv = [...'
WRITE(23,888) yv
WRITE(23,*) '];'



! Find polar values of v
CALL polar(u,ubar(:,1),v(:,1),ubar(:,ny+1),v(:,ny+1))

DO j = 2, ny
  jm = j - 1
  DO i = 1, nx
   im = MODULO(i - 2,nx) + 1
   xi(i,j) = ((v(i,j) - v(im,j))/dx - (u(i,j)*cosp(j) - u(i,jm)*cosp(jm))/dy)/(cosv(j)*rearth)
   mbar = 0.25*(phi(i,j) + phi(im,j) + phi(i,jm) + phi(im,jm))
   q(i,j) = (twoomega*singeolatz(i,j) + xi(i,j)) / mbar
  ENDDO
ENDDO
! South pole
temp = SUM(u(:,1))*cosp(1)*8/(dy*dy*nx*rearth)
xi(:,1) = -temp
temp = SUM(phi(:,1))/nx
q(:,1) = (twoomega*singeolatz(1,1)+xi(:,1))/temp
! North pole
temp = SUM(u(:,ny))*cosp(ny)*8/(dy*dy*nx*rearth)
xi(:,ny+1) = temp
temp = SUM(phi(:,ny))/nx
q(:,ny+1) = (twoomega*singeolatz(1,ny+1) + xi(:,ny+1))/temp

! Compute divergence
CALL divergence(u,v,div)


! WRITE(23,*)   'ytitle = ''ueqm  step ',ystep,''';'
! WRITE(23,*)   'q = [ ...'
! WRITE(23,888) ueqm
! WRITE(23,*)   ' ];'
! WRITE(23,*)   'q = reshape(q,nx,ny);'
! WRITE(23,*)   'figure(1)'
! WRITE(23,*)   'contour(xu,yu,q'',11,''k'')'
! WRITE(23,*)   'qx = max(max(q));'
! WRITE(23,*)   'qn = min(min(q));'
! WRITE(23,*)   'title ([ytitle,''  Min = '',num2str(qn),''  Max = '',num2str(qx) ])'
! WRITE(23,*)   'pause'

WRITE(23,*)   'ytitle = ''u  step ',ystep,''';'
WRITE(23,*)   'q = [ ...'
WRITE(23,888) u
WRITE(23,*)   ' ];'
WRITE(23,*)   'q = reshape(q,nx,ny);'
WRITE(23,*)   'figure(1)'
WRITE(23,*)   'contour(xu,yu,q'',11,''k'')'
WRITE(23,*)   'qx = max(max(q));'
WRITE(23,*)   'qn = min(min(q));'
WRITE(23,*)   'title ([ytitle,''  Min = '',num2str(qn),''  Max = '',num2str(qx) ])'
WRITE(23,*)   'pause'

! WRITE(23,*)   'ytitle = ''v  step ',ystep,''';'
! WRITE(23,*)   'q = [ ...'
! WRITE(23,888) v
! WRITE(23,*)   ' ];'
! WRITE(23,*)   'q = reshape(q,nx,ny+1);'
! WRITE(23,*)   'figure(1)'
! WRITE(23,*)   'contour(xv,yv,q'',11,''k'')'
! WRITE(23,*)   'qx = max(max(q));'
! WRITE(23,*)   'qn = min(min(q));'
! WRITE(23,*)   'title ([ytitle,''  Min = '',num2str(qn),''  Max = '',num2str(qx) ])'
! WRITE(23,*)   'pause'

! WRITE(23,*)   'ytitle = ''phieqm  step ',ystep,''';'
! WRITE(23,*)   'q = [ ...'
! WRITE(23,888) phieqm
! WRITE(23,*)   ' ];'
! WRITE(23,*)   'q = reshape(q,nx,ny);'
! WRITE(23,*)   'figure(1)'
! WRITE(23,*)   'contour(xv,yu,q'',11,''k'')'
! WRITE(23,*)   'qx = max(max(q));'
! WRITE(23,*)   'qn = min(min(q));'
! WRITE(23,*)   'title ([ytitle,''  Min = '',num2str(qn),''  Max = '',num2str(qx) ])'
! WRITE(23,*)   'pause'

WRITE(23,*)   'ytitle = ''phi  step ',ystep,''';'
WRITE(23,*)   'q = [ ...'
WRITE(23,888) phi
WRITE(23,*)   ' ];'
WRITE(23,*)   'q = reshape(q,nx,ny);'
WRITE(23,*)   'figure(1)'
WRITE(23,*)   'contour(xv,yu,q'',11,''k'')'
WRITE(23,*)   'qx = max(max(q));'
WRITE(23,*)   'qn = min(min(q));'
WRITE(23,*)   'title ([ytitle,''  Min = '',num2str(qn),''  Max = '',num2str(qx) ])'
WRITE(23,*)   'pause'

WRITE(23,*)   'ytitle = ''div  step ',ystep,''';'
WRITE(23,*)   'q = [ ...'
WRITE(23,888) div
WRITE(23,*)   ' ];'
WRITE(23,*)   'q = reshape(q,nx,ny);'
WRITE(23,*)   'figure(1)'
WRITE(23,*)   'subplot(2,1,1)'
WRITE(23,*)   'contour(xv,yu,q'',[-0.3e-4,-0.015e-4],''k'')'
WRITE(23,*)   'axis equal'
WRITE(23,*)   'axis ([0,2*pi,-0.5*pi,0.5*pi])'
WRITE(23,*)   'qx = max(max(q));'
WRITE(23,*)   'qn = min(min(q));'
WRITE(23,*)   'title ([ytitle,''  Min = '',num2str(qn,''%8.2e''),''  Max = '',num2str(qx,''%8.2e'') ])'
!!WRITE(23,*)   'title ([ytitle])'
WRITE(23,*)   'pause'

!WRITE(23,*)   'ytitle = ''phi + phis  step ',ystep,''';'
!WRITE(23,*)   'q = [ ...'
!WRITE(23,888) phi + phis
!WRITE(23,*)   ' ];'
!WRITE(23,*)   'q = reshape(q,nx,ny);'
!WRITE(23,*)   'figure(1)'
!WRITE(23,*)   'v = [5:0.5:12]*1E3;'   ! For Held-Suarez
!WRITE(23,*)   'contour(xv,yu,q'',v,''k'')'
!! WRITE(23,*)   'contour(xv,yu,q'',11,''k'')'
!WRITE(23,*)   'qx = max(max(q));'
!WRITE(23,*)   'qn = min(min(q));'
!WRITE(23,*)   'title ([ytitle,''  Min = '',num2str(qn),''  Max = '',num2str(qx) ])'
!WRITE(23,*)   'pause'

!WRITE(23,*)   'ytitle = ''xi  step ',ystep,''';'
!WRITE(23,*)   'q = [ ...'
!WRITE(23,888) xi
!WRITE(23,*)   ' ];'
!WRITE(23,*)   'q = reshape(q,nx,ny+1);'
!WRITE(23,*)   'figure(1)'
!! WRITE(23,*)   'v = [-8:1:8]*2E-5;'   ! For Galewsky
!WRITE(23,*)   'v = [-8:1:8]*5E-6;'   ! For Held-Suarez
!WRITE(23,*)   'contour(xu,yv,q'',v,''k'',''linewidth'',1)'
!WRITE(23,*)   'hold on'
!WRITE(23,*)   'v = [0,0];'
!WRITE(23,*)   'contour(xu,yv,q'',v,''k'',''linewidth'',2)'
!! WRITE(23,*)   'axis([0,2*pi,pi/18,4*pi/9])'
!WRITE(23,*)   'qx = max(max(q));'
!WRITE(23,*)   'qn = min(min(q));'
!WRITE(23,*)   'title ([ytitle,''  Min = '',num2str(qn),''  Max = '',num2str(qx) ])'
!WRITE(23,*)   'pause'
!WRITE(23,*)   'hold off'

WRITE(23,*)   'ytitle = ''PV  step ',ystep,''';'
WRITE(23,*)   'q = [ ...'
WRITE(23,888) q
WRITE(23,*)   ' ];'
WRITE(23,*)   'q = reshape(q,nx,ny+1);'
WRITE(23,*)   'figure(1)'
!WRITE(23,*)   'subplot(2,1,1)'
WRITE(23,*)   'subplot(1,1,1)'
!WRITE(23,*)   'v = [-3:0.3:3]*1E-8;'   ! For Held-Suarez
!WRITE(23,*)   'contour(xu,yv,q'',v,''k'')'
WRITE(23,*)   'contour(xu,yv,q'',21,''k'')'
WRITE(23,*)   'axis equal'
WRITE(23,*)   'axis ([0,2*pi,-0.5*pi,0.5*pi])'
WRITE(23,*)   'qx = max(max(q));'
WRITE(23,*)   'qn = min(min(q));'
! WRITE(23,*)   'title ([ytitle,''  Min = '',num2str(qn),''  Max = '',num2str(qx) ])'
WRITE(23,*)   'title ([ytitle,''  Min = '',num2str(qn,''%8.2e''),''  Max = '',num2str(qx,''%8.2e'') ])'
!WRITE(23,*)   'pause'

IF (lpdc) THEN

! Borrow div array to diagnose convection flag
div = 0.0D0
DO j = 1, ny
  DO i = 1, nx
    IF ((phi(i,j) - phimem(i,j) > thrconv).and.(philat(i,j) > 0.0D0)) div(i,j) = 1.0D0
  ENDDO
ENDDO

WRITE(23,*)   'ytitle = ''Convection  step ',ystep,''';'
WRITE(23,*)   'q = [ ...'
WRITE(23,888) div
WRITE(23,*)   ' ];'
WRITE(23,*)   'q = reshape(q,nx,ny);'
WRITE(23,*)   'figure(1)'
WRITE(23,*)   'subplot(2,1,2)'
WRITE(23,*)   'contour(xv,yu,q'',[0.1,0.1],''k'')'
WRITE(23,*)   'axis equal'
WRITE(23,*)   'axis ([0,2*pi,-0.5*pi,0.5*pi])'
WRITE(23,*)   'qx = max(max(q));'
WRITE(23,*)   'qn = min(min(q));'
WRITE(23,*)   'title ([ytitle,''  Min = '',num2str(qn),''  Max = '',num2str(qx) ])'
! WRITE(23,*)   'title ([ytitle])'
!WRITE(23,*)   'pause'

!WRITE(23,*)   'ytitle = ''philat  step ',ystep,''';'
!WRITE(23,*)   'q = [ ...'
!WRITE(23,888) philat
!WRITE(23,*)   ' ];'
!WRITE(23,*)   'q = reshape(q,nx,ny);'
!WRITE(23,*)   'figure(1)'
!WRITE(23,*)   'contour(xv,yu,q'',11,''k'')'
!WRITE(23,*)   'qx = max(max(q));'
!WRITE(23,*)   'qn = min(min(q));'
!WRITE(23,*)   'title ([ytitle,''  Min = '',num2str(qn),''  Max = '',num2str(qx) ])'
!WRITE(23,*)   'pause'

!WRITE(23,*)   'ytitle = ''phimem  step ',ystep,''';'
!WRITE(23,*)   'q = [ ...'
!WRITE(23,888) phimem
!WRITE(23,*)   ' ];'
!WRITE(23,*)   'q = reshape(q,nx,ny);'
!WRITE(23,*)   'figure(1)'
!WRITE(23,*)   'contour(xv,yu,q'',11,''k'')'
!WRITE(23,*)   'qx = max(max(q));'
!WRITE(23,*)   'qn = min(min(q));'
!WRITE(23,*)   'title ([ytitle,''  Min = '',num2str(qn),''  Max = '',num2str(qx) ])'
!WRITE(23,*)   'pause'

!WRITE(23,*)   'ytitle = ''phi - phimem  step ',ystep,''';'
!WRITE(23,*)   'q = [ ...'
!WRITE(23,888) phi - phimem
!WRITE(23,*)   ' ];'
!WRITE(23,*)   'q = reshape(q,nx,ny);'
!WRITE(23,*)   'figure(1)'
!WRITE(23,*)   'contour(xv,yu,q'',11,''k'')'
!WRITE(23,*)   'qx = max(max(q));'
!WRITE(23,*)   'qn = min(min(q));'
!WRITE(23,*)   'title ([ytitle,''  Min = '',num2str(qn),''  Max = '',num2str(qx) ])'
!WRITE(23,*)   'pause'

ENDIF




CLOSE(23)

888 FORMAT(E16.4)


END SUBROUTINE outstatem

! ==========================================================

SUBROUTINE outaccm(istep)

! Output accumulated diagnostics for matlab plotting

USE state
USE force
USE accumdiag

IMPLICIT NONE
INTEGER, INTENT(IN) :: istep
CHARACTER*8 :: ystep
CHARACTER*15 :: yname
CHARACTER*24 :: ytitle


INTEGER :: j
REAL*8 :: phibart(nx,ny), phistdd(nx,ny), phizm(ny), phi2zm(ny), phizsd(ny), &
          xibart(nx,ny+1), xistdd(nx,ny+1), xizm(ny+1), xi2zm(ny+1), xizsd(ny+1)


IF (naccum > 0) THEN

! Open file and write grid information
WRITE(ystep,'(I8.8)') istep
yname = 'accum'//ystep//'.m'
OPEN(23,FILE=yname)
WRITE(23,*) 'nx = ',nx,';'
WRITE(23,*) 'ny = ',ny,';'
WRITE(23,*) 'xu = [...'
WRITE(23,888) xu
WRITE(23,*) '];'
WRITE(23,*) 'yu = [...'
WRITE(23,888) yu
WRITE(23,*) '];'
WRITE(23,*) 'xv = [...'
WRITE(23,888) xv
WRITE(23,*) '];'
WRITE(23,*) 'yv = [...'
WRITE(23,888) yv
WRITE(23,*) '];'

! Time mean phi
phibart = accumphi/naccum
DO j = 1, ny
  phizm(j) = SUM(phibart(:,j))/nx
  phi2zm(j) = SUM(phibart(:,j)*phibart(:,j))
ENDDO
phizsd = SQRT((phi2zm - nx*phizm*phizm)/(nx - 1))
WRITE(23,*)   'ytitle = ''Time mean \Phi ',ystep,''';'
WRITE(23,*)   'q = [ ...'
WRITE(23,888) phibart
WRITE(23,*)   ' ];'
WRITE(23,*)   'q = reshape(q,nx,ny);'
WRITE(23,*)   'figure(1)'
WRITE(23,*)   'contour(xv,yu,q'',11,''k'')'
WRITE(23,*)   'qx = max(max(q));'
WRITE(23,*)   'qn = min(min(q));'
WRITE(23,*)   'title ([ytitle,''  Min = '',num2str(qn),''  Max = '',num2str(qx) ])'
WRITE(23,*)   'qzm = [ ...'
WRITE(23,888) phizm
WRITE(23,*)   ' ];'
WRITE(23,*)   'qsd = [ ...'
WRITE(23,888) phizsd
WRITE(23,*)   ' ];'
WRITE(23,*)   'figure(2)'
WRITE(23,*)   'plot(qzm,yu,''k'',qzm-qsd,yu,''k--'',qzm+qsd,yu,''k--'')'
WRITE(23,*)   'ylim([',yv(1),',',yv(ny+1),'])'
WRITE(23,*)   'title (ytitle)'
WRITE(23,*)   'pause'

! Standard deviation of phi
phistdd = SQRT((accumphi2 - naccum*phibart*phibart)/(naccum - 1))
DO j = 1, ny
  phizm(j) = SUM(phistdd(:,j))/nx
  phi2zm(j) = SUM(phistdd(:,j)*phistdd(:,j))
ENDDO
phizsd = SQRT((phi2zm - nx*phizm*phizm)/(nx - 1))
WRITE(23,*)   'ytitle = ''Std dev \Phi ',ystep,''';'
WRITE(23,*)   'q = [ ...'
WRITE(23,888) phistdd
WRITE(23,*)   ' ];'
WRITE(23,*)   'q = reshape(q,nx,ny);'
WRITE(23,*)   'figure(1)'
WRITE(23,*)   'contour(xv,yu,q'',11,''k'')'
WRITE(23,*)   'qx = max(max(q));'
WRITE(23,*)   'qn = min(min(q));'
WRITE(23,*)   'title ([ytitle,''  Min = '',num2str(qn),''  Max = '',num2str(qx) ])'
WRITE(23,*)   'qzm = [ ...'
WRITE(23,888) phizm
WRITE(23,*)   ' ];'
WRITE(23,*)   'qsd = [ ...'
WRITE(23,888) phizsd
WRITE(23,*)   ' ];'
WRITE(23,*)   'figure(2)'
WRITE(23,*)   'plot(qzm,yu,''k'',qzm-qsd,yu,''k--'',qzm+qsd,yu,''k--'')'
WRITE(23,*)   'ylim([',yv(1),',',yv(ny+1),'])'
WRITE(23,*)   'title (ytitle)'
WRITE(23,*)   'pause'

! Time mean div
phibart = accumdiv/naccum
DO j = 1, ny
  phizm(j) = SUM(phibart(:,j))/nx
  phi2zm(j) = SUM(phibart(:,j)*phibart(:,j))
ENDDO
phizsd = SQRT((phi2zm - nx*phizm*phizm)/(nx - 1))
WRITE(23,*)   'ytitle = ''Time mean div ',ystep,''';'
WRITE(23,*)   'q = [ ...'
WRITE(23,888) phibart
WRITE(23,*)   ' ];'
WRITE(23,*)   'q = reshape(q,nx,ny);'
WRITE(23,*)   'figure(1)'
WRITE(23,*)   'contour(xv,yu,q'',11,''k'')'
WRITE(23,*)   'qx = max(max(q));'
WRITE(23,*)   'qn = min(min(q));'
WRITE(23,*)   'title ([ytitle,''  Min = '',num2str(qn),''  Max = '',num2str(qx) ])'
WRITE(23,*)   'qzm = [ ...'
WRITE(23,888) phizm
WRITE(23,*)   ' ];'
WRITE(23,*)   'qsd = [ ...'
WRITE(23,888) phizsd
WRITE(23,*)   ' ];'
WRITE(23,*)   'figure(2)'
WRITE(23,*)   'plot(qzm,yu,''k'',qzm-qsd,yu,''k--'',qzm+qsd,yu,''k--'')'
WRITE(23,*)   'ylim([',yv(1),',',yv(ny+1),'])'
WRITE(23,*)   'title (ytitle)'
WRITE(23,*)   'pause'

! Standard deviation of div
phistdd = SQRT((accumdiv2 - naccum*phibart*phibart)/(naccum - 1))
DO j = 1, ny
  phizm(j) = SUM(phistdd(:,j))/nx
  phi2zm(j) = SUM(phistdd(:,j)*phistdd(:,j))
ENDDO
phizsd = SQRT((phi2zm - nx*phizm*phizm)/(nx - 1))
WRITE(23,*)   'ytitle = ''Std dev div ',ystep,''';'
WRITE(23,*)   'q = [ ...'
WRITE(23,888) phistdd
WRITE(23,*)   ' ];'
WRITE(23,*)   'q = reshape(q,nx,ny);'
WRITE(23,*)   'figure(1)'
WRITE(23,*)   'contour(xv,yu,q'',11,''k'')'
WRITE(23,*)   'qx = max(max(q));'
WRITE(23,*)   'qn = min(min(q));'
WRITE(23,*)   'title ([ytitle,''  Min = '',num2str(qn),''  Max = '',num2str(qx) ])'
WRITE(23,*)   'qzm = [ ...'
WRITE(23,888) phizm
WRITE(23,*)   ' ];'
WRITE(23,*)   'qsd = [ ...'
WRITE(23,888) phizsd
WRITE(23,*)   ' ];'
WRITE(23,*)   'figure(2)'
WRITE(23,*)   'plot(qzm,yu,''k'',qzm-qsd,yu,''k--'',qzm+qsd,yu,''k--'')'
WRITE(23,*)   'ylim([',yv(1),',',yv(ny+1),'])'
WRITE(23,*)   'title (ytitle)'
WRITE(23,*)   'pause'

! Time mean xi
xibart = accumxi/naccum
DO j = 1, ny+1
  xizm(j) = SUM(xibart(:,j))/nx
  xi2zm(j) = SUM(xibart(:,j)*xibart(:,j))
ENDDO
xizsd = SQRT((xi2zm - nx*xizm*xizm)/(nx - 1))
WRITE(23,*)   'ytitle = ''Time mean \xi ',ystep,''';'
WRITE(23,*)   'q = [ ...'
WRITE(23,888) xibart
WRITE(23,*)   ' ];'
WRITE(23,*)   'q = reshape(q,nx,ny+1);'
WRITE(23,*)   'figure(1)'
WRITE(23,*)   'contour(xu,yv,q'',11,''k'')'
WRITE(23,*)   'qx = max(max(q));'
WRITE(23,*)   'qn = min(min(q));'
WRITE(23,*)   'title ([ytitle,''  Min = '',num2str(qn),''  Max = '',num2str(qx) ])'
WRITE(23,*)   'qzm = [ ...'
WRITE(23,888) xizm
WRITE(23,*)   ' ];'
WRITE(23,*)   'qsd = [ ...'
WRITE(23,888) xizsd
WRITE(23,*)   ' ];'
WRITE(23,*)   'figure(2)'
WRITE(23,*)   'plot(qzm,yv,''k'',qzm-qsd,yv,''k--'',qzm+qsd,yv,''k--'')'
WRITE(23,*)   'ylim([',yv(1),',',yv(ny+1),'])'
WRITE(23,*)   'title (ytitle)'
WRITE(23,*)   'pause'

! Standard deviation of xi
xistdd = SQRT((accumxi2 - naccum*xibart*xibart)/(naccum - 1))
DO j = 1, ny+1
  xizm(j) = SUM(xistdd(:,j))/nx
  xi2zm(j) = SUM(xistdd(:,j)*xistdd(:,j))
ENDDO
xizsd = SQRT((xi2zm - nx*xizm*xizm)/(nx - 1))
WRITE(23,*)   'ytitle = ''Std dev \xi ',ystep,''';'
WRITE(23,*)   'q = [ ...'
WRITE(23,888) xistdd
WRITE(23,*)   ' ];'
WRITE(23,*)   'q = reshape(q,nx,ny+1);'
WRITE(23,*)   'figure(1)'
WRITE(23,*)   'contour(xu,yv,q'',11,''k'')'
WRITE(23,*)   'qx = max(max(q));'
WRITE(23,*)   'qn = min(min(q));'
WRITE(23,*)   'title ([ytitle,''  Min = '',num2str(qn),''  Max = '',num2str(qx) ])'
WRITE(23,*)   'qzm = [ ...'
WRITE(23,888) xizm
WRITE(23,*)   ' ];'
WRITE(23,*)   'qsd = [ ...'
WRITE(23,888) xizsd
WRITE(23,*)   ' ];'
WRITE(23,*)   'figure(2)'
WRITE(23,*)   'plot(qzm,yv,''k'',qzm-qsd,yv,''k--'',qzm+qsd,yv,''k--'')'
WRITE(23,*)   'ylim([',yv(1),',',yv(ny+1),'])'
WRITE(23,*)   'title (ytitle)'
WRITE(23,*)   'pause'

IF (lpdc) THEN

! Time mean convection
phibart = accumcon/naccum
DO j = 1, ny
  phizm(j) = SUM(phibart(:,j))/nx
  phi2zm(j) = SUM(phibart(:,j)*phibart(:,j))
ENDDO
phizsd = SQRT((phi2zm - nx*phizm*phizm)/(nx - 1))
WRITE(23,*)   'ytitle = ''Time mean convection ',ystep,''';'
WRITE(23,*)   'q = [ ...'
WRITE(23,888) phibart
WRITE(23,*)   ' ];'
WRITE(23,*)   'q = reshape(q,nx,ny);'
WRITE(23,*)   'figure(1)'
WRITE(23,*)   'contour(xv,yu,q'',11,''k'')'
WRITE(23,*)   'qx = max(max(q));'
WRITE(23,*)   'qn = min(min(q));'
WRITE(23,*)   'title ([ytitle,''  Min = '',num2str(qn),''  Max = '',num2str(qx) ])'
WRITE(23,*)   'qzm = [ ...'
WRITE(23,888) phizm
WRITE(23,*)   ' ];'
WRITE(23,*)   'qsd = [ ...'
WRITE(23,888) phizsd
WRITE(23,*)   ' ];'
WRITE(23,*)   'figure(2)'
WRITE(23,*)   'plot(qzm,yu,''k'',qzm-qsd,yu,''k--'',qzm+qsd,yu,''k--'')'
WRITE(23,*)   'ylim([',yv(1),',',yv(ny+1),'])'
WRITE(23,*)   'title (ytitle)'
WRITE(23,*)   'pause'

ENDIF

! Zonal mean u
DO j = 1, ny
  phizm(j) = SUM(u(:,j))/nx
ENDDO
WRITE(23,*)   'ytitle = ''Zonal mean u ',ystep,''';'
WRITE(23,*)   'qzm = [ ...'
WRITE(23,888) phizm
WRITE(23,*)   ' ];'
WRITE(23,*)   'figure(2)'
WRITE(23,*)   'plot(qzm,yu,''k'')'
WRITE(23,*)   'ylim([',yv(1),',',yv(ny+1),'])'
WRITE(23,*)   'title (ytitle)'
WRITE(23,*)   'pause'


CLOSE(23)

ENDIF

888 FORMAT(E16.4)


END SUBROUTINE outaccm

! ==========================================================

SUBROUTINE outthresh(istep)

! To compute some diagnostics of phi - phimem

USE state
USE thresholds

IMPLICIT NONE

INTEGER, INTENT(IN) :: istep
INTEGER :: i, j, k
REAL*8 :: aa, dphi, totarea



tcount(:) = 0
tarea(:) = 0.0D0
totarea = 0.0d0

DO j = 1, ny
  aa = area(j)
  DO i = 1, nx
    dphi = phi(i,j) - phimem(i,j)
    totarea = totarea + aa
    DO k = 1, 10
      IF (dphi > tvalues(k)) THEN
        tcount(k) = tcount(k) + 1
        tarea(k) = tarea(k) + aa
      ENDIF
    ENDDO
  ENDDO
ENDDO

print *,'Total count ',nx*ny,'    total area ',totarea
print *,'Threshold     Count     Area '
DO k = 1, 10
  WRITE(6,'(F8.1,I9,''      '',E12.5)') tvalues(k), tcount(k), tarea(k)
ENDDO


END SUBROUTINE outthresh

! ==========================================================

SUBROUTINE dump(f,ytitle,nx,ny)

IMPLICIT NONE

INTEGER, INTENT(IN) :: nx, ny
REAL*8, INTENT(IN) :: f(nx,ny)
REAL :: fstar4(nx,ny)
CHARACTER*(*) :: ytitle

! Convert to single precision to reduce size of output
! and improve readability!
fstar4 = f
WRITE(33,*) nx,ny
WRITE(33,*) ytitle
WRITE(33,*) fstar4

END SUBROUTINE dump

! ===============================================================

SUBROUTINE dumpm(f,ytitle,nx,ny)

! Dump for matlab plotting

IMPLICIT NONE

INTEGER, INTENT(IN) :: nx, ny
REAL*8, INTENT(IN) :: f(nx,ny)
REAL :: fstar4(nx,ny)
CHARACTER*(*) :: ytitle

! Convert to single precision to reduce size of output
! and improve readability!
fstar4 = f

WRITE(34,*) 'nx = ',nx,';'
WRITE(34,*) 'ny = ',ny,';'
WRITE(34,*)  ytitle//' = [ ...'
WRITE(34,888) fstar4
WRITE(34,*) ' ];'
WRITE(34,*) '[ min('//ytitle//') max('//ytitle//') ]'
WRITE(34,*) 'z = reshape('//ytitle//',nx,ny);'
WRITE(34,*) 'contour(z'')'
WRITE(34,*) 'title('''//ytitle//''')'
WRITE(34,*) 'pause'

888 FORMAT(E16.4)

END SUBROUTINE dumpm

! ===============================================================

SUBROUTINE outcon

! To finalize and output convergence test diagnostics

USE contest

OPEN(24,FILE='contest.dat')
WRITE(24,*) rmsu
WRITE(24,*) rmsv
WRITE(24,*) rmsphi


END SUBROUTINE outcon


! ======================================================

SUBROUTINE divergence(u,v,div)

! To compute the divergence of a vector field on the C-grid

USE grid
USE constants

IMPLICIT NONE

REAL*8, INTENT(IN) :: u(nx,ny), v(nx,ny+1)
REAL*8, INTENT(OUT) :: div(nx,ny)

INTEGER :: i, j, ip, jp
REAL*8 :: cm, cp, c0, cdx, cdy, rdx, rdy


rdx = rearth*dx
rdy = rearth*dy

DO j = 1, ny
  jp = j + 1
  cm = cosv(j)
  cp = cosv(jp)
  c0 = cosp(j)
  cdx = rdx*c0
  cdy = rdy*c0
  DO i = 1, nx
    ip = i + 1
    IF (i == nx) ip = 1
    div(i,j) = (u(ip,j) - u(i,j))/cdx &
             + (cp*v(i,jp) - cm*v(i,j))/cdy
  ENDDO
ENDDO


END SUBROUTINE divergence

! ======================================================

SUBROUTINE up_outer

! Update terms in outer loop

USE alldata
USE contest

IMPLICIT NONE

INTEGER :: i, j, jp, ng, im, jm, &
           nfixiter, iter

REAL*8 :: div(nx,ny), tempu(nx,ny+1), tempv(nx,ny), &
          sina, cosa, sind, cosd, sinad, cosad, sdl, cdl, &
	  m11, m12, m21, m22, den, tempphi(nx,ny), &
	  divd(nx,ny), &
	  adad(nx,ny), temp1(nx,ny), acorr(nx,ny), &
	  chi(nx,ny), nu, mass1, mass2, r2, cjr2, dyn, &
	  dabar, dxbar, xcorr(nx), dyc, totda, alpha, fac
	  

! Useful constant
r2 = rearth*rearth

! Interpolate current time level terms in momentum equations
! to departure points 
CALL lagrangeu(xdu,ydu,ru0,rud)
CALL lagrangev(xdv,ydv,rv0,rvd)
CALL lagrangeu(xdu,ydu,rv0bar,tempv)
CALL lagrangev(xdv,ydv,ru0bar,tempu)


! Rotate rud, rvd to arrival point local coordinate system

! u points
DO j = 1, ny
  sina = sinp(j)
  cosa = cosp(j)
  DO i = 1, nx
    sind = SIN(ydu(i,j))
    cosd = COS(ydu(i,j))
    sinad = sina*sind
    cosad = cosa*cosd
    sdl = SIN(xu(i) - xdu(i,j))
    cdl = COS(xu(i) - xdu(i,j))
    den = 1.0 + sinad + cosad*cdl
    m11 = (cosad + (1.0 + sinad)*cdl) / den
    m12 = (sina + sind)*sdl / den
    rud(i,j) = rud(i,j)*m11 + tempv(i,j)*m12
  ENDDO
ENDDO
! v points
DO j = 2, ny
  sina = sinv(j)
  cosa = cosv(j)
  DO i = 1, nx
    sind = SIN(ydv(i,j))
    cosd = COS(ydv(i,j))
    sinad = sina*sind
    cosad = cosa*cosd
    sdl = SIN(xv(i) - xdv(i,j))
    cdl = COS(xv(i) - xdv(i,j))
    den = 1.0 + sinad + cosad*cdl
    m22 = (cosad + (1.0 + sinad)*cdl) / den
    m21 = -(sina + sind)*sdl / den   
    rvd(i,j) = tempu(i,j)*m21 + rvd(i,j)*m22
  ENDDO
ENDDO   
! Polar v points
rvd(:,1) = 0.0
rvd(:,ny+1) = 0.0


! Add in orography contribution
rud = rud - adt*dphisdx
rvd = rvd - adt*dphisdy


! Arrival point divergence
CALL divergence(u,v,div)


! Begin alternatives

IF (ischeme == 1) THEN

  ! Semi-Lagrangian a la ENDGame
  
  ! Interpolate current time level terms to departure points
  CALL lagrangep(xdp,ydp,rphi0,rphid)

  ! Total terms from Phi equation
  rphi = rphid + adt*(phiref - phi)*div

ELSEIF (ischeme == 2) THEN

  ! Semi-Lagrangian as close as possible to SLICE
  
  ! Interpolate current time level terms to departure points
  CALL lagrangep(xdp,ydp,rphi0,rphid)

  ! Convert to SLICE-like departure point/cell values
  rphid = rphid / (1.0 + adt*div)

  ! Total terms from Phi equation
  rphi = rphid + adt*phiref*div
  
ELSEIF (ischeme == 3) THEN

  ! SLICE
  
  ! Interpolate current time level terms to departure points
  CALL lagrangep(xdp,ydp,rphi0,tempphi)

  ! Convert to SLICE-like departure point/cell values
  tempphi = tempphi / (1.0 + hdt*div)

  ! Departure point divergence
  CALL lagrangep(xdp,ydp,div0,divd)
  
  ! Current time level divergence associated with modified
  ! winds
  CALL divergence(u0mod,v0mod,div0mod)
  
  ! Arrival point divergence associated with modified winds
  CALL divergence(umod,vmod,divmod)
    
  IF (areafix == 7) THEN
  
    ! Improve estimate of departure area
    
    ! First estimate departure cell divergence
    CALL slice2da(xdumod,ydumod,xdvmod,ydvmod,aread,div0mod,divd)

    ! now step cell area
    DO j = 1, ny
      aread(:,j) = (1.0 - hdt*(divmod(:,j) + divd(:,j)))*area(j)
    ENDDO

    ! Finally advect mass using desired departure areas.
    CALL slice2da(xdumod,ydumod,xdvmod,ydvmod,aread,phi0,rphid)
    
  ELSE
  
    CALL slice2d(xdumod,ydumod,xdvmod,ydvmod,phi0,rphid)
    
  ENDIF


  ! Conservatively merge the SLICE values with the
  ! semi-Lagrangian values
  CALL merge(rphid,tempphi,jmods,jmodn)


  ! Total terms from Phi equation
  rphi = rphid + hdt*phiref*div


ELSE

  PRINT *,'ischeme = ',ischeme,' not implemented '
  STOP
  
ENDIF

! End of alternatives



END SUBROUTINE up_outer

! ======================================================

SUBROUTINE merge(rphi1,rphi2,jmods,jmodn)

! To merge values predicted by SLICE with values predicted
! by semi-Lagrangian advection, using the latter near the
! poles and the former elsewhere.

USE grid

IMPLICIT NONE

REAL*8, INTENT(INOUT) :: rphi1(nx,ny)
REAL*8, INTENT(IN) :: rphi2(nx,ny)
INTEGER, INTENT(IN) :: jmods, jmodn

INTEGER :: j, jj
REAL*8 :: fac, mass1, mass2, w(3)


! Weights for merging
w(1) = 5.0d0/32.0d0
w(2) = 0.5d0
w(3) = 27.0d0/32.0d0


! Southern hemisphere

! Find mass in modified region, correct SL values to
! conserve mass, and overwrite SLICE values
mass1 = 0.0d0
mass2 = 0.0d0
DO j = 1, jmods
  mass1 = mass1 + area(j)*SUM(rphi1(:,j))
  mass2 = mass2 + area(j)*SUM(rphi2(:,j))
ENDDO
fac = mass1/mass2
rphi1(:,1:jmods) = fac*rphi2(:,1:jmods)

! Merging region
DO jj = 1, 3
  j = jmods + jj
  ! Correct mass of SL values
  mass1 = SUM(rphi1(:,j))
  mass2 = SUM(rphi2(:,j))
  fac = mass1/mass2
  ! And merge with SLICE values
  rphi1(:,j) = w(jj)*rphi1(:,j) + fac*(1.0 - w(jj))*rphi2(:,j)
ENDDO


! Northern hemisphere

! Find mass in modified region, correct SL values to
! conserve mass, and overwrite SLICE values
mass1 = 0.0d0
mass2 = 0.0d0
DO j = jmodn, ny
  mass1 = mass1 + area(j)*SUM(rphi1(:,j))
  mass2 = mass2 + area(j)*SUM(rphi2(:,j))
ENDDO
fac = mass1/mass2
rphi1(:,jmodn:ny) = fac*rphi2(:,jmodn:ny)

! Merging region
DO jj = 1, 3
  j = jmodn - jj
  ! Correct mass of SL values
  mass1 = SUM(rphi1(:,j))
  mass2 = SUM(rphi2(:,j))
  fac = mass1/mass2
  ! And merge with SLICE values
  rphi1(:,j) = w(jj)*rphi1(:,j) + fac*(1.0 - w(jj))*rphi2(:,j)
ENDDO


END SUBROUTINE merge

! ======================================================

SUBROUTINE mergearea(area1,area2,jmods,jmodn)

! To merge departure areas predicted by SLICE with those predicted
! by semi-Lagrangian advection, using the latter near the
! poles and the former elsewhere.

USE grid

IMPLICIT NONE

REAL*8, INTENT(INOUT) :: area1(nx,ny)
REAL*8, INTENT(IN) :: area2(nx,ny)
INTEGER, INTENT(IN) :: jmods, jmodn

INTEGER :: j, jj
REAL*8 :: fac, a1, a2, w(3)


! Weights for merging
w(1) = 5.0d0/32.0d0
w(2) = 0.5d0
w(3) = 27.0d0/32.0d0


! Southern hemisphere

! Find area in modified region, correct SL values to
! conserve area, and overwrite SLICE values
a1 = 0.0d0
a2 = 0.0d0
DO j = 1, jmods
  a1 = a1 + SUM(area1(:,j))
  a2 = a2 + SUM(area2(:,j))
ENDDO
fac = a1/a2
area1(:,1:jmods) = fac*area2(:,1:jmods)

! Merging region
DO jj = 1, 3
  j = jmods + jj
  ! Correct area of SL values
  a1 = SUM(area1(:,j))
  a2 = SUM(area2(:,j))
  fac = a1/a2
  ! And merge with SLICE values
  area1(:,j) = w(jj)*area1(:,j) + fac*(1.0 - w(jj))*area2(:,j)
ENDDO


! Northern hemisphere

! Find area in modified region, correct SL values to
! conserve area, and overwrite SLICE values
a1 = 0.0d0
a2 = 0.0d0
DO j = jmodn, ny
  a1 = a1 + SUM(area1(:,j))
  a2 = a2 + SUM(area2(:,j))
ENDDO
fac = a1/a2
area1(:,jmodn:ny) = fac*area2(:,jmodn:ny)

! Merging region
DO jj = 1, 3
  j = jmodn - jj
  ! Correct mass of SL values
  a1 = SUM(area1(:,j))
  a2 = SUM(area2(:,j))
  fac = a1/a2
  ! And merge with SLICE values
  area1(:,j) = w(jj)*area1(:,j) + fac*(1.0 - w(jj))*area2(:,j)
ENDDO


END SUBROUTINE mergearea

! ======================================================


SUBROUTINE up_inner

! Update terms in inner loop

USE alldata
USE contest

IMPLICIT NONE

REAL*8 :: div(nx,ny)


! Update Coriolis terms and add to departure point
! momentum terms
!CALL coriolis(nx,ny,u,v,phi,fu,fv)
CALL coriolis(u,v,phi,fu,fv)
ru = rud + adt*fv
rv = rvd - adt*fu

! Compute divergence
CALL divergence(ru,rv,div)

! RHS of Helmholtz problem
rhs = rphi - phiref*adt*div


END SUBROUTINE up_inner

! ======================================================

SUBROUTINE backsub

! After solving Helmholtz, backsubstitute to find
! updated velocities

USE alldata
USE contest

IMPLICIT NONE

INTEGER :: i, j, im, jm

REAL*8 :: dphidx(nx,ny), dphidy(nx,ny+1), rdx, cdx, rdy

rdx = rearth*dx
rdy = rearth*dy

! Work out dphi/dx and dphi/dy

! At u-points
DO j = 1, ny
  cdx = cosp(j)*rdx
  DO i = 1, nx
    im = i - 1
    IF (i == 1) im = nx
    dphidx(i,j) = (phi(i,j) - phi(im,j))/cdx
  ENDDO
ENDDO

! At v-points
DO j = 2, ny
  jm = j - 1
  DO i = 1, nx 
    dphidy(i,j) = (phi(i,j) - phi(i,jm))/rdy
  ENDDO
ENDDO

! Set polar values to zero
dphidy(:,1) = 0.0
dphidy(:,ny+1) = 0.0

! Backsubstitute
u = ru - adt*dphidx
v = rv - adt*dphidy


END SUBROUTINE backsub

! ======================================================

SUBROUTINE mgsolve(phi,rr,nu,ng)

! Multigrid solver for elliptic equation
!
! Delsq phi - nu phi = rr
!
! using full multigrid algorithm

USE grid
USE constants

IMPLICIT NONE

! Numbers of iterations on coarsest grid and other grids
INTEGER, PARAMETER :: niterc = 10, niter = 2, npass = 2

INTEGER, INTENT(IN) :: ng
REAL*8, INTENT(IN) :: rr(nx,ny), nu
REAL*8, INTENT(OUT) :: phi(nx,ny)

INTEGER :: nnx(ng), nny(ng), igrid, igridm, jgrid, jgridm, iter, &
           nnxj, nnyj, ipass

REAL*8 :: ff(nx,ny,ng),temp1(nx,ny), rf(nx,ny,ng), &
          aa(ny,ng), bb(ny,ng), cc(ny,ng), &
          ddx, ddy, ccp(ny,ng), ccv(ny+1,ng), r2

integer :: j, j2


! Useful constant
r2 = rearth*rearth

! Map coefficients to each grid in the hierarchy
nnx(ng) = nx
nny(ng) = ny
ddx = dx
ddy = dy
ccp(:,ng) = cosp
ccv(:,ng) = cosv

DO j = 1, ny
  aa(j,ng) = (ccv(j+1,ng)/ccp(j,ng))/(r2*ddy*ddy)
  bb(j,ng) = 1.0/(r2*ccp(j,ng)*ccp(j,ng)*ddx*ddx)
  cc(j,ng) = (ccv(j,ng)/ccp(j,ng))/(r2*ddy*ddy)
ENDDO

DO igrid = ng, 2, -1
  igridm = igrid - 1
  nnxj = nnx(igrid)
  nnyj = nny(igrid)
  nnx(igridm) = nnxj/2
  nny(igridm) = nnyj/2
  ddx = 2.0*ddx
  ddy = 2.0*ddy

  DO j = 1, nny(igridm)
    j2 = 2*j
    ccp(j,igridm) = ccv(j2,igrid)
    ccv(j,igridm) = ccv(j2-1,igrid)
  ENDDO
  ccv(nny(igridm)+1,igridm) = 0.0
  
  DO j = 1, nny(igridm)
    aa(j,igridm) = (ccv(j+1,igridm)/ccp(j,igridm))/(r2*ddy*ddy)
    bb(j,igridm) = 1.0/(r2*ccp(j,igridm)*ccp(j,igridm)*ddx*ddx)
    cc(j,igridm) = (ccv(j,igridm)/ccp(j,igridm))/(r2*ddy*ddy)
  ENDDO
ENDDO


! Initialize solution to zero
phi = 0.0d0


DO ipass = 1, npass

! Initialize rhs as residual using latest estimate
IF (ipass == 1) THEN
  ! No need to do the calculation
  rf(:,:,ng) = rr
ELSE
  CALL residual(phi,rr,rf(:,:,ng), &
		  aa(:,ng),bb(:,ng),cc(:,ng),nu,nx,ny,nx,ny)
ENDIF

! Initialize solution to zero
ff = 0.0

! Inject right hand side to each grid in the hierarchy
DO igrid = ng, 2, -1
  igridm = igrid - 1
  nnxj = nnx(igrid)
  nnyj = nny(igrid)
  CALL inject(rf(1,1,igrid),rf(1,1,igridm),ccp(:,igrid),ccp(:,igridm),nx,ny,nnxj,nnyj)
ENDDO


! Iterate to convergence on coarsest grid
nnxj = nnx(1)
nnyj = nny(1)
ff(1:nnxj,1:nnyj,1) = 0.0d0
DO iter = 1, niterc
  CALL relax(ff(1,1,1),rf(1,1,1), &
  aa(:,1),bb(:,1),cc(:,1),nu,nx,ny,nnxj,nnyj)
ENDDO

! Sequence of growing V-cycles
DO igrid = 2, ng

  igridm = igrid - 1
  nnxj = nnx(igrid)
  nnyj = nny(igrid)

  ! Accurately prolong solution to grid igrid
  ! and execute one V-cycle starting from grid igrid

  !Accurately prolong
  CALL prolong2(ff(1,1,igridm),ff(1,1,igrid),nx,ny,nnxj,nnyj)

  ! Descending part of V-cycle
  DO jgrid = igrid, 2, -1
    
    jgridm = jgrid - 1
    nnxj = nnx(jgrid)
    nnyj = nny(jgrid)

    ! Relax on grid jgrid
    DO iter = 1, niter
      CALL relax(ff(1,1,jgrid),rf(1,1,jgrid), &
      aa(:,jgrid),bb(:,jgrid),cc(:,jgrid),nu,nx,ny,nnxj,nnyj)
    ENDDO

    ! Calculate residual on jgrid
    temp1 = 0.0d0
    CALL residual(ff(1,1,jgrid),rf(1,1,jgrid),temp1, &
		  aa(:,jgrid),bb(:,jgrid),cc(:,jgrid),nu,nx,ny,nnxj,nnyj)
    
    ! Inject residual to jgrid-1
    CALL inject(temp1,rf(1,1,jgridm),ccp(:,igrid),ccp(:,igridm),nx,ny,nnxj,nnyj)
   
    ! Set correction first guess to zero on grid jgrid-1
    ff(1:nnx(jgridm),1:nny(jgridm),jgridm) = 0.0d0
    
  ENDDO
  
  ! Relax to convergence on grid 1
  nnxj = nnx(1)
  nnyj = nny(1)
  DO iter = 1, niterc
    CALL relax(ff(1,1,1),rf(1,1,1), &
    aa(:,1),bb(:,1),cc(:,1),nu,nx,ny,nnxj,nnyj)
  ENDDO
  
  ! Ascending part of V-cycle
  DO jgrid = 2, igrid
    
    jgridm = jgrid - 1
    nnxj = nnx(jgrid)
    nnyj = nny(jgrid)
  
    ! Prolong correction to grid jgrid
    CALL prolong(ff(1,1,jgridm),temp1,nx,ny,nnxj,nnyj)
    
    ! Add correction to solution on jgrid
    ff(1:nnxj,1:nnyj,jgrid) = ff(1:nnxj,1:nnyj,jgrid) &
                            + temp1(1:nnxj,1:nnyj)

    ! Relax on grid jgrid
    DO iter = 1, niter
      CALL relax(ff(1,1,jgrid),rf(1,1,jgrid), &
      aa(:,jgrid),bb(:,jgrid),cc(:,jgrid),nu,nx,ny,nnxj,nnyj)
    ENDDO
  ENDDO
  
ENDDO

! Add correction to phi
phi = phi+ ff(1:nx,1:ny,ng)

ENDDO



END SUBROUTINE mgsolve

! ==========================================================

SUBROUTINE inject(ff,cf,ccf,ccc,nx,ny,nnx,nny)

! Inject data from a fine grid to a coarser grid
! using full area weighting for phi

IMPLICIT NONE

INTEGER, INTENT(IN) :: nx, ny, nnx, nny
REAL*8, INTENT(IN) :: ff(nx,ny),ccf(ny),ccc(ny)
REAL*8, INTENT(OUT) :: cf(nx,ny)
INTEGER :: i, i2, j, j2, i2m, i2p, j2m, j2p


DO j = 1, nny/2
  j2 = j + j
  j2m = MODULO(j2-2,nny) + 1
  j2p = MODULO(j2,nny) + 1
  DO i = 1, nnx/2
    i2 = i + i
    i2m = MODULO(i2-2,nnx) + 1
    i2p = MODULO(i2,nnx) + 1
    
    ! Basic version
    cf(i,j) = 0.25*( &
             (ff(i2m,j2m) + ff(i2,j2m))*ccf(j2m) &
	   + (ff(i2m,j2 ) + ff(i2,j2 ))*ccf(j2 ) &
	                                         )/ccc(j)

  ENDDO
ENDDO

END SUBROUTINE inject

! ==========================================================

SUBROUTINE prolong(cf,ff,nx,ny,nnx,nny)

! Prolong phi field from a coarse grid to a fine grid:
! Cheap version using linear fitting

IMPLICIT NONE

INTEGER, INTENT(IN) :: nx, ny, nnx, nny
REAL*8, INTENT(IN) :: cf(nx,ny)
REAL*8, INTENT(OUT) :: ff(nx,ny)
INTEGER :: i, i2, j, j2, im, jm, ip, i2m, jp, j2m, hnnx, hnny

hnnx = nnx/2
hnny = nny/2

DO j = 1, hnny
  j2 = j + j
  jm = MODULO(j-2,hnny) + 1
  jp = MODULO(j,hnny) + 1
  j2m = MODULO(j2-2,nny) + 1
  DO i = 1, hnnx
    i2 = i + i
    im = MODULO(i-2,hnnx) + 1
    ip = MODULO(i,hnnx) + 1
    i2m = MODULO(i2-2,nnx) + 1
    ff(i2m,j2m) = 0.0625*(cf(im,jm) + 3*cf(i,jm) + 3*cf(im,j) + 9*cf(i,j))
    ff(i2m,j2 ) = 0.0625*(cf(im,jp) + 3*cf(i,jp) + 3*cf(im,j) + 9*cf(i,j))
    ff(i2 ,j2m) = 0.0625*(cf(ip,jm) + 3*cf(i,jm) + 3*cf(ip,j) + 9*cf(i,j))
    ff(i2 ,j2 ) = 0.0625*(cf(ip,jp) + 3*cf(i,jp) + 3*cf(ip,j) + 9*cf(i,j))
  ENDDO
ENDDO

END SUBROUTINE prolong

! ==========================================================

SUBROUTINE prolong2(cf,ff,nx,ny,nnx,nny)

! Prolong phi field from a coarse grid to a fine grid:
! Accurate version using cubic fitting

IMPLICIT NONE

REAL*8, PARAMETER :: a = -0.1318, b = 0.8439, c = 0.4575, d = -0.1696

INTEGER, INTENT(IN) :: nx, ny, nnx, nny
REAL*8, INTENT(IN) :: cf(nx,ny)
REAL*8, INTENT(OUT) :: ff(nx,ny)
INTEGER :: hnnx, hnny, i, im, imm, ip, ipp, j, jm, jmm, jp, jpp, &
           i2, i2m, j2, j2m

hnnx = nnx/2
hnny = nny/2

DO j = 1, hnny
  j2 = j + j
  jm = MODULO(j-2,hnny) + 1
  jmm = MODULO(jm-2,hnny) + 1
  jp = MODULO(j,hnny) + 1
  jpp = MODULO(jp,hnny) + 1
  j2m = MODULO(j2-2,nny) + 1
  DO i = 1, hnnx
    i2 = i + i
    im = MODULO(i-2,hnnx) + 1
    imm = MODULO(im-2,hnnx) + 1
    ip = MODULO(i,hnnx) + 1
    ipp = MODULO(ip,hnnx) + 1
    i2m = MODULO(i2-2,nnx) + 1
    
    ff(i2 ,j2 ) = a*(a*cf(im,jm)  + b*cf(im,j)  + c*cf(im,jp)  + d*cf(im,jpp)) &
                + b*(a*cf(i,jm)   + b*cf(i,j)   + c*cf(i,jp)   + d*cf(i,jpp) ) &
		+ c*(a*cf(ip,jm)  + b*cf(ip,j)  + c*cf(ip,jp)  + d*cf(ip,jpp)) &
		+ d*(a*cf(ipp,jm) + b*cf(ipp,j) + c*cf(ipp,jp) + d*cf(ipp,jpp))
		
    ff(i2m,j2 ) = a*(a*cf(ip,jm)  + b*cf(ip,j)  + c*cf(ip,jp)  + d*cf(ip,jpp)) &
                + b*(a*cf(i,jm)   + b*cf(i,j)   + c*cf(i,jp)   + d*cf(i,jpp) ) &
		+ c*(a*cf(im,jm)  + b*cf(im,j)  + c*cf(im,jp)  + d*cf(im,jpp)) &
		+ d*(a*cf(imm,jm) + b*cf(imm,j) + c*cf(imm,jp) + d*cf(imm,jpp))
    
    ff(i2 ,j2m) = a*(a*cf(im,jp)  + b*cf(im,j)  + c*cf(im,jm)  + d*cf(im,jmm)) &
                + b*(a*cf(i,jp)   + b*cf(i,j)   + c*cf(i,jm)   + d*cf(i,jmm) ) &
		+ c*(a*cf(ip,jp)  + b*cf(ip,j)  + c*cf(ip,jm)  + d*cf(ip,jmm)) &
		+ d*(a*cf(ipp,jp) + b*cf(ipp,j) + c*cf(ipp,jm) + d*cf(ipp,jmm))
		
    ff(i2m,j2m) = a*(a*cf(ip,jp)  + b*cf(ip,j)  + c*cf(ip,jm)  + d*cf(ip,jmm)) &
                + b*(a*cf(i,jp)   + b*cf(i,j)   + c*cf(i,jm)   + d*cf(i,jmm) ) &
		+ c*(a*cf(im,jp)  + b*cf(im,j)  + c*cf(im,jm)  + d*cf(im,jmm)) &
		+ d*(a*cf(imm,jp) + b*cf(imm,j) + c*cf(imm,jm) + d*cf(imm,jmm))
    
  ENDDO
ENDDO


END SUBROUTINE prolong2

! ==========================================================

SUBROUTINE relax(ff,rf,a,b,c,nu,nx,ny,nnx,nny)

! Red-black relaxation for elliptic problem

IMPLICIT NONE

INTEGER, INTENT(IN) :: nx, ny, nnx, nny
REAL*8, INTENT(IN) :: rf(nx,ny), a(ny), b(ny), c(ny), nu
REAL*8, INTENT(INOUT) :: ff(nx,ny)

INTEGER :: i, j, im, jm, ip, jp, iparity, istart
INTEGER :: irb = 3 ! 1 for red-black, 0 for Gauss-Seidel, 2 for both
                   ! 3 for row by row simultaneous relaxation
REAL*8 :: c0, xa(nnx), xb(nnx), xc(nnx), xr(nnx), x(nnx)


IF (irb == 1 .OR. irb == 2) THEN

! Red-black relaxation of phi

DO iparity=1, 2

  ! First row
  j = 1
  jp = MODULO(j,nny) + 1
  c0 = a(j) + 2.0*b(j) + nu
  istart = MODULO(j + iparity,2) + 1
  DO i = istart, nnx, 2
    im = MODULO(i - 2,nnx) + 1
    ip = MODULO(i,nnx) + 1
    ff(i,j) = ( b(j)*(ff(ip,j) + ff(im,j)) + a(j)*ff(i,jp)    &
                -rf(i,j) ) &
	      / c0
  ENDDO

  ! Middle rows
  DO j = 2, nny - 1
    jm = MODULO(j - 2,nny) + 1
    jp = MODULO(j,nny) + 1
    c0 = a(j) + 2.0*b(j) + c(j) + nu
    istart = MODULO(j + iparity,2) + 1
    DO i = istart, nnx, 2
      im = MODULO(i - 2,nnx) + 1
      ip = MODULO(i,nnx) + 1
      ff(i,j) = ( b(j)*(ff(ip,j) + ff(im,j)) + a(j)*ff(i,jp) + c(j)*ff(i,jm)    &
                  -rf(i,j) ) &
	        / c0     
    ENDDO
  ENDDO
  
  ! Last row
  j = nny
  jm = MODULO(j - 2,nny) + 1
  c0 = 2.0*b(j) + c(j) + nu
  istart = MODULO(j + iparity,2) + 1
  DO i = istart, nnx, 2
    im = MODULO(i - 2,nnx) + 1
    ip = MODULO(i,nnx) + 1
    ff(i,j) = ( b(j)*(ff(ip,j) + ff(im,j)) + c(j)*ff(i,jm)    &
                -rf(i,j) ) &
	      / c0     
  ENDDO
    
ENDDO 

ENDIF


IF (irb == 0 .OR. irb == 2) THEN

! Gauss-Seidel

  ! First row
  j = 1
  jp = MODULO(j,nny) + 1
  c0 = a(j) + 2.0*b(j) + nu
  DO i = 1, nnx
    im = MODULO(i - 2,nnx) + 1
    ip = MODULO(i,nnx) + 1
    ff(i,j) = ( b(j)*(ff(ip,j) + ff(im,j)) + a(j)*ff(i,jp)    &
                -rf(i,j) ) &
	      / c0
  ENDDO

  ! Middle rows
  DO j = 2, nny - 1
    jm = MODULO(j - 2,nny) + 1
    jp = MODULO(j,nny) + 1
    c0 = a(j) + 2.0*b(j) + c(j) + nu
    DO i = 1, nnx
      im = MODULO(i - 2,nnx) + 1
      ip = MODULO(i,nnx) + 1
      ff(i,j) = ( b(j)*(ff(ip,j) + ff(im,j)) + a(j)*ff(i,jp) + c(j)*ff(i,jm)    &
                  -rf(i,j) ) &
	        / c0     
    ENDDO
  ENDDO
  
  ! Last row
  j = nny
  jm = MODULO(j - 2,nny) + 1
  c0 = 2.0*b(j) + c(j) + nu
  DO i = 1, nnx
    im = MODULO(i - 2,nnx) + 1
    ip = MODULO(i,nnx) + 1
    ff(i,j) = ( b(j)*(ff(ip,j) + ff(im,j)) + c(j)*ff(i,jm)    &
                -rf(i,j) ) &
	      / c0     
  ENDDO

ENDIF

IF (irb == 3) THEN

! Row by row simultaneous relaxation, alternating rows

xa(1:nnx) = b(1)
xc(1:nnx) = b(1)
xb(1:nnx) = -(2*b(1) + a(1) + c(1) + nu)
xr = rf(1:nnx,1) - a(1)*ff(1:nnx,2)
CALL trisolve(x,xa,xb,xc,xr,nnx)
ff(1:nnx,1) = x

DO j = 3, nny-1, 2
  xa(1:nnx) = b(j)
  xc(1:nnx) = b(j)
  xb(1:nnx) = -(2*b(j) + a(j) + c(j) + nu)
  xr = rf(1:nnx,j) - a(j)*ff(1:nnx,j+1) - c(j)*ff(1:nnx,j-1)
  CALL trisolve(x,xa,xb,xc,xr,nnx)
  ff(1:nnx,j) = x
ENDDO

DO j = 2, nny-2, 2
  xa(1:nnx) = b(j)
  xc(1:nnx) = b(j)
  xb(1:nnx) = -(2*b(j) + a(j) + c(j) + nu)
  xr = rf(1:nnx,j) - a(j)*ff(1:nnx,j+1) - c(j)*ff(1:nnx,j-1)
  CALL trisolve(x,xa,xb,xc,xr,nnx)
  ff(1:nnx,j) = x
ENDDO

xa(1:nnx) = b(nny)
xc(1:nnx) = b(nny)
xb(1:nnx) = -(2*b(nny) + a(nny) + c(nny) + nu)
xr = rf(1:nnx,nny) - c(nny)*ff(1:nnx,nny-1)
CALL trisolve(x,xa,xb,xc,xr,nnx)
ff(1:nnx,nny) = x

ENDIF


END SUBROUTINE relax

! ==========================================================

SUBROUTINE residual(ff,rf,resf,a,b,c,nu,nx,ny,nnx,nny)

! Calculate residual for elliptic problem

IMPLICIT NONE

INTEGER, INTENT(IN) :: nx, ny, nnx, nny
REAL*8, INTENT(IN) :: ff(nx,ny), rf(nx,ny), a(ny), b(ny), c(ny), nu 
REAL*8, INTENT(OUT) :: resf(nx,ny)

INTEGER :: i, j, im, jm, ip, jp
REAL*8 :: c0


! Residual in phi equation

! First row
j = 1
jp = MODULO(j,nny) + 1
c0 = a(j) + 2.0*b(j) + nu
DO i = 1, nnx
  im = MODULO(i - 2,nnx) + 1
  ip = MODULO(i,nnx) + 1
  resf(i,j) = ( rf(i,j) &
             - b(j)*(ff(ip,j) + ff(im,j)) - a(j)*ff(i,jp)   &
	     + c0*ff(i,j) )
ENDDO

! Middle rows
DO j = 2, nny - 1
  jm = MODULO(j - 2,nny) + 1
  jp = MODULO(j,nny) + 1
  c0 = a(j) + 2.0*b(j) + c(j) + nu
  DO i = 1, nnx
    im = MODULO(i - 2,nnx) + 1
    ip = MODULO(i,nnx) + 1
    resf(i,j) = ( rf(i,j) &
               - b(j)*(ff(ip,j) + ff(im,j)) - a(j)*ff(i,jp) - c(j)*ff(i,jm)   &
	       + c0*ff(i,j) )
  ENDDO
ENDDO
  
! Last row
j = nny
jm = MODULO(j - 2,nny) + 1
c0 = 2.0*b(j) + c(j) + nu
DO i = 1, nnx
  im = MODULO(i - 2,nnx) + 1
  ip = MODULO(i,nnx) + 1
  resf(i,j) = ( rf(i,j) &
             - b(j)*(ff(ip,j) + ff(im,j)) - c(j)*ff(i,jm)   &
	     + c0*ff(i,j) )
ENDDO


END SUBROUTINE residual

! ==========================================================

SUBROUTINE lagrangeu(xd,yd,q,qnew)

! 2D Cubic Lagrange interpolation for the u-points on the
! Spherical C-grid

USE grid
USE util

IMPLICIT NONE

REAL*8, INTENT(IN) :: xd(nx,ny), yd(nx,ny), &
		      q(nx,ny)
REAL*8, INTENT(OUT) :: qnew(nx,ny)
INTEGER :: i, j, id, idm, idp, idpp, jd, jdm, jdp, jdpp, &
           i1, i1m, i1p, i1pp, hnx
REAL*8 :: xdd, ydd, d12, d13, d14, d23, d24, d34, &
          denx1, denx2, denx3, denx4, deny1, deny2, deny3, deny4, &
	  dd1, dd2, dd3, dd4, &
	  fac1, fac2, fac3, fac4, facy1, facy2, facy3, facy4, &
	  q1, q2, q3, q4, qqm, qq, qqp, qqpp, &
	  flip, x(nx), y(ny)

! Regular gridded values are at u points
x = xu
y = yu

! Handy value for handling poles
hnx = nx/2

! Interpolation factors can be pre-computed on uniform grid

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

    ! Find indices of departure point and stencil 
    xdd = xd(i,j)
    ydd = yd(i,j)
  
    id =   MODULO(FLOOR((xdd-x(1))/dx),nx) + 1
    idm =  MODULO(id-2,nx) + 1
    idp =  MODULO(id,nx) + 1
    idpp = MODULO(idp,nx) + 1
    jd =   FLOOR((ydd-y(1))/dy) + 1
    jdm =  jd - 1
    jdp =  jd + 1
    jdpp = jdp + 1
    
    ! Factors for x-interpolation
    dd1 = near(xdd - x(idm), twopi)
    dd2 = near(xdd - x(id), twopi)
    dd3 = near(xdd - x(idp), twopi)
    dd4 = near(xdd - x(idpp), twopi)    
    fac1 = dd2*dd3*dd4/denx1
    fac2 = dd1*dd3*dd4/denx2
    fac3 = dd1*dd2*dd4/denx3
    fac4 = dd1*dd2*dd3/denx4
    
    ! Factors for y-interpolation
    dd2 = ydd - (y(1) + jdm*dy)
    dd1 = dd2 + dy
    dd3 = dd2 - dy
    dd4 = dd2 - 2*dy
    facy1 = dd2*dd3*dd4/deny1
    facy2 = dd1*dd3*dd4/deny2
    facy3 = dd1*dd2*dd4/deny3
    facy4 = dd1*dd2*dd3/deny4    
    
    ! Interpolate at four rows
    ! First
    i1m = idm
    i1 = id
    i1p = idp
    i1pp = idpp
    flip = 1.0
    IF (jd .le. 1) THEN
      jdm = 1 - jdm
      i1m  = MODULO(i1m - hnx - 1, nx) + 1
      i1   = MODULO(i1 - hnx - 1, nx) + 1
      i1p  = MODULO(i1p - hnx - 1, nx) + 1
      i1pp = MODULO(i1pp - hnx - 1, nx) + 1
      flip = -1.0
    ENDIF
    q1 = q(i1m ,jdm)
    q2 = q(i1  ,jdm)
    q3 = q(i1p ,jdm)
    q4 = q(i1pp,jdm)    
    qqm = (q1*fac1 + q2*fac2 + q3*fac3 + q4*fac4)*flip

    ! Second
    i1m = idm
    i1 = id
    i1p = idp
    i1pp = idpp
    flip = 1.0
    IF (jd .eq. 0) THEN
      jd = 1
      i1m  = MODULO(i1m - hnx - 1, nx) + 1
      i1   = MODULO(i1 - hnx - 1, nx) + 1
      i1p  = MODULO(i1p - hnx - 1, nx) + 1
      i1pp = MODULO(i1pp - hnx - 1, nx) + 1
      flip = -1.0
    ENDIF
    q1 = q(i1m ,jd)
    q2 = q(i1  ,jd)
    q3 = q(i1p ,jd)
    q4 = q(i1pp,jd)    
    qq = (q1*fac1 + q2*fac2 + q3*fac3 + q4*fac4)*flip

    ! Third
    i1m = idm
    i1 = id
    i1p = idp
    i1pp = idpp
    flip = 1.0
    IF (jd .eq. ny) THEN
      jdp = ny
      i1m  = MODULO(i1m - hnx - 1, nx) + 1
      i1   = MODULO(i1 - hnx - 1, nx) + 1
      i1p  = MODULO(i1p - hnx - 1, nx) + 1
      i1pp = MODULO(i1pp - hnx - 1, nx) + 1
      flip = -1.0
    ENDIF
    q1 = q(i1m ,jdp)
    q2 = q(i1  ,jdp)
    q3 = q(i1p ,jdp)
    q4 = q(i1pp,jdp)    
    qqp = (q1*fac1 + q2*fac2 + q3*fac3 + q4*fac4)*flip

    ! Fourth
    i1m = idm
    i1 = id
    i1p = idp
    i1pp = idpp
    flip = 1.0
    IF (jd .ge. ny - 1) THEN
      jdpp = 2*ny + 1 - jdpp
      i1m  = MODULO(i1m - hnx - 1, nx) + 1
      i1   = MODULO(i1 - hnx - 1, nx) + 1
      i1p  = MODULO(i1p - hnx - 1, nx) + 1
      i1pp = MODULO(i1pp - hnx - 1, nx) + 1
      flip = -1.0
    ENDIF
    q1 = q(i1m ,jdpp)
    q2 = q(i1  ,jdpp)
    q3 = q(i1p ,jdpp)
    q4 = q(i1pp,jdpp)
    qqpp = (q1*fac1 + q2*fac2 + q3*fac3 + q4*fac4)*flip

    ! Interpolate in y
    qnew(i,j) = qqm*facy1 + qq*facy2 + qqp*facy3 + qqpp*facy4

  ENDDO
ENDDO


END SUBROUTINE lagrangeu

! ==========================================================

SUBROUTINE lagrangev(xd,yd,q,qnew)

! 2D Cubic Lagrange interpolation for the v-points on the
! Spherical C-grid

! Polar values are used in interpolation, but new polar
! values are not interpolated as they are not needed later.
! Array incies still run from 1 to ny+1 to avoid confusion

USE grid
USE util

IMPLICIT NONE

REAL*8, INTENT(IN) :: xd(nx,ny+1), yd(nx,ny+1), &
		      q(nx,ny+1)
REAL*8, INTENT(OUT) :: qnew(nx,ny+1)
INTEGER :: i, j, id, idm, idp, idpp, jd, jdm, jdp, jdpp, &
           i1, i1m, i1p, i1pp, hnx
REAL*8 :: xdd, ydd, d12, d13, d14, d23, d24, d34, &
          denx1, denx2, denx3, denx4, deny1, deny2, deny3, deny4, &
	  dd1, dd2, dd3, dd4, &
	  fac1, fac2, fac3, fac4, facy1, facy2, facy3, facy4, &
	  q1, q2, q3, q4, qqm, qq, qqp, qqpp, &
	  flip, x(nx), y(ny+1)


! Regular gridded values are at v points
x = xv
y = yv

! Handy value for handling poles
hnx = nx/2

! Interpolation factors can be pre-computed on uniform grid

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

    ! Find indices of departure point and stencil 
    xdd = xd(i,j)
    ydd = yd(i,j)
  
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

    
    ! Factors for x-interpolation    
    dd1 = near(xdd - x(idm), twopi)
    dd2 = near(xdd - x(id), twopi)
    dd3 = near(xdd - x(idp), twopi)
    dd4 = near(xdd - x(idpp), twopi)    
    fac1 = dd2*dd3*dd4/denx1
    fac2 = dd1*dd3*dd4/denx2
    fac3 = dd1*dd2*dd4/denx3
    fac4 = dd1*dd2*dd3/denx4
    
    ! Factors for y-interpolation  
    dd2 = ydd - y(jd)
    dd1 = dd2 + dy
    dd3 = dd2 - dy
    dd4 = dd2 - 2*dy   
    facy1 = dd2*dd3*dd4/deny1
    facy2 = dd1*dd3*dd4/deny2
    facy3 = dd1*dd2*dd4/deny3
    facy4 = dd1*dd2*dd3/deny4    
    
    ! Interpolate at four rows
    ! First
    i1m = idm
    i1 = id
    i1p = idp
    i1pp = idpp
    flip = 1.0
    IF (jd .eq. 1) THEN
      jdm = 2
      i1m  = MODULO(i1m - hnx - 1, nx) + 1
      i1   = MODULO(i1 - hnx - 1, nx) + 1
      i1p  = MODULO(i1p - hnx - 1, nx) + 1
      i1pp = MODULO(i1pp - hnx - 1, nx) + 1
      flip = -1.0
    ENDIF
    q1 = q(i1m ,jdm)
    q2 = q(i1  ,jdm)
    q3 = q(i1p ,jdm)
    q4 = q(i1pp,jdm)    
    qqm = (q1*fac1 + q2*fac2 + q3*fac3 + q4*fac4)*flip

    ! Second
    i1m = idm
    i1 = id
    i1p = idp
    i1pp = idpp
    q1 = q(i1m ,jd)
    q2 = q(i1  ,jd)
    q3 = q(i1p ,jd)
    q4 = q(i1pp,jd)    
    qq = q1*fac1 + q2*fac2 + q3*fac3 + q4*fac4

    ! Third
    i1m = idm
    i1 = id
    i1p = idp
    i1pp = idpp
    q1 = q(i1m ,jdp)
    q2 = q(i1  ,jdp)
    q3 = q(i1p ,jdp)
    q4 = q(i1pp,jdp)    
    qqp = q1*fac1 + q2*fac2 + q3*fac3 + q4*fac4

    ! Fourth
    i1m = idm
    i1 = id
    i1p = idp
    i1pp = idpp
    flip = 1.0
    IF (jd .eq. ny) THEN
      jdpp = ny
      i1m  = MODULO(i1m - hnx - 1, nx) + 1
      i1   = MODULO(i1 - hnx - 1, nx) + 1
      i1p  = MODULO(i1p - hnx - 1, nx) + 1
      i1pp = MODULO(i1pp - hnx - 1, nx) + 1
      flip = -1.0
    ENDIF
    q1 = q(i1m ,jdpp)
    q2 = q(i1  ,jdpp)
    q3 = q(i1p ,jdpp)
    q4 = q(i1pp,jdpp)
    qqpp = (q1*fac1 + q2*fac2 + q3*fac3 + q4*fac4)*flip

    ! Interpolate in y
    qnew(i,j) = qqm*facy1 + qq*facy2 + qqp*facy3 + qqpp*facy4

  ENDDO
ENDDO


END SUBROUTINE lagrangev

! ==========================================================

SUBROUTINE gradphis

! Calculate the components of the gradient of the orography

USE alldata

IMPLICIT NONE

INTEGER :: i, im, j, jm
REAL*8 :: rdx, rdy, cdx

rdx = rearth*dx
rdy = rearth*dy

! At u-points
DO j = 1, ny
  cdx = cosp(j)*rdx
  DO i = 1, nx
    im = i - 1
    IF (i == 1) im = nx
    dphisdx(i,j) = (phis(i,j) - phis(im,j))/cdx
  ENDDO
ENDDO

! At v-points
DO j = 2, ny
  jm = j - 1
  DO i = 1, nx 
    dphisdy(i,j) = (phis(i,j) - phis(i,jm))/rdy
  ENDDO
ENDDO

! Set polar valuesto zero
dphisdy(:,1) = 0.0
dphisdy(:,ny+1) = 0.0


END SUBROUTINE gradphis

! ==========================================================

SUBROUTINE uvatphi(nx,ny,u,v,up,vp)

! To average u and v to phi points
! on the C-grid

IMPLICIT NONE

INTEGER, INTENT(IN) :: nx, ny
REAL*8, INTENT(IN) :: u(nx,ny), v(nx,ny+1)
REAL*8, INTENT(OUT) :: up(nx,ny), vp(nx,ny)

INTEGER :: i, ip, j, jp

! Assume polar values of v have been defined by a previous
! call to polar, e.g. within cgridave

DO j = 1, ny
  jp = j+1
  DO i = 1, nx
    ip = i+1
    IF ( i == nx ) ip = 1
    up(i,j) = 0.5*(u(i,j) + u(ip,j))
    vp(i,j) = 0.5*(v(i,j) + v(i,jp))
  ENDDO
ENDDO

END SUBROUTINE uvatphi

! ========================================================

SUBROUTINE departurefgp

! First guess for departure point calculation at phi points.
! Use velocity at old time level.

USE alldata

IMPLICIT NONE
INTEGER :: i, j
REAL*8 :: sina, cosa, x, y, r, sind, dlambda

DO  j = 1, ny
  sina = sinp(j)
  cosa = cosp(j)
  DO i = 1, nx
    ! Displacement in local Cartesian system
    x = -u0p(i,j)*dt
    y = -v0p(i,j)*dt
    ! Project back to spherical coordinate system
    r = SQRT(x*x + y*y + rearth*rearth)
    sind = (y*cosa + rearth*sina)/r
    ydp(i,j) = ASIN(sind)
    dlambda = ATAN2(x,rearth*cosa - y*sina)
    xdp(i,j) = MODULO(xp(i) + dlambda, twopi)
  ENDDO
ENDDO


END SUBROUTINE departurefgp

! ========================================================

SUBROUTINE departurep

! Calculate departure points for phi points for a given estimate
! of the current and new u and v


USE alldata

IMPLICIT NONE

INTEGER :: ndepit = 2, idepit, i, j, k, l, kp, lp, k1, k1p, k2, k2p, &
           hnx

REAL*8 :: a1, a2, b1, b2, ud, vd, sina, cosa, sind, cosd, rk, rl, &
          flip1, flip2, sinad, cosad, sdl, cdl, den, urot, vrot, &
	  m11, m12, m21, m22, x, y, r, dlambda


! ---------------------------------------------------------

! Handy quantity for polar interpolation
hnx = nx/2

! Loop over iterations
DO idepit = 1, ndepit

  ! phi-point departure points
  DO  j = 1, ny
    sina = sinp(j)
    cosa = cosp(j)
    DO i = 1, nx
    
      ! Trig factors at estimated departure point
      sind = SIN(ydp(i,j))
      cosd = COS(ydp(i,j))
    
      ! Determine departure cell index based on latest estimate
      ! xdp, ydp
      rk = (xdp(i,j) - xp(1))/dx
      k = FLOOR(rk)
      a2 = rk - k
      a1 = 1.0 - a2
      k = MODULO(k, nx) + 1
      rl = (ydp(i,j) - yp(1))/dy
      l = FLOOR(rl)
      b2 = rl - l
      b1 = 1.0 - b2
      l = l + 1
      kp = k+1
      IF (k == nx) kp = 1
      lp = l+1
      
      ! Tricks to handle polar case
      ! Here we interpolate across the pole.
      ! (An alternative would be to use the polar values of u and v
      ! and interpolate between nearest u latitude and the pole.)
      k1 = k
      k1p = kp
      k2 = k
      k2p = kp
      flip1 = 1.0
      flip2 = 1.0
      IF (l == 0) THEN ! South pole
        l = 1
        k1 = MODULO(k1 + hnx - 1, nx) + 1
        k1p = MODULO(k1p + hnx - 1, nx ) + 1
        flip1 = -1.0
      ELSEIF(l == ny) THEN ! North pole
        lp = ny
        k2 = MODULO(k2 + hnx - 1, nx) + 1
        k2p = MODULO(k2p + hnx - 1, nx ) + 1
        flip2 = -1.0
      ENDIF
      
      ! Linearly interpolate  velocity to estimated departure point
      ud = (a1*b1*u0p(k1,l) &
          + a2*b1*u0p(k1p,l))*flip1 &
	 + (a1*b2*u0p(k2,lp) &
	  + a2*b2*u0p(k2p,lp))*flip2
      vd = (a1*b1*v0p(k1,l) &
          + a2*b1*v0p(k1p,l))*flip1 &
	 + (a1*b2*v0p(k2,lp) &
	  + a2*b2*v0p(k2p,lp))*flip2
	 
      ! Rotate to arrival point Cartesian system
      sinad = sina*sind
      cosad = cosa*cosd
      sdl = SIN(xp(i) - xdp(i,j))
      cdl = COS(xp(i) - xdp(i,j))
      den = 1.0 + sinad + cosad*cdl
      m11 = (cosad + (1.0 + sinad)*cdl) / den
      m12 = (sina + sind)*sdl / den
      m21 = -m12
      m22 = m11
      urot = ud*m11 + vd*m12
      vrot = ud*m21 + vd*m22 

      ! Hence calculate better estimate of departure point
      ! in arrival point Cartesian system
      x = -hdt*(up(i,j) + urot)
      y = -hdt*(vp(i,j) + vrot)

      ! Project back to spherical coordinate system
      r = SQRT(x*x + y*y + rearth*rearth)
      sind = (y*cosa + rearth*sina)/r
      ydp(i,j) = ASIN(sind)
      dlambda = ATAN2(x,rearth*cosa - y*sina)
      xdp(i,j) = MODULO(xp(i) + dlambda, twopi)
    ENDDO
  ENDDO

ENDDO

END SUBROUTINE departurep

! ========================================================

SUBROUTINE lagrangep(xd,yd,q,qnew)

! 2D Cubic Lagrange interpolation for the phi-points on the
! Spherical C-grid

USE grid
USE util

IMPLICIT NONE

REAL*8, INTENT(IN) :: xd(nx,ny), yd(nx,ny), &
		      q(nx,ny)
REAL*8, INTENT(OUT) :: qnew(nx,ny)
INTEGER :: i, j, id, idm, idp, idpp, jd, jdm, jdp, jdpp, &
           i1, i1m, i1p, i1pp, hnx
REAL*8 :: xdd, ydd, d12, d13, d14, d23, d24, d34, &
          denx1, denx2, denx3, denx4, deny1, deny2, deny3, deny4, &
	  dd1, dd2, dd3, dd4, &
	  fac1, fac2, fac3, fac4, facy1, facy2, facy3, facy4, &
	  q1, q2, q3, q4, qqm, qq, qqp, qqpp, &
	  x(nx), y(ny)


! Regular gridded values are at phi points
x = xp
y = yp

! Handy value for handling poles
hnx = nx/2

! Interpolation factors can be pre-computed on uniform grid

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

    ! Find indices of departure point and stencil 
    xdd = xd(i,j)
    ydd = yd(i,j)
  
    id =   MODULO(FLOOR((xdd-x(1))/dx),nx) + 1
    idm =  MODULO(id-2,nx) + 1
    idp =  MODULO(id,nx) + 1
    idpp = MODULO(idp,nx) + 1
    jd =   FLOOR((ydd-y(1))/dy) + 1
    jdm =  jd - 1
    jdp =  jd + 1
    jdpp = jdp + 1
    
    ! Factors for x-interpolation
    dd1 = near(xdd - x(idm), twopi)
    dd2 = near(xdd - x(id), twopi)
    dd3 = near(xdd - x(idp), twopi)
    dd4 = near(xdd - x(idpp), twopi)    
    fac1 = dd2*dd3*dd4/denx1
    fac2 = dd1*dd3*dd4/denx2
    fac3 = dd1*dd2*dd4/denx3
    fac4 = dd1*dd2*dd3/denx4
    
    ! Factors for y-interpolation
    dd2 = ydd - (y(1) + jdm*dy)
    dd1 = dd2 + dy
    dd3 = dd2 - dy
    dd4 = dd2 - 2*dy
    facy1 = dd2*dd3*dd4/deny1
    facy2 = dd1*dd3*dd4/deny2
    facy3 = dd1*dd2*dd4/deny3
    facy4 = dd1*dd2*dd3/deny4        
    
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
    q1 = q(i1m ,jdm)
    q2 = q(i1  ,jdm)
    q3 = q(i1p ,jdm)
    q4 = q(i1pp,jdm)    
    qqm = (q1*fac1 + q2*fac2 + q3*fac3 + q4*fac4)

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
    q1 = q(i1m ,jd)
    q2 = q(i1  ,jd)
    q3 = q(i1p ,jd)
    q4 = q(i1pp,jd)    
    qq = (q1*fac1 + q2*fac2 + q3*fac3 + q4*fac4)

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
    q1 = q(i1m ,jdp)
    q2 = q(i1  ,jdp)
    q3 = q(i1p ,jdp)
    q4 = q(i1pp,jdp)    
    qqp = (q1*fac1 + q2*fac2 + q3*fac3 + q4*fac4)

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
    q1 = q(i1m ,jdpp)
    q2 = q(i1  ,jdpp)
    q3 = q(i1p ,jdpp)
    q4 = q(i1pp,jdpp)
    qqpp = (q1*fac1 + q2*fac2 + q3*fac3 + q4*fac4)

    ! Interpolate in y
    qnew(i,j) = qqm*facy1 + qq*facy2 + qqp*facy3 + qqpp*facy4

  ENDDO
ENDDO


END SUBROUTINE lagrangep

! ==========================================================

SUBROUTINE slphieqn

! Calulate terms needed in the phi equation at departure points
! when semi-Lagrangian advection of phi is used

USE alldata

IMPLICIT NONE


! Assume polar values of v0 have been defined by a call to polar,
! e.g. within cgridave

! First find flow divergence
CALL divergence(u0,v0,div0)

! Hence compute required term
rphi0 = phi0*(1.0 - bdt*div0)


END SUBROUTINE slphieqn

! ==========================================================

SUBROUTINE modifydep

! Compute modified departure points such that the
! poles don't move, smoothly merging with true
! departure points at lower latitudes

USE alldata

IMPLICIT NONE

INTEGER :: jmod, i, j, jp
REAL*8 :: margin, ydpole, ymod, w1, w2, frac, corr, xxa, xxd, yya, yyd


! width of buffer zone
margin = 4*dy

! Copy departure points
xdumod = xdu
ydumod = ydu
xdvmod = xdv
ydvmod = ydv


! Southern hemisphere

! Estimate departure latitude of the pole
ydpole = MAXVAL(ydu(:,1))

! Add buffer zone and calculate index to this latitude
ymod = ydpole + margin
jmod = CEILING((ymod+piby2)/dy) + 1
jmods = jmod

! Modify departure points
DO j = 1, jmod
  jp = j + 1
  DO i = 1,nx
    ! Transform to cartesian coodinates based on pole,
    ! merge, and transform back
    frac = MAX((ymod - yu(j))/(ymod + piby2), 0.0d0)
    w1 = COS(piby2*frac)**2
    w2 = 1.0 - w1
    xxa = (yu(j) + piby2)*COS(xu(i))
    yya = (yu(j) + piby2)*SIN(xu(i))
    xxd = (ydumod(i,j) + piby2)*COS(xdumod(i,j))
    yyd = (ydumod(i,j) + piby2)*SIN(xdumod(i,j))
    xxd = w1*xxd + w2*xxa
    yyd = w1*yyd + w2*yya
    xdumod(i,j) = MODULO(ATAN2(yyd,xxd), twopi)
    ydumod(i,j)  = SQRT(xxd*xxd + yyd*yyd) - piby2
    frac = MAX((ymod - yv(jp))/(ymod + piby2), 0.0d0)
    w1 = COS(piby2*frac)**2
    w2 = 1.0 - w1
    xxa = (yv(jp) + piby2)*COS(xv(i))
    yya = (yv(jp) + piby2)*SIN(xv(i))
    xxd = (ydvmod(i,jp) + piby2)*COS(xdvmod(i,jp))
    yyd = (ydvmod(i,jp) + piby2)*SIN(xdvmod(i,jp))
    xxd = w1*xxd + w2*xxa
    yyd = w1*yyd + w2*yya
    xdvmod(i,jp) = MODULO(ATAN2(yyd,xxd),twopi)
    ydvmod(i,jp)  = SQRT(xxd*xxd + yyd*yyd) - piby2
  ENDDO
ENDDO


! Northern hemisphere

! Estimate departure latitude of the pole
ydpole = MINVAL(ydu(:,ny))

! Add buffer zone and calculate index to this latitude
ymod = ydpole - margin
jmod = FLOOR((ymod+piby2)/dy)
jmodn = jmod


! Modify departure points
DO j = jmod, ny
  DO i = 1,nx
    frac = MAX((ymod - yu(j))/(ymod - piby2), 0.0d0)
    w1 = COS(piby2*frac)**2
    w2 = 1.0 - w1
    xxa = (piby2 - yu(j))*COS(xu(i))
    yya = (piby2 - yu(j))*SIN(xu(i))
    xxd = (piby2 - ydumod(i,j))*COS(xdumod(i,j))
    yyd = (piby2 - ydumod(i,j))*SIN(xdumod(i,j))
    xxd = w1*xxd + w2*xxa
    yyd = w1*yyd + w2*yya
    xdumod(i,j) = MODULO(ATAN2(yyd,xxd),twopi)
    ydumod(i,j)  = piby2 - SQRT(xxd*xxd + yyd*yyd)
    frac = MAX((ymod - yv(j))/(ymod - piby2), 0.0d0)
    w1 = COS(piby2*frac)**2
    w2 = 1.0 - w1
    xxa = (piby2 - yv(j))*COS(xv(i))
    yya = (piby2 - yv(j))*SIN(xv(i))
    xxd = (piby2 - ydvmod(i,j))*COS(xdvmod(i,j))
    yyd = (piby2 - ydvmod(i,j))*SIN(xdvmod(i,j))
    xxd = w1*xxd + w2*xxa
    yyd = w1*yyd + w2*yya
    xdvmod(i,j) = MODULO(ATAN2(yyd,xxd),twopi)
    ydvmod(i,j)  = piby2 - SQRT(xxd*xxd + yyd*yyd)
  ENDDO
ENDDO


! Assign polar departure points to equal arrival points
xdvmod(:,1) = xv(:)
ydvmod(:,1) = -piby2
xdvmod(:,ny+1) = xv(:)
ydvmod(:,ny+1) = piby2


END SUBROUTINE modifydep

! ==========================================================

SUBROUTINE modifywind

! Compute modified wind field such that the
! poles don't move, smoothly merging with true
! wind field at lower latitudes

USE alldata

IMPLICIT NONE

INTEGER :: jmod, i, j, jp
REAL*8 :: margin, ydpole, ymod, w1, frac


! width of buffer zone
margin = 4*dy

! Copy winds
umod = u
vmod = v


! Southern hemisphere

! Estimate departure latitude of the pole
ydpole = MAXVAL(ydu(:,1))

! Add buffer zone and calculate index to this latitude
ymod = ydpole + margin
jmod = CEILING((ymod+piby2)/dy) + 1
jmods = jmod

! Modify winds
DO j = 1, jmod
  jp = j + 1
  DO i = 1,nx
    
    frac = MAX((ymod - yu(j))/(ymod + piby2), 0.0d0)
    w1 = COS(piby2*frac)**2
    umod(i,j) = w1*u(i,j)
    frac = MAX((ymod - yv(jp))/(ymod + piby2), 0.0d0)
    w1 = COS(piby2*frac)**2
    vmod(i,jp) = w1*v(i,jp)
    
  ENDDO
ENDDO


! Northern hemisphere

! Estimate departure latitude of the pole
ydpole = MINVAL(ydu(:,ny))

! Add buffer zone and calculate index to this latitude
ymod = ydpole - margin
jmod = FLOOR((ymod+piby2)/dy)
jmodn = jmod


! Modify departure points
DO j = jmod, ny
  DO i = 1,nx
  
    frac = MAX((ymod - yu(j))/(ymod - piby2), 0.0d0)
    w1 = COS(piby2*frac)**2
    umod(i,j) = w1*u(i,j)
    frac = MAX((ymod - yv(j))/(ymod - piby2), 0.0d0)
    w1 = COS(piby2*frac)**2
    vmod(i,j) = w1*v(i,j)

  ENDDO
ENDDO


! Assign polar winds to equal zero
vmod(:,1) = 0.0
vmod(:,ny+1) = 0.0



END SUBROUTINE modifywind

! ==========================================================

SUBROUTINE dumpgrid

! Output the grid coordinates in a simple format for use in generating
! reference solutions

USE grid

IMPLICIT NONE
INTEGER :: i, j
REAL*8 :: long, lat
CHARACTER*30 :: ygridcoords
CHARACTER*5 :: ynlon, ynlat

! -----------------------------------------------------------------------

WRITE(ynlon,'(I5.5)') nx
WRITE(ynlat,'(I5.5)') ny
IF (ABS(rotgrid) < 0.00001D0) THEN
  ygridcoords = 'gridcoords_ll__'//ynlon//'x'//ynlat//'.dat'
ELSE
  ygridcoords = 'gridcoords_llr_'//ynlon//'x'//ynlat//'.dat'
ENDIF
OPEN(88,FILE=ygridcoords,FORM='UNFORMATTED')

WRITE(88) nx*ny
DO j = 1, ny
  DO i = 1, nx
    WRITE(88) geolonp(i,j), geolatp(i,j)
  ENDDO
ENDDO

CLOSE(88)

! -----------------------------------------------------------------------

END SUBROUTINE dumpgrid

! =======================================================================

SUBROUTINE forcing

! Add forcing terms for a forced-dissipative shallow water analogue
! of Held-Suarez

USE timeinfo
USE constants
USE state
USE work
USE force
USE errdiag

IMPLICIT NONE
INTEGER :: i, im, j, jm
LOGICAL :: lhitzero
REAL*8 :: fac, rdx, rdy, cdx, dphi, dphifac, div(nx,ny), pp, ppbound


! -----------------------------------------------------------------------

! Relaxation of phi for dry dynamics
fac = 1.0d0/(1.0d0 + dt/tauh)
phi = (phi + (dt/tauh)*phieqm)*fac

! Relaxation of u and v for dry dynamics
fac = 1.0d0/(1.0d0 + dt/tauv)
u = (u + (dt/tauv)*ueqm)*fac
v = (v + (dt/tauv)*veqm)*fac

print *,'Forcing terms added'


IF (lpdc) THEN

  ! phimem provides memory of recent phi values
  ! relax towards phi
  fac = 1.0d0/(1.0d0 + dt/taum)
  phimem = (phimem + (dt/taum)*phi)*fac

  ! relax `latent energy' towards equilibrium
  fac = 1.0d0/(1.0d0 + dt/taul)
  philat = (philat + (dt/taul)*phileq)*fac

  phiphys = 0.0D0

  print *,'Max phi-phim = ',MAXVAL(phi-phimem), &
          '  Range phil = ',MINVAL(philat),MAXVAL(philat)

  ! Detect where convection is happening
  ccount = 0
  DO j = 1, ny
    DO i = 1, nx
      dphi = phi(i,j) - phimem(i,j)
      IF ((dphi > thrconv).and.(philat(i,j) > 0.0D0)) THEN
        ccount = ccount + 1

        ! Convection is occurring
        ! First guess for phiphys field
        pp = thrconv - dphi - betal*philat(i,j)

        ! Sink of latent energy; bound phiphys if needed
        lhitzero = .FALSE.
        dphifac = (phi(i,j) - phi0(i,j))/phi0(i,j)
        IF (dphifac > 0.0D0) THEN
          ppbound = -philat(i,j)/dphifac
          IF (pp <= ppbound) THEN
            lhitzero = .TRUE.
            pp = ppbound
          ENDIF
        ENDIF
        philat(i,j) = philat(i,j) + MIN(0.0D0,pp*dphifac)

        ! If convection has taken philat to zero then set it
        ! to a negative value to allow time to `recharge'
        IF (lhitzero) philat(i,j) = philrecharge

        ! And finalize phiphys
        phiphys(i,j) = pp

      ENDIF
    ENDDO
  ENDDO
  ! ccount is saved in module errdiag for output to diagnostics file
  print *,'Convection at ',ccount,' grid points'

  ! Now take grad of phiphys field and use to increment momentum
  ! using simple forward Euler step

  rdx = rearth*dx/dt
  rdy = rearth*dy/dt

  ! At u-points
  DO j = 1, ny
    cdx = cosp(j)*rdx
    DO i = 1, nx
      im = i - 1
      IF (i == 1) im = nx
      u(i,j) = u(i,j) - (phiphys(i,j) - phiphys(im,j))/cdx
    ENDDO
  ENDDO

  ! At v-points
  DO j = 2, ny
    jm = j - 1
    DO i = 1, nx 
      v(i,j) = v(i,j) - (phiphys(i,j) - phiphys(i,jm))/rdy
    ENDDO
  ENDDO

  print *,'Physics terms added'


ENDIF


! -----------------------------------------------------------------------

END SUBROUTINE forcing

! =======================================================================

SUBROUTINE writeref

USE timeinfo
USE state

IMPLICIT NONE

INTEGER, PARAMETER :: &
ngref = 16,            &  ! Number of grids on which to dump reference solution
nreftime = 5              ! Number of times at which reference solution is required

INTEGER :: ilist, igref, nface, if0
CHARACTER*10 :: ytime
CHARACTER*31 :: ygc(ngref)
CHARACTER*22 :: yrefpre(ngref)
REAL*8 :: reftime(nreftime)
REAL*8, ALLOCATABLE :: flong(:), flat(:), href(:)

! ----------------------------------------------------------

! List of grid coordinate files
ygc( 1) = 'gridcoords_hex_0000000642.dat '
ygc( 2) = 'gridcoords_hex_0000002562.dat '
ygc( 3) = 'gridcoords_hex_0000010242.dat '
ygc( 4) = 'gridcoords_hex_0000040962.dat '
ygc( 5) = 'gridcoords_cube_0000000864.dat'
ygc( 6) = 'gridcoords_cube_0000003456.dat'
ygc( 7) = 'gridcoords_cube_0000013824.dat'
ygc( 8) = 'gridcoords_cube_0000055296.dat'
ygc( 9) = 'gridcoords_ll__00040x00020.dat'
ygc(10) = 'gridcoords_ll__00080x00040.dat'
ygc(11) = 'gridcoords_ll__00160x00080.dat'
ygc(12) = 'gridcoords_ll__00320x00160.dat'
ygc(13) = 'gridcoords_llr_00040x00020.dat'
ygc(14) = 'gridcoords_llr_00080x00040.dat'
ygc(15) = 'gridcoords_llr_00160x00080.dat'
ygc(16) = 'gridcoords_llr_00320x00160.dat'

! List of reference solution file prefixes
yrefpre( 1) = 'TC5ref_hex__0000000642'
yrefpre( 2) = 'TC5ref_hex__0000002562'
yrefpre( 3) = 'TC5ref_hex__0000010242'
yrefpre( 4) = 'TC5ref_hex__0000040962'
yrefpre( 5) = 'TC5ref_cube_0000000864'
yrefpre( 6) = 'TC5ref_cube_0000003456'
yrefpre( 7) = 'TC5ref_cube_0000013824'
yrefpre( 8) = 'TC5ref_cube_0000055296'
yrefpre( 9) = 'TC5ref_ll__00040x00020'
yrefpre(10) = 'TC5ref_ll__00080x00040'
yrefpre(11) = 'TC5ref_ll__00160x00080'
yrefpre(12) = 'TC5ref_ll__00320x00160'
yrefpre(13) = 'TC5ref_llr_00040x00020'
yrefpre(14) = 'TC5ref_llr_00080x00040'
yrefpre(15) = 'TC5ref_llr_00160x00080'
yrefpre(16) = 'TC5ref_llr_00320x00160'

! List of times at which reference solution is required
reftime(1) = 3600.0d0              ! 1 hour
reftime(2) = 86400.0d0             ! 1 day
reftime(3) = 5.0d0*86400.0         ! 5 days
reftime(4) = 10.0D0*86400.0        ! 10 days
reftime(5) = 15.0d0*86400.0        ! 15 days

DO ilist = 1, nreftime

  IF (istep == NINT(reftime(ilist)/dt)) THEN

    ! Reference solution is required at this time
    WRITE(ytime,'(I10.10)') NINT(reftime(ilist))

    DO igref = 1, ngref

      ! Read in grid on which reference solution is required
      print *,'opening ',ygc(igref)
      OPEN(48,FILE=ygc(igref),FORM='UNFORMATTED')
      READ(48) nface
      print *,'nface = ',nface
      ALLOCATE(flong(nface),flat(nface),href(nface))
      DO if0 = 1, nface
        READ(48) flong(if0), flat(if0)
      ENDDO
      print *,'coords read '
      CLOSE(48)

      ! Interpolate reference solution to this grid
      CALL interpref(flong,flat,href,nface)

      ! Output reference solution
      print *,'Creating reference solution file ',yrefpre(igref)//'_'//ytime//'.dat'
      OPEN(49,FILE=yrefpre(igref)//'_'//ytime//'.dat',FORM='UNFORMATTED')
      DO if0 = 1, nface
        WRITE(49) href(if0)
      ENDDO
      CLOSE(49)

      DEALLOCATE(flong,flat,href)

    ENDDO

  ENDIF

ENDDO

! ---------------------------------------------------------

END SUBROUTINE writeref

! ==========================================================

SUBROUTINE interpref(xd,yd,href,nface)

! 2D Cubic Lagrange interpolation to compute the surface height
! at a given list of coordinates - used for generating reference solutions
! for other models

USE grid
USE state
USE util

IMPLICIT NONE

INTEGER, INTENT(IN) :: nface
REAL*8, INTENT(IN) :: xd(nface), yd(nface)
REAL*8, INTENT(OUT) :: href(nface)
INTEGER :: i, j, id, idm, idp, idpp, jd, jdm, jdp, jdpp, &
           i1, i1m, i1p, i1pp, hnx
REAL*8 :: xdd, ydd, d12, d13, d14, d23, d24, d34, &
          denx1, denx2, denx3, denx4, deny1, deny2, deny3, deny4, &
	  dd1, dd2, dd3, dd4, &
	  fac1, fac2, fac3, fac4, facy1, facy2, facy3, facy4, &
	  q1, q2, q3, q4, qqm, qq, qqp, qqpp, &
	  x(nx), y(ny), q(nx,ny)


! Construct surface height
q = (phi + phis)/9.80616d0

! Regular gridded values are at phi points
x = xp
y = yp

! Handy value for handling poles
hnx = nx/2

! Interpolation factors can be pre-computed on uniform grid

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


  DO i = 1, nface

    ! Find indices of departure point and stencil 
    xdd = xd(i)
    ydd = yd(i)
  
    id =   MODULO(FLOOR((xdd-x(1))/dx),nx) + 1
    idm =  MODULO(id-2,nx) + 1
    idp =  MODULO(id,nx) + 1
    idpp = MODULO(idp,nx) + 1
    jd =   FLOOR((ydd-y(1))/dy) + 1
    jdm =  jd - 1
    jdp =  jd + 1
    jdpp = jdp + 1
    
    ! Factors for x-interpolation
    dd1 = near(xdd - x(idm), twopi)
    dd2 = near(xdd - x(id), twopi)
    dd3 = near(xdd - x(idp), twopi)
    dd4 = near(xdd - x(idpp), twopi)    
    fac1 = dd2*dd3*dd4/denx1
    fac2 = dd1*dd3*dd4/denx2
    fac3 = dd1*dd2*dd4/denx3
    fac4 = dd1*dd2*dd3/denx4
    
    ! Factors for y-interpolation
    dd2 = ydd - (y(1) + jdm*dy)
    dd1 = dd2 + dy
    dd3 = dd2 - dy
    dd4 = dd2 - 2*dy
    facy1 = dd2*dd3*dd4/deny1
    facy2 = dd1*dd3*dd4/deny2
    facy3 = dd1*dd2*dd4/deny3
    facy4 = dd1*dd2*dd3/deny4        
    
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
    q1 = q(i1m ,jdm)
    q2 = q(i1  ,jdm)
    q3 = q(i1p ,jdm)
    q4 = q(i1pp,jdm)    
    qqm = (q1*fac1 + q2*fac2 + q3*fac3 + q4*fac4)

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
    q1 = q(i1m ,jd)
    q2 = q(i1  ,jd)
    q3 = q(i1p ,jd)
    q4 = q(i1pp,jd)    
    qq = (q1*fac1 + q2*fac2 + q3*fac3 + q4*fac4)

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
    q1 = q(i1m ,jdp)
    q2 = q(i1  ,jdp)
    q3 = q(i1p ,jdp)
    q4 = q(i1pp,jdp)    
    qqp = (q1*fac1 + q2*fac2 + q3*fac3 + q4*fac4)

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
    q1 = q(i1m ,jdpp)
    q2 = q(i1  ,jdpp)
    q3 = q(i1p ,jdpp)
    q4 = q(i1pp,jdpp)
    qqpp = (q1*fac1 + q2*fac2 + q3*fac3 + q4*fac4)

    ! Interpolate in y
    href(i) = qqm*facy1 + qq*facy2 + qqp*facy3 + qqpp*facy4

  ENDDO



END SUBROUTINE interpref

! ==========================================================

SUBROUTINE diffref

! Compute and output the surface height difference from a reference solution

USE state
USE timeinfo
USE constants
IMPLICIT NONE

INTEGER, PARAMETER :: nreftime = 5, chanerrout = 86
INTEGER :: reftime(nreftime), ilist, i, j
REAL*8 :: h(nx,ny), href(nx,ny), l1, l2, linf, tota
CHARACTER*11 :: yrefpre
CHARACTER*10 :: ytime
CHARACTER*11 :: yres
CHARACTER*17 :: yname
CHARACTER*5 :: ynx, yny

! -----------------------------------------------------------------------

! Reference solution file prefixes
IF (ABS(rotgrid) < 0.00001D0) THEN
  yrefpre = 'TC5ref_ll__'
ELSE
  yrefpre = 'TC5ref_llr_'
ENDIF

! Resolution (for building file name)
WRITE(ynx,'(I5.5)') nx
WRITE(yny,'(I5.5)') ny
yres = ynx//'x'//yny

! List of times at which reference solution is required
reftime(1) = NINT(3600.0d0)              ! 1 hour
reftime(2) = NINT(86400.0d0)             ! 1 day
reftime(3) = NINT(5.0d0*86400.0)         ! 5 days
reftime(4) = NINT(10.0D0*86400.0)        ! 10 days
reftime(5) = NINT(15.0d0*86400.0)        ! 15 days

! Surface height
h = (phi + phis)/gravity

DO ilist = 1, nreftime

  IF (NINT(istep*dt) == reftime(ilist)) THEN

    ! Reference solution is required at this time
    WRITE(ytime,'(I10.10)') reftime(ilist)
    ! Read reference solution
    OPEN(87,FILE=yrefpre//yres//'_'//ytime//'.dat',FORM='UNFORMATTED')
    DO j = 1, ny
      DO i = 1, nx
        READ(87) href(i,j)
      ENDDO
    ENDDO
    CLOSE(87)

    ! File for Matlab primal grid output
    yname = 'err1_'//ytime//'.m'
    OPEN(chanerrout,FILE=yname)

    ! Write header information
    WRITE(chanerrout,*)   'nx = ',nx,';'
    WRITE(chanerrout,*)   'ny = ',ny,';'
    WRITE(chanerrout,*)   'long = [ ...'
    WRITE(chanerrout,888) xp
    WRITE(chanerrout,*)   ' ];'
    WRITE(chanerrout,*)   'lat = [ ...'
    WRITE(chanerrout,888) yp
    WRITE(chanerrout,*)   ' ];'

    ! Output
    WRITE(chanerrout,*)   'ytitle = ''Reference h time ',ytime,''';'
    WRITE(chanerrout,*)   'q = [ ...'
    WRITE(chanerrout,888) href
    WRITE(chanerrout,*)   ' ];'
    WRITE(chanerrout,*)   'q = reshape(q,nx,ny);'
    WRITE(chanerrout,*)   'figure(1)'
    WRITE(chanerrout,*)   'contour(long,lat,q'',11,''k'')'
    WRITE(chanerrout,*)   'qx = max(max(q));'
    WRITE(chanerrout,*)   'qn = min(min(q));'
    WRITE(chanerrout,*)   'title ([ytitle,''  Min = '',num2str(qn),''  Max = '',num2str(qx) ])'
    WRITE(chanerrout,*)   'ytitle = ''Solution h time ',ytime,''';'
    WRITE(chanerrout,*)   'q = [ ...'
    WRITE(chanerrout,888) h
    WRITE(chanerrout,*)   ' ];'
    WRITE(chanerrout,*)   'q = reshape(q,nx,ny);'
    WRITE(chanerrout,*)   'figure(2)'
    WRITE(chanerrout,*)   'contour(long,lat,q'',11,''k'')'
    WRITE(chanerrout,*)   'qx = max(max(q));'
    WRITE(chanerrout,*)   'qn = min(min(q));'
    WRITE(chanerrout,*)   'title ([ytitle,''  Min = '',num2str(qn),''  Max = '',num2str(qx) ])'
    h = h - href
    WRITE(chanerrout,*)   'ytitle = ''h error time ',ytime,''';'
    WRITE(chanerrout,*)   'q = [ ...'
    WRITE(chanerrout,888) h
    WRITE(chanerrout,*)   ' ];'
    WRITE(chanerrout,*)   'q = reshape(q,nx,ny);'
    WRITE(chanerrout,*)   'figure(3)'
    WRITE(chanerrout,*)   'contour(long,lat,q'',11,''k'')'
    WRITE(chanerrout,*)   'qx = max(max(q));'
    WRITE(chanerrout,*)   'qn = min(min(q));'
    WRITE(chanerrout,*)   'title ([ytitle,''  Min = '',num2str(qn),''  Max = '',num2str(qx) ])'
    CLOSE(chanerrout)

    ! Error norms
    tota = 0.0d0
    l1 = 0.0d0
    l2 = 0.0d0
    DO j = 1, ny
      tota = tota + area(j)
      l1 = l1 + area(j)*SUM(ABS(h(:,j)))
      l2 = l2 + area(j)*SUM(h(:,j)*h(:,j))
    ENDDO
    tota = tota*nx
    l1 = l1/tota
    l2 = SQRT(l2/tota)
    linf = MAXVAL(ABS(h))
    print *,'h error: L1   = ',l1 ,'  L2   = ',l2 ,'  Linf = ',linf

  ENDIF

ENDDO

888 FORMAT(E16.4)
889 FORMAT(I8)


! -----------------------------------------------------------------------

END SUBROUTINE diffref

! =======================================================================

SUBROUTINE testdep

! Test departure point calculation

USE alldata

IMPLICIT NONE

INTEGER :: i, j

! Old velocity field specified in initial
! Copy into new velocity field
! Calculate C-grid average velocities,
! and find first guess for departure points
CALL prelim


! Full calculation of departure points
CALL departure
CALL departurep


! Modify departure points so that poles don't move
! for use with SLICE
CALL modifydep
xdu = xdumod
ydu = ydumod
xdv = xdvmod
ydv = ydvmod




WRITE(22,*) nx, ny
DO j = 1, ny
  DO i = 1, nx
    WRITE(22,*) xdu(i,j), ydu(i,j), xu(i), yu(j)
  ENDDO
ENDDO
DO j = 2, ny
  DO i = 1, nx
    WRITE(22,*) xdv(i,j), ydv(i,j), xv(i), yv(j)
  ENDDO
ENDDO
DO j = 1, ny
  DO i = 1, nx
    WRITE(22,*) xdp(i,j), ydp(i,j), xp(i), yp(j)
  ENDDO
ENDDO


END SUBROUTINE testdep

! =========================================================

SUBROUTINE testrotate

! Test SL advection / interpolation

USE alldata

IMPLICIT NONE

INTEGER :: i, j

REAL*8 :: tempu(nx,ny+1), tempv(nx,ny), sina, cosa, sinad, cosad, &
          sdl, cdl, den, m11, m12, m21, m22, sind, cosd


! Old velocity field specified in initial
! Copy into new velocity field
! Calculate C-grid average velocities,
! and find first guess for departure points
CALL prelim

! Full calculation of departure points
CALL departure
CALL departurep

! Interpolate current time level terms in momentum equations
! to departure points 
CALL lagrangeu(xdu,ydu,u0,rud)
CALL lagrangev(xdv,ydv,v0,rvd)

! Rotate rud, rvd to arrival point local coordinate system
print *,' need to interpolate ave fields, not ave interpolated'
CALL cgridave(nx,ny,rud,rvd,tempu,tempv)

call dump(u0,' u ',nx,ny)
call dump(rud,' rdu before',nx,ny)
call dump(tempu,' tempu ',nx,ny+1)
call dump(v0,' v ',nx,ny+1)
call dump(rvd,' rdv before',nx,ny+1)
call dump(tempv,' tempv ',nx,ny)

! u points
DO j = 1, ny
  sina = sinp(j)
  cosa = cosp(j)
  DO i = 1, nx
    sind = SIN(ydu(i,j))
    cosd = COS(ydu(i,j))
    sinad = sina*sind
    cosad = cosa*cosd
    sdl = SIN(xu(i) - xdu(i,j))
    cdl = COS(xu(i) - xdu(i,j))
    den = 1.0 + sinad + cosad*cdl
    m11 = (cosad + (1.0 + sinad)*cdl) / den
    m12 = (sina + sind)*sdl / den
    rud(i,j) = rud(i,j)*m11 + tempv(i,j)*m12
  ENDDO
ENDDO
! v points
DO j = 2, ny
  sina = sinv(j)
  cosa = cosv(j)
  DO i = 1, nx
    sind = SIN(ydv(i,j))
    cosd = COS(ydv(i,j))
    sinad = sina*sind
    cosad = cosa*cosd
    sdl = SIN(xv(i) - xdv(i,j))
    cdl = COS(xv(i) - xdv(i,j))
    den = 1.0 + sinad + cosad*cdl
    m22 = (cosad + (1.0 + sinad)*cdl) / den
    m21 = -(sina + sind)*sdl / den
    if (j == 18 .and. i == 1) then
      print *,'j = ',j
      print *,'sina = ',sina,' cosa = ',cosa
      print *,'sind = ',sind,' cosd = ',cosd
      print *,'sdl = ',sdl,' cdl = ',cdl
      print *,'den = ',den,' m22 m21 = ',m22,m21
      print *,' tempu = ',tempu(i,j),' rvd = ',rvd(i,j)
    endif
    rvd(i,j) = tempu(i,j)*m21 + rvd(i,j)*m22
    if (j == 18 .and. i == 1) then
      print *,' rvd = ',rvd(i,j)
    endif
  ENDDO
ENDDO   
! Polar v points
rvd(:,1) = 0.0
rvd(:,ny+1) = 0.0

call dump(rud,' rdu after ',nx,ny)
call dump(rud - u0,' u diff ',nx,ny)

call dump(rvd,' rdv after ',nx,ny+1)
call dump(rvd - v0,' v diff ',nx,ny+1)



END SUBROUTINE testrotate

! =========================================================


