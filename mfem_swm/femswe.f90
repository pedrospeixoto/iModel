MODULE runtype

  ! Default values are set here; they may be overridden
  ! via the namelist

  IMPLICIT NONE

  ! File containing grid information
  !CHARACTER*31 :: ygridfile = 'gridopermap_hex_0000010242.dat'
  !CHARACTER*31 :: ygridfile = 'gridopermap_cube_0000013824.dat'
  CHARACTER*128 :: ygridfile = 'gridopermap_cube_0000221184.dat'

  ! File to output grid coordinates (for use in generating reference
  ! solutions).
  ! CHARACTER*30 :: ygridcoords = 'gridcoords_hex_0000010242.dat '
  !CHARACTER*30 :: ygridcoords = 'gridcoords_cube_0000013824.dat'
  !CHARACTER*30 :: ygridcoords = 'gridcoords_cube_0000000216.dat'
  CHARACTER*128 :: ygridcoords = 'gridcoords_cube_0000221184.dat'

  ! Run identifier for naming output files
  CHARACTER*6 :: runid = '000001'

  ! True for a restart run
  LOGICAL :: lrestart = .FALSE.

  ! Name of restart file
  CHARACTER*128 :: yresfile = 'run000002_restart_00000720.dat'

  ! File containing grid information
  CHARACTER*128 :: parameterFile = 'swenml.in'

  !Initial condition
  INTEGER :: ic=0
  CHARACTER*16 :: icname="0"

  INTEGER :: suppressOutput = 0

END MODULE runtype

! ================================================================

MODULE constants

  ! Various physical and geometrical constants

  IMPLICIT NONE

  ! Pi
  REAL*8, PARAMETER :: pi = 3.14159265358979323d0, piby2 = 0.5d0*pi

  ! Earth's radius
  REAL*8 :: rearth = 6371220.0d0

  ! Gravitational acceleration
  REAL*8 :: gravity = 9.80616d0

  ! Rotation rate of the earth
  REAL*8 :: rotatn = 7.29212d-5

  ! Dual cell integrals of planetary vorticity
  REAL*8, ALLOCATABLE :: planvort2(:)

  !Degrees to radians coversion (multiply to obtain conversion)
  REAL*8, PARAMETER :: deg2rad = pi / 180.0d0

  !Radians to Degrees coversion (multiply to obtain conversion)
  REAL*8, PARAMETER :: rad2deg = 1.0d0/deg2rad

END MODULE constants

! ================================================================

MODULE state

  ! Variables describing the model state.
  ! Note on naming convention: an appended 2 indicates a variable
  ! integrated over a cell or dual cell area (a discrete 2-form);
  ! an appended 1 indicates a variable integrated along a
  ! primal or dual edge (a discrete 1-form). Names with neither
  ! indicate point values or length or area averages.

  USE grid
  IMPLICIT NONE

  ! OROGRAPHY
  ! orog2             Primal cell area integral of orography
  REAL*8, ALLOCATABLE :: orog2(:)

  ! PROGNOSTIC VARIABLES
  ! phi2              Primal cell area integral of geopotential
  ! u1                Flux across primal edge
  ! c2                Area integral of tracer concentration
  REAL*8, ALLOCATABLE :: phi2(:), u1(:), c2(:)


  ! VARIABLES USED IN CHECKING THEORETICAL PROPERTIES
  ! xphibar2          Prognostic dual cell integral of geopotential
  ! xzeta2            Prognostic pv-like tracer
  REAL*8, ALLOCATABLE :: xphibar2(:), xzeta2(:)


END MODULE state

! ================================================================

MODULE work

  ! Variables used in various parts of the calculation, 
  ! particularly the advection scheme, which it is therefore
  ! useful to store

  IMPLICIT NONE

  ! v1                Circulation along dual edge
  ! vbar1             Time averaged circulation along dual edge
  ! ubar1             Time averaged flux across primal edge
  ! vperpbar1         Time averaged flux across dual edge
  ! uperpbar1         Time averaged circulation along primal edge
  ! mf1               Time integrated mass flux across primal edge
  ! mfperp1           Time integrated mass flux across dual edge
  ! divfac            Factor related to divergence used in modifying swept areas
  REAL*8, ALLOCATABLE :: v1(:), vbar1(:), ubar1(:), vperpbar1(:), uperpbar1(:), &
       mf1(:), mfperp1(:), divfac(:)


  REAL*8, ALLOCATABLE :: b2(:), pva(:), u1_temp(:), &
       b0(:), gb1(:), div2(:), zeta2(:), &
       pv(:), qperp1(:), hqperp1(:), divmf(:), phirefe(:), &
       b2_temp(:), phi_inc(:), u_inc(:), &
       rphi2(:), ru1(:), rhs(:), &
       phi2_new(:), u1_new(:), mru(:), &
       phibar2(:), mu(:), wu(:), wf(:)


END MODULE work

! ================================================================

MODULE errdiag

  ! Variables used in computing error diagnostics

  IMPLICIT NONE

  ! phi2_init           Initial geopotential
  ! u1_init             Initial u
  REAL*8, ALLOCATABLE :: phi2_init(:), u1_init(:), &
       phiexac(:), pvexac(:), &
       phierr(:), pverr(:), u1_err(:)

END MODULE errdiag

! ================================================================

MODULE helmcoeff

  ! Reference state values on the grid hierarchy needed for
  ! multigrid helmholtz solver

  USE grid
  IMPLICIT NONE

  ! phiref         Primal grid cell values of reference geopotential
  !                on finest grid
  ! nusq           Cell edge values of reference geopotential
  !                times alpha_pg*alpha_v*dt*dt
  !                on all grids
  ! helmdiag       Diagonal coefficient of Helmholtz operator
  !                on all grids
  ! underrel       Under-relaxation parameter on all grids
  REAL*8, ALLOCATABLE :: phiref(:), nusq(:,:), helmdiag(:,:), underrel(:)


END MODULE helmcoeff

! ================================================================

MODULE advection

  ! Fields related to advection scheme
  ! Default values are set here; they may be overridden
  ! via the namelist

  USE grid
  IMPLICIT NONE

  ! degree         Degree of polynomial fit
  ! monotone       Use limiter if true (NOT YET IMPLEMENTED)
  ! nmonomial      Number of monomials to build polynomials of
  !                given degree
  ! ngauss         Number of Gauss points needed in each direction
  !                to integrate fit over swept area
  ! nstenadvf      Number of cells in stencil for given face
  ! nexadvf        Number of cells in stencil to be fitted exactly
  ! nadvfmx        Maximum number of cells in advection stencil
  !                for given face
  ! stenadvf       Stencil for constructing polynomial fit for
  !                given face
  ! intmonf        Integrals of monomials over all cells of stencil
  !                of given face
  ! nstenadvv      Number of dual cells in stencil for given vertex
  ! nexadvv        Number of dual cells in stencil to be fitted exactly
  ! nadvvmx        Maximum number of dual cells in advection stencil
  !                for given vertex
  ! stenadvv       Stencil for constructing polynomial fit for
  !                given vertex
  ! intmonv        Integrals of monomials over all dual cells of stencil
  !                of given vertex
  ! xgauss         Gauss points for the interval [0,1]
  ! wgauss         Gauss weights for the interval [0,1]

  INTEGER :: degree = 2
  LOGICAL :: monotone = .FALSE.
  INTEGER :: nmonomial, ngauss
  INTEGER :: nadvfmx, nadvvmx
  INTEGER, ALLOCATABLE :: nstenadvf(:), nexadvf(:), &
       stenadvf(:,:), &
       nstenadvv(:), nexadvv(:), &
       stenadvv(:,:)
  REAL*8, ALLOCATABLE :: intmonf(:,:,:), intmonv(:,:,:)
  REAL*8 :: xgauss(3), wgauss(3)


END MODULE advection

! ================================================================

MODULE timestep

  ! Information related to timestepping
  ! Default values are set here; they may be overridden
  ! via the namelist

  IMPLICIT NONE

  ! Time step
  REAL*8 :: dt = 300.0d0

  ! Off-centring parameters
  ! alpha_v, beta_v         For the velocity used for advection
  ! alpha_pg, beta_pg       For the pressure gradient term               
  REAL*8 :: alpha_v = 0.5d0, alpha_pg = 0.5d0
  REAL*8 :: beta_v, beta_pg

  ! Number of iterations in nonlinear solver
  INTEGER :: niter = 4

  ! Length of integration (in steps)
  INTEGER :: nstop = 1

  ! Frequency of output dumps (in steps)
  INTEGER :: noutput = 1

  ! Frequency of restart dumps (in steps)
  INTEGER :: nrestart = 1

  ! Frequency of diagnostics
  INTEGER :: ndiagnostics = 1

  ! Current time step
  INTEGER :: istep

  ! Current time
  REAL*8 :: time


END MODULE timestep

! ================================================================

MODULE channels

  ! Tidy list of all I/O channels in one place to avoid accidental
  ! overuse of any channel number

  IMPLICIT NONE

  INTEGER, PARAMETER :: channml = 20          ! For reading namelists
  INTEGER, PARAMETER :: changrid = 25         ! Grid information
  INTEGER, PARAMETER :: chanerr = 42          ! Time series of basic error measures
  INTEGER, PARAMETER :: chandiag = 43         ! Time series of basic diagnostics
  INTEGER, PARAMETER :: chanrefgrd = 26       ! To dump grid coordinates for reference solutions
  INTEGER, PARAMETER :: chanrefin = 27        ! To read reference solution
  INTEGER, PARAMETER :: chanerrout = 28       ! To write difference from reference solution
  INTEGER, PARAMETER :: chanresin = 50        ! Input channel for restart run
  INTEGER, PARAMETER :: chanresout = 60       ! Restart dumps
  INTEGER, PARAMETER :: chandumpm1 = 80       ! Quick look dump file for Matlab (primal grid fields)
  INTEGER, PARAMETER :: chandumpm2 = 81       ! Quick look dump file for Matlab (dual grid fields)

END MODULE channels

! ================================================================

PROGRAM femswe

  USE runtype

  IMPLICIT NONE


  IF (IARGC() < 1) THEN
     WRITE (*,*) "Usage: ./femswe [parameter file] [suppress output]"
     CALL EXIT(-1)
  END IF

  IF (IARGC() > 1) THEN
     suppressOutput = 1
  END IF

  CALL GETARG(1, parameterFile)
  WRITE(*,*) "Using parameter file "//parameterFile

  ! ----------------------------------------------------------------

  CALL preliminary
  PRINT *,'Done preliminary'

  ! Dump the grid coordinates for use in generating
  ! reference solutions
  ! CALL dumpgrid
  ! print *,'Done dumpgrid'
  ! STOP

  ! Test the various exterior derivative and Hodge star operators
  ! CALL testop
  ! print *,'Done testop'

  ! Test the multigrid solver
  ! CALL testmg
  ! print *,'Done testmg'

  ! Generate system matrix for linearized SWEs on the f-sphere
  ! CALL fsphere
  ! print *,'Done fsphere'

  ! Generate system matrix for more general basic state on rotating sphere
  ! CALL nmodes
  ! print *,'Done nmodes'

  ! Test advection scheme
  ! CALL testadv
  ! print *,'Done testadv'

  ! Integrate in time
  CALL integrate
  PRINT *,'Done integration'

  ! ----------------------------------------------------------------

END PROGRAM femswe

! ================================================================

SUBROUTINE preliminary

  ! Preliminary calculations and setting up

  USE constants
  USE grid
  USE runtype
  USE timestep

  IMPLICIT NONE
  INTEGER :: iv0
  REAL*8 :: lat

  ! ----------------------------------------------------------------

  ! Read namelist information
  CALL readnml
  PRINT *,'Namelists read'

  ! ----------------------------------------------------------------

  ! Read in the grid data
  CALL readgrid
  PRINT *,'Done readgrid'

  ! ----------------------------------------------------------------

  ! Allocate array space now that resolution is known
  CALL allocateall
  PRINT *,'Done allocateall'

  ! ----------------------------------------------------------------

  ! Set up dual cell integrals of planetary vorticity
  PRINT *,'Should do proper area integral of f in preliminary'
  DO iv0 = 1, nvert(ngrids)
     lat = vlat(iv0,ngrids)
     planvort2(iv0) = 2.0d0*rotatn*SIN(lat)*varea(iv0,ngrids)
  ENDDO
  PRINT *,'Planetary vorticity field set up'

  ! ----------------------------------------------------------------

  ! Build a lumped version of the J, M and H matrices
  CALL buildjlump
  CALL buildmlump
  CALL buildhlump

  ! And build a local operator to approximate the inverse of M
  CALL buildxminv

  PRINT *,'Lumped J, M and H and approximate M-inverse created'

  ! ----------------------------------------------------------------

  ! Set up information needed by advection scheme
  CALL setupadv
  PRINT *,'Done setupadv'

  ! ----------------------------------------------------------------

  ! Set up orography
  CALL setorog
  PRINT *,'Done setorog'

  ! ----------------------------------------------------------------

  ! Set up initial or restart conditions
  IF (lrestart) THEN
     CALL readrestart
     PRINT *,'Done readrestart'
  ELSE
     CALL setini
     PRINT *,'Done setini'
  ENDIF

  ! ----------------------------------------------------------------

  ! Set initial time
  time = istep*dt

  ! Set beta values
  beta_v = 1.0d0 - alpha_v
  beta_pg = 1.0d0 - alpha_pg

  ! ----------------------------------------------------------------
  !  Setup table for lat lon grid conversion - P. Peixoto
  CALL setuplltable

END SUBROUTINE preliminary

! ================================================================

SUBROUTINE integrate

  ! Perform the time integration
  USE constants
  USE runtype
  USE grid
  USE state
  USE errdiag
  USE timestep

  IMPLICIT NONE
  INTEGER :: icount

  INTEGER :: clock_start, clock_end, clock_rate, itmp
  REAL*8 :: elapsed_time, u_maxerr, nonlin_alpha, phi_maxerror
  INTEGER :: nf, ne, nv, i

  REAL*8, ALLOCATABLE :: phiforce(:), uforce(:), tmp(:)
  REAL*8, ALLOCATABLE :: psi(:), utemp(:), v1(:)

  !REAL*8, ALLOCATABLE :: psigg(:), hgg(:)

  ALLOCATE(psi(nvertx), utemp(nedgex), v1(nedgex))
  ALLOCATE(phiforce(nface(ngrids)), uforce(nedge(ngrids)), tmp(nedge(ngrids)) )

  nf = nface(ngrids)
  ne = nedge(ngrids)
  nv = nvert(ngrids)

  ! ----------------------------------------------------------------

  CALL outputll
  CALL diagnostics


  CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate) ! Find the rate
  CALL SYSTEM_CLOCK(COUNT=clock_start) ! Start timing

  nonlin_alpha=1.0
  DO icount = 1, nstop
     istep = istep + 1
     CALL step
     ! PRINT *,'Calling convergence test diagnostics'
     ! CALL contest
     ! For testing impact of grid scale vorticity noise
     ! IF (istep == 720) CALL addcompmode
     time = time + dt
     IF (MODULO(istep,noutput) == 0) THEN
        !CALL output
        CALL outputll
     ENDIF
     IF (MODULO(istep,ndiagnostics) == 0) THEN
        CALL diagnostics
     ENDIF

    IF(ic==8.or.ic==9)THEN
        !Calculate error
        phi_maxerror = maxval(abs((phi2 - phi2_init)/farea(:,ngrids)))
        !print*, phi_maxerror
        !Check if error too large
        IF(phi_maxerror > 5)THEN
            PRINT*, 'Error too large, stopping here', phi_maxerror
            EXIT
        END IF

        !Hollingsworth test with linearizations
	    IF(ic==9)THEN
	 	    IF(istep==1)THEN

                !Forcing
                uforce=u1_init-u1
                phiforce=phi2_init-phi2

                !Perturb phi and u locally

                !phi2(nface(ngrids)/2)  = 1.00001*phi2(nface(ngrids)/2)
                if(nface(ngrids)>1000)then
                    i=1000
                else
                    i=18
                end if

                phi2(i)  = 1.01*phi2(nface(ngrids)/2)
                !print*, fnxtf(18,1,ngrids)
                !print*, fnxtf(18,2,ngrids)
                !print*, fnxtf(18,3,ngrids)
                phi2(fnxtf(i,1,ngrids))=1.001*phi2(fnxtf(18,1,ngrids))
                phi2(fnxtf(i,2,ngrids))=1.001*phi2(fnxtf(18,1,ngrids))
                phi2(fnxtf(i,3,ngrids))=1.001*phi2(fnxtf(18,1,ngrids))
                !u1(nedge(ngrids)/2)  = 1.00001*u1(nedge(ngrids)/2)
                !u1(nedge(ngrids)/2+1)  = 1.00001*u1(nedge(ngrids)/2+1)
                !u1(nedge(ngrids)/2-1)  = 1.00001*u1(nedge(ngrids)/2+1)
                CALL outputll
                CALL diagnostics
            END IF

            u1=u1+uforce
            phi2=phi2+phiforce

            !u_maxerr = maxval(abs((u1 - u1_init)/ldist(:,ngrids)))
            u_maxerr = SQRT((SUM((u1-u1_init)*(u1-u1_init)*ldist(:,ngrids)))/SUM(ldist(:,ngrids)))
            !h_l2 = SQRT((SUM(h*h*farea(:,ngrids)))/SUM(farea(:,ngrids)))

            if(u_maxerr>0.0000000000001) nonlin_alpha=0.00001/u_maxerr

            print '(a15, 2e16.8)', " Alpha / maxuerr:", nonlin_alpha, u_maxerr
            !PRINT *,' '
            !print '(i6,a2,1f6.2,a2,a8,1e16.8,a2,1e16.8)', istep,' (',time, 's)', &
            !    ' Alpha=',  nonlin_alpha, 'maxuerr=', u_maxerr

            u1=u1_init+nonlin_alpha*(u1-u1_init)
            phi2=phi2_init+nonlin_alpha*(phi2-phi2_init)


         END IF
     END IF

     !IF (MODULO(istep,nrestart) == 0) THEN
     !  CALL writerestart
     !ENDIF

     ! Compute and output difference from reference solution if
     ! required at this time
     ! CALL diffref
     PRINT *,'Done step ',istep,'   Time = ',time, 's (', time/24/60/60, ' dys)'
     PRINT *,' '
  ENDDO
  call outputll
  CALL SYSTEM_CLOCK(COUNT=clock_end)
  elapsed_time=REAL(clock_end-clock_start)/REAL(clock_rate)

  WRITE(*, '(''Simulation time: '',F10.2)') elapsed_time
  WRITE(*, '(''Cells per second: '',F10.2)') (REAL(nface(ngrids)*nstop) / elapsed_time)


  ! ----------------------------------------------------------------

END SUBROUTINE integrate

! ================================================================

SUBROUTINE readnml

  ! Read namelists

  USE runtype
  USE advection
  USE timestep
  USE channels
  IMPLICIT NONE


  NAMELIST /rundata/ ygridfile, ygridcoords, runid, lrestart, yresfile, ic
  NAMELIST /advectdata/ degree, monotone
  NAMELIST /timedata/ dt, niter, alpha_v, alpha_pg, nstop, noutput, nrestart, ndiagnostics

  OPEN(channml,FILE=parameterFile,DELIM='APOSTROPHE')
  READ(channml,rundata)
  READ(channml,advectdata)
  READ(channml,timedata)
  CLOSE(channml)


  ! Compute quantities derived from namelist value
  nmonomial = (degree + 1)*(degree + 2)/2
  ngauss = 1 + degree/2     ! Integer arithmetic!
  IF (ngauss > 3) THEN
     PRINT *,'*** Increase the dimension of xgauss and wgauss to ',ngauss
     PRINT *,'in MODULE advection ***'
     STOP
  ENDIF


END SUBROUTINE readnml

! ================================================================

SUBROUTINE setorog

  ! Set orography field

  USE constants
  USE grid
  USE state
  USE runtype

  IMPLICIT NONE

  INTEGER :: if0
  REAL*8 :: phis0, latc, longc, rr0, rr, rr2, long, lat

  ! ----------------------------------------------------------------

  ! Zero orography
  orog2 = 0.0d0*farea(:,ngrids)

  if(ic==5)then
  ! Williamson et al. test case 5
  phis0 = 2000.0d0*gravity
  rr0 = pi/9.0d0
  latc = pi/6.0d0
  longc = 0.5d0*pi   ! Should really be 1.5*pi but this
  ! makes plotting easier

  DO if0 = 1, nface(ngrids)
     ! Orography at dual vertices
     ! long = flong(if0,ngrids)
     ! lat = flat(if0,ngrids)
     ! Or orography at cell centroid
     CALL centroid(if0,long,lat)
     rr2 = (long - longc)**2 + (lat - latc)**2
     rr = SQRT(MIN(rr0*rr0,rr2))
     orog2(if0) = phis0*(1.0 - rr/rr0)*farea(if0,ngrids)
  ENDDO

  end if



  ! ----------------------------------------------------------------

END SUBROUTINE setorog

! ================================================================

SUBROUTINE setini

  ! Initialize prognostic fields

  USE constants
  USE runtype
  USE grid
  USE state
  USE errdiag
  USE timestep
  IMPLICIT NONE

  INTEGER :: nf, ne, nv, if0, iv0, ie0, j, jy, nygg
  REAL*8 :: u00, phi00, alpha, ca, sa, clat, slat, slon

  REAL*8, ALLOCATABLE :: psi(:), utemp(:), v1(:)

  REAL*8, ALLOCATABLE :: psigg(:), hgg(:)
  REAL*8 :: l1, l2, phiref, &
       lat0, lat1, lat2, en, umen, dygg, psi1, psi2, &
       beta, totvol, totarea, den, &
       cc1, cc2, cc3, cc4, uu1, uu2, e1, e2, hpert, &
       long, lat, x1, y1, z1, x2, y2, z2, xc, yc, zc, &
       xe, ye, ze, mag, dx, dy, dz, nx, ny, nz
  !INTEGER :: ic
  CHARACTER(len=12)::atmp !temporary char for file names , added by PXT


  ALLOCATE(psi(nvertx), utemp(nedgex), v1(nedgex))


  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(iv0) SHARED(nvertx, psi)
  DO iv0=1, nvertx
     psi(iv0) = 0
  END DO
  !$OMP END PARALLEL DO

  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(ie0) SHARED(nedgex, utemp, v1)
  DO ie0=1, nedgex
     utemp(ie0) = 0
     v1(ie0) = 0
  END DO
  !$OMP END PARALLEL DO

  ! ----------------------------------------------------------------

  nf = nface(ngrids)
  ne = nedge(ngrids)
  nv = nvert(ngrids)


  ! ic = 1   Resting, constant phi
  ! ic = 2   Williamson et al. test case 2
  ! ic = 4   Solid body rotation for normal mode calculation
  ! ic = 5   Williamson et al. test case 5
  ! ic = 7   Galewsky et al. barotropically unstable jet - specify stream fn
  ! ic = 71   Galewsky et al. barotropically unstable jet - specify U - P.Peixoto: tested and seems to have problems - avoid!
  ! ic = 8   Like case 2 but with orography balancing u
  ! ic = 9   Like case 2 but with orography balancing u - just for linear analysis - power method

  !ic = 3

  ! Construct ic element of filename
  print*,"Initial condition: ", ic
  WRITE(icname,'(i2)') ic
  icname="fem_ic"//trim(adjustl(trim(icname)))

  IF (ic ==1) THEN
     ! Resting state with uniform phi
     phi2 = 1.0d5*farea(:,ngrids)
     u1 = 0.0d0

  ELSEIF (ic == 2 .OR. ic == 5 .OR. ic == 4) THEN

     ! Balanced zonal flow at angle alpha
     IF (ic == 2) THEN
        ! Test case 2
        u00 = 2.0d0*pi*rearth/(12.0d0*86400.0d0)
        phi00 = 2.94d4
     ELSEIF (ic == 5) THEN
        ! Test case 5
        u00 = 20.0d0
        phi00 = 5960.0d0*gravity
     ELSEIF (ic == 4) THEN
        u00 = 0.0d0
        phi00 = 1.0d8
     ENDIF
     alpha = 0.0d0*pi    ! Rotation rate should be zero if alpha <> 0
     WRITE(*,*) pi
     WRITE(*,*) alpha
     WRITE(*,*) COS(alpha)
     ca = COS(alpha)
     sa = SIN(alpha)

     DO iv0 = 1, nv
        clat = COS(vlat(iv0,ngrids))
        slat = SIN(vlat(iv0,ngrids))
        slon = SIN(vlong(iv0,ngrids))
        ! Minus the stream function
        psi(iv0) = u00*rearth*(ca*slat + sa*clat*slon)
     ENDDO

     DO if0 = 1, nf
        ! Initialize at computational point...
        ! long = flong(if0,ngrids)
        ! lat = flat(if0,ngrids)
        ! ...or initialize at centroid
        CALL centroid(if0,long,lat)
        clat = COS(lat)
        slat = SIN(lat)
        slon = SIN(long)
        ! Geopotential
        phi2(if0) = phi00 - (rearth*2.0d0*rotatn*u00 + u00*u00)*0.5d0* &
             (ca*slat + sa*clat*slon)**2
        phi2(if0) = phi2(if0)*farea(if0,ngrids)
     ENDDO

     ! Non-divergent velocity field
     ! U = - D_1 (stream function); sign taken care of above
     CALL Dprimal1(psi,u1,ngrids,nv,ne)
     PRINT *,'u00 = ',u00
     PRINT *,'Max u = ',MAXVAL(u1/ldist(:,ngrids))
     PRINT *,'Should do a proper L2 projection of psi in setini'

   ! Correct depth to allow for orography, if any
     phi2 = phi2 - orog2  !For tc 5

  ELSEIF ((ic == 7) .OR. (ic == 71)) THEN

     ! Galewsky et al test case
     nygg = 2*FLOOR(SQRT(REAL(nf)))
     ALLOCATE(hgg(nygg+1), psigg(nygg+1))
     u00 = 80.0
     lat0 = pi/7.0
     lat1 = pi/2.0 - lat0
     en = exp(-4/(lat1 - lat0)**2)
     umen = u00/en
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
        IF (den .LT. 0.0D0) THEN
           uu1 = umen*exp(1.0D0/den)
        ELSE
           uu1 = 0.0D0
        ENDIF
        l2 = (j-1)*dygg - piby2
        den = (l2 - lat0)*(l2 - lat1)
        IF (den .LT. 0.0D0) THEN
           uu2 = umen*exp(1.0D0/den)
        ELSE
           uu2 = 0.0D0
        ENDIF
        psigg(j) = psigg(j-1) - 0.5d0*(uu1 + uu2)*dygg
        uu1 = uu1*(2.0d0*rotatn*SIN(l1) + TAN(l1)*uu1/rearth)
        uu2 = uu2*(2.0d0*rotatn*SIN(l2) + TAN(l2)*uu2/rearth)
        hgg(j) = hgg(j-1) - rearth*0.5d0*(uu1 + uu2)*dygg
        totarea = totarea + COS(l2)*dygg
        totvol = totvol + hgg(j)*COS(l2)*dygg
     ENDDO

     psigg(nygg+1) = psigg(nygg)
     hgg(nygg+1) = hgg(nygg)
     totvol = totvol/(totarea*gravity)
     hgg = hgg + (1.0D4 - totvol)*gravity

     ! Now assign h as a function of geographical latitude
     ! using interpolation from tabulated values
     totvol = 0.00
     totarea = 0.0D0
     ! PARALLELIZATION WARNING: THERE IS A LOOP DEPENDENCY to totvol
     DO if0 = 1, nf
        ! Initialize at computational point...
        ! long = flong(if0,ngrids)
        ! lat = flat(if0,ngrids)
        ! ...or initialize at centroid
        CALL centroid(if0,long,lat)
        l1 = lat + piby2
        jy = FLOOR(l1/dygg) + 1
        beta = (l1 - (jy - 1)*dygg)/dygg
        IF (jy == 1 .OR. jy == nygg) THEN
           ! Linear interpolation
           cc2 = 1.0D0 - beta
           cc3 = beta
           phi2(if0) = (cc2*hgg(jy) + cc3*hgg(jy+1))*farea(if0,ngrids)
        ELSE
           ! Cubic interpolation
           cc1 = -beta*(beta - 1.0D0)*(beta - 2.0D0)/6.0D0
           cc2 = 0.5D0*(beta + 1.0D0)*(beta - 1.0D0)*(beta - 2.0D0)
           cc3 = -0.5D0*(beta + 1.0D0)*beta*(beta - 2.0D0)
           cc4 = (beta + 1.0D0)*beta*(beta - 1.0D0)/6.0D0
           phi2(if0) = (cc1*hgg(jy-1) + cc2*hgg(jy) + cc3*hgg(jy+1) + cc4*hgg(jy+2))*farea(if0,ngrids)
        ENDIF
        totarea = totarea + farea(if0,ngrids)
        totvol = totvol + phi2(if0)
     ENDDO
     ! Now calculate velocity components by interpolating
     ! stream function to each vertex

     DO iv0 = 1, nv
        l1 = vlat(iv0,ngrids) + piby2
        jy = FLOOR(l1/dygg) + 1
        beta = (l1 - (jy - 1)*dygg)/dygg
        IF (jy == 1 .OR. jy == nygg) THEN
           ! Linear interpolation
           cc2 = 1.0D0 - beta
           cc3 = beta
           psi(iv0) = cc2*psigg(jy) + cc3*psigg(jy+1)
        ELSE
           ! Cubic interpolation
           cc1 = -beta*(beta - 1.0D0)*(beta - 2.0D0)/6.0D0
           cc2 = 0.5D0*(beta + 1.0D0)*(beta - 1.0D0)*(beta - 2.0D0)
           cc3 = -0.5D0*(beta + 1.0D0)*beta*(beta - 2.0D0)
           cc4 = (beta + 1.0D0)*beta*(beta - 1.0D0)/6.0D0
           psi(iv0) = cc1*psigg(jy-1) + cc2*psigg(jy) + cc3*psigg(jy+1) + cc4*psigg(jy+2)
        ENDIF
     ENDDO

     psi = psi*rearth
     DEALLOCATE(hgg, psigg)

     ! Geopotential perturbation
     alpha = 1.0D0/3.0D0
     beta = 1.0D0/15.0D0
     hpert = 120.0D0
     lat2 = 0.5D0*piby2

     DO if0 = 1, nf
        ! l2 = flat(if0,ngrids)
        ! l1 = flong(if0,ngrids)
        CALL centroid(if0,l1,l2)
        clat = COS(l2)
        IF (l1 > pi) l1 = l1 - 2.0d0*pi
        e1 = EXP(-(l1/alpha)**2)
        e2 = EXP(-((lat2 - l2)/beta)**2)
        phi2(if0) = phi2(if0) + gravity*hpert*clat*e1*e2*farea(if0,ngrids)
     ENDDO

     IF (ic == 7) THEN
        ! Non-divergent velocity field
        ! U = - D_1 (stream function)
        CALL Dprimal1(psi,utemp,ngrids,nv,ne)
        u1 = -utemp
        PRINT *,'Should do a proper L2 projection of psi in setini'
     ELSE

        DO ie0 = 1, ne
           ! Find location and latitude of primal edge midpoint and
           ! length times normal of primal edge
           iv0 = vofe(1,ie0,ngrids)
           long = vlong(if0,ngrids)
           lat = vlat(if0,ngrids)
           CALL ll2xyz(long,lat,x1,y1,z1)
           iv0 = vofe(2,ie0,ngrids)
           long = vlong(if0,ngrids)
           lat = vlat(if0,ngrids)
           CALL ll2xyz(long,lat,x2,y2,z2)
           xc = 0.5d0*(x1 + x2)
           yc = 0.5d0*(y1 + y2)
           zc = 0.5d0*(z1 + z2)
           mag = 1.0d0/(SQRT(xc*xc + yc*yc + zc*zc))
           xc = xc*mag
           yc = yc*mag
           zc = zc*mag
           dx = x2 - x1
           dy = y2 - y1
           dz = z2 - z1
           nx = -(yc*dz - zc*dy)
           ny = -(zc*dx - xc*dz)
           nz = -(xc*dy - yc*dx)
           CALL xyz2ll(xc,yc,zc,long,lat)
           ! Unit eastward vector
           xe = -yc
           ye = xc
           ze = 0.0d0
           mag = 1.0d0/(SQRT(xe*xe + ye*ye + ze*ze))
           xe = xe*mag
           ye = ye*mag
           ze = ze*mag
           den = (lat - lat0)*(lat - lat1)
           IF (den .LT. 0.0D0) THEN
              uu2 = umen*exp(1.0D0/den)
           ELSE
              uu2 = 0.0D0
           ENDIF
           u1(ie0) = uu2*rearth*(xe*nx + ye*ny + ze*nz)
        ENDDO

     ENDIF
     PRINT *,'u00 = ',u00
     PRINT *,'Max u = ',MAXVAL(u1/ldist(:,ngrids))

  ELSEIF (ic == 8 .or. ic == 9) THEN !Hollingsworth test

     u00 = 2.0d0*pi*rearth/(12.0d0*86400.0d0)
     phi00 = 2.94d4
	 phiref= 1.0d0*gravity
	 !phiref= 0.1d0*gravity

	 WRITE(atmp,'(f4.2)') phiref/gravity
     icname=trim(icname)//"_href"//trim(adjustl(trim(atmp)))

     alpha = 0.0d0*pi    ! Rotation rate should be zero if alpha <> 0
     !WRITE(*,*) pi
     !WRITE(*,*) alpha
     !WRITE(*,*) COS(alpha)
     ca = COS(alpha)
     sa = SIN(alpha)

     DO iv0 = 1, nv
        clat = COS(vlat(iv0,ngrids))
        slat = SIN(vlat(iv0,ngrids))
        slon = SIN(vlong(iv0,ngrids))
        ! Minus the stream function
        psi(iv0) = u00*rearth*(ca*slat + sa*clat*slon)
     ENDDO

     DO if0 = 1, nf
        ! Initialize at computational point...
        ! long = flong(if0,ngrids)
        ! lat = flat(if0,ngrids)
        ! ...or initialize at centroid
        CALL centroid(if0,long,lat)
        clat = COS(lat)
        slat = SIN(lat)
        slon = SIN(long)
        ! Geopotential
        !phi2(if0) = phi00 - (rearth*2.0d0*rotatn*u00 + u00*u00)*0.5d0* &
        !     (ca*slat + sa*clat*slon)**2
        !phi2(if0) = phi2(if0)*farea(if0,ngrids)
        phi2(if0) = phiref*farea(if0,ngrids)
        orog2(if0)= phi00 - (rearth*2.0d0*rotatn*u00 + u00*u00)*0.5d0* &
             (ca*slat + sa*clat*slon)**2
        orog2(if0) = orog2(if0)*farea(if0,ngrids)
     ENDDO

     ! Non-divergent velocity field
     ! U = - D_1 (stream function); sign taken care of above
     CALL Dprimal1(psi,u1,ngrids,nv,ne)
     PRINT *,'u00 = ',u00
     PRINT *,'Max u = ',MAXVAL(u1/ldist(:,ngrids))
  ENDIF


  ! Save initial state for error diagnostics
  phi2_init = phi2
  u1_init = u1

  ! Initialize dual grid diagnostic tracers
  CALL operR(phi2,xzeta2,ngrids,nf,nv)   ! Borrow xzeta2 array temporarily
  CALL HodgeJinv_TEST(xzeta2,xphibar2,ngrids,nv,-20)
  CALL massM(u1,utemp,ngrids,ne)
  CALL HodgeHinv(utemp,v1,ngrids,ne,-20)
  CALL Ddual2(v1,xzeta2,ngrids,ne,nv)
  xzeta2 = xzeta2 + planvort2


  ! Initialize step count
  istep = 0

  DEALLOCATE(psi, utemp, v1)
  ! ----------------------------------------------------------------

END SUBROUTINE setini


SUBROUTINE HodgeJinv_TEST(f1,f2,igrid,nv,niter)

  ! Apply the inverse Hodge star operator J^{-1} that maps from
  ! E_p to V_d on grid igrid.
  !
  ! If niter >= 0 then f2 is assumed to be an initial estimate
  ! for the solution, and niter further iterations are taken.
  ! If niter < 0 then f2 is initialized and -niter
  ! iterations are taken.
  ! If J is diagonal then there is no need to iterate.

  USE grid
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: niter
  INTEGER, INTENT(IN) :: igrid, nv
  REAL*8, INTENT(IN) :: f1(nv)
  REAL*8, INTENT(INOUT) :: f2(nv)
  INTEGER :: iv1, iter, miter
  REAL*8,ALLOCATABLE :: temp(:)
  REAL*8 :: relax = 1.4d0 ! relax = 1.4 is good for ijlump = 3 on hex and cube grids


  ALLOCATE(temp(nv))
  ! ----------------------------------------------------------------

  miter = ABS(niter)

  IF (niter < 0 .OR. njsmx == 1) THEN
     ! First guess based on lumped J
     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(iv1) SHARED(nv, f2, f1, jlump, igrid)
     DO iv1 = 1, nv
        f2(iv1) = f1(iv1)/jlump(iv1,igrid)
     ENDDO
     !$OMP END PARALLEL DO
  ENDIF


  IF (njsmx > 1) THEN
     ! J is not diagonal, so use Jacobi iteration to invert
     DO iter = 1, miter
        CALL HodgeJ(f2,temp,igrid,nv)
        !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(iv1) SHARED(nv, f2, temp, relax, f1, jlump, igrid)
        DO iv1 = 1, nv
           f2(iv1) = f2(iv1) + (relax*(f1(iv1) - temp(iv1)))/jlump(iv1,igrid)
        ENDDO
        !$OMP END PARALLEL DO
     ENDDO
  ENDIF

  ! ----------------------------------------------------------------


  DEALLOCATE(temp)

END SUBROUTINE HodgeJinv_TEST
! ================================================================

SUBROUTINE step

  ! Calculations for one time step

  USE constants
  USE state
  USE runtype
  USE work
  USE helmcoeff
  USE timestep
  IMPLICIT NONE

  INTEGER :: nf, ne, nv, iter, ie0, if1, if2, ninvit
  REAL*8 :: temp

  ! ----------------------------------------------------------------

  nf = nface(ngrids)
  ne = nedge(ngrids)
  nv = nvert(ngrids)

  ! Set up Coefficients for Helmholtz problem
  !$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(phiref, phi2)
  phiref = phi2
  !$OMP END PARALLEL WORKSHARE
  CALL buildhelm
  !$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(phiref, phi2, farea, ngrids)
  phiref = phi2/farea(:,ngrids)
  !$OMP END PARALLEL WORKSHARE
  CALL cell2edge(phiref,phirefe,ngrids,nf,ne)

  ! ----------------------------------------------------------------

  !print *,'*** Tidy up use of temporary arrays ***'

  ! Number of iterations for operator inverses (outside main loop)
  ! (5 to 7 is probably enough)
  ninvit = -10

  ! First guess for new time level fields
  !$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(phi2_new, phi2, u1_new, u1)
  phi2_new = phi2
  u1_new = u1
  !$OMP END PARALLEL WORKSHARE

  ! Calculate old values of dual edge circulations
  CALL massM(u1,mu,ngrids,ne) 
  CALL HodgeHinv(mu,v1,ngrids,ne,ninvit)

  ! Compute divergence - needed to modify swept areas in
  ! routine primaladvflx
  CALL Dprimal2(u1,div2,ngrids,ne,nf)
  !$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(divfac, div2, farea, ngrids, beta_v, dt)
  divfac = div2/farea(:,ngrids)
  divfac = 1.0d0/(1.0d0 + beta_v*dt*divfac)
  !$OMP END PARALLEL WORKSHARE

  ! Old values of cell integrals of KE
  CALL operT(u1,b2,ngrids,ne,nf)

  ! Add to geopotential
  !$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(b2, phi2, beta_pg, orog2)
  ! Temporal weighting, and include orography
  b2 = (phi2 + 0.5d0*b2)*beta_pg + orog2
  !$OMP END PARALLEL WORKSHARE


  ! Construct old dual cell integrals of vorticity, and hence PV
  ! First relative vorticity
  CALL Ddual2(v1,zeta2,ngrids,ne,nv)
  ! Add planetary contribution

  !$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(zeta2, planvort2)
  zeta2 = zeta2 + planvort2
  !$OMP END PARALLEL WORKSHARE

  ! Dual cell geopotential
  CALL operR(phi2,pv,ngrids,nf,nv)  ! Borrow pv array temporarily
  CALL HodgeJinv(pv,phibar2,ngrids,nv,ninvit)

  !$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pva, zeta2, phibar2, varea, ngrids,pv)
  ! Vertex PV
  pv = zeta2/phibar2
  ! Finally the dual cell area integrals of PV
  pva = pv*varea(:,ngrids)
  !$OMP END PARALLEL WORKSHARE

  ! -----------------------------------------------------------------

  ! Begin iterative solution of nonlinear problem
  DO iter = 1, niter

     ! Number of iterations for operator inverses
     IF (iter == 1) THEN
        ninvit = -2
     ELSE
        ninvit = 2
     ENDIF

     ! Calculate new values of M times U
     CALL massM(u1_new,mu,ngrids,ne)

     ! New values of cell integrals of KE
     CALL operT(u1_new,b2_temp,ngrids,ne,nf)

     !$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(b2_temp, phi2_new, b2, alpha_pg, dt)
     ! Add to geopotential
     b2_temp = (b2 + alpha_pg*(phi2_new + 0.5d0*b2_temp))*dt
     !$OMP END PARALLEL WORKSHARE

     ! And compute M times gradient
     CALL massL(b2_temp,b0,ngrids,nf)
     CALL Ddual1(b0,gb1,ngrids,nf,ne)

     ! Construct time averaged velocities
     !$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(ubar1, beta_v, u1, alpha_v, u1_new)
     ubar1 = beta_v*u1 + alpha_v*u1_new
     !$OMP END PARALLEL WORKSHARE

     CALL massM(ubar1,mu,ngrids,ne)
     ! Corresponding v field
     CALL HodgeHinv(mu,vbar1,ngrids,ne,ninvit)
     ! Perpendicular components
     CALL operW(ubar1,wu,ngrids,ne)
     CALL massMinv(wu,uperpbar1,ngrids,ne,ninvit)
     CALL HodgeHinv(wu,vperpbar1,ngrids,ne,ninvit)

     ! Compute primal grid mass fluxes
     CALL primaladvflx(phi2,mf1)

     ! Construct dual grid mass fluxes
     CALL operW(mf1,wf,ngrids,ne)
     CALL HodgeHinv(wf,mfperp1,ngrids,ne,ninvit)

     ! Compute dual grid PV fluxes
     CALL dualadvflx(pva,qperp1)
     CALL HodgeH(qperp1,hqperp1,ngrids,ne)

     ! Divergence of mass flux
     CALL Dprimal2(mf1,divmf,ngrids,ne,nf)

     ! Residual in mass equation
     !$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(rphi2, phi2_new, phi2, divmf, u_inc, u1_new, u1)
     rphi2 = phi2_new - phi2 + divmf

     ! Residual in momentum equation
     u_inc = u1_new - u1    ! Borrow array u_inc for temporary use
     !$OMP END PARALLEL WORKSHARE

     CALL massM(u_inc,ru1,ngrids,ne)
     !$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(ru1, hqperp1, gb1)
     ru1 = ru1 - hqperp1 + gb1
     !$OMP END PARALLEL WORKSHARE

     IF (suppressOutput == 0) THEN
        PRINT *,'RMS rphi = ',SQRT(SUM((rphi2/farea(:,ngrids))**2)/nf) &
             ,'    MAX rphi = ',MAXVAL(ABS(rphi2/farea(:,ngrids)))
        PRINT *,'RMS ru   = ',SQRT(SUM((ru1/ldist(:,ngrids))**2)/ne) &
             ,'    MAX ru   = ',MAXVAL(ABS(ru1/ldist(:,ngrids)))
     END IF
     ! Build RHS of Helmholtz problem
     CALL approxMinv(ru1,mru,ngrids,ne)    ! Use approximation to M-inverse
     !$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(mru, phirefe)
     mru = phirefe*mru
     !$OMP END PARALLEL WORKSHARE
     CALL Dprimal2(mru,rhs,ngrids,ne,nf)
     !$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(rphi2, rhs, alpha_v, dt)
     rhs =  rphi2 - rhs*alpha_v*dt
     !$OMP END PARALLEL WORKSHARE

     ! Solve Helmholtz problem
     CALL mgsolve(phi_inc,rhs,ngrids)

     ! Backsubstitute for u_inc
     CALL massL(phi_inc,b0,ngrids,nf)
     CALL Ddual1(b0,gb1,ngrids,nf,ne)
     !$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(mru, ru1, gb1, alpha_pg, dt)
     mru = -ru1 - gb1*alpha_pg*dt
     !$OMP END PARALLEL WORKSHARE
     CALL approxMinv(mru,u_inc,ngrids,ne)

     ! Update prognostic variables
     !$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(phi2_new, phi_inc, u1_new, u_inc)
     phi2_new = phi2_new + phi_inc
     u1_new = u1_new + u_inc
     !$OMP END PARALLEL WORKSHARE

  ENDDO

  ! Update diagnostic tracers using last estimate of
  ! dual mass flux
  CALL Ddual2(mfperp1,pva,ngrids,ne,nv)  ! Borrow pva array temporarily

  !$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pv, xzeta2, xphibar2, phibar2, pva, ngrids, varea)
  xphibar2 = xphibar2 + pva
  pv = xzeta2/phibar2
  pva = pv*varea(:,ngrids)
  !$OMP END PARALLEL WORKSHARE

  CALL dualadvflx(pva,qperp1)
  CALL Ddual2(qperp1,pva,ngrids,ne,nv)

  !$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(xzeta2, pva, phi2, phi2_new, u1, u1_new)
  xzeta2 = xzeta2 + pva
  ! Rename prognostic variables ready for next step
  phi2 = phi2_new
  u1 = u1_new
  !$OMP END PARALLEL WORKSHARE


  ! ----------------------------------------------------------------

END SUBROUTINE step

! ================================================================

SUBROUTINE contest

  ! Calculations for one time step
  ! with convergence diagnostics

  USE constants
  USE state
  USE work
  USE helmcoeff
  USE timestep
  IMPLICIT NONE

  INTEGER :: nf, ne, nv, iter, ie0, if1, if2, ninvit
  REAL*8 :: temp, maxrphi, maxru, maxephi, maxeu
  REAL*8, ALLOCATABLE :: phitru(:), utru(:)

  ALLOCATE(phitru(nfacex), utru(nedgex))
  ! ----------------------------------------------------------------

  nf = nface(ngrids)
  ne = nedge(ngrids)
  nv = nvert(ngrids)

  ! Set up Coefficients for Helmholtz problem
  phiref = phi2
  CALL buildhelm
  phiref = phi2/farea(:,ngrids)
  CALL cell2edge(phiref,phirefe,ngrids,nf,ne)

  ! ----------------------------------------------------------------

  !print *,'*** Tidy up use of temporary arrays ***'

  ! First guess for new time level fields
  phi2_new = phi2
  u1_new = u1

  ! Calculate old values of dual edge circulations
  CALL massM(u1,mu,ngrids,ne) 
  CALL HodgeHinv(mu,v1,ngrids,ne,-20)

  ! Compute divergence - needed to modify swept areas in
  ! routine primaladvflx
  CALL Dprimal2(u1,div2,ngrids,ne,nf)
  divfac = div2/farea(:,ngrids)
  divfac = 1.0d0/(1.0d0 + beta_v*dt*divfac)

  ! Old values of cell integrals of KE
  CALL operT(u1,b2,ngrids,ne,nf)

  ! Add to geopotential
  b2 = phi2 + 0.5d0*b2

  ! Temporal weighting, and include orography
  b2 = b2*beta_pg + orog2

  ! Construct old dual cell integrals of vorticity, and hence PV
  ! First relative vorticity
  CALL Ddual2(v1,zeta2,ngrids,ne,nv)
  ! Add planetary contribution
  zeta2 = zeta2 + planvort2
  ! Dual cell geopotential
  CALL operR(phi2,pv,ngrids,nf,nv)  ! Borrow pv array temporarily
  CALL HodgeJinv(pv,phibar2,ngrids,nv,-20)

  ! Vertex PV
  pv = zeta2/phibar2
  ! Finally the dual cell area integrals of PV
  pva = pv*varea(:,ngrids)

  ! -----------------------------------------------------------------

  ! Begin iterative solution of nonlinear problem
  DO iter = 1, 20

     PRINT *,'Iter = ',iter
     ! Number of iterations for operator inverses
     IF (iter == 1) THEN
        ninvit = -20
     ELSE
        ninvit = 20
     ENDIF

     ! Calculate new values of M times U
     CALL massM(u1_new,mu,ngrids,ne)

     ! New values of cell integrals of KE
     CALL operT(u1_new,b2_temp,ngrids,ne,nf)

     ! Add to geopotential
     b2_temp = phi2_new + 0.5d0*b2_temp

     ! Time averaged value, times dt
     b2_temp = (b2 + alpha_pg*b2_temp)*dt

     ! And compute M times gradient
     CALL massL(b2_temp,b0,ngrids,nf)
     CALL Ddual1(b0,gb1,ngrids,nf,ne)

     ! Construct time averaged velocities
     ubar1 = beta_v*u1 + alpha_v*u1_new
     CALL massM(ubar1,mu,ngrids,ne)
     ! Corresponding v field
     CALL HodgeHinv(mu,vbar1,ngrids,ne,ninvit)
     ! Perpendicular components
     CALL operW(ubar1,wu,ngrids,ne)
     CALL massMinv(wu,uperpbar1,ngrids,ne,ninvit)
     CALL HodgeHinv(wu,vperpbar1,ngrids,ne,ninvit)

     ! Compute primal grid mass fluxes
     CALL primaladvflx(phi2,mf1)

     ! Construct dual grid mass fluxes
     CALL operW(mf1,wf,ngrids,ne)
     CALL HodgeHinv(wf,mfperp1,ngrids,ne,ninvit)

     ! Compute dual grid PV fluxes
     CALL dualadvflx(pva,qperp1)
     CALL HodgeH(qperp1,hqperp1,ngrids,ne)

     ! Divergence of mass flux
     CALL Dprimal2(mf1,divmf,ngrids,ne,nf)

     ! Residual in mass equation
     rphi2 = phi2_new - phi2 + divmf

     ! Residual in momentum equation
     u_inc = u1_new - u1    ! Borrow array u_inc for temporary use
     CALL massM(u_inc,ru1,ngrids,ne)
     ru1 = ru1 - hqperp1 + gb1

     ! Build RHS of Helmholtz problem
     CALL approxMinv(ru1,mru,ngrids,ne)    ! Use approximation to M-inverse
     mru = phirefe*mru
     CALL Dprimal2(mru,rhs,ngrids,ne,nf)
     rhs =  rphi2 - rhs*alpha_v*dt

     ! Solve Helmholtz problem
     CALL mgsolve(phi_inc,rhs,ngrids)

     ! Backsubstitute for u_inc
     CALL massL(phi_inc,b0,ngrids,nf)
     CALL Ddual1(b0,gb1,ngrids,nf,ne)
     mru = -ru1 - gb1*alpha_pg*dt
     CALL approxMinv(mru,u_inc,ngrids,ne)

     ! Update prognostic variables
     phi2_new = phi2_new + phi_inc
     u1_new = u1_new + u_inc

  ENDDO

  ! Save converged solution
  phitru = phi2_new
  utru = u1_new


  ! Reset first guess for new time level fields
  phi2_new = phi2
  u1_new = u1

  ! Diagnose max residuals and errors
  maxephi = MAXVAL(ABS((phi2_new - phitru)/farea(:,ngrids)))
  maxeu = MAXVAL(ABS((u1_new - utru)/ldist(:,ngrids)))
  PRINT *,'    iter = ',0
  PRINT *,'    MAX ephi = ',maxephi,'    MAX eu = ',maxeu
  WRITE(66,'(4e12.3)') maxephi, maxeu

  ! -----------------------------------------------------------------

  ! Redo calculations on step n fields

  ! Calculate old values of dual edge circulations
  CALL massM(u1,mu,ngrids,ne) 
  CALL HodgeHinv(mu,v1,ngrids,ne,-20)

  ! Compute divergence - needed to modify swept areas in
  ! routine primaladvflx
  CALL Dprimal2(u1,div2,ngrids,ne,nf)
  divfac = div2/farea(:,ngrids)
  divfac = 1.0d0/(1.0d0 + beta_v*dt*divfac)

  ! Old values of cell integrals of KE
  CALL operT(u1,b2,ngrids,ne,nf)

  ! Add to geopotential
  b2 = phi2 + 0.5d0*b2

  ! Temporal weighting, and include orography
  b2 = b2*beta_pg + orog2

  ! Construct old dual cell integrals of vorticity, and hence PV
  ! First relative vorticity
  CALL Ddual2(v1,zeta2,ngrids,ne,nv)
  ! Add planetary contribution
  zeta2 = zeta2 + planvort2
  ! Dual cell geopotential
  CALL operR(phi2,pv,ngrids,nf,nv)  ! Borrow pv array temporarily
  CALL HodgeJinv(pv,phibar2,ngrids,nv,-20)

  ! Vertex PV
  pv = zeta2/phibar2
  ! Finally the dual cell area integrals of PV
  pva = pv*varea(:,ngrids)


  ! -----------------------------------------------------------------

  ! Now go again with diagnostics
  ! Begin iterative solution of nonlinear problem
  DO iter = 1, niter

     ! Number of iterations for operator inverses
     IF (iter == 1) THEN
        ninvit = -20
     ELSE
        ninvit = 20
     ENDIF

     ! Calculate new values of M times U
     CALL massM(u1_new,mu,ngrids,ne)

     ! New values of cell integrals of KE
     CALL operT(u1_new,b2_temp,ngrids,ne,nf)

     ! Add to geopotential
     b2_temp = phi2_new + 0.5d0*b2_temp

     ! Time averaged value, times dt
     b2_temp = (b2 + alpha_pg*b2_temp)*dt

     ! And compute M times gradient
     CALL massL(b2_temp,b0,ngrids,nf)
     CALL Ddual1(b0,gb1,ngrids,nf,ne)

     ! Construct time averaged velocities
     ubar1 = beta_v*u1 + alpha_v*u1_new
     CALL massM(ubar1,mu,ngrids,ne)
     ! Corresponding v field
     CALL HodgeHinv(mu,vbar1,ngrids,ne,ninvit)
     ! Perpendicular components
     CALL operW(ubar1,wu,ngrids,ne)
     CALL massMinv(wu,uperpbar1,ngrids,ne,ninvit)
     CALL HodgeHinv(wu,vperpbar1,ngrids,ne,ninvit)

     ! Compute primal grid mass fluxes
     CALL primaladvflx(phi2,mf1)

     ! Construct dual grid mass fluxes
     CALL operW(mf1,wf,ngrids,ne)
     CALL HodgeHinv(wf,mfperp1,ngrids,ne,ninvit)

     ! Compute dual grid PV fluxes
     CALL dualadvflx(pva,qperp1)
     CALL HodgeH(qperp1,hqperp1,ngrids,ne)

     ! Divergence of mass flux
     CALL Dprimal2(mf1,divmf,ngrids,ne,nf)

     ! Residual in mass equation
     rphi2 = phi2_new - phi2 + divmf

     ! Residual in momentum equation
     u_inc = u1_new - u1    ! Borrow array u_inc for temporary use
     CALL massM(u_inc,ru1,ngrids,ne)
     ru1 = ru1 - hqperp1 + gb1

     ! Build RHS of Helmholtz problem
     CALL approxMinv(ru1,mru,ngrids,ne)    ! Use approximation to M-inverse
     mru = phirefe*mru
     CALL Dprimal2(mru,rhs,ngrids,ne,nf)
     rhs =  rphi2 - rhs*alpha_v*dt

     ! Solve Helmholtz problem
     CALL mgsolve(phi_inc,rhs,ngrids)

     ! Backsubstitute for u_inc
     CALL massL(phi_inc,b0,ngrids,nf)
     CALL Ddual1(b0,gb1,ngrids,nf,ne)
     mru = -ru1 - gb1*alpha_pg*dt
     CALL approxMinv(mru,u_inc,ngrids,ne)

     ! Update prognostic variables
     phi2_new = phi2_new + phi_inc
     u1_new = u1_new + u_inc

     ! Diagnose max residuals and errors
     maxrphi = MAXVAL(ABS(rphi2/farea(:,ngrids)))
     maxru = MAXVAL(ABS(ru1/ldist(:,ngrids)))
     maxephi = MAXVAL(ABS((phi2_new - phitru)/farea(:,ngrids)))
     maxeu = MAXVAL(ABS((u1_new - utru)/ldist(:,ngrids)))
     PRINT *,'    iter = ',iter
     PRINT *,'    MAX rphi = ',maxrphi,'    MAX ru = ',maxru
     PRINT *,'    MAX ephi = ',maxephi,'    MAX eu = ',maxeu
     WRITE(66,'(4e12.3)') maxephi, maxeu

  ENDDO



  ! Update diagnostic tracers using last estimate of
  ! dual mass flux
  CALL Ddual2(mfperp1,pva,ngrids,ne,nv)  ! Borrow pva array temporarily
  xphibar2 = xphibar2 + pva
  pv = xzeta2/phibar2
  pva = pv*varea(:,ngrids)
  CALL dualadvflx(pva,qperp1)
  CALL Ddual2(qperp1,pva,ngrids,ne,nv)
  xzeta2 = xzeta2 + pva

  ! Rename prognostic variables ready for next step
  phi2 = phi2_new
  u1 = u1_new

  DEALLOCATE(phitru, utru)

  ! ----------------------------------------------------------------

END SUBROUTINE contest

! ================================================================

SUBROUTINE setupadv

  ! Set up information needed for advection scheme

  USE advection
  IMPLICIT NONE

  ! ----------------------------------------------------------------

  ! Build stencils for advection
  CALL buildadvsten
  PRINT *,'  Done buildadvsten'

  ! Build integrals of monomials over advection stencils
  CALL buildintmon
  PRINT *,'  Done buildintmon'

  ! Build coefficients used to construct polynomial fits
  ! for advection scheme
  CALL buildadvcoeff
  PRINT *,'  Done buildadvcoeff'

  ! Set up Gauss points for swept area integrals
  CALL setupgauss
  PRINT *,'  Done setupgauss'

  ! ----------------------------------------------------------------

END SUBROUTINE setupadv

! ================================================================

SUBROUTINE allocateall

  ! Allocate array space that will be needed throughout the code

  USE constants
  USE state
  USE work
  USE errdiag
  USE helmcoeff
  USE advection

  IMPLICIT NONE

  INTEGER :: igrid
  INTEGER :: if0, iv0, ie0

  ! ----------------------------------------------------------------

  ! Arrays in module constants
  ALLOCATE(planvort2(nvert(ngrids)))

  ! Arrays in module state
  ALLOCATE(orog2(nface(ngrids)))
  ALLOCATE(phi2(nface(ngrids)), u1(nedge(ngrids)), c2(nface(ngrids)))
  ALLOCATE(xphibar2(nvert(ngrids)), xzeta2(nvert(ngrids)))


  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(if0) SHARED(orog2, phi2, c2, nface, ngrids)
  DO if0=1, nface(ngrids)
     orog2(if0) = 0
     phi2(if0) = 0
     c2(if0) = 0
  END DO
  !$OMP END PARALLEL DO
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(ie0) SHARED(u1, nedge, ngrids)
  DO ie0=1, nedge(ngrids)
     u1(ie0) = 0
  END DO
  !$OMP END PARALLEL DO
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(iv0) SHARED(planvort2, xphibar2, xzeta2, nvert, ngrids)
  DO iv0=1, nvert(ngrids)
     planvort2(iv0) = 0
     xphibar2(iv0) = 0
     xzeta2(iv0) = 0
  END DO
  !$OMP END PARALLEL DO

  ! Arrays in module work
  ALLOCATE(v1(nedge(ngrids)), vbar1(nedge(ngrids)), ubar1(nedge(ngrids)), &
       vperpbar1(nedge(ngrids)), uperpbar1(nedge(ngrids)), &
       mf1(nedge(ngrids)), mfperp1(nedge(ngrids)), &
       divfac(nface(ngrids)))

  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(if0) SHARED(mfperp1, nface, ngrids)
  DO if0=1, nface(ngrids)
     mfperp1(if0) = 0
  END DO
  !$OMP END PARALLEL DO
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(ie0) &
  !$OMP    SHARED(v1, vbar1, ubar1, vperpbar1, uperpbar1, mf1, mfperp1, nedge, ngrids)
  DO ie0=1, nedge(ngrids)
     v1(ie0) = 0
     vbar1(ie0) = 0
     ubar1(ie0) = 0
     vperpbar1(ie0) = 0
     uperpbar1(ie0) = 0
     mf1(ie0) = 0
     mfperp1(ie0) = 0
  END DO
  !$OMP END PARALLEL DO


  ALLOCATE(b2(nfacex), pva(nvertx), u1_temp(nedgex), &
       b0(nfacex), gb1(nedgex), div2(nfacex), zeta2(nvertx), &
       pv(nvertx), qperp1(nedgex), hqperp1(nedgex), divmf(nfacex), phirefe(nedgex), &
       b2_temp(nfacex), phi_inc(nfacex), u_inc(nedgex), &
       rphi2(nfacex), ru1(nedgex), rhs(nfacex), &
       phi2_new(nfacex), u1_new(nedgex), mru(nedgex), &
       phibar2(nvertx), mu(nedgex), wu(nedgex), wf(nedgex))

  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(if0) &
  !$OMP    SHARED(b2, b0, div2, divmf, b2_temp, phi_inc, rphi2, rhs, phi2_new, nface, ngrids)
  DO if0=1, nface(ngrids)
     b2(if0) = 0
     b0(if0) = 0
     div2(if0) = 0
     divmf(if0) = 0
     b2_temp(if0) = 0
     phi_inc(if0) = 0
     rphi2(if0) = 0
     rhs(if0) = 0
     phi2_new(if0) = 0
  END DO
  !$OMP END PARALLEL DO
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(ie0) SHARED(u1_temp, gb1, qperp1, hqperp1, phirefe,&
  !$OMP         u_inc, u1_new, mru, mu, wu, wf, nedge, ngrids)
  DO ie0=1, nedge(ngrids)
     u1_temp(ie0) = 0
     gb1(ie0) = 0
     qperp1(ie0) = 0
     hqperp1(ie0) = 0
     phirefe(ie0) = 0
     u_inc(ie0) = 0
     u1_new(ie0) = 0
     mru(ie0) = 0
     mu(ie0) = 0
     wu(ie0) = 0
     wf(ie0) = 0
  END DO
  !$OMP END PARALLEL DO
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(iv0) SHARED(pva, zeta2, pv, phibar2, nvert, ngrids)
  DO iv0=1, nvert(ngrids)
     pva(iv0) = 0
     zeta2(iv0) = 0
     pv(iv0) = 0
     phibar2(iv0) = 0
  END DO
  !$OMP END PARALLEL DO


  ! Arrays in module errdiag
  ALLOCATE(phi2_init(nface(ngrids)), u1_init(nedge(ngrids)))
  ALLOCATE(phiexac(nface(ngrids)), pvexac(nvert(ngrids)))
  ALLOCATE(phierr(nface(ngrids)), pverr(nvert(ngrids)), u1_err(nedge(ngrids)))



  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(if0) SHARED(phi2_init, phiexac, phierr, nface, ngrids)
  DO if0=1, nface(ngrids)
     phi2_init(if0) = 0
     phiexac(if0) = 0
     phierr(if0) = 0
  END DO
  !$OMP END PARALLEL DO
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(ie0) SHARED(u1_init, u1_err, nedge, ngrids)
  DO ie0=1, nedge(ngrids)
     u1_init(ie0) = 0
     u1_err(ie0) = 0
  END DO
  !$OMP END PARALLEL DO
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(iv0) SHARED(pvexac, pverr, nvert, ngrids)
  DO iv0=1, nvert(ngrids)
     pvexac(iv0) = 0
     pverr(iv0) = 0
  END DO
  !$OMP END PARALLEL DO


  ! Arrays in module helmholtz
  ALLOCATE(phiref(nface(ngrids)))
  ALLOCATE(nusq(nedgex,ngrids))
  ALLOCATE(helmdiag(nfacex,ngrids))
  ALLOCATE(underrel(ngrids))


  DO igrid=1, ngrids
     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(if0) SHARED(phiref, helmdiag, nface, igrid, ngrids)
     DO if0=1, nface(ngrids)
        phiref(if0) = 0
        helmdiag(if0,igrid) = 0
     END DO
     !$OMP END PARALLEL DO
     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(ie0) SHARED(nusq, nedge, igrid, ngrids)
     DO ie0=1, nedge(ngrids)
        nusq(ie0,igrid) = 0
     END DO
     !$OMP END PARALLEL DO


     underrel(igrid) = 0
  END DO

  ! Arrays in module advection
  ! We don't know the maximum stencil size in general until
  ! we have built the stencils. To avoid convoluted code,
  ! set maximum stencil size here for known grids and degrees
  ! of fit and check later once stencils are built.
  IF (nefmx == 3) THEN
     ! Triangular primal grid
     IF (degree == 0) THEN
        ! At most one cell for piecewise constant
        nadvfmx = 1
     ELSEIF (degree == 1) THEN
        ! At most 4 cells for piecewise linear
        nadvfmx = 4
     ELSEIF (degree == 2) THEN
        ! At most 10 cells for piecewise quadratic
        nadvfmx = 10
     ELSEIF (degree == 3) THEN
        ! At most 10 cells for piecewise cubic
        nadvfmx = 10
     ELSEIF (degree == 4) THEN
        ! At most 22 cells for piecewise quartic
        nadvfmx = 22
     ELSE
        PRINT *,'Configuration not known.'
        PRINT *,'Unable to set nadvfmx in subroutine allocateall.'
        STOP
     ENDIF
  ELSEIF (nefmx == 4) THEN
     ! Quadrilateral primal grid
     IF (degree == 0) THEN
        ! At most one cell for piecewise constant
        nadvfmx = 1
     ELSEIF (degree == 1) THEN
        ! At most 5 cells for piecewise linear
        nadvfmx = 5
     ELSEIF (degree == 2) THEN
        ! At most 9 cells for piecewise quadratic
        nadvfmx = 9
     ELSEIF (degree == 3) THEN
        ! At most 21 cells for piecewise cubic
        nadvfmx = 21
     ELSEIF (degree == 4) THEN
        ! At most 23 cells for piecewise quartic
        nadvfmx = 23
     ELSE
        PRINT *,'Configuration not known.'
        PRINT *,'Unable to set nadvfmx in subroutine allocateall.'
        STOP
     ENDIF
  ELSEIF (nefmx == 6) THEN
     ! Hexagonal primal grid
     IF (degree == 0) THEN
        ! At most one cell for piecewise constant
        nadvfmx = 1
     ELSEIF (degree == 1) THEN
        ! At most 7 cells for piecewise linear
        nadvfmx = 7
     ELSEIF (degree == 2) THEN
        ! At most 7 cells for piecewise quadratic
        nadvfmx = 7
     ELSEIF (degree == 3) THEN
        ! At most 13 cells for piecewise cubic
        nadvfmx = 13
     ELSEIF (degree == 4) THEN
        ! At most 19 cells for piecewise quartic
        nadvfmx = 19
     ELSE
        PRINT *,'Configuration not known.'
        PRINT *,'Unable to set nadvfmx in subroutine allocateall.'
        STOP
     ENDIF
  ELSE
     PRINT *,'Configuration not known.'
     PRINT *,'Unable to set nadvfmx in subroutine allocateall.'
     STOP
  ENDIF

  IF (nevmx == 3) THEN
     ! Triangular dual grid
     IF (degree == 0) THEN
        ! At most one cell for piecewise constant
        nadvvmx = 1
     ELSEIF (degree == 1) THEN
        ! At most 4 cells for piecewise linear
        nadvvmx = 4
     ELSEIF (degree == 2) THEN
        ! At most 10 cells for piecewise quadratic
        nadvvmx = 10
     ELSEIF (degree == 3) THEN
        ! At most 10 cells for piecewise cubic
        nadvvmx = 10
     ELSEIF (degree == 4) THEN
        ! At most 22 cells for piecewise quartic
        nadvvmx = 22
     ELSE
        PRINT *,'Configuration not known.'
        PRINT *,'Unable to set nadvvmx in subroutine allocateall.'
        STOP
     ENDIF
  ELSEIF (nevmx == 4) THEN
     ! Quadrilateral dual grid
     IF (degree == 0) THEN
        ! At most one cell for piecewise constant
        nadvvmx = 1
     ELSEIF (degree == 1) THEN
        ! At most 5 cells for piecewise linear
        nadvvmx = 5
     ELSEIF (degree == 2) THEN
        ! At most 9 cells for piecewise quadratic
        nadvvmx = 9
     ELSEIF (degree == 3) THEN
        ! At most 21 cells for piecewise cubic
        nadvvmx = 21
     ELSEIF (degree == 4) THEN
        ! At most 21 cells for piecewise quartic
        nadvvmx = 21
     ELSE
        PRINT *,'Configuration not known.'
        PRINT *,'Unable to set nadvvmx in subroutine allocateall.'
        STOP
     ENDIF
  ELSEIF (nevmx == 6) THEN
     ! Hexagonal grid
     IF (degree == 0) THEN
        ! At most one cell for piecewise constant
        nadvvmx = 1
     ELSEIF (degree == 1) THEN
        ! At most 7 cells for piecewise linear
        nadvvmx = 7
     ELSEIF (degree == 2) THEN
        ! At most 7 cells for piecewise quadratic
        nadvvmx = 7
     ELSEIF (degree == 3) THEN
        ! At most 13 cells for piecewise cubic
        nadvvmx = 13
     ELSEIF (degree == 4) THEN
        ! At most 19 cells for piecewise quartic
        nadvvmx = 19
     ELSE
        PRINT *,'Configuration not known.'
        PRINT *,'Unable to set nadvvmx in subroutine allocateall.'
        STOP
     ENDIF
  ELSE
     PRINT *,'Configuration not known.'
     PRINT *,'Unable to set nadvvmx in subroutine allocateall.'
     STOP
  ENDIF

  ALLOCATE(nstenadvf(nfacex), nexadvf(nfacex), &
       stenadvf(nadvfmx,nfacex), intmonf(nadvfmx,nfacex,nmonomial))
  ALLOCATE(nstenadvv(nvertx), nexadvv(nvertx), &
       stenadvv(nadvvmx,nvertx), intmonv(nadvvmx,nvertx,nmonomial))


  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(if0) SHARED(nstenadvf, nexadvf, stenadvf, intmonf, nface, ngrids)
  DO if0=1, nface(ngrids)
     nstenadvf(if0) = 0
     nexadvf(if0) = 0
     stenadvf(:,if0) = 0
     intmonf(:,if0,:) = 0
  END DO
  !$OMP END PARALLEL DO
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(iv0) SHARED(nstenadvv, nexadvv, stenadvv, intmonv, nvert, ngrids)
  DO iv0=1, nvert(ngrids)
     nstenadvv(iv0) = 0
     nexadvv(iv0) = 0
     stenadvv(:,iv0) = 0
     intmonv(:,iv0,:) = 0
  END DO
  !$OMP END PARALLEL DO



  ! ----------------------------------------------------------------

END SUBROUTINE allocateall

! ================================================================

SUBROUTINE buildhelm

  ! Build the cell edge values related to phiref needed for
  ! the multigrid Helmholtz solver at all resolutions

  USE grid
  USE helmcoeff
  USE constants
  USE timestep

  IMPLICIT NONE
  INTEGER :: nf2, ne2, nf1, igrid, if1, if2, if3, ie1, ie2, &
       ixd2, ixm, ixd1, ixl
  REAL*8 :: const, temp, &
       cd2, cm, cd1, cl

  REAL*8, ALLOCATABLE :: temp1(:), temp2(:)
  ! real*8 :: phi00, mu2

  ALLOCATE(temp1(nfacex), temp2(nfacex))

  ! ---------------------------------------------------------------

  ! Useful constant
  const = alpha_pg*alpha_v*dt*dt

  ! Under-relaxation parameter. Should be the best compromise
  ! between damping small-scale error when Laplacian dominates
  ! and damping large-scale error when undifferentiated term
  ! dominates. Depends on grid structure and grid resolution.
  ! *** There is scope to optimize here ***
  DO igrid = 1, ngrids
     IF (nefmx == 6) THEN
        ! Hexagonal grid
        ! mu2 = const*phi00/fareamin(igrid)
        ! underrel(igrid) = 1.0d0
        underrel(igrid) = 0.8d0
        ! underrel(igrid) = (4.0d0*mu2 + 1.0d0)/(6.0d0*mu2 + 1.0d0)
     ELSEIF (nefmx == 4) THEN
        ! Cubic grid
        ! mu2 = const*phi00/fareamin(igrid)
        ! underrel(igrid) = 1.0d0
        underrel(igrid) = 0.8d0
        ! underrel(igrid) = (4.0d0*mu2 + 1.0d0)/(6.0d0*mu2 + 1.0d0)
     ELSE
        PRINT *,'Choose a sensible value for underrel in buildhelm'
        STOP
     ENDIF
  ENDDO


  ! Construct cell edge values of wave Courant number

  ! Finest grid
  nf2 = nface(ngrids)
  ne2 = nedge(ngrids)


  !$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(temp2, const, phiref, nf2, temp1, ngrids, farea)
  temp2 = const*phiref
  temp1(1:nf2) = temp2(1:nf2)/farea(1:nf2,ngrids)
  !$OMP END PARALLEL WORKSHARE

  CALL cell2edge(temp1,nusq(1,ngrids),ngrids,nf2,ne2)

  ! Coarser grids
  DO igrid = ngrids-1, 1, -1
     nf1 = nface(igrid+1)
     !$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(temp1, nf2, nf1, temp2)
     temp1(1:nf1) = temp2(1:nf1)
     !$OMP END PARALLEL WORKSHARE
     nf2 = nface(igrid)
     ne2 = nedge(igrid)
     CALL restrict(temp1,nf1,temp2,nf2,igrid)
!!!!!!!!!!!!!!!!!!!!!
     !CHANGED DUE TO BUG: ngrids to igrid in the next line
!!!!!!!!!!!!!!!!!!!!!
     !$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(temp1, nf2, temp2, farea, igrid)
     temp1(1:nf2) = temp2(1:nf2)/farea(1:nf2,igrid)
     !$OMP END PARALLEL WORKSHARE
     CALL cell2edge(temp1,nusq(1,igrid),igrid,nf2,ne2)
  ENDDO


  ! Extract diagonal coefficient of Helmholtz operator
  DO igrid = 1, ngrids
     ! Loop over cells
     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(temp, cl, ixm, ixd1, ixl, ixd2, ie1, cd1, cd2, ie2, if1, if2, if3, cm) &
     !$OMP SHARED(neoff, eoff, eoffin, nxminvsten, xminvsten, xminv, nusq, fnxte, nlsten, lsten, lmass, helmdiag, igrid, nface)
     DO if1 = 1, nface(igrid)

        temp = -1.0d0
        ! Loop over edges of if1 involved in Dprimal2 operator
        DO ixd2 = 1, neoff(if1,igrid)
           ie1 = eoff(ixd2,if1,igrid)
           cd2 = -eoffin(ixd2,if1,igrid)
           ! Loop over edges involved in approximate M-inverse operator
           DO ixm = 1, nxminvsten(ie1,igrid)
              ie2 = xminvsten(ixm,ie1,igrid)
              cm = nusq(ie1,igrid)*xminv(ixm,ie1,igrid)
              ! Loop over cells involved in Ddual1 operator
              cd1 = 1.0d0
              DO ixd1 = 1, 2
                 if2 = fnxte(ixd1,ie2,igrid)
                 cd1 = -cd1
                 ! Loop over cells in L operator
                 DO ixl = 1, nlsten(if2,igrid)
                    if3 = lsten(ixl,if2,igrid)
                    cl = lmass(ixl,if2,igrid)
                    IF (if3 == if1) temp = temp + cd2*cm*cd1*cl
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        helmdiag(if1,igrid) = temp

     ENDDO
     !$OMP END PARALLEL DO
  ENDDO


  DEALLOCATE(temp1, temp2)

  ! ---------------------------------------------------------------

END SUBROUTINE buildhelm

! ================================================================

SUBROUTINE buildadvsten

  ! Build stencils for advection scheme

  USE advection
  IMPLICIT NONE

  INTEGER :: if0, if1, if2, nlist, list(100), ixs, ixs2, ixn, ixl, &
       iv0, iv1, iv2, ie1
  LOGICAL :: alreadysten, alreadylist, lfound

  ! ---------------------------------------------------------------

  ! PRIMAL GRID

  PRINT *,'*** nexadvf frozen at 1 ***'

  ! Loop over faces
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)&
  !$OMP PRIVATE(if0, if1, if2, nlist, ixs, ixs2, list, ixn, ixl, iv0, iv1, iv2, ie1, alreadysten, alreadylist, lfound)&
  !$OMP SHARED(nstenadvf, stenadvf, nexadvf, nmonomial, neoff, fnxtf, ngrids, nface)
  DO if0 = 1, nface(ngrids)

     ! Start with the face itself
     nstenadvf(if0) = 1
     stenadvf(1,if0) = if0
     ! And demand that this cell should be fitted exactly
     nexadvf(if0) = 1

     ! Iteratively expand stencil until the number of cells
     ! is at least the number of monomials
     DO WHILE (nstenadvf(if0) < nmonomial)

        ! Current cells must be fitted exactly
        ! THIS OPTION SUSPENDED PENDING FURTHER INVESTIGATION
        ! SEE REPORT OF 15 MARCH 2012
        ! nexadvf(if0) = nstenadvf(if0)

        ! And look for some new ones
        nlist = 0
        list = 0

        ! Find neighbours of cells currently in the stencil
        DO ixs = 1, nstenadvf(if0)
           if1 = stenadvf(ixs,if0)
           DO ixn = 1, neoff(if1,ngrids)
              if2 = fnxtf(if1,ixn,ngrids)

              ! Is it already in the stencil?
              alreadysten = .FALSE.
              DO ixs2 = 1, nstenadvf(if0)
                 IF (if2 == stenadvf(ixs2,if0)) alreadysten = .TRUE.
              ENDDO

              IF (.NOT. alreadysten) THEN
                 ! If it's already on the temporary list, flag the fact that
                 ! we've seen it more than once
                 alreadylist = .FALSE.
                 DO ixl = 1, nlist
                    IF (if2 == ABS(list(ixl))) THEN
                       ! It's already on the list; make a note, and store as
                       ! a positive number 
                       alreadylist = .TRUE.
                       list(ixl) = if2
                    ENDIF
                 ENDDO
                 IF (.NOT. alreadylist) THEN
                    ! It's not already on the list; add it as a negative number
                    ! to indicate this is the first time
                    nlist = nlist + 1
                    list(nlist) = -if2
                 ENDIF
              ENDIF

           ENDDO
        ENDDO

        ! If we found any that are neighbours more than once then take them
        lfound = .FALSE.
        DO ixl = 1, nlist
           IF (list(ixl) > 0) THEN
              nstenadvf(if0) = nstenadvf(if0) + 1
              stenadvf(nstenadvf(if0),if0) = list(ixl)
              lfound = .TRUE.
           ENDIF
        ENDDO

        ! Otherwise, take those that are neighbours just once
        IF (.NOT. lfound) THEN
           DO ixl = 1, nlist
              nstenadvf(if0) = nstenadvf(if0) + 1
              stenadvf(nstenadvf(if0),if0) = ABS(list(ixl))
           ENDDO
        ENDIF

     ENDDO

  ENDDO
  !$OMP END PARALLEL DO

  IF (MAXVAL(nstenadvf) .NE. nadvfmx) THEN
     PRINT *,'nadvfmx = ',nadvfmx,' but MAXVAL(nstenadvf) = ',MAXVAL(nstenadvf)
     STOP
  ENDIF


  ! DUAL GRID

  PRINT *,'*** nexadvv frozen at 1 ***'

  ! Loop over dual cells
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)&
  !$OMP PRIVATE(if0, if1, if2, nlist, ixs, ixs2, list, ixn, ixl, iv0, iv1, iv2, ie1, alreadysten, alreadylist, lfound)&
  !$OMP SHARED(nstenadvv, stenadvv, nexadvv, nmonomial, neoff, fnxtf, ngrids, neofv, eofv, vofe, nvert)
  DO iv0 = 1, nvert(ngrids)

     ! Start with the dual cell itself
     nstenadvv(iv0) = 1
     stenadvv(1,iv0) = iv0
     ! And demand that this cell should be fitted exactly
     nexadvv(iv0) = 1

     ! Iteratively expand stencil until the number of cells
     ! is at least the number of monomials
     DO WHILE (nstenadvv(iv0) < nmonomial)

        ! Current cells must be fitted exactly
        ! THIS OPTION SUSPENDED PENDING FURTHER INVESTIGATION
        ! SEE REPORT OF 15 MARCH 2012
        ! nexadvv(iv0) = nstenadvv(iv0)

        ! And look for some new ones
        nlist = 0
        list = 0

        ! Find neighbours of cells currently in the stencil
        DO ixs = 1, nstenadvv(iv0)
           iv1 = stenadvv(ixs,iv0)
           DO ixn = 1, neofv(iv1,ngrids)
              ! We don't store the vertices neighbouring a given vertex,
              ! so look along edges of the vertex to find neighbours
              ie1 = eofv(ixn,iv1,ngrids)
              iv2 = vofe(1,ie1,ngrids)
              IF (iv2 == iv1) iv2 = vofe(2,ie1,ngrids)

              ! Is it already in the stencil?
              alreadysten = .FALSE.
              DO ixs2 = 1, nstenadvv(iv0)
                 IF (iv2 == stenadvv(ixs2,iv0)) alreadysten = .TRUE.
              ENDDO

              IF (.NOT. alreadysten) THEN
                 ! If it's already on the temporary list, flag the fact that
                 ! we've seen it more than once
                 alreadylist = .FALSE.
                 DO ixl = 1, nlist
                    IF (iv2 == ABS(list(ixl))) THEN
                       ! It's already on the list; make a note, and store as
                       ! a positive number 
                       alreadylist = .TRUE.
                       list(ixl) = iv2
                    ENDIF
                 ENDDO
                 IF (.NOT. alreadylist) THEN
                    ! It's not already on the list; add it as a negative number
                    ! to indicate this is the first time
                    nlist = nlist + 1
                    list(nlist) = -iv2
                 ENDIF
              ENDIF

           ENDDO
        ENDDO

        ! If we found any that are neighbours more than once then take them
        lfound = .FALSE.
        DO ixl = 1, nlist
           IF (list(ixl) > 0) THEN
              nstenadvv(iv0) = nstenadvv(iv0) + 1
              stenadvv(nstenadvv(iv0),iv0) = list(ixl)
              lfound = .TRUE.
           ENDIF
        ENDDO

        ! Otherwise, take those that are neighbours just once
        IF (.NOT. lfound) THEN
           DO ixl = 1, nlist
              nstenadvv(iv0) = nstenadvv(iv0) + 1
              stenadvv(nstenadvv(iv0),iv0) = ABS(list(ixl))
           ENDDO
        ENDIF

     ENDDO

  ENDDO
  !$OMP END PARALLEL DO

  IF (MAXVAL(nstenadvv) .NE. nadvvmx) THEN
     PRINT *,'nadvvmx = ',nadvvmx,' but MAXVAL(nstenadvv) = ',MAXVAL(nstenadvv)
     STOP
  ENDIF


  ! ---------------------------------------------------------------

END SUBROUTINE buildadvsten

! ================================================================

SUBROUTINE buildintmon

  ! Build integrals of monomials over all cells in advection stencils

  ! *** Note: it may be better to build these on the fly rather than
  ! store and retrieve them. ***

  USE constants
  USE advection
  IMPLICIT NONE

  INTEGER :: if0, if1, if2, ixs, ixt, ie1, iv0, iv1, iv2, igauss, &
       px, py, m

  REAL*8 :: long, lat, x0, y0, z0, x1, y1, z1, xn1, yn1, zn1, &
       xn2, yn2, zn2, mag, area, wgt, p(3), xg, yg, zg, &
       s, sinth, costh, xx, yy, &
       fn, x2, y2, z2, x3, y3, z3, x4, y4, z4

  ! ---------------------------------------------------------------

  ! PRIMAL GRID


  ! Loop over all faces
!!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
!!$OMP PRIVATE(if1, if2, ixs, ixt, ie1, iv0, iv1, iv2, igauss, px, py, &
!!$OMP   m, long, lat, x0, y0, z0, x1, y1, z1, xn1, yn1, zn1, xn2, yn2, zn2, &
!!$OMP   mag, area, wgt, p, xg, yg, zg, s, sinth, costh, xx, yy, fn, x2, y2, &
!!$OMP   z2, x3, y3, z3, x4, y4, z4, rearth) &
!!$OMP SHARED(nstenadvf, stenadvf, nexadvf, intmonf, nmonomial, eoff, neoff, &
!!$OMP   vlong, vlat, fnxtf, ngrids, neofv, eofv, vofe, flong, flat, nface)
  DO if0 = 1, nface(ngrids)
     ! Initialize to zero
     intmonf(:,if0,:) = 0

     ! Position vector of centre of face if0
     long = flong(if0,ngrids)
     lat = flat(if0,ngrids)
     CALL ll2xyz(long,lat,x0,y0,z0)

     ! Find direction of first neighbour to establish axes of
     ! Local coordinate system
     if1 = fnxtf(if0,1,ngrids)
     ! Position vector of face if1
     long = flong(if1,ngrids)
     lat = flat(if1,ngrids)
     CALL ll2xyz(long,lat,x1,y1,z1)

     ! Unit normal to plane containing points 0 and 1
     xn1 = y0*z1 - z0*y1
     yn1 = z0*x1 - x0*z1
     zn1 = x0*y1 - y0*x1
     mag = SQRT(xn1*xn1 + yn1*yn1 + zn1*zn1)
     xn1 = xn1/mag
     yn1 = yn1/mag
     zn1 = zn1/mag

     ! Loop over all cells in stencils
     DO ixs = 1, nstenadvf(if0)

        if2 = stenadvf(ixs,if0)

        ! Centre of this cell
        long = flong(if2,ngrids)
        lat = flat(if2,ngrids)
        CALL ll2xyz(long,lat,x2,y2,z2)

        ! Subdivide cell into triangles
        ! (This is overkill if the cell is itself a triangle)
        DO ixt = 1, neoff(if2,ngrids)
           ie1 = eoff(ixt,if2,ngrids)
           iv1 = vofe(1,ie1,ngrids)
           long = vlong(iv1,ngrids)
           lat = vlat(iv1,ngrids)
           CALL ll2xyz(long,lat,x3,y3,z3)
           iv2 = vofe(2,ie1,ngrids)
           long = vlong(iv2,ngrids)
           lat = vlat(iv2,ngrids)
           CALL ll2xyz(long,lat,x4,y4,z4)
           ! Area of this triangle
           CALL starea2(x2,y2,z2,x3,y3,z3,x4,y4,z4,area)
           wgt = area*rearth*rearth/3.0d0
           ! Loop over Gauss points for this triangle
           DO igauss = 1, 3
              p(1:3) = 1.0d0/6.0d0
              p(igauss) = 4.0d0/6.0d0
              xg = p(1)*x2 + p(2)*x3 + p(3)*x4
              yg = p(1)*y2 + p(2)*y3 + p(3)*y4
              zg = p(1)*z2 + p(2)*z3 + p(3)*z4
              mag = SQRT(xg*xg + yg*yg + zg*zg)
              xg = xg/mag
              yg = yg/mag
              zg = zg/mag
              ! Express Gauss point in terms of face if0's local coordinates
              ! Distance
              xn2 = y0*zg - z0*yg
              yn2 = z0*xg - x0*zg
              zn2 = x0*yg - y0*xg
              mag = SQRT(xn2*xn2 + yn2*yn2 + zn2*zn2)
              s = ASIN(mag)*rearth
              ! Unit normal
              xn2 = xn2/mag
              yn2 = yn2/mag
              zn2 = zn2/mag
              ! Angle relative to local x-axis
              sinth = x0*(yn1*zn2 - zn1*yn2) &
                   + y0*(zn1*xn2 - xn1*zn2) &
                   + z0*(xn1*yn2 - yn1*xn2)
              costh = xn1*xn2 + yn1*yn2 + zn1*zn2
              ! Finally obtain local coordinates
              xx = s*costh
              yy = s*sinth
              ! Loop over monomials
              px = 0
              py = 0
              DO m = 1, nmonomial
                 fn = (xx**px)*(yy**py)
                 intmonf(ixs,if0,m) = intmonf(ixs,if0,m) + wgt*fn
                 px = px - 1
                 py = py + 1
                 IF (px < 0) THEN
                    px = py
                    py = 0
                 ENDIF
              ENDDO      ! End loop over monomials
           ENDDO        ! End loop over Gauss points
        ENDDO          ! End loop over sub-triangles

     ENDDO            ! End loop over cells in stencil

  ENDDO              ! End loop over faces
!!$OMP END PARALLEL DO


  ! ---------------------------------------------------------------

  ! DUAL GRID

  ! Initialize to zero
  intmonv = 0.0d0

  ! Loop over all dual
  DO iv0 = 1, nvert(ngrids)

     ! Position vector of centre of dual cell iv0
     long = vlong(iv0,ngrids)
     lat = vlat(iv0,ngrids)
     CALL ll2xyz(long,lat,x0,y0,z0)

     ! Find direction of first neighbour to establish axes of
     ! Local coordinate system
     ie1 = eofv(1,iv0,ngrids)
     iv1 = vofe(1,ie1,ngrids)
     IF (iv1 == iv0) iv1 = vofe(2,ie1,ngrids)
     ! Position vector of dual cell iv1
     long = vlong(iv1,ngrids)
     lat = vlat(iv1,ngrids)
     CALL ll2xyz(long,lat,x1,y1,z1)

     ! Unit normal to plane containing points 0 and 1
     xn1 = y0*z1 - z0*y1
     yn1 = z0*x1 - x0*z1
     zn1 = x0*y1 - y0*x1
     mag = SQRT(xn1*xn1 + yn1*yn1 + zn1*zn1)
     xn1 = xn1/mag
     yn1 = yn1/mag
     zn1 = zn1/mag

     ! Loop over all dual cells in stencils
     DO ixs = 1, nstenadvv(iv0)

        iv2 = stenadvv(ixs,iv0)

        ! Centre of this dual cell
        long = vlong(iv2,ngrids)
        lat = vlat(iv2,ngrids)
        CALL ll2xyz(long,lat,x2,y2,z2)

        ! Subdivide cell into triangles
        ! (This is overkill if the cell is itself a triangle)
        DO ixt = 1, neofv(iv2,ngrids)
           ie1 = eofv(ixt,iv2,ngrids)
           if1 = fnxte(1,ie1,ngrids)
           long = flong(if1,ngrids)
           lat = flat(if1,ngrids)
           CALL ll2xyz(long,lat,x3,y3,z3)
           if2 = fnxte(2,ie1,ngrids)
           long = flong(if2,ngrids)
           lat = flat(if2,ngrids)
           CALL ll2xyz(long,lat,x4,y4,z4)
           ! Area of this triangle
           CALL starea2(x2,y2,z2,x3,y3,z3,x4,y4,z4,area)
           wgt = area*rearth*rearth/3.0d0
           ! Loop over Gauss points for this triangle
           DO igauss = 1, 3
              p(1:3) = 1.0d0/6.0d0
              p(igauss) = 4.0d0/6.0d0
              xg = p(1)*x2 + p(2)*x3 + p(3)*x4
              yg = p(1)*y2 + p(2)*y3 + p(3)*y4
              zg = p(1)*z2 + p(2)*z3 + p(3)*z4
              mag = SQRT(xg*xg + yg*yg + zg*zg)
              xg = xg/mag
              yg = yg/mag
              zg = zg/mag
              ! Express Gauss point in terms of dual cell iv0's local coordinates
              ! Distance
              xn2 = y0*zg - z0*yg
              yn2 = z0*xg - x0*zg
              zn2 = x0*yg - y0*xg
              mag = SQRT(xn2*xn2 + yn2*yn2 + zn2*zn2)
              s = ASIN(mag)*rearth
              ! Unit normal
              xn2 = xn2/mag
              yn2 = yn2/mag
              zn2 = zn2/mag
              ! Angle relative to local x-axis
              sinth = x0*(yn1*zn2 - zn1*yn2) &
                   + y0*(zn1*xn2 - xn1*zn2) &
                   + z0*(xn1*yn2 - yn1*xn2)
              costh = xn1*xn2 + yn1*yn2 + zn1*zn2
              ! Finally obtain local coordinates
              xx = s*costh
              yy = s*sinth
              ! Loop over monomials
              px = 0
              py = 0
              DO m = 1, nmonomial
                 fn = (xx**px)*(yy**py)
                 intmonv(ixs,iv0,m) = intmonv(ixs,iv0,m) + wgt*fn
                 px = px - 1
                 py = py + 1
                 IF (px < 0) THEN
                    px = py
                    py = 0
                 ENDIF
              ENDDO      ! End loop over monomials
           ENDDO        ! End loop over Gauss points
        ENDDO          ! End loop over sub-triangles

     ENDDO            ! End loop over dual cells in stencil

  ENDDO              ! End loop over dual cells



  ! ---------------------------------------------------------------

END SUBROUTINE buildintmon

! ================================================================

SUBROUTINE buildadvcoeff

  ! Build the coefficients used to construct polynomial fits for
  ! advection scheme. The polynomial fit should be exact for
  ! nexadvf or nexadvv cells in the stencil, and the mean square
  ! residual over the other cells should be minimized. This leads to
  ! to a linear problem for the coefficients of the polynomial fit
  ! and some Lagrange multipliers. Finding the inverse of the matrix
  ! involved allows us to save the matrix that gives the polynomial
  ! coefficients in terms of the cell integrals of the advected quantity.

  USE advection
  IMPLICIT NONE

  INTEGER :: if0, iv0, i, ns, ne, nm
  REAL*8, ALLOCATABLE :: l(:,:), lt(:,:), ltl(:,:), &
       m(:,:), minv(:,:), r(:,:), q(:,:)

  ! ---------------------------------------------------------------

  ! PRIMAL GRID

  ! Loop over faces
  DO if0 = 1, nface(ngrids)

     ! Determine the stencil size and the number of cells fitted
     ! exactly; hence allocate space for matrices
     ns = nstenadvf(if0)
     ne = nexadvf(if0)
     nm = nmonomial + ne
     ALLOCATE(l(ns,nmonomial), lt(nmonomial,ns), ltl(nmonomial,nmonomial))
     ALLOCATE(m(nm,nm), minv(nm,nm), r(nm,ns), q(nm,ns))

     ! Extract the matrix of integrals of monomials
     l = intmonf(1:ns,if0,1:nmonomial)
     lt = TRANSPOSE(l)
     ltl = MATMUL(lt,l)

     ! Build matrix for linear system
     m(1:nmonomial,1:nmonomial) = ltl
     m(1:nmonomial,nmonomial+1:nm) = lt(1:nmonomial,1:ne)
     m(nmonomial+1:nm,1:nmonomial) = l(1:ne,1:nmonomial)
     m(nmonomial+1:nm,nmonomial+1:nm) = 0.0d0

     ! Invert matrix
     CALL matinv(m,minv,nm)

     ! Matrix relating cell integrals to RHS of linear system
     r(1:nmonomial,1:ns) = lt(1:nmonomial,1:ns)
     r(nmonomial+1:nm,1:ns) = 0.0d0
     DO i = 1,ne
        r(nmonomial+i,i) = 1.0d0
     ENDDO

     ! Matrix giving polynomial coefficients in terms of cell integrals
     q = MATMUL(minv,r)

     ! Save the part we need
     lt = q(1:nmonomial,1:ns)
     l = TRANSPOSE(lt)
     intmonf(1:ns,if0,1:nmonomial) = l

     DEALLOCATE(l, lt, ltl, m, minv, r, q)

  ENDDO


  ! DUAL GRID

  ! Loop over dual cells
  DO iv0 = 1, nvert(ngrids)

     ! Determine the stencil size and the number of cells fitted
     ! exactly; hence allocate space for matrices
     ns = nstenadvv(iv0)
     ne = nexadvv(iv0)
     nm = nmonomial + ne
     ALLOCATE(l(ns,nmonomial), lt(nmonomial,ns), ltl(nmonomial,nmonomial))
     ALLOCATE(m(nm,nm), minv(nm,nm), r(nm,ns), q(nm,ns))

     ! Extract the matrix of integrals of monomials
     l = intmonv(1:ns,iv0,1:nmonomial)
     lt = TRANSPOSE(l)
     ltl = MATMUL(lt,l)

     ! Build matrix for linear system
     m(1:nmonomial,1:nmonomial) = ltl
     m(1:nmonomial,nmonomial+1:nm) = lt(1:nmonomial,1:ne)
     m(nmonomial+1:nm,1:nmonomial) = l(1:ne,1:nmonomial)
     m(nmonomial+1:nm,nmonomial+1:nm) = 0.0d0

     ! Invert matrix
     CALL matinv(m,minv,nm)

     ! Matrix relating cell integrals to RHS of linear system
     r(1:nmonomial,1:ns) = lt(1:nmonomial,1:ns)
     r(nmonomial+1:nm,1:ns) = 0.0d0
     DO i = 1,ne
        r(nmonomial+i,i) = 1.0d0
     ENDDO

     ! Matrix giving polynomial coefficients in terms of cell integrals
     q = MATMUL(minv,r)

     ! Save the part we need
     lt = q(1:nmonomial,1:ns)
     l = TRANSPOSE(lt)
     intmonv(1:ns,iv0,1:nmonomial) = l

     DEALLOCATE(l, lt, ltl, m, minv, r, q)

  ENDDO


  ! ---------------------------------------------------------------

END SUBROUTINE buildadvcoeff

! ================================================================

SUBROUTINE setupgauss

  ! Initialize Gauss points and weights for integrating
  ! over swept areas in advection scheme

  USE advection
  IMPLICIT NONE
  REAL*8 :: rr3, r3by5

  ! ---------------------------------------------------------------

  ! Only lower order versions required so simply assign
  ! explicitly

  IF (ngauss == 1) THEN
     xgauss(1) = 0.5d0
     wgauss(1) = 1.0d0
  ELSEIF (ngauss == 2) THEN
     rr3 = 1.0d0/SQRT(3.0D0)
     xgauss(1) = 0.5d0*(1.0d0 - rr3)
     xgauss(2) = 0.5d0*(1.0d0 + rr3)
     wgauss(1) = 0.5d0
     wgauss(2) = 0.5d0
  ELSEIF (ngauss == 3) THEN
     r3by5 = SQRT(3.0d0/5.0d0)
     xgauss(1) = 0.5d0*(1.0d0 - r3by5)
     xgauss(2) = 0.5d0
     xgauss(3) = 0.5d0*(1.0d0 + r3by5)
     wgauss(1) = 5.0d0/18.0d0
     wgauss(2) = 4.0d0/9.0d0
     wgauss(3) = 5.0d0/18.0d0
  ELSE
     PRINT *,'Need to set up Gauss points and weights for degree = ',degree
     STOP
  ENDIF


  ! ---------------------------------------------------------------

END SUBROUTINE setupgauss

! ================================================================

SUBROUTINE primaladvflx(chi2,flx1)

  ! Compute advective flux of input field.
  ! chi2 is the area integral over primal cells of the
  ! density or concentration of the field to be advected.
  ! flx1 is the net flux across primal cell edges over one time step.

  ! It is assumed that the fields u1 and uperp1 are already
  ! available.

  USE constants
  USE state
  USE work
  USE advection
  USE timestep
  IMPLICIT NONE

  REAL*8, INTENT(IN) :: chi2(nfacex)
  REAL*8, INTENT(OUT) :: flx1(nedgex)
  INTEGER :: nf, ne, ie0, if0, iv1, iv2, m, px, py, ixs, &
       if1, i, j, ns
  REAL*8 :: long, lat, x0, y0, z0, &
       x1, y1, z1, x2, y2, z2, xn1, yn1, zn1, xn2, yn2, zn2, &
       s, sinth, costh, xx1, yy1, xx2, yy2, mag, dx, dy, &
       tx, ty, nx, ny, wx, wy, udt, vdt, xg, yg, aswept, &
       flx, poly, fn, xgi, xgj, wgi, wgj, a(nmonomial), p, &
       sgn, temp
  REAL*8 :: cn, abscn, maxc

  ! ---------------------------------------------------------------

  nf = nface(ngrids)
  ne = nedge(ngrids)

  maxc = 0.0d0

  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
  !$OMP PRIVATE(  &
  !$OMP    nf, ie0, if0, iv1, iv2, m, px, py, ixs, if1, i, j, ns, abscn, &
  !$OMP    long, lat, x0, y0, z0, x1, y1, z1, x2, y2, z2, xn1, yn1, zn1, xn2, yn2, zn2, s, sinth, costh, xx1, yy1, xx2, yy2, &
  !$OMP    mag, dx, dy, tx, ty, nx, ny, wx, wy, udt, vdt, xg, yg, aswept, flx, poly, fn, xgi, xgj, wgi, wgj, a, p, sgn, temp, cn)&
  !$OMP SHARED(ne, ubar1, fnxte, vofe, xgauss, uperpbar1, wgauss, vlat, rearth, vlong, ldist, ngauss, flx1, flong, flat, nmonomial, &
  !$OMP    nstenadvf, stenadvf, chi2, farea, intmonf, divfac, ngrids, fnxtf, dt, maxc)
  DO ie0 = 1, ne

     ! Decide which is the upwind cell
     ! and orient the ends of the edge
     IF (ubar1(ie0) > 0.0d0) THEN
        if0 = fnxte(1,ie0,ngrids)
        iv1 = vofe(1,ie0,ngrids)
        iv2 = vofe(2,ie0,ngrids)
        sgn = 1.0d0
     ELSE
        if0 = fnxte(2,ie0,ngrids)
        iv1 = vofe(2,ie0,ngrids)
        iv2 = vofe(1,ie0,ngrids)
        sgn = -1.0d0
     ENDIF

     ! *** just for diagnostics ***
     cn = ubar1(ie0)*dt/farea(if0,ngrids)
     abscn = ABS(cn)

     !$OMP ATOMIC UPDATE
     maxc = MAX(abscn,maxc)

     ! Build the polynomial fit for cell if0
     ns = nstenadvf(if0)
     a = 0.0d0
     DO ixs = 1, ns
        if1 = stenadvf(ixs,if0)
        p = chi2(if1)
        DO m = 1, nmonomial
           a(m) = a(m) + intmonf(ixs,if0,m)*p
        ENDDO
     ENDDO

     ! Centre of upwind cell
     long = flong(if0,ngrids)
     lat = flat(if0,ngrids)
     CALL ll2xyz(long,lat,x0,y0,z0)

     ! Find direction of first neighbour to establish axes of
     ! local coordinate system
     if1 = fnxtf(if0,1,ngrids)
     ! Position vector of face if1
     long = flong(if1,ngrids)
     lat = flat(if1,ngrids)
     CALL ll2xyz(long,lat,x1,y1,z1)

     ! Unit normal to plane containing points 0 and 1
     xn1 = y0*z1 - z0*y1
     yn1 = z0*x1 - x0*z1
     zn1 = x0*y1 - y0*x1
     mag = SQRT(xn1*xn1 + yn1*yn1 + zn1*zn1)
     xn1 = xn1/mag
     yn1 = yn1/mag
     zn1 = zn1/mag

     ! Position vectors of edge ends,
     ! expressed in local coordinates.
     ! Note: these could be precomputed and stored if that
     ! was more efficient.
     ! First one:
     long = vlong(iv1,ngrids)
     lat = vlat(iv1,ngrids)
     CALL ll2xyz(long,lat,x2,y2,z2)
     ! Distance
     xn2 = y0*z2 - z0*y2
     yn2 = z0*x2 - x0*z2
     zn2 = x0*y2 - y0*x2
     mag = SQRT(xn2*xn2 + yn2*yn2 + zn2*zn2)
     s = ASIN(mag)*rearth
     ! Unit normal
     xn2 = xn2/mag
     yn2 = yn2/mag
     zn2 = zn2/mag
     ! Angle relative to local x-axis
     sinth = x0*(yn1*zn2 - zn1*yn2) &
          + y0*(zn1*xn2 - xn1*zn2) &
          + z0*(xn1*yn2 - yn1*xn2)
     costh = xn1*xn2 + yn1*yn2 + zn1*zn2
     ! Finally obtain local coordinates
     xx1 = s*costh
     yy1 = s*sinth
     ! Second one:
     long = vlong(iv2,ngrids)
     lat = vlat(iv2,ngrids)
     CALL ll2xyz(long,lat,x2,y2,z2)
     ! Distance
     xn2 = y0*z2 - z0*y2
     yn2 = z0*x2 - x0*z2
     zn2 = x0*y2 - y0*x2
     mag = SQRT(xn2*xn2 + yn2*yn2 + zn2*zn2)
     s = ASIN(mag)*rearth
     ! Unit normal
     xn2 = xn2/mag
     yn2 = yn2/mag
     zn2 = zn2/mag
     ! Angle relative to local x-axis
     sinth = x0*(yn1*zn2 - zn1*yn2) &
          + y0*(zn1*xn2 - xn1*zn2) &
          + z0*(xn1*yn2 - yn1*xn2)
     costh = xn1*xn2 + yn1*yn2 + zn1*zn2
     ! Finally obtain local coordinates
     xx2 = s*costh
     yy2 = s*sinth

     ! Unit normal and tangent expressed in local coordinates
     dx = xx2 - xx1
     dy = yy2 - yy1
     mag = SQRT(dx*dx + dy*dy)
     tx = dx/mag
     ty = dy/mag
     nx = ty
     ny = -tx

     ! Swept area
     aswept = ubar1(ie0)*dt*divfac(if0)

     ! Displacement in normal and tangential directions
     temp = sgn*dt/ldist(ie0,ngrids)
     udt = ubar1(ie0)*temp
     vdt = uperpbar1(ie0)*temp

     ! Displacement parallel to velocity
     wx = udt*nx + vdt*tx
     wy = udt*ny + vdt*ty

     ! Integrate over swept area
     flx = 0.0d0
     ! Loop over gauss points
     DO i = 1, ngauss
        xgi = xgauss(i)
        wgi = wgauss(i)
        DO j = 1, ngauss
           xgj = xgauss(j)
           wgj = wgauss(j)

           ! Gauss point in local coordinates
           xg = xx1 - xgi*wx + xgj*dx
           yg = yy1 - xgi*wy + xgj*dy

           ! Evaluate polynomial fit
           ! Loop over monomials
           poly = 0.0d0
           px = 0
           py = 0
           DO m = 1, nmonomial
              fn = (xg**px)*(yg**py)
              poly = poly + a(m)*fn
              px = px - 1
              py = py + 1
              IF (px < 0) THEN
                 px = py
                 py = 0
              ENDIF
           ENDDO
           flx = flx + poly*wgi*wgj
        ENDDO
     ENDDO
     flx1(ie0) = flx*aswept

  ENDDO
  !$OMP END PARALLEL DO

  ! print *,'Max Courant No. (primal) = ',maxc

  ! ---------------------------------------------------------------

END SUBROUTINE primaladvflx

! ================================================================

SUBROUTINE dualadvflx(chi2,flx1)

  ! Compute advective flux of input field.
  ! chi2 is the area integral over dual cells of the
  ! mixing ratio of the field to be advected.
  ! flx1 is the net flux across dual cell edges over one time step.

  ! It is assumed that the fields v1, vperp1 and mfperp1 are already
  ! available.

  USE constants
  USE state
  USE work
  USE advection
  USE timestep
  IMPLICIT NONE

  REAL*8, INTENT(IN) :: chi2(nvertx)
  REAL*8, INTENT(OUT) :: flx1(nedgex)
  INTEGER :: nv, ne, ie0, iv0, if1, if2, m, px, py, ixs, &
       iv1, i, j, ns, ie1
  REAL*8 :: long, lat, x0, y0, z0, &
       x1, y1, z1, x2, y2, z2, xn1, yn1, zn1, xn2, yn2, zn2, &
       s, sinth, costh, xx1, yy1, xx2, yy2, mag, dx, dy, &
       tx, ty, nx, ny, wx, wy, udt, vdt, xg, yg, mswept, &
       flx, poly, fn, xgi, xgj, wgi, wgj, a(nmonomial), p, &
       sgn, temp
  REAL*8 :: cn, abscn, maxc

  ! ---------------------------------------------------------------

  nv = nvert(ngrids)
  ne = nedge(ngrids)

  maxc = 0.0d0

  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
  !$OMP PRIVATE(  &
  !$OMP    nv, iv0, if1, if2, m, px, py, ixs, iv1, i, j, ns, ie1, &
  !$OMP    long, lat, x0, y0, z0, &
  !$OMP    x1, y1, z1, x2, y2, z2, xn1, yn1, zn1, xn2, yn2, zn2, s, sinth, costh, xx1, yy1, xx2, yy2, mag, dx, dy, &
  !$OMP    tx, ty, nx, ny, wx, wy, udt, vdt, xg, yg, mswept, flx, poly, fn, xgi, xgj, wgi, wgj, a, p, sgn, temp, &
  !$OMP    cn, abscn)&
  !$OMP SHARED(ne, ubar1, fnxte, vofe, xgauss, uperpbar1, wgauss, vlat, rearth, vlong, ldist, ngauss, flx1, flong, flat, nmonomial, &
  !$OMP    nstenadvf, stenadvf, chi2, farea, intmonf, divfac, ngrids, fnxtf, dt, maxc, mfperp1, &
  !$OMP    vperpbar1, jstar, nstenadvv, stenadvv, eofv, ddist, intmonv, vbar1)
  DO ie0 = 1, ne

     ! Decide which is the upwind dual cell
     ! and orient the ends of the edge
     IF (mfperp1(ie0) > 0.0d0) THEN
        iv0 = vofe(1,ie0,ngrids)
        if1 = fnxte(2,ie0,ngrids)
        if2 = fnxte(1,ie0,ngrids)
        sgn = 1.0d0
     ELSE
        iv0 = vofe(2,ie0,ngrids)
        if1 = fnxte(1,ie0,ngrids)
        if2 = fnxte(2,ie0,ngrids)
        sgn = -1.0d0
     ENDIF

     ! *** just for diagnostics ***
     cn = vperpbar1(ie0)*dt*jstar(1,iv0,ngrids)
     abscn = ABS(cn)

     !$OMP ATOMIC UPDATE
     maxc = MAX(abscn,maxc)

     ! Build the polynomial fit for dual cell iv0
     ns = nstenadvv(iv0)
     a = 0.0d0
     DO ixs = 1, ns
        iv1 = stenadvv(ixs,iv0)
        p = chi2(iv1)
        DO m = 1, nmonomial
           a(m) = a(m) + intmonv(ixs,iv0,m)*p
        ENDDO
     ENDDO

     ! Centre of upwind cell
     long = vlong(iv0,ngrids)
     lat = vlat(iv0,ngrids)
     CALL ll2xyz(long,lat,x0,y0,z0)

     ! Find direction of first neighbour to establish axes of
     ! local coordinate system
     ie1 = eofv(1,iv0,ngrids)
     iv1 = vofe(1,ie1,ngrids)
     IF (iv1 == iv0) iv1 = vofe(2,ie1,ngrids)
     ! Position vector of dual cell iv1
     long = vlong(iv1,ngrids)
     lat = vlat(iv1,ngrids)
     CALL ll2xyz(long,lat,x1,y1,z1)

     ! Unit normal to plane containing points 0 and 1
     xn1 = y0*z1 - z0*y1
     yn1 = z0*x1 - x0*z1
     zn1 = x0*y1 - y0*x1
     mag = SQRT(xn1*xn1 + yn1*yn1 + zn1*zn1)
     xn1 = xn1/mag
     yn1 = yn1/mag
     zn1 = zn1/mag

     ! Position vectors of edge ends,
     ! expressed in local coordinates.
     ! Note: these could be precomputed and stored if that
     ! was more efficient.
     ! First one:
     long = flong(if1,ngrids)
     lat = flat(if1,ngrids)
     CALL ll2xyz(long,lat,x2,y2,z2)
     ! Distance
     xn2 = y0*z2 - z0*y2
     yn2 = z0*x2 - x0*z2
     zn2 = x0*y2 - y0*x2
     mag = SQRT(xn2*xn2 + yn2*yn2 + zn2*zn2)
     s = ASIN(mag)*rearth
     ! Unit normal
     xn2 = xn2/mag
     yn2 = yn2/mag
     zn2 = zn2/mag
     ! Angle relative to local x-axis
     sinth = x0*(yn1*zn2 - zn1*yn2) &
          + y0*(zn1*xn2 - xn1*zn2) &
          + z0*(xn1*yn2 - yn1*xn2)
     costh = xn1*xn2 + yn1*yn2 + zn1*zn2
     ! Finally obtain local coordinates
     xx1 = s*costh
     yy1 = s*sinth
     ! Second one:
     long = flong(if2,ngrids)
     lat = flat(if2,ngrids)
     CALL ll2xyz(long,lat,x2,y2,z2)
     ! Distance
     xn2 = y0*z2 - z0*y2
     yn2 = z0*x2 - x0*z2
     zn2 = x0*y2 - y0*x2
     mag = SQRT(xn2*xn2 + yn2*yn2 + zn2*zn2)
     s = ASIN(mag)*rearth
     ! Unit normal
     xn2 = xn2/mag
     yn2 = yn2/mag
     zn2 = zn2/mag
     ! Angle relative to local x-axis
     sinth = x0*(yn1*zn2 - zn1*yn2) &
          + y0*(zn1*xn2 - xn1*zn2) &
          + z0*(xn1*yn2 - yn1*xn2)
     costh = xn1*xn2 + yn1*yn2 + zn1*zn2
     ! Finally obtain local coordinates
     xx2 = s*costh
     yy2 = s*sinth

     ! Unit normal and tangent expressed in local coordinates
     dx = xx2 - xx1
     dy = yy2 - yy1
     mag = SQRT(dx*dx + dy*dy)
     tx = dx/mag
     ty = dy/mag
     nx = ty
     ny = -tx

     ! Swept mass
     mswept = mfperp1(ie0)

     ! Displacement in normal and tangential directions
     temp = sgn*dt/ddist(ie0,ngrids)
     udt = vperpbar1(ie0)*temp
     vdt = -vbar1(ie0)*temp

     ! Displacement parallel to velocity
     wx = udt*nx + vdt*tx
     wy = udt*ny + vdt*ty

     ! Integrate over swept area
     flx = 0.0d0
     ! Loop over gauss points
     DO i = 1, ngauss
        xgi = xgauss(i)
        wgi = wgauss(i)
        DO j = 1, ngauss
           xgj = xgauss(j)
           wgj = wgauss(j)

           ! Gauss point in local coordinates
           xg = xx1 - xgi*wx + xgj*dx
           yg = yy1 - xgi*wy + xgj*dy

           ! Evaluate polynomial fit
           ! Loop over monomials
           poly = 0.0d0
           px = 0
           py = 0
           DO m = 1, nmonomial
              fn = (xg**px)*(yg**py)
              poly = poly + a(m)*fn
              px = px - 1
              py = py + 1
              IF (px < 0) THEN
                 px = py
                 py = 0
              ENDIF
           ENDDO
           flx = flx + poly*wgi*wgj
        ENDDO
     ENDDO
     flx1(ie0) = flx*mswept

  ENDDO
  !$OMP END PARALLEL DO

  !print *,'Max Courant No. (dual) =   ',maxc

  ! ---------------------------------------------------------------

END SUBROUTINE dualadvflx

! ================================================================

SUBROUTINE Dprimal1(f,df,igrid,nv,ne)

  ! To compute the exterior derivative df of the field f
  ! on primal grid number igrid. f comprises pointwise values
  ! at vertices; df comprises integrals of the derivative
  ! of f along primal cell edges.

  USE grid
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: igrid, nv, ne
  REAL*8, INTENT(IN) :: f(nv)
  REAL*8, INTENT(OUT) :: df(ne)
  INTEGER :: ie1, iv1, iv2

  ! ----------------------------------------------------------------

  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(ie1, iv1, iv2) SHARED(ne, vofe, igrid, df, f)
  DO ie1 = 1, ne
     iv1 = vofe(1,ie1,igrid)
     iv2 = vofe(2,ie1,igrid)
     df(ie1) = f(iv2) - f(iv1)
  ENDDO
  !$OMP END PARALLEL DO

  ! ----------------------------------------------------------------

END SUBROUTINE Dprimal1

! ================================================================

SUBROUTINE Dprimal2(f,df,igrid,ne,nf)

  ! To compute the exterior derivative df of the field f
  ! on primal grid number igrid. f comprises integrals along
  ! primal edges; df comprises integrals of the derivative
  ! over primal cells.

  USE grid
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: igrid, ne, nf
  REAL*8, INTENT(IN) :: f(ne)
  REAL*8, INTENT(OUT) :: df(nf)
  INTEGER :: if1, ix, ie1
  REAL*8 :: temp

  ! ----------------------------------------------------------------

  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(if1, ix, ie1, temp) SHARED(nf, neoff, eoff, eoffin, f, df, igrid)
  DO if1 = 1, nf
     temp = 0.0d0
     DO ix = 1, neoff(if1,igrid)
        ie1 = eoff(ix,if1,igrid)
        ! WARNING: LOOP DEPENDENCY in inner loop
        temp = temp - f(ie1)*eoffin(ix,if1,igrid)
     ENDDO
     df(if1) = temp
  ENDDO
  !$OMP END PARALLEL DO

  ! ----------------------------------------------------------------

END SUBROUTINE Dprimal2

! ================================================================

SUBROUTINE Ddual1(f,df,igrid,nf,ne)

  ! To compute the exterior derivative df of the field f
  ! on dual grid number igrid. f comprises pointwise values
  ! at face centres; df comprises integrals of the derivative
  ! of f along dual cell edges.

  USE grid
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: igrid, nf, ne
  REAL*8, INTENT(IN) :: f(nf)
  REAL*8, INTENT(OUT) :: df(ne)
  INTEGER :: ie1, if1, if2

  ! ----------------------------------------------------------------

  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(if1, if2) SHARED(fnxte, f, df, igrid, ne)
  DO ie1 = 1, ne
     if1 = fnxte(1,ie1,igrid)
     if2 = fnxte(2,ie1,igrid)
     df(ie1) = f(if2) - f(if1)
  ENDDO
  !$OMP END PARALLEL DO

  ! ----------------------------------------------------------------

END SUBROUTINE Ddual1

! ================================================================

SUBROUTINE Ddual2(f,df,igrid,ne,nv)

  ! To compute the exterior derivative df of the field f
  ! on dual grid number igrid. f comprises integrals along
  ! dual edges; df comprises integrals of the derivative
  ! over dual cells.

  USE grid
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: igrid, ne, nv
  REAL*8, INTENT(IN) :: f(ne)
  REAL*8, INTENT(OUT) :: df(nv)
  INTEGER :: iv1, ix, ie1
  REAL*8 :: temp

  ! ----------------------------------------------------------------

  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(ix, ie1, temp) SHARED(neofv, eofv, eofvin, f, df, igrid, nv)
  DO iv1 = 1, nv
     temp = 0.0d0
     DO ix = 1, neofv(iv1,igrid)
        ie1 = eofv(ix,iv1,igrid)
        temp = temp + f(ie1)*eofvin(ix,iv1,igrid)
     ENDDO
     df(iv1) = temp
  ENDDO
  !$OMP END PARALLEL DO

  ! ----------------------------------------------------------------

END SUBROUTINE Ddual2

! ================================================================

SUBROUTINE buildjlump

  ! Build a lumped version of the j matrix

  USE grid
  IMPLICIT NONE

  ! Two choices for lumped J
  ! ijlump = 1     Diagonal of J
  ! ijlump = 2     Exact answer when applied to constant input vector
  ! ijlump = 3     Exact answer when applied to the discrete representation
  !                  of a constant scalar (input vector proportional to varea)

  INTEGER, PARAMETER :: ijlump = 3
  INTEGER :: igrid, iv1, nv
  REAL*8, ALLOCATABLE :: p1(:), jp1(:)

  ! ----------------------------------------------------------------

  DO igrid = 1, ngrids

     IF (ijlump == 1) THEN

        DO iv1 = 1, nvert(igrid)
           jlump(iv1,igrid) = jstar(1,iv1,igrid)
        ENDDO

     ELSEIF (ijlump == 2) THEN

        nv = nvert(igrid)
        ALLOCATE(p1(nv), jp1(nv))
        p1 = 1.0d0
        CALL HodgeJ(p1,jp1,igrid,nv)
        jlump(1:nv,igrid) = jp1
        DEALLOCATE(p1,jp1)

     ELSEIF (ijlump == 3) THEN

        nv = nvert(igrid)
        ALLOCATE(p1(nv), jp1(nv))
        p1 = varea(1:nv,igrid)
        CALL HodgeJ(p1,jp1,igrid,nv)
        jlump(1:nv,igrid) = jp1/varea(1:nv,igrid)

        DEALLOCATE(p1,jp1)

     ELSE

        PRINT *,'Option ijlump = ',ijlump,' not available in subroutine buildjlump'
        STOP

     ENDIF

  ENDDO

  ! ----------------------------------------------------------------

END SUBROUTINE buildjlump

! ================================================================

SUBROUTINE buildmlump

  ! Build a lumped version of the M matrix

  USE grid
  IMPLICIT NONE

  ! Three choices for lumped M
  ! imlump = 1     Diagonal of M
  ! imlump = 2     Row sum of absolute values of M
  ! imlump = 3     Best fit to solid body rotation normal to each edge

  INTEGER, PARAMETER :: imlump = 3
  INTEGER :: igrid, ie1, ie2, ix, iv0, ne, nv, iv1, iv2
  REAL*8 :: temp, slon, slat, clon, clat, a1, a2, a3, num, den, &
       b1, b2, b3, c1, c2, c3, x, y, z
  REAL*8, ALLOCATABLE :: psi1(:), psi2(:), psi3(:),    &
       v1(:), v2(:), v3(:), vv(:),   &
       mv1(:), mv2(:), mv3(:)

  ! ----------------------------------------------------------------

  DO igrid = 1, ngrids

     IF (imlump == 1) THEN

        DO ie1 = 1, nedge(igrid)
           mlump(ie1,igrid) = mmass(1,ie1,igrid)
        ENDDO

     ELSEIF (imlump == 2) THEN

        DO ie1 = 1, nedge(igrid)
           temp = 0.0d0
           DO ix = 1, nmsten(ie1,igrid)
              ie2 = msten(ix,ie1,igrid)
              temp = temp + ABS(mmass(ix,ie1,igrid))
           ENDDO
           mlump(ie1,igrid) = temp
        ENDDO

     ELSEIF (imlump == 3) THEN

        nv = nvert(igrid)
        ne = nedge(igrid)
        ALLOCATE(psi1(nv), psi2(nv), psi3(nv),    &
             v1(ne), v2(ne), v3(ne), vv(ne),  &
             mv1(ne), mv2(ne), mv3(ne))
        ! Set up three solid body rotation fields
        DO 	iv0 = 1, nv
           slon = SIN(vlong(iv0,igrid))
           clon = COS(vlong(iv0,igrid))
           slat = SIN(vlat(iv0,igrid))
           clat = COS(vlat(iv0,igrid))
           psi1(iv0) = clat*clon
           psi2(iv0) = clat*slon
           psi3(iv0) = slat
        ENDDO
        CALL Dprimal1(psi1,v1,igrid,nv,ne)
        CALL Dprimal1(psi2,v2,igrid,nv,ne)
        CALL Dprimal1(psi3,v3,igrid,nv,ne)
        CALL massM(v1,mv1,igrid,ne)
        CALL massM(v2,mv2,igrid,ne)
        CALL massM(v3,mv3,igrid,ne)
        ! Now loop over edges
        DO ie1 = 1, ne
           ! Velocity field that maximizes v(ie1) is
           ! v = a1*v1 + a2*v2 + a3*v3
           ! with
           a1 = v1(ie1)
           a2 = v2(ie1)
           a3 = v3(ie1)
           den = SQRT(a1*a1 + a2*a2 + a3*a3)
           a1 = a1/den
           a2 = a2/den
           a3 = a3/den
           ! Demand that lumped matrix agrees with full M for this v
           num = a1*mv1(ie1) + a2*mv2(ie1) + a3*mv3(ie1)
           den = a1*v1(ie1)  + a2*v2(ie1)  + a3*v3(ie1)
           mlump(ie1,igrid) = num/den
           !      if (ie1 == 7) then
           !        print *,'edge ',ie1
           !        ! Find Cartesian coords of middle of edge
           !        iv1 = vofe(ie1,1,igrid)
           !        iv2 = vofe(ie1,2,igrid)
           !        x = 0.5d0*(psi1(iv1) + psi1(iv2))
           !        y = 0.5d0*(psi2(iv1) + psi2(iv2))
           !        z = 0.5d0*(psi3(iv1) + psi3(iv2))
           !        den = SQRT(x*x + y*y + z*z)
           !        ! Coefficients for SBR with axis at (x,y,z)
           !        b1 = x/den
           !        b2 = y/den
           !        b3 = z/den
           !        ! Coefficients for SBR tangential to edge ie1
           !        c1 = a2*b3 - a3*b2
           !        c2 = a3*b1 - a1*b3
           !        c3 = a1*b2 - a2*b1
           !        ! Flow normal to edge ie1
           !        vv = a1*v1 + a2*v2 + a3*v3
           !        print *,'Normal SBR'
           !        DO ix = 1, nmsten(ie1,igrid)
           !          ie2 = msten(ie1,ix,igrid)
           !          print *,'ix = ',ix,' ie2 = ',ie2,' M =',mmass(ie1,ix,igrid),' vv = ',vv(ie2)
           !        ENDDO
           !        print *,' M v = ',a1*mv1(ie1) + a2*mv2(ie1) + a3*mv3(ie1)
           !        print *,' Mlump v = ',mlump(ie1,igrid)*vv(ie1)
           !        ! Flow with axis on edge ie1
           !        vv = b1*v1 + b2*v2 + b3*v3
           !        print *,'zero SBR'
           !        DO ix = 1, nmsten(ie1,igrid)
           !          ie2 = msten(ie1,ix,igrid)
           !          print *,'ix = ',ix,' ie2 = ',ie2,' M =',mmass(ie1,ix,igrid),' vv = ',vv(ie2)
           !        ENDDO
           !        print *,' M v = ',b1*mv1(ie1) + b2*mv2(ie1) + b3*mv3(ie1)
           !        print *,' Mlump v = ',mlump(ie1,igrid)*vv(ie1)
           !        ! Flow tangential to edge ie1
           !        vv = c1*v1 + c2*v2 + c3*v3
           !        print *,'Tangential SBR'
           !        DO ix = 1, nmsten(ie1,igrid)
           !          ie2 = msten(ie1,ix,igrid)
           !          print *,'ix = ',ix,' ie2 = ',ie2,' M =',mmass(ie1,ix,igrid),' vv = ',vv(ie2)
           !        ENDDO
           !        print *,' M v = ',c1*mv1(ie1) + c2*mv2(ie1) + c3*mv3(ie1)
           !        print *,' Mlump v = ',mlump(ie1,igrid)*vv(ie1)
           !
           !        ! print *,'for v1  mv1 = ',mv1(ie1),'    mlump*v1 = ',mlump(ie1,igrid)*v1(ie1)
           !        ! print *,'for v2  mv2 = ',mv2(ie1),'    mlump*v2 = ',mlump(ie1,igrid)*v2(ie1)
           !        ! print *,'for v3  mv3 = ',mv3(ie1),'    mlump*v3 = ',mlump(ie1,igrid)*v3(ie1)
           !      endif
        ENDDO
        DEALLOCATE(psi1,psi2,psi3,v1,v2,v3,vv,mv1,mv2,mv3)

     ELSE

        PRINT *,'Option imlump = ',imlump,' not available in subroutine buildmlump'
        STOP

     ENDIF

  ENDDO

  ! ----------------------------------------------------------------

END SUBROUTINE buildmlump

! ================================================================

SUBROUTINE buildhlump

  ! Build a lumped version of the H matrix

  USE grid
  IMPLICIT NONE

  ! Three choices for lumped H
  ! ihlump = 1     Diagonal of H
  ! ihlump = 2     Row sum of absolute values of H
  ! ihlump = 3     Best fit to global divergent flow on dual grid
  ! ihlump = 4     Best fit to solid body rotation normal to each edge

  INTEGER, PARAMETER :: ihlump = 3
  INTEGER :: igrid, ie1, ie2, ix, iv0, if0, ne, nv, nf
  REAL*8 :: temp, slon, slat, clon, clat, a1, a2, a3, num, den
  REAL*8, ALLOCATABLE :: psi1(:), psi2(:), psi3(:),  &
       v1(:), v2(:), v3(:),        &
       hv1(:), hv2(:), hv3(:)

  ! ----------------------------------------------------------------

  DO igrid = 1, ngrids

     IF (ihlump == 1) THEN

        DO ie1 = 1, nedge(igrid)
           hlump(ie1,igrid) = hstar(1,ie1,igrid)
        ENDDO

     ELSEIF (ihlump == 2) THEN

        DO ie1 = 1, nedge(igrid)
           temp = 0.0d0
           DO ix = 1, nhsten(ie1,igrid)
              ie2 = hsten(ix,ie1,igrid)
              temp = temp + ABS(hstar(ix,ie1,igrid))
           ENDDO
           hlump(ie1,igrid) = temp
        ENDDO

     ELSEIF (ihlump == 3) THEN

        nf = nface(igrid)
        ne = nedge(igrid)
        ALLOCATE(psi1(nf), psi2(nf), psi3(nf),  &
             v1(ne), v2(ne), v3(ne),        &
             hv1(ne), hv2(ne), hv3(ne))
        ! Set up three global divergent fields
        DO 	if0 = 1, nf
           slon = SIN(flong(if0,igrid))
           clon = COS(flong(if0,igrid))
           slat = SIN(flat(if0,igrid))
           clat = COS(flat(if0,igrid))
           psi1(if0) = slat
           psi2(if0) = clat*slon
           psi3(if0) = clat*clon
        ENDDO
        CALL Ddual1(psi1,v1,igrid,nf,ne)
        CALL Ddual1(psi2,v2,igrid,nf,ne)
        CALL Ddual1(psi3,v3,igrid,nf,ne)
        CALL HodgeH(v1,hv1,igrid,ne)
        CALL HodgeH(v2,hv2,igrid,ne)
        CALL HodgeH(v3,hv3,igrid,ne)
        ! Now loop over edges
        DO ie1 = 1, ne
           ! Velocity field that maximizes v(ie1) is
           ! v = a1*v1 + a2*v2 + a3*v3
           ! with
           a1 = v1(ie1)
           a2 = v2(ie1)
           a3 = v3(ie1)
           ! Demand that lumped matrix agrees with full H for this v
           num = a1*hv1(ie1) + a2*hv2(ie1) + a3*hv3(ie1)
           den = a1*v1(ie1)  + a2*v2(ie1)  + a3*v3(ie1)
           hlump(ie1,igrid) = num/den
        ENDDO
        DEALLOCATE(psi1,psi2,psi3,v1,v2,v3,hv1,hv2,hv3)

     ELSEIF (ihlump == 4) THEN

        nv = nvert(igrid)
        ne = nedge(igrid)
        ALLOCATE(psi1(nv), psi2(nv), psi3(nv),  &
             v1(ne), v2(ne), v3(ne),        &
             hv1(ne), hv2(ne), hv3(ne))
        ! Set up three solid body rotation fields
        DO 	iv0 = 1, nv
           slon = SIN(vlong(iv0,igrid))
           clon = COS(vlong(iv0,igrid))
           slat = SIN(vlat(iv0,igrid))
           clat = COS(vlat(iv0,igrid))
           psi1(iv0) = slat
           psi2(iv0) = clat*slon
           psi3(iv0) = clat*clon
        ENDDO
        CALL Dprimal1(psi1,v1,igrid,nv,ne)
        CALL Dprimal1(psi2,v2,igrid,nv,ne)
        CALL Dprimal1(psi3,v3,igrid,nv,ne)
        CALL HodgeH(v1,hv1,igrid,ne)
        CALL HodgeH(v2,hv2,igrid,ne)
        CALL HodgeH(v3,hv3,igrid,ne)
        ! Now loop over edges
        DO ie1 = 1, ne
           ! Velocity field that maximizes v(ie1) is
           ! v = a1*v1 + a2*v2 + a3*v3
           ! with
           a1 = v1(ie1)
           a2 = v2(ie1)
           a3 = v3(ie1)
           ! Demand that lumped matrix agrees with full H for this v
           num = a1*hv1(ie1) + a2*hv2(ie1) + a3*hv3(ie1)
           den = a1*v1(ie1)  + a2*v2(ie1)  + a3*v3(ie1)
           hlump(ie1,igrid) = num/den
        ENDDO
        DEALLOCATE(psi1,psi2,psi3,v1,v2,v3,hv1,hv2,hv3)

     ELSE

        PRINT *,'Option ihlump = ',ihlump,' not available in subroutine buildhlump'
        STOP

     ENDIF

  ENDDO

  ! ----------------------------------------------------------------

END SUBROUTINE buildhlump

! ================================================================

SUBROUTINE buildxminv

  ! Determine stencil and coefficients for a local approximation
  ! to the inverse of M.
  !
  ! Two options are coded. The first is simply the inverse of
  ! the diagonal mass lumped approximation to M. The second is
  ! based on a single underrelaxed Jacobi iteration to the
  ! inverse of M.

  USE grid
  IMPLICIT NONE

  LOGICAL :: llump
  INTEGER :: igrid, ie0, ix, ie1
  REAL*8 :: temp, diag
  REAL*8 :: relax

  ! ----------------------------------------------------------------

  ! Underrelaxation coefficient depends on grid
  ! 0.9 and 1.4 when using mlump in Jacobi
  ! Check for consistency with massMinv
  IF (nefmx == 4) THEN
     ! Use sparse approximate inverse on cubed sphere
     llump = .FALSE.
     relax = 0.9
  ELSE
     ! Use diagonal approximate inverse on hexagonal-icosahedral grid
     llump = .TRUE.
     relax = 1.4
  ENDIF


  IF (llump) THEN

     ! Diagonal operator: stencil size is 1
     nxmisx = 1
     ALLOCATE(nxminvsten(nedgex,ngrids), xminvsten(nxmisx,nedgex,ngrids), &
          xminv(nxmisx,nedgex,ngrids))
     DO igrid = 1, ngrids
        DO ie0 = 1, nedge(igrid)
           ! Stencil for edge ie0 is ie0 itself, and coeff is
           ! inverse of the diagonal term of the lumped matrix
           nxminvsten(ie0,igrid) = 1
           xminvsten(1,ie0,igrid) = ie0
           xminv(1,ie0,igrid) = 1.0d0/mlump(ie0,igrid)
        ENDDO
     ENDDO

  ELSE

     ! Stencil is the same as for M itself
     nxmisx = nmsmx
     ALLOCATE(nxminvsten(nedgex,ngrids), xminvsten(nxmisx,nedgex,ngrids), &
          xminv(nxmisx,nedgex,ngrids))
     DO igrid = 1, ngrids
        DO ie0 = 1, nedge(igrid)
           ! Stencil for edge ie0 is the same as the stencil for M
           nxminvsten(ie0,igrid) = nmsten(ie0,igrid)
           DO ix = 1, nmsten(ie0,igrid)
              ie1 = msten(ix,ie0,igrid)
              xminvsten(ix,ie0,igrid) = ie1
              IF (ie1 == ie0) THEN
                 diag = 1.0d0 + relax
              ELSE
                 diag = 0.0d0
              ENDIF
              temp = mmass(ix,ie0,igrid)/mlump(ie1,igrid)
              xminv(ix,ie0,igrid) = (diag - relax*temp)/mlump(ie0,igrid)
           ENDDO
        ENDDO
     ENDDO

  ENDIF


  ! ----------------------------------------------------------------

END SUBROUTINE buildxminv

! ================================================================

SUBROUTINE massL(f1,f2,igrid,nf)

  ! Apply the mass matrix L to field f1 to obtain the result f2

  USE grid
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: igrid, nf
  REAL*8, INTENT(IN) :: f1(nf)
  REAL*8, INTENT(OUT) :: f2(nf)
  INTEGER :: if1, if2, ix
  REAL*8 :: temp

  ! ----------------------------------------------------------------

  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(temp, if2) SHARED(nlsten, lsten, f1, f2, lmass, igrid, nf)
  DO if1 = 1, nf
     temp = 0.0d0
     DO ix = 1, nlsten(if1,igrid)
        if2 = lsten(ix,if1,igrid)
        temp = temp + f1(if2)*lmass(ix,if1,igrid)
     ENDDO
     f2(if1) = temp
  ENDDO
  !$OMP END PARALLEL DO


  ! ----------------------------------------------------------------

END SUBROUTINE massL

! ================================================================

SUBROUTINE massM(f1,f2,igrid,ne)

  ! Apply the mass matrix M to field f1 to obtain field f2

  USE grid
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: igrid, ne
  REAL*8, INTENT(IN) :: f1(ne)
  REAL*8, INTENT(OUT) :: f2(ne)
  INTEGER :: ie1, ie2, ix
  REAL*8 :: temp

  ! ----------------------------------------------------------------


  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(temp, ie2, ix) SHARED(msten, f1, f2, nmsten, igrid, mmass, ne)
  DO ie1 = 1, ne
     temp = 0.0d0
     DO ix = 1, nmsten(ie1,igrid)
        ie2 = msten(ix,ie1,igrid)
        temp = temp + f1(ie2)*mmass(ix,ie1,igrid)
     ENDDO
     f2(ie1) = temp
  ENDDO
  !$OMP END PARALLEL DO

  ! ----------------------------------------------------------------

END SUBROUTINE massM

! ================================================================

SUBROUTINE approxMinv(f1,f2,igrid,ne)

  ! Apply an approximate inverse of the mass matrix M
  ! to field f1 to obtain field f2

  USE grid
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: igrid, ne
  REAL*8, INTENT(IN) :: f1(ne)
  REAL*8, INTENT(OUT) :: f2(ne)
  INTEGER :: ie1, ie2, ix
  REAL*8 :: temp

  ! ----------------------------------------------------------------

  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(temp, ie2) SHARED(ne, nxminvsten, xminvsten, f1, xminv, f2, igrid)
  DO ie1 = 1, ne
     temp = 0.0d0
     DO ix = 1, nxminvsten(ie1,igrid)
        ie2 = xminvsten(ix,ie1,igrid)
        temp = temp + f1(ie2)*xminv(ix,ie1,igrid)
     ENDDO
     f2(ie1) = temp
  ENDDO
  !$OMP END PARALLEL DO

  ! ----------------------------------------------------------------

END SUBROUTINE approxMinv

! ================================================================

SUBROUTINE HodgeJ(f1,f2,igrid,nv)

  ! Apply the Hodge star J operator that converts dual face
  ! integrals f1 to vertex values f2 on grid igrid.

  USE grid
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: igrid, nv
  REAL*8, INTENT(IN) :: f1(nv)
  REAL*8, INTENT(OUT) :: f2(nv)
  INTEGER :: iv1, iv2, ix
  REAL*8 :: temp

  ! ----------------------------------------------------------------

  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(iv1, iv2, temp) SHARED(nv, njsten, jsten, jstar, f1, f2, igrid)
  DO iv1 = 1, nv
     temp = 0.0d0
     DO ix = 1, njsten(iv1,igrid)
        iv2 = jsten(ix,iv1,igrid)
        temp = temp + f1(iv2)*jstar(ix,iv1,igrid)
     ENDDO
     f2(iv1) = temp
  ENDDO
  !$OMP END PARALLEL DO

  ! ----------------------------------------------------------------

END SUBROUTINE HodgeJ

! ================================================================

SUBROUTINE HodgeH(f1,f2,igrid,ne)

  ! Apply the Hodge star H operator that converts dual edge
  ! integrals (circulations) f1 to primal edge integrals (fluxes) f2
  ! on grid igrid.

  USE grid
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: igrid, ne
  REAL*8, INTENT(IN) :: f1(ne)
  REAL*8, INTENT(OUT) :: f2(ne)
  INTEGER :: ie1, ie2, ix
  REAL*8 :: temp

  ! ----------------------------------------------------------------

  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(ie1, ie2, ix, temp) SHARED(ne, nhsten, hsten, hstar, f1, f2, igrid)
  DO ie1 = 1, ne
     temp = 0.0d0
     DO ix = 1, nhsten(ie1,igrid)
        ie2 = hsten(ix,ie1,igrid)
        temp = temp + f1(ie2)*hstar(ix,ie1,igrid)
     ENDDO
     f2(ie1) = temp
  ENDDO
  !$OMP END PARALLEL DO

  ! ----------------------------------------------------------------

END SUBROUTINE HodgeH

! ================================================================

SUBROUTINE massLinv(f1,f2,igrid,nf,niter)

  ! Apply the inverse of the mass matrix L to field f1 to obtain
  ! the result f2.

  ! If niter >= 0 then f2 is assumed to be an initial estimate
  ! for the solution, and niter further iterations are taken.
  ! If niter < 0 then f2 is initialized and -niter
  ! iterations are taken.
  ! If L is diagonal then there is no need to iterate.

  USE grid
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: niter
  INTEGER, INTENT(IN) :: igrid, nf
  REAL*8, INTENT(IN) :: f1(nf)
  REAL*8, INTENT(INOUT) :: f2(nf)
  INTEGER :: if1, iter, miter
  REAL*8, ALLOCATABLE :: temp(:)

  ALLOCATE(temp(nf))
  ! ----------------------------------------------------------------

  miter = ABS(niter)

  IF (niter < 0 .OR. nlsmx == 1) THEN
     ! First guess based on diagonal L
     DO if1 = 1, nf
        f2(if1) = f1(if1)/lmass(1,if1,igrid)
     ENDDO
  ENDIF

  IF (nlsmx > 1) THEN
     ! L is not diagonal, so use Jacobi iteration to invert
     DO iter = 1, miter
        CALL massL(f2,temp,igrid,nf)
        !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(nf, f1, f2, temp, lmass, igrid)
        DO if1 = 1, nf
           f2(if1) = f2(if1) + (f1(if1) - temp(if1))/lmass(1,if1,igrid)
        ENDDO
        !$OMP END PARALLEL DO
     ENDDO
  ENDIF

  DEALLOCATE(temp)

  ! ----------------------------------------------------------------

END SUBROUTINE massLinv

! ================================================================

SUBROUTINE massMinv(f1,f2,igrid,ne,niter)

  ! Apply the inverse of the mass matrix M to the field f1
  ! to obtain the result f2

  ! If niter >= 0 then f2 is assumed to be an initial estimate
  ! for the solution, and niter further iterations are taken.
  ! If niter < 0 then f2 is initialized and -niter
  ! iterations are taken.
  ! If M is diagonal then there is no need to iterate.

  USE grid
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: niter
  INTEGER, INTENT(IN) :: igrid, ne
  REAL*8, INTENT(IN) :: f1(ne)
  REAL*8, INTENT(INOUT) :: f2(ne)
  INTEGER :: ie1, iter, miter
  REAL*8, ALLOCATABLE :: temp(:)
  REAL*8 :: relax

  ALLOCATE(temp(ne))

  ! ----------------------------------------------------------------

  ! Underrelaxation coefficient depends on grid
  IF (nefmx == 4) THEN
     relax = 0.9
  ELSE
     relax = 1.4
  ENDIF

  miter = ABS(niter)

  IF (niter < 0 .OR. nmsmx == 1) THEN
     ! First guess based on lumped M
     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(ie1) SHARED(ne, f2, f1, mlump, igrid)
     DO ie1 = 1, ne
        f2(ie1) = f1(ie1)/mlump(ie1,igrid)
     ENDDO
     !$OMP END PARALLEL DO
  ENDIF

  IF (nmsmx > 1) THEN
     ! M is not diagonal, so use Jacobi iteration to invert
     DO iter = 1, miter
        CALL massM(f2,temp,igrid,ne)

        !$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(temp, relax, f1)
        temp = relax*(f1 - temp)
        !$OMP END PARALLEL WORKSHARE

        !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(ie1) SHARED(ne, f2, f1, temp, mlump, igrid)
        DO ie1 = 1, ne
           f2(ie1) = f2(ie1) + temp(ie1)/mlump(ie1,igrid)
        ENDDO
        !$OMP END PARALLEL DO
     ENDDO
  ENDIF

  DEALLOCATE(temp)

  ! ----------------------------------------------------------------

END SUBROUTINE massMinv

! ================================================================

SUBROUTINE HodgeJinv(f1,f2,igrid,nv,niter)

  ! Apply the inverse Hodge star operator J^{-1} that maps from
  ! E_p to V_d on grid igrid.
  !
  ! If niter >= 0 then f2 is assumed to be an initial estimate
  ! for the solution, and niter further iterations are taken.
  ! If niter < 0 then f2 is initialized and -niter
  ! iterations are taken.
  ! If J is diagonal then there is no need to iterate.

  USE grid
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: niter
  INTEGER, INTENT(IN) :: igrid, nv
  REAL*8, INTENT(IN) :: f1(nv)
  REAL*8, INTENT(INOUT) :: f2(nv)
  INTEGER :: iv1, iter, miter
  REAL*8, ALLOCATABLE :: temp(:)
  REAL*8 :: relax = 1.4d0 ! relax = 1.4 is good for ijlump = 3 on hex and cube grids 

  ALLOCATE(temp(nv))

  ! ----------------------------------------------------------------

  miter = ABS(niter)

  IF (niter < 0 .OR. njsmx == 1) THEN
     ! First guess based on lumped J
     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(iv1) SHARED(nv, f2, f1, jlump, igrid)
     DO iv1 = 1, nv
        f2(iv1) = f1(iv1)/jlump(iv1,igrid)
     ENDDO
     !$OMP END PARALLEL DO
  ENDIF

  IF (njsmx > 1) THEN
     ! J is not diagonal, so use Jacobi iteration to invert
     DO iter = 1, miter
        CALL HodgeJ(f2,temp,igrid,nv)
        !$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(temp, relax, f1)
        temp = relax*(f1 - temp)
        !$OMP END PARALLEL WORKSHARE

        !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(iv1) SHARED(nv, f2, temp, jlump, igrid)
        DO iv1 = 1, nv
           f2(iv1) = f2(iv1) + temp(iv1)/jlump(iv1,igrid)
        ENDDO
        !$OMP END PARALLEL DO
     ENDDO
  ENDIF

  DEALLOCATE(temp)
  ! ----------------------------------------------------------------

END SUBROUTINE HodgeJinv

! ================================================================

SUBROUTINE HodgeHinv(f1,f2,igrid,ne,niter)

  ! Apply the inverse Hodge star operator H^{-1} that maps from
  ! S_p to S_d on grid igrid.
  !
  ! If niter >= 0 then f2 is assumed to be an initial estimate
  ! for the solution, and niter further iterations are taken.
  ! If niter < 0 then f2 is initialized and -niter
  ! iterations are taken.
  ! If H is diagonal then there is no need to iterate.

  USE grid
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: niter
  INTEGER, INTENT(IN) :: igrid, ne
  REAL*8, INTENT(IN) :: f1(ne)
  REAL*8, INTENT(INOUT) :: f2(ne)
  INTEGER :: ie1, iter, miter
  REAL*8, ALLOCATABLE :: temp(:)
  REAL*8 :: relax = 1.4d0 ! relax = 1.4 is good for ihlump = 3 on hex and cube grids 

  ALLOCATE(temp(ne))

  ! ----------------------------------------------------------------

  miter = ABS(niter)

  IF (niter < 0 .OR. nhsmx == 1) THEN
     ! First guess based on diagonal H

     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(ie1) SHARED(f1, f2, igrid, hlump, ne)
     DO ie1 = 1, ne
        f2(ie1) = f1(ie1)/hlump(ie1,igrid)
     ENDDO
     !$OMP END PARALLEL DO
  ENDIF

  IF (nhsmx > 1) THEN
     ! H is not diagonal, so use Jacobi iteration to invert
     DO iter = 1, miter
        CALL HodgeH(f2,temp,igrid,ne)

        !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(ie1) SHARED(f2, ne, relax, f1, temp, hlump, igrid)
        DO ie1 = 1, ne
           f2(ie1) = f2(ie1) + relax*(f1(ie1) - temp(ie1))/hlump(ie1,igrid)
        ENDDO
        !$OMP END PARALLEL DO
     ENDDO
  ENDIF

  ! ----------------------------------------------------------------
  DEALLOCATE(temp)

END SUBROUTINE HodgeHinv

! ================================================================

SUBROUTINE operW_original(f1,f2,igrid,ne)

  ! Apply the W operator:
  ! given fluxes f1 across primal edges, construct
  ! the rotated fluxes across dual edges f2, on grid igrid.

  ! This is the original formulation, building W from R
  ! a la TRiSK. It probably requires an MPI reduce operation
  ! so is likely to be inefficient.

  USE grid
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: igrid, ne
  REAL*8, INTENT(IN) :: f1(ne)
  REAL*8, INTENT(OUT) :: f2(ne)
  INTEGER :: if1, ie1, ie2, ix, ix1, ix2, ixv, ne1
  REAL*8 :: ss, w

  ! ----------------------------------------------------------------

  ! Initialize to zero
  f2 = 0.0d0

  ! Loop over faces
  DO if1 = 1, nface(igrid)
     ne1 = neoff(if1,igrid)
     ! For each edge of this face
     DO ix1 = 1, ne1
        ss = -0.5
        ie1 = eoff(ix1,if1,igrid)
        ! Find the contribution to f2 from every other
        ! edge of this face
        DO ix = 0, ne1 - 2
           ixv = MOD(ix1 + ix - 1,ne1) + 1
           ix2 = MOD(ix1 + ix,ne1) + 1
           ie2 = eoff(ix2,if1,igrid)
           ss = ss + rcoeff(if1,ixv,igrid)
           w = -ss*eoffin(ix1,if1,igrid)*eoffin(ix2,if1,igrid)
           f2(ie1) = f2(ie1) + w*f1(ie2)
        ENDDO
     ENDDO
  ENDDO

  ! ----------------------------------------------------------------

END SUBROUTINE operW_original

! ================================================================

SUBROUTINE operR_original(f1,f2,igrid,nf,nv)

  ! Apply the R operator:
  ! map from V_p to E_p

  ! This is the original formulation. It loops over `source'
  ! entities rather than target entities and so will require
  ! an MPI reduce.

  USE grid
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: igrid, nf, nv
  REAL*8, INTENT(IN) :: f1(nf)
  REAL*8, INTENT(OUT) :: f2(nv)
  INTEGER :: if1, iv1, ix1, ne1

  ! ----------------------------------------------------------------

  ! Initialize to zero
  f2 = 0.0d0

  ! Loop over faces
  DO if1 = 1, nface(igrid)
     ne1 = neoff(if1,igrid)
     ! Share out this face's contributions to its surrounding vertices
     DO ix1 = 1, ne1
        iv1 = rsten(if1,ix1,igrid)
        f2(iv1) = f2(iv1) + f1(if1)*rcoeff(if1,ix1,igrid)
     ENDDO
  ENDDO

  ! ----------------------------------------------------------------

END SUBROUTINE operR_original

! ================================================================

SUBROUTINE operW(f1,f2,igrid,ne)

  ! Apply the W operator:
  ! given edge integrals of normal components f1 on primal edges,
  ! construct edge integrals of normal components of perpendicular
  ! field f2, on grid igrid.

  ! This formulation uses pre-build stencil and coefficients to
  ! avoid the need for MPI reduce. It is mathematically equivalent
  ! to the original formulation.

  USE grid
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: igrid, ne
  REAL*8, INTENT(IN) :: f1(ne)
  REAL*8, INTENT(OUT) :: f2(ne)
  INTEGER :: ie0, ne1, ix1, ie1
  REAL*8 :: temp

  ! ----------------------------------------------------------------

  ! Loop over vertices
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(temp, ie1, ne1) SHARED(nwsten, igrid, wsten, f1, f2, wcoeff, nedge)
  DO ie0 = 1, nedge(igrid)
     ne1 = nwsten(ie0,igrid)
     ! Collect contributions from stencil
     temp = 0.0d0
     DO ix1 = 1, ne1
        ie1 = wsten(ix1,ie0,igrid)
        temp = temp + f1(ie1)*wcoeff(ix1,ie0,igrid)
     ENDDO
     f2(ie0) = temp
  ENDDO
  !$OMP END PARALLEL DO

  ! ----------------------------------------------------------------

END SUBROUTINE operW

! ================================================================

SUBROUTINE operR(f1,f2,igrid,nf,nv)

  ! Apply the R operator:
  ! given face integrals f1 on primal faces, map to dual cell
  ! integrals f2, on grid igrid.

  ! This formulation stores the coefficients in the transpose of
  ! the original formulation to avoid an MPI reduce. It is
  ! mathematically equivalent to the original formulation.

  USE grid
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: igrid, nf, nv
  REAL*8, INTENT(IN) :: f1(nf)
  REAL*8, INTENT(OUT) :: f2(nv)
  INTEGER :: iv0, if1, ix1, ne1
  REAL*8 :: temp

  ! ----------------------------------------------------------------

  ! Loop over vertices
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(temp, ne1, if1) SHARED(nrxsten, rxsten, rxcoeff, f1, f2, igrid, nvert)
  DO iv0 = 1, nvert(igrid)
     ne1 = nrxsten(iv0,igrid)
     ! Collect contributions from surrounding faces
     temp = 0.0d0
     DO ix1 = 1, ne1
        if1 = rxsten(ix1,iv0,igrid)
        temp = temp + f1(if1)*rxcoeff(ix1,iv0,igrid)
     ENDDO
     f2(iv0) = temp
  ENDDO
  !$OMP END PARALLEL DO

  ! ----------------------------------------------------------------

END SUBROUTINE operR

! ================================================================

SUBROUTINE operT(f1,f2,igrid,ne,nf)

  ! Apply the T operator:
  ! compute cell integrals of 2 x kinetic energy from edge integrals
  ! of normal fluxes

  USE grid
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: igrid, ne, nf
  REAL*8, INTENT(IN) :: f1(ne)
  REAL*8, INTENT(OUT) :: f2(nf)
  INTEGER :: if1, ix1, ix2, ne1, ie1, ie2
  REAL*8 :: temp

  ! ----------------------------------------------------------------

  ! Loop over faces
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(ne1, temp, ie1, ix1, ix2, ie2) &
  !$OMP    SHARED(ntsten, tsten, tcoeff, f1, igrid, f2, nface)
  DO if1 = 1, nface(igrid)
     ne1 = ntsten(if1,igrid)
     temp = 0.0d0
     ! Loop over all pairs of edges of this cell
     DO ix1 = 1, ne1
        ie1 = tsten(ix1,if1,igrid)
        DO ix2 = 1, ne1
           ie2 = tsten(ix2,if1,igrid)
           temp = temp + tcoeff(ix2,ix1,if1,igrid)*f1(ie1)*f1(ie2)
        ENDDO
     ENDDO
     f2(if1) = temp
  ENDDO
  !$OMP END PARALLEL DO

  ! ----------------------------------------------------------------

END SUBROUTINE operT

! ================================================================

SUBROUTINE graddiv(f1,f2,igrid,ne,nf)

  ! Calculate the grad-div component of the vector Laplacian
  ! operator applied to fluxes

  USE grid
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: igrid, ne, nf
  REAL*8, INTENT(IN) :: f1(ne)
  REAL*8, INTENT(OUT) :: f2(ne)
  REAL*8, ALLOCATABLE :: temp1(:), temp2(:), temp3(:)

  ALLOCATE(temp1(nf), temp2(nf), temp3(ne))

  ! ----------------------------------------------------------------

  ! Build up from simpler operators
  CALL Dprimal2(f1,temp1,igrid,ne,nf)
  CALL massL(temp1,temp2,igrid,nf)
  CALL Ddual1(temp2,temp3,igrid,nf,ne)
  CALL massMinv(temp3,f2,igrid,ne,-20)

  ! ----------------------------------------------------------------

  DEALLOCATE(temp1, temp2, temp3)

END SUBROUTINE graddiv

! ================================================================

!SUBROUTINE curlcurl(f1,f2,igrid,ne,nv)
!
!! Calculate the curlcurl component of the vector Laplacian
!! operator applied to circulations
!
!USE grid
!IMPLICIT NONE
!
!INTEGER, INTENT(IN) :: igrid, ne, nv
!REAL*8, INTENT(IN) :: f1(ne)
!REAL*8, INTENT(OUT) :: f2(ne)
!REAL*8 :: temp1(nv), temp2(nv), temp3(ne)
!
!! ----------------------------------------------------------------
!
!! Build up from simpler operators
!CALL Ddual2(f1,temp1,igrid,ne,nv)
!CALL HodgeJ(temp1,temp2,igrid,nv)
!CALL Dprimal1(temp2,temp3,igrid,nv,ne)
!CALL HodgeHinv(temp3,f2,igrid,ne)
!f2 = -f2
!
!print *,'We need the inverse of the E_p space mass matrix'
!print *,'to define curl of curl'
!print *,'STOP in subroutine curlcurl'
!stop
!
!! ----------------------------------------------------------------
!
!
!END SUBROUTINE curlcurl
!
! ================================================================

SUBROUTINE cell2edge(f1,f2,igrid,nf,ne)

  ! Simple interpolation from cell centre point values to
  ! edges

  USE grid

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: igrid, nf, ne
  REAL*8, INTENT(IN) :: f1(nf)
  REAL*8, INTENT(OUT) :: f2(ne)
  INTEGER :: ie1, if1, if2


  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(if1, if2) SHARED(fnxte, igrid, f1, f2, ne)
  DO ie1 = 1, ne
     if1 = fnxte(1,ie1,igrid)
     if2 = fnxte(2,ie1,igrid)
     f2(ie1) = 0.5d0*(f1(if1) + f1(if2))
  ENDDO
  !$OMP END PARALLEL DO


END SUBROUTINE cell2edge

! ================================================================

SUBROUTINE restrict(f1,nf1,f2,nf2,igrid)

  ! To perform the restriction operation needed for a multigrid solver.
  ! Restrict field f1 from grid igrid + 1 to grid igrid and put the
  ! result in field f2. f1 and f2 are assumed to be area integrals
  ! (discrete 2-forms).

  USE grid

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nf1, nf2, igrid
  REAL*8, INTENT(IN) :: f1(nf1)
  REAL*8, INTENT(OUT) :: f2(nf2)
  INTEGER :: if1, if2, ix
  REAL*8 :: wgt

  ! Safety check
  IF (nf2 .NE. nface(igrid)) THEN
     PRINT *,'Wrong size array in subroutine restrict'
     STOP
  ENDIF

  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(ix, if1, wgt) SHARED(nres, ressten, reswgt, f2, f1, igrid, nf2)
  DO if2 = 1, nf2
     f2(if2) = 0.0d0
     DO ix = 1, nres(if2,igrid)
        if1 = ressten(ix,if2,igrid)
        wgt = reswgt(ix,if2,igrid)
        f2(if2) = f2(if2) + wgt*f1(if1)
     ENDDO
  ENDDO
  !$OMP END PARALLEL DO

  ! ----------------------------------------------------------------

END SUBROUTINE restrict

! ================================================================

SUBROUTINE prolong(f2,nf2,f1,nf1,igrid)

  ! To perform the prolongation operation needed for a multigrid solver.
  ! Prolong field f2 from grid igrid to grid igrid + 1 and put the
  ! result in field f1. f1 and f2 are assumed to be area integrals
  ! (discrete 2-forms); so f2 must be converted to point values
  ! (zero-form) first and f1 must be converted back to an area-integral
  ! (2-form) at the end. The prolongation operator is the adjoint of
  ! the restriction operator, so uses the same stencil and weights.

  ! *** Currently this operator uses a loop over source entities,
  ! implying a need for an MPI reduce. The stencil and weights could
  ! be stored in transpose form to allow a loop over target entities
  ! and hence avoid an MPI reduce ***

  USE grid
  USE state

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nf1, nf2, igrid
  REAL*8, INTENT(IN) :: f2(nf2)
  REAL*8, INTENT(OUT) :: f1(nf1)
  INTEGER :: if1, if2, ix, igridp
  INTEGER :: if3, ixx
  INTEGER :: igridr
  REAL*8 :: wgt
  !REAL*8 :: f2if2, temp1(nf1)

  ! Safety check
  IF (nf2 .NE. nface(igrid)) THEN
     PRINT *,'Wrong size array in subroutine prolong'
     STOP
  ENDIF

  IF (.TRUE.) THEN
     igridp = igrid+1
     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(ix, if2, wgt) &
     !$OMP    SHARED(nf1, f1, f2, ninj, injsten, injwgt, farea, igrid, igridp)
     DO if1 = 1, nf1
        f1(if1) = 0.0d0

        DO ix = 1, ninj(if1,igridp)
           if2 = injsten(ix,if1,igridp)
           wgt = injwgt(ix,if1,igridp)

           f1(if1) = f1(if1) + wgt*(f2(if2)/farea(if2, igrid))
        END DO

        f1(if1) = f1(if1)*farea(if1,igridp)
     END DO
     !$OMP END PARALLEL DO
  END IF


  !IF (.FALSE.) THEN
  !    igridp = igrid + 1
  !    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(temp1)
  !    temp1 = 0.0d0
  !    !$OMP END PARALLEL WORKSHARE

  !    DO if2 = 1, nf2
  !      f2if2 = f2(if2)/farea(if2,igrid)

  !      DO ix = 1, nres(if2,igrid)
  !        if1 = ressten(if2,ix,igrid)
  !        wgt = reswgt(if2,ix,igrid)

  !        temp1(if1) = temp1(if1) + wgt*f2if2
  !      ENDDO
  !    ENDDO

  !    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(f1, nf1, temp1, farea, igridp)
  !    f1(1:nf1) = temp1(1:nf1)*farea(1:nf1,igridp)
  !    !$OMP END PARALLEL WORKSHARE
  !END IF

  ! ----------------------------------------------------------------

END SUBROUTINE prolong

! ================================================================

SUBROUTINE helmholtz(f,hf,igrid,nf,ne)

  ! To apply the Helmholtz operator to the input field f,
  ! on grid igrid, the result appearing in the output field hf.
  ! Note f and hf are area integrals (2-forms).

  USE helmcoeff

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: igrid, nf, ne
  REAL*8, INTENT(IN) :: f(nf)
  REAL*8, INTENT(OUT) :: hf(nf)
  REAL*8, ALLOCATABLE :: temp1(:), temp2(:), temp3(:)

  ALLOCATE(temp1(nf), temp2(ne), temp3(ne))

  CALL massL(f,temp1,igrid,nf)
  CALL Ddual1(temp1,temp2,igrid,nf,ne)
  CALL approxMinv(temp2,temp3,igrid,ne)
  !$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(temp2, temp3, nusq, ne, igrid)
  temp2 = temp3*nusq(1:ne,igrid)
  !$OMP END PARALLEL WORKSHARE
  CALL Dprimal2(temp2,hf,igrid,ne,nf)

  !$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(hf, f)
  hf = hf - f
  !$OMP END PARALLEL WORKSHARE

  DEALLOCATE(temp1, temp2, temp3)

END SUBROUTINE helmholtz

! ================================================================

SUBROUTINE residual(f,rhs,res,igrid,nf,ne)

  ! Compute the residual res in the helmholtz equation on grid igrid
  ! when f is the input field and rhs is the right hand side. Note that
  ! f, rhs and res are area integrals (2-forms).

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: igrid, nf, ne
  REAL*8, INTENT(IN) :: f(nf), rhs(nf)
  REAL*8, INTENT(OUT) :: res(nf)


  CALL helmholtz(f,res,igrid,nf,ne)
  !$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(res, rhs)
  res = rhs - res
  !$OMP END PARALLEL WORKSHARE


END SUBROUTINE residual

! ================================================================

SUBROUTINE relax(f,rhs,igrid,nf,ne,niter)

  ! To carry out niter Jacobi relaxation iterations for the multigrid
  ! solver on grid igrid

  USE helmcoeff

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: igrid, nf, ne, niter
  REAL*8, INTENT(IN) :: rhs(nf)
  REAL*8, INTENT(INOUT) :: f(nf)
  REAL*8, ALLOCATABLE :: res(:), finc(:)
  REAL*8 :: u
  INTEGER :: iter

  ALLOCATE(res(nf), finc(nf))

  u = underrel(igrid)
  DO iter = 1, niter
     CALL residual(f,rhs,res,igrid,nf,ne)
     !$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(finc, res, helmdiag, nf, igrid, f, u)
     finc = res/helmdiag(1:nf,igrid)
     f = f + u*finc
     !$OMP END PARALLEL WORKSHARE
  ENDDO

  DEALLOCATE(res, finc)

END SUBROUTINE relax

! ================================================================

SUBROUTINE mgsolve(phi,rr,ng)

  ! Multigrid solver for elliptic equation
  !
  ! Dprimal2 nusq H Ddual1 I phi - phi = RHS
  !
  ! using full multigrid algorithm.
  ! Coefficients are contained in module helmholtz.

  USE grid
  USE state

  IMPLICIT NONE

  ! Numbers of iterations on coarsest grid and other grids
  INTEGER, PARAMETER :: niterc = 10, niter = 2, npass = 1

  INTEGER, INTENT(IN) :: ng
  REAL*8, INTENT(IN) :: rr(nfacex)
  REAL*8, INTENT(OUT) :: phi(nfacex)

  INTEGER :: ipass, nf1, nf2, ne1, ne2, igrid, igridp, jgrid, jgridp, iter
  REAL*8, ALLOCATABLE :: ff(:,:), rf(:,:)
  REAL*8, ALLOCATABLE :: temp1(:)

  ALLOCATE(temp1(nfacex))
  ! ------------------------------------------------------------------------

  ! Allocate space on all grids
  ALLOCATE(ff(nfacex,ngrids),rf(nfacex,ngrids))

  ! One pass should be enough. Warn user if npass is set to
  ! some other value for testing
  IF (npass > 1) PRINT *,'mgsolve: npass = ',npass

  ! ------------------------------------------------------------------------

  ! Initialize solution to zero
  !$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(phi)
  phi = 0.0d0
  !$OMP END PARALLEL WORKSHARE

  ! For diagnostics
  !nf1 = nface(ngrids)
  !ne1 = nedge(ngrids)
  !CALL residual(phi,rr,temp1,ngrids,nf1,ne1)
  !print *,'Pass ',0,'  RMS residual = ',SQRT(SUM(temp1*temp1)/nf1)

  DO ipass = 1, npass

     ! Initialize rhs as residual using latest estimate
     IF (ipass == 1) THEN
        ! No need to do the calculation
        !$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(rf,rr,ngrids)
        rf(:,ngrids) = rr(:)
        !$OMP END PARALLEL WORKSHARE
     ELSE
        nf1 = nface(ngrids)
        ne1 = nedge(ngrids)
        CALL residual(phi,rr,rf(1,ngrids),ngrids,nf1,ne1)
     ENDIF

     ! Initialize correction to solution on all grids to zero
     ff = 0.0d0

     ! Restrict right hand side to each grid in the hierarchy
     DO igrid = ngrids-1, ngrids-ng+1, -1
        igridp = igrid + 1
        nf1 = nface(igridp)
        nf2 = nface(igrid)
        CALL restrict(rf(1,igridp),nf1,rf(1,igrid),nf2,igrid)
     ENDDO

     ! Iterate to convergence on coarsest grid
     igrid = ngrids-ng+1
     nf1 = nface(igrid)
     ne1 = nedge(igrid)
     !$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(ff, nf1, igrid)
     ff(1:nf1,igrid) = 0.0d0
     !$OMP END PARALLEL WORKSHARE
     CALL relax(ff(:,igrid),rf(:,igrid),igrid,nf1,ne1,niterc)

     ! Sequence of growing V-cycles
     DO igridp = ngrids-ng+2, ngrids

        igrid = igridp - 1
        nf1 = nface(igridp)
        ne1 = nedge(igridp)
        nf2 = nface(igrid)
        ne2 = nedge(igrid)

        ! Prolong solution to grid igridp
        ! and execute one V-cycle starting from grid igridp

        ! Prolong
        CALL prolong(ff(1,igrid),nf2,ff(1,igridp),nf1,igrid)

        ! Descending part of V-cycle
        DO jgrid = igrid, ngrids-ng+1, -1

           jgridp = jgrid + 1
           nf1 = nface(jgridp)
           ne1 = nedge(jgridp)
           nf2 = nface(jgrid)
           ne2 = nedge(jgrid)

           ! Relax on grid jgridp
           CALL relax(ff(1,jgridp),rf(1,jgridp),jgridp,nf1,ne1,niter)

           ! Calculate residual on jgridp
           CALL residual(ff(1,jgridp),rf(1,jgridp),temp1,jgridp,nf1,ne1)

           ! Restrict residual to jgrid
           CALL restrict(temp1,nf1,rf(1,jgrid),nf2,jgrid)

           ! Set correction first guess to zero on grid jgrid-1
           !$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(ff, nf2, jgrid)
           ff(1:nf2,jgrid) = 0.0d0
           !$OMP END PARALLEL WORKSHARE

        ENDDO

        ! Iterate to convergence on coarsest grid
        jgrid = ngrids-ng+1
        nf1 = nface(jgrid)
        ne1 = nedge(jgrid)
        ff(1:nf1,jgrid) = 0.0d0
        CALL relax(ff(1,jgrid),rf(1,jgrid),jgrid,nf1,ne1,niterc)

        ! Ascending part of V-cycle
        DO jgrid = ngrids-ng+1, igrid

           jgridp = jgrid + 1
           igrid = igrid - 1
           nf1 = nface(jgridp)
           ne1 = nedge(jgridp)
           nf2 = nface(jgrid)
           ne2 = nedge(jgrid)

           ! Prolong correction to jgridp
           CALL prolong(ff(1,jgrid),nf2,temp1,nf1,jgrid)

           ! Add correction to solution on jgridp
           !$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(ff, nf1, temp1, jgridp)
           ff(1:nf1,jgridp) = ff(1:nf1,jgridp) + temp1(1:nf1)
           !$OMP END PARALLEL WORKSHARE

           ! Relax on grid jgridp
           CALL relax(ff(1,jgridp),rf(1,jgridp),jgridp,nf1,ne1,niter)

        ENDDO

     ENDDO

     ! Add correction to phi
     !$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(phi, ff, ngrids)
     phi = phi + ff(:,ngrids)
     !$OMP END PARALLEL WORKSHARE


     ! For diagnostics
     ! nf1 = nface(ngrids)
     ! ne1 = nedge(ngrids)
     ! CALL residual(phi,rr,temp1,ngrids,nf1,ne1)
     ! print *,'Pass ',ipass,'  RMS residual = ',SQRT(SUM(temp1*temp1)/nf1)

  ENDDO


  DEALLOCATE(ff,rf)

  DEALLOCATE(temp1)
  ! ----------------------------------------------------------------

END SUBROUTINE mgsolve

! ================================================================

SUBROUTINE readgrid

  ! To allocate array space for the grid information in module grid
  ! and to read the information from file

  USE runtype
  USE constants
  USE grid
  USE channels

  IMPLICIT NONE

  INTEGER :: if0, ie0, iv0, igrid, ix, ixx, if1, if2, iv1, iv2, ie1

  ! ----------------------------------------------------------------

  ! Open file for reading
  OPEN(changrid,FILE=ygridfile,FORM='UNFORMATTED')

  ! First read ngrids
  READ(changrid) ngrids


  ! Allocate nface, nedge, nvert
  ALLOCATE(nface(ngrids), nedge(ngrids), nvert(ngrids))

  ! Read numbers of faces, edges and vertices on each grid
  READ(changrid) nface
  READ(changrid) nedge
  READ(changrid) nvert

  ! Find maximum values in order to allocated subsequent arrays
  nfacex = MAXVAL(nface)
  nedgex = MAXVAL(nedge)
  nvertx = MAXVAL(nvert)

  ! Allocate neoff, neofv
  ALLOCATE(neoff(nfacex,ngrids), neofv(nvertx,ngrids))

  ! Initialize to account for first touch policy
  ! We hope that the compiler is not removing these lines
  DO igrid= 1, ngrids
     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(if0) SHARED(nfacex, neoff, igrid)
     DO if0=1, nfacex
        neoff(if0,igrid) = 0
     END DO
     !$OMP END PARALLEL DO
     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(iv0) SHARED(nvertx, neofv, igrid)
     DO iv0=1, nvertx
        neofv(iv0,igrid) = 0
     END DO
     !$OMP END PARALLEL DO
  END DO

  ! Read the numbers of edges of each face and vertex on each grid
  neoff = 0
  neofv = 0
  READ(changrid) ((neoff(if0,igrid),          &
       if0 = 1, nface(igrid)), &
       igrid = 1, ngrids)
  READ(changrid) ((neofv(iv0,igrid),          &
       iv0 = 1, nvert(igrid)), &
       igrid = 1, ngrids)

  ! Find maximum values in order to allocate subsequent arrays
  nefmx = MAXVAL(neoff)
  nevmx = MAXVAL(neofv)


  ! Allocate connectivity arrays arrays
  ALLOCATE(fnxtf(nfacex,nefmx,ngrids), eoff(nefmx,nfacex,ngrids), &
       voff(nfacex,nefmx,ngrids),  fnxte(2,nedgex,ngrids),    &
       vofe(2,nedgex,ngrids),      fofv(nvertx,nevmx,ngrids), &
       eofv(nevmx,nvertx,ngrids))


  DO igrid = 1, ngrids
     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(if0) SHARED(nfacex, fnxtf, eoff, voff, igrid)
     DO if0=1, nfacex
        fnxtf(if0,:,igrid) = 0
        eoff(:,if0,igrid) = 0
        voff(if0,:,igrid) = 0
     END DO
     !$OMP END PARALLEL DO
     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(ie0) SHARED(nedgex, fnxte, vofe, igrid)
     DO ie0=1, nedgex
        fnxte(:,ie0,igrid) = 0
        vofe(:,ie0,igrid) = 0
     END DO
     !$OMP END PARALLEL DO
     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(iv0) SHARED(nvertx, fofv, eofv, igrid)
     DO iv0=1, nvertx
        fofv(iv0,:,igrid) = 0
        eofv(:,iv0,igrid) = 0
     END DO
     !$OMP END PARALLEL DO
  END DO

  ! Read the connectivity arrays
  READ(changrid) (((fnxtf(if0,ix,igrid),          &
       if0 = 1, nface(igrid)),    &
       ix = 1, nefmx),            &
       igrid = 1, ngrids)
  READ(changrid) (((eoff(ix,if0,igrid),           &
       ix = 1, nefmx),            &
       if0 = 1, nface(igrid)),    &
       igrid = 1, ngrids)
  READ(changrid) (((voff(if0,ix,igrid),           &
       if0 = 1, nface(igrid)),    &
       ix = 1, nefmx),            &
       igrid = 1, ngrids)
  READ(changrid) (((fnxte(ix,ie0,igrid),          &
       ix = 1, 2),                &
       ie0 = 1, nedge(igrid)),    &
       igrid = 1, ngrids)
  READ(changrid) (((vofe(ix,ie0,igrid),           &
       ix = 1, 2),    &
       ie0 = 1, nedge(igrid)),                &
       igrid = 1, ngrids)
  READ(changrid) (((fofv(iv0,ix,igrid),           &
       iv0 = 1, nvert(igrid)),    &
       ix = 1, nevmx),            &
       igrid = 1, ngrids)
  READ(changrid) (((eofv(ix,iv0,igrid),           &
       ix = 1, nevmx),            &
       iv0 = 1, nvert(igrid)),    &
       igrid = 1, ngrids)


  ! Allocate the geometrical information arrays
  ALLOCATE(flong(nfacex,ngrids), flat(nfacex,ngrids),  &
       vlong(nvertx,ngrids), vlat(nvertx,ngrids),  &
       farea(nfacex,ngrids), varea(nvertx,ngrids), &
       ldist(nedgex,ngrids), ddist(nedgex,ngrids), &
       fareamin(ngrids))


  DO igrid = 1, ngrids
     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(if0) SHARED(nfacex, flong, flat, farea, igrid)
     DO if0=1, nfacex
        flong(if0,igrid) = 0
        flat(if0,igrid) = 0
        farea(if0,igrid) = 0
     END DO
     !$OMP END PARALLEL DO
     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(ie0) SHARED(nedgex, ldist, ddist, igrid)
     DO ie0=1, nedgex
        ldist(ie0,igrid) = 0
        ddist(ie0,igrid) = 0
     END DO
     !$OMP END PARALLEL DO
     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(ie0) SHARED(nvertx, vlong, vlat, varea, igrid)
     DO iv0=1, nvertx
        vlong(iv0,igrid) = 0
        vlat(iv0,igrid) = 0
        varea(iv0,igrid) = 0
     END DO
     !$OMP END PARALLEL DO
  END DO

  ! Read the geometrical information arrays
  READ(changrid) ((flong(if0,igrid),              &
       if0 = 1, nface(igrid)),    &
       igrid = 1, ngrids)
  READ(changrid) ((flat(if0,igrid),               &
       if0 = 1, nface(igrid)),    &
       igrid = 1, ngrids)
  READ(changrid) ((vlong(iv0,igrid),              &
       iv0 = 1, nvert(igrid)),    &
       igrid = 1, ngrids)
  READ(changrid) ((vlat(iv0,igrid),               &
       iv0 = 1, nvert(igrid)),    &
       igrid = 1, ngrids)
  READ(changrid) ((farea(if0,igrid),              &
       if0 = 1, nface(igrid)),    &
       igrid = 1, ngrids)
  READ(changrid) ((varea(iv0,igrid),              &
       iv0 = 1, nvert(igrid)),    &
       igrid = 1, ngrids)
  READ(changrid) ((ldist(ie0,igrid),              &
       ie0 = 1, nedge(igrid)),    &
       igrid = 1, ngrids)
  READ(changrid) ((ddist(ie0,igrid),              &
       ie0 = 1, nedge(igrid)),    &
       igrid = 1, ngrids)

  ! Dimensionalize
  farea = farea*rearth*rearth
  varea = varea*rearth*rearth
  ldist = ldist*rearth
  ddist = ddist*rearth

  ! Determine smallest face area on each grid
  DO igrid = 1, ngrids
     fareamin(igrid) = MINVAL(farea(1:nface(igrid),igrid))
  ENDDO


  ! Allocate arrays for size of operator stencils
  ALLOCATE(nlsten(nfacex,ngrids), nmsten(nedgex,ngrids), &
       njsten(nvertx,ngrids), nhsten(nedgex,ngrids), &
       nrsten(nfacex,ngrids), nrxsten(nvertx,ngrids), &
       nwsten(nedgex,ngrids), ntsten(nfacex,ngrids))



  DO igrid = 1, ngrids
     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(if0) SHARED(nfacex, nlsten, nrsten, ntsten, igrid)
     DO if0=1, nfacex
        nlsten(if0,igrid) = 0
        nrsten(if0,igrid) = 0
        ntsten(if0,igrid) = 0
     END DO
     !$OMP END PARALLEL DO
     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(ie0) SHARED(nedgex, nmsten, nhsten, nwsten, igrid)
     DO ie0=1, nedgex
        nmsten(ie0,igrid) = 0
        nhsten(ie0,igrid) = 0
        nwsten(ie0,igrid) = 0
     END DO
     !$OMP END PARALLEL DO
     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(iv0) SHARED(nvertx, njsten, nrxsten, igrid)
     DO iv0=1, nvertx
        njsten(iv0,igrid) = 0
        nrxsten(iv0,igrid) = 0
     END DO
     !$OMP END PARALLEL DO
  END DO


  ! Read the sizes of the operator stencils on each grid
  nlsten = 0
  nmsten = 0
  njsten = 0
  nhsten = 0
  nrsten = 0
  nrxsten = 0
  nwsten = 0
  ntsten = 0
  READ(changrid) ((nlsten(if0,igrid),             &
       if0 = 1, nface(igrid)),    &
       igrid = 1, ngrids)
  READ(changrid) ((nmsten(ie0,igrid),             &
       ie0 = 1, nedge(igrid)),    &
       igrid = 1, ngrids)
  READ(changrid) ((njsten(iv0,igrid),             &
       iv0 = 1, nvert(igrid)),    &
       igrid = 1, ngrids)
  READ(changrid) ((nhsten(ie0,igrid),             &
       ie0 = 1, nedge(igrid)),    &
       igrid = 1, ngrids)
  READ(changrid) ((nrsten(if0,igrid),             &
       if0 = 1, nface(igrid)),    &
       igrid = 1, ngrids)
  READ(changrid) ((nrxsten(iv0,igrid),            &
       iv0 = 1, nvert(igrid)),    &
       igrid = 1, ngrids)
  READ(changrid) ((nwsten(ie0,igrid),             &
       ie0 = 1, nedge(igrid)),    &
       igrid = 1, ngrids)
  READ(changrid) ((ntsten(if0,igrid),             &
       if0 = 1, nface(igrid)),    &
       igrid = 1, ngrids)

  ! Find maximum values in order to allocate subsequent arrays
  nlsmx = MAXVAL(nlsten)
  nmsmx = MAXVAL(nmsten)
  njsmx = MAXVAL(njsten)
  nhsmx = MAXVAL(nhsten)
  nrsmx = MAXVAL(nrsten)
  nrxsmx = MAXVAL(nrxsten)
  nwsmx = MAXVAL(nwsten)
  ntsmx = MAXVAL(ntsten)

  PRINT *,'Maximum stencil sizes:'
  PRINT *,'massL ...  ',nlsmx
  PRINT *,'massM ...  ',nmsmx
  PRINT *,'HodgeJ ... ',njsmx
  PRINT *,'HodgeH ... ',nhsmx
  PRINT *,'operR ...  ',nrxsmx
  PRINT *,'operW ...  ',nwsmx
  PRINT *,'operT ...  ',ntsmx
  PRINT *,' '


  ! Allocate arrays for operator stencils and coefficients
  ALLOCATE(lsten(nlsmx,nfacex,ngrids), msten(nmsmx,nedgex,ngrids), &
       jsten(njsmx,nvertx,ngrids), hsten(nhsmx,nedgex,ngrids), &
       rsten(nfacex,nrsmx,ngrids), rxsten(nrxsmx,nvertx,ngrids), &
       wsten(nwsmx,nedgex,ngrids), tsten(ntsmx,nfacex,ngrids))
  ALLOCATE(lmass(nlsmx,nfacex,ngrids), mmass(nmsmx,nedgex,ngrids), &
       jstar(njsmx,nvertx,ngrids), hstar(nhsmx,nedgex,ngrids), &
       rcoeff(nfacex,nrsmx,ngrids), rxcoeff(nrxsmx,nvertx,ngrids), &
       wcoeff(nwsmx,nedgex,ngrids), tcoeff(ntsmx,ntsmx,nfacex,ngrids), &
       jlump(nvertx,ngrids), mlump(nedgex,ngrids), hlump(nedgex,ngrids))


  DO igrid = 1, ngrids
     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(if0) SHARED(nfacex, lsten, rsten, tsten, lmass, rcoeff, tcoeff, igrid)
     DO if0=1, nfacex
        lsten(:,if0,igrid) = 0
        rsten(if0,:,igrid) = 0
        tsten(:,if0,igrid) = 0

        lmass(:,if0,igrid) = 0
        rcoeff(if0,:,igrid) = 0
        tcoeff(:,:,if0,igrid) = 0
     END DO
     !$OMP END PARALLEL DO
     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(ie0) &
     !$OMP    SHARED(nedgex, msten, hsten, wsten, mmass, hstar, wcoeff, mlump, hlump, igrid)
     DO ie0=1, nedgex
        msten(:,ie0,igrid) = 0
        hsten(:,ie0,igrid) = 0
        wsten(:,ie0,igrid) = 0

        mmass(:,ie0,igrid) = 0
        hstar(:,ie0,igrid) = 0
        wcoeff(:,ie0,igrid) = 0
        mlump(ie0,igrid) = 0
        hlump(ie0,igrid) = 0
     END DO
     !$OMP END PARALLEL DO
     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(iv0) &
     !$OMP    SHARED(nvertx, jsten, rxsten, jstar, rxcoeff, jlump, igrid)
     DO iv0=1, nvertx
        jsten(:,iv0,igrid) = 0
        rxsten(:,iv0,igrid) = 0

        jstar(:,iv0,igrid) = 0
        rxcoeff(:,iv0,igrid) = 0
        jlump(iv0,igrid) = 0
     END DO
     !$OMP END PARALLEL DO
  END DO

  ! Read the operator stencils and coefficients
  READ(changrid) (((lsten(ix,if0,igrid),          &
       ix = 1, nlsmx),    &
       if0 = 1, nface(igrid)),            &
       igrid = 1, ngrids)
  READ(changrid) (((msten(ix,ie0,igrid),          &
       ix = 1, nmsmx),    &
       ie0 = 1, nedge(igrid)),            &
       igrid = 1, ngrids)
  READ(changrid) (((jsten(ix,iv0,igrid),          &
       ix = 1, njsmx),    &
       iv0 = 1, nvert(igrid)),            &
       igrid = 1, ngrids)
  READ(changrid) (((hsten(ix,ie0,igrid),          &
       ix = 1, nhsmx),    &
       ie0 = 1, nedge(igrid)),            &
       igrid = 1, ngrids)
  READ(changrid) (((rsten(if0,ix,igrid),          &
       if0 = 1, nface(igrid)),    &
       ix = 1, nrsmx),            &
       igrid = 1, ngrids)
  READ(changrid) (((rxsten(ix,iv0,igrid),         &
       ix = 1, nrxsmx),    &
       iv0 = 1, nvert(igrid)),           &
       igrid = 1, ngrids)
  READ(changrid) (((wsten(ix,ie0,igrid),          &
       ix = 1, nwsmx),    &
       ie0 = 1, nedge(igrid)),            &
       igrid = 1, ngrids)
  READ(changrid) (((tsten(ix,if0,igrid),          &
       ix = 1, ntsmx),            &
       if0 = 1, nface(igrid)),    &
       igrid = 1, ngrids)
  READ(changrid) (((lmass(ix,if0,igrid),          &
       ix = 1, nlsmx),    &
       if0 = 1, nface(igrid)),            &
       igrid = 1, ngrids)
  READ(changrid) (((mmass(ix,ie0,igrid),          &
       ix = 1, nmsmx),    &
       ie0 = 1, nedge(igrid)),            &
       igrid = 1, ngrids)
  READ(changrid) (((jstar(ix,iv0,igrid),          &
       ix = 1, njsmx),    &
       iv0 = 1, nvert(igrid)),            &
       igrid = 1, ngrids)
  READ(changrid) (((hstar(ix,ie0,igrid),          &
       ix = 1, nhsmx),    &
       ie0 = 1, nedge(igrid)),            &
       igrid = 1, ngrids)
  READ(changrid) (((rcoeff(if0,ix,igrid),         &
       if0 = 1, nface(igrid)),    &
       ix = 1, nrsmx),            &
       igrid = 1, ngrids)
  READ(changrid) (((rxcoeff(ix,iv0,igrid),        &
       ix = 1, nrxsmx),    &
       iv0 = 1, nvert(igrid)),           &
       igrid = 1, ngrids)
  READ(changrid) (((wcoeff(ix,ie0,igrid),         &
       ix = 1, nwsmx),    &
       ie0 = 1, nedge(igrid)),            &
       igrid = 1, ngrids)
  READ(changrid) ((((tcoeff(ixx,ix,if0,igrid),    &
       ixx = 1, ntsmx),    &
       ix = 1, ntsmx),            &
       if0 = 1, nface(igrid)),           &
       igrid = 1, ngrids)

  ! Dimensionalize
  lmass = lmass/(rearth*rearth)

  ! Construct the tables eoffin and eofvin
  ALLOCATE(eoffin(nefmx,nfacex,ngrids), eofvin(nevmx,nvertx,ngrids))



  DO igrid = 1, ngrids

     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(if1, ie1, if2) SHARED(nface, neoff, eoff, eoffin, fnxte, igrid)
     DO if1 = 1, nface(igrid)
        DO ix = 1, neoff(if1,igrid)
           ie1 = eoff(ix,if1,igrid)
           if2 = fnxte(1,ie1,igrid)
           IF (if1 == if2) THEN
              ! This edge points out of face if1
              eoffin(ix,if1,igrid) = -1.0d0
           ELSE
              ! This edge points into face if1
              eoffin(ix,if1,igrid) = 1.0d0
           ENDIF
        ENDDO
     ENDDO
     !$OMP END PARALLEL DO

     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(iv1, ie1, iv2) SHARED(nvert, neofv, eofv, vofe, eofvin, igrid)
     DO iv1 = 1, nvert(igrid)
        DO ix = 1, neofv(iv1,igrid)
           ie1 = eofv(ix,iv1,igrid)
           iv2 = vofe(1,ie1,igrid)
           IF (iv1 == iv2) THEN
              ! This edge points away from vertex iv1
              eofvin(ix,iv1,igrid) = -1.0d0
           ELSE
              ! This edge points towards vertex iv1
              eofvin(ix,iv1,igrid) = 1.0d0
           ENDIF
        ENDDO
     ENDDO
     !$OMP END PARALLEL DO

  ENDDO


  ! Allocate array for size of restriction stencil
  ALLOCATE(nres(nfacex,ngrids-1))


  DO igrid = 1, ngrids-1
     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(if0) SHARED(nfacex, nres, igrid)
     DO if0=1, nfacex
        nres(if0, igrid) = 0
     END DO
     !$OMP END PARALLEL DO
  END DO

  ! Read the size of the restriction stencil on each grid
  nres = 0
  READ(changrid) ((nres(if0,igrid),              &
       if0 = 1, nface(igrid)),    &
       igrid = 1, ngrids-1)

  ! Find maximum value in order to allocate subsequent arrays
  nresmx = MAXVAL(nres)

  ! Allocate arrays for restriction stencils and weights
  ALLOCATE(ressten(nresmx,nfacex,ngrids-1))
  ALLOCATE(reswgt(nresmx,nfacex,ngrids-1))

  DO igrid = 1, ngrids-1
     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(if0) SHARED(nfacex, ressten, reswgt, igrid)
     DO if0=1, nfacex
        ressten(:, if0, igrid) = 0
        reswgt(:, if0, igrid) = 0
     END DO
     !$OMP END PARALLEL DO
  END DO


  ! Read the restriction stencil and weights
  READ(changrid) (((ressten(ix,if0,igrid),       &
       ix = 1, nresmx),    &
       if0 = 1, nface(igrid)),           &
       igrid = 1, ngrids-1)
  READ(changrid) (((reswgt(ix,if0,igrid),        &
       ix = 1, nresmx),    &
       if0 = 1, nface(igrid)),           &
       igrid = 1, ngrids-1)



  ! Allocate array for size of injection stencil
  ALLOCATE(ninj(nfacex,2:ngrids))


  DO igrid = 2, ngrids
     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(if0) SHARED(nfacex, ninj, igrid)
     DO if0=1, nfacex
        ninj(if0, igrid) = 0
     END DO
     !$OMP END PARALLEL DO
  END DO


  ! Read the size of the injection stencil on each grid
  ninj = 0
  READ(changrid) ((ninj(if0,igrid),              &
       if0 = 1, nface(igrid)),    &
       igrid = 2, ngrids)

  ! Find maximum value in order to allocate subsequent arrays
  ninjmx = MAXVAL(ninj)

  ! Allocate arrays for injection stencils and weights
  ALLOCATE(injsten(ninjmx,nfacex,2:ngrids))
  ALLOCATE(injwgt(ninjmx,nfacex,2:ngrids))


  DO igrid = 2, ngrids
     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(if0) SHARED(nfacex, injsten, injwgt, igrid)
     DO if0=1, nfacex
        injsten(:, if0, igrid) = 0
        injwgt(:, if0, igrid) = 0
     END DO
     !$OMP END PARALLEL DO
  END DO

  ! Read the injection stencil and weights
  READ(changrid) (((injsten(ix,if0,igrid),       &
       ix = 1, ninjmx),    &
       if0 = 1, nface(igrid)),           &
       igrid = 2, ngrids)
  READ(changrid) (((injwgt(ix,if0,igrid),        &
       ix = 1, ninjmx),    &
       if0 = 1, nface(igrid)),           &
       igrid = 2, ngrids)


  READ(changrid) SFCIndexAvailable

  IF (SFCIndexAvailable == 1) THEN
     ALLOCATE(fNewFaceId(nfacex, ngrids), fNewFaceIdInverse(nfacex, ngrids), &
          fNewEdgeId(nedgex, ngrids), fNewEdgeIdInverse(nedgex, ngrids), &
          fNewVertId(nvertx, ngrids), fNewVertIdInverse(nvertx, ngrids))

     WRITE(*,*) "Using SFC information for MG"
     READ(changrid) ((fNewFaceId(if0,igrid),          &
          if0 = 1, nface(igrid)),    &
          igrid = 1, ngrids)

     READ(changrid) ((fNewFaceIdInverse(if0,igrid),          &
          if0 = 1, nface(igrid)),    &
          igrid = 1, ngrids)

     READ(changrid) ((fNewEdgeId(if0,igrid),          &
          if0 = 1, nedge(igrid)),    &
          igrid = 1, ngrids)

     READ(changrid) ((fNewEdgeIdInverse(if0,igrid),          &
          if0 = 1, nedge(igrid)),    &
          igrid = 1, ngrids)

     READ(changrid) ((fNewVertId(if0,igrid),          &
          if0 = 1, nvert(igrid)),    &
          igrid = 1, ngrids)

     READ(changrid) ((fNewVertIdInverse(if0,igrid),          &
          if0 = 1, nvert(igrid)),    &
          igrid = 1, ngrids)

     WRITE(*,*) "Done reading SFC information"
  END IF



  ! -------------------------------------------------------------------

END SUBROUTINE readgrid

!     ===============================================================
!
SUBROUTINE LL2XYZ(LONG,LAT,X,Y,Z)
  !
  !     To convert longitude and latitude to cartesian coordinates
  !     on the unit sphere
  !
  IMPLICIT NONE
  !
  REAL*8 LONG,LAT,X,Y,Z,CLN,SLN,CLT,SLT
  !
  !     ------------------------------------------------------------------
  !
  SLN=SIN(LONG)
  CLN=COS(LONG)
  SLT=SIN(LAT)
  CLT=COS(LAT)
  !
  X=CLN*CLT
  Y=SLN*CLT
  Z=SLT
  !
  !     ------------------------------------------------------------------
  !
  RETURN
END SUBROUTINE LL2XYZ
!
!     ==================================================================
!
SUBROUTINE XYZ2LL(X,Y,Z,LONG,LAT)
  !
  !     To convert cartesian coordinates to longitude and latitude
  !
  IMPLICIT NONE
  !   
  REAL*8 X,Y,Z,LONG,LAT,PI,TLN,TLT,R
  !
  !     -------------------------------------------------------------------
  !
  PI=4.0D0*ATAN(1.0D0)
  !
  IF (X.EQ.0.0D0) THEN
     IF (Y.GE.0.0D0) THEN
        LONG=0.5D0*PI
     ELSE
        LONG=1.5D0*PI
     ENDIF
  ELSE
     TLN=Y/X
     LONG=ATAN(TLN)
     IF (X.LT.0.0D0) THEN
        LONG=LONG+PI
     ENDIF
     IF (LONG.LT.0.0D0) THEN
        LONG=LONG+2.0D0*PI
     ENDIF
  ENDIF
  !
  R=SQRT(X*X+Y*Y)
  IF (R.EQ.0.0D0) THEN
     IF (Z.GT.0.0D0) THEN
        LAT=0.5D0*PI
     ELSE
        LAT=-0.5D0*PI
     ENDIF
  ELSE
     TLT=Z/R
     LAT=ATAN(TLT)
  ENDIF
  !
  !     --------------------------------------------------------------------
  !
  RETURN
END SUBROUTINE XYZ2LL
!
!     ====================================================================
!
SUBROUTINE STAREA2(X0,Y0,Z0,X1,Y1,Z1,X2,Y2,Z2,AREA)
  !
  !     Calculate the area of the spherical triangle whose corners
  !     have Cartesian coordinates (X0,Y0,Z0), (X1,Y1,Z1), (X2,Y2,Z2)
  !     The formula below is more robust to roundoff error than the
  !     better known sum of angle - PI formula
  !
  IMPLICIT NONE
  !
  REAL*8 X0,Y0,Z0,X1,Y1,Z1,X2,Y2,Z2, &
       D0,D1,D2,S,T0,T1,T2,T3,AREA
  !
  !
  !     Distances between pairs of points
  CALL SPDIST(X0,Y0,Z0,X1,Y1,Z1,D2)
  CALL SPDIST(X1,Y1,Z1,X2,Y2,Z2,D0)
  CALL SPDIST(X2,Y2,Z2,X0,Y0,Z0,D1)
  !
  !     Half perimeter
  S=0.5D0*(D0+D1+D2)
  !
  !     Tangents
  T0 = TAN(0.5D0*(S-D0))
  T1 = TAN(0.5D0*(S-D1))
  T2 = TAN(0.5D0*(S-D2))
  T3 = TAN(0.5D0*S)
  !
  !     Area
  AREA = 4.0D0*ATAN(SQRT(T0*T1*T2*T3))
  !
  RETURN
END SUBROUTINE STAREA2
!
!     ===================================================================
!
SUBROUTINE SPDIST(X1,Y1,Z1,X2,Y2,Z2,S)
  !
  !     Calculate the spherical distance S between two points with Cartesian
  !     coordinates (X1,Y1,Z1), (X2,Y2,Z2) on the unit sphere

  IMPLICIT NONE

  REAL*8 X1, Y1, Z1, X2, Y2, Z2, S, DX, DY, DZ, AD


  DX = X2 - X1
  DY = Y2 - Y1
  DZ = Z2 - Z1
  AD = SQRT(DX*DX + DY*DY + DZ*DZ)
  S = 2.0D0*ASIN(0.5D0*AD)


  RETURN
END SUBROUTINE SPDIST
!
!     ===================================================================

SUBROUTINE centroid(if0,long,lat)

  ! Find the centroid of cell if0 on grid ngrids

  USE grid
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: if0
  REAL*8, INTENT(OUT) :: long, lat
  INTEGER :: ixe, ie1, iv1, iv2
  REAL*8 :: long1, lat1, x0, y0, z0, x1, y1, z1, x2, y2, z2, &
       xc, yc, zc, a, aby3, mag

  ! -----------------------------------------------------------------------

  ! Coordinates of `centre' of face (i.e. dual vertex)
  long1 = flong(if0,ngrids)
  lat1 = flat(if0,ngrids)
  CALL ll2xyz(long1,lat1,x0,y0,z0)

  ! Loop over edges in turn and calculate area of triangle
  ! formed by the edge and the centre of the face
  ! Hence find area of face and centroid
  xc = 0.0d0
  yc = 0.0d0
  zc = 0.0d0
  DO ixe = 1, neoff(if0,ngrids)
     ie1 = eoff(ixe,if0,ngrids)
     iv1 = vofe(1,ie1,ngrids)
     iv2 = vofe(2,ie1,ngrids)
     long1 = vlong(iv1,ngrids)
     lat1 = vlat(iv1,ngrids)
     CALL ll2xyz(long1,lat1,x1,y1,z1)
     long1 = vlong(iv2,ngrids)
     lat1 = vlat(iv2,ngrids)
     CALL ll2xyz(long1,lat1,x2,y2,z2)
     CALL starea2(x0,y0,z0,x1,y1,z1,x2,y2,z2,a)
     aby3 = a/3.0d0
     xc = xc + (x0 + x1 + x2)*aby3
     yc = yc + (y0 + y1 + y2)*aby3
     zc = zc + (z0 + z1 + z2)*aby3
  ENDDO
  mag = SQRT(xc*xc + yc*yc + zc*zc)
  xc = xc/mag
  yc = yc/mag
  zc = zc/mag
  CALL xyz2ll(xc,yc,zc,long,lat)

  ! -----------------------------------------------------------------------

END SUBROUTINE centroid

!     ===================================================================

SUBROUTINE matinv(a,ainv,n)

  ! Invert an n x n matrix a by Gaussian elimination and put the
  ! result in matrix ainv. No fancy pivoting or checking for
  ! divide by zero!! Upon exit a should be the n x n identity matrix.

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  REAL*8, INTENT(INOUT) :: a(n,n)
  REAL*8, INTENT(OUT) :: ainv(n,n)
  INTEGER :: i, j
  REAL*8 :: temp

  ! -----------------------------------------------------------------------

  ! Initialize ainv to the identify matrix
  ainv = 0.0d0
  DO i = 1, n
     ainv(i,i) = 1.0d0
  ENDDO

  ! Eliminate subdiagonal terms
  DO j = 1, n
     temp = a(j,j)
     a(j,j:n) = a(j,j:n)/temp
     ainv(j,:) = ainv(j,:)/temp
     DO i = j+1, n
        temp = a(i,j)
        a(i,j:n) = a(i,j:n) - temp*a(j,j:n)
        ainv(i,:) = ainv(i,:) - temp*ainv(j,:)
     ENDDO
  ENDDO

  ! Eliminate above diagonal terms
  DO j = n, 2, -1
     DO i = 1, j-1
        temp = a(i,j)
        a(i,j) = a(i,j) - temp*a(j,j)   ! Should be zero if we've done it all right!
        ainv(i,:) = ainv(i,:) - temp*ainv(j,:)
     ENDDO
  ENDDO


END SUBROUTINE matinv

! =======================================================================

SUBROUTINE opendumpfiles

  ! Initialize input output tasks:
  ! Open files and write headers

  USE channels
  USE grid
  USE timestep
  IMPLICIT NONE

  REAL*4, ALLOCATABLE :: flongstar4(:), flatstar4(:), &
       vlongstar4(:), vlatstar4(:)
  CHARACTER*8 :: ytime
  CHARACTER*54 :: yname

  ! -----------------------------------------------------------------------

  ALLOCATE(flongstar4(nfacex), flatstar4(nfacex), &
       vlongstar4(nvertx), vlatstar4(nvertx))

  ! Construct timestep element of filename
  WRITE(ytime,'(I8.8)') istep

  ! File for Matlab primal grid output
  yname = 'dump/dump1_'//ytime//'.m'
  OPEN(chandumpm1,FILE=yname)

  ! Convert to single precision to reduce size of output
  ! and improve readability!
  flongstar4 = flong(:,ngrids)
  flatstar4 = flat(:,ngrids)
  ! Write header information
  WRITE(chandumpm1,*)   'npt = ',nfacex,';'
  WRITE(chandumpm1,*)   'nr = ',nvertx,';'
  WRITE(chandumpm1,*)   'nprmx = ',nevmx,';'
  WRITE(chandumpm1,*)   'rlist = [ ...'
  WRITE(chandumpm1,889) fofv(:,:,ngrids)
  WRITE(chandumpm1,*)   ' ];'
  WRITE(chandumpm1,*)   'rlist = reshape(rlist,nr,nprmx);'
  WRITE(chandumpm1,*)   'long = [ ...'
  WRITE(chandumpm1,888) flongstar4
  WRITE(chandumpm1,*)   ' ];'
  WRITE(chandumpm1,*)   'lat = [ ...'
  WRITE(chandumpm1,888) flatstar4
  WRITE(chandumpm1,*)   ' ];'


  ! File for Matlab primal grid output
  yname = 'dump/dump2_'//ytime//'.m'
  OPEN(chandumpm2,FILE=yname)

  ! Convert to single precision to reduce size of output
  ! and improve readability!
  vlongstar4 = vlong(:,ngrids)
  vlatstar4 = vlat(:,ngrids)
  ! Write header information
  WRITE(chandumpm2,*)   'npt = ',nvertx,';'
  WRITE(chandumpm2,*)   'nr = ',nfacex,';'
  WRITE(chandumpm2,*)   'nprmx = ',nefmx,';'
  WRITE(chandumpm2,*)   'rlist = [ ...'
  WRITE(chandumpm2,889) voff(:,:,ngrids)
  WRITE(chandumpm2,*)   ' ];'
  WRITE(chandumpm2,*)   'rlist = reshape(rlist,nr,nprmx);'
  WRITE(chandumpm2,*)   'long = [ ...'
  WRITE(chandumpm2,888) vlongstar4
  WRITE(chandumpm2,*)   ' ];'
  WRITE(chandumpm2,*)   'lat = [ ...'
  WRITE(chandumpm2,888) vlatstar4
  WRITE(chandumpm2,*)   ' ];'

888 FORMAT(E16.4)
889 FORMAT(I8)

  DEALLOCATE(flongstar4, flatstar4, vlongstar4, vlatstar4)



  ! -------------------------------------------------------------------

END SUBROUTINE opendumpfiles

! ===================================================================

SUBROUTINE closedumpfiles

  ! Close Matlab output files

  USE channels
  IMPLICIT NONE

  ! -------------------------------------------------------------------

  CLOSE(chandumpm1)
  CLOSE(chandumpm2)

  ! -------------------------------------------------------------------

END SUBROUTINE closedumpfiles

! ===================================================================

SUBROUTINE dumpm(q,ytitle,npt,ygrid)

  ! Dump a 2D field on an unstructured grid ready for reading and
  ! contouring by Matlab

  ! q(npt)                 The data to be output
  ! ytitle                 Title of data
  ! ygrid                  'primal' or 'dual  '

  USE channels
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: npt
  REAL*8, INTENT(IN) :: q(npt)
  CHARACTER*(*), INTENT(IN) :: ytitle
  CHARACTER*6, INTENT(IN) :: ygrid
  REAL*4, ALLOCATABLE :: qstar4(:)
  INTEGER chanm

  ALLOCATE(qstar4(npt))

  ! -----------------------------------------------------------------------

  ! Convert to single precision to reduce size of output
  ! and improve readability! (Shouldn't matter with formatted output)
  qstar4 = q

  IF (ygrid == 'primal') THEN
     chanm = chandumpm1
  ELSEIF (ygrid == 'dual  ') THEN
     chanm = chandumpm2
  ELSE
     PRINT *,'Invalid option: ygrid = ',ygrid,' in subroutine dumpm'
     STOP
  ENDIF

  ! Output
  WRITE(chanm,*)   'ytitle = ''',ytitle,''';'
  WRITE(chanm,*)   'q = [ ...'
  WRITE(chanm,888) qstar4
  WRITE(chanm,*)   ' ];'
  WRITE(chanm,*)   'jtcontour'

888 FORMAT(E16.4)

  ! -----------------------------------------------------------------------

  DEALLOCATE(qstar4)

END SUBROUTINE dumpm

! =======================================================================

SUBROUTINE output

  ! Standard output fields

  USE constants
  USE grid
  USE state
  USE errdiag
  USE timestep
  IMPLICIT NONE

  INTEGER :: nf, ne, nv
  REAL*8, ALLOCATABLE :: phi(:), zeta2(:), phiv2(:), pv(:), v1(:), mu(:)

  ALLOCATE(phi(nfacex), zeta2(nvertx), phiv2(nvertx), pv(nvertx), &
       v1(nedgex), mu(nedgex))

  ! -----------------------------------------------------------------------

  nf = nface(ngrids)
  ne = nedge(ngrids)
  nv = nvert(ngrids)


  ! Open output files and write header information
  CALL opendumpfiles

  ! Surface height
  ! phi = (phi2 + orog2)/(gravity*farea(:,ngrids))
  ! CALL dumpm(phi,'h',nf,'primal')
  ! Surface geopotential
  phi = (phi2 + orog2)/(farea(:,ngrids))
  CALL dumpm(phi,'phi',nf,'primal')
  ! phierr = (phi2 - phi2_init)/farea(:,ngrids)
  ! CALL dumpm(phierr,'phierr',nf,'primal')

  ! Vorticity and potential vorticity
  CALL massM(u1,mu,ngrids,ne)
  CALL HodgeHinv(mu,v1,ngrids,ne,-20)
  CALL Ddual2(v1,zeta2,ngrids,ne,nv)
  CALL operR(phi2,pv,ngrids,nf,nv)  ! Borrow pv array temporarily
  CALL HodgeJinv(pv,phiv2,ngrids,nv,-20)
  pv = (zeta2 + planvort2)/phiv2
  zeta2 = zeta2/varea(:,ngrids)
  CALL dumpm(zeta2,'Vorticity',nv,'dual  ')
  !CALL dumpm(pv,'PV',nv,'dual  ')
  !pv = xzeta2/phiv2
  !CALL dumpm(pv,'xPV',nv,'dual  ')

  ! Close output files
  CALL closedumpfiles
  PRINT*, "Data dumped in dump folder!"

  DEALLOCATE(phi, zeta2, phiv2, pv, v1, mu)

  ! -----------------------------------------------------------------------

END SUBROUTINE output

! =======================================================================

SUBROUTINE diffref

  ! Compute and output the surface height difference from a reference solution

  USE state
  USE timestep
  USE constants
  USE channels
  IMPLICIT NONE

  INTEGER, PARAMETER :: nreftime = 5
  INTEGER :: reftime(nreftime), ilist, if0
  REAL*8, ALLOCATABLE :: h(:), href(:)
  REAL*8 :: l1, l2, linf
  CHARACTER*12 :: yrefpre
  CHARACTER*10 :: yres, ytime
  CHARACTER*54 :: yname

  ALLOCATE(h(nfacex), href(nfacex))
  ! -----------------------------------------------------------------------

  ! Reference solution file prefixes
  IF (nefmx == 4) THEN
     yrefpre = 'TC5ref_cube_'
  ELSE
     yrefpre = 'TC5ref_hex__'
  ENDIF

  ! Resolution (for building file name)
  WRITE(yres,'(I10.10)') nface(ngrids)

  ! List of times at which reference solution is required
  reftime(1) = NINT(3600.0d0)              ! 1 hour
  reftime(2) = NINT(86400.0d0)             ! 1 day
  reftime(3) = NINT(5.0d0*86400.0)         ! 5 days
  reftime(4) = NINT(10.0D0*86400.0)        ! 10 days
  reftime(5) = NINT(15.0d0*86400.0)        ! 15 days

  ! Surface height
  h = (phi2 + orog2)/(gravity*farea(:,ngrids))

  DO ilist = 1, nreftime

     IF (NINT(istep*dt) == reftime(ilist)) THEN

        ! Reference solution is required at this time
        WRITE(ytime,'(I10.10)') reftime(ilist)
        ! Read reference solution
        OPEN(chanrefin,FILE=yrefpre//yres//'_'//ytime//'.dat',FORM='UNFORMATTED')
        DO if0 = 1, nface(ngrids)
           READ(chanrefin) href(if0)
        ENDDO
        CLOSE(chanrefin)

        ! File for Matlab primal grid output
        yname = 'dump/err1_'//ytime//'.m'
        OPEN(chanerrout,FILE=yname)

        ! Write header information
        WRITE(chanerrout,*)   'npt = ',nfacex,';'
        WRITE(chanerrout,*)   'nr = ',nvertx,';'
        WRITE(chanerrout,*)   'nprmx = ',nevmx,';'
        WRITE(chanerrout,*)   'rlist = [ ...'
        WRITE(chanerrout,889) fofv(:,:,ngrids)
        WRITE(chanerrout,*)   ' ];'
        WRITE(chanerrout,*)   'rlist = reshape(rlist,nr,nprmx);'
        WRITE(chanerrout,*)   'long = [ ...'
        WRITE(chanerrout,888) flong(:,ngrids)
        WRITE(chanerrout,*)   ' ];'
        WRITE(chanerrout,*)   'lat = [ ...'
        WRITE(chanerrout,888) flat(:,ngrids)
        WRITE(chanerrout,*)   ' ];'

        ! Output
        WRITE(chanerrout,*)   'ytitle = ''Reference h time ',ytime,''';'
        WRITE(chanerrout,*)   'q = [ ...'
        WRITE(chanerrout,888) href
        WRITE(chanerrout,*)   ' ];'
        WRITE(chanerrout,*)   'figure(1)'
        WRITE(chanerrout,*)   'jtcontour'
        WRITE(chanerrout,*)   'ytitle = ''Solution h time ',ytime,''';'
        WRITE(chanerrout,*)   'q = [ ...'
        WRITE(chanerrout,888) h
        WRITE(chanerrout,*)   ' ];'
        WRITE(chanerrout,*)   'figure(2)'
        WRITE(chanerrout,*)   'jtcontour'
        h = h - href
        WRITE(chanerrout,*)   'ytitle = ''h error time ',ytime,''';'
        WRITE(chanerrout,*)   'q = [ ...'
        WRITE(chanerrout,888) h
        WRITE(chanerrout,*)   ' ];'
        WRITE(chanerrout,*)   'figure(3)'
        WRITE(chanerrout,*)   'jtcontour'

        CLOSE(chanerrout)

        ! Error norms
        l1 = SUM(ABS(h)*farea(:,ngrids))/SUM(farea(:,ngrids))
        l2 = SQRT((SUM(h*h*farea(:,ngrids)))/SUM(farea(:,ngrids)))
        linf = MAXVAL(ABS(h))
        PRINT *,'h error: L1   = ',l1 ,'  L2   = ',l2 ,'  Linf = ',linf

     ENDIF

  ENDDO

888 FORMAT(E16.4)
889 FORMAT(I8)


  DEALLOCATE(h, href)
  ! -----------------------------------------------------------------------

END SUBROUTINE diffref

! =======================================================================

SUBROUTINE diagnostics

  ! To compute and output some basic diagnostics

  USE channels
  USE constants
  USE timestep
  USE state
  USE work
  USE errdiag
  IMPLICIT NONE

  INTEGER :: nf, ne, nv, ie0, if1, if2
  REAL*8 :: l1, l2, linf, meanphis, ape, ke, pz, temp, mass, &
       phibarerr, xpverr, l1u, l2u, linfu

  ! -----------------------------------------------------------------------

  nf = nface(ngrids)
  ne = nedge(ngrids)
  nv = nvert(ngrids)


  !! For test case 2 or other steady test cases
  phierr = (phi2 - phi2_init)/farea(:,ngrids)
  l1 = SUM(ABS(phierr)*farea(:,ngrids))/SUM(farea(:,ngrids))
  l2 = SQRT((SUM(phierr*phierr*farea(:,ngrids)))/SUM(farea(:,ngrids)))
  linf = MAXVAL(ABS(phierr))
  u1_err = (u1 - u1_init)/ldist(:,ngrids)
  l1u = SUM(ABS(u1_err)*ldist(:,ngrids)*ddist(:,ngrids)) &
       /SUM(ldist(:,ngrids)*ddist(:,ngrids))
  l2u = SQRT((SUM(u1_err*u1_err*ldist(:,ngrids)*ddist(:,ngrids))) &
       /SUM(ldist(:,ngrids)*ddist(:,ngrids)))
  linfu = MAXVAL(ABS(u1_err))
  print *,'Phierr: L1   = ',l1 ,'  L2   = ',l2 ,'  Linf = ',linf
  !print *,'Uerr:   L1   = ',l1u,'  L2   = ',l2u,'  Linf = ',linfu
  ! Write to file
  WRITE(chanerr,*) istep, l1, l2, linf


  ! Compute dual grid geopotential and compare with tracer
  CALL operR(phi2,pv,ngrids,nf,nv)  ! Borrow pv array for temporary use
  CALL HodgeJinv(pv,phibar2,ngrids,nv,-20)
  phibarerr = MAXVAL(ABS((phibar2 - xphibar2)/varea(:,ngrids))) &
       /MAXVAL(ABS(phibar2/varea(:,ngrids)))
  !PRINT *,'Normalized max difference between dual phi and tracer: ',phibarerr


  ! Compute PV and compare with tracer
  ! Relative vorticity
  CALL massM(u1,u1_temp,ngrids,ne)  ! Borrow u1_temp array for temporary use
  CALL HodgeHinv(u1_temp,v1,ngrids,ne,-20)
  CALL Ddual2(v1,zeta2,ngrids,ne,nv)
  ! Add planetary contribution
  zeta2 = zeta2 + planvort2
  xpverr = MAXVAL(ABS((zeta2 - xzeta2)/phibar2)) &
       /MAXVAL(ABS(zeta2/phibar2))
  !PRINT *,'Normalized max difference between dual PV and tracer: ', xpverr


  ! Total mass
  mass = SUM(phi2)
  !PRINT *,'Total mass = ',mass


  ! Available energy
  ! First find mean surface geopotential
  meanphis = SUM(phi2 + orog2)/SUM(farea(:,ngrids))
  ! Available potential energy
  b2 = phi2 + orog2 - meanphis*farea(:,ngrids)
  b0 = b2/farea(:,ngrids)
  ape = 0.5d0*SUM(b0*b2)
  ! Kinetic energy
  ! Cell integrals of KE
  CALL operT(u1,b2,ngrids,ne,nf)
  b2 = 0.5d0*b2
  b0 = b2/farea(:,ngrids)
  ke = SUM(b0*phi2)
  !PRINT *,'APE = ',ape,'    KE = ',ke,'    Total = ',ape + ke


  ! Potential enstrophy
  pv = zeta2/phibar2
  pz = 0.5d0*SUM(zeta2*pv)
  !PRINT *,'Potential enstrophy: ',pz

  ! Write to file
  WRITE(chandiag,888) istep, phibarerr, xpverr, mass, ape, ke, pz
888 FORMAT(I8,6E20.12)


  ! -----------------------------------------------------------------------

END SUBROUTINE diagnostics

! =======================================================================

SUBROUTINE dumpgrid

  ! Output the grid coordinates in a simple format for use in generating
  ! reference solutions

  USE runtype
  USE channels
  USE grid

  IMPLICIT NONE
  INTEGER :: if0
  REAL*8 :: long, lat

  ! -----------------------------------------------------------------------

  OPEN(chanrefgrd,FILE=ygridcoords,FORM='UNFORMATTED')

  WRITE(chanrefgrd) nface(ngrids)
  DO if0 = 1, nface(ngrids)
     ! Dual vertices
     ! long = flong(if0,ngrids)
     ! lat = flat(if0,ngrids)
     ! Or cell centroids
     CALL centroid(if0,long,lat)
     WRITE(chanrefgrd) long, lat
  ENDDO

  CLOSE(chanrefgrd)

  ! -----------------------------------------------------------------------

END SUBROUTINE dumpgrid

! =======================================================================

SUBROUTINE writerestart

  USE grid
  USE state
  USE channels
  USE runtype
  USE timestep
  IMPLICIT NONE
  CHARACTER*8 :: ystep
  CHARACTER*30 :: yresout
  INTEGER :: if0, ie0, iv0
  INTEGER :: outputMode

  ! -----------------------------------------------------------------------

  ! Construct name of restart file
  WRITE(ystep,'(I8.8)') istep
  yresout = 'run'//runid//'_restart_'//ystep//'.dat'



  ! 1: binary, SFC reordered (default)
  ! 2: raw
  ! 3: ascii, SFC reordered
  outputMode = 1

  IF (outputMode == 1) THEN
     ! Binary output

     ! Open file and write state
     OPEN(chanresout,FILE=yresout,FORM='UNFORMATTED')
     !OPEN(chanresout,FILE=yresout)
     WRITE(chanresout) istep, nfacex, nedgex, nvertx

     IF (SFCIndexAvailable == 0) THEN
        DO if0 = 1, nface(ngrids)
           WRITE(chanresout) phi2(if0)      ! cell
        END DO
        DO ie0 = 1, nedge(ngrids)
           WRITE(chanresout) u1(ie0)        ! edge
        END DO
        DO iv0 = 1, nvert(ngrids)
           WRITE(chanresout) xphibar2(iv0)  ! vertices
        END DO
        DO iv0 = 1, nvert(ngrids)
           WRITE(chanresout) xzeta2(iv0)    ! vertices
        END DO
     ELSE
        DO if0 = 1, nface(ngrids)
           WRITE(chanresout) phi2(fNewFaceIdInverse(if0, ngrids))
        END DO
        DO ie0 = 1, nedge(ngrids)
           WRITE(chanresout) u1(fNewEdgeIdInverse(ie0,ngrids))        ! edge
        END DO
        DO iv0 = 1, nvert(ngrids)
           WRITE(chanresout) xphibar2(fNewVertIdInverse(iv0,ngrids))  ! vertices
        END DO
        DO iv0 = 1, nvert(ngrids)
           WRITE(chanresout) xzeta2(fNewVertIdInverse(iv0,ngrids))    ! vertices
        END DO
     ENDIF

  END IF

  IF (outputMode == 2) THEN
     ! Binary output

     ! Open file and write state
     OPEN(chanresout,FILE=yresout,FORM='UNFORMATTED')
     !OPEN(chanresout,FILE=yresout)
     WRITE(chanresout) istep, nfacex, nedgex, nvertx

     WRITE(chanresout) phi2      ! cell
     WRITE(chanresout) u1        ! edge
     WRITE(chanresout) xphibar2  ! vertices
     WRITE(chanresout) xzeta2    ! vertices
  END IF


  IF (outputMode == 3) THEN
     ! Binary output
     ! ASCII text output, easier to compare results

     ! Open file and write state
     OPEN(chanresout,FILE=yresout)
     !OPEN(chanresout,FILE=yresout)
     WRITE(chanresout,*) istep, nfacex, nedgex, nvertx

     IF (SFCIndexAvailable == 0) THEN
        DO if0 = 1, nface(ngrids)
           WRITE(chanresout,*) phi2(if0)      ! cell
        END DO
        DO ie0 = 1, nedge(ngrids)
           WRITE(chanresout,*) u1(ie0)        ! edge
        END DO
        DO iv0 = 1, nvert(ngrids)
           WRITE(chanresout,*) xphibar2(iv0)  ! vertices
        END DO
        DO iv0 = 1, nvert(ngrids)
           WRITE(chanresout,*) xzeta2(iv0)    ! vertices
        END DO
     ELSE !NOCH FALSCH!
        DO if0 = 1, nface(ngrids)
           WRITE(chanresout,*) phi2(fNewFaceIdInverse(if0, ngrids))
        END DO
        DO ie0 = 1, nedge(ngrids)
           WRITE(chanresout,*) u1(ie0)        ! edge
        END DO
        DO iv0 = 1, nvert(ngrids)
           WRITE(chanresout,*) xphibar2(iv0)  ! vertices
        END DO
        DO iv0 = 1, nvert(ngrids)
           WRITE(chanresout,*) xzeta2(iv0)    ! vertices
        END DO
     ENDIF
  END IF


  ! Close file
  CLOSE(chanresout)


  ! -----------------------------------------------------------------------

END SUBROUTINE writerestart

! =======================================================================

SUBROUTINE readrestart

  USE state
  USE channels
  USE runtype
  USE timestep
  IMPLICIT NONE
  INTEGER :: nf, ne, nv

  ! -----------------------------------------------------------------------

  ! Open file and read and check data
  OPEN(chanresin,FILE=yresfile,FORM='UNFORMATTED')
  READ(chanresin) istep, nf, ne, nv
  IF ((nf.NE.nfacex).OR.(ne.NE.nedgex).OR.(nv.NE.nvertx)) THEN
     PRINT *,'Error in readrestart:'
     PRINT *,'nf = ',nf,' nfacex = ',nfacex
     PRINT *,'ne = ',ne,' nedgex = ',nedgex
     PRINT *,'nv = ',nv,' nvertx = ',nvertx
     STOP
  ENDIF
  IF (SFCIndexAvailable == 0) THEN
     WRITE(*,*) "TODO: reorder input along the SFC"
     CALL EXIT(1)
  END IF
  READ(chanresin) phi2
  READ(chanresin) u1
  READ(chanresin) xphibar2
  READ(chanresin) xzeta2

  ! Close file
  CLOSE(chanresin)

  ! -----------------------------------------------------------------------

END SUBROUTINE readrestart

! =======================================================================

SUBROUTINE addcompmode

  ! Construct a grid scale pattern in a stream function
  ! and use this to perturb the velocity field.

  USE state
  IMPLICIT NONE
  INTEGER :: nf, ne, nv, iter
  REAL*8, ALLOCATABLE :: psi(:), du(:), dv(:), dz(:), dpsi(:)
  REAL*8 :: amp

  ALLOCATE(psi(nvertx), du(nedgex), dv(nedgex), dz(nvertx), dpsi(nvertx))

  ! -----------------------------------------------------------------------

  nf = nface(ngrids)
  ne = nedge(ngrids)
  nv = nvert(ngrids)

  ! Initialize to zero
  psi = 0.0d0

  ! Plant seeds:
  ! Vertex near (pi,0)
  psi(13110) = 1.0d0
  ! Vertex near (0,pi)
  psi(16593) = 1.0d0

  ! Grow noise by antidiffusion
  DO iter = 1, 20
     CALL Dprimal1(psi,du,ngrids,nv,ne)
     CALL massM(du,dv,ngrids,ne)
     CALL Ddual2(dv,dz,ngrids,ne,nv)
     dz = -dz
     ! Last step is a hack because we don't have the E_p mass matrix
     ! needed to do a proper scalar Laplacian in E_p
     dpsi = dz/varea(:,ngrids)
     PRINT *,'maxval dpsi = ',maxval(dpsi)
     psi = psi - 2.0E9*dpsi
     psi = MIN(psi,1.0d0)
     psi = MAX(psi,-1.0d0)
  ENDDO

  ! Now compute velocity increments
  CALL Dprimal1(psi,du,ngrids,nv,ne)

  ! Normalize to desired amplitude
  amp = MAXVAL(ABS(du)/ldist(:,ngrids))
  du = du/amp

  ! And add to velocity field
  u1 = u1 + du

  DEALLOCATE(psi, du, dv, dz, dpsi)

  ! -----------------------------------------------------------------------

END SUBROUTINE addcompmode

! =======================================================================

SUBROUTINE testop

  ! To test various exterior derivative and Hodge star operators

  USE grid
  USE constants
  IMPLICIT NONE

  INTEGER :: igrid, nf, ne, nv, if1, iv1
  REAL*8, ALLOCATABLE :: ff1(:), ff2(:), ff3(:), &
       fe1(:), fe2(:), fe3(:), &
       fv1(:), fv2(:), fv3(:)
  REAL*8 :: long, lat, mxerr

  ! Which properties to check ?
  LOGICAL :: ldata = .FALSE.        ! Test map from primal to dual points
  LOGICAL :: ldd = .FALSE.          ! DD f = 0 on primal and dual grids
  LOGICAL :: linv = .FALSE.         ! Inverses of L, J, M and H operators
  LOGICAL :: llaplacian = .TRUE.    ! Accuracy of scalar and vector Laplacians
  LOGICAL :: lpd = .FALSE.          ! Mimetic and primal-dual-related properties

  INTEGER :: ix, ixmx

  ! -------------------------------------------------------------------

  ! Which grids shall we test?
  DO igrid = 1, ngrids
     ! DO igrid = 2, 2

     nf = nface(igrid)
     ne = nedge(igrid)
     nv = nvert(igrid)

     PRINT *,' '
     PRINT *,'--------------------------'
     PRINT *,' '
     PRINT *,'GRID ',igrid,'  nf = ',nf,'  ne = ',ne,'  nv = ',nv,'  dof = ',nf + ne
     PRINT *,' '

     PRINT *,'Maximum ddist = ',MAXVAL(ddist(1:ne,igrid))
     PRINT *,'max l / min l = ',MAXVAL(ldist(1:ne,igrid))/MINVAL(ldist(1:ne,igrid)), &
          '  max d / min d = ',MAXVAL(ddist(1:ne,igrid))/MINVAL(ddist(1:ne,igrid)), &
          '  max A / min A = ',MAXVAL(farea(1:nf,igrid))/MINVAL(farea(1:nf,igrid))
     PRINT *,' '

     ALLOCATE(ff1(nf),ff2(nf),ff3(nf))
     ALLOCATE(fe1(ne),fe2(ne),fe3(ne))
     ALLOCATE(fv1(nv),fv2(nv),fv3(nv))

     ! Set up test data
     DO if1 = 1, nf
        long = flong(if1,igrid)
        lat = flat(if1,igrid)
        ! ff1(if1) = SIN(lat)
        ff1(if1) = COS(lat)*SIN(long)
        ! ff1(if1) = 1.0d0
     ENDDO

     DO iv1 = 1, nv
        long = vlong(iv1,igrid)
        lat = vlat(iv1,igrid)
        ! fv1(iv1) = SIN(lat)
        fv1(iv1) = COS(lat)*SIN(long)
        ! fv1(iv1) = 1.0d0
     ENDDO


     IF (ldata) THEN

        ! Check ff1 and fv1 are consistent
        ff2 = ff1*farea(:,igrid)
        CALL operR(ff2,fv2,igrid,nf,nv)
        CALL HodgeJinv(fv2,fv3,igrid,nv,-20)
        fv2 = fv3/varea(:,igrid)
        PRINT *,'Check dual grid data consistent with primal'
        PRINT *,'fv1 = ',fv1(1:20)
        PRINT *,'fv2 = ',fv2(1:20)
        PRINT *,'fv2 should approx equal fv1 '
        PRINT *,'Biggest diff ',MAXVAL(ABS(fv2 - fv1))
        PRINT *,'  '

     ENDIF

     IF (ldd) THEN

        ! Check DD = 0
        CALL Ddual1(ff1,fe1,igrid,nf,ne)
        CALL Ddual2(fe1,fv2,igrid,ne,nv)
        PRINT *,'ff1 = ',ff1(1:40)
        PRINT *,'DD ff1 = ',fv2(1:40)
        PRINT *,'DD ff1 should exactly equal zero'
        PRINT *,' '

        CALL Dprimal1(fv1,fe1,igrid,nv,ne)
        CALL Dprimal2(fe1,ff2,igrid,ne,nf)
        PRINT *,'fv1 = ',fv1(1:40)
        PRINT *,'DD fv1 = ',ff2(1:40)
        PRINT *,'DD fv1 should exactly equal zero'
        PRINT *,' '

     ENDIF

     IF (linv) THEN

        ! Check inverse of L
        CALL massLinv(ff1,ff2,igrid,nf,-20)
        CALL massL(ff2,ff3,igrid,nf)
        PRINT *,'Inverse of L'
        PRINT *,'ff1 = ',ff1(1:40)
        PRINT *,'ff3 = ',ff3(1:40)
        PRINT *,'ff3 should (almost) exactly equal ff1'
        PRINT *,'Biggest diff ',MAXVAL(ABS(ff3 - ff1))
        PRINT *,' '

        ! Check inverse of J
        CALL HodgeJinv(fv1,fv2,igrid,nv,-20)
        CALL HodgeJ(fv2,fv3,igrid,nv)
        PRINT *,'Inverse of J'
        PRINT *,'fv1 = ',fv1(1:40)
        PRINT *,'fv3 = ',fv3(1:40)
        PRINT *,'fv3 should (almost) exactly equal fv1'
        PRINT *,'Biggest diff ',MAXVAL(ABS(fv3 - fv1))
        PRINT *,' '

        ! Check inverse of M
        CALL Ddual1(ff1,fe1,igrid,nf,ne)
        CALL massM(fe1,fe2,igrid,ne)
        CALL massMinv(fe2,fe3,igrid,ne,-20)
        PRINT *,'Inverse of M'
        PRINT *,'fe1 = ',fe1(1:40)
        PRINT *,'fe3 = ',fe3(1:40)
        PRINT *,'fe3 should (almost) exactly equal fe1'
        PRINT *,'Biggest diff ',MAXVAL(ABS(fe3 - fe1))
        PRINT *,' '
        CALL approxMinv(fe2,fe3,igrid,ne)
        PRINT *,'Approximate inverse of M'
        PRINT *,'fe4 = ',fe3(1:40)
        PRINT *,'fe4 should be an approximation to fe1'
        PRINT *,'Biggest diff ',MAXVAL(ABS(fe3 - fe1))
        PRINT *,' '

        ! Check inverse of H
        CALL Ddual1(ff1,fe1,igrid,nf,ne)
        CALL HodgeH(fe1,fe2,igrid,ne)
        CALL HodgeHinv(fe2,fe3,igrid,ne,-20)
        PRINT *,'Inverse of H'
        PRINT *,'fe1 = ',fe1(1:40)
        PRINT *,'fe3 = ',fe3(1:40)
        PRINT *,'fe3 should (almost) exactly equal fe1'
        PRINT *,'Biggest diff ',MAXVAL(ABS(fe3 - fe1))
        PRINT *,' '

     ENDIF

     IF (llaplacian) THEN

        ! Check accuracy of lumped M compared to M
        CALL Dprimal1(fv1,fe2,igrid,nv,ne)

        fe3 = mlump(1:ne,igrid)*fe2
        PRINT *,'mlump*u = ',fe3(1:20)
        CALL massM(fe2,fe1,igrid,ne)
        PRINT *,'    M*u = ',fe1(1:20)
        fe3 = fe3 - fe1
        PRINT *,' error  = ',fe3(1:40)
        PRINT *,'max error = ',MAXVAL(ABS(fe3)) 

        CALL massM(fe2,fe3,igrid,ne)
        PRINT *,'fe2 M fe2 = ',SUM(fe3*fe2)
        PRINT *,'Exact answer = ',8.0d0*pi/3.0d0,'    Err = ',SUM(fe3*fe2) - 8.0d0*pi/3.0d0
        PRINT *,' '
        fe3 = fe2*mlump(1:ne,igrid)
        PRINT *,'fe2 Mlump fe2 = ',SUM(fe3*fe2)
        PRINT *,'Exact answer = ',8.0d0*pi/3.0d0,'    Err = ',SUM(fe3*fe2) - 8.0d0*pi/3.0d0
        PRINT *,' '
        CALL Ddual1(ff1,fe2,igrid,nf,ne)
        CALL massMinv(fe2,fe3,igrid,ne,-20)
        PRINT *,'fe2 Minv fe2 = ',SUM(fe3*fe2)
        PRINT *,'Exact answer = ',8.0d0*pi/3.0d0,'    Err = ',SUM(fe3*fe2) - 8.0d0*pi/3.0d0
        PRINT *,' '

        ! Check Laplacian of scalar (in Vp)
        ff2 = ff1*farea(:,igrid)
        CALL massL(ff2,ff3,igrid,nf)
        CALL Ddual1(ff3,fe2,igrid,nf,ne)
        CALL massMinv(fe2,fe3,igrid,ne,-20)
        CALL Dprimal2(fe3,ff2,igrid,ne,nf)
        ff3 = ff2/farea(:,igrid)
        ff3 = ff3*rearth*rearth
        ff2 = -2.0d0*ff1  
        PRINT *,'Laplacian of scalar (primal)'
        PRINT *,'True:   ff2 = ',ff2(1:40)
        PRINT *,'Approx: ff3 = ',ff3(1:40)
        PRINT *,'ff3 should approximately equal ff2'
        ff3 = ff3 - ff2
        PRINT *,'Err = ',ff3(1:40)
        PRINT *,'Max error = ',MAXVAL(ABS(ff3))
        PRINT *,'RMS error = ',SQRT((SUM(ff3*ff3))/nf)
        mxerr = 0.0d0
        DO ix = 1, nf
           IF (ABS(ff3(ix)) > mxerr) THEN
              ixmx = ix
              mxerr = ABS(ff3(ix))
           ENDIF
        ENDDO
        PRINT *,'max err ',mxerr,' occurs in cell ',ixmx
        PRINT *,' '

        ! Check Laplacian of a scalar (in Ep)
        CALL Dprimal1(fv1,fe1,igrid,nv,ne)
        CALL massM(fe1,fe2,igrid,ne)
        CALL Ddual2(fe2,fv2,igrid,ne,nv)
        CALL HodgeJinv(fv2,fv3,igrid,nv,-20)
        fv3 = -fv3*rearth*rearth/varea(:,igrid)
        fv2 = -2.0d0*fv1
        PRINT *,'Laplacian of scalar (dual)'
        PRINT *,'fv2 = ',fv2(1:40)
        PRINT *,'fv3 = ',fv3(1:40)
        PRINT *,'fv3 should approximately equal fv2'
        fv2 = fv3 - fv2
        PRINT *,'Err = ',fv2(1:40)
        PRINT *,'Max error = ',MAXVAL(ABS(fv2))
        PRINT *,'RMS error = ',SQRT((SUM(fv2*fv2))/nv) 
        PRINT *,' '

        ! Check Laplacian of vector
        !  CALL Ddual1(ff1,fe1,igrid,nf,ne)
        !  CALL Dprimal1(fv1,fe2,igrid,nv,ne)
        !  CALL HodgeHinv(fe2,fe3,igrid,ne)
        !  fe1 = fe1 - fe3
        !  CALL graddiv(fe1,fe2,igrid,ne,nf)
        !  CALL curlcurl(fe1,fe3,igrid,ne,nv)
        !  fe2 = fe2 + fe3
        !  fe3 = -2.0*fe1
        !  fe2 = fe2*rearth*rearth
        !  print *,'Laplacian of vector'
        !!  print *,'fe3 = ',fe3(1:40)
        !!  print *,'fe2 = ',fe2(1:40)
        !!  print *,'fe2 should approximately equal fe3'
        !  fe2 = fe2 - fe3
        !!  print *,'Err = ',fe2(1:40)
        !  print *,'Max error = ',MAXVAL(ABS(fe2))
        !  print *,'RMS error = ',SQRT((SUM(fe2*fe2))/ne) 
        !  print *,' '

     ENDIF

     ! Check re-ordered W calculation agrees with original
     !  CALL Dprimal1(fv1,fe1,igrid,nv,ne)
     !  CALL operW_original(fe1,fe2,igrid,ne)
     !  CALL operW(fe1,fe3,igrid,ne)
     !  print *,'fe2 = ',fe2(1:20)
     !  print *,'fe3 = ',fe3(1:20)
     !  print *,'Diff = ',fe3(1:20) - fe2(1:20)
     !  print *,' '


     IF (lpd) THEN

        ! Check TRiSK property
        CALL Dprimal1(fv1,fe1,igrid,nv,ne)
        CALL operW(fe1,fe2,igrid,ne)
        CALL Ddual2(fe2,fv2,igrid,ne,nv)
        PRINT *,'fv1 = ',fv1(1:40)
        PRINT *,'DWD fv1 = ',fv2(1:40)
        PRINT *,'DWD fv1 should exactly equal zero'
        PRINT *,'Max Abs value = ',MAXVAL(ABS(fv2))
        PRINT *,' '

        CALL Ddual1(ff1,fe1,igrid,nf,ne)
        CALL HodgeH(fe1,fe2,igrid,ne)
        CALL Dprimal2(fe2,ff2,igrid,ne,nf)
        CALL operR(ff2,fv2,igrid,nf,nv)
        CALL operW(fe2,fe1,igrid,ne)
        CALL Ddual2(fe1,fv3,igrid,ne,nv)
        PRINT *,'fv2 = RDHD ff1 = ',fv2(1:40)
        PRINT *,'fv3 = DWHD ff1 = ',fv3(1:40)
        PRINT *,'fv3 should exactly equal -fv2'
        PRINT *,'Max Abs diff = ',MAXVAL(ABS(fv3 + fv2))
        PRINT *,' '

        ! Check antisymmetry of W
        CALL Dprimal1(fv1,fe1,igrid,nv,ne)
        CALL operW(fe1,fe2,igrid,ne)
        fe2 = fe1*fe2
        PRINT *,'fe1 W fe1 = ',SUM(fe2),'   should be zero'
        PRINT *,' '

        ! Check dual-primal map commutes with derivative
        CALL Ddual1(ff1,fe1,igrid,nf,ne)
        CALL HodgeH(fe1,fe2,igrid,ne)
        CALL Ddual2(fe2,fv2,igrid,ne,nv)
        PRINT *,'ff1 = ',ff1(1:40)
        PRINT *,'DHD ff1 = ',fv2(1:40)
        PRINT *,'DHD ff1 should exactly equal zero'
        PRINT *,'Max Abs value = ',MAXVAL(ABS(fv2))
        PRINT *,' '

        CALL Dprimal1(fv1,fe2,igrid,nv,ne)
        CALL massM(fe2,fe1,igrid,ne)
        CALL Ddual2(fe1,fv2,igrid,ne,nv)
        CALL HodgeJ(fv2,fv3,igrid,nv)
        CALL HodgeH(fe1,fe2,igrid,ne)
        CALL Ddual2(fe2,fv2,igrid,ne,nv)
        PRINT *,'JD fe1 = ',fv3(1:40)
        PRINT *,'DH fe1 = ',fv2(1:40)
        PRINT *,'JD fe1 should exactly equal DH fe1'
        PRINT *,'Max Abs diff = ',MAXVAL(ABS(fv3 - fv2))
        PRINT *,' '

        CALL Dprimal1(fv1*rearth,fe1,igrid,nv,ne)
        CALL operW(fe1,fe2,igrid,ne)
        CALL Ddual1(ff1*rearth,fe1,igrid,nf,ne)
        CALL HodgeH(fe1,fe3,igrid,ne)
        fe2 =  fe2/ldist(:,igrid)
        fe3 = -fe3/ldist(:,igrid)
        PRINT *,' W D psi(v) = ',fe2(1:40)
        PRINT *,'-H D psi(f) = ',fe3(1:40)
        fe2 = fe2 - fe3
        PRINT *,'Accuracy of W: rotational flow'
        PRINT *,'Max error = ',MAXVAL(ABS(fe2))
        PRINT *,'RMS error = ',SQRT(SUM(fe2*fe2)/ne)
        PRINT *,' '

     ENDIF

     DEALLOCATE(ff1,ff2,ff3,fe1,fe2,fe3,fv1,fv2,fv3)

  ENDDO

  ! -------------------------------------------------------------------

END SUBROUTINE testop

!     ================================================================

SUBROUTINE fsphere

  ! To construct system matrix for linearized equations on an f-sphere
  ! for subsequent calculation of eigenvalues

  USE grid
  IMPLICIT NONE

  INTEGER :: igrid, nf, ne, nv, nmat, if1, ie1, i, j
  REAL*8, ALLOCATABLE :: phi(:), u(:), vperp(:), temp(:), a(:,:), &
       phidot(:), vdot(:), udot(:)
  REAL*8 :: f0 = 1.4584d-4, phi0 = 1.0d5

  ! -------------------------------------------------------------------

  ! Which grid shall we test?
  igrid = 2

  nf = nface(igrid)
  ne = nedge(igrid)
  nv = nvert(igrid)
  nmat = nf + ne
  PRINT *,'igrid = ',igrid,' nmat = ',nmat

  ALLOCATE(phi(nf),u(ne),vperp(ne),temp(nf),a(nmat,nmat),phidot(nf),vdot(ne),udot(ne))

  ! Set each input variable in turn to 1
  DO if1 = 1, nf
     phi = 0.0d0
     u = 0.0d0
     phi(if1) = 1.0d0

     ! Compute u and phi tendencies
     CALL Dprimal2(u,temp,igrid,ne,nf)
     phidot = -phi0*temp

     CALL massL(phi,temp,igrid,nf)
     CALL Ddual1(temp,vdot,igrid,nf,ne)
     vdot = -vdot

     CALL operW(u,vperp,igrid,ne)
     vdot = vdot - f0*vperp

     CALL massMinv(vdot,udot,ngrids,ne,-20)

     ! And save as column of system matrix
     a(1:nf,if1) = phidot
     a(nf+1:nmat,if1) = udot
  ENDDO

  ! Set each input variable in turn to 1
  DO ie1 = 1, ne
     phi = 0.0d0
     u = 0.0d0
     u(ie1) = 1.0d0

     ! Compute u and phi tendencies
     CALL Dprimal2(u,temp,igrid,ne,nf)
     phidot = -phi0*temp

     CALL massL(phi,temp,igrid,nf)
     CALL Ddual1(temp,vdot,igrid,nf,ne)
     vdot = -vdot

     CALL operW(u,vperp,igrid,ne)
     vdot = vdot - f0*vperp

     CALL massMinv(vdot,udot,ngrids,ne,-20)

     ! And save as column of system matrix
     a(1:nf,nf+ie1) = phidot
     a(nf+1:nmat,nf+ie1) = udot
  ENDDO


  ! Write out matrix for use by Matlab
  OPEN(24,FILE='aaa.m')
  WRITE(24,*) 'nmat = ',nmat
  WRITE(24,*) 'a = ['
  DO j = 1, nmat
     DO i = 1, nmat
        WRITE(24,*) a(i,j)
     ENDDO
  ENDDO
  WRITE(24,*) '];'


  DEALLOCATE(phi,u,vperp,temp,a,phidot,vdot,udot)

  ! -------------------------------------------------------------------

END SUBROUTINE fsphere

! ====================================================================

SUBROUTINE nmodes

  ! To construct system matrix for linearized equations on a rotating sphere
  ! for subsequent calculation of eigenvalues. A solid body rotation basic
  ! state should be set up in initial.

  USE grid
  USE state
  USE errdiag
  IMPLICIT NONE

  INTEGER :: nf, ne, nv, nmat, if1, ie1, i, j
  REAL*8, ALLOCATABLE :: u1exac(:), a(:,:)
  REAL*8 :: pert

  ! -------------------------------------------------------------------

  nf = nface(ngrids)
  ne = nedge(ngrids)
  nv = nvert(ngrids)
  nmat = nf + ne
  PRINT *,'grid = ',ngrids,' nmat = ',nmat

  ALLOCATE(u1exac(ne),a(nmat,nmat))


  ! First compute effect of time stepping the basic state (because
  ! it won't be perfectly steady)
  phi2 = phi2_init
  u1 = u1_init
  CALL step
  phiexac = phi2
  u1exac = u1


  ! Perturb each input variable in turn
  ! First the mass field
  DO if1 = 1, nf
     phi2 = phi2_init
     u1 = u1_init
     pert = farea(if1,ngrids)
     phi2(if1) = phi2(if1) + pert

     ! Time step perturbed state
     CALL step

     ! And save as column of system matrix
     a(1:nf,if1) = (phi2 - phiexac)/pert
     a(nf+1:nmat,if1) = (u1 - u1exac)/pert
  ENDDO

  ! Now the velocity field
  DO ie1 = 1, ne
     phi2 = phi2_init
     u1 = u1_init
     pert = ldist(ie1,ngrids)
     u1(ie1) = u1(ie1) + pert

     ! Time step perturbed state
     CALL step

     ! And save as column of system matrix
     a(1:nf,nf+ie1) = (phi2 - phiexac)/pert
     a(nf+1:nmat,nf+ie1) = (u1 - u1exac)/pert
  ENDDO


  ! Write out matrix for use by Matlab
  OPEN(24,FILE='aaa.m')
  WRITE(24,*) 'nmat = ',nmat
  WRITE(24,*) 'a = ['
  DO j = 1, nmat
     DO i = 1, nmat
        WRITE(24,*) a(i,j)
     ENDDO
  ENDDO
  WRITE(24,*) '];'


  DEALLOCATE(u1exac,a)

  ! -------------------------------------------------------------------

END SUBROUTINE nmodes

! ====================================================================

SUBROUTINE testmg

  ! To test the multigrid Helmholtz solver

  USE grid
  USE helmcoeff    ! Just for testing
  IMPLICIT NONE

  INTEGER :: nf, ne, nv, if1, iv1
  REAL*8, ALLOCATABLE :: ff1(:), ff2(:), ff3(:)
  REAL*8 :: long, lat

  ! -------------------------------------------------------------------

  nf = nface(ngrids)
  ne = nedge(ngrids)
  nv = nvert(ngrids)

  PRINT *,' '
  PRINT *,'--------------------------'
  PRINT *,' '
  PRINT *,'Testing mgsolve '
  PRINT *,' '

  ALLOCATE(ff1(nf),ff2(nf),ff3(nf))


  ! Build coefficients used in Helmholtz operator on all grids
  phiref = 1.0d5*farea(:,ngrids)
  CALL buildhelm


  ! Set up test data
  ! Large-scale part
  DO if1 = 1, nf
     long = flong(if1,ngrids)
     lat = flat(if1,ngrids)
     ! ff2(if1) = SIN(lat)
     ff2(if1) = COS(lat)*SIN(long)
  ENDDO
  ! Convert to area integrals
  ff1 = ff2/farea(:,ngrids)
  ! Plus small-scale part
  ff1(10) = 10.0d0*ff1(10)
  PRINT *,'Original field ff1 =     ',ff1(1:40)
  PRINT *,' '

  CALL helmholtz(ff1,ff2,ngrids,nf,ne)
  PRINT *,'ff2 = Helm(ff1) =        ',ff2(1:40)
  PRINT *,' '

  CALL mgsolve(ff3,ff2,ngrids)
  PRINT *,'Soln of Helmholtz ff3 = ', ff3(1:40)
  PRINT *,' '


END SUBROUTINE testmg

! ====================================================================

SUBROUTINE testadv

  USE channels
  USE constants
  USE state
  USE work
  USE errdiag
  USE timestep
  USE advection
  IMPLICIT NONE

  INTEGER :: if0, iv0, nf, ne, nv
  REAL*8 :: long, lat, &
       xc, yc, zc, x0, y0, z0, u00, r, r0, h0, alpha, ca, sa, &
       clat, slat, slon, &
       lonrot, l1, l2, linf

  REAL*8, ALLOCATABLE :: psi(:), flx1(:), dphi2(:), phibar(:), dzeta2(:)


  !real*8 :: div2(nfacex), zeta2(nvertx) ! Now in module work
  REAL*8 :: q, qmax
  INTEGER :: ivmx, ifmx


  ALLOCATE(psi(nvertx), flx1(nedgex), dphi2(nfacex), phibar(nvertx), dzeta2(nvertx))

  ! --------------------------------------------------------------------

  nf = nface(ngrids)
  ne = nedge(ngrids)
  nv = nvert(ngrids)

  ! Constant geopotential
  !phi2 = 1.0d0*farea(:,ngrids)

  ! Dual grid geopotential
  !CALL operR(phi2,xphibar2,ngrids,nf,nv)
  !CALL HodgeJinv(xphibar2,phibar,ngrids,nv,-20)
  !phibar = phibar/varea(:,ngrids)

  ! Cosine bell: Williamson et al 1992
  long = 1.0d0*pi
  lat = 0.0d0
  CALL ll2xyz(long,lat,xc,yc,zc)
  r0 = 1.0d0/3.0d0
  h0 = 1000.0d0
  qmax = 0.0d0

  DO if0 = 1, nf
     long = flong(if0,ngrids)
     lat = flat(if0,ngrids)
     CALL ll2xyz(long,lat,x0,y0,z0)
     CALL spdist(x0,y0,z0,xc,yc,zc,r)
     r = r/r0
     IF (r < 1.0d0) THEN
        q = 0.5d0*h0*(1.0d0 + COS(pi*r))
     ELSE
        q = 0.0d0
     ENDIF
     phiexac(if0) = q
     phi2(if0) = q*farea(if0,ngrids)
     IF (q > qmax) THEN
        qmax = q
        ifmx = if0
     ENDIF
  ENDDO
  PRINT *,'Max phi = ',qmax,' in cell ',ifmx

  ! Cosine bell on dual grid
  !DO iv0 = 1, nv
  !  long = vlong(iv0,ngrids)
  !  lat = vlat(iv0,ngrids)
  !  CALL ll2xyz(long,lat,x0,y0,z0)
  !  CALL spdist(x0,y0,z0,xc,yc,zc,r)
  !  r = r/r0
  !  IF (r < 1.0d0) THEN
  !    q = 0.5d0*h0*(1.0d0 + COS(pi*r))
  !  ELSE
  !    q = 0.0d0
  !  ENDIF
  !  pvexac(iv0) = q
  !  xzeta2(iv0) = q
  !  if (q > qmax) then
  !    qmax = q
  !    ivmx = iv0
  !  endif
  !ENDDO
  !print *,'Max pv = ',qmax,' in dual cell ',ivmx

  ! Construct absolute vorticity
  !zeta2 = xzeta2*xphibar2

  !print *,'max xphibar2 = ',MAXVAL(xphibar2)
  !print *,'Max zeta2 = ',MAXVAL(zeta2)

  dphi2 = phi2/farea(:,ngrids)
  PRINT *,'min and max phi ',minval(dphi2), maxval(dphi2)
  CALL output

  ! Stream function: solid body rotation at angle alpha
  u00 = 2.0d0*pi*rearth/(12.0d0*86400.0d0)
  !alpha = 0.0d0
  alpha = 0.25d0*pi
  !alpha = 0.5*pi
  ca = COS(alpha)
  sa = SIN(alpha)
  DO iv0 = 1, nv
     clat = COS(vlat(iv0,ngrids))
     slat = SIN(vlat(iv0,ngrids))
     slon = SIN(vlong(iv0,ngrids))
     psi(iv0) = u00*rearth*(ca*slat + sa*clat*slon)
  ENDDO

  ! Non-divergent velocity field
  ! U is - D_1 (psi); sign of psi taken care of above
  CALL Dprimal1(psi,ubar1,ngrids,nv,ne)
  CALL massM(ubar1,mu,ngrids,ne)
  ! Corresponding v field
  CALL HodgeHinv(mu,vbar1,ngrids,ne,-20)
  ! Perpendicular components
  CALL operW(ubar1,wu,ngrids,ne)
  CALL massMinv(wu,uperpbar1,ngrids,ne,-20)
  CALL HodgeHinv(wu,vperpbar1,ngrids,ne,-20)

  ! Compute divergence - needed to modify swept areas in
  ! routine primaladvflx
  ! (Should really use old velocity, but since div = 0 here
  ! it shouldn't matter.)
  CALL Dprimal2(ubar1,div2,ngrids,ne,nf)
  divfac = div2/farea(:,ngrids)
  divfac = 1.0d0/(1.0d0 + beta_v*dt*divfac)


  ! Primal grid mass flux
  ! CALL primaladvflx(phi2,mf1)

  ! Dual grid mass flux
  ! CALL operW(mf1,wu,ngrids,ne)
  ! CALL HodgeHinv(wu,mfperp1,ngrids,ne,-20)

  ! Loop over time steps
  DO istep = 1, 12*48

     ! Compute advective fluxes
     CALL primaladvflx(phi2,flx1)

     ! Divergence of flux
     CALL Dprimal2(flx1,dphi2,ngrids,ne,nf)

     ! Update phi2
     phi2 = phi2 - dphi2

     dphi2 = phi2/farea(:,ngrids)
     PRINT *,'Step ',istep,'  Range of phi ',MINVAL(dphi2), MAXVAL(dphi2)
     PRINT *,'Total mass = ',SUM(phi2)

     ! Construct area integral of PV
     !pv = zeta2/phibar

     ! Compute advective fluxes
     !CALL dualadvflx(pv,flx1)

     ! Minus the divergence of flux
     !CALL Ddual2(flx1,dzeta2,ngrids,ne,nv)

     ! Update zeta2
     !zeta2 = zeta2 + dzeta2

     !pv = zeta2/xphibar2
     !print *,'Step ',istep,'  Range of pv ',MINVAL(pv), MAXVAL(pv)
     !print *,'Total pv = ',SUM(zeta2)

     time = time + dt

     ! Compute errors

     ! Centre of bell (initial bell should be at (pi,0))
     lonrot = pi + u00*time/rearth
     xc =  COS(lonrot)
     yc =  SIN(lonrot)*ca
     zc = -SIN(lonrot)*sa
     CALL xyz2ll(xc,yc,zc,long,lat)

     DO if0 = 1, nf
        long = flong(if0,ngrids)
        lat = flat(if0,ngrids)
        CALL ll2xyz(long,lat,x0,y0,z0)
        CALL spdist(x0,y0,z0,xc,yc,zc,r)
        r = r/r0
        IF (r < 1.0d0) THEN
           q = 0.5d0*h0*(1.0d0 + COS(pi*r))
        ELSE
           q = 0.0d0
        ENDIF
        phiexac(if0) = q
        IF (q > qmax) THEN
           qmax = q
           ifmx = if0
        ENDIF
     ENDDO

     ! Cosine bell on dual grid
     !DO iv0 = 1, nv
     !  long = vlong(iv0,ngrids)
     !  lat = vlat(iv0,ngrids)
     !  CALL ll2xyz(long,lat,x0,y0,z0)
     !  CALL spdist(x0,y0,z0,xc,yc,zc,r)
     !  r = r/r0
     !  IF (r < 1.0d0) THEN
     !    q = 0.5d0*h0*(1.0d0 + COS(pi*r))
     !  ELSE
     !    q = 0.0d0
     !  ENDIF
     !  pvexac(iv0) = q
     !  if (q > qmax) then
     !    qmax = q
     !    ivmx = iv0
     !  endif
     !ENDDO

     phierr = dphi2 - phiexac
     l1 = SUM(ABS(phierr)*farea(:,ngrids))/SUM(farea(:,ngrids))
     l2 = SQRT((SUM(phierr*phierr*farea(:,ngrids)))/SUM(farea(:,ngrids)))
     linf = MAXVAL(ABS(phierr))
     PRINT *,'L1   = ',l1
     PRINT *,'L2   = ',l2
     PRINT *,'Linf = ',linf
     WRITE(chanerr,*) istep,l1,l2,linf

     !pverr = pv - pvexac
     !print *,'L1   = ',SUM(ABS(pverr)*varea(:,ngrids))/SUM(varea(:,ngrids))
     !print *,'L2   = ',SQRT((SUM(pverr*pverr*varea(:,ngrids)))/SUM(varea(:,ngrids)))
     !print *,'Linf = ',MAXVAL(ABS(pverr))

     IF (modulo(istep,144) == 0) THEN
        CALL output
     ENDIF

  ENDDO

  ! Compute errors
  PRINT *,'L1   = ',SUM(ABS(phierr)*farea(:,ngrids))/SUM(farea(:,ngrids))
  PRINT *,'L2   = ',SQRT((SUM(phierr*phierr*farea(:,ngrids)))/SUM(farea(:,ngrids)))
  PRINT *,'Linf = ',MAXVAL(ABS(phierr))

  !pverr = pv - pvexac
  !print *,'L1   = ',SUM(ABS(pverr)*varea(:,ngrids))/SUM(varea(:,ngrids))
  !print *,'L2   = ',SQRT((SUM(pverr*pverr*varea(:,ngrids)))/SUM(varea(:,ngrids)))
  !print *,'Linf = ',MAXVAL(ABS(pverr))


  DEALLOCATE(psi, flx1, dphi2, phibar, dzeta2)

END SUBROUTINE testadv

! ====================================================================
!   OUTPUT ROUTINES ADDED BY P. PEIXOTO
!   - Create a map between regular lat-lon grid and the geodesic/cube grid
!   - Save map to disk
!   - Load map is exists, else create
!   - Interpolate scalar fields to regular lon lat files.
!====================================================================

SUBROUTINE setuplltable
  !Setup lat lon table

  USE grid
  USE runtype

  IMPLICIT NONE

  !Auxiliar local variables
  INTEGER:: i, j, k, nlon, nlat, iunit
  REAL*8 :: lon, lat, dlat, dlon, tlon, tlat

  !File name for output
  CHARACTER (LEN=256):: filename
  LOGICAL:: ifile

  !Global table dimensions - if changed, all pre saved files need to be precomputed
  lltable_nlat=720
  lltable_nlon=1440

  nlat=lltable_nlat
  nlon=lltable_nlon

  !Allocate space for table
  ALLOCATE(lltable_primal(nlon, nlat))
  ALLOCATE(lltable_dual(nlon, nlat))

  !Filename for lat lon table
  IF(ygridfile(1:3)=="grd")THEN
     filename="grd/lltable_"//trim(ygridfile(5:len(trim(ygridfile))))
  ELSE
     filename="lltable_"//trim(ygridfile)
  END IF

  !Check if it exists
  INQUIRE(FILE=filename, EXIST=ifile)

  IF(ifile)THEN
     !Load table
     CALL getunit(iunit)
     OPEN(iunit,FILE=filename, STATUS='old',  FORM='unformatted')
     READ(iunit) lltable_primal
     READ(iunit) lltable_dual
     CLOSE(iunit)
     PRINT*, "Table loaded : " , trim(filename)
  ELSE
     !Create table
	 print*, "Creating lat-lon table"
	 print*, "    This might take a while, but is done only once per grid"
	 print*, "    and the table is then saved for future use."
     dlat=180.0/real(nlat, 8)
     dlon=360.0/real(nlon, 8)
     !Pixel registration mode (GMT) (at midpoint of cell)
     tlat=-90+dlat/2.0
     DO j=1,nlat
        tlon=-180+dlon/2.0
        DO i=1,nlon
           lltable_primal(i,j)=getnearfaceind(tlon , tlat)
           lltable_dual(i,j)=getnearvertexind(tlon , tlat, lltable_primal(i,j))
           tlon=tlon+dlon
        END DO
        tlat=tlat+dlat
        if(mod(j, 8)==0)then
        	PRINT*, "Lat: ", tlat
        end if
     END DO

     PRINT*, "Created reference table to lat-long grid"

     PRINT*, "Saving table in file: ", filename
     CALL getunit(iunit)
     OPEN(iunit,FILE=filename, STATUS='replace', FORM='unformatted')
     WRITE(iunit) lltable_primal
     WRITE(iunit) lltable_dual
     CLOSE(iunit)

  END IF

  !call plot_longlat



CONTAINS

  FUNCTION getnearfaceind(lon , lat)
    ! Get the index of the nearest face of the point (lon, lat)

    USE grid
    USE constants

	integer*4::faceind !Nearest face index to this lat lon
    REAL*8:: lon, lat, tlat, tlon
    REAL*8:: d, dnear
    REAL*8:: x, y, z, x0, y0, z0
    INTEGER*4:: i, iv, inear
    INTEGER*4:: getnearfaceind

    !Convert -180-180 long to 0-360 longitudes in radians
    IF(lon<0)THEN
       tlon=(lon+360.0d0)*deg2rad
    ELSE
       tlon=(lon)*deg2rad
    END IF
    tlat=lat*deg2rad
    CALL LL2XYZ(tlon, tlat, x0, y0, z0)

    dnear=10000000000.00
    inear=0
    DO i=1,nface(ngrids)
       CALL LL2XYZ(flong(i, ngrids), flat(i, ngrids), x, y, z)
       CALL SPDIST(X0,Y0,Z0,X,Y,Z,d)
       IF(d<dnear)THEN
          dnear=d
          inear=i
       END IF
    END DO

    getnearfaceind=inear

  END FUNCTION getnearfaceind

  FUNCTION getnearvertexind(lon , lat, iface)
    ! Get the index of the nearest vertex of the point (lon, lat)

    USE grid
    USE constants

    REAL*8:: lon, lat, tlat, tlon
    REAL*8:: d, dnear
    REAL*8:: x, y, z, x0, y0, z0
    INTEGER*4:: i, iv, inear, iface
    INTEGER*4:: getnearvertexind

    !Convert -180-180 long to 0-360 longitudes in radians
    IF(lon<0)THEN
       tlon=(lon+360.0d0)*deg2rad
    ELSE
       tlon=(lon)*deg2rad
    END IF
    tlat=lat*deg2rad
    CALL LL2XYZ(tlon, tlat, x0, y0, z0)

    dnear=10000000000.00
    inear=0
    DO i=1,neoff(iface, ngrids)
		iv=voff(iface, i, ngrids)
       CALL LL2XYZ(vlong(iv, ngrids), vlat(iv, ngrids), x, y, z)
       CALL SPDIST(X0,Y0,Z0,X,Y,Z,d)
       IF(d<dnear)THEN
          dnear=d
          inear=iv
       END IF
    END DO

    getnearvertexind=inear

  END FUNCTION getnearvertexind

END SUBROUTINE setuplltable

SUBROUTINE outputll
  !Output data in latitude longitude table
  ! Standard output fields

  USE runtype
  USE constants
  USE grid
  USE state
  USE errdiag
  USE timestep
  IMPLICIT NONE

  INTEGER :: nf, ne, nv
  REAL*8, ALLOCATABLE :: phi(:), zeta2(:), phiv2(:), pv(:), v1(:), mu(:)
  CHARACTER*12 :: ytime
  CHARACTER*54 :: yname

  ! Construct timestep element of filename
  WRITE(ytime,*) nint(time)
  ytime=trim(adjustl(trim(ytime)))

  ALLOCATE(phi(nfacex), zeta2(nvertx), phiv2(nvertx), pv(nvertx), &
       v1(nedgex), mu(nedgex))

  ! -----------------------------------------------------------------------

  nf = nface(ngrids)
  ne = nedge(ngrids)
  nv = nvert(ngrids)

  ! Surface height
  ! phi = (phi2 + orog2)/(gravity*farea(:,ngrids))
  ! CALL dumpm(phi,'h',nf,'primal')

  ! Surface geopotential
  !phi = (phi2 + orog2)/(farea(:,ngrids))
  phi = (phi2)/(farea(:,ngrids))
  !CALL dumpm(phi,'phi',nf,'primal')
  CALL plot_longlat(phi/gravity,trim(icname)//'_h_t'//trim(ytime),nf,'primal')

  if(ic==2.or.ic==8.or.ic==9)then
    phierr = (phi2 - phi2_init)/farea(:,ngrids)
    CALL plot_longlat(phierr/gravity,trim(icname)//'_herr_t'//trim(ytime),nf,'primal')
  ! CALL dumpm(phierr,'phierr',nf,'primal')
  end if

  ! Vorticity and potential vorticity
  CALL massM(u1,mu,ngrids,ne)
  CALL HodgeHinv(mu,v1,ngrids,ne,-20)
  CALL Ddual2(v1,zeta2,ngrids,ne,nv)
  CALL operR(phi2,pv,ngrids,nf,nv)  ! Borrow pv array temporarily
  CALL HodgeJinv(pv,phiv2,ngrids,nv,-20)
  pv = (zeta2 + planvort2)/phiv2
  zeta2 = zeta2/varea(:,ngrids)
  CALL plot_longlat(zeta2,trim(icname)//'_vort_t'//trim(ytime),nv,'dual  ')
  !CALL plot_longlat(pv,trim(icname)//'_pv_t'//trim(ytime),nv,'dual  ')
  !CALL dumpm(pv,'PV',nv,'dual  ')
  !pv = xzeta2/phiv2
  !CALL dumpm(pv,'xPV',nv,'dual  ')

  !PRINT*, "Lat-lon data dumped in dump folder!"

  DEALLOCATE(phi, zeta2, phiv2, pv, v1, mu)

  ! -----------------------------------------------------------------------

END SUBROUTINE outputll

SUBROUTINE plot_longlat(q, ytitle, npt, ygrid)
  ! Plot primal field variable using lat lon regular grid
  ! P. Peixoto

  USE grid
  USE runtype

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: npt
  REAL*8, INTENT(IN) :: q(npt)
  CHARACTER*(*), INTENT(IN) :: ytitle
  CHARACTER*6, INTENT(IN) :: ygrid

  !Auxiliar local variables
  INTEGER*4:: i, j, k, nlon, nlat, iunit
  REAL*8 :: lon, lat, dlat, dlon, tlon, tlat

  !Buffer for data output (faster to write)
  ! Single precision output
  REAL*4, ALLOCATABLE :: buffer(:,:)

  !File name for output
  CHARACTER (LEN=256):: filename

  nlat=lltable_nlat
  nlon=lltable_nlon
  dlat=180.0/real(nlat, 8)
  dlon=360.0/real(nlon, 8)

  !Buffer for binary plotting - faster
  ALLOCATE(buffer(3, 1:nlat*nlon))

  IF (ygrid == 'primal') THEN
     !Write to buffer
     k=1
     !Pixel registration mode (GMT) (at midpoint of cell)
     tlat=-90+dlat/2.0
     DO j=1,nlat
        tlon=-180+dlon/2.0
        DO i=1,nlon
           buffer(1, k)=tlon
           buffer(2, k)=tlat
           buffer(3, k)=q(lltable_primal(i,j))
           k=k+1
           !write(iunit) tlon, tlat, varunif(i,j)
           tlon=tlon+dlon
        END DO
        tlat=tlat+dlat
     END DO
  ELSEIF (ygrid == 'dual  ') THEN
     !Write to buffer
     k=1
     !Pixel registration mode (GMT) (at midpoint of cell)
     tlat=-90+dlat/2.0
     DO j=1,nlat
        tlon=-180+dlon/2.0
        DO i=1,nlon
           buffer(1, k)=tlon
           buffer(2, k)=tlat
           buffer(3, k)=q(lltable_dual(i,j))
           k=k+1
           !write(iunit) tlon, tlat, varunif(i,j)
           tlon=tlon+dlon
        END DO
        tlat=tlat+dlat
     END DO
  ELSE
     PRINT *,'Invalid option: ygrid = ',ygrid,' in subroutine plot_longlat'
     STOP
  ENDIF

  !Filename for lat lon field
  IF(ygridfile(1:3)=="grd")THEN
     filename="dump/"//trim(ytitle)//"_"//trim(ygridfile(5+12:len(trim(ygridfile))))
  ELSE
     filename=trim(ytitle)//"_"//trim(ygridfile(12:len(trim(ygridfile))))
  END IF

  !Write values on file
  CALL getunit(iunit)
  OPEN(iunit,FILE=filename, STATUS='replace', ACCESS='stream', FORM='unformatted')
  !Write whole block to file (much faster)
  WRITE(iunit) buffer
  CLOSE(iunit)
  PRINT*, "Long-lat output written in : ", trim(filename)

  RETURN
END SUBROUTINE plot_longlat



SUBROUTINE getunit ( iunit )
  !----------------------------------------------------------
  ! GETUNIT returns a free FORTRAN unit number.
  !
  !    A "free" FORTRAN unit number is an integer between 1 and 99 which
  !    is not currently associated with an I/O device.  A free FORTRAN unit
  !    number is needed in order to open a file with the OPEN command.
  !
  !    If IUNIT = 0, then no free FORTRAN unit could be found, although
  !    all 99 units were checked (except for units 5, 6 and 9, which
  !    are commonly reserved for console I/O).
  !
  !    Otherwise, IUNIT is an integer between 1 and 99, representing a
  !    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
  !    are special, and will never return those values.
  !
  !    John Burkardt
  !    18 September 2005
  !----------------------------------------------------------------------------
  INTEGER :: i
  INTEGER :: ios
  INTEGER :: iunit
  LOGICAL:: lopen

  iunit = 0
  DO i = 11, 99
     IF ( i /= 5 .AND. i /= 6 .AND. i /= 9 ) THEN
        INQUIRE ( UNIT = i, OPENED = lopen, IOSTAT = ios )
        IF ( ios == 0 ) THEN
           IF ( .NOT. lopen ) THEN
              iunit = i
              RETURN
           END IF
        END IF
     END IF
  END DO

  RETURN
END SUBROUTINE getunit
