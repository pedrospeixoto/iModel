MODULE channels

! Tidy list of all I/O channels in one place to avoid accidental
! overuse of any channel number


INTEGER, PARAMETER :: changin = 25          ! Input grid information
INTEGER, PARAMETER :: changout = 26         ! Output grid and operator information


END MODULE channels

! ===============================================================

PROGRAM buildop

! To take basic grid information from grid generator program
! and construct mass matrices and other operators needed for
! a mimetic primal-dual mixed finite element shallow water model.
!
! The algrithm takes the primal and dual meshes and constructs
! a supermesh of triangles. The basis elements on the primal and
! dual meshes are then constructed by gluing together the lowest
! order elements on the supermesh.
!
! Most of the code is applicable to arbitrary polygonal grids.
! The multigrid restriction operator is currently only coded
! for a hexagonal icosahedral grid and a cubed sphere grid.
!
!
! John Thuburn July/August/September/October 2012
!
! ---------------------------------------------------------------

! Don't need this
USE grid

IMPLICIT NONE

INTEGER :: igrid

! ---------------------------------------------------------------


! Read input grid data
CALL readgrid
PRINT *,'Input data read'

! Construct some basic geometrical information
CALL buildgeom
PRINT *,'Done buildgeom'

! Information is assumed to be in a certain order in arrays
! defining connectivity. This may already hold for the output
! of the grid generator, but check and correct if necessary.
CALL ordering
PRINT *,'Done ordering'


! Loop over hierarchy of grids
DO igrid = 1, ngrids

  PRINT *,' '
  PRINT *,'Grid ',igrid

  ! Build compound elements on grid igrid
  CALL buildcomp(igrid)
  PRINT *,'Done buildcomp'

  ! Build matrix operators on grid igrid
  CALL buildmat(igrid)
  PRINT *,'Done buildmat'

ENDDO
PRINT *,' '

! Build injection and restriction operator
CALL buildInjectionAndRestriction
PRINT *,'Multigrid restriction operator built'

! Write everything out
CALL writegrid
PRINT *,'Data written'


! ---------------------------------------------------------------

END PROGRAM buildop

! ===============================================================

SUBROUTINE ordering

! Some of the later routines assume the grid connectivity data
! obeys certain ordering conventions. Some of these may be
! imposed by the grid generation program, but check and fix
! here if necessary.

USE grid

IMPLICIT NONE

INTEGER :: igrid, if0, ne, ix1, if1, ixmin, ifmin, ix2, if2, &
           if3, iv0, ie1, if21, if22, iv11, iv12, iv21, iv22, &
           ie2, ie3, iemin, if11, if12, ie0, iv1, iv2
REAL*8 :: pi, long, lat, x0(3), x1(3), d1x, d1y, d1z, thetamin, &
          x2(3), d2x, d2y, d2z, cs, sn, theta
LOGICAL :: lfound

! ---------------------------------------------------------------

pi = 4.0d0*ATAN(1.0D0)


! Loop over all grids in the hierarchy
DO igrid = 1, ngrids

  ! Loop over all faces on this grid
  DO if0 = 1, nface(igrid)

    ! Coordinates of face if0
    long = flong(if0,igrid)
    lat = flat(if0,igrid)
    CALL ll2xyz(long,lat,x0)

    ! Number of edges/neighbours
    ne = neoff(if0,igrid)

    ! Sort FNXTF into anticlockwise order 
    DO ix1 = 1, ne - 2

      ! Coordinates of IX1'th neighbour
      if1 = fnxtf(if0,ix1,igrid)
      long = flong(if1,igrid)
      lat = flat(if1,igrid)
      CALL ll2xyz(long,lat,x1)
      d1x = x1(1) - x0(1)
      d1y = x1(2) - x0(2)
      d1z = x1(3) - x0(3)

      ! Find next neighbour (anticlockwise)
      thetamin = pi
      ixmin = 0
      ifmin = 0

      DO ix2 = ix1 + 1, ne

        ! Coordinates of IX2'th neighbour
        if2 = fnxtf(if0,ix2,igrid)
        long = flong(if2,igrid)
        lat = flat(if2,igrid)
        CALL ll2xyz(long,lat,x2)
        d2x=x2(1) - x0(1)
        d2y=x2(2) - x0(2)
        d2z=x2(3) - x0(3)
        cs = d1x*d2x + d1y*d2y + d1z*d2z
        sn = x0(1)*(d1y*d2z - d1z*d2y) &
           + x0(2)*(d1z*d2x - d1x*d2z) &
           + x0(3)*(d1x*d2y - d1y*d2x)
        theta = ATAN2(sn,cs)
        IF ((theta < thetamin) .AND. (theta > 0.0d0)) THEN
          ixmin = ix2
          ifmin = if2
          thetamin = theta
        ENDIF
      ENDDO

!     The face in position IXMIN belongs in position IX1+1 so swap them
      if3 = fnxtf(if0,ix1+1,igrid)
      fnxtf(if0,ix1+1,igrid) = ifmin
      fnxtf(if0,ixmin,igrid) = if3

    ENDDO

    ! Now sort EOFF to correspond to FNXTF
    DO ix1 = 1, ne

      if1 = fnxtf(if0,ix1,igrid)
      ix2 = ix1 - 1
      lfound = .FALSE.
      DO WHILE (.NOT. lfound)
         ix2 = ix2 + 1
         ie1 = eoff(ix2,if0,igrid)
         if21 = fnxte(1,ie1,igrid)
         if22 = fnxte(2,ie1,igrid)
         IF ((if21 + if22) == (if0 + if1)) lfound = .TRUE.
      ENDDO

      ! Edge IE2 corresponds to face IF1
      eoff(ix2,if0,igrid) = eoff(ix1,if0,igrid)
      eoff(ix1,if0,igrid) = ie1

    ENDDO

    ! Order VOFF so that the k'th vertex is between the
    ! k'th and (k+1)'th edges in EOFF
    DO ix1 = 1, ne
      ix2 = ix1 + 1
      IF (ix2 > ne) ix2 = 1
      ie1 = eoff(ix1,if0,igrid)
      ie2 = eoff(ix2,if0,igrid)
      ! Find the common vertex of IE1 and IE2
      iv11 = vofe(1,ie1,igrid)
      iv12 = vofe(2,ie1,igrid)
      iv21 = vofe(1,ie2,igrid)
      iv22 = vofe(2,ie2,igrid)
      IF ((iv11 == iv21) .OR. (iv11 == iv22)) THEN
        iv0 = iv11
      ELSEIF ((iv12 == iv21) .OR. (iv12 == iv22)) THEN
        iv0 = iv12
      ELSE
        PRINT *,'Common vertex not found'
        STOP
      ENDIF
      voff(if0,ix1,igrid) = iv0
    ENDDO

  ENDDO

  ! Loop over all vertices on this grid
  DO iv0 = 1, nvert(igrid)

    ! Coordinates of vertex iv0
    long = vlong(iv0,igrid)
    lat = vlat(iv0,igrid)
    CALL ll2xyz(long,lat,x0)

    ! Number of edges / adjacent faces
    ne = neofv(iv0,igrid)

    ! Sort EOFV into anticlockwise order
    DO ix1 = 1, ne - 2

      ! Coordinates of IX1'th edge
      ie1 = eofv(ix1,iv0,igrid)
      long = elong(ie1,igrid)
      lat = elat(ie1,igrid)
      CALL ll2xyz(long,lat,x1)
      d1x = x1(1) - x0(1)
      d1y = x1(2) - x0(2)
      d1z = x1(3) - x0(3)

!     Find next neighbour (anticlockwise)
      thetamin = pi
      ixmin = 0
      ifmin = 0

      DO ix2 = ix1 + 1, ne

        ! Coordinates of IX2'th neighbour
        ie2 = eofv(ix2,iv0,igrid)
        long = elong(ie2,igrid)
        lat = elat(ie2,igrid)
        CALL ll2xyz(long,lat,x2)
        d2x=x2(1) - x0(1)
        d2y=x2(2) - x0(2)
        d2z=x2(3) - x0(3)
        cs = d1x*d2x + d1y*d2y + d1z*d2z
        sn = x0(1)*(d1y*d2z - d1z*d2y) &
           + x0(2)*(d1z*d2x - d1x*d2z) &
           + x0(3)*(d1x*d2y - d1y*d2x)
        theta = ATAN2(sn,cs)
        IF ((theta < thetamin) .AND. (theta > 0.0d0)) THEN
          ixmin = ix2
          iemin = ie2
          thetamin = theta
        ENDIF
      ENDDO

!     The edge in position IXMIN belongs in position IX1+1 so swap them
      ie3 = eofv(ix1+1,iv0,igrid)
      eofv(ix1+1,iv0,igrid) = iemin
      eofv(ixmin,iv0,igrid) = ie3

    ENDDO

    ! Order FOFV so that the k'th face is between the
    ! k'th and (k+1)'th edges in EOFV
    DO ix1 = 1, ne
      ix2 = ix1 + 1
      IF (ix2 > ne) ix2 = 1
      ie1 = eofv(ix1,iv0,igrid)
      ie2 = eofv(ix2,iv0,igrid)
      ! Find the common face of IE1 and IE2
      if11 = fnxte(1,ie1,igrid)
      if12 = fnxte(2,ie1,igrid)
      if21 = fnxte(1,ie2,igrid)
      if22 = fnxte(2,ie2,igrid)
      IF ((if11 == if21) .OR. (if11 == if22)) THEN
        if0 = if11
      ELSEIF ((if12 == if21) .OR. (if12 == if22)) THEN
        if0 = if12
      ELSE
        PRINT *,'Common face not found'
        PRINT *,'Grid ',igrid,' vertex ',iv0
        STOP
      ENDIF
      fofv(iv0,ix1,igrid) = if0
    ENDDO

  ENDDO


  ! Loop over all edges on this grid
  DO ie0 = 1, nedge(igrid)

    ! Sort VOFE so that VOFE(1) -> VOFE(2) (tangent vector)
    ! is 90 degrees anticlockwise of FNXTE(1) -> FNXTE(2) (normal vector)
    if1 = fnxte(1,ie0,igrid)
    if2 = fnxte(2,ie0,igrid)
    long = flong(if1,igrid)
    lat = flat(if1,igrid)
    CALL ll2xyz(long,lat,x0)
    long = flong(if2,igrid)
    lat = flat(if2,igrid)
    CALL ll2xyz(long,lat,x1)
    d1x = x1(1) - x0(1)
    d1y = x1(2) - x0(2)
    d1z = x1(3) - x0(3)
    iv1 = vofe(1,ie0,igrid)
    iv2 = vofe(2,ie0,igrid)
    long = vlong(iv1,igrid)
    lat = vlat(iv1,igrid)
    CALL ll2xyz(long,lat,x0)
    long = vlong(iv2,igrid)
    lat = vlat(iv2,igrid)
    CALL ll2xyz(long,lat,x1)
    d2x = x1(1) - x0(1)
    d2y = x1(2) - x0(2)
    d2z = x1(3) - x0(3)
    sn = x0(1)*(d1y*d2z - d1z*d2y) &
       + x0(2)*(d1z*d2x - d1x*d2z) &
       + x0(3)*(d1x*d2y - d1y*d2x)
    IF (sn < 0.0d0) THEN
      ! Swap the two vertices
      vofe(1,ie0,igrid) = iv2
      vofe(2,ie0,igrid) = iv1
    ENDIF

  ENDDO

ENDDO


! ---------------------------------------------------------------

END SUBROUTINE ordering

! ===============================================================

SUBROUTINE buildgeom

! Build some basic geometrical information:
! locations of edge crossing points;
! lengths of primal and dual edges;
! cell areas.

USE grid

IMPLICIT NONE
INTEGER :: igrid, ie0, if1, if2, iv1, iv2, if0, ie1, ix
REAL*8 :: long, lat, x1(3), x2(3), y1(3), y2(3), &
          n1(3), n2(3), r1(3), mag, l1sq, l2sq, l3sq, &
          area, a

! ---------------------------------------------------------------

! Loop over grids
DO igrid = 1, ngrids

  ! Loop over edges
  DO ie0 = 1, nedge(igrid)

    ! Locate cell centres either side (i.e dual vertices)
    if1 = fnxte(1,ie0,igrid)
    long = flong(if1,igrid)
    lat = flat(if1,igrid)
    CALL ll2xyz(long,lat,x1)
    if2 = fnxte(2,ie0,igrid)
    long = flong(if2,igrid)
    lat = flat(if2,igrid)
    CALL ll2xyz(long,lat,x2)

    ! And ends of edge (i.e. primal vertices)
    iv1 = vofe(1,ie0,igrid)
    long = vlong(iv1,igrid)
    lat = vlat(iv1,igrid)
    CALL ll2xyz(long,lat,y1)
    iv2 = vofe(2,ie0,igrid)
    long = vlong(iv2,igrid)
    lat = vlat(iv2,igrid)
    CALL ll2xyz(long,lat,y2)

    ! Normal to plane of dual edge
    n1(1) = x1(2)*x2(3) - x1(3)*x2(2)
    n1(2) = x1(3)*x2(1) - x1(1)*x2(3)
    n1(3) = x1(1)*x2(2) - x1(2)*x2(1)
    ! Normal to plane of primal edge
    n2(1) = y1(2)*y2(3) - y1(3)*y2(2)
    n2(2) = y1(3)*y2(1) - y1(1)*y2(3)
    n2(3) = y1(1)*y2(2) - y1(2)*y2(1)
    ! Hence radial vector of crossing point
    r1(1) = n1(2)*n2(3) - n1(3)*n2(2)
    r1(2) = n1(3)*n2(1) - n1(1)*n2(3)
    r1(3) = n1(1)*n2(2) - n1(2)*n2(1)
    ! Normalize to unit sphere
    mag = SQRT(r1(1)*r1(1) + r1(2)*r1(2) + r1(3)*r1(3))
    r1 = r1/mag
    ! Convert back to lat long and save
    CALL xyz2ll(r1,long,lat)
    elong(ie0,igrid) = long
    elat(ie0,igrid) = lat

    ! Dual edge is now composed of two straight edge pieces
    l1sq = (r1(1) - x1(1))**2 + (r1(2) - x1(2))**2 + (r1(3) - x1(3))**2
    l2sq = (r1(1) - x2(1))**2 + (r1(2) - x2(2))**2 + (r1(3) - x2(3))**2
    ddist(ie0,igrid) = SQRT(l1sq) + SQRT(l2sq)

    ! Primal edge is now composed of two straight edge pieces
    l1sq = (r1(1) - y1(1))**2 + (r1(2) - y1(2))**2 + (r1(3) - y1(3))**2
    l2sq = (r1(1) - y2(1))**2 + (r1(2) - y2(2))**2 + (r1(3) - y2(3))**2
    ldist(ie0,igrid) = SQRT(l1sq) + SQRT(l2sq)

  ENDDO

  ! Intialize dual cell areas to zero and
  ! accumulate as we touch each vertex
  varea(:,igrid) = 0.0d0
  ! Loop over primal cells
  DO if1 = 1, nface(igrid)
    IF (SFCIndexAvailable == 0) THEN
      if0 = if1
    ELSE
      if0 = fNewFaceIdInverse(if1, igrid)
    END IF

    ! Initialize area to zero and locate call centre
    area = 0.0d0
    long = flong(if0,igrid)
    lat = flat(if0,igrid)
    CALL ll2xyz(long,lat,r1)
    ! Accumulate contributions to cell area from supermesh cells
    ! Loop over primal cell edges
    DO ix = 1, neoff(if0,igrid)
      ! There are two supermesh cells incident on each
      ! primal cell edge
      ie1 = eoff(ix,if0,igrid)
      long = elong(ie1,igrid)
      lat = elat(ie1,igrid)
      CALL ll2xyz(long,lat,x1)
      ! First supermesh cell
      iv1 = vofe(1,ie1,igrid)
      long = vlong(iv1,igrid)
      lat = vlat(iv1,igrid)
      CALL ll2xyz(long,lat,x2)
      CALL triangle(r1,x1,x2,l1sq,l2sq,l3sq,a)
      area = area + a
      varea(iv1,igrid) = varea(iv1,igrid) + a
      ! Second supermesh cell
      iv1 = vofe(2,ie1,igrid)
      long = vlong(iv1,igrid)
      lat = vlat(iv1,igrid)
      CALL ll2xyz(long,lat,x2)
      CALL triangle(r1,x1,x2,l1sq,l2sq,l3sq,a)
      area = area + a
      varea(iv1,igrid) = varea(iv1,igrid) + a
    ENDDO
    farea(if0,igrid) = area

  ENDDO

ENDDO


! Construct the tables eoffin and eofvin
ALLOCATE(eoffin(nefmx,nfacex,ngrids), eofvin(nevmx,nvertx,ngrids))
DO igrid = 1, ngrids

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

ENDDO



! ---------------------------------------------------------------

END SUBROUTINE buildgeom

! ===============================================================

SUBROUTINE buildcomp(igrid)

! Construct the information about `internal' degrees of freedom
! needed to define the compound elements on grid igrid.

IMPLICIT NONE
INTEGER, INTENT(IN) :: igrid

! ---------------------------------------------------------------

! Compound P0 elements on primal cells (space Vp)
CALL buildvp(igrid)
PRINT *,'  Done buildvp'

! Compound RT0 elements on primal cells (space Sp)
CALL buildsp(igrid)
PRINT *,'  Done buildsp'

! Compound P1 elements on primal cells (space Ep)
CALL buildep(igrid)
PRINT *,'  Done buildep'

! Compound P0 elements on dual cells (space Vd)
CALL buildvd(igrid)
PRINT *,'  Done buildvd'

! Compound N0 elements on dual cells (space Sd)
CALL buildsd(igrid)
PRINT *,'  Done buildsd'

! ---------------------------------------------------------------

END SUBROUTINE buildcomp

! ===============================================================

SUBROUTINE buildmat(igrid)

! Build the matrix representations, in terms of stencils
! and ceofficients, of the various operators and mass matrices
! needed.

IMPLICIT NONE
INTEGER, INTENT(IN) :: igrid

! ---------------------------------------------------------------

! Mass matrix L
CALL buildL(igrid)
PRINT *,'  Done buildL'

! Mass matrix M
CALL buildM(igrid)
PRINT *,'  Done buildM'

! Vp to Ep transfer operator R
CALL buildR(igrid)
PRINT *,'  Done buildR'

! W operator for constructing perp of vectors
! (maps E_p tp E_p)
CALL buildW(igrid)
PRINT *,'  Done buildW'

! Vd to Ep transfer operator J
CALL buildJ(igrid)
PRINT *,'  Done buildJ'

! Sd to Sp transfer operator H
CALL buildH(igrid)
PRINT *,'  Done buildH'

! ---------------------------------------------------------------

END SUBROUTINE buildmat

! ===============================================================

SUBROUTINE buildvp(igrid)

! Find the internal degrees of freedom that define the
! compound P0 elements on primal cells (space Vp).
!
! The compound element is constant over each cell with value
! equal to the inverse of the cell area. The micro elements
! are themselves constant over each supermesh cell with value
! equal to the inverse of the supermesh cell area. The coefficient
! of the micro element is thus the supermesh cell area divided by
! the compound cell area.
!
! There are 2*neoff micro elements in one compound element,
! and they are assumed to be ordered anticlockwise, starting with
! those incident on the first edge of the compound element.

USE grid

IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid
INTEGER :: nf, if0, ie1, iv1, ix, ixm
REAL*8 :: long, lat, x0(3), x1(3), x2(3), l1sq, l2sq, l3sq, a

! ---------------------------------------------------------------

! Number of faces on this grid
nf = nface(igrid)

! Number of coefficients to define element - depends on
! shape of primal cell.
ncvp(1:nf) = 2*neoff(:nf,igrid)


! Loop over primal cells
DO if0 = 1, nf

  ! Cell centre
  long = flong(if0,igrid)
  lat = flat(if0,igrid)
  CALL ll2xyz(long,lat,x0)

  ! Loop over edges of cell
  DO ix = 1, neoff(if0,igrid)
    ixm = ix - 1
    IF (ixm < 1) ixm = neoff(if0,igrid)
    ie1 = eoff(ix,if0,igrid)

    ! Edge crossing point
    long = elong(ie1,igrid)
    lat = elat(ie1,igrid)
    CALL ll2xyz(long,lat,x1)

    ! First micro element 
    iv1 = voff(if0,ixm,igrid)
    long = vlong(iv1,igrid)
    lat = vlat(iv1,igrid)
    CALL ll2xyz(long,lat,x2)
    CALL triangle(x0,x1,x2,l1sq,l2sq,l3sq,a)
    cvp(if0,2*ix - 1) = a/farea(if0,igrid)

    ! Second micro element
    iv1 = voff(if0,ix,igrid)
    long = vlong(iv1,igrid)
    lat = vlat(iv1,igrid)
    CALL ll2xyz(long,lat,x2)
    CALL triangle(x0,x1,x2,l1sq,l2sq,l3sq,a)
    cvp(if0,2*ix) = a/farea(if0,igrid)

  ENDDO

ENDDO


! ---------------------------------------------------------------

END SUBROUTINE buildvp

! ===============================================================

SUBROUTINE buildsp(igrid)

! Find the internal degrees of freedom that define the
! compound RT0 elements on primal cells (space Sp).
!
! Each compound element is associated with a unit normal flux
! at some primal edge. The first two degrees of freedom are
! the fluxes at the two micro edges that make up the primal edge.
! The sign convention is the same as the sign convention for
! the primal edge.
! The remaining degrees of freedom are the micro edge fluxes
! internal to the primal cells either side of the primal edge.
! These dofs are assumed to be ordered anticlockwise starting
! with the microedge that touches the middle of the primal edge,
! for the first primal cell incident on the primal edge, then for
! the second primal cell incident on the primal edge. The sign convention
! is that clockwise circulation about a primal cell centre is positive.
!
! The compound element has constant divergence in the primal cell
! each side of the primal edge. This constraint fixes all but two
! degrees of freedom. These last two are fixed by demanding that
! a certain measure of the vorticity in each compound cell should
! vanish (this is the usual mixed finite element vorticity defined
! by integration by parts against a test function on the supermesh).

USE grid

IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid
INTEGER :: ne, ie0, if1, if2, ixf, ixv, iv1, nef1, ixc, &
           ie1, ixe, ix, ixv1, ixv2, ixe1, offset(2), &
           i1, i2, ixcp
REAL*8 :: long, lat, x0(3), x1(3), l1, sgnflx, &
          l1sq, l2sq, l3sq, a, atot, x2(3), x3(3), &
          sum1, sum2, c00, sg, contrib

! ---------------------------------------------------------------

! Number of edges on this grid
ne = nedge(igrid)

! Number of coefficients to define element - depends on
! shape of cells either side of edge. Computed below.


! Loop over edges
DO ie0 = 1, ne

  ! Find number of coefficients for this edge
  if1 = fnxte(1,ie0,igrid)
  if2 = fnxte(2,ie0,igrid)
  ncsp(ie0) = 2*(neoff(if1,igrid) + neoff(if2,igrid) + 1)

  ! Save related offsets
  offset(1) = 2
  offset(2) = 2 + 2*neoff(if1,igrid)

  ! Find coordinates of edge ie0
  long = elong(ie0,igrid)
  lat = elat(ie0,igrid)
  CALL ll2xyz(long,lat,x0)

  ! Find the lengths of the two micro edges that make up
  ! edge ie0, and hence determine the first two coefficients
  DO ixv = 1, 2
    iv1 = vofe(ixv,ie0,igrid)
    long = vlong(iv1,igrid)
    lat = vlat(iv1,igrid)
    CALL ll2xyz(long,lat,x1)
    l1 = SQRT((x1(1) - x0(1))**2 + (x1(2) - x0(2))**2 + (x1(3) - x0(3))**2)
    csp(ie0,ixv) = l1/ldist(ie0,igrid)
  ENDDO

  ! Loop over faces next to edge ie0
  DO ixf = 1, 2
    if1 = fnxte(ixf,ie0,igrid)

    ! How many edges does this cell have
    nef1 = neoff(if1,igrid)

    ! Is edge flux in or out of cell if1 ?
    IF (ixf == 1) THEN
      sgnflx =  1.0d0 ! outward
    ELSE
      sgnflx = -1.0d0 ! inward
    ENDIF

    ! Area of cell if1
    atot = farea(if1,igrid)

    ! Find coordinates of centre of face if1
    long = flong(if1,igrid)
    lat = flat(if1,igrid)
    CALL ll2xyz(long,lat,x1)

    ! Which edge of cell if1 is ie0 ?
    ixe = 1
    DO WHILE (ie0 .NE. eoff(ixe,if1,igrid))
      ixe = ixe + 1
    ENDDO

    ! Reset pointer to dofs
    ixc = offset(ixf)

    ! Internal fluxes
    ! First guess to get the divergence right
    ! First two internal fluxes
    ixc = ixc + 1
    csp(ie0,ixc) = 0.0d0
    ixv = 3 - ixf
    iv1 = vofe(ixv,ie0,igrid)
    long = vlong(iv1,igrid)
    lat = vlat(iv1,igrid)
    CALL ll2xyz(long,lat,x2)
    CALL triangle(x1,x0,x2,l1sq,l2sq,l3sq,a)
    ixc = ixc + 1
    csp(ie0,ixc) = sgnflx*(csp(ie0,ixv) - a/atot)

    ! Loop to find remaining internal fluxes
    ixe1 = ixe
    DO ix = 1, nef1 - 1
      ixv1 = ixe1
      ixe1 = ixe1 + 1
      IF (ixe1 > nef1) ixe1 = 1
      ixv2 = ixe1
      ie1 = eoff(ixe1,if1,igrid)
      long = elong(ie1,igrid)
      lat = elat(ie1,igrid)
      CALL ll2xyz(long,lat,x3)
      iv1 = voff(if1,ixv1,igrid)
      long = vlong(iv1,igrid)
      lat = vlat(iv1,igrid)
      CALL ll2xyz(long,lat,x2)
      CALL triangle(x1,x2,x3,l1sq,l2sq,l3sq,a)
      ixc = ixc + 1
      csp(ie0,ixc) = csp(ie0,ixc - 1) - sgnflx*a/atot
      iv1 = voff(if1,ixv2,igrid)
      long = vlong(iv1,igrid)
      lat = vlat(iv1,igrid)
      CALL ll2xyz(long,lat,x2)
      CALL triangle(x1,x2,x3,l1sq,l2sq,l3sq,a)
      ixc = ixc + 1
      csp(ie0,ixc) = csp(ie0,ixc - 1) - sgnflx*a/atot
    ENDDO

    ! Now correct to make vorticity vanish
    ! Accumulate integrals over micro elements
    sum1 = 0.0d0
    sum2 = 0.0d0
    ! First the contributions from the flux across edge ie0
    sg = 1.0d0
    DO ixv = 1, 2
      iv1 = vofe(ixv,ie0,igrid)
      long = vlong(iv1,igrid)
      lat = vlat(iv1,igrid)
      CALL ll2xyz(long,lat,x2)
      CALL triangle(x1,x0,x2,l1sq,l2sq,l3sq,a)
      sum2 = sum2 + sg*csp(ie0,ixv)*(l2sq - l3sq)/(12.0d0*a)
      sg = -sg
    ENDDO
    ! Now the contributions from all the internal fluxes
    ! Loop over vertices with two micro elements on each vertex
    ixe1 = ixe
    ixv1 = ixe
    ixc = offset(ixf)
    DO ix = 1, nef1
      
      ! Find the vertex
      iv1 = voff(if1,ixv1,igrid)
      long = vlong(iv1,igrid)
      lat = vlat(iv1,igrid)
      CALL ll2xyz(long,lat,x2)
      ! Find the preceding edge
      ie1 = eoff(ixe1,if1,igrid)
      long = elong(ie1,igrid)
      lat = elat(ie1,igrid)
      CALL ll2xyz(long,lat,x3)
      CALL triangle(x1,x2,x3,l1sq,l2sq,l3sq,a)
      ixc = ixc + 1
      ixcp = ixc + 1
      sum1 = sum1 + l1sq/(4.0d0*a)
      contrib = ( csp(ie0,ixc )*(3.0d0*l1sq + l3sq - l2sq)   &
                + csp(ie0,ixcp)*(3.0d0*l1sq + l2sq - l3sq) ) &
                                               /(24.0d0*a)
      sum2 = sum2 + contrib
      ! Find the following edge
      ixe1 = ixe1 + 1
      IF (ixe1 > nef1) ixe1 = 1
      ie1 = eoff(ixe1,if1,igrid)
      long = elong(ie1,igrid)
      lat = elat(ie1,igrid)
      CALL ll2xyz(long,lat,x3)
      CALL triangle(x1,x2,x3,l1sq,l2sq,l3sq,a)
      ixc = ixc + 1
      ixcp = ixc + 1
      IF (ix == nef1) ixcp = offset(ixf) + 1
      sum1 = sum1 + l1sq/(4.0d0*a)
      contrib = ( csp(ie0,ixcp)*(3.0d0*l1sq + l3sq - l2sq)   &
                + csp(ie0,ixc )*(3.0d0*l1sq + l2sq - l3sq) ) &
                                               /(24.0d0*a)
      sum2 = sum2 + contrib

      ixv1 = ixv1 + 1
      IF (ixv1 > nef1) ixv1 = 1
    ENDDO

    ! Correction to remove vorticity
    c00 = - sum2/sum1
    i1 = offset(ixf) + 1
    i2 = offset(ixf) + 2*nef1
    csp(ie0,i1:i2) = csp(ie0,i1:i2) + c00

  ENDDO

ENDDO


! ---------------------------------------------------------------

END SUBROUTINE buildsp

! ===============================================================

SUBROUTINE buildep(igrid)

! Find the internal degrees of freedom that define the
! compound P1 elements on primal cells (space Ep).
!
! The compound element is piecewise linear. It takes the value
! 1 at the central primal grid vertex, zero at all other
! primal grid vertices, and varies linearly along primal edges.
! To define the compound element we need to compute
! (i) the values at the supermesh vertices adjacent to the
! vertex with value 1, and
! (ii) the values at the supermesh vertices at the centres of
! of the primal cells incident on the central primal vertex.
! Thus there are 2*neofv degrees of freedom to be determined.
!
! The degrees of freedom (i) are easily obtained by linear
! interpolation. The degrees of freedom (ii) are fixed by
! demanding that k X grad of the compound element has zero
! vorticity within each face, where the vorticity is defined
! integration by parts, as for the Sp elements.

USE grid

IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid
INTEGER :: nv, iv0, ix, ixp, ie1, if1, ixv, ixvp, iv1, iv2, &
           nef1, nev0
REAL*8 :: long, lat, x0(3), x1(3), x2(3), x3(3), l1, &
          l1sq, l2sq, l3sq, a, sum1, sum2

! ---------------------------------------------------------------

! Number of vertices on this grid
nv = nvert(igrid)

! Number of coefficients to define element - depends on
! shape of primal cell.
ncep(1:nv) = 2*neofv(:nv,igrid)


! Loop over primal vertices
DO iv0 = 1, nv

  ! Locate the vertex
  long = vlong(iv0,igrid)
  lat = vlat(iv0,igrid)
  CALL ll2xyz(long,lat,x0)

  ! Loop over edges around this vertex and set the
  ! first neofv degres of freedom
  nev0 = neofv(iv0,igrid)
  DO ix = 1, nev0

    ! Locate the supermesh vertex on this edge
    ie1 = eofv(ix,iv0,igrid)
    long = elong(ie1,igrid)
    lat = elat(ie1,igrid)
    CALL ll2xyz(long,lat,x1)
    ! calculate distance of supermesh vertex from central vertex
    ! and hence the coefficient
    l1 = SQRT((x1(1) - x0(1))**2 + (x1(2) - x0(2))**2 + (x1(3) - x0(3))**2)
    cep(iv0,ix) = 1.0d0 - l1/ldist(ie1,igrid)

  ENDDO

  ! Now loop over the faces around this vertex and build the compound
  ! element in that face
  DO ix = 1, nev0

    ! Locate face centre
    if1 = fofv(iv0,ix,igrid)
    long = flong(if1,igrid)
    lat = flat(if1,igrid)
    CALL ll2xyz(long,lat,x1)

    ! Two edges of this face (four supermesh edges/cells)
    ! contribute to an integral
    ! First edge
    ie1 = eofv(ix,iv0,igrid)
    long = elong(ie1,igrid)
    lat = elat(ie1,igrid)
    CALL ll2xyz(long,lat,x2)
    ! Find the vertex of this edge that isn't iv0
    iv2 = vofe(1,ie1,igrid) + vofe(2,ie1,igrid) - iv0
    long = vlong(iv2,igrid)
    lat = vlat(iv2,igrid)
    CALL ll2xyz(long,lat,x3)
    CALL triangle(x2,x3,x1,l1sq,l2sq,l3sq,a)
    sum2 = cep(iv0,ix)*(l1sq - l2sq + l3sq)/(8.0d0*a)
    CALL triangle(x0,x2,x1,l1sq,l2sq,l3sq,a)
    sum2 = sum2 + (l1sq - l2sq + l3sq)/(8.0d0*a) &
                + cep(iv0,ix)*(-l1sq + l2sq + l3sq)/(8.0d0*a)
    ! Second edge
    ixp = ix + 1
    IF (ixp > nev0) ixp = 1
    ie1 = eofv(ixp,iv0,igrid)
    long = elong(ie1,igrid)
    lat = elat(ie1,igrid)
    CALL ll2xyz(long,lat,x2)
    ! Find the vertex of this edge that isn't iv0
    iv2 = vofe(1,ie1,igrid) + vofe(2,ie1,igrid) - iv0
    long = vlong(iv2,igrid)
    lat = vlat(iv2,igrid)
    CALL ll2xyz(long,lat,x3)
    CALL triangle(x3,x2,x1,l1sq,l2sq,l3sq,a)
    sum2 = sum2 + cep(iv0,ixp)*(-l1sq + l2sq + l3sq)/(8.0d0*a)
    CALL triangle(x2,x0,x1,l1sq,l2sq,l3sq,a)
    sum2 = sum2 + (-l1sq + l2sq + l3sq)/(8.0d0*a) &
                + cep(iv0,ixp)*(l1sq - l2sq + l3sq)/(8.0d0*a)

    ! Now construct a second integral summing over all microelements
    sum1 = 0.0d0
    nef1 = neoff(if1,igrid)
    DO ixv = 1, nef1
      
      ! Find the vertex
      iv1 = voff(if1,ixv,igrid)
      long = vlong(iv1,igrid)
      lat = vlat(iv1,igrid)
      CALL ll2xyz(long,lat,x2)
      ! Find the preceding edge
      ie1 = eoff(ixv,if1,igrid)
      long = elong(ie1,igrid)
      lat = elat(ie1,igrid)
      CALL ll2xyz(long,lat,x3)
      CALL triangle(x1,x2,x3,l1sq,l2sq,l3sq,a)
      sum1 = sum1 + l1sq/(4.0d0*a)
      ! Find the following edge
      ixvp = ixv + 1
      IF (ixvp > nef1) ixvp = 1
      ie1 = eoff(ixvp,if1,igrid)
      long = elong(ie1,igrid)
      lat = elat(ie1,igrid)
      CALL ll2xyz(long,lat,x3)
      CALL triangle(x1,x2,x3,l1sq,l2sq,l3sq,a)
      sum1 = sum1 + l1sq/(4.0d0*a)

    ENDDO

    ! Finish calculation of degree of freedom
    cep(iv0,ix + nev0) = sum2/sum1

  ENDDO

ENDDO


! ---------------------------------------------------------------

END SUBROUTINE buildep

! ===============================================================

SUBROUTINE buildvd(igrid)

! Find the internal degrees of freedom that define the
! compound P0 elements on dual cells (space Vd).
!
! The compound element is constant over each dual cell with value
! equal to the inverse of the dual cell area. The micro elements
! are themselves constant over each supermesh cell with value
! equal to the inverse of the supermesh cell area. The coefficient
! of the micro element is thus the supermesh cell area divided by
! the compound dual cell area.
!
! There are 2*neofv micro elements in one compound element,
! and they are assumed to be ordered anticlockwise, starting with
! those incident on the first edge of the compound element.

USE grid

IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid
INTEGER :: nv, iv0, ie1, if1, ix, ixm
REAL*8 :: long, lat, x0(3), x1(3), x2(3), l1sq, l2sq, l3sq, a

! ---------------------------------------------------------------

! Number of dual cells (vertices) on this grid
nv = nvert(igrid)

! Number of coefficients to define element - depends on
! shape of dual cell.
ncvd(1:nv) = 2*neofv(:nv,igrid)


! Loop over dual cells
DO iv0 = 1, nv

  ! Dual cell centre
  long = vlong(iv0,igrid)
  lat = vlat(iv0,igrid)
  CALL ll2xyz(long,lat,x0)

  ! Loop over edges of dual cell
  DO ix = 1, neofv(iv0,igrid)
    ixm = ix - 1
    IF (ixm < 1) ixm = neofv(iv0,igrid)
    ie1 = eofv(ix,iv0,igrid)

    ! Edge crossing point
    long = elong(ie1,igrid)
    lat = elat(ie1,igrid)
    CALL ll2xyz(long,lat,x1)

    ! First micro element 
    if1 = fofv(iv0,ixm,igrid)
    long = flong(if1,igrid)
    lat = flat(if1,igrid)
    CALL ll2xyz(long,lat,x2)
    CALL triangle(x0,x1,x2,l1sq,l2sq,l3sq,a)
    cvd(iv0,2*ix - 1) = a/varea(iv0,igrid)

    ! Second micro element
    if1 = fofv(iv0,ix,igrid)
    long = flong(if1,igrid)
    lat = flat(if1,igrid)
    CALL ll2xyz(long,lat,x2)
    CALL triangle(x0,x1,x2,l1sq,l2sq,l3sq,a)
    cvd(iv0,2*ix) = a/varea(iv0,igrid)

  ENDDO

ENDDO


! ---------------------------------------------------------------

END SUBROUTINE buildvd

! ===============================================================

SUBROUTINE buildsd(igrid)

! Find the internal degrees of freedom that define the
! compound N0 elements on dual cells (space Sd).
!
! Each compound element is associated with a unit tangential circulation
! at some dual edge. The first two degrees of freedom are
! the circulations at the two micro edges that make up the dual edge.
! The sign convention is the same as the sign convention for
! the dual edge.
! The remaining degrees of freedom are the micro edge circulations
! internal to the dual cells either side of the dual edge.
! These dofs are assumed to be ordered anticlockwise starting
! with the microedge that touches the middle of the dual edge,
! for the first dual cell incident on the dual edge, then for
! the second dual cell incident on the dual edge. The sign convention
! is that outwards from the dual cell centre is positive.
!
! The compound element has constant vorticity in the dual cell
! each side of the dual edge. This constraint fixes all but two
! degrees of freedom. These last two are fixed by demanding that
! a certain measure of the divergence in each compound cell should
! vanish (this is the usual mixed finite element divergence defined
! by integration by parts against a test function on the supermesh).

USE grid

IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid
INTEGER :: ne, ie0, iv1, iv2, offset(2), ixf, if1, ixv, nev1, &
           ixe, ixc, ixe1, ixf1, ixf2, ixcp, ix, i1, i2, ie1
REAL*8 :: long, lat, x0(3), x1(3), x2(3), x3(3), l1, sgncrc, atot, &
          l1sq, l2sq, l3sq, a, sum1, sum2, contrib, c00, sg

! ---------------------------------------------------------------

! Number of edges on this grid
ne = nedge(igrid)

! Number of coefficients to define element - depends on
! shape of cells either side of edge. Computed below.


! Loop over edges
DO ie0 = 1, ne

  ! Find number of coefficients for this edge
  iv1 = vofe(1,ie0,igrid)
  iv2 = vofe(2,ie0,igrid)
  ncsd(ie0) = 2*(neofv(iv1,igrid) + neofv(iv2,igrid) + 1)

  ! Save related offsets
  offset(1) = 2
  offset(2) = 2 + 2*neofv(iv1,igrid)

  ! Find coordinates of edge ie0
  long = elong(ie0,igrid)
  lat = elat(ie0,igrid)
  CALL ll2xyz(long,lat,x0)

  ! Find the lengths of the two micro edges that make up
  ! edge ie0, and hence determine the first two coefficients
  DO ixf = 1, 2
    if1 = fnxte(ixf,ie0,igrid)
    long = flong(if1,igrid)
    lat = flat(if1,igrid)
    CALL ll2xyz(long,lat,x1)
    l1 = SQRT((x1(1) - x0(1))**2 + (x1(2) - x0(2))**2 + (x1(3) - x0(3))**2)
    csd(ie0,ixf) = l1/ddist(ie0,igrid)
  ENDDO
 
  ! Loop over dual cells next to edge ie0
  DO ixv = 1, 2
    iv1 = vofe(ixv,ie0,igrid)

    ! How many edges does this dual cell have
    nev1 = neofv(iv1,igrid)

    ! Is circulation positive or negative for dual cell iv1 ?
    IF (ixv == 1) THEN
      sgncrc = -1.0d0
    ELSE
      sgncrc = 1.0d0
    ENDIF

    ! Area of dual cell iv1
    atot = varea(iv1,igrid)

    ! Find coordinates of centre of dual cell iv1
    long = vlong(iv1,igrid)
    lat = vlat(iv1,igrid)
    CALL ll2xyz(long,lat,x1)

    ! Which edge of dual cell iv1 is ie0 ?
    ixe = 1
    DO WHILE (ie0 .NE. eofv(ixe,iv1,igrid))
      ixe = ixe + 1
    ENDDO

    ! Reset pointer to dofs
    ixc = offset(ixv)

    ! Internal circulations
    ! First guess to get the vorticity right
    ! First two internal circulations
    ixc = ixc + 1
    csd(ie0,ixc) = 0.0d0
    ixf = ixv
    if1 = fnxte(ixf,ie0,igrid)
    long = flong(if1,igrid)
    lat = flat(if1,igrid)
    CALL ll2xyz(long,lat,x2)
    CALL triangle(x1,x0,x2,l1sq,l2sq,l3sq,a)
    ixc = ixc + 1
    csd(ie0,ixc) = sgncrc*(csd(ie0,ixf) - a/atot)

    ! Loop to find remaining internal circulations
    ixe1 = ixe
    DO ix = 1, nev1 - 1
      ixf1 = ixe1
      ixe1 = ixe1 + 1
      IF (ixe1 > nev1) ixe1 = 1
      ixf2 = ixe1
      ie1 = eofv(ixe1,iv1,igrid)
      long = elong(ie1,igrid)
      lat = elat(ie1,igrid)
      CALL ll2xyz(long,lat,x3)
      if1 = fofv(iv1,ixf1,igrid)
      long = flong(if1,igrid)
      lat = flat(if1,igrid)
      CALL ll2xyz(long,lat,x2)
      CALL triangle(x1,x2,x3,l1sq,l2sq,l3sq,a)
      ixc = ixc + 1
      csd(ie0,ixc) = csd(ie0,ixc - 1) - sgncrc*a/atot
      if1 = fofv(iv1,ixf2,igrid)
      long = flong(if1,igrid)
      lat = flat(if1,igrid)
      CALL ll2xyz(long,lat,x2)
      CALL triangle(x1,x2,x3,l1sq,l2sq,l3sq,a)
      ixc = ixc + 1
      csd(ie0,ixc) = csd(ie0,ixc - 1) - sgncrc*a/atot
    ENDDO

    ! Now correct to make vorticity vanish
    ! Accumulate integrals over micro elements
    sum1 = 0.0d0
    sum2 = 0.0d0
    ! First the contributions from the circulations at edge ie0
    sg = -1.0d0
    DO ixf = 1, 2
      if1 = fnxte(ixf,ie0,igrid)
      long = flong(if1,igrid)
      lat = flat(if1,igrid)
      CALL ll2xyz(long,lat,x2)
      CALL triangle(x1,x0,x2,l1sq,l2sq,l3sq,a)
      sum2 = sum2 + sg*csd(ie0,ixf)*(l2sq - l3sq)/(12.0d0*a)
      sg = -sg
    ENDDO

    ! Now the contributions from all the internal circulations
    ! Loop over dual vertices with two micro elements on each vertex
    ixe1 = ixe
    ixf1 = ixe
    ixc = offset(ixv)
    DO ix = 1, nev1
      
      ! Find the dual vertex
      if1 = fofv(iv1,ixf1,igrid)
      long = flong(if1,igrid)
      lat = flat(if1,igrid)
      CALL ll2xyz(long,lat,x2)
      ! Find the preceding edge
      ie1 = eofv(ixe1,iv1,igrid)
      long = elong(ie1,igrid)
      lat = elat(ie1,igrid)
      CALL ll2xyz(long,lat,x3)
      CALL triangle(x1,x2,x3,l1sq,l2sq,l3sq,a)
      ixc = ixc + 1
      ixcp = ixc + 1
      sum1 = sum1 + l1sq/(4.0d0*a)
      contrib = ( csd(ie0,ixc )*(3.0d0*l1sq + l3sq - l2sq)   &
                + csd(ie0,ixcp)*(3.0d0*l1sq + l2sq - l3sq) ) &
                                               /(24.0d0*a)
      sum2 = sum2 + contrib
      ! Find the following edge
      ixe1 = ixe1 + 1
      IF (ixe1 > nev1) ixe1 = 1
      ie1 = eofv(ixe1,iv1,igrid)
      long = elong(ie1,igrid)
      lat = elat(ie1,igrid)
      CALL ll2xyz(long,lat,x3)
      CALL triangle(x1,x2,x3,l1sq,l2sq,l3sq,a)
      ixc = ixc + 1
      ixcp = ixc + 1
      IF (ix == nev1) ixcp = offset(ixv) + 1
      sum1 = sum1 + l1sq/(4.0d0*a)
      contrib = ( csd(ie0,ixcp)*(3.0d0*l1sq + l3sq - l2sq)   &
                + csd(ie0,ixc )*(3.0d0*l1sq + l2sq - l3sq) ) &
                                               /(24.0d0*a)
      sum2 = sum2 + contrib

      ixf1 = ixf1 + 1
      IF (ixf1 > nev1) ixf1 = 1
    ENDDO

    ! Correction to remove vorticity
    c00 = - sum2/sum1
    i1 = offset(ixv) + 1
    i2 = offset(ixv) + 2*nev1
    csd(ie0,i1:i2) = csd(ie0,i1:i2) + c00

  ENDDO

ENDDO


! ---------------------------------------------------------------

END SUBROUTINE buildsd

! ===============================================================

SUBROUTINE buildL(igrid)

! To build the stencil and coefficients for the Vp mass matrix L
!
! The stencil is a single cell. The coefficient is just the reciprocal
! of the cell area, but let's go through the motions of summing the
! integrals over microelements to check.

USE grid

IMPLICIT NONE
INTEGER, INTENT(IN) :: igrid
INTEGER :: nf, if0, ix, ixm, ie1, iv1
REAL*8 :: long, lat, x0(3), x1(3), x2(3), c, &
          l1sq, l2sq, l3sq, a, sum1

! ---------------------------------------------------------------

nf = nface(igrid)

! Loop over primal cells
DO if0 = 1, nf

  ! Stencil is just the cell itself
  nlsten(if0,igrid) = 1
  lsten(1,if0,igrid) = if0

  sum1 = 0.0d0

  ! Cell centre
  long = flong(if0,igrid)
  lat = flat(if0,igrid)
  CALL ll2xyz(long,lat,x0)

  ! Loop over edges of cell
  DO ix = 1, neoff(if0,igrid)
    ixm = ix - 1
    IF (ixm < 1) ixm = neoff(if0,igrid)
    ie1 = eoff(ix,if0,igrid)

    ! Edge crossing point
    long = elong(ie1,igrid)
    lat = elat(ie1,igrid)
    CALL ll2xyz(long,lat,x1)

    ! First micro element 
    iv1 = voff(if0,ixm,igrid)
    long = vlong(iv1,igrid)
    lat = vlat(iv1,igrid)
    CALL ll2xyz(long,lat,x2)
    CALL triangle(x0,x1,x2,l1sq,l2sq,l3sq,a)
    c = cvp(if0,2*ix - 1)
    sum1 = sum1 + c*c/a

    ! Second micro element
    iv1 = voff(if0,ix,igrid)
    long = vlong(iv1,igrid)
    lat = vlat(iv1,igrid)
    CALL ll2xyz(long,lat,x2)
    CALL triangle(x0,x1,x2,l1sq,l2sq,l3sq,a)
    c = cvp(if0,2*ix)
    sum1 = sum1 + c*c/a
  ENDDO
  lmass(1,if0,igrid) = sum1

ENDDO


! ---------------------------------------------------------------

END SUBROUTINE buildL

! ===============================================================

SUBROUTINE buildM(igrid)

! To build the stencil and coefficients for the Sp mass matrix M
!
! The stencil for edge ie0 comprises all the edges of the two
! primal cells either side of edge ie0
!
! Also build the stencil and coefficients for the T operator
!
! The stencil for cell if0 is all the edges of cell if0

USE grid

IMPLICIT NONE
INTEGER, INTENT(IN) :: igrid
INTEGER :: ne, ie0, ixsten, ixf, if1, ix0, ixe, nef1, off1, &
           ixf2, ie2, ix3, ie3, ixm, ixc, ixp, ixm2, ixc2, ixp2, &
           disp, iv2, ixv
REAL*8 :: long, lat, x0(3), x1(3), x2(3), x3(3), c, &
          l1sq, l2sq, l3sq, a, sum1, &
          a1(2), a2(2), b1(2*nefmx), b2(2*nefmx)

! ---------------------------------------------------------------

ntsten(:,igrid) = neoff(:,igrid)
tsten(:,:,igrid) = eoff(:,:,igrid)

ne = nedge(igrid)

! Loop over edges
DO ie0 = 1, ne

  ! Locate the edge crossing point
  long = elong(ie0,igrid)
  lat = elat(ie0,igrid)
  CALL ll2xyz(long,lat,x0)

  ixsten = 1

  ! Make sure first stencil edge is ie0 itself
  msten(1,ie0,igrid) = ie0
  mmass(1,ie0,igrid) = 0.0d0

  ! Loop over the two primal cells either side
  DO ixf = 1, 2
    if1 = fnxte(ixf,ie0,igrid)
    nef1 = neoff(if1,igrid)

    ! Extract coefficients of ie0 basis function in this cell
    IF (ixf == 1) THEN
      a1(1) = csp(ie0,1)
      a1(2) = csp(ie0,2)
      b1(1:2*nef1) = csp(ie0,3:2*(nef1 + 1))
    ELSE
      off1 = 2*neoff(fnxte(1,ie0,igrid),igrid)
      a1(1) = -csp(ie0,2)
      a1(2) = -csp(ie0,1)
      b1(1:2*nef1) = csp(ie0,3 + off1:2*(nef1 + 1) + off1)
    ENDIF

    ! Locate the primal cell centre
    long = flong(if1,igrid)
    lat = flat(if1,igrid)
    CALL ll2xyz(long,lat,x1)

    ! Which edge of if1 is ie0 ?
    ix0 = 1
    DO WHILE (eoff(ix0,if1,igrid) .NE. ie0)
      ix0 = ix0 + 1
    ENDDO

    ! Loop over edges of if1
    DO 	ixe = 1, nef1
      ie2 = eoff(ixe,if1,igrid)
      
      ! Which face of ie2 is if1 ?
      ixf2 = 1
      IF (fnxte(2,ie2,igrid) == if1) ixf2 = 2

      ! Locate the edge crossing point
      long = elong(ie2,igrid)
      lat = elat(ie2,igrid)
      CALL ll2xyz(long,lat,x3)

      ! Extract coefficients of ie2 basis function in this cell
      IF (ixf2 == 1) THEN
        a2(1) = csp(ie2,1)
        a2(2) = csp(ie2,2)
        b2(1:2*nef1) = csp(ie2,3:2*(nef1 + 1))
      ELSE
        off1 = 2*neoff(fnxte(1,ie2,igrid),igrid)
        a2(1) = -csp(ie2,2)
        a2(2) = -csp(ie2,1)
        b2(1:2*nef1) = csp(ie2,3 + off1:2*(nef1 + 1) + off1)
      ENDIF

      ! Displacement of edge ie2 from edge ie0 around cell if1
      disp = ixe - ix0
      IF (disp < 0) disp = disp + nef1

      ! Integrate the product of the two basis functions

      IF (disp > 0) THEN
        ! Initialize sum to zero
        sum1 = 0.0d0
      ELSE
        ! We need the product of the exterior flux contributions
        ixv = ix0 - 1
        IF (ixv < 1) ixv = nef1
        iv2 = voff(if1,ixv,igrid)
        long = vlong(iv2,igrid)
        lat = vlat(iv2,igrid)
        CALL ll2xyz(long,lat,x2)
        CALL triangle(x1,x2,x0,l1sq,l2sq,l3sq,a)
        sum1 = a1(1)*a2(1)*(3.0d0*(l2sq + l3sq) - l1sq)/(48.0d0*a)
        ixv = ix0
        iv2 = voff(if1,ixv,igrid)
        long = vlong(iv2,igrid)
        lat = vlat(iv2,igrid)
        CALL ll2xyz(long,lat,x2)
        CALL triangle(x1,x0,x2,l1sq,l2sq,l3sq,a)
        sum1 = sum1 + a1(2)*a2(2)*(3.0d0*(l2sq + l3sq) - l1sq)/(48.0d0*a)
      ENDIF

      ! Now cross products of exterior and interior fluxes
      ! Indices of edges of micro elements on ie0 relative to ie2
      ixm = 2*(nef1 - disp)
      ixc = ixm + 1
      IF (ixc > 2*nef1) ixc = 1
      ixp = ixc + 1
      ixv = ix0 - 1
      IF (ixv < 1) ixv = nef1
      iv2 = voff(if1,ixv,igrid)
      long = vlong(iv2,igrid)
      lat = vlat(iv2,igrid)
      CALL ll2xyz(long,lat,x2)
      CALL triangle(x1,x2,x0,l1sq,l2sq,l3sq,a)
      sum1 = sum1 + a1(1)*( - b2(ixc)*(l1sq + l2sq - 3.0d0*l3sq)    &
                            + b2(ixm)*(l1sq + l3sq - 3.0d0*l2sq) )  &
                                                     /(48.0d0*a)
      ixv = ix0
      iv2 = voff(if1,ixv,igrid)
      long = vlong(iv2,igrid)
      lat = vlat(iv2,igrid)
      CALL ll2xyz(long,lat,x2)
      CALL triangle(x1,x0,x2,l1sq,l2sq,l3sq,a)
      sum1 = sum1 + a1(2)*( - b2(ixp)*(l1sq + l2sq - 3.0d0*l3sq)    &
                            + b2(ixc)*(l1sq + l3sq - 3.0d0*l2sq) )  &
                                                     /(48.0d0*a)
      ! Indices of edges of micro elements on ie2 relative to ie0
      ixm = 2*disp
      ixc = ixm + 1
      ixp = ixc + 1
      IF (ixm == 0) ixm = 2*nef1
      ixv = ixe - 1
      IF (ixv < 1) ixv = nef1
      iv2 = voff(if1,ixv,igrid)
      long = vlong(iv2,igrid)
      lat = vlat(iv2,igrid)
      CALL ll2xyz(long,lat,x2)
      CALL triangle(x1,x2,x3,l1sq,l2sq,l3sq,a)
      sum1 = sum1 + a2(1)*( - b1(ixc)*(l1sq + l2sq - 3.0d0*l3sq)    &
                            + b1(ixm)*(l1sq + l3sq - 3.0d0*l2sq) )  &
                                                     /(48.0d0*a)
      ixv = ixe
      iv2 = voff(if1,ixv,igrid)
      long = vlong(iv2,igrid)
      lat = vlat(iv2,igrid)
      CALL ll2xyz(long,lat,x2)
      CALL triangle(x1,x3,x2,l1sq,l2sq,l3sq,a)
      sum1 = sum1 + a2(2)*( - b1(ixp)*(l1sq + l2sq - 3.0d0*l3sq)    &
                            + b1(ixc)*(l1sq + l3sq - 3.0d0*l2sq) )  &
                                                     /(48.0d0*a)

      ! And finally products of interior flux contributions
      ! Loop over edges of if1
      DO ix3 = 1, nef1

        ! Indices of elements on edge ie3 relative to ie0
        disp = ix3 - ix0
        IF (disp < 0) disp = disp + nef1
        ixm = 2*disp
        ixc = ixm + 1
        ixp = ixc + 1
        IF (ixm == 0) ixm = 2*nef1
        ! Indices of elements on edge ie3 relative to ie2
        disp = ix3 - ixe
        IF (disp < 0) disp = disp + nef1
        ixm2 = 2*disp
        ixc2 = ixm2 + 1
        ixp2 = ixc2 + 1
        IF (ixm2 == 0) ixm2 = 2*nef1

        ! Locate edge ie3
        ie3 = eoff(ix3,if1,igrid)
        long = elong(ie3,igrid)
        lat = elat(ie3,igrid)
        CALL ll2xyz(long,lat,x3)

        ! Two micro elements incident on edge ie3
        ixv = ix3 - 1
        IF (ixv < 1) ixv = nef1
        iv2 = voff(if1,ixv,igrid)
        long = vlong(iv2,igrid)
        lat = vlat(iv2,igrid)
        CALL ll2xyz(long,lat,x2)
        CALL triangle(x1,x2,x3,l1sq,l2sq,l3sq,a)
        sum1 = sum1 + ( b1(ixm)*b2(ixm2)*(3.0d0*(l1sq + l2sq) - l3sq)        &
                      + b1(ixc)*b2(ixc2)*(3.0d0*(l1sq + l3sq) - l2sq)        &
                     - (b1(ixm)*b2(ixc2) + b1(ixc)*b2(ixm2))                 &
                                        *(l2sq + l3sq - 3.0d0*l1sq)    )     &
                                                       /(48.0d0*a)
        ixv = ix3
        iv2 = voff(if1,ixv,igrid)
        long = vlong(iv2,igrid)
        lat = vlat(iv2,igrid)
        CALL ll2xyz(long,lat,x2)
        CALL triangle(x1,x3,x2,l1sq,l2sq,l3sq,a)
        sum1 = sum1 + ( b1(ixc)*b2(ixc2)*(3.0d0*(l1sq + l2sq) - l3sq)        &
                      + b1(ixp)*b2(ixp2)*(3.0d0*(l1sq + l3sq) - l2sq)        &
                     - (b1(ixc)*b2(ixp2) + b1(ixp)*b2(ixc2))                 &
                                        *(l2sq + l3sq - 3.0d0*l1sq)    )     &
                                                       /(48.0d0*a)

      ENDDO

      ! Finally save the result
      IF (ie2 == ie0) THEN
        mmass(1,ie0,igrid) = mmass(1,ie0,igrid) + sum1
      ELSE
        ixsten = ixsten + 1
        msten(ixsten,ie0,igrid) = ie2
        mmass(ixsten,ie0,igrid) = sum1
      ENDIF
      tcoeff(ixe,ix0,if1,igrid) = sum1

    ENDDO

  ENDDO
  nmsten(ie0,igrid) = ixsten

ENDDO


! ---------------------------------------------------------------

END SUBROUTINE buildM

! ===============================================================

SUBROUTINE buildR(igrid)

! To build the R operator to transfer from Vp to Ep
!
! The operator is built in two forms, one suitable for
! looping over faces (which can be used for the TRiSK
! implementation of W), one suitable for looping over
! vertices (avoids reduce_all).

USE grid

IMPLICIT NONE
INTEGER, INTENT(IN) :: igrid
INTEGER :: nv, nf, iv0, nev1, if1, ixv0, ixf, ixfp, nef1, &
           ixv, ixvm, iv1, ixe, ie1, ne1
REAL*8 :: v1, v2, vc, long, lat, x1(3), x2(3), x3(3), &
          sum1, l1sq, l2sq, l3sq, a, aby3

! ---------------------------------------------------------------

nv = nvert(igrid)

nrsten(:,igrid) = neoff(:,igrid)

! Loop over vertices
DO iv0 = 1, nv

  nev1 = neofv(iv0,igrid)

  ! Loop over the faces around iv0
  DO ixf = 1, nev1

    if1 = fofv(iv0,ixf,igrid)

    ! Locate the face centre
    long = flong(if1,igrid)
    lat = flat(if1,igrid)
    CALL ll2xyz(long,lat,x1)

    ! Pick out degrees of freedom defining Ep basis element
    ! in this cell
    ixfp = ixf + 1
    if (ixfp > nev1) ixfp = 1
    v1 = cep(iv0,ixfp)
    v2 = cep(iv0,ixf)
    vc = cep(iv0,ixf + nev1)

    ! Which vertex of if1 is iv0 ?
    ixv0 = 1
    DO WHILE (voff(if1,ixv0,igrid) .NE. iv0)
      ixv0 = ixv0 + 1
    ENDDO

    ! Now loop over the micro elements in if1, integrating
    ! the product of the basis functions
    nef1 = neoff(if1,igrid)
    sum1 = 0.0d0
    DO ixe = 1, nef1
      ie1 = eoff(ixe,if1,igrid)
      ! Locate edge crossing point
      long = elong(ie1,igrid)
      lat = elat(ie1,igrid)
      CALL ll2xyz(long,lat,x2)
      ! Two micro elements on this edge
      ixv = ixe
      ixvm = ixv - 1
      IF (ixvm < 1) ixvm = nef1
      ! First one
      iv1 = voff(if1,ixvm,igrid)
      long = vlong(iv1,igrid)
      lat = vlat(iv1,igrid)
      CALL ll2xyz(long,lat,x3)
      CALL triangle(x1,x3,x2,l1sq,l2sq,l3sq,a)
      aby3 = a/3.0d0
      sum1 = sum1 + vc*aby3
      IF (ixv  == ixv0) sum1 = sum1 + v1*aby3
      IF (ixvm == ixv0) sum1 = sum1 + (1.0d0 + v2)*aby3
      ! Second one
      iv1 = voff(if1,ixv,igrid)
      long = vlong(iv1,igrid)
      lat = vlat(iv1,igrid)
      CALL ll2xyz(long,lat,x3)
      CALL triangle(x1,x2,x3,l1sq,l2sq,l3sq,a)
      aby3 = a/3.0d0
      sum1 = sum1 + vc*aby3
      IF (ixv  == ixv0) sum1 = sum1 + (1.0d0 + v1)*aby3
      IF (ixvm == ixv0) sum1 = sum1 + v2*aby3
    ENDDO
    rsten(if1,ixv0,igrid) = iv0
    rcoeff(if1,ixv0,igrid) = sum1/farea(if1,igrid)

  ENDDO

ENDDO




! Re-tabulate the coefficients of the R operator in a form
! suitable for looping over vertices rather than faces

nf = nface(igrid)

nrxsten(:,igrid) = neofv(:,igrid)

DO if1 = 1, nf

  ne1 = neoff(if1,igrid)
  DO ixv = 1, ne1

    iv1 = voff(if1,ixv,igrid)

    ! What is the index of face if1 relative to vertex iv1 ?
    ixf = 1
    DO WHILE (fofv(iv1,ixf,igrid) .NE. if1)
      ixf = ixf + 1
    ENDDO
    
    ! Save the stencil face and coefficient
    rxsten(ixf,iv1,igrid) = if1
    rxcoeff(ixf,iv1,igrid) = rcoeff(if1,ixv,igrid)

  ENDDO

ENDDO



! ---------------------------------------------------------------

END SUBROUTINE buildR

! ===============================================================

SUBROUTINE buildW(igrid)

! To construct the stencil and coefficients for the W operator,
! which maps from Ep to Ep.

USE grid

IMPLICIT NONE
INTEGER, INTENT(IN) :: igrid
INTEGER :: ie0, ixf, if1, ne1, ix1, ix, ixv, ix2, ie2, &
           isten, ne
REAL*8 :: ss, w


! ---------------------------------------------------------------

! The W operator is built from the R operator a la TRiSK
! (avoids having to do more integrals, which would be another
! way of doing it)

ne = nedge(igrid)

! Loop over edges
DO ie0 = 1, ne

  isten = 0

  ! Loop over the faces either side of this edge
  DO ixf = 1, 2
    if1 = fnxte(ixf,ie0,igrid)
    ne1 = neoff(if1,igrid)
    ss = -0.5d0

    ! Which edge of face if1 is edge ie0?
    ix1 = 1
    DO WHILE (eoff(ix1,if1,igrid) .NE. ie0)
      ix1 = ix1 + 1
    ENDDO

    ! Find the contribution from every other
    ! edge of this face
    DO ix = 0, ne1 - 2
      ixv = MOD(ix1 + ix - 1,ne1) + 1
      ix2 = MOD(ix1 + ix,ne1) + 1
      ie2 = eoff(ix2,if1,igrid)
      ss = ss + rcoeff(if1,ixv,igrid)
      w = -ss*eoffin(ix1,if1,igrid)*eoffin(ix2,if1,igrid)
      isten = isten + 1
      wsten(isten,ie0,igrid) = ie2
      wcoeff(isten,ie0,igrid) = w
    ENDDO
  ENDDO
  nwsten(ie0,igrid) = isten
ENDDO

! ---------------------------------------------------------------

END SUBROUTINE buildW

! ===============================================================

SUBROUTINE buildJ(igrid)

! To construct the stencil and coefficients for the J operator,
! which maps from Vd to Ep.

USE grid

IMPLICIT NONE
INTEGER, INTENT(IN) :: igrid
INTEGER :: nv, iv0, nev1, ixf, nef1, ns, ixfp, ix, iv1, &
           ixv0, ixvm, ixvp, ixe, ixep, ie1, ixv, if1
REAL*8 :: v1, v2, vc, long, lat, x1(3), x2(3), x3(3), &
          l1sq, l2sq, l3sq, a, aby3, sum1

! ---------------------------------------------------------------

nv = nvert(igrid)

! Loop over vertices
DO iv0 = 1, nv

  nev1 = neofv(iv0,igrid)  

  ! Initialize to make sure iv0 itself is first in the stencil
  jsten(:,iv0,igrid) = 0
  jsten(1,iv0,igrid) = iv0
  ns = 1
  jstar(:,iv0,igrid) = 0.0d0

  ! Loop over the faces comprising the element in Ep
  DO ixf = 1, nev1

    if1 = fofv(iv0,ixf,igrid)
    nef1 = neoff(if1,igrid)

    ! Locate the face centre
    long = flong(if1,igrid)
    lat = flat(if1,igrid)
    CALL ll2xyz(long,lat,x1)

    ! Which vertex of if1 is iv0 ?
    ixv0 = 1
    DO WHILE (voff(if1,ixv0,igrid) .NE. iv0)
      ixv0 = ixv0 + 1
    ENDDO
    ixvm = ixv0 - 1
    IF (ixvm < 1) ixvm = nef1
    ixvp = ixv0 + 1
    IF (ixvp > nef1) ixvp = 1

    ! Pick out degrees of freedom defining Ep basis element
    ! in this cell
    ixfp = ixf + 1
    if (ixfp > nev1) ixfp = 1
    v1 = cep(iv0,ixfp)
    v2 = cep(iv0,ixf)
    vc = cep(iv0,ixf + nev1)

    ! Loop over vertices of face if1
    DO ixv = 1, nef1
      ixe = ixv
      ixep = ixe + 1
      IF (ixep > nef1) ixep = 1
      iv1 = voff(if1,ixv,igrid)
      ! Locate the vertex iv1
      long = vlong(iv1,igrid)
      lat = vlat(iv1,igrid)
      CALL ll2xyz(long,lat,x2)
      ! Two micro elements are incident on vertex iv1 in cell if1
      ! First one
      ie1 = eoff(ixe,if1,igrid)
      long = elong(ie1,igrid)
      lat = elat(ie1,igrid)
      CALL ll2xyz(long,lat,x3)
      CALL triangle(x1,x3,x2,l1sq,l2sq,l3sq,a)
      aby3 = a/3.0d0
      sum1 = vc*aby3
      IF (ixv == ixv0) sum1 = sum1 + (1.0d0 + v1)*aby3
      IF (ixv == ixvp) sum1 = sum1 + v2*aby3
      ! Second one
      ie1 = eoff(ixep,if1,igrid)
      long = elong(ie1,igrid)
      lat = elat(ie1,igrid)
      CALL ll2xyz(long,lat,x3)
      CALL triangle(x1,x2,x3,l1sq,l2sq,l3sq,a)
      aby3 = a/3.0d0
      sum1 = sum1 + vc*aby3
      IF (ixv == ixvm) sum1 = sum1 + v1*aby3
      IF (ixv == ixv0) sum1 = sum1 + (1.0d0 + v2)*aby3
      ! Now include this contribution in the coefficient
      CALL findinlist(iv1,jsten(:,iv0,igrid),njsmx,ix)
      ns = MAX(ns,ix)
      jstar(ix,iv0,igrid) = jstar(ix,iv0,igrid) + sum1/varea(iv1,igrid)
    ENDDO

  ENDDO
  njsten(iv0,igrid) = ns

ENDDO


! ---------------------------------------------------------------

END SUBROUTINE buildJ

! ===============================================================

SUBROUTINE buildH(igrid)

! To construct the stencil and coefficients for the H operator,
! which maps from Sd to Sp.

Use grid

IMPLICIT NONE
INTEGER, INTENT(IN) :: igrid
INTEGER :: ne, ie0, ns, ixf, if1, nef1, iv1, nev1, &
           ie2, ixf2, ixe2, offset, offset2, ivfrst, d, &
           ip3, ip4, ip5, iq3, iq4, iq5, ixe0, ixv, ix
REAL*8 :: flx1, flx2, flx3, flx4, flx5, crc1, crc2, crc3, crc4, crc5, &
          extf1, extf2, extc1, extc2, sum1

! ---------------------------------------------------------------

ne = nedge(igrid)

! Loop over all edges
DO ie0 = 1, ne

  ! Initialize to make sure ie0 itself is first in the stencil
  hsten(:,ie0,igrid) = 0
  hsten(1,ie0,igrid) = ie0
  ns = 1
  hstar(:,ie0,igrid) = 0.0d0

  ! Two faces either side of edge ie0
  DO ixf = 1, 2

    if1 = fnxte(ixf,ie0,igrid)
    nef1 = neoff(if1,igrid)

    ! Which edge of if1 is ie0 ?
    ixe0 = 1
    DO WHILE (eoff(ixe0,if1,igrid) .NE. ie0)
      ixe0 = ixe0 + 1
    ENDDO

    ! Offset to coefficients for internal degrees of freedom,
    ! and orient external fluxes
    IF (ixf == 1) THEN
      offset = 2
      extf1 = csp(ie0,1)
      extf2 = csp(ie0,2)
    ELSE
      offset = 2 + 2*neoff(fnxte(1,ie0,igrid),igrid)
      extf1 = -csp(ie0,2)
      extf2 = -csp(ie0,1)
    ENDIF

    ! Loop over vertices of face if1
    DO ixv = 1, nef1

      iv1 = voff(if1,ixv,igrid)

      ! Pick out coefficients defining the part of the Sp
      ! element in this corner
      d = ixv - ixe0
      IF (d < 0) d = d + nef1
      IF (d == 0) THEN
        flx1 = extf2
      ELSE
        flx1 = 0.0d0
      ENDIF
      IF (d == nef1 - 1) THEN
        flx2 = extf1
      ELSE
        flx2 = 0.0d0
      ENDIF
      ip3 = 2*d + 1
      flx3 = csp(ie0,offset + ip3)
      ip4 = ip3 + 1
      flx4 = csp(ie0,offset + ip4)
      ip5 = ip4 + 1
      IF (ip5 > 2*nef1) ip5 = 1
      flx5 = csp(ie0,offset + ip5)

      ! Which face of iv1 is if1 ?
      ixf2 = 1
      DO WHILE (fofv(iv1,ixf2,igrid) .NE. if1)
        ixf2 = ixf2 + 1
      ENDDO

      ! Loop over all edges around vertex iv1 - these
      ! are the stencil edges
      nev1 = neofv(iv1,igrid)
      DO ixe2 = 1, nev1

        ie2 = eofv(ixe2,iv1,igrid)
        
        ! Is iv1 the first or second vertex of edge ie2 ?
        ivfrst = vofe(1,ie2,igrid)
        IF (iv1 == ivfrst) THEN
          offset2 = 2
          extc1 = csd(ie2,1)
          extc2 = csd(ie2,2)
        ELSE
          offset2 = 2 + 2*neofv(ivfrst,igrid)
          extc1 = -csd(ie2,2)
          extc2 = -csd(ie2,1)
        ENDIF

        ! Pick out coefficients defining the part of the Sd
        ! element in this corner
        d = ixf2 - ixe2
        IF (d < 0) d = d + nev1
        IF (d == 0) THEN
          crc2 = -extc1
        ELSE
          crc2 = 0.0d0
        ENDIF
        IF (d == nev1 - 1) THEN
          crc1 = -extc2
        ELSE
          crc1 = 0.0d0
        ENDIF
        iq3 = 2*d + 1
        crc3 = csd(ie2,offset2 + iq3)
        iq4 = iq3 + 1
        crc4 = csd(ie2,offset2 + iq4)
        iq5 = iq4 + 1
        IF (iq5 > 2*nev1) iq5 = 1
        crc5 = csd(ie2,offset2 + iq5)

        ! Integral over the two overlapping microelements
        sum1 = (-flx1*crc4 + flx1*crc1 + flx4*crc5 + flx4*crc1   &
               + flx3*crc5 + flx3*crc4 - flx2*crc2 - flx2*crc4   &
               - flx5*crc3 - flx5*crc4 - flx4*crc3 + flx4*crc2)  &
                                             / 6.0d0

        ! Update stencil and coefficients
        CALL findinlist(ie2,hsten(:,ie0,igrid),nhsmx,ix)
        ns = MAX(ns,ix)
        hstar(ix,ie0,igrid) = hstar(ix,ie0,igrid) + sum1

      ENDDO

    ENDDO

  ENDDO
  nhsten(ie0,igrid) = ns

ENDDO


! ---------------------------------------------------------------

END SUBROUTINE buildH

! ===============================================================

SUBROUTINE buildInjectionAndRestriction

! To build restriction operator for multigrid
! Note this code assumes the same grid cell numbering as
! the gengrid_hex.f and gengrid_cube.f90 grid generation programs.

USE grid

IMPLICIT NONE
INTEGER :: rif0, if0, if1, if2, ifSFC, igrid, igridp, nf, ix, ixp, &
           n2, n, n2p, np, i, j, &
           t1, s1, s2, s3, s4, p, ixx

! ---------------------------------------------------------------

! Allocate array for size of operator stencil
ALLOCATE(nres(nfacex,ngrids-1))
ALLOCATE(ninj(nfacex,2:ngrids))

nres = 0
ninj = 0

IF ((nedgex == 3*(nfacex - 2)) .AND. (nefmx == 6) .AND. (nevmx == 3)) THEN

  ! It looks like we're on a hexagonal icoshedral grid, so build the
  ! appropriate multigrid restriction operator

  ! Stencil is 6 or 7 cells on standard buckyball grid
  nresmx = 7
  ninjmx = 2

  ALLOCATE(ressten(nresmx,nfacex,ngrids-1))
  ALLOCATE(reswgt(nresmx,nfacex,ngrids-1))

  ! Allocate arrays for operator stencils and coefficients
  ALLOCATE(injsten(ninjmx,nfacex,2:ngrids))
  ALLOCATE(injwgt(ninjmx,nfacex,2:ngrids))

  ! Initialize to zero
  injsten = 0
  injwgt = 0.0d0

  ressten = 0
  reswgt = 0.0d0

  ! And define stencil and coefficients
  DO igrid = 1, ngrids-1
    igridp = igrid + 1

    DO ifSFC = 1, nface(igrid)
      ! Face if0 exists on grid igrid and grid igrid+1 and is
      ! the centre of the stencil

      ! SFC reordering?
      IF (SFCIndexAvailable == 0) THEN
        if0 = ifSFC
        if1 = if0
      ELSE
        ! use reordered face traversal to account for bit reproducability
        if0 = fNewFaceIdInverse(ifSFC,igrid)
        ! lookup old index on coarse grid and new index position on finer grid
        if1 = fNewFaceIdInverse(fNewFaceId(if0, igrid), igridp)
      END IF

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! RESTRICTION
      ressten(1,if0,igrid) = if1
      reswgt(1,if0,igrid) = 1.0d0
      ! The neighbours of if1 on grid igrid+1 are the other members
      ! of the stencil
      nf = neoff(if1,igridp)
      nres(if0,igrid) = 1 + nf

      IF (if1 .lt. 1 .OR. if1 .gt. nface(igridp)) THEN
        WRITE(*,*) if1
        WRITE(*,*) nface(igridp)
        WRITE(*,*) "Invalid face id (1) -> EXIT"
        CALL EXIT(1)
      END IF

      DO ix = 1, nf
        ixp = ix + 1
        ressten(ixp,if0,igrid) = fnxtf(if1,ix,igridp)
        reswgt(ixp,if0,igrid) = 0.5d0
      ENDDO


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! INJECTION
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (.TRUE.) THEN
          DO ix = 1, nres(if0,igrid)
            if2 = ressten(ix, if0, igrid)
            ninj(if2,igridp) = ninj(if2,igridp)+1

            ! setup values in injection stencil
            injsten(ninj(if2,igridp),if2,igridp) = if0
            ! weight
            injwgt(ninj(if2,igridp),if2,igridp) = reswgt(ix,if0,igrid)
          END DO
      END IF

      IF (.FALSE.) THEN
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! INJECTION
          ! Here, we assume that there is only one coarse cell
          ! overlapping with the current one
          !
          IF (injwgt(1,if1,igridp) .ne. 0.0d0) THEN
            WRITE(*,*) "Weight already set! There can be only a 1:1 correspondence for the direct matching cell"
            WRITE(*,*) injwgt(1,if1,igridp)
            CALL EXIT(1)
          END IF

          IF (ninj(if1,igridp) .ne. 0) THEN
            WRITE(*,*) "ninj already set"
          END IF

          injsten(1,if1,igridp) = if0
          injwgt(1,if1,igridp) = 1.0d0
          ninj(if1,igridp) = 1

          ! The neighbours of if1 on grid igrid+1 are the other members
          ! of the stencil
          nf = neoff(if1,igridp)

          DO ix = 1, nf
            ! adjacent face index
            if2 = fnxtf(if1,ix,igridp)

            ! increment to account for additional stencil value
            ninj(if2,igridp) = ninj(if2,igridp)+1

            ! setup values in injection stencil
            injsten(ninj(if2,igridp),if2,igridp) = if0
            injwgt(ninj(if2,igridp),if2,igridp) = 0.5d0
          ENDDO
      END IF
    ENDDO
  ENDDO

! ---------------------------------------------------------------

ELSEIF ((nedgex == 2*nfacex) .AND. (nefmx == 4) .AND. (nevmx == 4)) THEN

  ! It looks like we're on a cubed sphere grid, so build the
  ! appropriate multigrid restriction operator

  ! Stencil is always 4 cells
  nresmx = 4
  ninjmx = 2

  ! Allocate arrays for operator stencils and coefficients
  ALLOCATE(injsten(ninjmx,nfacex,2:ngrids))
  ALLOCATE(injwgt(ninjmx,nfacex,2:ngrids))

  ALLOCATE(ressten(nresmx,nfacex,ngrids-1))
  ALLOCATE(reswgt(nresmx,nfacex,ngrids-1))

  ! Initialize to zero
  injsten = 0
  injwgt = 0.0d0

  ressten = 0
  reswgt = 0.0d0

  ! And define stencil and coefficients
  DO igrid = 1, ngrids-1
    igridp = igrid + 1
    ! Number of cells per face and cells per edge on grid igrid
    n2 = nface(igrid)/6
    n = NINT(SQRT(1.0d0*n2))
    ! And on grid igridp
    n2p = nface(igridp)/6
    np = NINT(SQRT(1.0d0*n2p))
    ! Loop over cells on a panel
    DO j = 1, n
      DO i = 1, n
        t1 = (j-1)*n + i
        s1 = (2*j - 2)*np + (2*i - 1)
        s2 = s1 + 1
        s3 = s2 + np
        s4 = s3 - 1
        ! Loop over panels
        DO p = 1, 6
          if0 = t1 + (p - 1)*n2
          IF (SFCIndexAvailable == 0) THEN
            nres(if0,igrid) = 4
            ressten(1,if0,igrid) = s1 + (p - 1)*n2p
            ressten(2,if0,igrid) = s2 + (p - 1)*n2p
            ressten(3,if0,igrid) = s3 + (p - 1)*n2p
            ressten(4,if0,igrid) = s4 + (p - 1)*n2p
          ELSE
            nres(fNewFaceIdInverse(if0, igrid),igrid) = 4
            ressten(1,fNewFaceIdInverse(if0, igrid),igrid) = fNewFaceIdInverse(s1 + (p-1)*n2p, igridp)
            ressten(2,fNewFaceIdInverse(if0, igrid),igrid) = fNewFaceIdInverse(s2 + (p-1)*n2p, igridp)
            ressten(3,fNewFaceIdInverse(if0, igrid),igrid) = fNewFaceIdInverse(s3 + (p-1)*n2p, igridp)
            ressten(4,fNewFaceIdInverse(if0, igrid),igrid) = fNewFaceIdInverse(s4 + (p-1)*n2p, igridp)
          END IF
          
	  IF (SFCIndexAvailable == 0) THEN

	  ELSE
		if0 = fNewFaceIdInverse(if0, igrid)
	  ENDIF
	  reswgt(1:4,if0,igrid) = 1.0d0

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! INJECTION
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          DO ix = 1, nres(if0,igrid)
            if2 = ressten(ix, if0, igrid)
            ixx = ninj(if2,igridp)+1

            ninj(if2,igridp) = ixx

            ! setup values in injection stencil
            injsten(ixx,if2,igridp) = if0
            ! weight
            injwgt(ixx,if2,igridp) = reswgt(ix,if0,igrid)
          END DO

        ENDDO
      ENDDO
    ENDDO
  ENDDO

! ---------------------------------------------------------------

ELSE

  PRINT *,' '
  PRINT *,'**********'
  PRINT *,'I do not recognize this grid.'
  PRINT *,'Please programme me to build the multigrid restriction'
  PRINT *,'operator for this grid. Thank you.'
  PRINT *,'**********'
  PRINT *,' '
  STOP

ENDIF



!*********************
! VALIDATE restriction stencil by finding correspondenes in injection stencil
!*********************

IF (.TRUE.) THEN
  WRITE(*,*) "validating injection stencil"
  DO igrid = 1, ngrids-1
    igridp = igrid + 1

    DO if2 = 1, nface(igridp)
      DO ix = 1, ninj(if2,igridp)

        if1 = injsten(ix, if2, igridp)

        ! search for corresponding entry in injection stencil
        DO ixx = 1, nres(if1, igrid)
          IF (ressten(ixx, if1, igrid) .eq. if2) THEN
            ! Validate weight
            IF (reswgt(ixx, if1, igrid) .ne. injwgt(ix, if2, igridp)) THEN
              WRITE(*,*) "VALIDATION FOR INJECTION WEIGHT FAILED"
              WRITE(*,*) injwgt(ixx, if1, igridp)
              WRITE(*,*) reswgt(ix, if2, igridp)
              CALL EXIT(1)
            END IF
            EXIT
          END IF
        END DO
        IF (ixx .eq. nres(if1, igrid)+1) THEN
            WRITE(*,*) "INJECTION VALIDATION FAILED (entry not found)"
            CALL EXIT(1)
        END IF

      ENDDO
    ENDDO
  END DO
  WRITE(*,*) "DONE"
END IF


! Update maximum value
ninjmx = MAXVAL(ninj)
nresmx = MAXVAL(nres)


END SUBROUTINE buildInjectionAndRestriction

! ===============================================================
!
SUBROUTINE triangle(x1,x2,x3,l1sq,l2sq,l3sq,area)

! Compute the squares of the lengths of the sides and the area of
! a triangle defined by three points x1, x2 and x3 in 3D Euclidean
! space. Side 1 is opposite vertex 1, etc.

IMPLICIT NONE
REAL*8, INTENT(IN) :: x1(3), x2(3), x3(3)
REAL*8, INTENT(OUT) :: l1sq, l2sq, l3sq, area
REAL*8 :: a, b, c

! ---------------------------------------------------------------

! Squares of side by Pythagoras
l1sq = (x2(1) - x3(1))**2 + (x2(2) - x3(2))**2 + (x2(3) - x3(3))**2
l2sq = (x3(1) - x1(1))**2 + (x3(2) - x1(2))**2 + (x3(3) - x1(3))**2
l3sq = (x1(1) - x2(1))**2 + (x1(2) - x2(2))**2 + (x1(3) - x2(3))**2

! Area by Heron's formula
a = SQRT(l1sq)
b = SQRT(l2sq)
c = SQRT(l3sq)
area = 0.25d0*SQRT((a + b + c)*(b + c - a)*(c + a - b)*(a + b - c))

! ---------------------------------------------------------------

END SUBROUTINE triangle

! ===============================================================

SUBROUTINE findinlist(i,list,nlx,ix)

! Find the index ix to the entry i in the input list.
! If i is not already in the list then it is added into
! the first empty (i.e. zero) position.

IMPLICIT NONE
INTEGER, INTENT(IN) :: i, nlx
INTEGER, INTENT(INOUT) :: list(nlx)
INTEGER, INTENT(OUT) :: ix


ix = 1
DO WHILE ((list(ix) .NE. i) .AND. (list(ix) .NE. 0))
  IF (ix == nlx) THEN
    PRINT *,'search past end of list in SUBROUTINE findinlist.'
    PRINT *,'i = ',i
    PRINT *,'list = ',list
    STOP
  ENDIF
  ix = ix + 1
ENDDO
list(ix) = i


END SUBROUTINE findinlist

! ===============================================================
!

      SUBROUTINE LL2XYZ(LONG,LAT,X)

!
!     To convert longitude and latitude to cartesian coordinates
!     on the unit sphere
!
      IMPLICIT NONE
!
      REAL*8 LONG,LAT,X(3),CLN,SLN,CLT,SLT
!
!     ------------------------------------------------------------------
!
      SLN=SIN(LONG)
      CLN=COS(LONG)
      SLT=SIN(LAT)
      CLT=COS(LAT)
!
      X(1)=CLN*CLT
      X(2)=SLN*CLT
      X(3)=SLT
!
!     ------------------------------------------------------------------
!
      RETURN
      END
!
!     ==================================================================
!
      SUBROUTINE XYZ2LL(X,LONG,LAT)
!
!     To convert cartesian coordinates to longitude and latitude
!
      IMPLICIT NONE
!   
      REAL*8 X(3),LONG,LAT,PI,TLN,TLT,R
!
!     -------------------------------------------------------------------
!
      PI=4.0D0*ATAN(1.0D0)
!
      IF (X(1).EQ.0.0D0) THEN
        IF (x(2).GE.0.0D0) THEN
          LONG=0.5D0*PI
        ELSE
          LONG=1.5D0*PI
        ENDIF
      ELSE
        TLN=X(2)/X(1)
        LONG=ATAN(TLN)
        IF (X(1).LT.0.0D0) THEN
          LONG=LONG+PI
        ENDIF
        IF (LONG.LT.0.0D0) THEN
          LONG=LONG+2.0D0*PI
        ENDIF
      ENDIF
!
      R=SQRT(X(1)*X(1)+X(2)*X(2))
      IF (R.EQ.0.0D0) THEN
        IF (X(3).GT.0.0D0) THEN
          LAT=0.5D0*PI
        ELSE
          LAT=-0.5D0*PI
        ENDIF
      ELSE
        TLT=X(3)/R
        LAT=ATAN(TLT)
      ENDIF
!
!     --------------------------------------------------------------------
!
      RETURN
      END
!
!     ====================================================================
!
SUBROUTINE readgrid

! To allocate array space for the grid information in module grid
! and to read the information from file

USE grid
USE channels

IMPLICIT NONE

INTEGER :: if0, ie0, iv0, igrid, ix
CHARACTER*127 :: ygridfile

! ----------------------------------------------------------------

CHARACTER*127 :: gridType
CHARACTER*127 :: argGridCells
INTEGER :: gridCells

IF (IARGC() < 2) THEN
  WRITE (*,*) "Usage: ./buildop_fem [cube/hex] [cells]"
  WRITE (*,*) "    or    "
  WRITE (*,*) "Usage: ./buildop_fem [cube/hex] 0 [gridfilename]"
  WRITE (*,*) "       "
  WRITE (*,*) "You may add as 4th argument the name for output."
  CALL EXIT(-1)
END IF

CALL GETARG(1, gridType)
gridType = TRIM(gridType)

IF (gridType == "hex") THEN
  igtype=1
ELSE IF (gridType == "cube") THEN
  igtype = 2
ELSE
  WRITE(*,*) "Unknown grid type"
  CALL EXIT(-1)
END IF

CALL GETARG(2, argGridCells)
READ (argGridCells, '(i10)') nfacex

IF (igtype == 1) THEN
  WRITE(ygridfile,'(''gridmap_hex_'',I10.10,''.dat'')') nfacex
ELSEIF (igtype == 2) THEN
  WRITE(ygridfile,'(''gridmap_cube_'',I10.10,''.dat'')') nfacex
ELSEIF (igtype == 3) THEN
  ! Special case for testing non-orthog version of hex grid
  WRITE(ygridfile,'(''gridmap_hexb_'',I10.10,''.dat'')') nfacex
ELSE
  PRINT *,'Unknown grid type'
  STOP
ENDIF
ygridfile="grd/"//trim(ygridfile)

IF (IARGC() > 2) THEN
    CALL GETARG(3, ygridfile)
    ygridfile=TRIM(ygridfile)
    WRITE(*,*) "using "//trim(ygridfile)//" as input"
END IF

! Open file for reading
OPEN(changin,FILE=ygridfile,FORM='UNFORMATTED')

! First read ngrids
READ(changin) ngrids

! Allocate nface, nedge, nvert
ALLOCATE(nface(ngrids), nedge(ngrids), nvert(ngrids))

! Read numbers of faces, edges and vertices on each grid
READ(changin) nface
READ(changin) nedge
READ(changin) nvert

! Find maximum values in order to allocated subsequent arrays
nfacex = MAXVAL(nface)
nedgex = MAXVAL(nedge)
nvertx = MAXVAL(nvert)

! Allocate neoff, neofv
ALLOCATE(neoff(nfacex,ngrids), neofv(nvertx,ngrids))
neoff = 0
neofv = 0

! Read the numbers of edges of each face and vertex on each grid
READ(changin) ((neoff(if0,igrid),           &
                    if0 = 1, nface(igrid)), &
                    igrid = 1, ngrids)
READ(changin) ((neofv(iv0,igrid),           &
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

! Read the connectivity arrays
READ(changin) (((fnxtf(if0,ix,igrid),           &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, nefmx),            &
                     igrid = 1, ngrids)
READ(changin) (((eoff(ix,if0,igrid),            &
                     ix = 1, nefmx),            &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)
READ(changin) (((voff(if0,ix,igrid),            &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, nefmx),            &
                     igrid = 1, ngrids)
READ(changin) (((fnxte(ix,ie0,igrid),           &
                     ix = 1, 2),                &
                     ie0 = 1, nedge(igrid)),    &
                     igrid = 1, ngrids)
READ(changin) (((vofe(ix,ie0,igrid),            &
                     ix = 1, 2),                &
                     ie0 = 1, nedge(igrid)),    &
                     igrid = 1, ngrids)
READ(changin) (((fofv(iv0,ix,igrid),            &
                     iv0 = 1, nvert(igrid)),    &
                     ix = 1, nevmx),            &
                     igrid = 1, ngrids)
READ(changin) (((eofv(ix,iv0,igrid),            &
                     ix = 1, nevmx),            &
                     iv0 = 1, nvert(igrid)),    &
                     igrid = 1, ngrids)

! Allocate the geometrical information arrays
ALLOCATE(flong(nfacex,ngrids), flat(nfacex,ngrids), &
         vlong(nvertx,ngrids), vlat(nvertx,ngrids), &
         elong(nedgex,ngrids), elat(nedgex,ngrids), &
         farea(nfacex,ngrids), varea(nvertx,ngrids), &
         ldist(nedgex,ngrids), ddist(nedgex,ngrids))

! Read the geometrical information arrays
READ(changin) ((flong(if0,igrid),               &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)
READ(changin) ((flat(if0,igrid),                &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)
READ(changin) ((vlong(iv0,igrid),               &
                     iv0 = 1, nvert(igrid)),    &
                     igrid = 1, ngrids)
READ(changin) ((vlat(iv0,igrid),                &
                     iv0 = 1, nvert(igrid)),    &
                     igrid = 1, ngrids)
READ(changin) ((farea(if0,igrid),               &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)
READ(changin) ((ldist(ie0,igrid),               &
                     ie0 = 1, nedge(igrid)),    &
                     igrid = 1, ngrids)
READ(changin) ((ddist(ie0,igrid),               &
                     ie0 = 1, nedge(igrid)),    &
                     igrid = 1, ngrids)

READ(changin) SFCIndexAvailable

IF (SFCIndexAvailable == 1) THEN
    ALLOCATE(fNewFaceId(nfacex, ngrids), fNewFaceIdInverse(nfacex, ngrids), &
		fNewEdgeId(nedgex, ngrids), fNewEdgeIdInverse(nedgex, ngrids), &
		fNewVertId(nvertx, ngrids), fNewVertIdInverse(nvertx, ngrids))

     WRITE(*,*) "Using SFC information for MG"
     READ(changin) ((fNewFaceId(if0,igrid),          &
          if0 = 1, nface(igrid)),    &
          igrid = 1, ngrids)

     READ(changin) ((fNewFaceIdInverse(if0,igrid),          &
          if0 = 1, nface(igrid)),    &
          igrid = 1, ngrids)

     READ(changin) ((fNewEdgeId(if0,igrid),          &
          if0 = 1, nedge(igrid)),    &
          igrid = 1, ngrids)

     READ(changin) ((fNewEdgeIdInverse(if0,igrid),          &
          if0 = 1, nedge(igrid)),    &
          igrid = 1, ngrids)

     READ(changin) ((fNewVertId(if0,igrid),          &
          if0 = 1, nvert(igrid)),    &
          igrid = 1, ngrids)

     READ(changin) ((fNewVertIdInverse(if0,igrid),          &
          if0 = 1, nvert(igrid)),    &
          igrid = 1, ngrids)

     WRITE(*,*) "Done reading SFC information"
END IF

! (elong and elat are computed later, and farea, ldist and ddist are
! recomputed for a faceted grid.)


!     ---------------------------------------------------------------

! Allocate space for defining compound elements
ncvpmx = 2*nefmx
ncspmx = 2 + 4*nefmx
ncepmx = 2*nefmx
ncvdmx = 2*nefmx
ncsdmx = 2 + 4*nevmx
ALLOCATE(ncvp(nfacex), cvp(nfacex,ncvpmx), &
         ncsp(nedgex), csp(nedgex,ncspmx), &
         ncep(nvertx), cep(nvertx,ncepmx), &
         ncvd(nvertx), cvd(nvertx,ncvdmx), &
         ncsd(nedgex), csd(nedgex,ncsdmx))

!     ---------------------------------------------------------------

! Allocate space for operator stencils and coefficients
nlsmx = 1
nmsmx = 2*nefmx - 1
njsmx = nevmx*(nefmx - 2) + 1
nhsmx = 2*(nevmx - 1)*(nefmx - 1) - 1
nrsmx = nefmx
nrxsmx = nevmx
nwsmx = 2*(nefmx - 1)
ntsmx = nefmx
ALLOCATE(nlsten(nfacex,ngrids),lsten(nlsmx,nfacex,ngrids),lmass(nlsmx,nfacex,ngrids), &
         nmsten(nedgex,ngrids),msten(nmsmx,nedgex,ngrids),mmass(nmsmx,nedgex,ngrids), &
         njsten(nvertx,ngrids),jsten(njsmx,nvertx,ngrids),jstar(njsmx,nvertx,ngrids), &
         nhsten(nedgex,ngrids),hsten(nhsmx,nedgex,ngrids),hstar(nhsmx,nedgex,ngrids), &
         nrsten(nfacex,ngrids),rsten(nfacex,nrsmx,ngrids),rcoeff(nfacex,nrsmx,ngrids), &
         nrxsten(nvertx,ngrids),rxsten(nrxsmx,nvertx,ngrids),rxcoeff(nrxsmx,nvertx,ngrids), &
         nwsten(nedgex,ngrids),wsten(nwsmx,nedgex,ngrids),wcoeff(nwsmx,nedgex,ngrids), &
         ntsten(nfacex,ngrids),tsten(ntsmx,nfacex,ngrids),tcoeff(ntsmx,ntsmx,nfacex,ngrids))

!     ---------------------------------------------------------------

END SUBROUTINE readgrid

!     ===============================================================

SUBROUTINE writegrid

! To write the grid and operator information to file

USE grid
USE channels

IMPLICIT NONE

INTEGER :: if0, ie0, iv0, igrid, ix, ixx
CHARACTER*127 :: ygridfile

! ----------------------------------------------------------------

! Construct filename
IF (igtype == 1) THEN
  WRITE(ygridfile,'(''gridopermap_hex_'',I10.10,''.dat'')') nfacex
ELSEIF (igtype == 2) THEN
  WRITE(ygridfile,'(''gridopermap_cube_'',I10.10,''.dat'')') nfacex
ELSEIF (igtype == 3) THEN
  ! Special case for testing non-orthog version of hex grid
  WRITE(ygridfile,'(''gridopermap_hexb_'',I10.10,''.dat'')') nfacex
ENDIF

IF (IARGC() > 3) THEN
    CALL GETARG(4, ygridfile)
    ygridfile=TRIM(ygridfile)
    !WRITE(*,*) "using "//trim(ygridfile)//" as output"
END IF

ygridfile="grd/"//trim(ygridfile)  !P. Peixoto

! Open file for writing
OPEN(changout,FILE=ygridfile,FORM='UNFORMATTED')
WRITE(*,*) "using "//trim(ygridfile)//" as output"

! First write ngrids
WRITE(changout) ngrids

! Write numbers of faces, edges and vertices on each grid
WRITE(changout) nface
WRITE(changout) nedge
WRITE(changout) nvert

! Write the numbers of edges of each face and vertex on each grid
WRITE(changout) ((neoff(if0,igrid),         &
                    if0 = 1, nface(igrid)), &
                    igrid = 1, ngrids)
WRITE(changout) ((neofv(iv0,igrid),         &
                    iv0 = 1, nvert(igrid)), &
                    igrid = 1, ngrids)


! Write the connectivity arrays
WRITE(changout) (((fnxtf(if0,ix,igrid),         &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, nefmx),            &
                     igrid = 1, ngrids)
WRITE(changout) (((eoff(ix,if0,igrid),          &
                     ix = 1, nefmx),            &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) (((voff(if0,ix,igrid),          &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, nefmx),            &
                     igrid = 1, ngrids)
WRITE(changout) (((fnxte(ix,ie0,igrid),         &
                     ix = 1, 2),    &
                     ie0 = 1, nedge(igrid)),                &
                     igrid = 1, ngrids)
WRITE(changout) (((vofe(ix,ie0,igrid),          &
                     ix = 1, 2),    &
                     ie0 = 1, nedge(igrid)),                &
                     igrid = 1, ngrids)
WRITE(changout) (((fofv(iv0,ix,igrid),          &
                     iv0 = 1, nvert(igrid)),    &
                     ix = 1, nevmx),            &
                     igrid = 1, ngrids)
WRITE(changout) (((eofv(ix,iv0,igrid),          &
                     ix = 1, nevmx),            &
                     iv0 = 1, nvert(igrid)),    &
                     igrid = 1, ngrids)


! Write the geometrical information arrays
WRITE(changout) ((flong(if0,igrid),             &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) ((flat(if0,igrid),              &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) ((vlong(iv0,igrid),             &
                     iv0 = 1, nvert(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) ((vlat(iv0,igrid),              &
                     iv0 = 1, nvert(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) ((farea(if0,igrid),             &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) ((varea(iv0,igrid),             &
                     iv0 = 1, nvert(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) ((ldist(ie0,igrid),             &
                     ie0 = 1, nedge(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) ((ddist(ie0,igrid),             &
                     ie0 = 1, nedge(igrid)),    &
                     igrid = 1, ngrids)


! Write the sizes of the operator stencils on each grid
WRITE(changout) ((nlsten(if0,igrid),            &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) ((nmsten(ie0,igrid),            &
                     ie0 = 1, nedge(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) ((njsten(iv0,igrid),            &
                     iv0 = 1, nvert(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) ((nhsten(ie0,igrid),            &
                     ie0 = 1, nedge(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) ((nrsten(if0,igrid),            &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) ((nrxsten(iv0,igrid),           &
                     iv0 = 1, nvert(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) ((nwsten(ie0,igrid),            &
                     ie0 = 1, nedge(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) ((ntsten(if0,igrid),            &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)


! Write the operator stencils and coefficients
WRITE(changout) (((lsten(ix,if0,igrid),         &
                     ix = 1, nlsmx),    &
                     if0 = 1, nface(igrid)),            &
                     igrid = 1, ngrids)
WRITE(changout) (((msten(ix,ie0,igrid),         &
                     ix = 1, nmsmx),    &
                     ie0 = 1, nedge(igrid)),            &
                     igrid = 1, ngrids)
WRITE(changout) (((jsten(ix,iv0,igrid),         &
                     ix = 1, njsmx),    &
                     iv0 = 1, nvert(igrid)),            &
                     igrid = 1, ngrids)
WRITE(changout) (((hsten(ix,ie0,igrid),         &
                     ix = 1, nhsmx),    &
                     ie0 = 1, nedge(igrid)),            &
                     igrid = 1, ngrids)
WRITE(changout) (((rsten(if0,ix,igrid),         &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, nrsmx),            &
                     igrid = 1, ngrids)
WRITE(changout) (((rxsten(ix,iv0,igrid),        &
                     ix = 1, nrxsmx),    &
                     iv0 = 1, nvert(igrid)),           &
                     igrid = 1, ngrids)
WRITE(changout) (((wsten(ix,ie0,igrid),         &
                     ix = 1, nwsmx),    &
                     ie0 = 1, nedge(igrid)),            &
                     igrid = 1, ngrids)
WRITE(changout) (((tsten(ix,if0,igrid),         &
                     ix = 1, ntsmx),            &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) (((lmass(ix,if0,igrid),         &
                     ix = 1, nlsmx),    &
                     if0 = 1, nface(igrid)),            &
                     igrid = 1, ngrids)
WRITE(changout) (((mmass(ix,ie0,igrid),         &
                     ix = 1, nmsmx),    &
                     ie0 = 1, nedge(igrid)),            &
                     igrid = 1, ngrids)
WRITE(changout) (((jstar(ix,iv0,igrid),         &
                     ix = 1, njsmx),    &
                     iv0 = 1, nvert(igrid)),            &
                     igrid = 1, ngrids)
WRITE(changout) (((hstar(ix,ie0,igrid),         &
                     ix = 1, nhsmx),    &
                     ie0 = 1, nedge(igrid)),            &
                     igrid = 1, ngrids)
WRITE(changout) (((rcoeff(if0,ix,igrid),        &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, nrsmx),            &
                     igrid = 1, ngrids)
WRITE(changout) (((rxcoeff(ix,iv0,igrid),       &
                     ix = 1, nrxsmx),    &
                     iv0 = 1, nvert(igrid)),           &
                     igrid = 1, ngrids)
WRITE(changout) (((wcoeff(ix,ie0,igrid),        &
                     ix = 1, nwsmx),    &
                     ie0 = 1, nedge(igrid)),            &
                     igrid = 1, ngrids)
WRITE(changout) ((((tcoeff(ixx,ix,if0,igrid),   &
                     ixx = 1, ntsmx),    &
                     ix = 1, ntsmx),            &
                     if0 = 1, nface(igrid)),           &
                     igrid = 1, ngrids)


! Write the size of the restriction stencil
WRITE(changout) ((nres(if0,igrid),              &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids-1)

! Write the restriction stencil and coefficients
WRITE(changout) (((ressten(ix,if0,igrid),       &
                     ix = 1, nresmx),    &
                     if0 = 1, nface(igrid)),           &
                     igrid = 1, ngrids-1)
WRITE(changout) (((reswgt(ix,if0,igrid),        &
                     ix = 1, nresmx),    &
                     if0 = 1, nface(igrid)),           &
                     igrid = 1, ngrids-1)




! Write the size of the injtriction stencil
WRITE(changout) ((ninj(if0,igrid),              &
                     if0 = 1, nface(igrid)),    &
                     igrid = 2, ngrids)

! Write the injtriction stencil and coefficients
WRITE(changout) (((injsten(ix,if0,igrid),       &
                     ix = 1, ninjmx),    &
                     if0 = 1, nface(igrid)),           &
                     igrid = 2, ngrids)
WRITE(changout) (((injwgt(ix,if0,igrid),        &
                     ix = 1, ninjmx),    &
                     if0 = 1, nface(igrid)),           &
                     igrid = 2, ngrids)


WRITE(changout) SFCIndexAvailable

IF (SFCIndexAvailable == 1) THEN
    WRITE(changout) ((fNewFaceId(if0,igrid),          &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)

    WRITE(changout) ((fNewFaceIdInverse(if0,igrid),          &
                         if0 = 1, nface(igrid)),    &
                         igrid = 1, ngrids)

    WRITE(changout) ((fNewEdgeId(if0,igrid),          &
                     if0 = 1, nedge(igrid)),    &
                     igrid = 1, ngrids)

    WRITE(changout) ((fNewEdgeIdInverse(if0,igrid),          &
                         if0 = 1, nedge(igrid)),    &
                         igrid = 1, ngrids)

    WRITE(changout) ((fNewVertId(if0,igrid),          &
                     if0 = 1, nvert(igrid)),    &
                     igrid = 1, ngrids)

    WRITE(changout) ((fNewVertIdInverse(if0,igrid),          &
                         if0 = 1, nvert(igrid)),    &
                         igrid = 1, ngrids)
END IF

!     ---------------------------------------------------------------

END SUBROUTINE writegrid

!     ===============================================================
