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


END MODULE constants

MODULE runtype

! Default values are set here; they may be overridden
! via the namelist

IMPLICIT NONE

! grid type (cube/hex)
CHARACTER*31 :: gridType = 'cube'

! File containing grid information
! CHARACTER*31 :: ygridfile = 'gridopermap_hex_0000010242.dat'
CHARACTER*127 :: ygridcoords = 'gridmap_cube_0000221184.datxxxx'
CHARACTER*127 :: outygridcoords = 'gridmap_cube_0000221184.datxxxx'

! Run identifier for naming output files
CHARACTER*6 :: runid = '000001'

! True for a restart run
LOGICAL :: lrestart = .FALSE.

! Name of restart file
CHARACTER*30 :: yresfile = 'run000002_restart_00000720.datxxxx'


END MODULE runtype


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

INTEGER, PARAMETER :: chanaout = 90       ! Output channel for adjacency information
INTEGER, PARAMETER :: chanAdjGnuplotOUT = 91       ! Output channel for gnuplot file plotting adjacency information

END MODULE channels


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
    OPEN(changrid,FILE=ygridcoords,FORM='UNFORMATTED')

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
                         ix = 1, 2),    &
                         ie0 = 1, nedge(igrid)),                &
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
    ALLOCATE(flong(nfacex,ngrids), flat(nfacex,ngrids), &
             vlong(nvertx,ngrids), vlat(nvertx,ngrids), &
             elong(nedgex,ngrids), elat(nedgex,ngrids), &
             farea(nfacex,ngrids), varea(nvertx,ngrids), &
             ldist(nedgex,ngrids), ddist(nedgex,ngrids))


    ! Read the geometrical information arrays
    READ(changrid) ((flong(if0,igrid),               &
                         if0 = 1, nface(igrid)),    &
                         igrid = 1, ngrids)
    READ(changrid) ((flat(if0,igrid),                &
                         if0 = 1, nface(igrid)),    &
                         igrid = 1, ngrids)
    READ(changrid) ((vlong(iv0,igrid),               &
                         iv0 = 1, nvert(igrid)),    &
                         igrid = 1, ngrids)
    READ(changrid) ((vlat(iv0,igrid),                &
                         iv0 = 1, nvert(igrid)),    &
                         igrid = 1, ngrids)
    READ(changrid) ((farea(if0,igrid),               &
                         if0 = 1, nface(igrid)),    &
                         igrid = 1, ngrids)
    READ(changrid) ((ldist(ie0,igrid),               &
                         ie0 = 1, nedge(igrid)),    &
                         igrid = 1, ngrids)
    READ(changrid) ((ddist(ie0,igrid),               &
                         ie0 = 1, nedge(igrid)),    &
                         igrid = 1, ngrids)

    READ(changrid) SFCIndexAvailable

    IF (SFCIndexAvailable == 1) THEN
        WRITE(*,*) "SFC optimization already applied!"
        CALL EXIT(1)
    END IF

    CLOSE(changrid)
END SUBROUTINE readgrid





SUBROUTINE writegrid

    USE runtype
    USE constants
    USE grid
    USE channels

    IMPLICIT NONE

    INTEGER :: igrid, i, j, ixv, p1, p2, pp, jr, iv, iv0, n, n2, &
               ie0, ie1, ie2, if0, if1, if2, if3, iv1, iv2, ix1, ix2, &
               ixmin, ifmin, if21, if22, iv11, iv12, iv21, iv22, &
               ismooth

    ! Output gridmap file
    OPEN(changrid,FILE=outygridcoords,FORM='UNFORMATTED')

    ! WRITE(changrid,*) 'GRIDMAP for NGRIDS=',NGRIDS
    WRITE(changrid) ngrids
    WRITE(changrid) nface
    WRITE(changrid) nedge
    WRITE(changrid) nvert

    ! WRITE(changrid,*) 'Number of edges of each face - all grids'
    WRITE(changrid) ((neoff(if1,igrid),            &
                   if1 = 1, nface(igrid)),   &
                   igrid = 1, ngrids)
    ! WRITE(changrid,*) 'Number of edges of each vertex - all grids'
    WRITE(changrid) ((neofv(iv1,igrid),            &
                   iv1 = 1, nvert(igrid)),   &
                   igrid=1, ngrids)
    ! WRITE(changrid,*) 'Faces next to each face - all grids'
    WRITE(changrid) (((fnxtf(if1,if2,igrid),       &
                   if1 = 1, nface(igrid)),   &
                   if2 = 1, nefmx),              &
                   igrid = 1, ngrids)
    ! WRITE(changrid,*) 'Edges of each face - all grids'
    WRITE(changrid) (((eoff(ie1,if1,igrid),        &
                   ie1 = 1, nefmx),                &
                   if1 = 1, nface(igrid)),         &
                   igrid = 1, ngrids)
    ! WRITE(changrid,*) 'Vertices of each face - all grids'
    WRITE(changrid) (((voff(if1,iv1,igrid),        &
                   if1 = 1, nface(igrid)),   &
                   iv1 = 1, nefmx),              &
                   igrid = 1, ngrids)
    ! WRITE(changrid,*) 'Faces next to each edge - all grids'
    WRITE(changrid) (((fnxte(if2,ie1,igrid),       &
                   if2 = 1, 2),   &
                   ie1 = 1, nedge(igrid)),              &
                   igrid = 1, ngrids)
    ! WRITE(changrid,*) 'Vertices of each edge - all grids'
    WRITE(changrid) (((vofe(iv2,ie1,igrid),        &
                   iv2 = 1, 2),   &
                   ie1 = 1, nedge(igrid)),              &
                   igrid = 1, ngrids)
    ! WRITE(changrid,*) 'Faces around each vertex - all grids'
    WRITE(changrid) (((fofv(iv1,if2,igrid),        &
                   iv1 = 1, nvert(igrid)),   &
                   if2 = 1, nevmx),              &
                   igrid = 1, ngrids)
    ! WRITE(changrid,*) 'Edges around each vertex - all grids'
    WRITE(changrid) (((eofv(ie1,iv1,igrid),        &
                   ie1 = 1, nevmx),   &
                   iv1 = 1, nvert(igrid)),              &
                   igrid = 1, ngrids)


    ! WRITE(changrid,*) 'Longitudes of faces - all grids'
    WRITE(changrid) ((flong(if1,igrid),            &
                   if1 = 1, nface(igrid)),   &
                   igrid = 1, ngrids)
    ! WRITE(changrid,*) 'Latitudes of faces - all grids'
    WRITE(changrid) ((flat(if1,igrid),             &
                   if1 = 1, nface(igrid)),   &
                   igrid = 1, ngrids)
    ! WRITE(changrid,*) 'Longitudes of vertices - all grids'
    WRITE(changrid) ((vlong(iv1,igrid),            &
                   iv1 = 1, nvert(igrid)),   &
                   igrid = 1, ngrids)
    ! WRITE(changrid,*) 'Latitudes of vertices - all grids'
    WRITE(changrid) ((vlat(iv1,igrid),             &
                   iv1 = 1, nvert(igrid)),   &
                   igrid = 1, ngrids)
    ! WRITE(changrid,*) 'Areas of faces - all grids'
    WRITE(changrid) ((farea(if1,igrid),            &
                  if1 = 1, nface(igrid)),    &
                  igrid = 1, ngrids)
    ! WRITE(changrid,*) 'Lengths of edges - all grids'
    WRITE(changrid) ((ldist(ie1,igrid),            &
                  ie1 = 1, nedge(igrid)),    &
                  igrid = 1, ngrids)
    ! WRITE(changrid,*) 'Distance between faces across edges - all grids'
    WRITE(changrid) ((ddist(ie1,igrid),            &
                  ie1 = 1, nedge(igrid)),    &
                  igrid = 1, ngrids)

    if1 = 1
    WRITE(changrid) if1

    WRITE(changrid) ((fNewFaceId(if0,igrid),          &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)

    WRITE(changrid) ((fNewFaceIdInverse(if0,igrid),          &
                         if0 = 1, nface(igrid)),    &
                         igrid = 1, ngrids)

    WRITE(changrid) ((fNewEdgeId(if0,igrid),          &
                     if0 = 1, nedge(igrid)),    &
                     igrid = 1, ngrids)

    WRITE(changrid) ((fNewEdgeIdInverse(if0,igrid),          &
                         if0 = 1, nedge(igrid)),    &
                         igrid = 1, ngrids)

    WRITE(changrid) ((fNewVertId(if0,igrid),          &
                     if0 = 1, nvert(igrid)),    &
                     igrid = 1, ngrids)

    WRITE(changrid) ((fNewVertIdInverse(if0,igrid),          &
                         if0 = 1, nvert(igrid)),    &
                         igrid = 1, ngrids)

    CLOSE(changrid)

END SUBROUTINE writegrid



SUBROUTINE LL2XYZ(LONG,LAT,X)
!
!     To convert longitude and latitude to cartesian coordinates
!     on the unit sphere
!
      IMPLICIT NONE
      REAL*8 LONG,LAT,X(3),CLN,SLN,CLT,SLT

      SLN=SIN(LONG)
      CLN=COS(LONG)
      SLT=SIN(LAT)
      CLT=COS(LAT)

      X(1)=CLN*CLT
      X(2)=SLN*CLT
      X(3)=SLT

      RETURN
END



SUBROUTINE computeSFCIndex(X, o_scalar)
    ! compute SFC index based on Z-curve
    ! assume that X(i) \in [0;1]

    IMPLICIT NONE
    REAL*8 :: X(3)
    REAL*8 :: o_scalar
    REAL*8 :: scalarIntervalSize
    INTEGER :: i, d
    REAL*8 :: spaceIntervalSize

    ! return 1D index \in [0;1]
    o_scalar = 0

    ! current search interval size for 1D index
    scalarIntervalSize = 0.5

    ! space interval size in each dimension.
    ! this is useful if iterating over each dimension
    spaceIntervalSize = 0.5

    ! we shift and scale the search space to simplify the operations to

    ! iterate over all possible bit representations
    DO i = 1, 63
        ! iterate over all dimensions
        DO d = 1, 3
            IF (X(d) >= 0.5) THEN
                o_scalar = o_scalar + scalarIntervalSize
                X(d) = X(d) - 0.5
            END IF
            scalarIntervalSize = scalarIntervalSize * 0.5
        END DO
        X = X*2.0
    END DO
END

SUBROUTINE computeSFCIndexH3d(X, o_scalar)

  IMPLICIT NONE
  REAL*8 :: X(3)
  REAL*8 :: o_scalar

  interface
    function encode3d(a,b,c)
      REAL*8 :: encode3d, a, b, c
    end function
  end interface
 
  o_scalar = encode3d(%val(X(1)),%val(X(2)),%val(X(3)))
  !write(*,*) X(1), X(2), X(3), o_scalar

END SUBROUTINE computeSFCIndexH3d

SUBROUTINE computeSFCIndexH2d(x1, x2, o_scalar)

  IMPLICIT NONE
  REAL*8 :: x1, x2
  REAL*8 :: o_scalar

  interface
    function encode2d(a,b)
      REAL*8 :: encode2d, a, b
    end function
  end interface
 
  o_scalar = encode2d(%val(x1),%val(x2))
  !write(*,*) x1, x2, o_scalar

END SUBROUTINE computeSFCIndexH2d

SUBROUTINE sfcOptimize(face, edge, vert)

    USE runtype
    USE grid

    IMPLICIT NONE

    INTEGER, INTENT (IN) :: face, edge, vert

    REAL*8 :: sfcScalar
    REAL*8 :: x0(3)
    INTEGER :: igrid
    INTEGER :: nf, ne, nv
    INTEGER :: if0, if1, if2, oldf0
    INTEGER :: ie0, ie1, ie2, newe0, olde0
    INTEGER :: iv0, iv1, iv2, oldv0, newv0
    INTEGER :: ixn
    INTEGER :: i

    INTEGER, ALLOCATABLE :: SEED(:)
    INTEGER :: K
    REAL, ALLOCATABLE :: rnd(:)

    REAL*8, ALLOCATABLE :: sfcIndex(:), sfcSorting(:)

    INTEGER, ALLOCATABLE :: tmp_neoff(:), tmp_neofv(:)

    INTEGER, ALLOCATABLE :: tmp_fnxtf(:,:), tmp_eoff(:,:), tmp_voff(:,:), &
                            tmp_fnxte(:,:), tmp_vofe(:,:), &
                            tmp_fofv(:,:), tmp_eofv(:,:)

    REAL*8, ALLOCATABLE :: tmp_flong(:), tmp_flat(:), &
                           tmp_vlong(:), tmp_vlat(:), &
                           tmp_farea(:), tmp_varea(:), &
                           tmp_ldist(:), tmp_ddist(:)

    INTEGER, ALLOCATABLE :: visited_vert(:), visited_edge(:), faceIdSorting(:)

    REAL*8 :: long, lat, x1(3), x2(3), y1(3), y2(3), &
              n1(3), n2(3), r1(3), mag

    ! FOR HILBERT 2d
    INTEGER :: rightx, leftx, topx, bottomx, frontx, backx, nextStart, nextEnd
    INTEGER, ALLOCATABLE :: right(:), left(:), top(:), bottom(:), front(:), back(:)
    REAL*8, ALLOCATABLE :: rightSFC(:), leftSFC(:), topSFC(:), bottomSFC(:), frontSFC(:), backSFC(:)

    ALLOCATE(tmp_neoff(nfacex), tmp_neofv(nvertx))

    ! Allocate connectivity arrays arrays
    ALLOCATE(tmp_fnxtf(nfacex,nefmx), tmp_eoff(nefmx,nfacex), &
             tmp_voff(nfacex,nefmx),  tmp_fnxte(2,nedgex),    &
             tmp_vofe(2,nedgex),      tmp_fofv(nvertx,nevmx), &
             tmp_eofv(nevmx,nvertx))


    ! Allocate the geometrical information arrays
    ALLOCATE(tmp_flong(nfacex), tmp_flat(nfacex),  &
             tmp_vlong(nvertx), tmp_vlat(nvertx),  &
             tmp_farea(nfacex), tmp_varea(nvertx), &
             tmp_ldist(nedgex), tmp_ddist(nedgex)   &
        )

    ALLOCATE(fNewFaceId(nfacex, ngrids), fNewFaceIdInverse(nfacex, ngrids), &
		fNewEdgeId(nedgex, ngrids), fNewEdgeIdInverse(nedgex, ngrids), &
		fNewVertId(nvertx, ngrids), fNewVertIdInverse(nvertx, ngrids))

    ! Iterate over all grid levels
    DO igrid = 1, ngrids

        nf = nface(igrid)
        ne = nedge(igrid)
        nv = nvert(igrid)

        WRITE (*,*)
        WRITE (*,*) "Cells on current level: ", nf

IF (face == 0) THEN
	DO if0 = 1, nf
		fNewFaceIdInverse(if0, igrid) = if0
		fNewFaceId(if0, igrid) = if0
	END DO
ELSE

IF (face == 3) THEN

ALLOCATE(right(nf), left(nf), top(nf), bottom(nf), front(nf), back(nf))
ALLOCATE(rightSFC(nf), leftSFC(nf), topSFC(nf), bottomSFC(nf), frontSFC(nf), backSFC(nf))

rightx = 0
leftx = 0
topx = 0
bottomx = 0
frontx = 0
backx = 0

right(:) = 0
left(:) = 0
top(:) = 0
bottom(:) = 0
front(:) = 0
back(:) = 0

DO if0 = 1, nf
	CALL ll2xyz(flong(if0,igrid),flat(if0,igrid),x0)

	IF (ABS(x0(1)) > ABS(x0(2))) THEN
		IF (ABS(x0(1)) > ABS(x0(3))) THEN
			CALL computeSFCIndexH2d((x0(2)+1.0)*0.5, (x0(3)+1.0)*0.5, sfcScalar)
			IF (x0(1) > 0) THEN
				rightx = rightx + 1
				right(rightx) = if0
				rightSFC(rightx) = sfcScalar
			ELSE
				leftx = leftx + 1
				left(leftx) = if0
				leftSFC(leftx) = sfcScalar
			ENDIF	
		ELSE
			CALL computeSFCIndexH2d((x0(1)+1.0)*0.5, (x0(2)+1.0)*0.5, sfcScalar)
			IF (x0(3) > 0) THEN
				frontx = frontx + 1
				front(frontx) = if0
				frontSFC(frontx) = sfcScalar
			ELSE
				backx = backx + 1
				back(backx) = if0
				backSFC(backx) = sfcScalar
			ENDIF		
		END IF
	ELSE
		IF (ABS(x0(2)) > ABS(x0(3))) THEN
			CALL computeSFCIndexH2d((x0(1)+1.0)*0.5, (x0(3)+1.0)*0.5, sfcScalar)
			IF (x0(2) > 0) THEN
				topx = topx + 1
				top(topx) = if0
				topSFC(topx) = sfcScalar
			ELSE
				bottomx = bottomx + 1
				bottom(bottomx) = if0
				bottomSFC(bottomx) = sfcScalar
			ENDIF		
		ELSE
			CALL computeSFCIndexH2d((x0(1)+1.0)*0.5, (x0(2)+1.0)*0.5, sfcScalar)
			IF (x0(3) > 0) THEN
				frontx = frontx + 1
				front(frontx) = if0
				frontSFC(frontx) = sfcScalar
			ELSE
				backx = backx + 1
				back(backx) = if0
				backSFC(backx) = sfcScalar
			ENDIF			
		END IF
	ENDIF

END DO

WRITE(*,*) "Sorting according to SFC indices on each side"
ALLOCATE(sfcSorting((rightx+1)/2), faceIdSorting((rightx+1)/2))
CALL MergeSort(rightSFC, rightx, sfcSorting, right(1:rightx), faceIdSorting)
DEALLOCATE(sfcSorting, faceIdSorting)
ALLOCATE(sfcSorting((leftx+1)/2), faceIdSorting((leftx+1)/2))
CALL MergeSort(leftSFC, leftx, sfcSorting, left(1:leftx), faceIdSorting)
DEALLOCATE(sfcSorting, faceIdSorting)
ALLOCATE(sfcSorting((topx+1)/2), faceIdSorting((topx+1)/2))
CALL MergeSort(topSFC, topx, sfcSorting, top(1:topx), faceIdSorting)
DEALLOCATE(sfcSorting, faceIdSorting)
ALLOCATE(sfcSorting((bottomx+1)/2), faceIdSorting((bottomx+1)/2))
CALL MergeSort(bottomSFC, bottomx, sfcSorting, bottom(1:bottomx), faceIdSorting)
DEALLOCATE(sfcSorting, faceIdSorting)
ALLOCATE(sfcSorting((frontx+1)/2), faceIdSorting((frontx+1)/2))
CALL MergeSort(frontSFC, frontx, sfcSorting, front(1:frontx), faceIdSorting)
DEALLOCATE(sfcSorting, faceIdSorting)
ALLOCATE(sfcSorting((backx+1)/2), faceIdSorting((backx+1)/2))
CALL MergeSort(backSFC, backx, sfcSorting, back(1:backx), faceIdSorting)
DEALLOCATE(sfcSorting, faceIdSorting)

WRITE (*,*) left(:)

nextStart = 1
nextEnd = backx
fNewFaceId(nextStart:nextEnd,igrid) = back(1:backx)
nextStart = nextEnd + 1
nextEnd = nextEnd + topx
fNewFaceId(nextStart:nextEnd,igrid) = top(1:topx)
nextStart = nextEnd + 1
nextEnd = nextEnd + rightx
fNewFaceId(nextStart:nextEnd,igrid) = right(1:rightx)
nextStart = nextEnd + 1
nextEnd = nextEnd + frontx
fNewFaceId(nextStart:nextEnd,igrid) = front(1:frontx)
nextStart = nextEnd + 1
nextEnd = nextEnd + bottomx
fNewFaceId(nextStart:nextEnd,igrid) = bottom(1:bottomx)
nextStart = nextEnd + 1
nextEnd = nextEnd + leftx
fNewFaceId(nextStart:nextEnd,igrid) = left(1:leftx)

DEALLOCATE(right, left, top, bottom, front, back, rightSFC, leftSFC, topSFC, bottomSFC, frontSFC, backSFC)
ENDIF


IF ((face == 1) .OR. (face == 2) .OR. (face == 4)) THEN
        ALLOCATE(sfcIndex(nf),rnd(nf),SEED(nf))
!        ALLOCATE(faceId(nf))

IF (face == 4) THEN
	
	K = nf
	SEED(1:K) = 42 + igrid
	CALL RANDOM_SEED
	CALL RANDOM_SEED(SIZE=K)
	CALL RANDOM_SEED(PUT=SEED(1:K))

	CALL RANDOM_NUMBER(rnd)
	!WRITE(*,*) ' RND : ', rnd

	sfcIndex(:) = rnd(:)
	DO if0 = 1, nf
            fNewFaceId(if0, igrid) = if0
	END DO
ELSE
        WRITE(*,*) "Computing SFC 1D indices"
        ! Iterate over all cells of the current grid level igrid
        DO if0 = 1, nf
            ! initialize face id for reverse lookup after ordering cells along SFC
            fNewFaceId(if0, igrid) = if0

            ! compute 3D position
            CALL ll2xyz(flong(if0,igrid),flat(if0,igrid),x0)

            ! compute SFC index
	    IF (face == 1) THEN
		CALL computeSFCIndex((x0+1.0)*0.5, sfcScalar)
	    ENDIF
	    IF (face == 2) THEN
		CALL computeSFCIndexH3d((x0+1.0)*0.5, sfcScalar)
	    ENDIF

            ! store SFC scalar
            sfcIndex(if0) = sfcScalar
        END DO
END IF
	!WRITE(*,*) " array: ", sfcIndex
        WRITE(*,*) "Sorting according to SFC indices"
	ALLOCATE(sfcSorting((nf+1)/2), faceIdSorting((nf+1)/2))
	CALL MergeSort(sfcIndex, nf, sfcSorting, fNewFaceId(1:nf,igrid), faceIdSorting)
	!WRITE(*,*) " array: ", sfcIndex
	DEALLOCATE(sfcSorting, faceIdSorting)
	        DEALLOCATE(sfcIndex,rnd,SEED)
        ! SORT LIST
        ! TODO: optimize me, e.g. by using a merge sort!
!        DO if0 = 1, nf
!            DO if1 = if0+1, nf
!                IF (sfcIndex(if0) > sfcIndex(if1)) THEN
!                    sfcScalar = sfcIndex(if0)
!                    sfcIndex(if0) = sfcIndex(if1)
!                    sfcIndex(if1) = sfcScalar

!                    i = fNewFaceId(if0, igrid)
!                    fNewFaceId(if0, igrid) = fNewFaceId(if1, igrid)
!                    fNewFaceId(if1, igrid) = i
!                END IF
!            END DO
!        END DO
END IF

        WRITE(*,*) "Updating indices"

            DO if0 = 1, nf
                fNewFaceIdInverse(fNewFaceId(if0, igrid), igrid) = if0
            END DO


            DO if0 = 1, nf
                oldf0 = fNewFaceId(if0, igrid)

                ! look the data in the old storage
                tmp_neoff(if0) = neoff(oldf0, igrid)

                tmp_flong(if0) = flong(oldf0, igrid)
                tmp_flat(if0) = flat(oldf0, igrid)
                tmp_farea(if0) = farea(oldf0, igrid)

                ! iterate over all adjacent cells, given by the number of adjacent edges
                DO ixn = 1, neoff(oldf0,igrid)
                    ! adjacent cell id
                    ! 1st index: cell to search adjacent cells for
                    ! 2nd index: adjacent cell
                    ! 3rd index: grid level
                    tmp_fnxtf(if0, ixn) = fNewFaceIdInverse(fnxtf(oldf0,ixn,igrid), igrid)
                END DO
                tmp_eoff(:, if0) = eoff(:,oldf0,igrid)
                tmp_voff(if0, :) = voff(oldf0,:,igrid)
            END DO

            DO ie0 = 1, ne
                tmp_fnxte(1,ie0) = fNewFaceIdInverse(fnxte(1, ie0, igrid), igrid)
                tmp_fnxte(2,ie0) = fNewFaceIdInverse(fnxte(2, ie0, igrid), igrid)
            END DO

            DO iv0 = 1, nv
                ! iterate over cells adjacent to vertex
                DO ixn = 1, neofv(iv0,igrid)
                    tmp_fofv(iv0,ixn) = fNewFaceIdInverse(fofv(iv0, ixn, igrid), igrid)
                END DO
            END DO

            neoff(:,igrid) = tmp_neoff

            flong(:,igrid) = tmp_flong
            flat(:,igrid) = tmp_flat
            farea(:,igrid) = tmp_farea

            fnxtf(:,:,igrid) = tmp_fnxtf
            eoff(:,:,igrid) = tmp_eoff
            voff(:,:,igrid) = tmp_voff

            fnxte(:,:,igrid) = tmp_fnxte
            fofv(:,:,igrid) = tmp_fofv

END IF	

	!!! EDGE
	WRITE (*,*)
        WRITE (*,*) "Edges on current level: ", ne

IF (edge == 0) THEN
	DO ie0 = 1, ne
		fNewEdgeIdInverse(ie0, igrid) = ie0
		fNewEdgeId(ie0, igrid) = ie0
	END DO
ELSE


IF (edge == 3) THEN

ALLOCATE(right(ne), left(ne), top(ne), bottom(ne), front(ne), back(ne))
ALLOCATE(rightSFC(ne), leftSFC(ne), topSFC(ne), bottomSFC(ne), frontSFC(ne), backSFC(ne))

rightx = 0
leftx = 0
topx = 0
bottomx = 0
frontx = 0
backx = 0

right(:) = 0
left(:) = 0
top(:) = 0
bottom(:) = 0
front(:) = 0
back(:) = 0

DO ie0 = 1, ne

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

	IF (ABS(r1(1)) > ABS(r1(2))) THEN
		IF (ABS(r1(1)) > ABS(r1(3))) THEN
			CALL computeSFCIndexH2d((r1(2)+1.0)*0.5, (r1(3)+1.0)*0.5, sfcScalar)
			IF (r1(1) > 0) THEN
				rightx = rightx + 1
				right(rightx) = ie0
				rightSFC(rightx) = sfcScalar
			ELSE
				leftx = leftx + 1
				left(leftx) = ie0
				leftSFC(leftx) = sfcScalar
			ENDIF	
		ELSE
			CALL computeSFCIndexH2d((r1(1)+1.0)*0.5, (r1(2)+1.0)*0.5, sfcScalar)
			IF (r1(3) > 0) THEN
				frontx = frontx + 1
				front(frontx) = ie0
				frontSFC(frontx) = sfcScalar
			ELSE
				backx = backx + 1
				back(backx) = ie0
				backSFC(backx) = sfcScalar
			ENDIF		
		END IF
	ELSE
		IF (ABS(r1(2)) > ABS(r1(3))) THEN
			CALL computeSFCIndexH2d((r1(1)+1.0)*0.5, (r1(3)+1.0)*0.5, sfcScalar)
			IF (r1(2) > 0) THEN
				topx = topx + 1
				top(topx) = ie0
				topSFC(topx) = sfcScalar
			ELSE
				bottomx = bottomx + 1
				bottom(bottomx) = ie0
				bottomSFC(bottomx) = sfcScalar
			ENDIF		
		ELSE
			CALL computeSFCIndexH2d((r1(1)+1.0)*0.5, (r1(2)+1.0)*0.5, sfcScalar)
			IF (r1(3) > 0) THEN
				frontx = frontx + 1
				front(frontx) = ie0
				frontSFC(frontx) = sfcScalar
			ELSE
				backx = backx + 1
				back(backx) = ie0
				backSFC(backx) = sfcScalar
			ENDIF			
		END IF
	ENDIF

END DO

WRITE(*,*) "Sorting according to SFC indices on each side"
ALLOCATE(sfcSorting((rightx+1)/2), faceIdSorting((rightx+1)/2))
CALL MergeSort(rightSFC, rightx, sfcSorting, right(1:rightx), faceIdSorting)
DEALLOCATE(sfcSorting, faceIdSorting)
ALLOCATE(sfcSorting((leftx+1)/2), faceIdSorting((leftx+1)/2))
CALL MergeSort(leftSFC, leftx, sfcSorting, left(1:leftx), faceIdSorting)
DEALLOCATE(sfcSorting, faceIdSorting)
ALLOCATE(sfcSorting((topx+1)/2), faceIdSorting((topx+1)/2))
CALL MergeSort(topSFC, topx, sfcSorting, top(1:topx), faceIdSorting)
DEALLOCATE(sfcSorting, faceIdSorting)
ALLOCATE(sfcSorting((bottomx+1)/2), faceIdSorting((bottomx+1)/2))
CALL MergeSort(bottomSFC, bottomx, sfcSorting, bottom(1:bottomx), faceIdSorting)
DEALLOCATE(sfcSorting, faceIdSorting)
ALLOCATE(sfcSorting((frontx+1)/2), faceIdSorting((frontx+1)/2))
CALL MergeSort(frontSFC, frontx, sfcSorting, front(1:frontx), faceIdSorting)
DEALLOCATE(sfcSorting, faceIdSorting)
ALLOCATE(sfcSorting((backx+1)/2), faceIdSorting((backx+1)/2))
CALL MergeSort(backSFC, backx, sfcSorting, back(1:backx), faceIdSorting)
DEALLOCATE(sfcSorting, faceIdSorting)

nextStart = 1
nextEnd = backx
fNewEdgeId(nextStart:nextEnd,igrid) = back(1:backx)
nextStart = nextEnd + 1
nextEnd = nextEnd + topx
fNewEdgeId(nextStart:nextEnd,igrid) = top(1:topx)
nextStart = nextEnd + 1
nextEnd = nextEnd + rightx
fNewEdgeId(nextStart:nextEnd,igrid) = right(1:rightx)
nextStart = nextEnd + 1
nextEnd = nextEnd + frontx
fNewEdgeId(nextStart:nextEnd,igrid) = front(1:frontx)
nextStart = nextEnd + 1
nextEnd = nextEnd + bottomx
fNewEdgeId(nextStart:nextEnd,igrid) = bottom(1:bottomx)
nextStart = nextEnd + 1
nextEnd = nextEnd + leftx
fNewEdgeId(nextStart:nextEnd,igrid) = left(1:leftx)

DEALLOCATE(right, left, top, bottom, front, back, rightSFC, leftSFC, topSFC, bottomSFC, frontSFC, backSFC)
ENDIF


ALLOCATE(visited_edge(ne), sfcIndex(ne), rnd(ne), SEED(ne))
IF ((edge == 1) .OR. (edge == 2) .OR. (edge == 4)) THEN

IF (edge == 4) THEN
	
	K = ne
	SEED(1:K) = 12 + igrid
	CALL RANDOM_SEED
	CALL RANDOM_SEED(SIZE=K)
	CALL RANDOM_SEED(PUT=SEED(1:K))

	CALL RANDOM_NUMBER(rnd)
	!WRITE(*,*) ' RND : ', rnd

	sfcIndex(:) = rnd(:)
	DO ie0 = 1, ne
            fNewEdgeId(ie0, igrid) = ie0
	END DO
ELSE

  ! Loop over edges
  DO ie0 = 1, ne

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
    
            ! initialize face id for reverse lookup after ordering cells along SFC
            fNewEdgeId(ie0, igrid) = ie0

            ! compute SFC index
	    IF (edge == 1) THEN
			CALL computeSFCIndex((r1+1.0)*0.5, sfcScalar)
	    ENDIF
	    IF (edge == 2) THEN
		CALL computeSFCIndexH3d((r1+1.0)*0.5, sfcScalar)
	    ENDIF

            ! store SFC scalar
            sfcIndex(ie0) = sfcScalar
    END DO
END IF
	
	WRITE(*,*) "Sorting according to SFC indices"
	ALLOCATE(sfcSorting((ne+1)/2), faceIdSorting((ne+1)/2))
	CALL MergeSort(sfcIndex, ne, sfcSorting, fNewEdgeId(1:ne,igrid), faceIdSorting)
	!WRITE(*,*) " array: ", sfcIndex
	DEALLOCATE(sfcSorting, faceIdSorting)
ENDIF
IF (edge == 5) THEN

	visited_edge = -1
	newe0 = 1
	DO if0 = 1, nf
		DO ixn = 1, neoff(if0, igrid)
		
		ie0 = eoff(ixn, if0, igrid)
		IF (visited_edge(ie0) == -1) THEN
			visited_edge(ie0) = 1
			fNewEdgeId(newe0,igrid) = ie0
			newe0 = newe0 + 1
		ENDIF

		END DO
	END DO
END IF
        WRITE(*,*) "Updating indices for edges"

            DO ie0 = 1, ne
                fNewEdgeIdInverse(fNewEdgeId(ie0, igrid), igrid) = ie0
            END DO

            DO ie0 = 1, ne
                olde0 = fNewEdgeId(ie0, igrid)

		tmp_fnxte(:,ie0) = fnxte(:,olde0,igrid)
		tmp_vofe(:,ie0) = vofe(:,olde0,igrid)
            END DO
	   
 	    DO if0 = 1, nf
		DO ixn = 1, neoff(if0,igrid)
                    tmp_eoff(ixn,if0) = fNewEdgeIdInverse(eoff(ixn,if0,igrid), igrid)
                END DO
	    END DO
	    DO iv0 = 1, nv
		DO ixn = 1, neofv(iv0, igrid)
                    tmp_eofv(ixn, iv0) = fNewEdgeIdInverse(eofv(ixn, iv0, igrid), igrid)
                END DO
	    END DO

	eofv(:,:,igrid) = tmp_eofv		
	eoff(:,:,igrid) = tmp_eoff
	vofe(:,:,igrid) = tmp_vofe
	fnxte(:,:,igrid) = tmp_fnxte

	DEALLOCATE(visited_edge, sfcIndex, rnd, SEED)
END IF

	!!! VERTEX
	WRITE (*,*)
        WRITE (*,*) "Vertex on current level: ", nv

IF (vert == 0) THEN
	DO iv0 = 1, nv
		fNewVertIdInverse(iv0, igrid) = iv0
		fNewVertId(iv0, igrid) = iv0
	END DO
ELSE
        ALLOCATE(sfcIndex(nv), visited_vert(nv), rnd(nv), SEED(nv))
IF ((vert == 1) .OR. (vert == 2) .OR. (vert == 4)) THEN

IF (vert == 4) THEN
	
	K = nv
	SEED(1:K) = 29 + igrid
	CALL RANDOM_SEED
	CALL RANDOM_SEED(SIZE=K)
	CALL RANDOM_SEED(PUT=SEED(1:K))

	CALL RANDOM_NUMBER(rnd)
	!WRITE(*,*) ' RND : ', rnd

	sfcIndex(:) = rnd(:)
	DO iv0 = 1, nv
            fNewVertId(iv0, igrid) = iv0
	END DO
ELSE

WRITE(*,*) "Computing SFC 1D indices for Vertex"
        ! Iterate over all vertex of the current grid level igrid
        DO iv0 = 1, nv
            ! initialize vertex id for reverse lookup after ordering cells along SFC
            fNewVertId(iv0, igrid) = iv0
            ! compute 3D position
            CALL ll2xyz(vlong(iv0,igrid),vlat(iv0,igrid),x0)

            ! compute SFC index
	    IF (vert == 1) THEN
		CALL computeSFCIndex((x0+1.0)*0.5, sfcScalar)
	    ENDIF
	    IF (vert == 2) THEN
		CALL computeSFCIndexH3d((x0+1.0)*0.5, sfcScalar)
	    ENDIF

            ! store SFC scalar
            sfcIndex(iv0) = sfcScalar
        END DO
END IF

	WRITE(*,*) "Sorting according to SFC indices"
	ALLOCATE(sfcSorting((nv+1)/2), faceIdSorting((nv+1)/2))
	CALL MergeSort(sfcIndex, nv, sfcSorting, fNewVertId(1:nv,igrid), faceIdSorting)
	!WRITE(*,*) " array: ", sfcIndex
	DEALLOCATE(sfcSorting, faceIdSorting)
ENDIF

IF (vert == 3) THEN

ALLOCATE(right(nv), left(nv), top(nv), bottom(nv), front(nv), back(nv))
ALLOCATE(rightSFC(nv), leftSFC(nv), topSFC(nv), bottomSFC(nv), frontSFC(nv), backSFC(nv))

rightx = 0
leftx = 0
topx = 0
bottomx = 0
frontx = 0
backx = 0

DO iv0 = 1, nv
	CALL ll2xyz(vlong(iv0,igrid),vlat(iv0,igrid),x0)

	IF (ABS(x0(1)) > ABS(x0(2))) THEN
		IF (ABS(x0(1)) > ABS(x0(3))) THEN
			CALL computeSFCIndexH2d((x0(2)+1.0)*0.5, (x0(3)+1.0)*0.5, sfcScalar)
			IF (x0(1) > 0) THEN
				rightx = rightx + 1
				right(rightx) = iv0
				rightSFC(rightx) = sfcScalar
			ELSE
				leftx = leftx + 1
				left(leftx) = iv0
				leftSFC(leftx) = sfcScalar
			ENDIF	
		ELSE
			CALL computeSFCIndexH2d((x0(1)+1.0)*0.5, (x0(2)+1.0)*0.5, sfcScalar)
			IF (x0(3) > 0) THEN
				frontx = frontx + 1
				front(frontx) = iv0
				frontSFC(frontx) = sfcScalar
			ELSE
				backx = backx + 1
				back(backx) = iv0
				backSFC(backx) = sfcScalar
			ENDIF		
		END IF
	ELSE
		IF (ABS(x0(2)) > ABS(x0(3))) THEN
			CALL computeSFCIndexH2d((x0(1)+1.0)*0.5, (x0(3)+1.0)*0.5, sfcScalar)
			IF (x0(2) > 0) THEN
				topx = topx + 1
				top(topx) = iv0
				topSFC(topx) = sfcScalar
			ELSE
				bottomx = bottomx + 1
				bottom(bottomx) = iv0
				bottomSFC(bottomx) = sfcScalar
			ENDIF		
		ELSE
			CALL computeSFCIndexH2d((x0(1)+1.0)*0.5, (x0(2)+1.0)*0.5, sfcScalar)
			IF (x0(3) > 0) THEN
				frontx = frontx + 1
				front(frontx) = iv0
				frontSFC(frontx) = sfcScalar
			ELSE
				backx = backx + 1
				back(backx) = iv0
				backSFC(backx) = sfcScalar
			ENDIF			
		END IF
	ENDIF

END DO

WRITE(*,*) "Sorting according to SFC indices on each side"
ALLOCATE(sfcSorting((rightx+1)/2), faceIdSorting((rightx+1)/2))
CALL MergeSort(rightSFC, rightx, sfcSorting, right(1:rightx), faceIdSorting)
DEALLOCATE(sfcSorting, faceIdSorting)
ALLOCATE(sfcSorting((leftx+1)/2), faceIdSorting((leftx+1)/2))
CALL MergeSort(leftSFC, leftx, sfcSorting, left(1:leftx), faceIdSorting)
DEALLOCATE(sfcSorting, faceIdSorting)
ALLOCATE(sfcSorting((topx+1)/2), faceIdSorting((topx+1)/2))
CALL MergeSort(topSFC, topx, sfcSorting, top(1:topx), faceIdSorting)
DEALLOCATE(sfcSorting, faceIdSorting)
ALLOCATE(sfcSorting((bottomx+1)/2), faceIdSorting((bottomx+1)/2))
CALL MergeSort(bottomSFC, bottomx, sfcSorting, bottom(1:bottomx), faceIdSorting)
DEALLOCATE(sfcSorting, faceIdSorting)
ALLOCATE(sfcSorting((frontx+1)/2), faceIdSorting((frontx+1)/2))
CALL MergeSort(frontSFC, frontx, sfcSorting, front(1:frontx), faceIdSorting)
DEALLOCATE(sfcSorting, faceIdSorting)
ALLOCATE(sfcSorting((backx+1)/2), faceIdSorting((backx+1)/2))
CALL MergeSort(backSFC, backx, sfcSorting, back(1:backx), faceIdSorting)
DEALLOCATE(sfcSorting, faceIdSorting)

nextStart = 1
nextEnd = backx
fNewVertId(nextStart:nextEnd,igrid) = back(1:backx)
nextStart = nextEnd + 1
nextEnd = nextEnd + topx
fNewVertId(nextStart:nextEnd,igrid) = top(1:topx)
nextStart = nextEnd + 1
nextEnd = nextEnd + rightx
fNewVertId(nextStart:nextEnd,igrid) = right(1:rightx)
nextStart = nextEnd + 1
nextEnd = nextEnd + frontx
fNewVertId(nextStart:nextEnd,igrid) = front(1:frontx)
nextStart = nextEnd + 1
nextEnd = nextEnd + bottomx
fNewVertId(nextStart:nextEnd,igrid) = bottom(1:bottomx)
nextStart = nextEnd + 1
nextEnd = nextEnd + leftx
fNewVertId(nextStart:nextEnd,igrid) = left(1:leftx)


DEALLOCATE(right, left, top, bottom, front, back, rightSFC, leftSFC, topSFC, bottomSFC, frontSFC, backSFC)
ENDIF

IF (vert == 5) THEN
	  
	visited_vert = -1
	newv0 = 1
	  DO if0 = 1, nf
		DO ixn = 1, neoff(if0, igrid)
		
		iv0 = voff(if0, ixn, igrid)
		IF (visited_vert(iv0) == -1) THEN
			visited_vert(iv0) = 1
			fNewVertId(newv0,igrid) = iv0
			newv0 = newv0 + 1
		ENDIF

		END DO
	  END DO
ENDIF

        WRITE(*,*) "Updating indices"

            DO iv0 = 1, nv
                fNewVertIdInverse(fNewVertId(iv0, igrid), igrid) = iv0
            END DO
        WRITE(*,*) "HEREd"
            DO iv0 = 1, nv
                oldv0 = fNewVertId(iv0, igrid)

                ! look the data in the old storage
                tmp_neofv(iv0) = neofv(oldv0, igrid)

                tmp_vlong(iv0) = vlong(oldv0, igrid)
                tmp_vlat(iv0) = vlat(oldv0, igrid)
                tmp_varea(iv0) = varea(oldv0, igrid)

		tmp_fofv(iv0,:) = fofv(oldv0,:, igrid)
		tmp_eofv(:,iv0) = eofv(:,oldv0, igrid)
            END DO
	   
 	    DO if0 = 1, nf
		DO ixn = 1, neoff(if0,igrid)
                    tmp_voff(if0, ixn) = fNewVertIdInverse(voff(if0, ixn, igrid), igrid)
                END DO
	    END DO

	    DO ie0 = 1, ne
		DO ixn = 1, 2
                    tmp_vofe(ixn, ie0) = fNewVertIdInverse(vofe(ixn, ie0, igrid), igrid)
                END DO
	    END DO

	neofv(:,igrid) = tmp_neofv		
	
	vlong(:,igrid) = tmp_vlong
	vlat(:,igrid) = tmp_vlat
	varea(:,igrid) = tmp_varea

	fofv(:,:,igrid) = tmp_fofv
	eofv(:,:,igrid) = tmp_eofv

	voff(:,:,igrid) = tmp_voff
	vofe(:,:,igrid) = tmp_vofe	
	
	DEALLOCATE(sfcIndex, visited_vert, rnd, SEED)
ENDIF

END DO


    ! Deallocate the geometrical information arrays
    DEALLOCATE( &
            tmp_flong, tmp_flat,  &
            tmp_vlong, tmp_vlat,  &
            tmp_farea, tmp_varea, &
            tmp_ldist, tmp_ddist    &
        )


    ! Deallocate connectivity arrays arrays
    DEALLOCATE( &
            tmp_fnxtf,   tmp_eoff, &
            tmp_voff,      tmp_fnxte,    &
            tmp_vofe,      tmp_fofv, &
            tmp_eofv)


    DEALLOCATE(tmp_neoff, tmp_neofv)

END SUBROUTINE sfcOptimize

! Merge and Mergesort taken from http://rosettacode.org/wiki/Sorting_algorithms/Merge_sort#Fortran
! with modifications to fit this problem
subroutine Merge(A,NA,B,NB,C,NC,idA,idB,idC)
 
   integer, intent(in) :: NA,NB,NC         ! Normal usage: NA+NB = NC
   REAL*8, intent(in out) :: A(NA)        ! B overlays C(NA+1:NC)
   REAL*8, intent(in)     :: B(NB)
   REAL*8, intent(in out) :: C(NC)
 
   INTEGER, intent(in out) :: idA(NA)        ! B overlays C(NA+1:NC)
   INTEGER, intent(in)     :: idB(NB)
   INTEGER, intent(in out) :: idC(NC)

   integer :: I,J,K
 
   I = 1; J = 1; K = 1;
   do while(I <= NA .and. J <= NB)
      if (A(I) <= B(J)) then
         C(K) = A(I)
	 idC(K) = idA(I)
         I = I+1
      else
         C(K) = B(J)
	 idC(K) = idB(J)
         J = J+1
      endif
      K = K + 1
   enddo
   do while (I <= NA)
      C(K) = A(I)
      idC(K) = idA(I)
      I = I + 1
      K = K + 1
   enddo
   return
 
end subroutine merge
 
recursive subroutine MergeSort(A,N,T,idA,idT)
 
   integer, intent(in) :: N
   REAL*8, dimension(N), intent(in out) :: A
   REAL*8, dimension((N+1)/2), intent (out) :: T
 
   INTEGER, dimension(N), intent(in out) :: idA
   INTEGER, dimension((N+1)/2), intent (out) :: idT

   integer :: NA,NB, idV
   REAL*8 :: V	 


   if (N < 2) return
   if (N == 2) then
      if (A(1) > A(2)) then
         V = A(1)
         A(1) = A(2)
         A(2) = V
	 idV = idA(1)
	 idA(1) = idA(2)
	 idA(2) = idV
      endif
      return
   endif      
   NA=(N+1)/2
   NB=N-NA
 
   call MergeSort(A,NA,T,idA,idT)
   call MergeSort(A(NA+1),NB,T,idA(NA+1),idT)
 
   if (A(NA) > A(NA+1)) then
      T(1:NA)=A(1:NA)
      idT(1:NA)=idA(1:NA)
      call Merge(T,NA,A(NA+1),NB,A,N,idT,idA(NA+1),idA)
   endif
   return
 
end subroutine MergeSort

PROGRAM sfc_optimize

    USE runtype

    IMPLICIT NONE

    CHARACTER*127 :: argGridCells, argFace, argEdge, argVert
    INTEGER :: gridCells, face, edge, vert


    IF (IARGC() < 2) THEN
        WRITE (*,*) "Usage: ./sfc_optimize [cube/hex] [cells] [face_type] [edge_type] [vert_type]"
        CALL EXIT(-1)
    END IF

    CALL GETARG(1, gridType)
    gridType = TRIM(gridType)
    CALL GETARG(2, argGridCells)
    READ (argGridCells, '(i10)') gridCells

    WRITE (ygridcoords, '(''gridmap_'',A,''_'',I10.10,''.dat'')') TRIM(gridType), gridCells
    WRITE (*,*) ygridcoords

    WRITE (*,*) "Using grid file: "//ygridcoords

    CALL GETARG(3, argFace)
    READ (argFace, '(i1)') face
    CALL GETARG(4, argEdge)
    READ (argEdge, '(i1)') edge
    CALL GETARG(5, argVert)
    READ (argVert, '(i1)') vert

    outygridcoords = ygridcoords
    IF (IARGC() > 5) THEN
       CALL GETARG(6, outygridcoords)
    END IF

    WRITE (*,*) "Using grid output file: "//outygridcoords


    WRITE (*,*) "Reading grid data"
    CALL readgrid

    WRITE (*,*) "SFC Optimization of cell data"
    CALL sfcOptimize(face, edge, vert)

    WRITE (*,*) "Writing grid data"
    CALL writegrid

END PROGRAM
