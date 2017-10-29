
MODULE grid

IMPLICIT NONE


! Grid type (just used to determine file names)
! igtype       1 for hex, 2 for cube
INTEGER :: igtype


! COUNTING

! Number of grids in multigrid hierarchy
INTEGER :: ngrids

! nface        Number of faces on each grid
! nedge        Number of edges on each grid
! nvert        Number of vertices on each edge 
INTEGER, ALLOCATABLE :: nface(:), nedge(:), nvert(:)
INTEGER :: nfacex, nedgex, nvertx

! neoff        Number of edges and vertices of each face on each grid
! neofv        Number of edges of each vertex on each grid
INTEGER, ALLOCATABLE :: neoff(:,:), neofv(:,:)
INTEGER :: nefmx, nevmx


! CONNECTIVITY

! fnxtf        Faces next to each face on each grid
! eoff         Edges of each face on each grid
! voff         Vertices of each face on each grid
! fnxte        Faces either side of each edge on each grid
! vofe         Vertices at the ends of each edge on each grid
! fofv         Faces around each vertex on each grid
! eofv         Edges incident on each vertex on each grid
INTEGER, ALLOCATABLE :: fnxtf(:,:,:), eoff(:,:,:), voff(:,:,:), &
                        fnxte(:,:,:), vofe(:,:,:), &
                        fofv(:,:,:), eofv(:,:,:)

! Conventions
!
! 1. fnxtf and eoff are ordered anticlockwise and such that
! the i'th edge lies between the face in question and its i'th
! neighbour.
!
! 2. voff are ordered anticlockwise such that the k'th vertex
! is the common vertex of the k'th and (k+1)'th edge.
!
! 3. eofv are ordered anticlockwise.
!
! 4. fofv are ordered anticlockwise such that the k'th face lies
! between the k'th and (k+1)'th edge.
!
! 5. The positive normal direction n points from
! fnxte(:,1,:) -> fnxte(:,2,:)
! and the positive tangential direction t points from
! vofe(:,1,:) -> vofe(:,2,:)
! such that t = k x n (where k is vertical).
!
! These conventions are imposed by SUBROUTINE ordering in case they
! are not satisfied by the input data.


! eoffin(f,j,:)   Indicates whether the normal at the j'th edge is
!                 inward or outward relative to face f.
! eofvin(v,j,:)   Indicates whether the tangent at the j'th edge is
!                 inward or outward relative to vertex v.
INTEGER, ALLOCATABLE :: eoffin(:,:,:), eofvin(:,:,:)


! COORDINATES AND GEOMETRICAL INFORMATION

! flong        Longitude of faces on each grid
! flat         Latitude of faces on each grid
! vlong        Longitude of vertices on each grid
! vlat         Latitude of vertices on each grid
! elong        Longitude of edge crossing points on each grid
! elat         Latitude of edge crossing points on each grid
! farea        Area of primal faces on each grid
! varea        Area of dual faces on each grid
! ldist        Primal edge length, i.e. distance between neighbouring vertices
! ddist        Dual edge length, i.e. distance between neighbouring face centres
REAL*8, ALLOCATABLE :: flong(:,:), flat(:,:), &
                       vlong(:,:), vlat(:,:), &
                       elong(:,:), elat(:,:), &
                       farea(:,:), varea(:,:), &
                       ldist(:,:), ddist(:,:)

! Conventions
!
! Latitude and longitude in radians.
! Area and lengths are for the unit sphere.


! HODGE STAR, MASS MATRIX, AND RELATED OPERATORS

! nlsten        Number of faces in stencil for L mass matrix
! lsten         Stencil for L mass matrix
! lmass         Coefficients for L mass matrix
! nmsten        Number of faces in stencil for M mass matrix
! msten         Stencil for M mass matrix
! mmass         Coefficients for M mass matrix
! njsten        Number of vertices in stencil for J operator
! jsten         Stencil for J operator
! jstar         Coefficients for J operator
! nhsten        Number of edges in stencil for H operator
! hsten         Stencil for H operator
! hstar         Coefficients for H operator
! nrsten        Number of vertices in stencil for R operator (= neoff)
! rsten         Stencil for R operator (= voff)
! rcoeff        Coefficients for R operator
! nrxsten       Number of faces in stencil for R transpose operator (= neofv)
! rxsten        Stencil for R transpose operator (= fofv)
! rxcoeff       Coefficients for R transpose operator
! nwsten        Number of edges in stencil for W operator
! wsten         Stencil for W operator
! wcoeff        Coefficients for W operator
! ntsten        Number of edges in stencel for T operator
! tsten         Stencil for T operator
! tcoeff        Coefficients for T operator
INTEGER, ALLOCATABLE :: nlsten(:,:), nmsten(:,:), njsten(:,:), &
                        nhsten(:,:), nrsten(:,:), nrxsten(:,:), &
                        nwsten(:,:), ntsten(:,:)
INTEGER, ALLOCATABLE :: lsten(:,:,:), msten(:,:,:), jsten(:,:,:), &
                        hsten(:,:,:), rsten(:,:,:), rxsten(:,:,:), &
                        wsten(:,:,:), tsten(:,:,:)
REAL*8, ALLOCATABLE :: lmass(:,:,:), mmass(:,:,:), jstar(:,:,:), &
                       hstar(:,:,:), rcoeff(:,:,:), rxcoeff(:,:,:), &
                       wcoeff(:,:,:), tcoeff(:,:,:,:)
INTEGER :: nlsmx, nmsmx, njsmx, nhsmx, nrsmx, nrxsmx, nwsmx, ntsmx


! RESTRICTION AND PROLONGATION OPERATORS FOR MULTIGRID

! ninj           Number of faces in stencil for injection operator
! injsten        Stencil for injection operator
! injwgt         Weights for injection operator
INTEGER, ALLOCATABLE :: ninj(:,:), injsten(:,:,:)
REAL*8, ALLOCATABLE :: injwgt(:,:,:)
INTEGER :: ninjmx


! INFORMATION DEFINING COMPOUND ELEMENTS

! ncvp            Number of internal dofs to define a compound
!                 P0 element in space Vp.
! cvp             Dofs to define a compound element in space Vp.
! ncsp            Number of internal dofs to define a compound
!                 RT0 element in space Sp.
! csp             Dofs to define a compound element in space Sp.
! ncep            Number of internal dofs to define a compound
!                 P1 element in space Ep.
! cep             Dofs to define a compound element in space Ep.
! ncvd            Number of internal dofs to define a compound
!                 P0 element in space Vd.
! cvd             Dofs to define a compound element in space Vp.
! ncsd            Number of internal dofs to define a compound
!                 N0 element in space Sd.
! csd             Dofs to define a compound element in space Sd.
INTEGER, ALLOCATABLE :: ncvp(:), ncsp(:), ncep(:), &
                        ncvd(:), ncsd(:)
REAL*8, ALLOCATABLE :: cvp(:,:), csp(:,:), cep(:,:), &
                       cvd(:,:), csd(:,:)
INTEGER :: ncvpmx, ncspmx, ncepmx, ncvdmx, ncsdmx


END MODULE grid




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
CHARACTER*31 :: ygridfile = 'gridopermap_cube_0000221184.dat'

! File to output grid coordinates (for use in generating reference
! solutions).
! CHARACTER*30 :: ygridcoords = 'gridcoords_hex_0000010242.dat '
CHARACTER*30 :: ygridcoords = 'gridcoords_cube_0000221184.dat'

! Run identifier for naming output files
CHARACTER*6 :: runid = '000001'

! True for a restart run
LOGICAL :: lrestart = .FALSE.

! Name of restart file
CHARACTER*30 :: yresfile = 'run000002_restart_00000720.dat'


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
    ALLOCATE(fnxtf(nfacex,nefmx,ngrids), eoff(nfacex,nefmx,ngrids), &
             voff(nfacex,nefmx,ngrids),  fnxte(nedgex,2,ngrids),    &
             vofe(nedgex,2,ngrids),      fofv(nvertx,nevmx,ngrids), &
             eofv(nvertx,nevmx,ngrids))

    ! Read the connectivity arrays
    READ(changrid) (((fnxtf(if0,ix,igrid),          &
                         if0 = 1, nface(igrid)),    &
                         ix = 1, nefmx),            &
                         igrid = 1, ngrids)
    READ(changrid) (((eoff(if0,ix,igrid),           &
                         if0 = 1, nface(igrid)),    &
                         ix = 1, nefmx),            &
                         igrid = 1, ngrids)
    READ(changrid) (((voff(if0,ix,igrid),           &
                         if0 = 1, nface(igrid)),    &
                         ix = 1, nefmx),            &
                         igrid = 1, ngrids)
    READ(changrid) (((fnxte(ie0,ix,igrid),          &
                         ie0 = 1, nedge(igrid)),    &
                         ix = 1, 2),                &
                         igrid = 1, ngrids)
    READ(changrid) (((vofe(ie0,ix,igrid),           &
                         ie0 = 1, nedge(igrid)),    &
                         ix = 1, 2),                &
                         igrid = 1, ngrids)
    READ(changrid) (((fofv(iv0,ix,igrid),           &
                         iv0 = 1, nvert(igrid)),    &
                         ix = 1, nevmx),            &
                         igrid = 1, ngrids)
    READ(changrid) (((eofv(iv0,ix,igrid),           &
                         iv0 = 1, nvert(igrid)),    &
                         ix = 1, nevmx),            &
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

    CLOSE(changrid)
END SUBROUTINE readgrid





SUBROUTINE writeCellAdjacency

    USE runtype
    USE constants
    USE grid
    USE channels
    USE runtype

    IMPLICIT NONE

    ! identifyer string in the format
    ! [cube/hex]_[cells]
    CHARACTER*127 :: gridTypeAndCells

    CHARACTER*127 :: adjacencyFile
    CHARACTER*127 :: adjacencyFileGnuplot
    CHARACTER*127 :: adjacencyFileGnuplotOutput

    INTEGER :: igrid
    INTEGER :: nf, ne, nv
    INTEGER :: if0, if1, if2
    INTEGER :: ie2
    INTEGER :: iv2
    INTEGER :: ixn
    ! ----------------------------------------------------------------


    WRITE (gridTypeAndCells, '(A,''_'',I10.10)') TRIM(gridType), nface(ngrids)

    ! GNUPLOT: instruction file to generate png images
    adjacencyFileGnuplot = 'adjacencies_'//TRIM(gridTypeAndCells)//'.gnuplot'
    OPEN(chanAdjGnuplotOUT, FILE=adjacencyFileGnuplot)
    ! GNUPLOT: write to PNG files
    WRITE(chanAdjGnuplotOUT, *) "set terminal png size 1024,800"

    ! Iterate over all grid levels
    DO igrid = 1, ngrids

        nf = nface(igrid)
        ne = nedge(igrid)
        nv = nvert(igrid)

        WRITE (*,*)
        WRITE (*,*) "Cells on current level: ", nf

        !
        ! CELLS to CELLS
        !
            ! ADJ DATA: open file to write adjacency data to
            if1 = nface(ngrids)
            WRITE (adjacencyFile, '(''adjacencies_'',A,''_cell2cell_'',I10.10,''_level'',I2.2,''.dat'')') &
                   TRIM(gridType), if1, igrid
            WRITE (*,*) "Output file: "//TRIM(adjacencyFile)

            OPEN(chanaout, FILE=adjacencyFile)
            ! Iterate over all cells of the current grid level igrid
            DO if1 = 1, nf
                ! iterate over all adjacent cells, given by the number of adjacent edges
                DO ixn = 1, neoff(if1,igrid)
                    ! adjacent cell id
                    ! 1st index: cell to search adjacent cells for
                    ! 2nd index: adjacent cell
                    ! 3rd index: grid level
                    if2 = fnxtf(if1,ixn,igrid)

                    IF (if2 > nf) THEN
                        WRITE (*,*) "Invalid cell id"
                        WRITE (*,*) if2
                    END IF

                    WRITE(chanaout,*) if1, if2
                END DO
            END DO
            CLOSE(chanaout)


            ! GNUPLOT: plot instructions
            if1 = nface(ngrids)
            WRITE(adjacencyFileGnuplotOutput, '(''adjacencies_'',A,''_cell2cell_'',I10.10,''_level'',I2.2,''.png'')') &
                  TRIM(gridType), if1, igrid
            WRITE(chanAdjGnuplotOUT, '(A)') 'set output "'//TRIM(adjacencyFileGnuplotOutput)//'"'
            WRITE(chanAdjGnuplotOUT, '(A)') 'plot "'//TRIM(adjacencyFile)//'" with dots'
            WRITE(chanAdjGnuplotOUT, *)


        !
        ! CELLS to EDGES
        !
            ! ADJ DATA: open file to write adjacency data to
            if1 = nface(ngrids)
            WRITE (adjacencyFile, '(''adjacencies_'',A,''_cell2edge_'',I10.10,''_level'',I2.2,''.dat'')') &
                   TRIM(gridType), if1, igrid
            WRITE (*,*) "Output file: "//TRIM(adjacencyFile)

            OPEN(chanaout, FILE=adjacencyFile)
            ! Iterate over all cells of the current grid level igrid
            DO if1 = 1, nf
                ! iterate over all adjacent cells, given by the number of adjacent edges
                DO ixn = 1, neoff(if1,igrid)
                    ! adjacent edge id
                    ie2 = eoff(if1,ixn,igrid)

                    IF (ie2 > ne) THEN
                        WRITE (*,*) "Invalid edge id"
                        WRITE (*,*) ie2
                    END IF

                    WRITE(chanaout,*) if1, ie2
                END DO
            END DO
            CLOSE(chanaout)


            ! GNUPLOT: plot instructions
            if1 = nface(ngrids)
            WRITE(adjacencyFileGnuplotOutput, '(''adjacencies_'',A,''_cell2edge_'',I10.10,''_level'',I2.2,''.png'')') &
                  TRIM(gridType), if1, igrid
            WRITE(chanAdjGnuplotOUT, '(A)') 'set output "'//TRIM(adjacencyFileGnuplotOutput)//'"'
            WRITE(chanAdjGnuplotOUT, '(A)') 'plot "'//TRIM(adjacencyFile)//'" with dots'
            WRITE(chanAdjGnuplotOUT, *)



        !
        ! CELLS to VERTICES
        !
            ! ADJ DATA: open file to write adjacency data to
            if1 = nface(ngrids)
            WRITE(adjacencyFile, '(''adjacencies_'',A,''_cell2vertex_'',I10.10,''_level'',I2.2,''.dat'')') &
                  TRIM(gridType), if1, igrid
            WRITE(*,*) "Output file: "//TRIM(adjacencyFile)

            OPEN(chanaout, FILE=adjacencyFile)
            ! Iterate over all cells of the current grid level igrid
            DO if1 = 1, nf
                ! iterate over all adjacent cells, given by the number of adjacent edges
                DO ixn = 1, neoff(if1,igrid)
                    ! adjacent vertex id
                    iv2 = voff(if1,ixn,igrid)
                    IF (iv2 > nv) THEN
                        WRITE (*,*) "Invalid vertex id"
                        WRITE (*,*) iv2
                    END IF
                    WRITE(chanaout,*) if1, iv2
                END DO
            END DO
            CLOSE(chanaout)


            ! GNUPLOT: plot instructions
            if1 = nface(ngrids)
            WRITE(adjacencyFileGnuplotOutput, '(''adjacencies_'',A,''_cell2vertex_'',I10.10,''_level'',I2.2,''.png'')')&
                   TRIM(gridType), if1, igrid
            WRITE(chanAdjGnuplotOUT, '(A)') 'set output "'//TRIM(adjacencyFileGnuplotOutput)//'"'
            WRITE(chanAdjGnuplotOUT, '(A)') 'plot "'//TRIM(adjacencyFile)//'" with dots'
            WRITE(chanAdjGnuplotOUT, *)

    END DO

    CLOSE(chanAdjGnuplotOUT)

END SUBROUTINE


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

SUBROUTINE writeCellTraversal

    USE runtype
    USE constants
    USE grid
    USE channels
    USE runtype

    IMPLICIT NONE

    ! identifyer string in the format
    ! [cube/hex]_[cells]
    CHARACTER*127 :: gridTypeAndCells

    CHARACTER*127 :: adjacencyFile
    CHARACTER*127 :: adjacencyFileGnuplot
    CHARACTER*127 :: adjacencyFileGnuplotOutput

    INTEGER :: igrid
    INTEGER :: nf, ne, nv
    INTEGER :: if0, if1, if2
    INTEGER :: ie2
    INTEGER :: iv2
    INTEGER :: ixn
    REAL*8 :: x0(3)
    ! ----------------------------------------------------------------


    WRITE(gridTypeAndCells, '(A,''_'',I10.10)') TRIM(gridType), nface(ngrids)

    ! GNUPLOT: instruction file to generate png images
    adjacencyFileGnuplot = 'adjacencies_traversal_'//TRIM(gridTypeAndCells)//'.gnuplot'
    OPEN(chanAdjGnuplotOUT, FILE=adjacencyFileGnuplot)
    ! GNUPLOT: write to PNG files
    WRITE(chanAdjGnuplotOUT, '(A)') "set terminal png size 1024,800"

    ! Iterate over all grid levels
    DO igrid = 1, ngrids

        nf = nface(igrid)
        ne = nedge(igrid)
        nv = nvert(igrid)

        WRITE (*,*)
        WRITE (*,*) "Cells on current level: ", nf

        !
        ! CELLS to CELLS
        !
            ! ADJ DATA: open file to write adjacency data to
            if1 = nface(ngrids)
            WRITE (adjacencyFile, '(''adjacencies_traversal_'',A,''_cell_'',I10.10,''_level'',I2.2,''.dat'')')&
                   TRIM(gridType), if1, igrid
            WRITE (*,*) "Output file: "//TRIM(adjacencyFile)

            OPEN(chanaout, FILE=adjacencyFile)
            ! Iterate over all cells of the current grid level igrid
            DO if0 = 1, nf
                ! Coordinates of face if0
                CALL ll2xyz(flong(if0,igrid),flat(if0,igrid),x0)
                WRITE(chanaout,*) x0
            END DO
            CLOSE(chanaout)


            ! GNUPLOT: plot instructions

            if1 = nface(ngrids)
            WRITE(adjacencyFileGnuplotOutput, '(''adjacencies_traversal_'',A,''_cell_'',I10.10,''_level'',I2.2,''.png'')') &
                  TRIM(gridType), if1, igrid
            WRITE(chanAdjGnuplotOUT, '(A)') 'set output "'//TRIM(adjacencyFileGnuplotOutput)//'"'
            WRITE(chanAdjGnuplotOUT, '(A)') 'splot "'//TRIM(adjacencyFile)//'" u 1:2:3 with lines'
            WRITE(chanAdjGnuplotOUT, '(A)')


        !
        ! CELLS to EDGES
        !
            ! ADJ DATA: open file to write adjacency data to
!            WRITE (adjacencyFile, '(''adjacencies_'',A,''_cell2edge_'',I10.10,''_level'',I2.2,''.dat'')') TRIM(gridType), nface(ngrids), igrid
!            WRITE (*,*) "Output file: "//TRIM(adjacencyFile)
!
!            OPEN(chanaout, FILE=adjacencyFile)
!            ! Iterate over all cells of the current grid level igrid
!            DO if1 = 1, nf
!                ! iterate over all adjacent cells, given by the number of adjacent edges
!                DO ixn = 1, neoff(if1,igrid)
!                    ! adjacent edge id
!                    ie2 = eoff(if1,ixn,igrid)
!
!                    IF (ie2 > ne) THEN
!                        WRITE (*,*) "Invalid edge id"
!                        WRITE (*,*) ie2
!                    END IF
!
!                    WRITE(chanaout,*) if1, ie2
!                END DO
!            END DO
!            CLOSE(chanaout)


            ! GNUPLOT: plot instructions
!            WRITE (adjacencyFileGnuplotOutput, '(''adjacencies_'',A,''_cell2edge_'',I10.10,''_level'',I2.2,''.png'')') TRIM(gridType), nface(ngrids), igrid
!            WRITE(chanAdjGnuplotOUT, *) 'set output "'//TRIM(adjacencyFileGnuplotOutput)//'"'
!            WRITE(chanAdjGnuplotOUT, *) 'plot "'//TRIM(adjacencyFile)//'" with dots'
!            WRITE(chanAdjGnuplotOUT, *)



        !
        ! CELLS to VERTICES
        !
!            ! ADJ DATA: open file to write adjacency data to
!            WRITE (adjacencyFile, '(''adjacencies_'',A,''_cell2vertex_'',I10.10,''_level'',I2.2,''.dat'')') TRIM(gridType), nface(ngrids), igrid
!            WRITE (*,*) "Output file: "//TRIM(adjacencyFile)
!
!            OPEN(chanaout, FILE=adjacencyFile)
!            ! Iterate over all cells of the current grid level igrid
!            DO if1 = 1, nf
!                ! iterate over all adjacent cells, given by the number of adjacent edges
!                DO ixn = 1, neoff(if1,igrid)
!                    ! adjacent vertex id
!                    iv2 = voff(if1,ixn,igrid)
!                    IF (iv2 > nv) THEN
!                        WRITE (*,*) "Invalid vertex id"
!                        WRITE (*,*) iv2
!                    END IF
!                    WRITE(chanaout,*) if1, iv2
!                END DO
!            END DO
!            CLOSE(chanaout)


!            ! GNUPLOT: plot instructions
!            WRITE (adjacencyFileGnuplotOutput, '(''adjacencies_'',A,''_cell2vertex_'',I10.10,''_level'',I2.2,''.png'')') TRIM(gridType), nface(ngrids), igrid
!            WRITE(chanAdjGnuplotOUT, *) 'set output "'//TRIM(adjacencyFileGnuplotOutput)//'"'
!            WRITE(chanAdjGnuplotOUT, *) 'plot "'//TRIM(adjacencyFile)//'" with dots'
!            WRITE(chanAdjGnuplotOUT, *)

    END DO

    CLOSE(chanAdjGnuplotOUT)

END SUBROUTINE


PROGRAM adjacency_matrices

    USE runtype

    IMPLICIT NONE

    CHARACTER*127 :: argGridCells
    INTEGER :: gridCells

    IF (IARGC() < 2) THEN
        WRITE (*,*) "Usage: ./adjacency_matrices [cube/hex] [cells]"
        CALL EXIT(-1)
    END IF

    CALL GETARG(1, gridType)
    gridType = TRIM(gridType)
    CALL GETARG(2, argGridCells)
    READ (argGridCells, '(i10)') gridCells

    WRITE (ygridcoords, '(''gridmap_'',A,''_'',I10.10,''.dat'')') TRIM(gridType), gridCells
    WRITE (*,*) ygridcoords

    WRITE (ygridfile, '(''gridopermap_'',A,''_'',I10.10,''.dat'')') TRIM(gridType), gridCells
    WRITE (*,*) ygridfile

    WRITE (*,*) "Using grid file: "//ygridfile
    WRITE (*,*) "Using grid operator file: "//ygridcoords




    WRITE (*,*) "reading grid..."
    CALL readgrid


    WRITE (*,*) "writing cell adjacency..."
    CALL writeCellAdjacency

    WRITE (*,*) "writing cell traversal..."
    CALL writeCellTraversal

END PROGRAM
