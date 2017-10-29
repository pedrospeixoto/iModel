MODULE grid

! Based on grid of femswe, with additional variables of other files

IMPLICIT NONE

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
! farea        Area of faces on each grid
! varea        Area of dual faces on each grid
! ldist        Primal edge length, i.e. distance between neighbouring face centres
! ddist        Dual edge length, i.e. distance between neighbouring vertices
! fareamin     Minimum face area on each grid
REAL*8, ALLOCATABLE :: flong(:,:), flat(:,:), &
                       vlong(:,:), vlat(:,:), &
                       farea(:,:), varea(:,:), &
                       ldist(:,:), ddist(:,:), &
                       fareamin(:)

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
! jlump         Coefficients of lumped J matrix
! mlump         Coefficients of lumped M matrix
! hlump         Coefficients of lumped H matrix
! nxminvten     Number of edges in stencil for approximate inverse of M
! xminvsten     Stencil for approximate inverse of M
! xminv         Coefficients for approximate inverse of M
INTEGER, ALLOCATABLE :: nlsten(:,:), nmsten(:,:), njsten(:,:), &
                        nhsten(:,:), nrsten(:,:), nrxsten(:,:), &
                        nwsten(:,:), ntsten(:,:), nxminvsten(:,:)
INTEGER, ALLOCATABLE :: lsten(:,:,:), msten(:,:,:), jsten(:,:,:), &
                        hsten(:,:,:), rsten(:,:,:), rxsten(:,:,:), &
                        wsten(:,:,:), tsten(:,:,:), xminvsten(:,:,:)
REAL*8, ALLOCATABLE :: lmass(:,:,:), mmass(:,:,:), jstar(:,:,:), &
                       hstar(:,:,:), rcoeff(:,:,:), rxcoeff(:,:,:), &
                       wcoeff(:,:,:), tcoeff(:,:,:,:), jlump(:,:), &
                       mlump(:,:), hlump(:,:), xminv(:,:,:)
INTEGER :: nlsmx, nmsmx, njsmx, nhsmx, nrsmx, nrxsmx, nwsmx, ntsmx, nxmisx


! RESTRICTION AND PROLONGATION OPERATORS FOR MULTIGRID

! nres           Number of faces in stencil for restriction operator
! ressten        Stencil for restriction operator
! reswgt         Weights for restriction operator
INTEGER, ALLOCATABLE :: nres(:,:), ressten(:,:,:)
REAL*8, ALLOCATABLE :: reswgt(:,:,:)
INTEGER :: nresmx

INTEGER, ALLOCATABLE :: ninj(:,:), injsten(:,:,:)
REAL*8, ALLOCATABLE :: injwgt(:,:,:)
INTEGER :: ninjmx


! 0: no SFC optimization applied
! 1: SFC optimization applied and SFC Index information appended
INTEGER :: SFCIndexAvailable

INTEGER, ALLOCATABLE :: fNewFaceId(:,:), fNewFaceIdInverse(:,:), &
			fNewVertId(:,:), fNewVertIdInverse(:,:), &
			fNewEdgeId(:,:), fNewEdgeIdInverse(:,:) 

! Additional variables from grid_hex

INTEGER, ALLOCATABLE :: feofe(:,:,:)

REAL*8, ALLOCATABLE :: gdist(:,:), coeff(:,:,:)

! Additional variables from grid_cube

! n x n cells on each panel
! Smallest and largest n
INTEGER :: n0
INTEGER :: nl
INTEGER :: nx2

! Number of smoothing iterations. Must be at least 1 for
! consistency of H operator.
INTEGER :: nsmooth = 1

! Additional variables from grid_buildop

! Grid type (just used to determine file names)
! igtype       1 for hex, 2 for cube
INTEGER :: igtype

REAL*8, ALLOCATABLE :: elong(:,:), elat(:,:)

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


! Variables for lat-lon grid table - P. Peixoto
INTEGER*4, ALLOCATABLE :: lltable_primal(:,:), lltable_dual(:,:)
INTEGER :: lltable_nlat, lltable_nlon

END MODULE grid
