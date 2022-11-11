module datastruct
  !====================================================================
  !
  ! Data structures for a geodesic grid with primal 
  !  (Delaney triangular) and dual (Voronoi) grids
  !
  ! Commnets: 
  !   - To obtain better memory alignement, the data is ordered
  !    in each type by size, from large numeric arrays to smaller
  !    and then strings. Source:
  !   software.intel.com/sites/products/documentation/hpc/compilerpro/
  !     en-us/fortran/lin/compiler_f/optaps/fortran/optaps_prg_algn_f.htm
  !  
  ! P. Peixoto - 10-03-2011
  !=====================================================================

  !Use global constants and kinds
  use constants, only: i4, r8

  !-------------------------------------------------
  ! Simple 3d cartesian vector
  !------------------------------------------------
  type vector
     real (r8), dimension(1:3) :: v
  end type vector

  !-------------------------------------------------
  ! Simple string of 8 chars, to be used when creating array of strings
  !------------------------------------------------
  type string8
     character (len=8) :: s
  end type string8

  !-------------------------------------------------
  ! Coeficients of a 2nd order polinomial
  ! c(1:5) refers to a x^2, xy, y^2, x, y
  ! Also has the rotation parameters
  !------------------------------------------------
  type polcoef
     !Polinomial coeficients
     real (r8), allocatable :: c(:)

     ! Rotation parameters
     real(r8):: cx, sx, cy, sy

     !Kind of coeficients stored
     ! "spol2" = 2nd order scalar polinomial
     ! "vpol1" = 1st order vector polonomial
     character (len=8) :: polkind

  end type polcoef

  !-------------------------------------------------
  !General sphere point structure
  ! This point is not necessarily a grid point 
  !-------------------------------------------------
  type point_structure

     !Vector form of cartesian coordinates
     real (r8), dimension(1:3) :: p

     !Spherical/geographic coordinates of node in radians
     !  lat in [-pi/2, pi/2] , lon in [-pi, pi[   
     real (r8) :: lat, lon

     !Barycentric coordinates relative to a triangle kct
     !  Could be used for other purposes
     real (r8), allocatable :: b(:)

     !Barycentric coordinates relative to a triangle kt
     !  Could be used for other purposes
     integer(i4) :: kt

  end type point_structure

  !----------------------------------------------------------------
  !Triangle vertice/node structure
  ! This has a bit more information then regular sphere point structure
  !------------------------------------------------------------------
  type vertice_structure

     !List of distances to neighbour nodes.
     ! The distances are actually the -cos(angle), where angle is the angle between
     ! the 2 nodes, that is, -dot_product(node, neighbour)
     ! this distance is a increasing function of the real distance, belongs to [-1,1].
     ! It is used in interpolation procedures. Size: 1:nnb
     real (r8), allocatable :: nbd(:)

     !List of geodesical distances to neighbour nodes (in radians).
     !  Size: 1:nnb
     real (r8), allocatable :: nbdg(:)

     !List of neighbour nodes/vertices (list of indexes). Size: 1:nnb
     ! counterclockwise-ordered sequences of neighboring nodes 
     integer (i4), allocatable :: nb(:)

     !List of edges that has this point as vertice.
     !Size: 1:nnb
     integer (i4), allocatable :: ed(:)

     !List of triangle that has this point as vertice.
     !Size: 1:nnb
     integer (i4), allocatable :: tr(:)

     !Vector form of cartesian coordinates
     real (r8), dimension(1:3) :: p

     !Spherical/geographic coordinates of node in radians
     !  lat in [-pi/2, pi/2] , lon in [-pi, pi[   
     real (r8) :: lon, lat

     !Number of neighbour nodes (conected through edges)
     integer (i4):: nnb

  end type vertice_structure

  !--------------------------------------------------------
  !Structure for edges
  !--------------------------------------------------------
  type edge_structure

     !Midpoint of the edge / intersecting point
     !  between triangle and hexagonal edges
     type(point_structure) :: c

     !Tangent unit vector to edge placed at the midpoint, 
     ! it must be tangent to the sphere, and is 
     ! defined in the deriction from v(1) to v(2)
     real (r8), dimension(1:3) :: tg

     !Normal unit vector to edge, placed at the midpoint
     ! it must be tangent to the sphere, and
     ! it is defined so that  tg x nr = c, 
     ! that is, given the tangent vector, the normal vector is the one
     ! obtained rotation the tangent vector counter-clockwise 90 degrees
     real (r8), dimension(1:3) :: nr

     !End point vertices (in case of triangle edges),
     !  or tringle circumcenters (in case of voronoi cell
     !Only the index of vertice/triangle is used
     integer (i4), dimension(1:2) :: v

     !Triangles/hexagons sharing this edge 
     ! In the case of triangles, only their index is stored
     ! In the case of hexagons, vertice indexes are stored
     integer (i4), dimension(1:2) :: sh

     !Projected length- Euclidian distance in R3
     ! considers the unit sphere
     real (r8) :: lenp

     !Geodesical length considering the unit sphere
     real (r8) :: leng

     !Area of edge volume by tiles -- not automatic loaded!!
     ! d_e - geodesic length of voronoi edge e
     ! l_e - geodesic length of triangle edge e
     ! A_t= d_e*l_e
     ! NOT saved with mesh
     real (r8) :: areat

  end type edge_structure

  !--------------------------------------------------------
  !Structure for triangles
  !--------------------------------------------------------
  type triangle_structure
     !Internal angles
     real (r8), dimension(1:3) :: angles

     !Circumcenter point 
     type(point_structure) :: c

     !Barycenter Point
     type(point_structure) :: b

     !Circumradius of planar triangle of unit sphere
     real (r8) :: radp

     !Circumradius of geodesical triangle
     real (r8) :: radg

     !Area of planar triangle of unit sphere
     real (r8) :: areap

     !Area of geodesical triangle of unit sphere
     real (r8) :: areag

     !Area by tiles
     ! d_e - geodesic distance from circuncenter to tr edge midpoint
     ! l_e - geodesic length of edge e
     ! A_t= sum_edges d_e*l_e/2
     ! NOT saved with mesh
     real (r8) :: areat

     !Generating vertices (index only)
     !  Counter clockwise ordering
     integer (i4), dimension(1:3) :: v

     !Edges generating triangle (index only)
     !  The ith edge is conecting v(i) and v(i+1)
     integer (i4), dimension(1:3) :: ed

     !Neighbour triangles (index only)
     integer (i4), dimension(1:3) :: nb

     ! Vector correction sign for tangent vectors: 
     !   For each edge, set 1 if the direction of the tangent
     !   vector, set counterclockwisely with respect to the triangle,
     !   is the same as the edge tangent vector, or else
     !   store -1     
     integer (i4), dimension(1:3) :: tg

     !Vector correction sign for normal vectors
     ! For each edge, set 1 if the direction of the outer normal
     !  vector is the same as the edge normal vector, or else,
     !  store -1. (Obs: nr(1:3)=-tg(1:3)
     integer (i4), dimension(1:3) :: nr

     ! Voronoi normal vector correction relative to
     ! triangle CCw direction
     ! For each edge of the triangle,
     ! set 1 if the direction of the normal vector of Hx edge
     ! is given counter clockwiselly with respect to the triangle,
     ! or else store -1
     !  =dsign(1., dot_product(mesh%edhx()%nr,mesh%ed()%tg))*mesh%tr()%tg()
     !integer(i4), dimension(1:3) :: hxnrccw

     !Dual-primal cell intersection area
     ! Area of triangle intersection with Voronoi/hexagonal cell
     ! Areas relative to cells given in %v(1:n)
     ! NOT saved with mesh
     real (r8), dimension(1:3) :: trhx_area

  end type triangle_structure

  !--------------------------------------------------------
  !Structure for voronoi cells (pentagons/hexagons)
  !  Each cell is centered at on grid vertice/node
  !  so that the vertice structure is shared with this one 
  !--------------------------------------------------------
  type hexagon_structure

     ! Tangent conterclockwise vector correction
     ! For each hx edge, set 1 if the direction of the tangent
     !  vector is given counter clockwiselly with respect to the polygon,
     ! or else store -1
     integer(i4), allocatable :: tg(:)

     ! Normal outputing vector correction
     ! For each hx edge, set 1 if the direction of the outer normal
     !  vector is the same as the edge normal vector, or else,
     !  store -1
     integer(i4), allocatable :: nr(:)

     ! Tringle edge tangent vector correction relative to Voronoi cells
     ! For each edge (triangle edge with respect to the voronoi cell),
     ! set 1 if the direction of the tangent points outward, or else,
     !  store -1 (tirangle tangent outwards
     !  =dot_product(mesh%edhx()%nr,mesh%ed()%tg)*mesh%hx()%nr
     integer(i4), allocatable :: ttgout(:)

     ! Tringle edge normal vector correction relative to Voronoi cells
     ! For each edge (triangle edge with respect to the voronoi cell),
     ! set 1 if the direction of the normal vector of tr edge
     ! is given counter clockwiselly with respect to the polygon,
     ! or else store -1
     !  =dot_product(mesh%edhx()%tg,mesh%ed()%nr)*mesh%hx()%tg
     integer(i4), allocatable :: tnrccw(:)

     !Barycenter Point
     type(point_structure) :: b     

     !Area of planar polygon
     real (r8) :: areap

     !Area of geodesical polygon
     real (r8) :: areag

     !Area by tiles
     ! d_e - geodesic length of voronoi edge e
     ! l_e - geodesic length of triangle edge e
     ! A_t= sum_edges d_e*l_e/4
     ! NOT saved with mesh
     real (r8) :: areat

     !Alignment index
     real (r8) :: align

     !Weights of ith edge of hexagon relative to jth
     !TRISK - (see Thuburn et all 2009)
     ! NOT saved with mesh
     real (r8), allocatable :: trskw(:, :)

     !Dual-primal cell intersection area
     ! Area of triangle intersection with Voronoi/hexagonal cell
     ! Areas relative to cells given in %v(1:n)
     ! Divided by cell areas (normalized)!!!
     ! NOT saved with mesh
     real (r8), allocatable :: hxtr_area(:)

      !Weights for each vertex of the cell such
      ! that div remap is consistent
      ! and they sum 1 in triangular cell
      ! NOT saved with mesh
     real (r8), allocatable :: hxtr_remapw(:)

     !Sum of weights for each vertex of the cell such
     ! that div remap is consistent
     ! and they sum 1 in triangular cell
     ! NOT saved with mesh
     real (r8) :: hxtr_remapwsum

  end type hexagon_structure

  !--------------------------------------------------------
  !Structure for regular grid rectangles/quadrilaterals
  !--------------------------------------------------------
  type quadrilateral_structure

     !List of triangles intersecting this rectangle
     ! Size: 1:ntr
     integer (i4), allocatable :: tr(:)

     !Regular grid longitude in radians
     ! lon in [-pi, pi[
     real (r8) :: lon

     !Number of intersecting triangles
     integer (i4) :: ntr

  end type quadrilateral_structure

  !--------------------------------------------------------
  !List of by quads for fixed latitudes (varying longitudes)
  !--------------------------------------------------------
  type quad_lon_structure

     !List of quadrilaterals
     type(quadrilateral_structure), allocatable :: qd_lon(:)

     !Regular grid latitude in radians
     !  lat in [-pi/2, pi/2]
     real (r8) :: lat

     !Number of longitudes
     integer (i4) :: nlon

     !Regular grid spacing between points in radians
     ! This is obtained as: dlon=2pi/nlon
     real (r8) :: dlon

  end type quad_lon_structure

  !---------------------------------
  ! Global grid structure
  !---------------------------------
  type grid_structure

     !List of vertices. Size: 1:nv
     type(vertice_structure), allocatable :: v(:)

     !List of triangle edges. Size: 1:ne
     type(edge_structure), allocatable :: ed(:)

     !List of triangles. Size: 1:nt
     type(triangle_structure), allocatable :: tr(:)

     !List of voronoi cell edges, they have the same index of the 
     !  triangle edges intersecting it (1 to 1 relation)
     ! Size: 1:ne
     type(edge_structure), allocatable :: edhx(:)

     !List of pentagonal/hexagonal voronoi polygons. Size: 1:nv
     !  This list should be used toghether with vertices list.
     type(hexagon_structure), allocatable :: hx(:)

     !List of quadrilaterals. Size: (1:nlon, 1:nlat)
     !type(quadrilateral_structure), allocatable :: qd(:,:)
     type(quad_lon_structure), allocatable :: qd_lat(:)

     !Regular grid spacing between points in radians
     ! This is obtained as: dlat=pi/nlat ; dlon=2pi/nlon
     real (r8) :: dlat !dlon,

     !Number latitudes and longitudes in regular grid
     !  these variables are used for tabular searching scheme
     integer(i4):: nlat !, nlon

     !Minimum/Maximum angular distance between vertice points in radians
     real(r8):: minvdist, maxvdist, meanvdist

     !Minimum/Maximum angular distance between circumcenters in radians
     real(r8):: mincdist, maxcdist, meancdist

     !Minimum/Maximum triangle geodesical areas 
     real(r8):: mintrarea, maxtrarea, meantrarea

     !Minimum/Maximum triangle geodesical internal angles
     real(r8) :: mintrangle, maxtrangle, meantrangle

     !Minimum/Maximum voronoi cell geodesical areas 
     real(r8) :: minhxarea, maxhxarea, meanhxarea

     !Optimization parameter (dependes on the method)
     real(r8) :: optpar

     !Maximum number of neighbours for all nodes (triangle vertices)
     integer(i4):: maxvnb

     !Maximum number of triangles intersecting a regular qrid rectangle
     integer(i4):: maxtrsqint

     !Number of bisection (g-level)
     integer (i4) :: glevel

     !Flag for loadable grid (1- yes, 0- no)
     integer (i4) :: loadable

     !Number of triangle vertices, that is, 
     ! number of nodes on generated grid
     integer(i4)::nv

     !Number of triangles
     integer (i4)::nt 

     !Number of triangle edges
     integer (i4):: ne

     !Mesh kind
     !icos (icosahedral)
     !octg (octahedral)
     !read (read from file)
     !rand (random points - be careful...)
     character(len=16) :: kind

     !Mesh Positions:
     !  eqs (equator symmetric)
     !  pol (north pole point)
     !  ran (all random points)
     !  ref (random with local mesh refinement - needs scvt optimization)
     character(len=16) :: pos

     ! 'nopt'  -> No optimization
     ! 'scvt'   -> Spherical centroidal voronoi
     !            Nodes are the mass centers of the voronoi cells
     ! 'sprg'  -> Spring dynamics optimization
     character(len=16) :: optm

     !Hierarchical construction/optimization flag
     !  0-No;
     !  1-Yes;
     !  2 - Optimizes only the new points in hierarchy;
     !  3 - Do not optimize nodes that belong to primary icosahedral edges
     integer (i4) :: hrchy

     !Triangulation read from file
     !  0-No;
     !  1-Yes;
     ! mesh%kind must be 'read'
     ! and file to be read must have the same name as
     ! the files with nodes but extention .nbd (neighbours)
     !  If =1, it is not necessary to do the delaunay triangulation
     integer (i4) :: deltriread

     !Grid name, used for file names as outputs
     ! name=kind//pos//optm//hrchy//n
     ! or file name if nodes are read from file
     character (len=128) :: name
     character (len=256) :: filename

     real (r8), allocatable  :: densf_table(:,:)

  end type grid_structure

  !---------------------------------------------------------
  ! Variable for scalar values on geodesic grid
  !---------------------------------------------------------
  type scalar_field

     ! Values array, ordered in the same sequence as the
     !   original nodes/circumcenter/edges mid points
     !   indexation
     real (r8), allocatable  :: f(:)

     ! Gradient vectors on the respective nodes
     ! The gradient calculation routines are only done for
     !  position = 0
     type(vector), allocatable  :: g(:)

     ! 2nd Order polinomial least saure fit on node
     ! or 1st order vector field recosntruction least square fit
     type(polcoef), allocatable :: pol(:)

     !Number of values on the mesh
     integer (i4) :: n

     ! Position of the values relative to a mesh
     !   0 - Nodes - Triangle vertices (voronoi centers)
     !   1 - Triangle circumcenter (voronoi vertices)
     !   2 - Triangle's edge midpoint
     !   3 - Hexagon's edges midpoint
     !   4 - Voronoi cell barycenter
     !   5 - Triangle barycenter
     !   6 - Edges - Intersection point of tr x hx, reference hx
     ! Obs: 6 is essencially the same as 2, but when
     !   vector components are put at 6 the refence are the hx
     !   whereas in 2, the reference are the triangles
     integer (i4) :: pos

     !Variable name - long name - detailed name
     ! This is used to for filenames of this variable
     character (len=256) :: name

  end type scalar_field

  !---------------------------------------------------------
  ! Variables for mesh vector values
  !---------------------------------------------------------

  type vector_field_cart
     !Cartesian vector variable

     ! Vector in cartesian coord on 'pos'
     type(vector), allocatable  :: p(:)

     ! Position of the values relative to a mesh
     !   0 - Nodes - Triangle's vertices (voronoi centers)
     !   1 - Triangle's circumcenters (voronoi vertices)
     !   2 - Triangle's edge mid points
     !   3 - Hexagonal's edge mid points
     !   4 - Voronoi cell barycenter
     !   5 - Triangle barycenter
     !   6 - Edges - Intersection point of tr x hx
     integer (i4) :: pos

     !Number of values on the mesh
     integer (i4) :: n

     !Variable name - long name - detailed name
     ! This is used to for filenames of this variable
     character (len=256) :: name

  end type vector_field_cart

  type vector_variable_edge
     ! Vector defined on edges

     ! Normal vector component value
     ! Ordered in the same sequence as the
     !   original edges indexation
     real (r8), allocatable  :: nr(:)

     ! Tangent vector component value
     ! Ordered in the same sequence as the
     !   original edges indexation
     real (r8), allocatable  :: tg(:)

     ! Position of the values relative to a mesh
     !   0 - Triangle's edge mid points
     !   1 - Voronoi cell edge mid point
     integer (i4) :: pos

     !Number of edges on the mesh
     integer (i4) :: n

     !Variable name - long name - detailed name
     ! This is used to for filenames of this variable
     character (len=256) :: name

  end type vector_variable_edge

  type vector_variable_uv
     ! West-East vector value
     ! Ordered in the same sequence as the
     !   original nodes/circumcenter/edges mid points
     !   indexation
     real (r8), allocatable  :: u(:)

     ! South-North vector value
     ! Ordered in the same sequence as the
     !   original nodes/circumcenter/edges mid points
     !   indexation
     real (r8), allocatable  :: v(:)

     ! Position of the values relative to a mesh
     !   0 - Nodes - Triangle vertices (voronoi centers)
     !   1 - Triangle circumcenter (voronoi vertices)
     !   2 - Triangle's edge mid points
     integer (i4) :: pos

     !Number of values on the mesh
     integer (i4) :: n

     !Variable name - long name - detailed name
     ! This is used to for filenames of this variable
     character (len=256) :: name

  end type vector_variable_uv

  !---------------------------------------------------------
  ! Structures for RBF matrices
  !---------------------------------------------------------

  type rbf_matrix_structure
     !Radial Basis Function matrix for a given point

     !Cholesky decomposition (LDL') of rbf matrix
     !It has a lower triangular matrix L and 
     !  a diagonal D stored in the diagonal of L
     real(r8), allocatable :: L(:,:)

     !Right hand side of RBF system
     real(r8), allocatable :: b(:)

     ! RBF weights for interpolation
     real(r8), allocatable :: w(:)

     ! Polinomial to be added
     ! For scalar rbf, dim(pol)=1
     ! For vec recon, dim(pol)=3
     real(r8), allocatable :: pol(:)

     ! Estimative of the condition number
     real(r8) :: condnum, dmin, dmax

     !Nodes/triangles indexes involved in the rbf matrix
     integer (i4) , allocatable :: pts(:)

     !RBF width (shape parameter)
     real (r8) :: h

     !Number of points involved in the rbf matrix
     integer (i4) :: n

  end type rbf_matrix_structure

  !---------------------------------------------------------
  ! Structures for Natural Neighbour Coordinates
  !---------------------------------------------------------

  type general_coords
     !Natural coordinate vectors

     !Number of natural neighbours or near neighbours
     integer (i4) :: n

     ! Position of the values relative to a mesh
     !   0 - Nodes - Triangle's vertices (voronoi centers)
     !   1 - Triangle's circumcenters (voronoi vertices)
     !   2 - Triangle's edge mid points
     !   3 - Hexagonal's edge mid points
     !   4 - Voronoi cell barycenter
     !   5 - Triangle barycenter
     !   6 - Edges - Intersection point of tr x hx
     integer (i4) :: pos

     !Neighbour indexes (vertice index, or edge index)
     integer (i4), allocatable :: v(:)

     !Coordinate value, or weight relative to each vertice
     real (r8), allocatable :: w(:)

  end type general_coords

  !----------------------------------------------------------
  ! Structure for multi method vector reconstruction
  !----------------------------------------------------------
  type vectorinterpol_methods
     !Methods involved in vector reconstruction

     !Basic reconstruction method
     character (len=16):: recon

     !Vector interpolation method
     character (len=16):: interp

     !Higher order reconstruction
     character (len=16):: horecon

     !Full name
     character (len=24):: name

     !Alignment index cut off value
     real (r8) :: alignlimit

     !Percentage of aligned cells
     real (r8) :: alignpercent

     !Logic for hybrid method (horeconq/=0)
     logical :: hyb

     !Logic for reconstruction for barycenters
     logical :: massc

     !Logical variable to reconstruc only on ill aligned cells
     ! Do nothing on aligned cells
     logical :: nonalignedonly

     !Radial basis function parameter
     real (r8) :: rbf_par

  end type vectorinterpol_methods


end module datastruct
