module interpack
  !========================================================================
  !
  !  INTERPACK
  !
  !	Interpolation/Reconstruction pack for the geodesic spheric grids
  !
  !	Pedro da Silva Peixoto (pedrosp@ime.usp.br)
  !	March 2013
  !
  !========================================================================

  !Global constants
  use constants, only: &
       datadir, &
       deg2rad, &
       eps, &
       i4, &
       r8, &
       r4, &
       rad2deg

  !Mesh and variables data structures
  use datastruct, only: &
       general_coords, &
       grid_structure, &
       scalar_field, &
       rbf_matrix_structure, &
       vectorinterpol_methods, &
       vector, &
       vector_field_cart, &
       vector_variable_edge, &
       vector_variable_uv

  !Spherical unstructured mesh routines
  use smeshpack, only: &
       arclen, &
       bar_coord, &
       bar_coord_tr, &
       choleskydecomp, &
       choleskysolve, &
       condnumest, &
       convert_vec_cart2sph, &
       cross_product, datadir, &
       gcircarcintersec, &
       getnearhxedges, &
       getneartredge, &
       getnearhxedge, &
       getnearnode, &
       getnearnodes, &
       getneartrc, &
       getothernodeonedge, &
       gettr, &
       getunit, &
       hxtr_intersec_areas, &
       insmallvoredtr, &
       modint, &
       norm, &
       ortogonalarc, &
       planarpjhexagarea, &
       positive, &
       proj_vec_sphere, &
       solve3x3, &
       solvelintri, &
       sph2cart, &
       sphpolarea, &
       sphtriarea

  implicit none

contains

  !========================================================================
  !    GENERAL INTERPOLATION PROCEDURES
  !========================================================================

  function scalar_interpol (p, var, mesh, kinterp, rbf_mat, monot)
    !-----------------------------------------------------------
    !  SCALAR INTERPOLATION
    !
    !   Given the triangulation of a set of nodes on the unit
    ! sphere (mesh), along with data values at the nodes/vertices/triangle centers
    ! (in var), computes the value at a point p
    ! using 'kinterp' interpolation.
    !-------------------------------------------------------
    !Interpolation point
    real (r8), intent(in) :: p(1:3)

    !Variable for interpolation procedure
    type(scalar_field), intent(inout) :: var

    !Mesh in which the variable belongs
    type(grid_structure), intent(in) :: mesh

    !++++++++++ Kind of interpolations +++++++++++++++++++
    !*** var%pos=0 (For values given on triangle vertices) ***
    ! none     = Just get triangle
    ! neartrv  = Nearest node  (default)
    ! lintrv   = Linear (barycentric coords) - default
    ! wach     = Wachspress coorinates interpol (linear)
    ! lsqhx    = Quadratic Least Squares
    ! hermtrv  = Cubic (R.Renka C1 Hermite)
    ! rbftr    = RadialBasisFunction TR (3pts)
    ! rbfetr   = RadialBasisFunction ETR(11,12pts)
    ! rbfhx    = RadialBasisFunction HX (6,7pts)
    ! rbf***pc = RadialBasisFunction with constant polynomial
    ! natlap   = Natural Coordinate with Laplace
    ! natsib   = Natural Coordinate with Sibson
    ! natfar   = Natural Coordinate with Farin
    ! lmshep   = Local Modified Shepard Method
    ! qdshep   = Quadratic Shepard Method
    !*** var%pos=1 (For values given on triangle circumcenters) ***
    ! none     = Just get triangle
    ! neartrc  = Nearest triangle
    ! wach     = Wachspress interpolation  (default)
    !*** var%pos=2 (For values given on triangle edge midpoints) ***
    ! none     = Just get triangle
    ! neartred = Nearest triangle edge
    ! p1nc     = P1 Non corforming element  (default)
    !*** var%pos=3 (For values given on voronoi edge midpoints) ***
    ! none     = Just get triangle
    ! nearhxed = Nearest Voronoi cell edge
    ! wach     = Modified Wachspress interpolation (linear) - default
    ! lmshep   = Local Modified Shepard Method
    ! p1nc     = P1 Non corforming element - may do extrapolation
    !*** var%pos=4 (For values given on Voronoi centroids) ***
    ! neartrv  = Nearest node  (default)
    !*** var%pos=5 (For values given on triangle centroids) ***
    ! neartrc  = Nearest triangle  (default)
    character (len=8), optional :: kinterp

    !Aux kind of interpolation
    character (len=8):: kindinterpol

    !RBF matrix vector (used only on rbf interpolations)
    type(rbf_matrix_structure), optional :: rbf_mat(:)

    !Monotonicity flag
    integer(i4), optional :: monot
    integer(i4) :: monottmp

    !Interpolated value
    real (r8):: scalar_interpol

    !RBF/NAT Stencil
    character (len=4):: stencil

    !Radius of spherical cap (shepard method)
    real (r8):: r

    !Indexes
    integer (i4):: i
    integer (i4):: k

    !Set defaulf kind of interpolation (in case kindinterp not given)
    ! The default is always a linear interpolation
    if(.not.present(kinterp))then
       select case(var%pos)
       case(0)  !Values on triangle vertices / voronoi centers
          kindinterpol="lintrv"
       case(1)  !Values on triangles circumcenters / voronoi vertices
          kindinterpol="wach"
       case(2, 6) !Values on triangle edge midpoints
          kindinterpol="p1nc"
       case(3) !Values on hexagon/voronoi edge midpoints
          kindinterpol="wach"
       case(4) !Values on hexagon/voronoi centroids
          kindinterpol="neartrv"
       case(5) !Values given on triangle centroids
          kindinterpol="neartrc"
       case default
          print*, "Error on scalar_interpol:"
          print*, "    Unknown variable position", var%pos
          stop
       end select
    else
       !Local kind of interpolation
       kindinterpol=kinterp
    end if

    !Set RBF stencil kind
    if(kindinterpol(1:3)=='rbf')then
       if(.not.present(rbf_mat))then
          print*, "Scalar_interpol error: rbf_mat must be given!"
          stop
       end if
       stencil(1:4)=kindinterpol(4:7)
       kindinterpol='rbf'
    end if

    !Set NAT stencil kind - Natural Neighbour interpolation
    if(kindinterpol(1:3)=='nat')then
       stencil=kindinterpol(4:6)
       kindinterpol='nat'
    end if

    if(present(monot))then
       monottmp=monot
    else
       monottmp=0
    end if

    scalar_interpol=0.

    !Position selection
    select case(var%pos)
    case(0)  !Values on triangle vertices / voronoi nodes
       !Type of interpolation selection
       select case(trim(kindinterpol))
       case('none', 'testgrad') !Don't do any interpolation, just get triangle
          k=gettr(p, mesh)
          scalar_interpol=0._r8
       case('neartrv') !Nearest node
          !Triangle vertice value
          i=getnearnode(p, mesh)
          scalar_interpol=var%f(i)
       case('lintrv') !Linear - Barycentric coords
          scalar_interpol=scinterpol_linear_trv(p, var, mesh)
       case('wach') !Wachspress coordinates
          scalar_interpol=scinterpol_wachspress(p, var, mesh)
       case('lsqhx') !Quadratic least squares approximation
          scalar_interpol=scinterpol_lsqhx_trv(p, var, mesh)
       case('hermtrv') !Hermite - Renka Alg
          scalar_interpol=scinterpol_hermite_trv(p, var, mesh, monottmp)
       case('rbf') !RBF
          scalar_interpol=scinterpol_rbf(p, var, mesh, stencil, rbf_mat )
       case('nat') !Natural Neighbour Interpolation
          scalar_interpol=scinterpol_natneib(p, var, mesh, stencil)
       case('lmshep', 'LMSHEP') !Local Modified Shepard Method
          !Set radius of spherical cap
          r=mesh%maxvdist*1.05_r8
          scalar_interpol=scinterpol_lmshep(p, var, mesh, r)
       case('qdshep', 'QDSHEP') !Quadratic Shepard Method
          scalar_interpol=scinterpol_qdshep(p, var, mesh)
       case default
          print*, "Error on scalar_interpol:"
          print*, "    Unknown kind of interpolation for this variable"
          print*, "    kindinterpol=", kindinterpol, "  var%pos: ", var%pos
          stop
       end select

    case(1)  !Values on triangles circumcenters / voronoi vertices
       !Type of interpolation selector
       select case(trim(kindinterpol))
       case('none') !Don't do any interpolation, just get triangle
          k=gettr(p, mesh)
          scalar_interpol=0._r8
       case('neartrc') !Nearest triangle value
          !Triangle center value
          k=gettr(p, mesh)
          scalar_interpol=var%f(k)
       case('wach') !Wachspress coordinates
          scalar_interpol=scinterpol_wachspress(p, var, mesh)
       case default
          print*, "Error on scalar_interpol:"
          print*, "    Unknown kind of interpolation for this variable"
          print*, "    kindinterpol=", kindinterpol, "  var%pos: ", var%pos
          stop
       end select

    case(2, 6) !Values on triangle edge midpoints
       !Type of interpolation selector
       select case(trim(kindinterpol))
       case('none') !Don't do any interpolation, just get triangle
          k=gettr(p, mesh)
          scalar_interpol=0._r8
       case('neartred') !Nearest triangle edge
          k=getneartredge(p, mesh)
          scalar_interpol=var%f(k)
       case('p1nc')
          scalar_interpol=scinterpol_p1nc(p, var, mesh)
       case default
          print*, "Error on scalar_interpol:"
          print*, "    Unknown kind of interpolation for this variable"
          print*, "    kindinterpol=", kindinterpol, "  var%pos: ", var%pos
          stop
       end select

    case(3) !Values on hexagon/voronoi edge midpoints
       !Type of interpolation selector
       select case(trim(kindinterpol))
       case('none') !Don't do any interpolation, just get triangle
          k=gettr(p, mesh)
          scalar_interpol=0._r8
       case('nearhxed') !Nearest triangle edge
          k=getnearhxedge(p, mesh)
          scalar_interpol=var%f(k)
       case('p1nc')
          scalar_interpol=scinterpol_p1nc(p, var, mesh)
       case('wach')
          scalar_interpol=scinterpol_wachspress(p, var, mesh)
       case('lmshep')
          r=mesh%maxcdist*1.05_r8
          scalar_interpol=scinterpol_lmshep_edhx(p, var, mesh, r)
       case default
          print*, "Error on scalar_interpol:"
          print*, "    Unknown kind of interpolation for this variable"
          print*, "    kindinterpol=", kindinterpol, "  var%pos: ", var%pos
          stop
       end select

    case(4)  !Values on Voronoi centroids
       !Type of interpolation selection
       select case(trim(kindinterpol))
       case('none') !Don't do any interpolation, just get triangle
          k=gettr(p, mesh)
          scalar_interpol=0._r8
       case('neartrv') !Nearest node
          !Triangle vertice value
          i=getnearnode(p, mesh)
          scalar_interpol=var%f(i)
       case default
          print*, "Error on scalar_interpol:"
          print*, "    Unknown kind of interpolation for this variable"
          print*, "    kindinterpol=", kindinterpol, "  var%pos: ", var%pos
          stop
       end select

    case(5)  !Values on triangle centroids
       !Type of interpolation selection
       select case(trim(kindinterpol))
       case('none') !Don't do any interpolation, just get triangle
          k=gettr(p, mesh)
          scalar_interpol=0._r8
       case('neartrc') !Nearest triangle
          k=gettr(p, mesh)
          scalar_interpol=var%f(k)
       case default
          print*, "Error on scalar_interpol:"
          print*, "    Unknown kind of interpolation for this variable"
          print*, "    kindinterpol=", kindinterpol, "  var%pos: ", var%pos
          stop
       end select

    case default
       print*, "Error on scalar_interpol: Unknown variable position"
       print*, "    kindinterpol=", kindinterpol, "  var%pos: ", var%pos
       stop
    end select

    return
  end function scalar_interpol


  function vector_interpol (p, var, mesh, kinterp)
    !-----------------------------------------------------------
    !  VECTOR INTERPOLATION
    !
    !   Given the triangulation of a set of nodes on the unit
    ! sphere (mesh), along with vectors (cartesian coords) at
    ! the nodes/vertices/triangle centers (in var),
    ! computes the value at a point (p)
    ! using 'kinterp' interpolation.
    !-------------------------------------------------------
    !Interpolation point
    real (r8), intent(in) :: p(1:3)

    !Input Variable with vector field
    type(vector_field_cart), intent(in) :: var

    !Mesh in which the variable belongs
    type(grid_structure), intent(in) :: mesh

    !Kind of interpolation
    ! none    = Just get triangle
    ! lintrv  = Linear (barycentric coords) - var%pos=trv (0)
    ! wach   =  Wachspress interpolation - var%pos=trcc (1) or hxed (3)
    ! p1nc   =  P1 non conforming element - var%pos=tred (2)
    character (len=8), optional :: kinterp

    !Aux kind of interpolation
    character (len=8):: kindinterpol

    !Interpolated value
    real (r8):: vector_interpol(1:3)

    !Index vars
    integer(i4):: i, k

    !Set defaulf kind of interpolation (in case kinterp not given)
    if(.not.present(kinterp))then
       select case(var%pos)
       case(0)  !Values on triangle vertices / voronoi centers
          kindinterpol="lintrv"
       case(1)  !Values on triangles circumcenters / voronoi vertices
          kindinterpol="wach"
       case(2) !Values on triangle edge midpoints
          kindinterpol="p1nc"
       case(3) !Values on hexagon/voronoi edge midpoints
          kindinterpol="wach"
       case(4)  !Values on voronoi centroids
          kindinterpol="linhxb"
       case(5)  !Values on tr centroids
          kindinterpol="wach"
       case default
          print*, "Error on vector_interpol:"
          print*, "    Unknown variable position", var%pos
          stop
       end select
    else
       !Local kind of interpolation
       kindinterpol=kinterp
    end if

    !print*, var%pos, trim(kindinterpol)

    select case(var%pos)
    case(0) !Values on triangle vertices (voronoi nodes)
       !Type of interpolation selector
       select case(trim(kindinterpol))
       case('none') !Don't do any interpolation, just get triangle
          k=gettr(p, mesh)
          vector_interpol=0._r8
       case('lintrv') !Linear interpolation (barycentric coords)
          vector_interpol=vecinterpol_linear_trv(p, var, mesh)
       case default
          print*, "VECTOR INTERPOL ERROR: Don't know how to do this interpolation:", &
               kindinterpol, " with this positioning: ", var%pos
          stop
       end select
    case(1) !Values on triangle circumcenters (hx vertices)
       !Type of interpolation selector
       select case(trim(kindinterpol))
       case('none') !Don't do any interpolation, just get triangle
          i=getnearnode(p, mesh)
          vector_interpol=0._r8
       case('wach') !Linear interpolation (wachspress coords)
          vector_interpol=vecinterpol_wachspress(p, var,  mesh)
       case default
          print*, "VECTOR INTERPOL ERROR: Don't know how to do this interpolation:", &
               kindinterpol, " with this positioning: ", var%pos
          stop
       end select

    case(2) !Values on triangle edge midpoints
       !Type of interpolation selector
       select case(trim(kindinterpol))
       case('none') !Don't do any interpolation, just get triangle
          k=gettr(p, mesh)
          vector_interpol=0._r8
       case('p1nc') !Linear interpolation (P1 non conforming element)
          vector_interpol=vecinterpol_p1nc(p, var,  mesh)
       case default
          print*, "VECTOR INTERPOL ERROR: Don't know how to do this interpolation:", &
               kindinterpol, " with this positioning: ", var%pos
          stop
       end select

    case(3) !Linear interpolation for values on hexagon edge midpoints
       !Type of interpolation selector
       select case(trim(kindinterpol))
       case('none') !Don't do any interpolation, just get triangle
          i=getnearnode(p, mesh)
          vector_interpol=0._r8
       case('p1nc') !P1 non conformal element (uses triangle)
          vector_interpol=vecinterpol_p1nc(p, var,  mesh)
       case('wach') !Linear interpolation (Modified Wachspress coordinates)
          vector_interpol=vecinterpol_wachspress(p, var,  mesh)
       case default
          print*, "VECTOR INTERPOL ERROR: Don't know how to do this interpolation:", &
               kindinterpol, " with this positioning: ", var%pos
          stop
       end select

    case(4) !Values on Voronoi centroids
       !Type of interpolation selector
       select case(trim(kindinterpol))
       case('none') !Don't do any interpolation, just get triangle
          k=gettr(p, mesh)
          vector_interpol=0._r8
       case('linhxb') !Linear interpolation using Voronoi centroids
          vector_interpol=vecinterpol_linear_hxb(p, var, mesh)
       case('lintrv') !Linear interpolation using Voronoi nodes
          vector_interpol=vecinterpol_linear_trv(p, var, mesh)
       case default
          print*, "VECTOR INTERPOL ERROR: Don't know how to do this interpolation:", &
               kindinterpol, " with this positioning: ", var%pos
          stop
       end select

    case(5) !Values on triangle centroids
       !Type of interpolation selector
       select case(trim(kindinterpol))
       case('none') !Don't do any interpolation, just get triangle
          i=getnearnode(p, mesh)
          vector_interpol=0._r8
       case('wach') !Linear interpolation (wachspress coords)
          vector_interpol=vecinterpol_wachspress(p, var,  mesh)
       case default
          print*, "VECTOR INTERPOL ERROR: Don't know how to do this interpolation:", &
               kindinterpol, " with this positioning: ", var%pos
          stop
       end select

    case default
       print*, "VECTOR INTERPOL ERROR: Don't know how to do this interpolation:", &
            kindinterpol, " with this positioning: ", var%pos
       stop
    end select

    return
  end function vector_interpol

  function vector_reconstruct (p, var, mesh, kinterp, rbf_mat, pindex)
    !-----------------------------------------------------------
    !  VECTOR RECONSTRUCTION
    !
    ! Given the triangulation of a set of nodes on the unit
    !  sphere (mesh), along with vector components (normal) at
    !  the edge midpoints (in var), computes the value at a point (p)
    !  using 'kinterp' interpolation.
    ! The normal components could be given on midpoints of
    !  triangles or hexagons (var$pos = 2 or 3)
    ! If the point p is a node (triangle vertex), the index may
    !  be informed via pindex - this is mandatory for edge reconstruction
    !  such as trisk
    !-------------------------------------------------------
    !Reconstruction point
    real (r8), intent(in) :: p(1:3)

    !Input Variable with vector comp on hexagon midpoints
    type(scalar_field), intent(inout) :: var

    !Mesh in which the variable belongs
    type(grid_structure) :: mesh

    !RBF matrix vector (used only on rbf interpolations)
    type(rbf_matrix_structure), optional :: rbf_mat(:)

    !Kind of reconstruction
    ! % none   = Just find the nearest node
    ! % rbf    = Radial Basis Function
    !     rbfhx  - rbf for 6 point Voronoi cell
    !     rbftr  - rbf for 3 point triangle
    !     rbf*pc- rbf* with constant polynomial
    ! % per    = Perot (2000) recosntruction
    !     perhx  - Perot for Hexagons
    !     pertr  - Perot for Triangles
    !     perpj  - Perot for Hexagons projecting to tang plane
    !     pered  - Perot vector for edges (tangent recon)
    ! % lsq    = Least Square Fit
    !     lsqhxe - 12 point hexagon stencil
    !     lsqtrc - 9 point triangle based stencil
    ! % rt0    = Raviart Thomas 0th order element reconstruction
    ! % dtred  = Dual triangle tangent reconstruction
    ! % wht    = Whitney Edge element reconstruction
    ! % kls    = Klausen reconstruction with Wachspress coords
    ! % trsk   = TRISK tangent vector reconstruction
    character (len=8), optional :: kinterp
    character (len=8):: kindinterpol

    !Index of the point (tr circumcenter, node or edge) for which the recontructin
    !  will be made (avoids having to search for it)
    integer (i4), optional:: pindex
    integer (i4):: pindtmp

    !Interpolation/reconstruction stencil
    character (len=4):: stencil

    !Reconstructed vector
    real (r8):: vector_reconstruct(1:3)

    !Auxiliar index
    integer(i4):: i

    !Set defaulf kind of interpolation (in case kinterp not given)
    if(.not.present(kinterp))then
       kindinterpol="perhx"
    else
       !Local kind of interpolation
       kindinterpol=kinterp
    end if

    if(present(pindex))then
       pindtmp=pindex
    else
       pindtmp=0
    end if

    !Set RBF stencil kind
    if(kindinterpol(1:3)=='rbf')then
       if(.not.present(rbf_mat))then
          print*, "Vector_recon error: rbf_mat must be given!"
          stop
       end if
       stencil(1:4)=kindinterpol(4:7)
       kindinterpol='rbf'
    end if

    if(kindinterpol(1:3)=='per' .and. kindinterpol(1:5) /='pered' )then
       stencil(1:4)=kindinterpol(4:6)
       kindinterpol='per'
    end if

    if(kindinterpol(1:3)=='lsq')then
       stencil(1:4)=kindinterpol(4:7)
       kindinterpol='lsq'
       !print*, "hi", trim(kindinterpol), trim(stencil)
    end if


    select case(var%pos)
    case(2) !Normal given at triangle edge midpoints
       select case(trim(kindinterpol))
       case("none", "nonetr", "nonehx") !Just find the nearst node
          i=getnearnode(p, mesh)
          vector_reconstruct=0._r8
       case("per") !Perot vector reconstruction
          vector_reconstruct=vecrecon_perot(p, var, stencil, mesh, pindtmp)
       case("rt0") !Raviart-Thimas 0th order reconstruction
          vector_reconstruct=vecrecon_rt0(p, var, mesh)
       case("lsq") !Quadratic least square reconstruction
          !   vector_reconstruct=vecrecon_lsq (p, var, mesh)
          print*, trim(kindinterpol), " on TR edges : not implemented"
          stop
       case("p1nc") !P1 non conforming element
          !   vector_reconstruct=vecrecon_p1nc (p, var, mesh)
          print*, trim(kindinterpol), " on TR edges : not implemented"
          stop
       case default
          print*, "vector_reconstruct ERROR: Don't know how to do this interpolation:", &
               kindinterpol, var%pos
          stop
       end select
    case(3) !Normal given at hexagon edge midpoints
       select case(trim(kindinterpol))
       case("none", "nonetr", "nonehx") !Just find the triangle
          i=getnearnode(p, mesh)
          vector_reconstruct=0._r8
       case("per") !Perot (2000) reconstruction
          vector_reconstruct=vecrecon_perot(p, var, stencil, mesh, pindtmp)
          !vector_reconstruct=vecrecon_perot(p, var, stencil, mesh)
       case("rbf") !Radial basis function reconstruction
          vector_reconstruct=vecrecon_rbf(p, var, stencil, rbf_mat, mesh)
       case("wht") !Whitney Edge element - Tangent comp
          vector_reconstruct=vecrecon_whitney(p, var, mesh)
       case("kls") !Klausen reconstruction - Normal components
          vector_reconstruct=vecrecon_klausen(p, var, mesh)
       case("lsq") !Least Square fit
          vector_reconstruct=vecrecon_lsq (p, var, stencil, mesh, pindtmp)
          !vector_reconstruct=vecrecon_lsq (p, var, stencil, mesh)
       case("trsk") !TRISK (Thuburn et al 2009) tangent vector reconstruction
          vector_reconstruct=vecrecon_trsk(pindtmp, var, mesh)
       case("pered") !Perot for edges tangent vector reconstruction
          vector_reconstruct=vecrecon_pered(pindtmp, var, mesh)
       case("dtred") !Dual triangle recons for edges (tangent vector reconstruction)
          vector_reconstruct=vecrecon_dtred(pindtmp, var, mesh)
       case default
          print*, "VECTOR RECON ERROR: Don't know how to do this interpolation:", &
               kindinterpol, var%pos
          stop
       end select
    case(6) !Normal given at hexagon edge at intersection with tr edge midpoint
       select case(trim(kindinterpol))
       case("trsk") !TRISK (Thuburn et al 2009) tangent vector reconstruction
          vector_reconstruct=vecrecon_trsk(pindtmp, var, mesh)
       case("pered") !Perot for edges tangent vector reconstruction
          vector_reconstruct=vecrecon_pered(pindtmp, var, mesh)
       case("dtred") !Dual triangle recons for edges (tangent vector reconstruction)
          vector_reconstruct=vecrecon_dtred(pindtmp, var, mesh)
       case default
          print*, "VECTOR RECON ERROR: Don't know how to do this interpolation:", &
               kindinterpol, var%pos
          stop
       end select
    case default
       print*, "VECTOR RECON ERROR: Don't know how to do this interpolation:", &
            kindinterpol, var%pos
       stop
    end select

    !Force vetor to be tangent to the sphere
    vector_reconstruct=proj_vec_sphere(vector_reconstruct, p)

    return
  end function vector_reconstruct

  function tgrecon_index (ed, kindinterpol, pos, mesh)
    !-----------------------------------------------------------
    !  tangent reconstruction index
    !
    ! Given the triangulation of a set of nodes on the unit
    !  sphere (mesh) and a recon method 'kinterp', calculates
    !  a consistency index for the weights
    !-------------------------------------------------------

    !Mesh in which the variable belongs
    type(grid_structure) :: mesh

    !Kind of reconstruction
    ! % pered  - Perot vector for edges (tangent recon)
    ! % dtred  = Dual triangle tangent reconstruction
    ! % trsk   = TRISK tangent vector reconstruction
    character (len=8):: kindinterpol

    !Index of edge
    integer (i4):: ed

    !Vector positions
    integer(i4):: pos

    !Reconstructed vector
    real (r8):: tgrecon_index

    select case(trim(kindinterpol))
    case("trsk") !TRISK (Thuburn et al 2009) tangent vector reconstruction
       tgrecon_index=trsk_order_index(ed,  pos, mesh)
    case("pered") !Perot for edges tangent vector reconstruction
       tgrecon_index=pered_order_index(ed,  mesh)
    case("dtred") !Dual triangle recons for edges (tangent vector reconstruction)
       tgrecon_index=trsk_order_index(ed,  pos, mesh)
    case default
       print*, "TGRECON INDEX ERROR: Don't know how to do this interpolation:", &
            kindinterpol
       stop
    end select


    return
  end function tgrecon_index


  !  function vec_tg_reconstruct (ed, var, mesh, kinterp, rbf_mat)
  !    !-----------------------------------------------------------
  !    !  VECTOR RECONSTRUCTION OF TANGENT VECTOR
  !    !
  !    ! Given the triangulation of a set of nodes on the unit
  !    !  sphere (mesh), along with vector components (normal) at
  !    !  the edge midpoints (in var), computes the value at the midpoint
  !    !  of edge "ed"  using 'kinterp' method.
  !    ! The normal components should be given on midpoints of
  !    !  hexagons (var$pos = 3)
  !    !-------------------------------------------------------
  !    !Reconstruction point
  !    integer(i4), intent(in):: ed
  !
  !    !Input Variable with vector comp on hexagon midpoints
  !    type(interpolable_variable), intent(inout) :: var
  !
  !    !Mesh in which the variable belongs
  !    type(grid_structure), intent(in) :: mesh
  !
  !    !RBF matrix vector (used only on rbf interpolations)
  !    type(rbf_matrix_structure), optional :: rbf_mat(:)
  !
  !    !Kind of reconstruction
  !    ! trisk - Thuburn et al method
  !    ! pertg - Perot's method
  !    ! rbftg - RBF
  !    character (len=8), optional :: kinterp
  !    character (len=8):: kindinterpol
  !
  !    !Interpolation/reconstruction stencil
  !    character (len=4):: stencil
  !
  !    !Reconstructed vector
  !    real (r8):: vec_tg_reconstruct
  !
  !    !Auxiliar index
  !    integer(i4):: i
  !
  !    !Set defaulf kind of interpolation (in case kinterp not given)
  !    if(.not.present(kinterp))then
  !       kindinterpol="trisk"
  !    else
  !       !Local kind of interpolation
  !       kindinterpol=kinterp
  !    end if
  !
  !    !Set RBF stencil kind
  !    if(kindinterpol(1:3)=='rbf')then
  !       if(.not.present(rbf_mat))then
  !          print*, "Vector_recon error: rbf_mat must be given!"
  !          stop
  !       end if
  !       stencil(1:4)=kindinterpol(4:7)
  !       kindinterpol='rbf'
  !    end if
  !
  !    if(kindinterpol(1:3)=='per')then
  !       stencil(1:4)=kindinterpol(4:6)
  !       kindinterpol='per'
  !    end if
  !
  !    if(kindinterpol(1:3)=='lsq')then
  !       stencil(1:4)=kindinterpol(4:7)
  !       kindinterpol='lsq'
  !       !print*, "hi", trim(kindinterpol), trim(stencil)
  !    end if
  !
  !    if(kindinterpol(1:3)=='lsq')then
  !       stencil(1:4)=kindinterpol(4:7)
  !       kindinterpol='lsq'
  !       !print*, "hi", trim(kindinterpol), trim(stencil)
  !    end if
  !
  !
  !    select case(var%pos)
  !    case(3) !Normal given at hexagon edge midpoints
  !       select case(trim(kindinterpol))
  !       case("none", "nonetr", "nonehx") !Just find the triangle
  !          !   i=getnearnode(p, mesh)
  !          !   vector_reconstruct=0._r8
  !       case("trisk") !Perot (2000) reconstruction
  !          vec_tg_reconstruct=vec_tg_recon_trisk(ed, var, mesh)
  !          !case("per") !Perot (2000) reconstruction
  !          !  vector_reconstruct=vecrecon_perot(p, var, stencil, mesh, pindtmp)
  !          !vector_reconstruct=vecrecon_perot(p, var, stencil, mesh)
  !          !case("rbf") !Radial basis function reconstruction
  !          !   vector_reconstruct=vecrecon_rbf(p, var, stencil, rbf_mat, mesh)
  !          !case("lsq") !Least Square fit
  !          !   vector_reconstruct=vecrecon_lsq (p, var, stencil, mesh, pindtmp)
  !          !vector_reconstruct=vecrecon_lsq (p, var, stencil, mesh)
  !       case default
  !          print*, "VECTOR RECON ERROR: Don't know how to do this interpolation:", &
  !               kindinterpol, var%pos
  !          stop
  !       end select
  !    case default
  !       print*, "VECTOR RECON ERROR: Don't know how to do this interpolation:", &
  !            kindinterpol, var%pos
  !       stop
  !    end select
  !
  !    !Force vetor to be tangent to the sphere
  !    !vec_tg_reconstruct=proj_vec_sphere(vec_tg_reconstruct, mesh%ed(ed)%c%p)
  !
  !    return
  !  end function vec_tg_reconstruct

  !subroutine getmultirecon(kinterp, recon, interp, interp_ho, massc)
  subroutine getmultirecon(kinterp, recon_mtd)
    !------------------------------------------------
    ! Get the multiple vector reconstruction method
    !  to be used
    !
    ! Recieves stringread and separates into 3 sub strings
    ! kinterp must have sub strings separated by "+" and a number
    !  1-reconstruction method
    !  2-inteprolation method
    !  3-higher order method (if hybrid)
    ! If +b is appended at the end, do reconstruction to cell
    !  mass centers
    !
    ! Output:
    ! recon_mtd%
    ! recon = kind of reconstruction remap
    ! interp= kind of local interpolation method
    ! horecon = High order interpolation/reconstruct method
    ! kinterp receives a clean name (no + and numbers)
    !------------------------------------------------

    !Kind of interpolation to be used
    character (len=64):: kinterp !Full interpolation name
    type(vectorinterpol_methods):: recon_mtd

    !character (len=16):: recon !Decomposed names
    !character (len=16):: interp !Decomposed names
    !character (len=16):: horecon !Decomposed names
    !logical, optional :: massc

    !Aux vars
    integer:: i
    integer:: ini
    integer:: lst
    integer:: k
    integer:: n
    character (len=1):: stmp
    character (len=1):: skind

    recon_mtd%recon=""
    recon_mtd%interp=""
    recon_mtd%horecon=""

    !Counters
    ini=1
    k=0
    n=len(trim(kinterp))
    lst=n
    !Loop over chars
    do i=1, n
       stmp=kinterp(i:i)
       !print*, i, stmp
       !If a separation char is detected
       if(trim(stmp)=="+")then
          !If a full method is found
          if(k>0)then
             lst=i-1
             select case(trim(skind))
             case("1") !Reconstruction method
                recon_mtd%recon=kinterp(ini:lst)
             case("2") !Interpolation method
                recon_mtd%interp=kinterp(ini:lst)
             case("3") !Higher order method
                recon_mtd%horecon=kinterp(ini:lst)
             end select
          end if
          k=k+1
          ini=i+2
          skind=kinterp(i+1:i+1)
       end if
    end do

    !if(present(massc))then
    if(kinterp(n:n)=="b")then
       recon_mtd%massc=.true.
       !kinterp=trim(kinterp)//"_b"
    else
       recon_mtd%massc=.false.
    end if
    !end if

    if(k==0)then !Single method to be used
       recon_mtd%recon=trim(kinterp)
       recon_mtd%interp=trim(kinterp)
       recon_mtd%name=trim(recon_mtd%recon)//trim(recon_mtd%horecon)
       kinterp=recon_mtd%name
    elseif(len(trim(recon_mtd%interp)) == 0)then !Use recon as interpolation
       !Store recon mtd string as the interpolation method
       recon_mtd%interp=recon_mtd%recon
       recon_mtd%name=trim(recon_mtd%recon)//trim(recon_mtd%horecon)
       kinterp=recon_mtd%name
    else
       recon_mtd%name=trim(recon_mtd%recon)//trim(recon_mtd%interp)//trim(recon_mtd%horecon)
       kinterp=recon_mtd%name
    end if

    !if(present(massc))then
    if(recon_mtd%massc)then
       kinterp=trim(kinterp)//"_b"
       recon_mtd%name=trim(recon_mtd%name)//"_b"
    end if
    !end if

    if(trim(recon_mtd%horecon)=="" .or. &
         trim(recon_mtd%horecon)==trim(recon_mtd%recon))then
       recon_mtd%hyb=.false.
    else
       recon_mtd%hyb=.true.
    end if

    return
  end subroutine getmultirecon

  subroutine vrec_remap (vn, vr, mesh, recon_mtd, rbf_mat2)
    !-----------------------------------------------------------
    !  VECTOR RECONSTRUCTION REMAPPING
    !
    ! Reconstructs a vector field on the whole grid
    !
    ! vn - Normal components could be given on midpoints of
    ! triangles or hexagons (vn%pos = 2 or 3)
    !
    ! vr - Reconstructed vector field, it will be positioned
    !     depending on the reconstruction wanted or
    !     and the recon_mtd%massc flag
    !
    ! recon_mtd%
    ! recon - kind of reconstruction to be used (see vector_reconstruct)
    ! horecon - Second kind of reconstruction to be used (optional)
    !   to be used associated with mask_ho which indicates true where
    !    this second method must be used
    !
    ! massc - true indicates that the reconstruction will be made to
    !       hx centroids, and not nodes, or triangle centroids,
    !       and not circumcenters
    ! hyb - true for hybrid schemes
    ! alignlimit - alignment index threshold for hybrid method
    ! nonalignedonly - Reconstruct only on nonaligned cells
    !-------------------------------------------------------

    !Input Variable with vector normal comp on hexagon midpoints
    type(scalar_field):: vn

    !Reconstructed vector field  in cartesian coord
    type(vector_field_cart):: vr

    !Mesh in which the variable belongs
    type(grid_structure), intent(in) :: mesh

    !RBF matrix vector (used only on rbf interpolations)
    type(rbf_matrix_structure), allocatable :: rbf_mat(:)
    type(rbf_matrix_structure), optional :: rbf_mat2(:)

    !Methods to be used for the reconstruction
    type(vectorinterpol_methods) :: recon_mtd

    !Interpolation point
    real (r8):: p(1:3)

    !RBF shape parameter
    !real (r8), optional:: rbf_par
    real (r8):: rbf_h

    !Auxiliar vars
    integer (i4):: i, j, k
    logical :: isaligned

    !RBF stencil and kind of interpolation
    character(len=4):: stencil

    !RBF parameters
    rbf_h=recon_mtd%rbf_par

    !Kind of reconstruction
    select case (trim(recon_mtd%recon))
    case ("none", "nonehx", "rbfhx", "perhx", "perpj", "lsqhxe", "kls", &
         "rbfhxpc" )
       if(recon_mtd%massc)then
          !Reconstruct vectors to Voronoi cell centroid
          vr%pos=4
       else
          !Reconstruct vectors to triangle vertices/hexagon centers
          vr%pos=0
       end if
       vr%n=mesh%nv
    case ("nonetr", "pertr", "lsqtrc", "rt0", "wht", &
         "rbftr", "rbfetr", &
         "rbftrpc", "rbfetrpc", "rbfetrp")
       if(recon_mtd%massc)then
          !Reconstruct vector to triangle barycenters
          vr%pos=5
       else
          !Reconstruct vector to triangle circumcenters /hexagon vertices
          vr%pos=1
       end if
       vr%n=mesh%nt
    case ("nonetg", "pertg", "rbftg", "trsk", "trskb")
       !if(recon_mtd%massc)then
       !   !Reconstruct to the midpoint of the Voronoi cell edge
       vr%pos=3
       !else
       !   !Reconstruct to the midpoint of the triangular cell edge
       !   vr%pos=2
       !end if
       vr%n=mesh%ne
    case default
       print*, "vrec_remap error: unknown method ", recon_mtd%recon
       stop
    end select

    !Allocate space for reconstructed vectors if needed
    if(.not.allocated(vr%p))then
       allocate(vr%p(1:vr%n))
    elseif(ubound(vr%p, 1) /= vr%n) then
       deallocate(vr%p)
       allocate(vr%p(1:vr%n))
    end if

    !If RBF, calculate the matrix structure (precalculation, depends only on mesh)
    if(recon_mtd%recon(1:3)=="rbf")then
       stencil=trim(recon_mtd%recon(4:7))
       if(.not.present(rbf_mat2))then
          call rbfvec_matrix_build( stencil, rbf_h, rbf_mat, mesh )
       end if
    end if

    !Reconstruct to cell centers, tr centers or edges
    do i=1, vr%n
       !Do reconstruction based on vecreconpos
       select case(vr%pos)
       case(0) !Reconstruc to tr vertices
          p=mesh%v(i)%p
       case(1) !Reconstruction to tr circuncenters
          p=mesh%tr(i)%c%p
       case(2) !Reconstruction to tr edge center
          p=mesh%ed(i)%c%p
       case(3) !Reconstruction to hx edge center
          p=mesh%edhx(i)%c%p
       case(4) !Reconstruct to Voronoi cell centroid
          p=mesh%hx(i)%b%p
       case(5) !Reconstruc to tr barycenter
          p=mesh%tr(i)%b%p
       case default
          print*, "Error on vrec_remap: unknown recon position ", vr%pos
          stop
       end select

       select case(vr%pos)
       case(0,4) !Voronoi cells
          if( recon_mtd%hyb ) then
             !Check for hybrid method
             isaligned= (mesh%hx(i)%align < recon_mtd%alignlimit)
          else
             isaligned=.true.
          end if
       case(1,5) !Triangles
          !Hybrid method for triangular cells
          if( recon_mtd%hyb ) then
             !Check for hybrid method
             !If any of the triangle vertices is
             !  nonaligned use higherorder method
             isaligned=.false.
             do j=1, 3
                k=mesh%tr(i)%v(j)
                isaligned = (isaligned .or. (mesh%hx(k)%align < recon_mtd%alignlimit))
             end do
          else
             isaligned=.true.
          end if
       case(2,3) !Edges
          isaligned=.true.
       end select

       !Check for hybrid method
       if(.not.present(rbf_mat2))then
          if( isaligned ) then
             vr%p(i)%v=vector_reconstruct (p, vn, mesh, recon_mtd%recon, rbf_mat, i)
          else
             vr%p(i)%v=vector_reconstruct (p, vn, mesh, recon_mtd%horecon, rbf_mat, i)
          end if
       else
          if( isaligned ) then
             vr%p(i)%v=vector_reconstruct (p, vn, mesh, recon_mtd%recon, rbf_mat2, i)
          else
             vr%p(i)%v=vector_reconstruct (p, vn, mesh, recon_mtd%horecon, rbf_mat2, i)
          end if
       end if
    end do

    return
  end subroutine vrec_remap

  !========================================================================
  !    SCALAR REMAPPING ROUTINES
  !========================================================================

  subroutine scalar_remap_geo2ll(nlon, nlat, var, mesh, varunif, kindinterpol, rbf_mat)
    !---------------------------------------------------------------------
    ! Scalar Remap Geodesic to Lon Lat
    !
    !   Given (nlat, nlon) sizes, a scalar variable (var) already containing
    !   gradient estimatives and/ or the rbf matrices, and a mesh,
    !   this routine will interpolate
    !   a uniform global (nlat, nlon) mesh
    !---------------------------------------------------------------------
    !Uniform grid sizes
    integer (i4), intent(in) :: nlon
    integer (i4), intent(in) :: nlat

    !Variable for interpolation procedure
    type(scalar_field), intent(inout) :: var

    !Mesh in which the variable belongs
    type(grid_structure), intent(in) :: mesh

    !Output variable for uniform grid
    real (r8), intent(out) :: varunif(1:nlon, 1:nlat)

    !Kind of interpolation -see scalar_interpol
    character (len=8):: kindinterpol

    !RBF matrix vector (used only on rbf interpolations)
    type(rbf_matrix_structure), optional :: rbf_mat(:)

    !Auxiliar varibles
    integer (i4):: i
    integer (i4):: j
    real (r8):: tlat
    real (r8):: tlon
    real (r8):: dlat
    real (r8):: dlon
    real (r8):: p(1:3)
    real (r8)::  fest

    dlat=180._r8/real(nlat, r8)
    dlon=360._r8/real(nlon, r8)
    !Pixel registration mode (GMT)
    tlat=-90+dlat/2._r8
    lat: do j=1,nlat
       tlon=-180+dlon/2._r8
       lon: do i=1,nlon
          call sph2cart(tlon*deg2rad, tlat*deg2rad, p(1), p(2), p(3))
          if(present(rbf_mat))then
             fest=scalar_interpol(p, var, mesh, kindinterpol, rbf_mat)
          else
             fest=scalar_interpol(p, var, mesh, kindinterpol)
          end if
          varunif(i,j)=fest
          tlon=tlon+dlon
       end do lon
       tlat=tlat+dlat
    end do lat

    return
  end subroutine scalar_remap_geo2ll

  subroutine scalar_remap_trc2trv(var_trc, var_trv, mesh)
    !----------------------------------------------------
    !Scalar remapping routine TRC 2 TRV
    !
    !   Interpolates values on triangles vertices using
    !    values on triangle's circumcenters
    !
    !   The interpolation used is Wachspress coordinate based
    !-----------------------------------------------------

    !Mesh structure
    type(grid_structure), intent(in) :: mesh

    !Variable with values on tr circumcenters
    type(scalar_field), intent(inout) :: var_trc

    !Variable with values on tr vertices
    type(scalar_field), intent(out) :: var_trv

    !Kind of interpolation
    character (len=8):: kinterp

    !Indexes
    integer (i4):: i

    ! Check if values on triangle center,
    if(var_trc%pos/=1)then
       print*, "scalar_remap_trc2trv error: scalar field not on triangle centers"
       stop
    end if

    !Build new variable attributes
    var_trv%pos=0
    var_trv%n=mesh%nv
    var_trv%name=var_trc%name
    allocate(var_trv%f(1:var_trv%n))

    !Set kind of interpolation
    kinterp = "wach"

    !Do remapping
    do i=1, mesh%nv
       !Simples average
       !var_trv%f(i)=0
       !do j=1, mesh%v(i)%nnb
       !   k=mesh%v(i)%tr(j)
       !   var_trv%f(i)=var_trv%f(i)+var_trc%f(k)
       !end do
       !var_trv%f(i)=var_trv%f(i)/mesh%v(i)%nnb
       !p=mesh%v(i)%p
       var_trv%f(i)=scalar_interpol (mesh%v(i)%p, var_trc, mesh, kinterp)
    end do

    return
  end subroutine scalar_remap_trc2trv

  !========================================================================
  !    VECTOR REMAPPING ROUTINES
  !========================================================================

  subroutine vec_remap_ed2trc(ved, vtrc, mesh)
    !----------------------------------------------------
    !Vector remapping routine ED2TRC
    !   Interpolates vectors on triangles circumcenters
    !      using the values on the edges midpoints (tr or hexag)
    !   The interpolation is linear and the resulting vector
    !      is projected to the tangent plane at the end of the
    !      interpolation
    !-----------------------------------------------------
    !Input Variable with vector on hexagon edges midpoints
    type(vector_field_cart), intent(in) :: ved

    !Output Variable with vector on triangle circumcenters
    type(vector_field_cart), intent(inout) :: vtrc

    !Mesh in which the variable belongs
    type(grid_structure), intent(in) :: mesh

    !Kind of interpolation
    character (len=8):: kinterp

    !Counters
    integer (i4):: i

    !Set kind of interpolation
    kinterp = "p1nc"

    vtrc%pos=1
    vtrc%n=mesh%nt
    if(.not.allocated(vtrc%p))then
       allocate(vtrc%p(1:vtrc%n))
    end if

    !For each triangle
    do i=1, mesh%nt
       vtrc%p(i)%v=vector_interpol(mesh%tr(i)%c%p, ved, mesh, kinterp)
    end do

    return
  end subroutine vec_remap_ed2trc


  subroutine vec_remap_trc2ed(vtrc, ved, mesh)
    !----------------------------------------------------
    !Vector remapping routine
    !   Interpolates vector on triangles edges midpoints using
    !    vectors on triangle's circumcenters
    !   "ved%pos" must be given implicitly in order to
    !   specify which kind of edge will be remapped (tr or hexag)
    !   ved%pos=2 : triangle edges
    !   ved%pos=3 : hexag edges
    !
    !   The interpolation is linear and the resulting vector
    !      is projected to the tangent plane at the end of the
    !      interpolation
    !-----------------------------------------------------

    !Input Variable with vector on triangle circumcenters
    type(vector_field_cart), intent(in) :: vtrc

    !Output Variable with vector on triangle edges midpoints
    type(vector_field_cart), intent(inout) :: ved

    !Mesh in which the variable belongs
    type(grid_structure), intent(in) :: mesh

    !Counters
    integer (i4):: i

    !Check if vector field is on triangle circumcenters
    if(vtrc%pos /= 1)then
       print*, "Error on vec_remap_trc2ed:"
       print*, "    Input vector field not on tr cc, v%pos", vtrc%pos
       stop
    end if

    !ved%pos must be given
    ved%n=mesh%ne
    if(.not.allocated(ved%p))then
       allocate(ved%p(1:ved%n))
    end if

    if(ved%pos/=2 .or. ved%pos/=3)then
       print*, "Warning on vec_remap_trc2ed:"
       print*, "    Unknown edge position, ved%pos=", ved%pos
       print*, "    using vector on hexagon edges"
       ved%pos=3
    end if

    do i=1, mesh%ne
       ved%p(i)%v=vecinterpol_linear_trc2ed(i, vtrc, ved%pos, mesh)
    end do

    return
  end subroutine vec_remap_trc2ed

  subroutine vec_remap_trv2trc(vtrv, vtrc, mesh)
    !----------------------------------------------------
    !Vector remapping routine
    !
    !   Interpolates vector on triangles circumcenters using
    !    vectors on triangle's vertices
    !
    !   The interpolation is linear and the resulting vector
    !      is projected to the tangent plane at the end of the
    !      interpolation
    !-----------------------------------------------------

    !Input Variable with vector on triangle vertices
    type(vector_field_cart), intent(in) :: vtrv

    !Output Variable with vector on triangle circumcenters
    type(vector_field_cart), intent(inout) :: vtrc

    !Mesh in which the variable belongs
    type(grid_structure), intent(in) :: mesh

    !Counters
    integer (i4):: k

    !Kind of interpolation
    character (len=8):: kinterp

    !Set kind of interpolation
    kinterp = "lintrv"

    !Check position of vector field
    if(vtrv%pos/=0)then
       print*, "Error on vec_remap_trv2trc:"
       print*, "    Vectors not on triangle vertices, vtrv%pos=", vtrv%pos
       stop
    end if

    !Set space for net vector field
    vtrc%n=mesh%nt
    vtrc%pos=1
    if(.not.allocated(vtrc%p))then
       allocate(vtrc%p(1:vtrc%n))
    end if

    !Interpolate for each triangle
    do k=1, mesh%nt
       vtrc%p(k)%v=vector_interpol (mesh%tr(k)%c%p, vtrv, mesh, kinterp)
    end do

    return
  end subroutine vec_remap_trv2trc

  subroutine vec_remap_trv2ed(vtrv, ved, mesh)
    !----------------------------------------------------
    !Vector remapping routine
    !
    !   Interpolates vector on triangle edges midpoint using
    !    vectors on triangle's vertices
    !   ONLY FOR TRIANGLE EDGES
    !
    !   The interpolation is linear and the resulting vector
    !      is projected to the tangent plane at the end of the
    !      interpolation
    !-----------------------------------------------------

    !Input Variable with vector on triangle vertices
    type(vector_field_cart), intent(in) :: vtrv

    !Output Variable with vector on triangle edges midpoints
    type(vector_field_cart), intent(inout) :: ved

    !Mesh in which the variable belongs
    type(grid_structure), intent(in) :: mesh

    !Counters
    integer (i4):: k

    !Check position of vector field
    if(vtrv%pos/=0)then
       print*, "Error on vec_remap_trv2ed:"
       print*, "    Vectors not on triangle vertices, vtrv%pos=", vtrv%pos
       stop
    end if

    !Set space for vector field
    ved%n=mesh%ne
    if(.not.allocated(ved%p))then
       allocate(ved%p(1:ved%n))
    end if
    ved%pos=2

    !Interpolated for each edge
    do k=1, mesh%ne
       ved%p(k)%v=vecinterpol_linear_trv2ed(k, vtrv, ved%pos, mesh)
    end do

    return
  end subroutine vec_remap_trv2ed



  !========================================================================
  !    LINEAR / TRIANGLE BARYCENTRIC COORDS INTERPOLATION
  !========================================================================

  function scinterpol_linear_trv (p, var, mesh)
    !-----------------------------------------------------------
    !  scinterpol_LINEAR_TRV
    !
    !   Given the triangulation of a set of nodes on the unit
    ! sphere (mesh), along with data values at the nodes/vertices
    ! (in var), computes the value at a point (p)
    ! of a continuous function which interpolates the data values.
    ! !The interpolatory function is linear on each underlying triangle
    ! (planar triangle with the same vertices as a spherical
    ! triangle).
    !
    ! On input:
    !       p         =  point in cart coords
    !       mesh      =  mesh structure
    !       var       =  interpolable_variable structure
    !                      with var%pos=0
    ! Returns the linear interpolated value
    !-------------------------------------------------------
    !Interpolation point
    real (r8), intent(in) :: p(1:3)

    !Variable for interpolation procedure
    type(scalar_field), intent(in) :: var

    !Mesh in which the variable belongs
    type(grid_structure), intent(in) :: mesh

    !Value at interpolated point
    real (r8):: f
    real (r8):: scinterpol_linear_trv

    !loop index, indexes for nodes, triangle index
    integer(i4):: i
    integer(i4):: v(1:3)
    integer(i4):: k

    !Cartesian and barycentric coordinates of the points
    real (r8):: b(1:3)

    ! Locate p with respect to the triangulation.
    k=gettr(p, mesh)

    !Get barycentric coords
    b=bar_coord_tr(p, k, mesh)

    !Set the value on point
    v(1:3)=mesh%tr(k)%v(1:3)

    f = 0._r8
    do i=1,3
       f=f+ b(i)*var%f(v(i))
    end do
    scinterpol_linear_trv=f

    return
  end function scinterpol_linear_trv

  function vecinterpol_linear_trv(p, v, mesh)
    !----------------------------------------------------
    !Vector linear interpolation function
    !      TRIANGLE VERTICE (HX CENTER) -> TRIANGLE POINT
    !
    !   Interpolates vector on a triangles point "p" using
    !    vectors on triangle's vertices (hexag centers)
    !
    !   The interpolation is linear and the resulting vector
    !      is projected to the tangent plane at the end of the
    !      interpolation
    !-----------------------------------------------------
    !Triangle point to be interpolated
    real(r8), dimension(1:3), intent(in) :: p

    !Input Variable with vector on triangle vertices
    type(vector_field_cart), intent(in) :: v

    !Output Variable with vector on triangle circumcenters
    real(r8), dimension(1:3) :: vecinterpol_linear_trv

    !Mesh in which the variable belongs
    type(grid_structure), intent(in) :: mesh

    !Temp. vector
    real (r8)::  vec(1:3)

    !Barycentric coords. (planar)
    real (r8):: b(1:3)

    !Counters
    integer (i4):: i
    integer (i4):: tr

    !Get triangle which point belongs
    tr=gettr(p, mesh)

    !Get barycentric coords of point
    b=bar_coord_tr(p, tr, mesh)

    !Sum the  contribuitions of vector components on vertices
    vec(1:3)=0.0_r8
    do i=1,3
       vec=vec + b(i)*v%p(mesh%tr(tr)%v(i))%v
    end do
    !print*, b
    !print*, sum(b)
    !Project vector to the sphere
    vecinterpol_linear_trv=proj_vec_sphere(vec, p)

    return
  end function vecinterpol_linear_trv

  function vecinterpol_linear_hxb(p, v, mesh)
    !----------------------------------------------------
    !Vector linear interpolation function
    !      TRIANGLE VERTICE (HX CENTER) -> TRIANGLE POINT
    !
    !   Interpolates vector on a triangles point "p" using
    !    vectors on triangle's vertices (hexag centers)
    !
    !   The interpolation is linear and the resulting vector
    !      is projected to the tangent plane at the end of the
    !      interpolation
    !-----------------------------------------------------
    !Triangle point to be interpolated
    real(r8), dimension(1:3), intent(in) :: p

    !Input Variable with vector on triangle vertices
    type(vector_field_cart), intent(in) :: v

    !Output Variable with vector on triangle circumcenters
    real(r8), dimension(1:3) :: vecinterpol_linear_hxb

    !Mesh in which the variable belongs
    type(grid_structure), intent(in) :: mesh

    !Temp. vector
    real (r8)::  vec(1:3)

    !Barycentric coords. (planar)
    real (r8), dimension(1:3):: b
    real (r8), dimension(1:3):: p1
    real (r8), dimension(1:3):: p2
    real (r8), dimension(1:3):: p3

    !Counters
    integer (i4):: i
    integer (i4):: tr

    !Get triangle which point belongs
    tr=gettr(p, mesh)

    !Get barycentric coords of point
    if(v%pos/=4)then
       print*, "vecinterpol_linear_hxb error: variable not on correct pos", v%pos
    end if

    p1=mesh%hx(mesh%tr(tr)%v(1))%b%p
    p2=mesh%hx(mesh%tr(tr)%v(2))%b%p
    p3=mesh%hx(mesh%tr(tr)%v(3))%b%p

    b=bar_coord(p, p1, p2, p3)

    !Sum the  contribuitions of vector components on vertices
    vec(1:3)=0.0_r8
    do i=1,3
       vec=vec + b(i)*v%p(mesh%tr(tr)%v(i))%v
    end do

    !Project vector to the sphere
    vecinterpol_linear_hxb=proj_vec_sphere(vec, p)

    return
  end function vecinterpol_linear_hxb

  function vecinterpol_linear_trc2ed(ed, vtrc, edpos, mesh)
    !----------------------------------------------------
    !Vector interpolation function
    !      TRIANGLE CIRCUMCENTER -> EDGE
    !
    !   Interpolates vector on triangles edges midpoints using
    !    vectors on triangle's circumcenters
    !   "edpos" must be given in order to
    !   specify which kind of edge will be interpoled (tr or hexag)
    !   edpos="tr" : triangle edges
    !   edpos="hx" : hexag edges
    !
    !   Returns
    !     - vecinterpol_linear_trc2ed: interpolated vector in cartesian coords
    !
    !   The interpolation is linear and the resulting vector
    !      is projected to the tangent plane at the end of the
    !      interpolation
    !-----------------------------------------------------

    !Input Variable with edge index
    integer (i4), intent(in) :: ed

    !Input Variable with vector on triangle circumcenters
    type(vector_field_cart), intent(in) :: vtrc

    !Output Variable with vector on triangle edges midpoints
    real(r8), dimension(1:3)  :: vecinterpol_linear_trc2ed

    !Edge type : triangle (2) or hexag (3)
    integer (i4), intent(in) :: edpos

    !Mesh in which the variable belongs
    type(grid_structure), intent(in) :: mesh

    !Temp. vectors
    real (r8), dimension(1:3) ::  vec
    real (r8), dimension(1:3) :: p
    real (r8):: d1
    real (r8):: d2

    !Counters
    integer (i4):: i
    integer (i4):: j
    integer (i4):: k

    !Check if vector field is on triangle circumcenters
    if(vtrc%pos /= 1)then
       print*, "Error on vec_interpol_trc2ed:"
       print*, "    Input vector field not on tr cc, vtrc%pos", vtrc%pos
       stop
    end if

    !Use "i" for edge index
    i=ed

    !Get triangles sharing this edge
    j=mesh%ed(i)%sh(1)
    k=mesh%ed(i)%sh(2)

    if(edpos==2)then !Triangle edge
       !Get edge midpoint
       p=mesh%ed(i)%c%p
    elseif(edpos==3)then !Hexagon edge
       !Get edge midpoint
       p=mesh%edhx(i)%c%p
    else
       print*, "Error on vec_interpol_trc2ed:"
       print*, "    Unknown edge position, edpos=", edpos
       stop
    end if

    !Set distances to edge midpoint,
    !  used as weights in linear interpol
    d1=arclen(mesh%tr(j)%c%p, p)
    d2=arclen(mesh%tr(k)%c%p, p)

    !IMPORTANT: In case os triangle edge, d1=d2
    !  For hx edge, this might not happen
    if(abs(d1+d2)<eps)then
       print*, "Error on vec_interpol_trc2ed:"
       print*, "  Degenerated triangle (cc too close to edge), ed=", i, d1, d2
    end if

    !Linear interpolating vector
    vec=(d1*vtrc%p(j)%v+d2*vtrc%p(k)%v)/(d1+d2)

    !Project to the sphere
    vecinterpol_linear_trc2ed=proj_vec_sphere(vec, p)

    return
  end function vecinterpol_linear_trc2ed

  function vecinterpol_linear_trv2ed(ed, vtrv, edpos, mesh)
    !----------------------------------------------------
    !Vector interpolation function
    !      TRIANGLE VERTICE -> EDGE (triangle or Hex)
    !
    !   Interpolates vector on triangles edges midpoints using
    !    vectors on triangle's circumcenters
    !   "edpos" must be given in order to
    !   specify which kind of edge will be interpoled (tr or hexag)
    !   edpos="2" : triangle edges
    !   edpos="3" : hexag edges
    !
    !   Returns
    !     - vecinterpol_linear_trv2ed: interpolated vector in cartesian coords
    !
    !   The interpolation is linear and the resulting vector
    !      is projected to the tangent plane at the end of the
    !      interpolation
    !-----------------------------------------------------

    !Input Variable with edge index
    integer (i4), intent(in) :: ed

    !Input Variable with vector on triangle circumcenters
    type(vector_field_cart), intent(in) :: vtrv

    !Output Variable with vector on triangle edges midpoints
    real(r8), dimension(1:3)  :: vecinterpol_linear_trv2ed

    !Edge type : triangle (2) or hexag (3)
    integer (i4), intent(in) :: edpos

    !Mesh in which the variable belongs
    type(grid_structure), intent(in) :: mesh

    !Temp. vectors
    real (r8), dimension(1:3) ::  vec
    real (r8), dimension(1:3) :: p


    !Kind of interpolation
    character (len=8):: kinterp

    !Counters
    integer (i4):: i
    integer (i4):: j
    integer (i4):: k

    !Check if vector field is on triangle circumcenters
    if(vtrv%pos /= 1)then
       print*, "Error on vec_interpol_trc2ed:"
       print*, "    Input vector field not on trv, vtrv%pos", vtrv%pos
       stop
    end if

    !Use "i" for edge index
    i=ed

    !Get  sharing this edge
    j=mesh%ed(i)%v(1)
    k=mesh%ed(i)%v(2)

    if(edpos==2)then !Triangle edge
       !Get edge midpoint
       p=mesh%ed(i)%c%p
       !Do just simple average
       vec=(0.5*vtrv%p(j)%v+0.5*vtrv%p(k)%v)
       vec=proj_vec_sphere(vec, p)
    elseif(edpos==3)then !Hexagon edge
       !Get edge midpoint
       p=mesh%edhx(i)%c%p
       !The point might not be on the edge, so do
       ! linear interpolation
       kinterp="lintrv"
       vec=vector_interpol(p, vtrv, mesh, kinterp)
    else
       print*, "Error on vec_interpol_trv2ed:"
       print*, "    Unknown edge position, edpos=", edpos
       stop
    end if

    !Return vector
    vecinterpol_linear_trv2ed=vec

    return
  end function vecinterpol_linear_trv2ed

  !========================================================================
  !    HERMITE - C1 - RENKAS INTERPOLATION
  !========================================================================

  function scinterpol_hermite_trv(p, var, mesh, monot)
    !-------------------------------------------------------------------
    ! scinterpol_HERMITE_TRV
    !
    !   Given a triangulation of a set of nodes/vertices on the unit
    ! sphere, along with data values and gradients at the nodes,
    ! this routine computes a value f(p), where f interpolates
    ! the nodal data and is once-continuously differentiable
    ! over the convex hull of the nodes.  Refer to function FVAL
    ! for further details.  Renkas method.
    !
    ! Input:
    !       p    = interpolation point
    !       mesh = a mesh structure
    !       var = variable to be interpolated, it should contain
    !           gradients estimatives in var%grad, and var%pos=0
    !
    ! Returns the interpolated value
    !-----------------------------------------------------------------------
    !Interpolation point
    real (r8), intent(in) :: p(1:3)

    !Variable for interpolation procedure
    type(scalar_field), intent(inout) :: var

    !Mesh in which the variable belongs
    type(grid_structure), intent(in) :: mesh

    !Value at interpolated point
    real (r8):: scinterpol_hermite_trv

    !Monotonicity flag
    integer(i4), optional :: monot
    integer(i4) :: monottmp

    !Indexes for nodes and triangle
    integer(i4):: v(1:3)
    integer(i4):: i, k

    !Cartesian and barycentric coordinates of the points
    real (r8)::  b(1:3)

    !Triangle nodes and gradients
    real (r8), dimension(1:3) :: p1
    real (r8), dimension(1:3) :: p2
    real (r8), dimension(1:3) :: p3
    real (r8), dimension(1:3) :: g1
    real (r8), dimension(1:3) :: g2
    real (r8), dimension(1:3) :: g3

    !Min and max values on triangle
    real (r8) :: minf, maxf

    !Check if gradients exist
    if(.not.allocated(var%g))then
       call gradcalc(var, mesh)
    end if

    if(present(monot))then
       monottmp=monot
    else
       monottmp=0
    end if

    ! Locate p with respect to the triangulation.
    k=gettr(p, mesh)

    !Get barycentric coords
    b=bar_coord_tr(p, k, mesh)

    !Save vertex indexes
    v=mesh%tr(k)%v

    ! p is contained in the triangle (v(1),v(2),v(3))
    p1=mesh%v(v(1))%p
    p2=mesh%v(v(2))%p
    p3=mesh%v(v(3))%p
    g1=var%g(v(1))%v
    g2=var%g(v(2))%v
    g3=var%g(v(3))%v

    !Do interpolation
    scinterpol_hermite_trv= fval (b, p1, p2, p3, &
         var%f(v(1)), var%f(v(2)), var%f(v(3)), &
         g1, g2, g3)

    if(monottmp==1)then
       minf=100000000.
       maxf=-10000000.
       do i=1,3
          maxf=max(maxf,var%f(v(i)))
          minf=min(minf,var%f(v(i)))
       end do
       if(scinterpol_hermite_trv>maxf)then
          scinterpol_hermite_trv=maxf
       elseif(scinterpol_hermite_trv<minf)then
          scinterpol_hermite_trv=minf
       end if
    end if

    return

  contains

    function fval (b, v1, v2, v3, f1, f2, f3, g1, g2,g3)
      !--------------------------------------------------------
      ! FVAL  from ssrfpack by R. Renka
      !
      !   Given data values and gradients at the three vertices of
      ! a spherical triangle containing a point p, this routine
      ! computes the value of f at p where f interpolates the ver-
      ! tex data.  Along the triangle sides, the interpolatory
      ! function f is the hermite interpolatory spline
      ! defined by the values and tangential gradient components
      ! at the endpoints, and the gradient component normal to the
      ! triangle side varies linearly with respect to arc-length
      ! between the normal gradient components at the endpoints.
      ! A first-order C-1 blending method is used on the underly-
      ! ing planar triangle.  Since values and gradients on an arc
      ! depend only on the vertex data, the method results in c-1
      ! continuity when used to interpolate over a triangulation.
      !
      !   The blending method consists of taking f(p) to be a
      ! weighted sum of the values at pp of three univariate her-
      ! mite interpolatory splines defined on the line
      ! segments which join the vertices to the opposite sides and
      ! pass through pp:  the central projection of p onto the
      ! underlying planar triangle.
      !
      ! Input:
      !       b = barycentric coordinates of pp with re-
      !           spect to the (planar) underlying triangle
      !           (v1,v2,v3), where pp is the central
      !           projection of p onto this triangle.
      !       v1,v2,v3 = cartesian coordinates of the vertices of
      !                  a spherical triangle containing p.
      !                  v3 left v1->v2 (Counter-clockwise)
      !       f1,f2,f3 = data values associated with the vertices.
      !       g1,g2,g3 = gradients associated with the vertices.
      !                  gi is orthogonal to vi for i = 1,2,3.
      !
      ! Output:
      !       fval = interpolated value at p.
      !--------------------------------------------------------------
      !Barycentric Coords
      real (r8), dimension (1:3), intent(in) ::  b
      !Triangle nodes
      real (r8), dimension (1:3), intent(in) ::  v1
      real (r8), dimension (1:3), intent(in) :: v2
      real (r8), dimension (1:3), intent(in) :: v3
      !Values and gradients at the nodes
      real (r8), intent(in) ::  f1
      real (r8), intent(in) :: f2
      real (r8), intent(in) :: f3
      real (r8), dimension (1:3), intent(in) ::  g1
      real (r8), dimension (1:3), intent(in) :: g2
      real (r8), dimension (1:3), intent(in) :: g3

      !Return interpolated value
      real (r8) :: fval

      !Directional derivative (scaled by distance)
      ! at u1, u2, or u3:  ds = (g,u1-v1)/u1n =
      ! -(g,v1)/u1n on side opposite v1, where g/
      ! u1n (plus an orthogonal component) is the
      ! projection of g onto the planar triangle
      real (r8):: ds

      !Directional derivatives (scaled by distance)
      ! at a vertex:  d1 = (g1,u1-v1) = (g1,u1)
      real (r8):: dv

      !Points on the boundary of the planar triangle
      ! and lying on the lines containing proj point
      ! and the vertices.  u1 is opposite v1, etc.
      ! uin=magnitudes of u1, u2, and u3)
      real (r8):: u1(3)
      real (r8):: u2(3)
      real (r8):: u3(3)
      real (r8):: u1n
      real (r8):: u2n
      real (r8):: u3n

      ! Central projections of u1, u2, and u3 onto
      !  the sphere and thus lying on an arc of the
      !  spherical triangle
      real (r8):: q1(3)
      real (r8):: q2(3)
      real (r8):: q3(3)

      ! f,g
      ! value and gradient at q1 q2, or q3 obtained
      ! by interpolation along one of the arcs of
      ! the spherical triangle
      real (r8):: f
      real (r8):: g(3)
      !Auxiliar variables
      real (r8):: c(1:3)
      real (r8):: sum
      real (r8):: s1
      real (r8):: s2
      real (r8):: s3
      real (r8):: val
      integer (i4):: i

      ! Compute weight functions c1, c2, and c3.
      ! c1 = 1 on the edge opposite v1 and c1 = 0 on the other edges.
      ! Similarly for c2 and c3.    c1+c2+c3 = 1.

      c(1) = b(2) * b(3)
      c(2) = b(3) * b(1)
      c(3) = b(1) * b(2)
      sum = c(1) + c(2) + c(3)
      if (sum <= eps/100) then
         ! p coincides with a vertex.
         fval = b(1) * f1 + b(2) * f2 + b(3) * f3
         return
      endif

      ! normalize c1, c2, and c3.
      c=c/sum

      !Compute u1,u2,u3
      s1 = b(2) + b(3)
      s2 = b(3) + b(1)
      s3 = b(1) + b(2)
      u1n = 0.
      u2n = 0.
      u3n = 0.
      do i = 1, 3
         u1 (i) = (b(2) * v2 (i) + b(3) * v3 (i) ) / s1
         u2 (i) = (b(3) * v3 (i) + b(1) * v1 (i) ) / s2
         u3 (i) = (b(1) * v1 (i) + b(2) * v2 (i) ) / s3
         u1n = u1n + u1 (i) * u1 (i)
         u2n = u2n + u2 (i) * u2 (i)
         u3n = u3n + u3 (i) * u3 (i)
      end do
      u1n = sqrt (u1n)
      u2n = sqrt (u2n)
      u3n = sqrt (u3n)

      ! Compute q1, q2, and q3.
      do i = 1, 3
         q1 (i) = u1 (i) / u1n
         q2 (i) = u2 (i) / u2n
         q3 (i) = u3 (i) / u3n
      end do

      val = 0.
      ! Contribution from side opposite v1:
      !  Compute value and gradient at q1 by interpolating
      !  between v2 and v3.
      call arcint (q1, v2, v3, f2, f3, g2, g3, f, g)
      dv = dot_product(g1, u1)
      ds = - dot_product(g, v1)/u1n
      val = val + c(1) * hval (b(1), f1, f, dv, ds)

      ! Contribution from side opposite v2:
      !   Compute value and gradient at q2 by interpolating
      !   between v3 and v1.
      call arcint (q2, v3, v1, f3, f1, g3, g1, f, g)
      dv = dot_product(g2, u2)
      ds = - dot_product(g, v2)/u2n
      val = val + c(2) * hval (b(2), f2, f, dv, ds)

      ! Contribution from side opposite v3:
      !   Compute interpolated value and gradient at q3
      !   by interpolating between v1 and v2.
      call arcint (q3, v1, v2, f1, f2, g1, g2, f, g)
      dv = dot_product(g3, u3)
      ds = -dot_product(g, v3)/u3n
      fval = val + c(3) * hval (b(3), f3, f, dv, ds)

      return
    end function fval

    subroutine arcint (p, p1, p2, f1, f2, g1, g2, f, g)
      !----------------------------------------------------------
      !  ARCINT from SSRFPACK by R. Renka
      !
      !   Given 3 points P, P1, and P2 lying on a common geodesic
      ! of the unit sphere with P between P1 and P2, along with
      ! data values and gradients at P1 and P2, this subroutine
      ! computes an interpolated value F and a gradient vector G
      ! AT P.  F and the tangential component of G are taken to be
      ! the value and derivative (with respect to arc-length) of
      ! a Hermite interpolatory spline defined by the end-
      ! point values and tangential gradient components.  The nor-
      ! mal component of G is obtained by linear interpolation of
      ! the normal components of the gradients at P1 and P2.
      !
      ! On input:
      !       P = Cartesian coordinates of a point lying on the
      !           arc defined by P1 and P2.  P(1)**2 + P(2)**2 +
      !           P(3)**2 = 1.
      !       P1,P2 = Coordinates of distinct points on the unit
      !               sphere defining an arc with length less than
      !               180 degrees.
      !       F1,F2 = Data values associated with P1 and P2,
      !               respectively.
      !       G1,G2 = Gradient vectors associated with P1 and P2.
      !               G1 and G2 are orthogonal to P1 and P2,
      !               respectively.
      !
      ! On output:
      !       F = Interpolated value at P.
      !       G = Interpolated gradient at P.
      !----------------------------------------------------------
      real (r8), dimension (1:3), intent(in) :: p
      real (r8), dimension (1:3), intent(in) :: p1
      real (r8), dimension (1:3), intent(in) :: p2
      real (r8), dimension (1:3), intent(in) :: g1
      real (r8), dimension (1:3), intent(in) :: g2
      real (r8), intent(in) :: f1
      real (r8), intent(in) :: f2
      real (r8), intent(out) :: f
      real (r8), intent(out) :: g(1:3)


      real (r8):: a
      real (r8):: al
      real (r8):: b1
      real (r8):: b2
      real (r8):: d1
      real (r8):: d2
      real (r8):: gn
      real (r8):: gt
      real (r8)::       &
           s
      real (r8):: tau1
      real (r8):: tau2
      real (r8):: unorm
      real (r8), dimension (1:3) :: un(1:3)

      ! Compute unit normal UN.
      un=cross_product(p1,p2)
      unorm = norm(un)
      if (unorm == 0.) then
         ! P1 X P2 = 0.  Print an error
         print*,"ARCINT ERROR: P1 x P2 = 0"
         print*, p1
         print*, p2
         stop
      end if
      un=un/unorm

      ! Compute tangential derivatives at the endpoints:
      !   tau1 = (g1,un x p1) = (g1,p2)/unorm and
      !   tau2 = (g2,un x p2) = -(g2,p1)/unorm.
      tau1 = dot_product(g1,p2)/unorm
      tau2 = -dot_product(g2, p1)/ unorm

      ! Compute arc-lengths a, al.
      a = arclen (p1, p2)
      if (a == 0.) then
         ! Print an error
         print*,"ARCINT ERROR: arclen = 0", a
         print*, p1
         print*, p2
         stop
      end if
      al = arclen (p1, p)

      ! Compute local coordinates, slope, and second differences.
      b2 = al / a
      b1 = 1. - b2
      s = (f2 - f1) / a
      d1 = s - tau1
      d2 = tau2 - s

      ! Hermite cubic interpolation
      f = f1 + al * (tau1 + b2 * (d1 + b1 * (d1 - d2) ) )

      !Compute gradient
      gt = tau1 + b2 * (d1 + d2 + 3. * b1 * (d1 - d2) )
      gn = b1 * (dot_product(un, g1))+ b2 * dot_product(un, g2)
      g=gt*cross_product(un,p)+gn*un

      return
    end subroutine arcint

    function hval (b, h1, h2, hp1, hp2)
      !------------------------------------------------------------
      ! HVAL from ssrfpack R. Renka
      !
      !   Given a line segment p1-p2 containing a point p, along
      ! with values and derivatives at the endpoints, this func-
      ! tion returns the value h(p), where h is the hermite inter-
      ! polatory spline defined by the endpoint data.
      !
      !       b = local coordinate of p with respect to p1-p2:
      !           p = b*p1 + (1-b)*p2, and thus b = d(p,p2)/
      !           d(p1,p2), where d(p1,p2) is the distance between
      !           p1 and p2.  b < 0 or b > 1 results in extrapola-
      !           tion.
      !
      !       h1,h2 = values interpolated at p1 and p2, respec-
      !               tively.
      !
      !       hp1,hp2 = products of d(p1,p2) with first order der-
      !                 ivatives at p1 and p2, respectively.  hp1
      !                 may, for example, be the scalar product of
      !                 p2-p1 with a gradient at p1.
      !-----------------------------------------------------------
      real (r8), intent(in) ::  b
      real (r8), intent(in) :: h1
      real (r8), intent(in) :: h2
      real (r8), intent(in) :: hp1
      real (r8), intent(in) :: hp2
      real (r8):: hval

      real (r8) :: b1
      real (r8) :: b2
      real (r8) :: d1
      real (r8) :: d2
      real (r8) :: s
      b1 = b
      b2 = 1. - b1

      ! compute slope s and second differences d1 and d2 scaled
      !   by the separation between p1 and p2.
      s = h2 - h1
      d1 = s - hp1
      d2 = hp2 - s

      ! hermite cubic interpolation:
      hval = h1 + b2 * (hp1 + b2 * (d1 + b1 * (d1 - d2) ) )

      return
    end function hval

  end function scinterpol_hermite_trv


  !=====================================================================================
  !   WACHSPRESS INTERPOLATION METHODS
  !=====================================================================================

  function scinterpol_wachspress(p, var,  mesh)
    !--------------------------------------------------------
    ! WACHSPRESS COORDINATE BASED INTERPOLATION
    ! See Gillette, Rand Bajaj 2011 adapted
    ! Uses the areas of triangles formed by point and
    !   voronoi vertices
    !--------------------------------------------------------
    !Interpolation point
    real (r8), intent(in) :: p(1:3)

    !Variable for interpolation procedure
    type(scalar_field), intent(in) :: var

    !Mesh structure
    type(grid_structure), intent(in) :: mesh

    !Interpoled value
    real (r8):: scinterpol_wachspress

    !Wach coords
    type(general_coords):: wachc

    !Index
    integer:: i

    !Interpolate acording to variable position
    select case(var%pos)
    case(0)  !Values on triangle vertices / voronoi centers
       wachc=wachspress_coords_trv(p, mesh)
       !Do linear interpolation
       !scinterpol_wachspress=scinterpol_linear_trv (p, var, mesh)
       !return
    case(1)  !Values on triangles circumcenters / voronoi vertices
       wachc=wachspress_coords_hxv(p, mesh)
    case(3) !Values on hexagon/voronoi edge midpoints
       wachc=wachspress_coords_hxed(p, mesh)
    case default
       print*, "Error on scinterpol_wachspress: Variable position problem", var%pos
       stop
    end select

    scinterpol_wachspress=0._r8
    do i=1, wachc%n
       scinterpol_wachspress=scinterpol_wachspress+&
            wachc%w(i)*var%f(wachc%v(i))
    end do

    return
  end function scinterpol_wachspress

  function vecinterpol_wachspress(p, var,  mesh)
    !--------------------------------------------------------
    ! WACHSPRESS COORDINATE BASED VECTOR INTERPOLATION
    ! See Gillette, Rand Bajaj 2011 adapted for vector
    ! Uses the areas of triangles formed by point and
    !   voronoi vertices
    !--------------------------------------------------------
    !Interpolation point
    real (r8), intent(in) :: p(1:3)

    !Vector field
    type(vector_field_cart), intent(in) :: var

    !Mesh structure
    type(grid_structure), intent(in) :: mesh

    !Interpoled value
    real (r8):: vecinterpol_wachspress(1:3)

    !Wach coords
    type(general_coords):: wachc

    !Index
    integer:: i

    !Interpolate acording to variable position
    select case(var%pos)
    case(0)  !Values on triangle vertices / voronoi centers
       !Do linear interpolation
       vecinterpol_wachspress=vecinterpol_linear_trv(p, var, mesh)
       return
    case(4)  !Values on Voronoi centroids
       !Do linear interpolation
       vecinterpol_wachspress=vecinterpol_linear_hxb(p, var, mesh)
       return
    case(1)  !Values on triangles circumcenters / voronoi vertices
       wachc=wachspress_coords_hxv(p, mesh)
    case(5)  !Values on triangles barycenters
       wachc=wachspress_coords_hxb(p, mesh)
    case(3) !Values on hexagon/voronoi edge midpoints
       wachc=wachspress_coords_hxed(p, mesh)
    case default
       print*, "Error on vecinterpol_wachspress: Variable position problem", var%pos
       stop
    end select

    vecinterpol_wachspress(1:3)=(/0._r8, 0._r8, 0._r8/)
    do i=1, wachc%n
       vecinterpol_wachspress=vecinterpol_wachspress+&
            wachc%w(i)*var%p(wachc%v(i))%v
    end do

    vecinterpol_wachspress=proj_vec_sphere(vecinterpol_wachspress, p)

    return
  end function vecinterpol_wachspress

  function wachspress_coords_trv(p, mesh, triang)
    !--------------------------------------------------------
    ! WACHSPRESS COORDINATES FOR NODES ON TR vertices
    !   that is, hx nodes
    !
    ! Equivalent to barycentric coordinates using sphecrical areas
    !--------------------------------------------------------
    !Interpolation point
    real (r8), intent(in) :: p(1:3)

    !Mesh structure
    type(grid_structure), intent(in) :: mesh

    !Triangle containing p - optional
    integer (i4), optional:: triang
    integer (i4):: tri

    !Wachspress coordinates
    type(general_coords):: wachspress_coords_trv

    !Aux variables
    real (r8), dimension(1:3) :: p1
    real (r8), dimension(1:3) :: p2
    real (r8), dimension(1:3) :: p3

    !Indexes
    integer (i4):: i

    !Initialize returning variable for
    ! values on tr vertices / hx nodes
    wachspress_coords_trv%pos=0

    if(present(triang))then
       tri=triang
    else
       !First get the triangle
       tri=gettr(p, mesh)
    end if


    wachspress_coords_trv%n=3
    allocate(wachspress_coords_trv%v(1:3))
    allocate(wachspress_coords_trv%w(1:3))

    !Vertices for coordinates
    wachspress_coords_trv%v(1:3)=mesh%tr(tri)%v(1:3)

    p1=p
    do i=1, 3
       !Set vertice points
       p2=mesh%v(mesh%tr(tri)%v(i))%p
       p3=mesh%v(mesh%tr(tri)%v(modint(i+1,3)))%p
       !Calculate Bi
       wachspress_coords_trv%w(modint(i-1,3))= &
            sphtriarea(p1, p2, p3)/mesh%tr(tri)%areag
    end do

    return
  end function wachspress_coords_trv

  function wachspress_coords_hxv(p, mesh, nearnode)
    !--------------------------------------------------------
    ! WACHSPRESS COORDINATES FOR NODES ON HX Vertices
    !   that is, triangle circumcenters
    !
    ! See Gillette, Rand, Bajaj 2011
    ! Uses the areas of triangles formed by point and
    !   voronoi cell vertices
    !--------------------------------------------------------
    !Interpolation point
    real (r8), intent(in) :: p(1:3)

    !Mesh structure
    type(grid_structure), intent(in) :: mesh

    !Nearest node - optional
    integer (i4), optional:: nearnode

    !Wachspress coordinates
    type(general_coords):: wachspress_coords_hxv

    !Aux variables
    real (r8), dimension(1:3) :: p1
    real (r8), dimension(1:3) :: p2
    real (r8), dimension(1:3) :: p3

    !Weights
    real (r8), allocatable :: w(:)
    real (r8), allocatable :: A(:)
    real (r8), allocatable :: B(:)

    !Sum of weights
    real (r8):: sumw

    !Indexes to edges related to each weight
    integer (i4), allocatable :: wtrs(:)

    !Number of weights
    integer (i4):: wn

    !Indexes
    integer (i4):: i
    integer (i4):: j
    integer (i4):: node

    !Initialize returning variable for
    ! values on hx vertices / tr cc
    wachspress_coords_hxv%pos=1

    if(present(nearnode))then
       node=nearnode
    else
       !First get the nearest node (voronoi cell)
       node=getnearnode(p, mesh)
    end if

    !Calculate wachpress coordinates relative to hexagon verts / trcc
    wn=mesh%v(node)%nnb
    wachspress_coords_hxv%n=wn
    allocate(w(1:wn))
    allocate(wachspress_coords_hxv%v(1:wn))
    allocate(wachspress_coords_hxv%w(1:wn))
    allocate(A(1:wn))
    allocate(B(1:wn))
    allocate(wtrs(1:wn))
    wachspress_coords_hxv%v(1:wn)=mesh%v(node)%tr(1:wn)
    wtrs=mesh%v(node)%tr
    do i=1, wn
       !Set vertice points
       p1=mesh%tr(wtrs(modint(i-1,wn)))%c%p
       p2=mesh%tr(wtrs((i)))%c%p
       p3=mesh%tr(wtrs(modint(i+1,wn)))%c%p
       !Calculate Bi
       B(i)=sphtriarea(p1, p2, p3)/mesh%hx(node)%areag

       !Calculate Ai s
       !These are areas of triangles given cc-wisely
       A(i)=sphtriarea(p2, p3, p)/mesh%hx(node)%areag
    end do

    !Calculate weights
    do i=1, wn
       !Method 1 - Problems on edges
       !A2=A(i)*A(modint(i-1,wn))
       !w(i)=B(i)/A2

       !Method 2 - No problems on edges
       w(i)=B(i)
       do j=1, wn
          if(j/=i .and. j/=modint(i-1,wn))then
             w(i)=w(i)*A(j)
          end if
       end do
    end do

    !Calculate and save the coordinates
    sumw=sum(w(1:wn))
    if(sumw<eps)then
       print*, "WACHSPRESS_COORDS warning: possible numerical precision problem"
    end if
    do i=1, wn
       wachspress_coords_hxv%w(i)=w(i)/sumw
    end do

    return
  end function wachspress_coords_hxv

  function wachspress_coords_hxb(p, mesh, nearnode)
    !--------------------------------------------------------
    ! WACHSPRESS COORDINATES FOR NODES ON HX Vertices
    !   which are the triangle barycenters
    !
    ! See Gillette, Rand, Bajaj 2011
    ! Uses the areas of triangles formed by point and
    !   voronoi cell vertices
    !--------------------------------------------------------
    !Interpolation point
    real (r8), intent(in) :: p(1:3)

    !Mesh structure
    type(grid_structure), intent(in) :: mesh

    !Nearest node - optional
    integer (i4), optional:: nearnode

    !Wachspress coordinates
    type(general_coords):: wachspress_coords_hxb

    !Aux variables
    real (r8), dimension(1:3) :: p1
    real (r8), dimension(1:3) :: p2
    real (r8), dimension(1:3) :: p3

    !Weights
    real (r8), allocatable :: w(:)
    real (r8), allocatable :: A(:)
    real (r8), allocatable :: B(:)

    !Sum of weights
    real (r8):: sumw

    !Indexes to edges related to each weight
    integer (i4), allocatable :: wtrs(:)

    !Number of weights
    integer (i4):: wn

    !Indexes
    integer (i4):: i
    integer (i4):: j
    integer (i4):: node

    !Initialize returning variable for
    ! values on hx vertices / tr cc
    wachspress_coords_hxb%pos=5

    if(present(nearnode))then
       node=nearnode
    else
       !First get the nearest node (voronoi cell)
       node=getnearnode(p, mesh)
    end if

    !Calculate wachpress coordinates relative to hexagon verts / trcc
    wn=mesh%v(node)%nnb
    wachspress_coords_hxb%n=wn
    allocate(w(1:wn))
    allocate(wachspress_coords_hxb%v(1:wn))
    allocate(wachspress_coords_hxb%w(1:wn))
    allocate(A(1:wn))
    allocate(B(1:wn))
    allocate(wtrs(1:wn))
    wachspress_coords_hxb%v(1:wn)=mesh%v(node)%tr(1:wn)
    wtrs=mesh%v(node)%tr
    do i=1, wn
       !Set vertice points
       p1=mesh%tr(wtrs(modint(i-1,wn)))%b%p
       p2=mesh%tr(wtrs((i)))%b%p
       p3=mesh%tr(wtrs(modint(i+1,wn)))%b%p
       !Calculate Bi
       B(i)=sphtriarea(p1, p2, p3)/mesh%hx(node)%areag

       !Calculate Ai s
       !These are areas of triangles given cc-wisely
       A(i)=sphtriarea(p2, p3, p)/mesh%hx(node)%areag
    end do

    !Calculate weights
    do i=1, wn
       !Method 1 - Problems on edges
       !A2=A(i)*A(modint(i-1,wn))
       !w(i)=B(i)/A2

       !Method 2 - No problems on edges
       w(i)=B(i)
       do j=1, wn
          if(j/=i .and. j/=modint(i-1,wn))then
             w(i)=w(i)*A(j)
          end if
       end do
    end do

    !Calculate and save the coordinates
    sumw=sum(w(1:wn))
    if(sumw<eps)then
       print*, "WACHSPRESS_COORDS warning: possible numerical precision problem"
    end if
    do i=1, wn
       wachspress_coords_hxb%w(i)=w(i)/sumw
    end do

    return
  end function wachspress_coords_hxb

  function wachspress_coords_hxed(p, mesh)
    !--------------------------------------------------------
    ! WACHSPRESS COORDINATES FOR NODES ON HX EDGE MIDPOINTS
    !
    ! See Gillette, Rand Bajaj 2011
    ! Uses the areas of triangles formed by point and
    !   voronoi vertices
    !--------------------------------------------------------
    !Interpolation point
    real (r8), intent(in) :: p(1:3)

    !Mesh structure
    type(grid_structure), intent(in) :: mesh

    !Wachspress coordinates
    type(general_coords):: wachspress_coords_hxed

    !Triangle edge index
    integer (i4), dimension(1:3) :: eds

    !Points
    real (r8), dimension(1:3) :: p1
    real (r8), dimension(1:3) :: p2
    real (r8), dimension(1:3) :: p3

    !Logical to test if inside minor triangle
    logical:: intr

    !Weights
    real (r8), allocatable :: w(:)
    real (r8), allocatable :: A(:)
    real (r8), allocatable :: B(:)

    !Sum of weights
    real (r8):: sumw

    !Indexes to edges related to each weight
    integer (i4), allocatable :: weds(:)

    !Number of weights
    integer (i4):: wn

    !Indexes
    integer (i4):: i
    integer (i4):: j
    integer (i4):: node
    integer (i4):: tr

    !Initialize returning variable for
    ! values on hx midpoints
    wachspress_coords_hxed%pos=3

    !First get the nearest node (voronoi cell)
    node=getnearnode(p, mesh)

    !Find which triangle circumcenter is closest to the point
    tr=getneartrc(p, mesh)

    !Get triangle edges
    do i=1, 3
       eds(i)=mesh%tr(tr)%ed(i)
    end do

    !Check if it is inside the triangle generated by
    !   the 3 edges near a tr circumcenter, that is
    !   check if it is outside the hexagon formed by
    !   edge centers
    intr=insmallvoredtr(p, eds, mesh)

    !Calculate the wachspress coordinates
    if(intr)then
       !Calculate wachpress coordinates relative to minor triangle
       !In this case it is exactly the barycentric coordinates
       ! for the spherical triangle of this minor tr
       wn=3
       wachspress_coords_hxed%n=3
       allocate(w(1:3))
       allocate(wachspress_coords_hxed%v(1:3))
       allocate(wachspress_coords_hxed%w(1:3))
       allocate(A(1:3))
       allocate(B(1:1))
       wachspress_coords_hxed%v(1:3)=eds(1:3)
       p1=mesh%edhx(eds(1))%c%p
       p2=mesh%edhx(eds(2))%c%p
       p3=mesh%edhx(eds(3))%c%p
       B(1)=sphtriarea(p1, p2, p3)

       !These A(:) represent the areas of
       ! the oposite triangle to vertice
       A(1)=sphtriarea(p, p2, p3)
       A(2)=sphtriarea(p1, p, p3)
       A(3)=sphtriarea(p1, p2, p)

       !Set weigths
       do i=1, 3
          w(i)=A(i)/B(1)
       end do

    else
       !Calculate wachpress coordinates relative to larger cell (hexagon)
       wn=mesh%v(node)%nnb
       wachspress_coords_hxed%n=wn
       allocate(w(1:wn))
       allocate(wachspress_coords_hxed%v(1:wn))
       allocate(wachspress_coords_hxed%w(1:wn))
       allocate(A(1:wn))
       allocate(B(1:wn))
       allocate(weds(1:wn))
       wachspress_coords_hxed%v(1:wn)=mesh%v(node)%ed(1:wn)
       weds=mesh%v(node)%ed
       do i=1, wn
          !Set edge points
          p1=mesh%edhx(weds(modint(i-1,wn)))%c%p
          p2=mesh%edhx(weds(i))%c%p
          p3=mesh%edhx(weds(modint(i+1,wn)))%c%p
          !Calculate Bi
          B(i)=sphtriarea(p1, p2, p3)/mesh%hx(node)%areag

          !Calculate and Ai s
          !These are areas of triangles given cc-wisely
          A(i)=sphtriarea(p2, p3, p)/mesh%hx(node)%areag
       end do
       !Calculate weights
       do i=1, wn
          !Method 1 - Problems on edges
          !A2=A(i)*A(modint(i-1,wn))
          !w(i)=B(i)/A2

          !Method 2 - No problems on edges
          w(i)=B(i)
          do j=1, wn
             if(j/=i .and. j/=modint(i-1,wn))then
                w(i)=w(i)*A(j)
             end if
          end do
       end do
    end if

    !Calculate and save coordinates
    sumw=sum(w(1:wn))
    if(sumw<eps)then
       print*, "wachspress_coords_edhx warning: possible numerical precision problem"
    end if
    do i=1, wn
       wachspress_coords_hxed%w(i)=w(i)/sumw
    end do

    return
  end function wachspress_coords_hxed


  function vecrecon_wach(p, var, mesh)
    !--------------------------------------------------------
    ! vecrecon_wach
    ! New vector reconst.
    ! Uses RT0 with 2 components on vertices and wachspress
    !--------------------------------------------------------
    !Reconstruction point
    real (r8), intent(in) :: p(1:3)

    !Variable with normal vector values
    type(scalar_field), intent(in) :: var

    !Mesh structure
    type(grid_structure), intent(in) :: mesh

    !Interpoled value
    real (r8):: vecrecon_wach(1:3)

    !Aux var
    real (r8):: u1
    real (r8):: u2
    real (r8), dimension(1:3) :: n1
    real (r8), dimension(1:3) :: n2
    real (r8), dimension(1:3) :: v
    real (r8), dimension(1:3) :: q

    !Wach coords
    type(general_coords):: wachc

    !Indexes and aux vars
    integer (i4):: i
    integer (i4):: j
    integer (i4):: k1
    integer (i4):: k2
    integer (i4):: tr

    !Do reconstruction assuming point in hexagon
    !Get nearest node to point
    i=getnearnode(p, mesh)

    !Calculate Wachspress coordinates of point
    wachc=wachspress_coords_hxv(p, mesh)

    !Do reconstruction
    vecrecon_wach=(/0._r8,0._r8,0._r8/)

    do j=1,mesh%v(i)%nnb
       tr=mesh%v(i)%tr(j)
       k1=mesh%v(i)%ed(j)
       k2=mesh%v(i)%ed(modint(j+1,mesh%v(i)%nnb) )
       n1=mesh%edhx(k1)%nr
       n2=mesh%edhx(k2)%nr
       q=mesh%tr(tr)%c%p

       !Normal components
       u1=var%f(k1)
       u2=var%f(k2)
       !Solve the 3x3 system
       v=solve3x3(n1, n2, q, (/u1, u2, 0._r8/))
       vecrecon_wach=vecrecon_wach+wachc%w(j)*v
    end do

    return
  end function vecrecon_wach

  function vecrecon_howach(p, var, mesh)
    !--------------------------------------------------------
    ! vecrecon_howach
    ! New vector reconst higher order
    ! Uses 3 nodes RT0 reconstruction for vertices
    ! and then wachspress
    !--------------------------------------------------------
    !Reconstruction point
    real (r8), intent(in) :: p(1:3)

    !Variable with normal vector values
    type(scalar_field), intent(in) :: var

    !Mesh structure
    type(grid_structure), intent(in) :: mesh

    !Interpoled value
    real (r8):: vecrecon_howach(1:3)

    !Aux var
    real (r8):: ncor
    real (r8), dimension(1:3) :: tg
    real (r8), dimension(1:3) :: vtg
    real (r8), dimension(1:3) :: v
    real (r8), dimension(1:3) :: r

    !Wach coords
    type(general_coords):: wachc

    !Indexes and aux vars
    integer (i4):: i
    integer (i4):: j
    integer (i4):: l
    integer (i4):: kt
    integer (i4):: ed

    !Do reconstruction assuming point in hexagon
    !Get nearest node to point
    i=getnearnode(p, mesh)

    !Calculate Wachspress coordinates of point
    wachc=wachspress_coords_hxv(p, mesh)

    !Do reconstruction
    vecrecon_howach=(/0._r8,0._r8,0._r8/)
    ! Perot vorticity method
    do j=1,mesh%v(i)%nnb
       kt=mesh%v(i)%tr(j)
       v=(/0._r8, 0._r8, 0._r8 /)
       !print*, "Node, neighb, triangle:", i, j, kt
       do l=1, 3
          ed=mesh%tr(kt)%ed(l)
          !r=mesh%v(mesh%tr(kt)%v(l))%p-mesh%edhx(ed)%c%p
          r=mesh%ed(ed)%c%p-mesh%tr(kt)%c%p
          !Counter-clock-wise tangent direction vector on triangle edge midpoint
          tg=mesh%tr(kt)%tg(l)*mesh%ed(ed)%tg
          !tg=mesh%ed(ed)%tg
          !Check of vector components points in the cc-wise direction
          vtg=mesh%edhx(ed)%nr
          ncor=dot_product(vtg, tg)
          !print*, "   ed:" , ed, "Tangent correction:", ncor
          v=v + ncor*(var%f(ed))*r*mesh%ed(ed)%leng
       end do
       v=cross_product(mesh%tr(kt)%c%p, v)
       v=v/mesh%tr(kt)%areag
       vecrecon_howach=vecrecon_howach+wachc%w(j)*v
       !print*, vecrecon_howach
    end do

    return
  end function vecrecon_howach

  !========================================================================
  !    P1NC INTERPOLATION
  !========================================================================

  function scinterpol_p1nc (p, var, mesh)
    !-----------------------------------------------------------
    !  scinterpol_p1nc
    !
    !   Given the triangulation of a set of nodes on the unit
    ! sphere (mesh), along with data values at the edges of tr or hx
    ! (in var), computes the value at a point (p)
    ! of a non conforming function which interpolates the data values.
    !
    ! On input:
    !       p         =  point in cart coords
    !       mesh      =  mesh structure
    !       var       =  interpolable_variable structure
    !                      with var%pos=2 or 3
    ! Returns the interpolated value
    !-------------------------------------------------------
    !Interpolation point
    real (r8), intent(in) :: p(1:3)

    !Variable for interpolation procedure
    type(scalar_field), intent(in) :: var

    !Mesh in which the variable belongs
    type(grid_structure), intent(in) :: mesh

    !Value at interpolated point
    real (r8):: f
    real (r8):: scinterpol_p1nc

    !loop index, indexes for nodes, triangle index
    integer(i4):: i
    integer(i4):: j
    integer(i4):: k
    integer(i4):: ed(1:3)

    !Cartesian and barycentric coordinates of the points
    real (r8), dimension(1:3) :: b
    real (r8), dimension(1:3) :: p1
    real (r8), dimension(1:3) :: p2
    real (r8), dimension(1:3) :: p3

    ! Locate p with respect to the triangulation.
    k=gettr(p, mesh)

    !Edges associated with bary coords
    !  The edge is oposite to corresponding vertice on K
    do i=1,3
       j=i !modint(i+1, 3)
       ed(i)=mesh%tr(k)%ed(j)
    end do

    !Check var position
    select case (var%pos)
    case(2) !Values given on triangle k edge midpoints
       !Use edge as triangle
       p1=mesh%ed(ed(1))%c%p
       p2=mesh%ed(ed(2))%c%p
       p3=mesh%ed(ed(3))%c%p

       !Use vertices as triangle
       !p1=mesh%v(mesh%tr(k)%v(1))%p
       !p2=mesh%v(mesh%tr(k)%v(2))%p
       !p3=mesh%v(mesh%tr(k)%v(3))%p

       !Get barycentric coords based on the vertices
       ! of the triangle
       !b=bar_coord_tr(p, k, mesh)

    case(3) !Triangle given on hx edge midpoints
       !Calculate vertices of the triangle that
       !  conects the hx edge midpoints
       p1=mesh%edhx(ed(1))%c%p
       p2=mesh%edhx(ed(2))%c%p
       p3=mesh%edhx(ed(3))%c%p

       !p1=gcircintersec(p, mesh%edhx(ed(2))%tg, mesh%edhx(ed(3))%tg)
       !p2=gcircintersec(p, mesh%edhx(ed(3))%tg, mesh%edhx(ed(1))%tg)
       !p3=gcircintersec(p, mesh%edhx(ed(1))%tg, mesh%edhx(ed(2))%tg)
       !print*, "Triangle ", k
       !print*, "points "
       !print*, p1
       !print*, p2
       !print*, p3

    case default
       print*, "Error on scinterpol_p1nc: Variable position problem", var%pos
       stop
    end select

    !Calculate the barycentric coord
    b=bar_coord(p, p1, p2, p3)
    !b=bar_coord_plan(p, p1, p2, p3)

    !Calculate interpolation
    f = 0.0
    do i=1,3
       f=f+ b(i)*var%f(ed(i))
       !f = f + (1.0-2.0*b(i))*var%f(ed(i))
    end do

    scinterpol_p1nc=f

    return
  end function scinterpol_p1nc


  function vecinterpol_p1nc (p, var, mesh)
    !-----------------------------------------------------------
    !  vecinterpol_p1nc
    !
    !   Given the triangulation of a set of nodes on the unit
    ! sphere (mesh), along with vector field at the edges of tr or hx
    ! (in var), computes the value at a point (p)
    ! of a non conforming function which interpolates the data values.
    !
    ! On input:
    !       p         =  point in cart coords
    !       mesh      =  mesh structure
    !       var       =  interpolable_variable structure
    !                      with var%pos=2 or 3
    ! Returns the interpolated value
    !-------------------------------------------------------
    !Interpolation point
    real (r8), intent(in) :: p(1:3)

    !Vector field
    type(vector_field_cart), intent(in) :: var

    !Mesh in which the variable belongs
    type(grid_structure), intent(in) :: mesh

    !Vector at interpolated point
    real (r8), dimension(1:3) :: vec
    real (r8), dimension(1:3) :: vecinterpol_p1nc

    !loop index, indexes for nodes, triangle index
    integer(i4):: i
    integer(i4):: j
    integer(i4):: k
    integer(i4):: ed(1:3)

    !Cartesian and barycentric coordinates of the points
    real (r8), dimension(1:3) :: b
    real (r8), dimension(1:3) :: p1
    real (r8), dimension(1:3) :: p2
    real (r8), dimension(1:3) :: p3

    ! Locate p with respect to the triangulation.
    k=gettr(p, mesh)

    !Edges associated with bary coords
    !  The edge is oposite to corresponding vertice on K
    do i=1,3
       j=i !modint(i+1, 3)
       ed(i)=mesh%tr(k)%ed(j)
    end do

    !Check var position
    select case (var%pos)
    case(2) !Values given on triangle k edge midpoints
       !Use edge as triangle
       p1=mesh%ed(ed(1))%c%p
       p2=mesh%ed(ed(2))%c%p
       p3=mesh%ed(ed(3))%c%p

       !Use vertices as triangle
       !p1=mesh%v(mesh%tr(k)%v(1))%p
       !p2=mesh%v(mesh%tr(k)%v(2))%p
       !p3=mesh%v(mesh%tr(k)%v(3))%p

       !Get barycentric coords based on the vertices
       ! of the triangle
       !b=bar_coord_tr(p, k, mesh)

    case(3) !Triangle given on hx edge midpoints
       !Calculate vertices of the triangle that
       !  conects the hx edge midpoints
       p1=mesh%edhx(ed(1))%c%p
       p2=mesh%edhx(ed(2))%c%p
       p3=mesh%edhx(ed(3))%c%p

       !p1=gcircintersec(p, mesh%edhx(ed(2))%tg, mesh%edhx(ed(3))%tg)
       !p2=gcircintersec(p, mesh%edhx(ed(3))%tg, mesh%edhx(ed(1))%tg)
       !p3=gcircintersec(p, mesh%edhx(ed(1))%tg, mesh%edhx(ed(2))%tg)
       !print*, "Triangle ", k
       !print*, "points "
       !print*, p1
       !print*, p2
       !print*, p3

    case default
       print*, "Error on vecinterpol_p1nc: Variable position problem", var%pos
       stop
    end select

    !Calculate the barycentric coord
    b=bar_coord(p, p1, p2, p3)
    !b=bar_coord_plan(p, p1, p2, p3, mesh)

    !Calculate interpolation
    vec = (/0.0, 0.0, 0.0 /)
    do i=1,3
       vec=vec+ b(i)*var%p(ed(i))%v
       !f = f + (1.0-2.0*b(i))*var%f(ed(i))
    end do

    vecinterpol_p1nc=proj_vec_sphere(vec, p)

    return
  end function vecinterpol_p1nc


  !=====================================================================================
  !    Vector reconstructions
  !=====================================================================================

  function vecrecon_perot(p, var, stencil, mesh, pindex)
    !--------------------------------------------------------
    ! vecrecon_perot
    ! Perot vector reconstruction - See Perot (2000)
    !--------------------------------------------------------
    !Reconstruction point
    real (r8), intent(in) :: p(1:3)

    !Variable with normal vector values
    type(scalar_field), intent(in) :: var

    !Stencil used
    character(len=4), intent(in) :: stencil

    !Mesh structure
    type(grid_structure), intent(in) :: mesh

    !Index of the nearest node to p (optional)
    !  or the triangle that p belongs
    integer(i4), optional:: pindex

    !Interpoled value
    real (r8):: vecrecon_perot(1:3)

    !Aux var
    real (r8):: ncor, edlen
    real (r8):: ppj(1:3) !Projeted point
    real (r8):: p1(1:3) !Edge point
    real (r8):: p2(1:3) !Edge point

    !Indexes and aux vars
    integer (i4):: i
    integer (i4):: j
    integer (i4):: k
    integer (i4):: ktr

    i=0

    !Get nearest node
    if(present(pindex) .and. (trim(stencil)=="hx" .or. trim(stencil)=="pj"))then
       i=pindex
       ktr=-1
    end if

    !Get triangle contaning node
    if(present(pindex) .and. trim(stencil)=="tr" )then
       i=-1
       ktr=pindex
    end if

    if(i==0 .or. ktr==0)then
       !Get triangle containing point
       ktr=gettr(p,mesh)
       !Get nearest node to point
       i=getnearnode(p, mesh, ktr)
    end if

    !Do reconstruction
    vecrecon_perot=(/0._r8,0._r8,0._r8/)

    if(var%pos==3)then  !Normal components given on hx edge midpoints
       if(trim(stencil)=="hx")then !Reconstruct based on hexagons
          !See Perot 2007 - Mimetic reconstruction of vectors - eq (2.2)
          !or Perot 2000 - eq (95)
          do j=1,mesh%v(i)%nnb
             k=mesh%v(i)%ed(j)
             ncor=real(mesh%hx(i)%nr(j), r8)
             vecrecon_perot=vecrecon_perot + &
                  ncor*(var%f(k))*(mesh%edhx(k)%c%p-p)*mesh%edhx(k)%leng
          end do
          vecrecon_perot=vecrecon_perot/mesh%hx(i)%areag

       elseif(trim(stencil)=="pj")then
          !Reconstruct based on hexagons
          ! but first project elements to tangent plane
          do j=1,mesh%v(i)%nnb
             k=mesh%v(i)%ed(j)
             ncor=real(mesh%hx(i)%nr(j), r8)
             ppj=mesh%edhx(k)%c%p/dot_product(mesh%edhx(k)%c%p, p)
             p1=mesh%tr(mesh%edhx(k)%v(1))%c%p
             p1=p1/dot_product(p1,p)
             p2=mesh%tr(mesh%edhx(k)%v(2))%c%p
             p2=p2/dot_product(p2,p)
             edlen=norm(p1-p2)
             vecrecon_perot=vecrecon_perot + &
                  ncor*(var%f(k))*(ppj-p)*edlen
          end do
          vecrecon_perot=vecrecon_perot/planarpjhexagarea(i, p, mesh)
          !print*, dot_product(vecrecon_perot, p)
       else  !Stencil='tr'
          !Reconstruct based on triangle around vertices of hexagons
          ! See Perot 2000 eq 89
          !print*, "Triangle:", ktr
          do j=1,3
             k=mesh%tr(ktr)%ed(j)
             !Because we are mixing HC and TR normals and tangents,
             !  the corection needs some adjustment
             ncor=real(mesh%tr(ktr)%tg(j), r8)*dot_product(mesh%ed(k)%tg, mesh%edhx(k)%nr)
             if(ncor > 0 ) then
                ncor=1._r8
             else
                ncor=-1._r8
             end if
             vecrecon_perot=vecrecon_perot + &
                  ncor*(var%f(k))*(mesh%ed(k)%c%p-p)*mesh%ed(k)%leng
          end do
          vecrecon_perot=cross_product(p,vecrecon_perot)/mesh%tr(ktr)%areag
          !vecrecon_perot=vecrecon_perot/mesh%tr(ktr)%areag
       end if
    else !Data given on triangle edges
       if(trim(stencil)=="tr")then !Reconstruct based on triangles to circumcenters
          !See Perot 2007 - Mimetic reconstruction of vectors - eq (2.2)
          !or Perot 2000 - eq (95)
          do j=1,3
             k=mesh%tr(ktr)%ed(j)
             ncor=real(mesh%tr(ktr)%nr(j), r8)
             vecrecon_perot=vecrecon_perot + &
                  ncor*(var%f(k))*(mesh%ed(k)%c%p-p)*mesh%ed(k)%leng
          end do
          vecrecon_perot=vecrecon_perot/mesh%tr(ktr)%areag
       else  !Stencil='hx'
          ! See Perot 2000 eq 89
          !print*, stencil, " node: ", i
          !print*, p-mesh%v(i)%p
          do j=1,mesh%v(i)%nnb
             k=mesh%v(i)%ed(j)

             !Because we are mixing TC and  normals and tangents,
             !  the corection needs some adjustment
             ncor=real(mesh%hx(i)%tg(j), r8)*dot_product(mesh%edhx(k)%tg, mesh%ed(k)%nr)
             if(ncor > 0 ) then
                ncor=1._r8
             else
                ncor=-1._r8
             end if
             !print*, "Ed:", k, ncor, var%f(k)
             vecrecon_perot=vecrecon_perot + &
                  ncor*(var%f(k))*(mesh%edhx(k)%c%p-p)*mesh%edhx(k)%leng
          end do
          vecrecon_perot=cross_product(p,vecrecon_perot)/mesh%hx(i)%areag
       end if
    end if
    return
  end function vecrecon_perot

  function vecrecon_pered(ed, var, mesh)
    !--------------------------------------------------------
    ! vecrecon_pered
    ! Vector reconstruction with Perots method for edges - tangent component
    !--------------------------------------------------------
    !Reconstruction point
    integer(i4), intent(in) :: ed

    !Variable with normal vector values
    type(scalar_field), intent(in) :: var

    !Mesh structure
    type(grid_structure) :: mesh

    !Reconstructed vector
    real (r8), dimension(1:3):: vecrecon_pered

    !Just tangent component of reconstructed vector
    real (r8):: utg
    real (r8):: ncor
    real (r8):: w
    !real (r8):: die
    !real (r8):: nrtg

    !Cells adjacent to edge ed
    integer (i4):: cell(1:2)

    !Indexes and aux vars
    integer (i4):: i
    integer (i4):: j
    integer (i4):: k
    integer (i4):: ed_celli
    !integer (i4):: signcor

    !Check for argument entry error
    if(ed==0)then
       print*, "vecrecon_trsk error: Edge index must be informed"
       stop
    end if

    !Initialize reconstruction vectors
    vecrecon_pered=0._r8
    utg=0._r8

    !Get cell surrounding edge
    cell(1:2)=mesh%edhx(ed)%sh(1:2)

    !Debug
    !if(ed/=93)then
    !  return
    !end if

    !Works for if(var%pos==3)then  !Normal components given on hx edge midpoints
    do i=1,2
       ed_celli=0
       do j=1, mesh%v(cell(i))%nnb
          !Find the index of edge ed within cell(i)
          if(ed==mesh%v(cell(i))%ed(j))then
             ed_celli=j
             exit
          end if
       end do
       if(ed_celli==0)then !check for errors
          print*, "vecrecon_pered error: edge not in cell"
          stop
       end if
       !For all edges of cell, add contribution to reconstruction
       do j=1, mesh%v(cell(i))%nnb
          !Edge's global index
          k=mesh%v(cell(i))%ed(j)
          ncor=real(mesh%hx(cell(i))%nr(j), r8)
          !Use the dual edge
          !ncor=dsign( 1._r8, ncor* dot_product(mesh%edhx(k)%nr, mesh%ed(k)%tg))
          w=ncor*dot_product(mesh%edhx(k)%c%p-mesh%v(cell(i))%p, mesh%edhx(ed)%tg)
          !w=ncor*dot_product(mesh%edhx(k)%c%p-mesh%v(cell(i))%p, mesh%edhx(ed)%tg)
          w=w*mesh%edhx(k)%leng/(2.*mesh%hx(cell(i))%areag)
          utg=utg+var%f(k)*w
       end do
    end do
    !utg=utg/mesh%ed(ed)%leng

    !Get reconstructed vector including normal component
    vecrecon_pered=utg*mesh%edhx(ed)%tg+var%f(ed)*mesh%edhx(ed)%nr

    return
  end function vecrecon_pered

  function pered_order_index(ed, mesh)
    !--------------------------------------------------------
    ! pered_order_index
    ! Calculated the consistency conditions for Perot edge method
    !  and returns an index
    !
    !--------------------------------------------------------
    !Reconstruction point
    integer(i4), intent(in) :: ed

    !Mesh structure
    type(grid_structure) :: mesh

    !Reconstructed vector
    real (r8) :: pered_order_index

    !Reconstructed tangent vector
    real (r8) :: vtg(1:3)

    !Cells adjacent to edge ed
    integer (i4):: cell(1:2)

    !Just tangent component of reconstructed vector
    !real (r8):: utg
    real (r8):: ncor
    real (r8):: w
    !real (r8):: die
    !real (r8):: nrtg

    !Indexes and aux vars
    integer (i4):: i
    integer (i4):: j
    integer (i4):: k
    integer (i4):: ed_celli
    !integer (i4):: signcor

    !Check for argument entry error
    if(ed==0)then
       print*, "vecrecon_trsk error: Edge index must be informed"
       stop
    end if

    !Initialize reconstruction vectors
    pered_order_index=-1._r8

    !Debug
    !if(ed/=93)then
    !  return
    !end if

    !Get cell surrounding edge
    cell(1:2)=mesh%edhx(ed)%sh(1:2)

    vtg=0
    do i=1,2
       ed_celli=0
       do j=1, mesh%v(cell(i))%nnb
          !Find the index of edge ed within cell(i)
          if(ed==mesh%v(cell(i))%ed(j))then
             ed_celli=j
             exit
          end if
       end do
       !print*
       !print*, "Cell:", cell(i), ed_celli
       if(ed_celli==0)then !check for errors
          print*, "pered_order_index error: edge not in cell"
          stop
       end if
       !For all edges of cell, add contribution to reconstruction
       do j=1, mesh%v(cell(i))%nnb
          !Edge's global index

          k=mesh%v(cell(i))%ed(j)
          ncor=real(mesh%hx(cell(i))%nr(j), r8)
          w=ncor*dot_product(mesh%edhx(k)%c%p-mesh%v(cell(i))%p, mesh%edhx(ed)%tg)
          w=w*mesh%edhx(k)%leng/(2.*mesh%hx(cell(i))%areag)

          !die=arclen(mesh%ed(k)%c%p, mesh%v(cell(i))%p)
          !print*, "Edge: ", j, k
          !ncor=1. !real(mesh%hx(cell(i))%nr(j), r8)
          !nrtg=dot_product(mesh%edhx(k)%nr*ncor, mesh%edhx(ed)%tg)
          !print*, "nrtg: ", nrtg
          !w=(mesh%edhx(k)%leng*die*nrtg)/(2.*mesh%hx(cell(i))%areag)
          !print*, w
          vtg=vtg+mesh%edhx(k)%nr*w
          !vtg=vtg+mesh%edhx(k)%nr*mesh%edhx(k)%leng*mesh%hx(cell(i))%trskw(ed_celli, j)
          !vtg=vtg+(mesh%edhx(k)%tg)*mesh%hx(cell(i))%tg(j)
          !print*, k
          !print*, (mesh%edhx(k)%tg)*mesh%hx(cell(i))%tg(j)
          !print*, vtg
          !print*, "Weight", k, ed, mesh%edhx(k)%leng*mesh%hx(cell(i))%trskw(ed_celli, j)/mesh%ed(ed)%leng
       end do
       !print*, norm(vtg)
       !print*, vtg
       !print*
    end do
    !vtg=vtg/mesh%ed(ed)%leng
    !vtg=proj_vec_sphere(vtg, mesh%edhx(ed)%c%p)
    !print*, norm(vtg)
    !print*, vtg
    !print*, mesh%ed(ed)%tg
    !print*,"----------------------------"
    !Get reconstructed vector including normal component
    pered_order_index=norm(vtg-mesh%edhx(ed)%tg)

    return
  end function pered_order_index

  function vecrecon_dtred(ed, var, mesh)
    !--------------------------------------------------------
    ! vecrecon_dtred
    ! Vector reconstruction with dual triangule formulation
    !--------------------------------------------------------
    !Reconstruction point
    integer(i4), intent(in) :: ed

    !Variable with normal vector values
    type(scalar_field), intent(in) :: var

    !Mesh structure
    type(grid_structure) :: mesh

    !Reconstructed vector
    real (r8), dimension(1:3):: vecrecon_dtred

    !Just tangent component of reconstructed vector
    real (r8):: utg
    !real (r8):: ncor
    !real (r8):: ncor_ed
    !real (r8):: w
    real (r8):: die(1:2)
    real (r8):: v1(1:3),  v2(1:3), v(1:3)

    !Cells adjacent to edge ed
    !integer (i4):: cell(1:2)

    !Indexes and aux vars
    !integer (i4):: i
    !integer (i4):: j
    !integer (i4):: k
    integer (i4):: ktr
    !integer (i4):: ed_celli

    !Check for argument entry error
    if(ed==0)then
       print*, "vecrecon_dtred error: Edge index must be informed"
       stop
    end if

    !Initialize reconstruction vectors
    vecrecon_dtred=0._r8
    utg=0._r8

    !if(ed .ne. 93)then
    !  return
    !end if

    !ncor_ed=dsign( 1._r8, dot_product(mesh%edhx(ed)%tg, mesh%ed(ed)%nr))

    !Works for if(var%pos==3)then  !Normal components given on hx edge midpoints
    !do i=1,2
    !       ktr=mesh%edhx(ed)%v(i)
    !       ed_celli=0
    !       do j=1, 3
    !          !Find the index of edge ed within cell(i)
    !          if(ed==mesh%tr(ktr)%ed(j))then
    !             ed_celli=j
    !             exit
    !          end if
    !       end do
    ktr=mesh%edhx(ed)%v(1)
    v1=vecrecon_perot(mesh%tr(ktr)%c%p, var, "tr  ", mesh, ktr)
    die(1)=norm(mesh%ed(ed)%c%p-mesh%tr(ktr)%c%p)
    ktr=mesh%edhx(ed)%v(2)
    v2=vecrecon_perot(mesh%tr(ktr)%c%p, var, "tr  ", mesh, ktr)
    die(2)=norm(mesh%ed(ed)%c%p-mesh%tr(ktr)%c%p)
    v=(die(1)*v1+die(2)*v2)/sum(die(1:2))

    !For all edges of triangle, add contribution to reconstruction
    !do j=1, 3
    !   !Edge's global index
    !   k=mesh%tr(ktr)%ed(j)
    !   !Reference vector
    !   v=cross_product(mesh%tr(ktr)%c%p,mesh%ed(k)%c%p-mesh%tr(ktr)%c%p)
    !   !Use the dual edge
    !   ncor=real(mesh%tr(ktr)%tg(j), r8)
    !   ncor=dsign( 1._r8, ncor* dot_product(mesh%edhx(k)%nr, mesh%ed(k)%tg))
    !   w=dot_product(v, mesh%ed(ed)%nr)
    !   w=w*mesh%ed(k)%leng
    !   !die=arclen(mesh%ed(k)%c%p,mesh%tr(ktr)%c%p)
    !w=die*mesh%ed(k)%leng* dot_product(ncor*mesh%ed(k)%tg, &
    !  real(mesh%tr(ktr)%tg(ed_celli), r8)*real(mesh%ed(ed)%nr))

    !w=ncor*dot_product(mesh%edhx(k)%c%p-mesh%v(cell(i))%p, mesh%edhx(ed)%tg)
    !   utg=utg+var%f(k)*w !*dsign( 1._r8, dot_product(mesh%edhx(k)%tg, mesh%ed(k)%nr))
    !end do

    !die(i)=norm(mesh%ed(ed)%c%p-mesh%tr(ktr)%c%p)
    !utg=die(i)*dot_product(v, mesh%edhx(ed)%tg) !*ncor_ed/(mesh%tr(ktr)%areag)
    !end do
    !utg=utg/sum(die(1:2))
    utg=dot_product(v, mesh%edhx(ed)%tg)
    !Get reconstructed vector including normal component
    vecrecon_dtred=utg*mesh%edhx(ed)%tg+var%f(ed)*mesh%edhx(ed)%nr

    return
  end function vecrecon_dtred

  function vecrecon_rt0(p, var, mesh)
    !--------------------------------------------------------
    ! vecrecon_rt0
    ! Raviart-Thomas Element Basis Reconstruction
    ! See Bahriawati and Carstensen 2005
    !--------------------------------------------------------
    !Reconstruction point
    real (r8), intent(in) :: p(1:3)

    !Variable with normal vector values
    type(scalar_field), intent(in) :: var

    !Mesh structure
    type(grid_structure), intent(in) :: mesh

    !Interpoled value
    real (r8):: vecrecon_rt0(1:3)

    !Aux var
    real (r8):: ncor
    real (r8):: tau(1:3)

    !Indexes and aux vars
    integer (i4):: i
    integer (i4):: j
    integer (i4):: k
    integer (i4):: l
    integer (i4):: ktr

    !Get triangle containing point
    ktr=gettr(p,mesh)

    !Get nearest node to point
    i=getnearnode(p, mesh, ktr)

    !Check if var is in the correct position
    if(var%pos/=2)then  !Normal components given not given on tr edge midpoints
       print*, "vecrecon_rt0 warning: vector field given in incorrect position", var%pos
       return
    else !Normal components given on tr edge midpoints
       !Do reconstruction
       vecrecon_rt0=(/0._r8,0._r8,0._r8/)
       !print*, "Triangle:", ktr
       do j=1,3
          !Triangle
          k=mesh%tr(ktr)%ed(j)
          !Outward unit vector correction
          ncor=real(mesh%tr(ktr)%nr(j), r8)
          !Oposite node to edge k
          l=mesh%tr(ktr)%v(modint(j-1,3))
          !Basis function for this edge
          !tau=(p-mesh%v(l)%p)/(2._r8*norm(mesh%ed(k)%c%p-mesh%v(l)%p))
          tau=(p-mesh%v(l)%p)*(mesh%ed(k)%leng)/(2._r8*mesh%tr(ktr)%areag)
          vecrecon_rt0=vecrecon_rt0 + ncor*(var%f(k))*tau
       end do
    end if

    return
  end function vecrecon_rt0

  function vecrecon_whitney(p, var, mesh)
    !--------------------------------------------------------
    ! vecrecon_whitney
    ! Edge Element - Whitney Basis Reconstruction
    ! See PDE and FE Methods - Solin 2005 - Eq 7.100
    ! E/ou http://www.iue.tuwien.ac.at/phd/nentchev/node40.html
    !--------------------------------------------------------
    !Reconstruction point
    real (r8), intent(in) :: p(1:3)

    !Variable with normal vector values
    type(scalar_field), intent(in) :: var

    !Mesh structure
    type(grid_structure), intent(in) :: mesh

    !Interpoled value
    real (r8):: vecrecon_whitney(1:3)

    !Normal correction
    real (r8):: ncor(1:3)

    !Basis functions
    real (r8):: theta1(1:3)
    real (r8):: theta2(1:3)
    real (r8):: theta3(1:3)

    !Normal vectors
    real (r8):: nr1(1:3)
    real (r8):: nr2(1:3)
    real (r8):: nr3(1:3)

    !Tangent vectors
    real (r8):: tg1(1:3)
    real (r8):: tg2(1:3)
    real (r8):: tg3(1:3)

    !Barycentric coords
    real (r8):: b(1:3)

    !Indexes and aux vars
    integer (i4):: i
    integer (i4):: l
    integer (i4):: ktr
    integer (i4):: ed

    vecrecon_whitney=(/0._r8,0._r8,0._r8/)

    !Check if var is in the correct position
    if(var%pos/=3)then  !Normal components given not given on hx edge midpoints
       print*, "vecrecon_whitney warning: vector field given in incorrect position", var%pos
       return
    end if

    !Get triangle containing point
    ktr=gettr(p,mesh)
    !Get nearest node to point
    i=getnearnode(p, mesh, ktr)

    !Barycentric coords
    b=bar_coord_tr(p, ktr, mesh)

    !Normal components given on hx edge midpoints
    !Reconctruction to TR circumcenter using tangential components
    !  given by normal components at hx edges

    !Calculate normal to triangle edges - unit outward
    nr1=mesh%ed(mesh%tr(ktr)%ed(1))%nr * mesh%tr(ktr)%nr(1)
    nr2=mesh%ed(mesh%tr(ktr)%ed(2))%nr * mesh%tr(ktr)%nr(2)
    nr3=mesh%ed(mesh%tr(ktr)%ed(3))%nr * mesh%tr(ktr)%nr(3)

    !Calculate tangents to triangle edges - unit ccwise
    tg1=mesh%ed(mesh%tr(ktr)%ed(1))%tg * mesh%tr(ktr)%tg(1)
    tg2=mesh%ed(mesh%tr(ktr)%ed(2))%tg * mesh%tr(ktr)%tg(2)
    tg3=mesh%ed(mesh%tr(ktr)%ed(3))%tg * mesh%tr(ktr)%tg(3)

    !Calculate basis functions
    theta1= ( b(3) * nr2 / dot_product(nr2,tg1) + &
         b(2) * nr3 / dot_product(nr3,tg1) )
    theta2= ( b(1) * nr3 / dot_product(nr3,tg2) + &
         b(3) * nr1 / dot_product(nr1,tg2) )
    theta3=   ( b(2) * nr1 / dot_product(nr1,tg3) + &
         b(1) * nr2 / dot_product(nr2,tg3) )

    !Because we are mixing HC normals and TR tangents,
    !  the corection needs some adjustment
    !print*
    do l=1,3
       ed=mesh%tr(ktr)%ed(l)
       ncor(l)=dsign( 1._r8, real(mesh%tr(ktr)%tg(l), r8)* &
            dot_product(mesh%edhx(ed)%nr, &
            mesh%ed(ed)%tg))
       !print*, ed, real(mesh%tr(ktr)%tg(l), r8), dot_product(mesh%edhx(ed)%nr, &
       ! mesh%ed(ed)%tg), ncor(l)
    end do
    !print*

    vecrecon_whitney=ncor(1)*(var%f(mesh%tr(ktr)%ed(1)))*theta1+ &
         ncor(2)*(var%f(mesh%tr(ktr)%ed(2)))*theta2+ &
         ncor(3)*(var%f(mesh%tr(ktr)%ed(3)))*theta3

    return
  end function vecrecon_whitney

  function vecrecon_klausen(p, var, mesh)
    !--------------------------------------------------------
    ! vecrecon_klausen
    ! Klausen, Rasmussen and Stephansen 2012
    !--------------------------------------------------------
    !Reconstruction point
    real (r8), intent(in) :: p(1:3)

    !Variable with normal vector values
    type(scalar_field), intent(in) :: var

    !Mesh structure
    type(grid_structure), intent(in) :: mesh

    !Interpoled value
    real (r8):: vecrecon_klausen(1:3)

    !Normal correction
    real (r8):: ncor

    !Basis functions
    real (r8):: psi(1:3)

    !Normalization terms
    real (r8):: delta1
    real (r8):: delta2

    !Tangent vectors
    type(vector):: tg(-1:1)

    !Edge indexes
    integer (i4):: ed(-1:1)

    !Wachspress coords
    type(general_coords):: wach

    !Indexes and aux vars
    integer (i4):: i
    integer (i4):: j

    vecrecon_klausen=(/0._r8,0._r8,0._r8/)

    !Check if var is in the correct position
    if(var%pos/=3)then  !Normal components given not given on hx edge midpoints
       print*, "vecrecon_klausen warning: vector field given in incorrect position", var%pos
       return
    end if

    !Get nearest node to point
    i=getnearnode(p, mesh)

    !Wachspress coords
    wach=wachspress_coords_hxv(p, mesh, i)

    !print*, "Node :", i

    !Normal components given on hx edge midpoints
    !Reconctruction to TR vert, HX center using normal HX edge components

    do j=1, mesh%v(i)%nnb

       !Set edge indexes
       ed(-1)=mesh%v(i)%ed(modint(j-1, mesh%v(i)%nnb))
       ed( 0)=mesh%v(i)%ed(j)
       ed(+1)=mesh%v(i)%ed(modint(j+1, mesh%v(i)%nnb))

       !print*, j
       !print*, ed
       !Set tangent vectors
       tg(-1)%v=mesh%edhx(ed(-1))%tg*mesh%hx(i)%tg(modint(j-1, mesh%v(i)%nnb))
       tg( 0)%v=mesh%edhx(ed( 0))%tg*mesh%hx(i)%tg(j)
       tg(+1)%v=mesh%edhx(ed(+1))%tg*mesh%hx(i)%tg(modint(j+1, mesh%v(i)%nnb))

       !Normalization terms
       delta1=norm(cross_product(tg(-1)%v, tg( 0)%v ))
       delta2=norm(cross_product(tg(+1)%v, tg( 0)%v ))

       !Basis function
       psi=wach%w(modint(j-1, mesh%v(i)%nnb))*tg(-1)%v/delta1 - &
            wach%w(j)*tg(1)%v/delta2
       !print*, wach%w(modint(j-1, mesh%v(i)%nnb)),    wach%w(j)*tg(1)%v
       !Correction term
       ncor=mesh%hx(i)%nr(j)

       !Update reconstructed vector
       vecrecon_klausen = vecrecon_klausen + ncor*psi*var%f(ed(0))

    end do

    return
  end function vecrecon_klausen

  function vecrecon_trsk(ed, var, mesh)
    !--------------------------------------------------------
    ! vecrecon_trisk
    ! Vector reconstruction - tangent component - Thuburn et all (2009)
    !   TRISK
    !--------------------------------------------------------
    !Reconstruction point
    integer(i4), intent(in) :: ed

    !Variable with normal vector values
    type(scalar_field), intent(in) :: var

    !Mesh structure
    type(grid_structure) :: mesh

    !Method to be used
    !  "b" weights calculated with barycentric coord
    !  "o" or " " original Trisk scheme
    !character, optional :: mtd_in
    !character :: mtd

    !Reconstructed vector
    real (r8), dimension(1:3):: vecrecon_trsk

    !Just tangent component of reconstructed vector
    real (r8):: utg

    !Cells adjacent to edge ed
    integer (i4):: cell(1:2)

    !Indexes and aux vars
    integer (i4):: i
    integer (i4):: j
    integer (i4):: k
    integer (i4):: ed_celli
    integer (i4):: signcor

    !Check for argument entry error
    if(ed==0)then
       print*, "vecrecon_trsk error: Edge index must be informed"
       stop
    end if


    !Initialize reconstruction vectors
    vecrecon_trsk=0._r8
    utg=0._r8

    !Get cell surrounding edge
    cell(1:2)=mesh%edhx(ed)%sh(1:2)

    !Calculate primal-dual cell area ratios, Div weights and TRSK weights
    if(.not. allocated(mesh%hx(mesh%nv)%trskw))then
       call calc_trisk_weights(mesh)
    end if
    !print*, var%pos
    !stop
    !Works for if(var%pos==3 or 6)then  !Normal components given on hx edge midpoints
    !See Thuburn et al 2009 eq 8
    do i=1,2
       ed_celli=0
       do j=1, mesh%v(cell(i))%nnb
          !Find the index of edge ed within cell(i)
          if(ed==mesh%v(cell(i))%ed(j))then
             ed_celli=j
             exit
          end if
       end do
       if(ed_celli==0)then !check for errors
          print*, "vecrecon_trisk error: edge not in cell"
          stop
       end if
       !For all edges of cell, add contribution to reconstruction
       do j=1, mesh%v(cell(i))%nnb
          !Edge's global index
          k=mesh%v(cell(i))%ed(j)
          !Apply sign corrections
          if(var%pos==3)then
             signcor=-1.*mesh%hx(cell(i))%nr(j)*mesh%hx(cell(i))%tg(ed_celli)
          elseif(var%pos==6)then
             signcor=-1.*mesh%hx(cell(i))%ttgout(j)*mesh%hx(cell(i))%tnrccw(ed_celli)
          endif
          utg=utg+var%f(k)*mesh%edhx(k)%leng*mesh%hx(cell(i))%trskw(ed_celli, j)*signcor
          !Correct vector direction
          !mesh%hx(cell)%trskw(i,j)=-mesh%hx(cell)%trskw(i,j)*mesh%hx(cell)%nr(j)
          !Correct so that the tangent on edge ed is ccw
          !mesh%hx(cell)%trskw(i,j)=mesh%hx(cell)%trskw(i,j)*mesh%hx(cell)%tg(i)
       end do
    end do
    utg=utg/mesh%ed(ed)%leng

    !Get reconstructed vector including normal component
    if(var%pos==3)then
       vecrecon_trsk=utg*mesh%edhx(ed)%tg+var%f(ed)*mesh%edhx(ed)%nr
    elseif(var%pos==6)then
       vecrecon_trsk=utg*mesh%ed(ed)%nr+var%f(ed)*mesh%ed(ed)%tg
    else
       vecrecon_trsk=0.
    endif

    return
  end function vecrecon_trsk

  subroutine calc_trisk_weights(mesh)
    !---------------------------------------------------------
    !Calculates TRISK weights for Voronoi cell i
    ! Returns weight in mesh
    !-----------------------------------------------------------

    !Mesh structure
    type(grid_structure), intent(inout) :: mesh

    !Number of edges in cell
    integer(i4)::n

    !Index for edges
    integer(i4)::cell, i, j, k
    integer(i4)::ed1, ed2

    !Method to be used
    !  "b" weights calculated with barycentric coord
    !  "o" or " " original Trisk scheme
    !character, optional :: mtd_in
    !character :: mtd

    !Area ratios
    real (r8) :: R(1:10)
    real (r8) :: Rsum
    real (r8) :: alpha
    !real (r8) :: bcoord(1:3)

    !Calculate areas of intersection between triangles and voronoi cells
    !  required for TRISK schemes, see Thuburn 2009
    do i=1, mesh%nv
       n=mesh%v(i)%nnb
       if(.not. allocated(mesh%hx(i)%hxtr_area))then
          allocate(mesh%hx(i)%hxtr_area(1:n))
       end if
       mesh%hx(i)%hxtr_area= hxtr_intersec_areas(i, mesh%v(i)%nnb, mesh)
       mesh%hx(i)%hxtr_area(1:n)=mesh%hx(i)%hxtr_area(1:n)/mesh%hx(i)%areag
       ! print*, i, allocated(mesh%hx(i)%hxtr_area), mesh%hx(i)%hxtr_area(1:mesh%v(i)%nnb)
    end do

    !Calculate barycentric coordinates of cell vertex (tr cc)
    ! to be used as weights (for modified method only)
    !do i=1, mesh%nv
    ! print*, " Cell ", i
    !   if(.not. allocated(mesh%hx(i)%hxtr_remapw))then
    !      allocate(mesh%hx(i)%hxtr_remapw(1:mesh%v(i)%nnb))
    !   end if
    !   mesh%hx(i)%hxtr_remapwsum=0._r8
    !   do j=1, mesh%v(i)%nnb
    !      k=mesh%v(i)%tr(j)
    !      bcoord=bar_coord_tr(mesh%tr(k)%c%p, k, mesh)
    !      !print*, "  Nb: ", j," Tr: ", k
    !      !print*, bcoord
    !      !print*, mesh%tr(k)%v(1:3)
    !      do l=1, 3
    !        if(mesh%tr(k)%v(l)==i)then
    !          mesh%hx(i)%hxtr_remapw(j)=bcoord(l)*mesh%tr(k)%areag/mesh%hx(i)%areag
    !        end if
    !      end do
    !      mesh%hx(i)%hxtr_remapwsum=mesh%hx(i)%hxtr_remapwsum+mesh%hx(i)%hxtr_remapw(j)
    !      !print*, mesh%hx(i)%hxtr_remapw(j)
    !   end do
    !   print*, "Cell: ", i, " Rsum:", mesh%hx(i)%hxtr_remapwsum
    !   !print*
    !end do

    !See Thuburn et al 2009 eq 33
    !or Weller 2012
    do cell=1, mesh%nv
       n=mesh%v(cell)%nnb
       if(.not. allocated(mesh%hx(cell)%trskw))then
          allocate(mesh%hx(cell)%trskw(1:n, 1:n))
       end if
       !Check which weight is to be used
       !if(mtd=="o")then
       R(1:n)=mesh%hx(cell)%hxtr_area(1:n)
       alpha=0.5_r8
       !print*, sum(R)
       !else
       !  R(1:n)=mesh%hx(cell)%hxtr_remapw(1:n)
       !  alpha=mesh%hx(cell)%hxtr_remapwsum /2._r8
       !end if
       !print*
       !print*, "Cell: ", cell
       !For each edge
       do i=1, n

          ed1=mesh%v(cell)%ed(i)
          !print*, "Edge i: ", i, ed1

          !Calculate the weight relative to an other edge
          do j=1,n
             Rsum=0._r8
             if(i==j)then !Null weight
                mesh%hx(cell)%trskw(i,i)=0._r8
                cycle
             else if(j<i)then
                do k=i, n
                   Rsum=Rsum + R(k)
                end do
                do k=1, j-1
                   Rsum=Rsum + R(k)
                end do
             else !(i<j)
                do k=i, j-1
                   Rsum=Rsum + R(k)
                end do
             end if
             !Rsum=Rsum/mesh%hx(cell)%areag
             ed2=mesh%v(cell)%ed(j)
             mesh%hx(cell)%trskw(i,j)=0._r8
             !Subtract 1/2
             mesh%hx(cell)%trskw(i,j)=(Rsum-alpha)
             !print*, " Edge j, w(i,j): ", j, mesh%hx(cell)%trskw(i,j), mesh%hx(cell)%nr(j)
             !print*, "  w(i,j): ", j, ed2, mesh%hx(cell)%trskw(i,j), mesh%hx(cell)%tg(i), mesh%hx(cell)%nr(j), Rsum
             !mesh%hx(cell)%trskw(j,i)=-mesh%hx(cell)%trskw(i,j)
          end do
       end do
    end do

  end subroutine calc_trisk_weights

  function trsk_order_index(ed, pos, mesh)
    !--------------------------------------------------------
    ! trsk_order_index
    ! Calculated the consistency conditions for trisk and returns an index
    !
    !--------------------------------------------------------
    !Reconstruction point
    integer(i4), intent(in) :: ed

    !Mesh structure
    type(grid_structure) :: mesh

    !Reconstructed vector
    real (r8) :: trsk_order_index

    !Reconstructed tangent vector
    real (r8) :: vtg(1:3)

    !Cells adjacent to edge ed
    integer (i4):: cell(1:2)

    !Position of normals (hx ed midpoint==3, tr x hx intersec ==6)
    integer (i4):: pos

    !Indexes and aux vars
    integer (i4):: i
    integer (i4):: j
    integer (i4):: k
    integer (i4):: ed_celli
    integer (i4):: signcor

    !Check for argument entry error
    if(ed==0)then
       print*, "vecrecon_trsk error: Edge index must be informed"
       stop
    end if

    !Initialize reconstruction vectors
    trsk_order_index=-1._r8

    !Get cell surrounding edge
    cell(1:2)=mesh%edhx(ed)%sh(1:2)

    !Calculate primal-dual cell area ratios, Div weights and TRSK weights
    if(.not. allocated(mesh%hx(mesh%nv)%trskw))then
       call calc_trisk_weights(mesh)
    end if
    !print*, "------------------"
    !print*, "Edge:", ed
    !print*, cell
    !Works for if(var%pos==3)then  !Normal components given on hx edge midpoints

    vtg=0
    do i=1,2

       ed_celli=0
       do j=1, mesh%v(cell(i))%nnb
          !Find the index of edge ed within cell(i)
          if(ed==mesh%v(cell(i))%ed(j))then
             ed_celli=j
             exit
          end if
       end do
       if(ed_celli==0)then !check for errors
          print*, "trsk_order_index error: edge not in cell"
          stop
       end if
       !For all edges of cell, add contribution to reconstruction
       do j=1, mesh%v(cell(i))%nnb
          !Edge's global index
          k=mesh%v(cell(i))%ed(j)
          !Apply sign corrections
          if(pos==3)then
             signcor=-1.*mesh%hx(cell(i))%nr(j)*mesh%hx(cell(i))%tg(ed_celli)
             vtg=vtg+mesh%edhx(k)%nr*mesh%edhx(k)%leng*mesh%hx(cell(i))%trskw(ed_celli, j)*signcor
          elseif(pos==6)then
             signcor=-1.*mesh%hx(cell(i))%ttgout(j)*mesh%hx(cell(i))%tnrccw(ed_celli)
             vtg=vtg+mesh%ed(k)%tg*mesh%edhx(k)%leng*mesh%hx(cell(i))%trskw(ed_celli, j)*signcor
          endif
          !utg=utg+var%f(k)*mesh%edhx(k)%leng*mesh%hx(cell(i))%trskw(ed_celli, j)*signcor

          !vtg=vtg+(mesh%edhx(k)%tg)*mesh%hx(cell(i))%tg(j)
          !print*, k
          !print*, (mesh%edhx(k)%tg)*mesh%hx(cell(i))%tg(j)
          !print*, vtg
          !print*, "Weitght", k, ed, mesh%edhx(k)%leng*mesh%hx(cell(i))%trskw(ed_celli, j)/mesh%ed(ed)%leng
       end do
       !print*, norm(vtg)
       !print*, vtg
       !print*
    end do
    vtg=vtg/mesh%ed(ed)%leng
    !vtg=proj_vec_sphere(vtg, mesh%edhx(ed)%c%p)
    !print*, norm(vtg)
    !print*, vtg
    !print*, mesh%ed(ed)%tg
    !print*,"----------------------------"
    !Get reconstructed vector including normal component
    if(pos==3)then
       trsk_order_index=norm(vtg-mesh%edhx(ed)%tg)
    elseif(pos==6)then
       trsk_order_index=norm(vtg-mesh%ed(ed)%nr)
    endif


    return
  end function trsk_order_index


  !=====================================================================================
  !    RADIAL BASIS FUNCTIONS ROUTINES
  !=====================================================================================

  function rbf(r, h)
    !---------------------------------------------------------
    !RBF
    !  Calculates radial basis function
    !-----------------------------------------------------------
    !Distance (radius) variable
    real(r8), intent(in) :: r

    !Mesh characteristic factor (ex: min dist between nodes)
    real(r8), intent(in) :: h

    !RBF value
    real (r8):: rbf

    !Shape parameter
    real (r8):: e

    !r must be positive or zero
    if(r<0)then
       print*, "RBF error: negative radius", r
       stop
    end if

    rbf=h

    !Inverse multiquadratic
    !e= 2 !1._r8/h
    !e= 1._r8/(h*8._r8)
    !e= 4._r8/(1._r8)
    !rbf=1._r8/dsqrt(1+(e*r)**2)

    !Gaussian
    !e= 1._r8/(h*16._r8)
    !e= 1._r8/h
    !e= 1._r8
    !e=1._r8
    e=h
    !e= 1._r8/(1._r8)
    rbf=exp(-(e*r)**2)

    ! Cubic Spline based on
    ! en.wikipedia.org/wiki/Spline_%28mathematics%29
    !x=abs(r/h)/8.
    !if(r<1)then
    !    rbf=(3._r8*x**3-6._r8*x**2+4._r8)/4._r8
    !elseif(r<2)then
    !    rbf=((2._r8-x)**3)/4._r8
    !else
    !    rbf=0._r8
    !end if

    return
  end function rbf

  subroutine rbf_matrix_build( stencil, rbf_h, rbf_mat, mesh )
    !---------------------------------------------------------
    !RBF_MATRIX_BUILD
    !  Builds the structure for rbf matrices
    !-----------------------------------------------------------
    !Stencil or type of rbf point structure
    ! TR   = Points on triangles - 3 points only
    ! ETR = (Extended TR) Points neighbour to triangle also included
    ! HX   = Hexagonal, pentagonal neighbour points
    character(len=4), intent(in) :: stencil

    !Shape parameter
    real (r8):: rbf_h

    !Mesh structure
    type(grid_structure), intent(in) :: mesh

    !RBF matrix vector
    type(rbf_matrix_structure), allocatable :: rbf_mat(:)

    !Rbf relational matrix
    real (r8), allocatable :: A(:,:)

    !Tmp vector
    integer (i4), allocatable :: ptmp(:)

    !Indexes and aux vars
    integer (i4):: i
    integer (i4):: j
    integer (i4):: jj
    integer (i4):: k
    integer (i4):: l
    integer (i4):: n
    logical:: add

    print*
    print*, "Constructing RBF matrix for stencil:", stencil
    print*

    !Check if already matrix allocated
    if(allocated(rbf_mat))then
       do l=1, size(rbf_mat)
          deallocate(rbf_mat(l)%L, rbf_mat(l)%b,rbf_mat(l)%w,rbf_mat(l)%pts)
       end do
       deallocate(rbf_mat)
    end if

    !Build matrix for selected stencil
    select case(trim(stencil))

       !Poins given on triangle vertices, only 3 points used for each rbf
    case("TR" , "tr", "trpc", "TRPC")

       !Allocate space for matrices
       !There will be number of triangles rbf matrices
       allocate(rbf_mat(1:mesh%nt))

       !Save temporary A matrix
       allocate(A(1:3,1:3))

       !Create matrix for all triangles
       do k=1,mesh%nt

          !Allocate matrices
          rbf_mat(k)%n=3
          allocate(rbf_mat(k)%L(1:rbf_mat(k)%n,1:rbf_mat(k)%n))
          allocate(rbf_mat(k)%pts(1:rbf_mat(k)%n))
          allocate(rbf_mat(k)%b(1:rbf_mat(k)%n))
          allocate(rbf_mat(k)%w(1:rbf_mat(k)%n))
          !Save rbf points
          rbf_mat(k)%pts(1:3)=mesh%tr(k)%v(1:3)
          rbf_mat(k)%h=rbf_h

          A=0._r8
          !Calculate distances
          do i=1, rbf_mat(k)%n
             do j=1, i
                !Geodesic distance
                A(i,j)= rbf( &
                     arclen(mesh%v(rbf_mat(k)%pts(i))%p, mesh%v(rbf_mat(k)%pts(j))%p), &
                     rbf_mat(k)%h)
                !Euclidian distance
                !A(i,j)=rbf( &
                !norm(mesh%v(rbf_mat(k)%pts(i))%p - mesh%v(rbf_mat(k)%pts(j))%p),
                !rbf_mat(k)%h)
             end do
          end do
          ! print*, k
          ! print '(3f12.6)', transpose(A)
          ! print*

          !Make the cholesky decomposition
          call choleskydecomp(A, 3, rbf_mat(k)%L)

          !Calculate estimative of condition number
          rbf_mat(k)%condnum=condnumest(rbf_mat(k)%L, rbf_mat(k)%n)
          ! print '(3f12.6)', transpose(rbf_mat(k)%L)
          ! print*
          !print *, "Cond num:", rbf_mat(k)%condnum
          !print*
          !print*
          !print*,"------------------------------------------------"
          !print*

       end do

       !Points given on triangle vertices, uses TR neighbours as well
    case("ETR", "etr", "ETRP", "etrp" )

       !Allocate space for matrices
       !There will be number of triangles rbf matrices
       allocate(rbf_mat(1:mesh%nt))

       !Save temporary A matrix
       n=3*mesh%maxvnb-4
       allocate(A(1:n,1:n))
       allocate(ptmp(1:n))
       !print*, n
       !Create matrix for all triangles
       do k=1,mesh%nt
          !print*, "triangle : ", k

          !Verify number of points to be included
          !Save the triangle edges as first on the list
          !ptmp(1:3)=mesh%tr(k)%ed(1:3)
          ptmp(1:3)=mesh%tr(k)%v(1:3)
          !print*, ptmp(1:3)
          !Go for the next nodes
          l=4
          !For each triangle vertice, get neighbours
          do i=1,3
             !For each neighbour check if it is already on the list
             ! if not, add it
             do j=1,mesh%v(ptmp(i))%nnb
                ptmp(l)=mesh%v(ptmp(i))%nb(j)
                !We initially supose it is not on the list,
                ! so it is "true" to add it
                add=.true.
                !Check if it is already on the list
                do jj=1,l-1
                   if(ptmp(l)==ptmp(jj)) add=.false.
                end do
                !If node is to be added on the list, go to the next one
                ! if not, the next one will overwrite this one
                if(add)l=l+1
             end do
          end do
          l=l-1

          !Allocate matrices
          rbf_mat(k)%n=l
          allocate(rbf_mat(k)%L(1:rbf_mat(k)%n,1:rbf_mat(k)%n))
          allocate(rbf_mat(k)%pts(1:rbf_mat(k)%n))
          allocate(rbf_mat(k)%b(1:rbf_mat(k)%n))
          allocate(rbf_mat(k)%w(1:rbf_mat(k)%n))
          !Save rbf points
          rbf_mat(k)%pts(1:l)=ptmp(1:l)
          rbf_mat(k)%h=rbf_h
          !print*, "num of points:", l
          !print*, ptmp(1:l)
          A=0._r8
          !Calculate distances
          do i=1, rbf_mat(k)%n
             do j=1, i
                !Geodesic distance
                A(i,j)= rbf( &
                     arclen(mesh%v(rbf_mat(k)%pts(i))%p, mesh%v(rbf_mat(k)%pts(j))%p), &
                     rbf_mat(k)%h)
                !Euclidian distance
                !A(i,j)=rbf( &
                !norm(mesh%v(rbf_mat(k)%pts(i))%p - mesh%v(rbf_mat(k)%pts(j))%p),
                !rbf_mat(k)%h)
             end do
          end do
          !print*, k
          !fmt="(XXf16.8)"
          !write(fmt(2:3),'(i2)') rbf_mat(k)%n
          !print*, fmt
          !print fmt, transpose(A(1:rbf_mat(k)%n,1:rbf_mat(k)%n))
          !print*

          !Make the cholesky decomposition
          call choleskydecomp(A(1:rbf_mat(k)%n,1:rbf_mat(k)%n), rbf_mat(k)%n, rbf_mat(k)%L)

          !Calculate estimative of condition number
          rbf_mat(k)%condnum=condnumest(rbf_mat(k)%L, rbf_mat(k)%n)
          !print fmt, transpose(rbf_mat(k)%L)
          !print*
          !print *, "Cond num:", rbf_mat(k)%condnum
          !print*
          !print*,"------------------------------------------------"
          !print*

       end do

       !Poins given on triangle vertices, but assuming hexagon as center
    case("HX" , "hx", "hxpc", "HXPC" )

       !Allocate space for matrices
       !There will be 'number of nodes' rbf matrices
       allocate(rbf_mat(1:mesh%nv))

       !Save temporary A matrix
       n=mesh%maxvnb+1
       allocate(A(1:n,1:n))
       allocate(ptmp(1:n))

       !Create matrix for all hexagons
       do k=1,mesh%nv
          !Get all hexagons/pentagon indexes
          ptmp(1)=k
          do l=1,mesh%v(k)%nnb
             ptmp(l+1)=mesh%v(k)%nb(l)
          end do
          !print*, k, ptmp(1:l)

          !Allocate matrices
          rbf_mat(k)%n=mesh%v(k)%nnb+1
          allocate(rbf_mat(k)%L(1:rbf_mat(k)%n,1:rbf_mat(k)%n))
          allocate(rbf_mat(k)%pts(1:rbf_mat(k)%n))
          allocate(rbf_mat(k)%b(1:rbf_mat(k)%n))
          allocate(rbf_mat(k)%w(1:rbf_mat(k)%n))

          !Save rbf points
          rbf_mat(k)%pts(1:rbf_mat(k)%n)=ptmp(1:rbf_mat(k)%n)
          rbf_mat(k)%h=rbf_h

          A=0._r8
          !Calculate distances
          do i=1, rbf_mat(k)%n
             do j=1, i
                !Geodesic distance
                A(i,j)= rbf( &
                     arclen(mesh%v(rbf_mat(k)%pts(i))%p, mesh%v(rbf_mat(k)%pts(j))%p), &
                     rbf_mat(k)%h)
                !Euclidian distance
                !A(i,j)=rbf( &
                !norm(mesh%v(rbf_mat(k)%pts(i))%p - mesh%v(rbf_mat(k)%pts(j))%p),
                !rbf_mat(k)%h)
             end do
          end do
          !print*, k
          !fmt="(XXf16.8)"
          !write(fmt(2:3),'(i2)') rbf_mat(k)%n
          !print*, fmt
          !print fmt, transpose(A(1:rbf_mat(k)%n,1:rbf_mat(k)%n))
          !print*

          !Make the cholesky decomposition
          call choleskydecomp(A(1:rbf_mat(k)%n,1:rbf_mat(k)%n),rbf_mat(k)%n , rbf_mat(k)%L)

          !Calculate estimative of condition number
          rbf_mat(k)%condnum=condnumest(rbf_mat(k)%L, rbf_mat(k)%n)

          !print fmt, transpose(rbf_mat(k)%L)
          !print*
          !print *, "Cond num:", rbf_mat(k)%condnum
          !print*
          !print*,"------------------------------------------------"
          !print*

       end do

    case default
       print*, "RBF_MATRIX_BUILD ERROR : Unknown RBF stencil. ", stencil
       stop
    end select

    return
  end subroutine rbf_matrix_build

  subroutine rbf_addpol(rbf_mat_local)
    !--------------------------------------------------------
    ! rbfvec_addpol
    ! Adds a polynomial to a RBF vector reconstruction
    ! Must receive rbf%mat_local%w with weights (A^(-1)b)
    !--------------------------------------------------------
    !RBF variable structure
    type(rbf_matrix_structure), intent(inout) :: rbf_mat_local

    !Mesh structure
    !type(grid_structure), intent(in) :: mesh

    !Vector with ones
    real (r8), allocatable :: c(:)

    !Auxiliar vector
    real (r8), allocatable :: w(:)

    !Auxiliar matrix
    !real (r8), dimension(1:3) :: l1, l2, l3, d

    !Indexes and aux vars
    integer (i4):: n

    n=rbf_mat_local%n

    !Initialize polynomial vector
    if(.not.allocated(rbf_mat_local%pol))then
       allocate(rbf_mat_local%pol(1:1))
    end if
    rbf_mat_local%pol=0

    allocate(c(1:n))
    allocate(w(1:n))

    !Build problem
    c(1:n)=1

    !Solve Aw=c
    call choleskysolve(rbf_mat_local%L, w, c, n)

    !print*, "pol"
    !print*, rbf_mat_local%pol
    rbf_mat_local%pol(1)= sum(rbf_mat_local%w(1:n))/sum(w(1:n))

    !print*, "Old weights"
    !print*, rbf_mat_local%w

    !Correct original rbf weights with pol weights
    rbf_mat_local%w=rbf_mat_local%w -rbf_mat_local%pol(1)*w

    !print*, "New weights"
    !print*, rbf_mat_local%w
    !print*
    !print*
    !read*,i
    return
  end subroutine rbf_addpol

  subroutine rbfvec_matrix_build( stencil, rbf_h, rbf_mat, mesh )
    !---------------------------------------------------------
    !RBFVEC_MATRIX_BUILD
    !  Builds the structure for rbf vector reconstruction matrices
    !-----------------------------------------------------------
    !Stencil or type of rbf point structure
    ! HX   = Hexagonal, pentagonal neighbour points
    character(len=4), intent(in) :: stencil

    !Shape parameter
    real (r8):: rbf_h

    !Mesh structure
    type(grid_structure), intent(in) :: mesh

    !Rbf relational matrix
    real (r8), allocatable :: A(:,:)

    !RBF matrix vector
    type(rbf_matrix_structure), intent(out), allocatable :: rbf_mat(:)

    !Tmp vector
    integer (i4), allocatable :: ptmp(:)

    !Indexes and aux vars
    integer (i4):: i
    integer (i4):: j
    integer (i4):: k
    integer (i4):: l
    integer (i4):: n
    integer (i4):: kn

    real (r8):: tmp

    print*
    print*, "Constructing RBF VEC matrix for stencil:", stencil
    print*

    !Build matrix for selected stencil
    select case(trim(stencil))
       !Poins given on hexagonal edge midpoints
    case("HX" , "hx", "HXPC", "hxpc" )

       !Allocate space for matrices
       !There will be 'number of nodes' rbf matrices
       allocate(rbf_mat(1:mesh%nv))

       !Save temporary A matrix
       n=mesh%maxvnb
       allocate(A(1:n,1:n))
       allocate(ptmp(1:n))

       !Create matrix for all hexagons
       do k=1,mesh%nv
          !Get all edge indexes
          !ptmp(1)=k
          do l=1,mesh%v(k)%nnb
             ptmp(l)=mesh%v(k)%ed(l)
          end do
          !print*, k, ptmp(1:l)

          !Allocate matrices
          rbf_mat(k)%n=mesh%v(k)%nnb
          allocate(rbf_mat(k)%L(1:rbf_mat(k)%n,1:rbf_mat(k)%n))
          allocate(rbf_mat(k)%pts(1:rbf_mat(k)%n))
          allocate(rbf_mat(k)%b(1:rbf_mat(k)%n))
          allocate(rbf_mat(k)%w(1:rbf_mat(k)%n))

          !Save rbf edges
          rbf_mat(k)%pts(1:rbf_mat(k)%n)=ptmp(1:rbf_mat(k)%n)
          rbf_mat(k)%h=rbf_h

          A=0._r8
          !Calculate distances
          do i=1, rbf_mat(k)%n
             do j=1, i
                !Geodesic distance
                A(i,j)= rbf( &
                     arclen(mesh%edhx(rbf_mat(k)%pts(i))%c%p, &
                     mesh%edhx(rbf_mat(k)%pts(j))%c%p), &
                     rbf_mat(k)%h)
                !Calculate normal components dot product
                tmp=dot_product(mesh%edhx(rbf_mat(k)%pts(i))%nr, &
                     mesh%edhx(rbf_mat(k)%pts(j))%nr)
                A(i,j)=A(i,j)*tmp
             end do
          end do
          !print*, k
          !fmt="(XXf16.8)"
          !write(fmt(2:3),'(i2)') rbf_mat(k)%n
          !print*, fmt
          !print fmt, transpose(A(1:rbf_mat(k)%n,1:rbf_mat(k)%n))
          !print*

          !Make the cholesky decomposition
          call choleskydecomp(A(1:rbf_mat(k)%n,1:rbf_mat(k)%n),rbf_mat(k)%n , rbf_mat(k)%L)

          !Calculate estimative of condition number
          rbf_mat(k)%condnum=condnumest(rbf_mat(k)%L, rbf_mat(k)%n)

          !print fmt, transpose(rbf_mat(k)%L)
          !print*
          !print *, "Cond num:", rbf_mat(k)%condnum
          !print*
          !print*,"------------------------------------------------"
          !print*

       end do

    case("TR" , "tr", "trpc", "TRPC" )

       !Allocate space for matrices
       !There will be 'number of triangles' rbf matrices
       allocate(rbf_mat(1:mesh%nt))

       !Save temporary A matrix
       n=3
       allocate(A(1:n,1:n))
       allocate(ptmp(1:n))

       !Create matrix for all triangles
       do k=1,mesh%nt
          !Get all edge indexes
          !ptmp(1)=k
          do l=1, n
             ptmp(l)=mesh%tr(k)%ed(l)
          end do
          !print*, k, ptmp(1:l)

          !Allocate matrices
          rbf_mat(k)%n=n
          allocate(rbf_mat(k)%L(1:rbf_mat(k)%n,1:rbf_mat(k)%n))
          allocate(rbf_mat(k)%pts(1:rbf_mat(k)%n))
          allocate(rbf_mat(k)%b(1:rbf_mat(k)%n))
          allocate(rbf_mat(k)%w(1:rbf_mat(k)%n))

          !Save rbf edges
          rbf_mat(k)%pts(1:rbf_mat(k)%n)=ptmp(1:rbf_mat(k)%n)
          rbf_mat(k)%h=rbf_h

          A=0._r8
          !Calculate distances
          do i=1, rbf_mat(k)%n
             do j=1, i
                !Geodesic distance
                A(i,j)= rbf( &
                     arclen(mesh%edhx(rbf_mat(k)%pts(i))%c%p, &
                     mesh%edhx(rbf_mat(k)%pts(j))%c%p), &
                     rbf_mat(k)%h)
                !Calculate normal components dot product
                tmp=dot_product(mesh%edhx(rbf_mat(k)%pts(i))%nr, &
                     mesh%edhx(rbf_mat(k)%pts(j))%nr)
                A(i,j)=A(i,j)*tmp
             end do
          end do
          !print*, k
          !fmt="(XXf16.8)"
          !write(fmt(2:3),'(i2)') rbf_mat(k)%n
          !print*, fmt
          !print fmt, transpose(A(1:rbf_mat(k)%n,1:rbf_mat(k)%n))
          !print*

          !Make the cholesky decomposition
          call choleskydecomp(A(1:rbf_mat(k)%n,1:rbf_mat(k)%n),rbf_mat(k)%n , rbf_mat(k)%L)

          !Calculate estimative of condition number
          rbf_mat(k)%condnum=condnumest(rbf_mat(k)%L, rbf_mat(k)%n)

          !print fmt, transpose(rbf_mat(k)%L)
          !print*
          !print *, "Cond num:", rbf_mat(k)%condnum
          !print*
          !print*,"------------------------------------------------"
          !print*

       end do

    case("ETR" , "etr", "etrp", "ETRP" )

       !Allocate space for matrices
       !There will be 'number of triangles' rbf matrices
       allocate(rbf_mat(1:mesh%nt))

       !Save temporary A matrix
       n=9
       allocate(A(1:n,1:n))
       allocate(ptmp(1:n+3))

       !Create matrix for all triangles
       do k=1,mesh%nt
          !print*, "triangle : ", k

          !For each neighbour triangle, get the edges
          j=0
          do i = 1, 3
             !Neighbour triangle
             kn=mesh%tr(k)%nb(i)
             do l=1, 3
                j=j+1
                !Save index of edge
                ptmp(j)=mesh%tr(kn)%ed(l)
             end do
          end do
          if(j/=n)then
             print*, "Error on rbfvec_matrix_build ", j
             stop
          end if

          !Allocate matrices
          rbf_mat(k)%n=n
          allocate(rbf_mat(k)%L(1:rbf_mat(k)%n,1:rbf_mat(k)%n))
          allocate(rbf_mat(k)%pts(1:rbf_mat(k)%n))
          allocate(rbf_mat(k)%b(1:rbf_mat(k)%n))
          allocate(rbf_mat(k)%w(1:rbf_mat(k)%n))

          !Save rbf edges
          rbf_mat(k)%pts(1:rbf_mat(k)%n)=ptmp(1:rbf_mat(k)%n)
          rbf_mat(k)%h=rbf_h

          A=0._r8
          !Calculate distances
          do i=1, rbf_mat(k)%n
             do j=1, i
                !Geodesic distance
                A(i,j)= rbf( &
                     arclen(mesh%edhx(rbf_mat(k)%pts(i))%c%p, &
                     mesh%edhx(rbf_mat(k)%pts(j))%c%p), &
                     rbf_mat(k)%h)
                !Calculate normal components dot product
                tmp=dot_product(mesh%edhx(rbf_mat(k)%pts(i))%nr, &
                     mesh%edhx(rbf_mat(k)%pts(j))%nr)
                A(i,j)=A(i,j)*tmp
             end do
          end do
          !print*, k
          !fmt="(XXf16.8)"
          !write(fmt(2:3),'(i2)') rbf_mat(k)%n
          !print*, fmt
          !print fmt, transpose(A(1:rbf_mat(k)%n,1:rbf_mat(k)%n))
          !print*

          !Make the cholesky decomposition
          call choleskydecomp(A(1:rbf_mat(k)%n,1:rbf_mat(k)%n),rbf_mat(k)%n , rbf_mat(k)%L)

          !Calculate estimative of condition number
          rbf_mat(k)%condnum=condnumest(rbf_mat(k)%L, rbf_mat(k)%n)

          !print fmt, transpose(rbf_mat(k)%L)
          !print*
          !print *, "Cond num:", rbf_mat(k)%condnum
          !print*
          !print*,"------------------------------------------------"
          !print*

       end do

    case default
       print*, "RBFVEC_MATRIX_BUILD ERROR : Unknown RBF stencil. ", stencil
       stop
    end select

    return
  end subroutine rbfvec_matrix_build

  subroutine rbfvec_addpol(rbf_mat_local, mesh)
    !--------------------------------------------------------
    ! rbfvec_addpol
    ! Adds a polynomial to a RBF vector reconstruction
    ! Must receive rbf%mat_local%w with weights (A^(-1)b)
    !--------------------------------------------------------
    !RBF variable structure
    type(rbf_matrix_structure), intent(inout) :: rbf_mat_local

    !Mesh structure
    type(grid_structure), intent(in) :: mesh

    !Vector of normal components
    real (r8), allocatable :: c1(:)
    real (r8), allocatable :: c2(:)
    real (r8), allocatable :: c3(:)

    !Auxiliar vectors
    real (r8), allocatable :: w1(:)
    real (r8), allocatable :: w2(:)
    real (r8), allocatable :: w3(:)

    !Auxiliar matrix
    real (r8), dimension(1:3) :: l1, l2, l3, d

    !Indexes and aux vars
    integer (i4):: i, n

    n=rbf_mat_local%n

    !Initialize polynomial vector
    if(.not.allocated(rbf_mat_local%pol))then
       allocate(rbf_mat_local%pol(1:3))
    end if
    rbf_mat_local%pol=0

    allocate(c1(1:n))
    allocate(c2(1:n))
    allocate(c3(1:n))
    allocate(w1(1:n))
    allocate(w2(1:n))
    allocate(w3(1:n))

    !Build problem
    do i=1,rbf_mat_local%n
       !Vector with information of normal vectors
       c1(i)=mesh%edhx(rbf_mat_local%pts(i))%nr(1)
       c2(i)=mesh%edhx(rbf_mat_local%pts(i))%nr(2)
       c3(i)=mesh%edhx(rbf_mat_local%pts(i))%nr(3)

       !d=arclen(p,mesh%edhx(rbf_mat(k)%pts(i))%c%p)
       !d=norm(p-mesh%v(rbf_mat(k)%pts(i))%p)
    end do

    call choleskysolve(rbf_mat_local%L, w1, c1, n)
    call choleskysolve(rbf_mat_local%L, w2, c2, n)
    call choleskysolve(rbf_mat_local%L, w3, c3, n)
    !print*, "c"
    !print*, c1
    !print*, c2
    !print*, c3
    !print*, "w"
    !print*, w1
    !print*, w2
    !print*, w3

    l1(1)=dot_product(c1, w1)
    l1(2)=dot_product(c1, w2)
    l1(3)=dot_product(c1, w3)
    l2(1)=dot_product(c2, w1)
    l2(2)=dot_product(c2, w2)
    l2(3)=dot_product(c2, w3)
    l3(1)=dot_product(c3, w1)
    l3(2)=dot_product(c3, w2)
    l3(3)=dot_product(c3, w3)

    !print*, "A"
    !print*, l1
    !print*, l2
    !print*, l3

    d(1)=dot_product(c1, rbf_mat_local%w)
    d(2)=dot_product(c2, rbf_mat_local%w)
    d(3)=dot_product(c3, rbf_mat_local%w)

    rbf_mat_local%pol=solve3x3(l1, l2, l3, d)
    !print*, "pol"
    !print*, rbf_mat_local%pol

    !print*, "Old weights"
    !print*, rbf_mat_local%w
    rbf_mat_local%w=rbf_mat_local%w &
         -rbf_mat_local%pol(1)*w1&
         -rbf_mat_local%pol(2)*w2&
         -rbf_mat_local%pol(3)*w3

    !print*, "New weights"
    !print*, rbf_mat_local%w
    !print*
    !print*
    !read*,i
    return
  end subroutine rbfvec_addpol


  function vecrecon_rbf(p, var, stencil, rbf_mat, mesh)
    !--------------------------------------------------------
    ! vecrecon_rbf
    ! Radial bais function vector reconstruction
    !--------------------------------------------------------
    !Reconstruction point
    real (r8), intent(in) :: p(1:3)

    !Variable with normal vector values
    type(scalar_field), intent(in) :: var

    !Stencil or type of rbf point structure
    ! HX   = Hexagonal, pentagonal edge midpoints
    character(len=4), intent(in) :: stencil

    !RBF matrix vector
    type(rbf_matrix_structure), intent(inout) :: rbf_mat(:)

    !Mesh structure
    type(grid_structure), intent(in) :: mesh

    !Interpoled value
    real (r8):: vecrecon_rbf(1:3)

    !Aux var
    real (r8):: d

    !Indexes and aux vars
    integer (i4):: i
    integer (i4):: k

    !Do reconstruction
    if(stencil(1:2)=='HX'.or. stencil(1:2)=='hx')then
       !Do reconstruction assuming point in hexagon
       !Get nearest node to point
       k=getnearnode(p, mesh)
    elseif(trim(stencil(1:2))=='TR'.or.trim(stencil(1:2))=='tr' .or. &
         trim(stencil(1:3))=='ETR'.or.trim(stencil(1:3))=='etr')then
       !Do reconstruction assuming point in triangle
       !Get triangle
       k=gettr(p, mesh)
    else
       print*, "RBF_INTERPOL ERROR : Unknown RBF stencil. ", stencil
       stop
    end if

    !Save RHS of the system
    do i=1,rbf_mat(k)%n
       rbf_mat(k)%b(i)=var%f(rbf_mat(k)%pts(i))
    end do

    !Calculate weights with no polynomial
    call choleskysolve(rbf_mat(k)%L, rbf_mat(k)%w, rbf_mat(k)%b, rbf_mat(k)%n)

    if(trim(stencil(3:4))=="PC".or.trim(stencil(3:4))=="pc".or.&
         trim(stencil(4:4))=="P".or.trim(stencil(4:4))=="p")then
       !Add polinomial to reconstruction
       call rbfvec_addpol(rbf_mat(k), mesh)
    end if

    !Do reconstruction
    vecrecon_rbf=(/0._r8,0._r8,0._r8/)
    !print*, "RBF"
    !print*, " TR V w  d rbf(d)"

    do i=1,rbf_mat(k)%n
       !print*, i, rbf_mat(k)%n
       !Calculate distance
       d=arclen(p,mesh%edhx(rbf_mat(k)%pts(i))%c%p)
       !d=norm(p-mesh%v(rbf_mat(k)%pts(i))%p)

       !Apply weights
       vecrecon_rbf=vecrecon_rbf+rbf_mat(k)%w(i)*rbf(d, rbf_mat(k)%h)* &
            mesh%edhx(rbf_mat(k)%pts(i))%nr
       !print*, k, rbf_mat(k)%pts(i), rbf_mat(k)%w(i),  d, rbf(d, rbf_mat(k)%h)
    end do
    if(trim(stencil(3:4))=="PC".or.trim(stencil(3:4))=="pc".or.&
         trim(stencil(4:4))=="P".or.trim(stencil(4:4))=="p")then
       !Add polinomial to reconstruction
       vecrecon_rbf=vecrecon_rbf+rbf_mat(k)%pol
    end if

    !print*, "Interpolated value:", vecrecon_rbf
    return

  end function vecrecon_rbf

  function scinterpol_rbf(p, var, mesh, stencil, rbf_mat)
    !--------------------------------------------------------
    ! RBF_INTERPOL
    ! Radial bais function interpolation
    !--------------------------------------------------------
    !Interpolation point
    real (r8), intent(in) :: p(1:3)

    !Variable for interpolation procedure
    type(scalar_field), intent(in) :: var

    !Stencil or type of rbf point structure
    ! TR   = Points on triangles - 3 points only
    ! EXTR = (Extended TR) Points neighbour to triangle also included
    ! HX   = Hexagonal, pentagonal neighbour points
    character(len=4), intent(in) :: stencil

    !RBF matrix vector
    type(rbf_matrix_structure), intent(inout) :: rbf_mat(:)

    !Mesh structure
    type(grid_structure), intent(in) :: mesh

    !Interpoled value
    real (r8):: scinterpol_rbf

    !Aux var
    real (r8):: d

    !Indexes and aux vars
    integer (i4):: i
    integer (i4):: k

    !Do interpolation
    if(trim(stencil(1:2))=='TR' .or. trim(stencil(1:3))=='ETR' .or. &
         trim(stencil(1:2))=='tr' .or. trim(stencil(1:3))=='etr')then
       !Do interpolation assuming triangle as center
       !Get triangle that point belongs
       k=gettr(p, mesh)

    elseif(trim(stencil(1:2))=='HX'.or.trim(stencil(1:2))=='hx')then
       !Do interpolation assuming node as center
       !Get nearest node to point
       k=getnearnode(p, mesh)
    else
       print*, "scinterpol_rbf ERROR : Unknown RBF stencil. ", stencil
       stop
    end if

    !Save RHS of the system
    do i=1,rbf_mat(k)%n
       rbf_mat(k)%b(i)=var%f(rbf_mat(k)%pts(i))
    end do

    !Calculate weights for this triangle's interpolation
    call choleskysolve(rbf_mat(k)%L, rbf_mat(k)%w, rbf_mat(k)%b, rbf_mat(k)%n)

    if(trim(stencil(3:4))=="PC".or.trim(stencil(3:4))=="pc" .or. &
         trim(stencil(4:4))=="P".or.trim(stencil(4:4))=="p")then
       !Add constant polinomial to reconstruction
       call rbf_addpol(rbf_mat(k))
    end if

    !Do interpolation
    scinterpol_rbf=0._r8
    !print*, "RBF"
    !print*, " TR V w  d rbf(d)"
    do i=1,rbf_mat(k)%n
       !Calculate distance
       d=arclen(p,mesh%v(rbf_mat(k)%pts(i))%p)
       !d=norm(p-mesh%v(rbf_mat(k)%pts(i))%p)

       !Apply weights
       scinterpol_rbf=scinterpol_rbf+rbf_mat(k)%w(i)*rbf(d, rbf_mat(k)%h)
       !print*, k, rbf_mat(k)%pts(i), rbf_mat(k)%w(i),  d, rbf(d, rbf_mat(k)%h)
    end do

    if(trim(stencil(3:4))=="PC".or.trim(stencil(3:4))=="pc" .or. &
         trim(stencil(4:4))=="P".or.trim(stencil(4:4))=="p")then
       !Add polinomial to reconstruction
       scinterpol_rbf=scinterpol_rbf+rbf_mat(k)%pol(1)
    end if

    !print*, "Interpoled value:", rbf_interpol
    return

  end function scinterpol_rbf

  !=====================================================================================
  !    SHEPARDS METHOD INTERPOLATIONS
  !=====================================================================================

  function scinterpol_lmshep(p, var, mesh, r)
    !--------------------------------------------------------
    ! LOCAL MODIFIED SHEPARD INTERPOLATION
    !
    ! Uses the inverse of the squared distanced
    ! See Reanka 1988 - Multivar Interpol of Large...
    !--------------------------------------------------------
    !Interpolation point
    real (r8), intent(in) :: p(1:3)

    !Variable for interpolation procedure
    type(scalar_field), intent(in) :: var

    !Mesh structure
    type(grid_structure), intent(in) :: mesh

    !Interpoled value
    real (r8):: scinterpol_lmshep

    !Radius of spherical cap
    real (r8), intent(in) :: r

    !index
    integer:: i
    integer:: n

    !Weights
    real (r8), allocatable :: w(:)
    real (r8):: sumw

    !Lists of nodes/trs and distances
    integer (i4), allocatable :: list(:)
    real (r8), allocatable :: listd(:)

    !Find nodes that are in the spherical cap
    n=0
    call getnearnodes(p, mesh, r, list, listd, n)

    if(n<=0)then
       print*, "LMSHEP_INTERP ERROR: No nodes/trs in cap of radius r=", r
       stop
    end if
    !print*, n
    !do i=1, n
    !    print*, list(i), listd(i)
    !end do

    if(listd(1)<eps*10)then
       !Very near to a node, consider the value of the node
       !   and return
       scinterpol_lmshep=var%f(list(1))
       return
    end if

    !Calculate Weights
    allocate(w(1:n))
    sumw=0._r8
    do i=1,n
       w(i)=(positive(r-listd(i))/(r*listd(i)))**2
       sumw=sumw+w(i)
    end do
    !print*
    !Do interpolation
    scinterpol_lmshep=0._r8
    do i=1,n
       scinterpol_lmshep=scinterpol_lmshep+w(i)*var%f(list(i))/sumw
       !   print*, scinterpol_lmshep, list(i), w(i)/sumw, var%f(list(i))
    end do
    !print*
    !print*, scinterpol_lmshep

    return
  end function scinterpol_lmshep

  function scinterpol_lmshep_edhx(p, var, mesh, r)
    !--------------------------------------------------------
    ! LOCAL MODIFIED SHEPARD INTERPOLATION FOR VALUES
    !   ON HEXGANAL EDGES MIDPOINTS
    !
    ! Uses the inverse of the squared distanced
    !--------------------------------------------------------
    !Interpolation point
    real (r8), intent(in) :: p(1:3)

    !Variable for interpolation procedure
    type(scalar_field), intent(in) :: var

    !Mesh structure
    type(grid_structure), intent(in) :: mesh

    !Interpoled value
    real (r8):: scinterpol_lmshep_edhx

    !Radius of spherical cap
    real (r8), intent(in) :: r

    !index
    integer:: i
    integer:: n

    !Weights
    real (r8), allocatable :: w(:)
    real (r8):: sumw

    !Lists of nodes/trs and distances
    integer (i4), allocatable :: list(:)
    real (r8), allocatable :: listd(:)

    !Find nodes that are in the spherical cap
    n=0
    !call getnearnodes(p, mesh, r, list, listd, n)
    call getnearhxedges(p, mesh, r, list, listd, n)

    if(n<=0)then
       print*, "LMSHEP_EDHX_INTERPOL ERROR: No nodes/trs in cap of radius r=", r
       stop
    end if
    !print*, n
    !do i=1, n
    !    print*, list(i), listd(i)
    !end do

    if(listd(1)<eps*10)then
       !Very near to a node, consider the value of the node
       !   and return
       scinterpol_lmshep_edhx=var%f(list(1))
       return
    end if

    !Calculate Weights
    allocate(w(1:n))
    sumw=0._r8
    do i=1,n
       w(i)=(positive(r-listd(i))/(r*listd(i)))**2
       sumw=sumw+w(i)
    end do
    !print*
    !Do interpolation
    scinterpol_lmshep_edhx=0._r8
    do i=1,n
       scinterpol_lmshep_edhx=scinterpol_lmshep_edhx+w(i)*var%f(list(i))/sumw
       !   print*, scinterpol_lmshep_edhx, list(i), w(i)/sumw, var%f(list(i))
    end do
    !print*
    !print*, scinterpol_lmshep_edhx

    return
  end function scinterpol_lmshep_edhx


  !function scinterpol_qdshep(p, var, r, mesh)
  function scinterpol_qdshep(p, var, mesh)
    !--------------------------------------------------------
    ! QUADRATIC SHEPARD INTERPOLATION
    !
    ! Uses the inverse of the squared distanced
    !--------------------------------------------------------
    !Interpolation point
    real (r8), intent(in) :: p(1:3)

    !Variable for interpolation procedure
    type(scalar_field), intent(inout) :: var

    !Mesh structure
    type(grid_structure), intent(in) :: mesh

    !Interpoled value
    real (r8):: scinterpol_qdshep

    !Radius of spherical cap
    !real (r8), intent(in) :: r
    real (r8):: r

    !index
    integer:: i
    integer:: n

    !Weights
    real (r8), allocatable :: w(:)
    real (r8):: sumw

    !Lists of nodes/trs and distances
    integer (i4), allocatable :: list(:)
    real (r8), allocatable :: listd(:)

    !Find nodes that are in the spherical cap
    n=0
    r=mesh%maxvdist*1.05_r8
    call getnearnodes(p, mesh, r, list, listd, n)

    if(n<=0)then
       print*, "QDSHEP_INTERP ERROR: No nodes/trs in cap of radius r=", r
       stop
    end if
    !print*, "Nearnodes: (", n, ") index and distance"
    !do i=1, n
    !    print*, list(i), listd(i)
    !end do

    if(listd(1)<eps*10)then
       !Very near to a node, consider the value of the node
       !   and return
       scinterpol_qdshep=var%f(list(1))
       return
    end if

    !Calculate Weights
    allocate(w(1:n))
    sumw=0._r8
    do i=1,n
       w(i)=(positive(r-listd(i))/(r*listd(i)))**2
       sumw=sumw+w(i)
    end do

    !Do interpolation
    scinterpol_qdshep=0._r8
    do i=1,n
       scinterpol_qdshep=scinterpol_qdshep+w(i)*pol2lsq_eval_trv(p, var, mesh, list(i)) /sumw
       !print*, scinterpol_qdshep, list(i), w(i)/sumw, var%f(list(i)), pol2_eval_trv(p, list(i), mesh, var)
    end do
    !print*
    !print*, scinterpol_qdshep

    return
  end function scinterpol_qdshep

  !=====================================================================================
  !    NATURAL NEIGHBOUR COORDS INTERPOLATION
  !=====================================================================================

  function scinterpol_natneib(p, var, mesh, kinterp)
    !--------------------------------------------------------
    ! NATNEIB_INTERPOL
    !
    ! Calculates the natural neighbour interpolation based on
    !  natural neighbour coordinates given by kinterp = SIB or LAP
    !  standing respectively for Sibson and Laplace coordinates.
    !--------------------------------------------------------
    !Interpolation point
    real (r8), intent(in) :: p(1:3)

    !Variable for interpolation procedure
    type(scalar_field), intent(inout) :: var

    !Mesh structure
    type(grid_structure), intent(in) :: mesh

    !Kind of interpolation
    ! "LAP" => Laplacian coordinates
    ! "SIB" => Sibson C1 interpolation
    character (len=4), optional :: kinterp
    character (len=4):: kindinterpol

    !Interpoled value
    real (r8):: scinterpol_natneib

    !Natural coords
    type (general_coords):: lap
    type (general_coords):: sib

    !index
    integer:: i

    !Only interpolates for scalar variables
    ! with values on triangle vertices (voronoi centers)
    if(var%pos /= 0 )then
       print*, "NATNEIB_INTERPOL ERROR: Variable not defined on triangle vertices"
       stop
    end if

    !Set default interpolation kind as Sibson Coordinate interpolation
    if(.not.present(kinterp) )then
       kindinterpol="SIB"
    else
       kindinterpol=kinterp
    end if

    !Calculate coordinates
    call natural_coords(p, lap, sib, mesh)

    scinterpol_natneib=0._r8
    if(trim(kindinterpol)=="LAP"  .or. trim(kindinterpol)=="lap")then
       !Use Laplace
       do i=1, lap%n
          scinterpol_natneib=scinterpol_natneib+lap%w(i)*var%f(lap%v(i))
       end do
    elseif(trim(kindinterpol)=="FAR"  .or. trim(kindinterpol)=="far")then
       !Use Farin
       scinterpol_natneib=farin()
    else
       !Use Sibson
       do i=1, sib%n
          scinterpol_natneib=scinterpol_natneib+sib%w(i)*var%f(sib%v(i))
       end do
    end if

  contains

    function farin()
      !Calculates the Farin interpolant
      ! See Hiyoshi and Sugihara 2004

      real (r8):: farin
      real (r8):: z(1:sib%n , 1:sib%n)
      real (r8):: f
      real (r8):: v(1:3)
      real (r8):: sig
      real (r8):: sumf
      real (r8):: sumf2
      integer:: i
      integer:: j
      integer:: k

      !Check if gradients exist
      if(.not.allocated(var%g))then
         call gradcalc(var, mesh)
      end if

      !Calculate edge gradients
      farin=0._r8
      do i=1,sib%n
         do j=1, sib%n
            if(i/=j)then
               v=proj_vec_sphere(mesh%v(sib%v(j))%p-mesh%v(sib%v(i))%p, &
                    mesh%v(sib%v(i))%p)
               z(i,j)=dot_product(var%g(sib%v(i))%v, v)
            else
               z(i,j)=var%f(sib%v(i))
            endif
         end do
      end do

      !Compose polinomial
      farin=0._r8
      do i=1,sib%n
         do j=1, sib%n
            do k=1, sib%n
               sig=sib%w(i)*sib%w(j)*sib%w(k)
               if(i==j .and. j==k)then
                  !i=j=k
                  f=z(i,i)
               elseif(i==j .or. j==k .or. i==k)then
                  !2 index are equal, one diferent
                  if(i==j) f=z(i,i)+z(i,k)/3._r8
                  if(j==k) f=z(j,j)+z(j,i)/3._r8
                  if(i==k) f=z(k,k)+z(k,j)/3._r8
               else !i, j, k all diferent
                  sumf=(z(i,i)+z(j,j)+z(k,k))/3._r8
                  sumf2=(z(i,j)+z(i,k)+z(j,i)+z(j,k)+z(k,i)+z(k,j))/12._r8
                  f=sumf+sumf2
               endif
               farin=farin+f*sig
            end do
         end do
      end do
      return
    end function farin

  end function scinterpol_natneib

  subroutine natural_coords(p, lap, sib, mesh)
    !--------------------------------------------------------
    ! NATURAL_COORDS
    !
    ! Calculates the natural neighbour coordinate (laplacian and sibson)
    !  for p in a given mesh.
    ! See Bertram, Bobach, Umlauf 2006
    !--------------------------------------------------------
    !Interpolation point
    real (r8), intent(in) :: p(1:3)

    !Natural coords
    type (general_coords), intent(out) :: lap
    type (general_coords), intent(out) :: sib

    !Mesh structure
    type(grid_structure), intent(in) :: mesh

    !Coordinate auxiliar variables
    integer (i4):: node
    real (r8), allocatable :: laptmp(:)
    real (r8), allocatable :: sibtmp(:)
    real (r8):: sumlap
    real (r8):: sumsib

    !List of nodes to contribute on new voronoi region
    integer (i4), allocatable :: nodes(:)
    integer (i4):: nnodes

    !Index
    integer (i4):: i

    !Get nearest vertice node to point
    node = getnearnode(p, mesh)

    !Allocate temporary space
    allocate(nodes(1:2*mesh%v(node)%nnb))
    allocate(laptmp(1:2*mesh%v(node)%nnb))
    allocate(sibtmp(1:2*mesh%v(node)%nnb))

    !Calculate the contribuition of nearest node to new voronoi region
    nodes(1)=node
    nnodes=1
    call local_natural_coords(p, node, laptmp(1), sibtmp(1), nodes, nnodes, mesh)

    !Calculate the coordinate relative to each
    ! neighbour of 'node' in nodes(:) list
    do i=2, 2*mesh%v(node)%nnb
       if(nnodes==1)exit

       !Get coordinates for the point relative to nodes(i)
       call local_natural_coords(p, nodes(i), laptmp(i), sibtmp(i), nodes, nnodes, mesh)

       if(i>=nnodes)exit
    end do

    !Check for error
    if(i>=2*mesh%v(node)%nnb)then
       print*, "natural_coords Warning : too many neighbours!"
    end if

    !Allocate space for output
    lap%n=nnodes
    lap%pos=0
    allocate(lap%v(1:lap%n))
    allocate(lap%w(1:lap%n))

    sib%n=nnodes
    sib%pos=0
    allocate(sib%v(1:sib%n))
    allocate(sib%w(1:sib%n))

    !Save coordinates for output
    lap%v(1:nnodes)=nodes(1:nnodes)
    sumlap=sum(laptmp(1:nnodes))
    lap%w(1:nnodes)=laptmp(1:nnodes)/sumlap

    sib%v(1:nnodes)=nodes(1:nnodes)
    sumsib=sum(sibtmp(1:nnodes))
    sib%w(1:nnodes)=sibtmp(1:nnodes)/sumsib

    !Check if coordinates make sense
    do i=2, nnodes
       if( sib%w(i) > sib%w(1)*1.5_r8)then
          print*, "NATURAL_COORDS WARNING : Sibson coordinate problem"
       end if
    end do

    return

  contains

    subroutine local_natural_coords(p, node, lap, sib, nodes, nnodes, mesh)
      !--------------------------------------------------------
      ! LOCAL_NATURAL_COORDS
      !
      ! Calculates the contribuition of a certain 'node', that is one of the
      ! nearest nodes to the point 'p', to the
      ! natural neighbour coordinates of a point 'p' relative to 'mesh'
      ! Returns the laplacian and sibson coordinates
      ! Returns a list of neighbour nodes that may contribute to the
      !  coordinates (add to list)
      !--------------------------------------------------------
      !Interpolation point
      real (r8), dimension(1:3) , intent(in) :: p

      !Mesh node index for which the coordinates refer
      integer (i4), intent(in) :: node

      !Natural coords for this node
      real (r8), intent(out) :: lap
      real (r8), intent(out) :: sib

      !List of neighbour that may contribute to coordinates
      integer (i4), intent(inout) :: nodes(:)

      !Size of list
      integer (i4), intent(inout) :: nnodes

      !Mesh structure
      type(grid_structure), intent(in) :: mesh

      !Cartesian coordinate of points and other points/vectors
      real (r8), dimension(1:3) :: pb
      real (r8), dimension(1:3) :: p1
      real (r8), dimension(1:3) :: p2
      real (r8), dimension(1:3) :: normal

      !Intersection point
      real (r8), dimension(1:3) :: r
      real (r8), dimension(1:3) :: r1
      real (r8), dimension(1:3) :: r2
      real (r8):: normr

      !Intersection edges
      integer (i4):: ed
      integer (i4):: ed1
      integer (i4):: ed2
      integer (i4):: tr1
      integer (i4):: tr2
      integer (i4):: newnode

      !Auxiliar vars
      integer (i4):: i
      integer (i4):: j

      !Polygon of new voronoi region
      type(vector), allocatable :: pol1(:)
      type(vector), allocatable :: pol2(:)

      !Polygon indexes
      integer (i4):: ipol1
      integer (i4):: npol1
      integer (i4):: ipol2
      integer (i4):: npol2
      integer (i4):: nintersec
      integer (i4):: ispol1
      integer (i4):: ispol2

      !Point to node distance
      real (r8):: dpnode

      !Node addition flag
      logical:: add

      lap=0._r8
      sib=0._r8

      !Check if point is too close to the nearest node
      dpnode=arclen(p, mesh%v(node)%p)
      if( dpnode < eps )then
         lap=1._r8
         sib=1._r8
         return
      end if

      !Get bisection point and normal to bisectrix
      call ortogonalarc(p, mesh%v(node)%p, pb, normal)

      !Loop through edges to define intersection points
      nintersec=0
      do i=1, mesh%v(node)%nnb
         !For each edge, get its extremal points
         !  in counter-clockwise order
         tr1=mesh%v(node)%tr(modint(i-1,mesh%v(node)%nnb))
         p1=mesh%tr(tr1)%c%p

         tr2=mesh%v(node)%tr(i)
         p2=mesh%tr(tr2)%c%p

         !Test for intersection between bisectrix and edge
         r=gcircarcintersec(normal, p1, p2)
         normr=norm(r)
         !Check results of intersections
         if(normr > 0.5_r8 .and. normr<2) then
            !Only one point of intersection
            ed=i
            !Add one to intersection counter
            nintersec=nintersec+1
            if(nintersec==1)then
               r1=r
               ed1=ed
            elseif(nintersec>=2)then
               !First check if the point is not the one already found
               !That is, if the point is a vertex of the cell
               if(norm(r-r1)<eps*100)then
                  !Store this point in second position
                  r2=r
                  ed2=ed
                  !Check if found a diferent intersection point
               elseif(norm(r-r1)>=eps*100)then
                  !Save the point and it edge
                  r2=r
                  ed2=ed
                  exit
               end if
            end if
         elseif( normr > 2 ) then
            !There are infinite number of intersections
            !The bisectrix is one of the edges (arc in great circle)
            print*, "NATURALCOORDS ERROR: Bisectrix on voronoi edge in nearest node"
            print*, "p      pb     node"
            print*, p, pb, mesh%v(node)%p
            stop
         end if
      end do

      if(nintersec<2 .or. norm(r1-r2) <eps)then
         !print*, "NATURALCOORDSFIRST Warning : not intersecting 2 points"
         lap=0._r8
         sib=0._r8
         return
      end if

      !  Laplacian coordinate
      lap = arclen(r1, r2)/dpnode

      !Check if edges are correct
      if(ed2-ed1<= 0)then
         print*, "NATURALCOORDS ERROR: Intersecting edges incorrect"
         print*, "node, ed1, ed2"
         print*, node, ed1, ed2
         stop
      end if

      !Polygon 1 goes from ed1 to ed2
      npol1=ed2-ed1+2

      !Allocate space for polygons
      allocate(pol1(1:npol1))

      !Set polygon 1 points and check if it is the nearest polygon to 'p'
      ipol1=1
      pol1(ipol1)%v=r1
      ispol1=0
      do i=ed1, ed2-1
         !Check if point in nearer to p then to node
         ipol1=ipol1+1
         pol1(ipol1)%v=mesh%tr(mesh%v(node)%tr(i))%c%p
         if( arclen(pol1(ipol1)%v, p) < arclen (pol1(ipol1)%v, mesh%v(node)%p)-eps )then
            ispol1=ispol1+1
         end if
      end do
      ipol1=ipol1+1
      pol1(ipol1)%v=r2

      if(ispol1>0)then
         !Pol1 is the nearest polygon to point p

         !Add nodes to list
         do i=ed1, ed2

            !Find a possible node
            newnode=getothernodeonedge(mesh%v(node)%ed(i), node, mesh)

            !Check if it is on the list
            add=.true.
            do j=1,nnodes
               if(newnode==nodes(j))then
                  add=.false.
               end if
            end do

            !Add node
            if(add)then
               nnodes=nnodes+1
               nodes(nnodes)=newnode
            end if
         end do

         !Calculate the polygon area
         sib = sphpolarea(pol1, npol1)

         return

      end if


      !Polygon 2 goes from ed2 to nnb and then from 1 to ed1
      npol2=(mesh%v(node)%nnb-ed2 + 1) + ed1 +1

      !Allocate space for polygons
      allocate(pol2(1:npol2))

      !Set polygon 2 points and check if it is the nearest polygon to 'p'
      ipol2=1
      pol2(ipol2)%v=r2
      ispol2=0
      do i=ed2, mesh%v(node)%nnb
         ipol2=ipol2+1
         pol2(ipol2)%v=mesh%tr(mesh%v(node)%tr(i))%c%p
         if( arclen(pol2(ipol2)%v, p) < arclen (pol2(ipol2)%v, mesh%v(node)%p)-eps )then
            ispol2=ispol2+1
         end if
      end do
      do i=1, ed1-1
         ipol2=ipol2+1
         pol2(ipol2)%v=mesh%tr(mesh%v(node)%tr(i))%c%p
         if( arclen(pol2(ipol2)%v, p) < arclen (pol2(ipol2)%v, mesh%v(node)%p)-eps )then
            ispol2=ispol2+1
         end if
      end do
      ipol2=ipol2+1
      pol2(ipol2)%v=r1


      if (ispol2 > 0) then
         !Pol2 is the nearest polygon to point p

         !Add nodes to list
         do i=ed2, mesh%v(node)%nnb

            !Find a possible node
            newnode=getothernodeonedge(mesh%v(node)%ed(i), node, mesh)

            !Check if it is on the list
            add=.true.
            do j=1,nnodes
               if(newnode==nodes(j))then
                  add=.false.
               end if
            end do

            !Add node
            if(add)then
               nnodes=nnodes+1
               nodes(nnodes)=newnode
            end if
         end do
         do i=1, ed1

            !Find a possible node
            newnode=getothernodeonedge(mesh%v(node)%ed(i), node, mesh)

            !Check if it is on the list
            add=.true.
            do j=1,nnodes
               if(newnode==nodes(j))then
                  add=.false.
               end if
            end do

            !Add node
            if(add)then
               nnodes=nnodes+1
               nodes(nnodes)=newnode
            end if
         end do
         !Calculate the polygon area
         sib = sphpolarea(pol2, npol2)

         return
      end if

      return
    end subroutine local_natural_coords

  end subroutine natural_coords


  !=====================================================================================
  !   LEAST SQUARES POLINOMIAL FIT
  !=====================================================================================

  function scinterpol_lsqhx_trv (p, var, mesh)
    !-----------------------------------------------------------
    !  scinterpol_lsqhx_trv
    !
    !   Given the triangulation of a set of nodes on the unit
    ! sphere (mesh), along with data values at the nodes/vertices
    ! (in var), computes the value at a point (p)
    ! of a quadratic least squares fit.
    !
    ! On input:
    !       p         =  point in cart coords
    !       mesh      =  mesh structure
    !       var       =  interpolable_variable structure
    !                      with var%pos=0
    ! Returns the value at p
    !-------------------------------------------------------
    !Interpolation point
    real (r8), intent(in) :: p(1:3)

    !Variable for interpolation procedure
    type(scalar_field) :: var

    !Mesh in which the variable belongs
    type(grid_structure), intent(in) :: mesh

    !Value at interpolated point
    real (r8):: scinterpol_lsqhx_trv

    !Nearest node to p
    integer(i4):: node

    ! Locate p with respect to the triangulation.
    node=getnearnode(p, mesh)

    !Aply least squares approximation relative to node
    scinterpol_lsqhx_trv=pol2lsq_eval_trv (p, var, mesh, node)

    return
  end function scinterpol_lsqhx_trv


  subroutine precalcgradpol(var, mesh, kindinterpol)
    !--------------------------------------------------------------
    !   Precalcgradpol
    !
    ! Precalculates gradients or polinomial for interpolation if needed
    !---------------------------------------------------------------
    !Mesh
    type(grid_structure), intent(in) :: mesh

    !Scalar field
    type(scalar_field), intent(inout) :: var

    !Kind of interpolation
    character (len=8), optional:: kindinterpol
    character (len=8):: kindinterp

    if(present(kindinterpol))then
       kindinterp=kindinterpol
    else
       kindinterp="lsqhx"
    end if

    if(var%pos/=0)then
       print*, "Warning on precalcgradpol: "//&
            "   Variable not on tr vertices"
       return
    end if

    select case(trim(kindinterp))
    case('hermtrv') !Hermite - needs gradient as well
       call gradcalc(var, mesh)
    case default
       call pol2lsq_global(var, mesh)
    end select

    return
  end subroutine precalcgradpol

  subroutine gradcalc(var, mesh)
    !---------------------------------------------------------------------
    ! GRADCALC
    !
    !   Recieves an interpolable variable and a mesh
    !   Returns estimatives of the gradient vector on the entire mesh
    !---------------------------------------------------------------------
    type(scalar_field), intent(inout) :: var
    type(grid_structure), intent(in) :: mesh

    !Auxiliar vars
    real (r8) :: dx, dy
    integer (i4):: i

    !Allocate space if necessary
    if(.not.allocated(var%g))then
       allocate(var%g(1:var%n))
    end if

    !Calculate the gradient for each triangle node
    select case(var%pos)
    case(0) !Gradients on the triangle vertices
       !Calculate the polinomials
       call pol2lsq_global (var, mesh)
       !Calculate the gradients
       do i= 1, mesh%nv
          !Optimized SSRFPACK XY Least Squares
          !call gradxyls_trv (i, var, mesh)

          !Get gradient on the projected plane
          dx=var%pol(i)%c(4)
          dy=var%pol(i)%c(5)

          ! Rotate the gradient (DX,DY,0) back into the original
          !   coordinate system.
          call aplyrt (dx, dy, &
               var%pol(i)%cx, var%pol(i)%sx, var%pol(i)%cy, var%pol(i)%sy, &
               var%g(i)%v)
       end do

       !case(1) !Gradients on Triangle circumcenters
    case default
       print*, "Error on gradcalc: Var not on tr vertices"
       stop
    end select

    return
  end subroutine gradcalc

  function pol2lsq_eval_trv (p, var, mesh, node)
    !----------------------------------------------------------
    ! CALCULATES POLINOMIAL 2ND ORDER LEAST SQUARE FIT FOR A POINT
    !
    ! Uses polinomial coeficients estimation based on a quadratic least squares
    !  on the tangent plane (XY) for triangle vertices (nodes) to
    !  calculate aproximation for a given point p relative to the
    !  polinomial of the 'node'
    !
    !----------------------------------------------------------------
    !Approximation point
    real (r8), intent(in) :: p(1:3)
    real (r8):: q(1:3)
    real (r8):: x
    real (r8):: y

    !Polinomial for node
    integer (i4), intent(in) :: node

    !Mesh
    type(grid_structure), intent(in) :: mesh

    !Variable
    type(scalar_field), intent(inout) :: var

    !Returning value of approximation
    real (r8):: pol2lsq_eval_trv

    !Aux
    integer (i4):: i

    !If coeficients are not wet calculated, do so
    if(.not.allocated(var%pol))then
       allocate(var%pol(1:var%n))
       do i=1, var%n
          allocate(var%pol(i)%c(1:5))
       end do
       call pol2lsq_global (var, mesh)
    end if

    !Rotate point to pole (p -> q)
    call aplyr (p(1), p(2), p(3), &
         var%pol(node)%cx, var%pol(node)%sx, var%pol(node)%cy, var%pol(node)%sy, &
         q(1), q(2), q(3))
    x=q(1)
    y=q(2)

    pol2lsq_eval_trv=0._r8
    pol2lsq_eval_trv=x**2*var%pol(node)%c(1)+ &
         x*y*var%pol(node)%c(2)+ &
         y**2*var%pol(node)%c(3)+ &
         x*var%pol(node)%c(4)+ &
         y*var%pol(node)%c(5)+ &
         var%f(node)
    return
  end function pol2lsq_eval_trv

  subroutine pol2lsq_global (var, mesh)
    !---------------------------------------------------------------------
    !	POLINOM2LSFIT
    !
    !   Recieves an interpolable variable and a mesh
    !   Calculates the coeficients of a polinomial fit of degree 2 for each
    !    node in the mesh
    !---------------------------------------------------------------------
    type(scalar_field), intent(inout) :: var
    type(grid_structure), intent(in) :: mesh
    integer (i4):: i

    !Allocate space if necessary
    if(.not.allocated(var%pol))then
       allocate(var%pol(1:var%n))
       do i=1, var%n
          allocate(var%pol(i)%c(1:5))
       end do
    end if

    !Calculate the gradient for each triangle node
    select case(var%pos)
    case(0) !Variable on the triangle vertices
       do i= 1, mesh%nv
          !Fit polinomial for this node
          call pol2lsq_node (i, var, mesh)
       end do
    case(1) !Variable on Triangle circumcenters
       print*, "Error on polinom2lsfit: Only variables on triangle vertices implemented"
       stop
    end select

  end subroutine pol2lsq_global

  subroutine pol2lsq_node (node, var, mesh)
    !----------------------------------------------------------
    ! POLINOMIAL 2nd ORDER LEAST SQUARE FIT
    !
    !  Polinomial coeficients estimation based on a quadratic least squares
    !  on the tangent plane (XY) for triangle vertices (nodes)
    !
    !   Details:
    ! Given a triangulation of a set of nodes on the unit
    ! sphere (mesh) with their associated data values var%f,
    ! estimates a gradient vector at 'node'  as follows:
    ! 1) the coordinate system is rotated so that 'node'
    ! becomes the north pole
    ! 2) 'node' and a set of nearby nodes are projected
    ! orthogonally onto the tangent X-Y plane (in the new coordinate
    ! system)
    ! 3) A quadratic is fitted in a weighted least squares
    ! sense to the data values at the projected nodes such that
    ! the value (associated with 'node') at (0,0) is interpolated
    ! 4) The polinomial coeficients are saved in var, al weel as the
    !   rotation parameters
    !
    ! COEFICIENT ORDER: 1:5 => x^2, xy, y^2, x, y
    !         p(0,0)=f(node)
    !
    ! On input:
    !       node : central node for gradient calculation
    !       mesh : with at least 7 nodes
    !       var  : with var%f filled
    ! On output:
    !      var%pol%c is filled with the coeficients
    !      var%pol%cx/cy/sx/sy is filled with the rotaion parameters
    !
    !----------------------------------------------------------------
    integer (i4), intent(in) :: node
    type(grid_structure), intent(in) :: mesh
    type(scalar_field), intent(inout) :: var

    ! globally fixed minimum and maximal lengths
    integer, parameter ::  lmn=7
    integer, parameter :: lmx=30

    ! local minimum and maximal lengths
    integer:: lmin
    integer:: lmax

    !List of neighbour nodes
    integer (i4):: npts(1:lmx)

    !Floating point tolerance for ill conditioning detection
    real (r8), parameter :: dtol=0.001

    ! Sum of distances
    real (r8):: sum

    !Indexes
    integer (i4):: lnp
    integer (i4)::  np
    integer (i4):: i
    integer (i4):: j
    integer (i4):: k
    integer (i4):: l
    !integer :: ier

    !Function of distance, -cos(angle) between 2 nodes
    !It is also the negative Z component (in the rotated coordi-
    !nate system) of an element of NPTS
    !increasing function of the angular distance
    !between K and NP. DF lies in the interval (-1,1).
    real (r8):: df
    real (r8):: rf

    ! Root-mean-square distance (in the rotated
    ! coordinate system) between the origin and
    ! the nodes (other than K) in the least
    ! squares fit.  The first 3 columns of A**T
    ! are scaled by 1/AVSQ, the next 2 by 1/AV.
    ! AVSQ =AV*AV:  accumulated in SUM
    real (r8):: av
    real (r8):: avsq

    ! Inverse of a radius of influence R which
    ! enters into the weights, R = 1+RF unless all ele-
    ! ments of NPTS are used in the fit (LNP =
    ! LMAX+1), in which case R is the distance
    ! function associated with some point more
    ! distant from node than NPTS(LMAX)
    real(r8):: rin

    ! Components of a plane rotation about the X-
    ! axis and Y-axis. Then define a rotation
    ! of 'node' to the north pole
    real (r8)::  cx
    real (r8):: sx
    real (r8):: cy
    real (r8):: sy

    !Givens rotation parameters
    real (r8)::  c
    real (r8):: s

    !Coordinates of NP in the rotated coordinate
    !system unless ZP < 0, in which case
    !(XP,YP,0) lies on the equator
    real (r8):: xp
    real (r8):: yp
    real (r8):: zp

    ! Weight for the equation coresponding to a point
    ! WT = (R-D)/(R*D) = 1/D - RIN, where D = 1-ZP
    real (r8):: wt

    ! Polinomial coeficients
    real (r8):: coef(1:5)

    !Augmented Matrix for the least square problem
    real (r8):: a(1:6,1:6)

    !Auxiliar variables
    real (r8):: dmin

    !Check for wrong input parameters
    if (mesh%nv < 7 .or. node<1 .or. node>mesh%nv) then
       print*, "POLINOM2LSFIT ERROR : Invalid node or mesh"
       print*, " Size of mesh must be > 7, mesh size:", mesh%nv
       print*, " node must be in mesh, node: ", node
       return
    end if
    lmin = min (lmn, mesh%nv)
    lmax = min (lmx, mesh%nv)

    !Test for the values position
    if(var%pos>0)then
       print*, "POLINOM2LSFIT ERROR : Variable not on triangle vertices"
       print*, "pos:", var%pos
       stop
    end if

    !Test if the node has enough neighbours
    !print*, node, mesh%v(node)%nnb
    if(mesh%v(node)%nnb<5)then
       print*, "POLINOM2LSFIT:Not enough neighbours"
       print*, "node:", node, " nnb:",mesh%v(node)%nnb
       stop
    end if

    ! Compute NPTS, LNP, AVSQ, AV, and R.
    !   Set NPTS to the closest LMIN-1 nodes to K.
    !   DF contains
    !   the negative Z component (in the rotated coordinate
    !   system) of the new node on return from GETNP.
    !
    sum = 0.
    npts (1) = node
    !For each neighbour node
    do lnp=1, mesh%v(node)%nnb
       !Save index
       npts(lnp+1)=mesh%v(node)%nb(lnp)
       !Save distance
       df=mesh%v(node)%nbd(lnp)
       sum = sum + 1. - df * df
    end do

    !R is arbitrarily increased by 25 percent.
    !This is done because the R distance must be larger then
    ! the nodes distance, otherwise the weights would be null
    rf = df + .25*abs(df)

    !Average of euclidian distances
    avsq = sum / real(lnp - 2)
    av = sqrt (avsq)
    rin = 1. / (1. + rf)

    ! Construct the rotation coeficients
    call constr (mesh%v(node)%p(1),mesh%v(node)%p(2),mesh%v(node)%p(3), &
         cx, sx, cy, sy)

    ! Set up the first 5 equations of the augmented regression
    !   matrix (transposed) as the columns of A, and zero out
    !   the lower triangle (upper triangle of A) with Givens
    !   rotations.
    do  i = 1, 5
       np = npts (i + 1)
       !Apply the rotation for this neighbour
       call aplyr (mesh%v(np)%p(1),mesh%v(np)%p(2),mesh%v(np)%p(3), cx, sx, cy, sy, xp, yp, zp)
       !Define the weight of the neighbour
       wt = 1. / (1. - zp) - rin
       !Set the matrix column refering to the neighbour
       !The last line contains the transposed right hand side of 'Ax=b'
       a(1:6,i)=(/ (xp**2)*wt/avsq, xp*yp*wt/avsq, (yp**2)*wt/avsq, &
            xp*wt/av, yp*wt/av, (var%f(np)-var%f(node))*wt/)
       if (i>1) then
          !Apply a givens rotation to the 2 columns in order
          ! to triangulate the A matrix between
          ! columns i and j
          do  j = 1, i - 1
             l = 6 - j
             call givens (a (j, j), a (j, i), c, s)
             call rotate (l, c, s, a (j+1, j), a (j+1, i) )
          end do
       end if
    end do

    !Up to now we have in A the system for quadratic interpolation of
    !the values in the following format:
    !  X 0 0 0 0
    !  X X 0 0 0
    !  X X X 0 0
    !  X X X X 0
    !  Y Y Y Y Y
    ! where X are non zero values of A, and Y are values of b, and
    ! we assume we want to solve the system Ax=b

    ! If the 'node' has more than 5 neighbours
    ! then we have to add these to the Least Square problem
    ! See Lawson - Solving Least Square Problem for detais
    dmin=0.0
    k=7
    !At this point lnp is the number os neighbours + 1
    if(k<lnp+1)then
       do i=k,lnp
          np = npts (i)
          call addcolumn(np)
       end do
    end if

    ! Test the system for ill-conditioning
    call condition()

    !Solve the triangular 5x5 system
    coef=solvelin5x5tri(a, av, avsq)

    !Save the coeficients
    var%pol(node)%c=coef

    !Save the rotation parameters
    var%pol(node)%cx=cx
    var%pol(node)%cy=cy
    var%pol(node)%sx=sx
    var%pol(node)%sy=sy

    return

  contains
    function solvelin5x5tri(a, av, avsq)
      !---------------------------------------------
      !Solve a 5x5 Triangular linear system wiht weighted columns
      !-----------------------------------------------

      !Augmented Matrix for the least square problem
      real (r8), intent(in) :: a(1:6,1:6)

      ! Root-mean-square distance (in the rotated
      ! coordinate system) between the origin and
      ! the nodes (other than K) in the least
      ! squares fit.  The first 3 columns of A**T
      ! are scaled by 1/AVSQ, the next 2 by 1/AV.
      ! AVSQ =AV*AV:  accumulated in SUM
      real (r8):: av
      real (r8):: avsq

      !Solution
      real (r8), dimension(1:5) :: x
      real (r8), dimension(1:5) :: solvelin5x5tri
      real (r8):: tmp
      integer:: i
      integer:: j

      do i=5, 1, -1
         tmp=0
         do j=5, i+1, -1
            tmp=tmp+a(j,i)*x(j)
         end do
         x(i)=(a(6,i)-tmp)/a(i,i)

      end do
      solvelin5x5tri(1:3)=x(1:3)/avsq
      solvelin5x5tri(4:5)=x(4:5)/av
      return
    end function solvelin5x5tri

    subroutine addcolumn(np)
      !-------------------------------------------
      !Adds a column to the system with a new point
      !-------------------------------------------
      integer (i4), intent(in) :: np

      call aplyr (mesh%v(np)%p(1),mesh%v(np)%p(2),mesh%v(np)%p(3), &
           cx, sx, cy, sy, xp, yp, zp)
      wt = 1. / (1. - zp) - rin
      ! We use only the last column of A
      a(1:6,6)=(/ (xp**2)*wt/avsq, xp*yp*wt/avsq, (yp**2)*wt/avsq, &
           xp*wt/av, yp*wt/av, (var%f(np)-var%f(node))*wt/)
      !Triangulate the matrix
      do j = 1, 5
         l = 6 - j
         call givens (a (j, j), a (j, 6), c, s)
         call rotate (l, c, s, a (j+1, j), a (j+1, 6) )
      end do

      return
    end subroutine addcolumn

    subroutine condition()
      !---------------------------------------------------------
      !CONDITION
      ! Verifies if the matrix is well conditioned, if not, adds
      ! new nodes until reached weel condiotioning or maximum number
      ! of nodes. In this case it tries to stabilize the system
      !-----------------------------------------------------------
      dmin = min (abs (a (1, 1) ), abs (a (2, 2) ), abs (a (3, 3) ), &
           abs (a (4, 4) ), abs (a (5, 5) ) )
      !If the system is ill-conditioned stop the program
      if(dmin<dtol)then
         print*, "PLINOM2LSFIT_TRV ERROR: System ill conditioned"
         print*, "Min diagonal: ", dmin
         print*, "Node: ", node
         stop
      end if

    end subroutine condition

  end subroutine pol2lsq_node


  !=====================================================================================
  !    POLINOMIAL Least Square Vector Reconstruction
  !=====================================================================================

  function vecrecon_lsq (p, var, stencil, mesh, pindex)
    !----------------------------------------------------------
    ! Evaluate the reconstructed vector using a least square fit
    !
    !  using 12 hexagon midpoint normal vector components.
    !----------------------------------------------------------------
    !Approximation point
    real (r8), intent(in) :: p(1:3)

    !Variable
    type(scalar_field), intent(inout) :: var

    !Stencil used
    character(len=4), intent(in) :: stencil

    !Mesh
    type(grid_structure), intent(in) :: mesh

    !Index of the nearest node to p (optional)
    !  or the triangle that p belongs
    integer(i4), optional:: pindex

    !Node/Triangle for the center of reconstruction
    integer (i4):: k

    !Returning value of approximation
    real (r8):: vecrecon_lsq(1:3)

    if(var%pos/=3)then
       print*, "vecrecon_lsq warning: vector field given in incorrect position", var%pos
       return
    end if

    k=0

    !Get nearest node
    if(present(pindex))then
       k=pindex
    end if

    if(k==0)then !Get cell index
       select case (trim(stencil))
       case("hxe", "hx")
          !Reconstruct using cell centered least square polinomial
          !  with 12 edges

          !Get nearest node
          k=getnearnode(p, mesh)
          !print*, k
       case("trc")
          !Reconstruct using triangle as center of least square polinomial
          ! uses 9 edges

          !Get triangle
          k=gettr(p, mesh)

       case default
          print*, "vecrecon_lsq error: unknown stencil", trim(stencil)
          stop
       end select
    end if
    !print*, k
    call vecrecon_lsqfitpol (k, stencil, var, mesh)

    vecrecon_lsq=vecrecon_lsqeval(k, p, var)

    return

  end function vecrecon_lsq


  function vecrecon_lsqeval(i, p, var)
    !----------------------------------------------------------
    ! Evaluate the reconstructed vector using a least square fit
    !
    !  using 12 hexagon midpoint normal vector components.
    !----------------------------------------------------------------
    !Approximation point
    real (r8), intent(in) :: p(1:3)

    !Variable
    type(scalar_field), intent(inout) :: var

    !Rotation of the point
    real (r8):: x
    real (r8):: y
    real (r8):: z

    !Coeficients and rotated reconstructed vector
    real (r8):: a(1:2)
    real (r8):: b(1:2)
    real (r8):: c(1:2)
    real (r8):: nr(1:2)

    !Nod/Trianglee for the reference
    integer (i4), intent(in) :: i

    !Returning value of approximation
    real (r8):: vecrecon_lsqeval(1:3)

    vecrecon_lsqeval(1:3)=0._r8

    !print*, p
    !Rotate point using node as reference to pole (p -> (x,y,z))
    call aplyr (p(1), p(2), p(3), &
         var%pol(i)%cx, var%pol(i)%sx, var%pol(i)%cy, var%pol(i)%sy, &
         x, y, z)
    !print*, var%pol(node)%cx, var%pol(node)%sx, var%pol(node)%cy, var%pol(node)%sy
    !print*, x, y, z

    a=var%pol(i)%c(1:2)
    b=var%pol(i)%c(3:4)
    c=var%pol(i)%c(5:6)
    nr=a*x+b*y+c

    !Rotate vector to original position
    call aplyrt (nr(1), nr(2), &
         var%pol(i)%cx, var%pol(i)%sx, var%pol(i)%cy, var%pol(i)%sy, &
         vecrecon_lsqeval)
    return
  end function vecrecon_lsqeval

  subroutine vecrecon_lsqfitpol (kc, stencil, var, mesh)
    !----------------------------------------------------------
    ! POLINOMIAL 1st ORDER LEAST SQUARE vector reconstruction
    !
    ! On input:
    !       kc : central triangle/node index for reconstruction
    !            calculation
    !       stencil:
    !           - trc : 9 point centered on triangle
    !           - hxe : 12 point centered on node
    !       mesh : with at least 12 nodes
    !       var  : with var%f filled with normal components at edges
    ! On output:
    !      var%pol%c is filled with the coeficients
    !      var%pol%cx/cy/sx/sy is filled with the rotaion parameters
    !
    !----------------------------------------------------------------

    !Node or triangle index serving as center of reconstruction
    integer (i4), intent(in) :: kc

    !Stencil used
    character(len=4), intent(in) :: stencil

    !Mesh structure
    type(grid_structure), intent(in) :: mesh

    !Variable used in reconstruction
    type(scalar_field), intent(inout) :: var

    !Globally fixed minimum and maximal lengths
    integer, parameter :: lmx=30

    !List of neighbour nodes
    integer (i4):: eds(1:lmx)

    !Floating point tolerance for ill conditioning detection
    real (r8), parameter :: dtol=0.0001

    ! Distances
    real (r8):: d
    real (r8):: dsum
    real (r8):: rmax
    real (r8):: dmax

    !Normal vectors
    real (r8):: nx
    real (r8):: ny
    real (r8):: nz
    real (r8):: nr(1:3)

    !Indexes
    integer (i4):: i
    integer (i4):: j
    integer (i4):: k
    integer (i4):: l
    integer (i4):: n
    integer (i4):: ed
    !integer :: ier

    ! Weighting Average
    real (r8):: av

    ! Inverse of a radius of influence R which
    ! enters into the weights, R = 1+RF unless all ele-
    ! ments of NPTS are used in the fit (LNP =
    ! LMAX+1), in which case R is the distance
    ! function associated with some point more
    ! distant from node than NPTS(LMAX)
    real(r8):: rinv

    ! Components of a plane rotation about the X-
    ! axis and Y-axis. Then define a rotation
    ! of 'node' to the north pole
    real (r8)::  cx
    real (r8):: sx
    real (r8):: cy
    real (r8):: sy

    !Givens rotation parameters
    real (r8)::  c
    real (r8):: s

    !Coordinates of NP in the rotated coordinate
    !system unless ZP < 0, in which case
    !(XP,YP,0) lies on the equator
    real (r8):: xp
    real (r8):: yp
    real (r8):: zp

    ! Weight for the equation coresponding to a point
    ! WT = (R-D)/(R*D) = 1/D - RIN, where D = 1-ZP
    real (r8):: wt

    ! Polinomial coeficients
    real (r8):: coef(1:6)

    !Augmented Matrix for the least square problem
    real (r8):: a(1:7,1:7)

    !Reference point
    real (r8):: pref(1:3)

    !Auxiliar variables
    real (r8):: dmin

    !Test for the values position
    if(var%pos/=3)then
       print*, "vecrecon_lsfitpol ERROR : Variable not on hexagon edge midpoints"
       print*, "pos:", var%pos
       stop
    end if

    !Save the nearest edges to be used
    dsum = 0._r8
    dmax=-10._r8
    dmin=100000.0_r8


    select case(trim(stencil))
    case("trc") !9 point stencial centered on triangle

       !Alocate space if necessary
       if(.not.allocated(var%pol))then
          allocate(var%pol(1:mesh%nt))
       end if

       if(.not.allocated(var%pol(kc)%c))then
          allocate(var%pol(kc)%c(1:6))
          !else !Nothing to be done
          !   return
       end if

       !Check if index given is ok
       if(kc>mesh%nt .or. kc <1)then
          print*, "vecrecon_lsfitpol error: Index given is not correct", kc, mesh%nt
          stop
       end if

       !Reference point (center of polinomial)
       pref=mesh%tr(kc)%b%p

       !For each neighbour triangle, get the edges
       j=0
       do i = 1, 3
          k=mesh%tr(kc)%nb(i)
          do l=1, 3
             j=j+1
             !Save index of edge
             eds(j)=mesh%tr(k)%ed(l)
             !Calculate distance (-cos(angle))
             d=-dot_product(pref, mesh%edhx(eds(j))%c%p)
             dmax=max(d, dmax)
             dsum = dsum + 1. - d ** 2
          end do
       end do

    case("hxe", "hx") !hxe 12 stencil

       !Check if index given is ok
       if(kc>mesh%nv .or. kc <1)then
          print*, "vecrecon_lsfitpol error: Index given is not correct", kc, mesh%nv
          stop
       end if

       !Alocate space if necessary
       if(.not.allocated(var%pol))then
          allocate(var%pol(1:mesh%nv))
       end if

       !Check if coeficients already exist,
       ! if so, they will be overwritten
       if(.not.allocated(var%pol(kc)%c))then
          allocate(var%pol(kc)%c(1:6))
          !else !nothing to be done
          !   return
       end if

       pref=mesh%v(kc)%p

       !For each neighbour node
       do i = 1, mesh%v(kc)%nnb
          !Save index of edge
          eds(i)=mesh%v(kc)%ed(i)
          !Calculate distance (-cos(angle))
          d=-dot_product(mesh%v(kc)%p, mesh%edhx(eds(i))%c%p)
          dmax=max(d, dmax)
          dsum = dsum + 1. - d ** 2
       end do

       !Now add the surounding edges depending on the stencil
       if(trim(stencil)=="hxe" .or. mesh%v(kc)%nnb < 6 .or. dmin<dtol )then
          call addextraeqs()
       end if
    case default
       print*, "vecrecon_lsqfitpol error: unknown stencil", trim(stencil)
       stop
    end select

    !Number of edges used
    n=j

    !Average of distances
    av = dsqrt (dsum / real(n, r8))

    !The radius of action of the reconstruction is arbitrarily
    !  increased by 25 percent. This is done because the R
    !  distance must be larger then
    ! the nodes distance, otherwise weights would be null
    rmax = dmax + .25*abs(dmax)
    rinv = 1. / (1. + rmax)
    !print*, node, n
    !print*, "EDS:", eds
    !print*, "AV DMAX RMAX RINV", av, dmax, rmax, rinv

    ! Construct the rotation coeficients
    call constr(pref(1), pref(2), pref(3), cx, sx, cy, sy)

    !Build 6x6 least square system
    call buildsys()

    !do i=1,n+1
    !  print '(7f12.4)', a(i,1:7)
    !end do

    ! Test the system for ill-conditioning
    dmin = min (abs (a (1, 1) ), abs (a (2, 2) ), abs (a (3, 3) ), &
         abs (a (4, 4) ), abs (a (5, 5) ), abs (a (6, 6) ) )

    !If the system is ill-conditioned stop the program
    if(dmin<dtol .and. trim(stencil)=="hx")then !Add more nodes
       print*, "Warning on vecrecon_lsqfitpol: System ill conditioned - Adding more edges", &
            " Node:", kc, trim(stencil), " -> hxe ", dmin
       call addextraeqs()
       call buildsys()
    end if

    !Up to now we have in A the system for least sq of
    !the values in the following format:
    !  X 0 0 0 0
    !  X X 0 0 0
    !  X X X 0 0
    !  X X X X 0
    !  Y Y Y Y Y
    ! where X are non zero values of A, and Y are values of b, and
    ! we assume we want to solve the system Ax=b

    ! Now we add the remmaining edges using tha last column
    do i=7, n
       ed = eds (i)
       !Apply the rotation for this point
       call aplyr (mesh%edhx(ed)%c%p(1), mesh%edhx(ed)%c%p(2), &
            mesh%edhx(ed)%c%p(3), cx, sx, cy, sy, xp, yp, zp)

       !Define the weight of the edge
       wt = 1. / (1. - zp) - rinv

       !print*, ed, zp, wt

       !Get Normal vector component at this edge and rotate it
       nr=mesh%edhx(ed)%nr
       call aplyr (mesh%edhx(ed)%nr(1), mesh%edhx(ed)%nr(2), &
            mesh%edhx(ed)%nr(3), cx, sx, cy, sy, nx, ny, nz)

       !Ignore z normal component and correct norm of 2d vector
       nx=nx/dsqrt(1-nz**2)
       ny=ny/dsqrt(1-nz**2)

       ! We use only the last column of A
       a(1:7,7)=(/ (xp*nx)*wt/av, (xp*ny)*wt/av, (yp*nx)*wt/av, &
            (yp*ny)*wt/av, (nx)*wt, (ny)*wt, var%f(ed)*wt /)

       !Triangulate the matrix
       do j = 1, 6
          l = 7 - j
          call givens (a (j, j), a (j, 7), c, s)
          call rotate (l, c, s, a (j+1, j), a (j+1, 7) )
       end do
    end do

    !do i=1,n+1
    !  print '(7f12.4)', a(i,1:7)
    !end do

    ! Test the system for ill-conditioning
    dmin = min (abs (a (1, 1) ), abs (a (2, 2) ), abs (a (3, 3) ), &
         abs (a (4, 4) ), abs (a (5, 5) ), abs (a (6, 6) ) )
    !If the system is ill-conditioned stop the program
    if(dmin<dtol)then
       print*, "vecrecon_lsqfitpol ERROR: System ill conditioned"
       print*, "Min diagonal: ", dmin
       print*, "Node/Triangle: ", kc, stencil
       stop
    end if

    !Solve the triangular 6x6 system
    coef=solvelintri(a, 6)

    coef(1:4)=coef(1:4)/av
    !print*, coef

    !Save the coeficients
    var%pol(kc)%c=coef

    !Save the rotation parameters
    var%pol(kc)%cx=cx
    var%pol(kc)%cy=cy
    var%pol(kc)%sx=sx
    var%pol(kc)%sy=sy

    return
  contains

    subroutine addextraeqs()
      !---------------------------------------------------
      ! Add 6 extra edges to the least squares
      !---------------------------------------------------
      j=n
      !For each neighbour node
      do i=1, mesh%v(kc)%nnb
         !Edge counter
         j=i+mesh%v(kc)%nnb
         !Get the triangle
         k=mesh%v(kc)%tr(i)
         !Get the edge index that is missing from the list
         eds(j)=0
         do l=1,3
            if(mesh%tr(k)%ed(l) /= eds(i) .and. &
                 mesh%tr(k)%ed(l) /= eds(modint(i+1, mesh%v(kc)%nnb)) ) then
               !Add edge to list
               eds(j)=mesh%tr(k)%ed(l)
               exit
            end if
         end do
         if(eds(j)==0)then
            print*, "vecrecon_lsfit ERROR: could not add surrounding edge"
            stop
         end if
         !Calculate distance (-cos(angle))
         d=-dot_product(mesh%v(kc)%p, mesh%edhx(eds(j))%c%p)
         dmax=max(d, dmax)
         dsum = dsum + 1. - d ** 2
      end do

      !set n
      n=j

      !Average of distances
      av = dsqrt (dsum / real(n, r8))

      !The radius
      rmax = dmax + .25*abs(dmax)
      rinv = 1. / (1. + rmax)

    end subroutine addextraeqs


    subroutine buildsys()
      ! Set up the first 6 equations of the augmented regression
      !   matrix (transposed) as the columns of A, and zero out
      !   the lower triangle (upper triangle of A) with Givens
      !   rotations.
      do  i = 1, 6
         ed = eds (i)
         !Apply the rotation for this point
         call aplyr (mesh%edhx(ed)%c%p(1), mesh%edhx(ed)%c%p(2), &
              mesh%edhx(ed)%c%p(3), cx, sx, cy, sy, xp, yp, zp)

         !Define the weight of the edge
         wt = 1. / (1. - zp) - rinv

         !print*, ed, zp, wt

         !Get Normal vector component at this edge and rotate it
         nr=mesh%edhx(ed)%nr
         call aplyr (mesh%edhx(ed)%nr(1), mesh%edhx(ed)%nr(2), &
              mesh%edhx(ed)%nr(3), cx, sx, cy, sy, nx, ny, nz)

         !Ignore z normal component and correct norm of 2d vector
         nx=nx/dsqrt(1-nz**2)
         ny=ny/dsqrt(1-nz**2)

         !Set the matrix column refering to the edge
         !The last line contains the transposed right hand side of 'Ax=b'
         a(1:7,i)=(/ (xp*nx)*wt/av, (xp*ny)*wt/av, (yp*nx)*wt/av, &
              (yp*ny)*wt/av, (nx)*wt, (ny)*wt, (var%f(ed))*wt/)
         if (i>1) then
            !Apply a givens rotation to the 2 columns in order
            ! to triangulate the A matrix between
            ! columns i and j
            do  j = 1, i - 1
               l = 7 - j
               call givens (a (j, j), a (j, i), c, s)
               call rotate (l, c, s, a (j+1, j), a (j+1, i) )
            end do
         end if
      end do

    end subroutine buildsys

  end subroutine vecrecon_lsqfitpol

  !========================================================================
  !    AUXILIAR ROUTINES FOR LEAST SQUARES SOLVER
  !========================================================================

  subroutine aplyr (x, y, z, cx, sx, cy, sy, xp, yp, zp)
    !-------------------------------------------------------
    !  APLYR
    !
    !   This subroutine applies the rotation r defined by sub-
    ! routine constr to the unit vector (x y z)**t, i,e. (x,y,z)
    ! is rotated to (xp,yp,zp).  if (xp,yp,zp) lies in the
    ! southern hemisphere (zp < 0), (xp,yp) are set to the
    ! coordinates of the nearest point of the equator, zp re-
    ! maining unchanged.
    !
    ! On input:
    !       x,y,z = coordinates of a point on the unit sphere.
    !       cx,sx,cy,sy = elements of the rotation defined by
    !                     subroutine constr.
    ! On output:
    !       xp,yp,zp = coordinates of the rotated point on the
    !                  sphere unless zp < 0, in which case
    !                  (xp,yp,0) is the closest point of the
    !                  equator to the rotated point.  storage
    !                  for xp, yp, and zp may coincide with
    !                  storage for x, y, and z, respectively,
    !                  if the latter need not be saved.
    !---------------------------------------------------------
    real (r8) :: x
    real (r8) :: y
    real (r8) :: z
    real (r8) :: cx
    real (r8) :: sx
    real (r8) :: cy
    real (r8) :: sy
    real (r8) :: xp
    real (r8) :: yp
    real (r8) :: zp
    real (r8) :: t

    t = sx * y + cx * z
    yp = cx * y - sx * z
    zp = sy * x + cy * t
    xp = cy * x - sy * t
    !if (zp>=0.) return

    ! move (xp,yp,zp) to the equator.
    !t = sqrt (xp * xp + yp * yp)
    !if (t > 0.) then
    !   xp = xp / t
    !   yp = yp / t
    !else
    !   ! move the south pole to an arbitrary point of the equator.
    !   xp = 1.
    !   yp = 0.
    !end if
    return
  end subroutine aplyr

  subroutine aplyrt (g1p, g2p, cx, sx, cy, sy, g)
    !----------------------------------------------------------
    ! APLYRT
    !   This subroutine applies the inverse (transpose) of the
    ! rotation defined by subroutine constr to the vector
    ! (g1p g2p 0)**t, i.e., the gradient (g1p,g2p,0) in the rot-
    ! ated coordinate system is mapped to (g1,g2,g3) in the
    ! original coordinate system.
    !
    ! On input:
    !       g1p,g2p = x and y components, respectively, of the
    !                 gradient in the rotated coordinate system.
    !       cx,sx,cy,sy = elements of the rotation r constructed
    !                     by subroutine constr.
    ! On output:
    !       g = x, y, and z components (in that order) of the
    !           inverse rotation applied to (g1p,g2p,0) --
    !           gradient in the original coordinate system.
    !---------------------------------------------------
    real (r8) :: g1p
    real (r8) :: g2p
    real (r8) :: cx
    real (r8) :: sx
    real (r8) :: cy
    real (r8) :: sy
    real (r8) :: g (3)
    real (r8) :: t
    t = sy * g1p
    g (1) = cy * g1p
    g (2) = cx * g2p - sx * t
    g (3) = - sx * g2p - cx * t
    return
  end subroutine aplyrt

  subroutine constr (xk, yk, zk, cx, sx, cy, sy)
    !--------------------------------------------------------
    ! CONSTR
    !   This subroutine constructs the elements of a 3 by 3
    ! orthogonal matrix r which rotates a point (xk,yk,zk) on
    ! the unit sphere to the north pole, i.e.,
    !
    !      (xk)     (cy  0 -sy)   (1   0   0)   (xk)     (0)
    !  r * (yk)  =  ( 0  1   0) * (0  cx -sx) * (yk)  =  (0)
    !      (zk)     (sy  0  cy)   (0  sx  cx)   (zk)     (1)
    !
    ! On input:
    !       xk,yk,zk = components of a unit vector to be
    !                  rotated to (0,0,1).
    ! On output:
    !       cx,sx,cy,sy = elements of r:  cx,sx define a rota-
    !                     tion about the x-axis and cy,sy define
    !                     a rotation about the y-axis.
    !---------------------------------------------------------
    real (r8) :: xk
    real (r8) :: yk
    real (r8) :: zk
    real (r8) :: cx
    real (r8) :: sx
    real (r8) :: cy
    real (r8) :: sy
    cy = sqrt (yk * yk + zk * zk)
    sy = xk
    if (cy/=0.) then
       cx = zk / cy
       sx = yk / cy
    else
       ! (xk,yk,zk) lies on the x-axis.
       cx = 1.
       sx = 0.
    endif
    return
  end subroutine constr

  subroutine givens (a, b, c, s)
    !--------------------------------------------------------
    !  GIVENS
    !
    !   This subroutine constructs the givens plane rotation,
    !
    !           ( c  s)
    !       g = (     ) , where c*c + s*s = 1,
    !           (-s  c)
    !
    ! Which zeros the second component of the vector (a,b)**t
    ! (transposed).  subroutine rotate may be called to apply
    ! the transformation to a 2 by n matrix.
    !
    !   This routine is identical to subroutine srotg from the
    ! linpack blas (basic linear algebra subroutines). It was
    ! obtained from R Renka's ssrfpack
    !
    ! On input:
    !       a,b = components of the vector defining the rota-
    !             tion.  these are overwritten by values r
    !             and z (described below) which define c and s.
    !
    ! on output:
    !       a = signed euclidean norm r of the input vector:
    !           r = +/-sqrt(a*a + b*b)
    !       b = value z such that:
    !             c = sqrt(1-z*z) and s=z if abs(z) .le. 1, and
    !             c = 1/z and s = sqrt(1-c*c) if abs(z) > 1.
    !       c = +/-(a/r) or 1 if r = 0.
    !       s = +/-(b/r) or 0 if r = 0.
    !------------------------------------------------------------
    real (r8) :: a
    real (r8) :: b
    real (r8) :: c
    real (r8) :: s
    real (r8) :: aa
    real (r8) :: bb
    real (r8) :: r
    real (r8) :: u
    real (r8) :: v
    aa = a
    bb = b
    if (abs (aa) > abs (bb) ) then
       u = aa + aa
       v = bb / u
       r = sqrt (.25 + v * v) * u
       c = aa / r
       s = v * (c + c)
       !r has the sign of a, c > 0, s has sign(a)*sign(b)
       b = s
       a = r
    elseif (bb/=0.) then
       ! abs(a) .le. abs(b).
       u = bb + bb
       v = aa / u
       ! store r in a.
       a = sqrt (.25 + v * v) * u
       s = bb / a
       c = v * (s + s)
       ! note that r has the sign of b, s > 0, and c has
       !   sign(a)*sign(b).
       b = 1.
       if (c/=0.) b = 1. / c
    else
       ! a = b = 0.
       c = 1.
       s = 0.
    endif
    return
  end subroutine givens

  subroutine rotate (n, c, s, x, y)
    !-----------------------------------------------------------
    !  ROTATE
    !                                                ( C  S)
    !   This subroutine applies the Givens rotation  (     )  to
    !                                                (-S  C)
    !                    (X(1) ... X(N))
    ! the 2 by N matrix  (             ) .
    !                    (Y(1) ... Y(N))
    !
    !   This routine is identical to Subroutine SROT from the
    ! LINPACK BLAS (Basic Linear Algebra Subroutines).
    !
    !       N = Number of columns to be rotated.
    !       C,S = Elements of the Givens rotation.  Refer to
    !             Subroutine GIVENS.
    !       X,Y = Arrays containing the in/rotated vectors
    !------------------------------------------------------
    integer:: n
    real (r8):: c
    real (r8):: s
    real (r8):: x (n)
    real (r8):: y (n)
    integer:: i
    real (r8):: xi
    real (r8):: yi
    do i = 1, n
       xi = x (i)
       yi = y (i)
       x (i) = c * xi + s * yi
       y (i) = - s * xi + c * yi
    end do
    return
  end subroutine rotate


  !========================================================================
  !    MAP PLOTTING ROUTINES
  !========================================================================

  subroutine plot_scalarfield_sphgrid(nlon, nlat, var, filename)
    !----------------------------------------------
    ! PLOT_SCALARFIELD_shpgrid
    !   Writes in 'filename' a uniform mesh with
    !   nlat latitudes and nlon longitudes
    !   using linear interpolation
    !   This is intended to be used for GMT plots
    !---------------------------------------------
    !Uniform grid sizes
    integer (i4), intent(in)  :: nlat
    integer (i4), intent(in)  :: nlon

    !Variable for uniform grid
    real (r8), intent(in) :: var(1:nlon, 1:nlat)
    real*4, allocatable :: buffer(:,:)

    !File name for output
    character (len=256), optional :: filename

    !Aux variable
    real(r8):: dlon
    real(r8):: dlat
    real(r8):: tlon
    real(r8):: tlat
    integer (i4):: i
    integer (i4):: j
    integer (i4):: k
    integer (i4):: iunit

    if(.not.present(filename))then
       filename="scalar"
    end if
    print*, " Plotting scalarfield ( lonlat ): "
    filename=trim(datadir)//trim(filename)//".dat"
    write(*, '(a)') "   "//trim(filename)

    allocate(buffer(3, 1:nlat*nlon))

    !Pixel registration mode (GMT)
    k=1
    dlat=180._r8/real(nlat, r8)
    dlon=360._r8/real(nlon, r8)
    tlat=-90+dlat/2._r8
    lat: do j=1,nlat
       tlon=-180+dlon/2._r8
       lon: do i=1,nlon
          buffer(1, k)=tlon
          buffer(2, k)=tlat
          buffer(3, k)=var(i,j)
          k=k+1
          !write(iunit) tlon, tlat, var(i,j)
          tlon=tlon+dlon
       end do lon
       tlat=tlat+dlat
    end do lat

    !Write values on file
    call getunit(iunit)

    open(iunit,file=filename, status='replace', access='stream', form='unformatted')

    write (iunit) buffer

    close(iunit)

    return
  end subroutine plot_scalarfield_sphgrid

  subroutine plot_scalarfield(var, mesh, kindinterp)
    !----------------------------------------------
    ! PLOT_SCALARFIELD
    !   Writes in 'var%name file' a uniform mesh with
    !   nlat latitudes and nlon longitudes
    !   using interpolation
    !   This is intended to be used only for GMT plots
    !---------------------------------------------

    !Mesh structure
    type(grid_structure), intent(in) :: mesh

    !Variable to be plotted
    type(scalar_field), intent(inout) :: var

    !Uniform grid sizes
    !integer (i4), parameter  :: nlat=180, nlon=360
    !integer (i4), parameter  :: nlat=360, nlon=720
    !integer (i4), parameter  :: nlat=720
    !integer (i4), parameter  :: nlon=1440
    integer (i4), parameter  :: nlat=720
    integer (i4), parameter  :: nlon=1440

    !Kind of interpolation -see scalar_interpol
    !If not given, uses nearest node interpolation
    character (len=8), optional :: kindinterp
    character (len=8):: kindinterpol

    !Aux variable
    type(scalar_field):: var0

    !Variable for uniform grid
    real (r8), allocatable :: varunif(:, :)
    real (r4), allocatable :: buffer(:,:)
    real (r8):: dlon
    real (r8):: dlat
    real (r8):: tlon
    real (r8):: tlat
    integer (i4):: i
    integer (i4):: j
    integer (i4):: k
    integer (i4):: iunit

    !File name for output
    character (len=256):: filename

    allocate(varunif(1:nlon, 1:nlat))

    if(present(kindinterp))then
       kindinterpol=kindinterp
    else
       select case (var%pos)
       case(0, 4) !Values on Voronoi cells
          kindinterpol='neartrv'
       case(1,5) !Values on triangles
          kindinterpol='neartrc'
       case(2, 6) !Values on tr edges
          kindinterpol='neartred'
       case(3) !Values on edges
          kindinterpol='nearhxed'
       case default
          print*, "plot_scalarfield error: variable position not plottable", var%pos
          stop
       end select
    end if

    if( len_trim(var%name)==0) then
       filename=trim(datadir)//"scalar_"//trim(mesh%name)//".dat"
    else
       filename=trim(datadir)//trim(var%name)//"_"//trim(mesh%name)//".dat"
    end if
    print*, " Plotting scalarfield ( with ", trim(kindinterpol), " ): "
    write(*, '(a)') "   "//trim(filename)
    !print*, "   "//trim(var%name)
!    print*, "    "//trim(filename)

    select case(var%pos)
    case(1, 5) !Values on triangles
       !If values on triangle center, and kindinterpol not 'neartrc'
       !   interpolate with linear
       !   to triangle vertex points (this must be improved)
       if(trim(kindinterpol)=='neartrc')then
          call scalar_remap_geo2ll(nlon, nlat, var, mesh, varunif, kindinterpol)
       else
          call scalar_remap_trc2trv(var, var0, mesh)
          ! Remap to lat lon grid
          kindinterpol='neartrv'
          call scalar_remap_geo2ll(nlon, nlat, var0, mesh, varunif, kindinterpol)
       end if

    case(0, 4) !Values on Voronoi cells
       ! Remap to lat lon grid
       call scalar_remap_geo2ll(nlon, nlat, var, mesh, varunif, kindinterpol)

    case(2, 3, 6) !Values on edges
       ! Remap to lat lon grid
       call scalar_remap_geo2ll(nlon, nlat, var, mesh, varunif, kindinterpol)

    case default
       print*, "plot_scalarfield error: variable position not plottable", var%pos
       stop
    end select

    allocate(buffer(3, 1:nlat*nlon))
    !Pixel registration mode (GMT)

    !Copy data to a buffer
    k=1
    dlat=180._r8/real(nlat, r8)
    dlon=360._r8/real(nlon, r8)
    !Pixel registration mode (GMT) (at midpoint of cell)
    tlat=-90+dlat/2._r8
    lat: do j=1,nlat
       tlon=-180+dlon/2._r8
       lon: do i=1,nlon
          buffer(1, k)=tlon
          buffer(2, k)=tlat
          buffer(3, k)=varunif(i,j)
          k=k+1
          !write(iunit) tlon, tlat, varunif(i,j)
          tlon=tlon+dlon
       end do lon
       tlat=tlat+dlat
    end do lat

    !Write values on file
    call getunit(iunit)
    !print*, "-----------------------------"
    !print*, trim(filename)
    !print*
    open(iunit,file=filename, status='replace', access='stream', form='unformatted')

    !Write whole block to file (much faster)
    write(iunit) buffer

    close(iunit)

    return
  end subroutine plot_scalarfield

  subroutine plot_vectorfield_sphgrid(nlon, nlat, var, filename)
    !----------------------------------------------
    ! PLOT_VECTORFIELD_shpgrid
    !   Writes in 'filename' a uniform mesh with
    !   nlat latitudes and nlon longitudes
    !   using linear interpolation
    !   This is intended to be used for GMT plots
    !---------------------------------------------
    !Uniform grid sizes
    integer (i4), intent(in)  :: nlat
    integer (i4), intent(in)  :: nlon

    !Variable for uniform grid
    type(vector), intent(in) :: var(1:nlon, 1:nlat)
    real (r8), allocatable :: buffer(:,:)

    !File name for output
    character (len=256), optional :: filename

    !Aux variable
    real(r8):: dlon
    real(r8):: dlat
    real(r8):: tlon
    real(r8):: tlat
    real(r8):: vlon
    real(r8):: vlat
    real(r8):: p(1:3)
    real (r8):: length
    real (r8):: angle
    integer (i4):: i
    integer (i4):: j
    integer (i4):: k
    integer (i4):: iunit

    if(.not.present(filename))then
       filename="vector"
    end if
    print*, " Plotting vectorfield ( lonlat ): "
    filename=trim(datadir)//trim(filename)//".dat"
    write(*, '(a)') "   "//trim(filename)
    !print*, "   "//trim(filename)

    allocate(buffer(4, 1:nlat*nlon))

    !Pixel registration mode (GMT)
    k=1
    dlat=180._r8/real(nlat, r8)
    dlon=360._r8/real(nlon, r8)
    tlat=-90+dlat/2._r8
    lat: do j=1,nlat
       tlon=-180+dlon/2._r8
       lon: do i=1,nlon
          buffer(1, k)=tlon
          buffer(2, k)=tlat
          call sph2cart(tlon*deg2rad, tlat*deg2rad, p(1), p(2), p(3))
          call convert_vec_cart2sph(p, var(i,j)%v, vlon, vlat)
          angle=datan2(vlon, vlat)*rad2deg
          buffer(3, k)=angle
          length=dsqrt(vlat**2+vlon**2)
          buffer(4, k)=length
          k=k+1
          tlon=tlon+dlon
       end do lon
       tlat=tlat+dlat
    end do lat

    !Write values on file
    call getunit(iunit)

    open(iunit,file=filename, status='replace')
    write(iunit, *) "#    Longitude  Latitude  Direction Length "
    write (iunit,'(2f12.6, 2f30.15)') buffer

    close(iunit)

    return
  end subroutine plot_vectorfield_sphgrid

  subroutine plot_cart_vectorfield(vec, mesh)
    !----------------------------------------------
    ! PLOT_CART_VECTORFIELD
    !   Writes in 'filename' the R3 vectors
    !   Converts to west-east/south-north vector type
    !   to be read by gmt
    !   Plots only up to 'num' of the vectors given
    !   so that the graph is not too poluted
    !---------------------------------------------
    type(grid_structure), intent(in) :: mesh
    type(vector_field_cart), intent(in) :: vec
    character (len=256):: filename

    !Auxiliar varibles
    real (r8):: tlat
    real (r8):: tlon
    real (r8):: angle
    real (r8):: length
    integer (i4):: iunit
    integer (i4):: i

    call getunit(iunit)
    if( len_trim(vec%name)==0) then
       filename="vector_"//"_"//trim(mesh%name)
    else
       filename=trim(vec%name)//"_"//trim(mesh%name)
    end if

    print*, " Plotting vectorfield ( ptwise vector ): "
    filename=trim(datadir)//trim(filename)//".dat"
    write(*, '(a)') "   "//trim(filename)
    !print*, "   "//trim(filename)

    !print*, vec%pos, vec%name
    open(iunit,file=filename, status='replace')

    write(iunit, *) "#    Longitude  Latitude  Direction  Length "
    select case(vec%pos)
    case(0) !TR vertices
       do i=1, min(10242, mesh%nv)
          call convert_vec_cart2sph(mesh%v(i)%p, vec%p(i)%v, tlon, tlat)
          length=dsqrt(tlat**2+tlon**2)
          if(length>eps)then
             angle=datan2(tlon, tlat)*rad2deg
             write(iunit,'(2f12.6, 2f30.15)') mesh%v(i)%lon*rad2deg, mesh%v(i)%lat*rad2deg, &
                  angle, length
          end if
       end do
    case(1) !TR centers
       do i=1, min(mesh%nt, 20480)
          call convert_vec_cart2sph(mesh%tr(i)%c%p, vec%p(i)%v, tlon, tlat)

          length=dsqrt(tlat**2+tlon**2)
          if(length>0)then
             angle=datan2(tlon, tlat)*rad2deg
             write(iunit,'(2f12.6, 2f30.15)') mesh%tr(i)%c%lon*rad2deg, mesh%tr(i)%c%lat*rad2deg, &
                  angle, length
          end if
       end do
    case(2, 6) !TR edges
       do i=1, min(mesh%ne, 30720)
          call convert_vec_cart2sph(mesh%ed(i)%c%p, vec%p(i)%v, tlon, tlat)
          length= norm(vec%p(i)%v) !dsqrt(tlat**2+tlon**2)
          if(length>eps)then
             angle=datan2(tlon, tlat)*rad2deg
             write(iunit,'(2f12.6, 2f30.15)') mesh%ed(i)%c%lon*rad2deg, mesh%ed(i)%c%lat*rad2deg, &
                  angle, length
          end if
       end do
    case(3) !HX edges
       do i=1, min(mesh%ne, 30720)
          call convert_vec_cart2sph(mesh%edhx(i)%c%p, vec%p(i)%v, tlon, tlat)
          length= norm(vec%p(i)%v) !dsqrt(tlat**2+tlon**2)
          if(length>eps)then
             angle=datan2(tlon, tlat)*rad2deg
             write(iunit,'(2f12.6, 2f30.15)') mesh%edhx(i)%c%lon*rad2deg, mesh%edhx(i)%c%lat*rad2deg, &
                  angle, length
          end if
       end do
    case(4) !Voronoi centroids
       do i=1, min(10242, mesh%nv)
          call convert_vec_cart2sph(mesh%hx(i)%b%p, vec%p(i)%v, tlon, tlat)
          length=dsqrt(tlat**2+tlon**2)
          if(length>eps)then
             angle=datan2(tlon, tlat)*rad2deg
             write(iunit,'(2f12.6, 2f30.15)') mesh%v(i)%lon*rad2deg, mesh%v(i)%lat*rad2deg, &
                  angle, length
          end if
       end do
    case(5) !TR centroids
       do i=1, min(mesh%nt, 20480)
          call convert_vec_cart2sph(mesh%tr(i)%b%p, vec%p(i)%v, tlon, tlat)
          length=dsqrt(tlat**2+tlon**2)
          if(length>eps)then
             angle=datan2(tlon, tlat)*rad2deg
             write(iunit,'(2f12.6, 2f30.15)') mesh%tr(i)%c%lon*rad2deg, mesh%tr(i)%c%lat*rad2deg, &
                  angle, length
          end if
       end do
    case default
       print*, "PLOT_CART_VECTORFIELD Warning: Vector not positioned in mesh"
       stop
    end select
    close(iunit)

    return
  end subroutine plot_cart_vectorfield

  subroutine plot_uv_vectorfield(vec, mesh)
    !----------------------------------------------
    ! PLOT_UV_VECTORFIELD
    !   Writes in 'filename' the vector field given in
    !    west-east/south-north type
    !---------------------------------------------
    type(grid_structure), intent(in) :: mesh
    type(vector_variable_uv), intent(in) :: vec
    integer(i4), parameter :: num=40962
    integer(i4):: step
    character (len=256):: filename

    !Auxiliar varibles
    integer (i4):: iunit
    integer (i4):: i

    print*, " Plotting: ", vec%name
    call getunit(iunit)
    if( len_trim(vec%name)==0) then
       filename=trim(datadir)//"vector_trv_"// &
            trim(mesh%name)//".dat"
    else
       filename=trim(datadir)//trim(vec%name)//"_"// &
            trim(mesh%name)//".dat"
    end if
    open(iunit,file=filename, status='replace')

    write(iunit, *) "#    Longitude  Latitude  Direction  Magnitude "
    select case(vec%pos)
    case(0)
       if(mesh%nv<num)then
          step=1
       else
          step=mesh%nv/num
       end if
       do i=1,mesh%nv
          if(mod(i, step)==0)then
             write(iunit,'(2f12.6, 2f30.15)') mesh%v(i)%lon*rad2deg, mesh%v(i)%lat*rad2deg, &
                  datan2(vec%u(i), vec%v(i))*rad2deg, dsqrt(vec%u(i)**2 + vec%v(i)**2)
          end if
       end do
    case(1)
       if(mesh%nt<num)then
          step=1
       else
          step=mesh%nt/num
       end if
       do i=1,mesh%nt
          if(mod(i, step)==0)then
             write(iunit,'(2f12.6, 2f32.16)') mesh%tr(i)%c%lon*rad2deg, mesh%tr(i)%c%lat*rad2deg, &
                  datan2(vec%u(i), vec%v(i))*rad2deg, dsqrt(vec%u(i)**2 + vec%v(i)**2)
          end if
       end do
    case(2)
    end select
    close(iunit)

    return
  end subroutine plot_uv_vectorfield

  subroutine plot_edge_vectorfield(vec, mesh)
    !----------------------------------------------
    ! PLOTVECTORFIELD
    !   Writes in 'filename' the Edge vectors
    !   Converts to west-east/south-north vector type
    !   to be read by gmt
    !   Plots only up to 'num' of the vectors given
    !   so that the graph is not too poluted
    !---------------------------------------------
    type(grid_structure), intent(in) :: mesh
    type(vector_variable_edge), intent(in) :: vec
    integer(i4), parameter ::  num=40962
    integer(i4):: step
    character (len=256):: filename

    !Auxiliar varibles
    real (r8):: tlat
    real (r8):: tlon
    real (r8):: v(1:3)
    integer (i4):: iunit
    integer (i4):: i

    print*, " Plotting: ", vec%name

    call getunit(iunit)
    if( len_trim(vec%name)==0) then
       filename=trim(datadir)//"vector_trv_"// &
            trim(mesh%name)//".dat"
    else
       filename=trim(datadir)//trim(vec%name)//"_"// &
            trim(mesh%name)//".dat"
    end if
    open(iunit,file=filename, status='replace')

    write(iunit, *) "#    Longitude  Latitude  Direction Length "
    if(mesh%ne<num)then
       step=1
    else
       step=mesh%ne/num
    end if
    do i=1,mesh%ne
       if(mod(i, step)==0)then
          select case(vec%pos)
          case(0)
             v=vec%nr(i)*mesh%ed(i)%nr+vec%tg(i)*mesh%ed(i)%tg
             call convert_vec_cart2sph(mesh%ed(i)%c%p, v, tlon, tlat)
             write(iunit,'(2f12.6, 2f30.15)') mesh%ed(i)%c%lon*rad2deg, &
                  mesh%ed(i)%c%lat*rad2deg, &
                  datan2(tlon, tlat)*rad2deg, dsqrt(tlat**2+tlon**2)
          case(1)
             v=vec%nr(i)*mesh%edhx(i)%nr+vec%tg(i)*mesh%edhx(i)%tg
             call convert_vec_cart2sph(mesh%edhx(i)%c%p, v, tlon, tlat)
             write(iunit,'(2f12.6, 2f30.15)') mesh%edhx(i)%c%lon*rad2deg, &
                  mesh%edhx(i)%c%lat*rad2deg, &
                  datan2(tlon, tlat)*rad2deg, dsqrt(tlat**2+tlon**2)
          end select
       end if
    end do
    close(iunit)

    return
  end subroutine plot_edge_vectorfield

  subroutine plot_grad_vectorfield(var, mesh)
    !----------------------------------------------
    ! PLOT_GRAD_VECTORFIELD
    ! Plots a gradient vector field given implicitely
    ! inside the interpolable varible (var)
    !---------------------------------------------
    type(grid_structure), intent(in) :: mesh
    type(scalar_field), intent(in) :: var

    !Temporary vector var
    type(vector_field_cart) :: vec

    vec%n=var%n
    vec%pos=var%pos
    vec%name=trim(var%name)//"_gradest"
    allocate(vec%p(1:var%n))
    vec%p=var%g

    call plot_cart_vectorfield(vec, mesh)
    return

  end subroutine plot_grad_vectorfield

end module interpack
