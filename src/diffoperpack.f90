module diffoperpack
  !========================================================================
  !
  !  DIFFOPERPACK
  !  
  !
  !	Pack with discrete differential operators routines
  !      for the geodesic sphere
  !
  !	Pedro da Silva Peixoto (pedrosp@ime.usp.br)
  !	Nov 2011
  !
  !========================================================================

  !Global kinds and constants
  use constants, only: &
       i4, &
       r8

  !Mesh and fields structures/types
  use datastruct, only: &
       grid_structure, &
       scalar_field, &
       rbf_matrix_structure, &
       vector_field_cart, &
       vector_variable_edge, &
       vectorinterpol_methods

  !Spherical mesh pack
  use smeshpack, only: &
       arclen, &
       cross_product, &
       insidetrmesh, &
       modint, &
       proj_vec_sphere

  !Interpolation pack
  use interpack, only: &
       rbfvec_matrix_build, &
       vec_remap_trc2ed, &
       vec_remap_trv2ed, &
       vec_remap_trv2trc, &
       vector_reconstruct

  implicit none

contains 

  !========================================================================
  !    Discrete integral operators
  !========================================================================

  function int_simpson(s, f)
    !-----------------------------------------------
    !  Simpsons integral rule for 3 points on an arc
    !   not necessarily equally spaced
    !   s -> relative position between the points
    !   f -> values on points
    !-----------------------------------------------

    !Points positions on arc
    real (r8), intent(in) :: s(1:3)
    !Values of the function
    real (r8), intent(in) :: f(1:3)
    !Points positions 
    real (r8):: int_simpson
    real (r8):: l
    real (r8):: l2
    real (r8):: l3

    !Aux
    real (r8):: g(1:3)
    real (r8):: a
    real (r8):: b
    real (r8):: c

    l=s(3)-s(1)
    l2=(s(3)*s(3)-s(1)*s(1))/2
    l3=( s(3)*s(3)*s(3)-s(1)*s(1)*s(1) )/3

    g(1)=f(1)/((s(1)-s(2))*(s(1)-s(3)))
    g(2)=f(2)/((s(2)-s(1))*(s(2)-s(3)))
    g(3)=f(3)/((s(3)-s(1))*(s(3)-s(2)))

    !int_{x=s1..s3} (a x^2 + b x + c)
    a=sum(g)
    b=-(g(1)*(s(2)+s(3))+ g(2)*(s(1)+s(3)) + g(3)*(s(1)+s(2)))
    c= g(1)*s(2)*s(3) + g(2)*s(1)*s(3) + g(3)*s(1)*s(2)
    !print*, "s", s
    !print*, "f", f
    !print*, "a,b,c",a, b, c
    !print*, "ls", l, l2, l3
    int_simpson=a*l3+b*l2+c*l

    return
  end function int_simpson

  function int_trp(s, f)
    !-----------------------------------------------
    !  Trapezoidal integral rule for 2 points on an arc
    !   s -> relative position between the points
    !   f -> values on points
    !-----------------------------------------------

    !Points positions on arc
    real (r8), intent(in) :: s(1:2)
    !Values of the function
    real (r8), intent(in) :: f(1:2)
    !Points positions 
    real (r8):: int_trp

    int_trp=(s(2)-s(1))*(f(1)+f(2))/2

    return
  end function int_trp


  !========================================================================
  !    Discrete DIVERGENCE operators
  !========================================================================

  function div_cell_Cgrid(i, u, mesh)
    !-----------------------------------------------------------
    ! Given a the normal components of a vector field defined on  edges
    !  calculates the numerical divergence based on Gauss Theorem
    ! Valid for C grids only
    !-----------------------------------------------------------
    !   i -> node if u%pos=3, tr if u%pos=2
    !   u -> scalarfield defined on Triangle or hexagon edges (normal comp only)
    !   mesh -> mesh structure
    !
    !-----------------------------------------------------------
    !Node for divergence calculation
    integer (i4), intent(in) :: i

    !Input variable with vectors
    type(scalar_field), intent(in) :: u

    real (r8):: div_cell_Cgrid

    !Mesh in which the variable belongs
    type(grid_structure), intent(in) :: mesh

    !Logical var for considering or not divergente average
    !That is, division or not by the cell area
    !logical, intent(in) :: average

    !Temp. vars
    real (r8):: divtmp
    real (r8):: signcor
    real (r8):: nrv
    !real (r8), allocatable :: vec(:)

    !Counters
    integer (i4):: j
    integer (i4):: k

    !Check if vectors on hexag/pent edges

    select case(u%pos)
    case(2)
       !For all edges forming the triangle
       divtmp=0_r8
       do j=1, 3
          !Get edge index
          k=mesh%tr(i)%ed(j)
          !print*, "Edge:", j
          !Get edge outer normal related to the hexagon
          signcor=real(mesh%tr(i)%nr(j), r8)
          !Calculate normal component
          nrv=signcor*u%f(k)
          !Calculate numerical integration
          divtmp=divtmp+nrv*mesh%ed(k)%leng
       end do
       !if(average)then
       !Averaged divergence (division by cell area)
       divtmp=divtmp/mesh%tr(i)%areag
       !end if
    case(3)
       !For all edges forming the hexagon
       divtmp=0._r8
       do j=1, mesh%v(i)%nnb
          !Get edge index
          k=mesh%v(i)%ed(j)
          !print*, "Edge:", j
          !Get edge outer normal related to the hexagon
          signcor=real(mesh%hx(i)%nr(j), r8)
          !Calculate normal component
          nrv=signcor*u%f(k) !dot_product(v%p(k)%v, mesh%edhx(k)%nr)
          !Calculate numerical integration
          divtmp=divtmp+nrv*mesh%edhx(k)%leng
       end do
       !if(average)then
       !Averaged divergence (division by cell area)
       divtmp=divtmp/mesh%hx(i)%areag
       !end if
    case(6)
       !For all edges forming the hexagon
       divtmp=0._r8
       do j=1, mesh%v(i)%nnb
          !Get edge index
          k=mesh%v(i)%ed(j)
          !print*, "Edge:", j
          !Get edge outer normal related to the hexagon
          signcor=real(mesh%hx(i)%ttgout(j), r8)
          !Calculate normal component
          nrv=signcor*u%f(k) !dot_product(v%p(k)%v, mesh%edhx(k)%nr)
          !Calculate numerical integration
          divtmp=divtmp+nrv*mesh%edhx(k)%leng
       end do
       !if(average)then
       !Averaged divergence (division by cell area)
       divtmp=divtmp/mesh%hx(i)%areag
       !end if
    case default
       print*, "div_cell ERROR: Vectors not on hx or tr edges (pos=2,3, 6).", u%pos
       stop
    end select

    !Assign returnable variable the result
    div_cell_Cgrid=divtmp

    return
  end function div_cell_Cgrid

  function grad_edge_Cgrid(ed, h, pos, mtd, mesh)
    !-----------------------------------------------------------
    ! Given a scalar field h calculates the gradient at edge "ed"
    ! with orietation depending on the position "pos" wanted
    !
    ! If pos=3 - normal to edhx edge is used as reference
    !   (velocities at midpoint of Voronoi cell edge)
    ! If pos=6 - tangent to ed is used as reference
    !   (velocities at intersection of triangle and Voronoi edges)
    !
    ! MTD - method to be used
    ! If mtd="trsk" - use simple differences of values on cells that
    !     share the edge (2nd order only if pos=6)
    ! If mtd="bary" - use gradient of barycentric coordinate
    !     (2nd order with pos = 3 and 6)
    !     *may affect mimetic properties
    !
    ! Valid for C grids only
    !-----------------------------------------------------------
    !Node for divergence calculation
    integer (i4), intent(in) :: ed

    !Position to put gradient
    integer (i4), intent(in) :: pos

    !Method to be used
    character (len=4), intent(in) :: mtd

    !Input variable with vectors
    type(scalar_field), intent(in) :: h

    real (r8):: grad_edge_Cgrid

    !Mesh in which the variable belongs
    type(grid_structure), intent(in) :: mesh

    !Logical var for considering or not divergente average
    !That is, division or not by the cell area
    !logical, intent(in) :: average

    !Temp. vars
    real (r8), dimension(1:3):: gradtmp
    real (r8), dimension(1:3):: vectmp
    real (r8):: signcor
    !real (r8):: nrv
    !real (r8), allocatable :: vec(:)
    integer(i4):: i, k

    !Check if vectors on hexag/pent edges

    select case(pos)
    case(3)
       if(trim(mtd)=="trsk")then
          signcor=dsign( 1._r8, dot_product(mesh%edhx(ed)%nr, &
               mesh%v(mesh%edhx(ed)%sh(2))%p-mesh%v(mesh%edhx(ed)%sh(1))%p ))
          grad_edge_Cgrid=signcor*(h%f(mesh%edhx(ed)%sh(2))-h%f(mesh%edhx(ed)%sh(1)))
          grad_edge_Cgrid=grad_edge_Cgrid/mesh%ed(ed)%leng

       elseif(trim(mtd)=="bary")then
          !Find triangle that cointains the midpoint of the edge
          if(insidetrmesh(mesh%edhx(ed)%c%p, mesh%edhx(ed)%v(1) , mesh))then
             k=mesh%edhx(ed)%v(1)
          elseif(insidetrmesh(mesh%edhx(ed)%c%p, mesh%edhx(ed)%v(2) , mesh))then
             k=mesh%edhx(ed)%v(2)
          else
             print*, "grad_edge_Cgrid error:"
             print*, "Could not fit edge center into triangles"
             stop
          end if
          !print*, ed, k

          !Calculate vector gradient of field using barycentric coords
          gradtmp=0.
          do i=1, 3
             !get vertices
             !j=mesh%tr(k)%v(i)
             !get edge indexes
             !m=modint(i+1, 3)
             !l=mesh%tr(k)%ed(m)
             !vectmp=-mesh%ed(l)%nr*mesh%tr(k)%nr(m)*mesh%ed(l)%leng*h%f(j)
             vectmp=mesh%v(mesh%tr(k)%v(modint(i+2, 3)))%p-mesh%v(mesh%tr(k)%v(modint(i+1, 3)))%p
             gradtmp=gradtmp+h%f(mesh%tr(k)%v(i))*&
                  cross_product(mesh%edhx(ed)%c%p, vectmp)/mesh%tr(k)%areap/2._r8
             !print*, mesh%tr(k)%v(modint(i+2, 3)), mesh%tr(k)%v(modint(i+1, 3))
          end do
          gradtmp=proj_vec_sphere(gradtmp, mesh%edhx(ed)%c%p)
          grad_edge_Cgrid=dot_product(gradtmp, mesh%edhx(ed)%nr)
       end if

    case(6)
       ! Obs: the tangent of a triangle edge (which is used as normal
       !   of the voronoi cell edges) is always defined such
       !   that it point from its vertex 1 to its vertex 2
       !   Therefore, the gradient needs no correction term (n_ei)
       !   NOT normalized to earth radius
       grad_edge_Cgrid=(h%f(mesh%ed(ed)%v(2))-h%f(mesh%ed(ed)%v(1)))
       grad_edge_Cgrid=grad_edge_Cgrid/mesh%ed(ed)%leng
    case default
       print*, "grad_edge_Cgrid ERROR: position not set for HX C gridon hx or tr edges (pos=3, 6)."
       stop
    end select


    return
  end function grad_edge_Cgrid

  subroutine div_mesh_hx_trp(v, div, mesh)
    !-----------------------------------------------------------
    ! Given a vector variable defined on hexagonal edges
    !  calculates the numerical divergence based on Gauss Theorem
    !-----------------------------------------------------------
    !Variable with vectors on hexagonal corners (triang. cc)
    type(vector_field_cart), intent(in) :: v

    !Output variable for divergence in scalar variable defined on
    ! hexagon centers
    type(scalar_field), intent(inout) :: div

    !Mesh in which the variable belongs
    type(grid_structure), intent(in) :: mesh

    !Temp. vars
    real (r8):: divtmp
    real (r8):: sign
    real (r8), dimension(1:3) :: nr(1:3)
    real (r8), dimension(1:3) :: f(1:2)
    real (r8), dimension(1:3) :: s(1:2)

    !Counters
    integer (i4):: i
    integer (i4):: j
    integer (i4):: ed
    integer (i4):: tr(1:2)

    div%pos=0
    div%n=mesh%nv
    if(.not.allocated(div%f))then
       allocate(div%f(1:div%n))
    end if

    !call vec_edhx2trc_interpol(ved, v, mesh)

    !For all hexagons
    do i=1, mesh%nv
       !do i=8,8
       divtmp=0
       do j=1, mesh%v(i)%nnb 

          !Get edge index
          ed=mesh%v(i)%ed(j)
          !print*,"ed: ", ed

          !Get corners (triang circumcenters)
          tr=mesh%edhx(ed)%v
          !print*,"tr: ", tr(1:2)

          !Sign correction for this hx edge
          sign=real(mesh%hx(i)%nr(j), r8)
          !print*,"sign: ", sign

          !Get normal to edge with correction
          nr=sign*mesh%edhx(ed)%nr
          !print*,"Normal: ", nr

          !Get first corner vector dot prod normal
          s(1)=0
          f(1)=dot_product(v%p(tr(1))%v,  nr)
          !print*, "V1:", v%p(tr(1))%v
          !print*, "s1, f1:", s(1), f(1)

          !Get other corner vector dot prod normal
          s(2)=mesh%edhx(ed)%leng
          f(2)=dot_product(v%p(tr(2))%v,  nr)

          divtmp=divtmp+int_trp(s, f)
          !divtmp=divtmp+(s(3)/6)*(f(1)+4*f(2)+f(3))
          !print*,"simp:", int_simpson(s, f)
       end do
       div%f(i)=divtmp/mesh%hx(i)%areag
       !print*, "Hx:", i, "Div:", divtmp
    end do

  end subroutine div_mesh_hx_trp

  subroutine div_mesh_hx_simp(ved, v, div, mesh)
    !-----------------------------------------------------------
    ! Given a vector variable defined on hexagonal edges
    !  calculates the numerical divergence based on Gauss Theorem
    !-----------------------------------------------------------
    !Input variable with vectors on hexagonal edges midpoints
    type(vector_variable_edge), intent(in) :: ved

    !Output variable for divergence in scalar variable defined on
    ! hexagon centers
    type(scalar_field), intent(inout) :: div

    !Mesh in which the variable belongs
    type(grid_structure), intent(in) :: mesh

    !Variable with vectors on hexagonal corners (triang. cc)
    type(vector_field_cart):: v

    !Temp. vars
    real (r8):: divtmp
    real (r8):: signc
    real (r8), dimension(1:3) :: nr(1:3)
    real (r8), dimension(1:3) :: f(1:3)
    real (r8), dimension(1:3) :: s(1:3)

    !Counters
    integer (i4):: i
    integer (i4):: j
    integer (i4):: ed
    integer (i4):: tr(1:2)

    div%pos=0
    div%n=mesh%nv
    if(.not.allocated(div%f))then
       allocate(div%f(1:div%n))
    end if

    !call vec_edhx2trc_interpol(ved, v, mesh)

    !For all hexagons
    do i=1, mesh%nv
       !do i=8,8
       divtmp=0
       do j=1, mesh%v(i)%nnb 

          !Get edge index
          ed=mesh%v(i)%ed(j)
          !print*,"ed: ", ed

          !Get corners (triang circumcenters)
          tr=mesh%edhx(ed)%v
          !print*,"tr: ", tr(1:2)

          !Sign correction for this hx edge
          signc=real(mesh%hx(i)%nr(j), r8)
          !print*,"sign: ", sign

          !Get normal to edge with correction
          nr=signc*mesh%edhx(ed)%nr
          !print*,"Normal: ", nr

          !Get first corner vector dot prod normal
          s(1)=0
          f(1)=dot_product(v%p(tr(1))%v,  nr)
          !print*, "V1:", v%p(tr(1))%v
          !print*, "s1, f1:", s(1), f(1)

          !Get edge's midpoint vector
          s(2)=arclen(mesh%tr(tr(1))%c%p, mesh%edhx(ed)%c%p)
          f(2)=signc*ved%nr(ed)
          !print*, "s2, f2:", s(2), f(2)

          !Get other corner vector dot prod normal
          s(3)=s(2)+arclen(mesh%edhx(ed)%c%p, mesh%tr(tr(2))%c%p)
          !print*, s(3)-s(2)-s(2)
          f(3)=dot_product(v%p(tr(2))%v,  nr)
          !print*, "V3:", v%p(tr(2))%v
          !print*, "s3, f:", s(3), f(3)

          divtmp=divtmp+int_simpson(s, f)
          !divtmp=divtmp+(s(3)/6)*(f(1)+4*f(2)+f(3))
          !print*,"simp:", int_simpson(s, f)
       end do
       div%f(i)=divtmp/mesh%hx(i)%areag
       !print*, "Hx:", i, "Div:", divtmp
    end do

  end subroutine div_mesh_hx_simp

  !  subroutine divergence(v, div, mesh, average)
  !    !-----------------------------------------------------------
  !    ! Given a vector variable defined on either triangles' edges
  !    ! or hexagonals/pentagonals' edges calculates the numerical
  !    ! divergence average or integral based on Divergence Theorem
  !    !
  !    ! Input
  !    !   v -> vector variable defined on edges
  !    !   v%pos=0 -> triangle edges
  !    !   v%pos=1 -> hexag edges
  !    !   mesh -> mesh structure
  !    !
  !    !
  !    ! average = .true. -> Divides divergence cell integral by cell area
  !    ! average = .false. -> Don't divide by cell area
  !    !
  !    ! Output
  !    !   div  -> divergence on triangle or hexagon centers
  !    !  div%pos=1 -> triangle centers
  !    !  div%pos=0 -> hexagon centers
  !    !-----------------------------------------------------------
  !    !Input variable with vectors on triangles edge's midpoints
  !    type(vector_variable_edge), intent(in) :: v
  !
  !    !Output variable for divergence in scalar variable defined on
  !    ! triangle centers
  !    type(interpolable_variable), intent(inout) :: div
  !
  !    !Mesh in which the variable belongs
  !    type(grid_structure), intent(in) :: mesh
  !
  !    !Logical var for considering or not divergente average
  !    !That is, division or not by the cell area
  !    logical, intent(in):: average
  !
  !    !Temp. vars
  !    real (r8):: divtmp
  !    real (r8):: signc
  !    real (r8):: a
  !    real (r8):: nrm
  !    !real (r8), allocatable :: vec(:)
  !
  !    !Counters
  !    integer (i4):: i
  !    integer (i4):: j
  !    integer (i4):: k
  !
  !
  !    if(v%pos==0)then
  !       !Vectors on triangles' edges
  !
  !       !Divergence positioned on triangle centres
  !       div%pos=1
  !       div%n=mesh%nt
  !
  !       !Allocate memory for divergence
  !       if(.not.allocated(div%f))then
  !          allocate(div%f(1:div%n))
  !       end if
  !
  !       !For all triangles
  !       do i=1, mesh%nt
  !          !For all edges forming the triangle
  !          divtmp=0_r8
  !          do j=1, 3
  !             !Get edge index
  !             k=mesh%tr(i)%ed(j)
  !             !Get edge outer normal related to the triangle
  !             signc=real(mesh%tr(i)%nr(j), r8)
  !             !Get the norm of the edge vector
  !             nrm=dsqrt((v%nr(k))**2+(v%tg(k))**2)
  !             !Calculate sin(alpha), angle between v%nr and edge tangent
  !             a=signc*v%nr(k)/nrm
  !             !Calculate numerical integral term
  !             divtmp=divtmp+signc*v%nr(k)*mesh%ed(k)%leng
  !          end do
  !          if(average)then
  !             !Averaged divergence (division by cell area)
  !             divtmp=divtmp/mesh%tr(i)%areag
  !          end if
  !          !Assign returnable variable the result
  !          div%f(i)=divtmp
  !
  !       end do
  !    else
  !       !Vectors on hexag/pent edges
  !
  !       !Divergence on triangle nodes/hexagon centers
  !       div%pos=0
  !       div%n=mesh%nv
  !
  !       !Allocate memory for divergence
  !       if(.not.allocated(div%f))then
  !          allocate(div%f(1:div%n))
  !       end if
  !       !For all hexagons
  !       do i=1, mesh%nv
  !          !For all edges forming the hexagon
  !          !print*, "---- Hexagon:", hx , " ----"
  !          divtmp=0_r8
  !          do j=1, mesh%v(i)%nnb
  !             !Get edge index
  !             k=mesh%v(i)%ed(j)
  !             !print*, "Edge:", j
  !             !Get the norm of the edge vector
  !             nrm=dsqrt((v%nr(k))**2+(v%tg(k))**2)
  !             !Get edge outer normal related to the hexagon
  !             signc=real(mesh%hx(i)%nr(j), r8)
  !             !Calculate sin(alpha), angle between v%nr and edge tangent
  !             a=signc*v%nr(k)/nrm
  !             !Calculate numerical integration
  !             divtmp=divtmp+signc*v%nr(k)*mesh%edhx(k)%leng
  !          end do
  !          if(average)then
  !             !Averaged divergence (division by cell area)
  !             divtmp=divtmp/mesh%hx(i)%areag
  !          end if
  !          !Assign returnable variable the result
  !          div%f(i)=divtmp
  !       end do
  !    end if
  !
  !    return
  !  end subroutine divergence

  function div_mesh_fullvec(v, mesh, cellopt, average)
    !-----------------------------------------------------------
    ! Given a vector variable (cartesian) calculates the numerical 
    ! divergence average or integral based on Divergence Theorem
    !
    ! Input
    !   v -> vector variable defined 
    !   mesh -> mesh structure
    !   cellopt -> Cell to use: "TR" or "HX"
    !
    ! Output
    !   div  -> divergence estimates
    !-----------------------------------------------------------
    !Input variable with vectors
    type(vector_field_cart), intent(in) :: v

    !Output variable for divergence in scalar variable defined on
    type(scalar_field):: div_mesh_fullvec

    !Mesh in which the variable belongs
    type(grid_structure), intent(in) :: mesh

    !Cell to be used
    character (len=2), optional :: cellopt
    character (len=2):: cell

    !Logical var for considering or not divergente average
    !That is, division or not by the cell area
    ! Default is .true.
    logical, optional :: average
    logical:: av

    !Velocity field - aux
    type(vector_field_cart):: vtrc
    type(vector_field_cart):: ved

    !Set default values for optionals
    if(.not.present(cellopt))then
       cell="hx"
    else
       cell=cellopt
    end if

    if(.not.present(average))then
       av=.true.
    else
       av=average
    end if

    if(cell=="HX" .or. cell=="hx")then

       !Calculate divergence based on the position of the given vectors

       select case (v%pos)
       case(0)  !Vectors on tr vertices
          ! HA - Estimate Divergence NICAM method
          !Interpolate vector field on triangle circumcenters
          call vec_remap_trv2trc(v, vtrc, mesh)
          !Interpolate vector field to edges midpoints
          ved%pos=3 !Inform that vectors must be interpoalated to hx edge
          call vec_remap_trc2ed(vtrc, ved, mesh)
          !Estimate numerical divergence HA - NICAM Method
          call div_mesh_ed_fullvec(ved, div_mesh_fullvec, mesh, av)

       case (1) !HB !Vector on Tr circumcenters, or hx vertices
          !Interpolate vector field to edges midpoints
          ved%pos=3 !Inform that vectors must be interpoalated to hx edge
          call vec_remap_trc2ed(v, ved, mesh)
          !Estimate numerical divergence 
          call div_mesh_ed_fullvec(ved, div_mesh_fullvec, mesh, av)

       case (3, 6) !HC !Vectors on hx edges
          !Directly estimate numerical divergence 
          call div_mesh_ed_fullvec(v, div_mesh_fullvec, mesh, av)

       case default
          print*, "div_mesh error - Don't know how to calculate with this vector position:", v%pos, cell
          stop
       end select
    elseif(cell=="TR" .or. cell=="tr")then
       !Calculate divergence based on the position of the given vectors

       select case (v%pos)
       case(0)  ! TB !Vectors on tr vertices
          !Interpolate vector field to edges midpoints
          call vec_remap_trv2ed(v, ved, mesh)

          !Estimate numerical divergence 
          call div_mesh_ed_fullvec(ved, div_mesh_fullvec, mesh, av)

       case (1) ! TA !Vector on Tr circumcenters, or hx vertices
          !Interpolate vector field to tr vertices (nodes)
          ved%pos=2
          call vec_remap_trc2ed(v, ved, mesh)
          !Estimate numerical divergence 
          call div_mesh_ed_fullvec(ved, div_mesh_fullvec, mesh, av)

       case (2) !TC !Vectors on TR edges
          !Directly estimate numerical divergence 
          call div_mesh_ed_fullvec(v, div_mesh_fullvec, mesh, av)

       case default
          print*, "div_mesh error - Don't know how to calculate with this vector position:", v%pos, cell
          stop
       end select

    else
       print*, "div_mesh error - Cell has not been properly given",  cell
       stop
    end if

    return
  end function div_mesh_fullvec

  function divho_mesh_fullvec(ved, mesh, discinterp_mtd, cellopt, average)
    !-----------------------------------------------------------
    ! Given a vector variable (cartesian) calculates the numerical
    ! divergence average or integral based on Divergence Theorem
    !
    ! Input
    !   ved -> vector field
    !   mesh -> mesh structure
    !   cellopt -> Cell to use: "HX"
    !
    ! Output
    !   divho_mesh  -> divergence estimates
    !-----------------------------------------------------------
    !Input variable with scalars on Voronoi cell edge's midpoints
    type(vector_field_cart), intent(in) :: ved

    !Output variable for divergence in scalar variable defined on
    type(scalar_field):: divho_mesh_fullvec

    !Mesh in which the variable belongs
    type(grid_structure), intent(in) :: mesh

    !Inteprolation/Reconsgtruction method used for
    !     divergence/laplacian discretizations
    type(vectorinterpol_methods) :: discinterp_mtd

    !Normal components
    type(scalar_field):: vecncomp

    !Cell to be used
    character (len=2), optional :: cellopt
    character (len=2):: cell

    !Logical var for considering or not divergente average
    !That is, division or not by the cell area
    ! Default is .true.
    logical, optional :: average
    logical:: av

    !Velocity field - aux
    type(vector_field_cart):: vtrc

    !Aux point
    real(r8), dimension(1:3)::p

    !Index
    integer (i4):: i
    integer (i4):: j
    integer (i4):: k
    integer (i4):: l
    integer (i4):: n

    !Condition numbers
    real (r8):: condnummin
    real (r8):: condnummax

    !RBF matrix vector
    type(rbf_matrix_structure), allocatable :: rbf_mat(:)

    !RBF stencil and kind of interpolation
    character(len=4):: stencil

    !Set default values for optionals
    if(.not.present(cellopt))then
       cell="hx"
    else
       cell=cellopt
    end if

    if(.not.present(average))then
       av=.true.
    else
       av=average
    end if

    if(discinterp_mtd%recon(1:3)=="rbf")then
       !If RBF, calculate the matrix structure (precalculation, depends only on mesh)
       condnummin=100000000.
       condnummax=0.
       stencil=trim(discinterp_mtd%recon(4:6))
       !Build rbf matrices
       call rbfvec_matrix_build( stencil, discinterp_mtd%rbf_par, rbf_mat, mesh)
       !Calculate condition number estimates
       n=size(rbf_mat)
       print*
       print*, "Length of RBF matrices vector:", n
       do i=1,n
          condnummin=min(condnummin,rbf_mat(i)%condnum)
          condnummax=max(condnummax,rbf_mat(i)%condnum)
       end do
       print*, "Condition Number Estimates (min, max):"
       print "(2e12.3)", condnummin, condnummax
       print*, "Shape parameter:"
       print*, discinterp_mtd%rbf_par
       print*
    end if

    if(cell=="HX" .or. cell=="hx")then

       !Calculate divergence based on the position of the given vectors
       select case (ved%pos)
       case (3) !HC !Vectors on hx edges

          !Reconstruct vector field to triangle circuncenters/
          !  Voronoi corners - based on alignement index
          vtrc%pos=1
          vtrc%n=mesh%nt
          allocate(vtrc%p(1:vtrc%n))

          !Divergence on triangle nodes/hexagon centers
          divho_mesh_fullvec%pos=0
          divho_mesh_fullvec%n=mesh%nv
          if(.not.allocated(divho_mesh_fullvec%f))then
             allocate(divho_mesh_fullvec%f(1:divho_mesh_fullvec%n))
          end if

          !Vectors normal comp on HX edges
          vecncomp%pos=3
          vecncomp%n=mesh%ne
          allocate(vecncomp%f(1:vecncomp%n))
          do i=1, mesh%ne
             vecncomp%f(i)=dot_product(ved%p(i)%v, mesh%edhx(i)%nr)
          end do
          l=0
          do i=1, mesh%nv
             !Check if cell is aligned
             if(mesh%hx(i)%align<discinterp_mtd%alignlimit)then
                !Directly estimate numerical divergence
                divho_mesh_fullvec%f(i)=div_cell_fullvec(i, ved, mesh, average)
             else
                l=l+1
                !Reconstruct to surrounding trcs
                do j=1, mesh%v(i)%nnb
                   k=mesh%v(i)%tr(j)
                   p=mesh%tr(k)%c%p
                   vtrc%p(k)%v=vector_reconstruct (p, &
                        vecncomp, mesh, discinterp_mtd%recon, rbf_mat, pindex=k)
                   !vtrc%p(k)%v=(/-p(2), p(1), 0/)
                end do
                !Calculate higher order div
                divho_mesh_fullvec%f(i)=divho_cell_fullvec(i, ved, vtrc, mesh, average)

             end if
          end do
       case default
          print*, "divho_mesh error - Don't know how to calculate with"//&
               " this vector position:", ved%pos, cell
          stop
       end select
    else
       print*, "divho_mesh error - Cell has not been properly given",  cell
       stop
    end if
    print*, "Used higher order method in ", l, " cells"
    return
  end function divho_mesh_fullvec


  subroutine div_mesh_ed_fullvec(v, div, mesh, average)
    !-----------------------------------------------------------
    ! Given a vector variable defined on either triangles' edges
    ! or hexagonals/pentagonals' edges calculates the numerical 
    ! divergence average or integral based on Divergence Theorem
    !
    ! Input
    !   v -> cartesian vector variable defined on tr or hx edges
    !   v%pos=2 -> triangle edges
    !   v%pos=3 -> hexag edges
    !   mesh -> mesh structure
    !   
    ! average = .true. -> Divides divergence cell integral by cell area
    ! average = .false. -> Don't divide by cell area
    !
    ! Output
    !   div  -> divergence on triangle or hexagon centers 
    !  div%pos=1 -> triangle centers
    !  div%pos=0 -> hexagon centers
    !-----------------------------------------------------------
    !Input variable with vectors on triangles edge's midpoints
    type(vector_field_cart), intent(in) :: v

    !Output variable for divergence in scalar variable defined on
    ! triangle centers
    type(scalar_field), intent(inout) :: div

    !Mesh in which the variable belongs
    type(grid_structure), intent(in) :: mesh

    !Logical var for considering or not divergente average
    !That is, division or not by the cell area
    logical, intent(in):: average

    !Counters
    integer (i4):: i

    if(v%pos==2)then
       !Vectors on triangles' edges

       !Divergence positioned on triangle centres
       div%pos=1
       div%n=mesh%nt

       !Allocate memory for divergence
       if(.not.allocated(div%f))then
          allocate(div%f(1:div%n))
       end if

       !For all triangles
       do i=1, mesh%nt
          !Assign returnable variable the result
          div%f(i)=div_cell_fullvec(i, v, mesh, average)
       end do
    else
       !Vectors on hexag/pent edges

       !Divergence on triangle nodes/hexagon centers
       div%pos=0
       div%n=mesh%nv

       !Allocate memory for divergence
       if(.not.allocated(div%f))then
          allocate(div%f(1:div%n))
       end if
       !For all hexagons
       do i=1, mesh%nv
          !Assign returnable variable the result
          div%f(i)=div_cell_fullvec(i, v, mesh, average)
       end do
    end if

    return
  end subroutine div_mesh_ed_fullvec

  function div_cell_fullvec(i, v, mesh, average)
    !-----------------------------------------------------------
    ! Given a cartesian vector variable the numerical 
    ! divergence average or integral based on Divergence Theorem
    !
    !   i -> node if v%pos=3, tr if v%pos=2
    !   v -> vector variable defined on Triangle or hexagon edges !!!
    !   mesh -> mesh structure
    !
    !-----------------------------------------------------------
    !Node for divergence calculation
    integer (i4), intent(in) :: i

    !Input variable with vectors 
    type(vector_field_cart), intent(in) :: v

    real (r8):: div_cell_fullvec

    !Mesh in which the variable belongs
    type(grid_structure), intent(in) :: mesh

    !Logical var for considering or not divergente average
    !That is, division or not by the cell area
    logical, intent(in) :: average

    !Temp. vars
    real (r8):: divtmp
    real (r8):: signcor
    real (r8):: nrv
    !real (r8), allocatable :: vec(:)

    !Counters
    integer (i4):: j
    integer (i4):: k

    !Check if vectors on hexag/pent edges

    select case(v%pos)
    case(2)
       !For all edges forming the triangle
       divtmp=0_r8
       do j=1, 3
          !Get edge index
          k=mesh%tr(i)%ed(j)
          !print*, "Edge:", j
          !Get edge outer normal related to the hexagon
          signcor=real(mesh%tr(i)%nr(j), r8)
          !Calculate normal component
          nrv=signcor*dot_product(v%p(k)%v, mesh%ed(k)%nr)
          !Calculate numerical integration
          divtmp=divtmp+nrv*mesh%ed(k)%leng
       end do
       if(average)then
          !Averaged divergence (division by cell area)
          divtmp=divtmp/mesh%tr(i)%areag
       end if
    case(3)
       !For all edges forming the hexagon
       divtmp=0._r8
       do j=1, mesh%v(i)%nnb 
          !Get edge index
          k=mesh%v(i)%ed(j)
          !print*, "Edge:", j
          !Get edge outer normal related to the hexagon
          signcor=real(mesh%hx(i)%nr(j), r8)
          !Calculate normal component
          nrv=signcor*dot_product(v%p(k)%v, mesh%edhx(k)%nr)
          !Calculate numerical integration
          divtmp=divtmp+nrv*mesh%edhx(k)%leng
       end do
       if(average)then
          !Averaged divergence (division by cell area)
          divtmp=divtmp/mesh%hx(i)%areag
       end if
    case(6)
       !For all edges forming the hexagon
       divtmp=0._r8
       do j=1, mesh%v(i)%nnb
          !Get edge index
          k=mesh%v(i)%ed(j)
          !print*, "Edge:", j
          !Get edge outer normal related to the hexagon
          signcor=real(mesh%hx(i)%ttgout(j), r8)
          !Calculate normal component
          nrv=signcor*dot_product(v%p(k)%v, mesh%ed(k)%tg)
          !Calculate numerical integration
          divtmp=divtmp+nrv*mesh%edhx(k)%leng
       end do

       if(average)then
          !Averaged divergence (division by cell area)
          divtmp=divtmp/mesh%hx(i)%areag
       end if
    case default
       print*, "div_cell ERROR: Vectors not on hx or tr edges (pos=2,3, 6).", v%pos
       stop
    end select

    !Assign returnable variable the result
    div_cell_fullvec=divtmp    

    return
  end function div_cell_fullvec

  function divho_cell_fullvec(i, ved, vtrc, mesh, average)
    !-----------------------------------------------------------
    ! Given a cartesian vector variable the numerical
    ! divergence average or integral based on Divergence Theorem
    !
    !   i -> node (v%pos=3)
    !   ved -> vector variable defined on hexagon edges !!!
    !   vtrc -> vector variable defined on triangle circuncenters !!!
    !   mesh -> mesh structure
    !
    !-----------------------------------------------------------
    !Node for divergence calculation
    integer (i4), intent(in) :: i

    !Input variable with vectors
    type(vector_field_cart), intent(in) :: ved
    type(vector_field_cart), intent(in) :: vtrc

    real (r8):: divho_cell_fullvec

    !Mesh in which the variable belongs
    type(grid_structure), intent(in) :: mesh

    !Logical var for considering or not divergente average
    !That is, division or not by the cell area
    logical, intent(in) :: average

    !Counters
    integer (i4):: j
    !integer (i4):: k
    integer (i4):: ed
    integer (i4):: tr(1:2)

    !Temp. vars
    real (r8):: divtmp
    real (r8):: signc
    !real (r8):: nrv
    real (r8), dimension(1:3) :: nr(1:3)
    real (r8), dimension(1:3) :: f(1:3)
    real (r8), dimension(1:3) :: s(1:3)


    !Check if vectors on hexag/pent edges
    if(ved%pos==3)then
       !For all edges forming the hexagon
       divtmp=0._r8
       do j=1, mesh%v(i)%nnb

          !Get edge index
          ed=mesh%v(i)%ed(j)
          !print*,"ed: ", ed

          !Get corners (triang circumcenters)
          tr=mesh%edhx(ed)%v
          !print*,"tr: ", tr(1:2)

          !Sign correction for this hx edge
          signc=real(mesh%hx(i)%nr(j), r8)
          !print*,"sign: ", sign

          !Get normal to edge with correction
          nr=signc*mesh%edhx(ed)%nr
          !print*,"Normal: ", nr

          !Get first corner vector dot prod normal
          s(1)=0
          f(1)=dot_product(vtrc%p(tr(1))%v,  nr)
          !print*, "V1:", v%p(tr(1))%v
          !print*, "s1, f1:", s(1), f(1)

          !Get edge's midpoint vector
          s(2)=arclen(mesh%tr(tr(1))%c%p, mesh%edhx(ed)%c%p)
          f(2)=dot_product(ved%p(ed)%v, nr)
          !print*, "s2, f2:", s(2), f(2)

          !Get other corner vector dot prod normal
          s(3)=s(2)+arclen(mesh%edhx(ed)%c%p, mesh%tr(tr(2))%c%p)
          !print*, s(3)-s(2)-s(2)
          f(3)=dot_product(vtrc%p(tr(2))%v,  nr)
          !print*, "V3:", v%p(tr(2))%v
          !print*, "s3, f:", s(3), f(3)

          divtmp=divtmp+int_simpson(s, f)
          !divtmp=divtmp+(s(3)/6)*(f(1)+4*f(2)+f(3))
          !print*,"simp:", int_simpson(s, f)
       end do
       if(average)then
          !Averaged divergence (division by cell area)
          divtmp=divtmp/mesh%hx(i)%areag
       end if
    else
       print*, "div_cell ERROR: Vectors not on hx edges (pos=3).", ved%pos
       stop
    end if

    !Assign returnable variable the result
    divho_cell_fullvec=divtmp

    return
  end function divho_cell_fullvec

end module diffoperpack
