module smeshpack
  !========================================================================
  !
  !  SMESHPACK
  !
  !
  !	Pack for mesh generation, structuring and handling for the sphere
  !
  !     Includes: 
  !     - Mesh generation of icosahedral, ocatahedral and random nodes mesh
  !     - Mesh optimizations: spring dynamics and scvt including local refinement
  !     - Spherical geometry tools
  !     - Mesh restructuring
  !     - Mesh search tools
  !     - Lots of others ...
  !
  !	Pedro da Silva Peixoto (pedrosp@ime.usp.br)
  !	March 2013 - last updated
  !
  !========================================================================

  !Global constants
  use constants, only: &
       datadir, &
       deg2rad, &
       eps, &
       eps2, &
       griddir, &
       i4, &
       pardir, &
       pi, &
       pi2, &
       pio2, &
       r8, &
       r16, &
       rad2deg, &
       showonscreen, &
       simulcase

  !Data structures
  use datastruct, only: &
       grid_structure, &
       point_structure, &
       vector, &
       vector_field_cart

  !Spherical triangulation pack
  use stripackm, only: &
       trmesh

  implicit none

contains 

  !=================================================================================
  !    MAIN ROUTINES FOR THE GEODESIC MESH CONSTRUCTION
  !=================================================================================

  subroutine meshbuild(mesh)
    !---------------------------------------------------------------------
    !
    !	MESHBUILD
    !
    !	Creates or loads a mesh structure with parameters given implicitly by
    !
    !	mesh%kind="icos", "octg", "rand", "read" -> Mesh basis
    !          -icos: Icosahedron
    !		       -octg: Octahedron
    !          -rand: Random points - be careful ...
    !          -read: Read cartesian coords of nodes from file
    !
    !   mesh%pos = "eqs", "pol", "ran", "ref"   -> Mesh position
    !          -eqs: Equator Symmetric
    !          -pol: Nodes on poles
    !          -ran: Random nodes
    !          -ref: Random nodes with local refinement
    !
    !   mesh%optm = "nopt", "sprg", "scvt", "salt", "hr95" -> Mesh optimization
    !		       -nopt: No optimization
    !          -sprg: Spring Dynamics
    !          -scvt: Spherical Centroidal Voronoi Tesselation
    !          -hr95: HR95 edge alignment using Miura 2005 algorithm (has bug)
    !          -salt: Spherical Aligend Tesselation (Not a good idea)
    !
    !	mesh%nv  -> number of nodes
    !
    !	mesh%loadable -> 1 or 0, where 1 indicates to load and 0 not to load
    !          the grid structure from files
    !
    !---------------------------------------------------------------------
    type(grid_structure), intent(inout) :: mesh

    !Auxiliar variables
    character (len=256):: header
    logical:: ifile
    logical:: opt
    integer (i4):: nv
    integer (i4):: glevel
    integer (i4):: i
    integer (i4):: n
    integer (i4):: level
    integer (i4):: levelload
    integer (i4):: stat
    real(r8), allocatable :: x(:)
    real(r8), allocatable :: y(:)
    real(r8), allocatable :: z(:)

    !Adjust the total number of grid point if needed
    nv=mesh%nv
    call nodesizeadj(nv, glevel, mesh%kind)

    !Give a name to the grid
    mesh%glevel=glevel
    mesh%nv=nv
    call namegrid(mesh)

    !Zero flag for delaunay triangulation
    ! 0 => nodes must be triangulated
    mesh%deltriread=0

    !Check if file to be read exists
    if(trim(mesh%kind)=="read")then
       header=trim(griddir)//trim(mesh%name)
       inquire(file=header, exist=ifile)
       if(.not.ifile)then
          print*, "meshbuild error: cannot find mesh data file to read"
          print*, trim(mesh%name)
          stop
       end if
       !print*, "Mesh before:", mesh%name, len(trim(mesh%name))
       mesh%filename=trim(mesh%name)
       n=len(trim(mesh%name))-4
       if(trim(mesh%optm) /= "nopt")then
          mesh%name=trim(mesh%name(1:n))//"_"//trim(mesh%optm)
       else
          mesh%name=trim(mesh%name(1:n))
       end if
       !print*, "After: ", mesh%name
    end if

    !Set header filename and check existence
    header=trim(griddir)//trim(mesh%name)//"_header.dat"
    inquire(file=header, exist=ifile)

    if(mesh%loadable==1 .and. ifile) then
       !------------------------------------------------
       !Load the whole mesh
       !------------------------------------------------
       call meshload(mesh, header)
       ! Check if read correctly
       !call meshprint(mesh)
       !Save gmt mesh plotting files
       !call meshgmt(mesh)
    else
       !------------------------------------------------
       !Generate mesh
       !------------------------------------------------
       if(trim(mesh%kind)=="icos" .and. mesh%hrchy>0)then
          !Hierachical grid

          !Allocate space for all points
          allocate(mesh%v(1:nv))
          allocate(x(1:nv),y(1:nv),z(1:nv))

          !Check for sub grids that already exist
          levelload=-1
          mesh%nv=12
          if(mesh%loadable==1)then
             !Check the highest level that already exists
             do level=glevel, 0, -1
                mesh%glevel=level
                !Set header filename and check existence
                call namegrid(mesh)
                header=trim(griddir)//trim(mesh%name)//"_header.dat"
                inquire(file=header, exist=ifile)
                if(ifile)then
                   mesh%nv=10*2**(2*level)+2
                   levelload=level
                   exit
                end if
             end do
          end if
          !print*, mesh%nv, mesh%glevel, header

          if(levelload>=0)then !Load a mesh
             !Load mesh
             call meshload(mesh, header)

             !Save a local copy of node positions
             do i=1,mesh%nv
                x(i)=mesh%v(i)%p(1)
                y(i)=mesh%v(i)%p(2)
                z(i)=mesh%v(i)%p(3)
                !print*, mesh%v(i)%p
             end do
             !print*, mesh%nv, mesh%glevel

          else !Create basic icosahedral grid
             levelload=0
             mesh%glevel=0
             mesh%nv=12 !Number of nodes

             print*, "------------------------------------------------------"
             print*, "    Basic mesh is icosahedral with ", mesh%nv, " nodes"
             print*, "------------------------------------------------------"
             print*

             !Get basic coords
             print*, " Mesh position: ", mesh%pos
             call icosanodes0(x(1:12), y(1:12), z(1:12), mesh%pos)


             do i=1,mesh%nv
                mesh%v(i)%p(1)=x(i)
                mesh%v(i)%p(2)=y(i)
                mesh%v(i)%p(3)=z(i)
                !print*, mesh%v(i)%p
             end do

             !Do a complete mesh structuring
             opt=.false.
             call meshstruct(mesh, opt, stat)

             call meshstore(mesh, header)
          end if

          !Construct higher levels on top of the basic grid
          do level=levelload+1, glevel

             !Calculates refinements puting points on edge's midpoints
             call meshrefin(mesh%nv, x, y, z, nv)
             !print*, mesh%nv, mesh%glevel, nv

             !Save level number
             mesh%glevel=level
             call namegrid(mesh)

             print*, "------------------------------------------------"
             print*, "Building mesh with ", mesh%nv, " nodes"
             print*, "------------------------------------------------"
             print*

             do i=1,mesh%nv
                mesh%v(i)%p(1)=x(i)
                mesh%v(i)%p(2)=y(i)
                mesh%v(i)%p(3)=z(i)
                !print*, mesh%v(i)%p
             end do

             !Perform possible optimization
             !if(.not.(trim(mesh%optm)=="nopt") .and. level==glevel)then
             if(.not.(trim(mesh%optm)=="nopt") )then
                !Check multiple optimizations
                select case(trim(mesh%optm))
                case("scvt+hr95")
                   !SCVT optimization
                   mesh%optm="scvt"
                   print*, "Applying SCVT mesh optimization..."
                   call optimizegrid(mesh)

                   !Optimize grid with Edge Alignment
                   mesh%optm="hr95"
                   print*, "Applying Edge Alignment mesh optimization..."
                   call optimizegrid(mesh)

                   mesh%optm="scvt+hr95"
                case default
                   !Optimize grid
                   call optimizegrid(mesh)
                end select
             end if

             !Do a complete mesh structuring
             opt=.false.
             call meshstruct(mesh, opt, stat)

             !Store mesh
             header=trim(griddir)//trim(mesh%name)//"_header.dat"
             call meshstore(mesh, header)
          end do

       else

          !Generate full grid and then apply possible optimization
          select case(trim(mesh%kind))
          case("rand") !Random nodes
             call randnodes(mesh)
          case("icos") !Icosahedral nodes
             call icosanodes(mesh)
          case("octg") !Octahedral nodes
             call octgnodes(mesh)
          case("read") !Read nodes from file
             call meshread(mesh)
          case default
             print*, "MESH BUILD ERROR: Invalid mesh kind : ", mesh%kind
             stop
             return
          end select

          !------------------------------------------------
          !Structure the whole mesh
          !------------------------------------------------
          mesh%nt = 2 * mesh%nv - 4  !Number of triangles
          mesh%ne = 3 * mesh%nv - 6  !Number of edges

          !Print some mesh caracteristics on screen
          call printmesh(mesh)

          !In case loadability wanted
          if(mesh%loadable==1) print '(a)', "   But no grid files to load"
          print*

          !Perform possible optimization
          if(.not.(trim(mesh%optm)=="nopt"))then
             !Check multiple optimizations
             select case(trim(mesh%optm))
             case("scvt+hr95")
                !SCVT optimization
                mesh%optm="scvt"
                print*, "Applying SCVT mesh optimization..."
                call optimizegrid(mesh)

                !Optimize grid with Edge Alignment
                mesh%optm="hr95"
                print*, "Applying HR95 Miura's mesh optimization..."
                call optimizegrid(mesh)

                mesh%optm="scvt+hr95"
             case default
                !Optimize grid
                call optimizegrid(mesh)
             end select
          end if

          !Do a complete mesh structuring
          opt=.false.
          call meshstruct(mesh, opt, stat)

          call meshstore(mesh, header)

       end if

    end if

    print*, "Mesh created or loaded: ", mesh%name
    print*,"-------------------------------------------------"
    print*

    return
  end subroutine meshbuild


  subroutine meshstruct(mesh, opt, stat)
    !--------------------------------------------------------------------
    ! MESHSTRUCT
    !
    !  Subroutine for mesh structuring
    !     This routines creates all needed relations between
    !     vertices, triangles, edges, etc ...
    !
    !  OPT = TRUE - Does only some important parts of mesh structuring
    !        FALSE - Structures and calculates all mesh properties
    !  MESH = Receives mesh parameters implicitly, and returns
    !        all mesh structures also implicitly
    !        Must receive:
    !        - mesh%nv/ - number of vertices
    !        - mesh%v   - vertices positions
    !---------------------------------------------------------------------
    !Mesh structure
    type(grid_structure), intent(inout) :: mesh

    !Logical flag to calculate remaining things only on last
    ! structuring step.
    ! opt=true : do only necessary things
    ! opt=false : do all
    logical, intent(in) :: opt

    !Flag for error/warnings
    ! 0 - no problems found
    ! 1 - mesh need restructuring
    integer (i4), intent(out) :: stat

    !Delaney Triangulation List (from stripack)
    !    n1 = list(lend(node)) is the first neghbour of node
    !    n2 = list(lptr(lend(node))) is the next neghbour
    !    n3 = list(lptr(lptr(lend(node)))) is the next and so on until
    !         nk=n1
    integer(i4), dimension(:), allocatable:: list
    integer(i4), dimension(:), allocatable:: lptr
    integer(i4), dimension(:), allocatable:: lend
    integer(i4), dimension(:), allocatable:: near
    integer(i4), dimension(:), allocatable:: next
    integer(i4), dimension(:), allocatable:: lnbtmp
    real(r8), dimension(:), allocatable:: dist
    integer(i4):: lnew

    !Flags
    logical:: isintr

    !Nodes data
    real (r8), dimension(:),allocatable :: x
    real (r8), dimension(:),allocatable :: y
    real (r8), dimension(:),allocatable :: z

    !Triangle circumcenter and vertices coordinate
    real (r8), dimension (1:3) :: c
    real (r8), dimension (1:3) :: p1
    real (r8), dimension (1:3) :: p2
    real (r8), dimension (1:3) :: p3

    !Triangle areas for the unit sphere (min, max, average, temporary)
    real (r8):: area
    real (r8):: angles(1:3)

    !Counts of neighbours, nodes, triangles, ...
    integer (i4):: nb
    integer (i4):: lp
    integer (i4):: lpl
    integer (i4):: minind
    integer (i4):: ind(1:3)
    integer (i4):: n
    integer (i4):: i
    integer (i4):: j
    integer (i4):: k
    integer (i4):: l
    integer (i4):: m
    integer (i4):: i1
    integer (i4):: i2
    integer (i4):: i3
    integer (i4):: j1
    integer (i4):: j2
    integer (i4):: iwarn1
    integer (i4):: iwarn2
    integer (i4):: iwarn3

    !Others
    real (r8)::  disttmp
    real (r8):: costheta
    real (r8):: nr(1:3)
    real (r8):: dp
    integer:: ier
    integer:: status

    !Initialize flags or erros and warnings
    ier=0
    iwarn1=0
    iwarn2=0
    iwarn3=0
    stat=0

    !------------------------------------------
    ! Set temporary cartesian coordinates of the nodes
    !------------------------------------------

    n=mesh%nv
    mesh%nt = 2 * mesh%nv - 4  !Number of triangles
    mesh%ne = 3 * mesh%nv - 6  !Number of edges
    !print*, n, mesh%nt, mesh%ne

    if(mesh%deltriread==0)then
       !Temporary lists
       allocate(x(1:n),stat=ier)
       allocate(y(1:n),stat=ier)
       allocate(z(1:n),stat=ier)
       if(ier/=0) stop "MESHSTRUCT ERROR: Allocation problem on nodes!"
       do i=1, n
          !Save a temporary list
          x(i)=mesh%v(i)%p(1)
          y(i)=mesh%v(i)%p(2)
          z(i)=mesh%v(i)%p(3)
          !For debuging
          !print'(i8,3f24.14)', i, x(i), y(i), z(i)
       end do
    end if
    !print*

    !------------------------------------------------
    ! Convert nodes from cartesian to spherical coords
    !  Pre-processing feature
    !-----------------------------------------------
    if(.not.opt)then !Do only in case of full structuring
       do i=1,mesh%nv
          !Convert to spherical/geographic coords
          call cart2sph( mesh%v(i)%p(1), mesh%v(i)%p(2), mesh%v(i)%p(3), &
               mesh%v(i)%lon, mesh%v(i)%lat)
       end do
    end if

    !------------------------------------------------
    ! Genarate Delaney triangulation with stripack
    !-------------------------------------------------
    if(.not. opt .and. showonscreen .and. mesh%deltriread==0)then
       print*,"Generating Delauney triangulation...."
       print*
    end if
    ier=0

    if(mesh%deltriread==0)then
       !Node's neighbourhood list
       allocate(list(3*mesh%nt),stat=ier)
       allocate(lptr(3*mesh%nt),stat=ier)
       allocate(lend(mesh%nv),stat=ier)
       allocate(near(mesh%nv),stat=ier)
       allocate(next(mesh%nv),stat=ier)
       allocate(dist(mesh%nv),stat=ier)
       if(ier/=0) stop "MESHSTRUCT ERROR: Allocation problem on neighbours list!"

       !Stripacks triangulation function
       call trmesh ( n, x, y, z, list, lptr, lend, lnew, near, next, dist, ier)
       if ( ier == -2 ) then
          print*,'STRIPACK TRMESH ERROR: ', ier
          print*,'  The first three nodes are collinear.'
          print*,'  Consider re-ordering the nodes   '
          stop
       else if ( ier > 0 ) then
          print*, 'STRIPACK TRMESH ERROR: ', ier
          print*, '  Duplicate nodes encountered.'
          stop
       end if
       deallocate(near)
       deallocate(next)
       deallocate(dist)
       deallocate(x)
       deallocate(y)
       deallocate(z)
    end if

    !------------------------------------------------------
    ! Restructure triangulation data for our mesh structure
    !------------------------------------------------------
    if(.not.opt .and. mesh%deltriread==0)then
       print*, "Reestructuring triangulation ..."
       print*
    end if


    if(mesh%deltriread==0)then
       !Temporary neighbour list
       allocate(lnbtmp(1:n))

       !For each node on the mesh
       do i=1, mesh%nv
          !Get the position on the list for the last node neighbour
          lpl = lend(i)
          !Get the first neighbour index
          lp=lptr(lpl)
          nb=0
          !Get all neighbours
          do j=1, n
             lnbtmp(j)=list(lp)
             nb=nb+1
             if(lp==lpl) exit
             lp=lptr(lp)
          end do

          !Allocate (or reallocate) neighbourhood lists
          if(allocated(mesh%v(i)%nb))then
             !Check if number of neighbours changed
             ! and check if the neighbours changed
             if(mesh%v(i)%nnb /= nb .or. &
                  (maxval(mesh%v(i)%nb(1:mesh%v(i)%nnb)) /= maxval(lnbtmp(1:nb))) .or. &
                  (minval(mesh%v(i)%nb(1:mesh%v(i)%nnb)) /= minval(lnbtmp(1:nb))) )then
                !Update warnings flag
                iwarn2=iwarn2+1
                !if(iwarn2==1 .and. (.not. trim(mesh%optm) == "scvt"))then
                !   print*, "Meshstruct warning: Mesh changed on node", i
                !   !print*, "  Neib changed :" !, mesh%v(i)%nnb, nb
                !   print*, "    Old Nb:", mesh%v(i)%nb(1:mesh%v(i)%nnb)
                !   print*, "    New Nb:", lnbtmp(1:nb)
                !end if
             end if

             !Reallocate in case the number of neighbours changed
             if(mesh%v(i)%nnb /= nb )then
                !Update warnings flag
                iwarn1=iwarn1+1
                !Print a Warning in case it is happening a change in
                !  number of neighbours
                !if(iwarn1==1.and. (.not. trim(mesh%optm) == "scvt"))then
                !   print*, "WARNING: Mesh changed on node", i, &
                !        "   Num of neighbours changed:", mesh%v(i)%nnb, nb
                !end if
                !Resize array (loses data)
                call reallocint(mesh%v(i)%nb, nb)
             end if
          else
             !If it is the first time structuring the mesh
             allocate(mesh%v(i)%nb(1:nb))
          end if

          !Save neighbourhood structure in mesh
          do j=1, nb !For each neighbour
             !Save neighbour index
             mesh%v(i)%nb(j)=lnbtmp(j)
          end do
          mesh%v(i)%nnb=nb
       end do
       deallocate(lnbtmp)
       deallocate(list)
       deallocate(lptr)
       deallocate(lend)

       !Print warnings if needed
       if(iwarn1>0)then
          print*, "Meshstruct warning: Num nnb changes= ", iwarn1
          stat=1
       end if
       if(iwarn2>0 .and. (.not. trim(mesh%optm) == "scvt"))then
          print*, "Meshstruct warning: Num nb struc changes = ", iwarn2
          stat=1
       end if
    end if

    !For each node on the mesh calculate other informations
    ! such as distance metrics
    disttmp=0.
    mesh%minvdist=100000.
    mesh%maxvdist=0.
    mesh%meanvdist=0.
    mesh%maxvnb=0
    do i=1,mesh%nv
       !print*, mesh%v(i)%nnb, mesh%v(i)%nb(1:mesh%v(i)%nnb)
       !Allocate (or reallocate) neighbourhood lists
       call reallocreal(mesh%v(i)%nbd,mesh%v(i)%nnb)
       call reallocreal(mesh%v(i)%nbdg,mesh%v(i)%nnb)

       !Update max number of neighbours
       mesh%maxvnb=max(mesh%v(i)%nnb,  mesh%maxvnb)

       do j=1, mesh%v(i)%nnb !For each neighbour

          !Calculate distance to neighbours
          ! This distance function is actually the -cos(angle), where angle is the angle between
          ! the 2 nodes, that is, -dot_product(node, neighbour) this distance is a increasing
          ! function of the distance, belongs to [-1,1]. Done for interpolation purpuses.
          k=mesh%v(i)%nb(j)
          costheta=dot_product(mesh%v(i)%p, mesh%v(k)%p)
          mesh%v(i)%nbd(j)=-costheta

          !Calculates geodesical distance in radians
          mesh%v(i)%nbdg(j)=arclen(mesh%v(i)%p,mesh%v(k)%p)

          !Store max, min distances between 2 nodes
          !disttmp=arcdistll(mesh%v(i)%lon,mesh%v(i)%lat, &
          !     mesh%v(k)%lon,mesh%v(k)%lat)
          mesh%minvdist=min(mesh%minvdist,mesh%v(i)%nbdg(j))
          mesh%maxvdist=max(mesh%maxvdist,mesh%v(i)%nbdg(j))
          mesh%meanvdist=mesh%meanvdist+mesh%v(i)%nbdg(j)
       end do
    end do
    mesh%meanvdist=(mesh%meanvdist/mesh%ne)/2._r8

    if(.not.opt .and. showonscreen )then
       !Print statistics
       print '(a,i4)'," Maximum number of neighbours per node: ", mesh%maxvnb
       print '(a,f8.4)'," Minimum arc distance (in degrees) between 2 nodes: ", &
            mesh%minvdist*rad2deg
       print '(a,f8.4)'," Maximum arc distance (in degrees) between 2 nodes: ", &
            mesh%maxvdist*rad2deg
       print '(a,f8.4)'," Mean arc distance (in degrees) between 2 nodes:    ", &
            mesh%meanvdist*rad2deg
       print '(a,f8.4)'," Min/Max arc distance :                             ", &
            mesh%minvdist/mesh%maxvdist
       print*
    end if

    !---------------------------------------------------
    ! Define edge indexes and properties
    !---------------------------------------------------
    if(.not.opt .and. showonscreen )then
       print*, "Defining edges ..."
       print*
    end if

    !Allocate space
    if(.not.allocated(mesh%ed))then
       allocate(mesh%ed(1:mesh%ne))
    else
       !If the new array has a different size, reallocate
       !print*, ubound(mesh%ed, 1), mesh%ne
       if(ubound(mesh%ed, 1)/=mesh%ne)then
          print*, "Meshstruct warning: Num edges changed:", &
               ubound(mesh%ed, 1), mesh%ne
          deallocate(mesh%ed)
          allocate(mesh%ed(1:mesh%ne), stat=status)
       end if
    end if

    l=1 !Edge counter
    do i=1, mesh%nv
       !Allocate (or reallocate) edges surrounding node
       call reallocint(mesh%v(i)%ed, mesh%v(i)%nnb)
       !For each neighbour
       do j=1, mesh%v(i)%nnb
          k=mesh%v(i)%nb(j)
          if(k>i)then !This is the case where an edge is to be setted
             !print*, "Edge:", l, "Nodes:", i, k
             !Save the edge index in the vertice structure
             mesh%v(i)%ed(j)=l
             !The two endpoints of the edge
             mesh%ed(l)%v(1)=i
             mesh%ed(l)%v(2)=k

             !The edge's midpoint
             mesh%ed(l)%c=midpoint_ed(i, k, mesh)
             !Vector normal to edge - points to the left of the vector
             !  v(i) --> v(k)
             mesh%ed(l)%nr=cross_product(mesh%v(i)%p, mesh%v(k)%p)
             mesh%ed(l)%nr=mesh%ed(l)%nr/norm(mesh%ed(l)%nr)
             !Tangent vector - Points in the direction from v(i) --> v(k)
             mesh%ed(l)%tg=cross_product(mesh%ed(l)%nr, mesh%ed(l)%c%p)
             mesh%ed(l)%tg=mesh%ed(l)%tg/norm(mesh%ed(l)%tg)

             !Edge euclidian length (straight line)
             mesh%ed(l)%lenp=norm(mesh%v(i)%p-mesh%v(k)%p)
             !Edge geodesical length (curved line)
             mesh%ed(l)%leng=arclen(mesh%v(i)%p,mesh%v(k)%p)
             !Zero triangles index that share the edge (for future use)
             mesh%ed(l)%sh(1:2)=0
             l=l+1
          else !This is the case where an edge has already been set
             !We need to define it on the correct position
             do m=1, mesh%v(k)%nnb
                if(mesh%v(k)%nb(m)==i) then
                   mesh%v(i)%ed(j)=mesh%v(k)%ed(m)
                   exit
                end if
             end do
          end if
       end do
    end do
    if(l-1/=mesh%ne)then
       print*, "MESHSTRUCT ERROR: Calculated number of edges differ from iterated",  mesh%ne, l-1
       stop
    end if

    !---------------------------------------------------
    !  -Triangle vertices
    !  -Triangle neighbour triangles
    !  -Triangles sharing edge
    !  -Triangle vs nodes relation
    !  -Set triangle outer normal and tangent vector indicators
    !---------------------------------------------------
    if(.not.opt .and. showonscreen )then
       print*, "Defining triangles and their relations ..."
       print*
    end if
    !Allocate space
    if(.not.allocated(mesh%tr))then
       allocate(mesh%tr(1:mesh%nt))
    else
       !If the new array has a different size, reallocate
       if(ubound(mesh%tr, 1)/=mesh%nt)then
          print*, "Meshstruct warning: Num triang changed:", &
               ubound(mesh%tr, 1), mesh%nt
          deallocate(mesh%tr)
          allocate(mesh%tr(1:mesh%nt), stat=status)
       end if
    end if

    !Define triangles indexes
    k = 1     !Index used to count triangle
    !For each node on the mesh
    do i=1,mesh%nv
       !print*, "Vertice:", i
       nb=mesh%v(i)%nnb

       !Allocate/reallocate list of neighbour triangles
       call reallocint(mesh%v(i)%tr, nb)

       !Set first node as the actual node
       ind(1)=i
       !Get all node neighbours from list
       do j=1, nb
          ind(2)=mesh%v(i)%nb(j)
          ind(3)=mesh%v(i)%nb(modint(j+1, nb))

          minind=minval(ind(1:3))
          !print*, "   Neib:", j, "Indexes", ind(1:3), "minind", minind
          !Create a triangle if it doesn't exist
          if( minind == i)then
             !print*, "       New TR:", k
             !Set triangle for this node
             mesh%v(i)%tr(j)=k
             !Save nodes forming this triangle
             mesh%tr(k)%v(1:3)=ind(1:3)
             !Update counter
             k=k+1
          else !If it already exists
             !We need to define it based on previous settings
             do l=1, mesh%v(minind)%nnb
                !There are 2 triangles sharing the edge i -- minind
                !We check which is it based on the positions of i1, i2, i3
                if(mesh%v(minind)%nb(l)==i.and.minind==ind(3)) then
                   mesh%v(i)%tr(j)=mesh%v(minind)%tr(l)
                   !print*, "      Old TR", mesh%v(minind)%tr(l)
                   exit
                end if
                if(mesh%v(minind)%nb(l)==i.and.minind==ind(2)) then
                   m=modint(l-1, mesh%v(minind)%nnb)
                   mesh%v(i)%tr(j)=mesh%v(minind)%tr(m)
                   !print*, "      Old TR", mesh%v(minind)%tr(m)
                   exit
                end if
             end do
          end if
       end do
    end do

    !Create a list of properties of triangles
    do k=1, mesh%nt
       !Define properties of this triangle
       do l=1,3
          m=modint( l+1, 3)
          !Set triangles' edge indexes
          mesh%tr(k)%ed(l)=getedge(mesh%tr(k)%v(l), mesh%tr(k)%v(m), mesh)
          !Store the triangles sharing an edge
          if(mesh%ed(mesh%tr(k)%ed(l))%sh(1)==0)then
             mesh%ed(mesh%tr(k)%ed(l))%sh(1)=k
          else
             mesh%ed(mesh%tr(k)%ed(l))%sh(2)=k
          end if
          !Normal and tangent components to triangle's edge
          nr=-cross_product(mesh%v(mesh%tr(k)%v(l))%p, &
               mesh%v(mesh%tr(k)%v(m))%p)
          nr=nr/norm(nr)
          dp=dot_product(nr,mesh%ed((mesh%tr(k)%ed(l)))%nr)
          if(abs(abs(dp)-1_r8)>eps*100)then
             !If the vectors are not colinear
             print*, "MESHSTRUCT ERROR: triangle's normal and edge's", &
                  " normal vectors are not colinear"
             print*, "Dot product:", dp
             stop
          end if
          !Define sign correction
          if(dp>0)then
             mesh%tr(k)%nr(l)=+1
             mesh%tr(k)%tg(l)=-1
          else
             mesh%tr(k)%nr(l)=-1
             mesh%tr(k)%tg(l)=+1
          end if
       end do
    end do


    !Debug
    !do i=1, mesh%nv
    !  print*, i
    !  print*, mesh%v(i)%tr(1:mesh%v(i)%nnb)
    !  print*, mesh%v(i)%ed(1:mesh%v(i)%nnb)
    !end do

    if(.not.opt)then
       !Create list of triangles neighbour triangles
       do k=1, mesh%nt   !For each triangle
          do i=1,3     !For each edge
             l=mesh%tr(k)%ed(i)
             !Check triangles that share the edge
             !  and get the one diferent from 'k'
             if(mesh%ed(l)%sh(1)==k)then
                mesh%tr(k)%nb(i)=mesh%ed(l)%sh(2)
             else
                mesh%tr(k)%nb(i)=mesh%ed(l)%sh(1)
             end if
          end do
       end do
    end if

    !----------------------------------------------------------
    ! Calculate triangle's circumcenters, barycenter, circumradius and areas
    !----------------------------------------------------------
    mesh%mintrarea=(1.e20)
    mesh%maxtrarea=0.
    mesh%meantrarea=0.
    mesh%mintrangle=(1.e20)
    mesh%maxtrangle=0.
    mesh%meantrangle=0.

    do k=1,mesh%nt
       i1=mesh%tr(k)%v(1)
       i2=mesh%tr(k)%v(2)
       i3=mesh%tr(k)%v(3)
       !Calculate circumcenter
       c=trcircumc ( mesh%v(i1)%p,  mesh%v(i2)%p,  mesh%v(i3)%p)
       mesh%tr(k)%c%p=c

       !Test if it is inside the triangle - warn if not
       isintr=insidetr( c, mesh%v(i1)%p,  mesh%v(i2)%p,  mesh%v(i3)%p)
       if((.not.isintr))then
          iwarn3=iwarn3+1
          if(iwarn3==1.and. (.not. trim(mesh%optm) == "scvt"))then
             print*, "MESHSTRUCTURE Warning: Circumcenter outside triangle"
             print*, "     tr:", k, " ier of inside:", ier
          end if
          !Stop if using SPRG or SALT
          !if( trim(mesh%optm)=='sprg' .or. trim(mesh%optm)=='salt') then
          !opt=.false.  !To end optimization
          !return
          !end if
       end if

       if(.not.opt)then
          !Convert to spherical coordinates and store of circumcenter
          call cart2sph( c(1), c(2), c(3), mesh%tr(k)%c%lon, mesh%tr(k)%c%lat)
       end if

       !Calculate barycenter
       c=trbarycenter(k, mesh)
       mesh%tr(k)%b%p=c

       if(.not.opt)then
          !Convert to spherical coordinates and store
          call cart2sph( c(1), c(2), c(3), mesh%tr(k)%b%lon, mesh%tr(k)%b%lat)
       end if

       !Calculate geodesical/spherical circumradius in radians of barycenter
       mesh%tr(k)%radg=arclen(mesh%tr(k)%c%p, mesh%v(mesh%tr(k)%v(1))%p)

       !Calculate planar circumradius in unit sphere metrics
       mesh%tr(k)%radp=norm(mesh%tr(k)%c%p-mesh%v(mesh%tr(k)%v(1))%p)

       !Calculate geodesical/spherical areas for the unit sphere
       p1=mesh%v(mesh%tr(k)%v(1))%p
       p2=mesh%v(mesh%tr(k)%v(2))%p
       p3=mesh%v(mesh%tr(k)%v(3))%p
       area=sphtriarea(p1, p2, p3)
       mesh%tr(k)%areag=area
       mesh%mintrarea=min(mesh%mintrarea, area)
       mesh%maxtrarea=max(mesh%maxtrarea, area)
       mesh%meantrarea=mesh%meantrarea+area

       if(.not.opt)then
          !Calculate planar area for points in the unit sphere
          area = norm(cross_product( &
               mesh%v(mesh%tr(k)%v(2))%p-mesh%v(mesh%tr(k)%v(1))%p, &
               mesh%v(mesh%tr(k)%v(3))%p-mesh%v(mesh%tr(k)%v(1))%p))/2
          mesh%tr(k)%areap=area
       end if

       !Calculate the internal angles of the triangle
       angles=sphtriangles(p1, p2, p3)
       mesh%tr(k)%angles=angles
       !print '(i4, 3f10.2)', k, angles*rad2deg
       mesh%mintrangle=min(mesh%mintrangle, minval(angles(1:3)))
       mesh%maxtrangle=max(mesh%maxtrangle, maxval(angles(1:3)))
       mesh%meantrangle=mesh%meantrangle+sum(angles(1:3))
    end do
    mesh%meantrarea=mesh%meantrarea/mesh%nt
    mesh%meantrangle=mesh%meantrangle/mesh%nt/3

    !Print warnings if needed
    if(iwarn3>0)then
       print*, "Meshstruct warning: Num circumc", &
            " out of triangles=", iwarn3
       stat=1
    end if

    if(.not.opt .and. showonscreen )then
       !Print area info
       print*, "Areas of spherical triangles (unit sphere)"
       print*, "    Minimum    Maximum    Average     Min/Max"
       print '(3f12.6, f12.8)', mesh%mintrarea, mesh%maxtrarea, mesh%meantrarea, mesh%mintrarea/mesh%maxtrarea
       print*
    end if

    if(.not.opt .and. showonscreen)then
       print*, "Internal angles of spherical triangles (degrees)"
       print*, "    Minimum    Maximum    Average     Min/Max"
       print '(3f12.6, f12.8)', mesh%mintrangle*rad2deg, mesh%maxtrangle*rad2deg, &
            mesh%meantrangle*rad2deg, mesh%mintrangle/mesh%maxtrangle
       print*
    end if

    !-------------------------------------------------------
    ! Create edges for hexagons/pentagons (Voronoi cells)
    !-------------------------------------------------------
    if(.not.opt  .and. showonscreen)then
       print*, "Defining Voronoi structure ... "
       print*
    end if

    !Allocate voronoi edges array
    if(.not.allocated(mesh%edhx))then
       allocate(mesh%edhx(1:mesh%ne))
    else
       !If the new array has a different size, reallocate
       !print*, ubound(mesh%ed, 1), mesh%ne
       if(ubound(mesh%edhx, 1)/=mesh%ne)then
          deallocate(mesh%edhx)
          allocate(mesh%edhx(1:mesh%ne), stat=status)
       end if
    end if

    mesh%mincdist=1000000._r8
    mesh%maxcdist=0._r8
    mesh%meancdist=0._r8
    do i=1, mesh%ne
       !Define vertices as triangle index (circumcenters)
       mesh%edhx(i)%v=mesh%ed(i)%sh

       !Define neighbours as hexagons (vertice indexes)
       mesh%edhx(i)%sh=mesh%ed(i)%v

       !The hx edge meanpoint in not necessarily
       !  the crossing point between ed(i) and edhx(i)
       !  the crossing point is exactly the midpoint of ed(i)
       i1=mesh%edhx(i)%v(1)
       i2=mesh%edhx(i)%v(2)
       mesh%edhx(i)%c=midpoint_edhx(i1,i2,mesh)

       !Vector normal to edge
       mesh%edhx(i)%nr=cross_product(mesh%tr(i1)%c%p, mesh%tr(i2)%c%p)
       mesh%edhx(i)%nr=mesh%edhx(i)%nr/norm(mesh%edhx(i)%nr)
       !Tangent vector - Points in the direction from tr(i1) --> tr(i2)
       mesh%edhx(i)%tg=cross_product(mesh%edhx(i)%nr, mesh%edhx(i)%c%p)
       mesh%edhx(i)%tg=mesh%edhx(i)%tg/norm(mesh%edhx(i)%tg)

       !Edge euclidian length (straight line)
       mesh%edhx(i)%lenp=norm(mesh%tr(i1)%c%p-mesh%tr(i2)%c%p)

       !Edge geodesical length (curved line)
       mesh%edhx(i)%leng=arclen(mesh%tr(i1)%c%p, mesh%tr(i2)%c%p)

       !Mesh Distance measures
       mesh%mincdist=min(mesh%mincdist, mesh%edhx(i)%leng)
       mesh%maxcdist=max(mesh%maxcdist, mesh%edhx(i)%leng)
       mesh%meancdist=mesh%meancdist + mesh%edhx(i)%leng

    end do
    mesh%meancdist=mesh%meancdist/mesh%ne

    !-------------------------------------------------------
    ! Create hexagonal cell structure
    !------------------------------------------------------
    mesh%minhxarea=1000000._r8
    mesh%maxhxarea=0._r8
    mesh%meanhxarea=0._r8

    !Allocate Voronoi cell array
    if(.not.allocated(mesh%hx))then
       allocate(mesh%hx(1:mesh%nv))
    else
       !If the new array has a different size, reallocate
       !print*, ubound(mesh%ed, 1), mesh%ne
       if(ubound(mesh%hx, 1)/=mesh%nv)then
          deallocate(mesh%hx)
          allocate(mesh%hx(1:mesh%nv), stat=status)
       end if
    end if

    !For each voronoi cell
    do i=1, mesh%nv
       if(.not.opt)then
          !Allocate space for this cell
          nb=mesh%v(i)%nnb
          call reallocint(mesh%hx(i)%tg, nb)
          call reallocint(mesh%hx(i)%nr, nb)
          call reallocint(mesh%hx(i)%ttgout, nb)
          call reallocint(mesh%hx(i)%tnrccw, nb)

          !For the coresponding neighbour edges/triangles
          do j=1,nb
             !Get triangle center points for this edge
             if(j==1) then
                j1=mesh%v(i)%tr(nb)
                j2=mesh%v(i)%tr(1)
                k=mesh%v(i)%ed(1)
             else
                j1=mesh%v(i)%tr(j-1)
                j2=mesh%v(i)%tr(j)
                k=mesh%v(i)%ed(j-1)
             end if
             !Normal and tangent components
             nr=-cross_product(mesh%tr(j1)%c%p, mesh%tr(j2)%c%p)
             nr=nr/norm(nr)
             dp=dot_product(nr,mesh%edhx(mesh%v(i)%ed(j))%nr)
             if(abs(abs(dp)-1_r8)>eps*100)then
                !If the vectors are not colinear
                print*, "MESHSTRUCT ERROR: Normal to hexagon edge and edge vectors not colinear"
                print*, "Cell: ", i
                print*, "Tr:",  j1, j2
                print*, "Ed:",  k
                print*, "Dot product:", dp
                stop
             end if

             !Define sign correction
             if(dp>0)then
                mesh%hx(i)%nr(j)=+1
                mesh%hx(i)%tg(j)=-1
             else
                mesh%hx(i)%nr(j)=-1
                mesh%hx(i)%tg(j)=+1
             end if
             !Calculate correction relative to trinagle normal and tangents
             !=dot_product(mesh%edhx()%nr,mesh%ed()%tg)*mesh%hx()%nr
             !Calculate outward normal to edge product to tr tangent
             dp=dot_product(mesh%ed(mesh%v(i)%ed(j))%tg, &
                  mesh%edhx(mesh%v(i)%ed(j))%nr*mesh%hx(i)%nr(j))
             mesh%hx(i)%ttgout(j)=dsign( 1._r8, dp)
             !Calculate CCW tg of edge product to tr normal
             dp=dot_product(mesh%ed(mesh%v(i)%ed(j))%nr, &
                  mesh%edhx(mesh%v(i)%ed(j))%tg*mesh%hx(i)%tg(j))
             mesh%hx(i)%tnrccw(j)=dsign( 1._r8, dp)
          end do
          !print*, i
          !print*, "nb: ", mesh%v(i)%nb(1:nb)
          !print*, "ttgout: ", mesh%hx(i)%ttgout(1:nb)
          !print*, "tnrccw: ", mesh%hx(i)%tnrccw(1:nb)

          !Calculate aligment index
          mesh%hx(i)%align=alignind(i, mesh)
       end if

       !Calculate Barycenter (mass center)
       c=vorbarycenter(i, mesh)
       mesh%hx(i)%b%p=c
       if(.not.opt)then
          !Convert to spherical coordinates and store
          call cart2sph( c(1), c(2), c(3), mesh%hx(i)%b%lon, mesh%hx(i)%b%lat)
       end if

       !Calculate Areas
       mesh%hx(i)%areap=planarhexagarea(i, mesh)
       mesh%hx(i)%areag=sphhxarea(i, mesh)

       !Mesh Area measures
       mesh%minhxarea=min(mesh%minhxarea, mesh%hx(i)%areag)
       mesh%maxhxarea=max(mesh%maxhxarea, mesh%hx(i)%areag)
       mesh%meanhxarea=mesh%meanhxarea + mesh%hx(i)%areag

    end do

    mesh%meanhxarea=mesh%meanhxarea/mesh%nv

    if(.not.opt  .and. showonscreen)then
       print*, "Areas of spherical voronoi cells (unit sphere)"
       print*, "    Minimum    Maximum    Average     Min/Max"
       print '(3f12.6, f12.8)', mesh%minhxarea, mesh%maxhxarea, mesh%meanhxarea, mesh%minhxarea/mesh%maxhxarea
       print*
    end if


    !--------------------------------------------------------
    ! Geodesic to regular grid convertion tool - search table
    !-------------------------------------------------------
    if(.not.opt)then
       print*, "Creating geodesic to regular grid tool (Search table) ... "

       !Calculate the regular grid spacing
       if(trim(mesh%kind)=='rand' .and. trim(mesh%optm)=='nopt')then
          print *
          print *, "Mesh with random nodes: Search table may not be trustworthy"
          print *
          mesh%dlat=mesh%meanvdist/10._r8
       else
          ! Forces to be less then 3 times smaller than the smallest triangle
          if(trim(mesh%pos)=='ref')then
             mesh%dlat=mesh%minvdist/3._r8
          else
             mesh%dlat=mesh%minvdist/3.5_r8
          end if
       end if
       mesh%nlat=int (pi/mesh%dlat, i4)
       mesh%nlat=int(real(mesh%nlat, r8)/real(10.0,r8))*10_i4+10_i4
       mesh%dlat=pi/real(mesh%nlat, r8)
       print*, "nlat=",mesh%nlat," dlat=",mesh%dlat*rad2deg, " degrees"

       call geo2reg(mesh)

       print*,"Mesh structured."
       print*
    end if

    return
  end subroutine meshstruct

  subroutine icosanodes0( x, y, z, pos)
    !-----------------------------------------------------------
    ! ICOSANODES0
    !   Calculates the basic icasoedral point
    !     It positions the icosahedron symmetrical within hemespheres (eqs),
    !   or with nodes on poles (pol), depending on the variable mesh%pos
    !------------------------------------------------------------
    !Number of nodes
    integer (i4), parameter :: n=12

    !Array with nodes
    real(r8):: x(1:n)
    real(r8):: y(1:n)
    real(r8):: z(1:n)

    !Auxiliar variables
    real (r8):: radius
    real (r8):: sqrt5
    real (r8):: phi
    real (r8):: ratio
    real (r8):: a
    real (r8):: b
    real (r8):: lon
    real (r8):: lat
    real (r8):: dlon
    integer:: i
    character (len=3):: pos
    character (len=3):: pos0

    if(trim(pos)=="ref")then
       pos0="pol"
    else
       pos0=pos
    end if

    ! Creating basic icosahedral shape
    select case(pos0)
    case("eqs") !Equator symmetric
       radius = 1.0_r8
       sqrt5 = dsqrt(5.0_r8)
       phi = (1.0_r8+sqrt5)*0.5_r8
       ratio = dsqrt(10.0_r8+(2.0_r8*sqrt5))/(4.0*phi)

       a = (radius/ratio)*0.5_r8
       b = (radius/ratio)/(2.0_r8*phi)

       !The order of the nodes is very important, do not change as
       ! it may cause stripack to enter an infinite loop
       x(1) = 0; y(1) = b; z(1) =-a
       x(2) = 0; y(2) =-b; z(2) =-a
       x(3) = a; y(3) = 0; z(3) =-b
       x(4) =-a; y(4) = 0; z(4) =-b
       x(5) = b; y(5) = a; z(5) = 0
       x(6) =-b; y(6) = a; z(6) = 0
       x(7) =-a; y(7) = 0; z(7) = b
       x(8) = a; y(8) = 0; z(8) = b
       x(9) = b; y(9) =-a; z(9) = 0
       x(10) =-b; y(10) =-a; z(10) = 0
       x(11)= 0; y(11)= b; z(11)= a
       x(12)= 0; y(12)=-b; z(12)= a

    case("pol") !Point on poles

       !North pole
       x(1) = 0; y(1) = 0; z(1) =1
       !Latitude for North Hemisphere nodes
       lat= datan (real(0.5,8))
       !Longitude variation between nodes
       dlon=2/5.*pi
       lon=-pi
       do i=1,5
          call sph2cart(lon, lat, x(1+i), y(1+i), z(1+i))
          lon=lon+dlon
       end do
       !South pole
       x(7) = 0; y(7) = 0; z(7) =-1
       !Latitude for South Hemisphere nodes
       lat= -lat
       lon= -4*pi/5.
       do i=1,5
          call sph2cart(lon, lat, x(7+i), y(7+i), z(7+i))
          lon=lon+dlon
       end do
    end select

    !Writes in mesh variable the nodes coordinates (old)
    !call getunit(ixyzunit)
    !open(ixyzunit, file=axyzgrid,status='replace')
    !write(ixyzunit,'(i8)') n
    !do i=1,n
    !   mesh%v(i)%p(1)=x(i)
    !   mesh%v(i)%p(2)=y(i)
    !   mesh%v(i)%p(3)=z(i)
    !write(ixyzunit,'(i8,3f32.16)') i,x(i),y(i),z(i)
    !end do
    !close(ixyzunit)

    return
  end subroutine icosanodes0

  subroutine icosanodes(mesh)
    !-----------------------------------------------------------
    !	ICOSANODES
    !	  Calculates the icosahedral nodes on the sphere and writes
    !	it on a mesh structure variable.
    !     It positions the icosahedron symmetrical within hemespheres (eqs),
    !   or with nodes on poles (pol), depending on the variable mesh%pos
    !------------------------------------------------------------

    !Grid structure variable
    type(grid_structure), intent(inout) :: mesh   !Input and output

    !Filename
    !character (len=60), intent(in):: axyzgrid

    !Number of nodes
    integer (i4):: n

    !Auxiliar variables
    integer:: n1
    integer:: i
    integer:: j
    integer:: ier
    character (len=3):: pos
    real(r8), allocatable :: x(:)
    real(r8), allocatable :: y(:)
    real(r8), allocatable :: z(:)

    pos=trim(mesh%pos)
    !Check for the position of the mesh
    select case(mesh%pos)
    case("eqs") !Equator symmetric
       print*, "Generating mesh using hemesphere symmetrical icosahedral ..."
    case("pol") !Point on poles
       print*, "Generating mesh using icosahedral nodes on poles..."
    case("ref") !Refined grid - use points on poles
       pos="pol"
       if(trim(mesh%optm(1:4))=="scvt")then
          print *, "Generating mesh using icosahedral nodes on poles, but with local refinement"
       else
          print*, "  Generating mesh using icosahedral nodes on poles..."
          print *, "  (Ignored locally refined positioning, because SCVT was not set)"
          mesh%pos="pol"
       end if
    case default
       print*, "ERROR ICOSANODES: Cannot generate initial icosahedral mesh"
       print*, "Icosahedral mesh position undefined. ", mesh%pos
       stop
    end select

    !Alocate space
    n=mesh%nv
    ier=0
    allocate(x(1:n),y(1:n),z(1:n), stat=ier)
    allocate(mesh%v(1:n),stat=ier)
    if(ier>0) then
       print*,"ICOSANODES ERROR: allocation problem. Size of arrays:", n
       stop
    end if

    call icosanodes0(x(1:12), y(1:12), z(1:12), mesh%pos)
    n1=12
    !Calculates refinements puting points on edge's midpoints
    do j=1,mesh%glevel
       call meshrefin(n1,x,y,z,n)
    end do

    !Save nodes
    do i=1,n
       mesh%v(i)%p(1)=x(i)
       mesh%v(i)%p(2)=y(i)
       mesh%v(i)%p(3)=z(i)
       !write(ixyzunit,'(i8,3f32.16)') i,x(i),y(i),z(i)
    end do
    !close(ixyzunit)

    return
  end subroutine icosanodes


  subroutine octgnodes(mesh)
    !-----------------------------------------------------------
    !	OCTGNODES
    !	  Calculates the octahedron nodes on the sphere and writes
    !	it on a mesh structure variable.
    !     It positions the octahedron symmetrical within hemespheres (eqs) even
    !   if flag is set to 'pol'
    !------------------------------------------------------------

    !Grid structure variable
    type(grid_structure), intent(inout) :: mesh   !Input and output

    !Filename
    !character (len=60), intent(in):: axyzgrid

    !Number of nodes
    integer (i4):: n

    !Auxiliar variables
    real (r8):: lon
    real (r8):: lat
    real (r8):: dlon
    integer:: n1
    integer:: i
    integer:: j
    integer:: ier
    character (len=3):: pos
    real(r8), allocatable :: x(:)
    real(r8), allocatable :: y(:)
    real(r8), allocatable :: z(:)

    pos=trim(mesh%pos)
    !Check mesh position - Octahedral mesh do not depend on this flag
    if(trim(mesh%pos)=='eqs' .or. trim(mesh%pos)=='pol')then
       print*, "Generating mesh using hemisphere symmetrical octahedron,"
       print*, "    with point on poles..."
    elseif(trim(mesh%pos)=="ref")then !Refined grid - use points on poles
       pos="pol"
       if(trim(mesh%optm)=="scvt")then
          print*, "Generating mesh using hemisphere symmetrical octahedron,"
          print*, "    with point on poles, and with local refinement"
       else
          print*, "Generating mesh using hemisphere symmetrical octahedron,"
          print*, "    with point on poles..."
          print *, "  (Ignored locally refined positioning, because SCVT was not set)"
          mesh%pos="pol"
       end if
    else
       print*, "ERROR OCTGNODES: Cannot generate initial octahedron mesh"
       print*, "Octahedron mesh position undefined."
       stop
    end if

    !Allocate space
    n=mesh%nv
    ier=0
    allocate(x(1:n),y(1:n),z(1:n), stat=ier)
    allocate(mesh%v(1:n),stat=ier)
    if(ier>0) then
       print*,"OCTGNODES ERROR: allocation problem. Size of arrays:", n
       stop
    end if

    ! Creating basic octogonal shape
    n1=6
    select case(pos)
    case("eqs", "pol")
       !North pole
       x(1) = 0; y(1) = 0; z(1) =1
       !Equator nodes
       lat= 0.0_r8
       !Longitude variation between nodes
       dlon=0.5_r8*pi
       lon=-pi
       do i=1,4
          call sph2cart(lon, lat, x(1+i), y(1+i), z(1+i))
          lon=lon+dlon
       end do
       !South pole
       x(6) = 0; y(6) = 0; z(6) =-1
    end select

    !Calculates refinements
    do j=1,mesh%glevel
       call meshrefin(n1,x,y,z,n)
    end do

    !Writes in mesh variable the nodes coordinates
    !call getunit(ixyzunit)
    !open(ixyzunit, file=axyzgrid,status='replace')
    !write(ixyzunit,'(i8)') n
    do i=1,n
       mesh%v(i)%p(1)=x(i)
       mesh%v(i)%p(2)=y(i)
       mesh%v(i)%p(3)=z(i)
       !write(ixyzunit,'(i8,3f32.16)') i,x(i),y(i),z(i)
    end do
    !close(ixyzunit)

    return
  end subroutine octgnodes

  subroutine nodesizeadj(n, nbisect, kind)
    !---------------------------------------------------------------
    ! nodesizeadj
    !  Calculates the number of nodes for an icosahedral/octahedral mesh
    !  Adjusts the user input to possible icosahedron size
    !
    !  Input   :
    !     n - proposed number of nodes
    !     kind - icos or octg grid kind
    !  Returns :
    !     n - adjusted number of nodes (if n is suitable, returns n)
    !     nbisect - Number of refinements on icosahedron
    !---------------------------------------------------------------

    integer (i4):: n
    integer (i4):: nbisect
    character(len=16):: kind

    !Parameter for maximum number of nodes (Aprox. 42 million nodes)
    integer (i4), parameter :: maxsize=12

    !Auxiliar variables
    integer (i4):: sizes(maxsize)
    integer:: i

    if(trim(kind)=='rand' .or. trim(kind)=='read')then
       return
    end if

    !Ajust n to valid number of nodes for bisected octahedron
    print*," Adjusting number of mesh points ..."
    print*," Proposed number of vertices: ", n
    do i=1,maxsize
       if(trim(kind)=='icos')then
          sizes(i)=5*2**(2*i-1)+2
       elseif(trim(kind)=='octg')then
          if(i==1)then !octahedron
             sizes(1)=6
          else
             sizes(i)=4*sizes(i-1)-6
          end if
       else
          print*, "NODESIZEADJ ERROR: unknown kind of grid ", trim(kind)
       end if
       nbisect=i-1
       !print*,i, sizes(i)
       if(n<=sizes(i))then
          n=sizes(i)
          exit
       end if
    end do
    print*, " Adjusted number of vertices:",n
    print*

    if(n>sizes(nbisect+1))then
       print*, "NODESIZEADJ ERROR: number of nodes too large!"
       print*, " If you want more than 42 million nodes please change code parameters."
       stop
    end if

    return
  end subroutine nodesizeadj

  subroutine glevelnodeindex(glevel, ni, nf, kind)
    !---------------------------------------------------------------
    ! glevelnodeindex
    !
    !  Return the initial (ni) and final (nf) node index for the glevel wanted
    !  Kind indicates the type of grid
    !---------------------------------------------------------------

    integer (i4), intent(in) :: glevel
    integer (i4), intent(out) ::  ni
    integer (i4), intent(out) :: nf
    character(len=16):: kind
    integer:: i

    !Calculate indexes for this level
    i=glevel+1
    if(trim(kind)=='icos')then
       ni=5*2**(2*(i-1)-1)+2+1

       nf=5*2**(2*i-1)+2
    else
       print*, "glevelnodeindex warning: not implmented for this kind of grid : ", trim(kind)
    end if

    return
  end subroutine glevelnodeindex

  subroutine meshrefin(n1,x,y,z,n)
    !------------------------------------------------------------
    !   MESHREFIN
    !	Creates 1 level of refinement on mesh
    !   Includes nodes on the medium points lying on each edge
    !
    !   n1 -> original number of point on input. The new number of
    !         points are set in n1 on output
    !   x, y, z -> arrays of cartesian coords of n1 nodeson input
    !   n  -> Array sizes n>=n1 and must be pre-calculated correcly
    !
    !   The output is given in x, y, z, in the sequence
    !         imediately after the original given points
    !------------------------------------------------------------
    integer(i4), intent(in) :: n  !Total size of arrays
    integer(i4), intent(inout) :: n1 !Actual used size of arrays

    real(r8), intent(inout) ::  x(1:n)
    real(r8), intent(inout) :: y(1:n)
    real(r8), intent(inout) :: z(1:n)

    !Triangulation arrays
    integer(i4), allocatable :: list(:)
    integer(i4), allocatable :: lptr(:)
    integer(i4), allocatable :: lend(:)
    integer(i4), allocatable :: iwk(:)
    real(r8), allocatable :: dist(:)
    integer (i4):: lnew
    integer(i4), allocatable :: neigh(:,:)
    integer(i4), allocatable :: neisz(:)

    !Indexes
    integer(i4):: node
    integer(i4):: lpl
    integer(i4):: k
    integer(i4):: lp
    integer(i4):: ier
    integer(i4):: i
    integer(i4):: j

    !Auxiliar variable
    real (r8):: r

    allocate(list(1:3*(2*n1-4)), lptr(1:3*(2*n1-4)), lend(1:n1))
    allocate(iwk(1:3*(2*n1-4)))
    allocate(dist(1:n1))
    allocate(neigh(1:n1,1:30), neisz(1:n1))

    !Create a triangulation with given points
    call trmesh (n1,x,y,z,list,lptr,lend,lnew,iwk,iwk(n1+1),dist,ier)

    !Create a list of neighbours for each point
    !dir$ loop count min(8)
    do  node = 1,n1
       lpl = lend(node)
       lp = lpl
       do k=1,20
          lp = lptr(lp)
          neigh(node,k) = list(lp)
          if(lp == lpl)then
             exit
          end if
       end do
       neisz(node) = k
       !print*,"Node : ",node," neighb:", neigh(node,1:neisz(node))
    enddo

    !  Add points into the middle of edges
    k=n1+1
    do i = 1,n1
       do j = 1,neisz(i)
          if (neigh(i,j)>i) then
             x(k) = 0.5*(x(i)+x(neigh(i,j)))
             y(k) = 0.5*(y(i)+y(neigh(i,j)))
             z(k) = 0.5*(z(i)+z(neigh(i,j)))
             r = dsqrt(x(k)**2+y(k)**2+z(k)**2)
             x(k) = x(k)/r
             y(k) = y(k)/r
             z(k) = z(k)/r
             k=k+1
          endif
       end do
    enddo
    !New number of points (original + added ones)
    n1=k-1
    !print *,"Number of Final Generators = ",n1

    return
  end subroutine meshrefin

  subroutine randnodes(mesh)
    !-----------------------------------------------------------
    ! RANDNODES
    !   Calculates a ramdom number of points on the sphere
    !     If mesh%pos=pol, sets points on north and south pole
    !        mesh%pos=eqs, set equator symmetry
    !        mesh%pos=ran, set all nodes randomly
    !        mesh%pos=ref, random nodes with local refinement
    !------------------------------------------------------------

    !Grid structure variable
    type(grid_structure), intent(inout) :: mesh   !Input and output

    !Number of nodes
    integer (i4):: n
    integer (i4):: m

    !Auxiliar variables
    integer:: i

    !Number of nodes
    n=mesh%nv
    if(mod(n,2) == 1 .and. trim(mesh%pos) == "eqs") then
       print*, "RANDNODES warning: odd number of nodes and symmetric mesh wanted"
       print*, "    Adjusting number of nodes + 1 "
       n=n+1
       mesh%nv=n
    end if
    allocate(mesh%v(1:n))

    !  Set the seed for random number generator
    m=n
    call setrandomseed(m)

    ! Creating basic icosahedral shape
    select case(trim(mesh%pos))
    case("eqs") !Equator symmetric
       print*, "Generating mesh using hemesphere symmetrical random nodes ..."

       !North pole
       mesh%v(1)%p(1)= 0._r8
       mesh%v(1)%p(2)= 0._r8
       mesh%v(1)%p(3)= 1._r8

       !Set each pair of random nodes
       do i=2, n/2
          mesh%v(i)%p=sphrandom_denspt(mesh%pos)
          mesh%v(i+n/2-1)%p=mesh%v(i)%p
          mesh%v(i+n/2-1)%p(3)=-mesh%v(i+n/2-1)%p(3)
          ! print*, i, mesh%v(i)%p
          ! print*, i+1, mesh%v(i+1)%p
       end do
       !South pole
       mesh%v(n)%p=-mesh%v(1)%p

    case("pol") !Nodes on poles
       print*, "Generating mesh using random nodes but with nodes on poles..."

       !North pole
       mesh%v(1)%p(1)= 0._r8
       mesh%v(1)%p(2)= 0._r8
       mesh%v(1)%p(3)= 1._r8

       !Set each random node
       do i=2,n-1
          mesh%v(i)%p=sphrandom_denspt(mesh%pos)
       end do
       !South pole
       mesh%v(n)%p=-mesh%v(1)%p

    case default !Nodes random
       print*, "Generating mesh using random nodes ..."
       do i=1,n
          mesh%v(i)%p=sphrandom_denspt(mesh%pos)
       end do

    end select

    return
  end subroutine randnodes


  !================================================================================
  !    MESH OPTIMIZATION ROUTINES
  !================================================================================

  subroutine optimizegrid(mesh)
    !--------------------------------------------------------------
    !  OPTIMIZEGRID
    !   Optimizes grid based on mesh%optm
    !   Implemented:
    !   sprg - Spring Dynamics - tested ok
    !   scvt - Centroidal Voronoi Tesselation - tested ok
    !   hr95 - HR1995 Miura's algorithm (problematic)
    !   salt - Spherical Aligned Tesselation (problematic)
    !--------------------------------------------------------------

    !Original mesh
    type(grid_structure), intent(inout) :: mesh

    !Primary icosahedral grid
    type(grid_structure):: icos0

    !List of cartesian coordinates of nodes (temporary)
    type(vector_field_cart):: vec !, nodes
    type(vector_field_cart):: nodesmin !, nodes

    !List of nodes with info on how nodes will be optimized
    integer (i4), allocatable:: nodeslist(:)

    !Energy varibles
    real (r8):: energy
    real (r8):: energytmp

    !Aux vars
    real (r8):: dm
    !real (r8):: dt
    real (r8):: epsdif
    real (r8):: epsmin
    real (r8):: inienergy

    !Temp point variable
    real (r8)::  p(1:3)

    !Temp vector var
    real (r8)::  np(1:3)

    !Flag for otimization
    logical:: opt

    !Indexes
    integer (i4):: i
    integer (i4):: j
    integer (i4):: k
    integer (i4):: ni
    integer (i4):: nf
    integer (i4):: maxiter
    integer (i4):: miniter
    integer (i4):: counter
    integer (i4):: upstop
    integer (i4):: stat

    !File variables
    integer:: iunit
    integer:: iunit2
    character (len=100):: filename
    !character (len=100):: method
    logical:: ifile
    logical:: iopen

    print*, "Optimizing grid ..."
    print*

    if(mesh%nv<14)then
       print*, "Mesh with ", mesh%nv," points: "
       print*, "   No optimization needed"
       return
    end if

    !Error counter
    counter=0

    !Allocate space for temporary nodes
    allocate(nodesmin%p(1:mesh%nv))
    nodesmin%pos=0
    nodesmin%n=mesh%nv
    nodesmin%name="nodesmin_"//mesh%name

    !Save a copy of the grid points, as local minimum
    do j=1, mesh%nv
       nodesmin%p(j)%v=mesh%v(j)%p
    end do

    !Structure only important mesh characteristics
    if(mesh%hrchy>0)then
       print*, "Restructuring grid: ignore meshstruct warnings"
    end if
    opt=.true.
    call meshstruct(mesh, opt, stat)
    stat=0
    print*

    !Allocate space for tmp array of forcing terms
    allocate(vec%p(1:mesh%nv))
    vec%pos=0
    vec%n=mesh%nv
    vec%name="vec_"//mesh%name
    miniter=1
    select case (trim(mesh%optm))
    case ("scvt")  !Spherical Centroidal Voronoi Tesselation
       !No optimization parameter needed
       if(trim(mesh%pos)=="ref")then
          !Stop criteria - difference between energies
          if(mesh%nv<45000)then
             maxiter=100000 !6000 !2000
             epsdif=mesh%meanvdist*eps/10000.0_r8
             epsmin=mesh%meanvdist*eps*1.0_r8 !/100.0_r8
             miniter=40000
          elseif(mesh%nv>450000)then
             maxiter=40000 !6000 !2000
             epsdif=mesh%meanvdist*eps/1000.0_r8
             epsmin=mesh%meanvdist*eps*10.0_r8 !/100.0_r8
             miniter=2000
          else
             maxiter=40000 !6000 !2000
             epsdif=mesh%meanvdist*eps/1000.0_r8
             epsmin=mesh%meanvdist*eps*10.0_r8 !/100.0_r8
             miniter=5000
          end if
       else
          !Stop criteria - difference between energies
          epsdif=mesh%meanvdist*eps/10000.0_r8
          !Stop criteria - minimum energy
          epsmin=mesh%meanvdist*eps/100.0_r8
          maxiter=4000 !6000 !2000
       end if
       !print*, epsmin, mesh%meanvdist
       !Maximum number of iterations
       !maxiter=10000 !6000 !2000
       !Do not stop it energy grows
       upstop=0

    case ("sprg") !Spring Dynamics optimization

       !Relaxation parameter - Miura's scheme
       !  or Tomita's dt for ODE solver
       !mesh%optpar=0.016_r8  !Tomita
       !mesh%optpar=0.04     !OLAM
       !mesh%optpar=0.3     !Miura

       !If energy used
       !mesh%optpar=0.3_r8   !Miura
       !print*, "Using Miura's algorithm "
       !Stop criteria - difference between energies
       !epsdif=mesh%meanvdist*eps*0.0001 !100_r8
       !Stop criteria - minimum energy
       !epsmin=mesh%meanvdist*0.000001_r8

       !If disrtotion used - prefered due to finer grids
       mesh%optpar=0.1_r8   !Miura
       !Stop criteria - difference between energies
       epsdif=mesh%meanvdist*eps/10 !100_r8
       !Stop criteria - minimum energy
       epsmin=mesh%meanvdist*0.01_r8

       print*, "Optimization parameter:", mesh%optpar

       !Maximum number of iterations
       maxiter=4000 !6000 !2000

       !Stop it energy grows (local minimum)
       upstop=1

    case ("salt") !Spherical Aligned Tesselation
       !mesh%optpar=0.005_r8
       !For 10242 or finer grids
       !mesh%optpar=0.02_r8
       !For 2000
       mesh%optpar=0.0005_r8
       !For 642
       !mesh%optpar=0.00005_r8
       !For 160
       !mesh%optpar=0.0002_r8
       print*, "Optimization parameter:", mesh%optpar
       epsdif=0.00000001_r8 !mesh%meanvdist*eps*100_r8
       epsmin=0.0001_r8 !mesh%meanvdist*0.001_r8
       !Maximum number of iterations
       maxiter=10000 !6000 !2000
       !Do not stop it energy grows
       upstop=0

    case("hr95") !Edge alignement optimization
       !mesh%optpar=0.05_r8
       mesh%optpar=0.3_r8
       print*, "Optimization parameter:", mesh%optpar
       !Stop criteria - difference between energies
       epsdif=0._r8 !mesh%meanvdist*eps/10.0_r8
       !Stop criteria - minimum energy
       epsmin=mesh%meanvdist*eps*100.0_r8
       !Maximum number of iterations
       maxiter=5000 !6000 !2000
       !Do not stop it energy grows
       upstop=1
       !Select nodes that will be updated
    end select

    !print*, trim(mesh%name), mesh%glevel
    !Write log file
    filename=trim(datadir)//trim(mesh%name)//"_energy"//".txt"
    inquire(file=filename, exist=ifile, opened=iopen, number=iunit2)
    if(iopen)then
       iunit=iunit2
       write(iunit,*)
       write(iunit,*) "n=",mesh%nv, " optm:", trim(mesh%optm)
       write(iunit,*) "energy"
    else
       call getunit(iunit)
       open(iunit,file=filename, status='replace', position='append')
       write(iunit,*) "n=",mesh%nv, " optm:", trim(mesh%optm)
       write(iunit,*) "energy"
    end if

    !Energy variables
    energy=0._r8
    energytmp=0._r8
    inienergy=0._r8

    !Points to be moved
    !From
    ni=1
    !To
    nf=mesh%nv

    !Hierachical optimization (replaces ni, nf)
    ! optimizes only new node points
    if(mesh%hrchy==2 .or. mesh%hrchy==4)then
       !Optimize only new hierachy nodes
       call glevelnodeindex(mesh%glevel, ni, nf, mesh%kind)
       print*, "Level:", mesh%glevel, " Nodes optimized (from to):", ni, nf
    end if

    if(mesh%hrchy==3 .or. mesh%hrchy==4)then

       !Build or load primary icosahedral grid
       icos0%kind=mesh%kind
       icos0%pos=mesh%pos
       icos0%optm="nopt"
       icos0%nv=12
       icos0%loadable=1
       print*
       print*, "Setting primary icosahedral grid for optimization purposes"
       print*, "-----------------------------------------------------------"
       i=0
       if(showonscreen)then
          showonscreen=.false.
          i=1
       end if
       call meshbuild(icos0)
       if(i==1)then
          showonscreen=.true.
       end if

       allocate(nodeslist(1:mesh%nv))
       do i=1, mesh%nv
          !All nodes should be fully optimized
          nodeslist(i)=0
       end do
       !Assume that the nodes that lie on the edges of the primary icosahedral
       !  can only move on the edges, to mantain symmetry
       call icos0edlist(mesh, icos0, nodeslist)
    end if

    !Modify mesh
    do i=1, maxiter

       !Zero energy variable
       energy=0._r8

       !Select what kind of optimization is to be done
       select case(trim(mesh%optm))
       case("scvt") !Spherical Centroidal Voronoi Tesselation
          !Calculate new nodes as barycenters of voronoi cells
          do k=ni, nf
             !In case of locally refined grid
             !  uses a density function (dens_f)
             if(trim(mesh%pos)=="ref")then
                p=vorbarycenterdens(k, mesh)
                dm=arclen(p, mesh%v(k)%p)
                mesh%v(k)%p=p
                energy=energy+dm*dm
             else
                !Uniform density Centroidal Vor Tesselation
                dm=arclen(mesh%hx(k)%b%p, mesh%v(k)%p)
                mesh%v(k)%p=mesh%hx(k)%b%p
                energy=max(energy,dm)
             end if
             !energy=max(energy,arclen(mesh%hx(j)%b%p, mesh%v(j)%p))
             !energy=energy+dm*dm
          end do
          if(trim(mesh%pos)=="ref")then
             energy = dsqrt(energy/(1.0_r8*(nf-ni+1)))
          end if
       case("sprg") !Spring Dynamics

          !Calculate vector forcing terms for each nodes
          call springforcing(vec, mesh, energy)

          !In case of 10 icosahedral symmmety wanted
          if(mesh%hrchy==3)then
             call adjustforsymmetry()
          end if

          !Apply movement on nodes
          do k=ni, nf
             mesh%v(k)%p=mesh%v(k)%p+mesh%optpar*vec%p(k)%v
             mesh%v(k)%p=mesh%v(k)%p/norm(mesh%v(k)%p)
          end do

       case("salt") !Spherical Aligned Tesselation
          !Calculate forcing terms
          call saltforcing(vec, mesh, energy)
          !Move nodes
          !do j=1, sizenlist
          do k=ni, nf
             !k=nodeslist(j)
             mesh%v(k)%p=mesh%v(k)%p+mesh%optpar*vec%p(k)%v
             mesh%v(k)%p=mesh%v(k)%p/norm(mesh%v(k)%p)
          end do

       case("hr95") !Edge alignement optimization
          !Calculate forcing terms
          call hr95forcing(vec, mesh, energy)

          !In case of 10 icosahedral symmmety wanted
          if(mesh%hrchy==3)then
             call adjustforsymmetry()
          end if

          !Move nodes
          do k=ni, nf
             !mesh%v(k)%p=mesh%v(k)%p+mesh%optpar*norm(vec%p(k)%v) *vec%p(k)%v/(mesh%meanvdist) **2
             mesh%v(k)%p=mesh%v(k)%p+mesh%optpar*vec%p(k)%v
             mesh%v(k)%p=mesh%v(k)%p/norm(mesh%v(k)%p)
          end do
       end select

       !Structure/Restructure mesh (only important parts)
       if(trim(mesh%optm)=="scvt" .or. &
            trim(mesh%optm)=="sprg" )then
          call meshstruct(mesh, opt, stat)
       end if

       print*, "Iter, Energy: ", i, energy
       write(iunit,*) i, energy

       !Save initial energy
       if(i==1)then
          inienergy=energy
       end if

       !Save local minimun
       if(energy<=energytmp .and. i>1)then
          do j=ni, nf
             nodesmin%p(j)%v=mesh%v(j)%p
          end do
       elseif(energy>energytmp .and. i>1) then
          if(upstop==1)then
             print*, "OPTIMIZEGRID: Energy grew up!"
             print*, " Using last local minimum."
             print*
             !Use the local minimum saved
             do j=ni, nf
                mesh%v(j)%p=nodesmin%p(j)%v
             end do
             return
          end if
       end if

       !Check if grid maintained it's structure
       if(stat>0)then
          if(trim(mesh%optm) == "scvt") then !Centroidal Voronoi
             !In this case, re-estructuring the mesh is not problematic
             !print*, "SCVT: Neighbourhood relations changed"
          else !In case of SPRG and SALT
             print*, "Warning OPTIMIZEGRID: Mesh needed to be re-structured. "
             print*, " Using last local minimum."
             print*
             !Use the local minimum saved
             do j=ni, nf
                mesh%v(j)%p=nodesmin%p(j)%v
             end do
             return
          end if
       end if

       !Stop if error limit reached
       if(energy<epsmin)then
          print *, " Energy minimum limit reached "
          print '(a, f26.18, a,i8,a,i8,a,f8.4,a)', " Energy: ", energy, &
               " Iterations:", i, " of max ", maxiter, &
               "(", 100.*real(i)/real(maxiter), "% )"
          print*, " Criteria: energy < ", epsmin
          print*
          return
       end if

       !Stop moving too little
       if(abs(energy-energytmp)< epsdif .and. i > miniter )then
          print *, " Two steps with little change reached "
          print '(a, f26.18, a,i8,a,i8,a,f8.4,a)', " Energy: ", energy, &
               " Iterations:", i, " of max ", maxiter, &
               "(", 100.*real(i)/real(maxiter), "% )"
          print*, " Criteria: difference in steps < ", epsdif
          print*, "  Using last local minimum."
          print*
          !Use the local minimum saved
          do j=ni, nf
             mesh%v(j)%p=nodesmin%p(j)%v
          end do
          return
       end if

       !Stop if energy grew too much
       if(energy > inienergy*2.1_r8 )then
          print *, " Energy too high relative to initial energy  "
          print*, " Criteria: energy > ", inienergy*2.1_r8
          print*, "    Using last local minimum."
          do j=ni, nf
             mesh%v(j)%p=nodesmin%p(j)%v
          end do
          print*
          return
       end if
       energytmp=energy

    end do !optimization iterations

    !Max iterations reached
    print*,  " OPTIMIZEGRID WARNING: Max iterations reached - Energy, maxiter: ", energy, maxiter

    !Use the local minimum saved
    print*, "    Using last local minimun."
    do j=ni, nf
       mesh%v(j)%p=nodesmin%p(j)%v
    end do

    close(iunit)

    return
  contains

    subroutine adjustforsymmetry()
      !Correct forcings based on the nodeslist
      do k=1, mesh%nv
         !print*, k,nodeslist(k)
         !Check if node belongs to primary icosahedral edge
         if(nodeslist(k)>0)then !Primary edge points
            !Normal to plane which forcing vector must be projected
            np=cross_product(icos0%v(icos0%ed(nodeslist(k))%v(1))%p, &
                 icos0%v(icos0%ed(nodeslist(k))%v(2))%p)
            np=np/norm(np)
            !Project forcing vector to plane, so that new point will belong to
            !  primary icosahedral edge
            vec%p(i)%v=vec%p(i)%v-dot_product(vec%p(i)%v,np)*np
         elseif(nodeslist(k)<0)then !Pentagon points
            vec%p(i)%v(1:3)=0._r8
         end if
      end do

    end subroutine adjustforsymmetry

  end subroutine optimizegrid

  subroutine springforcing(vec, mesh, energy)
    !--------------------------------------------------------------
    !  SPRINGFORCING
    !   Calculates the spring dynamics optimized forcing vectors
    !   Uses Miuras 2005 iterative algorithm
    !--------------------------------------------------------------

    !Original mesh
    type(grid_structure), intent(in) :: mesh

    !List of forcing vectors
    type(vector_field_cart), intent(inout) :: vec

    !Energy variable
    real (r8):: energy

    !Other energy varibles
    real (r8):: energytmp
    real (r8):: edistor
    real (r8):: emean
    real (r8):: eforce

    !Grid/optimization parameters
    real (r8):: lambda
    real (r8):: beta
    real (r8):: meanl

    real (r8):: nbdist
    !real (r8):: dt

    !Auxiliar vars
    real (r8):: ev(1:3)
    real (r8):: force(1:3)
    real (r8):: distorv(1:3)

    !Indexes
    integer (i4):: j
    integer (i4):: l

    if(mesh%glevel==0)then
       energy=0.0
       return
    end if

    !Tomita's original mean lenght
    lambda=(2._r8*pi)/(10.*(2.**(mesh%glevel-1)))

    !Best beta
    beta=1.2_r8 !1.2_r8

    !Natural spring lenght
    meanl=lambda*beta

    !Forcing vectors
    if(.not.allocated(vec%p))then
       !Allocate space for tmp array of forcing terms
       allocate(vec%p(1:mesh%nv))
       vec%pos=0
       vec%n=mesh%nv
       vec%name="vec"
       !Zero forcing vectors
       do j=1, mesh%nv
          vec%p(j)%v(1:3)=0._r8
       end do
    end if

    !Zero vars
    energy=0.0_r8
    emean=0.0_r8
    edistor=0.0_r8
    eforce=0.0_r8

    !Calculate forcing terms based on each mesh hexagon
    do j=1, mesh%nv
       distorv(1:3)=0._r8
       force(1:3)=0._r8
       energytmp=0.0_r8
       do l=1, mesh%v(j)%nnb

          !Distance to neighbour
          nbdist=arclen(mesh%v(j)%p,mesh%v(mesh%v(j)%nb(l))%p)

          !Unit vector towards the neighbour
          ev=mesh%v(mesh%v(j)%nb(l))%p-mesh%v(j)%p

          !Forcing vector
          force=force+(nbdist-meanl)*ev/nbdist

          !Forcing energy - Miura
          emean=emean+((nbdist-meanl)**2)/2

       end do

       !Forcing energy - Tomita/Miura
       eforce=max(eforce, norm(force))

       !Calculate distortion
       ! Mesh needs to be restructured every step to use this
       edistor=max(edistor, distortionhx(j, mesh))

       !Forcing term Miura
       vec%p(j)%v=force

    end do

    !Select energy to be used
    !enerfy=eforce
    !if(mesh%glevel<7)then
    !  energy=emean
    !else
    energy=edistor
    !end if
    !print*, mesh%glevel, energy, emean, edistor

    return

  end subroutine springforcing

  subroutine saltforcing(vec, mesh, energy)
    !--------------------------------------------------------------
    !  SALTFORCING
    !   Calculates the aligned optimized forcing vectors
    !--------------------------------------------------------------
    !Original mesh
    type(grid_structure), intent(in) :: mesh

    !List of forcing vectors
    type(vector_field_cart), intent(inout) :: vec

    !Energy varible
    real (r8):: energy

    !Temp vec variable
    real (r8):: v1(1:3)
    real (r8):: v2(1:3)
    real (r8):: v(1:3)

    !Indexes
    integer (i4):: j
    integer (i4):: l
    integer (i4):: m

    if(.not.allocated(vec%p))then
       !Allocate space for tmp array of forcing terms
       allocate(vec%p(1:mesh%nv))
       vec%pos=0
       vec%n=mesh%nv
       vec%name="vec"
    end if

    !Zero forcing vectors
    do j=1, mesh%nv
       vec%p(j)%v(1:3)=0._r8
    end do

    energy=0._r8

    !Calculate forcing terms based on each mesh hexagon
    !  Pentagons are forced indirectly by the surrounding hexagons
    do j=1, mesh%nv
       if(mesh%v(j)%nnb==6)then
          do l=1, 3
             m=modint(l+4, 6)
             v1=mesh%v(mesh%v(j)%nb(m))%p-mesh%v(mesh%v(j)%nb(l))%p
             v2=mesh%v(mesh%v(j)%nb(l+3))%p-mesh%v(mesh%v(j)%nb(l+1))%p
             !Main forcing vector
             v=v2-v1
             !Individual forcing vector
             vec%p(mesh%v(j)%nb(l))%v   = vec%p(mesh%v(j)%nb(l))%v    &
                  - proj_vec_sphere(v, mesh%v(mesh%v(j)%nb(l))%p)
             vec%p(mesh%v(j)%nb(l+1))%v = vec%p(mesh%v(j)%nb(l+1))%v  &
                  + proj_vec_sphere(v, mesh%v(mesh%v(j)%nb(l+3))%p)
             vec%p(mesh%v(j)%nb(l+3))%v = vec%p(mesh%v(j)%nb(l+3))%v  &
                  - proj_vec_sphere(v, mesh%v(mesh%v(j)%nb(l+3))%p)
             vec%p(mesh%v(j)%nb(m))%v   = vec%p(mesh%v(j)%nb(m))%v    &
                  + proj_vec_sphere(v, mesh%v(mesh%v(j)%nb(m))%p)
             !Forcing energy
             !energy=max(energy, norm(v))
          end do
       end if
       energy=energy+ alignind(j, mesh)**2
    end do
    energy=dsqrt(energy)/mesh%nv

    !Make forcings tangent to the sphere
    do j=1, mesh%nv
       vec%p(j)%v=proj_vec_sphere(vec%p(j)%v, mesh%v(j)%p)
    end do

    return
  end subroutine saltforcing

  subroutine hr95forcing(vec, mesh, energy)
    !--------------------------------------------------------------
    !  EDGE ALIGN FORCING
    !   Calculates the aligned triangle/Voronoi edges forcing vectors
    !--------------------------------------------------------------
    !Original mesh
    type(grid_structure), intent(in) :: mesh

    !List of forcing vectors
    type(vector_field_cart), intent(inout) :: vec

    !Energy varible, mesh edge miss alignement
    real (r8):: energy
    real (r8):: r

    !Temp vec variable
    real (r8):: v(1:3)
    real (r8):: pmed(1:3)
    real (r8):: pmedhx(1:3)
    real (r8), dimension(1:3) :: ptrc1
    real (r8), dimension(1:3) :: ptrc2

    !Indexes
    integer (i4):: j
    integer (i4):: l
    integer (i4):: k
    integer (i4):: k1
    integer (i4):: k2
    integer (i4):: i1
    integer (i4):: i2
    integer (i4):: i3


    if(.not.allocated(vec%p))then
       !Allocate space for tmp array of forcing terms
       allocate(vec%p(1:mesh%nv))
       vec%pos=0
       vec%n=mesh%nv
       vec%name="vec"
    end if

    energy=0._r8

    !Calculate forcing terms based on each mesh hexagon
    do j=1, mesh%nv
       vec%p(j)%v(1:3)=(/ 0._r8, 0._r8, 0._r8 /)
       !print*, j, mesh%v(j)%tr(1:mesh%v(j)%nnb)
       !if(j==29)then
       ! print*, j, "nbtr:", mesh%v(j)%tr(1:mesh%v(j)%nnb)
       !end if
       do l=1, mesh%v(j)%nnb

          !Forcing vector
          k=mesh%v(j)%ed(l)
          !r=arclen(mesh%edhx(k)%c%p,mesh%ed(k)%c%p)/mesh%edhx(k)%leng
          !v=(mesh%edhx(k)%c%p-mesh%ed(k)%c%p)

          k1=mesh%v(j)%tr(modint(l-1, mesh%v(j)%nnb))
          i1=mesh%tr(k1)%v(1)
          i2=mesh%tr(k1)%v(2)
          i3=mesh%tr(k1)%v(3)
          !Calculate the circumcenter of triangles
          ptrc1=trcircumc ( mesh%v(i1)%p,  mesh%v(i2)%p,  mesh%v(i3)%p)

          k2=mesh%v(j)%tr(l)
          i1=mesh%tr(k2)%v(1)
          i2=mesh%tr(k2)%v(2)
          i3=mesh%tr(k2)%v(3)
          ptrc2=trcircumc ( mesh%v(i1)%p,  mesh%v(i2)%p,  mesh%v(i3)%p)

          !Calculate mean point of edges
          pmedhx=(ptrc1+ptrc2)/2._r8
          pmedhx=pmedhx/norm(pmedhx)
          pmed=(mesh%v( mesh%v(j)%nb(l) )%p+mesh%v(j)%p)/2._r8
          pmed=pmed/norm(pmed)

          v=pmedhx-pmed
          !v=-v
          r=arclen(pmedhx, pmed)/arclen(ptrc1, ptrc2)

          !r=arclen(pmedhx,pmed)/arclen(mesh%tr(k2)%c%p,mesh%tr(k1)%c%p)
          !print*, j, l, r, norm(v)
          !v=v/norm(v)
          if(j>182 .and. j <186)then
             print '(4i6, 3f20.14)', j, mesh%v(j)%nb(l), k1, k2, &
                  arclen(pmedhx,pmed), arclen(mesh%edhx(k)%c%p, mesh%ed(k)%c%p), &
                  arclen(mesh%v( mesh%v(j)%nb(l) )%p,mesh%v(j)%p)
          end if
          !print*, arclen(pmedhx,pmed), arclen(mesh%edhx(k)%c%p, mesh%ed(k)%c%p)
          !if(dot_product(v, (mesh%edhx(k)%c%p-mesh%ed(k)%c%p)) > 0) then
          ! dcor=1._r8
          !else
          ! dcor=-1._r8
          !end if
          !vec%p(j)%v=vec%p(j)%v+dcor*r*v
          vec%p(j)%v=vec%p(j)%v+v !(r)*v/norm(v)



          !Energy
          !energy=max(energy, r ) ! + (r)**4
          energy = energy + (r)**4
       end do
    end do
    !energy=energy !/mesh%nv

    return
  end subroutine hr95forcing

  !================================================================================
  !    DENSITY FUNCTION - TO BE USED IN CASE OF LOCALLY REFINED GRIDS
  !================================================================================

  function dens_f(p, refin)
    !--------------------------------------------------------------------
    ! dens_f
    !     Define the density function for a given point on the sphere
    !       to be used in case of locally refined grids
    !
    !     Is to used together with SCVT optimization
    !     Densities must be in [0,1]
    !     Need manual adjustments to define refined region
    !
    !     refin = "ref" gives refined mesh, otherwise uniform
    !        (refin may be passed as mesh%pos)
    !--------------------------------------------------------------------
    real (kind=8), intent(in) :: p(1:3)
    character(len=3), intent(in) :: refin
    real (kind=8):: dens_f

    real (kind=8):: lat
    real (kind=8):: lon
    real (kind=8):: latc
    real (kind=8):: lonc
    real (kind=8):: gammas
    real (kind=8):: epsilons
    real (kind=8):: dists
    real (kind=8):: maxdist
    real (kind=8):: sx

    !Set this flag to choose density function
    if(refin/="ref")then
       dens_f=1
    else
       call cart2sph(p(1), p(2), p(3), lon, lat)

       !See Voronoi Tessellations and their application to climate and global modeling
       ! by Lili Ju, Todd Ringler and Max Gunzburger
       !Center of refined region
       latc=-50.0*pi/180._r8
       lonc=-30.0*pi/180._r8
       !Distance to center
       dists=dsqrt((lat-latc)**2+(lon-lonc)**2)
       !Density function parameters
       gammas=3._r8
       epsilons=pi/12._r8
       maxdist=pi/6._r8
       !Distance to center metric
       sx=(dists-maxdist)/epsilons
       !Set density
       if(dists<=maxdist)then
          !Point close to the center
          dens_f=gammas**4
       elseif((dists<=maxdist+epsilons) .and. (dists>maxdist))then
          !Point in the transition
          dens_f=((1._r8-sx)*gammas+sx)**4
       else
          !Point far from the center
          dens_f=1
       end if
       !Normalization - Make it in [0,1]
       dens_f=dens_f/gammas**4

    end if
    return

  end function dens_f

  !================================================================================
  !    GEODESIC TO REGULAR GRID ROUTINES
  !================================================================================

  subroutine geo2reg(mesh)
    !----------------------------------------------------------
    !GEO2REG
    !
    ! Generates a relationship list between the triangulation
    !   and a regular mesh
    ! Saves the results in files for future loads
    !--------------------------------------------------------------
    type(grid_structure), intent(inout) :: mesh   !Input and output

    !Aux variables
    !integer (i4):: nlat
    real (r8):: sqlat
    real (r8):: sqlon
    real (r8):: sqlatc
    real (r8):: sqlonc
    real (r8):: sqcp(1:3)

    !Temporary list for intersecting triangles
    !   It is set to a maximum of 20 intersection of
    !    triangles with a rectangle/square
    !   For icosahedral meshes the maximum is 6
    integer (i4), dimension(1:20) :: trlist
    !integer (i4), allocatable :: trlist(:)

    !Sizes for each square and triangle
    real(r8):: nodesqdist
    real(r8):: sqsize
    real(r8):: trsize
    real(r8):: trsqdist
    real(r8):: mindist
    real(r8):: trsqdistmin
    integer(i4):: intersec
    integer(i4):: isqtr
    integer(i4)::itmp

    !Auxiliar variables
    integer (i4):: i
    integer (i4):: j
    integer (i4):: k
    integer (i4):: l
    integer (i4):: ll
    integer (i4):: kt
    integer (i4):: kk
    integer (i4):: node
    integer (i4):: ercount
    logical:: possibleint

    !OpenMP variables
    !integer :: n_threads, id_thread, chunk,  OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM

    !print*
    !print*,"Creating geodesic to regular grid structure..."
    !print*

    !Temp vars
    !nlat=mesh%nlat

    !Print progress of calculations
    if(mesh%nv>11000 .and. mesh%nlat>100 )then
       write(*, '(a)') "   Progress (5% per mark): "
       write(*,'(a2, $)') "  "
       !write(*, '(a)', advance='no') "    "
    end if

    !Allocate Latitudes in table
    if(.not.allocated(mesh%qd_lat))then
       allocate(mesh%qd_lat(1:mesh%nlat))
    else
       !If the new array has a different size, reallocate
       if(ubound(mesh%qd_lat,1)/=mesh%nlat)then
          deallocate(mesh%qd_lat)
          allocate(mesh%qd_lat(1:mesh%nlat))
       end if
    end if

    !Allocate longitudes in table
    do j=1, mesh%nlat
       mesh%qd_lat(j)%nlon=nlon_redgrd(j,mesh%nlat)
       allocate(mesh%qd_lat(j)%qd_lon(1:mesh%qd_lat(j)%nlon))
       mesh%qd_lat(j)%dlon=pi2/real(mesh%qd_lat(j)%nlon, r8)
    end do

    mesh%maxtrsqint=0

    !Count errors
    ercount=0

    !OPENMP PARALLEL REGION
    !$omp parallel  shared(mesh, ercount)  &
    !$omp default(private)

    !allocate(trlist(1:20))

    !OPENMP Debug
    !n_threads=OMP_GET_NUM_THREADS()
    !id_thread=OMP_GET_THREAD_NUM()
    !print*, "Thread ",id_thread+1, " of ", n_threads, " started..."

    !OPENMP PARALLEL DO REGION
    !$omp do schedule(dynamic, 10)
    do j=1,mesh%nlat
       !do j=5,5 !26,26 !31,31
       !print*, "Thread ",id_thread+1, " j ", j, ercount
       !Square corner lat
       sqlat=-pio2+(j-1)*mesh%dlat
       !Square center
       sqlatc=sqlat + mesh%dlat/2
       !Save in mesh
       mesh%qd_lat(j)%lat=sqlat
       !sqlon=-180.0_r8*deg2rad
       !print*, j, nlon_redgrd(j, mesh%nlat), mesh%nlat
       do i=1,mesh%qd_lat(j)%nlon
          !do i=3,3 !211,211
          !Square corner lon
          sqlon=-pi+(i-1)*mesh%qd_lat(j)%dlon
          !Square center
          sqlonc=sqlon + mesh%qd_lat(j)%dlon/2
          !Save in mesh
          mesh%qd_lat(j)%qd_lon(i)%lon=sqlon

          !Find nearest node to square center
          call sph2cart(sqlonc, sqlatc, sqcp(1), sqcp(2), sqcp(3))
          node=nearnodesearch(sqcp, node, mesh)

          !zero counters and flags
          intersec=0
          trsqdistmin=10000._r8
          mindist=100000._r8
          nodesqdist=10000._r8
          itmp=0

          !Check if node is in square (slightly enlarged square)
          if(mesh%v(node)%lon>=sqlon-mesh%qd_lat(j)%dlon/3._r8 .and.  &
               mesh%v(node)%lon <= sqlon+mesh%qd_lat(j)%dlon+mesh%qd_lat(j)%dlon/3._r8 &
               .and. mesh%v(node)%lat>=sqlat-mesh%dlat/3._r8 &
               .and. mesh%v(node)%lat <= sqlat+mesh%dlat+mesh%dlat/3._r8)then
             !In this case, all triangles around the node intersect the square
             do l=1, mesh%v(node)%nnb
                k=mesh%v(node)%tr(l)
                trlist(l)=k
                !Get triangle with closest barycenter
                trsqdist=arclen(mesh%tr(k)%b%p, sqcp)
                if(trsqdist<=trsqdistmin)then
                   trsqdistmin=trsqdist
                   !kt is the closest triangle to the center of the square
                   kt=k
                   !index of the closest triangle position in trlist
                   ll=l
                end if
             end do
             intersec=mesh%v(node)%nnb
             !Set closest triangle as number 1 in trlist
             !  changing it with the first
             kk=trlist(1)
             trlist(1)=kt
             trlist(ll)=kk
          else
             ! The node is not in the square

             !Get nearest triangle (barycenter)
             do l=1, mesh%v(node)%nnb
                k=mesh%v(node)%tr(l)
                !Get distance from square center to triangle barycenter !circumcenter
                !trsqdist=arclen(mesh%tr(k)%c%p, sqcp)
                trsqdist=arclen(mesh%tr(k)%b%p, sqcp)
                if(trsqdist<=trsqdistmin)then
                   trsqdistmin=trsqdist
                   kt=k
                end if
             end do
             ! kt will be the first element on the triangle list
             !   because it is probably the one with most intersection
             !   with the rectangle
             intersec=intersec+1
             trlist(intersec)=kt

             ! Check the surronding triangles of "node"
             !Max arc dist from square center to nodes
             sqsize=max(arcdistll(sqlonc, sqlatc, sqlon, sqlat), &
                  arcdistll(sqlonc, sqlatc, sqlon+mesh%qd_lat(j)%dlon, &
                  sqlat+mesh%dlat))
             !The square is not convex on the sphere with geodesic metric, so
             ! the square radius must be slightly increase
             sqsize=1.2_r8*sqsize

             !For all triangles surrounding the node
             nodetrs: do l=1, mesh%v(node)%nnb
                !Get triangle index
                k=mesh%v(node)%tr(l)
                isqtr=0
                if(k==kt)then
                   !Check if sqaure is totaly in the closest triangle
                   if(sqtriintersec(i, j, k, mesh)==2)then
                      !square is totally in kt triangle
                      intersec=1
                      trlist(intersec)=kt
                      exit nodetrs
                   else
                      !kt is already on the list, as the closest triangle,
                      ! so it doesn't need to be re-checked for intersection
                      cycle
                   end if
                end if
                !Get distance from square center to triangle circumcenter
                trsqdist=arclen(mesh%tr(k)%c%p, sqcp)
                !Get triangle circumradius + the distance to the barycenter
                !  This is an enlarged triangles circumcircle, to avoid round off errors
                trsize=mesh%tr(k)%radg+arclen(mesh%tr(k)%b%p, mesh%tr(k)%c%p)
                !See if the squares circumcircle intersects the triangle's circumcircle
                possibleint= trsqdist < trsize+sqsize
                !If there is a possible intersection
                if(possibleint) then
                   !Test for intersection
                   !print*, mesh%tr(k)%v
                   if(abs(sqlat) < 60.*deg2rad .and. &
                        abs(sqlon) < 176.*deg2rad)then
                      !If square in the tropics, the square is well defined
                      !  so check for intersection with triangle
                      isqtr=sqtriintersec(i, j, k, mesh)
                   else
                      !If square near the pole, then just assume that
                      ! the intersection exists
                      isqtr=1
                   end if

                   !In case of existing intersection
                   if(isqtr>0)then
                      !Add triangle to list
                      intersec=intersec+1
                      trlist(intersec)=k
                   end if
                end if
             end do nodetrs
             !print *, id_thread+1, i, j, kt, intersec

          end if

          !As the rectangles are much smaller then the triangles, the
          !  intersections obtained up until now are generaly enough
          !But some special cases it might fail, so in order to
          !  guarantee the list, do some further checking

          !Check neighbour triangles
          do kk=1,intersec
             kt=trlist(kk)
             neib: do l=1, 3
                k=mesh%tr(kt)%nb(l)
                !Check if k is already on the list
                do ll=1,intersec
                   if(k==trlist(ll))then
                      cycle neib
                   end if
                end do
                isqtr=0
                !Get distance from square center to triangle center
                trsqdist=arclen(mesh%tr(k)%c%p,sqcp)
                !Get triangle circumradius
                !trsize=mesh%tr(k)%radg
                trsize=mesh%tr(k)%radg+arclen(mesh%tr(k)%b%p, mesh%tr(k)%c%p)
                !See if the square's circumcircle intersects the triangle's circumcircle
                possibleint=trsqdist<trsize+sqsize
                !If there is intersection, there might be an intersection between
                !  square and triangle
                if(possibleint) then
                   !Test for intersection
                   if(abs(sqlat) < 60.*deg2rad .and. &
                        abs(sqlon) < 176.*deg2rad)then
                      !If square in the tropics, the square is well defined
                      !  so check for intersection with triangle
                      isqtr=sqtriintersec(i, j, k, mesh)
                   else
                      !If square near the pole, then just assume that
                      ! the intersection exists
                      isqtr=1
                   end if

                   !In case of existing intersection, add triangle to list
                   if(isqtr>0)then
                      if(isqtr==2)then
                         !Found that the square is totally inside triangle k
                         !  Put is as first tr
                         intersec=intersec+1
                         trlist(intersec)=trlist(1)
                         trlist(1)=k
                         !exit !stop searching
                      else !just add the triangle to the list
                         intersec=intersec+1
                         trlist(intersec)=k
                      end if
                   end if
                end if
             end do neib
          end do

          if(intersec==0 )then
             if(ercount==0)then
                print*, "Sqlon  sqlat  sqcenterlon sqcenterlat  node nodelon nodelat "
                print '(4f10.4,i6,2f10.4)', sqlon*rad2deg, sqlat*rad2deg, sqlonc*rad2deg, sqlatc*rad2deg,node, &
                     mesh%v(node)%lon*rad2deg, mesh%v(node)%lat*rad2deg
                print*, "GEO2REG ERROR: Uniform mesh square without triangle intersection!!"
                ercount=ercount+1
             else
                ercount=ercount+1
             end if
          end if

          !Update the maximum of intersections for each square
          ! In icosahedral meshes this must be less then 6 if squares suficiently small
          mesh%maxtrsqint=max(intersec,mesh%maxtrsqint)
          !print*,i,j,intersec, maxintersec
          !maxsqsize=max(maxsqsize,sqsize)

          !Save intersections in grid structure
          mesh%qd_lat(j)%qd_lon(i)%ntr=intersec
          !print *, id_thread+1, i, j, intersec

          if(intersec==0)then
             if(.not.allocated(mesh%qd_lat(j)%qd_lon(i)%tr))then
                allocate(mesh%qd_lat(j)%qd_lon(i)%tr(1:intersec+1))
             end if
             mesh%qd_lat(j)%qd_lon(i)%tr(1)=0
          else
             if(.not.allocated(mesh%qd_lat(j)%qd_lon(i)%tr))then
                allocate(mesh%qd_lat(j)%qd_lon(i)%tr(1:intersec))
             elseif(ubound(mesh%qd_lat(j)%qd_lon(i)%tr, 1)/=intersec)then
                deallocate(mesh%qd_lat(j)%qd_lon(i)%tr)
                allocate(mesh%qd_lat(j)%qd_lon(i)%tr(1:intersec))
             end if
             mesh%qd_lat(j)%qd_lon(i)%tr(1:intersec)=trlist(1:intersec)
          end if

       end do
       !Print progress of calculations
       if(mesh%nv>11000 .and. mesh%nlat>100 )then
          if(mod(j,mesh%nlat/20)==0)then
             write(*,'(a1, $)') "*" !print*, bar(1:12)
          end if
       end if

    end do
    !$omp end do

    !Openmp debug
    !print*, "Thread ",id_thread+1, " of ", n_threads, " ended!"
    !!omp barrier
    !deallocate(trlist)

    !$omp end parallel 

    !Error check
    if(ercount>0)then
       print*
       print*, "  Number of errors encountered: ", ercount
       print*, "  Avoid using the table based triangle search!!!"
    end if

    print*
    print*,"Geodesic to regular grid structure created."
    print*

    return

  end subroutine geo2reg

  function sqtriintersec(ilon, jlat, kt, mesh)
    !----------------------------------------------------------
    !
    !	SQTRIINTERSEC
    !
    !	Verifies the intersection between geodesic squares 
    !		and triangles for a given geodesic mesh
    !   Uses the triangle and square edges to test for intersection.
    !   Garantees intersection, if it does not exist, it will return zero.
    !
    !	Recieves the (ilat,ilon) as indexes of a square
    !		and the number of the triangle (kt)
    !	Returns 0 if there is no intersection
    !	Returns 1 or 2 if there is intersection
    !   Returns 2 if square totally inside triangle.
    !
    !
    ! Warning: If the square is near the pole, then it is too narrow
    !   and the result may be incorrect
    !--------------------------------------------------------------	
    !Input Variables
    integer (i4), intent(in) :: ilon
    integer(i4), intent(in) :: jlat
    real (r8) :: sqlon
    real (r8) :: sqlat
    integer(i4), intent(in) :: kt
    type(grid_structure), intent(in) :: mesh

    !Output Variables
    integer:: sqtriintersec

    !Aux variables
    real(r8):: sqlatv(4)
    real(r8):: sqlonv(4)
    real(r8):: sqxyz(4,3)
    !real(r8):: trp(3,3)
    real (r8), dimension(1:3):: p1
    real (r8), dimension(1:3):: p2
    real (r8), dimension(1:3):: q1
    real (r8), dimension(1:3):: q2
    !real (r8) :: sqcenterlat, sqcenterlon, sqcenter(1:3)
    logical:: isintr
    integer(i4):: i
    integer(i4):: j
    integer(i4):: sqcount

    sqtriintersec=0
    !ind(1:3)=(/1,2,3/)

    !Left bottom square corner
    sqlat=mesh%qd_lat(jlat)%lat
    sqlon=mesh%qd_lat(jlat)%qd_lon(ilon)%lon

    !Square nodes, enlarge the square to ensure
    ! that the non convexity of the
    ! square does not create problems
    sqlatv(1)=sqlat-mesh%dlat/3._r8
    sqlonv(1)=sqlon-mesh%qd_lat(jlat)%dlon/3._r8
    sqlatv(2)=sqlat-mesh%dlat/3._r8
    sqlonv(2)=(sqlon+mesh%qd_lat(jlat)%dlon)+mesh%qd_lat(jlat)%dlon/3._r8
    sqlatv(3)=(sqlat+mesh%dlat)+mesh%dlat/3._r8
    sqlonv(3)=(sqlon+mesh%qd_lat(jlat)%dlon)+mesh%qd_lat(jlat)%dlon/3._r8
    sqlatv(4)=(sqlat+mesh%dlat)+mesh%dlat/3._r8
    sqlonv(4)=sqlon-mesh%qd_lat(jlat)%dlon/3._r8

    !Ensure that square is in boundaries -180, 180 ; -90;90
    do i=1,4
       if(sqlatv(i)<=-pio2)then
          sqlatv(i)=-pio2
       elseif(sqlatv(i)>=pio2)then
          sqlatv(i)=pio2
       endif
       if(sqlonv(i)<=-pi)then
          sqlonv(i)=-pi
       elseif(sqlonv(i)>=pi)then
          sqlonv(i)=pi
       endif
    end do

    !Convert to cartesian coords
    !(vectorial diretive optimization)
    !dir$ ivdep
    do i=1,4
       call sph2cart(sqlonv(i), sqlatv(i), sqxyz(i,1), sqxyz(i,2), sqxyz(i,3))
    end do

    !Check if the square is totaly inside the triangle

    !For each sq node
    sqcount=0
    do i=1,4
       !Check if it is inside the triangle
       !isintr=insidetr( sqxyz(i,1:3), &
       !   mesh%v(mesh%tr(kt)%v(1))%p,  &
       !   mesh%v(mesh%tr(kt)%v(2))%p,  &
       !   mesh%v(mesh%tr(kt)%v(3))%p )
       isintr=insidetrmesh(sqxyz(i,1:3), kt, mesh)
       !print '(a,2f8.3,a,i4,2f8.3,a,l,i2)', "SQPOINT:", &
       !	tlat(i)*rad2deg, tlon(i)*rad2deg, " inside TR:", kt,  &
       !	mesh%clat(kt)*rad2deg, mesh%clon(kt)*rad2deg, " inside? ", isintr, ier
       if(isintr) then
          sqtriintersec=1
          sqcount=sqcount+1
       end if
    end do

    if(sqcount>=4) then
       !The square is totally in the triangle
       sqtriintersec=2
       return
    end if

    if(sqtriintersec==1)then
       !There is a square point in the triangle
       !print*,"Point inside!!"
       return
    end if

    !Test for intersecting edges
    !For each square edge
    do i=1,4
       p1=sqxyz(i,1:3)
       p2=sqxyz(modint(i+1,4),1:3)
       !For each triangle edge
       do j=1,3
          q1=mesh%v(mesh%tr(kt)%v(j))%p
          q2=mesh%v(mesh%tr(kt)%v(modint(j+1,3)))%p
          !Check if there is intersection between edges
          if(arcintersec(p1,p2,q1,q2))then
             sqtriintersec=1
             return
          end if
       end do
    end do

    return

  end function sqtriintersec

  !====================================================================
  !   SEARCH TOOLS
  !====================================================================

  function nearnodesearch(p, inode, mesh)
    !----------------------------------------------------------
    !	NEARNODESEARCH
    !
    !	Given a point 'p=(x,y,z)' on the sphere  returns 
    !		the nearest node in the mesh to this point
    !           using linear search methods
    !           
    !   'inode' is the initial near point guess to begin search
    !--------------------------------------------------------------	
    !Input variables
    real (r8), intent(in) :: p(1:3)
    integer (i4), intent(inout) :: inode
    type(grid_structure), intent(in) :: mesh

    !Output variable
    integer (i4):: nearnodesearch

    !Auxiliar Variables
    integer (i4):: i
    integer (i4):: j
    integer (i4):: l
    integer (i4):: nearnode
    integer (i4):: reached
    integer (i4):: newnode
    real (r8):: dmin
    real (r8):: dtmp

    reached=0
    if(inode>1.and. inode<mesh%nv)then
       nearnode=inode
    else
       nearnode=1
    end if
    dmin=arclen(mesh%v(nearnode)%p, p)

    !Iterate for at most all mesh%nv nodes
    iter: do i=1, mesh%nv 
       !Check distances on 'nearnode' neighbours
       newnode=0
       do j=1, mesh%v(nearnode)%nnb
          l=mesh%v(nearnode)%nb(j)
          dtmp=arclen(mesh%v(l)%p, p)
          !If this neighbor is closer to p, save it
          if(dtmp<dmin)then
             newnode=l
             dmin=dtmp
          end if
       end do
       if(newnode==0)then
          !In this case, the nearnode is the nearest node to p
          !  so we can stop the search
          exit iter
       else
          !In this case, a new node is to be set,
          !   and the process is restarted
          nearnode=newnode
       end if
    end do iter
    if(dmin>mesh%maxvdist)then
       print*, "NEARNODESEARCH WARNING: node found too distant from target"
       print*, "dist:",  dmin, " maxvdist for mesh: ", mesh%maxvdist
    end if
    nearnodesearch=nearnode

    return

  end function nearnodesearch


  subroutine getnearnodes(p, mesh, r, listv, listd, n)
    !----------------------------------------------------------
    !	GETNEARNODES
    !
    !	Given a point p returns 
    !		a list of the nearest nodes in the mesh to this point
    !           using the regular to geodesic structure using r or n
    !    - 'r' is a radius that defines the region of interest
    !         all nodes from mesh in this region will be returned.
    !    - The routines returns listv as an array of node indexes, and sets
    !       n to the size of the list.
    !    - 'n' is the number of nearest nodes to points to be found, if
    !       the number of nodes are to be defined by 'r', then call the routine
    !       with n=0. At the end, 'n' will contain the number os node index 
    !       on the list.
    !    - If r=0,  returns the n nearest nodes.
    !      If n=0 and r=0, returns the nearst node to point
    !    - The list returned is in order of proximity
    !
    !   Uses lat/lon search table
    !--------------------------------------------------------------	
    !Maximum number of elements on the list
    integer (i4), parameter :: nmax=500

    !Input variables
    real (r8), intent(in) :: p(1:3)

    !Mesh
    type(grid_structure), intent(in) :: mesh

    !Output list of node index
    integer (i4), allocatable :: listv(:)

    !Output list of node distances
    real (r8), allocatable :: listd(:)

    !Temporary lists
    integer (i4):: listvtmp(1:nmax)
    real (r8):: listdtmp(1:nmax)

    !Radius of region
    real (r8), intent(in) :: r
    real (r8):: rmax

    !Number of nodes
    integer (i4), intent(inout) :: n

    !Auxiliar Variables
    integer (i4):: i
    integer (i4):: ntmp
    integer (i4):: node
    integer (i4):: node1
    real (r8):: d


    !Check for consistent input
    if(n>nmax)then
       print*, "GETNEARNODES WARNING: Too many neighbours wanted", n
       print*, "  Using the maximum: ", nmax
       n=nmax
    end if

    !Get nearest node
    node1=getnearnode(p, mesh)

    !Start the temporary list
    listvtmp(1)=node1
    listdtmp(1)=arclen(mesh%v(node1)%p, p)
    ntmp=1

    !Change radius in case of radius given too small
    rmax=r
    if(rmax<eps)then
       rmax=999._r8
    end if

    !Get nearest nodes to point
    do i=2, nmax

       !Get next nearest node
       node=nextnearestnode(p, mesh, i-1, listvtmp, d)
       !print*, i, node, d
       !Check if radius reached or maximum number of nodes
       if(d>rmax .or. i==mesh%nv/2 .or. d>pio2)then
          n=i-1 !The last node is not required in the list since d>r
          exit
       elseif(i==nmax .or. i==n )then
          n=i
          !Add last node to list
          listvtmp(i)=node
          listdtmp(i)=d
          exit
       end if

       !Add new node to list
       listvtmp(i)=node
       listdtmp(i)=d

    end do

    if((.not.allocated(listv)))then
       allocate(listv(1:n))
       allocate(listd(1:n))
    elseif(size(listv)<n)then
       print*, "NEARNODES WARNING: Received allocated array with"
       print*, "   too little space. Returning up until ", size(listv)
       listd(1:size(listv))=listdtmp(1:size(listv))
       listv(1:size(listv))=listvtmp(1:size(listv))
       return
    end if

    listd(1:n)=listdtmp(1:n)
    listv(1:n)=listvtmp(1:n)

    return

  end subroutine getnearnodes


  function nextnearestnode(p, mesh, n, listv, dist)
    !----------------------------------------------------------
    !	NEXTNEARESTNODE
    !
    !	Given a point p returns 
    !		the index of the next nearest node for a given list 'listv'
    !       of the 'n' already nearest nodes in order
    !
    !--------------------------------------------------------------	
    !Input variables
    real (r8):: p(1:3)

    !Mesh
    type(grid_structure), intent(in) :: mesh

    !Number of nodes already in list
    integer (i4), intent(in) :: n

    !Output array
    integer (i4), intent(in) :: listv(1:n)

    !Distance to point
    real (r8), intent(out) :: dist

    !Next nerest node
    integer (i4):: nextnearestnode

    !Auxiliar Variables
    integer (i4):: i
    integer (i4):: j
    integer (i4):: k

    integer (i4):: node
    integer (i4):: nodenb
    integer (i4):: nearestnode
    real (r8):: dmin
    real (r8):: d

    !Check if list is reliable
    if(listv(n)<=0 .or. listv(n)>mesh%nv )then
       print*, "NEXTNEARESTNODE ERROR: Listv incorect", listv(1:n)
       stop
    end if

    !  For all members of the list check for the nearest neighbour that
    !    is closer to the the given point
    dmin=999._r8
    listmembers: do i = 1, n
       !Get node index
       node=listv(i)
       !For all its neighbours
       neighbours: do j = 1, mesh%v(node)%nnb
          !Get node neighbour index
          nodenb=mesh%v(node)%nb(j)

          !Check if already on the list
          do k=1, n
             if(nodenb==listv(k))then
                cycle neighbours
             end if
          end do

          !Calculate distance to point
          d = arclen(mesh%v(nodenb)%p, p)
          if(d < dmin)then
             !New nearest node!
             nearestnode=nodenb
             dmin=d
          end if
       end do neighbours
    end do listmembers

    if(dmin > 2._r8)then
       !No new node to include
       nextnearestnode=0
       dist=dmin
    else
       !Return new node
       nextnearestnode=nearestnode
       dist=dmin
    end if

    return

  end function nextnearestnode

  subroutine getnearhxedges(p, mesh, r, listed, listd, n)
    !----------------------------------------------------------
    !	GETNEARHXEDGES
    !
    !	Given a point p returns 
    !		a list of the nearest voronoi edges midpoints in the mesh to this point
    !           using the regular to geodesic structure using r or n
    !    - 'r' is a radius that defines the region of interest
    !         all nodes from mesh in this region will be returned.
    !    - The routines returns listv as an array of node indexes, and sets
    !       n to the size of the list.
    !    - 'n' is the number of nearest nodes to points to be found, if
    !       the number of nodes are to be defined by 'r', then call the routine
    !       with n=0. At the end, 'n' will contain the number os node index 
    !       on the list.
    !    - If r=0,  returns the n nearest nodes.
    !      If n=0 and r=0, returns the nearst node to point
    !    - The list returned is in order of proximity
    !
    !   Uses lat/lon search table
    !--------------------------------------------------------------	
    !Maximum number of elements on the list
    integer (i4), parameter :: maxn=500

    !Input variables
    real (r8), intent(in) :: p(1:3)

    !Mesh
    type(grid_structure), intent(in) :: mesh

    !Output list of edge index
    integer (i4), allocatable :: listed(:)

    !Output list of node distances
    real (r8), allocatable :: listd(:)

    !Temporary lists
    integer (i4):: listedtmp(1:maxn)
    real (r8):: listdtmp(1:maxn)

    !Radius of region
    real (r8), intent(in) :: r
    real (r8):: rmax

    !Number of nodes
    integer (i4), intent(inout) :: n

    !Auxiliar Variables
    integer (i4):: i
    integer (i4):: ntmp
    integer (i4):: ed
    integer (i4):: ed1
    real (r8):: d

    !Check for consistent input
    if(n>maxn)then
       print*, "GETNEARNODES WARNING: Too many neighbours wanted", n
       print*, "  Using the maximum: ", maxn
       n=maxn
    end if

    !Get nearest node
    ed1=getnearhxedge(p, mesh)

    !Start the temporary list
    listedtmp(1)=ed1
    listdtmp(1)=arclen(mesh%edhx(ed1)%c%p, p)
    ntmp=1

    !Change radius in case of radius given too small
    rmax=r
    if(rmax<eps)then
       rmax=999._r8
    end if

    !Get nearest edges to point
    do i=2, maxn

       !Get next nearest node
       ed=nextnearesthxedge(p, mesh, i-1, listedtmp, d)
       !print*, i, node, d
       !Check if radius reached or maximum number of edges
       if(d>rmax .or. i==n .or. i==mesh%ne/2 .or. d>pio2)then
          n=i-1
          exit
       elseif(i==maxn)then
          n=i
          exit
       end if

       !Add new node to list
       listedtmp(i)=ed
       listdtmp(i)=d

    end do

    if((.not.allocated(listed)))then
       allocate(listed(1:n))
       allocate(listd(1:n))
    elseif(size(listed)<n)then
       print*, "NEAREDGES WARNING: Received allocated array with"
       print*, "   too little space. Returning up until ", size(listed)
       listd(1:size(listed))=listdtmp(1:size(listed))
       listed(1:size(listed))=listedtmp(1:size(listed))
       return
    end if

    listd(1:n)=listdtmp(1:n)
    listed(1:n)=listedtmp(1:n)

    return

  end subroutine getnearhxedges

  function nextnearesthxedge(p, mesh, n, listed, dist)
    !----------------------------------------------------------
    !	NEXTNEARESTHXEDGE
    !
    !	Given a point p returns 
    !		the index of the next nearest hx edge for a given list 'listv'
    !       of the 'n' already nearest nodes in order
    !
    !--------------------------------------------------------------	
    !Input variables
    real (r8), intent(in) :: p(1:3)

    !Mesh
    type(grid_structure), intent(in) :: mesh

    !Number of nodes already in list
    integer (i4), intent(in) :: n

    !Output array
    integer (i4), intent(in) :: listed(1:n)

    !Distance to point
    real (r8), intent(out) :: dist

    !Next nerest node
    integer (i4):: nextnearesthxedge

    !Auxiliar Variables
    integer (i4):: i
    integer (i4):: j
    integer (i4):: k
    integer (i4):: l
    integer (i4):: kt

    integer (i4):: ed
    integer (i4):: ednb
    integer (i4):: nearested
    real (r8):: dmin
    real (r8):: d


    !Check if list is reliable
    if(listed(n)<=0 .or. listed(n)>mesh%ne )then
       print*, "NEXTNEARESTHXEDGE ERROR: Listv incorect", listed(1:n)
       stop
    end if

    !  For all members of the list check for the nearest neighbour edge that
    !    is closer to the the given point
    dmin=9999._r8
    listmembers: do i = 1, n
       !Get node index
       ed=listed(i)


       !Check the other edges of the triangle for proximity
       trs: do k=1,2
          !Triangles that contain edge
          kt=mesh%edhx(ed)%v(k)
          treds: do j = 1, 3
             !Get edge index
             ednb=mesh%tr(kt)%ed(j)

             !Check if already on the list
             do l=1, n
                if(ednb==listed(l))then
                   cycle treds
                end if
             end do

             !Calculate distance to point
             d = arclen(mesh%edhx(ednb)%c%p, p)
             if(d < dmin)then
                !New nearest edge!
                nearested=ednb
                dmin=d
             end if
          end do treds
       end do trs
    end do listmembers

    if(dmin > 2._r8)then
       !No new node to include
       nextnearesthxedge=0
       dist=dmin
    else
       !Return new node
       nextnearesthxedge=nearested
       dist=dmin
    end if

    return

  end function nextnearesthxedge

  function getothernodeonedge(ed, node, mesh)
    !----------------------------------------------------------
    !	GETOTHERNODEONEDGE
    !
    !	Given an edge index, and a node index, returns the other node
    !    on the edge, diferent from 'node'
    !
    !--------------------------------------------------------------	
    !Edge index
    integer (i4), intent(in) :: ed
    !Node index
    integer (i4), intent(in) :: node
    !Mesh
    type(grid_structure), intent(in) :: mesh
    !Other node
    integer (i4):: getothernodeonedge

    if(mesh%ed(ed)%v(1)==node)then
       getothernodeonedge=mesh%ed(ed)%v(2)
    else
       getothernodeonedge=mesh%ed(ed)%v(1)
    end if

    return
  end function getothernodeonedge

  function getnearnode(p, mesh, tr)
    !----------------------------------------------------------
    !	GETNEARNODE
    !
    !	Given a point p returns 
    !		the nearest node in the mesh to this point
    !           using the regular to geodesic structure
    !   Uses lat/lon search table
    !--------------------------------------------------------------	
    !Input variables
    real (r8), intent(in) :: p(1:3)
    type(grid_structure), intent(in) :: mesh

    !Triangle already known
    integer (i4), optional :: tr

    !Output variable
    integer (i4):: getnearnode

    !Auxiliar Variables
    integer (i4):: kt
    integer (i4):: i
    integer (i4):: j
    integer (i4):: nearnode
    real (r8):: d
    real (r8):: dmin

    if(.not.present(tr))then
       kt=gettr(p, mesh)
    else
       kt=tr
    end if

    dmin=10000000._r8
    nearnode=0
    do i=1,3
       j=mesh%tr(kt)%v(i)
       d = norm(mesh%v(j)%p-p)
       if(d<dmin)then
          nearnode=j
          dmin=d
       end if
    end do

    getnearnode=nearnode
    return

  end function getnearnode

  function getnearhxedge(p, mesh)
    !----------------------------------------------------------
    !	GETNEAREDGE
    !
    !	Given a point p returns 
    !		the index of the hexagonal edges that has 
    !       midpoint closest to the point
    !   Uses lat/lon search table
    !--------------------------------------------------------------	
    !Input variables
    real (r8), intent(in) :: p(1:3)
    type(grid_structure), intent(in) :: mesh

    !Output variable
    integer (i4):: getnearhxedge

    !Auxiliar Variables
    integer (i4):: kt
    integer (i4):: i
    integer (i4):: j
    integer (i4):: nearedge
    integer (i4):: node
    real (r8):: d
    real (r8):: dcc
    real (r8):: dmin


    kt=gettr(p,mesh)

    dcc = arclen(mesh%tr(kt)%c%p, p) 
    if(dcc < mesh%mincdist/2._r8)then
       !The point is close to a triangle circumcenter
       !Check the triangles edges
       dmin=100000._r8
       do i=1,3
          j=mesh%tr(kt)%ed(i)
          d = arclen(mesh%edhx(j)%c%p, p)
          !Keep the nearest edge up until now
          if(d<dmin)then
             nearedge=j
             dmin=d
          end if
       end do

    else
       !The point is close to a node
       node=getnearnode(p, mesh)
       !Check the voronoi cell edges
       dmin=100000._r8
       do i=1,mesh%v(node)%nnb
          j=mesh%v(node)%ed(i)
          d = arclen(mesh%edhx(j)%c%p, p)
          !Keep the nearest edge up until now
          if(d<dmin)then
             nearedge=j
             dmin=d
          end if
       end do

    end if

    getnearhxedge=nearedge
    return

  end function getnearhxedge

  function getneartredge(p, mesh)
    !----------------------------------------------------------
    ! GETNEARTREDGE
    !
    ! Given a point p returns
    !   the index of the triangle edges that has
    !       midpoint closest to the point
    !   Uses lat/lon search table
    !--------------------------------------------------------------
    !Input variables
    real (r8), intent(in) :: p(1:3)
    type(grid_structure), intent(in) :: mesh

    !Output variable
    integer (i4):: getneartredge

    !Auxiliar Variables
    integer (i4):: kt
    integer (i4):: i
    integer (i4):: j
    integer (i4):: nearedge
    integer (i4):: node
    real (r8):: d
    real (r8):: dcc
    real (r8):: dmin


    kt=gettr(p,mesh)

    dcc = arclen(mesh%tr(kt)%c%p, p)
    if(dcc < mesh%mincdist/2._r8)then
       !The point is close to a triangle circumcenter
       !Check the triangles edges
       dmin=100000._r8
       do i=1,3
          j=mesh%tr(kt)%ed(i)
          d = arclen(mesh%ed(j)%c%p, p)
          !Keep the nearest edge up until now
          if(d<dmin)then
             nearedge=j
             dmin=d
          end if
       end do

    else
       !The point is close to a node
       node=getnearnode(p, mesh)
       !Check the voronoi cell edges
       dmin=100000._r8
       do i=1,mesh%v(node)%nnb
          j=mesh%v(node)%ed(i)
          d = arclen(mesh%ed(j)%c%p, p)
          !Keep the nearest edge up until now
          if(d<dmin)then
             nearedge=j
             dmin=d
          end if
       end do

    end if

    getneartredge=nearedge
    return

  end function getneartredge

  function getneartrc(p, mesh)
    !----------------------------------------------------------
    !	GETNEARTRC
    !
    !	Given a point p returns 
    !		the triangle index that has the circumcenter nearest 
    !       to the point
    !   Uses lat/lon search table
    !--------------------------------------------------------------	
    !Input variables
    real (r8), intent(in) :: p(1:3)
    type(grid_structure), intent(in) :: mesh

    !Output variable
    integer (i4):: getneartrc

    !Auxiliar Variables
    integer (i4):: i
    integer (i4):: kt
    integer (i4):: k
    integer (i4):: neartr
    real (r8):: d
    real (r8):: dmin

    kt=gettr(p, mesh)

    !Define initially the triangle that contains the point
    !  as the closest one
    dmin=arclen(p, mesh%tr(kt)%c%p)
    neartr=kt

    !Check if the the neighbour trs have cc
    ! closer to point, this happens when triangles
    ! are a bit distorced
    do i=1,3
       k=mesh%tr(kt)%nb(i)
       d = arclen(p, mesh%tr(k)%c%p)
       if(d<dmin)then
          neartr=k
          dmin=d
       end if
    end do

    getneartrc=neartr
    return

  end function getneartrc

  function insmallvoredtr(p, eds, mesh)
    !-------------------------------------------------------------
    !	INSMALLVOREDTR
    !
    !	Given a point in cartesian coords and 3 edge indexes
    !   tests if the point belongs to the triangle formed by
    !   the 3 voronoi edges midpoints.
    !   The eds indexes must be in counter clock wise order
    !--------------------------------------------------------------	
    !Point
    real (r8), intent(in) :: p(1:3)

    !Edges' indexes
    integer (i4), dimension(1:3), intent(in) :: eds

    !Mesh structure
    type(grid_structure), intent(in) :: mesh

    !Logical return value
    logical:: insmallvoredtr
    logical:: intr

    !Check if the point is inside the triangle
    intr=.false.
    intr=insidetr(p, mesh%edhx(eds(1))%c%p, &
         mesh%edhx(eds(2))%c%p, &
         mesh%edhx(eds(3))%c%p)

    !Return 
    insmallvoredtr=intr

    return
  end function insmallvoredtr


  function gettr(p, mesh )
    !-------------------------------------------------------------
    !	GETTR
    !
    !	Given a point p  returns 
    !		the triangle it belongs (kt) for a given mesh
    !           using the regular to geodesic structure
    !
    ! This routine depende on the "inside" funciton of Stripack
    !   which must be improved. For this reason several repated checks
    !   are done to ensure tha correct trinagles is obtained,
    !   altough this makes it more expensive.
    !
    ! You are encouraged to develop a fail proof point in triangle
    !   check routine to solve this out!
    !--------------------------------------------------------------	
    !Input variables
    real (r8), intent(in):: p(1:3)
    type(grid_structure), intent(in) :: mesh

    !Output variable
    integer(i4):: gettr

    !Auxiliar variables
    real (r8):: lon
    real (r8):: lat
    real (r8):: near
    real (r8):: neartmp
    real (r8):: latij
    real (r8):: lonij
    integer(i4):: kt
    integer(i4):: i
    integer(i4):: j
    integer(i4):: k
    integer(i4):: intersec
    logical:: isintr
    real (r8), dimension (1:3) :: q

    gettr=0

    !Get search tabel index
    call cart2sph ( p(1), p(2), p(3), lon, lat )
    call convllij( lon, lat, i, j, mesh)
    !print*, "Checking: ", lon*rad2deg, lat*rad2deg

    !Check is square table is good
    if (mesh%qd_lat(j)%qd_lon(i)%ntr==0) then
       print*, "GETTR Error: No triangle for this square. Check search table matrix"
       print*, "Lon:", lon, i
       print*, "Lat:", lat, j
       stop
    end if


    !Check the triangle list for this square
    do k=1, mesh%qd_lat(j)%qd_lon(i)%ntr
       gettr=mesh%qd_lat(j)%qd_lon(i)%tr(k)
       ! Check if point in triangle
       isintr=insidetrmesh(p, gettr, mesh)
       if(isintr) then !Found the triangle
          !print*, "Tr:", gettr, " Nodes: ", mesh%trnodes(gettr,:),&
          ! "Circumcen (lat, lon): ", mesh%clat(gettr)*rad2deg, mesh%clon(gettr)*rad2deg
          return
       end if
    end do

    !print*, "gettr warning: Could not fit point into any triangle. Doing further check..."
    !print "(a8,2f8.2)", "  Point:", lon*rad2deg, lat*rad2deg
    !print*, "  Triangles in sq:", mesh%qd_lat(j)%qd_lon(i)%tr(1:mesh%qd_lat(j)%qd_lon(i)%ntr)
    !print "(a10,2i8,2f8.2)", "  Square:", i, j, &
    !     -180._r8+mesh%qd_lat(j)%dlon*rad2deg*(real(i-1)), &
    !     -90._r8+mesh%dlat*rad2deg*(real(j-1))

    !stop
    !If the point could not be put into any triangle, the point may be at the boundaries of a 
    ! triangle or too near a node, so we do some extra checking

    !Now get the nearest triangle (relative to its barycenter) on the list to point
    near=100000.0
    do k=1, mesh%qd_lat(j)%qd_lon(i)%ntr
       kt=mesh%qd_lat(j)%qd_lon(i)%tr(k)
       neartmp=arclen(mesh%tr(kt)%b%p, p)
       if(neartmp<=near)then
          gettr=kt
          near=neartmp
       end if
    end do
    !print*, "  Nearest triangle: ", gettr

    !print*, lon*rad2deg, lat*rad2deg, mesh%tr(gettr)%b%lon*rad2deg, &
    !      mesh%tr(gettr)%b%lat*rad2deg

    !Check if the point is too near to the triangle node
    do k=1, 3
       kt=mesh%tr(gettr)%v(k)
       q=p-mesh%v(kt)%p
       if(norm(q)<eps/1000._r8)then
          !Return the nearest triangle
          return
       end if
    end do

    !Now check if point is inside one of the surrounding triangles
    ! this is the case when the point is at an edge
    do k=1, 3
       kt=mesh%tr(gettr)%nb(k)
       isintr=insidetrmesh(p, kt, mesh)
       if(isintr) then
          !Found a triangle that has the point
          gettr=kt
          !print*, "Found triangle as neighbour of nearest tr: ", gettr
          return
       end if
    end do
    !read*, k
    !If it still doesn't find the triangle, it will do a full check
    ! and will do some extra tests to make sure its is correct
    print*, "gettr warning: Could not fit point into any triangle. Doing further check..."
    print "(a9,2f8.2)", "   Point:", lon*rad2deg, lat*rad2deg
    print*, "  Triangles in sq:", mesh%qd_lat(j)%qd_lon(i)%tr(1:mesh%qd_lat(j)%qd_lon(i)%ntr)
    print "(a10,2i8,2f8.2)", "  Square:", i, j, &
         -180._r8+mesh%qd_lat(j)%dlon*rad2deg*(real(i-1)), &
         -90._r8+mesh%dlat*rad2deg*(real(j-1))
    print*, "  Best guess is nearest tr:", gettr

    !stop
    !If the mesh is randomly generated it might have too many
    !    issues, so due further testings, else return
    !if(.not. (trim(mesh%kind)=="rand"))then
    !   return
    !end if

    !print*," GETTR Warning, point not in any of the triangles of the regular grid list!!"
    !print*, "lon, lat", lon*rad2deg, lat*rad2deg, "trs:", mesh%qd_lat(j)%qd_lon(i)%tr(1:mesh%qd_lat(j)%qd_lon(i)%ntr)
    !print*

    !In the awful case where everything fails, to ensure that the function will
    !return the correct triangle, go through all triangles and test them
    !This may happen if you are using a very distorted mesh (random points)
    do k=1, mesh%nt
       isintr=insidetrmesh(p, k, mesh)
       if(isintr) then
          !print*, "Tr:", gettr, " Nodes: ", mesh%trnodes(gettr,:),&
          ! "Circumcen (lat, lon): ", mesh%clat(gettr)*rad2deg, mesh%clon(gettr)*rad2deg
          gettr=k
          print*, "  Full check returned tr: ", k
          print*
          return
       end if
    end do

    !In order to verify if the list of regular to triangular list is correct
    ! do a simple check if the square intersects the obtained triangle
    call convijll(i, j, lonij, latij, mesh)
    intersec=sqtriintersec(i, j, gettr, mesh)
    if(intersec==0.and.trim(mesh%kind)/="rand")then
       print '(a,2f32.16,a,2f14.4)',"GETTR Warning: Point",  lon*rad2deg, &
            lat*rad2deg," in sqr ", lonij*rad2deg, latij*rad2deg
       print '(a,i8, a, 2f14.4,a,i4)', " in triangle ", gettr, " with CC ", &
            mesh%tr(gettr)%c%lon*rad2deg, mesh%tr(gettr)%c%lat*rad2deg, &
            " but Sq does not intersec Tri ", intersec
       print*
    end if

    !Nothing worked
    print*, "GETTR ERROR: Could not find triangle for this point: ", p
    stop
    return
  end function gettr

  function getedge(v1, v2, mesh)
    !----------------------------------------------------------------------
    ! GETEDGE
    !   Given 2 triangle vertice indexes and a mesh structure
    !   return the index of the edge that has these vertices as endpoints
    !   To be used only after mesh is fully structured
    !----------------------------------------------------------------------
    integer (i4), intent(in) :: v1
    integer (i4), intent(in) :: v2
    type(grid_structure), intent(in) :: mesh

    integer (i4):: getedge

    integer (i4):: i
    integer (i4):: nb

    getedge=0
    do i=1, mesh%v(v1)%nnb !For all v1 neighbours
       nb=mesh%v(v1)%nb(i)
       !print*, v1, i, nb, v2
       if(nb==v2)then
          getedge=mesh%v(v1)%ed(i)
          return
       end if
    end do
    !The vertices are not neighbours
    return
  end function getedge

  function getedindexonhx(ed, cell, mesh)
    !----------------------------------------------------------------------
    ! getedindexonhx
    !   Given an edge number (global index) and a cell number (global)
    !   returns the local index of the edge with respect to cell
    !   returns 0 if edge not in cell
    !   To be used only after mesh is fully structured
    !----------------------------------------------------------------------
    integer (i4), intent(in) :: ed
    integer (i4), intent(in) :: cell
    type(grid_structure), intent(in) :: mesh

    integer (i4):: getedindexonhx

    integer (i4):: i

    getedindexonhx=0 !index of cell with respect to cell
    do i=1, mesh%v(cell)%nnb
       !Find the index of edge ed within cell(i)
       if(ed==mesh%v(cell)%ed(i))then
          getedindexonhx=i
          exit
       end if
    end do

    return
  end function getedindexonhx

  !========================================================================
  !    GEOMETRIC TOOLS FOR THE SPHERE
  !========================================================================

  subroutine sph2cart (lon, lat, x, y, z )
    !------------------------------------------------------------------------------------
    ! SPH2CART 
    !
    !     Transforms geographical coordinates (lat,lon) to Cartesian coordinates.
    !     Similar to stripack's TRANS
    !
    !    Input: LAT, latitudes of the node in radians [-pi/2,pi/2]
    !           LON, longitudes of the nodes in radians [-pi,pi]
    !
    !    Output:  X, Y, Z, the coordinates in the range -1 to 1. 
    !                    X**2 + Y**2 + Z**2 = 1 
    !---------------------------------------------------------------------
    real (r8), intent(in) :: lon
    real (r8), intent(in) :: lat
    real (r8), intent(out) :: x
    real (r8), intent(out) :: y
    real (r8), intent(out) :: z
    real (r8):: coslat

    coslat = dcos (lat)
    x = coslat * dcos (lon)
    y = coslat * dsin (lon)
    z = dsin (lat)

    return
  end subroutine sph2cart

  subroutine cart2sph ( x, y, z, lon, lat )
    !---------------------------------------------------------------------
    ! CART2SPH 
    !     Transforms  Cartesian coordinates to geographical (lat,lon) coordinates .
    !     Similar to stripack's SCOORD
    !
    !    Input:  X, Y, Z, the coordinates in the range -1 to 1. 
    !
    !    Output, LON, longitude of the node in radians [-pi,pi].
    !                       LON=0 if point lies on the Z-axis.  
    !    Output, LAT, latitude of the node in radians [-pi/2,pi/2].
    !                       LAT=0 if   X**2 + Y**2 + Z**2 = 0.
    !------------------------------------------------------------------------------------
    real    (r8), intent(in) :: x
    real    (r8), intent(in) :: y
    real    (r8), intent(in) :: z
    real    (r8), intent(out) :: lat
    real    (r8), intent(out) :: lon
    real    (r8):: pnrm

    pnrm = dsqrt ( x **2 + y**2 + z**2 )
    if ( x /= 0.0D+00 .or. y /= 0.0D+00 ) then
       lon = datan2 ( y, x )
    else
       lon = 0.0D+00
    end if

    if ( pnrm /= 0.0D+00 ) then
       lat = dasin ( z / pnrm )
    else
       print*, "CART2SPH Warning: Point not in the unit sphere. Norm= ", pnrm
       lat = 0.0D+00
    end if

    return
  end subroutine cart2sph

  subroutine convert_vec_cart2sph(p, v , vlon, vlat)
    !---------------------------------------------------------------------
    !	CONVERT_VEC_CART2SPH
    !
    !   Recieves a point p=(x,y,z) and a vector (v) at this point
    !   Returns the vector at this point in (lon,lat) reference
    !      vlon=West-East direction
    !      vlat=South-North direction
    !
    !   The vector must be tangent to the sphere, that is, v.(x,y,z)=0
    !   This is done in order to plot the vectorfield on GMT
    !---------------------------------------------------------------------
    !Point cartesian coords
    real(r8), intent(in) :: p(1:3)
    !Cartesian vector on point
    real(r8), intent(in) :: v(1:3)
    !Spherical coord vector on point
    real(r8), intent(out) :: vlat
    real(r8), intent(out) :: vlon
    !Auxiliar variables
    real(r8):: r
    real(r8):: rho
    real(r8):: rvec(1:3)
    real(r8):: latvec(1:3)
    real(r8):: lonvec(1:3)
    real(r8):: zero
    real(r8):: test

    zero=0
    r=dsqrt(p(1)**2+p(2)**2+p(3)**2)
    rho=dsqrt(p(1)**2+p(2)**2)

    !Case where the point is in the north or south pole
    if(rho<10*eps)then
       !print*, "Pole:", v
       vlon=v(2)
       vlat=v(1)
       return
    end if

    rvec=(/p(1),p(2),p(3)/)
    rvec=rvec/r

    latvec=(/-p(1)*p(3),-p(2)*p(3),rho**2/)
    latvec=latvec/rho

    lonvec=(/-p(2), p(1), zero /)
    lonvec=lonvec/rho

    test=dot_product(v,rvec)
    if(abs(test)>10e-5)then
       print *,"CONVERT_VEC_CART2SPH Warning: Vector not tangent to sphere."
       print '(a,3f10.6)',"Vector:",v(1:3)
       print '(a,3f10.6)',"Point:",rvec(1:3)
       print '(a,f10.6)',"Dot Product:", test
       stop
    end if
    vlat=dot_product(v,latvec)
    vlon=dot_product(v,lonvec)

    return    
  end subroutine convert_vec_cart2sph

  subroutine convert_vec_sph2cart(vlon, vlat, p, v)
    !---------------------------------------------------------------------
    !	CONVERT_VEC_SPH2CART
    !
    !   Recieves a point p=(x,y,z) and a vector at this point in lon, lat system
    !   in radians (vlon, vlat), ie, 
    !      vlon=West-East direction
    !      vlat=South-North direction
    !   Returns the vector at this point in cartesian reference (v)
    !---------------------------------------------------------------------
    !Point cartesian coords
    real(r8), intent(in) :: p(1:3)
    real(r8), intent(in) :: vlon
    real(r8), intent(in) :: vlat
    !Cartesian vector on point
    real(r8), intent(out) :: v(1:3)
    !Auxiliar variables
    real(r8):: r
    real(r8):: rho

    r=dsqrt(p(1)**2+p(2)**2+p(3)**2)
    rho=dsqrt(p(1)**2+p(2)**2)

    !Case where the point is in the north or south pole
    if(rho==0)then
       v(1)=vlat
       v(2)=vlon
       v(3)=0
       return
    else    
       !The minus sign in vlat is due to the diference between spherical coords and
       !   geographical coords
       v(1)=-vlon*(p(2)/rho) - (vlat*p(1)*p(3))/(rho)
       v(2)=vlon*(p(1)/rho) - (vlat*p(2)*p(3))/(rho)
       v(3)=vlat*rho
    end if

    return    
  end subroutine convert_vec_sph2cart

  function proj_vec_sphere(v, p)
    !-----------------------------------------------------------
    !  Projects a vector 'v' on the plane tangent to a sphere
    !   Uses the the tangent plane relative to the unit sphere's
    !   point 'p', in cartesian coords
    !-----------------------------------------------------------
    real (r8), intent(in), dimension(1:3) :: v
    real (r8), intent(in), dimension(1:3) :: p
    real (r8), dimension(1:3)  :: proj_vec_sphere

    proj_vec_sphere(1:3)=&
         v(1:3)-dot_product(v,p)*p(1:3)/norm(p)

    return
  end function proj_vec_sphere

  function nlon_redgrd(j, nlat)
    !-----------------------------------------------------------
    !  nlon_redgrd - Number of longitudes in reduced grid
    !
    !  Given a number of latitudes (nlat) and a latitude index (j)
    !  varying from south pole (1) to near northo pole (nlat+1)
    !  return the number of longitudes for this latitude
    !-----------------------------------------------------------
    integer (i4), intent(in) :: j
    integer (i4), intent(in) :: nlat
    integer (i4) :: nlon_redgrd
    real (r8) :: lat

    lat=-pio2+(real(j,r8)-0.5)*pi/(real(nlat,r8))
    nlon_redgrd=ceiling(2*nlat*dcos(lat))
    !nlon_redgrd=2*nlat
    return
  end function nlon_redgrd


  subroutine convijll(i, j, lon, lat, mesh)
    !----------------------------------------------------------
    !	CONVIJLL
    !
    !	Converts the (i,j) indexation to (lat,lon) in radians
    !      according to a given mesh specifications
    !      Lat [-pi/2, pi/2] , Lon [-pi, pi[   
    !--------------------------------------------------------------	
    !Input variables
    integer(i4), intent(in) :: i
    integer(i4), intent(in) ::j
    type(grid_structure), intent(in) :: mesh

    !Output variables
    real (r8), intent(out) :: lat
    real (r8), intent(out) :: lon

    lat=-90.*deg2rad+(j-1)*mesh%dlat

    lon=-180.*deg2rad+(i-1)*mesh%qd_lat(j)%dlon

    if(abs(lat)>90.*deg2rad.or.abs(lon)>180.*deg2rad)then
       print*, "Convijll Warning:(Lat, Lon) out of range !"
    end if

    return

  end subroutine convijll

  subroutine convllij(lon, lat, i, j, mesh)
    !----------------------------------------------------------
    !	CONVLLIJ
    !
    !	Converts the (lon,lat) in radians indexation 
    !           to (i,j) according to
    !		a given mesh specifications
    !           Lat [-pi/2, pi/2] , lon [-pi, pi[   
    !--------------------------------------------------------------	
    !Input variables
    real (r8), intent(in) :: lat
    real (r8), intent(in) :: lon
    type(grid_structure), intent(in) :: mesh

    !Output variables
    integer(i4), intent(out) :: i
    integer(i4), intent(out) ::j

    !Aux variable
    real (r8):: tlat
    real (r8):: tlon

    tlat=lat
    tlon=lon

    !Check bounds
    if(abs(tlat)>90*deg2rad+eps .or. abs(tlon) >180.*deg2rad + eps )then
       print*, "Convllij Error:(Lon, Lat) out of range !"
       print*, tlon, tlat
       stop
    end if

    !Calculate latitude index (j)
    if(tlat>90.*deg2rad-eps) then
       j=mesh%nlat
    else
       j=int((tlat+90.*deg2rad)/mesh%dlat+eps)+1
    end if

    if(tlon>180.*deg2rad-eps) then
       i=1
    end if

    if(j>mesh%nlat) then
       !North pole
       j=int((tlat+90.*deg2rad)/mesh%dlat+eps)
       !print*, "Convllij Warning:(Lon, Lat) conversion might have error!"
       !print*, "Lon, lat ", tlon, tlat
       !print*, "i  , j ",  i, j
       !print*, mesh%dlat, (tlat+90.*deg2rad)/mesh%dlat
       !print*, mesh%qd_lat(j)%nlon, mesh%qd_lat(j)%dlon
    elseif(j<0)then
       print*, "Convllij Error:(Lon, Lat) out of range"
       stop
    else
       !Longitude index (i)
       i=int((tlon+180.*deg2rad)/mesh%qd_lat(j)%dlon+eps)+1
    endif

    if(i>mesh%qd_lat(j)%nlon) then
       i=mesh%qd_lat(j)%nlon
       print*, "Convllij Warning:(Lon, Lat) conversion might have error!"
       print*, "Lon, lat ", tlon, tlat
       print*, "i  , j ",  i, j
       print*, mesh%qd_lat(j)%dlon, (tlon+180.*deg2rad)/mesh%qd_lat(j)%dlon
    end if

    return
  end subroutine convllij

  function arcdistll(lon1, lat1, lon2 , lat2)
    !----------------------------------------------------
    ! ARCDISTLL
    !   Arc length between p1 and p2 points on the sphere
    !     Receives latitude longitude coordinates in radians 
    !	  Returns the angle between points in radians
    !     Lat [-pi/2, pi/2] , lon [-pi, pi[   
    !
    !	Uses Vincenty formula
    !-----------------------------------------------------
    real (r8), intent(in) :: lat1
    real (r8), intent(in) :: lon1
    real (r8), intent(in) :: lat2
    real (r8), intent(in) :: lon2
    real (r8):: arcdistll
    real (r8):: dlat
    real (r8):: dlon

    dlat=lat1-lat2
    dlon=lon1-lon2

    arcdistll=0

    arcdistll= datan2(dsqrt( (dcos(lat2)*dsin(dlon))**2 + &
         (dcos(lat1)*dsin(lat2)-dsin(lat1)*dcos(lat2)*dcos(dlon))**2 ),&
         (dsin(lat1)*dsin(lat2)+dcos(lat1)*dcos(lat2)*dcos(dlon)) )
    !See wikipedia - Vincenty formula for the sphere
    arcdistll=abs(arcdistll)
    return

  end function arcdistll

  function arcdistxyz(p1, p2)
    !----------------------------------------------------
    ! ARCDISTXYZ - Prefer to use ARCLEN
    !    Arc length between p1 and p2 points on the sphere
    !     Receives cartesian coords. of points
    !     Returns the angle between points in radians
    !     It if preferable to use the routine ARCLEN
    !      to avoid non standard inverse cossine 
    !-----------------------------------------------------
    real (r8), dimension(3), intent(in) :: p1
    real (r8), dimension(3), intent(in) :: p2
    real (r8):: arcdistxyz
    real (r8):: p1norm
    real (r8):: p2norm

    arcdistxyz=0
    p1norm=dsqrt(dot_product(p1,p1))
    p2norm=dsqrt(dot_product(p2,p2))

    if(p1norm==0.or.p2norm==0) then
       print*
       print*,"ARCDIST Warning: Distance calculation error"
       print*,"p1:", p1, " p2:", p2
       return
    end if
    arcdistxyz=dacos(dot_product(p1,p2)/(p1norm*p2norm))

    if((p1norm >= p2norm+100*eps).or. &
         (p1norm <= p2norm-100*eps)) then
       print*, "ARCDIST WARNING: Vectors with different norms (dif):",p1norm-p2norm
       print*, "ARCDIST WARNING: Vectors are not on the same sphere! "
       print*
       return
    end if

    return
  end function arcdistxyz

  function arclen(p, q)
    !-----------------------------------------------------------
    ! ARCLEN from ssrfpack by r. renka            
    !
    !   This function computes the arc-length (angle in radians)            
    !    between a pair of points on the unit sphere. It is similar to
    !    arcdistxyz, but uses a calculation to avoid inverse cosine function
    !    p,q = arrays of length 3 containing the x, y, and z             
    !             coordinates (in that order) of points on the              
    !             unit sphere.                                                 
    !   Returns the angle in radians between the unit vectors              
    !       p and q.  0 <= arclen <= pi.                       
    !---------------------------------------------------------
    real (r8), intent(in) :: p(3)
    real (r8), intent(in) :: q(3)
    real (r8):: arclen

    !Dot product
    real (r8):: d

    d = dot_product(p+q, p+q)

    if (d==0.) then 
       ! p and q are separated by 180 degrees.
       arclen = pi
    elseif (d>=4.) then 
       ! p and q coincide.
       arclen = 0._r8
    else 
       arclen = 2._r8 * datan (dsqrt ( (4._r8 - d) / d) )
    endif

    return 
  end function arclen

  function trcircumc(p1, p2, p3)
    !------------------------------------------
    !  SPHERICAL TRIANGLE CIRCUMCENTER
    ! Calculates the circumcentre of the triangle
    !  given by points p1, p2, p3
    !--------------------------------------------

    !Cartesian coordinates of the nodes of spherical triangle
    real (r8),intent(in) :: p1(1:3)
    real (r8),intent(in) :: p2(1:3)
    real (r8),intent(in) :: p3(1:3)

    !Cross product
    real (r8):: cp(1:3)

    !Spherical triangle area
    real (r8):: trcircumc(1:3)

    !Check for degeneration
    if(norm(p1-p2)<eps/100 .or. norm(p2-p3)<eps/100 .or. norm(p1-p3)<eps/100)then
       trcircumc=0._r8
       return
    end if

    cp=cross_product(p2-p1, p3-p1)
    trcircumc=cp/norm(cp)

    return
  end function trcircumc

  function triarea(p1, p2, p3)
    !------------------------------------------
    !  TRIAREA
    ! Calculates the area of
    !    a planar triangle on the unit sphere
    !    given the cartesian coords of nodes
    !--------------------------------------------
    !Cartesian coordinates of the nodes of spherical triangle
    real (r8),intent(in) :: p1(1:3)
    real (r8),intent(in) :: p2(1:3)
    real (r8),intent(in) :: p3(1:3)

    !Area
    real (r8):: triarea

    triarea= norm(cross_product( p2-p1, p3-p1))/2

    return
  end function triarea


  function sphtriarea(p1, p2, p3)
    !------------------------------------------
    !  SPHTRIAREA
    ! Calculates the area of 
    !    a spherical triangle on the unit sphere
    !    given the cartesian coords of nodes
    !--------------------------------------------

    !Cartesian coordinates of the nodes of spherical triangle
    real (r8),intent(in) :: p1(1:3)
    real (r8),intent(in) :: p2(1:3)
    real (r8),intent(in) :: p3(1:3)

    !Variables for spherical triangle area calculation
    ! a, b, c are lengths of the the 3 sides
    ! s is the semi-perimiter s=(a+b+c)/2
    ! e is the spherical excess of triangle
    !   e = a + b + c - 180, but we will useing L'Huilier's Theorem
    real (r8):: a
    real (r8):: b
    real (r8):: c
    real (r8):: s
    real (r8):: e
    real (r8):: tmp

    !Spherical triangle area
    real (r8):: sphtriarea

    !Check for degeneration
    if(norm(p1-p2)<eps/100 .or. norm(p2-p3)<eps/100 .or. norm(p1-p3)<eps/100)then
       sphtriarea=0._r8
       return
    end if

    !Calculate the sides length's and the semiperimiter
    s=0.
    a=arclen(p1,p2) !arcdistll(lon(1),lat(1), lon(2),lat(2))
    b=arclen(p2,p3)  !arcdistll(lon(2),lat(2), lon(3),lat(3))
    c=arclen(p3,p1)  !arcdistll(lon(3),lat(3), lon(1),lat(1))
    s=(a+b+c)/2

    !Calculate spherical triangle excess using L'Huilier's Theorem
    tmp=dtan(s/2)*dtan((s-a)/2)*dtan((s-b)/2)*dtan((s-c)/2)
    !Round off error might give almost zero negative numbers => assume zero
    if(tmp<0)then
       e=0.
    else
       e = 4*datan(dsqrt(tmp))
    end if
    sphtriarea=e

    return
  end function sphtriarea


  function sphhxarea(hx, mesh)
    !-------------------------------------------------
    !  SPHHXAREA
    ! Calculates the area of 
    !    a spherical voronoi/hexagonal/pentagonal cell
    !    on the unit sphere
    !-------------------------------------------------
    integer(i4), intent(in) :: hx
    type(grid_structure), intent(in) :: mesh

    real (r8):: sphhxarea

    integer(i4):: i
    integer(i4):: i1
    integer(i4):: i2
    integer(i4):: nb
    real (r8), dimension(1:3) :: p1
    real (r8), dimension(1:3) :: p2
    real (r8), dimension(1:3) :: p3

    nb=mesh%v(hx)%nnb

    p1=mesh%v(hx)%p

    sphhxarea=0
    do i=1,nb !For the coresponding neighbour edges/triangles
       !Get triangle center points for this edge
       if(i==1) then
          i1=mesh%v(hx)%tr(nb)
          i2=mesh%v(hx)%tr(1)
       else
          i1=mesh%v(hx)%tr(i-1)
          i2=mesh%v(hx)%tr(i)
       end if
       p2=mesh%tr(i1)%c%p
       p3=mesh%tr(i2)%c%p
       sphhxarea=sphhxarea+sphtriarea(p1, p2, p3)
    end do

    return
  end function sphhxarea

  function sphpolarea(p, n)
    !-------------------------------------------------
    !  SPHPOLAREA
    ! Calculates the area of 
    !    a spherical polygon formed by (p(1,:), p(2,:), ..., p(n,:))
    !    points ordered counterclockwisely
    ! They may degenerated (have equal points)
    !-------------------------------------------------

    !Polygon points (cartesian coords on the sphere)
    type(vector), dimension(1:n), intent(in) :: p
    type(vector), dimension(1:n) :: q

    !Geodesical area
    real (r8):: sphpolarea

    !Number of vertices
    integer (i4), intent(in) :: n
    integer (i4):: m

    !Euclidian Distance
    real (r8):: d

    !Indexes
    integer (i4):: i
    integer (i4):: j

    !Check for degenerancy
    m=n
    q(1)%v=p(1)%v
    j=2
    do i=2, n
       d=norm(p(i-1)%v-p(i)%v)
       !print*, i, j, d
       if(d<eps)then
          m=m-1
       else
          q(j)%v=p(i)%v
          j=j+1
       end if
    end do

    sphpolarea=0
    do i=2,m-1
       sphpolarea=sphpolarea+sphtriarea(q(1)%v, q(i)%v, q(i+1)%v)
       !print*, sphtriarea(q(1)%v, q(i)%v, q(i+1)%v)
    end do

    return
  end function sphpolarea

  function planarhexagarea(hx, mesh)
    !-------------------------------------------------
    ! PLANARHEXAGAREA
    ! Calculates the area of 
    !    a planar voronoi/hexagonal/pentagonal cell
    !    based on the underlying triangles
    !-------------------------------------------------
    integer(i4), intent(in) :: hx
    type(grid_structure), intent(in) :: mesh

    real (r8):: planarhexagarea

    integer(i4):: i
    integer(i4):: i1
    integer(i4):: i2
    integer(i4):: nb
    real (r8), dimension(1:3) :: p1
    real (r8), dimension(1:3) :: p2
    real (r8), dimension(1:3) :: p3
    real (r8):: area

    nb=mesh%v(hx)%nnb

    p1=mesh%v(hx)%p

    area=0
    do i=1,nb !For the coresponding neighbour edges/triangles
       !Get triangle center points for this edge
       if(i==1) then
          i1=mesh%v(hx)%tr(nb)
          i2=mesh%v(hx)%tr(1)
       else
          i1=mesh%v(hx)%tr(i-1)
          i2=mesh%v(hx)%tr(i)
       end if
       p2=mesh%tr(i1)%c%p
       p3=mesh%tr(i2)%c%p

       area = area+ abs(det(p1, p2, p3)/2)

    end do
    planarhexagarea=area

    return
  end function planarhexagarea

  function planarpjhexagarea(hx, p, mesh)
    !-------------------------------------------------
    ! PLANARPJHEXAGAREA
    ! Calculates the area of
    !    a planar voronoi/hexagonal/pentagonal cell
    !    which is the polygon radially projected
    !    to the tangent plane of "p" in the polygon
    !-------------------------------------------------
    integer(i4), intent(in) :: hx
    real (r8), intent(in):: p(1:3)
    type(grid_structure), intent(in) :: mesh

    real (r8):: planarpjhexagarea

    integer(i4):: i
    integer(i4):: i1
    integer(i4):: i2
    integer(i4):: nb
    real (r8), dimension(1:3) :: p1
    real (r8), dimension(1:3) :: p2
    real (r8), dimension(1:3) :: p3
    real (r8):: area

    nb=mesh%v(hx)%nnb

    p1=p !mesh%v(hx)%p

    area=0
    do i=1,nb !For the coresponding neighbour edges/triangles
       !Get triangle center points for this edge
       if(i==1) then
          i1=mesh%v(hx)%tr(nb)
          i2=mesh%v(hx)%tr(1)
       else
          i1=mesh%v(hx)%tr(i-1)
          i2=mesh%v(hx)%tr(i)
       end if
       p2=mesh%tr(i1)%c%p
       p2=p2/dot_product(p2,p1)
       p3=mesh%tr(i2)%c%p
       p3=p3/dot_product(p3,p1)

       area = area+ abs(det(p1, p2, p3)/2)

    end do
    planarpjhexagarea=area

    return
  end function planarpjhexagarea

  function sphtriangles(p1, p2, p3)
    !------------------------------------------
    !  SPHTRIANGLES
    ! Calculates the angles of 
    !    a spherical triangle on the unit sphere
    !--------------------------------------------

    !Latitudes and longitudes of nodes on the sphere
    real (r8),intent(in) :: p1(1:3)
    real (r8),intent(in) :: p2(1:3)
    real (r8),intent(in) :: p3(1:3)

    !Variables for spherical triangle area calculation
    ! a, b, c are lengths of the the 3 sides
    real (r8):: a
    real (r8):: b
    real (r8):: c

    !Spherical triangle area
    real (r8):: sphtriangles(1:3)

    !Calculate the sides length's and the semiperimiter
    a=arclen(p1,p2) !arcdistll(lon(1),lat(1), lon(2),lat(2))
    b=arclen(p2,p3)  !arcdistll(lon(2),lat(2), lon(3),lat(3))
    c=arclen(p3,p1)  !arcdistll(lon(3),lat(3), lon(1),lat(1))
    if((pi/2-a)<eps.or.(pi/2-b)<eps.or.(pi/2-c)<eps)then
       print*, "SPHTRIANGLES WARNING: Triangle with internal angle too large."
       print*, "Angles:", a, b, c
    end if
    !Using spherical law of cosine 
    sphtriangles(1)=dacos( (dcos(b)-dcos(a)*dcos(c)) / (dsin(a)*dsin(c)) )
    sphtriangles(2)=dacos( (dcos(c)-dcos(a)*dcos(b)) / (dsin(a)*dsin(b)) )
    sphtriangles(3)=dacos( (dcos(a)-dcos(b)*dcos(c)) / (dsin(b)*dsin(c)) )

    return
  end function sphtriangles

  function hxtr_intersec_areas(cell, n, mesh)
    !---------------------------------------------------------
    !Dual-primal cell area for "cell"
    ! Area of triangle intersection of Voronoi/hexagonal cell
    ! Receives index of cell and it's number of vertices "n"
    ! Returns areas in mesh
    !-----------------------------------------------------------
    !Voronoi cell index
    integer(i4), intent(in):: cell

    !Mesh structure
    type(grid_structure), intent(in) :: mesh

    !Calculated areas
    real(r8), dimension(1:n)::hxtr_intersec_areas

    !Number of triangles intersecting cell
    integer(i4), intent(in)::n

    !Index cell for vertices
    integer(i4):: j

    !Vertices of polygons
    !Cell node (central point)
    real(r8), dimension(1:3) :: pc
    !Edge point 1
    real(r8), dimension(1:3) :: ped1
    !Cell vertex point
    real(r8), dimension(1:3) :: pv
    !Edge point 2
    real(r8), dimension(1:3) :: ped2

    !Areas
    real(r8) :: area
    real(r8) :: area1
    real(r8) :: area2

    !Indexes for points
    integer(i4) :: ed1
    integer(i4) :: ed2
    integer(i4) :: tr

    !print*, "Cell:", cell
    !n=mesh%v(cell)%nnb
    !if(.not. allocated(mesh%hx(cell)%hxtr_area))then
    !  allocate(mesh%hx(cell)%hxtr_area(1:n))
    !end if
    !mesh%hx(cell)%hxtr_area(1:n)=0.

    !Central point
    pc=mesh%v(cell)%p

    !print*, "Cell: ", cell
    !For cell vertex
    do j=1, n
       !Edges around vertex j
       ed1=mesh%v(cell)%ed(j)
       ed2=mesh%v(cell)%ed(modint(j+1,n))
       ped1=mesh%edhx(ed1)%c%p
       ped2=mesh%edhx(ed2)%c%p
       !Triangle centered at vertex j
       tr=mesh%v(cell)%tr(j)
       !Vertex j - circumcenter of tr
       pv=mesh%tr(tr)%c%p
       !Calculate first area
       area1=sphtriarea(pc, ped1, pv)
       !Correct sign, in case circuncenter out of tr
       if(det(pc, ped1, pv)<0)then
          area1=-area1
       end if
       !Calculate 2nd area
       area2=sphtriarea(pc, pv, ped2)
       !Correct sign, in case circuncenter out of tr
       if(det(pc, pv, ped2)<0)then
          area2=-area2
       end if
       area=area1+area2
       !print*, j, area1 + area2
       hxtr_intersec_areas(j)=area
       !print*, "Area of vertex ", j, area, mesh%hx(cell)%areag
       !print*, mesh%hx(cell)%hxtr_area(j)
    end do
    !print*, "Sum Areas and error:", sum(hxtr_intersec_areas(1:n)), &
    !    sum(hxtr_intersec_areas(1:n))-mesh%hx(cell)%areag
    !read(*,*) j

    return
  end function hxtr_intersec_areas

  subroutine trhx_intersec_areas(mesh)
    !---------------------------------------------------------
    !Dual-primal cell area for all trinalge in mesh
    ! Area of triangle intersection with Voronoi/hexagonal
    ! Receives the mesh and fills in mesh%tr%trhx_area
    !-----------------------------------------------------------

    !Mesh structure
    type(grid_structure), intent(inout) :: mesh

    !Vertices of polygons
    !Cell node (central point)
    real(r8), dimension(1:3) :: pc
    !Edge point 1
    real(r8), dimension(1:3) :: ped1
    !Cell vertex point
    real(r8), dimension(1:3) :: pv
    !Edge point 2
    real(r8), dimension(1:3) :: ped2

    !Areas
    real(r8) :: area1
    real(r8) :: area2

    !Indexes for points
    integer(i4) :: ed1
    integer(i4) :: ed2
    integer(i4) :: j
    integer(i4) :: k


    !For all triangles
    do k=1, mesh%nt

       !Central point
       pc=mesh%tr(k)%c%p
       !print*
       !print*, "Triangle:", k
       !For tr vertex
       do j=1, 3
          !Edges around vertex j
          ed1=mesh%tr(k)%ed(modint(j-1,3))
          ed2=mesh%tr(k)%ed(j)
          ped1=mesh%ed(ed1)%c%p
          ped2=mesh%ed(ed2)%c%p
          ! print*, ed1, ed2, mesh%tr(k)%v(j)
          !Vertex j - hx node
          pv=mesh%v(mesh%tr(k)%v(j))%p

          !Calculate first area
          area1=sphtriarea(pc, ped1, pv)
          !Correct sign, in case circuncenter out of tr
          if(det(pc, ped1, pv)<0)then
             area1=-area1
          end if
          !Calculate 2nd area
          area2=sphtriarea(pc, pv, ped2)
          !Correct sign, in case circuncenter out of tr
          if(det(pc, pv, ped2)<0)then
             area2=-area2
          end if
          !print*, j, area1, area2
          mesh%tr(k)%trhx_area(j)=area1+area2
          !print*, "Area of vertex ", j, area, mesh%hx(cell)%areag
          !print*, mesh%hx(cell)%hxtr_area(j)
       end do

       !print*, "Sum Areas and error:", sum(mesh%tr(k)%trhx_area(1:3)), &
       !    sum(mesh%tr(k)%trhx_area(1:3))-mesh%tr(k)%areag
       !    if(k==153)then
       !read(*,*) l
       !end if
    end do
    return
  end subroutine trhx_intersec_areas


  subroutine calc_tiled_areas(mesh)
    !-------------------------------------------------
    !  SPHHXAREA
    ! Calculates the area of
    !    a spherical voronoi/hexagonal/pentagonal cell
    !    on the unit sphere
    !-------------------------------------------------
    type(grid_structure), intent(inout) :: mesh

    real (r8):: d

    integer(i4):: i, j, ed
    integer(i4):: nb

    !for all edges calculate tiled edge volume area
    do i=1,mesh%ne
      mesh%ed(i)%areat=mesh%ed(i)%leng*mesh%edhx(i)%leng
      mesh%edhx(i)%areat=mesh%ed(i)%areat
    end do

    !for all triangle calculate tiled triangle area
    do i=1,mesh%nt
      mesh%tr(i)%areat=0.0
      do j=1,3
        ed=mesh%tr(i)%ed(j)
        d=arclen(mesh%tr(i)%c%p, mesh%ed(ed)%c%p)
        mesh%tr(i)%areat=mesh%tr(i)%areat + d*mesh%ed(ed)%leng/2.0
      end do
    end do

    !for all hx/voronoi cells calculate tiled area
    do i=1,mesh%nv
      mesh%hx(i)%areat=0.0
      do j=1,mesh%v(i)%nnb
        ed=mesh%v(i)%ed(j)
        mesh%hx(i)%areat=mesh%hx(i)%areat + mesh%ed(ed)%areat/4.0
      end do
    end do

    return
  end subroutine calc_tiled_areas


  function bar_coord_tr(p, tr, mesh)
    !----------------------------------------------------------
    !	BARYCENTRIC COORDINATES for a Triangle in mesh
    !
    !	Given a point in cartesian coords and a triangle that 
    !   contains this point, returns the barycentric coords
    !   of the point. (p must be in k)
    !   The barycentric coordinates refer to the planar triangle
    !   and because the point inside the triangle is in fact
    !   in the spherical triangle, we use linearity formula
    !   and the projected value to calculate the coords.
    !--------------------------------------------------------------	    
    !Point inside triangle tr
    real (r8), intent(in) :: p(1:3)

    !Triangle containing the point p
    integer (i4), intent(in) :: tr

    !Mesh
    type(grid_structure), intent(in) :: mesh

    !Barycentric coordinates
    real (r8):: bar_coord_tr(1:3)
    real (r8):: b(1:3)

    !Tr nodes values
    real (r8):: p1(1:3)
    real (r8):: p2(1:3)
    real (r8):: p3(1:3)

    ! 2* Area of the triangle
    real (r8):: area

    p1=mesh%v(mesh%tr(tr)%v(1))%p
    p2=mesh%v(mesh%tr(tr)%v(2))%p
    p3=mesh%v(mesh%tr(tr)%v(3))%p

    b(1)=abs(det(p,p2,p3))
    b(2)=abs(det(p,p3,p1))
    b(3)=abs(det(p,p1,p2))
    area=sum(b(1:3))

    bar_coord_tr=b/area

    return
  end function bar_coord_tr

  function bar_coord(p, p1, p2, p3)
    !----------------------------------------------------------
    ! BARYCENTRIC COORDINATES given 3 points on the sphere
    !
    ! Given a point in cartesian coords and 3 points (p1, p2, p3)
    ! that form a triangle, and are given in the counter-clockwise
    ! ordering, returns the barycentric coords
    !   of the point.
    !
    !   The barycentric coordinates refer to the planar triangle
    !--------------------------------------------------------------
    !Point inside triangle tr and triangle vertices
    real (r8), intent(in), dimension(1:3) :: p
    real (r8), intent(in), dimension(1:3) :: p1
    real (r8), intent(in), dimension(1:3) :: p2
    real (r8), intent(in), dimension(1:3) :: p3

    !Barycentric coordinates
    real (r8):: bar_coord(1:3)
    real (r8):: b(1:3)

    ! 2* Area of the triangle
    real (r8):: area

    b(1)=(det(p,p2,p3))
    b(2)=(det(p,p3,p1))
    b(3)=(det(p,p1,p2))
    area=sum(b(1:3)) !abs(det(p1,p2,p3))

    bar_coord=b/area

    return
  end function bar_coord

  function bar_coord_plan(p, p1, p2, p3)
    !----------------------------------------------------------
    ! BARYCENTRIC COORDINATES given 3 points on the sphere
    !
    ! Given a point in cartesian coords and 3 points (p1, p2, p3)
    ! that form a triangle returns the barycentric coords
    !   of the point. (p is not required to be in the triangle)
    !
    !   The barycentric coordinates refer to the planar triangle
    !   The point is projected to the plane and then the coordinate
    !    is calculated on the plane
    !  Should give exactly the same results as bar_coord
    !--------------------------------------------------------------
    !Point inside triangle tr and triangle vertices
    real (r8), intent(in), dimension(1:3) :: p
    real (r8), intent(in), dimension(1:3) :: p1
    real (r8), intent(in), dimension(1:3) :: p2
    real (r8), intent(in), dimension(1:3) :: p3

    !Barycentric coordinates
    real (r8):: bar_coord_plan(1:3)
    real (r8):: b(1:3)

    !Auxiliar points
    real (r8), dimension(1:3) :: q
    real (r8), dimension(1:3) :: nq

    ! 2* Area of the triangle
    real (r8):: area

    !Get normal to plane
    nq=cross_product(p2-p1, p3-p1)
    nq=nq/norm(nq)

    !Project point onto normal
    q=dot_product(p1,nq)*p/dot_product(p, nq)

    !Calculate the planar coordinates
    b(1)=(det(q,p2,p3))
    b(2)=(det(q,p3,p1))
    b(3)=(det(q,p1,p2))
    area=(det(p1,p2,p3))

    bar_coord_plan=b/area

    return
  end function bar_coord_plan

  function insidetrmesh(p, tr, mesh)
    !----------------------------------------------------------
    !	INSIDETRMESH
    !
    !	Checks if 'p' is inside the geodesical triangle of index tr
    ! from the mesh
    ! The algorithm checks if the point left of each edge
    !--------------------------------------------------------------
    !Point
    real (r8), intent(in) :: p(1:3)

    !Triangle index to be tested
    integer (i4), intent(in) :: tr

    !Mesh
    type(grid_structure), intent(in) :: mesh

    !Return true if point in triangle
    logical:: insidetrmesh

    !Aux
    integer (i4):: l
    real (r8):: vnn
    real (r8)::  determ
    real (r8), dimension(1:3) :: vn
    real (r8), dimension(1:3) :: p1
    real (r8), dimension(1:3) :: p2
    real (r8):: determ2
    real (r8):: tol
    real (r8):: tol2

    insidetrmesh=.false.
    tol=eps/10000._r8
    tol2=eps2*100._r8
    !print'(f48.32)', tol
    !For every edge
    !Check if point is at left of edge
    do l=1, 3
       p1=mesh%v(mesh%tr(tr)%v(l))%p
       p2=mesh%v(mesh%tr(tr)%v(modint(l+1,3)))%p
       !Check if point is a node
       if(norm(p-p1)<tol)then
          insidetrmesh=.true.
          return
       end if

       ! Get normal to plane formed by P1, 0, P2
       !  that is, vn=P1 x P2
       !Calculate the compoenent of p in the direction of n
       ! that is, np=p.n
       !This is equivalent to calculating
       !  the determinant of p, p1, p2
       !  and if this is positive, p is at left of p1->p2
       vn=cross_product(p1, p2)
       vnn=norm(vn)
       if(vnn<tol)then
          print*, "Warning on insidetrmesh: arc points p1 and p2 too close together", vnn
          stop
       end if
       vn=vn/vnn
       determ=dot_product(vn, p)
       if(determ < -tol )then !Determinant is negative, then
          !Point not in triangle
          return
       elseif(determ > tol)then
          !Point may be in triangle, check other edges
          cycle
       end if
       ! To be robust, this determinant has to be calculated with higher precision
       !   if it is too near zero
       if(abs(determ)  <= tol )then !Determinant is small, use quad precision
          determ2=robdet(p, p1, p2)
          !print'(2e48.32)', determ2, determ
          if(abs(determ2) < tol2 .or. determ*determ2 < 0 )then !On the edge
             !Point is on the edge, assume it may belong to this triangle
             cycle
          elseif(determ2<-tol2)then
             !Point close to edge, but outside triangle
             return
          else
             ! Point close to edge, maybe inside triangle
             cycle
          end if
       end if
    end do

    !If left <=0  to all edge, then the point
    !  is inside the triangle, or on the edge
    insidetrmesh=.true.

    return
  end function insidetrmesh

  function insidetr(p, p1, p2, p3)
    !----------------------------------------------------------
    ! INSIDETR
    !
    ! Checks if 'p' is inside the geodesical triangle formed bu p1, p2, p3
    !  The vertices of the triangle must be given ccwisely
    ! The algorithm checks if the point left of each edge
    !--------------------------------------------------------------
    !Point
    real (r8), intent(in), dimension(1:3) :: p
    real (r8), intent(in), dimension(1:3) :: p1
    real (r8), intent(in), dimension(1:3) :: p2
    real (r8), intent(in), dimension(1:3) :: p3

    !Return true if point in triangle
    logical:: insidetr

    !Aux
    real (r8):: left

    insidetr=.false.

    !For every edge
    !Check if point is at left of edge
    left=det(p, p1, p2)
    if(left< -eps )then
       !Point not in triangle
       return
    end if

    left=det(p, p2, p3)
    if(left< -eps )then
       !Point not in triangle
       return
    end if

    left=det(p, p3, p1)
    if(left< -eps )then
       !Point not in triangle
       return
    end if

    !If left <=0  to all edge, then the point
    !  is inside the triangle, or on the edge
    insidetr=.true.

    return
  end function insidetr


  function trbarycenter(tr, mesh)
    !----------------------------------------------------------
    !	BARYCENTER
    !
    !	Given a tringle 'tr' return its barycenter
    !--------------------------------------------------------------	    
    !Triangles barycenter
    real (r8):: trbarycenter(1:3)

    !Triangle index
    integer (i4), intent(in) :: tr

    !Mesh
    type(grid_structure), intent(in) :: mesh

    integer:: i

    !Tr nodes values
    real (r8):: p(1:3)

    p=(/0._r8,0._r8,0._r8/)
    do i=1,3
       p=p+mesh%v(mesh%tr(tr)%v(i))%p
    end do
    p=p/norm(p)

    trbarycenter=p

    return
  end function trbarycenter

  function vorbarycenter(hx, mesh)
    !----------------------------------------------------------
    !	VORBARYCENTER
    !
    !	Given a voronoi cell 'hx' return its barycenter
    !   The cell is usualy an hexagon or pentagon
    !--------------------------------------------------------------	    
    !Voronoi cell barycenter
    real (r8):: vorbarycenter(1:3)

    !Voronoi index
    integer (i4), intent(in) :: hx

    !Mesh
    type(grid_structure), intent(in) :: mesh

    integer:: i
    integer:: k1
    integer:: k2
    !Tr nodes values
    real (r8):: p(1:3)
    real (r8):: p1(1:3)
    real (r8):: p2(1:3)
    real (r8):: p3(1:3)
    real (r8):: area

    !Aux var for barycenter
    p=(/0._r8,0._r8,0._r8/)

    !Center node
    p1=mesh%v(hx)%p

    !For all neighbour traingles
    do i=1,mesh%v(hx)%nnb
       k1=mesh%v(hx)%tr(i)
       k2=mesh%v(hx)%tr(modint(i+1, mesh%v(hx)%nnb))
       p2=mesh%tr(k1)%c%p
       p3=mesh%tr(k2)%c%p
       area=sphtriarea(p1, p2, p3)
       p=p+area*(p1+p2+p3)
       !print*, "    area:      ", area
       !print*, "    barycenter:", (p1+p2+p3)/3
    end do
    vorbarycenter=p/norm(p)
    !print*, vorbarycenter
    !print*

    return
  end function vorbarycenter

  function vorbarycenterdens(hx, mesh)
    !----------------------------------------------------------
    !	VORBARYCENTERDENS
    !
    !	Given a voronoi cell 'hx' return its barycenter using
    !   a density function (dens_f)
    !   The cell is usualy an hexagon or pentagon
    !--------------------------------------------------------------
    !Voronoi cell barycenter with density
    real (r8):: vorbarycenterdens(1:3)

    !Voronoi index
    integer (i4), intent(in) :: hx

    !Mesh
    type(grid_structure), intent(in) :: mesh

    !Aux variables
    integer:: i
    integer:: k1
    integer:: k2
    real (r8), dimension(1:3) :: p
    real (r8), dimension(1:3) :: p1
    real (r8), dimension(1:3) :: p2
    real (r8), dimension(1:3) :: p3
    real (r8), dimension(1:3) :: dfw
    real (r8):: area
    real (r8):: wsum

    p=(/0._r8,0._r8,0._r8/)

    !print*, "Voronoi Cell:" , hx
    !print*, mesh%v(hx)%p

    !Center point
    p1=mesh%v(hx)%p
    dfw(1)=dens_f(p1, mesh%pos)

    !Sum of density weights with triangle areas
    wsum=0._r8

    !For all triangles surounding the node
    do i=1,mesh%v(hx)%nnb
       !Trinagle indexes
       k1=mesh%v(hx)%tr(i)
       k2=mesh%v(hx)%tr(modint(i+1, mesh%v(hx)%nnb))

       !Triangle vertices and associated densities
       p2=mesh%tr(k1)%c%p
       p3=mesh%tr(k2)%c%p
       dfw(2)=dens_f(p2, mesh%pos)
       dfw(3)=dens_f(p3, mesh%pos)

       !Triangle geodesical area
       area=sphtriarea(p1, p2, p3)

       !Weighted barycenter of planar triangle
       p=p+area*(p1*dfw(1)+p2*dfw(2)+p3*dfw(3))

       !Sum weighted factors
       wsum=wsum+area*sum(dfw(1:3))

    end do
    !Normalize weights
    p=p/wsum
    !Project on the sphere
    p=p/norm(p)

    !Set output variable
    vorbarycenterdens=p
    !print*, vorbarycenter
    !print*

    return
  end function vorbarycenterdens

  function arcintersec(p1, p2, q1, q2, intpt)
    !-------------------------------------------------------------
    !  ARCINTERSEC
    !   
    !  Tests for existence of intersection between 2 geodesic arcs
    !  given by p1->p2 and q1->q2. Returns a logical .true. in case of
    !  intersection. An optional returning argument is the intersection point 'intpt'
    !-------------------------------------------------------------
    real(r8), dimension(1:3), intent(in) :: p1
    real(r8), dimension(1:3), intent(in) :: p2
    real(r8), dimension(1:3), intent(in) :: q1
    real(r8), dimension(1:3), intent(in) :: q2

    !intpt -> intersectin point
    real(r8), dimension(1:3), optional :: intpt
    real(r8), dimension(1:3) :: r

    !np -> normal to (p1, p2, zero) plane
    !nq -> normal to (q1, q2, zero) plane
    !nl -> np x nq
    real(r8), dimension(1:3) :: np
    real(r8), dimension(1:3) :: nq
    real(r8), dimension(1:3) :: nl

    !arcpr
    real (r8):: arcpr
    real (r8):: arcqr
    real (r8):: arcpq
    real (r8):: normnl
    real (r8), parameter :: eps=10e-6
    logical:: arcintersec
    logical:: arcintp

    !Initate variables
    arcintersec=.false.
    if(present(intpt))then
       intpt=0._r8
    endif
    r=0._r8

    !Calculate normal to planes formed by
    !  the arc points and the origin
    np=cross_product(p1,p2)
    np=np/norm(np)
    nq=cross_product(q1,q2)
    nq=nq/norm(nq)

    !Plane intersection line/vector
    !This vector velongs to the intersection of the 2 planes 
    ! defined by np and nq 
    nl=cross_product(np,nq)

    !print*
    !print '(6f8.2)',p1,p2
    !print '(6f8.2)', q1, q2
    !print*,norm(nl)

    normnl=norm(nl)
    if(normnl<eps)then
       !Arcs are on the same great circle

       !Check if q1 is in arc (p1,p2)
       arcpq=arclen(p1,q1)+arclen(p2,q1)
       arcpr=arclen(p1,p2)
       if(arcpq<=arcpr+eps)then
          arcintersec=.true.
          return
       end if

       !Check if q2 is in arc (p1,p2)
       arcpq=arclen(p1,q2)+arclen(p2,q2)
       if(arcpq<=arcpr+eps)then
          arcintersec=.true.
          return
       end if

       !Check if p1 is in arc (q1,q2)
       arcpq=arclen(q1,p1)+arclen(q2,p1)
       arcpr=arclen(q1,q2)
       if(arcpq<=arcpr+eps)then
          arcintersec=.true.
          return
       end if

       !Check if p2 is in arc (q1,q2)
       arcpq=arclen(q1,p2)+arclen(q2,p2)
       if(arcpq<=arcpr+eps)then
          arcintersec=.true.
          return
       end if

    else
       !Arcs are on diferent great circles

       !The intersection point must be in the intersection of
       ! the planes defined by np and nq and also belong to the sphere
       !So this point is +/- nl/norm(nl)
       r=nl/norm(nl)

       arcintp=.false.
       !Check if r is in arc (p1,p2)
       arcpr=arclen(p1,r)+arclen(p2,r)
       if(arcpr<=arclen(p1,p2)+eps)then
          arcintp=.true.
       end if
       if(.not.arcintp)then
          !Check the other point (-r)
          arcpr=arclen(p1,-r)+arclen(p2,-r)
          if(arcpr<=arclen(p1,p2)+eps)then
             arcintp=.true.
             r=-r
          else
             !Neither points (+/-r) are in (p1,p2) arc
             !  There cannot be intersection
             return
          end if
       endif

       !Check if the intersec point of (p1,p2) is in arc (q1,q2)
       arcqr=arclen(q1,r)+arclen(q2,r)
       if(arcintp .and. arcqr<=arclen(q1,q2)+eps)then
          arcintersec=.true.
       end if

    end if

    if(present(intpt).and. arcintersec)then
       intpt=r
    endif

    return
  end function arcintersec

  function gcircarcintersec(n, p1, p2)
    !-------------------------------------------------------------
    !  GCIRCARCINTERSEC
    !   
    !  Calculates the intersection between a great circle, 
    !  given by its normal component (n), and a geodesic arc
    !  given by p1->p2. Returns the intersection point.
    !
    !  If no point of intersection exists, r=(0, 0, 0)
    !  If arc on great circle, r=(9,0,0)
    !-------------------------------------------------------------
    real(r8), dimension(1:3), intent(in) :: n
    real(r8), dimension(1:3), intent(in) :: p1
    real(r8), dimension(1:3), intent(in) :: p2

    !intpt -> intersectin point
    real(r8), dimension(1:3) :: gcircarcintersec
    real(r8), dimension(1:3) :: r

    !np -> normal to (p1, p2, zero) plane
    !nl -> np x nq
    real(r8), dimension(1:3) :: np
    real(r8), dimension(1:3) :: nl

    !arcpr
    real (r8):: arcpr
    real (r8):: arcp1r
    real (r8):: arcp2r
    real (r8):: arcpq
    real (r8):: normnl
    !real (r8), parameter :: eps=10e-6
    logical:: lcircarcintersec

    !Initate variables
    lcircarcintersec=.false.
    gcircarcintersec=0._r8
    r=0._r8

    !Calculate normal to plane formed by
    !  the arc points and the origin
    np=cross_product(p1,p2)

    !Plane intersection line/vector
    !This vector velongs to the intersection of the 2 planes 
    ! defined by np and n
    nl=cross_product(np,n)
    normnl=norm(nl)

    if(normnl<eps)then
       !Arc is on the great circle

       !There are infinite the intersection points
       ! intersection point is returned as (9, 0, 0)
       lcircarcintersec=.true.
       r(1)=9

    else
       !Arc not on the great circles

       !The intersection point must be in the intersection of
       ! the planes defined by np and nq and also belong to the sphere
       !So this point is +/- nl/norm(nl)
       r=nl/norm(nl)

       !Check if r is in arc (p1,p2)
       arcp1r=arclen(p1,r)
       arcp2r=arclen(p2,r)
       arcpr=arclen(p1,r)+arclen(p2,r)
       arcpq=arclen(p1,p2)
       if(arcp1r<=arcpq+eps .and. arcp2r<=arcpq+eps )then
          lcircarcintersec=.true.
       end if
       if(.not.lcircarcintersec)then
          !Check the other point (-r)
          arcp1r=arclen(p1,-r)
          arcp2r=arclen(p2,-r)
          arcpr=arclen(p1,-r)+arclen(p2,-r)
          if(arcp1r<=arcpq+eps .and. arcp2r<=arcpq+eps )then
             lcircarcintersec=.true.
             r=-r
          else
             !Neither points (+/-r) are in (p1,p2) arc
             !  There cannot be intersection
             return
          end if
       endif
    end if

    if(lcircarcintersec)then
       gcircarcintersec=r
    endif

    return
  end function gcircarcintersec

  function gcircintersec(p, np, nq)
    !-------------------------------------------------------------
    !  GCIRCINTERSEC
    !
    !  Calculates the intersection between 2 great circle,
    !  given by its normal component (np, nq) and
    !  returns the nearest intersection to the point p
    !-------------------------------------------------------------
    real(r8), dimension(1:3), intent(in) :: p
    real(r8), dimension(1:3), intent(in) :: np
    real(r8), dimension(1:3), intent(in) :: nq

    !intpt -> intersectin point
    real(r8), dimension(1:3) :: gcircintersec
    real(r8), dimension(1:3) :: r1
    real(r8), dimension(1:3) :: r2

    !nl -> np x nq
    real(r8), dimension(1:3) ::  nl

    !Initate variables
    gcircintersec=0._r8

    !Plane intersection line/vector
    !This vector velongs to the intersection of the 2 planes
    ! defined by np and n
    nl=cross_product(np,nq)

    !The intersection point must be in the intersection of
    ! the planes defined by np and nq and also belong to the sphere
    !So this point is +/- nl/norm(nl)
    r1=nl/norm(nl)
    r2=-nl/norm(nl)

    if(norm(r1-p) < norm(r2-p) ) then
       gcircintersec=r1
    else
       gcircintersec=r2
    endif

    return
  end function gcircintersec

  function ptingarc(p, p1, p2, limit)
    !-------------------------------------------------------------
    !  ptingarc
    !
    !  Evaluates if a a point p is in a geodesic arc formed
    !    by p1-p2
    !
    !  Returns true if p in arc, or false if not
    !
    !  Limit : maximum angle (rad) to be considered still in the arc
    !      ex: 10e-7
    !-------------------------------------------------------------
    real(r8), dimension(1:3), intent(in) :: p
    real(r8), dimension(1:3), intent(in) :: p1
    real(r8), dimension(1:3), intent(in) :: p2
    real(r8), intent(in) :: limit
    !real(r8) :: limitcos, limitsin

    !Logical output
    logical:: ptingarc

    !Normal - v1 x v2
    real(r8), dimension(1:3) :: nl
    real(r8):: normnl

    !Dot product between plane normal and point
    real(r8):: dp

    !Projection of point on plane
    !real(r8), dimension(1:3) :: proj

    !Angle
    real(r8):: d
    real(r8):: d1
    real(r8):: d2

    !Default return value
    ptingarc=.false.

    !Normal to plane formed by p1, p2
    nl=cross_product(p1,p2)
    normnl=norm(nl)
    !limitsin=dsin(limit)

    if(normnl<eps)then
       print*, "ptingarc warning: Points are too close or are the same"
       print*, "P1:", p1
       print*, "P2:", p2
       print*, "Norm dif:", normnl
       return
    end if

    !Calculate plane p1, 0, p2, normal
    nl=nl/normnl
    !Project the point on the plane
    dp=dot_product(nl, p)
    !proj=p-dp*nl
    !Put the projection on the sphere
    !proj=proj/norm(proj)

    !print*, abs(dp), limit, limitcos

    !Check if point is in great circle
    if(abs(dp)>eps)then
       !print*, "Point not in great circle"
       return
    end if

    !Calculate distances to points
    d1=arclen(p,p1)
    d2=arclen(p,p2)
    d=arclen(p1,p2)

    !print*, dp
    !print*, d1
    !print*, d2
    !print*, d

    !Test if point in arc
    if(d1+d2<d+limit)then
       ptingarc=.true.
    end if

    return
  end function ptingarc


  function dist2lines(p1, v1, p2, v2)
    !-------------------------------------------------------------
    !  DIST2LINES
    !
    !  Calculates the distânce between 2 striaght lines in R3
    !
    !  Line 1 passes in point p1 and is parallel to v1
    !  Line 2 passes in point p2 and is parallel to v2
    !-------------------------------------------------------------
    real(r8), dimension(1:3), intent(in) :: p1
    real(r8), dimension(1:3), intent(in) :: p2
    real(r8), dimension(1:3), intent(in) :: v1
    real(r8), dimension(1:3), intent(in) :: v2

    !distance
    real(r8):: dist2lines

    !Normal - v1 x v2
    real(r8), dimension(1:3) :: nl
    real(r8):: normnl

    !p1p2
    real(r8), dimension(1:3) :: p1p2

    !Initialize variable
    dist2lines=0._r8

    !Normal to plane
    nl=cross_product(v1,v2)
    normnl=norm(nl)

    !P2-P1
    p1p2=p2-p1

    !Calculate:
    !  Based on http://pages.pacificcoast.net/~cazelais/251/distance.pdf
    ! or equivalent: http://www.easycalculation.com/analytical/shortest-distance-between-lines.php
    if( normnl < eps ) then
       dist2lines=norm(cross_product(p1p2, v1))/norm(v1)
    else
       dist2lines=abs(dot_product(p1p2, nl))/normnl
    endif

    return
  end function dist2lines

  subroutine icos0edlist(mesh, icos0, nodeslist)
    !----------------------------------------------------------
    ! icos0edlist - List of node that belong to primary icosahedral grid
    !
    ! Creates a list describing the position of the nodes
    ! relative to the primary icosahedral
    !
    ! For each node in the list, the following values are set:
    !  edindex of icos0: if nodes is on an edge of the primmary icosahedral
    !  0 :  if the node is interior to the primmary icosahedral triangle
    ! -1 :  if node is one of the primary/secondary icosahedral vertices
    !--------------------------------------------------------------
    !Input variables
    type(grid_structure), intent(in) :: mesh

    !Icos level 0 mesh
    type(grid_structure), intent(in) :: icos0

    !List (output)
    integer (i4), allocatable:: nodeslist(:)

    !Aux
    integer(i4):: i
    integer(i4):: j


    real (r8), dimension(1:3):: p
    real (r8), dimension(1:3):: p1
    real (r8), dimension(1:3):: p2
    logical:: inarc

    if(.not.allocated(nodeslist))then
       allocate(nodeslist(1:mesh%nv))
    end if

    if(trim(mesh%kind)/="icos")then
       print*, "icos0edlist warning: mesh not icosahedral"
       print*, "  returned all nodes with value 0 (not on the edge of icos0)"
       do i=1, mesh%nv
          nodeslist(i)=0
       end do
       return
    end if

    !Create the list for output
    do i=1, mesh%nv
       if(i<43)then !Do not move primary and secondary nodes
          nodeslist(i)=-1
       end if
       p=mesh%v(i)%p
       inarc=.false.
       !print*, i
       do j=1, 30
          p1=icos0%v(icos0%ed(j)%v(1))%p
          p2=icos0%v(icos0%ed(j)%v(2))%p
          !print*, "Testing:", i, j
          inarc=ptingarc(p, p1, p2, mesh%minvdist/2._r8)
          !print*, inarc
          !print*, "Testing:", i, " with edge:", edv(j,1:2)
          if(inarc)then
             !print*, "Node:", i, " in icos0 edge:", edv(j,1:2)
             nodeslist(i)=j
             exit
          end if
       end do
       if(.not.inarc)then
          nodeslist(i)=0
          !print*, "Point:", i, " is interior icos0!"
       endif
    end do
    !do i=1, n
    !print*, intlist(i)
    !end do

    return
  end subroutine icos0edlist

  function midpoint_ed(v1, v2, mesh)
    !--------------------------------------------------------------
    ! MIDPOINT 
    !  Calculates the midpoint of a geodesic arc and some other
    !   informations about it
    !  Recieves 2 vertices indexes (v1, v2) and a mesh
    !  Returns a 'point_structure' kind
    !-------------------------------------------------------------
    type(grid_structure), intent(in) :: mesh
    integer(i4), intent(in)          :: v1
    integer(i4), intent(in)          :: v2
    type(point_structure):: midpoint_ed

    real(r8):: p(1:3)

    !Position
    p=(mesh%v(v1)%p+mesh%v(v2)%p)/2
    p=p/norm(p)
    midpoint_ed%p=p
    !Trigonometric properties
    call cart2sph(p(1), p(2), p(3), midpoint_ed%lon, midpoint_ed%lat)

    return
  end function midpoint_ed

  function midpoint_edhx(t1, t2, mesh)
    !--------------------------------------------------------------
    ! MIDPOINT 
    !  Calculates the midpoint of a geodesic arc and some other
    !   informations about it
    !  Recieves 2 triangle indexes (t1, t2) and a mesh
    !  Returns a 'point_structure' kind
    !-------------------------------------------------------------
    type(grid_structure), intent(in) :: mesh
    integer(i4), intent(in)          :: t1
    integer(i4), intent(in)          :: t2
    type(point_structure):: midpoint_edhx

    real(r8):: p(1:3)

    !Position
    p=(mesh%tr(t1)%c%p+mesh%tr(t2)%c%p)/2
    p=p/norm(p)
    midpoint_edhx%p=p
    !Trigonometric properties
    call cart2sph(p(1), p(2), p(3), midpoint_edhx%lon, midpoint_edhx%lat)

    return
  end function midpoint_edhx

  function distortionhx(cell, mesh)
    !--------------------------------------------------------------
    !  DISTORTION of Voronoi Cells
    !   Calculates the distortion of a Voronoi cell
    !
    !--------------------------------------------------------------

    !Original mesh
    type(grid_structure), intent(in) :: mesh

    !Cell index
    integer (i4), intent(in) :: cell

    !Output value
    real (r8):: distortionhx

    !Temp vars
    real (r8):: tmpmin
    real (r8):: tmpmax
    real (r8):: tmpmean
    real (r8):: tmp
    real (r8):: vtmp(1:20)
    integer (i4):: i
    integer (i4):: j
    integer (i4):: k
    integer (i4):: l

    !Calculate distortion index
    i=cell
    tmpmin=1000000.
    tmpmax=0.
    tmpmean=0.
    do j=1,mesh%v(i)%nnb
       k=mesh%v(i)%tr(j)
       l=mesh%v(i)%tr(modint(j+1, mesh%v(i)%nnb))
       vtmp(j)=arclen(mesh%tr(k)%c%p, mesh%tr(l)%c%p)
       tmpmin=min(tmpmin, vtmp(j))
       tmpmax=max(tmpmax, vtmp(j))
       tmpmean=tmpmean+vtmp(j)**2
       !tmpmean=tmpmean+vtmp(j) !**2
    end do

    !Internal product adequate mean length
    tmpmean=dsqrt(tmpmean/mesh%v(i)%nnb)

    !Miura mean length
    !tmpmean=tmpmean/mesh%v(i)%nnb

    !Tomita mean length
    !tmpmean=dsqrt(2._r8*mesh%hx(i)%areag/(3._r8*dsqrt(3._r8)))

    !Calculate distortion
    tmp=0.
    do j=1,mesh%v(i)%nnb
       tmp=tmp+(vtmp(j)-tmpmean)**2
    end do
    distortionhx=dsqrt(tmp/mesh%v(i)%nnb)/tmpmean
    !distortionhx=tmpmin/tmpmax

    return
  end function distortionhx

  function alignind(i, mesh)
    !--------------------------------------------------------------
    ! ALIGNEMENT INDEX
    !  Calculates the parallelism/alignement index for voronoi cells
    ! returns 1 if cell is not an hexagon/square
    !-------------------------------------------------------------
    type(grid_structure), intent(in) :: mesh
    integer(i4), intent(in)          :: i
    real(r8):: alignind

    integer (i4):: i0
    integer (i4):: i1
    integer (i4):: i2
    integer (i4):: i3
    integer (i4):: j

    real (r8):: tmp1
    real (r8):: tmp2
    real (r8):: tmp
    real (r8):: meanl

    tmp=0.
    meanl=0.
    !Check if cell is a rectangle or haxagon
    !The index is not prepered for octagonal cells
    if(mesh%v(i)%nnb==6)then
       do j=1,3
          !get suposed paralel edges' index
          ! i1 -> i2 // i3 -> i4
          i1=mesh%v(i)%tr(j)
          i2=mesh%v(i)%tr(j+1)
          i3=mesh%v(i)%tr(j+3)
          i0=mesh%v(i)%tr(modint(j+4,6))

          !Calculate lenghts of pared edges to verify regularity
          tmp1=arclen(mesh%tr(i1)%c%p, mesh%tr(i2)%c%p)
          tmp2=arclen(mesh%tr(i3)%c%p, mesh%tr(i0)%c%p)

          !Store edge lengths to calculate mean edge lenght
          meanl=tmp1+tmp2

          !Check for diferences in lenghts
          tmp=tmp+abs(tmp1-tmp2)

          !Calculate the distances between 2 possible aligned/parallel edges
          tmp1=arclen(mesh%tr(i1)%c%p, mesh%tr(i0)%c%p)
          tmp2=arclen(mesh%tr(i2)%c%p, mesh%tr(i3)%c%p)

          !Set the diference between distances:
          !  paralel edges with same lenght result is zero
          tmp=tmp+abs(tmp1-tmp2)
       end do
       tmp=tmp/6
       meanl=meanl/6
       tmp=tmp/meanl

    elseif(mesh%v(i)%nnb==4)then !In case of squares
       do j=1,2
          !get suposed paralel edges' index
          ! i1 -> i2 // i3 -> i4
          i1=mesh%v(i)%tr(j)
          i2=mesh%v(i)%tr(j+1)
          i3=mesh%v(i)%tr(j+2)
          i0=mesh%v(i)%tr(modint(j+3,4))

          !Calculate lenghts of pared edges to verify regularity
          tmp1=arclen(mesh%tr(i1)%c%p, mesh%tr(i2)%c%p)
          tmp2=arclen(mesh%tr(i3)%c%p, mesh%tr(i0)%c%p)

          !Store edge lengths to calculate mean edge lenght
          meanl=tmp1+tmp2

          !Check for diferences in lenghts
          tmp=tmp+abs(tmp1-tmp2)

          !Calculate the distances between 2 possible aligned/parallel edges
          tmp1=arclen(mesh%tr(i1)%c%p, mesh%tr(i0)%c%p)
          tmp2=arclen(mesh%tr(i2)%c%p, mesh%tr(i3)%c%p)

          !Set the diference between distances:
          !  paralel edges with same lenght result is zero
          tmp=tmp+abs(tmp1-tmp2)
       end do
       tmp=tmp/6
       meanl=meanl/6
       tmp=tmp/meanl
    else !Polygon not rectangle or hexagon set the value of 1
       tmp=1
    end if
    alignind=tmp

    return
  end function alignind

  function alignindlimit(percent, mesh)
    !--------------------------------------------------------------
    ! ALIGNEMENT INDEX LIMIT
    !  Given a percentage amount of badly aligned cells
    !  returns the approximate alignement index to be used as cut off
    !
    !  Sorts the cells in order by alignement index and returns
    !  the index that separates the "percent" amount
    !  of badly aligned cells
    !
    !  Alingment index must be given within mes structure
    !-------------------------------------------------------------
    !Grid structure
    type(grid_structure):: mesh

    !Percentage of non alinged cells
    real(r8):: percent

    !Alignment index threshold
    real(r8):: alignindlimit

    !Vector used to sort indexes
    real (r8), allocatable :: alignvec(:)
    real (r8), allocatable :: alignsort(:)

    !Counters
    integer (i4):: i
    integer (i4):: j

    !Temporary alignement index vector
    allocate(alignvec(1:mesh%nv))
    do i=1, mesh%nv
       !alignvec(i)=alignind(i,mesh)
       alignvec(i)=mesh%hx(i)%align
    end do

    !Test using % of most aligned/nonaligned cells
    allocate(alignsort(1:mesh%nv))
    alignsort=alignvec
    call Qsort(alignsort(1:mesh%nv))

    j=int ( (100.0-percent*100.0)*mesh%nv/100 , i4)
    print '(a,f6.2,a)', " Asked for   ", real(percent*100.0,4), &
         "% bad alignment index cut off "
    alignindlimit=  alignsort(j) - eps/10000. !mesh%meanvdist*mesh%meanvdist/10
    print'(a,f16.8)', " Alignment index limit:", alignindlimit

    j=0
    do i=1, mesh%nv
       if(alignsort(i)<alignindlimit)then
          j=j+1
       end if
    end do
    percent=real(mesh%nv-j,4)/real(mesh%nv,4)
    print '(a,i8,a,i8,a,f7.2,a )', " Considering ", int(mesh%nv-j,4), &
         " cells of ", int(mesh%nv,4), " (", percent*100 , &
         "% ) as badly aligned"
    print*

    return
  end function alignindlimit



  !===============================================================================================
  !    NUMERICAL LINEAR ALGEBRA
  !===============================================================================================

  function cross_product(a,b)
    !-----------------------------------------------------------------------
    !  CROSS_PRODUCT
    !
    !  Returns the right-handed vector cross product of two 3-vectors:  
    !				C = A x B.
    !-----------------------------------------------------------------------
    implicit none

    real (r8), intent(in):: a(1:3)
    real (r8), intent(in):: b(1:3)
    real (r8):: cross_product(1:3)

    cross_product(1) = a(2)*b(3) - a(3)*b(2)                                    
    cross_product(2) = a(3)*b(1) - a(1)*b(3)
    cross_product(3) = a(1)*b(2) - a(2)*b(1)

    return
  end function cross_product

  function det(p1, p2, p3)
    !-----------------------------------------------------------------------
    !  DET
    !
    !  Returns the determinant of the matrix made of the 3 points p1, p2, p3
    !   as columns
    !-----------------------------------------------------------------------
    real (r8), intent(in) :: p1(1:3)
    real (r8), intent(in) :: p2(1:3)
    real (r8), intent(in) :: p3(1:3)
    real (r8):: det

    det=dot_product(cross_product(p1,p2),p3)

    return
  end function det

  function robdet(p1, p2, p3)
    !-----------------------------------------------------------------------
    !  ROBDET
    !
    !  Returns the robust determinant of the matrix made of the 3 points p1, p2, p3
    !   as columns - The robust part is to ensure that the sign is correct if
    !   it is near zero.
    !-----------------------------------------------------------------------
    real (r8), intent(in) :: p1(1:3)
    real (r8), intent(in) :: p2(1:3)
    real (r8), intent(in) :: p3(1:3)
    real (r8):: robdet

    real (r16):: a(1:6)
    real (r16):: robdettmp
    real(r16), dimension(1:3) :: q1
    real(r16), dimension(1:3) :: q2
    real(r16), dimension(1:3) :: q3

    q1=p1
    q2=p2
    q3=p3
    a(1)=q1(1)*q2(2)*q3(3)
    a(2)=q2(1)*q3(2)*q1(3)
    a(3)=q3(1)*q1(2)*q2(3)
    a(4)=-q3(1)*q2(2)*q1(3)
    a(5)=-q1(1)*q3(2)*q2(3)
    a(6)=-q2(1)*q1(2)*q3(3)

    robdettmp=a(1)+a(2)+a(3)+a(4)+a(5)+a(6)
    robdet=real(robdettmp, r8)

    return
  end function robdet

  function solve3x3(l1, l2, l3, b)
    !-----------------------------------------------------------------------
    !  3x3 linear system solution
    !  li = line i
    !  b - right hand side vector
    !-----------------------------------------------------------------------
    real (r8), intent(in) :: l1(1:3)
    real (r8), intent(in) :: l2(1:3)
    real (r8), intent(in) :: l3(1:3)
    real (r8), intent(in) :: b(1:3)
    real (r8):: solve3x3(1:3) !Solution
    real (r8):: detdiv !Matrix determinant
    real (r8):: v(1:3) !Auxiliar vector

    detdiv=det(l1,l2,l3)

    !Check for null determinant
    if(abs(detdiv)<eps)then
       print*, "    solve3x3 error: null determinant"
       stop
    end if

    v(1) = (l1(2) * l2(3) * b(3) - l1(2) * b(2) * l3(3) + l1(3) * l3(2) * &
         b(2) - l1(3) * b(3) * l2(2) + b(1) * l2(2) * l3(3) - b(1) * l2(3) * l3(2))
    v(2) = (-l1(1) * l2(3) * b(3) + l1(1) * b(2) * l3(3) + l2(1) * l1(3) * &
         b(3) - l2(1) * b(1) * l3(3) + l2(3) * l3(1) * b(1) - b(2) * l3(1) * l1(3))
    v(3) =  (l3(1) * l1(2) * b(2) + l3(2) * l2(1) * b(1) - l3(2) * l1(1) * b(2) - &
         l3(1) * b(1) * l2(2) - b(3) * l2(1) * l1(2) + b(3) * l1(1) * l2(2))

    solve3x3=v/detdiv

    return
  end function solve3x3

  function solvelintri(a, n)
    !---------------------------------------------
    !Solve a nxn Triangular linear system with weighted columns
    !The A matrix be of format (n+1) x (n+1)
    !  X 0 0 0 0
    !  X X 0 0 0
    !  X X X 0 0
    !  X X X X 0
    !  Y Y Y Y Y   -> Right Hand Side
    !-----------------------------------------------
    !Matrix size
    integer (i4):: n

    !Augmented Matrix for the least square problem
    real (r8), intent(in) :: a(1:n+1,1:n+1)

    !Solution
    real (r8):: solvelintri(1:n)

    !Aux vars
    real (r8), allocatable :: x(:)
    real (r8):: tmp
    integer:: i
    integer:: j

    allocate(x(1:n))

    do i=n, 1, -1
       tmp=0
       do j=n, i+1, -1
          tmp=tmp+a(j,i)*x(j)
       end do
       x(i)=(a(n+1,i)-tmp)/a(i,i)
    end do
    solvelintri(1:n)=x(1:n)

    return
  end function solvelintri

  function rot_point (p, theta)
    !-------------------------------------------------------
    !  ROT_POINT
    !
    !   This subroutine applies a rotation 'gimble like'
    ! around x, y, z-axis respectively using angles theta(x, y, z),
    ! of the point p=(x, y, z)
    !
    ! On input:
    !       p=(x,y,z) = coordinates of a point on the unit sphere.
    !       theta_? = angles of rotation in radians
    !
    ! On output:
    !       pr=(xr,yr,zr) = rotated point
    !---------------------------------------------------------
    real (r8):: p(1:3)
    real (r8):: theta(1:3)
    real (r8):: pr(1:3, 1)
    real (r8):: rot_point(1:3)
    real(r8):: R(1:3, 1:3)

    !Save rotated point
    pr(1:3,1)=p(1:3)

    !Rotate around x axis
    R(1,1)=1; R(1,2)=0;             R(1,3) =0
    R(2,1)=0; R(2,2)=dcos(theta(1)); R(2,3) = -dsin(theta(1))
    R(3,1)=0; R(3,2)=dsin(theta(1)); R(3,3) = dcos(theta(1))
    pr=matmul(R, pr)

    !Rotate around y axis
    R(1,1)=dcos(theta(2));  R(1,2)=0; R(1,3) = dsin(theta(2))
    R(2,1)=0;              R(2,2)=1; R(2,3) =0
    R(3,1)=-dsin(theta(2)); R(3,2)=0; R(3,3) = dcos(theta(2))
    pr=matmul(R, pr)

    !Rotate around z axis
    R(1,1)=dcos(theta(3));  R(1,2)=-dsin(theta(3)); R(1,3) = 0
    R(2,1)=dsin(theta(3)); R(2,2)=dcos(theta(3)); R(2,3) = 0
    R(3,1)=0;              R(3,2)=0; R(3,3) = 1
    pr=matmul(R, pr)

    rot_point(1:3)=pr(1:3, 1)
    return
  end function rot_point

  subroutine ortogonalarc(p1, p2, p, n)
    !-------------------------------------------------------------
    !  ORTOGONALARC
    !
    !  Given an arc (p1,p2), computes the bisection point p and the
    !  normal n relative to the arc that crosses perpendicularly
    !  (p1,p2) at p.
    !-------------------------------------------------------------
    !Arc points
    real(r8), dimension(1:3), intent(in) :: p1
    real(r8), dimension(1:3), intent(in) :: p2
    !Midpoint, normal relative to arc ortogonal to (p1,p2)
    ! Obviously n is tangent to (p1,p2)
    real(r8), dimension(1:3), intent(out) :: p
    real(r8), dimension(1:3), intent(out) :: n

    p=(p1+p2)/2._r8
    p=p/norm(p)

    n=proj_vec_sphere(p-p1, p)
    n=n/norm(n)

    return
  end subroutine ortogonalarc

  subroutine choleskydecomp(A, n, L)
    !-----------------------------------------
    ! Cholesky Decomposition (A=L D L')
    !    A is assumed to be symetric
    !    A must be positive definite
    !    A needs only values on lower triangular part (i>=j)
    !    D is diagonal
    !    L is lower triangular with ones in diagonal
    ! As output, D is places in the diagonal os L, only for
    !   storing purposes
    !
    ! Wikipedia version
    !------------------------------------------
    !Dimension of the matrix
    integer(i4), intent(in) :: n

    !Matrix to be decomposed
    real(r8), intent(inout) :: A(1:n,1:n)

    !Decomposed Matrix
    !  Contains lower triangular values and
    !   on its diagonal the D matriz values
    real(r8), intent(out) :: L(1:n,1:n)

    !Diagonal Matrix
    real(r8):: D(1:n)
    real(r8):: D2(1:n,1:n)
    real(r8):: D3(1:n,1:n)
    real(r8)::  T(1:n,1:n)

    !indexes
    integer (i4):: i
    integer (i4):: j
    integer (i4):: k

    !Check ?
    logical:: check
    real(r8):: error
    character(len=16):: fmt

    !Do the LDL' decomp
    L(1:n,1:n)=0._r8
    D(1:n)=0._r8
    do j=1,n
       D(j)=A(j,j)
       do k=1,j-1
          D(j)=D(j)-D(k)*L(j,k)**2
       end do
       if(D(j)<-eps/1000)then
          !print*, "Warningcholeskydecomp:"
          !print*, "   Matrix not positive definite, diagonal up to now="
          !print*, D(1:j)
          !print*
          !elseif(abs(D(j))<eps/1000000)then
       elseif(abs(D(j))==0)then
          !print*, "Error on choleskydecomp:"
          !print*, "   Matrix not positive definite, diagonal with zero"
          !print*, D(1:j)
          !print*
          !stop
       end if

       L(j,j)=1._r8
       do i=j+1,n
          L(i,j)=A(i,j)
          do k=1,j-1
             L(i,j)=L(i,j)-L(i,k)*L(j,k)*D(k)
          end do
          if(abs(D(j))==0)then
             L(i,j)=0
          else
             L(i,j)=L(i,j)/D(j)
          end if
       end do
    end do

    !Debugging checks
    check=.false.
    if(check)then
       !Fill in symetric part of the matrix A (received blank)
       do i=1, n
          do j=i+1, n
             A(i,j)=A(j,i)
          end do
       end do
       fmt="(xxf16.8)"
       write(fmt(2:3),'(i2)') n
       print*, "A"
       print fmt,transpose(A)
       print*
       print*,"L"
       print fmt,transpose(L)
       print*
       print *, "D"
       print fmt,D
       print*
       do i=1,n
          do j=1,n
             if(i==j)then
                D2(i,j)=D(i)
                D3(i,j)=dsqrt(D(i))
             else
                D2(i,j)=0._r8
                D3(i,j)=0._r8
             end if
          end do
       end do
       T=matmul(matmul(L, D2), transpose(L))-A
       error=maxval(T(1:n,1:n))
       print*, "T=LDL'-A"
       print fmt, transpose(T)
       print*, "Max error (max T):", error
       if(abs(error)>eps)then
          print*,"Error in choleskydecomp:"
          print*,"   Decomposition with too much error, error=", error
          stop
       end if
    end if

    !Put diagonal D in the diagonal of L
    do i=1, n
       L(i,i)=D(i)
    end do

    return
  end subroutine choleskydecomp

  function condnumest(L, n)
    !---------------------------------------------------------
    ! MATRI_CONDNUMEST
    !Calculate the matrix condition number lower estimate
    !
    !Uses maximum value diveded by minimun value of the diagonal
    !  matrix generated by a cholesky LDLt decomposition
    !
    !Receives a matrix decomposed by cholesky LDLt, that is,
    !  a lower triangular matrix L and a diagonal D stored
    !  in the diagonal of L, and tha size 'n' of the matrix
    !
    !Returns tha max/min value of the diagonal
    !---------------------------------------------------------
    !Condition number estimative
    real (r8):: condnumest

    !Size the matrix
    integer (i4):: n

    !Matrix
    real (r8):: L(1:n, 1:n)

    !Auxiliar variables
    integer (i4):: i
    real (r8):: dmin
    real (r8):: dmax

    condnumest=0

    !Get max and min of diagonal
    dmin=100000.
    dmax=0.
    do i=1, n
       dmin=min(abs(L(i,i)), dmin)
       dmax=max(abs(L(i,i)), dmax)
    end do

    if(dmin==0)then
       print*, "Warning condnumest: dmin=0"
       condnumest=9999999999999.
    else
       condnumest=dmax/dmin
    end if

    return
  end function condnumest

  subroutine choleskysolve(L, x, b, n)
    !--------------------------------------------------
    !  solve cholesky decomposed system (L D L')
    !
    !  L is lower diagonal and has D on its diagonal
    !  x is the solution of L D L' x = b
    !
    !--------------------------------------------------
    !Dimension of the matrix
    integer(i4), intent(in) :: n

    !Lower triangular matrix
    real(r8), intent(inout) :: L(1:n,1:n)

    !RHS vetor
    real(r8), intent(in) :: b(1:n)

    !Solution
    real(r8), intent(out) :: x(1:n)

    !auxiliar Matrix
    real(r8):: y(1:n)
    real(r8):: D(1:n,1:n)
    real(r8)::  T(1:n,1:n)

    !indexes
    integer (i4):: i
    integer (i4):: j

    !Check ?
    logical:: check
    real(r8):: res(1:n)
    real(r8):: normres
    character(len=16):: fmt

    !Solve Ly=b
    do i=1,n
       y(i)=b(i)
       do j=1,i-1
          y(i)=y(i)-y(j)*L(i,j)
       end do
    end do

    !Solve Dz=y, and put z-> y
    do i=1,n
       if(L(i,i)==0)then
          y(i)=0
       else
          y(i)=y(i)/L(i,i)
       end if
    end do

    !Solve L'x=y
    do i=n,1, -1
       x(i)=y(i)
       do j=n,i+1,-1
          x(i)=x(i)-x(j)*L(j,i)
       end do
    end do

    !Debugging checks
    check= .false. !.true.
    if(check)then
       print*, "CholeskySolve Debug"
       fmt="(xxf16.8)"
       write(fmt(2:3),'(i1)') n
       print*, "x:"
       print fmt, x
       D(1:n,1:n)=0._r8
       do i=1,n
          D(i,i)=L(i,i)
          L(i,i)=1._r8
       end do
       T=matmul(matmul(L,D), transpose(L))
       res=b-matmul(T,x)
       print*,"Residual:"
       print fmt, res
       normres=norm(res)
       print*, "Residual norm : ", normres
       if(normres>eps)then
          print*,"Warning in choleskysolve:"
          print*,"   Residual too large, norm(residual)=", normres
       end if
    end if

    return
  end subroutine choleskysolve


  !===============================================================================================
  !    ERROR NORMS
  !===============================================================================================

  function positive(x)
    !-----------------------------------------
    ! POSITIVE
    ! If x negative return 0, else return x
    !----------------------------------------------
    real(r8), intent(in) :: x
    real(r8):: positive

    if(x<0)then
       positive=0._r8
    else
       positive=x
    end if

    return
  end function positive

  function norm(p)
    !-----------------------------------------
    ! NORM
    ! Calculates the euclidian norm of a vector
    !----------------------------------------------
    real(r8), intent(in) :: p(:)
    real(r8):: norm

    norm=dot_product( p, p)
    norm=dsqrt(norm)

    return
  end function norm

  subroutine normalize(p)
    !-----------------------------------------
    ! NORMALIZE
    ! Normalizes a vector
    !----------------------------------------------
    real(r8), intent(inout):: p(:)
    real(r8):: normp

    normp=norm(p) !sqrt(dot_product ( p(1:3), p(1:3) ))
    p=p/normp

    return
  end subroutine normalize

  function distance(x,y)
    !------------------------------------------------------------
    ! DISTANCE
    !	Calculates the distance between the points x e y in R3
    !
    !	Pedro Peixoto - Dec 2010
    !---------------------------------------------------------------

    real (r8), dimension(3)::x
    real (r8), dimension(3)::y
    real (r8), dimension(3)::z
    real (r8):: distance

    z=x-y
    distance=norm(z)

    return
  end function distance

  function error_norm_max(f, g, n)
    !-------------------------------------------
    !Calculates the maximum absolute value of
    !  f-g vector both having size 1:n
    !-------------------------------------------
    integer (i4), intent(in) :: n
    real (r8), dimension(1:n), intent(in) :: f
    real (r8), dimension(1:n), intent(in) :: g
    real (r8):: error_norm_max

    error_norm_max=maxval(abs(f(1:n)-g(1:n)))


    return
  end function error_norm_max

  function error_norm_max_rel(f, g, n)
    !-------------------------------------------
    !Calculates the maximum absolute value of
    !  f-g vector divided by the max of g
    !   both having size 1:n
    !-------------------------------------------
    integer (i4), intent(in) :: n
    real (r8), dimension(1:n), intent(in) :: f
    real (r8), dimension(1:n), intent(in) :: g
    real (r8):: error_norm_max_rel
    real (r8):: maxfg
    real (r8):: maxg

    maxfg=maxval(abs(f(1:n)-g(1:n)))
    maxg=maxval(abs(g(1:n)))
    error_norm_max_rel=maxfg/maxg

    return
  end function error_norm_max_rel

  function error_norm_2(f, g, n)
    !-------------------------------------------
    !Calculates the square root of the
    !  the sum of the squares of (f-g) vector
    !  both having size 1:n
    !-------------------------------------------
    integer (i4), intent(in) :: n
    real (r8), dimension(1:n), intent(in) :: f
    real (r8), dimension(1:n), intent(in) :: g
    real (r8):: error_norm_2

    real (r8):: sum_sq

    sum_sq=dot_product((f-g),(f-g))

    error_norm_2=dsqrt(sum_sq/real(n, r8))

    return
  end function error_norm_2

  function error_norm_2_rel(f, g)
    !-------------------------------------------
    !Calculates relative L2 error
    !  That is, the square root of the
    !  the sum of the squares of (f-g) vector
    !  divided by the norm 2 of g
    !  both having size 1:n
    !-------------------------------------------
    real (r8),  intent(in) :: f(:)
    real (r8),  intent(in) :: g(:)
    real (r8):: error_norm_2_rel

    real (r8):: sum_sq
    real (r8):: sum_sq_g

    sum_sq=dot_product((f-g),(f-g))
    sum_sq_g=dot_product(g,g)

    if(abs(sum_sq_g) == 0 )then
       print*, "error_norm_2_rel error: division by zero"
       stop
    end if
    error_norm_2_rel=dsqrt(sum_sq/sum_sq_g)

    return
  end function error_norm_2_rel

  function error_norm_1(f, g, n)
    !-------------------------------------------
    !Calculates the
    !  the sum of the absolute values of (f-g) vector
    !  both having size 1:n
    !-------------------------------------------
    integer (i4), intent(in) :: n
    real (r8), dimension(1:n), intent(in) :: f
    real (r8), dimension(1:n), intent(in) :: g
    real (r8):: error_norm_1

    error_norm_1=sum(abs(f-g))

    return
  end function error_norm_1

  !===============================================================================================
  !    ARRAY CONTROL ROUTINES - ALLOCATION
  !===============================================================================================

  function modint(i, n)
    !--------------------------------------
    ! MODINT
    !   Returns i if 1<=i<=n
    !   Returns i-n if n< i <2n
    !   Returns i+n if -n < i <= 0
    !   Returns i otherwise
    !----------------------------------------
    integer (i4), intent(in) :: i
    integer (i4), intent(in) :: n
    integer (i4):: modint

    modint=i

    if( i> n .and. i < 2*n )then
       modint = i-n
    elseif (-n < i .and. i <=0 )then
       modint = i + n
    end if

    return
  end function modint

  subroutine reallocint(array, nnew)
    !--------------------------------------
    ! reallocint
    !   Changes the size of an array of integers
    !   Loses all data in it
    !----------------------------------------
    integer (i4), dimension(:), allocatable, intent(inout) :: array
    integer (i4), intent(in) :: nnew

    integer (i4):: nold
    integer (i4):: i
    integer (i4):: status

    if(.not.allocated(array))then
       allocate(array(1:nnew), stat=status)
       if(status>0)then
          print *, "ERROR on reallocint: Allocation problem ", status
          stop
       end if
    else
       nold=ubound(array, 1)
       i=lbound(array, 1)
       !Check for dimension problems
       if(i/=1)then
          print *, "Warningreallocint: Array lbound not one ", i
          print *, "  please check array."
       end if
       !If the new array has a different size, reallocate
       if(nold/=nnew)then
          deallocate(array)
          allocate(array(1:nnew))
       end if
    end if
    return

  end subroutine reallocint

  subroutine reallocreal(array, nnew)
    !--------------------------------------
    ! reallocreal
    !   Changes the size of an array of reals
    !   Loses all data in it
    !----------------------------------------
    real (r8), dimension(:), allocatable, intent(inout) :: array
    integer (i4), intent(in) :: nnew

    integer (i4):: nold
    integer (i4):: i
    integer (i4):: status

    if(.not.allocated(array))then
       allocate(array(1:nnew), stat=status)
       if(status>0)then
          print *, "ERROR on reallocreal: Allocation problem ", status
          stop
       end if
    else
       nold=ubound(array, 1)
       i=lbound(array, 1)
       !Check for dimension problems
       if(i/=1)then
          print *, "Warning reallocreal: Array lbound not one ", i
          print *, "  please check array."
       end if
       !If the new array has a different size, reallocate
       if(nold/=nnew)then
          deallocate(array)
          allocate(array(1:nnew), stat=status)
          if(status>0)then
             print *, "ERROR on reallocreal: Allocation problem ", status
             stop
          end if
       end if
    end if
    return

  end subroutine reallocreal


  !===============================================================================================
  !    IN/OUTPUT ROUTINES
  !===============================================================================================

  subroutine printheader()
    !---------------------------------------
    !Simple header printing routine
    !------------------------------------
    print*,"-----------------------------------------------------------"
    print*,"  Numerical Analysis on a Geodesical Icosahedral Sphere    "
    print*,"                                                           "
    print*
    print*,"  Pedro Peixoto and collaborators                          "
    print*,"    Oct  2018                                              "
    print*,"-----------------------------------------------------------"
    print*
  end subroutine printheader

  subroutine printending()
    !---------------------------------------
    !Simple header printing routine
    !------------------------------------

    print*,"-----------------------------------------------------------"
    print*,"  End of program  "
    print*,"-----------------------------------------------------------"
    print*  

  end subroutine printending

  subroutine getparameters(mesh)
    !---------------------------------------------------
    ! GETPARAMETERS
    !    Reads mesh parameters from file named "mesh.par"
    !    Saves parameters on mesh structure
    !--------------------------------------------------
    type(grid_structure):: mesh
    character (len=60):: filename
    character (len=300):: buffer
    integer (i4):: fileunit
    integer:: i
    integer:: n

    !Variables for line argument reading
    character(len=60):: argv
    integer:: iargc
    integer:: nargs

    !Standard definition of the mesh caracteristics    
    mesh%kind="icos" !'scvt' or 'ico' or 'icop'
    mesh%nv=100
    mesh%loadable=0

    !Standard mesh parameters file
    filename=trim(pardir)//"mesh.par"

    n=0
    nargs=iargc()
    select case (nargs)
    case(0) !If no arguments are given, read file "mesh.par"
    case(1) !If a file is given, use it
       call getarg(1, argv)
       write(filename,'(a)') argv
    case(2) !If a file is given, and a number of mesh points
       call getarg(1, argv)
       write(filename,'(a)') argv
       !This number of grid points will overlap the file given one
       call getarg(2, argv)
       read (argv, *) n
    end select
    print*,"Mesh parameters: ", trim(filename)
    print*
    call getunit(fileunit)

    !A parameters file must exist 
    open(fileunit,file=filename,status='old')

    read(fileunit,*)  buffer
    read(fileunit,*)  buffer
    read(fileunit,*)  mesh%nv
    read(fileunit,*)  buffer
    read(fileunit,*)  mesh%kind
    read(fileunit,*)  buffer
    read(fileunit,*)  mesh%pos
    read(fileunit,*)  buffer
    read(fileunit,*)  mesh%optm
    read(fileunit,*)  buffer
    read(fileunit,*)  mesh%loadable
    read(fileunit,*)  buffer
    read(fileunit,*)  i  !showonscreen
    read(fileunit,*)  buffer
    read(fileunit,*)  simulcase
    read(fileunit,*)  buffer
    read(fileunit,*)  mesh%name
    read(fileunit,*)  buffer
    read(fileunit,*)  mesh%hrchy

    close(fileunit)

    if(i==1) showonscreen=.true.
    if(n>0) mesh%nv=n

    return
  end subroutine getparameters

  !===============================================================================================
  !    READ, WRITE, PRINT MESH - SAVING AND LOADING A MESH
  !===============================================================================================

  subroutine namegrid(mesh)
    !---------------------------------------------------------------------
    !
    ! NAMEGRID
    !
    ! Sets a name for the grid bases on type, number of vertices, ...
    ! Name is returned in mesh%name
    !---------------------------------------------------------------------
    type(grid_structure), intent(inout) :: mesh

    character (len=60):: an
    character (len=60):: ah

    !Set number of noveds/level
    select case(trim(mesh%kind))
    case('icos', 'octg')
       !Set name based on glevel
       write(an,'(i8)') mesh%glevel
    case default
       !Set named based on number of nodes
       write(an,'(i8)') mesh%nv
    end select

    !Set hierarchycal flag
    if(mesh%glevel==0)then
       write(ah,'(i2)') 0
    else
       write(ah,'(i2)') mesh%hrchy
    end if

    if(trim(mesh%kind)=="read")then
       mesh%name=trim(mesh%name)
    elseif(trim(mesh%optm)/="nopt")then
       mesh%name=trim(mesh%kind)//"_"//trim(mesh%pos)//"_"// &
            trim(mesh%optm)//"_h"//trim(adjustl(ah))//"_"//trim(adjustl(an))
    elseif(mesh%glevel==0)then
       mesh%name=trim(mesh%kind)//"_"//trim(mesh%pos)//"_"//trim(adjustl(an))
    else
       mesh%name=trim(mesh%kind)//"_"//trim(mesh%pos)//"_"// &
            trim(mesh%optm)//"_"//trim(adjustl(an))
    end if

  end subroutine namegrid

  subroutine printmesh(mesh)
    !-------------------------------------------
    !PRINTMESH
    ! Print main mesh caracteristics on screen
    ! Do not make confusion with meshprint routine, that
    !  writes all mesh structure onto txt files
    !-------------------------------------------
    type(grid_structure):: mesh

    print*
    select case(trim(mesh%kind))
    case("rand")
       !"SCVT (Spherical Centroidal Voronoi Tesselation)"
       print'(a)', " Mesh        : Random Density Based Point"
    case("icos")
       print'(a)', " Mesh        : Icosahedral points "
    case("octg")
       print'(a)', " Mesh        : Octogonal points "
    case("read")
       print'(a)', " Mesh        : Read from file "//trim(mesh%filename)
       return
    case default
       print'(a)', " PRINTMESH ERROR: Invalid mesh kind : ", mesh%kind
       stop
    end select

    select case(trim(mesh%pos))
    case("eqs")
       print'(a)', " Position    : EQS (Hemesphere Symmetry)"
    case("pol")
       print'(a)', " Position    : POL (Nodes on South and North Poles)"
    case("ran")
       print'(a)', " Position    : RAN (All nodes random)"
    case("ref")
       print'(a)', " Position    : REF (Random nodes with local refinement)"
    case default
       print'(a)', " PRINTMESH ERROR: Invalid mesh position : ", trim(mesh%pos)
       stop
    end select

    select case(trim(mesh%optm))
    case("nopt")
       print '(a)', " Optimization: No mesh optimization"
    case("scvt")
       print '(a)', " Optimization: Spherical Centroidal Voronoi Tesselation"
    case("salt")
       print '(a)', " Optimization: Spherical Aligned Tesselation "
    case("sprg")
       print '(a)', " Optimization: Spring Dynamics "
    case("hr95")
       print '(a)', " Optimization: HR1995 optimization using Miura's algorithm "
    case default
       print '(a)', " Optimization:  "//trim(mesh%optm)
    end select

    print '(a,i8)', " G-level     :", mesh%glevel
    print '(a,i8)', " Vertices    :", mesh%nv
    print '(a,i8)', " Edges       :", mesh%ne
    print '(a,i8)', " Triangles   :", mesh%nt
    print '(a,l8)', " Loadable    ?", mesh%loadable==1

    return
  end subroutine printmesh

  subroutine meshstore(mesh, header)
    !Store mesh in disk
    type(grid_structure), intent(in) :: mesh
    character (len=128)::  header

    !Save mesh data for future uses
    print*, "Saving generated grid in binary format (.dat)"
    print*, "Directory: ", trim(griddir)
    print*
    call meshwrite(mesh, header)

    !Write readable mesh, for debuging purpuses
    !print*, "Saving generated grid in text format (.txt)"
    !print*, "Directory: ", trim(griddir)
    !print*
    !call meshprint(mesh)

    !Save gmt mesh plotting files
    print*, "Saving generated grid in gmt format (.gmt)"
    print*, "Directory: ", trim(griddir)
    print*
    call meshgmt(mesh)

  end subroutine meshstore

  subroutine meshread(mesh)
    !--------------------------------------------------------------
    !meshread
    !  Subroutine that reads the mesh nodes from a file
    !  The file must be in grid/
    !  The file name is given implicitily via mesh%readpath
    !--------------------------------------------------------------
    type(grid_structure), intent(inout) :: mesh

    !Names for files and numbers
    character (len=128):: filename
    character (len=128):: tmp
    character (len=3):: ext
    !character (len=256):: buffer

    !I/O unit
    integer (i4):: iunit

    !Auxiliar vars
    integer (i4):: i
    integer (i4):: j
    integer (i4):: n
    integer (i4):: m
    integer (i4):: ios
    real (r8):: lon
    real (r8):: lat
    logical:: ifile
    integer (i4), dimension(1:6):: nbtmp, nbtmp2

    !print*
    !print*, "Reading mesh nodes from file..."
    !print*

    !Read vertice list structure
    !---------------------------------------
    !filename=trim(griddir)//trim(mesh%name)//".txt"
    filename=trim(griddir)//trim(mesh%filename)

    call getunit(iunit)
    inquire(file = filename, exist = ifile)
    if(ifile) then
       open(iunit, file=filename, status='old')
       print*, "Reading nodes:", trim(filename)

    else
       print*, "meshread ERROR: ", trim(filename), " not found"
       print*, "Please set 'loadable=0' in the parameters file"
       stop
    end if

    m=len(trim(mesh%filename))
    ext=mesh%filename(m-2:m)

    if(ext=="xyz")then
       !Read number of nodes
       read(iunit, *, IOSTAT=ios) n
       if(ios>0 .or. n<2)then
          print*, "Mesh read error: Could not get number of nodes on first line"
          stop
       end if
    elseif(ext=="gmt")then
       !Find out the number of nodes
       n = 0
       do
          read(iunit,*,IOSTAT=ios) tmp
          if (ios /= 0) exit
          !print*, tmp
          n = n + 1
       end do
       rewind(iunit)
    else
       print*, "Mesh read error: Unknown kind of file format. Extension: ", ext
       stop
    end if

    mesh%nv=n
    if(.not.allocated(mesh%v)) allocate(mesh%v(1:n))
    print*, "Number of nodes : ", n

    !Read coordinates
    if(ext=="xyz")then
       do i=1, mesh%nv
          read(iunit, *) mesh%v(i)%p(1), mesh%v(i)%p(2), mesh%v(i)%p(3)
          !print*, mesh%v(i)%p(1), mesh%v(i)%p(2), mesh%v(i)%p(3)
       end do
    else
       do i=1, mesh%nv
          read(iunit, *) lon, lat
          lon=lon*deg2rad
          lat=lat*deg2rad
          call sph2cart(lon, lat, mesh%v(i)%p(1), mesh%v(i)%p(2), mesh%v(i)%p(3))
          !print*, mesh%v(i)%p(1), mesh%v(i)%p(2), mesh%v(i)%p(3)
       end do
    end if

    close(iunit)

    !Read delaunay triangulation structure
    ! if it exists
    !- read neighbours of each vertex
    !---------------------------------------
    !filename=trim(griddir)//trim(mesh%name)//".txt"
    filename=trim(griddir)//mesh%filename(1:m-4)

    filename=trim(filename)//".ngb"

    call getunit(iunit)
    inquire(file = filename, exist = ifile)
    if(ifile) then
       open(iunit, file=filename, status='old')
       print*, "Reading Delaunay structure: ", trim(filename)
       print*, " Obs: Each line of the file must have exactly 6 integers (neighbours indexes)"
       mesh%deltriread=1
    else
       print*, "Files for delaunay structure not found"
       print*, filename
       print*, "Please wait, as it will be automatically generated in the next steps."
       mesh%deltriread=0
       return
    end if

    !Read neighbour indexes
    do i=1, mesh%nv
       read(iunit, *) nbtmp(1:6)
       n=0
       do j=1, 6
          if(nbtmp(j)>0)then
             n=n+1
             nbtmp2(n)=nbtmp(j)
          end if
       end do
       mesh%v(i)%nnb=n
       allocate(mesh%v(i)%nb(1:n))
       mesh%v(i)%nb(1:n)=nbtmp2(1:n)
       !print*,  mesh%v(i)%nb(1:n)
    end do

    close(iunit)

    !print*, "Mesh read from file."
    print*
  end subroutine meshread

  subroutine meshload(mesh, header)
    !--------------------------------------------------------------
    !meshload
    !  Subroutine that reads the mesh structure from grid files
    !--------------------------------------------------------------
    type(grid_structure), intent(inout) :: mesh
    character (len=128), intent(in)::  header

    !Names for files and numbers
    character (len=128):: filename
    character (len=16):: buffer

    !I/O unit
    integer (i4):: iunit

    !Auxiliar vars
    integer (i4):: i
    integer (i4):: j
    integer (i4):: k
    integer (i4):: l
    !integer (i4):: m
    integer (i4):: n
    integer (i4):: nnb

    print*
    print*, "Loading mesh..."
    print*

    !Read Header file
    !----------------------------------------
    call getunit(iunit)
    open(iunit, file=header, status='old', form='unformatted')
    read(iunit) buffer
    if(trim(mesh%kind)/=trim(buffer))then
       if(trim(buffer)/="read")then
          print*, "Warning loading mesh: mesh kind should be", trim(mesh%kind)
          print*, "   but header gives ", trim(buffer)
          !stop
       else
          print*, "Loading mesh that was originally read from file"
          print*, " Ignoring some elements of header file"
       end if
    end if
    read(iunit) buffer
    !if(trim(mesh%pos)/=trim(buffer))then
    !   print*, "Warning on loading mesh: mesh position should be", trim(mesh%pos)
    !   print*, "   but header gives ", trim(buffer)
    !   stop
    !end if
    read(iunit) buffer !mesh%optm
    read(iunit) i !mesh%glevel
    read(iunit) mesh%nv, mesh%nt, mesh%ne
    read(iunit) mesh%nlat, mesh%dlat
    read(iunit) mesh%maxvnb, mesh%maxtrsqint, &
         mesh%minvdist, mesh%maxvdist, mesh%meanvdist, &
         mesh%mincdist, mesh%maxcdist,mesh%meancdist
    read(iunit) mesh%mintrarea, mesh%maxtrarea, mesh%meantrarea, &
         mesh%minhxarea, mesh%maxhxarea, mesh%meanhxarea
    read(iunit) mesh%mintrangle, mesh%maxtrangle, mesh%meantrangle

    read(iunit) mesh%name

    close(iunit)
    !print*
    !print*, mesh%kind
    !print*, mesh%pos
    !print*, mesh%optm
    !print*, mesh%glevel
    !if(trim(mesh%kind)/="icos".and.trim(mesh%kind)/="octg".and. &
    !     trim(mesh%kind)/="rand")then
    !   print*, "meshload ERROR: Header file incorect, or missread. Mesh: ",mesh%kind
    !   stop
    !end if

    !Print mesh caracteristics
    !------------------------------------

    if(showonscreen)then
       call printmesh(mesh)
       print '(a, i6)',   " nlat:", mesh%nlat
       print '(a, f8.4)', " dlat:", mesh%dlat
       print '(a, i6)', " Max Quad Triangle intersection: ",  mesh%maxtrsqint
       print '(a, i6)', " Max Vertices Neighbours:        ",  mesh%maxvnb
       print '(a, 3f8.4)', " Min Max Mean Vertice Distances (degrees): ",  &
            mesh%minvdist*rad2deg, mesh%maxvdist*rad2deg, mesh%meanvdist*rad2deg
       print*
    end if

    !Read vertice list structure
    !---------------------------------------
    filename=trim(griddir)//trim(mesh%name)//"_vert_coord.dat"
    call openfile(filename, iunit)

    k=0
    read(iunit) n
    if(mesh%nv==n)then
       if(.not.allocated(mesh%v)) then
          allocate(mesh%v(1:n))
       end if
    else
       print*, "meshload ERROR: Inconsistent header and data file number of vertices "
       print*, "Header: ",  mesh%nv, " Datafile: ", n
       stop
    end if
    !Read coordinates
    read(iunit) (k, mesh%v(i)%p(1), mesh%v(i)%p(2), mesh%v(i)%p(3), &
         mesh%v(i)%lon, mesh%v(i)%lat, &
         mesh%v(i)%nnb, i=1, mesh%nv)
    close(iunit)

    filename=trim(griddir)//trim(mesh%name)//"_vert_nb.dat"
    call openfile(filename, iunit)
    !Read neighbour lists
    do i=1, mesh%nv
       call reallocint(mesh%v(i)%nb, mesh%v(i)%nnb)
       call reallocreal(mesh%v(i)%nbd, mesh%v(i)%nnb)
       call reallocreal(mesh%v(i)%nbdg, mesh%v(i)%nnb)
       call reallocint(mesh%v(i)%ed, mesh%v(i)%nnb)
       call reallocint(mesh%v(i)%tr, mesh%v(i)%nnb)

       read(iunit) (k, l, nnb, mesh%v(i)%nb(j), mesh%v(i)%nbd(j), &
            mesh%v(i)%nbdg(j), mesh%v(i)%ed(j), mesh%v(i)%tr(j), &
            j=1,mesh%v(i)%nnb)
       if(l < nnb .or. nnb/=mesh%v(i)%nnb) then
          print*, "meshload ERROR: vertice nb list not well structured"
          print*, "Header: ",mesh%v(i)%nnb, " Datafile: ", l, nnb
          stop
       end if
    end do
    close(iunit)

    !Read triangular edge list structure
    !---------------------------------------

    filename=trim(griddir)//trim(mesh%name)//"_ed_cc.dat"
    call openfile(filename, iunit)
    read(iunit) n
    if(mesh%ne==n)then
       allocate(mesh%ed(1:n))
    else
       print*, "meshload ERROR: Inconsistent header and data file number of edges "
       print*, "Header: ",  mesh%ne, " Datafile: ", n
       stop
    end if

    !Read center/midpoint coordinates
    read(iunit) (k, mesh%ed(i)%c%p(1), mesh%ed(i)%c%p(2), mesh%ed(i)%c%p(3), &
         mesh%ed(i)%c%lon, mesh%ed(i)%c%lat, i=1, mesh%ne)
    close(iunit)

    filename=trim(griddir)//trim(mesh%name)//"_ed_nb.dat"
    call openfile(filename, iunit)
    !Read other caracteristics
    read(iunit) (k, mesh%ed(i)%v(1), mesh%ed(i)%v(2), &
         mesh%ed(i)%sh(1), mesh%ed(i)%sh(2), &
         mesh%ed(i)%tg(1), mesh%ed(i)%tg(2), mesh%ed(i)%tg(3), &
         mesh%ed(i)%nr(1), mesh%ed(i)%nr(2), mesh%ed(i)%nr(3), &
         mesh%ed(i)%lenp, mesh%ed(i)%leng,i=1, mesh%ne)
    close(iunit)

    !Read hexagonal edge list structure
    !---------------------------------------

    filename=trim(griddir)//trim(mesh%name)//"_edhx_cc.dat"
    call openfile(filename, iunit)
    read(iunit) n
    if(mesh%ne==n)then
       allocate(mesh%edhx(1:n))
    else
       print*, "meshload ERROR: Inconsistent header and data file number of edges "
       print*, "Header: ",  mesh%ne, " Datafile: ", n
       stop
    end if

    !Read center/midpoint coordinates
    read(iunit) (k, mesh%edhx(i)%c%p(1), mesh%edhx(i)%c%p(2), mesh%edhx(i)%c%p(3), &
         mesh%edhx(i)%c%lon, mesh%edhx(i)%c%lat, i=1, mesh%ne)
    close(iunit)

    filename=trim(griddir)//trim(mesh%name)//"_edhx_nb.dat"
    call openfile(filename, iunit)
    !Read other caracteristics
    read(iunit) (k, mesh%edhx(i)%v(1), mesh%edhx(i)%v(2), &
         mesh%edhx(i)%sh(1), mesh%edhx(i)%sh(2), &
         mesh%edhx(i)%tg(1), mesh%edhx(i)%tg(2), mesh%edhx(i)%tg(3), &
         mesh%edhx(i)%nr(1), mesh%edhx(i)%nr(2), mesh%edhx(i)%nr(3), &
         mesh%edhx(i)%lenp, mesh%edhx(i)%leng, i=1, mesh%ne)
    close(iunit)


    !Read triangle list structure
    !---------------------------------------

    filename=trim(griddir)//trim(mesh%name)//"_tr_cc.dat"
    call openfile(filename, iunit)

    read(iunit) n
    if(mesh%nt==n)then
       allocate(mesh%tr(1:n))
    else
       print*, "meshload ERROR: Inconsistent header and data file number of triangles"
       print*, "Header: ", mesh%nt, " Datafile: ", n
    end if

    !Read circumcenter coordinates
    read(iunit) ( j, mesh%tr(i)%c%p(1), mesh%tr(i)%c%p(2), mesh%tr(i)%c%p(3), &
         mesh%tr(i)%c%lon, mesh%tr(i)%c%lat, &
         i=1, mesh%nt)
    !Read barycenters
    read(iunit) ( j, mesh%tr(i)%b%p(1), mesh%tr(i)%b%p(2), mesh%tr(i)%b%p(3), &
         mesh%tr(i)%b%lon, mesh%tr(i)%b%lat, &
         i=1, mesh%nt)
    close(iunit)

    filename=trim(griddir)//trim(mesh%name)//"_tr_nb.dat"
    call openfile(filename, iunit)
    !Read other caracteristics
    read(iunit) (j, mesh%tr(i)%v(1), mesh%tr(i)%v(2), mesh%tr(i)%v(3), &
         mesh%tr(i)%ed(1), mesh%tr(i)%ed(2), mesh%tr(i)%ed(3), &
         mesh%tr(i)%nb(1), mesh%tr(i)%nb(2), mesh%tr(i)%nb(3), &
         mesh%tr(i)%tg(1), mesh%tr(i)%tg(2), mesh%tr(i)%tg(3), &
         mesh%tr(i)%nr(1), mesh%tr(i)%nr(2), mesh%tr(i)%nr(3), &
         mesh%tr(i)%radp, mesh%tr(i)%radg, mesh%tr(i)%areap, mesh%tr(i)%areag, &
         mesh%tr(i)%angles(1), mesh%tr(i)%angles(2), mesh%tr(i)%angles(3), &
         i=1, mesh%nt)
    close(iunit)

    !Read hexagon list structure
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_hx.dat"
    call openfile(filename, iunit)

    !Read the number of hexagons
    read(iunit) n
    if(mesh%nv==n)then
       allocate(mesh%hx(1:n))
    else
       print*, "meshload ERROR: Inconsistent header and data file number of hexagons"
       print*, "Header: ", mesh%nv, " Datafile: ", n
       stop
    end if

    !Read data
    do i=1, mesh%nv
       !First read the barycenters
       read(iunit) l, mesh%hx(i)%b%p(1),mesh%hx(i)%b%p(2),mesh%hx(i)%b%p(3), &
            mesh%hx(i)%b%lon, mesh%hx(i)%b%lat

       !Number of edges
       l=mesh%v(i)%nnb
       !Read vector basis for each cell edge
       allocate(mesh%hx(i)%tg(1:l))
       allocate(mesh%hx(i)%nr(1:l))
       allocate(mesh%hx(i)%ttgout(1:l))
       allocate(mesh%hx(i)%tnrccw(1:l))
       read(iunit) (mesh%hx(i)%tg(j), mesh%hx(i)%nr(j), &
            mesh%hx(i)%ttgout(j), mesh%hx(i)%tnrccw(j), &
            j=1, mesh%v(i)%nnb)

       !Hx properties
       read(iunit) mesh%hx(i)%areap, mesh%hx(i)%areag, mesh%hx(i)%align
    end do
    close(iunit)

    !Read regular grid structure
    !---------------------------------------
    filename=trim(griddir)//trim(mesh%name)//"_quad.dat"
    call openfile(filename, iunit)

    read(iunit) n
    if(mesh%nlat==n)then
       allocate(mesh%qd_lat(1:n))
    else
       print*, "meshload ERROR: Inconsistent header and data file for regular grid"
       print*, "Header nlat: ", mesh%nlat, " Datafile nlat: ", n
    end if

    do j=1, mesh%nlat
       read(iunit)  mesh%qd_lat(j)%nlon
       read(iunit)  mesh%qd_lat(j)%dlon
       read(iunit)  mesh%qd_lat(j)%lat
       allocate(mesh%qd_lat(j)%qd_lon(1:mesh%qd_lat(j)%nlon))
       do i=1, mesh%qd_lat(j)%nlon
          read(iunit) mesh%qd_lat(j)%qd_lon(i)%lon, mesh%qd_lat(j)%qd_lon(i)%ntr
          allocate(mesh%qd_lat(j)%qd_lon(i)%tr(1: mesh%qd_lat(j)%qd_lon(i)%ntr))
          !write(atmp,'(i4)') mesh%qd_lat(j)%qd_lon(i)%ntr
          !fmt1='('//trim(adjustl(atmp))//'i8)'
          !fmt1='(<mesh%qd_lat(j)%qd_lon(i)%ntr> i8)'
          read(iunit) mesh%qd_lat(j)%qd_lon(i)%tr(1:mesh%qd_lat(j)%qd_lon(i)%ntr)
       end do
    end do
    close(iunit)

    print*, "Mesh loaded."
    print*
  contains

    subroutine openfile(filename, iunit)
      !Name for file
      character (len=128), intent(in):: filename
      !I/O unit
      integer (i4), intent(out) :: iunit
      !Auxiliar vars
      logical:: ifile

      call getunit(iunit)
      inquire(file = filename, exist = ifile)
      if(ifile) then
         open(iunit, file=filename, status='old', form='unformatted')
      else
         print*, "meshload ERROR: ", trim(filename), " not found"
         print*, "Please set 'loadable=0' in the parameters file"
         stop
      end if
      return
    end subroutine openfile

  end subroutine meshload

  subroutine meshwrite(mesh, header)
    !--------------------------------------------------------------
    ! MESHWRITE
    !  Subroutine that writes on files the mesh structure
    !  It is written on binary files .dat
    !--------------------------------------------------------------
    type(grid_structure), intent(in) :: mesh
    character (len=128), intent(in)::  header

    !Names for files and numbers
    character (len=128):: filename

    !I/O unit
    integer (i4):: iunit

    !Auxiliar vars
    integer (i4):: i
    integer (i4):: j

    ! Write Header file
    !-----------------------------------------------------
    call getunit(iunit)
    open(iunit, file=header, status='replace', form='unformatted')
    write(iunit) mesh%kind
    write(iunit) mesh%pos
    write(iunit) mesh%optm
    write(iunit) mesh%glevel
    write(iunit) mesh%nv, mesh%nt, mesh%ne
    write(iunit) mesh%nlat, mesh%dlat
    write(iunit) mesh%maxvnb, mesh%maxtrsqint, &
         mesh%minvdist, mesh%maxvdist, mesh%meanvdist, &
         mesh%mincdist, mesh%maxcdist, mesh%meancdist
    write(iunit) mesh%mintrarea, mesh%maxtrarea, mesh%meantrarea, &
         mesh%minhxarea, mesh%maxhxarea, mesh%meanhxarea
    write(iunit) mesh%mintrangle, mesh%maxtrangle, mesh%meantrangle
    write(iunit) mesh%name
    close(iunit)

    !Write vertice list structure
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_vert_coord.dat"
    open(iunit, file=filename, status='replace', form='unformatted')
    !Write number of vertices
    write(iunit) mesh%nv
    !Write coordinates
    write(iunit) (i, mesh%v(i)%p(1), mesh%v(i)%p(2), mesh%v(i)%p(3), &
         mesh%v(i)%lon, mesh%v(i)%lat, &
         mesh%v(i)%nnb,  i=1, mesh%nv)
    close(iunit)

    !Write vertices in ASCII format (for future loads)
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//".xyz"
    open(iunit, file=filename, status='replace')
    !Write number of vertices
    write(iunit, *) mesh%nv
    !Write coordinates
    do i=1, mesh%nv
       write(iunit, *) mesh%v(i)%p(1:3)
    end do
    close(iunit)

    !Write neighbour info
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_vert_nb.dat"
    open(iunit, file=filename, status='replace', form='unformatted')
    do i=1, mesh%nv
       !Write neighbour lists
       write(iunit) (i, j, mesh%v(i)%nnb, mesh%v(i)%nb(j), mesh%v(i)%nbd(j), &
            mesh%v(i)%nbdg(j), mesh%v(i)%ed(j), mesh%v(i)%tr(j),  &
            j=1, mesh%v(i)%nnb)
    end do
    close(iunit)

    !Write triangle's edges list structure
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_ed_cc.dat"
    open(iunit, file=filename, status='replace', form='unformatted')

    !Write number of edges
    write(iunit) mesh%ne
    !Write coordinates
    write(iunit) (i, mesh%ed(i)%c%p(1), mesh%ed(i)%c%p(2), mesh%ed(i)%c%p(3), &
         mesh%ed(i)%c%lon, mesh%ed(i)%c%lat, i=1, mesh%ne)
    close(iunit)

    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_ed_nb.dat"
    open(iunit, file=filename, status='replace', form='unformatted')
    !Write more info
    write(iunit) (i, mesh%ed(i)%v(1), mesh%ed(i)%v(2), &
         mesh%ed(i)%sh(1), mesh%ed(i)%sh(2), &
         mesh%ed(i)%tg(1), mesh%ed(i)%tg(2), mesh%ed(i)%tg(3), &
         mesh%ed(i)%nr(1), mesh%ed(i)%nr(2), mesh%ed(i)%nr(3), &
         mesh%ed(i)%lenp, mesh%ed(i)%leng, i=1, mesh%ne)
    close(iunit)

    !Write hexagonal edge list structure
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_edhx_cc.dat"
    open(iunit, file=filename, status='replace', form='unformatted')

    !Write number of edges
    write(iunit) mesh%ne
    !Write coordinates
    write(iunit) (i, mesh%edhx(i)%c%p(1), mesh%edhx(i)%c%p(2), mesh%edhx(i)%c%p(3), &
         mesh%edhx(i)%c%lon, mesh%edhx(i)%c%lat, i=1, mesh%ne)
    close(iunit)

    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_edhx_nb.dat"
    open(iunit, file=filename, status='replace', form='unformatted')
    !Write more info
    write(iunit) (i, mesh%edhx(i)%v(1), mesh%edhx(i)%v(2), &
         mesh%edhx(i)%sh(1), mesh%edhx(i)%sh(2), &
         mesh%edhx(i)%tg(1), mesh%edhx(i)%tg(2), mesh%edhx(i)%tg(3), &
         mesh%edhx(i)%nr(1), mesh%edhx(i)%nr(2), mesh%edhx(i)%nr(3), &
         mesh%edhx(i)%lenp, mesh%edhx(i)%leng, i=1, mesh%ne)
    close(iunit)

    !Write triangle list structure
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_tr_cc.dat"
    open(iunit, file=filename, status='replace', form='unformatted')

    !Write number of triangles
    write(iunit) mesh%nt
    !Write coordinates of circumcenters
    write(iunit) (i, mesh%tr(i)%c%p(1), mesh%tr(i)%c%p(2), mesh%tr(i)%c%p(3), &
         mesh%tr(i)%c%lon, mesh%tr(i)%c%lat, i=1, mesh%nt)
    !Write coordinates of barycenter
    write(iunit) (i, mesh%tr(i)%b%p(1), mesh%tr(i)%b%p(2), mesh%tr(i)%b%p(3), &
         mesh%tr(i)%b%lon, mesh%tr(i)%b%lat, i=1, mesh%nt)
    close(iunit)

    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_tr_nb.dat"
    open(iunit, file=filename, status='replace', form='unformatted')
    !Write more info
    write(iunit) (i, mesh%tr(i)%v(1), mesh%tr(i)%v(2), mesh%tr(i)%v(3), &
         mesh%tr(i)%ed(1), mesh%tr(i)%ed(2), mesh%tr(i)%ed(3), &
         mesh%tr(i)%nb(1), mesh%tr(i)%nb(2), mesh%tr(i)%nb(3), &
         mesh%tr(i)%tg(1), mesh%tr(i)%tg(2), mesh%tr(i)%tg(3), &
         mesh%tr(i)%nr(1), mesh%tr(i)%nr(2), mesh%tr(i)%nr(3), &
         mesh%tr(i)%radp, mesh%tr(i)%radg, mesh%tr(i)%areap, mesh%tr(i)%areag, &
         mesh%tr(i)%angles(1), mesh%tr(i)%angles(2), mesh%tr(i)%angles(3), &
         i=1, mesh%nt)
    close(iunit)

    !Write hexagon/pentagon (vornoi cells) list structure
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_hx.dat"
    open(iunit, file=filename, status='replace', form='unformatted')

    !Write number of voronoi cells
    write(iunit) mesh%nv
    !Write coordinates
    do i=1, mesh%nv
       !First write the barycenters
       write(iunit) i, mesh%hx(i)%b%p(1), mesh%hx(i)%b%p(2), mesh%hx(i)%b%p(3), &
            mesh%hx(i)%b%lon, mesh%hx(i)%b%lat

       !Vector basis for edges
       write(iunit) (mesh%hx(i)%tg(j), mesh%hx(i)%nr(j), &
            mesh%hx(i)%ttgout(j), mesh%hx(i)%tnrccw(j), &
            j=1, mesh%v(i)%nnb)

       !Hx properties
       write(iunit) mesh%hx(i)%areap, mesh%hx(i)%areag, mesh%hx(i)%align
    end do
    close(iunit)

    !Write regular grid structure
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_quad.dat"
    open(iunit, file=filename, status='replace',form='unformatted')
    write(iunit) mesh%nlat
    do j=1, mesh%nlat
       write(iunit) mesh%qd_lat(j)%nlon
       write(iunit) mesh%qd_lat(j)%dlon
       write(iunit) mesh%qd_lat(j)%lat
       do i=1, mesh%qd_lat(j)%nlon
          !write(atmp,'(i4)') mesh%qd_lat(j)%qd_lon(i)%ntr
          !fmt='(2i8, 2f32.16, i8,'//trim(adjustl(atmp))//'i8)'
          !Write triangle lists
          write(iunit) mesh%qd_lat(j)%qd_lon(i)%lon, mesh%qd_lat(j)%qd_lon(i)%ntr
          write(iunit) mesh%qd_lat(j)%qd_lon(i)%tr(1:mesh%qd_lat(j)%qd_lon(i)%ntr)
       end do
    end do
    !write(iunit, *) "      i       j              lon                                lat", &
    !     "                  ntr      trs  "
    close(iunit)

  end subroutine meshwrite

  subroutine meshprint(mesh)
    !--------------------------------------------------------------
    !  Subroutine that prints on files the mesh structure
    !  on human readable format
    !--------------------------------------------------------------
    type(grid_structure), intent(in) :: mesh

    !Names for failes and numbers
    character (len=100):: filename
    character (len=100):: fmt
    character (len=100):: atmp

    !I/O unit
    integer (i4):: iunit
    integer (i4):: iunit1
    integer (i4):: iunit2

    !Auxiliar vars
    integer (i4):: i
    integer (i4):: j

    ! Print Header file
    !-----------------------------------------------------
    filename=trim(griddir)//trim(mesh%name)//"_header.txt"
    call getunit(iunit)
    open(iunit, file=filename, status='replace')
    write(iunit, *) "Header file for a geodesic mesh"
    write(iunit, *) "Kind : ", trim(mesh%kind)
    write(iunit, *) "Position :", trim(mesh%pos)
    write(iunit, *) "Optimizations :",  trim(mesh%optm)
    write(iunit, *) "Number of bisections (glevel) : ", mesh%glevel
    write(iunit, *) "Geodesic grid: Vertices, Triangles, Edges"
    write(iunit, *) mesh%nv, mesh%nt, mesh%ne
    write(iunit, *) "Regular grid:  nlat,  dlat"
    write(iunit, *)  mesh%nlat, mesh%dlat
    write(iunit, *) "Info: maxvnb, maxtrsqint, minvdist, ", &
         "maxvdist, meanvdist, mincdist, maxcdist, meancdist"
    write(iunit, '(2i4, 6f12.6)') mesh%maxvnb, mesh%maxtrsqint, &
         mesh%minvdist, mesh%maxvdist, mesh%meanvdist, &
         mesh%mincdist, mesh%maxcdist,mesh%meancdist
    write(iunit, *) "Info: mintrarea, maxtrarea, meantrarea, ",&
         "minhxarea, maxhxarea, meanhxarea"
    write(iunit, '(6f12.6)') mesh%mintrarea, mesh%maxtrarea, mesh%meantrarea, &
         mesh%minhxarea, mesh%maxhxarea, mesh%meanhxarea
    write(iunit, *) "Info: mintrangle, maxtrangle, mentrangle"
    write(iunit, '(3f12.6)') mesh%mintrangle, mesh%maxtrangle, mesh%meantrangle
    write(iunit, *) "Info: name : ",  trim(mesh%name)
    close(iunit)

    !Print vertice list structure
    !---------------------------------------
    call getunit(iunit1)
    filename=trim(griddir)//trim(mesh%name)//"_vert_coord.txt"
    open(iunit1, file=filename, status='replace')

    call getunit(iunit2)
    filename=trim(griddir)//trim(mesh%name)//"_vert_nb.txt"
    open(iunit2, file=filename, status='replace')

    write(iunit1, *) mesh%nv
    write(iunit1, *) "      v         "//&
         "     x             y            z       "//&
         "   lon             lat          nnb"

    write(iunit2, *) "     v           i          nnb      "//&
         "    nb         nbd        nbdg         ed         tr "

    do i=1, mesh%nv
       !Write coordinates
       write(iunit1, '(i8, 5f16.8, i4)') i, &
            mesh%v(i)%p(1), mesh%v(i)%p(2), mesh%v(i)%p(3), &
            mesh%v(i)%lon, mesh%v(i)%lat, &
            mesh%v(i)%nnb
       !Write neighbour lists
       do j=1, mesh%v(i)%nnb
          write(iunit2, '(4i8, 2f16.8, 2i8)') i, j, mesh%v(i)%nnb, mesh%v(i)%nb(j), &
               mesh%v(i)%nbd(j), mesh%v(i)%nbdg(j), mesh%v(i)%ed(j), mesh%v(i)%tr(j)
       end do
    end do
    close(iunit1)
    close(iunit2)

    !Print triangular edge list structure
    !---------------------------------------
    call getunit(iunit1)
    filename=trim(griddir)//trim(mesh%name)//"_ed_cc.txt"
    open(iunit1, file=filename, status='replace')

    call getunit(iunit2)
    filename=trim(griddir)//trim(mesh%name)//"_ed_nb.txt"
    open(iunit2, file=filename, status='replace')

    write(iunit1, *) mesh%ne
    write(iunit1, *) "    ed   "// &
         "      ccx               ccy             ccz  "// &
         "     cclon                cclat  "
    write(iunit2, *) "       ed         v(1)         v(2)    "// &
         "    tr(1)       tr(2)         tgx"// &
         "        tgy              tgz          "// &
         "        nrx          nry            nrz "// &
         "           lenp             leng       "

    do i=1, mesh%ne
       !Write coordinates
       write(iunit1, '(i8, 5f18.6)') i, &
            mesh%ed(i)%c%p(1), mesh%ed(i)%c%p(2), mesh%ed(i)%c%p(3), &
            mesh%ed(i)%c%lon, mesh%ed(i)%c%lat
       write(iunit2, '(5i8,8f16.8)') i, &
            mesh%ed(i)%v(1), mesh%ed(i)%v(2), &
            mesh%ed(i)%sh(1), mesh%ed(i)%sh(2), &
            mesh%ed(i)%tg(1), mesh%ed(i)%tg(2), mesh%ed(i)%tg(3), &
            mesh%ed(i)%nr(1), mesh%ed(i)%nr(2), mesh%ed(i)%nr(3), &
            mesh%ed(i)%lenp, mesh%ed(i)%leng
    end do
    close(iunit1)
    close(iunit2)

    !Print hexagonal edge list structure
    !---------------------------------------
    call getunit(iunit1)
    filename=trim(griddir)//trim(mesh%name)//"_edhx_cc.txt"
    open(iunit1, file=filename, status='replace')

    call getunit(iunit2)
    filename=trim(griddir)//trim(mesh%name)//"_edhx_nb.txt"
    open(iunit2, file=filename, status='replace')

    write(iunit1, *) mesh%ne
    write(iunit1, *) "         ed   "// &
         "      ccx                       ccy                          ccz  "// &
         "                   cclon                        cclat    "
    write(iunit2, *) "    ed      v(1)         v(2)    "// &
         "    sh(1)       sh(2)       tgx          tgy          tgz"// &
         "     nrx         nry            nrz      lenp         leng   "

    do i=1, mesh%ne
       !Write coordinates
       write(iunit1, '(i8,5f16.8)') i, &
            mesh%edhx(i)%c%p(1), mesh%edhx(i)%c%p(2), mesh%edhx(i)%c%p(3), &
            mesh%edhx(i)%c%lon, mesh%edhx(i)%c%lat
       write(iunit2, '(5i8,8f16.8)') i, mesh%edhx(i)%v(1), mesh%edhx(i)%v(2), &
            mesh%edhx(i)%sh(1), mesh%edhx(i)%sh(2), &
            mesh%edhx(i)%tg(1), mesh%edhx(i)%tg(2), mesh%edhx(i)%tg(3), &
            mesh%edhx(i)%nr(1), mesh%edhx(i)%nr(2), mesh%edhx(i)%nr(3), &
            mesh%edhx(i)%lenp, mesh%edhx(i)%leng
    end do
    close(iunit1)
    close(iunit2)

    !Print triangle list structure
    !---------------------------------------
    call getunit(iunit1)
    filename=trim(griddir)//trim(mesh%name)//"_tr_cc.txt"
    open(iunit1, file=filename, status='replace')

    call getunit(iunit2)
    filename=trim(griddir)//trim(mesh%name)//"_tr_nb.txt"
    open(iunit2, file=filename, status='replace')

    write(iunit1, *) mesh%nt
    write(iunit1, *) "CIRCUMCENTERS"
    write(iunit1, *) "  tr     ccx      ccy        ccz  "// &
         "     cclon        cclat    "
    write(iunit2, *) "         tr         v(1)         v(2)        v(3)  "// &
         "     ed(1)       ed(2)       ed(3)  "// &
         "     nb(1)       nb(2)       nb(3)  "// &
         "     tg(1)       tg(2)       tg(3)  "// &
         "     nr(1)       nr(2)       nr(3)  "// &
         "      radp                    radg                   "// &
         "      areap   areag            InternalAngles(deg)"

    do i=1, mesh%nt
       !Write coordinates
       write(iunit1, '(i8, 5f16.8)') i, &
            mesh%tr(i)%c%p(1), mesh%tr(i)%c%p(2), mesh%tr(i)%c%p(3), &
            mesh%tr(i)%c%lon, mesh%tr(i)%c%lat
       write(iunit2, '(16i8, 7f16.8)') i, &
            mesh%tr(i)%v(1), mesh%tr(i)%v(2), mesh%tr(i)%v(3), &
            mesh%tr(i)%ed(1), mesh%tr(i)%ed(2), mesh%tr(i)%ed(3), &
            mesh%tr(i)%nb(1), mesh%tr(i)%nb(2), mesh%tr(i)%nb(3), &
            mesh%tr(i)%tg(1), mesh%tr(i)%tg(2), mesh%tr(i)%tg(3), &
            mesh%tr(i)%nr(1), mesh%tr(i)%nr(2), mesh%tr(i)%nr(3), &
            mesh%tr(i)%radp, mesh%tr(i)%radg, mesh%tr(i)%areap, mesh%tr(i)%areag, &
            mesh%tr(i)%angles(1:3)*rad2deg
    end do

    write(iunit1, *) "BARYCENTERS"
    write(iunit1, *) "         tr   "// &
         "      ccx         ccy                  ccz  "// &
         "          cclon            cclat    "
    do i=1, mesh%nt
       !Write coordinates
       write(iunit1, '(i8, 5f16.8)') i, &
            mesh%tr(i)%b%p(1), mesh%tr(i)%b%p(2), mesh%tr(i)%b%p(3), &
            mesh%tr(i)%b%lon, mesh%tr(i)%b%lat
    end do

    close(iunit1)
    close(iunit2)

    !Print hexagon list structure
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_hx.txt"
    open(iunit, file=filename, status='replace')

    !Write number of triangles
    write(iunit, *) mesh%nv
    write(iunit, *) "         hx   "// &
         "      bcx         bcy          bcz  "// &
         "           bclon                bclat      "// &
         "       tg_sign_corection                            "// &
         "       nr_sig_corection                            "// &
         "  areap             areag         aligmentindex"
    !Write coordinates
    do i=1, mesh%nv
       write(atmp,'(i4)') 4*mesh%v(i)%nnb
       fmt='(i8, 5f16.8, '//trim(adjustl(atmp))//'i4, 3f16.8)'
       write(iunit, fmt) i, &
            mesh%hx(i)%b%p(1), mesh%hx(i)%b%p(2), mesh%hx(i)%b%p(3), &
            mesh%hx(i)%b%lon, mesh%hx(i)%b%lat, &
            mesh%hx(i)%tg(1:mesh%v(i)%nnb), &
            mesh%hx(i)%nr(1:mesh%v(i)%nnb), &
            mesh%hx(i)%ttgout(1:mesh%v(i)%nnb), &
            mesh%hx(i)%tnrccw(1:mesh%v(i)%nnb), &
            mesh%hx(i)%areap, mesh%hx(i)%areag, &
            mesh%hx(i)%align
    end do
    close(iunit)


    !Print regular grid structure
    !---------------------------------------
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_quad.txt"
    open(iunit1, file=filename, status='replace')

    write(iunit, *) mesh%nlat
    write(iunit, *) "      i       j              lon                                lat", &
         "                  ntr      trs  "

    do j=1, mesh%nlat
       do i=1, mesh%qd_lat(j)%nlon
          write(atmp,'(i4)') mesh%qd_lat(j)%qd_lon(i)%ntr
          fmt='(2i8, 2f32.16, i8,'//trim(adjustl(atmp))//'i8)'
          !Write triangle lists
          write(iunit, fmt) i, j, mesh%qd_lat(j)%qd_lon(i)%lon*rad2deg, &
               mesh%qd_lat(j)%lat*rad2deg, mesh%qd_lat(j)%qd_lon(i)%ntr, &
               mesh%qd_lat(j)%qd_lon(i)%tr(1:mesh%qd_lat(j)%qd_lon(i)%ntr)
       end do
    end do
    close(iunit)

  end subroutine meshprint

  subroutine meshgmt(mesh)
    !------------------------------------------------
    ! Write mesh information for GMT mesh plotting
    !    All GMT outputs are in degrees
    !------------------------------------------------
    type(grid_structure), intent(in) :: mesh

    !Names
    character (len=60):: filename

    !Units
    integer:: iunit

    !Aux
    real (r8):: vlon
    real (r8):: vlat
    real (r8):: v(1:3)
    integer (i4):: i

    !Nodes/vertices longitude and latitude
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_nodes.gmt"
    !open(iunit, file=filename, status='replace', access='stream', form='unformatted')
    open(iunit, file=filename, status='replace')
    !write(iunit, *) "# Triangle vertices (lon, lat)"
    do i=1, mesh%nv
       write(iunit,'(2f16.8)') mesh%v(i)%lon*rad2deg, mesh%v(i)%lat*rad2deg
       !write(iunit) mesh%v(i)%lon*rad2deg, mesh%v(i)%lat*rad2deg
    end do
    close(iunit)

    !Triangle's circumcenters
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_trcc.gmt"
    open(iunit, file=filename, status='replace')
    !open(iunit, file=filename, status='replace', access='stream', form='unformatted')
    !write(iunit, *) "# Triangle circumcenters (lon, lat)"
    do i=1, mesh%nt
       write(iunit,'(2f16.8)') mesh%tr(i)%c%lon*rad2deg, mesh%tr(i)%c%lat*rad2deg
    end do
    close(iunit)

    !Triangular Edges
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_ed.gmt"
    open(iunit, file=filename, status='replace')
    write(iunit, '(a27)') "# Triangle Edges (lon, lat)"
    do i=1, mesh%ne
       write(iunit, '(a8,i12)') "> Edge: ", i
       write(iunit,'(2f16.8)') mesh%v(mesh%ed(i)%v(1))%lon*rad2deg, &
            mesh%v(mesh%ed(i)%v(1))%lat*rad2deg
       write(iunit,'(2f16.8)') mesh%v(mesh%ed(i)%v(2))%lon*rad2deg, &
            mesh%v(mesh%ed(i)%v(2))%lat*rad2deg
    end do
    close(iunit)

    !Edges center/midpoint
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_edc.gmt"
    open(iunit, file=filename, status='replace')
    write(iunit, '(a33)') "# Edges center/midpoint(lon, lat)"
    do i=1, mesh%ne
       write(iunit,'(2f16.8)') mesh%ed(i)%c%lon*rad2deg, &
            mesh%ed(i)%c%lat*rad2deg
    end do
    close(iunit)

    !Edges center/midpoint and normal vectors
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_ednr.gmt"
    open(iunit, file=filename, status='replace')
    write(iunit, '(a77)') "# Edges center/midpoint(lon, lat) and "// &
         "Normal vec (degrees from north, length)"
    do i=1, mesh%ne
       call convert_vec_cart2sph(mesh%ed(i)%c%p, mesh%ed(i)%nr, vlon, vlat)
       !For gmt the output must be: lon, lat, direction(degrees from north), length
       write(iunit,'(4f16.8)') mesh%ed(i)%c%lon*rad2deg, &
            mesh%ed(i)%c%lat*rad2deg, datan2(vlon, vlat)*rad2deg, &
            dsqrt(vlat**2+vlon**2)/2
    end do
    close(iunit)

    !Edges center/midpoint and tangent vectors
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_edtg.gmt"
    open(iunit, file=filename, status='replace')
    write(iunit, '(a77)') "# Edges center/midpoint(lon, lat) and "// &
         "Tangent vec (degrees from north, length)"
    do i=1, mesh%ne
       call convert_vec_cart2sph(mesh%ed(i)%c%p, mesh%ed(i)%tg, vlon, vlat)
       write(iunit,'(4f16.8)') mesh%ed(i)%c%lon*rad2deg, &
            mesh%ed(i)%c%lat*rad2deg, datan2(vlon, vlat)*rad2deg, &
            dsqrt(vlat**2+vlon**2)/2
    end do
    close(iunit)

    !Hexagonal Edges
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_edhx.gmt"
    open(iunit, file=filename, status='replace')
    write(iunit, '(a28)') "# Hexagonal Edges (lon, lat)"
    do i=1, mesh%ne
       write(iunit, '(a8,i12)') "> Edge: ", i
       write(iunit,'(2f16.8)') mesh%tr(mesh%edhx(i)%v(1))%c%lon*rad2deg, &
            mesh%tr(mesh%edhx(i)%v(1))%c%lat*rad2deg
       write(iunit,'(2f16.8)') mesh%tr(mesh%edhx(i)%v(2))%c%lon*rad2deg, &
            mesh%tr(mesh%edhx(i)%v(2))%c%lat*rad2deg
    end do
    close(iunit)

    !Hx Edges center/midpoint and normal vectors
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_edhxnr.gmt"
    open(iunit, file=filename, status='replace')
    write(iunit, '(a80)') "# Hexagonal Edges center/midpoint(lon, lat) and "// &
         "Normal vec (degrees from north, length)"
    do i=1, mesh%ne
       v=mesh%edhx(i)%nr
       call convert_vec_cart2sph(mesh%edhx(i)%c%p, v, vlon, vlat)
       write(iunit,'(4f16.8)') mesh%edhx(i)%c%lon*rad2deg, &
            mesh%edhx(i)%c%lat*rad2deg, datan2(vlon, vlat)*rad2deg, &
            dsqrt(vlat**2+vlon**2)/2
    end do
    close(iunit)

    !Hx Edges center/midpoint and tangent vectors
    call getunit(iunit)
    filename=trim(griddir)//trim(mesh%name)//"_edhxtg.gmt"
    open(iunit, file=filename, status='replace')
    write(iunit, '(a80)') "# Hexagonal Edges center/midpoint(lon, lat) and "//&
         "Tangent vec (degrees from north, length)"
    do i=1, mesh%ne
       v=mesh%edhx(i)%tg
       call convert_vec_cart2sph(mesh%edhx(i)%c%p, v, vlon, vlat)
       write(iunit,'(4f16.8)') mesh%edhx(i)%c%lon*rad2deg, &
            mesh%edhx(i)%c%lat*rad2deg, datan2(vlon, vlat)*rad2deg, &
            dsqrt(vlat**2+vlon**2)/2
    end do
    close(iunit)

    return
  end subroutine meshgmt


  !===============================================================================================
  !    RANDOM NUMBERS ROUTINES
  !===============================================================================================

  function sphrandom_denspt(refin)
    !-------------------------------------------------------------------
    !  Generate a point on the unit sphere randomly with an associated
    !  density function (dens_f)
    !  Input: seed - integer, dens_max - maximum of density function
    !------------------------------------------------------------------
    character(len=3), intent(in) :: refin
    real (r8):: sphrandom_denspt(1:3)

    real (r8):: p(1:3)
    real (r8):: df
    real (r8):: u
    integer (i4):: i
    integer (i4):: itmax

    !Maximum of iterations
    itmax=1000

    do i=1, itmax
       !Generate a ramdom point on the sphere
       p=sphrandom_pt()
       !Generate a random number in [0,1]
       !  limit of probabilities
       u=randnum() !call random_number(u)
       !Get density associates with p
       ! probability of this point to occur
       df=dens_f(p, refin)
       !If point probability is larger then random limit, keep this point
       !Else, make another point and limit
       if (u <= df) then
          exit
       endif
    end do

    !Save point
    sphrandom_denspt=p

  end function sphrandom_denspt

  function sphrandom_pt()
    !-------------------------------------------------------------------
    !  Generate a point on the unit sphere randomly with uniform
    !  density. Receives and integer seed.
    !------------------------------------------------------------------
    real (r8):: p(1:3)
    real (r8):: sphrandom_pt(1:3)
    real (r8):: normp
    integer (i4):: i

    sphrandom_pt(1:3)=(/1._r8,1._r8,1._r8/)

    !Create ramdom point on the sphere
    do i=1, 10
       p(1)=2._r8*randnum()-1._r8
       p(2)=2._r8*randnum()-1._r8
       p(3)=2._r8*randnum()-1._r8
       normp=norm(p)
       !If the norm is too small, the try again (maximum of 10 times)
       if(normp > eps ) then
          sphrandom_pt=p/norm(p)
          return
       end if
    end do

    print*, "SPHRANDOMPT WARNING: After 10 tries, could not generate random point on the sphere"
    print*, "norm=", normp, "  Returning point (1,1,1)/sqrt(3)"
    sphrandom_pt=sphrandom_pt/norm(sphrandom_pt)

  end function sphrandom_pt

  function random(my_seed)
    !--------------------------------------------------------
    ! RANDOM
    !  Random number generation, recieves any given seed
    !  and returns a number between 0,inclusive, and 1
    !  Use only when generation 1 random number only, else
    !   use setrandomseed and the use random_number instrinsic
    !--------------------------------------------------------
    real (r8):: random
    real (r8):: r
    integer (i4):: my_seed
    integer (i4):: rand_seed
    integer,allocatable :: seed(:)
    integer:: the_size
    integer:: j

    !Creates a  new seed based on clock and the seed given by my_seed
    call system_clock(rand_seed)
    rand_seed=abs(rand_seed+my_seed-my_seed*rand_seed*rand_seed/(my_seed+5))

    !Loads fortran seed array
    call random_seed(size=the_size) ! Checks seed size
    allocate(seed(the_size))        ! Allocate space
    do j=1,the_size                 ! Create seed array
       seed(j)=j+my_seed+abs(rand_seed)+(j-1)*rand_seed
    enddo
    call random_seed(put=seed)      ! Set seed
    deallocate(seed)                ! Deallocate space

    !Generates pseudo-random number
    call random_number(r)
    random=r
    return
  end function random

  subroutine setrandomseed(my_seed)
    !--------------------------------------------------------
    ! setrandomseed
    !  Random number seed generation, recieves any given seed
    !  and returns a number between 0,inclusive, and 1
    !--------------------------------------------------------

    integer (i4):: my_seed
    integer (i4):: rand_seed
    integer,allocatable :: seed(:)
    integer:: the_size
    integer:: j

    !Creates a  new seed based on clock and the seed given by my_seed
    call system_clock(rand_seed)
    rand_seed=abs(rand_seed+my_seed-my_seed*rand_seed*rand_seed/(my_seed+5))

    !Loads fortran seed array
    call random_seed(size=the_size) ! Checks seed size
    allocate(seed(the_size))        ! Allocate space
    do j=1,the_size                 ! Create seed array
       seed(j)=j+my_seed+abs(rand_seed)+(j-1)*rand_seed
    enddo
    call random_seed(put=seed)      ! Set seed
    deallocate(seed)                ! Deallocate space

    !Generates pseudo-random number
    !call random_number(r)
    !random=r
    return
  end subroutine setrandomseed

  function randnum()
    !--------------------------------------------------------
    ! RANDNUM
    !  Random number generation, recieves any given seed
    !  and returns a number between 0,inclusive, and 1
    !  It is a simple modification of the random_number instrinsic
    !   subroutine to a function
    !  A seed must have been previously set - use "setrandomseed"
    !--------------------------------------------------------
    real (r8):: randnum
    real (r8):: r

    !Generates pseudo-random number
    call random_number(r)
    randnum=r
    return
  end function randnum

  !===============================================================================================
  !    ARRAY ORDERING ROUTINES
  !===============================================================================================

  recursive subroutine Qsort(A)
    !--------------------------------------------------------
    ! Recursive Fortran 95 quicksort routine
    ! sorts real numbers into ascending numerical order
    ! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
    ! Based on algorithm from Cormen et al., Introduction to Algorithms,
    ! 1997 printing
    ! Made F conformant by Walt Brainerd
    !--------------------------------------------------------
    real(r8), intent(in out), dimension(:) :: A
    integer:: iq

    if(size(A) > 1) then
       call Partition(A, iq)
       call Qsort(A(:iq-1))
       call Qsort(A(iq:))
    endif

  contains

    subroutine Partition(A, marker)
      ! Part of quicksort
      real (r8), intent(in out), dimension(:) :: A
      integer, intent(out) :: marker
      integer:: i
      integer:: j
      real(r8):: temp
      real(r8):: x      ! pivot point
      x = A(1)
      i= 0
      j= size(A) + 1

      do
         j = j-1
         do
            if (A(j) <= x) exit
            j = j-1
         end do
         i = i+1
         do
            if (A(i) >= x) exit
            i = i+1
         end do
         if (i < j) then
            ! exchange A(i) and A(j)
            temp = A(i)
            A(i) = A(j)
            A(j) = temp
         elseif (i == j) then
            marker = i+1
            return
         else
            marker = i
            return
         endif
      end do

    end subroutine Partition
  end subroutine Qsort

  !===============================================================================================
  !    MISCELLANEOUS
  !===============================================================================================

  subroutine getunit ( iunit )
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
    integer ( i4 ):: i
    integer ( i4 ):: ios
    integer ( i4 ):: iunit
    logical:: lopen

    iunit = 0
    do i = 11, 99
       if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then
          inquire ( unit = i, opened = lopen, iostat = ios )
          if ( ios == 0 ) then
             if ( .not. lopen ) then
                iunit = i
                return
             end if
          end if
       end if
    end do

    return
  end subroutine getunit

end module smeshpack
