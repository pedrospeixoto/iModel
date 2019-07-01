module simulpack
  !=============================================================================
  !  SIMULPACK
  !
  !	Pack for several simulations and tests for SMESHPACK and INTERPACK
  !
  !	Pedro da Silva Peixoto (pedrosp@ime.usp.br)
  !	March 2013
  !=============================================================================

  !Global constants
  use constants, only: &
       datadir, &
       altdir, &
       deg2rad, &
       eps, &
       erad, &
       i4, &
       pardir, &
       pi, &
       r8, &
       rad2deg, &
       nlat_alt, &
       nlon_alt, &
       n_lat, &
       n_lon, &
       simulcase

  !Data structures
  use datastruct, only: &
       general_coords, &
       grid_structure, &
       scalar_field, &
       rbf_matrix_structure, &
       vectorinterpol_methods, &
       vector, &
       vector_field_cart

  !Spherical mesh routines
  use smeshpack, only: &
       alignind, &
       alignindlimit, &
       arcdistll, &
       arcintersec, &
       arclen, &
       calc_tiled_areas, &
       cart2sph, &
       choleskydecomp, &
       choleskysolve, &
       convert_vec_sph2cart, &
       convijll, &
       convllij, &
       distortionhx, &
       error_norm_2, &
       error_norm_max, &
       gcircarcintersec, &
       getnearnode, &
       getnearnodes, &
       gettr, &
       getunit, &
       gethxedgeconnection, &
       hxtr_intersec_areas, &
       modint, &
       norm, &
       ortogonalarc, &
       proj_vec_sphere, &
       sph2cart, &
       sphpolarea, &
       sqtriintersec, &
       vorbarycenter, &
       vorbarycenterdens

  !Interpolation pack
  use interpack, only: &
       calc_trisk_weights, &
       general_coords, &
       getmultirecon, &
       gradcalc, &
       natural_coords, &
       plot_cart_vectorfield, &
       plot_grad_vectorfield, &
       plot_scalarfield, &
       plot_scalarfield_sphgrid, &
       plot_vectorfield_sphgrid, &
       rbf_matrix_build, &
       rbf_matrix_structure, &
       rbfvec_matrix_build, &
       scalar_interpol, &
       tgrecon_index, &
       vector_interpol, &
       vector_reconstruct, &
       vrec_remap, &
       vecrecon_trsk, &
       vecrecon_trsk_tg

  !Differential operators pack
  use diffoperpack, only: &
       div_mesh_fullvec, &
       divho_mesh_fullvec


  !Intepolation for local refinement
  use refinter, only: &
    andes_density_table, &
    densftable

  implicit none

  !Staggering kind
  character (len=3):: stag

  !Function for test
  integer (i4):: testfunc

  !Logical for other tests
  logical:: discretcentroid !Discretizations to centroids
  logical:: testgrad !Test gradient
  logical:: plots !Perform plots

  !Methods to be used for reconstruction and interpolation
  type(vectorinterpol_methods) :: recon_mtd

  !Inteprolation/Reconsgtruction method used for
  !     divergence/laplacian discretizations
  type(vectorinterpol_methods) :: discinterp_mtd

  !Kind of interpolation to be used
  character (len=64):: kinterp

  !RBF shape parameter
  real (r8):: rbf_par

  !Alignment index cut off value
  real (r8) :: alignlimit

  !Percentage of aligned cells
  real (r8) :: alignpercent

  !Name for interpolation files
  character (len=128):: simulname

contains 

  subroutine getsimulpars(mesh)
    !---------------------------------------------------
    ! GETSIMULPARS
    !    Reads simulation parameters from file named "simul.par"
    !    Saves parameters on global variables
    !
    ! Pedro Peixoto Jan 2012
    !--------------------------------------------------

    !Mesh
    type(grid_structure) :: mesh

    !Buffer for strings
    character (len=300):: buffer

    !Alignement index type (value (val) or percentage (per))
    character (len=3):: aligntype

    !Alignement index threshold
    real (r8):: alignpar

    !Flag for plots
    integer:: iplots

    !Flag for reconstruction to cell centroids
    integer:: ireconcentroid

    !Flag for discretization to cell centroids
    integer:: idiscretcentroid

    !Flag for variable shape parameter
    integer(i4)::rbf_adj

    !Parameters file
    character (len=256):: filename

    !Unit for input file
    integer (i4):: fileunit
    !integer :: error

    !Couters
    integer(i4)::i
    integer(i4)::j
    character (len=16):: atmp

    !Standard definition of the interpolation tests
    testfunc=1      !Non-div case 1
    plots=.true.   !Do plots or not
    stag="HA"      !HA staggering

    !Standard interpolation parameters file
    filename=trim(pardir)//"simul.par"

    print*,"Simulation parameters (file): ", trim(filename)
    print*
    call getunit(fileunit)

    !A parameters file already must exist
    open(fileunit,file=filename,status='old')

    read(fileunit,*)  buffer
    read(fileunit,*)  buffer
    read(fileunit,*)  stag
    read(fileunit,*)  buffer
    read(fileunit,*)  kinterp
    read(fileunit,*)  buffer
    read(fileunit,*)  recon_mtd%recon
    read(fileunit,*)  buffer
    read(fileunit,*)  recon_mtd%horecon
    read(fileunit,*)  buffer
    read(fileunit,*)  discinterp_mtd%interp
    read(fileunit,*)  buffer
    read(fileunit,*)  ireconcentroid
    read(fileunit,*)  buffer
    read(fileunit,*)  idiscretcentroid
    read(fileunit,*)  buffer
    read(fileunit,*)  alignpar, aligntype
    read(fileunit,*)  buffer
    read(fileunit,*)  rbf_par, rbf_adj
    read(fileunit,*)  buffer
    read(fileunit,*)  testfunc
    read(fileunit,*)  buffer
    read(fileunit,*)  iplots

    close(fileunit)

    !Set ploting parameter
    if(iplots==1)then
       plots=.true.
    else
       plots=.false.
    end if

    !Set RBF shape parameter
    if(rbf_adj>0)then
       rbf_par=rbf_par/mesh%meanvdist
    end if

    !No higher order reconstruction method
    !    ("0" or "none" given)
    if(len(trim(recon_mtd%horecon))<2 .or. &
         trim(recon_mtd%horecon) == "none" )then
       recon_mtd%horecon=""
       recon_mtd%hyb=.false.
    else
       recon_mtd%hyb=.true.
    end if

    if(kinterp(1:1)=="0")then
      kinterp=""
    end if

    !Set interpolation for reconstructions
    recon_mtd%interp=trim(kinterp)
    recon_mtd%name=trim(recon_mtd%recon)//&
         trim(recon_mtd%horecon)//trim(recon_mtd%interp)

    !Set reconstruction position (centroid?)
    if(ireconcentroid==1)then
       recon_mtd%massc=.true.
       recon_mtd%name=trim(recon_mtd%name)//"_b"
    else
       recon_mtd%massc=.false.
    end if
    recon_mtd%rbf_par=rbf_par
    !recon_mtd%rbf_par=rbf_par/mesh%meanvdist
    !rbf_par=rbf_par/mesh%meanvdist

    discinterp_mtd%hyb=.false.
    !Set interpolation/ reconstruction method for discrete div/lap
    if(len(trim(discinterp_mtd%interp))<=2 .or. &
         trim(discinterp_mtd%interp)=="none")then
       discinterp_mtd%recon=""
       discinterp_mtd%interp=""
       discinterp_mtd%horecon=""
    else
       if(trim(stag)=="HC".or.trim(stag)=="TC")then
          discinterp_mtd%recon=discinterp_mtd%interp
          discinterp_mtd%hyb=.true.
       end if
    end if
    !Set interpolation name for discretizations
    discinterp_mtd%name=trim(discinterp_mtd%interp)
    discinterp_mtd%rbf_par=rbf_par
    ! discinterp_mtd%rbf_par=rbf_par/mesh%meanvdist

    !Set discretization position (centroid?)
    if(idiscretcentroid==1)then
       discretcentroid=.true.
       if(trim(discinterp_mtd%name)=="")then
          discinterp_mtd%name=trim(discinterp_mtd%name)//"b"
       else
          discinterp_mtd%name=trim(discinterp_mtd%name)//"_b"
       end if
    else
       discretcentroid=.false.
    end if
    discinterp_mtd%massc=discretcentroid

    !Set aligment index limit
    if(alignpar>=0)then
       if(trim(aligntype)=="per")then
          !Find index value relative to the percentage given
          alignlimit=alignindlimit(alignpar, mesh)
          alignpercent=alignpar
       else
          alignlimit=alignpar
          j=0
          do i=1, mesh%nv
             if(mesh%hx(i)%align<alignlimit+eps)then
                j=j+1
             end if
          end do
          alignpercent=real(j,r8)/real(mesh%nv,r8)
       end if
    else
       alignlimit=0
       alignpercent=0
    end if
    recon_mtd%alignlimit=alignlimit
    recon_mtd%alignpercent=alignpercent
    discinterp_mtd%alignlimit=alignlimit
    discinterp_mtd%alignpercent=alignpercent

    !Set gradient calculation parameter
    if(trim(kinterp)=="testgrad" )then
       testgrad=.true.
    else
       testgrad=.false.
    end if

    print*, "Staggering            : ", stag
    print*, "Function used         : ", testfunc
    print*, "Plots? (nplots)       : ", plots

    print*, "Interpolation method  : ", recon_mtd%interp
    print*, "Reconstruction method : ", recon_mtd%recon
    print*, "Higher order method   : ", recon_mtd%horecon
    if(recon_mtd%massc)then
       print*,"Reconstructing to mass centers"
    else
       print*,"Reconstructing to nodes or circumcenters"
    end if
    print*, "Interpolation for div mtd  : ", discinterp_mtd%interp
    if(discretcentroid)then
       print*,"Discretizing to mass centers"
    else
       print*,"Discretizing to nodes or circumcenters"
    end if
    print*, "Alignment index threshold  : ", alignlimit
    print*, "Percentage of aligned cells: ", alignpercent

    !Names depending on test cases
    select case(simulcase)
    case(3)
       simulname="meshprop"
    case(4)
       simulname="divtest"
       kinterp=trim(discinterp_mtd%name)
    case(5)
       simulname="laptest"
    case(6)
       simulname="sinterp"
    case(7)
       simulname="vinterp"
    case(8)
       simulname="vrecon"
       kinterp=trim(recon_mtd%name)
    case(12)
       simulname="vrecontg"
       kinterp=trim(recon_mtd%name)
    case(13)
       simulname="rot_test"

    case default
       simulname="simul_"
    end select

    if(simulcase/=3)then
       !Test function (vector or scalar field)
       write(atmp,'(i8)') int(testfunc)
       simulname=trim(adjustl(trim(simulname)))//"_f"//trim(adjustl(trim(atmp)))

       !Grid positioning, stagering
       simulname=trim(adjustl(trim(simulname)))//"_"//trim(adjustl(trim(stag)))

      if(trim(kinterp)/="")then
        !Interpolation methods
        simulname=trim(adjustl(trim(simulname)))//"_"//trim(kinterp)
       end if
    end if

    !Insert RBF shape parameter in name for simulations
    if(kinterp(1:3)=="rbf")then
       write(atmp,'(es8.1e1)') real(rbf_par, 4)
       !print*, real(rbf_par, 4)
       !print*, atmp
       !print*, trim(atmp)
       simulname=trim(adjustl(trim(simulname)))//"_rbf"//trim(adjustl(trim(atmp)))
    end if

    print*, "Simulation Name for Plots: ", trim(simulname)
    print*

    return
  end subroutine getsimulpars

  !======================================================================
  !    Mesh Distortions
  !======================================================================
  subroutine meshquality(mesh)
    !------------------------------------------------
    ! Mesh distorition testing routine
    !------------------------------------------------
    !Mesh
    type(grid_structure) :: mesh

    !Mesh properties variables
    ! tr=triangle
    ! hx=Voronoi cell (hexagons, pentagons)
    ! ed=edge
    type(scalar_field):: areatr
    type(scalar_field):: areahx
    type(scalar_field):: distortr
    type(scalar_field):: distorhx
    type(scalar_field):: offsettr
    type(scalar_field):: offsethx
    type(scalar_field):: alignindex
    type(scalar_field):: alignindex01
    type(scalar_field):: distshx
    type(scalar_field):: eddisp
    type(scalar_field):: edtrhxdist
    type(scalar_field):: diamhx
    type(scalar_field):: triskareas

    !Counters
    integer (i4):: i
    integer (i4):: j
    integer (i4):: j1
    integer (i4):: k
    integer (i4):: k1
    integer (i4):: l
    integer (i4):: imin
    integer (i4):: imax

    !Errors
    logical:: ifile
    integer:: iunit
    integer:: iunit2
    integer:: izero
    character (len=256):: filename, filename2

    !Field parameters and aux vars
    real (r8):: tmpmin
    real (r8):: tmpmax
    real (r8):: tmpmean

    !Aux
    real(r8):: tmp
    real(r8):: vtmp(1:30)
    real(r8):: p(1:3)

    !Min, max, mean of arrays
    real(r8):: mindistortr
    real(r8):: maxdistortr
    real(r8):: meandistortr
    real(r8):: mindistorhx
    real(r8):: maxdistorhx
    real(r8):: meandistorhx
    real(r8):: minoffsettr
    real(r8):: maxoffsettr
    real(r8):: meanoffsettr
    real(r8):: minoffsethx
    real(r8):: maxoffsethx
    real(r8):: meanoffsethx
    real(r8):: minalignindex
    real(r8):: maxalignindex
    real(r8):: meanalignindex
    real(r8):: mindistshx
    real(r8):: maxdistshx
    real(r8):: meandistshx
    real(r8):: mineddisp
    real(r8):: maxeddisp
    real(r8):: meaneddisp

    print*
    print*,"Mesh distortion testing "
    print*

    !-------------------------------------------
    !  Read parameters from file "simul.par"
    !------------------------------------------

    call getsimulpars(mesh)

    ! AREAS
    !----------------------

    !Store triangle areas in plotable variable
    areatr%n=mesh%nt
    areatr%pos=1
    allocate(areatr%f(1:areatr%n))
    do i=1,mesh%nt
       areatr%f(i)=mesh%tr(i)%areag !/mesh%maxtrarea
    end do

    !Store hexagon areas in plotable variable
    areahx%n=mesh%nv
    areahx%pos=0
    allocate(areahx%f(1:areahx%n))
    do i=1,mesh%nv
       areahx%f(i)=mesh%hx(i)%areag !/mesh%maxhxarea
       !print*, i, areahx%f(i)
    end do


    ! MASS CENTER OFFSET
    !----------------------

    !Offset between mass center and hexag centers
    offsethx%n=mesh%nv
    offsethx%pos=0

    if(trim(mesh%pos)=="readref_andes" .or. trim(mesh%pos)=="readref")then
       if(.not.allocated(mesh%densf_table))then
         allocate (mesh%densf_table(nlat_alt*nlon_alt, 3))
         call getunit(iunit2)

         if(trim(mesh%pos)=="readref_andes")then
           call andes_density_table(mesh%densf_table,iunit2)

         else !trim(mesh%pos)=="readref"
           filename2 = trim(altdir)//"densf_table.dat"
           call getunit(iunit2)
           call densftable(filename2, iunit2)
           call getunit(iunit2)
           open(iunit2, file=filename2,status='old')
             read(iunit2,*) ((mesh%densf_table(i,j), j=1,3), i=1,n_lat*n_lon)
           close(iunit2) 
         endif

       end if

         allocate(offsethx%f(1:offsethx%n))
         do i=1,mesh%nv
           p = vorbarycenterdens(i,mesh) 
           offsethx%f(i)=arclen(p,mesh%v(i)%p)
         end do
       deallocate (mesh%densf_table)
    else    
        allocate(offsethx%f(1:offsethx%n))
        do i=1,mesh%nv
          offsethx%f(i)=arclen(mesh%hx(i)%b%p,mesh%v(i)%p)
        end do
    end if


    minoffsethx=minval(offsethx%f(1:offsethx%n))
    maxoffsethx=maxval(offsethx%f(1:offsethx%n))
    meanoffsethx=sum(offsethx%f(1:offsethx%n))/offsethx%n

    !Offset between barycenter and circumcenter of triangles
    offsettr%n=mesh%nt
    offsettr%pos=1
    allocate(offsettr%f(1:offsettr%n))
    do i=1,mesh%nt
       !For each triangle
       offsettr%f(i)=arclen(mesh%tr(i)%b%p, mesh%tr(i)%c%p)
    end do

    minoffsettr=minval(offsettr%f(1:offsettr%n))
    maxoffsettr=maxval(offsettr%f(1:offsettr%n))
    meanoffsettr=sum(offsettr%f(1:offsettr%n))/offsettr%n


    ! DISTORTION
    !----------------------
    !Calculate distortion indexes for triangles
    distortr%n=mesh%nt
    distortr%pos=1
    allocate(distortr%f(1:distortr%n))
    do i=1,mesh%nt
       tmpmin=1000000.
       tmpmax=0.
       tmpmean=0.
       do j=1,3
          k=mesh%tr(i)%v(j)
          l=mesh%tr(i)%v(modint(j+1,3))
          tmp=arcdistll(mesh%v(k)%lon,mesh%v(k)%lat, &
               mesh%v(l)%lon,mesh%v(l)%lat)
          vtmp(j)=tmp
          tmpmin=min(tmpmin, tmp)
          tmpmax=max(tmpmax, tmp)
          tmpmean=tmpmean+tmp**2
       end do
       tmpmean=dsqrt(tmpmean/3)
       tmp=0.
       do j=1,3
          tmp=tmp+(vtmp(j)-tmpmean)**2
       end do
       distortr%f(i)=sqrt(tmp)/tmpmean
       !distortr%f(i)=tmpmin/tmpmax
    end do
    mindistortr=minval(distortr%f(1:distortr%n))
    maxdistortr=maxval(distortr%f(1:distortr%n))
    meandistortr=sum(distortr%f(1:distortr%n))/distortr%n

    !Calculate distortion indexes for hexagons
    !Calculate de distance index
    !Calculate mean Voronoi cell radius
    distorhx%n=mesh%nv
    distorhx%pos=0
    distshx%n=mesh%nv
    distshx%pos=0
    diamhx%n=mesh%nv
    diamhx%pos=0
    allocate(distorhx%f(1:distorhx%n))
    allocate(distshx%f(1:distshx%n))
    allocate(diamhx%f(1:diamhx%n))
    do i=1,mesh%nv

       !Calculate distortion index
       distorhx%f(i)=distortionhx(i, mesh)

       !Calculate distance index
       tmpmin=1000000.
       tmpmax=0.
       tmpmean=0.
       !Calculate Diameter
       diamhx%f(i)=0.
       do j=1,mesh%v(i)%nnb
          do k=1, mesh%v(i)%nnb
             !Triangle indexes (to get circumcenters)
             j1=mesh%v(i)%tr(j)
             k1=mesh%v(i)%tr(k)
             diamhx%f(i)=max(diamhx%f(i), &
                  arclen(mesh%tr(j1)%c%p, mesh%tr(k1)%c%p))
          end do
          vtmp(j)=mesh%v(i)%nbdg(j)
          tmpmin=min(tmpmin, tmp)
          tmpmax=max(tmpmax, tmp)
          tmpmean=tmpmean+vtmp(j)
       end do
       tmpmean=tmpmean/mesh%v(i)%nnb
       tmp=0.
       do j=1,mesh%v(i)%nnb
          tmp=tmp+(vtmp(j)-tmpmean)**2
       end do
       distshx%f(i)=dsqrt(tmp/mesh%v(i)%nnb)   
    end do

    mindistorhx=minval(distorhx%f(1:distorhx%n))
    maxdistorhx=maxval(distorhx%f(1:distorhx%n))
    meandistorhx=sum(distorhx%f(1:distorhx%n))/distorhx%n    
    mindistshx=minval(distshx%f(1:distshx%n))
    maxdistshx=maxval(distshx%f(1:distshx%n))
    meandistshx=sum(distshx%f(1:distshx%n))/distshx%n    


    ! Alignment/Paralelism index
    !-------------------------------------------------------
    alignindex%n=mesh%nv
    alignindex%pos=0
    alignindex01%n=mesh%nv
    alignindex01%pos=0
    imin=0
    imax=0
    izero=0
    maxalignindex=0._r8
    allocate(alignindex%f(1:alignindex%n))
    allocate(alignindex01%f(1:alignindex01%n))
    do i=1,mesh%nv
       if(mesh%v(i)%nnb==6.or.mesh%v(i)%nnb==4)then
          alignindex%f(i)=mesh%hx(i)%align
          !Set 1 for aligned cells and zero for others
          if(abs(alignindex%f(i))<=alignlimit)then
             izero=izero+1
             alignindex01%f(i)=1.
          else
             alignindex01%f(i)=0.
          end if
       else
          alignindex%f(i)=0.
       end if
       !if(maxalignindex < alignindex%f(i) ) then
       ! maxalignindex = alignindex%f(i)
       ! print*, i, alignindex%f(i), maxalignindex
       !end if
    end do
    minalignindex=minval(alignindex%f(1:alignindex%n))
    imin=minloc(alignindex%f(13:alignindex%n), dim=1)+12
    maxalignindex=maxval(alignindex%f(1:alignindex%n))
    imax=maxloc(alignindex%f(13:alignindex%n), dim=1)+12
    meanalignindex=sum(alignindex%f(1:alignindex%n))/alignindex%n

    !print*, "Alignment Index"
    !print*, "Almost zero:", izero, " of ", mesh%nv, " - ", izero, real(izero)/real(mesh%nv)*100,"%"
    ! print*, "Aligned HX - pt of min value:", imin, mesh%v(imin)%lon*rad2deg, mesh%v(imin)%lat*rad2deg
    ! print*, "Aligned HX - pt of max value:", imax, mesh%v(imax)%lon*rad2deg, mesh%v(imax)%lat*rad2deg
    !print*

    ! Hx and Tri Edge intersection - edge displacement
    !-------------------------------------------------------

    eddisp%n=mesh%nv
    eddisp%pos=0
    imin=0
    imax=0
    allocate(eddisp%f(1:eddisp%n))
    do i=1, mesh%nv
       eddisp%f(i)=0._r8
       do j=1,mesh%v(i)%nnb
          k=mesh%v(i)%ed(j)
          eddisp%f(i)=max(eddisp%f(i), (arclen(mesh%ed(k)%c%p, mesh%edhx(k)%c%p)/ &
               mesh%edhx(k)%leng))
       end do
    end do

    mineddisp=minval(eddisp%f(1:eddisp%n))
    imin=minloc(eddisp%f(1:eddisp%n), dim=1)
    maxeddisp=maxval(eddisp%f(1:eddisp%n))
    imax=maxloc(eddisp%f(1:eddisp%n), dim=1)
    meaneddisp=sum(eddisp%f(1:eddisp%n))/eddisp%n

    !Distance from tr edge center to Voronoi cell edge center
    edtrhxdist%n=mesh%ne
    edtrhxdist%pos=3
    allocate(edtrhxdist%f(1:edtrhxdist%n))
    do i=1, mesh%ne
       edtrhxdist%f(i)=arclen(mesh%ed(i)%c%p, mesh%edhx(i)%c%p)
    end do

    ! TRISK areas distortion - Hx and Tri areas
    !-------------------------------------------------------
    triskareas%n=mesh%nv
    triskareas%pos=0

    !Calculate areas of intersection between triangles and voronoi cells
    !  required for TRISK schemes, see Thuburn 2009
    do i=1, mesh%nv
       if(.not. allocated(mesh%hx(i)%hxtr_area))then
          allocate(mesh%hx(i)%hxtr_area(1:mesh%v(i)%nnb))
       end if
       mesh%hx(i)%hxtr_area= hxtr_intersec_areas(i, mesh%v(i)%nnb, mesh)
       ! print*, i, allocated(mesh%hx(i)%hxtr_area), mesh%hx(i)%hxtr_area(1:mesh%v(i)%nnb)
    end do

    !Set index max/min
    allocate(triskareas%f(1:triskareas%n))
    do i=1, mesh%nv
       !print*, i, allocated(mesh%hx(i)%hxtr_area), mesh%hx(i)%hxtr_area(1:mesh%v(i)%nnb)
       triskareas%f(i)=maxval(mesh%hx(i)%hxtr_area(1:mesh%v(i)%nnb))/&
            minval(mesh%hx(i)%hxtr_area(1:mesh%v(i)%nnb))
       !print*,triskareas%f(i)
    end do

    ! Save indexes
    !-------------------------------------------------------

    !Mesh Caracteristics for triangular mesh
    filename=trim(datadir)//trim(simulname)//"_indexes_tr.txt"
    call getunit(iunit)   

    inquire(file=filename, exist=ifile)
    if(ifile)then
       open(iunit,file=filename, status='old', position='append')
    else
       open(iunit,file=filename, status='replace')
       write(iunit, '(a)') &
            " nt ne opt mincdist  maxcdist meancdist "// &
            " minarea maxarea meanarea"// &
            " mindistor maxdistor meandistor"// &
            " minoffset maxoffset meanoffset"
    end if

    write(iunit, '(2i8, a12, 12f16.12)') mesh%nt, mesh%ne, trim(mesh%optm), &
         mesh%mincdist*rad2deg,  mesh%maxcdist*rad2deg,  mesh%meancdist*rad2deg, &
         mesh%mintrarea, mesh%maxtrarea, mesh%meantrarea, mindistortr, maxdistortr, meandistortr, &
         minoffsettr, maxoffsettr, meanoffsettr

    close(iunit)
    print*, "Results in : ", trim(filename)

    !Mesh Caracteristics for hexagonal mesh
    filename=trim(datadir)//trim(simulname)//"_indexes_hx.txt"
    call getunit(iunit)   

    inquire(file=filename, exist=ifile)
    if(ifile)then
       open(iunit,file=filename, status='old', position='append')
    else
       open(iunit,file=filename, status='replace')
       write(iunit, '(a)') &
            "      nv                    name                 "//&
            "    mindist        maxdist    "//&
            "      meandist        minarea         maxarea  "//&
            "      meanarea         mindistor      maxdistor      meandistor "//&
            "      minoffset      maxoffset        meanoffset      mindisthx "//&
            "      maxdisthx      meandisthx      minalignindex   maxalignindex"//&
            "   meanalignindex  n_aligned   mineddisp            maxeddisp  "//&
            "       meaneddisp"
    end if

    write(iunit, '(i8, a36, 18f16.12, i8, 3f20.12)') mesh%nv,  trim(mesh%name), &
         mesh%minvdist*rad2deg,  mesh%maxvdist*rad2deg,  mesh%meanvdist*rad2deg, &
         mesh%minhxarea, mesh%maxhxarea, mesh%meanhxarea, mindistorhx, maxdistorhx, meandistorhx, &
         minoffsethx, maxoffsethx, meanoffsethx, mindistshx, maxdistshx, meandistshx, &
         minalignindex, maxalignindex, meanalignindex, izero, mineddisp, maxeddisp, meaneddisp

    close(iunit)
    print*, "Results in : ", trim(filename)

    ! Plot fields 
    !-------------------------
    !plots=.true.
    if(plots) then
       print*, "Plotting variables ... "

       distshx%name=trim(simulname)//"_distshx"
       call plot_scalarfield(distshx, mesh)

       areahx%name=trim(simulname)//"_areahx"
       call plot_scalarfield(areahx, mesh)

       distorhx%name=trim(simulname)//"_distorhx"
       call plot_scalarfield(distorhx, mesh)

       offsethx%name=trim(simulname)//"_offsethx"
       call plot_scalarfield(offsethx, mesh)

       alignindex%name=trim(simulname)//"_alignindex"
       call plot_scalarfield(alignindex, mesh)

       alignindex01%name=trim(simulname)//"_alignindex01"
       call plot_scalarfield(alignindex01, mesh)

       areatr%name=trim(simulname)//"_areatr"
       call plot_scalarfield(areatr, mesh)

       distortr%name=trim(simulname)//"_distortr"
       call plot_scalarfield(distortr, mesh)

       offsettr%name=trim(simulname)//"_offsettr"
       call plot_scalarfield(offsettr, mesh)

       diamhx%name=trim(simulname)//"_diamhx"
       call plot_scalarfield(diamhx, mesh)

       eddisp%name=trim(simulname)//"_eddisp"
       call plot_scalarfield(eddisp, mesh)

       edtrhxdist%name=trim(simulname)//"_edtrhxdist"
       call plot_scalarfield(edtrhxdist, mesh)

       triskareas%name=trim(simulname)//"_triskareas"
       call plot_scalarfield(triskareas, mesh)
    end if

    return

  end subroutine meshquality

 subroutine meshquality_tiledareas(mesh)
    !------------------------------------------------
    ! Mesh distorition testing routine
    !------------------------------------------------
    !Mesh
    type(grid_structure) :: mesh

    !Mesh properties variables
    ! tr=triangle
    ! hx=Voronoi cell (hexagons, pentagons)
    ! ed=edge
    type(scalar_field):: areatr
    type(scalar_field):: areahx


    !Counters
    integer (i4):: i
    integer (i4):: j

    !Errors
    logical:: ifile
    integer:: iunit
    integer:: izero
    character (len=256):: filename

    !Min, max, mean of arrays
    real(r8):: min_trarea
    real(r8):: max_trarea
    real(r8):: mean_trarea
    real(r8):: min_hxarea
    real(r8):: max_hxarea
    real(r8):: mean_hxarea

    print*
    print*,"Mesh tiled area testing "
    print*

    !-------------------------------------------
    !  Read parameters from file "simul.par"
    !------------------------------------------

    call getsimulpars(mesh)
    call calc_tiled_areas(mesh)

    ! AREAS
    !----------------------

    !Store triangle areas in plotable variable
    areatr%n=mesh%nt
    areatr%pos=1
    allocate(areatr%f(1:areatr%n))
    min_trarea=100000000.0
    max_trarea=0.0
    mean_trarea=0.0
    do i=1,mesh%nt
       areatr%f(i)=mesh%tr(i)%areat !/mesh%maxtrarea
       min_trarea=min(min_trarea, areatr%f(i))
       max_trarea=max(max_trarea, areatr%f(i))
       mean_trarea=mean_trarea+areatr%f(i)
    end do
    mean_trarea=mean_trarea/mesh%nt

    !Store hexagon areas in plotable variable
    areahx%n=mesh%nv
    areahx%pos=0
    allocate(areahx%f(1:areahx%n))
    min_hxarea=100000000.0
    max_hxarea=0.0
    mean_hxarea=0.0
    do i=1,mesh%nv
       areahx%f(i)=mesh%hx(i)%areat !/mesh%maxhxarea
       min_hxarea=min(min_hxarea, areahx%f(i))
       max_hxarea=max(max_hxarea, areahx%f(i))
       mean_hxarea=mean_hxarea+areahx%f(i)
       !print*, i, areahx%f(i)
    end do
    mean_hxarea=mean_hxarea/mesh%nv

    ! Save indexes
    !-------------------------------------------------------

    !Mesh Caracteristics for tiled areas mesh
    filename=trim(datadir)//trim(simulname)//"_indexes_tiled_areas.txt"
    call getunit(iunit)

    inquire(file=filename, exist=ifile)
    if(ifile)then
       open(iunit,file=filename, status='old', position='append')
    else
       open(iunit,file=filename, status='replace')
       write(iunit, '(a)') &
            "      nv                    name                 "//&
            "      minarea_tr         maxarea_tr  meanarea_tr "//&
            "      minarea_hx         maxarea_hx  meanarea_hx "
    end if

    write(iunit, '(i8, a36, 18f16.12, i8, 3f20.12)') mesh%nv,  trim(mesh%name), &
         min_trarea, max_trarea, mean_trarea, &
         min_hxarea, max_hxarea, mean_hxarea

    close(iunit)
    print*, "Results in : ", trim(filename)

    ! Plot fields
    !-------------------------
    !plots=.true.
    if(plots) then
       print*, "Plotting variables ... "

       areahx%name=trim(simulname)//"_tareahx"
       call plot_scalarfield(areahx, mesh)

       areatr%name=trim(simulname)//"_tareatr"
       call plot_scalarfield(areatr, mesh)

    end if

    return

  end subroutine meshquality_tiledareas


  !======================================================================
  !    DIVERGENCE TESTS
  !======================================================================

  subroutine divergence_test(mesh)
    !-------------------------------------------------------------------
    ! Divergence_test
    !
    ! Test divergence numerical estimates and grid imprinting
    !--------------------------------------------------------------------

    !Mesh
    type(grid_structure) :: mesh

    !Scalar fields
    type(scalar_field):: div_ex
    type(scalar_field):: div_est
    type(scalar_field):: error
    type(scalar_field):: align

    !Velocity field
    type(vector_field_cart):: v_ex

    !Counters
    integer:: i

    !Errors
    real (r8):: error2
    real (r8):: errormax

    !Alignment index
    real (r8):: limit
    integer(i4):: countalin

    !File variables
    character (len=256):: filename
    logical:: ifile
    integer(i4):: iunit

    !Aux
    real(r8):: p(1:3)
    character (len=2):: cellopt
    logical:: average
    logical:: alignedpol

    print*
    print*,"Divergence Testing "
    print*

    !-------------------------------------------
    !  Read parameters from file "simul.par"
    !   and set default names for outputs
    !------------------------------------------

    call getsimulpars(mesh)


    print*,"Testing divergence discretization ..."

    ! Exact/analitical Values
    !-------------------------------------

    !Allocate variables

    select case(trim(stag))
    case("HA", "ha", "HB", "Hb", "HC", "hc", "HTC", "htc")

       !Scalars on hexagon centers (nodes) 
       cellopt="HX"
       div_ex%pos=0
       div_ex%n=mesh%nv
       div_ex%name="div_ex"
       allocate(div_ex%f(1:div_ex%n))

       align%n=mesh%nv
       align%name="alignind"
       align%pos=0
       allocate(align%f(1:mesh%nv))

       do i=1, mesh%nv
          !Divergence on point
          p=mesh%v(i)%p
          div_ex%f(i)=div_exact(p)

          !Divergence on barycenter
          p=mesh%hx(i)%b%p
          div_ex%f(i)=div_exact(p)

          !Alignment index
          align%f(i)=alignind(i,mesh)
       end do

       !Vector positions
       select case(trim(stag))
       case("HA", "ha")
          !Vectors on nodes
          v_ex%pos=0
          v_ex%n=mesh%nv
          allocate(v_ex%p(1:v_ex%n))
          do i=1, mesh%nv
             v_ex%p(i)%v=vecfield(p)
          end do
       case("HB", "hb")
          !Vectors on HX vertices
          v_ex%pos=2
          v_ex%n=mesh%nt
          allocate(v_ex%p(1:v_ex%n))
          do i=1, mesh%nt
             p=mesh%tr(i)%c%p
             v_ex%p(i)%v=vecfield(p)
          end do
       case("HC", "hc")
          !Vectors on HX edges
          v_ex%pos=3
          v_ex%n=mesh%ne
          allocate(v_ex%p(1:v_ex%n))
          do i=1, mesh%ne
             p=mesh%edhx(i)%c%p
             v_ex%p(i)%v=vecfield(p)
          end do
       case("HTC", "htc")
          !Vectors on HX/TR edge intersec
          v_ex%pos=6
          v_ex%n=mesh%ne
          allocate(v_ex%p(1:v_ex%n))
          do i=1, mesh%ne
             p=mesh%ed(i)%c%p
             v_ex%p(i)%v=vecfield(p)
          end do
       end select

    case("TA", "ta","TC", "tc")
       cellopt="TR"
       !Scalars on Tr circumcenters or barycenters
       div_ex%pos=1
       div_ex%n=mesh%nt
       allocate(div_ex%f(1:div_ex%n))

       !Set vector fields
       select case(trim(stag))
       case("TA", "ta")
          !Vector on TR centers
          v_ex%pos=1
          v_ex%n=mesh%nv
          allocate(v_ex%p(1:v_ex%n))
          do i=1, mesh%nt
             p=mesh%tr(i)%c%p
             div_ex%f(i)=div_exact(p)
             v_ex%p(i)%v=vecfield(p)
             !On barycenter
             p=mesh%tr(i)%b%p
             div_ex%f(i)=div_exact(p)
          end do
       case("TC", "tc")
          !Vectors on TR edges 
          v_ex%pos=2
          v_ex%n=mesh%ne
          allocate(v_ex%p(1:v_ex%n))
          do i=1, mesh%ne
             p=mesh%ed(i)%c%p
             v_ex%p(i)%v=vecfield(p)
          end do
       end select
    case default 
       print*, "DIVERGENCE TESTS error: unknown staggering", trim(stag)
       stop
    end select

    !Calculate divergence for the whole mesh
    average=.true.

    if(discinterp_mtd%hyb)then
       print*, "Using higher order method"
       div_est=divho_mesh_fullvec(v_ex, mesh, discinterp_mtd, cellopt, average)
    else
       print*, "Using usual method"
       div_est=div_mesh_fullvec(v_ex, mesh, cellopt, average)
    end if

    !Calculate Errors
    error=div_est
    error%name=trim(simulname)//"_error"
    error%f=div_ex%f-div_est%f

    !Global Errors (comment if you want only the aligned pol errors)
    error2=error_norm_2(div_est%f, div_ex%f, error%n)
    errormax=error_norm_max(div_est%f, div_ex%f, error%n)
    print*, "Error (max, L2): ", errormax, error2
    print*

    !Calculate error for aligned polygons (set .true. to use it)
    countalin=mesh%nv
    alignedpol=.false.
    if(trim(stag)=="HC".and. alignedpol )then

       countalin=0
       !Test using % of most aligned/nonaligned cells
       !allocate(alignsort(1:mesh%nv))
       !alignsort=align%f
       !call Qsort(alignsort(1:mesh%nv))

       !j=95*mesh%nv/100
       !j=mesh%nv-72
       !if(j<0)then
       !		j=1
       !end if
       !print*, "Using ", j, " of ", mesh%nv, " ( ", real(j)/real(mesh%nv)*100, " % ) "
       !limit=  alignsort(j) !mesh%meanvdist*mesh%meanvdist/10

       !Test using alignement index
       limit=  1._r8/100._r8
       !limit=  2._r8/3._r8
       print*, "Limit = ", limit
       do i=1, error%n
          if(dabs(align%f(i)) >= limit )then
             align%f(i)=0
          else
             align%f(i)=error%f(i)
             countalin=countalin+1
          end if
       end do

       !Test an especific region
       !tlat=-datan (real(0.5,8))
       !tlon=-2._r8*pi/5._r8
       !call sph2cart ( tlon, tlat, p(1), p(2), p(3))
       !i=getnearnode(p,mesh)
       !print*, "Node:", i, mesh%v(i)%nnb
       !do j=1, mesh%v(i)%nnb
       !	  k=mesh%v(i)%nb(j)
       !	  align%f(k)=error%f(k)
       !	  count=count+1
       !end do

       error2=dsqrt(dot_product(align%f,align%f)/real(countalin,r8))
       errormax=maxval(abs(align%f(1:mesh%nv)))
       print*, "Aligned count:", countalin, real(countalin)/real(mesh%nv)*100, "%"
       print*, "Aligned Error (max, L2): ", errormax, error2
    end if

    ! Save error estimates
    !-------------------------------------------------------

    !For hexagonal methods
    filename=trim(datadir)//"div_errors_"//trim(mesh%kind)// &
         "_"//trim(mesh%pos)//"_"//trim(mesh%optm)//".txt"
    call getunit(iunit)   

    inquire(file=filename, exist=ifile)
    if(ifile)then
       open(iunit,file=filename, status='old', position='append')
    else
       open(iunit,file=filename, status='replace')
       write(iunit, '(a)') &
            " n    distance stag testfunc errormax error2 alinged? numaligned"
    end if

    write(iunit, '(i8, f18.8, a12, i8, 2f22.12, l3, i10)') mesh%nv, mesh%meanvdist*rad2deg, &
         trim(stag), testfunc, errormax, error2, alignedpol, countalin

    close(iunit)

    ! Save point error estimates or table of error vs align
    !-------------------------------------------------------

    !For hexagonal methods
    !filename=trim(datadir)//"divalign_errors.txt"
    !call getunit(iunit)   

    !inquire(file=filename, exist=ifile)
    !ifile=.false.
    !if(ifile)then
    !   open(iunit,file=filename, status='old', position='append')
    !else
    !   open(iunit,file=filename, status='replace')
    !   !write(iunit, '(a)') &
    !   !     " n    distance  stag testfunc node error"
    !   write(iunit, '(a)') &
    !        " pt    align error"
    !end if

    !Choose node to print error
    !tlon= 0*deg2rad
    !tlat=-90*deg2rad
    !tlon=  -56.678556*deg2rad
    !tlat= 1.16503931150*deg2rad

    !call sph2cart(tlon, tlat, p(1), p(2), p(3))

    !i=getnearnode(p, mesh)
    !errormax=error%f(i)
    !write(iunit, '(i8, f18.8, a12, i8, i8, f22.12)') mesh%nv, mesh%meanvdist*rad2deg, &
    !     trim(stag), testfunc, i, errormax

    !do i=1, mesh%nv
    !  write(iunit, '(i8, 2f18.8)') i, align%f(i), error%f(i)
    !end do

    close(iunit)

    ! Plot fields 
    !-------------------------
    if(plots) then
       print*, "Plotting variables ... "

       div_ex%name=trim(simulname)//"_exact"
       call plot_scalarfield(div_ex, mesh)

       div_est%name=trim(simulname)//"_est"
       call plot_scalarfield(div_est, mesh)

       !error%f=error%f/(mesh%meanvdist**2)
       !error%f=abs(error%f)
       error%name=trim(simulname)//"_error"
       call plot_scalarfield(error, mesh)

       align%name=trim(simulname)//"_align"
       call plot_scalarfield(align, mesh)

       v_ex%name=trim(simulname)//"_vec"
       call plot_cart_vectorfield(v_ex, mesh)

    end if

    return

  end subroutine divergence_test

  !======================================================================
  !    LAPLACIAN TESTS
  !======================================================================
  subroutine laplacian_test(mesh)

    !Mesh
    type(grid_structure) :: mesh


    !Scalar fields
    type(scalar_field):: func
    type(scalar_field):: lap_ex
    type(scalar_field):: lap
    type(scalar_field):: error
    type(scalar_field):: ergrad
    type(scalar_field):: eddif

    !Counters
    integer:: i
    integer:: j
    integer:: k

    !Errors
    logical:: ifile
    integer:: iunit
    real (r8)::  error2
    real (r8):: errormax
    real (r8):: errorgrad2
    real (r8):: errorgradsup
    real (r8):: maxeddif
    character (len=256):: filename

    !Aux
    real(r8):: p(1:3)
    real(r8):: gradex
    real(r8):: gradest
    real(r8):: len


    print*
    print*,"Laplacian Testing "
    print*

    !-------------------------------------------
    !  Read parameters from file "simul.par"
    !------------------------------------------

    call getsimulpars(mesh)


    print*,"Setting up variables ..."

    !Scalars on hexagon centers (nodes) - Exact
    lap_ex%pos=0
    lap_ex%n=mesh%nv
    lap_ex%name="lap_ex"
    allocate(lap_ex%f(1:lap_ex%n))
    do i=1, mesh%nv
       !Laplacian on node
       p=mesh%v(i)%p
       !lap_ex%f(i)=lap_exact(p)

       !Laplacian on barycenter
       !p=mesh%hx(i)%b%p
       lap_ex%f(i)=lap_exact(p)
       !print "(i8, 4f16.8)", i, lap_ex%f(i), p
    end do

    !Edge displacement
    eddif%n=mesh%nv
    eddif%name="edgedif"
    eddif%pos=0
    allocate(eddif%f(1:mesh%nv))
    do i=1, mesh%nv
       eddif%f(i)=0._r8
       do j=1,mesh%v(i)%nnb
          k=mesh%v(i)%ed(j)
          eddif%f(i)=max(eddif%f(i), (arclen(mesh%ed(k)%c%p, mesh%edhx(k)%c%p)/ &
               mesh%edhx(k)%leng))
       end do
       !eddif%f(i)=eddif%f(i)/mesh%v(i)%nnb
    end do

    !Edge intersection maximum difference
    maxeddif=maxval(abs(eddif%f(1:eddif%n)))

    !Gradient error
    ergrad%n=mesh%nv
    ergrad%name="ergrad"
    ergrad%pos=0
    allocate(ergrad%f(1:mesh%nv))
    do i=1, mesh%nv
       ergrad%f(i)=0._r8
    end do

    !Laplacian error
    error%n=mesh%nv
    error%name="erlap"
    error%pos=0
    allocate(error%f(1:mesh%nv))

    !Test function - Scalars on hexagon centers (nodes)
    func%pos=0
    func%n=mesh%nv
    func%name="func"
    allocate(func%f(1:func%n))
    do i=1, mesh%nv
       !Function on node
       p=mesh%v(i)%p
       func%f(i)=f(p)
    end do

    !Calculate Numerical Laplacian
    lap%pos=0
    lap%n=mesh%nv
    lap%name="lap_est"
    allocate(lap%f(1:lap%n))
    do i=1, mesh%nv
       lap%f(i)=0._r8
       ergrad%f(i)=0._r8
       do j=1, mesh%v(i)%nnb
          !Edge index
          k=mesh%v(i)%ed(j)
          !Hexagonal edge midpoint
          p=mesh%edhx(k)%c%p
          !ExactGrad=ExactGradVector*NormalVectorEdge*CorrectionOutHx
          gradex=dot_product(df(p),mesh%edhx(k)%nr)*mesh%hx(i)%nr(j)
          !Estimated Gradient Normal component
          gradest=func%f(mesh%v(i)%nb(j))-func%f(i)
          len=mesh%ed(k)%leng
          !arclen(p, mesh%v(i)%p) + &
          !arclen(p, mesh%v(mesh%v(i)%nb(j))%p)
          gradest=gradest/ len !mesh%ed(k)%leng
          !Maximum gradient error for cell
          ergrad%f(i)=max(abs(gradest-gradex), ergrad%f(i))
          !Updade Laplacian
          lap%f(i)=lap%f(i)+gradest*mesh%edhx(k)%leng
       end do
       lap%f(i)=lap%f(i)/mesh%hx(i)%areag
       !Error in laplacian
       error%f(i)=lap_ex%f(i)-lap%f(i)
       !print"(i4, 3f16.8)", i, lap%f(i), lap_ex%f(i), error%f(i)
    end do

    !Global Errors

    error2=error_norm_2(lap%f, lap_ex%f, error%n)
    errormax=error_norm_max(lap%f, lap_ex%f, error%n)
    print*, "Error Lap (max, L2): ", errormax, error2

    !Global Errors for gradients
    errorgrad2=dsqrt(dot_product(ergrad%f,ergrad%f)/ergrad%n)
    errorgradsup=maxval(abs(ergrad%f(1:ergrad%n)))
    print*, "Error GRAD (max, L2): ", errorgradsup, errorgrad2
    print*

    ! Save error estimates
    !-------------------------------------------------------

    !For hexagonal methods
    filename=trim(datadir)//"lap_errors_"//trim(mesh%kind)// &
         "_"//trim(mesh%pos)//"_"//trim(mesh%optm)//".txt"
    call getunit(iunit)

    inquire(file=filename, exist=ifile)
    if(ifile)then
       open(iunit,file=filename, status='old', position='append')
    else
       open(iunit,file=filename, status='replace')
       write(iunit, '(a)') &
            " n    distance stag testfunc errormax error2 errorgradsup errorgrad2 maxeddif"
    end if

    write(iunit, '(i8, f18.8, a12, i8, 5f22.12)') mesh%nv, mesh%meanvdist*rad2deg, &
         trim(stag), testfunc, errormax, error2, errorgradsup, errorgrad2, maxeddif

    close(iunit)

    ! Plot fields
    !-------------------------
    if(plots) then
       print*, "Plotting variables ... "

       lap_ex%name=trim(simulname)//"_exact"
       call plot_scalarfield(lap_ex, mesh)

       lap%name=trim(simulname)//"_est"
       call plot_scalarfield(lap, mesh)

       error%name=trim(simulname)//"_error"
       call plot_scalarfield(error, mesh)

       eddif%name=trim(simulname)//"_eddif"
       call plot_scalarfield(eddif, mesh)

       ergrad%name=trim(simulname)//"_ergrad"
       call plot_scalarfield(ergrad, mesh)

    end if

    return

  end subroutine laplacian_test

 subroutine rotational_test(mesh)
    !-------------------------------------------------------------------
    ! rotational_test
    !
    ! Test rotational discretization
    !--------------------------------------------------------------------

    !Mesh
    type(grid_structure) :: mesh

    !Scalar fields
    type(scalar_field):: rot_ex
    type(scalar_field):: rot_est
    type(scalar_field):: error

    !Velocity field
    type(vector_field_cart):: v_ex

    !Counters
    integer:: i, j, k

    !Errors
    real (r8):: error2
    real (r8):: errormax
    real (r8):: signcor


    !File variables
    character (len=256):: filename
    logical:: ifile
    integer(i4):: iunit

    !Aux
    real(r8):: p(1:3)
    !character (len=2):: cellopt
    !logical:: average
    !logical:: alignedpol

    print*
    print*,"Rotational Testing "
    print*

    !-------------------------------------------
    !  Read parameters from file "simul.par"
    !   and set default names for outputs
    !------------------------------------------

    call getsimulpars(mesh)


    print*,"Testing rotational discretization ..."

    ! Exact/analitical Values
    !-------------------------------------

    !Allocate variables

    select case(trim(stag))
    case("HC", "hc")
       !Rotaional on triangle center
       ! data given on hexagon egdes
       !cellopt="HX"
       rot_ex%pos=1
       rot_ex%n=mesh%nt
       rot_ex%name="rot_ex"
       allocate(rot_ex%f(1:rot_ex%n))

       do i=1, mesh%nt
          !Divergence on point
          !p=mesh%tr(i)%c%p
          p=mesh%tr(i)%b%p
          rot_ex%f(i)=rot_exact(p)

       end do

        !Vectors on HX edges
          v_ex%pos=3
          v_ex%n=mesh%ne
          allocate(v_ex%p(1:v_ex%n))
          do i=1, mesh%ne
             p=mesh%edhx(i)%c%p
             v_ex%p(i)%v=vecfield(p) !df(p)
          end do

       rot_est%pos=1
       rot_est%n=mesh%nt
       rot_est%name="rot_est"
       allocate(rot_est%f(1:rot_est%n))

    case default
       print*, "Rotational tests error: unknown staggering", trim(stag)
       stop
    end select

    !Calculate rotational for the whole mesh

    !For all edges forming the triangle
    do i=1, mesh%nt
     rot_est%f(i)=0._r8
       do j=1, 3
          !Get edge index
          k=mesh%tr(i)%ed(j)
          !print*, "Edge:", j
          !Get edge outer normal related to the hexagon
          !signcor=real(mesh%tr(i)%nr(j), r8)

          signcor=dsign( 1._r8, real(mesh%tr(i)%tg(j), r8)* &
            dot_product(mesh%edhx(k)%nr, &
            mesh%ed(k)%tg))

          !Calculate numerical integration
          rot_est%f(i)=rot_est%f(i)+ &
              signcor*dot_product(v_ex%p(k)%v, mesh%edhx(k)%nr)*mesh%ed(k)%leng
       end do

          !Averaged divergence (division by cell area)
          rot_est%f(i)=rot_est%f(i)/mesh%tr(i)%areag

    end do

    !Calculate Errors
    error=rot_est
    error%name=trim(simulname)//"_error"
    error%f=rot_ex%f-rot_est%f

    !Global Errors (comment if you want only the aligned pol errors)
    error2=error_norm_2(rot_est%f, rot_ex%f, error%n)
    errormax=error_norm_max(rot_est%f, rot_ex%f, error%n)
    print*, "Error (max, L2): ", errormax, error2
    print*


    ! Save error estimates
    !-------------------------------------------------------

    !For hexagonal methods
    filename=trim(datadir)//"rot_errors_"//trim(mesh%kind)// &
         "_"//trim(mesh%pos)//"_"//trim(mesh%optm)//".txt"
    call getunit(iunit)

    inquire(file=filename, exist=ifile)
    if(ifile)then
       open(iunit,file=filename, status='old', position='append')
    else
       open(iunit,file=filename, status='replace')
       write(iunit, '(a)') &
            " n    distance stag testfunc errormax error2 "
    end if

    write(iunit, '(i8, f18.8, a12, i8, 2f22.12)') mesh%nv, mesh%meanvdist*rad2deg, &
         trim(stag), testfunc, errormax, error2

    close(iunit)


    ! Plot fields
    !-------------------------
    if(plots) then
       print*, "Plotting variables ... "

       rot_ex%name=trim(simulname)//"_exact"
       call plot_scalarfield(rot_ex, mesh)

       rot_est%name=trim(simulname)//"_est"
       call plot_scalarfield(rot_est, mesh)

       !error%f=error%f/(mesh%meanvdist**2)
       !error%f=abs(error%f)
       error%name=trim(simulname)//"_error"
       call plot_scalarfield(error, mesh)

       v_ex%name=trim(simulname)//"_vec"
       call plot_cart_vectorfield(v_ex, mesh)

    end if

    return

  end subroutine rotational_test

  !=====================================================================================
  !    TEST INTERPOLATION ROUTINES
  !=====================================================================================

  subroutine scalar_interpolation_test(mesh)
    !-----------------------------------------------
    !  TEST_SCALARINTERPOLATION
    !   Test scalar interpolation routines
    !-----------------------------------------------

    !Mesh structure
    type(grid_structure) :: mesh

    !Variable for interpolation
    type(scalar_field):: fdata   !Original data
    type(scalar_field):: frmap   !Remaped variable
    type(scalar_field):: frmaper !Remapping error

    !Gradient variables
    type(scalar_field):: gradnerr !Norm of errors
    type(vector_field_cart):: gradext   !Exat gradients
    type(vector_field_cart):: graderr   !Vector errors

    !Regular lon-lat grid variables
    real (r8), allocatable :: fexact(:,:) !Exact values
    real (r8), allocatable :: fest(:,:)   !Estimated values
    real (r8), allocatable :: ferror(:,:) !Error values

    !Indexes
    integer (i4):: i
    integer (i4):: j
    !integer (i4):: n
    integer (i4):: nlon !Total number of longitudes
    integer (i4):: nlat !Total number of latitudes

    !Files variables
    logical:: ifile                !Logical for file existence
    integer:: errorsunit           !Unit for error files
    character (len=256):: filename !Name for files

    !Errors
    real (r8):: error2       !Quadratic errors
    real (r8):: errormax     !Maximum errors
    real (r8):: graderror2   !Quad error for gradients
    real (r8):: graderrormax !Max error for gradients
    real (r8):: condnummin   !Min condition number estimate
    real (r8):: condnummax   !Max condition number estimate

    !Auxiliar variables
    real (r8):: p(1:3)  !General point in R3
    real (r8):: dlat    !Latitude step
    real (r8):: dlon    !Longitude step
    real (r8):: tlat    !Latitude
    real (r8):: tlon    !Longitude
    real (r8):: tlatmax !Latitude of max error
    real (r8):: tlonmax !Longitude of max error
    real (r8):: intval  !Interpolated value
    real (r8):: exval   !Exact value
    real (r8):: mass    ! Global mass of variable
    real (r8):: mass_remap! Global mass of variable after remaping

    !RBF stencil and kind of interpolation
    character(len=4):: stencil

    !RBF matrix vector
    type(rbf_matrix_structure), allocatable :: rbf_mat(:)

    !Time counting variables
    real(r8):: elapsed_time
    real(r8):: start
    real(r8):: finish

    print*
    print*,"Test scalar interpolation functions "
    print*,"------------------------------------"
    print*

    !-------------------------------------------
    !  Read parameters from file "simul.par"
    !------------------------------------------
    !Read variable positions on mesh
    ! interpolation method and test function
    call getsimulpars(mesh)

    !-----------------------------------------
    ! Define variable on mesh
    !----------------------------------------
    select case (trim(stag))
    case ("HA") !Scalars on triangle vertices
       !Data values
       fdata%pos=0
       fdata%n=mesh%nv
       allocate(fdata%f(1:fdata%n))

       !Remap to Triangle Circumcenters
       frmap%pos=1
       frmap%n=mesh%nt
       allocate(frmap%f(1:frmap%n))

       if(testgrad)then
          !Vector variables
          gradext%pos=0
          gradext%n=mesh%nv
          allocate(gradext%p(1:gradext%n))

          graderr%pos=0
          graderr%n=mesh%nv
          allocate(graderr%p(1:graderr%n))

          gradnerr%pos=0
          gradnerr%n=mesh%nv
          allocate(gradnerr%f(1:fdata%n))
       end if

       !Data and gradient value calculations
       mass=0.
       do i=1,mesh%nv
          p=mesh%v(i)%p
          fdata%f(i)=f(p)
          if(testgrad)then
             gradext%p(i)%v(1:3)=proj_vec_sphere(df(p), p)
          end if
          mass=mass+mesh%hx(i)%areag*fdata%f(i)
       end do

    case ("TA") !Scalars on triangle circumcenters
       !Data values
       fdata%pos=1
       fdata%n=mesh%nt
       allocate(fdata%f(1:fdata%n))

       !Remap to Triangle Vertices
       frmap%pos=0
       frmap%n=mesh%nv
       allocate(frmap%f(1:frmap%n))

       !Data value calculation
       do i=1,mesh%nt
          p=mesh%tr(i)%c%p
          fdata%f(i)=f(p)
       end do

    case ("TC") !Supose scalars on triangle edges
       !Data values
       fdata%pos=2
       fdata%n=mesh%ne
       allocate(fdata%f(1:fdata%n))

       !Remap to Triangle Circumcenters
       frmap%pos=1
       frmap%n=mesh%nt
       allocate(frmap%f(1:frmap%n))

       !Data value calculation
       do i=1,mesh%ne
          p=mesh%ed(i)%c%p
          fdata%f(i)=f(p)
       end do

    case ("HC") !Supose scalars on Hexagon edges
       !Data values
       fdata%pos=3
       fdata%n=mesh%ne
       allocate(fdata%f(1:fdata%n))

       !Remap to Triangle Vertices
       frmap%pos=0
       frmap%n=mesh%nv
       allocate(frmap%f(1:frmap%n))

       !Data value calculation
       do i=1,mesh%ne
          p=mesh%edhx(i)%c%p
          fdata%f(i)=f(p)
       end do
    end select

    !Create error variable for remapping
    frmaper=frmap

    !If RBF, calculate the matrix structure (precalculation, depends only on mesh)
    condnummin=100000000.
    condnummax=0.
    if(kinterp(1:3)=="RBF" .or. kinterp(1:3)=="rbf")then
       stencil=trim(kinterp(4:7))
       call rbf_matrix_build( stencil, rbf_par, rbf_mat, mesh )
       !Calculate condition number estimates
       print*
       print*, "Length of RBF matrices vector:", size(rbf_mat)
       do i=1,size(rbf_mat)
          condnummin=min(condnummin,rbf_mat(i)%condnum)
          condnummax=max(condnummax,rbf_mat(i)%condnum)
       end do
       print*, "Condition Number Estimates (min, max):", condnummin, condnummax
       print*
    end if

    !----------------------------------------
    ! Test  interpolation for a point
    !----------------------------------------
    print*
    print*,"Test point interpolation"
    print*

    !Point
    !tlon=-64.5
    !tlat=-13.5
    tlon=-2.5
    tlat= -2.5
    tlon=mesh%tr(min(154,mesh%nt) )%c%lon * rad2deg
    tlon=mesh%tr(min(154,mesh%nt) )%c%lat * rad2deg
    print*, " Point interpolation   (Lon,Lat): ", &
         real(tlon,4), real(tlat,4)

    call sph2cart(tlon*deg2rad, tlat*deg2rad, p(1), p(2), p(3))

    intval=scalar_interpol(p, fdata, mesh, kinterp, rbf_mat)

    exval=f(p)
    print *, " Exact  : ", exval
    print *, " Interp : ", intval
    print *, " Error  : ", exval-intval
    print*

    !-----------------------------------
    !Test remmaping
    !----------------------------------

    !Zero error variables
    error2=0.
    errormax=0.

    !start counting time
    call cpu_time(start)
    mass_remap=0.
    !Loop over remapping points
    do i=1, frmap%n
       if(frmap%pos==1)then !Remap to TR circumcenters
          p=mesh%tr(i)%c%p
       elseif(frmap%pos==0)then !Remap to TR Vertices
          p=mesh%v(i)%p
       else
          print*, "Scalar_interpolation_test error: Don't "//&
               "know how to do this remap - ", frmap%pos
          stop
       end if

       !Interpolate
       frmap%f(i)=scalar_interpol(p, fdata, mesh, kinterp, rbf_mat)
        mass_remap=mass_remap+frmap%f(i)*mesh%tr(i)%areag
       !Calculate error
       frmaper%f(i)=f(p)-frmap%f(i)
       error2=error2+(frmaper%f(i))**2
       errormax=max(errormax, abs(frmaper%f(i)) )
    end do

    !Stop counting time
    call cpu_time(finish)
    elapsed_time= finish-start !real(clock_end-clock_start, r8)/real(clock_rate,r8)

    !Calculate quadratic error
    error2= dsqrt(error2/real(frmap%n, r8))

    !Print Errors
    print*, "Remapping on grid error  (inf/ L2) : ", real(errormax,4), real(error2,4)
    print*, "Time : ", elapsed_time
    print*, "Mass (before, after remap)", mass, mass_remap, abs(mass-mass_remap)
    print*
    !File for errors
    filename=trim(datadir)//"scinterp_rmap_errors.txt"
    call getunit(errorsunit)

    !Check if file exists, to append info
    inquire(file=filename, exist=ifile)
    if(ifile)then
       open(errorsunit,file=filename, status='old', position='append')
    else
       open(errorsunit,file=filename, status='replace')
       write(errorsunit, '(a)') &
            "          n                meanvdist              "//&
            " meshname             kinterp        stag     func      interperror(inf)"//&
            "      interperror(L2)       time     "
    end if

    ! Write interpolation erros
    write(errorsunit, '(i12, f22.8, a32, 2a12, i8, 7e22.8)') mesh%nv,  &
         mesh%meanvdist*rad2deg, trim(mesh%name), trim(kinterp), &
         trim(stag), testfunc, errormax, error2, elapsed_time

    !Close file
    close(errorsunit)

    !Plot fields
    if(plots)then
       print*,"Writing grid remapping data on files for plotting"

       frmap%name=trim(simulname)//"_rmap"
       call plot_scalarfield(frmap, mesh)

       frmaper%name=trim(simulname)//"_rmaper"
       call plot_scalarfield(frmaper, mesh)
       print*
    end if

    !-----------------------------------------
    !Calculate errors in gradient estimatives
    !----------------------------------------
    graderror2=0.
    graderrormax=0.
    if(testgrad)then
       !Calculate gradients
       call gradcalc(fdata, mesh)

       !Calculate errors
       do i=1, fdata%n
          graderr%p(i)%v=fdata%g(i)%v-gradext%p(i)%v
          gradnerr%f(i)=norm(fdata%g(i)%v-gradext%p(i)%v)
          graderror2=graderror2+(gradnerr%f(i))**2
          graderrormax=max(graderrormax, gradnerr%f(i))
       end do
       graderror2=dsqrt(graderror2/fdata%n)

       print*, "GRAD errors: inf/ L2", graderrormax, graderror2

       !Plot fields
       if(plots)then
          print*
          print*,"Writing gradient data on files for plotting"

          fdata%name=trim(simulname)
          call plot_grad_vectorfield(fdata, mesh)

          gradext%name=trim(simulname)//"_gradexact"
          call plot_cart_vectorfield(gradext, mesh)

          graderr%name=trim(simulname)//"_graderror"
          call plot_cart_vectorfield(graderr, mesh)

          gradnerr%name=trim(simulname)//"_gradnerror"
          call plot_scalarfield(gradnerr, mesh, kinterp)
          print*
       end if
    end if

    !-----------------------------------------
    ! Do remmaping to lat-lon grid and evaluate errors
    !  for each interpolation type
    !----------------------------------------

    !Lat-lon grid sizes
    nlat=720 !180 !720
    nlon=1440 !360 !1440
    dlat=180._r8/real(nlat, r8)
    dlon=360._r8/real(nlon, r8)

    !Allocate space
    allocate(fexact(1:nlon,1:nlat))
    allocate(fest(1:nlon,1:nlat))
    allocate(ferror(1:nlon,1:nlat))

    !Zero variables
    error2=0.
    errormax=0.

    !Do a remapping to lat lon grid

    !Start counting time
    call cpu_time(start)

    !Pixel registration mode (for GMT graph)
    tlat=-90._r8+dlat/2._r8
    lat: do j=1,nlat
       tlon=-180._r8+dlon/2._r8
       lon: do i=1,nlon

          !Convert to cartesian coords
          call sph2cart(tlon*deg2rad,tlat*deg2rad, p(1), p(2), p(3))

          !Interpolate
          fest(i,j)=scalar_interpol(p, fdata, mesh, kinterp, rbf_mat)

          !Calculate errors
          fexact(i,j)=f(p)
          ferror(i,j)=fest(i,j) - fexact(i,j)
          error2=error2+(ferror(i,j))**2 
          !print*, tlon, tlat, ferror(i,j)

          !Set the point of maximum error (for debuging purpuses)
          if(errormax < abs(ferror(i,j)))then
             tlonmax=tlon
             tlatmax=tlat
          end if
          errormax=max(errormax, abs(ferror(i,j)))

          tlon=tlon+dlon
       end do lon
       tlat=tlat+dlat
    end do lat

    !Stop counting time
    call cpu_time(finish)
    elapsed_time= finish-start

    !Calculate quadratic error
    error2= dsqrt(error2/real(nlon*nlat, r8))

    !Print Errors
    print *, "Remap to lat/lon error (inf/ L2): ", real(errormax,4), real(error2,4)
    print *, "  Largest error in point (lon, lat):", real(tlonmax, 4), real(tlatmax, 4)
    print *, "  Time: ", elapsed_time
    print *

    !File for errors
    filename=trim(datadir)//"scinterp_rmapll_errors.txt"
    call getunit(errorsunit)

    !Check if file exists to append info
    inquire(file=filename, exist=ifile)
    if(ifile)then
       open(errorsunit,file=filename, status='old', position='append')
    else
       open(errorsunit,file=filename, status='replace')
       write(errorsunit, '(a)') &
            "          n                meanvdist              "//&
            " meshname             kinterp        stag     func      interperror(inf)"//&
            "      interperror(L2)       time                     graderror(inf) "//&
            "       graderror(L2)        condnummin            condnummax"
    end if

    ! Write interpolation erros
    write(errorsunit, '(i12, f22.8, a32, 2a12, i8, 7e22.8)') mesh%nv,  &
         mesh%meanvdist*rad2deg, trim(mesh%name), trim(kinterp), &
         trim(stag), testfunc, errormax, error2, elapsed_time, &
         graderrormax, graderror2, condnummin, condnummax

    !Plot interpolation maps
    if(plots)then
       print*,"Writing lont-lat remapping data on files for plotting"
       filename=trim(simulname)//"_"//trim(mesh%name)
       call plot_scalarfield_sphgrid(nlon, nlat, fest, filename)

       filename=trim(simulname)//"_error_"//trim(mesh%name)
       call plot_scalarfield_sphgrid(nlon, nlat, ferror, filename)

       filename=trim(simulname)//"_exact_"//trim(mesh%name)
       call plot_scalarfield_sphgrid(nlon, nlat, fexact, filename)
    end if

    print*
    print*,"Run plot.sh in gmt/ to see results"
    print*

    return
  end subroutine scalar_interpolation_test

  subroutine vector_interpolation_test(mesh)
    !-----------------------------------------------
    !  TEST_VECTORINTERPOLATION
    !   Test vector interpolation routines
    !-----------------------------------------------
    !Mesh structure
    type(grid_structure) :: mesh

    !Indexes
    integer (i4):: i
    integer (i4):: j
    integer (i4):: nlon
    integer (i4):: nlat

    !Units for files
    logical:: ifile

    integer:: errorsunit
    character (len=256):: filename

    !Errors
    real (r8)::  error2
    real (r8):: errormax

    !Interpolable vector variable
    type(vector_field_cart):: w

    !Regular grid estimates
    type(vector), allocatable :: fexact(:,:)
    type(vector), allocatable :: fest(:,:)
    type(vector), allocatable :: ferror(:,:)
    real (r8), allocatable :: ferrorabs(:,:)

    !Auxiliar variables
    real (r8), dimension(1:3):: p
    real (r8):: dlat
    real (r8):: dlon
    real (r8):: tlat
    real (r8):: tlon
    real (r8):: tlatmax
    real (r8):: tlonmax
    real (r8):: vecint(1:3)
    real (r8):: vecexact(1:3)

    !Time counting variables

    real(r8):: elapsed_time
    real(r8):: start
    real(r8):: finish

    print*
    print*,"Test vector interpolation functions "
    print*,"-----------------------------"
    print*

    !-------------------------------------------
    !  Read parameters from file "simul.par"
    !------------------------------------------
    call getsimulpars(mesh)

    !-----------------------------------------
    ! Define variable on mesh
    !----------------------------------------
    select case (trim(stag))
    case ("HA") !Vectors on triangle vertices
       w%pos=0
       w%n=mesh%nv
       allocate(w%p(1:w%n))
       !Exact vector field
       do i=1,mesh%nv
          p=mesh%v(i)%p
          w%p(i)%v=vecfield(p)
       end do
    case ("TC") !Vectors on triangle edges
       w%pos=2
       w%n=mesh%ne
       allocate(w%p(1:w%n))
       do i=1,mesh%ne
          p=mesh%ed(i)%c%p
          w%p(i)%v=vecfield(p)
       end do
    case ("HC") !Vectors on Hexagon edges
       w%pos=3
       w%n=mesh%ne
       allocate(w%p(1:w%n))
       !Exact function/gradient calculation
       do i=1,mesh%ne
          p=mesh%edhx(i)%c%p
          w%p(i)%v=vecfield(p)
       end do
    end select

    !-----------------------------------------
    ! Test  interpolation for a point
    !----------------------------------------
    !Point
    !tlon=-64.5
    !tlat=-13.5
    tlon=-36.0878869808453
    tlon=tlon-2.5
    tlat=58.2825557357298


    !tlon=mesh%tr(min(85,mesh%nt) )%c%lon * rad2deg
    !tlon=mesh%tr(min(85,mesh%nt) )%c%lat * rad2deg
    print*, " Point interpolation   (Lon,Lat): ", &
         real(tlon,4), real(tlat,4)
    call sph2cart(tlon*deg2rad, tlat*deg2rad, p(1), p(2), p(3))

    vecint=vector_interpol(p, w, mesh, kinterp)
    vecexact=vecfield(p)
    print *, "  Exact     :", real(vecexact,4)
    print *, "  Interpoled:", real(vecint, 4)
    print *, "  Error     :", real(norm(vecint-vecexact),4)
    print*

    !stop
    !-----------------------------------------
    ! Do remmaping to lat-lon grid and evaluate errors
    !  for each interpolation type
    !----------------------------------------

    !Lat-lon grid sizes
    nlat=720 !180
    nlon=1440 !360
    dlat=180._r8/real(nlat, r8)
    dlon=360._r8/real(nlon, r8)

    !Allocate space
    allocate(fexact(1:nlon,1:nlat))
    allocate(fest(1:nlon,1:nlat))
    allocate(ferror(1:nlon,1:nlat))
    allocate(ferrorabs(1:nlon,1:nlat))

    !File for errors
    filename=trim(datadir)//"vecinterp_errors.txt"
    call getunit(errorsunit)   

    inquire(file=filename, exist=ifile)
    if(ifile)then
       open(errorsunit,file=filename, status='old', position='append')
    else
       open(errorsunit,file=filename, status='replace')
       write(errorsunit, '(a)') "         n        meanvdist        kinterp   stag "//&
            "   func      interperror(inf)      interperror(L2)    time "
    end if

    !Zero variables
    !fest(1:nlon, 1:nlat)=0. 
    error2=0.
    errormax=0.

    !Do a remapping to lat lon grid
    call cpu_time(start)

    !Pixel registration mode (GMT)
    tlat=-90._r8+dlat/2._r8
    lat: do j=1,nlat
       tlon=-180._r8+dlon/2._r8
       lon: do i=1,nlon

          !Convert to cartesian coords
          call sph2cart(tlon*deg2rad,tlat*deg2rad, p(1), p(2), p(3))
          !Interpol
          fest(i,j)%v=vector_interpol(p, w, mesh, kinterp)
          !Calculate errors
          fexact(i,j)%v=vecfield(p)
          ferror(i,j)%v=fest(i,j)%v - fexact(i,j)%v
          ferrorabs(i,j)=norm(ferror(i,j)%v)

          error2=error2+(norm(ferror(i,j)%v))**2 
          !print*, tlon, tlat, ferror(i,j)
          !Set the point of maximum error (for debuging purpuses)
          if(errormax < norm(ferror(i,j)%v))then
             tlonmax=tlon
             tlatmax=tlat
          end if
          errormax=max(errormax, norm(ferror(i,j)%v))

          tlon=tlon+dlon
       end do lon
       tlat=tlat+dlat
    end do lat
    call cpu_time(finish)
    elapsed_time= finish-start !real(clock_end-clock_start, r8)/real(clock_rate,r8)

    error2= dsqrt(error2/real(nlon*nlat, r8))

    !Print Errors
    print*, "Interpol error (inf/ L2) :", errormax, error2
    print*, "  Point of largest error (lon, lat):", real(tlonmax,4), real(tlatmax,4)
    print*, "  Time: ", elapsed_time
    print*

    if(plots)then
       !Plot interpolation maps

       w%name=trim(simulname)//"_vec_"
       call plot_cart_vectorfield(w, mesh)

       filename=trim(simulname)//"_estim_"//trim(mesh%name)
       call plot_vectorfield_sphgrid(nlon, nlat, fest, filename)

       filename=trim(simulname)//"_error_"//trim(mesh%name)
       call plot_vectorfield_sphgrid(nlon, nlat, ferror, filename)

       filename=trim(simulname)//"_exact_"//trim(mesh%name)
       call plot_vectorfield_sphgrid(nlon, nlat, fexact, filename)

       filename=trim(simulname)//"_nerror_"//trim(mesh%name)
       call plot_scalarfield_sphgrid(nlon, nlat, ferrorabs, filename)
    end if

    ! Write interpolation erros
    write(errorsunit, '(i12, f22.8, 2a12, i8, 3e22.8)') mesh%nv,  &
         mesh%meanvdist*rad2deg, trim(kinterp), trim(stag), testfunc, errormax, error2, elapsed_time

    print*,"Run plot.sh to see results"
    print*
    return

  end subroutine vector_interpolation_test

  subroutine vector_reconstruction_test(mesh)
    !-----------------------------------------------
    !  TEST VECTOR RECONSTRUCTION
    !   Test vector reconstruction routines
    !-----------------------------------------------
    !Mesh structure
    type(grid_structure) :: mesh

    !Indexes
    integer (i4):: i
    integer (i4):: j
    integer (i4):: n
    integer (i4):: nlon
    integer (i4):: nlat

    !Units for files
    logical:: ifile
    integer:: errorsunit
    character (len=256):: filename

    !Errors
    real (r8):: error2
    real (r8):: errormax

    !Condition numbers
    real (r8):: condnummin
    real (r8):: condnummax

    !Interpolable vector variable
    type(vector_field_cart):: vecedfull   !Full vector field
    type(vector_field_cart):: vecnormal !Normal comp of vectors
    type(vector_field_cart):: vecrec    !Reconstructed vectors
    type(vector_field_cart):: vecerror  !Error vector
    type(vector_field_cart):: vecexact  !Exact vector field
    type(scalar_field):: vecnerror  !Magnitude of error field

    !Variable with normal components
    type(scalar_field):: vecncomp

    !Auxiliar variables
    real (r8), dimension(1:3):: p !Arbitrary R3 point
    real (r8):: dlat               !Latitude step
    real (r8):: dlon               !Longititude step
    real (r8):: tlat               !Latitutes
    real (r8):: tlon               !Longitudes
    real (r8):: tlatmax            !Latitude of max error
    real (r8):: tlonmax            !Longitude of max error

    integer (i4):: vecreconpos      !Position for vector reconstruction
    integer (i4):: vecpos           !Position for given vectors
    logical :: intisrecon           !Flag when interpol mtd is the recon mtd

    !Time counting variables
    real(r8):: time_remap
    real(r8):: time_latlon
    real(r8):: totaltime
    real(r8):: time_rbf
    real(r8):: start
    real(r8):: finish

    !Regular grid estimates
    type(vector), allocatable :: llexact(:,:)
    type(vector), allocatable :: llrecon(:,:)
    type(vector), allocatable :: llerror(:,:)
    real (r8),    allocatable :: llnerror(:,:)

    !RBF matrix vector
    type(rbf_matrix_structure), allocatable :: rbf_mat(:)

    !RBF stencil and kind of interpolation
    character(len=4):: stencil

    print*
    print*,"Test vector reconstruction"
    print*,"-----------------------------"
    print*

    !-------------------------------------------
    !  Read parameters from file "simul.par"
    !------------------------------------------
    !Read kinterp for kind of vector reconstruction to be used
    call getsimulpars(mesh)

    !-----------------------------------------
    ! Define variables on mesh
    !----------------------------------------
    select case (trim(stag))
    case ("TC") !Vectors on triangle edges
       vecpos=2
    case ("HC") !Vectors on Hexagon edges
       vecpos=3
    end select

    !Full vector
    vecedfull%pos=vecpos
    vecedfull%n=mesh%ne
    allocate(vecedfull%p(1:vecedfull%n))

    !Normal component as vector
    vecnormal%pos=vecpos
    vecnormal%n=mesh%ne
    allocate(vecnormal%p(1:vecnormal%n))

    !Normal component - scalar
    vecncomp%pos=vecpos
    vecncomp%n=mesh%ne
    allocate(vecncomp%f(1:vecncomp%n))

    select case (trim(stag))
    case ("TC") !Vectors on triangle edges
       !Calculate vector field and components
       do i=1,mesh%ne
          p=mesh%ed(i)%c%p
          vecedfull%p(i)%v=vecfield(p)
          vecncomp%f(i)=dot_product(vecedfull%p(i)%v, mesh%ed(i)%nr)
          vecnormal%p(i)%v=vecncomp%f(i)*mesh%ed(i)%nr
       end do
    case ("HC") !Vectors on Hexagon edges
       !Calculate vector field and components
       do i=1,mesh%ne
          p=mesh%edhx(i)%c%p
          vecedfull%p(i)%v=vecfield(p)
          vecncomp%f(i)=dot_product(vecedfull%p(i)%v, mesh%edhx(i)%nr)
          vecnormal%p(i)%v=vecncomp%f(i)*mesh%edhx(i)%nr
       end do
    end select

    !----------------------------------------
    !Alocate space for reconstructed vector field
    !----------------------------------------

    !Kind of reconstruction
    select case (trim(recon_mtd%recon))
    case ("rbfhx", "perhx", "perpj", "lsqhxe", "kls", "none", "nonehx", &
         "rbfhxpc")
       if(recon_mtd%massc)then
          !Reconstruct vectors to Voronoi cell centroid
          vecreconpos=4
       else
          !Reconstruct vectors to triangle vertices/hexagon centers
          vecreconpos=0
       end if
    case ("pertr", "lsqtrc", "rt0", "wht", "rbftr", "rbfetr", "nonetr", &
         "rbftrpc", "rbfetrpc", "rbfetrp")
       if(recon_mtd%massc)then
          !Reconstruct vector to triangle barycenters
          vecreconpos=5
       else
          !Reconstruct vector to triangle circumcenters /hexagon vertices
          vecreconpos=1
       end if
    case default
       print*, "Error on vector reconstruction test: unknown recon method ", trim(recon_mtd%recon)
       stop
       !  hexagon centers
       vecreconpos=0
    end select

    !Set number of vector reconstructions
    select case(vecreconpos)
    case(0, 4) !Reconstruction to triangle vertices/ Voronoi barycenters
       vecrec%n=mesh%nv
    case(1, 5) !Reconstruction to triangle circumcenters / barycenters
       vecrec%n=mesh%nt
    end select

    !Reconstructed vectors
    vecrec%pos=vecreconpos
    allocate(vecrec%p(1:vecrec%n))

    !Reconstruction error
    vecerror%pos=vecreconpos
    vecerror%n=vecrec%n
    allocate(vecerror%p(1:vecerror%n))

    !Exact vector field error
    vecexact%pos=vecreconpos
    vecexact%n=vecrec%n
    allocate(vecexact%p(1:vecexact%n))

    !Norm of the reconstruction error
    vecnerror%pos=vecreconpos
    vecnerror%n=vecrec%n
    allocate(vecnerror%f(1:vecnerror%n))

    !If RBF, calculate the matrix structure (precalculation, depends only on mesh)
    condnummin=100000000.
    condnummax=0.
    time_rbf=0.
    if(recon_mtd%recon(1:3)=="rbf")then
       stencil=trim(recon_mtd%recon(4:7))
       !Build rbf matrices
       call cpu_time(start)
       call rbfvec_matrix_build( stencil, rbf_par, rbf_mat, mesh)
       call cpu_time(finish)
       !time_remap = (finish-start)/5.
       time_rbf = (finish-start)
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
       print*, rbf_par
       print*
    end if

    !-----------------------------------------
    ! Test reconstruction
    !----------------------------------------
    print*
    print*,"Test reconstruction re-map on geodesic grid"
    print*

    !Reconstruct vector field
    call cpu_time(start)
    !do i=1,5
    call vrec_remap(vecncomp, vecrec, mesh, recon_mtd, rbf_mat)
    !end do
    call cpu_time(finish)
    !time_remap = (finish-start)/5.
    time_remap = (finish-start)
    !time_remap=time_remap-time_rbf

    !Calculate errors
    error2=0._r8
    errormax=0.0_r8
    do i=1,  vecrec%n
       !Do reconstruction based on vecreconpos
       select case(vecreconpos)
       case(0) !Reconstruc to tr vertices
          p=mesh%v(i)%p
       case(1) !Reconstruction to tr circuncenters
          p=mesh%tr(i)%c%p
       case(4) !Reconstruct to Voronoi cell centroid
          p=mesh%hx(i)%b%p
       case(5) !Reconstruc to tr barycenter
          p=mesh%tr(i)%b%p
       case default
          print*, "Error on vector reconstruction test: unknown recon position ", vecreconpos
          stop
       end select

       vecexact%p(i)%v=vecfield(p)
       vecerror%p(i)%v=vecexact%p(i)%v-vecrec%p(i)%v
       vecnerror%f(i)=norm(vecerror%p(i)%v)
       !print*, i, vecnerror%f(i)
       errormax=max(errormax, vecnerror%f(i))
       !errormax=max(errormax, dot_product(vecexact%p(i)%v,vecexact%p(i)%v)- &
       !         dot_product(vecrec%p(i)%v, vecrec%p(i)%v))
       error2=error2+vecnerror%f(i)**2
    end do
    !Error 2
    error2=dsqrt(error2/real(vecrec%n, r8))
    print *, "Max and L2 errors for remap: ", errormax, error2
    print*, "Elapsed time for remap: ", time_remap
    print*, "Elapsed time for rbf matrix build:", time_rbf

    !File for errors
    filename=trim(datadir)//"vreconremap_errors.txt"
    call getunit(errorsunit)

    inquire(file=filename, exist=ifile)
    if(ifile)then
       open(errorsunit,file=filename, status='old', position='append')
    else
       open(errorsunit,file=filename, status='replace')
       write(errorsunit, '(a)') "        n              meanvdist          "//&
            "          kinterp   stag    func          error(inf)          error(L2)"//&
            "             condnummin             condnummax            RBFshape "//&
            "        alignpercent             time                             mesh "
    end if
    ! Write errors
    write(errorsunit, '(i12, f22.8, a24, a8, i8, 7e22.8, a32)') mesh%nv,  &
         mesh%meanvdist*rad2deg, trim(recon_mtd%name), trim(stag), testfunc, errormax, error2, &
         condnummin, condnummax, rbf_par, recon_mtd%alignpercent, time_remap, trim(mesh%name)

    !-----------------------------------------
    ! Test reconstruction to lat-lon grid
    !----------------------------------------
    print*
    print*,"Test reconstruction re-map to lat-lon grid"
    print*

    !Lat-lon grid sizes
    nlat=180*4
    nlon=360*4
    dlat=180._r8/real(nlat, r8)
    dlon=360._r8/real(nlon, r8)

    !Allocate space
    allocate(llexact(1:nlon,1:nlat))
    allocate(llrecon(1:nlon,1:nlat))
    allocate(llerror(1:nlon,1:nlat))
    allocate(llnerror(1:nlon,1:nlat))

    !Zero variables
    !fest(1:nlon, 1:nlat)=0.
    error2=0.
    errormax=0.
    intisrecon=(trim(recon_mtd%interp)==trim(recon_mtd%recon))

    !Do a remapping to lat lon grid
    call cpu_time(start)

    !Pixel registration mode (GMT)
    tlat=-90._r8+dlat/2._r8
    do j=1,nlat
       tlon=-180._r8+dlon/2._r8
       do i=1,nlon

          !Convert to cartesian coords
          call sph2cart(tlon*deg2rad,tlat*deg2rad, p(1), p(2), p(3))

          !Interpol
          if(intisrecon)then
             llrecon(i,j)%v=vector_reconstruct(p, vecncomp, mesh, recon_mtd%interp, rbf_mat )
          else
             llrecon(i,j)%v=vector_interpol(p, vecrec, mesh, recon_mtd%interp)
          end if

          !Calculate errors
          llexact(i,j)%v=vecfield(p)
          llerror(i,j)%v=llrecon(i,j)%v - llexact(i,j)%v
          llnerror(i,j)=norm(llerror(i,j)%v)
          error2=error2+(norm(llerror(i,j)%v))**2
          !print*, tlon, tlat, ferror(i,j)

          !Set the point of maximum error (for debuging purpuses)
          if(errormax < norm(llerror(i,j)%v))then
             tlonmax=tlon
             tlatmax=tlat
          end if
          errormax=max(errormax, norm(llerror(i,j)%v))

          tlon=tlon+dlon
       end do
       tlat=tlat+dlat
    end do

    call cpu_time(finish)
    time_latlon =  finish-start

    !Error 2
    error2= dsqrt(error2/real(nlon*nlat, r8))

    !Print Errors
    print *, "Max and L2 errors: ",  errormax, error2
    print*, "  Point of largest error (lon, lat):", real(tlonmax,4), real(tlatmax,4)
    print*, "  Time for regrid to lat lon: ", time_latlon
    print*
    totaltime=time_remap+time_latlon
    print*, "  Total Time: ", totaltime
    print*

    !File for errors
    filename=trim(datadir)//"vreconll_errors.txt"
    call getunit(errorsunit)

    inquire(file=filename, exist=ifile)
    if(ifile)then
       open(errorsunit,file=filename, status='old', position='append')
    else
       open(errorsunit,file=filename, status='replace')
       write(errorsunit, '(a)') "         n                meanvdist "//&
            "                kinterp               stag    func        error(inf) "//&
            "             error(L2)           condnummin            condnummax    "//&
            "          rbf_shape              aligned%               time_latlon  "//&
            "          totaltime                       mesh"
    end if

    write(errorsunit, '(i12, f22.8, a32, a12, i8, 8e22.8, a32)') mesh%nv,  &
         mesh%meanvdist*rad2deg, trim(recon_mtd%name), trim(stag), testfunc, errormax, error2, &
         condnummin, condnummax, rbf_par, recon_mtd%alignpercent, time_latlon, totaltime, trim(mesh%name)

    if(plots)then

       !Plot vector remaps error
       vecnerror%name=trim(simulname)//"_vecnerror"
       call plot_scalarfield(vecnerror, mesh)

       !Plot vector field maps
       vecrec%name=trim(simulname)//"_vecrec"
       call plot_cart_vectorfield(vecrec, mesh)

       vecexact%name=trim(simulname)//"_vecexact"
       call plot_cart_vectorfield(vecexact, mesh)

       vecerror%name=trim(simulname)//"_vecerror"
       call plot_cart_vectorfield(vecerror, mesh)

       !vecedfull%name=trim(simulname)//"_vecfull"
       !call plot_cart_vectorfield(vecedfull, mesh)

       !vecnormal%name=trim(simulname)//"_vecnormal"
       !call plot_cart_vectorfield(vecnormal, mesh)

       !Plot lat-lon interpolation maps
       filename=trim(simulname)//"_llestim_"//trim(mesh%name)
       call plot_vectorfield_sphgrid(nlon, nlat, llrecon, filename)

       filename=trim(simulname)//"_llerror_"//trim(mesh%name)
       call plot_vectorfield_sphgrid(nlon, nlat, llerror, filename)

       filename=trim(simulname)//"_llexact_"//trim(mesh%name)
       call plot_vectorfield_sphgrid(nlon, nlat, llexact, filename)

       filename=trim(simulname)//"_llnerror_"//trim(mesh%name)
       call plot_scalarfield_sphgrid(nlon, nlat, llnerror, filename)

       print*,"Run plot.sh to see results"

    end if

    print*
    return

  end subroutine vector_reconstruction_test

  subroutine vec_tg_reconstruction_test(mesh)
    !-----------------------------------------------
    !  TEST VECTOR RECONSTRUCTION
    !   Test vector reconstruction routines
    !-----------------------------------------------
    !Mesh structure
    type(grid_structure) :: mesh

    !Indexes
    integer (i4):: i, j

    !Units for files
    logical:: ifile
    integer:: errorsunit
    character (len=256):: filename

    !Errors
    real (r8):: error2
    real (r8):: errormax
    real (r8):: indmax
    real (r8):: ind2
    real (r8):: energy
    real (r8):: tgcomp
    real (r8):: vecmax


    !Interpolable vector variable
    type(vector_field_cart):: vecedfull !Full vector field
    type(vector_field_cart):: vecnormal !Normal comp of vectors
    type(vector_field_cart):: vecrec    !Reconstructed vectors
    type(vector_field_cart):: vecerror  !Error vector
    type(vector_field_cart):: vecexact  !Exact vector field
    type(scalar_field):: vecnerror !Magnitude of error field
    type(scalar_field):: trskind !Consistency index


    !Variable with normal components
    type(scalar_field):: vecncomp
    type(scalar_field):: vectgcomp

    !Auxiliar variables
    real (r8), dimension(1:3):: p !Arbitrary R3 point

    !RBF matrix vector
    type(rbf_matrix_structure), allocatable :: rbf_mat(:)

    integer (i4):: vecreconpos      !Position for vector reconstruction
    integer (i4):: vecpos           !Position for data vectors
    integer (i4):: vecpos_out       !Position for out data vectors
    !logical :: intisrecon           !Flag when interpol mtd is the recon mtd

    !Time counting variables
    real(r8):: time_remap
    !real(r8):: start
    !real(r8):: finish


    !RBF stencil and kind of interpolation
    !character(len=4):: stencil

    print*
    print*,"Test tangent vector reconstruction"
    print*,"-----------------------------"
    print*

    !-------------------------------------------
    !  Read parameters from file "simul.par"
    !------------------------------------------
    !Read kinterp for kind of vector reconstruction to be used
    call getsimulpars(mesh)

    !-----------------------------------------
    ! Define variables on mesh
    !----------------------------------------
    select case (trim(stag))
    case ("HCT", "HWT") !Vectors on triangle edges instersec hx edges
       vecpos=6
    case ("HC", "HW") !Vectors on Hexagon edge midpoints
       vecpos=3
    case default
        print*, "vec_tg_reconstruction_test error: I only know how to "//&
          " reconstruct tg vectors from HC and HCT grids - for now!"
        stop
    end select

    vecpos_out=vecpos !same as data

    !Swap position for reconstructed vector for HW grids
    if(trim(stag)=="HW")then
        vecpos_out=6 !edge intersection
    elseif(trim(stag)=="HWT")then
        vecpos_out=3 !same as data
    end if

    !Full vector
    vecedfull%pos=vecpos
    vecedfull%n=mesh%ne
    allocate(vecedfull%p(1:vecedfull%n))

    !Normal component as vector
    vecnormal%pos=vecpos
    vecnormal%n=mesh%ne
    allocate(vecnormal%p(1:vecnormal%n))

    !Normal component - scalar
    vecncomp%pos=vecpos
    vecncomp%n=mesh%ne
    allocate(vecncomp%f(1:vecncomp%n))

    !Tangent component - scalar
    vectgcomp%pos=vecpos_out
    vectgcomp%n=mesh%ne
    allocate(vectgcomp%f(1:vectgcomp%n))


    select case (trim(stag))
    case ("HCT", "HWT") !Vectors on Hexagon edges
       !Calculate vector field and components
       vecmax=0._r8
       do i=1,mesh%ne
          p=mesh%ed(i)%c%p
          vecedfull%p(i)%v=vecfield(p)
          vecmax=max(vecmax,norm(vecedfull%p(i)%v))
          vecncomp%f(i)=dot_product(vecedfull%p(i)%v, mesh%ed(i)%tg)
          if(trim(stag)=="HWT")then !Reconstruct to intersection point
            vectgcomp%f(i)=dot_product(vecedfull%p(i)%v, mesh%edhx(i)%tg)
          else
            vectgcomp%f(i)=dot_product(vecedfull%p(i)%v, mesh%ed(i)%nr)
          end if
          vecnormal%p(i)%v=vecncomp%f(i)*mesh%ed(i)%tg
       end do
    case ("HC", "HW") !Vectors on Hexagon edges
       !Calculate vector field and components
       vecmax=0._r8
       do i=1,mesh%ne
          p=mesh%edhx(i)%c%p
          vecedfull%p(i)%v=vecfield(p)
          vecmax=max(vecmax,norm(vecedfull%p(i)%v))
          vecncomp%f(i)=dot_product(vecedfull%p(i)%v, mesh%edhx(i)%nr)
          !Tangent component - scalar
          if(trim(stag)=="HW")then !Reconstruct to intersection point
            p=mesh%ed(i)%c%p
            vectgcomp%f(i)=dot_product(vecfield(p), mesh%ed(i)%nr)
          else  !Reconstruct to midpoint
            vectgcomp%f(i)=dot_product(vecedfull%p(i)%v, mesh%edhx(i)%tg)
          end if
          vecnormal%p(i)%v=vecncomp%f(i)*mesh%edhx(i)%nr
       end do
    end select

    !----------------------------------------
    !Alocate space for reconstructed vector field
    !----------------------------------------

    vecreconpos=vecpos_out
    vecrec%n=mesh%ne

    !Reconstructed vectors
    vecrec%pos=vecreconpos
    allocate(vecrec%p(1:vecrec%n))

    !Consistency index for reconstruction
    trskind%pos=vecreconpos
    trskind%n=vecrec%n
    allocate(trskind%f(1:trskind%n))

    !Reconstruction error
    vecerror%pos=vecreconpos
    vecerror%n=vecrec%n
    allocate(vecerror%p(1:vecerror%n))

    !Exact vector field
    vecexact%pos=vecreconpos
    vecexact%n=vecrec%n
    allocate(vecexact%p(1:vecexact%n))

    !Norm of the reconstruction error
    vecnerror%pos=vecreconpos
    vecnerror%n=vecrec%n
    allocate(vecnerror%f(1:vecnerror%n))

    !-----------------------------------------
    ! Test reconstruction
    !----------------------------------------
    print*
    print*,"Test reconstruction re-map on geodesic grid"
    print*

    !Reconstruct vector field
    !call cpu_time(start)
    !do i=1,5
    !call vrec_remap(vecncomp, vecrec, mesh, recon_mtd)
    !end do
    !call cpu_time(finish)
    time_remap=0.
    !time_remap = (finish-start)/5.
    !time_remap = (finish-start)
    !time_remap=time_remap-time_rbf
    call calc_trisk_weights(mesh)

    !Calculate errors
    error2=0._r8
    ind2=0._r8
    errormax=0._r8
    indmax=0._r8
    energy=0._r8
    do i=1,  vecrec%n !For each edge
       if(vecpos_out==6)then !reconstruct to edge intersection
          p=mesh%ed(i)%c%p
       elseif(vecpos_out==3)then !Do reconstruction to hx edges midpoint
          p=mesh%edhx(i)%c%p
       endif
       vecexact%p(i)%v=vecfield(p)

       !Consistency index (does not depend on data, only on grid)
       trskind%f(i)=tgrecon_index(i, recon_mtd%recon, vecpos, mesh, vecreconpos)

       !Calculate reconstruction and errors
       !Calculate error only of the tg component
       if(trim(recon_mtd%recon)=="trsk")then !.and. trim(stag)=="HW")then
        vecrec%p(i)%v=vecrecon_trsk_tg (i, vecncomp, mesh, p) !Returns just tg vector
        if(vecpos_out==3)then
            vecerror%p(i)%v=vectgcomp%f(i)*mesh%edhx(i)%tg-vecrec%p(i)%v !
        elseif(vecpos_out==6)then
            vecerror%p(i)%v=vectgcomp%f(i)*mesh%ed(i)%nr-vecrec%p(i)%v !
        endif
        vecnerror%f(i)=norm(vecerror%p(i)%v)
       else !Calculate error based on full vector - should give same results if vecpos_out=vecpos
        vecrec%p(i)%v=vector_reconstruct (p, vecncomp, mesh, recon_mtd%recon, rbf_mat, i)
        vecerror%p(i)%v=vecedfull%p(i)%v-vecrec%p(i)%v
        vecnerror%f(i)=norm(vecerror%p(i)%v)
       endif
       !print*,vecnerror%f(i),   trskind%f(i)
       !print*, i, vecnerror%f(i)
       errormax=max(errormax, vecnerror%f(i))
       indmax=max(indmax, trskind%f(i))
       error2=error2+vecnerror%f(i)**2
       ind2=ind2+trskind%f(i)**2
       !dot_product(vecexact%p(i)%v, mesh%edhx(i)%tg) !
       if(vecpos_out==3)then
          tgcomp=dot_product(vecrec%p(i)%v, mesh%edhx(i)%tg)
       elseif(vecpos_out==6)then
          tgcomp=dot_product(vecrec%p(i)%v, mesh%ed(i)%nr)
       endif
       energy=energy+mesh%ed(i)%leng*mesh%edhx(i)%leng*vecncomp%f(i)*tgcomp/2
       !energy=energy+mesh%ed(i)%lenp*mesh%edhx(i)%lenp*vecncomp%f(i)*tgcomp/2
    end do
    !Normalize energy with respect to vector field
    energy=energy/(vecmax*vecmax)
    !print*, vecmax
    !Error 2
    error2=dsqrt(error2/real(vecrec%n, r8))
    ind2=dsqrt(ind2/real(vecrec%n, r8))
    print *, "Max and L2 errors for remap: ", errormax, error2
    print *, "Max and L2 norms of index : ", indmax, ind2
    print*, "Energy: ", energy
    !print*, "Elapsed time for remap: ", time_remap
    !print*, "Elapsed time for rbf matrix build:", time_rbf

    !File for errors
    filename=trim(datadir)//"vrecontgremap_errors.txt"
    call getunit(errorsunit)

    inquire(file=filename, exist=ifile)
    if(ifile)then
       open(errorsunit,file=filename, status='old', position='append')
    else
       open(errorsunit,file=filename, status='replace')
       write(errorsunit, '(a)') "        n              meanvdist          "//&
            "          kinterp   stag    func          error(inf)          error(L2)"//&
            "         indmax    ind2     energy    time                             mesh "
    end if
    ! Write errors
    write(errorsunit, '(i12, f22.8, a24, a8, i8, 6e22.8, a32)') mesh%nv,  &
         mesh%meanvdist*rad2deg, trim(recon_mtd%name), trim(stag), testfunc, errormax, error2, &
         indmax, ind2, energy, time_remap, trim(mesh%name)


    if(plots)then

       !Plot vector remaps error
       vecnerror%name=trim(simulname)//"_vecnerror"
       call plot_scalarfield(vecnerror, mesh)

       trskind%name=trim(simulname)//"_tgind"
       call plot_scalarfield(trskind, mesh)

       !Plot vector field maps
       vecrec%name=trim(simulname)//"_vecrec"
       call plot_cart_vectorfield(vecrec, mesh)

       vecexact%name=trim(simulname)//"_vecexact"
       call plot_cart_vectorfield(vecexact, mesh)

       vecerror%name=trim(simulname)//"_vecerror"
       call plot_cart_vectorfield(vecerror, mesh)

       !vecedfull%name=trim(simulname)//"_vecfull"
       !call plot_cart_vectorfield(vecedfull, mesh)

       !vecnormal%name=trim(simulname)//"_vecnormal"
       !call plot_cart_vectorfield(vecnormal, mesh)


       print*,"Run plot.sh to see results"

    end if

    print*
    return

  end subroutine vec_tg_reconstruction_test


  !-----------------------------------------------------------------------------------------
  !    FUNCTIONS FOR TESTS
  !-----------------------------------------------------------------------------------------

  function f(p)
    !-----------------------------------
    !  F - initial conditions for scalar fields
    !
    !   p is a point in cartesian coords
    !
    !  P. Peixoto - Feb2012
    !---------------------------------------------
    real (r8), intent(in) :: p(1:3)
    real (r8):: lon
    real (r8):: lat

    real (r8):: f

    !General parameters
    real (r8):: h1
    real (r8):: h2

    !Cosbell parameters
    real (r8):: b
    real (r8):: c
    real (r8):: r
    real (r8):: r1
    real (r8):: r2

    !Gaussian parameters
    real (r8), dimension(1:3) :: p1
    real (r8), dimension(1:3) :: p2
    real (r8):: b0

    !Center of bell, slot, hill
    real (r8):: lat1
    real (r8):: lon1
    real (r8):: lat2
    real (r8):: lon2

    !Wavelength of test function
    real (r8):: m
    real (r8):: n

    m=1
    n=1

    !convert to spherical coords
    call cart2sph ( p(1), p(2), p(3), lon, lat )

    if(testfunc>0 .and. testfunc <4)then
       !Deformational test functions
       lon1=-pi/(4._r8)
       lat1=0._r8
       lon2=pi/(4._r8)
       lat2=0._r8
    end if

    !scalar fields
    select case(testfunc)
    case(1) !Cosine bell
       b=0.1_r8
       c=1._r8 !0.9
       r=1._r8/2._r8
       r1= arcdistll(lon, lat, lon1, lat1)
       if(r1<r)then
          h1=(1._r8/2._r8)*(1._r8+dcos(pi*r1/r))
          f=b+c*h1
          return
       end if
       r2= arcdistll(lon, lat, lon2, lat2)
       if(r2<r)then
          h2=(1._r8/2._r8)*(1._r8+dcos(pi*r2/r))
          f=b+c*h2
          return
       end if
       f=b

    case(2) !Gaussian
       b0=5
       call sph2cart(lon1, lat1, p1(1), p1(2), p1(3))
       call sph2cart(lon2, lat2, p2(1), p2(2), p2(3))
       h1=dexp(-b0*norm(p-p1)**2)
       h2=dexp(-b0*norm(p-p2)**2)
       f=h1+h2

    case(3) !Non-smooth- Slotted cylinder
       b=0.1_r8
       c=1._r8 !0.9
       r=1._r8/2._r8
       f=b

       r1= arcdistll(lon, lat, lon1, lat1)
       if(r1<=r.and.abs(lon-lon1)>=r/(6._r8))then
          f=c
          return
       end if
       r2= arcdistll(lon, lat, lon2, lat2)
       if(r2<=r.and.abs(lon-lon2)>=r/(6._r8))then
          f=c
          return
       end if
       if(r1<=r.and.abs(lon-lon1) < r/(6._r8) .and. &
            (lat-lat1)<-(5._r8)*r/(12._r8))then
          f=c
          return
       end if
       if(r2<=r.and.abs(lon-lon2)<r/(6._r8).and. &
            (lat-lat2)>(5._r8)*r/(12._r8))then
          f=c
          return
       end if

    case(4) !Smooth Trigonometric
       f=((dcos(lat))**3)*((dsin(lon))**2)

    case(5) !Exponential
       f= (exp(p(1))+2*exp(p(2)+p(3)))/10

    case(6) !Linear in R3
       f=(1.0+2*p(1)+3*p(2)+4*p(3))/6.0

    case(7) !Constant 1
       f=1

    case(8) !Trigonom oscilation (beta of HR95)
       f=dcos(m * lon) * dcos(n * lat) ** 4
       !((dcos(lat))**4)*((dsin(lon))**7)

    case default
       print*, "F error: unknown scalar field (testfunc)", testfunc
       stop
    end select

    return
  end function f


  function df(p)
    !-------------------------------------------
    !  Grad of function "f"
    ! used to test interpolations
    ! in cartesian coords
    ! Returns the gradient in cartesian coords
    !------------------------------------------
    real (r8):: df(1:3)
    real (r8):: p(1:3)
    real (r8):: vlat
    real (r8):: vlon
    real (r8):: lat
    real (r8):: lon


    !Wavelength of test function
    real (r8):: m
    real (r8):: n

    m=1
    n=1

    call cart2sph (p(1), p(2), p(3), lon, lat)
    df=0.

    select case(testfunc)
    case(1) !Cosine bell

    case(2) !Gaussian

    case(3) !Non-smooth- Slotted cylinder
       df=0
    case(4) !Smooth Trigonometric
       !f=((dcos(lat))**3)*((dsin(lon))**2)
       vlon=2*dsin(lon)*dcos(lon)*((dcos(lat))**2)
       vlat=-3*((dcos(lat))**2)*((dsin(lon))**2)*dsin(lat)
       call convert_vec_sph2cart(vlon, vlat, p, df)
    case(5) !Exponential
       !f= (exp(p(1))+2*exp(p(2)+p(3)))/10
       df(1)=exp(p(1))/10
       df(2)=exp(p(2)+p(3))/5
       df(3)=exp(p(2)+p(3))/5
    case(6) !Linear in R3
       !f=(1.0+2*p(1)+3*p(2)+4*p(3))/6.0
       df(1)=2./6.
       df(2)=3./6.
       df(3)=4./6.
    case(7) !Constant
       df=0
    case(8) !Oscilating trig (beta of HR95)
       vlon = -dsin(m*lon)*m*dcos(n*lat)**4 / dcos(lat)
       vlat = -4._r8*dcos(m*lon)*dcos(n*lat)**3 * dsin(n*lat)* n
       !vlon=7*(dsin(lon)**6)*dcos(lon)*((dcos(lat))**3)
       !vlat=-4*((dcos(lat))**3)*((dsin(lon))**7)*dsin(lat)
       call convert_vec_sph2cart(vlon, vlat, p, df)
    case default
       print*, "DF error: unknown gradient for field (testfunc)", testfunc
       stop
    end select

    !Project vector to the sphere (if not already)
    df=proj_vec_sphere(df, p)

    return
  end function df

  function lap_exact(p)
    !-----------------------------------
    !  Laplacian of scalar field given in f()
    !
    !   p is a point in cartesian coords
    !---------------------------------------------
    real (r8), intent(in) :: p(1:3)
    real (r8):: lon
    real (r8):: lat

    real (r8):: lap_exact

    !Wavelength of test function
    real (r8):: m
    real (r8):: n

    !convert to spherical coords
    call cart2sph ( p(1), p(2), p(3), lon, lat )

    !Laplacian
    select case(testfunc)
    case(4) !Smooth Trigonometric
       lap_exact=dcos(lat) * (-5.0_r8 * dcos(lon) ** 2 + 7.0_r8 &
            - 12.0_r8 * dcos(lat)** 2 * dsin(lon) ** 2)
    case(7) !Constant 1
       lap_exact=0
    case (8) !Oscilat trig (beta of HR95)
       m=1
       n=1
       lap_exact = -dcos(m*lon)*dcos(n*lat)**2 * (m**2*dcos(n * lat) **2 &
            - 4._r8 * dsin(lat) * dcos(n*lat)*dsin(n*lat)*n*dcos(lat) &
            - 12._r8 * dcos(lat)** 2 * n** 2 + 16._r8*dcos(lat)**2 *dcos(n* lat)**2* n** 2) &
            / dcos(lat) ** 2
       !lap_exact= (dcos(lat)**2) * (dsin(lon) **5 )* &
       !    (33._r8 * dcos(lon) ** 2 + 9.0_r8 &
       !    - 20.0_r8 * dcos(lat)** 2 * dsin(lon) ** 2)
    case default
       print*, "Laplacian (lap_exact) error: unknown scalar field (testfunc) :", testfunc
       stop
    end select

    return
  end function lap_exact


  function vecfield(p)
    !-------------------------------------------
    !  Vector fields
    ! used to test vector interpolations
    ! in cartesian coords
    !------------------------------------------
    real (r8):: vecfield(1:3)
    real (r8):: vec(1:3)
    real (r8):: p(1:3)
    real (r8):: lat
    real (r8):: lon
    real (r8):: m
    real (r8):: n
    real (r8):: rho
    real (r8):: lonref
    real (r8):: latref

    real (r8):: omega_m
    real (r8):: omega_0
    real (r8):: a
    integer (i4):: wave_m

    !Convert to spherical coords
    call cart2sph (p(1), p(2), p(3), lon, lat)
    vec=0.

    select case(testfunc)
    case(1:3, 6) !Deformal test cases and trignometric case
       m=1
       n=1
       !m=3
       !n=3
       call convert_vec_sph2cart(u(lon, lat), v(lon, lat), p, vec)
    case(4) !Constant Parallel Field
       rho=dsqrt(p(1)**2+p(2)**2)
       if(rho<=eps)then
          vec=0.
       else
          vec(1)=-p(2)/rho
          vec(2)=p(1)/rho
          vec(3)=0._r8
       end if
    case(5) !Rotation around z-axis
       !vec(1)=-p(2)
       !vec(2)=p(1)
       !vec(3)=0._r8
       call convert_vec_sph2cart(u(lon, lat), v(lon, lat), p, vec)
    case(7) !Rotation with angle
       !New north pole
       latref=45._r8*deg2rad
       lonref=45._r8*deg2rad
       call convert_vec_sph2cart(u(lon, lat), v(lon, lat), p, vec)
       !print*, lon, lat, u(lon, lat), v(lon, lat)
    case(8) !Rossby-Haurwitz wave (Wan Thesis 2009)
       a=erad
       omega_m=7.848e-6_r8
       omega_0=7.848e-6_r8
       !wave_m=4_i4
       wave_m=8_i4
       !print*, a, omega_m, omega_0, m
       !print "(4f32.8)", lon, lat, u(lon, lat), v(lon, lat)
       call convert_vec_sph2cart(u(lon, lat), v(lon, lat), p, vec)
    case default
       print*, "Vecfield error: unknown vector field (testfunc)", testfunc
       stop
    end select

    !To ensure it is tangent to the sphere, project it
    vecfield=proj_vec_sphere(vec, p)

    return

  contains

    function u(lon, lat)
      !-----------------------------------
      !  U - velocity in West-East direction
      !   for a given testcase
      !
      !   Lat in [-pi/2,pi/2], Lon in [-pi,pi]
      !
      !  P. Peixoto - Feb2011
      !---------------------------------------------
      real (r8), intent(in) :: lon
      real (r8), intent(in) :: lat
      real (r8):: u

      !Auxiliar variables
      real (r8):: k

      u=0.
      ! Velocity for each testcase
      select case(testfunc)
      case(1) !Nondivergent case-1
         k=2.4
         u=k*(dsin((lon+pi)/2.)**2)*(dsin(2.*lat))

      case(2) !Nondivergent case-2
         k=2.
         u=k*(dsin((lon+pi))**2)*(dsin(2.*lat))

      case(3) !Divergent case-3
         k=1.
         u=-k*(dsin((lon+pi)/2._r8)**2)*(dsin(2.*lat))*(dcos(lat)**2)
      case(5)
         u=dcos(lat) !*dsin(lat)
      case(6) !Trigonometric (Tomita)
         if(abs(lat*rad2deg)>=90 - eps)then
            u=0
         else
            !u=-(dsin(lon)**2)*(dcos(lat)**3)
            u=-m*(dsin(lon)*dsin(m*lon)*dcos(n*lat)**4)/dcos(lat)
         end if

      case(7) !Rotated rotation (Ritchie 1987 eq 11)
         u=dcos(lat)*dsin(latref)-dcos(lon-lonref)*dsin(lat)*dcos(latref)

      case(8) !Rossby  Haurwitz wave
         u=a*omega_0*dcos(lat)+ &
              a*omega_m*(dcos(lat))**(wave_m-1)* &
              (wave_m*dsin(lat)**2 - dcos(lat)**2) * dcos(wave_m*lon)
      end select

      return
    end function u

    function v(lon, lat)
      !-----------------------------------
      !  V - velocity in South-North direction
      !   for a given testcase
      !
      !   Lat in [-pi/2,pi/2], Lon in [-pi,pi]
      !
      !  P. Peixoto - Feb2011
      !---------------------------------------------
      real (r8), intent(in) :: lon
      real (r8), intent(in) :: lat
      real (r8):: v

      !Auxiliar variables
      real (r8):: k

      ! Velocity for each testcase

      v=0.0
      select case(testfunc)
      case(1) !Nondivergent case-1
         k=2.4
         v=k*(dsin(lon+pi))*(dcos(lat))/2._r8
      case(2) !Nondivergent case-2
         k=2
         v=k*(dsin(2*(lon+pi)))*(dcos(lat))
      case(3) !Divergent case-3
         k=1
         v=(k/2._r8)*(dsin((lon+pi)))*(dcos(lat)**3)
      case(5)
        v=0._r8
      case(6) !Trigonometric (Tomita)
         !v=-4*(dsin(lon)**2)*dcos(lon)*(dcos(lat)**3)
         v=-4*n*(dcos(n*lat)**3)*dsin(n*lat)*dcos(m*lon)*dsin(lon)
      case(7) !Rotated rotation (Ritchie 1987 eq 12)
         v=dsin(lon-lonref)*dcos(latref)
      case(8) !Rossby  Haurwitz wave
         v = - a*omega_m*wave_m*((dcos(lat))**(wave_m-1))* &
              dsin(lat)*dsin(wave_m*lon)
      end select
      return
    end function v

  end function vecfield

  function div_exact(p)
    ! Exact divergence of velocity field

    real (r8):: p(1:3)
    real (r8):: lat
    real (r8):: lon

    real (r8):: div_exact
    real (r8):: m
    real (r8):: n
    real (r8):: k

    call cart2sph(p(1), p(2), p(3), lon, lat)

    !Calculate divergence for field with non null divergence
    select case(testfunc)
    case(3) !Deformal test cases divergence
       k=1.
       div_exact=k*dcos(lat)*dsin(2.*lat)*dsin((1./2.)*lon)*dcos((1./2.)*lon) &
            +2*k*dcos(lat)**2*dsin(lon)*dsin(lat)
    case(6) !Trigonometric (Tomita)
       m=1
       n=1
       !div_exact=(-2*(dcos(lat)**3)*dsin(lon)*dcos(lon)+ &
       !   16*(dsin(lon)**2)*dcos(lon)*(dcos(lat)**3)*dsin(lat))/dcos(lat)
       !div_exact=1._r8/dcos(lat)*(-dcos(lon)/dcos(lat)*dsin(m*lon)*m*dcos(n*lat)**4 &
       !     - dsin(lon)/dcos(lat)*dcos(m*lon)*m**2*dcos(n*lat)**4 &
       !     + 12._r8*dsin(lon)*dcos(m*lon)*dcos(n*lat)**2*dsin(n*lat)**2*n**2*dcos(lat) &
       !     - 4._r8*dsin(lon)*dcos(m*lon)*dcos(n*lat)**4*n**2*dcos(lat) &
       !     + 4._r8*dsin(lon)*dcos(m*lon)*dcos(n*lat)**3*dsin(n*lat)*n*dsin(lat))
       div_exact=(-dcos(lon) * dsin(m * lon) * m * dcos(n * lat) ** 4 / dcos(lat) - &
            dsin(lon) * dcos(m * lon) * m ** 2 * dcos(n * lat) ** 4 / dcos(lat) + &
            12.0 * dsin(lon) * dcos(m * lon) * dcos(n * lat) ** 2 * dsin(n *lat) ** 2 * n ** 2 * dcos(lat) - &
            4.0 * dsin(lon) * dcos(m * lon) * dcos(n * lat) ** 4 * n ** 2 * cos(lat) + &
            4.0 * dsin(lon) * dcos(m * lon) * dcos(n * lat) ** 3 * dsin(n * lat) * n * dsin(lat)) / dcos(lat)
    case(8)
      print*, "Exact RW wave divergence not implemented"
      stop
    case default
       div_exact=0.
    end select

    return
  end function div_exact

   function rot_exact(p)
    ! Exact rotational of velocity field

    real (r8):: p(1:3)
    real (r8):: lat
    real (r8):: lon

    real (r8):: rot_exact
    !real (r8):: m
    !real (r8):: n
    real (r8):: wave_m
    real (r8):: omega_m
    real (r8):: omega_0
    real (r8):: a

    call cart2sph(p(1), p(2), p(3), lon, lat)

    !Calculate divergence for field with non null divergence
    select case(testfunc)
    case(5) !Rotation about z azis
      rot_exact=2.*dsin(lat)
    case(8) !
        a=erad
       omega_m=7.848e-6_r8
       omega_0=omega_m
       !Assuming omega_0=omega_m
       !wave_m=4_i4
       wave_m=4_i4
        rot_exact=-a *omega_0* (-2 + (2 + 3* wave_m +wave_m**2) *dcos(lat)**wave_m &
          *dcos(wave_m*lon)) *dsin(lat)
    case default
       print*, "Vector field with curl/rotational not implemented"
       stop
    end select

    return
  end function rot_exact

  !----------------------------------------------------------------------------------------
  !   Tests for mesh structure and others
  !---------------------------------------------------------------------------------------

  subroutine test_geo2reg(mesh)
    !-----------------------------------------------
    !  TEST_GEO2REG
    !   Test geodesic to regular grid conversion tool 
    !-----------------------------------------------
    type(grid_structure), intent(in):: mesh
    integer (i4):: i
    integer (i4):: j
    integer (i4):: k
    integer (i4):: kt
    integer (i4):: itmp
    real (r8):: tlat
    real (r8):: tlon
    real (r8):: latij
    real (r8):: lonij
    real (r8), dimension(1:3) :: p


    print*
    print*, "Test geodesic to regular grid conversion tool ..."
    print*
    print*, "Set lon (degrees):"
    !read (*, *) tlon
    print*, "Set lat (degrees):"
    !read (*, *) tlat
    !Point
    !tlon=-64.5
    !tlat=-13.5

    !tlon=35.9121339133462_r8
    !tlat=-47.1678780596014_r8
    !tlon=-168.5
    !tlat=-86.5
    !tlon=-160.5
    !tlat=-89.5
    !tlon=-36.0878869808453
    !tlat=58.2825557357298
    !tlon=-23.8744234108428
    !tlat=34.7927455551923
    !tlon=-111.857582256335
    !tlat=28.0029729308683
    tlon=-161.875
    tlat=-53.875

    call sph2cart(tlon*deg2rad,tlat*deg2rad, p(1), p(2),p(3))
    call convllij(tlon*deg2rad, tlat*deg2rad, i, j, mesh)
    call convijll(i, j, lonij, latij, mesh)
    print*, "Point:"
    print*, "Lon:", tlon, " in ", i, lonij*rad2deg
    print*, "Lat:", tlat, " in ", j, latij*rad2deg
    !print*, " (x,y,z):", p
    print*
    print*, "Sq higher corner:", (lonij+mesh%qd_lat(j)%dlon)*rad2deg, &
         (latij+mesh%dlat)*rad2deg


    print*, "Intersecting triangles: "
    do k=1, mesh%qd_lat(j)%qd_lon(i)%ntr
       kt=mesh%qd_lat(j)%qd_lon(i)%tr(k)
       itmp=sqtriintersec(i, j, kt, mesh)
       print '(i8,a,2f8.3,a,i4)', kt, "  CC:", mesh%tr(kt)%c%lon*rad2deg, &
            mesh%tr(kt)%c%lat*rad2deg, " Really intersects? :", itmp
    end do

    print*
    print*,"Point in triangle:"

    kt=gettr(p, mesh)
    print*,kt
    print*

    print*, "Nearest node to point"
    k=getnearnode(p, mesh)
    print*,"node:",k," lon", mesh%v(k)%lon*rad2deg," lat",mesh%v(k)%lat*rad2deg


    !p=mesh%v(mesh%tr(kt)%v(1))%p+mesh%v(mesh%tr(kt)%v(2))%p
    !p=p/norm(p)
    !kt=gettr(p, mesh)
    !print*,kt
    !p(1)=p(1)+0.1_r8*eps*p(1)
    !p=p/norm(p)
    !kt=gettr(p, mesh)
    !print*,kt

    print*



    !print*, "Projection of point on the planar triangle"

    !kt=549
    !call sph2cart(tlon*deg2rad,tlat*deg2rad, p(1), p(2), p(3))

    !p1=mesh%v(mesh%tr(kt)%v(1))%p
    !p2=mesh%v(mesh%tr(kt)%v(2))%p
    !p3=mesh%v(mesh%tr(kt)%v(3))%p
    !test=insidetr(p, kt, mesh)

    !print*,p
    !print*,p1
    !print*,p2
    !print*,p3
    !print*
    !print*, test


    print*
    print*,"End of geodesic to regular grid conversion tool test."

    return
  end subroutine test_geo2reg

  subroutine test_trsearch(mesh)
    !-----------------------------------------------
    !  TEST_TRSEARCH
    !   Test triangle search routines
    !-----------------------------------------------
    type(grid_structure) :: mesh

    real (r8):: tlon
    real (r8):: tlat
    real (r8):: p(1:3)
    integer:: i
    integer:: j
    integer:: nst
    integer:: unit
    integer:: nmax
    character (len=60):: filename
    logical:: ifile

    !Time counting variables
    real(r8):: elapsed_time
    real(r8):: start
    real(r8):: finish
    nst=1
    nmax=1000000
    call cpu_time(start)
    do i=1,1000000
       tlon=real(i*j, r8)*360/nmax-180
       tlat=real(j, r8)*180/nmax-90
       !nst=1
       call sph2cart(tlon*deg2rad,tlat*deg2rad, p(1), p(2),p(3))
       !call trfind ( nst, p, mesh%n, mesh%x, mesh%y, mesh%z, mesh%list, &
       !     mesh%lptr, mesh%lend, b1, b2, b3, i1, i2, i3 )         
       j=gettr(p, mesh)
    end do
    call cpu_time(finish)
    elapsed_time= finish-start !real(clock_end-clock_start, r8)/real(clock_rate,r8)
    print*, "Search time: ", elapsed_time

    filename=trim(datadir)//"trsearchtimes.txt"
    call getunit(unit)   
    inquire(file=filename, exist=ifile)
    if(.not.ifile)then
       open(unit,file=filename, status='replace')
       write(unit, '(a)') "       n    time"
    else
       open(unit,file=filename, status='old', position='append')
    end if
    write(unit, '(i8, 5e16.8)') mesh%nv,  elapsed_time
    close(unit)
    return
  end subroutine test_trsearch

  subroutine choleskytest()
    !-----------------------------------------------
    !  CHOLESKY TEST
    !   Test cholesky decomposition
    !-----------------------------------------------
    integer, parameter :: n=8
    real (r8):: A(1:n,1:n)
    real (r8):: L(1:n,1:n)
    real (r8):: x(1:n)
    real (r8):: b(1:n)
    integer:: i
    integer:: j
    print*
    print*, "Cholesky Solver"
    print*


    do i=1,n
       do j=1,n
          A(i,j)= 1._r8/(real(i+j-1, r8))
       end do
       b(i)=1
    end do

    print*, "A"
    print '(8f8.4)', transpose(A)
    print*
    call choleskydecomp(A, n, L)
    call choleskysolve(L, x, b, n)
    print*
    print*, x

  end subroutine choleskytest

  subroutine rbftest(mesh)
    !-----------------------------------------------
    !  RBF TEST
    !   
    !-----------------------------------------------

    !RBF stencil type
    character (len=4):: stencil

    !RBF matrix vector
    type(rbf_matrix_structure), allocatable :: rbf_mat(:)

    !Mesh structure
    type(grid_structure) :: mesh


    stencil="HX"

    print*
    print*, "RBF tester"
    print*
    rbf_par=1.

    call rbf_matrix_build( stencil, rbf_par, rbf_mat, mesh )

  end subroutine rbftest


  subroutine natneibinterptest(mesh)
    !-----------------------------------------------
    !  NATURAL NEIGHBOUR INTERPOLATION TEST
    !   
    !-----------------------------------------------
    !Mesh structure
    type(grid_structure) :: mesh


    !Natural coords
    type (general_coords):: lap
    type (general_coords):: sib

    real (r8), dimension(1:3) :: p
    real (r8), dimension(1:3) :: p1
    real (r8), dimension(1:3) :: p2
    real (r8), dimension(1:3) :: q1
    real (r8), dimension(1:3) :: q2
    real (r8), dimension(1:3) :: r
    real (r8), dimension(1:3) :: rt
    real (r8), dimension(1:3) :: n
    type(vector), dimension(1:8) :: pts
    integer (i4), allocatable :: listv(:)
    real (r8), allocatable :: listd(:)
    real (r8):: area
    real (r8):: tlon
    real (r8):: tlat
    real (r8):: rad

    logical:: intersec

    integer:: i
    integer:: j
    integer:: k

    print*
    print*, "NAT NEIB INTERP TEST"
    print*

    !Point assumed to be in the 42 icosahedral mesh
    p1=mesh%v(4)%p
    p2=mesh%v(5)%p
    q1=mesh%v(25)%p
    q2=mesh%v(31)%p
    rt=(q1+q2)/2
    rt=rt/norm(rt)

    intersec=arcintersec(p1, p2, q1, q2, r)

    print*,intersec
    print*,r

    print*
    call ortogonalarc(p1, p2, p, n)
    print*, p
    print*,n
    r=gcircarcintersec(n, q1, q2)
    print*, r
    print*, rt
    print*
    print*, "-------------------"
    print*
    i=40
    print*, "Cell:", i
    print*, "Cell Area:", mesh%hx(i)%areag

    do j=1,mesh%v(i)%nnb 
       k=mesh%v(i)%tr(j)
       pts(j)%v=mesh%tr(k)%c%p
    end do

    area=sphpolarea(pts(1:mesh%v(i)%nnb), mesh%v(i)%nnb)
    print*, area
    print*
    print*

    print*,"Test point coordinate"

    !Point
    tlat=0
    tlon=1
    print*, "Lon    Lat"
    print*, tlon, tlat
    call sph2cart(tlon*deg2rad, tlat*deg2rad, p(1), p(2), p(3))
    call natural_coords(p, lap, sib, mesh)

    rad=10.5_r8*mesh%maxvdist
    k=550
    call getnearnodes(p, mesh, rad, listv, listd, k)
    print*, "NEARNODES:"
    print*, rad
    print*, k
    do i=1,k
       print*, i, listv(i), listd(i)
    end do
    return
  end subroutine natneibinterptest

  subroutine test_edgeconnections(mesh)
    !-----------------------------------------------
    !  TEST_EDGECONNECTIONS
    !   Test triangle search routines
    !-----------------------------------------------
    type(grid_structure) :: mesh

    real (r8):: tlon
    real (r8):: tlat
    real (r8):: p(1:3)
    integer:: i, j, connect
    character (len=60):: filename
    logical:: ifile
    integer:: unit

    print*, "Testing edge connection..."

    filename=trim(datadir)//"edgeconnections.txt"
    call getunit(unit)
    open(unit,file=filename, status='replace')
    write(unit, '(a)') "       edge1      edge2     connecting_edge"

    do i=1, mesh%ne
      do j=1,mesh%ne
        connect=gethxedgeconnection(i,j,mesh)
        if(connect>0)then
          write(unit, *) i, j, connect
        endif
      end do
    end do

    close(unit)
    print*, "See file: ", filename
    return
  end subroutine test_edgeconnections


end module simulpack
