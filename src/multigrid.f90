module multigrid
  !====================================================================
  ! Main module for Multigrid method
  ! Marline I. Silva - Sept 2014
  !=====================================================================
  !Global constants 
  use constants, only: &
    datadir, &
    i4, &
    pardir, &
    r4, &
    r8

  !Data structures
  use datastruct, only:  &
    grid_structure, &
    scalar_field

  use datastructmult, only: &
    interpolable_meshes, &
    multimesh, &
    multisol

  !Spherical mesh routines
  use smeshpack, only: &
    cart2sph, &
    error_norm_2, &
    error_norm_max, &
    getparameters, &
    getunit, &
    meshbuild

  !Interpolation pack
  use interpack, only: &
    plot_scalarfield

  implicit none

contains

  !=====================================================================
  !                     functions for tests
  !===================================================================== 

  function func(p)
    real (r8)::p(1:3)
    real(r8)::func
    real(r8)::lat
    real(r8)::lon
    real(r8)::m
    real(r8)::n

    m=1
    n=1

    !TESTE 1
    !func=(1+2*p(1)+3*p(2)+4*p(3))/6._r8

    !TESTE 2
    !func = (exp(p(1))+ 2*(exp(p(2)+p(3))))/10

    call cart2sph (p(1), p(2), p(3), lon, lat)
    !TESTE 3
    !func=((dcos(lat))**3)*((dsin(lon))**2)

    !TESTE 4
    func=dcos(m*lon)*dcos(n*lat)**4

    !TESTE 5
    !func=dsin(lat)

    !TESTE 6
    !func=dsin(3*lon)*((dcos(3*lat)))**4

    !TESTE 7
    !func=-(1/2)*((dcos(lat))**3)*dcos(3*lon)

    !TESTE 8
    !func=dcos(lat)*dcos(lon)

    return
  end function func

  !----------------------------------------------------------------
  !     Exact Laplacian of the preceding functions
  !----------------------------------------------------------------

  function lap_ext(p)
    real (r8)::p(1:3)
    real (r8):: lat
    real (r8)::lon
    real (r8)::m
    real (r8)::n
    real(r8)::lap_ext

    m=1
    n=1

    !Converte para coordenadas esféricas
    call cart2sph (p(1), p(2), p(3), lon, lat)

    !TESTE 1 
    !lap_ext= (((3*dsin(lon)+2*dcos(lon))*(1-dcos(lat)**2+dsin(lat)**2))/6._r8*dcos(lat)) &
    !- 4._r8*dsin(lat)/3._r8

    !TESTE 2
    !lap_ext=(1.0_r8/10.0_r8*dcos(lat))*(2.0_r8*exp(dcos(lat)*dsin(lon)+dsin(lat))*(dcos(lat) &
    !*((dcos(lon))**2)-dsin(lon))+exp(dcos(lat)*dcos(lon))*(dcos(lat)*((dsin(lon))**2) &
    !-dcos(lon))-2.0_r8*dsin(lat)*exp(dcos(lat)*dsin(lon)+dsin(lat))*(dcos(lat)-dsin(lon) &
    !*dsin(lat))+dcos(lon)*((dsin(lat))**2)*exp(dcos(lat)*dcos(lon))) +0.1*(exp(dcos(lat) &
    !*dcos(lon))*(((dcos(lon))**2)*((dsin(lat))**2)-dcos(lon)*dcos(lat))+2.0_r8*exp(dcos(lat) &
    !*dsin(lon)+dsin(lat))*((dcos(lat)-dsin(lon)*(dsin(lat)))**2)-(dsin(lon)*dcos(lat) &
    !+dsin(lat)))

    !TESTE 3 (-lap u + u=f)
    !lap_ext=-dcos(lat)*(-5.0_r8*dcos(lon)**2+7.0_r8-12.0_r8*dcos(lat)**2*dsin(lon)**2)  &
    !    +((dcos(lat))**3)*((dsin(lon))**2)

    !-------------------------------------------------------------
    !TESTE 4
    ! lap_ext = -dcos(m*lon)*dcos(n*lat)**2 * (m**2*dcos(n * lat) **2 &
    !    - 4._r8 * dsin(lat) * dcos(n*lat)*dsin(n*lat)*n*dcos(lat) &
    !  -12._r8 * dcos(lat)**2*n**2+16._r8*dcos(lat)**2 *dcos(n* lat)**2* n**2)&
    !/dcos(lat)**2


    !TESTE 4 (-lap u +u =f)
    lap_ext = -(-dcos(m*lon)*dcos(n*lat)**2 * (m**2*dcos(n * lat) **2 &
      - 4._r8 * dsin(lat) * dcos(n*lat)*dsin(n*lat)*n*dcos(lat) &
      -12._r8 * dcos(lat)**2*n**2+16._r8*dcos(lat)**2 *dcos(n* lat)**2* n**2)&
      /dcos(lat)**2)   +  dcos(m*lon)*dcos(n*lat)**4

    !------------------------------------------------------------
    !TESTE 5
    !lap_ext=-2*dsin(lat)

    !TESTE 5 (-lap u+u=f)
    !lap_ext=2*dsin(lat)+dsin(lat)
    !----------------------------------------------------

    !TESTE 6
    !lap_ext=(-9*dsin(3*lon)*((dcos(3*lat))**4)/dcos(lat)**2) &
    ! -6*dsin(3*lon)*(dcos(lat))**2*(1-2*dcos(2*lat))**2*(2*dcos(2*lat) &
    !-2*dcos(4*lat)+13*dcos(6*lat)-7)   


    !TESTE 7
    !lap_ext=0!6*(dcos(lat)**3)*(dcos(3*lon))


    !TESTE 8
    !lap_ext=0.0

    return
  end function lap_ext



  !------------------------------------------------------------
  !    Initial condition
  !------------------------------------------------------------

  function inic(p)
    real (r8)::p(1:3)
    real(r8)::inic
    real(r8)::lat
    real(r8)::lon

    call cart2sph (p(1), p(2), p(3), lon, lat)

    inic=dsin(35*lat)+(dcos(10*lon))

    return
  end function inic


  subroutine mesh_generation(meshvet)
    !--------------------------------------------------------------------
    !  MESh_GENERATION 
    !       Construction of meshes  
    !--------------------------------------------------------------------

    !List of meshes
    type(multimesh)::meshvet

    !Couter
    integer(i4)::i

    !Finer grid level
    integer(r4)::levelmf

    character (len=300):: buffer    !Buffer for strings
    character(len=60)::filename     !Parameters file
    integer::fileunit               !Unit for input file

    !Standard mesh parameters file
    filename=trim(pardir)//"multigrid.par"
    call getunit(fileunit)

    !A parameters file already must exist
    open(fileunit,file=filename,status='old')

    read(fileunit,*)  buffer
    read(fileunit,*)  buffer
    read(fileunit,*)  levelmf

    close(fileunit)

    !Allocate meshes
    allocate(meshvet%mesh(0:levelmf))

    do i=0,levelmf
      !Get parameters from file mesh.par
      call getparameters(meshvet%mesh(i))

      if(trim(meshvet%mesh(i)%kind)=="icos")then
        !Overwrite number of nodes
        meshvet%mesh(i)%nv=10*(2**(2*(i)))+2

      elseif(trim(meshvet%mesh(i)%kind)=="read")then
        write(meshvet%mesh(i)%name(10:10),'(i1)')i
      end if
      print*, i,meshvet%mesh(i)%name

      !meshbuild from module smeshpack
      call meshbuild(meshvet%mesh(i))

    enddo



  end subroutine mesh_generation


  subroutine relation_meshes(meshvet)
    !----------------------------------------------------------
    !  RELATION_MESHES
    !  
    !This subroutine relates cells of two consecutive levels
    !----------------------------------------------------------

    !List of meshes
    type(multimesh)::meshvet

    !meshe
    type(grid_structure)::mesh

    !grid level
    integer(i4)::niv

    !finer grid level
    integer(i4)::levelmf

    !levelmf-1
    integer(i4)::level

    !Counter
    integer(i4)::i
    integer(i4)::j
    integer(i4)::k

    !Index of fine mesh
    integer(i4)::indcel

    !The two endpoints of the edge of the triangle
    integer(i4)::i1
    integer(i4)::i2

    !neighboring cells the two endpoints of
    !the edge of the triangle
    integer(i4)::l1
    integer(i4)::l2

    character (len=300):: buffer  !Buffer for strings
    character(len=60)::filename   !Parameters file
    integer::fileunit             !Unit for input file

    !Standard mesh parameters file
    filename=trim(pardir)//"multigrid.par"
    call getunit(fileunit)
    !A parameters file must exist 
    open(fileunit,file=filename,status='old')

    read(fileunit,*)  buffer
    read(fileunit,*)  buffer
    read(fileunit,*)  levelmf
    close(fileunit)


    level=levelmf-1

    !fine mesh
    allocate(meshvet%mf(0:levelmf))

    !coarse mesh
    allocate(meshvet%mg(0:levelmf))

    mesh=meshvet%mesh(levelmf)
    do i=0, levelmf-1
      allocate(meshvet%mf(i+1)%indmf(1:mesh%ne))
      allocate(meshvet%mf(i+1)%cf(1:mesh%nv))
      allocate(meshvet%mg(i)%cng(1:mesh%nv))
      allocate(meshvet%mg(i)%cg(1:mesh%nv))
    enddo


    !For each grid level
    do niv=0,levelmf-1
      mesh=meshvet%mesh(niv)

      !for each edge
      do i=1,mesh%ne

        !Two endpoints of the edge
        i1=mesh%ed(i)%v(1)
        i2=mesh%ed(i)%v(2)

        !for each neighbor of i1
        do k=1,mesh%v(i1)%nnb

          l1=meshvet%mesh(niv+1)%v(i1)%nb(k)

          !for each neighbor of i2
          do j=1, mesh%v(i2)%nnb

            l2=meshvet%mesh(niv+1)%v(i2)%nb(j)

            if(l1==l2) then

              !index of fine cell between two neighboring cells
              indcel=l2

              meshvet%mf(niv+1)%indmf(i)%celm=indcel
              meshvet%mf(niv+1)%indmf(i)%cel(1)=i1
              meshvet%mf(niv+1)%indmf(i)%cel(2)=i2

            endif
          enddo
        enddo
      enddo
    enddo



    do niv=0,levelmf-1

      !for all cells in the coarse level
      do i=1, meshvet%mesh(niv)%nv

        meshvet%mg(niv)%cg(i)=i
      enddo

      !Given a cell in the fine mesh provides index in the coarse mesh

      !for all cells in the fine level
      do k=1, meshvet%mesh(niv+1)%nv
        if (k<=meshvet%mesh(niv)%nv) then
          meshvet%mf(niv+1)%cf(k)%cmg=meshvet%mg(niv)%cg(k)

        else if (k>meshvet%mesh(niv)%nv) then

          EXIT
        endif
      enddo



      !Given a cell in the coarse mesh provides index in the coarse fine
      do k=1, meshvet%mesh(niv)%nv
        meshvet%mg(niv)%cng(k)%cmf=meshvet%mg(niv)%cg(k)
      enddo
    enddo

  end subroutine relation_meshes


  subroutine laplacian(meshvet)
    !-------------------------------------------------------
    !LAPLACIAN
    !Test for laplacian for each point of the mesh
    ! For this test we use the subroutine "lap_cell"
    !-------------------------------------------------------
    !List of meshes
    type(multimesh), intent(in) :: meshvet
    !Mesh
    type(grid_structure)::mesh
    !Exact value of the laplacian on each grid point
    type(scalar_field)::lap_exact
    !Approximate value of the laplacian on each grid point
    type(scalar_field)::lap_ap
    !Function you want to calculate the Laplacian
    type(scalar_field)::fcalc
    !Error between the exact value and the approximate value of the Laplacian
    type(scalar_field)::erro


    integer(i4)::levelmf  !Finer grid level
    integer::niv          !Grid level
    integer::i            !Couter

    !Norm errors on each grid level
    real(r8),allocatable::nrm2(:)
    real(r8),allocatable::nrmmax(:)

    character (len=300):: buffer  !Buffer for strings
    character(len=60)::filename   !Parameters file
    integer::fileunit             !Unit for input file
    integer::iunit                !Aux variable


    !Standard mesh parameters file
    filename=trim(pardir)//"multigrid.par"
    call getunit(fileunit)

    !A parameters file must exist 
    open(fileunit,file=filename,status='old')

    read(fileunit,*)  buffer
    read(fileunit,*)  buffer
    read(fileunit,*)  levelmf
    close(fileunit)

    !Allocate vectors
    allocate(nrm2(0:levelmf))
    allocate(nrmmax(0:levelmf))

    !For each grid level
    do niv=0,levelmf   

      mesh=meshvet%mesh(niv)

      !Allocate exact value
      allocate(lap_exact%f(1:mesh%nv))

      !Allocate approximate value
      allocate(lap_ap%f(1:mesh%nv))

      !Allocate error
      allocate(erro%f(1:mesh%nv))

      !Allocation function to test
      allocate(fcalc%f(1:mesh%nv))

      !Calculate exact value of the laplacian on each grid point
      !Calculate function you want to calculate the Laplacian
      do i=1,mesh%nv
        lap_exact%f(i)=lap_ext(mesh%v(i)%p)
        fcalc%f(i)=func(mesh%v(i)%p)
      enddo

      !Calculate approximate value of the laplacian on each grid point
      call lap_cell(mesh,lap_ap,fcalc)


      !Calculate errors in laplacian estimatives
      nrm2(niv)=error_norm_2(lap_exact%f,lap_ap%f,mesh%nv)
      nrmmax(niv)=error_norm_max(lap_exact%f,lap_ap%f,mesh%nv)

      !Calculate error in laplacian estimative on each grid point
      do i=1,mesh%nv
        erro%f(i)=abs(lap_ap%f(i)-lap_exact%f(i))
      enddo

      !PLOTS

      !erro%name="ERRO_AP_LAP"
      !call plot_scalarfield(erro,meshvet%mesh(niv))

      !Deallocate
      deallocate(lap_exact%f,lap_ap%f,fcalc%f,erro%f)

    enddo

    !Write values on file
    filename=trim(datadir)//"NRM_ERRO_LAP.txt"

    call getunit ( iunit )
    open(iunit, file=filename, status='replace')
    write(iunit,*) "erro_nrm2, erro_nrm_max"
    do i=0,levelmf  
      write(iunit,"(2f16.8)") nrm2(i),nrmmax(i)
    enddo
    close (iunit)

    !Deallocate
    deallocate(nrm2,nrmmax)

    return
  end subroutine laplacian


  subroutine relaxac(meshvet)
    !--------------------------------------------------------
    !RELAXAC
    !
    ! Test for relaxation of the problem using the 
    !  method of Jacobi.
    ! For this test we use the subroutine "relaxac_mult_cell"
    !---------------------------------------------------------

    type(multimesh)::meshvet                 !List of meshes
    type(grid_structure)::mesh               !Mesh
    type(scalar_field)::fexact      !Exact value of the function on each grid point
    type(scalar_field)::fap         !Approximate value of the function on each grid point
    type(scalar_field)::resid       !Residue
    type(scalar_field)::faperro     !Error between the exact value and the approximate value of the function
    type(scalar_field)::g           !Independent term of equation
    type(scalar_field)::lap_ap      !Approximate value of the laplacian on each grid point
    type(interpolable_meshes)::fap_aux       !Auxiliary vector

    !Couter
    integer::i

    integer(i4)::levelmf  !Finer grid level
    integer::niv          !Grid level
    integer(i4)::numit    !Maximum number of iterations

    real (r8):: w         !Relaxation parameter
    real (r8)::nrm_res    !Norm residue
    real (r8)::TOL        !Tolerance


    !Norm errors on each grid level
    real(r8),allocatable::nrm2(:)
    real(r8),allocatable::nrmmax(:)
    real(r8),allocatable::nrmres(:)


    !A parameters file must exist
    character(len=60)::filename
    character (len=300):: buffer
    integer::fileunit
    integer::iunit

    !Standard mesh parameters file
    filename=trim(pardir)//"multigrid.par"
    call getunit(fileunit)
    !A parameters file must exist 
    open(fileunit,file=filename,status='old')

    read(fileunit,*)  buffer
    read(fileunit,*)  buffer
    read(fileunit,*)  levelmf
    read(fileunit,*)  buffer
    read(fileunit,*)  w
    close(fileunit)

    !Allocate vectors
    allocate(nrm2(0:levelmf))
    allocate(nrmmax(0:levelmf))
    allocate(nrmres(0:levelmf))


    TOL=0.00000001             !Tolerance
    nrm_res=100                !Norm residue
    numit=4                   !Maximum number of iterations



    !For each grid level
    do niv=0,levelmf

      mesh=meshvet%mesh(niv)

      !Approximate value of the function on each grid point
      allocate(fap%f(1:mesh%nv))

      !Exact value of the function on each grid point
      allocate(fexact%f(1:mesh%nv))

      !Residue
      allocate(resid%f(1:mesh%nv))

      !Independent term of equation
      allocate(g%f(1:mesh%nv))

      !Auxiliary vector
      allocate(fap_aux%f(1:mesh%nv))

      !Approximate value of the laplacian on each grid point
      allocate(lap_ap%f(1:mesh%nv))

      !Error between the exact value and the approximate value of the function
      allocate(faperro%f(1:mesh%nv))


      !For each grid point
      do i=1,mesh%nv
        fap%f(i)=0!inic(mesh%v(i)%p)                        !initial condition
        !initial condition
        fexact%f(i)=func(mesh%v(i)%p)    !Exact value
        g%f(i)=lap_ext(mesh%v(i)%p)      !Independent term
      enddo



      !Calculate approximate value of the function on each grid point
      call relaxac_mult_cell_GS(mesh,g,fap,numit,w)

      !Calculation of the residue after each iteration
      call lap_cell(mesh,lap_ap,fap)
      resid%f=g%f-lap_ap%f

      !Norm residue
      nrm_res=error_nrm_max(resid%f,mesh%nv)

      !if(nrm_res<=TOL) then
      !  print*, 'Tolerance was achieved'

      !Norm errors on each grid level
      ! nrm2(niv)=error_norm_2(fexact%f,fap%f,mesh%nv)
      !nrmmax(niv)=error_norm_max(fexact%f,fap%f,mesh%nv)
      !nrmres(niv)=error_nrm_max(resid%f,mesh%nv)
      ! print*,nrm2(niv),nrmmax(niv),nrmres(niv)

      !EXIT
      !endif


      !Calculate error in funtion estimative on each grid point
      !do i=1,meshvet%mesh(niv)%nv
      !  faperro%f(i)=abs(fap%f(i)-fexact%f(i))
      !enddo

      !PLOTS
      !faperro%name="ERRO_RELAXAC"
      !call plot_scalarfield(faperro,meshvet%mesh(niv))

      !Norm errors on each grid level

      nrm2(niv)=error_norm_2(fexact%f,fap%f,mesh%nv)
      nrmmax(niv)=error_norm_max(fexact%f,fap%f,mesh%nv)
      nrmres(niv)=nrm_res


      !Deallocate
      deallocate(fap%f,fexact%f,faperro%f,g%f,fap_aux%f,resid%f,lap_ap%f)
    enddo

    !Write values on file
    filename=trim(datadir)//"NRM_ERRO_RELAXAC.txt"
    call getunit ( iunit )
    open(iunit, file=filename, status='replace')
    write(iunit,*) "nrm2, nrmmax"
    do niv=0,levelmf
      write(iunit,"(3f32.16)") nrm2(niv), nrmmax(niv)
    enddo
    close (iunit)

    !Deallocate
    deallocate(nrm2,nrmmax,nrmres)

    return
  end subroutine relaxac

  subroutine interpol_linear(meshvet)
    !--------------------------------------------------------------------
    !  INTERPOL_LINEAR
    ! Test for interpolation linear.
    ! For this test we use the subroutine "interpol_cell"
    !-------------------------------------------------------------------

    !List of meshes
    type(multimesh),intent(in)::meshvet
    !Exact value of the function on each grid point
    type(scalar_field)::fexact
    !Interpolate value
    type(scalar_field)::fap_interpol
    !Auxiliar scalar field for interpolation
    type(scalar_field)::fap
    !Error between the exact value and the interpolate value of the function
    type(scalar_field)::faperro

    !Couter
    integer(i4)::i

    integer(i4)::niv                  !Grid level
    integer(i4)::levelmf              !Finer grid level

    !Norm errors on each grid level
    real(r8),allocatable::nrm2(:)
    real(r8),allocatable::nrmmax(:)

    !A parameters file must exist
    character(len=60)::filename
    character (len=300):: buffer
    integer::fileunit
    integer::iunit

    !Standard mesh parameters file
    filename=trim(pardir)//"multigrid.par"
    call getunit(fileunit)
    !A parameters file must exist 
    open(fileunit,file=filename,status='old')

    read(fileunit,*)  buffer
    read(fileunit,*)  buffer
    read(fileunit,*)  levelmf
    close(fileunit)

    !Allocate vectors
    allocate(nrm2(0:levelmf))
    allocate(nrmmax(0:levelmf))

    !For each grid level
    do niv=0,levelmf-1

      !Allocate
      allocate(fexact%f(1:meshvet%mesh(niv+1)%nv))
      allocate(fap_interpol%f(1:meshvet%mesh(niv+1)%nv))
      allocate(faperro%f(1:meshvet%mesh(niv+1)%nv))
      allocate(fap%f(1:meshvet%mesh(niv+1)%nv))


      !For each grid point
      do i=1, meshvet%mesh(niv)%nv
        !Exact value function
        fap%f(i)=func(meshvet%mesh(niv)%v(i)%p)
      enddo

      !Calculate interpolate value
      call interpol_cell(meshvet,niv,fap,fap_interpol)

      !For each grid point
      do i=1, meshvet%mesh(niv+1)%nv
        !Exact value
        fexact%f(i)=func(meshvet%mesh(niv+1)%v(i)%p)

        !Error between the exact value and the interpolate value
        faperro%f(i)=abs(fap%f(i)-fexact%f(i))
      enddo

      !Norm errors on each grid level
      nrm2(niv)=error_norm_2(fap_interpol%f,fexact%f,meshvet%mesh(niv+1)%nv)
      nrmmax(niv)=error_norm_max(fap_interpol%f,fexact%f,meshvet%mesh(niv+1)%nv)

      !PLOTS
      !faperro%name="PLOT_ERRO_INTERPOL"
      !call plot_scalarfield(faperro,meshvet%mesh(niv+1))

      !Deallocate
      deallocate(fap%f,fexact%f,faperro%f,fap_interpol%f)
    enddo

    !Write values on file
    filename=trim(datadir)//"NRM_ERRO_INTERPOL.txt"
    call getunit ( iunit )
    open(iunit, file=filename, status='replace')
    write(iunit,*) "erro_nrm2, erro_nrm_max"
    do i=0,levelmf-1
      write(iunit,"(2f16.8)") nrm2(i),nrmmax(i)
    enddo
    close (iunit)
  end subroutine interpol_linear


  subroutine multigrid_solver(meshvet)
    !=============================================================================
    ! multigrid_solver
    !
    ! For this multigrid solver we use the subroutines
    !
    ! "lap_cell" - for calculate laplacian
    ! "relaxac_mult_cell" - for relaxation
    ! "interpol_cell" - for interpolation
    ! "transfer_res" - for transfer residue
    ! "decomLU" - for resolution system coarse mesh
    !
    !=============================================================================
    type(multimesh)::meshvet                     !List of meshes
    type(grid_structure)::mesh                   !mesh
    type(multisol)::meshsol                      !approximate solution in each level 
    type(interpolable_meshes)::fap_aux           !Auxiliary vector
    type(scalar_field)::fap             !Approximate value of the function on each grid point
    type(scalar_field)::residuoa        !Residue before the cycle
    type(scalar_field):: residuod       !Residue after the cycle
    type(scalar_field)::g               !Independent term of equation
    type(scalar_field)::lapc            !Approximate value of the laplacian on each grid point
    type(scalar_field)::fap_interpol    !Interpolate value
    type(scalar_field)::fap_exact       !Exact function 
    type(scalar_field)::error           !error=fap_exact-fap

    integer::niv           !Grid level
    integer(i4)::levelmf   !Finer grid level

    !counters
    integer::i
    integer::j
    integer::k

    integer::contv       !Couter cycle
    integer(i4)::numvc   !Number cycle

    !number relaxation before the cycle
    integer(i4)::n1

    !number relaxation after the cycle
    integer(i4)::n2

    !gamma=1 v-ciclo/ gamma=2 W-ciclo
    integer::gamma

    !parameter of relaxation
    real(r8)::w

    !nrm residue before the cycle
    real(r8)::nrmg1

    !nrm residue after the cycle
    real(r8)::nrmg2

    !Smoothness factor - (nrmg2/nrmg1)**(1./(n1+n2))
    real(r8)::sfact

    real (r8)::Lchx      !Distance between two neighboring
    real(r8)::ledhx      !Edge length
    real(r8)::laptmp

    !Auxiliary vector
    real(r8), allocatable::c(:)

    !A parameters file must exist
    character (len=300):: buffer
    character(len=60)::filename
    integer::fileunit


    !Standard mesh parameters file
    filename=trim(pardir)//"multigrid.par"
    call getunit(fileunit)
    !A parameters file must exist 
    open(fileunit,file=filename,status='old')

    read(fileunit,*)  buffer
    read(fileunit,*)  buffer
    read(fileunit,*)  levelmf
    read(fileunit,*)  buffer
    read(fileunit,*)  w
    read(fileunit,*)  buffer
    read(fileunit,*)  numvc
    read(fileunit,*)  buffer
    read(fileunit,*)  n1
    read(fileunit,*)  n2
    close(fileunit)


    allocate(meshsol%solm(0:levelmf))
    allocate(c(0:levelmf))

    do niv=0,levelmf

      !mesh level "niv"
      mesh=meshvet%mesh(niv)

      !value function each grid level
      allocate(meshsol%solm(niv)%valf(1:mesh%nv))

      !value residue each grid level
      allocate(meshsol%solm(niv)%resf(1:mesh%nv))
    enddo

    mesh=meshvet%mesh(levelmf)

    allocate(lapc%f(1:mesh%nv))
    allocate(g%f(1:mesh%nv))
    allocate(fap%f(1:mesh%nv))
    allocate(fap_aux%f(1:mesh%nv))
    allocate(fap_interpol%f(1:mesh%nv))
    allocate(fap_exact%f(1:mesh%nv))
    allocate(error%f(1:mesh%nv))

    allocate(residuoa%f(1:mesh%nv))
    allocate(residuod%f(1:mesh%nv))

    do i=1,mesh%nv

      !Independent term
      g%f(i)=lap_ext(mesh%v(i)%p)

      !initial condition
      fap%f(i)=inic(mesh%v(i)%p)

      !Exact value
      fap_exact%f(i)=func(mesh%v(i)%p)

    enddo

    !--------------------------------------------------
    !PLOTS

    !Exact function
    !fap_exact%name="PLOT_FUNC_EXACT"
    !call plot_scalarfield(fap_exact,meshvet%mesh(levelmf))

    !inicial condition
    !fap%name="PLOT_INIC_COND"
    !call plot_scalarfield(fap,meshvet%mesh(levelmf))
    !-------------------------------------------------------  

    gamma=1 !gamma=1 vciclo,  gamma=2 wciclo
    contv=1

        !Begin V-ciclo

    10  do while(contv<=numvc)

      do  k=0,levelmf
        c(k)=0
      enddo

      k=levelmf

      mesh=meshvet%mesh(levelmf)

      !residue before the cycle
      call lap_cell(mesh,lapc,fap)

      !For each grid point
      do i=1,mesh%nv
        residuoa%f(i)=g%f(i)-lapc%f(i)
      enddo

      !----------------------------------------------------------------------
      !                (I) relax the values ​​in the fine mesh
      !----------------------------------------------------------------------
20    mesh=meshvet%mesh(k)

      !call relaxac_mult_cell(mesh,g,fap,fap_aux,n1,w)
      call relaxac_mult_cell_GS(mesh,g,fap,n1,w)

      c(k)=c(k)+1

      !Save the solution to use later
      meshsol%solm(k)%valf=fap%f

      !-----------------------------------------------------------------------
      !            (II) calculation of the residue
      !-----------------------------------------------------------------------

      call lap_cell(mesh,lapc,fap)
      do i=1,mesh%nv
        g%f(i)=g%f(i)-lapc%f(i)
      enddo

      !------------------------------------------------------------------------
      !  (III) Transference of the residue to the coarser grid
      !
      !Restriction by:
      !Guohua Zhou and Scoot R. Fulton, Fourier analysis of
      !multigrid methods on hexagonal grids,
      !Society for Industrial and Applied Mathematics. (31):1518-1538, 2009.
      !------------------------------------------------------------------------

      call transfer_res(meshvet,g,k)

      k=k-1

      mesh=meshvet%mesh(k)
      do i=1,mesh%nv
        meshsol%solm(k)%resf(i)=g%f(i)
      enddo


      if (k.NE.0) then
        do i=1,mesh%nv
          fap%f(i)=0.0
        enddo

        goto 20 !relax the values
      endif

      !--------------------------------------------------------
      !    (IV)   solve the problem in the coarse grid
      !-------------------------------------------------------

      ! matrix problem in coarser level of resolution

      if(k.EQ.0) then

        do i=1,mesh%nv
          fap%f(i)=0.0
        enddo

        !call relaxac_mult_cell(mesh,g,fap,fap_aux,15,w)
        call relaxac_mult_cell_GS(mesh,g,fap,15,w)
      endif

      !---------------------------------------------------------------------
      !            (V)nterpolate the correction for fine mesh
      !                 and add the solution obtained in the step I
      !---------------------------------------------------------------------

40    call interpol_cell(meshvet,k,fap,fap_interpol)

      k=k+1
      mesh=meshvet%mesh(k)
      do i=1,mesh%nv
        fap%f(i)=fap_interpol%f(i)+meshsol%solm(k)%valf(i)
      enddo

      !---------------------------------------------------------------------
      !             (VI) apply relaxations in the previous solution
      !---------------------------------------------------------------------

      if (k==levelmf) then
        do i=1,mesh%nv
          g%f(i)=lap_ext(mesh%v(i)%p)
        enddo
      else
        do i=1,mesh%nv
          g%f(i)=meshsol%solm(k)%resf(i)
        enddo
      endif


      !call relaxac_mult_cell(mesh,g,fap,fap_aux,n2,w)
      call relaxac_mult_cell_GS(mesh,g,fap,n2, w)

      if(k.NE.levelmf.AND.c(k)==gamma) then
        c(k)=0
        goto 40

      else if (k.NE.levelmf.AND.c(k).NE.gamma) then
        goto 20
      else if  (k.EQ.levelmf) then
        contv=contv+1

        !residue after the cycle V-CICLO
        call lap_cell(mesh,lapc,fap)

        !For each grid point
        do i=1,mesh%nv
          residuod%f(i)=g%f(i)-lapc%f(i)
        enddo

        !calculation of the convergence factor
        nrmg1=error_norm_g(residuoa%f,mesh%nv)
        nrmg2=error_norm_g(residuod%f,mesh%nv)

        sfact=(nrmg2/nrmg1)

        print*,'Ciclo',contv-1
        print*, 'convergence factor:', sfact

      endif

       
      if(contv==numvc) then
        goto 60
      endif

      goto 10

    enddo

60  print*,''

    do i=1,meshvet%mesh(levelmf)%nv
      error%f(i)=abs(fap%f(i)-fap_exact%f(i))
    enddo

    !PLOTS
    !Approximate solution
    !fap%name="PLOT_SOL"
    !call plot_scalarfield(fap,meshvet%mesh(levelmf))

    !Error
    !error%name="PLOT_ERRO"
    !call plot_scalarfield(error,meshvet%mesh(levelmf))

    print*,''
    print*,'nrm2, nrmmax'
    print*, error_norm_2(fap%f,fap_exact%f,meshvet%mesh(levelmf)%nv), error_norm_max(fap%f,fap_exact%f,meshvet%mesh(levelmf)%nv)


50  return

  end subroutine multigrid_solver



  !==========================================================================
  ! Subroutines for relaxation, interpolation and transfer residue
  !==========================================================================

  subroutine lap_cell(mesh,lapc,fap)
    !------------------------------------------------------------
    !LAP_CELL
    !calculation of the Laplacian each point of mesh.
    !------------------------------------------------------------

    type(grid_structure)::mesh            !mesh
    !value function
    type(scalar_field)::fap
    !Approximate value of the laplacian on each grid point
    type(scalar_field)::lapc

    !counter
    integer::i
    integer::j
    integer::k

    !number points of mesh
    integer::n

    real (r8):: ledhx    !Edge length
    real (r8)::Lchx      !Distance between two neighboring
    real (r8)::laptmp


    !number points of the mesh
    n=mesh%nv

    !For each grid point
    do  i=1,n

      laptmp=0_r8
      ledhx=0_r8
      Lchx=0_r8

      !For each grid neighboring of i
      do j=1, mesh%v(i)%nnb

        !Edge index
        k=mesh%v(i)%ed(j)

        !Edge length
        ledhx=mesh%edhx(k)%leng

        !Distance between two neighboring
        Lchx=mesh%ed(k)%leng

        laptmp=laptmp+(ledhx/Lchx)*(fap%f(mesh%v(i)%nb(j))-fap%f(i))

      enddo


      !para o problema -lap f+f =g
      lapc%f(i)=-(laptmp/mesh%hx(i)%areag)+fap%f(i)

       !para o problema lap f =g
       !lapc%f(i)=laptmp

    enddo
    return 
  end subroutine lap_cell


  subroutine relaxac_mult_cell(mesh,g,fap,fap_aux,numit, w)
    !------------------------------------------------------------
    !RELAXAC_MULT_CELL
    ! Relaxation Scheme - Weighted Jacobi
    !------------------------------------------------------------
    type(grid_structure)::mesh             !Mesh
    type(interpolable_meshes)::fap_aux     !Auxiliary vector
    !Approximate value of the function on each grid point
    type(scalar_field)::fap
    !Independent term of equation
    type(scalar_field)::g

    !parameter of relaxation
    real(r8)::w

    !Conter
    integer::k
    integer::i
    integer::j

    !Number of Iterations
    integer::numit

    !number of grid points
    integer::n

    real (r8):: ledhx     !Edge length
    real (r8)::Lchx       !Distance between two neighboring
    real (r8)::laptmp
    real (r8)::laptmp1

    !number points of the mesh
    n=mesh%nv

    !For each iteration
    do k=1,numit

      !For each grid point
      do i=1,n

        laptmp=0_r8
        laptmp1=0_r8
        ledhx=0_r8
        Lchx=0_r8

        !For each grid neighboring of i
        do j=1,mesh%v(i)%nnb

          !Edge length
          ledhx=mesh%edhx(mesh%v(i)%ed(j))%leng

          !Distance between two neighboring
          !arclen(mesh%v(mesh%v(i)%nb(j))%p,mesh%v(i)%p)
          Lchx=mesh%ed(mesh%v(i)%ed(j))%leng

          laptmp1=laptmp1+(ledhx/Lchx)
          laptmp=laptmp+(ledhx/Lchx)*(fap%f(mesh%v(i)%nb(j)))

        enddo

        laptmp=laptmp/mesh%hx(i)%areag
        laptmp1=laptmp1/mesh%hx(i)%areag

        !problem -lap f+f =g
        fap_aux%f(i)=(1-w)*(fap%f(i))+w*((laptmp+g%f(i))/(laptmp1+1))


         !problem lap f =g
         !fap_aux%f(i)=(1-w)*fap%f(i)+w*((laptmp-g%f(i))/(laptmp1))

      enddo

      fap%f=fap_aux%f

    enddo

    return 
  end subroutine relaxac_mult_cell

  subroutine relaxac_mult_cell_GS(mesh,g,fap,numit, w)
    !------------------------------------------------------------
    !RELAXAC_MULT_CELL
    ! Relaxation Scheme - Weighted Gauss-Seidel
    !------------------------------------------------------------
    type(grid_structure)::mesh  !Mesh
    type(scalar_field)::fap     !Approximate value of the function on each grid point
    type(scalar_field)::g       !Independent term of equation

    !parameter of relaxation
    real(r8)::w

    !Conter
    integer::k
    integer::i
    integer::j

    !Number of Iterations
    integer::numit

    !number of grid points
    integer::n

    real (r8):: ledhx     !Edge length
    real (r8)::Lchx       !Distance between two neighboring
    real (r8)::laptmp
    real (r8)::laptmp1

    !number points of the mesh
    n=mesh%nv

    !For each iteration
    do k=1,numit

      !For each grid point
      do i=1,n

        laptmp=0_r8
        laptmp1=0_r8
        ledhx=0_r8
        Lchx=0_r8

        !For each grid neighboring of i
        do j=1,mesh%v(i)%nnb

          !Edge length
          ledhx=mesh%edhx(mesh%v(i)%ed(j))%leng

          !Distance between two neighboring
          Lchx=mesh%ed(mesh%v(i)%ed(j))%leng
          laptmp1=laptmp1+(ledhx/Lchx)
          laptmp=laptmp+(ledhx/Lchx)*(fap%f(mesh%v(i)%nb(j)))
        enddo

        laptmp=laptmp/mesh%hx(i)%areag
        laptmp1=laptmp1/mesh%hx(i)%areag

        !para o problema lap f =g
        !fap%f(i)=(1-w)*fap%f(i)+w*((laptmp-g%f(i))/(laptmp1))


        !para o problema -lap f+f =g
        fap%f(i)=(1-w)*(fap%f(i))+w*((laptmp+g%f(i))/(laptmp1+1))


      enddo

    enddo

    return 
  end subroutine relaxac_mult_cell_GS


  subroutine interpol_cell(meshvet,niv,fap,fap_interpol)
    !---------------------------------------------------------------------
    !INTERPOL_CELL
    ! Linear Interpolation
    !The center of cell on the fine grid coincide with centers of cell on 
    !the coarser grid remain at the same value. For interpolated values
    !using linear interpolation between two cells neighboring. 
    !----------------------------------------------------------------------

    type(multimesh)::meshvet            !list of meshes
    type(scalar_field)::fap_interpol    !Interpolate value
    !Auxiliar interlation field
    type(scalar_field)::fap

    !neighboring cells the two endpoints of
    !the edge of the triangle
    integer(i4)::l1    
    integer(i4)::l2    

    integer::i          !counter
    integer(i4)::niv    !Grid level

    integer::n          !Number of points
    integer::na         !number of edges of the triangle

    n=meshvet%mesh(niv)%nv
    na=meshvet%mesh(niv)%ne

    !values ​​remain the same
    do i=1, n
      fap_interpol%f(i)=fap%f(i)
    enddo

    !interpolated values
    do i=1,na
      l1=meshvet%mf(niv+1)%indmf(i)%cel(1)
      l2=meshvet%mf(niv+1)%indmf(i)%cel(2)

      fap_interpol%f(meshvet%mf(niv+1)%indmf(i)%celm)=(fap%f(l1)+fap%f(l2))/2

    enddo

    return
  end subroutine interpol_cell


  subroutine transfer_res(meshvet,residuo,level)
    !-----------------------------------------------------------
    !TRANSFER_RES
    !Transfer the residue of the finer mesh to coarser grid
    !Full weighting
    !-----------------------------------------------------------

    type(multimesh)::meshvet                 !List of meshes
    type(scalar_field)::residuo     !Residue

    !counter
    integer::i
    integer::j
    integer::l

    integer::level        !lever grid

    !Auxiliary
    real(r8)::laptmp       
    real(r8)::laptmp1   

    integer::k

    do i=1,meshvet%mesh(level-1)%nv

      laptmp=0_r8
      laptmp1=0_r8


      !Injeção
      !residuo%f(i)=residuo%f(meshvet%mf(level)%cf(i)%cmg)


      !Full weighting
      do j=1, meshvet%mesh(level)%v(i)%nnb

        !neighbor index
        l=meshvet%mesh(level)%v(i)%nb(j)

        !k=meshvet%mg(level-1)%cng(l)%cmf

        laptmp=laptmp+(residuo%f(l)*meshvet%mesh(level)%hx(l)%areag)
        laptmp1= laptmp1+meshvet%mesh(level)%hx(l)%areag

      enddo

      residuo%f(i)=(meshvet%mesh(level)%hx(i)%areag*residuo%f(i)+0.5*laptmp) &
        /(meshvet%mesh(level)%hx(i)%areag+0.5*laptmp1)
    enddo

    return
  end subroutine transfer_res


  function error_nrm_max(f, n)
    !-------------------------------------------
    !Calcula o maximo valor absoluto
    !  do vetor f de tamanho n
    !-------------------------------------------
    integer (i4):: n
    real (r8), dimension(1:n), intent(in) :: f
    real (r8):: error_nrm_max

    error_nrm_max=maxval(abs(f(1:n)))

    return
  end function error_nrm_max


  function error_norm_g(f,n)
    !-------------------------------------------
    !Calcula a raiz quadrada da soma
    !  dos elementos do vetor f
    !-------------------------------------------
    integer (i4):: n
    real (r8), dimension(1:n), intent(in) :: f
    real (r8):: error_norm_g
    real (r8):: sum_sq
    integer::i

    sum_sq=0.0

    do i=1,n 
      sum_sq=sum_sq+(f(i)**2)
    enddo

    error_norm_g=dsqrt(sum_sq)

    return
  end function error_norm_g


end module multigrid

