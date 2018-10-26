module swm_operators
  !=============================================================================
  !  Global data for shallow water model
  !
  ! Pedro da Silva Peixoto (pedrosp@ime.usp.br)
  ! Oct 2018
  !=============================================================================

   !Use global constants and kinds
  use constants !Everything

  !Global variables for shallow water model
  use swm_data !Everything

  !Use main grid data structures
  use datastruct, only: &
    grid_structure, &
    scalar_field, &
    vectorinterpol_methods, &
    vector_field_cart

  !Use routines from the spherical mesh pack
  use smeshpack, only: &
    alignind, &
    alignindlimit, &
    arcdistll, &
    arclen, &
    bar_coord_tr, &
    calc_tiled_areas, &
    convert_vec_sph2cart, &
    cross_product, &
    error_norm_2, &
    error_norm_max, &
    error_norm_max_rel, &
    getunit, &
    getedindexonhx, &
    insidetrmesh, &
    norm, &
    sph2cart, &
    trhx_intersec_areas, &
    cart2sph, &
    proj_vec_sphere, &
    sphtriarea

  !Use interpolation routines
  use interpack, only: &
    calc_trisk_weights, &
    gradcalc, &
    plot_cart_vectorfield, &
    plot_scalarfield, &
    precalcgradpol, &
    scinterpol_linear_trv, &
    wachspress_coords_hxv, &
    scinterpol_wachspress, &
    trsk_order_index, &
    vector_interpol, &
    vector_field_cart, &
    vrec_remap, &
    vector_reconstruct

  !Use differential operator routines
  use diffoperpack, only: &
    div_cell_Cgrid, &
    grad_edge_Cgrid

  implicit none

contains

  subroutine scalar_hx2ed(fhx, fed, mesh)
    !---------------------------------------------------------------
    !Interpolate from cell centers to edges (linear interpolation)
    ! in: fhx - scalar field defined at cells
    ! out: fed - scalar field defined at edges (must already be allocated)
    !---------------------------------------------------------------

    type(grid_structure), intent(in) :: mesh
    type(scalar_field), intent(in):: fhx
    type(scalar_field), intent(inout) :: fed

    integer(i4)::  i, l

    !Quick check for dimensions
    if(fed%n/=mesh%ne)then
      print*, "ERROR in scalar_hx2ed: dimensions do not match", mesh%ne, fed%n
      stop
    end if


    if(useSinterpolTrisk .or. useSinterpolGass .or. useStagHTC)then
      !$omp parallel do default(none) &
      !$omp shared(mesh, fhx, fed) &
      !$omp schedule(static)
      do l=1,mesh%ne
        !Simple average
        fed%f(l)=(fhx%f(mesh%edhx(l)%sh(1))+fhx%f(mesh%edhx(l)%sh(2)))*0.5_r8
         ! h*u is obtained at edges and stored
         !uh%f(l)=u%f(l)*h_ed%f(l)
      end do
       !$omp end parallel do

    elseif(useSinterpolBary .and. useStagHC)then
      !$omp parallel do default(none) &
      !$omp shared(mesh, fhx, fed ) &
      !$omp schedule(static)
      do l=1,mesh%ne
        !Linear interpolation
        fed%f(l)=0._r8 !scinterpol_linear_trv(mesh%edhx(l)%c%p, h, mesh)
        do i=1,3
          fed%f(l)=fed%f(l)+ mesh%edhx(l)%c%b(i)*fhx%f(mesh%tr(mesh%edhx(l)%c%kt)%v(i))
        end do
         ! h*u is obtained at edges and stored
         !uh%f(l)=u%f(l)*h_ed%f(l)
      end do
       !$omp end parallel do

    else
      print*, "ERROR in scalar_hx2ed: don't know what to use"
    end if

    return

  end subroutine scalar_hx2ed

  subroutine scalar_tr2ed(ftr, fed, mesh)
    !---------------------------------------------------------------
    !Interpolate from triang centers to edges (linear interpolation)
    ! in: ftr - scalar field defined at triangles
    ! out: fed - scalar field defined at edges (must already be allocated)
    !---------------------------------------------------------------
    type(grid_structure), intent(in) :: mesh
    type(scalar_field), intent(in):: ftr
    type(scalar_field), intent(inout) :: fed

    integer(i4)::  i, l
    real(r8):: d1, d2

    !Quick check for dimensions
    if(fed%n/=mesh%ne)then
      print*, "ERROR in scalar_tr2ed: dimensions do not match", mesh%ne, fed%n
      stop
    end if

    !$omp parallel do &
    !$omp default(none) &
    !$omp shared(mesh,  fed, ftr) &
    !$omp shared(useSinterpolTrisk, useSinterpolBary) &
    !$omp shared(useStagHTC, useStagHC) &
    !$omp private(d1, d2) &
    !$omp schedule(static)
    do l=1,mesh%ne
      !Calculate PV on edges
      if(useStagHC.or. useSinterpolTrisk)then
        fed%f(l)=0.5_r8*(ftr%f(mesh%ed(l)%sh(1))+ftr%f(mesh%ed(l)%sh(2)))
      else
        d1=arclen(mesh%ed(l)%c%p,mesh%tr(mesh%ed(l)%sh(1))%c%p)
        d2=arclen(mesh%ed(l)%c%p,mesh%tr(mesh%ed(l)%sh(2))%c%p)
        !print*, d1/(d1+d2), d2/(d1+d2)
        fed%f(l)=(d2*ftr%f(mesh%ed(l)%sh(1))+d1*ftr%f(mesh%ed(l)%sh(2)))/(d1+d2)
      end if
    end do
    !$omp end parallel do

    return

  end subroutine scalar_tr2ed



  subroutine scalar_hx2trcc(fhx, ftr, mesh)
    !---------------------------------------------------------------
    !Interpolate from cell centers to triangle centers (constant or linear interpolation)
    ! in: fhx - scalar field defines at cells
    ! out: ftr - scalar field defined at triangles (must already be allocated)
    !---------------------------------------------------------------
    !Grid structure
    type(grid_structure), intent(in) :: mesh
    type(scalar_field), intent(in):: fhx
    type(scalar_field), intent(inout) :: ftr

    integer(i4):: i, k
    real(r8):: signcor

    !Quick check for dimensions
    if(ftr%n/=mesh%nt)then
      print*, "ERROR in scalar_hx2trcc: dimensions do not match", mesh%nt, ftr%n
      stop
    end if

    !$omp parallel do &
    !$omp default(none) &
    !$omp shared(mesh, fhx, ftr) &
    !$omp shared(useSinterpolTrisk, useSinterpolBary, useSinterpolGass) &
    !$omp private(i) &
    !$omp schedule(static)
    do k=1, mesh%nt
      !Interpolate the thickness to triangle cc
      ftr%f(k)=0.
      if(useSinterpolTrisk)then
        !Use area weighted interpolation - 1st order
        do i=1,3 !for each tr vertex
          ftr%f(k)=ftr%f(k)+mesh%tr(k)%trhx_area(i)*fhx%f(mesh%tr(k)%v(i))
        end do
        ftr%f(k)=ftr%f(k)/mesh%tr(k)%areag

      elseif(useSinterpolBary)then
        !Interpolate
        do i=1,3
          ftr%f(k)=ftr%f(k)+ mesh%tr(k)%c%b(i)*fhx%f(mesh%tr(k)%v(i))
        end do

      elseif(useSinterpolGass)then
        print*, "Error at scalar_hx2trcc: Gassmann's scheme needs edge values"
        stop
      end if
    enddo
    !$omp end parallel do

    return

  end subroutine scalar_hx2trcc

  subroutine scalar_ed2trcc(fed, ftr, mesh)
    !---------------------------------------------------------------
    !Interpolate from cell edges to triangle centers (constant or linear interpolation)
    ! in: fed - scalar field defines at edges
    ! out: ftr - scalar field defined at triangles (must already be allocated)
    !---------------------------------------------------------------
    !Grid structure
    type(grid_structure), intent(in) :: mesh
    type(scalar_field), intent(in):: fed
    type(scalar_field), intent(inout) :: ftr

    real(r8):: p1(1:3)
    real(r8):: p2(1:3)
    real(r8):: p3(1:3)

    integer(i4):: i, k, l, ed
    real(r8):: signcor, ed_area

    !Quick check for dimensions
    if(ftr%n/=mesh%nt)then
      print*, "ERROR in scalar_hx2trcc: dimensions do not match", mesh%nt, ftr%n
      stop
    end if

    if(.not.useSinterpolGass)then
      print*, "Warning in scalar_ed2trcc: this should only be required in Gassmann's scheme"
    endif

    !$omp parallel do &
    !$omp default(none) &
    !$omp shared(mesh, fed, ftr, useTiledAreas) &
    !$omp private(i, ed, p1, p2, p3, ed_area) &
    !$omp schedule(static)
    do k=1, mesh%nt
      !Interpolate the thickness to triangle cc
      ftr%f(k)=0.
      !Use area weighted interpolation based on edge values - 1st order
      do i=1,3 !for each tr edge
        !Calculate weight
        ed=mesh%tr(k)%ed(i)
        if(useTiledAreas)then ! Dist(trcc,tredmidpoint)*tred_length/2
          ed_area=0.5*arclen(mesh%tr(k)%c%p, mesh%ed(ed)%c%p)*mesh%ed(ed)%leng
        else !Use geodesic areas !Area of tr formed by trcc and two endpoints of edge
          p1=mesh%tr(k)%c%p
          p2=mesh%v(mesh%ed(ed)%v(1))%p
          p3=mesh%v(mesh%ed(ed)%v(2))%p
          ed_area=sphtriarea(p1, p2, p3)
        endif
        ftr%f(k)=ftr%f(k)+ed_area*fed%f(ed)
      end do
      if(useTiledAreas)then
        ftr%f(k)=ftr%f(k)/mesh%tr(k)%areat
      else
        ftr%f(k)=ftr%f(k)/mesh%tr(k)%areag
      end if
    enddo
    !$omp end parallel do

    return

  end subroutine scalar_ed2trcc

  subroutine vector_edhx2tr_perot(ued, vtr, mesh)
    !---------------------------------------------------------------
    !Reconstruct from cell edges to triangle centers (perot)
    ! in: ued - vector field defined at cell edges with normal components only
    ! out: vtr - vector field defined at triangles (must already be allocated)
    !---------------------------------------------------------------

    type(grid_structure), intent(in) :: mesh
    type(scalar_field), intent(in):: ued
    type(vector_field_cart), intent(inout) :: vtr

    integer(i4):: i, k, l, ed
    real(r8):: signcor, ed_area

    !Quick check for dimensions
    if(vtr%n/=mesh%nt)then
      print*, "ERROR in vector_ed2tr_perot: dimensions do not match", mesh%nt, vtr%n
      stop
    end if

    if(.not.(useCoriolisMtdDtred .or. test_lterror==1))then
      print*, "Warning at vector_ed2tr_perot: this should only be required if useCoriolisMtdDtred or test_lterror selected"
    end if

    !$omp parallel do &
    !$omp default(none) &
    !$omp shared(mesh, ued, vtr) &
    !$omp shared(useStagHTC, useStagHC, useTiledAreas) &
    !$omp private(l, signcor, ed) &
    !$omp schedule(static)
    do k=1, mesh%nt

      !Reconstruct vhq to triangle circumcenter (use Perot's method)
      !  Only used in dtred Coriolis method

      vtr%p(k)%v=0._r8
      do l=1,3
        ed=mesh%tr(k)%ed(l)
        if(useStagHTC)then
          signcor=dsign(1._r8, real(mesh%tr(k)%tg(l), r8))
          !print*, j, ed, signcor, uh%f(ed)
          vtr%p(k)%v=vtr%p(k)%v + &
            signcor*(ued%f(ed))*(mesh%ed(ed)%c%p-mesh%tr(k)%c%p)*mesh%ed(ed)%leng
        elseif(useStagHC)then
          signcor=dsign( 1._r8, real(mesh%tr(k)%tg(l), r8)* &
            dot_product(mesh%edhx(ed)%nr, mesh%ed(ed)%tg))
          !print*, j, ed, signcor, uh%f(ed)
          vtr%p(k)%v=vtr%p(k)%v + &
            signcor*(ued%f(ed))*(mesh%edhx(ed)%c%p-mesh%tr(k)%c%p)*mesh%ed(ed)%leng
        endif
         !print*,uhq_tr%p(k)%v, norm(uhq_tr%p(k)%v)
      end do
      vtr%p(k)%v=cross_product(mesh%tr(k)%c%p, vtr%p(k)%v)/mesh%tr(k)%areag
          !vtr%p(k)%v=vtr%p(k)%v !*q_tr%f(k)
    end do
    !$omp end parallel do

    return

  end subroutine vector_edhx2tr_perot

  subroutine vector_edtr2tr_perot(ued, vtr, mesh)
    !---------------------------------------------------------------
    !Reconstruct from cell edges to triangle centers (perot)
    ! in: ued - vector field defined at triangle edges with normal components only
    ! out: vtr - vector field defined at triangles (must already be allocated)
    !---------------------------------------------------------------

    type(grid_structure), intent(in) :: mesh
    type(scalar_field), intent(in):: ued
    type(vector_field_cart), intent(inout) :: vtr

    integer(i4):: i, k, l, ed
    real(r8):: signcor, ed_area

    !Quick check for dimensions
    if(vtr%n/=mesh%nt)then
      print*, "ERROR in vector_ed2tr_perot: dimensions do not match", mesh%nt, vtr%n
      stop
    end if

    !$omp parallel do &
    !$omp default(none) &
    !$omp shared(mesh, ued, vtr) &
    !$omp private(l, signcor, ed) &
    !$omp schedule(static)
    do k=1, mesh%nt
      vtr%p(k)%v=0._r8
      do l=1,3
        ed=mesh%tr(k)%ed(l)
        signcor=dsign(1._r8, real(mesh%tr(k)%nr(l), r8))
        !print*, j, ed, signcor, uh%f(ed)
        vtr%p(k)%v=vtr%p(k)%v + &
          signcor*(ued%f(ed))*(mesh%ed(ed)%c%p-mesh%tr(k)%c%p)*mesh%ed(ed)%leng
         !print*,ed, q_grad_ed%f(ed)
      end do
      vtr%p(k)%v=vtr%p(k)%v/mesh%tr(k)%areag !/erad !(TODO: check this erad)
       !print*, k, norm(vtr%p(k)%v), norm(vtr%p(k)%v)/erad
       !print*, k, q_grad_tr%p(k)%v
    end do
       !$omp end parallel do

    return

  end subroutine vector_edtr2tr_perot

  subroutine vort_tr(u, eta, mesh)
    !---------------------------------------------------------------
    !Calculate fields at triangles
    !   vorticity at triangle cc based on edge velocities
    !---------------------------------------------------------------
    type(grid_structure), intent(in) :: mesh
    type(scalar_field), intent(in):: u ! velocity at cell edges
    type(scalar_field), intent(inout):: eta !absolute vorticity at tr cc

    integer(i4):: k, l, ed
    real(r8):: signcor

    !OPENMP PARALLEL DO
    !$omp parallel do &
    !$omp default(none) &
    !$omp shared(mesh, u, eta) &
    !$omp shared(useStagHTC, useStagHC, useTiledAreas) &
    !$omp shared(fsphere, fcte) &
    !$omp private(l, ed, signcor) &
    !$omp schedule(static)
    do k=1, mesh%nt
      !print*, k, omp_get_thread_num()
      !First calculate relative vorticity
      eta%f(k)=0.
      !print*, k
      !loop over triangle edges
      do l=1, 3
        ed=mesh%tr(k)%ed(l)
        !if(stag=="HTC")then
        if(useStagHTC)then
          signcor=mesh%tr(k)%tg(l)
          eta%f(k)=eta%f(k)+u%f(ed)*mesh%ed(ed)%leng*signcor
           !print*, eta%f(k)
        elseif(useStagHC)then
          signcor=dsign( 1._r8, real(mesh%tr(k)%tg(l), r8)* &
            dot_product(mesh%edhx(ed)%nr, mesh%ed(ed)%tg))
          eta%f(k)=eta%f(k)+u%f(ed)*mesh%ed(ed)%leng*signcor
        end if
      end do
      !This is the relative vort.
      if(useTiledAreas)then
        eta%f(k)=eta%f(k)/mesh%tr(k)%areat/erad
      else
        eta%f(k)=eta%f(k)/mesh%tr(k)%areag/erad
      endif
      !Add coriolis term - this is the absolute vorticity
      if(fsphere==0)then !variable f
        eta%f(k)=eta%f(k) +2.*Omega*dsin(mesh%tr(k)%c%lat)
         !For the linearized equation, there is only the coriolis effect
      elseif(fsphere==1)then !Just put the linear term - f sphere
        ! Set to constant f
        eta%f(k)=  fcte !1.4584e-4_r8 !+2.*Omega*dsin(mesh%tr(k)%c%lat)
      elseif(fsphere==2)then !nonlinear with constant fm - f sphere
        eta%f(k)=eta%f(k)+  fcte
      end if
    end do
    !$omp end parallel do

    return

  end subroutine vort_tr

  subroutine ke_tr(u, ke, mesh)
    !---------------------------------------------------------------
    !Calculate kinetic energy at triangles
    ! in: u at cell edges
    ! out: ke at triangles
    !---------------------------------------------------------------
    type(grid_structure), intent(in) :: mesh
    type(scalar_field), intent(in):: u ! velocity at cell edges
    type(scalar_field), intent(inout):: ke !kinetic energy at triangles

    integer(i4):: k, l, ed
    real(r8)::ed_area, cell_area

    real(r8):: p1(1:3)
    real(r8):: p2(1:3)
    real(r8):: p3(1:3)

    !$omp parallel do &
    !$omp default(none) &
    !$omp shared(mesh, u, ke) &
    !$omp shared(useTiledAreas) &
    !$omp private(l, ed, cell_area, ed_area, p1, p2, p3) &
    !$omp schedule(static)
    do k=1, mesh%nt

      !Calculate K energy a la A. Gassman

      ke%f(k)=0._r8
      cell_area=0._r8
      do l=1, 3
        ed=mesh%tr(k)%ed(l)
        !Gassmann
        if(useTiledAreas)then ! Dist(trcc,tredmidpoint)*tred_length/2
          !Dubos/Gass
          ed_area=0.5*arclen(mesh%tr(k)%c%p, mesh%ed(ed)%c%p)*mesh%ed(ed)%leng
          !simplified version (assumes arclen(mesh%tr(k)%c%p, mesh%ed(ed)%c%p)=mesh%edhx(ed)%leng/2
          !ed_area=0.25*mesh%edhx(ed)%leng*mesh%ed(ed)%leng
          cell_area=cell_area+ed_area
        else !Use geodesic areas !Area of tr formed by trcc and two endpoints of edge
          p1=mesh%tr(k)%c%p
          p2=mesh%v(mesh%ed(ed)%v(1))%p
          p3=mesh%v(mesh%ed(ed)%v(2))%p
          ed_area=sphtriarea(p1, p2, p3)
        endif
        ke%f(k)=ke%f(k)+ed_area*u%f(ed)**2
      end do
      if(useTiledAreas)then
        ke%f(k)=ke%f(k)/cell_area
      else
        ke%f(k)=ke%f(k)/mesh%tr(k)%areag
      endif
    end do
    !$omp end parallel do

    return

  end subroutine ke_tr

  subroutine grad_tr2ed(ftr, fed, mesh)
    !---------------------------------------------------------------
    !Gradient using data of triang centers to aprox grad on edges
    ! in: ftr - scalar field defined at triangles
    ! out: fed - grad scalar field defined at edges (must already be allocated)
    !---------------------------------------------------------------
    type(grid_structure), intent(in) :: mesh
    type(scalar_field), intent(in):: ftr
    type(scalar_field), intent(inout) :: fed

    integer(i4)::  l, k1, k2
    real(r8):: signcor

    !Quick check for dimensions
    if(fed%n/=mesh%ne)then
      print*, "ERROR in scalar_tr2ed: dimensions do not match", mesh%ne, fed%n
      stop
    end if

    !$omp parallel do default(none) &
    !$omp shared(mesh, fed, ftr) &
    !$omp private(signcor, k1, k2) &
    !$omp schedule(static)
    do l=1,mesh%ne
      !Calculate gradient of scalar on tr edges
      !  always relative to the normal of the triangle edge - independent of HCT and HCm
      k1=mesh%ed(l)%sh(1)
      k2=mesh%ed(l)%sh(2)
      signcor=dsign( 1._r8, dot_product(mesh%ed(l)%nr, &
        mesh%tr(k2)%c%p-mesh%tr(k1)%c%p ))
      fed%f(l)=signcor*(ftr%f(k2)-ftr%f(k1))/mesh%edhx(l)%leng !/erad
    end do
       !$omp end parallel do

    return

  end subroutine grad_tr2ed

  subroutine scalar_elem_product(f1, f2, fprod)
    !---------------------------------------------------------------
    ! Element wise product of scalar fields
    !    fprod=f1*f2
    !---------------------------------------------------------------

    type(scalar_field), intent(in):: f1, f2
    type(scalar_field), intent(inout) :: fprod

    integer(i4)::l

    !Quick check for dimensions
    if(f1%n/=f2%n)then
      print*, "ERROR in scalar_elem_product: dimensions do not match", f1%n, f2%n
      stop
    end if
    if(f1%n/=fprod%n)then
      print*, "ERROR in scalar_elem_product: dimensions do not match", f1%n, fprod%n
      stop
    end if

    !$omp parallel do default(none) &
    !$omp shared(f1, f2, fprod) &
    !$omp schedule(static)
    do l=1,fprod%n
      fprod%f(l)=f1%f(l)*f2%f(l)
    end do
    !$omp end parallel do

    return

  end subroutine scalar_elem_product

  subroutine scalar_elem_divide(f1, f2, fdiv)
    !---------------------------------------------------------------
    ! Element wise division of scalar fields
    !    fprod=f1/f2
    !---------------------------------------------------------------

    type(scalar_field), intent(in):: f1, f2
    type(scalar_field), intent(inout) :: fdiv

    integer(i4)::l

    !Quick check for dimensions
    if(f1%n/=f2%n)then
      print*, "ERROR in scalar_elem_product: dimensions do not match", f1%n, f2%n
      stop
    end if
    if(f1%n/=fdiv%n)then
      print*, "ERROR in scalar_elem_product: dimensions do not match", f1%n, fdiv%n
      stop
    end if

    !$omp parallel do default(none) &
    !$omp shared(f1, f2, fdiv) &
    !$omp schedule(static)
    do l=1,fdiv%n
      if(abs(f2%f(l))<eps/100000.)then
        print*, "ERROR at scalar_elem_devide: Division by zero"
        stop "  This test case is not adequate for calculating the mommentum eq."
      else
        fdiv%f(l)=f1%f(l)/f2%f(l)
      endif
    end do
    !$omp end parallel do

    return

  end subroutine scalar_elem_divide

  subroutine vector_elem_product(v, f, vprod)
    !---------------------------------------------------------------
    ! Element wise product of vector and scalar field
    !    vprod=v*f
    !---------------------------------------------------------------

    type(scalar_field), intent(in):: f
    type(vector_field_cart), intent(in) :: v
    type(vector_field_cart), intent(inout) :: vprod

    integer(i4)::l

    !Quick check for dimensions
    if(v%n/=f%n)then
      print*, "ERROR in vector_elem_product: dimensions do not match", v%n, f%n
      stop
    end if
    if(f%n/=vprod%n)then
      print*, "ERROR in vector_elem_product: dimensions do not match", f%n, vprod%n
      stop
    end if

    !$omp parallel do default(none) &
    !$omp shared(v, f, vprod) &
    !$omp schedule(static)
    do l=1,vprod%n
      vprod%p(l)%v=v%p(l)%v*f%f(l)
    end do
    !$omp end parallel do

    return

  end subroutine vector_elem_product


  subroutine zero_vector(v)
    !---------------------------------------------------------------
    ! Zero a vector (omp paralel)
    !    v must already be allocated
    !---------------------------------------------------------------

    real(r8)::v(:)

    !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
    !$OMP SHARED(v)
    v=0.
    !$OMP END PARALLEL WORKSHARE

    return

  end subroutine zero_vector

end module swm_operators
