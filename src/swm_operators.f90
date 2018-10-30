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

  subroutine scalar_tr2hx(ftr, fhx, mesh)
    !---------------------------------------------------------------
    !Interpolate from triang centers to cell centers (linear interpolation)
    ! in: ftr - scalar field defined at triangles
    ! out: fhx - scalar field defined at cells (must already be allocated)
    !---------------------------------------------------------------
    type(grid_structure), intent(in) :: mesh
    type(scalar_field), intent(in):: ftr
    type(scalar_field), intent(inout) :: fhx

    integer(i4)::  i, j

    !Quick check for dimensions
    if(fhx%n/=mesh%nv)then
      print*, "ERROR in scalar_tr2hx: dimensions do not match", mesh%nv, fhx%n
      stop
    end if

    !$omp parallel do &
    !$omp default(none) &
    !$omp shared(mesh,  wachc_tr2v, fhx, ftr) &
    !$omp shared(useSinterpolTrisk, useSinterpolBary) &
    !$omp private(j) &
    !$omp schedule(static)
    do i=1,mesh%nv
        !Calculate PV at cell centers
        fhx%f(i)=0._r8
        if(useSinterpolTrisk)then
          !Weighted average, using areas - 1st order
          do j=1, mesh%v(i)%nnb
            fhx%f(i)=fhx%f(i)+&
              mesh%hx(i)%hxtr_area(j)*ftr%f(mesh%v(i)%tr(j))
          end do
        elseif(useSinterpolBary)then
          !Wachspress coordinates - 2nd order
          do j=1, wachc_tr2v(i)%n
            fhx%f(i)=fhx%f(i)+&
              wachc_tr2v(i)%w(j)*ftr%f(wachc_tr2v(i)%v(j))
          end do
           !print*, q_hx%f(i), scinterpol_wachspress(mesh%v(i)%p, q_tr, mesh)
        end if
    end do
    !$omp end parallel do

    return

  end subroutine scalar_tr2hx

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

 subroutine vector_edhx2hx_perot(ued, vhx, mesh)
    !---------------------------------------------------------------
    !Reconstruct from cell edges to hexagonal centers (perot)
    ! in: ued - vector field defined at cell edges with normal components only
    ! out: vhx - vector field defined at cell centers (must already be allocated)
    !---------------------------------------------------------------

    type(grid_structure), intent(in) :: mesh
    type(scalar_field), intent(in):: ued
    type(vector_field_cart), intent(inout) :: vhx

    integer(i4):: i, l, ed
    real(r8):: signcor, ed_area, vectmp(1:3)

    !Quick check for dimensions
    if(vhx%n/=mesh%nv)then
      print*, "ERROR in vector_edhx2hx_perot: dimensions do not match", mesh%nv, vhx%n
      stop
    end if

    !$omp parallel do &
    !$omp default(none) &
    !$omp shared(mesh, ued, vhx) &
    !$omp shared(useStagHTC, useStagHC, useTiledAreas) &
    !$omp private(l, signcor, vectmp, ed) &
    !$omp schedule(static)
    do i=1,mesh%nv
      !Reconstruct velocities to cells - use Perot's method
      vectmp=0._r8
      if(useStagHTC)then
        do l=1,mesh%v(i)%nnb
          ed=mesh%v(i)%ed(l)
          signcor=mesh%hx(i)%ttgout(l)
          vectmp=vectmp + &
            signcor*(ued%f(ed))*(mesh%ed(ed)%c%p-mesh%v(i)%p)*mesh%edhx(ed)%leng
        end do
      elseif(useStagHC)then
        do l=1,mesh%v(i)%nnb
          ed=mesh%v(i)%ed(l)
          signcor=mesh%hx(i)%nr(l)
          vectmp=vectmp + &
            signcor*(ued%f(ed))*(mesh%edhx(ed)%c%p-mesh%v(i)%p)*mesh%edhx(ed)%leng
        end do
      end if
      if(useTiledAreas)then
        vhx%p(i)%v=vectmp/mesh%hx(i)%areat
      else
        vhx%p(i)%v=vectmp/mesh%hx(i)%areag
      end if
      vhx%p(i)%v=proj_vec_sphere(vhx%p(i)%v, mesh%v(i)%p)
    end do
    !$omp end parallel do

    return

  end subroutine vector_edhx2hx_perot

  subroutine grad_ed(f, grad, mesh)
    !---------------------------------------------------------------
    !Calculate gradient on edges based on scalar at cells
    !---------------------------------------------------------------
    type(grid_structure), intent(in) :: mesh
    type(scalar_field), intent(in):: f ! scalar at cells
    type(scalar_field), intent(inout):: grad !gradient at edges

    integer(i4):: l
    real(r8):: signcor

    !$omp parallel do &
    !$omp default(none) &
    !$omp shared(mesh, f, grad) &
    !$omp shared(useStagHTC, useStagHC) &
    !$omp private(signcor) &
    !$omp schedule(static)
    do l=1, mesh%ne
      !print*, l
      !Calculate gradient of ghbk on edges
      if(useStagHC)then
        signcor=dsign( 1._r8, dot_product(mesh%edhx(l)%nr, &
          mesh%v(mesh%edhx(l)%sh(2))%p-mesh%v(mesh%edhx(l)%sh(1))%p ))
        grad%f(l)=signcor*(f%f(mesh%edhx(l)%sh(2))-f%f(mesh%edhx(l)%sh(1)))
      elseif(useStagHTC)then
        ! Obs: the tangent of a triangle edge (which is used as normal
        !   of the voronoi cell edges) is always defined such
        !   that it point from its vertex 1 to its vertex 2
        !   Therefore, the gradient needs no correction term (n_ei)
        grad%f(l)=(f%f(mesh%ed(l)%v(2))-f%f(mesh%ed(l)%v(1)))
      end if
      grad%f(l)=grad%f(l)/mesh%ed(l)%leng/erad
    end do
    !$omp end parallel do

    return

  end subroutine grad_ed


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

  subroutine div_hx(u, div, mesh)
    !---------------------------------------------------------------
    !Calculate divergence at voronoi cells (hexagons)
    !   based on edge normal velocities
    !---------------------------------------------------------------
    type(grid_structure), intent(in) :: mesh
    type(scalar_field), intent(in):: u ! velocity at cell edges
    type(scalar_field), intent(inout):: div !divergence - must be already allocated

    integer(i4):: i, j, k, l, ed
    real(r8):: signcor

    !$omp parallel do &
    !$omp default(none) &
    !$omp shared(mesh, u, div) &
    !$omp shared(useStagHTC, useStagHC, useTiledAreas ) &
    !$omp private(j, l, signcor) &
    !$omp schedule(static)
    do i=1,mesh%nv
      !Divergence of uh on unit sphere using div_cell_Cgrig
      !For all edges forming the hexagon
      div%f(i)=0._r8
      do j=1, mesh%v(i)%nnb
        !Get edge index
        l=mesh%v(i)%ed(j)
        !Get edge outer normal related to the hexagon
        if(useStagHTC)then
          signcor=real(mesh%hx(i)%ttgout(j), r8)
        elseif(useStagHC)then
          signcor=real(mesh%hx(i)%nr(j), r8)
        end if
        !Calculate numerical integration
        !lengths are multiplied by erad, and area by erad**2, therefore /erad
        div%f(i)=div%f(i)+signcor*u%f(l)*mesh%edhx(l)%leng
         !divu%f(i)=divu%f(i)+signcor*u%f(l)*mesh%edhx(l)%leng
      end do
      if(useTiledAreas)then
        div%f(i)=div%f(i)/mesh%hx(i)%areat/erad
      else
        div%f(i)=div%f(i)/mesh%hx(i)%areag/erad
      end if
      !divu%f(i)=divu%f(i)/mesh%hx(i)%areag/erad
      !divuh%f(i)=div_cell_Cgrid(i,uh,mesh)/erad
    end do
    !$omp end parallel do

    return

  end subroutine div_hx

  subroutine diffusion_hx(u, lapu, mesh)
    !---------------------------------------------------------------
    !Calculate diffusion (Laplacian) of velocities
    !---------------------------------------------------------------
    type(grid_structure), intent(in) :: mesh
    type(scalar_field), intent(in):: u ! velocity at cell edges
    type(scalar_field), intent(inout):: lapu !divergence - must be already allocated


    integer(i4):: i, j, k, l, ed, k1, k2
    real(r8):: signcor, lapu_n, lapu_t
    real(r8):: etarel1, etarel2

    !if(difus==0)then !nothing to do
    !  return
    !end if

    !Needs verification!!!

      !Add difusion
      do l=1, mesh%ne

        !Calculate normal part
        lapu_n=0.
        !Calculate gradient of div for this edge
        if(useStagHC)then
          signcor=dsign( 1._r8, dot_product(mesh%edhx(l)%nr, &
            mesh%v(mesh%edhx(l)%sh(2))%p-mesh%v(mesh%edhx(l)%sh(1))%p ))
          lapu_n=signcor*(divu%f(mesh%edhx(l)%sh(2))-divu%f(mesh%edhx(l)%sh(1)))
        elseif(useStagHTC)then
          ! Obs: the tangent of a triangle edge (which is used as normal
          !   of the voronoi cell edges) is always defined such
          !   that it point from its vertex 1 to its vertex 2
          !   Therefore, the gradient needs no correction term (n_ei)
          lapu_n=(divu%f(mesh%ed(l)%v(2))-divu%f(mesh%ed(l)%v(1)))
        end if
        lapu_n=lapu_n/mesh%ed(l)%leng/erad

        !Calculate tangent part
        lapu_t=0.

        !Calculate gradient of relative vorticity for this edge
        if(useStagHC)then
          k1=mesh%edhx(l)%v(1)
          k2=mesh%edhx(l)%v(2)
          !Get relative vorticity
          etarel1=eta%f(k1)-2.*Omega*dsin(mesh%tr(k1)%c%lat)
          etarel2=eta%f(k2)-2.*Omega*dsin(mesh%tr(k2)%c%lat)
          signcor=dsign( 1._r8, dot_product(mesh%edhx(l)%tg, &
            mesh%tr(k2)%c%p-mesh%tr(k1)%c%p ))
          lapu_t=signcor*(etarel2-etarel1)
        elseif(useStagHTC)then
          k1=mesh%ed(l)%sh(1)
          k2=mesh%ed(l)%sh(2)
          !Get relative vorticity
          etarel1=eta%f(k1)-2.*Omega*dsin(mesh%tr(k1)%c%lat)
          etarel2=eta%f(k2)-2.*Omega*dsin(mesh%tr(k2)%c%lat)
          signcor=dsign( 1._r8, dot_product(-mesh%ed(l)%nr, &
            mesh%tr(k2)%c%p-mesh%tr(k1)%c%p ))
          lapu_t=signcor*(etarel2-etarel1)
        end if
        lapu_t=lapu_t/mesh%ed(l)%leng/erad

        !Global vector laplacian
        lapu%f(l)=lapu_n-lapu_t
        !print*, l, lapu%f(l), lapu_n, lapu_t
        !print*, l, momeq(l), difus*lapu%f(l), momeq(l)+ difus*lapu%f(l), abs(momeq(l))>abs(momeq(l)+ difus*lapu%f(l))
        !momeq(l)=momeq(l) + difus*lapu%f(l)
      end do

    return

  end subroutine diffusion_hx

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

  subroutine ke_hx(u, vhx, ketr, uh, h, ke, mesh)
    !---------------------------------------------------------------
    !Calculate kinetic energy at cells
    ! in: u at cell edges
    ! in (optional) : kinetic energy at triangles - used in Gasman's 2013 scheme
    ! out: ke at cells
    !---------------------------------------------------------------
    type(grid_structure), intent(in) :: mesh
    type(scalar_field), intent(in):: u ! velocity at cell edges
    type(vector_field_cart), intent(in):: vhx !reconstructed velocity - used in PXT16 scheme
    type(scalar_field), intent(in):: ketr !kinetic energy at triangles - used in GASS13 scheme
    type(scalar_field), intent(in):: uh !uh field - used in melvin scheme
    type(scalar_field), intent(in):: h  !h field - used in melvin scheme
    type(scalar_field), intent(inout):: ke !kinetic energy at cells

    integer(i4):: i, j, k, l, ed
    real(r8):: ed_area, cell_area, signcor, vectmp(1:3)

    !$omp parallel do &
    !$omp default(none) &
    !$omp shared(mesh, u, vhx, ketr, ke, uh, h) &
    !$omp shared(gasscoef) &
    !$omp shared(useSinterpolTrisk, useSinterpolBary, useTiledAreas) &
    !$omp shared(useStagHTC, useStagHC, useReconmtdTrisk,useReconmtdPerhx ) &
    !$omp shared(useReconmtdGass, useReconmtdMelv) &
    !$omp private(l, j, k, signcor, vectmp, ed, cell_area, ed_area) &
    !$omp schedule(static)
    do i=1,mesh%nv
      ke%f(i)=0._r8
      !Calculate K energy
      if(useReconmtdTrisk .or. useReconmtdGass)then
        cell_area=0._r8
        do l=1, mesh%v(i)%nnb
          ed=mesh%v(i)%ed(l)
          ed_area=0.25*mesh%ed(ed)%leng*mesh%edhx(ed)%leng
          cell_area=cell_area+ed_area
          ke%f(i)=ke%f(i)+ed_area*u%f(ed)**2
        end do
        !Kin_energy%f(i)=Kin_energy%f(i)/mesh%hx(i)%areag
        ke%f(i)=ke%f(i)/cell_area
      elseif(useReconmtdPerhx)then
        ke%f(i)=dot_product( vhx%p(i)%v,  vhx%p(i)%v)/2._r8
      elseif(useReconmtdMelv)then
        if(useStagHTC)then !use trsk method
          cell_area=0._r8
          do l=1, mesh%v(i)%nnb
            ed=mesh%v(i)%ed(l)
            ed_area=0.25*mesh%ed(ed)%leng*mesh%edhx(ed)%leng
            cell_area=cell_area+ed_area
            ke%f(i)=ke%f(i)+ed_area*uh%f(ed)*u%f(ed)
          end do
          !Kin_energy%f(i)=Kin_energy%f(i)/mesh%hx(i)%areag
          ke%f(i)=ke%f(i)/h%f(i)/cell_area
        elseif(useStagHC)then !use perots method
          vectmp=0._r8
          do l=1,mesh%v(i)%nnb
            ed=mesh%v(i)%ed(l)
            signcor=mesh%hx(i)%nr(l)
            vectmp=vectmp + &
              signcor*(uh%f(ed))*(mesh%edhx(ed)%c%p-mesh%v(i)%p)*mesh%edhx(ed)%leng
          end do
          vectmp=vectmp/mesh%hx(i)%areag
          vectmp=proj_vec_sphere(vectmp, mesh%v(i)%p)
          ke%f(i)=dot_product( vectmp,  vhx%p(i)%v)/h%f(i)/2._r8
        end if
      end if

      !Calculate K energy a la A. Gassman
      if(useReconmtdGass)then
        !Kin_energy_tr%f(k)=0._r8
        ke%f(i)=ke%f(i)*gasscoef
        do j=1, mesh%v(i)%nnb
          k=mesh%v(i)%tr(j)
          ke%f(i)=ke%f(i)+&
            (1.0_r8-gasscoef)*mesh%hx(i)%hxtr_area(j)*ketr%f(k)
           !print*, mesh%hx(i)%hxtr_area(j), Kin_energy_tr%f(k),Kin_energy%f(i)
        end do
         !print*, Kin_energy%f(i)
      end if

    end do
    !$omp end parallel do

    return

  end subroutine ke_hx

 subroutine coriolis_ed(u, uh, eta_ed, q_ed, uhq_perp, mesh)
    !---------------------------------------------------------------
    !Calculate Coriolis term
    !---------------------------------------------------------------
    type(grid_structure), intent(in) :: mesh
    type(scalar_field), intent(in):: u, uh, eta_ed, q_ed ! scalar fields
    type(scalar_field), intent(inout):: uhq_perp !gradient at edges

    !Indexes
    integer(i4):: i, i1, i2 !For node values
    integer(i4):: k, k1, k2 !For triangles
    integer(i4):: l, l1, l2 !For edges
    integer(i4):: j     !Auxiliar
    integer(i4):: ed    !Edge index
    integer(i4):: tr   !Triangle index
    integer(i4):: node !Node index
    integer(i4):: edcelli !Edge relative to cell i

    !Time
    real(r8):: t, dt

    !Temporary scalars
    real(r8):: ed_area
    real(r8):: signcor
    real(r8):: qtmp
    real(r8):: cell_area
    real(r8):: vectmp(1:3)
    real(r8):: vectmp1(1:3)
    real(r8):: vectmp2(1:3)
    real(r8):: d1
    real(r8):: d2

    real(r8):: lambda
    ! eta_ct, h_ct, h0, grad_ct
    real(r8):: ratiotrsk

    !OPENMP PARALLEL DO
    !$omp parallel do &
    !$omp default(none) &
    !$omp shared(mesh, uh, u, h, q_ed, momeq, uhq_perp, vhq_tr, vhq_hx) &
    !$omp shared(isTrskindLow, uhq_perp_exact, eta_ed) &
    !$omp shared(useSinterpolTrisk, useSinterpolBary) &
    !$omp shared(useStagHTC, useStagHC, useCoriolisMtdHyb, useCoriolisMtdDtred) &
    !$omp shared(useCoriolisMtdPered, useCoriolisMtdTrisk, useCoriolisMtdExact, noPV) &
    !$omp private(edcelli, signcor, qtmp) &
    !$omp private(i, i1, i2, j, l1, l2, node, k, k1, k2, d1, d2) &
    !$omp private(vectmp, vectmp1, vectmp2) &
    !$omp schedule(static)
    do l=1, mesh%ne
      uhq_perp%f(l)=0._r8
      !trskind%f(l)=trsk_order_index(l, u%pos, mesh)
      if(useCoriolisMtdHyb)then
        if(.not. isTrskindLow(l))then
          !mtdtmp="pered"
          useCoriolisMtdTrisk=.false.
          useCoriolisMtdPered=.true.
        else
          !mtdtmp="trsk"
          useCoriolisMtdTrisk=.true.
          useCoriolisMtdPered=.false.
        end if
      end if
      !For cells that share the edge
      if(useCoriolisMtdTrisk)then
        do i=1,2 !cells sharing edge l
          node=mesh%edhx(l)%sh(i)
          edcelli=getedindexonhx(l, node, mesh)
          !For all edges of cell, add contribution to reconstruction
          do j=1, mesh%v(node)%nnb
            !Edge's global index
            k=mesh%v(node)%ed(j)
            !Apply sign corrections
            if(useStagHTC)then
              signcor=dsign(1._r8,1._r8*mesh%hx(node)%ttgout(j)*mesh%hx(node)%tnrccw(edcelli))
            elseif(useStagHC)then
              signcor=dsign(1._r8,-1._r8*mesh%hx(node)%nr(j)*mesh%hx(node)%tg(edcelli))
            endif
            if(noPV)then
              qtmp=0.5_r8*(eta_ed%f(l)+eta_ed%f(k))
              uhq_perp%f(l)= uhq_perp%f(l)+ &
                u%f(k)*qtmp*mesh%edhx(k)%leng*mesh%hx(node)%trskw(edcelli, j)*signcor
            else
              qtmp=0.5_r8*(q_ed%f(l)+q_ed%f(k))
              uhq_perp%f(l)= uhq_perp%f(l)+ &
                uh%f(k)*qtmp*mesh%edhx(k)%leng*mesh%hx(node)%trskw(edcelli, j)*signcor
            endif
          end do
        end do
        uhq_perp%f(l)=uhq_perp%f(l)/mesh%ed(l)%leng
         !print*, l
      elseif(useCoriolisMtdDtred)then !This method only works for layer model (with PV)
        !Use dual triangle reconstruction
        k1=mesh%ed(l)%sh(1)
        k2=mesh%ed(l)%sh(2)

        if(useStagHTC)then
          d1=arclen(mesh%ed(l)%c%p,mesh%tr(k1)%c%p)
          d2=arclen(mesh%ed(l)%c%p,mesh%tr(k2)%c%p)
          vectmp=(d2*vhq_tr%p(k1)%v+d1*vhq_tr%p(k2)%v)/(d1+d2)
          vectmp=proj_vec_sphere(vectmp, mesh%ed(l)%c%p)
          uhq_perp%f(l)=-dot_product(vectmp,mesh%ed(l)%nr)
        elseif(useStagHC)then
          !d1=arclen(mesh%edhx(l)%c%p,mesh%tr(k1)%c%p)
          !d2=arclen(mesh%edhx(l)%c%p,mesh%tr(k2)%c%p)
          ! d1/(d1+d2), d2/( d1+d2) should be 1/2
          !vectmp=(d2*vhq_tr%p(k1)%v+d1*vhq_tr%p(k2)%v)/(d1+d2)
          vectmp=0.5*vhq_tr%p(k1)%v+0.5*vhq_tr%p(k2)%v
          !vectmp=proj_vec_sphere(vectmp, mesh%edhx(l)%c%p)
          uhq_perp%f(l)=dot_product(vectmp,mesh%edhx(l)%tg)
        end if
      elseif(useCoriolisMtdPered)then !This method only works for layer model (with PV)
        !print*, l
        do i=1,2 !cells sharing edge l
          node=mesh%edhx(l)%sh(i)
          !Reconstruct velocity with pv to cell centers
          vectmp=0._r8
          if(useStagHTC)then
            do l1=1,mesh%v(node)%nnb
              l2=mesh%v(i)%ed(l1)
              qtmp=0.5_r8*(q_ed%f(l)+q_ed%f(l2))
              signcor=mesh%hx(node)%ttgout(l1)
              vectmp=vectmp + &
                signcor*(uh%f(l2))*qtmp*(mesh%ed(l2)%c%p-mesh%v(node)%p)*mesh%edhx(l2)%leng
            end do
          elseif(useStagHC)then
            do l1=1,mesh%v(node)%nnb
              !print*, l, i
              l2=mesh%v(node)%ed(l1)
              qtmp=0.5_r8*(q_ed%f(l)+q_ed%f(l2))
              signcor=mesh%hx(node)%nr(l1)
              vectmp=vectmp + &
                signcor*(uh%f(l2))*qtmp*(mesh%edhx(l2)%c%p-mesh%v(node)%p)*mesh%edhx(l2)%leng
            end do
          end if
          vectmp=vectmp/mesh%hx(node)%areag
          if(i==1)then
            vectmp1=proj_vec_sphere(vectmp, mesh%v(node)%p) !*h%f(node) !vectmp*h%f(node) !
          else
            vectmp2=proj_vec_sphere(vectmp, mesh%v(node)%p) !*h%f(node) !vectmp*h%f(node) !
          end if
        end do
        !print*, l, vectmp
        !print*
        !Interpolate to edges
        i1=mesh%edhx(l)%sh(1)
        i2=mesh%edhx(l)%sh(2)
        !Simples average works on both HC and HTC cases
        vectmp=0.5*vectmp1+0.5*vectmp2
        if(useStagHTC)then
          !d1=arclen(mesh%ed(l)%c%p,mesh%v(i1)%p)
          !d2=arclen(mesh%ed(l)%c%p,mesh%v(i2)%p)
          !vectmp=(d2*vectmp1+d1*vectmp2)/(d1+d2)
          vectmp=proj_vec_sphere(vectmp, mesh%ed(l)%c%p)
          uhq_perp%f(l)=-dot_product(vectmp,mesh%ed(l)%nr)
        elseif(useStagHC)then
          !d1=arclen(mesh%edhx(l)%c%p,mesh%v(k1)%p)
          !d2=arclen(mesh%edhx(l)%c%p,mesh%v(k2)%p)
          ! d1/(d1+d2), d2/( d1+d2) should be 1/2
          !print*, l, d1, d2, d1/(d1+d2), d2/(d1+d2)
          vectmp=proj_vec_sphere(vectmp, mesh%edhx(l)%c%p)
          uhq_perp%f(l)=dot_product(vectmp,mesh%edhx(l)%tg)
           !write(55,*) l, i1, i2, uhq_perp%f(l), vectmp
           !write(55,*)
        end if
         !print*, l, "pered"
      elseif(useCoriolisMtdExact)then
        uhq_perp%f(l)=uhq_perp_exact%f(l)
      endif

     ! momeq(l)=momeq(l) -uhq_perp%f(l)  !-uhq_perp_exact%f(l) !
    end do
    !$omp end parallel do

    return

  end subroutine coriolis_ed


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

  subroutine apvm(q_grad_tr, v_hx, pvspar, dt, q_ed, mesh)
    !---------------------------------------------------------------
    !Calculates Anticipated Vorticity Method into PV at edges
    !  see Ringler 2010 and Peixoto 2016
    ! in: PV gradient at triangles q_grad_tr
    ! in: reconstructed velocity at cells v_hx
    ! in: pvspar, dt, parameters used in method
    ! out: q_ed - updates pv on edges
    !---------------------------------------------------------------
    type(grid_structure), intent(in) :: mesh
    type(vector_field_cart), intent(in):: q_grad_tr
    type(vector_field_cart), intent(in):: v_hx
    type(scalar_field), intent(inout) :: q_ed
    real(r8), intent(in):: pvspar, dt

    integer(i4)::  l


    if(.not.useAPVM)then
      !don't add apvm
      return
    end if

    if(useCLUST)then
        print*, "ERROR in apvm: Clust needs debugging"
        stop
    end if

    !$omp parallel do &
    !$omp default(none) &
    !$omp shared(mesh, q_grad_tr, v_hx, pvspar, dt, q_ed) &
    !$omp schedule(static)
    do l=1, mesh%ne
        q_ed%f(l)=q_ed%f(l)-pvspar*dt*dot_product( &
          q_grad_tr%p(mesh%ed(l)%sh(1))%v +q_grad_tr%p(mesh%ed(l)%sh(2))%v, &
          v_hx%p(mesh%ed(l)%v(1))%v + v_hx%p(mesh%ed(l)%v(2))%v )
    end do
    !$omp end parallel do

    !Clust stuff
    !!$omp shared(useSinterpolTrisk, useSinterpolBary, useCoriolisMtdDtred) &
    !!$omp shared(useStagHTC, useStagHC, useAPVM, useCLUST, dt, bclust, pvspar) &
    !!$omp shared(test_lterror) &
    !!$omp private(signcor, d1, d2, k, lambda, vectmp, nvectmp, upwindq) &

        !Calculate parameter to check whether the flow is parallel to the edge
        !Total velocity at edges (based on Perots method)
        !vectmp=0.5_r8*(v_hx%p(mesh%ed(l)%v(1))%v+v_hx%p(mesh%ed(l)%v(2))%v)
        !nvectmp=norm(vectmp)

        !Use CLUST method
        !if(nvectmp>eps/1000000.)then
         ! lambda=abs(u%f(l))/nvectmp
         ! bclust=pvspar
         ! if(dot_product(vectmp, mesh%tr(mesh%ed(l)%sh(1))%c%p-mesh%ed(l)%c%p)<0)then
         !   !The upwind triangle point is
         !   k= mesh%ed(l)%sh(1)
         ! else
         !   k= mesh%ed(l)%sh(2)
         ! end if
         ! if(useStagHC)then
         !   upwindq=q_tr%f(k)- &
         !     dot_product(mesh%tr(mesh%ed(l)%sh(1))%c%p-mesh%edhx(l)%c%p, &
         !     q_grad_tr%p(k)%v)
         ! else
         !   upwindq=q_tr%f(k)- &
         !     dot_product(mesh%tr(mesh%ed(l)%sh(1))%c%p-mesh%ed(l)%c%p, &
         !     q_grad_tr%p(k)%v)
         ! end if
         ! !print*, l, lambda, (1-bclust*lambda), q_ed%f(l)
         ! q_ed%f(l)=(1-bclust*lambda)*q_ed%f(l)+lambda*bclust*upwindq
           !print*, l, lambda*bclust, upwindq, q_ed%f(l)
           !print*
        !end if


    return

  end subroutine apvm


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
