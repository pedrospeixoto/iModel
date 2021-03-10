!======================================================================
! Module for the Matsuno baroclinic wave (Paldor et al 2019)
! tc56 - Eastward inertial gravity wave (EIGW)
! tc57 - Rossby wave (RW)
!======================================================================
module swm_matsuno
    use swm_data
    use constants, only: &
    i4, &
    r8, &
    erad, &
    Omega, &
    rad2deg, &
    deg2rad, &
    pi, &
    piby2, &
    grav, &
    datadir, &
    sec2day

    !Use main grid data structures
    use datastruct, only: &
    grid_structure, &
    scalar_field

    !Use routines from the spherical mesh pack
    use smeshpack, only: &
    convert_vec_sph2cart, &
    convert_vec_cart2sph, &
    insidetrmesh, &
    sph2cart

    !Use routines from the interpolation pack
    use interpack, only: &
    plot_hovmoller_diagram


contains

  subroutine parameters_matsuno(w,ks,n,Am,H0,coef)
  !----------------------------------------------------------------------------------
  ! 
  !----------------------------------------------------------------------------------
    real(r8), intent(inout) :: w, ks, Am, H0
    complex(r8),intent(inout) :: coef(0:4)
    real(r8) :: eps ! Lamb parameter
    real(r8) :: w_nk(1:3), a, g, k
    complex(r8) :: D(0:4), imag
    complex(r8) :: c1, c2, c3, c4, c5, c6
    integer, intent(inout) :: n
    integer :: j

    !Parameters
    a = erad
    g = grav
    H0 = 30.d0
    eps = (2.d0*omega*a)**2/(g*H0)
    Am = 0.00001
    n = 1
    k = 5.d0
    ks = 5.d0/a

    !Computes the values of delta
    imag = (0.d0,1.d0)
    D(0) = 3.d0*(g*H0*ks*ks + 2.d0*omega*sqrt(g*H0)*(2.d0*n+1.d0)/a)
    D(4) = -54.d0*omega*g*H0*ks/a
  
    do j = 1,3
      D(j) = sqrt(D(4)**2-4.d0*D(0)**3)
      D(j) = (D(j) + D(4))*0.5d0
      D(j) = D(j)**(1.d0/3.d0)
      D(j) = D(j)*exp((2.d0*pi*j/3.d0)*imag)
    end do
  
    ! Computes the omegas
    do j = 1, 3
      w_nk(j) = -real((D(j) + D(0)/D(j)))/3.d0
    end do
  
    if(testcase==56)then
      w = maxval(w_nk) !EIG
    else if(testcase==57)then
      w = -minval(abs(w_nk)) ! Rossby
    end if


    !Others parameters
    c1 = eps**(0.25d0)/a
    coef(0) = c1

    c2 = g*H0*eps**(0.25d0)/((w*w-g*H0*ks*ks)*a)/imag
    c3 = -sqrt(real((n+1.d0)*0.5d0))*(w/sqrt(g*H0)+ks)
    coef(1) = c2*c3
    c4 = -sqrt(real(n*0.5d0))*(w/sqrt(g*H0)-ks)
    coef(2) = c2*c4
  
    c5 = -sqrt(real((n+1.d0)*0.5d0))*(w+sqrt(g*H0)*ks)
    coef(3) = c2*c5
    c6 = sqrt(real(n*0.5d0))*(w-sqrt(g*H0)*ks)
    coef(4) = c2*c6
    return
  end subroutine parameters_matsuno


  function h_matsuno(lon,lat,w,coef,n,ks,Am,H0,t)
  !----------------------------------------------------------------------------------
  !
  !----------------------------------------------------------------------------------
    real(r8) :: h_matsuno
    real(r8) :: phi
    real(r8), intent(in) :: lon, lat, w, ks, Am,H0,t
    complex(r8),intent(in) :: coef(0:4)
    real(r8) :: x,y, yy
    integer(i4), intent(in) :: n
    complex(r8):: ubar, vbar(-1:1), phibar
    complex(r8) :: imag

    x = lon*erad
    y = lat*erad
    imag = (0.d0,1.d0)

    ! Computes vbar
    yy = coef(0)*y
    vbar(-1) = Am*Hermite(yy,n-1)*exp(-0.5d0*yy**2)
    vbar(0) = Am*Hermite(yy,n)*exp(-0.5d0*yy**2)
    vbar(1) = Am*Hermite(yy,n+1)*exp(-0.5d0*yy**2)
    ! Computes u_n and h_n
    !ubar = coef(1)*vbar(1) + coef(2)*vbar(-1)
    phibar = coef(3)*vbar(1) + coef(4)*vbar(-1)

    ! Computes u, v and h
    !u = real(ubar*exp((ks*x-w*t)*imag))
    !v = real(vbar(0)*exp((ks*x-w*t)*imag))
    phi = real(phibar*exp((ks*x-w*t)*imag))
    h_matsuno = phi/grav + H0
    !print*,hbar
    return
  end function h_matsuno


  ! Set velocity fields for Matsuno test case
  subroutine velocity_matsuno(u,v,lon,lat,w,coef,n,ks,Am,t)
    real(r8), intent(out) :: u,v
    real(r8), intent(in) :: lon, lat, w, ks, Am,t
    complex(r8),intent(in) :: coef(0:4)
    real(r8) :: x,y, yy
    integer(i4), intent(in) :: n
    complex(r8):: ubar, vbar(-1:1), hbar
    complex(r8) :: imag

    x = lon*erad
    y = lat*erad
    imag = (0.d0,1.d0)

    ! Computes vbar
    yy = coef(0)*y
    vbar(-1) = Am*Hermite(yy,n-1)*exp(-0.5d0*yy**2)
    vbar(0) = Am*Hermite(yy,n)*exp(-0.5d0*yy**2)
    vbar(1) = Am*Hermite(yy,n+1)*exp(-0.5d0*yy**2)
    ! Computes u_n and h_n
    ubar = coef(1)*vbar(1) + coef(2)*vbar(-1)
    !hbar = coef(3)*vbar(1) + coef(4)*vbar(-1)

    ! Computes u, v and h
    u = real(ubar*exp((ks*x-w*t)*imag))
    v = real(vbar(0)*exp((ks*x-w*t)*imag))
  end subroutine velocity_matsuno

  function Hermite(x,n)
  !----------------------------------------------------------------------------------
  ! compute the normalized Hermite polynomial
  !----------------------------------------------------------------------------------
    real(r8):: Hermite
    real(r8), intent(in) :: x
    real(r8) :: H(-1:n)
    integer, intent(in) :: n
    integer :: i

    H(-1) = 0.d0
    H(0) = 1.d0/pi**(0.25d0)
    do i = 0, n-1
      H(i+1) = x*sqrt(2.d0/(i+1))*H(i) - sqrt(real(i)/real(i+1))*H(i-1)
    end do
    Hermite = H(n)
  end function Hermite


  subroutine exact_matsuno(t)
  !----------------------------------------------------------------------------------
  ! compute the exact solution
  !----------------------------------------------------------------------------------
    real(r8), intent(in) ::  t
    complex(r8) :: coef(0:4)
    real(r8) :: w, ks, Am, H0, lon, lat, u_aux, v_aux
    real(r8):: vectmp(1:3)
    integer :: i,j, n
    
    !Computes the parameters
    call parameters_matsuno(w,ks,n,Am,H0,coef)
    
    !Computes the height field
    do i=1, mesh%nv
      lon=mesh%v(i)%lon
      lat=mesh%v(i)%lat
      h_exact%f(i) = h_matsuno(lon,lat,w,coef,n,ks,Am,H0,t)
    end do

    !Set velocity field
    do l=1,mesh%ne
      if(useStagHC)then
        lat = mesh%edhx(l)%c%lat
        lon = mesh%edhx(l)%c%lon
        call velocity_matsuno(u_aux,v_aux,lon,lat,w,coef,n,ks,Am,t)
        call convert_vec_sph2cart(u_aux, v_aux, mesh%edhx(l)%c%p, vectmp)
        u_exact%f(l)=dot_product(vectmp,mesh%edhx(l)%nr)
      elseif(useStagHTC)then
        lat = mesh%ed(l)%c%lat
        lon = mesh%ed(l)%c%lon
        call velocity_matsuno(u_aux,v_aux,lon,lat,w,coef,n,ks,Am,t)
        call convert_vec_sph2cart(u_aux, v_aux, mesh%ed(l)%c%p, vectmp)
        u_exact%f(l)=dot_product(vectmp,mesh%ed(l)%tg)
      end if
    end do
  end subroutine exact_matsuno


  subroutine matsuno_analysis(k,time)
  !----------------------------------------------------------------------------------
  ! Plot the Hovmoller diagrams
  !----------------------------------------------------------------------------------
    integer, intent(in) :: k

    !Filename
    character (len=70) :: filename

    real(r8), intent(in) :: time
    
    real(r8) :: t

    !Matsuno parameters
    real(r8) :: w, ks, Am, H0
    complex(r8) :: coef(0:4)

    !Open files
    integer :: iunit
    logical::  ifile
    
    !Counter 
    integer :: i

    call parameters_matsuno(w,ks,n,Am,H0,coef)
    call exact_matsuno(time)

    !Wave period
    if(testcase==56)then !EIGW
      t = 1.9d0
    else if(testcase==57)then !RW
      t = 18.5d0
    endif

    if(time*sec2day>(99.d0)*t .and. time*sec2day<(100.d0)*t) then 
    !if(time*sec2day>(0.d0)*t .and. time*sec2day<(1.d0)*t) then     
      !Plot the Hovmoller diagrams
      h%f = h%f - 30._r8
      h%name=trim(swmname)//"_h"!//trim(adjustl(trim(atime)))
      call plot_hovmoller_diagram(h, mesh,0.d0,-65.d0,k)
      
      !Exact solution diagrams
      !Plot the Hovmoller diagrams
      h_exact%f = h_exact%f - 30._r8
      h_exact%name=trim(swmname)//"_h_exact"!//trim(adjustl(trim(atime)))
      call plot_hovmoller_diagram(h_exact, mesh,0.d0,-65.d0,k)
      
      h%f = h%f + 30._r8
      h_exact%f = h_exact%f + 30._r8  
    end if
    
  
  end subroutine


end module
