

! ========================================================

SUBROUTINE coriolis(u,v,phi,fu,fv)

! To evaluate the Coriolis terms on the C-grid,
! taking account of energy conservation and improved
! Rossby mode dispersion

USE grid
USE constants

IMPLICIT NONE

! INTEGER, INTENT(IN) :: nx, ny
REAL*8, INTENT(IN) :: u(nx,ny), v(nx,ny), phi(nx,ny)
!REAL*8, INTENT(OUT) :: fu(nx,ny), fv(nx,ny)
REAL*8, INTENT(OUT) :: fu(nx,ny+1), fv(nx,ny)
REAL*8 :: tempv(nx,ny+1), tempu(nx,ny), tempp(nx,ny)

INTEGER :: i, im, ip, j, jm, jp

! ------

! phi v cos(lat) at v points
DO j = 2, ny
  jm = j - 1
  DO i = 1, nx
    tempv(i,j) = 0.5*(phi(i,jm) + phi(i,j))*cosv(j)*v(i,j)
  ENDDO
ENDDO
! zero at polar latitudes
tempv(:,1) = 0.0
tempv(:,ny+1) = 0.0

! Average to phi points and times f / phi
DO j = 1, ny
  jp = j + 1
  DO i = 1, nx
    tempp(i,j) = 0.5*(tempv(i,j) + tempv(i,jp))*twoomega*singeolatp(i,j)/phi(i,j)
  ENDDO
ENDDO

! Average to u points and divide by cos(lat) to get
! fv at u points
DO j = 1, ny
  DO i = 1, nx
    im = i-1
    IF (im == 0) im = nx
    fv(i,j) = 0.5*(tempp(im,j) + tempp(i,j))/cosp(j)
  ENDDO
ENDDO

! ------

! phi u at u points
DO j = 1, ny
  DO i = 1, nx
    im = i - 1
    IF (im == 0) im = nx
    tempu(i,j) = 0.5*(phi(im,j) + phi(i,j))*u(i,j)
  ENDDO
ENDDO

! Average to phi points and times f / phi
DO j = 1, ny
  DO i = 1, nx
    ip = i + 1
    IF (i == nx) ip = 1
    tempp(i,j) = 0.5*(tempu(i,j) + tempu(ip,j))*twoomega*singeolatp(i,j)/phi(i,j)
  ENDDO
ENDDO

! Average to v points to get
! fu at v points
DO j = 2, ny
  jm = j - 1
  DO i = 1, nx
    fu(i,j) = 0.5*(tempp(i,jm) + tempp(i,j))
  ENDDO
ENDDO
! zero at polar latitudes
fu(:,1) = 0.0
fu(:,ny+1) = 0.0

! ------

END SUBROUTINE coriolis

! ========================================================

! ========================================================

SUBROUTINE newcoriolis(u,v,phi,fu,fv)

! To evaluate the Coriolis terms on the C-grid,
! taking account of energy conservation and improved
! stability (but probably impairing Rossby mode dispersion)

USE grid
USE constants

IMPLICIT NONE

! INTEGER, INTENT(IN) :: nx, ny
REAL*8, INTENT(IN) :: u(nx,ny), v(nx,ny), phi(nx,ny)
!REAL*8, INTENT(OUT) :: fu(nx,ny), fv(nx,ny)
REAL*8, INTENT(OUT) :: fu(nx,ny+1), fv(nx,ny)
REAL*8 :: tempv(nx,ny+1), tempu(nx,ny), tempz(nx,ny), tempz2(nx,ny+1)

INTEGER :: i, im, ip, j, jm, jp

! ------

! phi at vorticity points
DO j = 2, ny
  jm = j - 1
  DO i = 1, nx
    im = i - 1
    IF (im == 0) im = nx
    tempz2(i,j) = 0.25*(phi(im,jm) + phi(im,j) + phi(i,jm) + phi(i,j))
  ENDDO
ENDDO

! ------

! phi v cos(lat) at v points
DO j = 2, ny
  jm = j - 1
  DO i = 1, nx
    tempv(i,j) = 0.5*(phi(i,jm) + phi(i,j))*cosv(j)*v(i,j)
  ENDDO
ENDDO
! zero at polar latitudes
tempv(:,1) = 0.0
tempv(:,ny+1) = 0.0

! Average to vorticity points and times f / phi
DO j = 2, ny
  DO i = 1, nx
    im = i-1
    IF (im == 0) im = nx
    tempz(i,j) = 0.5*(tempv(im,j) + tempv(i,j))*twoomega*singeolatz(i,j)/tempz2(i,j)
  ENDDO
ENDDO
tempz(:,1) = 0.0d0
tempz(:,nyp) = 0.0d0

! Average to u points and divide by cos(lat) to get
! fv at u points
DO j = 1, ny
  jp = j + 1
  DO i = 1, nx
    fv(i,j) = 0.5*(tempz(i,j) + tempz(i,jp))/cosp(j)
  ENDDO
ENDDO

! ------

! phi u at u points
DO j = 1, ny
  DO i = 1, nx
    im = i - 1
    IF (im == 0) im = nx
    tempu(i,j) = 0.5*(phi(im,j) + phi(i,j))*u(i,j)
  ENDDO
ENDDO

! Average to vorticity points and times f / phi
DO j = 2, ny
  jm = j - 1
  DO i = 1, nx
    tempz(i,j) = 0.5*(tempu(i,jm) + tempu(i,j))*twoomega*singeolatz(i,j)/tempz2(i,j)
  ENDDO
ENDDO

! Average to v points to get
! fu at v points
DO j = 2, ny
  DO i = 1, nx
    ip = i + 1
    IF (i == nx) ip = 1
    fu(i,j) = 0.5*(tempz(i,j) + tempz(ip,j))
  ENDDO
ENDDO
! zero at polar latitudes
fu(:,1) = 0.0
fu(:,ny+1) = 0.0

! ------

END SUBROUTINE newcoriolis

! ========================================================

