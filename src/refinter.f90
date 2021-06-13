!================================================================================
! Refinter - module for topography based locally refined SCVT grids
! Generates smooth density functions for Lloyd's method
!================================================================================
module refinter
    !parameters
    use constants, only: &
        pi, &
        datadir, &
        pio2, &
        pi2, &
        r8, &
        rad2deg, &
        deg2rad, &
        n_lat, &
        n_lon, &
        latmin, &
        latmax, &
        lonmin, &
        lonmax, &
        altdir, &
        nlat_alt, &
        nlon_alt
    implicit none

contains
    !--------------------------------------------------------------------------
    ! Compute the density function defined in 'Voronoi Tessellations and their application to
    !  climate and global modeling' by Lili Ju, Todd Ringler and Max Gunzburger
    !--------------------------------------------------------------------------
    real (r8) function densf(x) result(dens_f)
        real (r8), dimension(1:2), intent(in) :: x
        real (r8), dimension(1:3) :: p
        real (r8), dimension(1:3) :: pc
        real (r8) :: coslat
        real (r8) :: lat
        real (r8) :: lon
        real (r8) :: latc
        real (r8) :: lonc
        real (r8) :: gammas
        real (r8) :: epsilons
        real (r8) :: dists
        real (r8) :: maxdist
        real (r8) :: sx

        !Center of refined region

        latc=-20.0*pi/180._r8
        lonc=-60.0*pi/180._r8
        coslat = dcos (latc)
        pc(1) = coslat * dcos (lonc)
        pc(2) = coslat * dsin (lonc)
        pc(3) = dsin (latc)

        lat = x(1)
        lon = x(2)
        coslat = dcos (lat)
        p(1) = coslat * dcos (lon)
        p(2) = coslat * dsin (lon)
        p(3) = dsin (lat)

        !Distance to center
        dists = dsqrt((p(1)-pc(1))**2+(p(2)-pc(2))**2+(p(3)-pc(3))**2)
        !dists = dsqrt((lat-latc)**2+(lon-lonc)**2)
        !Density function parameters
        gammas = 3._r8
        epsilons = pi/12._r8
        maxdist = pi/6._r8
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
        dens_f = dens_f/gammas**4
    end function densf



    !--------------------------------------------------------------------------
    ! Bilinear interpolation
    !--------------------------------------------------------------------------
    real (r8) function interpol_densf(p, densftable, xmin, xmax, ymin, ymax, nx, ny) result(ff)
        integer, intent(in) :: nx, ny
        real (r8), intent(in) :: xmin, xmax, ymin, ymax
        real (r8), dimension(1:2), intent(in) :: p
        real (r8), dimension(nx*ny,3), intent(in) :: densftable
        integer :: i, j
        integer :: i1, i2, j1, j2
        real (r8) :: f11, f12, f21, f22
        real (r8) :: x1, x2, y1, y2, x0, y0
        real (r8) :: dx, dy

        dx = (xmax-xmin)/(nx-1)
        dy = (ymax-ymin)/(ny-1)

        x0 = p(1)
        y0 = p(2)
        i1 = INT((x0 - xmin)/dx+1)
        i2 = CEILING((x0 - xmin)/dx+1)
        j1 = INT((y0 - ymin)/dy+1)
        j2 = CEILING((y0 - ymin)/dy+1)

        x1 = densftable((i1-1)*ny + 1, 1)
        y1 = densftable(j1,2)
        x2 = densftable((i2-1)*ny + 1,1)
        y2 = densftable(j2,2)

        f21 = densftable((i2-1)*ny + j1, 3)
        f12 = densftable((i1-1)*ny + j2, 3)
        f11 = densftable((i1-1)*ny + j1, 3)
        f22 = densftable((i2-1)*ny + j2, 3)


        if(x1 /= x2 .and. y1 /= y2) then
            ff = f11*(x2-x0)*(y2-y0)
            ff = ff + f21*(x0-x1)*(y2-y0)
            ff = ff + f12*(x2-x0)*(y0-y1)
            ff = ff + f22*(x0-x1)*(y0-y1)
            ff = ff/((x2-x1)*(y2-y1))
        else if (x1 == x2 .and. y1 /= y2 ) then
            ff = f11 + (f12-f11)*(y0-y1)/(y2-y1)        !linear interpolation
        else if (x1 /= x2 .and. y1 == y2 ) then
            ff = f11 + (f21-f11)*(x0-x1)/(x2-x1)    !linear interpolation
        else  
            ff = f11
        end if
	
    end function interpol_densf

    !--------------------------------------------------------------------------
    ! Load earth elevation data from ETOPO
    ! Etopo data must be in the directory /altitude/
    !--------------------------------------------------------------------------
    subroutine earth_elevation(iunit, altitude_table)
        integer :: iunit
        real (r8), dimension(nlat_alt*nlon_alt,3), intent(inout) :: altitude_table
        integer :: i, j
        character (len=100) :: filename
        character (len=100) :: buffer

        filename = trim(altdir)//"elevation_table_720x1440.dat"

        print*,'Loading ETOPO data...'
        open(iunit, file=filename)
          read(iunit,*) buffer
          read(iunit,*) buffer
          read(iunit,*) buffer
          read(iunit,*) ((altitude_table  (i,j), j=1,3), i=1,nlat_alt*nlon_alt)
        close(iunit)

    end subroutine earth_elevation



    !--------------------------------------------------------------------------
    ! Defines smooth Andes density function in lat-lon grid and saves in 
    ! altitude_table
    !--------------------------------------------------------------------------
    subroutine andes_density_table(altitude_table, iunit)
      real (r8), dimension(nlat_alt*nlon_alt,3), intent(inout) :: altitude_table
      integer, intent(in) :: iunit
      integer :: i, j, k, n, m, l
      integer :: old, new, aux
      real (r8) :: lat, lon, latc, lonc
      real (r8) :: dists, maxdist, gammas, epsilons
      real (r8), dimension(1:3) :: p
      real (r8), dimension(1:3) :: pc
      real (r8) :: coslat, a, b, c, sx, cutoff, alpha
      real (r8), dimension(:,:,:), allocatable  :: auxiliary_function
      real (r8), dimension(:,:,:), allocatable  :: topo
      real (r8) :: temp, dlat, dlon, lambda
      character (len=100) :: filename

      ! Load ETOPO data
      call earth_elevation(iunit, altitude_table)
      
      ! Allocation
      allocate (topo(nlat_alt, nlon_alt,0:1))
      allocate (auxiliary_function(nlat_alt, nlon_alt,0:1))
      
      ! Parameters
      ! Center of circular refinement region in lat-lon
      latc=-20._r8*deg2rad
      lonc=-70._r8*deg2rad
      
      ! Center of circular refinement region in R^3 coordinates
      coslat = dcos(latc)
      pc(1) = coslat*dcos(lonc)
      pc(2) = coslat*dsin(lonc)
      pc(3) = dsin(latc)

      !radius of circular refinement region
      !maxdist = 45._r8*deg2rad     
      maxdist = 28._r8*deg2rad
      
      !width of transition zone
      epsilons = 15._r8*deg2rad
      
      gammas = 2._r8
    
      !Cut-off altitude
      cutoff = 500._r8
    
      !ETOPO lat-lon grid space
      dlat = pi/(nlat_alt-1)
      dlon = pi2/(nlon_alt-1)

      ! Saves the topography data
      k = 1
      do i = 1, nlat_alt
        do j = 1, nlon_alt
          topo(i,j,0) = altitude_table(k,3)
          k = k + 1
        end do
      end do
      
      print*,'Generating smooth density function based on Andes topography...'

      ! Defines an auxiliary function the defines a circular refinement centered in pc
      auxiliary_function = 0._r8
      do i = 1, nlat_alt
        lat = -pio2 + (i-1)*dlat
        do j = 1, nlon_alt
          lon = -pi + (j-1)*dlon
          !Point in R^3 coordinates
          coslat = dcos(lat)
          p(1) = coslat*dcos(lon)
          p(2) = coslat*dsin(lon)
          p(3) = dsin(lat)
          
          !Evaluates the R^3 distance
          dists = dsqrt((p(1)-pc(1))**2+(p(2)-pc(2))**2+(p(3)-pc(3))**2)

          !Point in circular region
          if(dists<=maxdist)then
            auxiliary_function(i,j,0) = 1._r8
           !Point in transition zone
           else if(dists<=maxdist+epsilons)then
             sx = (maxdist+epsilons-dists)/epsilons
             auxiliary_function(i,j,0) = sx
          endif
        enddo
      enddo
      
      ! Applies moving average smoothing method in the auxiliary_function data
      ! Uses a box centered in (i,j) with l^2 elements
      l = 55
      do i = l, nlat_alt-l
        do j = l, nlon_alt-l
          temp =  0._r8
          do m = 0, l-1
            do n = 0, l-1
              temp = temp + auxiliary_function(i-(l-1)/2+m, j-(l-1)/2+n,0)
            end do
          end do
          auxiliary_function(i,j,1) = temp/l**2
        end do
      end do

      !Normalizes the auxiliary data in [0,1]
      auxiliary_function(:,:,1) = auxiliary_function(:,:,1)/maxval(auxiliary_function(:,:,1))

      !Defines the non smooth Andes data
      do i = 1, nlat_alt
        lat = -pio2 + (i-1)*dlat
        do j = 1, nlon_alt
          lon = -pi + (j-1)*dlon
        
          !Set topography equal to zero if it is less than 500m or if the point (lat,lon)
          !is not close enough to be considered lying in the Andes   
          if(topo(i,j,0)<=cutoff .or. lon*rad2deg>-62._r8 .or.lon*rad2deg<-82._r8 &
             .or.lat*rad2deg>15._r8 .or.lat*rad2deg<-55._r8)then
            topo(i,j,0) = 0._r8
            
          ! Removes Mount Roraima
          else if(lat*rad2deg > -5._r8 .and. lat*rad2deg < 10._r8 .and. lon*rad2deg > -68._r8 .and. lon*rad2deg < -50._r8 )then
            topo(i,j,0) = 0._r8
          end if   

          ! In this case, we consider the point (lat,lon) lying in the Andes
          ! and the topography is defined equal to 1
          if(topo(i,j,0)>0._r8)then
            topo(i,j,0) = 1._r8
          endif 
        enddo
      enddo

      ! Applies 500 iterations of Jacobi smoothing method in Andes topography data
      old = 0
      new = 1
      do k = 1, 500
        do i = 2, nlat_alt-1
          do j = 2, nlon_alt-1
            topo(i,j,new) = 0.25_r8*(topo(i+1,j,old)+topo(i-1,j,old)+topo(i,j+1,old)+topo(i,j-1,old))
          enddo
        enddo
        !Update
        aux = new
        new = old
        old = aux
      enddo

      !Normalizes the topography data in [0,1]
      topo(:,:,old) = topo(:,:,old)/maxval(topo(:,:,old))

      !Defines the density function in topo(:,:,old)    
      lambda = 0.8_r8 
      topo(:,:,old) = 1._r8/gammas**4 + (1._r8-1._r8/gammas**4)*(lambda*topo(:,:,old)+(1._r8-lambda)*auxiliary_function(:,:,1))

      filename = "topo.dat"
      open(4, file=filename)
        do i=1, nlat_alt
            write(4,*) topo(i,:,old)
        end do
      close(4)

      print*,(maxval(topo(:,:,old))/minval(topo(:,:,old)))**(0.25_r8),gammas
      !stop

      ! Update values in the table
      k = 1
      do i = 1, nlat_alt
        do j = 1, nlon_alt
            altitude_table(k,3) = topo(i,j,old)
            k = k + 1
        end do
      end do
      
      deallocate(topo,auxiliary_function)
    end subroutine andes_density_table

    !--------------------------------------------------------------------------
    ! Andes density function
    !--------------------------------------------------------------------------
    real (r8) function dens_andes(x, dens_table) result(densf)
        real (r8), dimension(1:2), intent(in) :: x
        real (r8), dimension(nlat_alt*nlon_alt,3), intent(in) :: dens_table
        real (r8) :: lat, lon

        lat = x(1)
        lon = x(2)
        lat = lat*rad2deg
        lon = lon*rad2deg
        densf = interpol_densf([lat,lon], dens_table, -90._r8, 90._r8, -180._r8, 180._r8, nlat_alt, nlon_alt)
    end function dens_andes

    !--------------------------------------------------------------------------
    ! Andes mountain topography - for shallow water tests
    ! - altitude_table stores the smooth Andes topography
    ! - Given x in lat/lon coordinates, this routine applies bilinear interpolation
    ! using altitude_table to estimate the topography value at x.
    !--------------------------------------------------------------------------
    real (r8) function smooth_andes_mountain(x, altitude_table) result(densf)
        real (r8), dimension(1:2), intent(in) :: x
        real (r8), dimension(nlat_alt*nlon_alt,3), intent(in) :: altitude_table
        real (r8) :: lat, lon
    
        lat = x(1)
        lon = x(2)
        lat = lat*rad2deg
        lon = lon*rad2deg

        densf = interpol_densf([lat,lon], altitude_table, -90._r8, 90._r8, -180._r8, 180._r8, nlat_alt, nlon_alt)
    end function smooth_andes_mountain

    !--------------------------------------------------------------------------
    ! Loads ETOPO data and generates a smooth Andes topography normalized
    ! in [0,1] in altitude_table applying Jacobi-like smoothing method. 
    ! ETOPO data must be in the directory /altitude/ 
    !--------------------------------------------------------------------------    
    subroutine generate_smooth_andes_mountain(altitude_table, iunit)
      real (r8), dimension(nlat_alt*nlon_alt,3), intent(inout) :: altitude_table
      integer, intent(in) :: iunit
      integer :: i, j, k
      integer :: old, new, aux
      real (r8) :: lat, lon
      real (r8) :: cutoff
      real (r8), dimension(:,:,:), allocatable  :: topo, topo2
      real (r8) :: dlat, dlon

      !Load ETOPO data
      call earth_elevation(iunit, altitude_table)
    
      !Allocate topography matrices
      allocate (topo(nlat_alt, nlon_alt,0:1))
    
      !Cut-off altitude
      cutoff = 500._r8
    
      !ETOPO lat-lon grid space
      dlat = pi/(nlat_alt-1)
      dlon = pi2/(nlon_alt-1)

      ! Saves the topography data
      k = 1
      do i = 1, nlat_alt
        do j = 1, nlon_alt
          topo(i,j,0) = altitude_table(k,3)
          k = k + 1
        end do
      end do

      print*,'Generating smooth Andes topography...'
      
      
      ! Define the non smooth Andes data
      do i = 1, nlat_alt
        lat = -pio2 + (i-1)*dlat
        do j = 1, nlon_alt
          lon = -pi + (j-1)*dlon
        
          !Set topography equal to zero if it is less than 500m or if the point (lat,lon)
          !is not close enough to be considered lying in the Andes   
          if(topo(i,j,0)<=cutoff .or. lon*rad2deg>-62._r8 .or.lon*rad2deg<-82._r8 & 
             .or.lat*rad2deg>15._r8 .or.lat*rad2deg<-55._r8)then
            topo(i,j,0) = 0._r8
            
          ! Remove Mount Roraima
          else if(lat*rad2deg > -5._r8 .and. lat*rad2deg < 10._r8 .and. lon*rad2deg > -68._r8 .and. lon*rad2deg < -50._r8 )then
            topo(i,j,0) = 0._r8
          end if   

          ! In this case, we considered the point lying (lat,lon) in the Andes
          ! and the topography defined equal to 1
          if(topo(i,j,0)>0._r8)then
            topo(i,j,0) = 1._r8
          endif 
        enddo
      enddo

      ! Apply 500 iterations of Jacobi smoothing method
      old = 0
      new = 1
      do k = 1, 500
        do i = 2, nlat_alt-1
          do j = 2, nlon_alt-1
            topo(i,j,new) = 0.25_r8*(topo(i+1,j,old)+topo(i-1,j,old)+topo(i,j+1,old)+topo(i,j-1,old))
          enddo
        enddo
        !Update
        aux = new
        new = old
        old = aux
      enddo

      !Normalize the topography data in [0,1] 
      topo(:,:,old) = topo(:,:,old)/maxval(topo(:,:,old))
    
      ! Update values in the table
      k = 1
      do i = 1, nlat_alt
        do j = 1, nlon_alt
          altitude_table(k,3) = topo(i,j,old)
          k = k + 1
        end do
      end do

      deallocate(topo)
    end subroutine generate_smooth_andes_mountain
    
    
    !--------------------------------------------------------------------------
    ! Creates the table for interpolation
    !--------------------------------------------------------------------------
    subroutine densftable(filename, iunit)
        character (len=100), intent(in) :: filename
        integer, intent(in) :: iunit
        integer :: nx, ny
        integer :: i, j
        real (r8) :: xmin, xmax, ymin, ymax
        real (r8) :: dx, dy
        real (r8), dimension(1:n_lat)  :: x
        real (r8), dimension(1:n_lon)  :: y
        real (r8), allocatable  :: table(:,:)

        allocate(table(n_lat*n_lon, 3))
        xmin = latmin
        xmax = latmax
        ymin = lonmin
        ymax = lonmax
        nx = n_lat
        ny = n_lon

        !Creates the grid
        dx = (xmax-xmin)/(nx-1)
        dy = (ymax-ymin)/(ny-1)
 
        do i = 1, nx
            x(i) = xmin + (i-1)*dx
        end do

        do j = 1, ny
            y(j) = ymin + (j-1)*dy
        end do

        !!!
        do i = 1, nx
            do j = 1, ny
                table((i-1)*ny + j, 1) = x(i)
                table((i-1)*ny + j, 2) = y(j)
            end do
        end do
  
        do i = 1, nx*ny
            table(i,3) = densf([table(i,1), table(i,2)]) !exact value
        end do

        !Output data to a file
        open(iunit, file=filename)
        do i=1, nx*ny
           write(iunit,*) table(i,1), table(i,2), table(i,3)
        end do
        close(iunit)
        deallocate(table)
    end subroutine densftable


end module refinter
