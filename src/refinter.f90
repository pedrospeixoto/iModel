module refinter
    !parameters
    use constants, only: &
        pi, &
        datadir, &
        pio2, &
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
        nx_andes, &
        ny_andes, &
        nlat_alt, &
        nlon_alt, &
        latmin_andes, &
        latmax_andes, &
        lonmin_andes, &
        lonmax_andes
    implicit none

contains
    !--------------------------------------------------------------------------
    !Compute the density function
    !--------------------------------------------------------------------------
    real (kind=8) function densf(x) result(dens_f)
        real (kind=8), dimension(1:2), intent(in) :: x
        real (kind=8), dimension(1:3) :: p
        real (kind=8), dimension(1:3) :: pc
        real (kind=8) :: coslat
        real (kind=8) :: lat
        real (kind=8) :: lon
        real (kind=8) :: latc
        real (kind=8) :: lonc
        real (kind=8) :: gammas
        real (kind=8) :: epsilons
        real (kind=8) :: dists
        real (kind=8) :: maxdist
        real (kind=8) :: sx


        !See Voronoi Tessellations and their application to climate and global modeling
        ! by Lili Ju, Todd Ringler and Max Gunzburger
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
    ! Creates the table for interpolation
    !--------------------------------------------------------------------------
    subroutine densftable(filename, iunit)
        character (len=100), intent(in) :: filename
        integer, intent(in) :: iunit
        integer :: nx, ny
        integer :: i, j
        real (kind=8) :: xmin, xmax, ymin, ymax
        real (kind=8) :: dx, dy
        real (kind=8), dimension(1:n_lat)  :: x
        real (kind=8), dimension(1:n_lon)  :: y
        real (kind=8), dimension(n_lat*n_lon, 3)  :: table

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
   
    end subroutine densftable


    !--------------------------------------------------------------------------
    ! Bilinear interpolation
    !--------------------------------------------------------------------------
    real (kind=8) function interpol_densf(p, densftable, xmin, xmax, ymin, ymax, nx, ny) result(ff)
        integer, intent(in) :: nx, ny
        real (kind=8), intent(in) :: xmin, xmax, ymin, ymax
        real (kind=8), dimension(1:2), intent(in) :: p
        real (kind=8), dimension(nx*ny,3), intent(in) :: densftable
        integer :: i, j
        integer :: i1, i2, j1, j2
        real (kind=8) :: f11, f12, f21, f22
        real (kind=8) :: x1, x2, y1, y2, x0, y0
        real (kind=8) :: dx, dy

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
    ! Earth elevation data
    !--------------------------------------------------------------------------
    subroutine earth_elevation(iunit, alt_table)
        integer :: iunit
        real (kind=8), dimension(nlat_alt*nlon_alt,3), intent(inout) :: alt_table
        integer :: i, j
        character (len=100) :: filename
        character (len=100) :: buffer

        !filename = trim(altdir)//"elevation_table_180x360.dat"
        !filename = trim(altdir)//"elevation_table_360x720.dat"
        filename = trim(altdir)//"elevation_table_720x1440.dat"

        open(iunit, file=filename)
        read(iunit,*) buffer
        read(iunit,*) buffer
        read(iunit,*) buffer
        read(iunit,*) ((alt_table  (i,j), j=1,3), i=1,nlat_alt*nlon_alt)
        close(iunit)

    end subroutine earth_elevation


    !--------------------------------------------------------------------------
    ! Moving average
    !--------------------------------------------------------------------------
    subroutine moving_average1(alt_table, gammas)
        real (kind=8), dimension(nlat_alt*nlon_alt,3), intent(inout) :: alt_table
        real (kind=8), intent(in) :: gammas
        integer :: i, j, k, m, n
        integer :: l, l2, l3, l4,l2aux,l3aux
        real (kind=8), dimension(:,:), allocatable  :: aux
        real (kind=8), dimension(:,:), allocatable  :: aux2
        real (kind=8) :: temp, alpha

        allocate (aux(nlat_alt, nlon_alt))
        allocate (aux2(nlat_alt, nlon_alt))
   
        k = 1
        do i = 1, nlat_alt
            do j = 1, nlon_alt
         
                aux(i,j) = alt_table(k,3)
                aux2(i,j) = alt_table(k,3)
                k = k + 1
            end do
        end do

        l = 50
        l2 = 90  !90
        l3 = 90 !90
        l4 = 20

        do i = l, nlat_alt-l
            do j = l, nlon_alt-l

                if(aux(i,j) == 500._r8/6320._r8)then

                    temp =  0._r8
                    do m = 0,l2-1
                        do n = 0,l2-1
                            temp = temp + aux(i-(l2-1)/2+m, j-(l2-1)/2+n)
                        end do
                    end do
                    temp = temp - aux(i,j)
                    aux2(i,j) = temp/(l2**2-1)
	
                end if
 
                alpha = 5500._r8/6320._r8
                !alpha = alpha**(0.25)
                if(aux(i,j) == alpha*1._r8/gammas**4)then
                    ! if(aux(i,j) == 1._r8/gammas**4)then

                    temp =  0._r8
                    do m = 0,l3-1
                        do n = 0,l3-1
                            temp = temp + aux(i-(l3-1)/2+m, j-(l3-1)/2+n)
                        end do
                    end do
                    temp = temp - aux(i,j)
                    aux2(i,j) = temp/(l3**2-1)

                end if


            !if(aux(i,j) > 500._r8/6320._r8)then

            !         temp =  0._r8
            !          do m = 0,l4-1
            !          do n = 0,l4-1
            !           temp = temp + aux(i-(l4-1)/2+m, j-(l4-1)/2+n)
            !        end do
            !     end do
            !    temp = temp - aux(i,j)
             !   aux2(i,j) = temp/(l4**2-1)
 
            !end if




            end do
        end do

        k = 1
        do i = 1, nlat_alt
            do j = 1, nlon_alt
                alt_table(k,3) = aux2(i,j)
                k = k + 1
            end do
        end do

        deallocate(aux)
        deallocate(aux2)

    end subroutine moving_average1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine moving_average2(alt_table, gammas)
        real (kind=8), dimension(nlat_alt*nlon_alt,3), intent(inout) :: alt_table
        real (kind=8), intent(in) :: gammas
        integer :: i, j, k, m, n
        integer :: l, l2
        real (kind=8), dimension(:,:), allocatable  :: aux
        real (kind=8), dimension(:,:), allocatable  :: aux2
        real (kind=8) :: temp


        allocate (aux(nlat_alt, nlon_alt))
        allocate (aux2(nlat_alt, nlon_alt))
   

        k = 1
        do i = 1, nlat_alt
            do j = 1, nlon_alt
         
                aux(i,j) = alt_table(k,3)
                aux2(i,j) = alt_table(k,3)
                k = k + 1
            end do
        end do

        l = 17
        do i = l, nlat_alt-l
            do j = l, nlon_alt-l


                if(aux(i,j) > 500._r8/6320._r8)then

                    temp =  0._r8
                    do m = 0,l-1
                        do n = 0,l-1
                            temp = temp + aux(i-(l-1)/2+m, j-(l-1)/2+n)
                        end do
                    end do
                    temp = temp - aux(i,j)
                    aux2(i,j) = temp/(l**2-1)

                end if

            end do
        end do

        k = 1
        do i = 1, nlat_alt
            do j = 1, nlon_alt
                alt_table(k,3) = aux2(i,j)
                k = k + 1
            end do
        end do

        deallocate(aux)
        deallocate(aux2)

    end subroutine moving_average2


    !--------------------------------------------------------------------------
    ! Andes density data for interpolation
    !--------------------------------------------------------------------------
    subroutine andes_density_table(alt_table)
        real (kind=8), dimension(nlat_alt*nlon_alt,3), intent(inout) :: alt_table
        integer :: i, j
        real (kind=8) :: maxx, minn, lat, lon, latc, lonc
        real (kind=8) :: dists, maxdist, gammas, epsilons
        real (kind=8), dimension(1:3) :: p
        real (kind=8), dimension(1:3) :: pc
        real (kind=8) :: coslat, alpha

        latc=-20.0*pi/180._r8
        lonc=-60.0*pi/180._r8
        maxdist = pi/6._r8
        epsilons = pi/12._r8
        alpha = 5500._r8/6320._r8
        !alpha = alpha**(0.25)
        minn = 0._r8
        maxx = 0._r8

        latc=-20.0*pi/180._r8
        lonc=-60.0*pi/180._r8
        coslat = dcos(latc)
        pc(1) = coslat*dcos(lonc)
        pc(2) = coslat*dsin(lonc)
        pc(3) = dsin(latc)



        do i = 1, nlat_alt*nlon_alt
            lat = alt_table(i,1)
            lon = alt_table(i,2)
            lat = lat*deg2rad
            lon = lon*deg2rad

            coslat = dcos(lat)
            p(1) = coslat*dcos(lon)
            p(2) = coslat*dsin(lon)
            p(3) = dsin (lat)

            !Distance to center
            dists = dsqrt((p(1)-pc(1))**2+(p(2)-pc(2))**2+(p(3)-pc(3))**2)
            ! dists = dsqrt((lat-latc)**2+(lon-lonc)**2)
      
            if(dists<=maxdist)then
                if(alt_table(i,3) < minn) then
                    minn = alt_table(i,3)

                else if(alt_table(i,3) > maxx) then
                    maxx = alt_table(i,3)
                end if
            end if
        end do


        gammas = 3._r8

        do i = 1, nlat_alt*nlon_alt
            lat = alt_table(i,1)
            lon = alt_table(i,2)
            lat = lat*deg2rad
            lon = lon*deg2rad

            coslat = dcos(lat)
            p(1) = coslat*dcos(lon)
            p(2) = coslat*dsin(lon)
            p(3) = dsin(lat)

            !Distance to center
            dists = dsqrt((p(1)-pc(1))**2+(p(2)-pc(2))**2+(p(3)-pc(3))**2)
            !dists = dsqrt((lat-latc)**2+(lon-lonc)**2)

            if(dists<=maxdist)then
                if(alt_table(i,3) <= 500._r8 )then
                    alt_table(i,3) = 500._r8
                end if
 
                alt_table(i,3) = alt_table(i,3)/maxx

            else
                alt_table(i,3) = alpha*1._r8/gammas**4
               !alt_table(i,3) = 1._r8/gammas**4
            end if

        end do


        call moving_average1(alt_table, gammas)
        call moving_average2(alt_table, gammas)

    end subroutine andes_density_table


    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------


    subroutine andes_density_table2(alt_table)
        real (kind=8), dimension(nlat_alt*nlon_alt,3), intent(inout) :: alt_table
        integer :: i, j
        real (kind=8) :: maxx, minn, lat, lon, latc, lonc
        real (kind=8) :: dists, maxdist, gammas, epsilons
        real (kind=8), dimension(1:3) :: p
        real (kind=8), dimension(1:3) :: pc
        real (kind=8) :: coslat

        latc=-20.0*pi/180._r8
        lonc=-60.0*pi/180._r8
        maxdist = pi/6._r8
        epsilons = pi/12._r8

        coslat = dcos(latc)
        pc(1) = coslat*dcos(lonc)
        pc(2) = coslat*dsin(lonc)
        pc(3) = dsin(latc)

        minn = 0._r8
        maxx = 0._r8


        do i = 1, nlat_alt*nlon_alt
            lat = alt_table(i,1)
            lon = alt_table(i,2)
            lat = lat*deg2rad
            lon = lon*deg2rad

            coslat = dcos(lat)
            p(1) = coslat*dcos(lon)
            p(2) = coslat*dsin(lon)
            p(3) = dsin(lat)
      
            dists = dsqrt((p(1)-pc(1))**2+(p(2)-pc(2))**2+(p(3)-pc(3))**2)

            if(dists<=maxdist)then
                if(alt_table(i,3) > maxx) then
                    maxx = alt_table(i,3)
                end if
            end if
        end do


        do i = 1, nlat_alt*nlon_alt
            if(alt_table(i,3) <= 200._r8) then
                alt_table(i,3) = 200._r8
            end if
            alt_table(i,3) = alt_table(i,3)/maxx
        end do
    end subroutine andes_density_table2


    subroutine andes_density_table3(alt_table)
        real (kind=8), dimension(nlat_alt*nlon_alt,3), intent(inout) :: alt_table
        integer :: i, j, k, n, m, l,l2,l3
        real (kind=8) :: maxx, lat, lon, latc, lonc
        real (kind=8) :: dists, maxdist, gammas, epsilons
        real (kind=8), dimension(1:3) :: p
        real (kind=8), dimension(1:3) :: pc
        real (kind=8) :: coslat, a, b, sx
        real (kind=8), dimension(:,:), allocatable  :: aux
        real (kind=8), dimension(:,:), allocatable  :: aux2
        real (kind=8) :: temp


        allocate (aux(nlat_alt, nlon_alt))
        allocate (aux2(nlat_alt, nlon_alt))

        latc=-20.0*pi/180._r8
        lonc=-60.0*pi/180._r8
        maxdist = pi/6._r8
        epsilons = pi/12._r8
        a = 500._r8/6320._r8
        b = 5500._r8/6320._r8
        a = a**(0.25)
        b = b**(0.25)
        maxx = 0._r8

        latc=-20.0*pi/180._r8
        lonc=-60.0*pi/180._r8
        coslat = dcos(latc)
        pc(1) = coslat*dcos(lonc)
        pc(2) = coslat*dsin(lonc)
        pc(3) = dsin(latc)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do i = 1, nlat_alt*nlon_alt
            lat = alt_table(i,1)
            lon = alt_table(i,2)
            lat = lat*deg2rad
            lon = lon*deg2rad

            coslat = dcos(lat)
            p(1) = coslat*dcos(lon)
            p(2) = coslat*dsin(lon)
            p(3) = dsin (lat)

            !Distance to center
            dists = dsqrt((p(1)-pc(1))**2+(p(2)-pc(2))**2+(p(3)-pc(3))**2)
            !dists = dsqrt((lat-latc)**2+(lon-lonc)**2)
      
            if(dists<=maxdist)then
                if(alt_table(i,3) > maxx) then
                    maxx = alt_table(i,3)
                end if
            end if
        end do


        gammas = 3._r8

        do i = 1, nlat_alt*nlon_alt
            lat = alt_table(i,1)
            lon = alt_table(i,2)
            lat = lat*deg2rad
            lon = lon*deg2rad

            coslat = dcos(lat)
            p(1) = coslat*dcos(lon)
            p(2) = coslat*dsin(lon)
            p(3) = dsin(lat)

            !Distance to center
            dists = dsqrt((p(1)-pc(1))**2+(p(2)-pc(2))**2+(p(3)-pc(3))**2)
            !dists = dsqrt((lat-latc)**2+(lon-lonc)**2)

            if(dists<=maxdist)then
                if(alt_table(i,3) <= 500._r8 )then
                    alt_table(i,3) = 500._r8
                end if
 
                alt_table(i,3) = alt_table(i,3)/maxx

            else if(dists<=maxdist+epsilons)then
                !Distance to center metric
                sx=(dists-maxdist)/epsilons
                alt_table(i,3) = ((1._r8-sx)*a*gammas+sx*b)**4/gammas**4
            else
                alt_table(i,3) = b**4/gammas**4
            end if
        end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        k = 1
        do i = 1, nlat_alt
            do j = 1, nlon_alt
                aux(i,j) = alt_table(k,3)
                aux2(i,j) = alt_table(k,3)
                k = k + 1
            end do
        end do

        l = 40
        l2 = 17
        l3 = 40

        do i = l3, nlat_alt-l3
            do j = l3, nlon_alt-l3
                if(aux(i,j) == 500._r8/6320._r8)then
                    temp =  0._r8
                    do m = 0,l-1
                        do n = 0,l-1
                            temp = temp + aux(i-(l-1)/2+m, j-(l-1)/2+n)
                        end do
                    end do
                    aux2(i,j) = temp/l**2
                end if
            end do
        end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        k = 1
        do i = 1, nlat_alt
            do j = 1, nlon_alt
                lat = alt_table(k,1)
                lon = alt_table(k,2)
                lat = lat*deg2rad
                lon = lon*deg2rad

                coslat = dcos(lat)
                p(1) = coslat*dcos(lon)
                p(2) = coslat*dsin(lon)
                p(3) = dsin(lat)

                !Distance to center
                dists = dsqrt((p(1)-pc(1))**2+(p(2)-pc(2))**2+(p(3)-pc(3))**2)
                if(dists<=maxdist)then
                    alt_table(k,3) = aux2(i,j)
                    aux(i,j) = aux2(i,j)
                end if
                k = k + 1
            end do
        end do


        do i = l3, nlat_alt-l3
            do j = l3, nlon_alt-l3
                if(aux(i,j) > 500._r8/6320._r8)then
                    temp =  0._r8
                    do m = 0,l2-1
                        do n = 0,l2-1
                            temp = temp + aux(i-(l2-1)/2+m, j-(l2-1)/2+n)
                        end do
                    end do
                    aux2(i,j) = temp/l2**2
                end if
            end do
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        k = 1
        do i = 1, nlat_alt
            do j = 1, nlon_alt
                lat = alt_table(k,1)
                lon = alt_table(k,2)
                lat = lat*deg2rad
                lon = lon*deg2rad

                coslat = dcos(lat)
                p(1) = coslat*dcos(lon)
                p(2) = coslat*dsin(lon)
                p(3) = dsin(lat)

                !Distance to center
                dists = dsqrt((p(1)-pc(1))**2+(p(2)-pc(2))**2+(p(3)-pc(3))**2)
                !dists = dsqrt((lat-latc)**2+(lon-lonc)**2)
                if(dists<=maxdist)then
                    alt_table(k,3) = aux2(i,j)
                end if
                k = k + 1
            end do
        end do

        deallocate(aux)
        deallocate(aux2)
    end subroutine andes_density_table3


    !--------------------------------------------------------------------------
    ! Andes density function
    !--------------------------------------------------------------------------
    real (kind=8) function dens_andes(x, dens_table) result(densf)
        real (kind=8), dimension(1:2), intent(in) :: x
        real (kind=8), dimension(nlat_alt*nlon_alt,3), intent(in) :: dens_table
        real (kind=8) :: lat, lon
    
        lat = x(1)
        lon = x(2)
        lat = lat*rad2deg
        lon = lon*rad2deg

        densf = interpol_densf([lat,lon], dens_table, -90._r8, 90._r8, -180._r8, 180._r8, nlat_alt, nlon_alt)

    end function dens_andes



    real (kind=8) function dens_andes1(x, dens_table) result(densf)
        real (kind=8), dimension(1:2), intent(in) :: x
        real (kind=8), dimension(nlat_alt*nlon_alt,3), intent(in) :: dens_table
        real (kind=8) :: lat, lon, latc, lonc
        real (kind=8) :: dists, maxdist, gammas, epsilons, coslat
        real (kind=8), dimension(1:3) :: p
        real (kind=8), dimension(1:3) :: pc

        latc=-20.0*pi/180._r8
        lonc=-60.0*pi/180._r8
        maxdist = pi/6._r8
        epsilons = pi/12._r8
        gammas = 3._r8

        lat = x(1)
        lon = x(2)
        coslat = dcos(lat)
        p(1) = coslat*dcos(lon)
        p(2) = coslat*dsin(lon)
        p(3) = dsin(lat)

        coslat = dcos(latc)
        pc(1) = coslat*dcos(lonc)
        pc(2) = coslat*dsin(lonc)
        pc(3) = dsin(latc)

        dists = dsqrt((p(1)-pc(1))**2+(p(2)-pc(2))**2+(p(3)-pc(3))**2)

        if(dists<=maxdist+epsilons)then
            lat = lat*rad2deg
            lon = lon*rad2deg
            densf = interpol_densf([lat,lon], dens_table, -90._r8, 90._r8, -180._r8, 180._r8, nlat_alt, nlon_alt)
        else
            densf = 1._r8/gammas**4
        end if

    end function dens_andes1

    !--------------------------------------------------------------------------
    ! Andes density function 2
    !--------------------------------------------------------------------------


    real (kind=8) function dens_andes2(x, dens_table) result(dens_f)
        real (kind=8), dimension(1:2), intent(in) :: x
        real (kind=8), dimension(1:3) :: p, pc
        real (kind=8), dimension(nlat_alt*nlon_alt,3), intent(in) :: dens_table
        real (kind=8) :: lat, lon,latc, lonc, coslat
        real (kind=8) :: epsilons, gammas, sx, maxdist, dists, a, b

        latc=-20.0*pi/180._r8
        lonc=-60.0*pi/180._r8
        lat = x(1)
        lon = x(2)


        !Density function parameters
        gammas = 3._r8
        epsilons = pi/12._r8
        maxdist = pi/6._r8
        a = 500._r8/6320._r8
        a = a**(0.25)
        b = 5500._r8/6320._r8
        b = b**(0.25)

        coslat = dcos(lat)
        p(1) = coslat*dcos(lon)
        p(2) = coslat*dsin(lon)
        p(3) = dsin(lat)

        coslat = dcos(latc)
        pc(1) = coslat*dcos(lonc)
        pc(2) = coslat*dsin(lonc)
        pc(3) = dsin(latc)

        dists = dsqrt((p(1)-pc(1))**2+(p(2)-pc(2))**2+(p(3)-pc(3))**2)
        !Distance to center metric
        sx=(dists-maxdist)/epsilons

        if(dists<=maxdist)then
            lat = lat*rad2deg
            lon = lon*rad2deg
            dens_f = interpol_densf([lat,lon], dens_table, -90._r8, 90._r8, -180._r8, 180._r8, nlat_alt, nlon_alt)
        elseif((dists<=maxdist+epsilons) .and. (dists>maxdist))then
            dens_f = ((1._r8-sx)*a*gammas+sx*b)**4/gammas**4
        else
            dens_f=b**4/gammas**4
        end if

    end function dens_andes2

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    subroutine andean_mountain_data(alt_table, iunit)
        real (kind=8), dimension(nlat_alt*nlon_alt,3), intent(inout) :: alt_table
        integer :: i, j, k, n, m, l,l2,l3
        integer :: iunit  !, iunit2
        character (len=100) :: filename
        real (kind=8) :: maxx, lat, lon, latc, lonc
        real (kind=8) :: dists, maxdist
        real (kind=8), dimension(1:3) :: p
        real (kind=8), dimension(1:3) :: pc
        real (kind=8) :: coslat, a, b, sx
        real (kind=8), dimension(:,:), allocatable  :: aux
        real (kind=8), dimension(:,:), allocatable  :: aux2
        real (kind=8) :: temp, corte

        allocate (aux(nlat_alt, nlon_alt))
        allocate (aux2(nlat_alt, nlon_alt))
        corte = 200.d0
        aux = corte
        aux2 = corte
        
        call earth_elevation(iunit, alt_table)
        latc=-20.0*pi/180._r8
        lonc=-60.0*pi/180._r8
        maxdist = pi/6._r8

        coslat = dcos(latc)
        pc(1) = coslat*dcos(lonc)
        pc(2) = coslat*dsin(lonc)
        pc(3) = dsin(latc)

        do i = 1, nlat_alt*nlon_alt
            lat = alt_table(i,1)
            lon = alt_table(i,2)
            lat = lat*deg2rad
            lon = lon*deg2rad

            coslat = dcos(lat)
            p(1) = coslat*dcos(lon)
            p(2) = coslat*dsin(lon)
            p(3) = dsin(lat)

            !Distance to center
            dists = dsqrt((p(1)-pc(1))**2+(p(2)-pc(2))**2+(p(3)-pc(3))**2)

            if(dists<=maxdist)then
                if(alt_table(i,3) <= corte .or. lon*rad2deg > -60.d0)then
                  alt_table(i,3) = corte
                  ! Remove Mount Roraima
                else if(lat*rad2deg > 0.d0 .and. lat*rad2deg < 10.d0 .and. lon*rad2deg > -68.d0 .and. lon*rad2deg < -58.d0 )then
                  alt_table(i,3) = corte
                end if
                
            else 
                alt_table(i,3) = corte
            end if
        end do

        !Output data to a file
        !filename = "topo.dat"
        !open(13, file=filename)
        !  do i=1, nlat_alt*nlon_alt
        !    write(13,*) alt_table(i,1), alt_table(i,2), alt_table(i,3)
        !  end do
        !close(13)
        !print*, minval(alt_table(:,3)),maxval(alt_table(:,3))

        k = 1
        do i = 1, nlat_alt
            do j = 1, nlon_alt
                aux(i,j) = alt_table(k,3)
                aux2(i,j) = alt_table(k,3)
                k = k + 1
            end do
        end do

        l = 40
        l2 = 17
        l3 = 40
        do i = l3, nlat_alt-l3
          do j = l3, nlon_alt-l3
            if(aux(i,j) == corte)then
              temp =  0._r8
              do m = 0,l-1
                do n = 0,l-1
                  temp = temp + aux(i-(l-1)/2+m, j-(l-1)/2+n)
                end do
              end do
              aux2(i,j) = temp/l**2
            end if
          end do
        end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        k = 1
        do i = 1, nlat_alt
          do j = 1, nlon_alt
            alt_table(k,3) = aux2(i,j)
            aux(i,j) = aux2(i,j)
            k = k + 1
          end do
        end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do i = l3, nlat_alt-l3
          do j = l3, nlon_alt-l3
            if(aux(i,j) > corte)then
              temp =  0._r8
              do m = 0,l2-1
                do n = 0,l2-1
                  temp = temp + aux(i-(l2-1)/2+m, j-(l2-1)/2+n)
                end do
              end do
              aux2(i,j) = temp/l2**2
            end if
          end do
        end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        k = 1
        do i = 1, nlat_alt
          do j = 1, nlon_alt
            alt_table(k,3) = aux2(i,j) - corte
            k = k + 1
          end do
        end do

        !print*, minval(alt_table(:,3)),maxval(alt_table(:,3))
        !Output data to a file
        !filename = "topo2.dat"
        !open(11, file=filename)
        !  do i=1, nlat_alt*nlon_alt
        !    write(11,*) alt_table(i,1), alt_table(i,2), alt_table(i,3)
        !  end do
        !close(11)

        deallocate(aux)
        deallocate(aux2)
    end subroutine andean_mountain_data

    !--------------------------------------------------------------------------
    ! Andes mountain
    !--------------------------------------------------------------------------
    real (kind=8) function smooth_andean_mountain(x, alt_table) result(densf)
        real (kind=8), dimension(1:2), intent(in) :: x
        real (kind=8), dimension(nlat_alt*nlon_alt,3), intent(in) :: alt_table
        real (kind=8) :: lat, lon
    
        lat = x(1)
        lon = x(2)
        lat = lat*rad2deg
        lon = lon*rad2deg

        densf = interpol_densf([lat,lon], alt_table, -90._r8, 90._r8, -180._r8, 180._r8, nlat_alt, nlon_alt)
    end function smooth_andean_mountain

end module refinter
