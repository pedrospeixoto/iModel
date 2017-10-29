!
!     ---------------------------------------------------------------
!
PROGRAM gengrid
  !
  !     Program to generate a cubed sphere grid, including cross-
  !     reference tables of adjacent faces, edges and vertices,
  !     coordinates of faces and vertices, and lengths of edges
  !     and areas of faces.
  !
  !             .....
  !            :  3  |
  !            :     |
  !             -----
  !      -----  .....  .....  -----
  !     |  5  :|  1  :|  2  :|  4  :
  !     |     :|     :|     :|     :
  !      .....  -----  -----  .....
  !             .....
  !            |  6  :
  !            |     :
  !             -----
  !
  !     Solid lines: left and bottom edges of panel
  !     Dotted lines: right and top edges of panel
  !
  !     John Thuburn Nov 2011
  !
  !     ---------------------------------------------------------------

  USE grid

  IMPLICIT NONE

  INTEGER :: igrid, i, j, ixv, p1, p2, pp, jr, iv, iv0, n, n2, &
       ie0, ie1, ie2, if0, if1, if2, if3, iv1, iv2, ix1, ix2, &
       ixmin, ifmin, if21, if22, iv11, iv12, iv21, iv22, &
       ismooth

  REAL*8 :: dlambda, lambda1, lambda2, t1, t2, long, lat, &
       x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, &
       x5, y5, z5, x6, y6, z6, piby4, xc, yc, zc, x0, y0, z0, &
       rmag, atri, aface, lmn, lmx, dmn, dmx, dav, cs, s, &
       theta, thetamin, sn, d1x, d1y, d1z, d2x, d2y, d2z

  CHARACTER (LEN=56) :: yname, atmp

  CHARACTER*56 :: ygridfile

  LOGICAL :: lfound

  ! Pi
  REAL*8, PARAMETER :: pi = 3.14159265358979323d0

  !Degrees to radians coversion (multiply to obtain conversion)
  REAL*8, PARAMETER :: deg2rad = pi / 180.0d0

  !Radians to Degrees coversion (multiply to obtain conversion)
  REAL*8, PARAMETER :: rad2deg = 1.0d0/deg2rad


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! PROGRAM ARGUMENTS


  CHARACTER*127 :: argLevels
  INTEGER :: levels

  IF (IARGC() < 1) THEN
     WRITE (*,*) "Usage: ./gengrid_cube [levels > 0]"
     CALL EXIT(-1)
  END IF

  CALL GETARG(1, argLevels)
  READ (argLevels, '(i10)') levels

  ! update variables in grid module
  ngrids = levels

  ! n x n cells on each panel
  ! Smallest and largest n

  if(ngrids > 0 )then
     n0 = 3
  else
     n0=1
     ngrids=1
  endif
  nl = n0*(2**(ngrids-1))

  nx2 = nl*nl

  nfacex = 6*nl*nl
  nedgex = 2*nfacex
  nvertx = nfacex + 2
  !print*, nl, nfacex, nedgex

  ALLOCATE(neoff(nfacex,ngrids), neofv(nvertx,ngrids), &
       nface(ngrids), nedge(ngrids), nvert(ngrids))
  ALLOCATE(fnxtf(nfacex,4,ngrids), eoff(4,nfacex,ngrids), &
       voff(nfacex,4,ngrids), fnxte(2,nedgex,ngrids), &
       vofe(2,nedgex,ngrids), fofv(nvertx,4,ngrids), &
       eofv(4,nvertx,ngrids))
  ALLOCATE(flong(nfacex,ngrids), flat(nfacex,ngrids), &
       vlong(nvertx,ngrids), vlat(nvertx,ngrids), &
       farea(nfacex,ngrids), &
       ldist(nedgex,ngrids), ddist(nedgex,ngrids))

  !     ---------------------------------------------------------------

  ! Constants

  !piby4 = ATAN(1.0d0)
  !pi = 4.0d0*piby4
  piby4=pi/real(4,8)

  !     ---------------------------------------------------------------
  !
  print*, "Basic Calculations"
  DO igrid = 1, ngrids
     print*, "Grid ", igrid, " of ", ngrids
     ! Size of panels on this grid
     n = n0*(2**(igrid-1))
     n2 = n*n
     dlambda = 0.5d0*pi/n

     nface(igrid) = 6*n2
     nedge(igrid) = 2*nface(igrid)
     nvert(igrid) = nface(igrid) + 2

     !
     ! Loop over vertices/faces of one panel
     DO j = 1, n
        lambda2 = (j-1)*dlambda - piby4
        t2 = TAN(lambda2)
        DO i = 1, n
           lambda1 = (i-1)*dlambda - piby4
           t1 = TAN(lambda1)

           !   Set up coordinates of vertices
           !   Panel 1
           !   Index of vertex
           ixv = (j-1)*n + i
           !   Cartesian coordinates of vertex
           x1 = 1.0d0/SQRT(1.0d0 + t1*t1 + t2*t2)
           y1 = x1*t1
           z1 = x1*t2
           !   Lat long coordinates of vertex
           CALL xyz2ll(x1,y1,z1,long,lat)
           vlong(ixv,igrid) = long
           vlat(ixv,igrid) = lat

           !   Panel 2
           !   Index of vertex
           ixv = ixv + n2
           !   Cartesian coordinates of vertex
           x2 = -y1
           y2 = x1
           z2 = z1
           !   Lat long coordinates of vertex
           CALL xyz2ll(x2,y2,z2,long,lat)
           vlong(ixv,igrid) = long
           vlat(ixv,igrid) = lat

           !   Panel 3
           !   Index of vertex
           ixv = ixv + n2
           !   Cartesian coordinates of vertex
           x3 = x2
           y3 = -z2
           z3 = y2
           !   Lat long coordinates of vertex
           CALL xyz2ll(x3,y3,z3,long,lat)
           vlong(ixv,igrid) = long
           vlat(ixv,igrid) = lat

           !   Panel 4
           !   Index of vertex
           ixv = ixv + n2
           !   Cartesian coordinates of vertex
           x4 = -z3
           y4 = y3
           z4 = x3
           !   Lat long coordinates of vertex
           CALL xyz2ll(x4,y4,z4,long,lat)
           vlong(ixv,igrid) = long
           vlat(ixv,igrid) = lat

           !   Panel 5
           !   Index of vertex
           ixv = ixv + n2
           !   Cartesian coordinates of vertex
           x5 = -y4
           y5 = x4
           z5 = z4
           !   Lat long coordinates of vertex
           CALL xyz2ll(x5,y5,z5,long,lat)
           vlong(ixv,igrid) = long
           vlat(ixv,igrid) = lat

           !   Panel 6
           !   Index of vertex
           ixv = ixv + n2
           !   Cartesian coordinates of vertex
           x6 = x5
           y6 = -z5
           z6 = y5
           !   Lat long coordinates of vertex
           CALL xyz2ll(x6,y6,z6,long,lat)
           vlong(ixv,igrid) = long
           vlat(ixv,igrid) = lat

           !   Set up incidence tables ignoring complications at
           !   panel edges
           DO p1 = 1, 6
              ixv = (p1 - 1)*n2 + (j-1)*n + i

              !     Edges of the face
              eoff(1,ixv,igrid) = 2*ixv - 1
              eoff(2,ixv,igrid) = 2*ixv
              eoff(3,ixv,igrid) = 2*ixv + 1
              eoff(4,ixv,igrid) = 2*ixv + 2*n
              !     Vertices of the face
              voff(ixv,1,igrid) = ixv
              voff(ixv,2,igrid) = ixv + 1
              voff(ixv,3,igrid) = ixv + n + 1
              voff(ixv,4,igrid) = ixv + n
              !     Faces neighboring this face
              fnxtf(ixv,1,igrid) = ixv - 1
              fnxtf(ixv,2,igrid) = ixv - n
              fnxtf(ixv,3,igrid) = ixv + 1
              fnxtf(ixv,4,igrid) = ixv + n
              !     Edges incident on the vertex
              eofv(1,ixv,igrid) = 2*ixv - 2
              eofv(2,ixv,igrid) = 2*(ixv-n) - 1
              eofv(3,ixv,igrid) = 2*ixv 
              eofv(4,ixv,igrid) = 2*ixv - 1

           ENDDO

        ENDDO
     ENDDO


     ! Now sort out complications at panel edges
     DO j = 1, n
        jr = n + 1 - j

        DO pp = 1, 3

           ! Odd numbered panels
           p1 = 2*pp - 1

           ! Left edge of panel p1 joins to top edge of panel p1 - 2
           ! Reverse order
           p2 = MODULO(p1 + 3, 6) + 1
           ixv = (p1 - 1)*n2 + n*(j - 1) + 1
           fnxtf(ixv,1,igrid) = p2*n2 - n + jr
           eofv(1,ixv,igrid) = 2*p2*n2 - 2*n - 1 + 2*(jr + 1)

           ! Bottom edge of panel p1 joins to top edge of panel p1 - 1
           p2 = MODULO(p1 + 4, 6) + 1
           ixv = (p1 - 1)*n2 + j
           fnxtf(ixv,2,igrid) = p2*n2 - n + j
           eofv(2,ixv,igrid) = 2*p2*n2 - 2*n - 1 + 2*j

           ! Right edge of panel p1 joins to left edge of panel p1 + 1
           p2 = MODULO(p1, 6) + 1
           ixv = (p1 - 1)*n2 + n*j
           eoff(3,ixv,igrid) = 2*(p2 - 1)*n2 + 2*(j - 1)*n + 1
           voff(ixv,2,igrid) = (p2 - 1)*n2 + (j - 1)*n + 1
           voff(ixv,3,igrid) = (p2 - 1)*n2 + j*n + 1
           fnxtf(ixv,3,igrid) = (p2 - 1)*n2 + (j - 1)*n + 1

           ! Top edge of panel p1 joins to left edge of panel p1 + 2
           ! Reverse order
           p2 = MODULO(p1 + 1, 6) + 1
           ixv = p1*n2 - n + j
           eoff(4,ixv,igrid) = 2*(p2 - 1)*n2 + 2*(jr - 1)*n + 1
           voff(ixv,3,igrid) = (p2 - 1)*n2 + (jr - 1)*n + 1
           voff(ixv,4,igrid) = (p2 - 1)*n2 + jr*n + 1
           fnxtf(ixv,4,igrid) = (p2 - 1)*n2 + (jr - 1)*n + 1

           ! Even numbered panels
           p1 = 2*pp

           ! Left edge of panel p1 joins to right edge of panel p1 - 1
           p2 = MODULO(p1 + 4, 6) + 1
           ixv = (p1 - 1)*n2 + n*(j - 1) + 1
           fnxtf(ixv,1,igrid) = (p2 - 1)*n2 + n*j
           eofv(1,ixv,igrid) = 2*(p2 - 1)*n2 + 2*n*j

           ! Bottom edge of panel p1 joins to right edge of panel p1 - 2
           ! Reverse order
           p2 = MODULO(p1 + 3, 6) + 1
           ixv = (p1 - 1)*n2 + j
           fnxtf(ixv,2,igrid) = (p2 - 1)*n2 + n*jr
           eofv(2,ixv,igrid) = 2*(p2 - 1)*n2 + 2*n*(jr + 1)

           ! Top edge of panel p1 joins to bottom edge of panel p1 + 1
           p2 = MODULO(p1, 6) + 1
           ixv = p1*n2 - n + j
           eoff(4,ixv,igrid) = 2*(p2 - 1)*n2 + 2*j
           voff(ixv,3,igrid) = (p2 - 1)*n2 + j + 1
           voff(ixv,4,igrid) = (p2 - 1)*n2 + j
           fnxtf(ixv,4,igrid) = (p2 - 1)*n2 + j

           ! Right edge of panel p1 joins to bottom edge of panel p1 + 2
           ! Reverse order
           p2 = MODULO(p1 + 1, 6) + 1
           ixv = (p1-1)*n2 + n*j
           eoff(3,ixv,igrid) = 2*(p2 - 1)*n2 + 2*jr
           voff(ixv,2,igrid) = (p2 - 1)*n2 + jr + 1
           voff(ixv,3,igrid) = (p2 - 1)*n2 + jr
           fnxtf(ixv,3,igrid) = (p2 - 1)*n2 + jr

        ENDDO

     ENDDO

     ! All faces have 4 edges and vertices
     neoff(1:nface(igrid),igrid) = 4

     ! Almost all vertices have 4 edges (exceptions dealt with below)
     neofv(1:nvert(igrid),igrid) = 4

     ! Vertices not correctly captured by the above
     ! Corner vertices only have 3 edges
     DO pp = 1, 3

        ! Bottom left of odd numbered panels
        p1 = 2*pp - 1
        ixv = (p1 - 1)*n2 + 1
        ! First edge needs to be deleted
        eofv(1,ixv,igrid) = eofv(2,ixv,igrid)
        eofv(2,ixv,igrid) = eofv(3,ixv,igrid)
        eofv(3,ixv,igrid) = eofv(4,ixv,igrid)
        eofv(4,ixv,igrid) = 0
        neofv(ixv,igrid) = 3

        ! Bottom left of even numbered panels
        p1 = 2*pp
        ixv = (p1 - 1)*n2 + 1
        ! Second edge needs to be deleted
        eofv(2,ixv,igrid) = eofv(3,ixv,igrid)
        eofv(3,ixv,igrid) = eofv(4,ixv,igrid)
        eofv(4,ixv,igrid) = 0
        neofv(ixv,igrid) = 3

     ENDDO

     ! Vertex 6*n2 + 1 is at top left of panels 1, 3, and 5
     iv = 6*n2 + 1
     lambda2 = piby4
     t2 = TAN(lambda2)
     lambda1 = - piby4
     t1 = TAN(lambda1)
     ! Cartesian coordinates of vertex
     x1 = 1.0d0/SQRT(1.0d0 + t1*t1 + t2*t2)
     y1 = x1*t1
     z1 = x1*t2
     ! Lat long coordinates of vertex
     CALL xyz2ll(x1,y1,z1,long,lat)
     vlong(iv,igrid) = long
     vlat(iv,igrid) = lat
     DO pp = 1, 3
        p1 = 2*pp - 1
        ixv = p1*n2 - n + 1
        voff(ixv,4,igrid) = iv
        eofv(pp,iv,igrid) = 2*p1*n2 - 2*n + 1
     ENDDO
     neofv(iv,igrid) = 3

     ! Vertex 6*n2 + 2 is at bottom right of panels 2, 4, and 6
     iv = 6*n2 + 2
     x1 = -x1
     y1 = -y1
     z1 = -z1
     ! Lat long coordinates of vertex
     CALL xyz2ll(x1,y1,z1,long,lat)
     vlong(iv,igrid) = long
     vlat(iv,igrid) = lat
     DO pp = 1, 3
        p1 = 2*pp
        ixv = (p1 - 1)*n2 + n
        voff(ixv,2,igrid) = iv
        eofv(pp,iv,igrid) = 2*(p1 - 1)*n2 + 2*n
     ENDDO
     neofv(iv,igrid) = 3


  ENDDO ! End of main loop over grids


  ! ---------------------------------------------------------------

  ! Now construct inverse tables
  ! First initialize entries to zero
  fnxte = 0
  vofe = 0
  fofv = 0

  print*, "Inverse tables"
  DO igrid = 1, ngrids
     print*, "Grid ", igrid, " of ", ngrids

     DO j = 1, nface(igrid)
        DO i = 1, 4
           ie1 = eoff(i,j,igrid)
           CALL addtab2(fnxte(1,1,igrid),ie1,j,2,nedgex)
           ixv = voff(j,i,igrid)
           CALL addtab(fofv(1,1,igrid),ixv,j,nvertx,4)
        ENDDO
     ENDDO

     DO j = 1, nvert(igrid)
        DO i = 1, 4
           ixv = eofv(i,j,igrid)
           IF (ixv > 0) THEN
              CALL addtab2(vofe(1,1,igrid),ixv,j,2,nedgex)
           ENDIF
        ENDDO
     ENDDO

  ENDDO


  ! Calculate geometrical quantities
  print*, "Geometrical quantities"
  DO igrid = 1, ngrids
     print*, "Grid ", igrid, " of ", ngrids
     !DO igrid = 1, ngrids

     ! Smoothing iterations
     DO ismooth = 1, nsmooth

        ! First locate face centres at barycentres of 
        ! surrounding vertices
        DO if1 = 1, nface(igrid)
           xc = 0.0d0
           yc = 0.0d0
           zc = 0.0d0
           DO i = 1, 4
              ixv = voff(if1,i,igrid)
              long = vlong(ixv,igrid)
              lat = vlat(ixv,igrid)
              CALL ll2xyz(long,lat,x1,y1,z1)
              xc = xc + x1
              yc = yc + y1
              zc = zc + z1
           ENDDO
           rmag = 1.0d0/SQRT(xc*xc + yc*yc + zc*zc)
           xc = xc*rmag
           yc = yc*rmag
           zc = zc*rmag
           CALL xyz2ll(xc,yc,zc,long,lat)
           flong(if1,igrid) = long
           flat(if1,igrid) = lat
        ENDDO

        ! Next relocate vertices at barycentres of 
        ! surrounding face centres - needed for H operator
        DO iv1 = 1, nvert(igrid)
           xc = 0.0d0
           yc = 0.0d0
           zc = 0.0d0
           DO i = 1, neofv(iv1,igrid)
              if1 = fofv(iv1,i,igrid)
              long = flong(if1,igrid)
              lat = flat(if1,igrid)
              CALL ll2xyz(long,lat,x1,y1,z1)
              xc = xc + x1
              yc = yc + y1
              zc = zc + z1
           ENDDO
           rmag = 1.0d0/SQRT(xc*xc + yc*yc + zc*zc)
           xc = xc*rmag
           yc = yc*rmag
           zc = zc*rmag
           CALL xyz2ll(xc,yc,zc,long,lat)
           vlong(iv1,igrid) = long
           vlat(iv1,igrid) = lat
        ENDDO

     ENDDO

     ! Tabulate areas
     DO if1 = 1, nface(igrid)
        long = flong(if1,igrid)
        lat = flat(if1,igrid)
        CALL ll2xyz(long,lat,x0,y0,z0)
        !   Compute face area
        aface = 0.0d0
        DO i = 1, 4
           ie1 = eoff(i,if1,igrid)
           ixv = vofe(1,ie1,igrid)
           long = vlong(ixv,igrid)
           lat = vlat(ixv,igrid)
           CALL ll2xyz(long,lat,x1,y1,z1)
           ixv = vofe(2,ie1,igrid)
           long = vlong(ixv,igrid)
           lat = vlat(ixv,igrid)
           CALL ll2xyz(long,lat,x2,y2,z2)
           CALL starea2(x0,y0,z0,x1,y1,z1,x2,y2,z2,atri)
           aface = aface + atri
        ENDDO
        farea(if1,igrid) = aface
     ENDDO

     ! Tabulate lengths of edges and distances between face centres
     ! across each edge
     lmn=5.0D0
     lmx=0.0D0
     dmn=5.0D0
     dmx=0.0D0
     dav=0.0D0
     DO ie0 = 1, nedge(igrid)
        !   Vertices at ends of this edge
        iv1 = vofe(1,ie0,igrid)
        iv2 = vofe(2,ie0,igrid)
        long = vlong(iv1,igrid)
        lat = vlat(iv1,igrid)
        CALL ll2xyz(long,lat,x1,y1,z1)
        long = vlong(iv2,igrid)
        lat = vlat(iv2,igrid)
        CALL ll2xyz(long,lat,x2,y2,z2)
        CALL spdist(x1,y1,z1,x2,y2,z2,s)
        ldist(ie0,igrid) = s
        lmn = MIN(lmn,ldist(ie0,igrid))
        lmx = MIN(lmx,ldist(ie0,igrid))
        !   Faces either side of this edge
        if1 = fnxte(1,ie0,igrid)
        if2 = fnxte(2,ie0,igrid)
        long = flong(if1,igrid)
        lat = flat(if1,igrid)
        CALL ll2xyz(long,lat,x1,y1,z1)
        long = flong(if2,igrid)
        lat = flat(if2,igrid)
        CALL ll2xyz(long,lat,x2,y2,z2)
        CALL spdist(x1,y1,z1,x2,y2,z2,s)
        ddist(ie0,igrid) = s
        dmn = MIN(dmn,ddist(ie0,igrid))
        dmx = MIN(dmx,ddist(ie0,igrid))
        dav = dav + ddist(ie0,igrid)/nedge(igrid)
     ENDDO

  ENDDO


  ! Sort FNXTF into anticlockwise order on each grid
  ! and sort EOFF to correspond to FNXTF
  ! Also sort fofv into anticlockwise order
  print*, "Sort FNXTF"
  DO igrid = 1, ngrids
     print*, "Grid ", igrid, " of ", ngrids
     !DO igrid = 1, ngrids

     DO if0 = 1, nface(igrid)
        !   Coordinates of face if0
        long = flong(if0,igrid)
        lat = flat(if0,igrid)
        CALL ll2xyz(long,lat,x0,y0,z0)
        DO ix1 = 1, 2
           !     Coordinates of IX1'th neighbour
           if1 = fnxtf(if0,ix1,igrid)
           long = flong(if1,igrid)
           lat = flat(if1,igrid)
           CALL ll2xyz(long,lat,x1,y1,z1)
           d1x = x1 - x0
           d1y = y1 - y0
           d1z = z1 - z0
           !     Find next neighbour (anticlockwise)
           thetamin = pi
           ixmin = 0
           ifmin = 0
           DO ix2 = ix1 + 1, 4
              !       Coordinates of IX2'th neighbour
              if2 = fnxtf(if0,ix2,igrid)
              long = flong(if2,igrid)
              lat = flat(if2,igrid)
              CALL ll2xyz(long,lat,x2,y2,z2)
              d2x=x2 - x0
              d2y=y2 - y0
              d2z=z2 - z0
              cs = d1x*d2x + d1y*d2y + d1z*d2z
              sn = x0*(d1y*d2z - d1z*d2y) &
                   + y0*(d1z*d2x - d1x*d2z) &
                   + z0*(d1x*d2y - d1y*d2x)
              theta = ATAN2(sn,cs)
              IF ((theta < thetamin) .AND. (theta > 0.0d0)) THEN
                 ixmin = ix2
                 ifmin = if2
                 thetamin = theta
              ENDIF
           ENDDO
           !     The face in position IXMIN belongs in position IX1+1 so swap them
           if3 = fnxtf(if0,ix1+1,igrid)
           fnxtf(if0,ix1+1,igrid) = ifmin
           fnxtf(if0,ixmin,igrid) = if3
        ENDDO

        DO ix1 = 1, 4
           if1 = fnxtf(if0,ix1,igrid)
           ix2 = ix1 - 1
           lfound = .FALSE.
           DO WHILE (.NOT. lfound)
              ix2 = ix2 + 1
              ie1 = eoff(ix2,if0,igrid)
              if21 = fnxte(1,ie1,igrid)
              if22 = fnxte(2,ie1,igrid)
              IF ((if21 + if22) == (if0 + if1)) lfound = .TRUE.
           ENDDO
           !     Edge IE2 corresponds to face IF1
           eoff(ix2,if0,igrid) = eoff(ix1,if0,igrid)
           eoff(ix1,if0,igrid) = ie1
        ENDDO

     ENDDO

     DO iv0 = 1, nvert(igrid)
        !   Coordinates of vertex iv0
        long = vlong(iv0,igrid)
        lat = vlat(iv0,igrid)
        CALL ll2xyz(long,lat,x0,y0,z0)
        DO ix1 = 1, neofv(iv0,igrid) - 2
           !     Coordinates of IX1'th face
           if1 = fofv(iv0,ix1,igrid)
           long = flong(if1,igrid)
           lat = flat(if1,igrid)
           CALL ll2xyz(long,lat,x1,y1,z1)
           d1x = x1 - x0
           d1y = y1 - y0
           d1z = z1 - z0
           !     Find next neighbour (anticlockwise)
           thetamin = pi
           ixmin = 0
           ifmin = 0
           DO ix2 = ix1 + 1, neofv(iv0,igrid)
              !       Coordinates of IX2'th neighbour
              if2 = fofv(iv0,ix2,igrid)
              long = flong(if2,igrid)
              lat = flat(if2,igrid)
              CALL ll2xyz(long,lat,x2,y2,z2)
              d2x=x2 - x0
              d2y=y2 - y0
              d2z=z2 - z0
              cs = d1x*d2x + d1y*d2y + d1z*d2z
              sn = x0*(d1y*d2z - d1z*d2y) &
                   + y0*(d1z*d2x - d1x*d2z) &
                   + z0*(d1x*d2y - d1y*d2x)
              theta = ATAN2(sn,cs)
              IF ((theta < thetamin) .AND. (theta > 0.0d0)) THEN
                 ixmin = ix2
                 ifmin = if2
                 thetamin = theta
              ENDIF
           ENDDO
           !     The face in position IXMIN belongs in position IX1+1 so swap them
           if3 = fofv(iv0,ix1+1,igrid)
           fofv(iv0,ix1+1,igrid) = ifmin
           fofv(iv0,ixmin,igrid) = if3
        ENDDO
     ENDDO

  ENDDO

  !
  ! Order VOFF so that the k'th vertex is between the
  ! k'th and (k+1)'th edges in EOFF
  print*, "Sort VOFF"
  DO igrid = 1, ngrids
     print*, "Grid ", igrid, " of ", ngrids
     !DO igrid = 1, ngrids
     DO if0 = 1, nface(igrid)
        DO ix1 = 1, neoff(if0,igrid)
           ix2 = ix1 + 1
           IF (ix2 > neoff(if0,igrid)) ix2 = 1
           ie1 = eoff(ix1,if0,igrid)
           ie2 = eoff(ix2,if0,igrid)
           ! Find the common vertex of IE1 and IE2
           iv11 = vofe(1,ie1,igrid)
           iv12 = vofe(2,ie1,igrid)
           iv21 = vofe(1,ie2,igrid)
           iv22 = vofe(2,ie2,igrid)
           IF ((iv11 == iv21) .OR. (iv11 == iv22)) THEN
              iv0 = iv11
           ELSEIF ((iv12 == iv21) .OR. (iv12 == iv22)) THEN
              iv0 = iv12
           ELSE
              PRINT *,'Common vertex not found'
              STOP
           ENDIF
           voff(if0,ix1,igrid) = iv0
        ENDDO
     ENDDO
  ENDDO


  ! Sort VOFE so that VOFE(1) -> VOFE(2) (tangent vector)
  ! is 90 degrees anticlockwise of FNXTE(1) -> FNXTE(2) (normal vector)
  print*, "Sort VOFE"
  DO igrid = 1, ngrids
     print*, "Grid ", igrid, " of ", ngrids
     !DO igrid = 1, ngrids
     DO ie0 = 1, nedge(igrid)
        if1 = fnxte(1,ie0,igrid)
        if2 = fnxte(2,ie0,igrid)
        long = flong(if1,igrid)
        lat = flat(if1,igrid)
        CALL ll2xyz(long,lat,x0,y0,z0)
        long = flong(if2,igrid)
        lat = flat(if2,igrid)
        CALL ll2xyz(long,lat,x1,y1,z1)
        d1x = x1 - x0
        d1y = y1 - y0
        d1z = z1 - z0
        iv1 = vofe(1,ie0,igrid)
        iv2 = vofe(2,ie0,igrid)
        long = vlong(iv1,igrid)
        lat = vlat(iv1,igrid)
        CALL ll2xyz(long,lat,x0,y0,z0)
        long = vlong(iv2,igrid)
        lat = vlat(iv2,igrid)
        CALL ll2xyz(long,lat,x1,y1,z1)
        d2x = x1 - x0
        d2y = y1 - y0
        d2z = z1 - z0
        sn = x0*(d1y*d2z - d1z*d2y) &
             + y0*(d1z*d2x - d1x*d2z) &
             + z0*(d1x*d2y - d1y*d2x)
        IF (sn < 0.0d0) THEN
           ! Swap the two vertices
           vofe(1,ie0,igrid) = iv2
           vofe(2,ie0,igrid) = iv1
        ENDIF
     ENDDO
  ENDDO


  ! Write out coordinates of edges for plotting
  !OPEN(44,FILE='primalgrid.dat',FORM='FORMATTED')
  !P.Peixoto edit to include grid level in name
  WRITE(ATMP,'(i3.3)') levels
  YNAME="grd/primalgrid_cube_"//trim(ATMP)//".dat"
  print*, "Grid file saved as:", trim(yname)
  OPEN(44,FILE=YNAME,FORM='FORMATTED')
  DO j = 1, nedgex
     iv = vofe(1,j,ngrids)
     long = vlong(iv,ngrids)
     lat = vlat(iv,ngrids)
     CALL ll2xyz(long,lat,x1,y1,z1)
     iv = vofe(2,j,ngrids)
     long = vlong(iv,ngrids)
     lat = vlat(iv,ngrids)
     CALL ll2xyz(long,lat,x2,y2,z2)
     WRITE(44,'(6e15.7)') x1,y1,z1,x2,y2,z2
  ENDDO
  CLOSE(44)

  if (levels <4) then
     !Output for plotting grid in GMT - vertices
     YNAME="grd/cube_"//trim(ATMP)//"_nodes.gmt"
     print*, "Grid file saved as:", trim(yname)
     OPEN(44,FILE=YNAME,FORM='FORMATTED')
     do j=1, nvert(ngrids)
        ! WRITE(22,*) 'Longitudes of vertices - all grids'
        long=vlong(j,ngrids)*rad2deg
        lat=vlat(j,ngrids)*rad2deg
        if(long>180)then
           long=long-360.0
        end if
        !print*, long, lat
        WRITE(44,'(2e15.7)') long, lat
     ENDDO
     CLOSE(44)

     !Write edges for gmt
     YNAME="grd/cube_"//trim(ATMP)//"_ed.gmt"
     print*, "Grid file saved as:", trim(yname)
     OPEN(44,FILE=YNAME)
     write(44, '(a24)') "# Primal Edges (lon, lat)"
     DO j = 1, nedgex
        write(44, '(a8,i10)') "> Edge: ", j
        iv = vofe(1,j,ngrids)
        long = vlong(iv,ngrids)*rad2deg
        lat = vlat(iv,ngrids)*rad2deg
        if(long>180)then
           long=long-360.0
        end if
        write(44,'(2f16.8)') long, lat
        !CALL ll2xyz(long,lat,x1,y1,z1)
        iv = vofe(2,j,ngrids)
        long = vlong(iv,ngrids)*rad2deg
        lat = vlat(iv,ngrids)*rad2deg
        if(long>180)then
           long=long-360.0
        end if
        write(44,'(2f16.8)') long, lat
        !CALL ll2xyz(long,lat,x2,y2,z2)
        !WRITE(44,'(6e15.7)') x1,y1,z1,x2,y2,z2
     ENDDO
     CLOSE(44)

  end if

  !OPEN(44,FILE='dualgrid.dat',FORM='FORMATTED')
  WRITE(ATMP,'(i3.3)') levels
  YNAME="grd/dualgrid_cube_"//trim(ATMP)//".dat"
  print*, "Grid file saved as:", trim(yname)
  OPEN(44,FILE=YNAME,FORM='FORMATTED')
  DO j = 1, nedgex
     if1 = fnxte(1,j,ngrids)
     long = flong(if1,ngrids)
     lat = flat(if1,ngrids)
     CALL ll2xyz(long,lat,x1,y1,z1)
     if1 = fnxte(2,j,ngrids)
     long = flong(if1,ngrids)
     lat = flat(if1,ngrids)
     CALL ll2xyz(long,lat,x2,y2,z2)
     WRITE(44,'(6e15.7)') x1,y1,z1,x2,y2,z2
  ENDDO
  CLOSE(44)

  if (levels <4) then
     !Output for plotting grid in GMT - faces
     YNAME="grd/cube_"//trim(ATMP)//"_trcc.gmt"
     print*, "Grid file saved as:", trim(yname)
     OPEN(44,FILE=YNAME,FORM='FORMATTED')
     do j=1, nface(ngrids)
        ! WRITE(22,*) 'Longitudes of vertices - all grids'
        long = flong(j,ngrids)*rad2deg
        lat = flat(j,ngrids)*rad2deg
        if(long>180)then
           long=long-360.0
        end if
        !print*, long, lat
        WRITE(44,'(2e15.7)') long, lat
     ENDDO
     CLOSE(44)

     !Write edges for gmt
     YNAME="grd/cube_"//trim(ATMP)//"_edhx.gmt"
     print*, "Grid file saved as:", trim(yname)
     OPEN(44,FILE=YNAME)
     write(44, '(a23)') "# Dual Edges (lon, lat)"
     DO j = 1, nedgex
        write(44, '(a8,i10)') "> Edge: ", j
        if1 = fnxte(1,j,ngrids)
     	long = flong(if1,ngrids)*rad2deg
     	lat = flat(if1,ngrids)*rad2deg
        if(long>180)then
           long=long-360.0
        end if
        write(44,'(2f16.8)') long, lat
        !CALL ll2xyz(long,lat,x1,y1,z1)
        if1 = fnxte(2,j,ngrids)
     	long = flong(if1,ngrids)*rad2deg
     	lat = flat(if1,ngrids)*rad2deg
        if(long>180)then
           long=long-360.0
        end if
        write(44,'(2f16.8)') long, lat
        !CALL ll2xyz(long,lat,x2,y2,z2)
        !WRITE(44,'(6e15.7)') x1,y1,z1,x2,y2,z2
     ENDDO
     CLOSE(44)

  end if


  ! Output gridmap file !Edit P. Peixoto grd/
  WRITE(ATMP,'(i10.10)') nfacex
  ygridfile="grd/gridmap_cube_"//trim(ATMP)//".dat"
  !WRITE(ygridfile,'(''grd/gridmap_cube_'',I10.10,''.dat'')') nfacex
  print*, "Grid file saved as:", trim(ygridfile)
  OPEN(22,FILE=ygridfile,FORM='UNFORMATTED')

  ! WRITE(22,*) 'GRIDMAP for NGRIDS=',NGRIDS
  WRITE(22) ngrids
  WRITE(22) nface
  WRITE(22) nedge
  WRITE(22) nvert
  ! WRITE(22,*) 'Number of edges of each face - all grids'
  WRITE(22) ((neoff(if1,igrid),            &
       if1 = 1, nface(igrid)),   &
       igrid = 1, ngrids)
  ! WRITE(22,*) 'Number of edges of each vertex - all grids'
  WRITE(22) ((neofv(iv1,igrid),            &
       iv1 = 1, nvert(igrid)),   &
       igrid=1, ngrids)
  ! WRITE(22,*) 'Faces next to each face - all grids'
  WRITE(22) (((fnxtf(if1,if2,igrid),       &
       if1 = 1, nface(igrid)),   &
       if2 = 1, 4),              &
       igrid = 1, ngrids)
  ! WRITE(22,*) 'Edges of each face - all grids'
  WRITE(22) (((eoff(ie1,if1,igrid),        &
       ie1 = 1, 4),   &
       if1 = 1, nface(igrid)),              &
       igrid = 1, ngrids)
  ! WRITE(22,*) 'Vertices of each face - all grids'
  WRITE(22) (((voff(if1,iv1,igrid),        &
       if1 = 1, nface(igrid)),   &
       iv1 = 1, 4),              &
       igrid = 1, ngrids)
  ! WRITE(22,*) 'Faces next to each edge - all grids'
  WRITE(22) (((fnxte(if2,ie1,igrid),       &
       if2 = 1, 2),              &
       ie1 = 1, nedge(igrid)),   &
       igrid = 1, ngrids)
  ! WRITE(22,*) 'Vertices of each edge - all grids'
  WRITE(22) (((vofe(iv2,ie1,igrid),        &
       iv2 = 1, 2),   &
       ie1 = 1, nedge(igrid)),              &
       igrid = 1, ngrids)
  ! WRITE(22,*) 'Faces around each vertex - all grids'
  WRITE(22) (((fofv(iv1,if2,igrid),        &
       iv1 = 1, nvert(igrid)),   &
       if2 = 1, 4),              &
       igrid = 1, ngrids)
  ! WRITE(22,*) 'Edges around each vertex - all grids'
  WRITE(22) (((eofv(ie1,iv1,igrid),        &
       ie1 = 1, 4),   &
       iv1 = 1, nvert(igrid)),              &
       igrid = 1, ngrids)
  ! WRITE(22,*) 'Longitudes of faces - all grids'
  WRITE(22) ((flong(if1,igrid),            &
       if1 = 1, nface(igrid)),   &
       igrid = 1, ngrids)
  ! WRITE(22,*) 'Latitudes of faces - all grids'
  WRITE(22) ((flat(if1,igrid),             &
       if1 = 1, nface(igrid)),   &
       igrid = 1, ngrids)
  ! WRITE(22,*) 'Longitudes of vertices - all grids'
  WRITE(22) ((vlong(iv1,igrid),            &
       iv1 = 1, nvert(igrid)),   &
       igrid = 1, ngrids)
  ! WRITE(22,*) 'Latitudes of vertices - all grids'
  WRITE(22) ((vlat(iv1,igrid),             &
       iv1 = 1, nvert(igrid)),   &
       igrid = 1, ngrids)
  ! WRITE(22,*) 'Areas of faces - all grids'
  WRITE(22) ((farea(if1,igrid),            &
       if1 = 1, nface(igrid)),    &
       igrid = 1, ngrids)
  ! WRITE(22,*) 'Lengths of edges - all grids'
  WRITE(22) ((ldist(ie1,igrid),            &
       ie1 = 1, nedge(igrid)),    &
       igrid = 1, ngrids)
  ! WRITE(22,*) 'Distance between faces across edges - all grids'
  WRITE(22) ((ddist(ie1,igrid),            &
       ie1 = 1, nedge(igrid)),    &
       igrid = 1, ngrids)

  ! WRITE zero value for 'SFC odering not initialized'
  WRITE(22) 0

  !     ---------------------------------------------------------------
  !
END PROGRAM gengrid
!
!     ==================================================================
!
SUBROUTINE LL2XYZ(LONG,LAT,X,Y,Z)
  !
  !     To convert longitude and latitude to cartesian coordinates
  !     on the unit sphere
  !
  IMPLICIT NONE
  !
  REAL*8 LONG,LAT,X,Y,Z,CLN,SLN,CLT,SLT
  !
  !     ------------------------------------------------------------------
  !
  SLN=SIN(LONG)
  CLN=COS(LONG)
  SLT=SIN(LAT)
  CLT=COS(LAT)
  !
  X=CLN*CLT
  Y=SLN*CLT
  Z=SLT
  !
  !     ------------------------------------------------------------------
  !
  RETURN
END SUBROUTINE LL2XYZ
!
!     ==================================================================
!
SUBROUTINE XYZ2LL(X,Y,Z,LONG,LAT)
  !
  !     To convert cartesian coordinates to longitude and latitude
  !
  IMPLICIT NONE
  !   
  REAL*8 X,Y,Z,LONG,LAT,PI,TLN,TLT,R
  !
  !     -------------------------------------------------------------------
  !
  PI=4.0D0*ATAN(1.0D0)
  !
  IF (X.EQ.0.0D0) THEN
     IF (Y.GE.0.0D0) THEN
        LONG=0.5D0*PI
     ELSE
        LONG=1.5D0*PI
     ENDIF
  ELSE
     TLN=Y/X
     LONG=ATAN(TLN)
     IF (X.LT.0.0D0) THEN
        LONG=LONG+PI
     ENDIF
     IF (LONG.LT.0.0D0) THEN
        LONG=LONG+2.0D0*PI
     ENDIF
  ENDIF
  !
  R=SQRT(X*X+Y*Y)
  IF (R.EQ.0.0D0) THEN
     IF (Z.GT.0.0D0) THEN
        LAT=0.5D0*PI
     ELSE
        LAT=-0.5D0*PI
     ENDIF
  ELSE
     TLT=Z/R
     LAT=ATAN(TLT)
  ENDIF
  !
  !     --------------------------------------------------------------------
  !
  RETURN
END SUBROUTINE XYZ2LL
!
!     ====================================================================
!
SUBROUTINE ADDTAB(TAB,INDEX,ENTRY,DIM1,DIM2)
  !
  !     Subroutine to add an entry to a table
  !
  INTEGER DIM1,DIM2,TAB(DIM1,DIM2),INDEX,ENTRY
  !
  !     --------------------------------------------------------------------
  !
  I=0
  !
100 CONTINUE
  I=I+1
  IF (I.GT.DIM2) THEN
     PRINT *,'**********'
     PRINT *,'TABLE FULL'
     PRINT *,'**********'
     print *,index,entry,dim1,dim2
     print *,tab(index,:)
     STOP
  ENDIF
  IF (TAB(INDEX,I).NE.0) GOTO 100
  TAB(INDEX,I)=ENTRY
  !
  !     ---------------------------------------------------------------------
  !
  RETURN
END SUBROUTINE ADDTAB

SUBROUTINE addtab2(tab,INDEX,ENTRY,dim1,dim2)

  !     Subroutine to add an entry to a table

  IMPLICIT NONE

  INTEGER, INTENT(IN)                  :: dim1
  INTEGER, INTENT(IN)                  :: dim2
  INTEGER, INTENT(IN OUT)              :: tab(dim1,dim2)
  INTEGER, INTENT(IN OUT)              :: INDEX
  INTEGER, INTENT(IN)                  :: ENTRY

  INTEGER :: i


  !     --------------------------------------------------------------------

  i=0

100 CONTINUE
  i=i+1
  IF (i > dim1) THEN
     PRINT *,'**********'
     PRINT *,'TABLE FULL'
     PRINT *,'**********'
     STOP
  END IF
  IF (tab(i,INDEX) /= 0) GO TO 100
  tab(i,INDEX)=ENTRY

  !     ---------------------------------------------------------------------

  RETURN
END SUBROUTINE addtab2


!
!     ====================================================================
!
SUBROUTINE STAREA2(X0,Y0,Z0,X1,Y1,Z1,X2,Y2,Z2,AREA)
  !
  !     Calculate the area of the spherical triangle whose corners
  !     have Cartesian coordinates (X0,Y0,Z0), (X1,Y1,Z1), (X2,Y2,Z2)
  !     The formula below is more robust to roundoff error than the
  !     better known sum of angle - PI formula
  !
  IMPLICIT NONE
  !
  REAL*8 X0,Y0,Z0,X1,Y1,Z1,X2,Y2,Z2, &
       D0,D1,D2,S,T0,T1,T2,T3,AREA
  !
  !
  !     Distances between pairs of points
  CALL SPDIST(X0,Y0,Z0,X1,Y1,Z1,D2)
  CALL SPDIST(X1,Y1,Z1,X2,Y2,Z2,D0)
  CALL SPDIST(X2,Y2,Z2,X0,Y0,Z0,D1)
  !
  !     Half perimeter
  S=0.5D0*(D0+D1+D2)
  !
  !     Tangents
  T0 = TAN(0.5D0*(S-D0))
  T1 = TAN(0.5D0*(S-D1))
  T2 = TAN(0.5D0*(S-D2))
  T3 = TAN(0.5D0*S)
  !
  !     Area
  AREA = 4.0D0*ATAN(SQRT(T0*T1*T2*T3))
  !
  RETURN
END SUBROUTINE STAREA2
!
!     ===================================================================
!
SUBROUTINE SPDIST(X1,Y1,Z1,X2,Y2,Z2,S)
  !
  !     Calculate the spherical distance S between two points with Cartesian
  !     coordinates (X1,Y1,Z1), (X2,Y2,Z2) on the unit sphere

  IMPLICIT NONE

  REAL*8 X1, Y1, Z1, X2, Y2, Z2, S, DX, DY, DZ, AD


  DX = X2 - X1
  DY = Y2 - Y1
  DZ = Z2 - Z1
  AD = SQRT(DX*DX + DY*DY + DZ*DZ)
  S = 2.0D0*ASIN(0.5D0*AD)


  RETURN
END SUBROUTINE SPDIST
!
!     ===================================================================

