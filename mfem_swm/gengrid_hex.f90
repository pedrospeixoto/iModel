
PROGRAM gengrid

  ! Code converted using TO_F90 by Alan Miller
  ! Date: 2013-11-20  Time: 15:37:58

  !     Program to generate a nested set of hexagonal-
  !     icosahedral grids on the sphere, including cross-
  !     reference tables of adjacent faces, edges and vertices,
  !     coordinates of faces and vertices, and lengths of edges
  !     and areas of faces.

  !     John Thuburn 1/November/1994

  !     Modified to write GRIDMAP file unformatted

  !     John Thuburn 2/May/1997

  !     Adapted to optimise on both Heikes-Randall criterion
  !     and centroidal criterion

  !     John Thuburn October 2011

  !     Modified on Apr 2015 by P. Peixoto - Stop criteria, output format

  !     Modified on Aug 2015 by P. Peixoto - New way of grid outputs
  
  !     ---------------------------------------------------------------

  USE grid

  IMPLICIT NONE


  REAL*8 weight1,weight2

  INTEGER :: if1,if2,ie1,iv1,if3,nf2,igrid,iedge,iv2,iv0,  &
       igridm,ie2,ie3,if11,if12,if21,if22,if31,if32,  &
       ivn,ien,ifn1,ifn2,ifn3,ix,ie0,if0,ifn,ie11,ie12,  &
       ie21,ie22,ie31,ie32,NE,nf,ifac,iter,emax,fmax, ix1,ix2,iv11,iv12,iv21,iv22

  REAL*8 pi,zz,phi0,long,lat,x,y,z,px(3),py(3),pz(3),  &
       d1x,d1y,d1z,d2x,d2y,d2z,mag,vx,vy,vz,  &
       x1,y1,z1,x2,y2,z2,x0,y0,z0,atot,aface,da,s12,s22,s32,  &
       amn,amx,sfac,emn,emx,dmn,dmx,cs,sn,gmn,gmx,x3,y3,z3,  &
       w1,w2,w3,tot,dav,cost,totcst,ftol,dd,maxdd,disp, mxdisp,&
       bdisp,s, tmp, totcstold
  CHARACTER (LEN=56) :: yname, atmp
  LOGICAL :: lswap,lprint

  COMMON /comwgt/ weight1,weight2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! PROGRAM ARGUMENTS


  CHARACTER*127 :: argLevels
  INTEGER :: levels

  IF (IARGC() < 1) THEN
     WRITE (*,*) "Usage: ./gengrid_hex [glevel from 0-10]"
     CALL EXIT(-1)
  END IF

  CALL GETARG(1, argLevels)
  READ (argLevels, '(i10)') levels

  ! update variables in grid module
  ngrids = levels+1

  nfacex=5*2**(2*ngrids-1)+2
  nedgex=15*2**(2*ngrids-1)
  nvertx=5*2**(2*ngrids)


  ALLOCATE(nface(ngrids),nedge(ngrids),nvert(ngrids))

  ALLOCATE(neoff(nfacex,ngrids), neofv(nvertx,ngrids))

  ALLOCATE(fnxtf(nfacex,6,ngrids),eoff(6,nfacex,ngrids),  &
       fnxte(2,nedgex,ngrids),feofe(nedgex,2,ngrids),  &
       vofe(2,nedgex,ngrids),voff(nfacex,6,ngrids),  &
       fofv(nvertx,3,ngrids),eofv(3,nvertx,ngrids))

  ALLOCATE(flong(nfacex,ngrids),flat(nfacex,ngrids), farea(nfacex,ngrids),  &
       vlong(nvertx,ngrids),vlat(nvertx,ngrids),  &
       ldist(nedgex,ngrids),ddist(nedgex,ngrids),  &
       gdist(nedgex,ngrids),coeff(nvertx,3,ngrids))

  !     ---------------------------------------------------------------

  !     Initialize all tables to zero
  DO  igrid=1,ngrids
     DO  if1=1,nfacex
        DO  if2=1,6
           fnxtf(if1,if2,igrid)=0
           eoff(if2,if1,igrid)=0
           voff(if1,if2,igrid)=0
        END DO
        neoff(if1,igrid)=0
     END DO
     DO  ie1=1,nedgex
        DO  if2=1,2
           fnxte(if2,ie1,igrid)=0
           feofe(ie1,if2,igrid)=0
           vofe(if2,ie1,igrid)=0
        END DO
     END DO
     DO  iv1=1,nvertx
        DO  if2=1,3
           fofv(iv1,if2,igrid)=0
           eofv(if2,iv1,igrid)=0
        END DO
        neofv(iv1,igrid)=0
     END DO
  END DO

  !     ----------------------------------------------------------------

  !     Cross-reference tables for dodecahedron

  igrid=1

  OPEN(80,FILE='dodecahedron.xref')

  nface(1)=12
  nedge(1)=30
  nvert(1)=20

  !     Faces next to each face
  READ(80,*)
  READ(80,*) ((fnxtf(if1,if2,1),if2=1,5),if1=1,12)

  !     Edges of each face
  READ(80,*)
  READ(80,*) ((eoff(ie1,if1,1),if1=1,12),ie1=1,5)

  !     Faces next to each edge
  READ(80,*)
  READ(80,*) ((fnxte(if1,ie1,1),if1=1,2),ie1=1,30)

  !     Faces at the ends of each edge
  READ(80,*)
  READ(80,*) ((feofe(ie1,if1,1),if1=1,2),ie1=1,30)

  !     Vertices of each edge
  READ(80,*)
  READ(80,*) ((vofe(iv1,ie1,1),iv1=1,2),ie1=1,30)

  !     Vertices of each face
  READ(80,*)
  READ(80,*) ((voff(if1,iv1,1),iv1=1,5),if1=1,12)

  !     Faces around each vertex
  READ(80,*)
  READ(80,*) ((fofv(iv1,if1,1),if1=1,3),iv1=1,20)

  !     Construct table for edges around each vertex
  DO  ie1=1,nedge(1)
     iv1=vofe(1,ie1,1)
     iv2=vofe(2,ie1,1)
     CALL addtab2(eofv(1,1,1),iv1,ie1,3,nvertx)
     CALL addtab2(eofv(1,1,1),iv2,ie1,3,nvertx)
  END DO

  !     -----------------------------------------------------------------

  !     Set up coordinates of face centres
  pi=4.0D0*ATAN(1.0D0)
  zz=2.0D0*COS(0.3D0*pi)
  phi0=2.0D0*ASIN(1.0D0/zz)-0.5D0*pi

  flong(1,1)=0.0D0
  flat(1,1)=0.5D0*pi
  DO  if1=1,5
     flong(1+if1,1)=(if1-1)*0.4D0*pi+0.2D0*pi
     flat(1+if1,1)=phi0
     flong(6+if1,1)=(if1-1)*0.4D0*pi
     flat(6+if1,1)=-phi0
  END DO
  flong(12,1)=0.0D0
  flat(12,1)=-0.5D0*pi

  !     ---------------------------------------------------------------

  !     Areas of faces
  DO  if1=1,12
     farea(if1,1)=pi/3.0D0
  END DO

  !     ---------------------------------------------------------------

  !     Coordinates of vertices
  DO  iv1=1,20
     !       First find cartesian coords of surrounding faces
     DO  if1=1,3
        if2=fofv(iv1,if1,1)
        long=flong(if2,1)
        lat=flat(if2,1)
        CALL ll2xyz(long,lat,x,y,z)
        px(if1)=x
        py(if1)=y
        pz(if1)=z
     END DO
     !       Two sides of the triangle joining the face centres
     d1x=px(2)-px(1)
     d1y=py(2)-py(1)
     d1z=pz(2)-pz(1)
     d2x=px(3)-px(1)
     d2y=py(3)-py(1)
     d2z=pz(3)-pz(1)
     !       Find the middle (it's an equilateral triangle)
     vx=(px(1)+px(2)+px(3))/3.0D0
     vy=(py(1)+py(2)+py(3))/3.0D0
     vz=(pz(1)+pz(2)+pz(3))/3.0D0
     !       Project back onto the sphere
     mag=SQRT(vx*vx+vy*vy+vz*vz)
     vx=vx/mag
     vy=vy/mag
     vz=vz/mag
     !       Convert back to latitude/longitude
     CALL xyz2ll(vx,vy,vz,long,lat)
     vlong(iv1,1)=long
     vlat(iv1,1)=lat
  END DO

  !     Tabulate lengths of edges and distances between face centres
  emn=1.0D0
  emx=0.0D0
  dmn=2.0D0
  dmx=0.0D0
  DO  ie0=1,nedge(igrid)
     !       Vertices at ends of this edge
     iv1=vofe(1,ie0,igrid)
     iv2=vofe(2,ie0,igrid)
     long=vlong(iv1,igrid)
     lat=vlat(iv1,igrid)
     CALL ll2xyz(long,lat,x1,y1,z1)
     long=vlong(iv2,igrid)
     lat=vlat(iv2,igrid)
     CALL ll2xyz(long,lat,x2,y2,z2)
     CALL spdist(x1,y1,z1,x2,y2,z2,s)
     ldist(ie0,igrid)=s
     emn=MIN(emn,ldist(ie0,igrid))
     emx=MAX(emx,ldist(ie0,igrid))
     !       Faces either side of this edge
     if1=fnxte(1,ie0,igrid)
     if2=fnxte(2,ie0,igrid)
     long=flong(if1,igrid)
     lat=flat(if1,igrid)
     CALL ll2xyz(long,lat,x1,y1,z1)
     long=flong(if2,igrid)
     lat=flat(if2,igrid)
     CALL ll2xyz(long,lat,x2,y2,z2)
     CALL spdist(x1,y1,z1,x2,y2,z2,s)
     ddist(ie0,igrid)=s
     dmn=MIN(dmn,ddist(ie0,igrid))
     dmx=MAX(dmx,ddist(ie0,igrid))
  END DO

  !     ------------------------------------------------------------------

  !     Generate the next grid in the series
  !     ====================================

  IF (ngrids == 1) GO TO 699

700 CONTINUE

  igridm=igrid
  igrid=igrid+1
  nface(igrid)=nface(igridm)+nedge(igridm)
  nedge(igrid)=4*nedge(igridm)
  nvert(igrid)=4*nvert(igridm)

  !     Old faces keep their coordinates
  DO  if1=1,nface(igridm)
     flong(if1,igrid)=flong(if1,igridm)
     flat(if1,igrid)=flat(if1,igridm)
  END DO

  !     Loop over old edges
  DO  ie0=1,nedge(igridm)

     !       Each old edge generates a new face
     ifn=nface(igridm)+ie0

     !       This is next to the two old faces either side of the old edge
     if1=fnxte(1,ie0,igridm)
     if2=fnxte(2,ie0,igridm)
     CALL addtab(fnxtf(1,1,igrid),if1,ifn,nfacex,6)
     CALL addtab(fnxtf(1,1,igrid),if2,ifn,nfacex,6)
     CALL addtab(fnxtf(1,1,igrid),ifn,if1,nfacex,6)
     CALL addtab(fnxtf(1,1,igrid),ifn,if2,nfacex,6)

     !       There is a new edge between IF1 and IFN and similarly
     !       between IF2 and IFN
     ie1=2*ie0-1
     ie2=2*ie0
     CALL addtab2(fnxte(1,1,igrid),ie1,if1,2,nedgex)
     CALL addtab2(fnxte(1,1,igrid),ie1,ifn,2,nedgex)
     CALL addtab2(eoff(1,1,igrid),if1,ie1,6,nfacex)
     CALL addtab2(eoff(1,1,igrid),ifn,ie1,6,nfacex)
     CALL addtab2(fnxte(1,1,igrid),ie2,if2,2,nedgex)
     CALL addtab2(fnxte(1,1,igrid),ie2,ifn,2,nedgex)
     CALL addtab2(eoff(1,1,igrid),if2,ie2,6,nfacex)
     CALL addtab2(eoff(1,1,igrid),ifn,ie2,6,nfacex)

     !       Coordinates of the new face
     !         Midpoint of ends of edge
     !          IV1=VOFE(IE0,1,IGRIDM)
     !          IV2=VOFE(IE0,2,IGRIDM)
     !          LONG=VLONG(IV1,IGRID)
     !          LAT=VLAT(IV1,IGRID)
     !          CALL LL2XYZ(LONG,LAT,X1,Y1,Z1)
     !          LONG=VLONG(IV2,IGRID)
     !          LAT=VLAT(IV2,IGRID)
     !          CALL LL2XYZ(LONG,LAT,X,Y,Z)
     !       Midpoint of old faces
     long=flong(if1,igrid)
     lat=flat(if1,igrid)
     CALL ll2xyz(long,lat,x1,y1,z1)
     long=flong(if2,igrid)
     lat=flat(if2,igrid)
     CALL ll2xyz(long,lat,x,y,z)
     !       Find the middle
     vx=0.5*(x1+x)
     vy=0.5*(y1+y)
     vz=0.5*(z1+z)
     !       Project back onto the sphere
     mag=SQRT(vx*vx+vy*vy+vz*vz)
     vx=vx/mag
     vy=vy/mag
     vz=vz/mag
     !       Convert back to latitude/longitude
     CALL xyz2ll(vx,vy,vz,long,lat)
     flong(ifn,igrid)=long
     flat(ifn,igrid)=lat

  END DO

  !     Loop over old vertices
  DO  iv0=1,nvert(igridm)

     !       Each vertex has three new faces
     ie1=eofv(1,iv0,igridm)
     ie2=eofv(2,iv0,igridm)
     ie3=eofv(3,iv0,igridm)
     ifn1=nface(igridm)+ie1
     ifn2=nface(igridm)+ie2
     ifn3=nface(igridm)+ie3
     CALL addtab(fofv(1,1,igrid),iv0,ifn1,nvertx,3)
     CALL addtab(fofv(1,1,igrid),iv0,ifn2,nvertx,3)
     CALL addtab(fofv(1,1,igrid),iv0,ifn3,nvertx,3)
     CALL addtab(voff(1,1,igrid),ifn1,iv0,nfacex,6)
     CALL addtab(voff(1,1,igrid),ifn2,iv0,nfacex,6)
     CALL addtab(voff(1,1,igrid),ifn3,iv0,nfacex,6)

     !       These faces are mutual neighbours
     CALL addtab(fnxtf(1,1,igrid),ifn1,ifn2,nfacex,6)
     CALL addtab(fnxtf(1,1,igrid),ifn1,ifn3,nfacex,6)
     CALL addtab(fnxtf(1,1,igrid),ifn2,ifn1,nfacex,6)
     CALL addtab(fnxtf(1,1,igrid),ifn2,ifn3,nfacex,6)
     CALL addtab(fnxtf(1,1,igrid),ifn3,ifn1,nfacex,6)
     CALL addtab(fnxtf(1,1,igrid),ifn3,ifn2,nfacex,6)

     !       Note old faces next to new faces and corresponding edges
     if11=fnxtf(ifn1,1,igrid)
     ie11=eoff(1,ifn1,igrid)
     if12=fnxtf(ifn1,2,igrid)
     ie12=eoff(2,ifn1,igrid)
     if21=fnxtf(ifn2,1,igrid)
     ie21=eoff(1,ifn2,igrid)
     if22=fnxtf(ifn2,2,igrid)
     ie22=eoff(2,ifn2,igrid)
     if31=fnxtf(ifn3,1,igrid)
     ie31=eoff(1,ifn3,igrid)
     if32=fnxtf(ifn3,2,igrid)
     ie32=eoff(2,ifn3,igrid)

     !       Each old vertex generates three new edges and three
     !       new vertices
     DO  ix=1,3

        if0=fofv(iv0,ix,igridm)
        ien=2*nedge(igridm)+3*(iv0-1)+ix
        ivn=nvert(igridm)+3*(iv0-1)+ix
        CALL addtab2(eofv(1,1,igrid),iv0,ien,3,nvertx)
        CALL addtab2(eofv(1,1,igrid),ivn,ien,3,nvertx)
        CALL addtab2(vofe(1,1,igrid),ien,iv0,2,nedgex)
        CALL addtab2(vofe(1,1,igrid),ien,ivn,2,nedgex)

        !         New vertex is a vertex of the old face
        CALL addtab(voff(1,1,igrid),if0,ivn,nfacex,6)
        CALL addtab(fofv(1,1,igrid),ivn,if0,nvertx,3)

        !         If the old face is a neighbour of a new face then
        !         the new vertex and edge also pertain to the new face
        IF ((if0 == if11).OR.(if0 == if12)) THEN
           CALL addtab(voff(1,1,igrid),ifn1,ivn,nfacex,6)
           CALL addtab(fofv(1,1,igrid),ivn,ifn1,nvertx,3)
           CALL addtab2(eoff(1,1,igrid),ifn1,ien,6,nfacex)
           CALL addtab2(fnxte(1,1,igrid),ien,ifn1,2,nedgex)
        END IF
        IF ((if0 == if21).OR.(if0 == if22)) THEN
           CALL addtab(voff(1,1,igrid),ifn2,ivn,nfacex,6)
           CALL addtab(fofv(1,1,igrid),ivn,ifn2,nvertx,3)
           CALL addtab2(eoff(1,1,igrid),ifn2,ien,6,nfacex)
           CALL addtab2(fnxte(1,1,igrid),ien,ifn2,2,nedgex)
        END IF
        IF ((if0 == if31).OR.(if0 == if32)) THEN
           CALL addtab(voff(1,1,igrid),ifn3,ivn,nfacex,6)
           CALL addtab(fofv(1,1,igrid),ivn,ifn3,nvertx,3)
           CALL addtab2(eoff(1,1,igrid),ifn3,ien,6,nfacex)
           CALL addtab2(fnxte(1,1,igrid),ien,ifn3,2,nedgex)
        END IF

        !         If the old face is a neighbour of a new face then
        !         the edge between the old and new faces pertains
        !         to the new vertex
        IF (if0 == if11) THEN
           CALL addtab2(eofv(1,1,igrid),ivn,ie11,3,nvertx)
           CALL addtab2(vofe(1,1,igrid),ie11,ivn,2,nedgex)
        END IF
        IF (if0 == if12) THEN
           CALL addtab2(eofv(1,1,igrid),ivn,ie12,3,nvertx)
           CALL addtab2(vofe(1,1,igrid),ie12,ivn,2,nedgex)
        END IF
        IF (if0 == if21) THEN
           CALL addtab2(eofv(1,1,igrid),ivn,ie21,3,nvertx)
           CALL addtab2(vofe(1,1,igrid),ie21,ivn,2,nedgex)
        END IF
        IF (if0 == if22) THEN
           CALL addtab2(eofv(1,1,igrid),ivn,ie22,3,nvertx)
           CALL addtab2(vofe(1,1,igrid),ie22,ivn,2,nedgex)
        END IF
        IF (if0 == if31) THEN
           CALL addtab2(eofv(1,1,igrid),ivn,ie31,3,nvertx)
           CALL addtab2(vofe(1,1,igrid),ie31,ivn,2,nedgex)
        END IF
        IF (if0 == if32) THEN
           CALL addtab2(eofv(1,1,igrid),ivn,ie32,3,nvertx)
           CALL addtab2(vofe(1,1,igrid),ie32,ivn,2,nedgex)
        END IF

     END DO

  END DO

  !     Calculate coordinates of vertices - Voronoi grid
  DO  iv1=1,nvert(igrid)
     CALL vrtcrd(iv1,igrid)
  END DO

  !     Now tweak face and vertex positions a la Heikes and Randall
  !     to minimise cost function related to inaccuracy in Laplacian
  !     Outer loop over iterations
  !      FTOL=1.0D-5
  ftol=1.0D-2
  totcst=0.0D0
  maxdd=0.0D0
  DO ie1=1,nedge(igrid)
     CALL penalt(ie1,igrid,cost,dd)
     totcst=totcst+cost
     IF (dd > maxdd) THEN
        maxdd=dd
        emax=ie1
     END IF
  END DO
  PRINT *,'Init. Heikes-Randall cost func ',totcst
  PRINT *,'MAXDD=                         ',maxdd,' Edge ',emax
  totcst=0.0D0
  maxdd=0.0D0
  DO if1=1,nface(igrid)
     CALL penalc(if1,igrid,cost,dd)
     totcst=totcst+cost
     IF (dd > maxdd) THEN
        maxdd=dd
        fmax=if1
     END IF
  END DO
  PRINT *,'Init. centroidal cost func     ',totcst
  PRINT *,'MAXDD=                         ',maxdd,' Face ',fmax

  !     Weights for HR (weight1) and centroidal (weight2) contributions
  !     to cost function
  weight1=1.0D0
  weight2=0.0D0
  !     Initial interval size for BRENT
  bdisp=2.0D0**(-igrid)

  !DO iter=1,40
  DO iter=1,2000 !Peixoto 
     lprint = .true.
     !       Loop over faces of grid
     mxdisp=0.0D0
     DO if2=13,nface(igrid)
        !         Minimize cost associated with face IF1
        CALL powell(if2,igrid,ftol,bdisp,disp,ngrids)
        mxdisp=MAX(mxdisp,disp)
     END DO
     IF (lprint) THEN
        !PRINT *,'Done iteration ',iter
        totcst=0.0D0
        maxdd=0.0D0
        DO ie1=1,nedge(igrid)
           CALL penalt(ie1,igrid,cost,dd)
           totcst=totcst+cost
           IF (dd > maxdd) THEN
              maxdd=dd
              emax=ie1
           END IF
        END DO
        print '(a6,i4,a6,i8,a12,d20.12,a12,d22.12)', 'GRID: ', &
             IGRID-1, ' ITER: ', ITER, &
             ' HR Cost: ', totcst, " MAXDD:", MAXDD
        !    PRINT *,'Heikes-Randall cost function ',totcst
        !    PRINT *,'MAXDD=                       ',maxdd, ' Edge ',emax


        !     Stop criteria  - P. Peixoto
        TMP=1.0D-4/(2**REAL(NGRIDS))
        IF (TOTCST.LT.TMP ) THEN
           print*, "Stopping here due to energy too small"
           print*, "Limit:", TMP
           exit
        ENDIF
        tmp=TOTCSTOLD-TOTCST
        IF (abs(tmp) < 1.0D-11 .and. iter>10) THEN
           print*, "Stopping here due to small changes in energy"
           print*, "Old: ", TOTCSTOLD
           print*, "Now: ", TOTCST
           exit
        ENDIF

        TOTCSTOLD=TOTCST
        totcst=0.0D0
        maxdd=0.0D0
        DO if1=1,nface(igrid)
           CALL penalc(if1,igrid,cost,dd)
           totcst=totcst+cost
           IF (dd > maxdd) THEN
              maxdd=dd
              fmax=if1
           END IF
        END DO

	!Print grids in intermediate steps for very fine grids
        IF ((IGRID>7.and. mod(iter,20)==0)) THEN
           !     Create file of grid coordinates only - x,y,z - P. Peixoto
           WRITE(ATMP,'(i3.3)') IGRID-1
           YNAME="grd/HR95JT_"//trim(ATMP)//"_evol.xyz"
           !WRITE(ATMP,'(i6.6)') iter
           !YNAME=trim(YNAME)//"_"//ATMP//".xyz"
           print*, "Grid file with nodes saved as:", yname
           OPEN(44,FILE=YNAME,FORM='FORMATTED')

           WRITE(44,'(i10)') NFACE(IGRID)
           DO if1 = 1, NFACE(IGRID)
              long = flong(if1,IGRID)
              lat = flat(if1,IGRID)
              CALL ll2xyz(long,lat,x1,y1,z1)
              WRITE(44,'(3e32.16)') x1,y1,z1
           END DO
           CLOSE(44)
        ENDIF

        !    PRINT *,'Centroidal cost function     ',totcst
        !    PRINT *,'MAXDD=                       ',maxdd, ' Face ',fmax
        !    PRINT *,'Max disp in BRENT = ',mxdisp
     END IF

     !       BRENT interval size for next iteration is max displacement
     !       from this iteration
     bdisp=mxdisp

  END DO

  !     Create file of grid coordinates only - x,y,z - P. Peixoto
  WRITE(ATMP,'(i3.3)') IGRID-1
  YNAME="grd/HR95JT_"//trim(ATMP)//".xyz"
  print*, "Grid file with nodes saved as:", trim(yname)
  OPEN(44,FILE=YNAME,FORM='FORMATTED')
  
  WRITE(44,'(i10)') NFACE(IGRID)
  DO if1 = 1, NFACE(IGRID)
     long = flong(if1,IGRID)
     lat = flat(if1,IGRID)
     CALL ll2xyz(long,lat,x1,y1,z1)
     WRITE(44,'(3e32.16)') x1,y1,z1
  END DO
  CLOSE(44)

  !     Tabulate lengths of edges and distances between face centres
  !     across each edge
  emn=5.0D0
  emx=0.0D0
  dmn=5.0D0
  dmx=0.0D0
  dav=0.0D0
  DO  ie0=1,nedge(igrid)
     !       Vertices at ends of this edge
     iv1=vofe(1,ie0,igrid)
     iv2=vofe(2,ie0,igrid)
     long=vlong(iv1,igrid)
     lat=vlat(iv1,igrid)
     CALL ll2xyz(long,lat,x1,y1,z1)
     long=vlong(iv2,igrid)
     lat=vlat(iv2,igrid)
     CALL ll2xyz(long,lat,x2,y2,z2)
     CALL spdist(x1,y1,z1,x2,y2,z2,s)
     ldist(ie0,igrid)=s
     emn=MIN(emn,ldist(ie0,igrid))
     emx=MAX(emx,ldist(ie0,igrid))
     !       Faces either side of this edge
     if1=fnxte(1,ie0,igrid)
     if2=fnxte(2,ie0,igrid)
     long=flong(if1,igrid)
     lat=flat(if1,igrid)
     CALL ll2xyz(long,lat,x1,y1,z1)
     long=flong(if2,igrid)
     lat=flat(if2,igrid)
     CALL ll2xyz(long,lat,x2,y2,z2)
     CALL spdist(x1,y1,z1,x2,y2,z2,s)
     ddist(ie0,igrid)=s
     dmn=MIN(dmn,ddist(ie0,igrid))
     dmx=MAX(dmx,ddist(ie0,igrid))
     dav=dav+ddist(ie0,igrid)/nedge(igrid)
  END DO

  PRINT *,' Done grid ',igrid-1

  PRINT *,'min side: ',emn,' max side: ',emx,' ratio: ',emn/emx
  PRINT *,'min dist: ',dmn,' max dist: ',dmx,' ratio: ',dmn/dmx
  PRINT *,'average side * rearth: ',dav*6371220.0D0



  IF (igrid < ngrids) GO TO 700
699 CONTINUE

  !     ------------------------------------------------------------------

  !     Calculate areas of grid cells on each grid

  DO  igrid=1,ngrids

     amn=4.0D0*pi
     amx=0.0D0
     atot=0.0D0
     DO  if1=1,nface(igrid)
        aface=0.0D0
        !       Coordinates of centre of face
        long=flong(if1,igrid)
        lat=flat(if1,igrid)
        CALL ll2xyz(long,lat,x0,y0,z0)
        !       Loop over edges in turn and calculate area of triangle
        !       formed by the edge and the centre of the face
        IF (if1 <= 12) THEN
           NE=5
        ELSE
           NE=6
        END IF
        DO  ie1=1,NE
           ie2=eoff(ie1,if1,igrid)
           iv1=vofe(1,ie2,igrid)
           iv2=vofe(2,ie2,igrid)
           long=vlong(iv1,igrid)
           lat=vlat(iv1,igrid)
           CALL ll2xyz(long,lat,x1,y1,z1)
           long=vlong(iv2,igrid)
           lat=vlat(iv2,igrid)
           CALL ll2xyz(long,lat,x2,y2,z2)
           CALL starea2(x0,y0,z0,x1,y1,z1,x2,y2,z2,da)
           aface=aface+da
        END DO
        farea(if1,igrid)=aface
        atot=atot+aface
     END DO

     PRINT *,'Grid ',igrid,' total area=',atot
     !     Now scale up areas so total is 4*pi
     !     (Not really necessary now that we use the correct spherical
     !     triangle area)
     sfac=4.0D0*pi/atot
     PRINT *,'fraction of sphere: ',1.0D0/sfac
     DO  if1=1,nface(igrid)
        farea(if1,igrid)=farea(if1,igrid)*sfac
        amn=MIN(amn,farea(if1,igrid))
        amx=MAX(amx,farea(if1,igrid))
     END DO

     PRINT *,'min area: ',amn,' max area: ',amx,' ratio: ',amn/amx
     PRINT *,'average area =',5.101D14/nface(igrid)

  END DO

  !     ------------------------------------------------------------------

  !     Every vertex has three edges
  DO  igrid=1,ngrids
     DO  iv0=1,nvert(igrid)
        neofv(iv0,igrid)=3
     END DO
  END DO

  !     Sort FNXTF into anticlockwise order on each grid
  !     and sort EOFF to correspond to FNXTF
  DO  igrid=1,ngrids
     DO  if0=1,nface(igrid)
        IF (if0 <= 12) THEN
           nf=5
        ELSE
           nf=6
        END IF
        neoff(if0,igrid)=nf
        !       Coordinates of face IF0
        long=flong(if0,igrid)
        lat=flat(if0,igrid)
        CALL ll2xyz(long,lat,x,y,z)
        DO  ix1=1,nf-2
           !         Coordinates of IX1'th neighbour
           if1=fnxtf(if0,ix1,igrid)
           long=flong(if1,igrid)
           lat=flat(if1,igrid)
           CALL ll2xyz(long,lat,x1,y1,z1)
           d1x=x1-x
           d1y=y1-y
           d1z=z1-z
           ix2=ix1
920        CONTINUE
           ix2=ix2+1
           !           Coordinates of IX2'th neighbour
           if2=fnxtf(if0,ix2,igrid)
           long=flong(if2,igrid)
           lat=flat(if2,igrid)
           CALL ll2xyz(long,lat,x2,y2,z2)
           d2x=x2-x
           d2y=y2-y
           d2z=z2-z
           cs=d1x*d2x+d1y*d2y+d1z*d2z
           IF (cs < 0.0D0) GO TO 920
           sn=x*(d1y*d2z-d1z*d2y) +y*(d1z*d2x-d1x*d2z)  &
                +z*(d1x*d2y-d1y*d2x)
           IF (sn < 0.0D0) GO TO 920
           !         IF2 belongs in position IX1+1 so swap them
           if3=fnxtf(if0,ix1+1,igrid)
           fnxtf(if0,ix1+1,igrid)=if2
           fnxtf(if0,ix2,igrid)=if3
        END DO
        DO  ix1=1,nf
           if1=fnxtf(if0,ix1,igrid)
           ix2=ix1-1
940        CONTINUE
           ix2=ix2+1
           ie2=eoff(ix2,if0,igrid)
           if21=fnxte(1,ie2,igrid)
           if22=fnxte(2,ie2,igrid)
           IF ((if21+if22) /= (if0+if1)) GO TO 940
           !         Edge IE2 corresponds to face IF1
           eoff(ix2,if0,igrid)=eoff(ix1,if0,igrid)
           eoff(ix1,if0,igrid)=ie2
        END DO
     END DO
  END DO

  !     Order VOFF so that the k'th vertex is between the
  !     k'th and (k+1)'th edges in EOFF
  DO  igrid=1,ngrids
     DO  if0=1,nface(igrid)
        DO  ix1=1,neoff(if0,igrid)
           ix2=ix1+1
           IF (ix2 > neoff(if0,igrid)) ix2=1
           ie1=eoff(ix1,if0,igrid)
           ie2=eoff(ix2,if0,igrid)
           ! Find the common vertex of IE1 and IE2
           iv11=vofe(1,ie1,igrid)
           iv12=vofe(2,ie1,igrid)
           iv21=vofe(1,ie2,igrid)
           iv22=vofe(2,ie2,igrid)
           IF ((iv11 == iv21).OR.(iv11 == iv22)) THEN
              iv0=iv11
           ELSE IF ((iv12 == iv21).OR.(iv12 == iv22)) THEN
              iv0=iv12
           ELSE
              PRINT *,'Common vertex not found'
              STOP
           END IF
           voff(if0,ix1,igrid) = iv0
        END DO
     END DO
  END DO

  !     Find faces at ends of edges for output grid and
  !     sort FEOFE so that FEOFE(1) -> FEOFE(2) is 90 deg anticlockwise
  !     of FNXTE(1) -> FNXTE(2) and sort VOFE so that VOFE(1) -> VOFE(2)
  !     is 90 degrees anticlockwise of FNXTE(1) -> FNXTE(2)
  DO  igrid=1,ngrids
     DO  ie0=1,nedge(igrid)
        if1=fnxte(1,ie0,igrid)
        if2=fnxte(2,ie0,igrid)
        long=flong(if1,igrid)
        lat=flat(if1,igrid)
        CALL ll2xyz(long,lat,x,y,z)
        long=flong(if2,igrid)
        lat=flat(if2,igrid)
        CALL ll2xyz(long,lat,x1,y1,z1)
        d1x=x1-x
        d1y=y1-y
        d1z=z1-z
        lswap=.false.
        DO  ix1=1,2
           iv1=vofe(ix1,ie0,igrid)
           DO  ix2=1,3
              if3=fofv(iv1,ix2,igrid)
              IF ((if3 /= if1).AND.(if3 /= if2)) THEN
                 long=flong(if3,igrid)
                 lat=flat(if3,igrid)
                 CALL ll2xyz(long,lat,x2,y2,z2)
                 d2x=x2-x
                 d2y=y2-y
                 d2z=z2-z
                 sn=x*(d1y*d2z-d1z*d2y) +y*(d1z*d2x-d1x*d2z)  &
                      +z*(d1x*d2y-d1y*d2x)
                 IF (sn > 0.0D0) THEN
                    feofe(ie0,2,igrid)=if3
                    IF (ix1 == 1) lswap=.true.
                 ELSE
                    feofe(ie0,1,igrid)=if3
                 END IF
              END IF
           END DO
        END DO
        IF (lswap) THEN
           iv1=vofe(1,ie0,igrid)
           vofe(1,ie0,igrid)=vofe(2,ie0,igrid)
           vofe(2,ie0,igrid)=iv1
        END IF
     END DO
  END DO


  !     Find distance between faces along each edge
  igrid=ngrids
  gmn=5.0D0
  gmx=0.0D0
  DO  ie0=1,nedge(ngrids)
     !       Faces either end of this edge
     if1=feofe(ie0,1,igrid)
     if2=feofe(ie0,2,igrid)
     long=flong(if1,igrid)
     lat=flat(if1,igrid)
     CALL ll2xyz(long,lat,x1,y1,z1)
     long=flong(if2,igrid)
     lat=flat(if2,igrid)
     CALL ll2xyz(long,lat,x2,y2,z2)
     CALL spdist(x1,y1,z1,x2,y2,z2,s)
     gdist(ie0,igrid)=s
     gmn=MIN(gmn,gdist(ie0,igrid))
     gmx=MAX(gmx,gdist(ie0,igrid))
  END DO

  PRINT *,'min area: ',amn,' max area: ',amx,' ratio: ',amn/amx
  PRINT *,'min side: ',emn,' max side: ',emx,' ratio: ',emn/emx
  PRINT *,'min dist: ',dmn,' max dist: ',dmx,' ratio: ',dmn/dmx
  PRINT *,'min dist2: ',gmn,' max dist2: ',gmx,' ratio: ',gmn/gmx

  !     ------------------------------------------------------------------

  !     Calculate coefficients for interpolating stream function and
  !     velocity potential to vertices.

  DO  igrid=1,ngrids

     DO  iv1=1,nvert(igrid)
        if1=fofv(iv1,1,igrid)
        if2=fofv(iv1,2,igrid)
        if3=fofv(iv1,3,igrid)
        long=vlong(iv1,igrid)
        lat=vlat(iv1,igrid)
        CALL ll2xyz(long,lat,x,y,z)
        long=flong(if1,igrid)
        lat=flat(if1,igrid)
        CALL ll2xyz(long,lat,x1,y1,z1)
        long=flong(if2,igrid)
        lat=flat(if2,igrid)
        CALL ll2xyz(long,lat,x2,y2,z2)
        long=flong(if3,igrid)
        lat=flat(if3,igrid)
        CALL ll2xyz(long,lat,x3,y3,z3)
        x1=x1-x
        y1=y1-y
        z1=z1-z
        x2=x2-x
        y2=y2-y
        z2=z2-z
        x3=x3-x
        y3=y3-y
        z3=z3-z
        w1=x*(y2*z3-z2*y3) +y*(z2*x3-x2*z3)  &
             +z*(x2*y3-y2*x3)
        w2=x*(y3*z1-z3*y1) +y*(z3*x1-x3*z1)  &
             +z*(x3*y1-y3*x1)
        w3=x*(y1*z2-z1*y2) +y*(z1*x2-x1*z2)  &
             +z*(x1*y2-y1*x2)
        w1=ABS(w1)
        w2=ABS(w2)
        w3=ABS(w3)
        tot=w1+w2+w3
        w1=w1/tot
        w2=w2/tot
        w3=w3/tot
        coeff(iv1,1,igrid)=w1
        coeff(iv1,2,igrid)=w2
        coeff(iv1,3,igrid)=w3
     END DO
     DO  iv1=nvert(igrid)+1,nvertx
        coeff(iv1,1,igrid)=0.0D0
        coeff(iv1,2,igrid)=0.0D0
        coeff(iv1,3,igrid)=0.0D0
     END DO

  END DO

  !     ------------------------------------------------------------------

  !     Create a file of cross reference tables and coordinates etc
  !     for use by the model

  !PRINT *,'Create a GRIDMAP file (0 or 1) ?'
  !READ (5,*) if0
  !IF (if0 == 1) THEN
  IF (.TRUE.) THEN
     WRITE(ATMP,'(i10.10)') nfacex
     yname="grd/gridmap_hex_"//trim(ATMP)//".dat"
     print*, "Grid file saved as:", trim(yname)
     OPEN(22,FILE=yname,FORM='UNFORMATTED')
     
     !WRITE(yname,'(''grd/gridmap_hex_'',I10.10,''.dat'')') nfacex
     !OPEN(22,FILE=yname,FORM='UNFORMATTED')

     PRINT *,'NFACE = ',nface
     PRINT *,'NEDGE = ',nedge
     PRINT *,'NVERT = ',nvert

     !        WRITE(22,*) 'GRIDMAP for NGRIDS=',NGRIDS
     WRITE(22) ngrids
     WRITE(22) nface
     WRITE(22) nedge
     WRITE(22) nvert
     !        WRITE(22,*) 'Number of edges of each face - all grids'
     WRITE(22) ((neoff(if1,igrid), if1=1,nface(igrid)),  &
          igrid=1,ngrids)
     !        WRITE(22,*) 'Number of edges of each vertex - all grids'
     WRITE(22) ((neofv(iv1,igrid), iv1=1,nvert(igrid)),  &
          igrid=1,ngrids)
     !        WRITE(22,*) 'Faces next to each face - all grids'
     WRITE(22) (((fnxtf(if1,if2,igrid), if1=1,nface(igrid)),  &
          if2=1,6), igrid=1,ngrids)
     !        WRITE(22,*) 'Edges of each face - all grids'
     WRITE(22) (((eoff(ie2,if1,igrid), ie2=1,6),  &
          if1=1,nface(igrid)), igrid=1,ngrids)
     !        WRITE(22,*) 'Vertices of each face - all grids'
     WRITE(22) (((voff(if1,iv1,igrid), if1=1,nface(igrid)),  &
          iv1=1,6), igrid=1,ngrids)
     !        WRITE(22,*) 'Faces next to each edge - all grids'
     WRITE(22) (((fnxte(if2,ie1,igrid), if2=1,2),  &
          ie1=1,nedge(igrid)), igrid=1,ngrids)
     !        WRITE(22,*) 'Vertices of each edge - all grids'
     WRITE(22) (((vofe(iv2,ie1,igrid), iv2=1,2),  &
          ie1=1,nedge(igrid)), igrid=1,ngrids)
     !        WRITE(22,*) 'Faces around each vertex - all grids'
     WRITE(22) (((fofv(iv1,if2,igrid), iv1=1,nvert(igrid)),  &
          if2=1,3), igrid=1,ngrids)
     !        WRITE(22,*) 'Edges around each vertex - all grids'
     WRITE(22) (((eofv(ie1,iv1,igrid), ie1=1,3),  &
          iv1=1,nvert(igrid)), igrid=1,ngrids)
     !        WRITE(22,*) 'Coefficients for interpolation - all grids'
     !        WRITE(22) (((COEFF(IV1,IF2,IGRID),
     !     :                 IV1=1,NVERT(IGRID)),
     !     :                 IGRID=1,NGRIDS),
     !     :                 IF2=1,3)
     !        WRITE(22,*) 'Longitudes of faces - all grids'
     WRITE(22) ((flong(if1,igrid), if1=1,nface(igrid)),  &
          igrid=1,ngrids)
     !        WRITE(22,*) 'Latitudes of faces - all grids'
     WRITE(22) ((flat(if1,igrid), if1=1,nface(igrid)),  &
          igrid=1,ngrids)
     !        WRITE(22,*) 'Longitudes of vertices - all grids'
     WRITE(22) ((vlong(iv1,igrid), iv1=1,nvert(igrid)),  &
          igrid=1,ngrids)
     !        WRITE(22,*) 'Latitudes of vertices - all grids'
     WRITE(22) ((vlat(iv1,igrid), iv1=1,nvert(igrid)),  &
          igrid=1,ngrids)
     !        WRITE(22,*) 'Areas of faces - all grids'
     WRITE(22) ((farea(if1,igrid), if1=1,nface(igrid)),  &
          igrid=1,ngrids)
     !        WRITE(22,*) 'Lengths of edges - all grids'
     WRITE(22) ((ldist(ie1,igrid), ie1=1,nedge(igrid)),  &
          igrid=1,ngrids)
     !        WRITE(22,*) 'Distance between faces across edges - all grids'
     WRITE(22) ((ddist(ie1,igrid), ie1=1,nedge(igrid)),  &
          igrid=1,ngrids)

     ! NO SFC Ordering, yet
     WRITE(22) 0

     CLOSE(22)

  END IF

  !     Create file of grid coordinates only, for plotting
!  OPEN(44,FILE='primalgrid.dat',FORM='FORMATTED') 
!P.Peixoto edit to include grid level
  WRITE(ATMP,'(i3.3)') NGRIDS-1
  YNAME="grd/primalgrid_"//trim(ATMP)//".dat"
  print*, "Grid file saved as:", trim(yname)
  OPEN(44,FILE=YNAME,FORM='FORMATTED')
  DO  ie1 = 1, nedgex
     iv1 = vofe(1,ie1,ngrids)
     long = vlong(iv1,ngrids)
     lat = vlat(iv1,ngrids)
     CALL ll2xyz(long,lat,x1,y1,z1)
     iv1 = vofe(2,ie1,ngrids)
     long = vlong(iv1,ngrids)
     lat = vlat(iv1,ngrids)
     CALL ll2xyz(long,lat,x2,y2,z2)
     WRITE(44,'(6e15.7)') x1,y1,z1,x2,y2,z2
  END DO
  CLOSE(44)
  
  !OPEN(44,FILE='dualgrid.dat',FORM='FORMATTED') 
  !P.Peixoto edit to include grid level
  WRITE(ATMP,'(i3.3)') NGRIDS-1
  YNAME="grd/dualgrid_"//trim(ATMP)//".dat"
  print*, "Grid file saved as:", trim(yname)
  OPEN(44,FILE=YNAME,FORM='FORMATTED')
  DO  ie1 = 1, nedgex
     if1 = fnxte(1,ie1,ngrids)
     long = flong(if1,ngrids)
     lat = flat(if1,ngrids)
     CALL ll2xyz(long,lat,x1,y1,z1)
     if1 = fnxte(2,ie1,ngrids)
     long = flong(if1,ngrids)
     lat = flat(if1,ngrids)
     CALL ll2xyz(long,lat,x2,y2,z2)
     WRITE(44,'(6e15.7)') x1,y1,z1,x2,y2,z2
  END DO
  CLOSE(44)

  !     Create file of grid coordinates only - x,y,z - P. Peixoto
  !WRITE(YNAME,'(''HR95JT_'',I3.3,''.xyz'')') NGRIDS-1
  !OPEN(44,FILE=YNAME,FORM='FORMATTED')
  
  !WRITE(44,'(i10)') nfacex
  !DO  if1 = 1, nfacex 
  !   long = flong(if1,ngrids)
  !   lat = flat(if1,ngrids)
  !   CALL ll2xyz(long,lat,x1,y1,z1)
  !   WRITE(44,'(3e32.16)') x1,y1,z1
  !END DO
  !CLOSE(44)
  
  !     -------------------------------------------------------------------

  !     Check grid cross-references by writing out coords for plotting

  GO TO 503
  
  !P.Peixoto edit to include grid level
  WRITE(ATMP,'(i3.3)') NGRIDS-1
  YNAME="grd/grid.coords_"//trim(ATMP)//".dat"
  print*, "Grid file with nodes saved as:", trim(yname)
  OPEN(82,FILE=YNAME,FORM='FORMATTED')
  !OPEN(82,FILE='grid.coords')

  DO  igrid=1,ngrids

     !     VLONG(IV1,IGRID) and VLAT(IV1,IGRID) contain the latitude
     !     and longitude of the IV1'th vertex on grid IGRID

     !     1) plot edges of dodecahedron
     WRITE(82,*) nedge(igrid)
     DO  iedge=1,nedge(igrid)
        iv1=vofe(1,iedge,igrid)
        iv2=vofe(2,iedge,igrid)
        long=vlong(iv1,igrid)
        lat=vlat(iv1,igrid)
        CALL ll2xyz(long,lat,x1,y1,z1)
        long=vlong(iv2,igrid)
        lat=vlat(iv2,igrid)
        CALL ll2xyz(long,lat,x2,y2,z2)
        WRITE (82,598) x1,y1,z1,x2,y2,z2
     END DO

  END DO

  !     Plot pairs of lines crossing each edge
  igrid=ngrids
  WRITE(82,*) nedge(igrid)
  DO  iedge=1,nedge(igrid)
     !        IV1=FNXTE(IEDGE,1,IGRID)
     !        IV2=FNXTE(IEDGE,2,IGRID)
     !        LONG=FLONG(IV1,IGRID)
     !        LAT=FLAT(IV1,IGRID)
     !        CALL LL2XYZ(LONG,LAT,X1,Y1,Z1)
     !        LONG=FLONG(IV2,IGRID)
     !        LAT=FLAT(IV2,IGRID)
     !        CALL LL2XYZ(LONG,LAT,X2,Y2,Z2)
     !        WRITE (82,598) X1,Y1,Z1,X2,Y2,Z2
     iv1=feofe(iedge,1,igrid)
     iv2=feofe(iedge,2,igrid)
     long=flong(iv1,igrid)
     lat=flat(iv1,igrid)
     CALL ll2xyz(long,lat,x1,y1,z1)
     long=flong(iv2,igrid)
     lat=flat(iv2,igrid)
     CALL ll2xyz(long,lat,x2,y2,z2)
     WRITE (82,598) x1,y1,z1,x2,y2,z2
  END DO

  DO  igrid=1,ngrids

     !     FLONG(IF1,IGRID) and FLAT(IF1,IGRID) contain the latitude
     !     and longitude of the IF1'th face centre on grid IGRID

     !     2) plot edges of icosahedron
     WRITE(82,*) nedge(igrid)
     DO  if1=1,nface(igrid)
        long=flong(if1,igrid)
        lat=flat(if1,igrid)
        CALL ll2xyz(long,lat,x1,y1,z1)
        IF (if1 <= 12) THEN
           nf2=5
        ELSE
           nf2=6
        END IF
        DO  if2=1,nf2
           if3=fnxtf(if1,if2,igrid)
           IF (if3 > if1) THEN
              long=flong(if3,igrid)
              lat=flat(if3,igrid)
              CALL ll2xyz(long,lat,x2,y2,z2)
              WRITE(82,598) x1,y1,z1,x2,y2,z2
           END IF
        END DO
     END DO
  END DO

598 FORMAT(6F9.4)

  CLOSE(82)
503 CONTINUE

  !     ------------------------------------------------------------------
  print*
  STOP
END PROGRAM gengrid

!     ==================================================================

SUBROUTINE ll2xyz(long,lat,x,y,z)

  !     To convert longitude and latitude to cartesian coordinates
  !     on the unit sphere

  IMPLICIT NONE

  REAL*8, INTENT(IN OUT)                   :: long
  REAL*8, INTENT(IN OUT)                   :: lat
  REAL*8, INTENT(OUT)                      :: x
  REAL*8, INTENT(OUT)                      :: y
  REAL*8, INTENT(OUT)                      :: z

  REAL*8  cln,sln,clt,slt

  !     ------------------------------------------------------------------

  sln=SIN(long)
  cln=COS(long)
  slt=SIN(lat)
  clt=COS(lat)

  x=cln*clt
  y=sln*clt
  z=slt

  !     ------------------------------------------------------------------

  RETURN
END SUBROUTINE ll2xyz

!     ==================================================================

SUBROUTINE xyz2ll(x,y,z,long,lat)

  !     To convert cartesian coordinates to longitude and latitude

  IMPLICIT NONE

  REAL*8, INTENT(IN)                       :: x
  REAL*8, INTENT(IN)                       :: y
  REAL*8, INTENT(IN)                       :: z
  REAL*8, INTENT(OUT)                      :: long
  REAL*8, INTENT(OUT)                      :: lat

  REAL*8  pi,tln,tlt,r

  !     -------------------------------------------------------------------

  pi=4.0D0*ATAN(1.0D0)

  IF (x == 0.0D0) THEN
     IF (y >= 0.0D0) THEN
        long=0.5D0*pi
     ELSE
        long=1.5D0*pi
     END IF
  ELSE
     tln=y/x
     long=ATAN(tln)
     IF (x < 0.0D0) THEN
        long=long+pi
     END IF
     IF (long < 0.0D0) THEN
        long=long+2.0D0*pi
     END IF
  END IF

  r=SQRT(x*x+y*y)
  IF (r == 0.0D0) THEN
     IF (z > 0.0D0) THEN
        lat=0.5D0*pi
     ELSE
        lat=-0.5D0*pi
     END IF
  ELSE
     tlt=z/r
     lat=ATAN(tlt)
  END IF

  !     --------------------------------------------------------------------

  RETURN
END SUBROUTINE xyz2ll

!     ====================================================================

SUBROUTINE addtab(tab,INDEX,ENTRY,dim1,dim2)

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
  IF (i > dim2) THEN
     PRINT *,'**********'
     PRINT *,'TABLE FULL'
     PRINT *,'**********'
     STOP
  END IF
  IF (tab(INDEX,i) /= 0) GO TO 100
  tab(INDEX,i)=ENTRY

  !     ---------------------------------------------------------------------

  RETURN
END SUBROUTINE addtab

!     =====================================================================

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

!     =====================================================================

SUBROUTINE brent(ax,bx,cx,tol,n,p,xit,xmin,fmin,if1,igrid,ngchck)

  USE grid

  !     Perform line minimization of the function F using Brent's method.
  !     Based on Numerical Recipes
  !     AX,BX,CX is assumed to bracket the minimum, i.e. F(BX)<F(AX),F(CX)

  IMPLICIT NONE

  REAL*8, INTENT(IN)                       :: ax
  REAL*8, INTENT(IN)                       :: bx
  REAL*8, INTENT(IN)                       :: cx
  REAL*8, INTENT(IN)                       :: tol
  INTEGER, INTENT(IN)                      :: n
  REAL*8, INTENT(IN OUT)                   :: p(n)
  REAL*8, INTENT(IN OUT)                   :: xit(n)
  REAL*8, INTENT(OUT)                      :: xmin
  REAL*8, INTENT(OUT)                      :: fmin
  INTEGER, INTENT(IN OUT)                  :: if1
  INTEGER, INTENT(IN OUT)                  :: igrid
  INTEGER, INTENT(IN)                      :: ngchck


  INTEGER :: iter,itmax, j
  REAL*8 cgold, a,b,x,u,v,w,fx,fu,fv,fw,tol1,tol2,  e,d,zeps,pp,q,r,etemp,xm

  !     Check resolution
  IF (ngrids /= ngchck) THEN
     PRINT *,'ERROR'
     PRINT *,'NGRIDS=',ngrids,' in routine BRENT but'
     PRINT *,'NGRIDS=',ngchck,' in the calling routine.'
     STOP
  END IF

  !     Maximum number of iterations
  !itmax=100
  itmax=200  !Peixoto modif

  !     Golden ratio
  cgold=0.3819660D0

  !     Protect against divide by zero
  zeps=1.0D-10

  !     Sort A, B, into ascending order
  IF (ax < bx) THEN
     a=ax
     b=cx
  ELSE
     a=cx
     b=ax
  END IF

  !     Initialize search points and function values
  x=bx
  w=bx
  v=bx
  flong(if1,igrid)=p(1)+x*xit(1)
  flat(if1,igrid)=p(2)+x*xit(2)
  CALL pen(if1,igrid,fx)
  fw=fx
  fv=fx

  e=0.0D0
  d=0.0D0

  !     MAIN LOOP
  iter=0
100 CONTINUE

  iter=iter+1
  xm=0.5D0*(a+b)

  !     Check for convergence
  !      TOL1=TOL*ABS(X)+ZEPS
  !     In the case where the min happens to fall at x=0 but P is
  !     non-zero, better to put a typical P value (1.0) instead of x
  tol1=tol+zeps
  tol2=2.0D0*tol1
  IF (ABS(x-xm) <= tol2-0.5D0*(b-a)) GO TO 900

  !     Construct a trial parabolic fit
  IF (ABS(e) > tol1) THEN
     r=(x-w)*(fx-fv)
     q=(x-v)*(fx-fw)
     pp=(x-v)*q-(x-w)*r
     q=2.0D0*(q-r)
     IF (q > 0.0D0) pp=-pp
     q=ABS(q)
     etemp=e
     e=d
     IF (ABS(pp) >= ABS(0.5D0*q*etemp).OR. pp <= q*(a-x).OR.  &
          pp >= q*(b-x)) THEN
        IF (x >= xm) THEN
           e=a-x
        ELSE
           e=b-x
        END IF
        d=cgold*e
     ELSE
        d=pp/q
        u=x+d
        IF (u-a < tol2.OR.b-u < tol2) d=SIGN(tol1,xm-x)
     END IF
  ELSE
     IF (x >= xm) THEN
        e=a-x
     ELSE
        e=b-x
     END IF
     d=cgold*e
  END IF

  IF (ABS(d) >= tol1) THEN
     u=x+d
  ELSE
     u=x+SIGN(tol1,d)
  END IF
  flong(if1,igrid)=p(1)+u*xit(1)
  flat(if1,igrid)=p(2)+u*xit(2)
  CALL pen(if1,igrid,fu)

  IF (fu <= fx) THEN
     IF (u >= x) THEN
        a=x
     ELSE
        b=x
     END IF
     v=w
     w=x
     x=u
     fv=fw
     fw=fx
     fx=fu
  ELSE
     IF (u < x) THEN
        a=u
     ELSE
        b=u
     END IF
     IF (fu <= fw.OR.w == x) THEN
        v=w
        w=u
        fv=fw
        fw=fu
     ELSE IF (fu <= fv.OR.v == x.OR.v == w) THEN
        v=u
        fv=fu
     END IF
  END IF

  IF (iter < itmax) GO TO 100

  !      PRINT *,'Maximum iterations exceeded in BRENT'
  !      STOP

900 CONTINUE
  xmin=x
  fmin=fx

  !     Save vector displacement and new position of min
  DO j=1,n
     xit(j)=xmin*xit(j)
     p(j)=p(j)+xit(j)
  END DO

  !      print *,'Brent: ',iter,' iterations '

  RETURN
END SUBROUTINE brent

!     =================================================================

SUBROUTINE powell(if1,igrid,ftol,bdisp,mxdisp,ngchck)

  USE grid

  !     Minimize a function in 2 dimensions using Powell's method.
  !     Following Numerical Recipes

  IMPLICIT NONE

  INTEGER, INTENT(IN OUT)                  :: if1
  INTEGER, INTENT(IN OUT)                  :: igrid
  REAL*8, INTENT(IN)                       :: ftol
  REAL*8, INTENT(IN)                       :: bdisp
  REAL*8, INTENT(OUT)                      :: mxdisp
  INTEGER, INTENT(IN)                      :: ngchck


  INTEGER, PARAMETER :: n=2


  INTEGER :: j,iter,itmax,ibig,i, iv1
  REAL*8 p(2),xi(2,2),pt(2),ptt(2),xit(2), fp,fptt,ax,bx,cx,tol,  &
       del,fret,t,fac,cgold,xmin


  !     Check resolution
  IF (ngrids /= ngchck) THEN
     PRINT *,'ERROR'
     PRINT *,'NGRIDS=',ngrids,' in routine POWELL but'
     PRINT *,'NGRIDS=',ngchck,' in the calling routine.'
     STOP
  END IF

  !     Golden ratio
  cgold=0.3819660D0

  mxdisp=0.0D0

  !     Brackets for 1D search
  !      FAC=2.0D0**(-IGRID)
  !      DEL=1.18D0*FAC
  !c      DEL=0.01D0*FAC
  !      AX=-DEL
  !      BX=0.0
  !C      BX=(2.0D0*CGOLD-1.0D0)*DEL
  !      CX=DEL
  ax=-bdisp
  bx=0.0D0
  cx=bdisp

  !     Tolerance for 1d search
  !      TOL=1.0D-8
  !tol=0.01D0*bdisp
  tol=0.001D0*bdisp  !Peixoto modif

  p(1)=flong(if1,igrid)
  p(2)=flat(if1,igrid)
  CALL pen(if1,igrid,fret) !Fret is total penalty cost

  !     Initialize search directions
  DO i=1,n
     DO j=1,n
        xi(j,i)=0.0D0
     END DO
  END DO
  xi(1,1)=1.0D0/COS(flat(if1,igrid))
  xi(2,2)=1.0D0

  !     Save initial point
  DO j=1,n
     pt(j)=p(j)
  END DO

  !itmax=6
  itmax=60 !Peixoto modif

  !     MAIN LOOP
  iter=0
100 CONTINUE
  iter=iter+1
  fp=fret
  ibig=0
  del=0.0D0

  !     Loop over all directions
  DO i=1,n

     !       Copy the direction
     DO j=1,n
        xit(j)=xi(j,i)
     END DO
     fptt=fret

     !       Line minimization along direction XIT from P
     CALL brent(ax,bx,cx,tol,n,p,xit,xmin,fret,if1,igrid,ngrids)
     mxdisp=MAX(mxdisp,xmin)

     IF (ABS(fptt-fret) > del) THEN
        del=ABS(fptt-fret)
        ibig=i
     END IF

  END DO

  !     Are we finished?
  IF (2.0D0*ABS(fp-fret) <= ftol*(ABS(fp)+ABS(fret))) GO TO 900

  !     Construct extrapolated point and average direction moved
  DO j=1,n
     ptt(j)=2.0D0*p(j)-pt(j)
     xit(j)=p(j)-pt(j)
     pt(j)=p(j)
  END DO

  !     Function value at extrapolated point
  flong(if1,igrid)=ptt(1)
  flat(if1,igrid)=ptt(2)
  CALL pen(if1,igrid,fptt)

  IF (fptt < fp) THEN
     t=2.0D0*(fp-2.0D0*fret+fptt)*SQRT(fp-fret-del)-del*SQRT(fp-fptt)
     IF (t < 0.0D0) THEN
        CALL brent(ax,bx,cx,tol,n,p,xit,xmin,fret,if1,igrid,ngrids)
        mxdisp=MAX(mxdisp,xmin)
        DO j=1,n
           xi(j,ibig)=xi(j,n)
           xi(j,n)=xit(j)
        END DO
     END IF
  END IF


  IF (iter < itmax) GO TO 100


900 CONTINUE

  !     Reset long and lat to best value
  flong(if1,igrid)=p(1)
  flat(if1,igrid)=p(2)

  !     And reset positions of affected vertices
  DO i=1,6
     iv1=voff(if1,i,igrid)
     CALL vrtcrd(iv1,igrid)
  END DO

  CALL pen(if1,igrid,fret)


  RETURN
END SUBROUTINE powell

!     =================================================================

SUBROUTINE vrtcrd(iv1,igrid)

  USE grid

  !     Find the coordinates of a vertex

  IMPLICIT NONE

  INTEGER, INTENT(IN OUT)                  :: iv1
  INTEGER, INTENT(IN OUT)                  :: igrid


  INTEGER :: if1
  REAL*8 long,lat,x,y,z,x1,y1,z1,d1x,d1y,d1z,d2x,d2y,d2z, vx,vy,vz,mag


  if1=fofv(iv1,1,igrid)
  long=flong(if1,igrid)
  lat=flat(if1,igrid)
  CALL ll2xyz(long,lat,x,y,z)
  if1=fofv(iv1,2,igrid)
  long=flong(if1,igrid)
  lat=flat(if1,igrid)
  CALL ll2xyz(long,lat,x1,y1,z1)
  d1x=x-x1
  d1y=y-y1
  d1z=z-z1
  if1=fofv(iv1,3,igrid)
  long=flong(if1,igrid)
  lat=flat(if1,igrid)
  CALL ll2xyz(long,lat,x1,y1,z1)
  d2x=x-x1
  d2y=y-y1
  d2z=z-z1
  vx=(d1y*d2z-d2y*d1z)
  vy=(d1z*d2x-d2z*d1x)
  vz=(d1x*d2y-d2x*d1y)
  !     Make sure it's in the right hemisphere
  IF ((vx*x+vy*y+vz*z) < 0.0D0) THEN
     vx=-vx
     vy=-vy
     vz=-vz
  END IF
  !     Project back onto the sphere
  mag=SQRT(vx*vx+vy*vy+vz*vz)
  vx=vx/mag
  vy=vy/mag
  vz=vz/mag
  !     Convert back to latitude/longitude
  CALL xyz2ll(vx,vy,vz,long,lat)
  vlong(iv1,igrid)=long
  vlat(iv1,igrid)=lat

  RETURN
END SUBROUTINE vrtcrd

!     =======================================================================

SUBROUTINE penalt(ie1,igrid,cost,dd)

  USE grid

  !     Find the contribution to the Heikes+Randall penalty function
  !     from edge IE1

  IMPLICIT NONE

  INTEGER, INTENT(IN OUT)                  :: ie1
  INTEGER, INTENT(IN OUT)                  :: igrid
  REAL*8, INTENT(OUT)                      :: cost
  REAL*8, INTENT(OUT)                      :: dd

  INTEGER :: iv1,if1
  REAL*8 long,lat,x1,y1,z1,x2,y2,z2,xm,ym,zm,x,y,z,dx,dy,dz,l2, r2, mag,d2


  !     First find the ends of the edge
  iv1=vofe(1,ie1,igrid)
  long=vlong(iv1,igrid)
  lat=vlat(iv1,igrid)
  CALL ll2xyz(long,lat,x1,y1,z1)
  iv1=vofe(2,ie1,igrid)
  long=vlong(iv1,igrid)
  lat=vlat(iv1,igrid)
  CALL ll2xyz(long,lat,x2,y2,z2)

  !     Calculate length
  dx=x1-x2
  dy=y1-y2
  dz=z1-z2
  l2=dx*dx+dy*dy+dz*dz

  !     And midpoint
  xm=0.5D0*(x1+x2)
  ym=0.5D0*(y1+y2)
  zm=0.5D0*(z1+z2)
  mag=1.0D0/SQRT(xm*xm+ym*ym+zm*zm)
  xm=xm*mag
  ym=ym*mag
  zm=zm*mag

  !     Find faces either side of edge
  if1=fnxte(1,ie1,igrid)
  long=flong(if1,igrid)
  lat=flat(if1,igrid)
  CALL ll2xyz(long,lat,x1,y1,z1)
  if1=fnxte(2,ie1,igrid)
  long=flong(if1,igrid)
  lat=flat(if1,igrid)
  CALL ll2xyz(long,lat,x2,y2,z2)

  !     Find midpoint
  x=0.5D0*(x1+x2)
  y=0.5D0*(y1+y2)
  z=0.5D0*(z1+z2)
  mag=1.0D0/SQRT(x*x+y*y+z*z)
  x=x*mag
  y=y*mag
  z=z*mag

  !     Contribution to penalty function
  dx=x-xm
  dy=y-ym
  dz=z-zm
  d2=dx*dx+dy*dy+dz*dz
  dd=SQRT(d2)
  r2=d2/l2
  cost=r2*r2
  !      print *,'Edge ',IE1,'  displacement ',DD

  RETURN
END SUBROUTINE penalt

!     =======================================================================

SUBROUTINE penalf1(if1,igrid,totcst)

  USE grid

  !     Find the contribution to the Heikes+Randall penalty function
  !     affected by face IF1

  IMPLICIT NONE

  INTEGER, INTENT(IN OUT)                  :: if1
  INTEGER, INTENT(IN OUT)                  :: igrid
  REAL*8, INTENT(OUT)                      :: totcst

  INTEGER :: iv1,ie1,je1,je2,ielist(12),i,j
  REAL*8 cost, dd


  IF (if1 <= 12) THEN
     PRINT *,'Attempt to find cost function PENALF1 for pentagon'
     STOP
  END IF

  !     Find vertices that are affected and update their coordinates
  !     (IF1 will never be a pentagon, so only consider hexagons)
  !     List edges that are affected by face IF1
  je1=0
  DO i=1,6
     iv1=voff(if1,i,igrid)
     CALL vrtcrd(iv1,igrid)
     DO j=1,3
        ie1=eofv(j,iv1,igrid)
        !         Check we haven't got this one already
        je2=1
300     CONTINUE
        IF (je2 > je1) GO TO 100
        IF (ielist(je2) == ie1) GO TO 200
        je2=je2+1
        GO TO 300
100     CONTINUE
        !         Haven't got it so add it to the list
        je1=je1+1
        ielist(je1)=ie1
200     CONTINUE
     END DO
  END DO

  !     Now add up contributions to penalty function
  totcst=0.0D0
  DO je1=1,12
     ie1=ielist(je1)
     CALL penalt(ie1,igrid,cost,dd)
     totcst=totcst+cost
  END DO

  RETURN
END SUBROUTINE penalf1

!     ===================================================================

SUBROUTINE penalc(if1,igrid,cost,dd)

  USE grid

  !     Find the contribution to the penalty function that measures departures
  !     from centroidal for face IF1
  !     PENALF1 should be called first to ensure that vertices are up to date

  IMPLICIT NONE

  INTEGER, INTENT(IN OUT)                  :: if1
  INTEGER, INTENT(IN OUT)                  :: igrid
  REAL*8, INTENT(OUT)                      :: cost
  REAL*8, INTENT(OUT)                      :: dd

  INTEGER :: iv1,iv2,ie1,ie2,NE
  REAL*8  long,lat,x0,y0,z0,x1,y1,z1,x2,y2,z2,aface,  &
       xc,yc,zc,r2,dx,dy,dz,mag,da,d2


  IF (if1 <= 12) THEN
     NE = 5
  ELSE
     NE = 6
  END IF

  !     Coordinates of centre of face
  long=flong(if1,igrid)
  lat=flat(if1,igrid)
  CALL ll2xyz(long,lat,x0,y0,z0)
  !     Loop over edges in turn and calculate area of triangle
  !     formed by the edge and the centre of the face
  !     Hence find area of face and centroid
  xc=0.0D0
  yc=0.0D0
  zc=0.0D0
  aface=0.0D0
  DO  ie1=1,NE
     ie2=eoff(ie1,if1,igrid)
     iv1=vofe(1,ie2,igrid)
     iv2=vofe(2,ie2,igrid)
     long=vlong(iv1,igrid)
     lat=vlat(iv1,igrid)
     CALL ll2xyz(long,lat,x1,y1,z1)
     long=vlong(iv2,igrid)
     lat=vlat(iv2,igrid)
     CALL ll2xyz(long,lat,x2,y2,z2)
     CALL starea2(x0,y0,z0,x1,y1,z1,x2,y2,z2,da)
     aface=aface+da
     xc=xc+(x0+x1+x2)*da/3.0D0
     yc=yc+(y0+y1+y2)*da/3.0D0
     zc=zc+(z0+z1+z2)*da/3.0D0
  END DO
  mag=SQRT(xc*xc+yc*yc+zc*zc)
  xc=xc/mag
  yc=yc/mag
  zc=zc/mag

  !     Contribution to penalty function
  dx=x0-xc
  dy=y0-yc
  dz=z0-zc
  d2=dx*dx+dy*dy+dz*dz
  dd=SQRT(d2)
  r2=d2/aface
  cost=r2*r2

  RETURN
END SUBROUTINE penalc

!     ===================================================================

SUBROUTINE penalf2(if1,igrid,totcst)

  USE grid

  !     Find the contribution to the centroidal penalty function
  !     affected by face IF1

  IMPLICIT NONE

  INTEGER, INTENT(IN)                      :: if1
  INTEGER, INTENT(IN OUT)                  :: igrid
  REAL*8, INTENT(OUT)                      :: totcst


  INTEGER :: if2,iflist(7),i
  REAL*8 cost, dd


  IF (if1 <= 12) THEN
     PRINT *,'Attempt to find cost function PENALF2 for pentagon'
     STOP
  END IF

  !     Find neighbouring faces, i.e. those whose contribution is
  !     affected by IF1
  !     (IF1 will never be a pentagon, so only consider hexagons)
  !     List faces that are affected by face IF1
  DO i=1,6
     if2=fnxtf(if1,i,igrid)
     iflist(i)=if2
  END DO
  iflist(7)=if1

  !     Now add up contributions to penalty function
  totcst=0.0D0
  DO i=1,7
     if2=iflist(i)
     CALL penalc(if2,igrid,cost,dd)
     totcst=totcst+cost
  END DO

  RETURN
END SUBROUTINE penalf2

!     ===================================================================
!C

SUBROUTINE pen(if1,igrid,totcst)


  USE grid

  !     Find the contribution to the penalty function affected by face IF1
  !     Weighted contributions from the Heikes-Randall penalty
  !     function and the centroidal penalty function

  IMPLICIT NONE

  INTEGER, INTENT(IN OUT)                  :: if1
  INTEGER, INTENT(IN OUT)                  :: igrid
  REAL*8, INTENT(OUT)                      :: totcst

  REAL*8 c1,c2, dd
  REAL*8 weight1,weight2

  COMMON /comwgt/ weight1,weight2


  CALL penalf1(if1,igrid,c1) !HR95 penalty cost for this cell
  CALL penalf2(if1,igrid,c2) !Peixoto - remove Centroidal optimization to reduce overhead
  totcst=weight1*c1 + weight2*c2

  RETURN
END SUBROUTINE pen

!     ===================================================================

SUBROUTINE starea2(x0,y0,z0,x1,y1,z1,x2,y2,z2,area)

  !     Calculate the area of the spherical triangle whose corners
  !     have Cartesian coordinates (X0,Y0,Z0), (X1,Y1,Z1), (X2,Y2,Z2)
  !     The formula below is more robust to roundoff error than the
  !     better known sum of angle - PI formula

  IMPLICIT NONE

  REAL*8, INTENT(IN OUT)                   :: x0
  REAL*8, INTENT(IN OUT)                   :: y0
  REAL*8, INTENT(IN OUT)                   :: z0
  REAL*8, INTENT(IN OUT)                   :: x1
  REAL*8, INTENT(IN OUT)                   :: y1
  REAL*8, INTENT(IN OUT)                   :: z1
  REAL*8, INTENT(IN OUT)                   :: x2
  REAL*8, INTENT(IN OUT)                   :: y2
  REAL*8, INTENT(IN OUT)                   :: z2
  REAL*8, INTENT(OUT)                      :: area

  REAL*8  d0,d1,d2,s,t0,t1,t2,t3


  !     Distances between pairs of points
  CALL spdist(x0,y0,z0,x1,y1,z1,d2)
  CALL spdist(x1,y1,z1,x2,y2,z2,d0)
  CALL spdist(x2,y2,z2,x0,y0,z0,d1)

  !     Half perimeter
  s=0.5D0*(d0+d1+d2)

  !     Tangents
  t0 = TAN(0.5D0*(s-d0))
  t1 = TAN(0.5D0*(s-d1))
  t2 = TAN(0.5D0*(s-d2))
  t3 = TAN(0.5D0*s)

  !     Area
  area = 4.0D0*ATAN(SQRT(t0*t1*t2*t3))

  RETURN
END SUBROUTINE starea2

!     ===================================================================

SUBROUTINE spdist(x1,y1,z1,x2,y2,z2,s)

  !     Calculate the spherical distance S between two points with Cartesian
  !     coordinates (X1,Y1,Z1), (X2,Y2,Z2) on the unit sphere

  IMPLICIT NONE

  REAL*8, INTENT(IN)                       :: x1
  REAL*8, INTENT(IN)                       :: y1
  REAL*8, INTENT(IN)                       :: z1
  REAL*8, INTENT(IN)                       :: x2
  REAL*8, INTENT(IN)                       :: y2
  REAL*8, INTENT(IN)                       :: z2
  REAL*8, INTENT(OUT)                      :: s

  REAL*8  dx, dy, dz, ad


  dx = x2 - x1
  dy = y2 - y1
  dz = z2 - z1
  ad = SQRT(dx*dx + dy*dy + dz*dz)
  s = 2.0D0*ASIN(0.5D0*ad)


  RETURN
END SUBROUTINE spdist

!     ===================================================================


