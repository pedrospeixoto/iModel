module constants
  !======================================================
  !
  !	Module containing parameters needed for several 
  !       other routines
  !
  !
  !   Pedro Peixoto - 03-03-2011
  !======================================================
  implicit none
  save 
  public

  !---------------------------------------------------
  !Kind attributions
  !---------------------------------------------------

  integer, parameter :: i2  = selected_real_kind(2,20)
  integer, parameter :: i4  = selected_real_kind(6,20)
  integer, parameter :: i8  = selected_real_kind(14,40)

  integer, parameter :: r4  = selected_real_kind(6,37)  
  integer, parameter :: r8  = selected_real_kind(12,100)
  integer, parameter :: r16 = max(r8,selected_real_kind(27,2400))

  !---------------------------------------------------
  ! General Parameters
  !---------------------------------------------------

  !Pi 
  real(r8), parameter :: pi   = 4._r8* 0.78539816339744830961566 !datan (1._r8)
  real(r8), parameter :: pi2  = 2._r8*pi
  real(r8), parameter :: pio2 = pi/2._r8
  real(r8), parameter :: piby2 = pi*0.5_r8
  real(r8), parameter :: pio4 = pi/4._r8

  !Degrees to radians coversion (multiply to obtain conversion)
  real(r8), parameter :: deg2rad = pi / 180._r8

  !Radians to Degrees coversion (multiply to obtain conversion)
  real(r8), parameter :: rad2deg = 1._r8/deg2rad

  !Very small real number (aprox 1.10e-7)
  !real(r8), parameter :: eps = epsilon(1.)
  real(r8), parameter :: eps = epsilon(1.)

  !Very very small real number (aprox 1.10e-16)
  real(r8), parameter :: eps2 = epsilon(pi)

  !---------------------------------------------------
  ! Physical Parameters
  !---------------------------------------------------

  ! Earth mean radius (meters)
  real(r8), parameter :: erad     = 6.37122e6_r8
  real(r8), parameter :: rearth     = 6.37122e6_r8
  real(r8), parameter :: eradi    = 1._r8/6.37122e6_r8
  real(r8), parameter :: unitspharea    = 4._r8*pi

  !Gravitational accerlaration of the Earth (m/s^2)
  real(r8), parameter  :: grav    = 9.80616_r8
  real(r8), parameter  :: gravity = 9.80616_r8
  real(r8), parameter  :: gravi   = 1._r8/9.80616_r8

  ! Angular velocity of the Earth (rot/s)
  real (r8), parameter :: omega   = 7.292e-5_r8
  real (r8), parameter :: rotatn   = 7.292e-5_r8

  !Days to seconds
  real (r8), parameter :: day2sec = 86400_r8
  real (r8), parameter :: sec2day = 1._r8/86400_r8

  ! Dry air gas constant [J/(kg K)]
  real(r8), parameter :: rdry     = 287.   

  ! Dry air spec heat at const P [J/(kg K)]
  real(r8), parameter :: cp       = 1004.  

  ! Dry air spec heat at const vol [J/(kg K)]
  real(r8), parameter :: cv       = 717.              

  ! Water vapor gas constant [J/(kg K)]
  real(r8), parameter :: rvap     = 461.               

  ! Reference pressure [Pa]
  real(r8), parameter :: p00      = 1.e5            

  ! 0 Celsius temperature [K]
  real(r8), parameter :: t00      = 273.15         

  !---------------------------------------------------
  ! Directories
  !--------------------------------------------------
  !Output data directory
  character (len=60), parameter::  datadir= "data/"
  !Grid data directory
  character (len=60), parameter::  griddir= "grid/"
  !Directory for ploting thing using GMT
  character (len=60), parameter::  gmtdir = "gmt/"
  !Sources directory
  character (len=60), parameter::  srcdir = "src/"
  !Parameter files  directory
  character (len=60), parameter::  pardir = "par/"
  !Reference solutions directory - data will be read from here
  character (len=60), parameter::  refdir = "ref/"

  !---------------------------------------------------
  ! Parameters files
  !-------------------------------------------------- 
  !character (len=60), parameter::  deformpar = "deform.par"
  !character (len=60), parameter::  meshpar = "mesh.par"
  !character (len=60), parameter::  simulpar = "simul.par"

  !---------------------------------------------------
  ! Global constant variables
  !---------------------------------------------------
  !Flag for verbose output
  logical :: showonscreen=.false.

  !Simulation to be done
  integer (i4) :: simulcase

  !---------------------------------------------------
  ! Density Interpolation Parameters
  !---------------------------------------------------

  integer, parameter :: n_lat = 4*180+1
  integer, parameter :: n_lon = 4*360+1
  real (kind=8), parameter :: latmin = -pio2
  real (kind=8), parameter :: latmax = pio2
  real (kind=8), parameter :: lonmin = -pi
  real (kind=8), parameter :: lonmax = pi

  !---------------------------------------------------
  ! Andes Density Interpolation Parameters
  !---------------------------------------------------

  real (kind=8), parameter :: latmin_andes = -50
  real (kind=8), parameter :: latmax_andes = 10
  real (kind=8), parameter :: lonmin_andes = -100
  real (kind=8), parameter :: lonmax_andes = -55 

  integer, parameter :: nlat_alt = 4*180+1 
  integer, parameter :: nlon_alt = 4*360+1
  integer, parameter :: nx_andes = (latmax_andes - latmin_andes)*4+1
  integer, parameter :: ny_andes = (lonmax_andes - lonmin_andes)*4+1

  character (len=60), parameter::  altdir = "altitude/"

end module constants
