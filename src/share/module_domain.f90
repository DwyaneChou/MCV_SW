MODULE domain
  use constants  , only: DOF, D2R
  use parameters , only: dt, dx, dy, xhalo, yhalo
  
  implicit none
  
  ! Index parameter
  integer :: ids      ! The starting index in the x-direction (Physical domain)
  integer :: ide      ! The ending index in the x-direction  (Physical domain)
  integer :: jds      ! The starting index in the y-direction  (Physical domain)
  integer :: jde      ! The ending index in the y-direction  (Physical domain)
                      
  integer :: ics      ! The starting index in the x-direction (Physical cell/element domain)
  integer :: ice      ! The ending index in the x-direction  (Physical cell/element domain)
  integer :: jcs      ! The starting index in the y-direction  (Physical cell/element domain)
  integer :: jce      ! The ending index in the y-direction  (Physical cell/element domain)
                      
  integer :: ims      ! The starting index in the x-direction (Memory domain)
  integer :: ime      ! The ending index in the x-direction  (Memory domain)
  integer :: jms      ! The starting index in the y-direction  (Memory domain)
  integer :: jme      ! The ending index in the y-direction  (Memory domain)
                      
  integer :: Nx       ! Element numbers in the x-direction
  integer :: Ny       ! Element numbers in the y-direction
  
  integer :: nPVx     ! Point-value numbers in the x-direction
  integer :: nPVy     ! Point-value numbers in the y-direction
  
  integer :: Nlambda  ! grid points in the lambda direction
  integer :: Ntheta   ! grid points in the theta direction
  
  integer, parameter :: Nf = 6           ! Number of cube faces
  
  real    :: x_min = -45.d0   !  start location of x-direction
  real    :: x_max =  45.d0   !  end location of x-direction
  real    :: y_min = -45.d0   !  start location of y-direction
  real    :: y_max =  45.d0   !  end location of y-direction
  
  ! MCV basic definiton
  real, dimension(:,:,:,:), allocatable :: PV   !  Point Values.
  real, dimension(:,:,:,:), allocatable :: VIA  !  Volume Integrated Average.
  
  ! Jacobian and Metric matrices
  real, dimension(:,:    ), allocatable :: jacobTransform !  jacobian of Transformation
  real, dimension(:,:,:  ), allocatable :: metricTensor   !  horizontal metrics Tensor
  real, dimension(:,:,:,:), allocatable :: compositeJacob !  composite Jacobian of transformation
  
  ! ghost cell location
  integer, dimension(:,:,:,:), allocatable :: ghostCellIndex
  real   , dimension(:,:,:,:), allocatable :: ghostCellPosition
  
  ! coordinate
  real, dimension(:), allocatable :: x
  real, dimension(:), allocatable :: y
  real, dimension(:), allocatable :: xc
  real, dimension(:), allocatable :: yc
  
  contains
  
  subroutine initDomain
    
    ! Check if dx and dy are avaliable to achieve the integral element number
    if( mod(x_max - x_min, dx) /= 0 )then
        stop '90 divide dx must be integer, choose another dx'
    end if
    
    if( mod(y_max - y_min, dy) /= 0 )then
        stop '90 divide dy must be integer, choose another dy'
    end if
    
    ! Calculate element numbers on x/y direction
    Nx = (x_max - x_min)/dx + 1
    Ny = (y_max - y_min)/dy + 1
    
    ! Calculate PV number on x/y direction
    nPVx = Nx * (DOF - 1) + 1
    nPVy = Ny * (DOF - 1) + 1
    
    ! Calculate grid number on longitude/latitude coordinate
    Nlambda = 360/dx + 1
    Ntheta  = 180/dy + 1
    
    ! Calculate starting and ending index for physical domain
    ids  = 1
    ide  = Nx
    jds  = 1
    jde  = Ny
  
    ! Calculate starting and ending index for element cell
    ics  = 1  - xhalo
    ice  = Nx + xhalo
    jcs  = 1  - yhalo
    jce  = Ny + yhalo
    
    ! Calculate starting and ending index for memory array
    ims  = 1  - xhalo * (DOF - 1)
    ime  = Nx + xhalo * (DOF - 1)
    jms  = 1  - yhalo * (DOF - 1)
    jme  = Ny + yhalo * (DOF - 1)
    
    ! Convert Degree to Radian
    dx    = dx    * D2R
    dy    = dt    * D2R
    x_min = x_min * D2R
    x_max = x_max * D2R
    y_min = y_min * D2R
    y_max = y_max * D2R
    
  end subroutine initDomain
  
END MODULE domain

