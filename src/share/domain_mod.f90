MODULE domain_mod
  use constants_mod  , only: DOF, D2R, radius
  use parameters_mod , only: dt, dx, dy, xhalo, yhalo, integral_scheme
  
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
  
  integer :: nIntegralSubSteps ! number of integral substeps in temporal integration scheme
  
  integer, parameter :: Nf = 6           ! Number of cube faces
  
  real    :: x_min = -45.   !  start location of x-direction
  real    :: x_max =  45.   !  end location of x-direction
  real    :: y_min = -45.   !  start location of y-direction
  real    :: y_max =  45.   !  end location of y-direction
  
  ! MCV basic definiton
  type cell
    real, dimension(DOF,DOF) :: PV         ! Point Values.
    real                     :: VIA        ! Volume Integrated Average.
    real, dimension(DOF)     :: VIA_x      ! VIA on x direction
    real, dimension(DOF)     :: VIA_y      ! VIA on y direction
  end type cell
  
  type fields
    type(cell), dimension(:,:,:), allocatable :: U               ! covariant u-wind component, the 1st dimension is the patch index, 
                                                                 !                             the 2nd dimension is the x direcion index, 
                                                                 !                             the 3rd dimension is the y direcion index
    type(cell), dimension(:,:,:), allocatable :: V               ! covariant v-wind component, dimension setting is the same as U
    type(cell), dimension(:,:,:), allocatable :: contraU         ! contravariant u-wind component, dimension setting is the same as U
    type(cell), dimension(:,:,:), allocatable :: contraV         ! contravariant v-wind component, dimension setting is the same as U
    type(cell), dimension(:,:,:), allocatable :: phi             ! geopotential height, dimension setting is the same as covariantU
    
    type(cell), dimension(:,:,:), allocatable :: zonalWind       ! contravariant u-wind component, dimension setting is the same as U
    type(cell), dimension(:,:,:), allocatable :: meridionalWind  ! contravariant v-wind component, dimension setting is the same as U
  end type fields
  
  type(fields), dimension(:), allocatable :: state ! allocated by n time points, which is used by temporal integration schemes
  type(fields), dimension(:), allocatable :: tend  ! allocated by n time points, which is used by temporal integration schemes
  
  type fitOnCell ! 1st order derivatives determined by polynominal fitting
    real, dimension(DOF) :: derivLeft_x
    real, dimension(DOF) :: derivRight_x
    real, dimension(DOF) :: derivTop_y
    real, dimension(DOF) :: derivBottom_y
  end type fitOnCell
  
  type(fitOnCell), dimension(:,:,:), allocatable :: fitU
  type(fitOnCell), dimension(:,:,:), allocatable :: fitV
  type(fitOnCell), dimension(:,:,:), allocatable :: fitPhi
  
  ! Jacobian and Metric matrices
  real, dimension(:,:    ), allocatable :: jacobTransform !  jacobian of Transformation
  real, dimension(:,:,:  ), allocatable :: metricTensor   !  horizontal metrics Tensor
  real, dimension(:,:,:,:), allocatable :: compositeJacob !  composite Jacobian of transformation
  
  ! ghost cell location
  integer, dimension(:,:,:,:), allocatable :: ghostCellIndex
  real   , dimension(:,:,:,:), allocatable :: ghostCellPosition
  
  ! coordinate
  type mesh_info
    real, dimension(DOF,DOF) :: alpha              ! central angle on x direction for each patch
    real, dimension(DOF,DOF) :: beta               ! central angle on y direction for each patch
    real, dimension(DOF,DOF) :: x                  ! length of the arc on x direction for each patch, x = radius * alpha
    real, dimension(DOF,DOF) :: y                  ! length of the arc on y direction for each patch, y = radius * beta
    real, dimension(DOF,DOF) :: lon                ! longitude on sphere coordinate
    real, dimension(DOF,DOF) :: lat                ! latitude on sphere coordinate
    
    real, dimension(2,2,DOF,DOF) :: metricTensor   ! horizontal metric Tensor, which transform covariant vectors to contravariant vectors
    real, dimension(2,2,DOF,DOF) :: jacobTransform ! jacobian of Transformation, which transform the contravariant vectors to zonal vector and meridional vector
  end type mesh_info
  
  type(mesh_info), dimension(:,:,:), allocatable :: mesh
  
  contains
  
  subroutine initDomain
    integer :: i,j
    
    ! Check if dx and dy are avaliable to achieve the integral element number
    if( (x_max - x_min)/dx - int((x_max - x_min)/dx) /= 0. )then
      stop '90 divide dx must be integer, choose another dy'
    end if
    
    if( (y_max - y_min)/dy - int((y_max - y_min)/dy) /= 0. )then
        stop '90 divide dy must be integer, choose another dy'
    end if
    
    ! Calculate element numbers on x/y direction
    Nx = int((x_max - x_min)/dx) + 1
    Ny = int((y_max - y_min)/dy) + 1
    
    ! Calculate PV number on x/y direction
    nPVx = Nx * (DOF - 1) + 1
    nPVy = Ny * (DOF - 1) + 1
    
    ! Calculate grid number on longitude/latitude coordinate
    Nlambda = nPVx
    Ntheta  = nPVy
    
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
    
    ! Allocate the data structure
    if(trim(adjustl(integral_scheme)) == 'RK4')then
      nIntegralSubSteps = 4
    else
      stop 'Unknown integral scheme, please select from RK4 ...'
    endif
    
    allocate( state(nIntegralSubSteps) )
    allocate( tend (nIntegralSubSteps) )
    
    do i = 1,nIntegralSubSteps
      allocate( state(i)%U             (Nf, ics:ice, jcs:jce) )
      allocate( state(i)%V             (Nf, ics:ice, jcs:jce) )
      allocate( state(i)%contraU       (Nf, ics:ice, jcs:jce) )
      allocate( state(i)%contraV       (Nf, ics:ice, jcs:jce) )
      allocate( state(i)%phi           (Nf, ics:ice, jcs:jce) )
      allocate( state(i)%zonalWind     (Nf, ics:ice, jcs:jce) )
      allocate( state(i)%meridionalWind(Nf, ics:ice, jcs:jce) )
      
      allocate( tend (i)%U      (Nf, ics:ice, jcs:jce) )
      allocate( tend (i)%V      (Nf, ics:ice, jcs:jce) )
      !allocate( tend (i)%contraU(Nf, ics:ice, jcs:jce) )
      !allocate( tend (i)%contraV(Nf, ics:ice, jcs:jce) )
      allocate( tend (i)%phi    (Nf, ics:ice, jcs:jce) )
    enddo
    
    allocate( fitU  (Nf, ics:ice, jcs:jce) )
    allocate( fitV  (Nf, ics:ice, jcs:jce) )
    allocate( fitPhi(Nf, ics:ice, jcs:jce) )
    
  end subroutine initDomain
  
END MODULE domain_mod

