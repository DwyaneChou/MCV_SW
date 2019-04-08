module parameters_mod
  use constants_mod
  implicit none
  
  ! Namelist parameters
  ! time_settings
  integer :: run_days
  integer :: run_hours
  integer :: run_minutes
  integer :: run_seconds
  real    :: dt               ! time step
  integer :: history_interval ! output interval in seconds
  
  character*200 :: integral_scheme
  
  ! Case select
  integer :: case_num
  
  ! Domain
  real    :: dx     !  grid-spacing in the x-direction
  real    :: dy     !  grid-spacing in the y-direction
  
  integer :: xhalo  !  halo number of x-diretion
  integer :: yhalo  !  halo number of y-diretion
  
  ! Index parameter
  integer :: ids      ! The starting index in the x-direction (Physical domain)
  integer :: ide      ! The ending index in the x-direction  (Physical domain)
  integer :: jds      ! The starting index in the y-direction  (Physical domain)
  integer :: jde      ! The ending index in the y-direction  (Physical domain)
                      
  integer :: ics      ! The starting index in the x-direction (Physical cell/element domain)
  integer :: ice      ! The ending index in the x-direction  (Physical cell/element domain)
  integer :: jcs      ! The starting index in the y-direction  (Physical cell/element domain)
  integer :: jce      ! The ending index in the y-direction  (Physical cell/element domain)
                      
  integer :: ips      ! The starting index in the x-direction (PV domain)
  integer :: ipe      ! The ending index in the x-direction  (PV domain)
  integer :: jps      ! The starting index in the y-direction  (PV domain)
  integer :: jpe      ! The ending index in the y-direction  (PV domain)
  
  integer :: ifs      ! The starting index of patch(face)
  integer :: ife      ! The ending index of patch(face)
                      
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
  
  ! Model run time control variables
  integer :: total_run_time   ! total run time for this model in seconds, this variable is determined by run_days, run_hours ...
  integer :: total_run_steps  ! total run steps for this model in seconds, this variable is determined by total_run_time and dt
  
  namelist /time_settings/ dt               ,&
                           run_days         ,&
                           run_hours        ,&
                           run_minutes      ,&
                           run_seconds      ,&
                           history_interval ,&
                           integral_scheme
  
  namelist /case_select/   case_num
  
  namelist /domain/        dx          ,&
                           dy          ,&
                           xhalo       ,&
                           yhalo
  
  contains
  
  subroutine readNamelist
    
    open(1, file = 'namelist.input',status='old')
    read(1, nml  = time_settings)
    read(1, nml  = case_select  )
    read(1, nml  = domain       )
    close(1)
    
  end subroutine readNamelist
  
  subroutine initParameters
    
    ! Setting default values
    dt    = 300.
    dx    = 2.
    dy    = 2.
    xhalo = 1
    yhalo = 1
    
    run_days         = 1
    run_hours        = 0
    run_minutes      = 0
    run_seconds      = 0
    history_interval = 3600
    integral_scheme  = 'RK4'
    
    ! Read namelist
    call readNamelist
    
    ! Calculate total run time in seconds
    total_run_time  = run_days * 86400 + run_hours * 3600 + run_minutes * 60 + run_seconds
    
    ! Calculate total run steps
    total_run_steps = ceiling(total_run_time/dt)
    
    ! Check if dx and dy are avaliable to achieve the integral element number
    if( (x_max - x_min)/dx - int((x_max - x_min)/dx) /= 0. )then
      stop '90 divide dx must be integer, choose another dx'
    end if
    
    if( (y_max - y_min)/dy - int((y_max - y_min)/dy) /= 0. )then
        stop '90 divide dy must be integer, choose another dy'
    end if
    
    ! Calculate element numbers on x/y direction
    Nx = int((x_max - x_min)/dx)
    Ny = int((y_max - y_min)/dy)
    
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
    ips  = 1    - xhalo * (DOF - 1)
    ipe  = nPVx + xhalo * (DOF - 1)
    jps  = 1    - yhalo * (DOF - 1)
    jpe  = nPVy + yhalo * (DOF - 1)
    
    ! Setting the starting patch index and ending patch index
    ifs = 1
    ife = Nf
    
    ! Convert Degree to Radian
    dx    = dx    * D2R
    dy    = dy    * D2R
    x_min = x_min * D2R
    x_max = x_max * D2R
    y_min = y_min * D2R
    y_max = y_max * D2R
    
    ! Setting the number of substeps in temporal integration scheme
    if(trim(adjustl(integral_scheme)) == 'RK4')then
      nIntegralSubSteps = 4
    else
      stop 'Unknown integral scheme, please select from RK4 ...'
    endif
    
  end subroutine initParameters
  
end module parameters_mod
    