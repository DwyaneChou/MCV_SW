module parameters_mod
  implicit none
  
  ! Namelist parameters
  ! time_settings
  integer :: run_days
  integer :: run_hours
  integer :: run_minutes
  integer :: run_seconds
  real    :: dt               ! time step
  integer :: history_interval ! output interval in seconds
  
  ! Domain
  real    :: dx     !  grid-spacing in the x-direction
  real    :: dy     !  grid-spacing in the y-direction
  
  integer :: xhalo  !  halo number of x-diretion
  integer :: yhalo  !  halo number of y-diretion
  
  ! Model run time control variables
  integer :: total_run_time   ! total run time for this model in seconds, this variable is determined by run_days, run_hours ...
  integer :: total_run_steps  ! total run steps for this model in seconds, this variable is determined by total_run_time and dt
  
  namelist /time_settings/ dt               ,&
                           run_days         ,&
                           run_hours        ,&
                           run_minutes      ,&
                           run_seconds      ,&
                           history_interval
  
  namelist /domain/        dx          ,&
                           dy          ,&
                           xhalo       ,&
                           yhalo
  
  contains
  
  subroutine readNamelist
    
    open(1, file = 'namelist.input',status='old')
    read(1, nml  = time_settings)
    read(1, nml  = domain       )
    close(1)
    
  end subroutine readNamelist
  
  subroutine initParameters
    
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
    
    call readNamelist
    
    total_run_time  = run_days * 86400 + run_hours * 3600 + run_minutes * 60 + run_seconds
    
    total_run_steps = ceiling(total_run_time/dt)
    
  end subroutine initParameters
  
end module parameters_mod
    