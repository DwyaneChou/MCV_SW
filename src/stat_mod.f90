MODULE stat_mod
  use parameters_mod
  
  implicit none
  
  ! MCV basic definiton
  type stat_field
    real, dimension(:,:,:), allocatable :: uP   ! covariant wind on x direction on points
    real, dimension(:,:,:), allocatable :: vP   ! covariant wind on y direction on points
    real, dimension(:,:,:), allocatable :: phiP ! geopotential height on points
    
    real, dimension(:,:,:), allocatable :: uC   ! covariant wind on x direction on cells
    real, dimension(:,:,:), allocatable :: vC   ! covariant wind on y direction on cells
    real, dimension(:,:,:), allocatable :: phiC ! geopotential height on cells
  end type stat_field
  
  type(stat_field), dimension(:), allocatable :: stat ! allocated by n time points, which is used by temporal integration schemes
  
  contains
  
  subroutine initStat
    integer :: iT
    
    allocate( stat(-nIntegralSubSteps:1) )
    
    do iT = -nIntegralSubSteps, 1
      allocate(stat(iT)%uP  (ips:ipe,jps:jpe,ifs:ife))
      allocate(stat(iT)%vP  (ips:ipe,jps:jpe,ifs:ife))
      allocate(stat(iT)%phiP(ips:ipe,jps:jpe,ifs:ife))
      
      allocate(stat(iT)%uC  (ics:ice,jcs:jce,ifs:ife))
      allocate(stat(iT)%vC  (ics:ice,jcs:jce,ifs:ife))
      allocate(stat(iT)%phiC(ics:ice,jcs:jce,ifs:ife))
    enddo
    
  end subroutine initStat
  
END MODULE stat_mod

