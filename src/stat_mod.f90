MODULE stat_mod
  use parameters_mod
  
  implicit none
  
  ! MCV basic definiton
  type stat_field
    real, dimension(:,:,:), allocatable :: u   ! covariant wind on x direction on points
    real, dimension(:,:,:), allocatable :: v   ! covariant wind on y direction on points
    real, dimension(:,:,:), allocatable :: phi ! geopotential height on points
  end type stat_field
  
  type(stat_field), dimension(:), allocatable :: stat ! allocated by n time points, which is used by temporal integration schemes
  
  contains
  
  subroutine initStat
    integer :: iT
    
    allocate( stat(-nIntegralSubSteps:1) )
    
    do iT = -nIntegralSubSteps, 1
      allocate(stat(iT)%u  (ips:ipe,jps:jpe,ifs:ife))
      allocate(stat(iT)%v  (ips:ipe,jps:jpe,ifs:ife))
      allocate(stat(iT)%phi(ips:ipe,jps:jpe,ifs:ife))
    enddo
    
  end subroutine initStat
  
END MODULE stat_mod

