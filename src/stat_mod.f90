MODULE stat_mod
  use parameters_mod
  
  implicit none
  
  ! MCV basic definiton
  type stat_field
    real, dimension(:,:,:), allocatable :: u                ! covariant wind on x direction on points
    real, dimension(:,:,:), allocatable :: v                ! covariant wind on y direction on points
    real, dimension(:,:,:), allocatable :: phi              ! geopotential height on points
    real, dimension(:,:,:), allocatable :: phiG             ! phi * sqrtG on points
    real, dimension(:,:,:), allocatable :: contraU
    real, dimension(:,:,:), allocatable :: contraV
    real, dimension(:,:,:), allocatable :: zonal_wind
    real, dimension(:,:,:), allocatable :: meridional_wind
    
    real, dimension(:,:,:), allocatable :: uC               ! covariant wind on x direction on cells
    real, dimension(:,:,:), allocatable :: vC               ! covariant wind on y direction on cells
    real, dimension(:,:,:), allocatable :: phiC             ! geopotential height on points
  end type stat_field
  
  type(stat_field), dimension(:), allocatable :: stat  ! allocated by n time points, which is used by temporal integration schemes
  type(stat_field)                            :: statC
  contains
  
  subroutine initStat
    integer :: iT
    
    allocate( stat(-nIntegralSubSteps:1) )
    
    do iT = -nIntegralSubSteps, 1
      allocate(stat(iT)%u              (ips:ipe,jps:jpe,ifs:ife))
      allocate(stat(iT)%v              (ips:ipe,jps:jpe,ifs:ife))
      allocate(stat(iT)%phi            (ips:ipe,jps:jpe,ifs:ife))
      allocate(stat(iT)%phiG           (ips:ipe,jps:jpe,ifs:ife))
      allocate(stat(iT)%contraU        (ips:ipe,jps:jpe,ifs:ife))
      allocate(stat(iT)%contraV        (ips:ipe,jps:jpe,ifs:ife))
      allocate(stat(iT)%zonal_wind     (ips:ipe,jps:jpe,ifs:ife))
      allocate(stat(iT)%meridional_wind(ips:ipe,jps:jpe,ifs:ife))
      
      allocate(stat(iT)%uC             (ics:ice,jcs:jce,ifs:ife))
      allocate(stat(iT)%vC             (ics:ice,jcs:jce,ifs:ife))
      allocate(stat(iT)%phiC           (ics:ice,jcs:jce,ifs:ife))
    enddo
    
  end subroutine initStat
  
  subroutine copyStat(stat_out,stat_in)
    type(stat_field),intent(out) :: stat_out
    type(stat_field),intent(in ) :: stat_in
  
    stat_out%u               = stat_in%u
    stat_out%v               = stat_in%v
    stat_out%phi             = stat_in%phi
    stat_out%phiG            = stat_in%phiG
    stat_out%contraU         = stat_in%contraU
    stat_out%contraV         = stat_in%contraV
    stat_out%zonal_wind      = stat_in%zonal_wind
    stat_out%meridional_wind = stat_in%meridional_wind
  end subroutine copyStat
  
END MODULE stat_mod

