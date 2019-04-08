MODULE domain_mod
  use constants_mod
  use parameters_mod
  
  implicit none
  
  ! MCV basic definiton
  type state_field
    real, dimension(:,:,:,:), allocatable :: PV        ! Point Values.
    real, dimension(:,:,:,:), allocatable :: VIA       ! Volume Integrated Average.
    real, dimension(:,:,:,:), allocatable :: VIAx      ! VIA on x direction
    real, dimension(:,:,:,:), allocatable :: VIAy      ! VIA on y direction
    real, dimension(:,:,:,:), allocatable :: dLeft     ! Left derivative on cell
    real, dimension(:,:,:,:), allocatable :: dRight    ! Right derivative on cell
    real, dimension(:,:,:,:), allocatable :: dTop      ! Top derivative on cell
    real, dimension(:,:,:,:), allocatable :: dBottom   ! Bottom derivative on cell
  end type state_field
  
  type tend_field
    real, dimension(:,:,:,:), allocatable :: PV         ! Point Values.
    real, dimension(:,:,:,:), allocatable :: VIA        ! Volume Integrated Average.
    real, dimension(:,:,:,:), allocatable :: VIAx       ! VIA on x direction
    real, dimension(:,:,:,:), allocatable :: VIAy       ! VIA on y direction
  end type tend_field
  
  type(state_field), dimension(:), allocatable :: state ! allocated by n time points, which is used by temporal integration schemes
  type(tend_field ), dimension(:), allocatable :: tend  ! allocated by n time points, which is used by temporal integration schemes
  
  contains
  
  subroutine initDomain
    integer :: iT
    
    allocate( state(-nIntegralSubSteps:1) )
    allocate( tend (-nIntegralSubSteps:1) )
    
    do iT = -nIntegralSubSteps, 1
      allocate( state(iT)%PV     (numVar, ips:ipe, jps:jpe, ifs:ife))
      allocate( state(iT)%VIA    (numVar, ics:ice, jcs:jce, ifs:ife))
      allocate( state(iT)%VIAx   (numVar, ics:ice, jps:jpe, ifs:ife))
      allocate( state(iT)%VIAy   (numVar, ips:ipe, jps:jpe, ifs:ife))
      allocate( state(iT)%dLeft  (numVar, ics:ice, jps:jpe, ifs:ife))
      allocate( state(iT)%dRight (numVar, ics:ice, jps:jpe, ifs:ife))
      allocate( state(iT)%dTop   (numVar, ips:ipe, jcs:jce, ifs:ife))
      allocate( state(iT)%dBottom(numVar, ips:ipe, jcs:jce, ifs:ife))
       
      allocate( tend (iT)%PV     (numVar, ips:ipe, jps:jpe, ifs:ife))
      allocate( tend (iT)%VIA    (numVar, ics:ice, jcs:jce, ifs:ife))
      allocate( tend (iT)%VIAx   (numVar, ics:ice, jps:jpe, ifs:ife))
      allocate( tend (iT)%VIAy   (numVar, ips:ipe, jps:jpe, ifs:ife))
    end do
    
  end subroutine initDomain
  
END MODULE domain_mod

