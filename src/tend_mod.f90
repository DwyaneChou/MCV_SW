MODULE tend_mod
  use parameters_mod
  
  implicit none
  
  ! MCV basic definiton
  type tend_field
    real, dimension(:,:,:), allocatable :: uP
    real, dimension(:,:,:), allocatable :: vP
    real, dimension(:,:,:), allocatable :: phiP
    
    real, dimension(:,:,:), allocatable :: uC
    real, dimension(:,:,:), allocatable :: vC
    real, dimension(:,:,:), allocatable :: phiC
  end type tend_field
  
  type(tend_field), dimension(:), allocatable :: tend ! allocated by n time points, which is used by temporal integration schemes
  
  contains
  
  subroutine initTend
    integer :: iT
    
    allocate( tend(-nIntegralSubSteps:1) )
    
    do iT = -nIntegralSubSteps, 1
      allocate(tend(iT)%uP  (ips:ipe,jps:jpe,Nf))
      allocate(tend(iT)%vP  (ips:ipe,jps:jpe,Nf))
      allocate(tend(iT)%phiP(ips:ipe,jps:jpe,Nf))
      
      allocate(tend(iT)%uC  (ics:ice,jcs:jce,Nf))
      allocate(tend(iT)%vC  (ics:ice,jcs:jce,Nf))
      allocate(tend(iT)%phiC(ics:ice,jcs:jce,Nf))
    enddo
    
  end subroutine initTend
  
END MODULE tend_mod

