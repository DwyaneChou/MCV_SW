MODULE tend_mod
  use parameters_mod
  
  implicit none
  
  ! MCV basic definiton
  type tend_field
    real, dimension(:,:,:), allocatable :: u
    real, dimension(:,:,:), allocatable :: v
    real, dimension(:,:,:), allocatable :: phiG
  end type tend_field
  
  type(tend_field), dimension(:), allocatable :: tend ! allocated by n time points, which is used by temporal integration schemes
  
  contains
  
  subroutine initTend
    integer :: iT
    
    allocate( tend(-nIntegralSubSteps:1) )
    
    do iT = -nIntegralSubSteps, 1
      allocate(tend(iT)%u   (ids:ide,jds:jde,ifs:ife))
      allocate(tend(iT)%v   (ids:ide,jds:jde,ifs:ife))
      allocate(tend(iT)%phiG(ids:ide,jds:jde,ifs:ife))
    enddo
    
  end subroutine initTend
  
END MODULE tend_mod

