module initialize
  use parameters, only: initParameters
  use domain    , only: initDomain
  implicit none
  
  contains
    
  subroutine initMCV
  
    call initParameters
    call initDomain
  
  end subroutine initMCV
    
end module initialize
    