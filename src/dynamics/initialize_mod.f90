module initialize_mod
  use parameters_mod, only: initParameters
  use domain_mod    , only: initDomain
  implicit none
  
  contains
    
  subroutine initMCV
  
    call initParameters
    call initDomain
  
  end subroutine initMCV
    
end module initialize_mod
    