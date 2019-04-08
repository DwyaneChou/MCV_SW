module initialize_mod
  use parameters_mod, only: initParameters
  use mesh_mod      , only: initMesh
  use ghost_mod     , only: initGhost
  use domain_mod    , only: initDomain
  use test_case_mod , only: initTestCase
  implicit none
  
  contains
    
  subroutine initMCV
  
    call initParameters
    call initMesh
    call initGhost
    call initDomain
    call initTestCase
  
  end subroutine initMCV
    
end module initialize_mod
    