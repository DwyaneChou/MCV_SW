module test_case_mod
  use parameters_mod
  use mesh_mod
  use ghost_mod
  use domain_mod
  use projection_mod
  implicit none
  
  contains
    
  subroutine initTestCase
    if(case_num == 2) call case2
  end subroutine initTestCase

  
  ! Global steady flow
  subroutine case2
    real    :: u0
    real    :: gh0 = 29400.
    real    :: gh
    real    :: u, v
    real    :: contraU, contraV
    real    :: coU, coV
    
    integer :: iPatch, iPV, jPV
    
    u0 = 2. * pi * radius / (12. * 86400.)
    
    do iPatch = ifs, ife
      do jPV = jps, jpe
        do iPV = ips, ipe
          gh = gh0 - (radius * OMEGA * u0 + u0**2 / 2.) * sin(PVmesh%lat(iPV,jPV,iPatch))**2
          
          u  = u0 * cos(PVmesh%lat(iPV,jPV,iPatch))
          
          call contravProjSphere2Plane(contraU, contraV, u      , 0.     , PVmesh%matrixIA(:,:,iPV,jPV,iPatch))
          call contrav2cov            (coU    , coV    , contraU, contraV, PVmesh%matrixG (:,:,iPV,jPV,iPatch))
          
          state(0)%PV(1,iPV,jPV,iPatch) = gh
          state(0)%PV(2,iPV,jPV,iPatch) = coU
          state(0)%PV(3,iPV,jPV,iPatch) = coV
        enddo
      enddo
    enddo
    
    !hs = 0.
  end subroutine case2
  
  ! Rossby-Haurwitz wave with wavenumber 4
  subroutine case6
  
  end subroutine case6
  
end module test_case_mod
    