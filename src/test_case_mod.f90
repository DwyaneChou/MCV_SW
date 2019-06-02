module test_case_mod
  use parameters_mod
  use mesh_mod
  use ghost_mod
  use stat_mod
  use projection_mod
  implicit none
  
  contains
    
  subroutine initTestCase
    if(case_num == 2) call case2
    
    print*,''
    print*,'max/min value of u   : ',maxval(stat(0)%u  ),minval(stat(0)%u  )
    print*,'max/min value of v   : ',maxval(stat(0)%v  ),minval(stat(0)%v  )
    print*,'max/min value of phi : ',maxval(stat(0)%phi),minval(stat(0)%phi)
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
          gh = gh0 - (radius * OMEGA * u0 + u0**2 / 2.) * sin(mesh%latP(iPV,jPV,iPatch))**2
          
          u  = u0 * cos(mesh%latP(iPV,jPV,iPatch))
          
          call contravProjSphere2Plane(contraU, contraV, u      , 0.     , mesh%matrixIA(:,:,iPV,jPV,iPatch))
          call contrav2cov            (coU    , coV    , contraU, contraV, mesh%matrixG (:,:,iPV,jPV,iPatch))
          
          stat(0)%phi(iPV,jPV,iPatch) = gh
          stat(0)%u  (iPV,jPV,iPatch) = coU
          stat(0)%v  (iPV,jPV,iPatch) = coV
        enddo
      enddo
    enddo
    
    mesh%dphisdx = 0.
    mesh%dphisdy = 0.
  end subroutine case2
  
  ! Rossby-Haurwitz wave with wavenumber 4
  subroutine case6
  
  end subroutine case6
  
end module test_case_mod
    