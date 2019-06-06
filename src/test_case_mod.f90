module test_case_mod
  use parameters_mod
  use mesh_mod
  use ghost_mod
  use stat_mod
  use projection_mod
  use spatial_operators_mod
  implicit none
  
  contains
    
  subroutine initTestCase
    if(case_num == 2) call case2(stat(0))
    
    print*,''
    print*,'max/min value of u   : ',maxval(stat(0)%u  ),minval(stat(0)%u  )
    print*,'max/min value of v   : ',maxval(stat(0)%v  ),minval(stat(0)%v  )
    print*,'max/min value of phi : ',maxval(stat(0)%phi),minval(stat(0)%phi)
  end subroutine initTestCase

  ! Global steady flow
  subroutine case2(stat)
    type(stat_field), intent(inout) :: stat
    real    :: u0
    real    :: gh0 = 29400.
    real    :: gh
    real    :: u, v
    
    integer :: i,j,iPatch
    
    stat%meridional_wind = 0.
    
    u0 = 2. * pi * radius / (12. * 86400.)
    do iPatch = ifs, ife
      do j = jps, jpe
        do i = ips, ipe
          stat%phi(i,j,iPatch) = gh0 - (radius * Omega * u0 + u0**2 / 2.) * sin(mesh%latP(i,j,iPatch))**2
          
          stat%zonal_wind(i,j,iPatch) = u0 * cos(mesh%latP(i,j,iPatch))
          
          call covProjSphere2Plane    (stat%u      (i,j,iPatch), stat%v      (i,j,iPatch), stat%zonal_wind(i,j,iPatch), stat%meridional_wind(i,j,iPatch), mesh%matrixIA(:,:,i,j,iPatch), mesh%matrixG(:,:,i,j,iPatch))
          call contravProjSphere2Plane(stat%contraU(i,j,iPatch), stat%contraV(i,j,iPatch), stat%zonal_wind(i,j,iPatch), stat%meridional_wind(i,j,iPatch), mesh%matrixIA(:,:,i,j,iPatch))
        enddo
      enddo
    enddo
    
    stat%phiG = stat%phi * mesh%sqrtG
    
    mesh%phi_s = 0.
  end subroutine case2
  
  ! Rossby-Haurwitz wave with wavenumber 4
  subroutine case6
  
  end subroutine case6
  
end module test_case_mod
    