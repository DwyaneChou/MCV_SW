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
    if(case_num == 6) call case6(stat(0))
    
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
  subroutine case6(stat)
    type(stat_field), intent(inout) :: stat
    real,parameter              :: omg  = 7.848d-6         ! angular velocity of RH wave
    real,parameter              :: R    = 4              ! wave number of RH wave
    real,parameter              :: h0   = 8000.          ! wave number of RH wave
    
    real,dimension(ips:ipe,jps:jpe,ifs:ife) :: u,v,phi                 ! working array
    real,dimension(ips:ipe,jps:jpe,ifs:ife) :: u1,u2,u3                ! working array
    real,dimension(ips:ipe,jps:jpe,ifs:ife) :: AA1,Ac,A21,A22,A23,Ah   ! working array
    real,dimension(ips:ipe,jps:jpe,ifs:ife) :: Bc,BB1,BB2,Bh           ! working array
    real,dimension(ips:ipe,jps:jpe,ifs:ife) :: CC,CC1,CC2,Ch           ! working array
    real,dimension(ips:ipe,jps:jpe,ifs:ife) :: sinlat                  ! working array
    real,dimension(ips:ipe,jps:jpe,ifs:ife) :: coslat                  ! working array
    real,dimension(ips:ipe,jps:jpe,ifs:ife) :: longitude               ! working array
    
    integer                                 :: i,j,iPatch              ! working variable
    
    sinlat    = sin(mesh%latP)
    coslat    = cos(mesh%latP)
    longitude = mesh%lonP
    
    do iPatch = ifs, ife
      do j = jps, jpe
        do i = ips, ipe
          u1(i,j,iPatch) = coslat(i,j,iPatch)
          u2(i,j,iPatch) = R*coslat(i,j,iPatch)**(R-1)*sinlat(i,j,iPatch)**2*cos(R*longitude(i,j,iPatch))
          u3(i,j,iPatch) = coslat(i,j,iPatch)**(R+1)*cos(R*longitude(i,j,iPatch))
          u (i,j,iPatch) = radius*omg*(u1(i,j,iPatch)+u2(i,j,iPatch)-u3(i,j,iPatch))
          
          v (i,j,iPatch) = -radius*omg*R*coslat(i,j,iPatch)**(R-1)*sinlat(i,j,iPatch)*sin(R*longitude(i,j,iPatch))
          
          AA1 (i,j,iPatch) = omg*0.5*(2.*Omega+omg)*coslat(i,j,iPatch)**2
          Ac  (i,j,iPatch) = 0.25*omg**2
          A21 (i,j,iPatch) = (R+1.)*coslat(i,j,iPatch)**(2.*R+2.)
          A22 (i,j,iPatch) = (2.*R**2-R-2.)*coslat(i,j,iPatch)**(2.*R)
          A23 (i,j,iPatch) = 2.*R**2*coslat(i,j,iPatch)**(2.*R-2)
          Ah  (i,j,iPatch) = AA1(i,j,iPatch)+Ac(i,j,iPatch)*(A21(i,j,iPatch)+A22(i,j,iPatch)-A23(i,j,iPatch))
          
          Bc  (i,j,iPatch) = 2.*(Omega+omg)*omg/((R+1)*(R+2))*coslat(i,j,iPatch)**R
          BB1 (i,j,iPatch) = R**2+2.*R+2.
          BB2 (i,j,iPatch) = (R+1.)**2.*coslat(i,j,iPatch)**2.
          Bh  (i,j,iPatch) = Bc(i,j,iPatch)*(BB1(i,j,iPatch)-BB2(i,j,iPatch))
          
          CC  (i,j,iPatch) = 0.25*omg**2*coslat(i,j,iPatch)**(2.*R)
          CC1 (i,j,iPatch) = (R+1.)*coslat(i,j,iPatch)**2;
          CC2 (i,j,iPatch) = R+2.
          Ch  (i,j,iPatch) = CC(i,j,iPatch)*(CC1(i,j,iPatch)-CC2(i,j,iPatch))
          
          phi (i,j,iPatch) = gravity*h0+radius**2*(Ah(i,j,iPatch) + Bh(i,j,iPatch)*cos(R*longitude(i,j,iPatch)) + Ch(i,j,iPatch)*cos(2.0*R*longitude(i,j,iPatch)))
        enddo
      enddo
    enddo
    
    stat%zonal_wind      = u
    stat%meridional_wind = v
    stat%phi             = phi
    stat%phiG            = stat%phi * mesh%sqrtG
    
    do iPatch = ifs, ife
      do j = jps, jpe
        do i = ips, ipe
          call covProjSphere2Plane    (stat%u      (i,j,iPatch), stat%v      (i,j,iPatch), stat%zonal_wind(i,j,iPatch), stat%meridional_wind(i,j,iPatch), mesh%matrixIA(:,:,i,j,iPatch), mesh%matrixG(:,:,i,j,iPatch))
          call contravProjSphere2Plane(stat%contraU(i,j,iPatch), stat%contraV(i,j,iPatch), stat%zonal_wind(i,j,iPatch), stat%meridional_wind(i,j,iPatch), mesh%matrixIA(:,:,i,j,iPatch))
        enddo
      enddo
    enddo
    
    mesh%phi_s = 0.
    
    !u  (0  ,:) = u (nx,:)
    !v  (0  ,:) = v (nx,:)
    !phi(0  ,:) = phi(nx,:)
    !u  (nx1,:) = u (1 ,:)
    !v  (nx1,:) = v (1 ,:)
    !phi(nx1,:) = phi(1 ,:)
    !
    !u  (:,0  ) = 0
    !u  (:,ny1) = 0
    !v  (:,0  ) = 0
    !v  (:,ny1) = 0
    !phi(:,0  ) = 0
    !phi(:,ny1) = 0
  end subroutine case6
  
end module test_case_mod
    