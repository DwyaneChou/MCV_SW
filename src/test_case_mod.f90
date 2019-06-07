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
  subroutine case6(stat)
    type(stat_field), intent(inout) :: stat
    real,parameter              :: omg  = 7.848d-6         ! angular velocity of RH wave
    real,parameter              :: R    = 4d0              ! wave number of RH wave
    real,parameter              :: h0   = 8000.d0          ! wave number of RH wave
    
    real,dimension(0:nx1,0:ny1) :: u1,u2,u3                ! working array
    real,dimension(0:nx1,0:ny1) :: AA1,Ac,A21,A22,A23,Ah   ! working array
    real,dimension(0:nx1,0:ny1) :: Bc,BB1,BB2,Bh           ! working array
    real,dimension(0:nx1,0:ny1) :: CC,CC1,CC2,Ch           ! working array
    real,dimension(0:nx1,0:ny1) :: coslat                  ! working array
    
    integer                     :: i,j                     ! working variable
    
    coslat       = cos_lat
    coslat(:,1 ) = 0.d0
    coslat(:,ny) = 0.d0
    
    do j=1,ny
        do i=1,nx
            u1(i,j) = coslat(i,j)
            u2(i,j) = R*coslat(i,j)**(R-1)*sin_lat(i,j)**2*dcos(R*longitude(i,j))
            u3(i,j) = coslat(i,j)**(R+1)*dcos(R*longitude(i,j))
            u (i,j) = a*omg*(u1(i,j)+u2(i,j)-u3(i,j))
            
            v (i,j) = -a*omg*R*coslat(i,j)**(R-1)*sin_lat(i,j)*dsin(R*longitude(i,j))
            
            AA1 (i,j) = omg*0.5d0*(2.d0*omg0+omg)*coslat(i,j)**2
            Ac  (i,j) = 0.25*omg**2
            A21 (i,j) = (R+1.d0)*coslat(i,j)**(2.d0*R+2.d0)
            A22 (i,j) = (2.d0*R**2-R-2.d0)*coslat(i,j)**(2.d0*R)
            A23 (i,j) = 2.d0*R**2*coslat(i,j)**(2.d0*R-2)
            Ah  (i,j) = AA1(i,j)+Ac(i,j)*(A21(i,j)+A22(i,j)-A23(i,j))
            
            Bc  (i,j) = 2.*(omg0+omg)*omg/((R+1)*(R+2))*coslat(i,j)**R
            BB1 (i,j) = R**2+2.d0*R+2.d0
            BB2 (i,j) = (R+1.d0)**2.*coslat(i,j)**2
            Bh  (i,j) = Bc(i,j)*(BB1(i,j)-BB2(i,j))
            
            CC  (i,j) = 0.25*omg**2*coslat(i,j)**(2.d0*R)
            CC1 (i,j) = (R+1.d0)*coslat(i,j)**2;
            CC2 (i,j) = R+2.d0
            Ch  (i,j) = CC(i,j)*(CC1(i,j)-CC2(i,j))
            
            wh (i,j) = g*h0+a**2*(Ah(i,j) + Bh(i,j)*dcos(R*longitude(i,j)) + Ch(i,j)*dcos(2.0*R*longitude(i,j)))
        enddo
    enddo
    
    u (0  ,:) = u (nx,:)
    v (0  ,:) = v (nx,:)
    wh(0  ,:) = wh(nx,:)
    u (nx1,:) = u (1 ,:)
    v (nx1,:) = v (1 ,:)
    wh(nx1,:) = wh(1 ,:)
    
    u (:,0  ) = 0
    u (:,ny1) = 0
    v (:,0  ) = 0
    v (:,ny1) = 0
    wh(:,0  ) = 0
    wh(:,ny1) = 0
  end subroutine case6
  
end module test_case_mod
    