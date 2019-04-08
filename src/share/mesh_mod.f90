MODULE mesh_mod
  use constants_mod
  use parameters_mod
  use projection_mod
  
  implicit none
  
  ! coordinate
  type mesh_info
    real, dimension(:,:,:    ), allocatable :: x        ! central angle on x direction for each patch, unit: degree
    real, dimension(:,:,:    ), allocatable :: y        ! central angle on y direction for each patch, unit: degree
    real, dimension(:,:,:    ), allocatable :: xi       ! length of the arc on x direction for each patch, xi  = radius * x
    real, dimension(:,:,:    ), allocatable :: eta      ! length of the arc on y direction for each patch, eta = radius * y
    real, dimension(:,:,:    ), allocatable :: lon      ! longitude on sphere coordinate
    real, dimension(:,:,:    ), allocatable :: lat      ! latitude on sphere coordinate
    real, dimension(:,:,:    ), allocatable :: sqrtG    ! jacobian of Transformation, sqrt(G)
    real, dimension(:,:,:,:,:), allocatable :: matrixG  ! horizontal metric Tensor, which transform covariant vectors to contravariant vectors
    real, dimension(:,:,:,:,:), allocatable :: matrixIG ! horizontal metric Tensor, which transform contravariant vectors to covariant vectors
    real, dimension(:,:,:,:,:), allocatable :: matrixA  ! horizontal metric Tensor, which transform 
    real, dimension(:,:,:,:,:), allocatable :: matrixIA ! horizontal metric Tensor, which transform 
    real, dimension(:,:,:    ), allocatable :: dsqrtGdx ! \partial \sqrt(G) / \partial x
    real, dimension(:,:,:    ), allocatable :: dsqrtGdy ! \partial \sqrt(G) / \partial y
    real, dimension(:,:,:    ), allocatable :: dG11dx   ! \partial \G11     / \partial x
    real, dimension(:,:,:    ), allocatable :: dG11dy   ! \partial \G11     / \partial y
    real, dimension(:,:,:    ), allocatable :: dG12dx   ! \partial \G12     / \partial x
    real, dimension(:,:,:    ), allocatable :: dG12dy   ! \partial \G12     / \partial y
    real, dimension(:,:,:    ), allocatable :: dG22dx   ! \partial \G22     / \partial x
    real, dimension(:,:,:    ), allocatable :: dG22dy   ! \partial \G22     / \partial y
    
    real, dimension(:,:,:    ), allocatable :: sinx    ! trigonometric function
    real, dimension(:,:,:    ), allocatable :: cosx    ! trigonometric function
    real, dimension(:,:,:    ), allocatable :: tanx    ! trigonometric function
    real, dimension(:,:,:    ), allocatable :: cotx    ! trigonometric function
    real, dimension(:,:,:    ), allocatable :: secx    ! trigonometric function
    real, dimension(:,:,:    ), allocatable :: cscx    ! trigonometric function
    
    real, dimension(:,:,:    ), allocatable :: siny    ! trigonometric function
    real, dimension(:,:,:    ), allocatable :: cosy    ! trigonometric function
    real, dimension(:,:,:    ), allocatable :: tany    ! trigonometric function
    real, dimension(:,:,:    ), allocatable :: coty    ! trigonometric function
    real, dimension(:,:,:    ), allocatable :: secy    ! trigonometric function
    real, dimension(:,:,:    ), allocatable :: cscy    ! trigonometric function
  end type mesh_info
  
  ! PV index on Cell
  real   , dimension(:,:,:,:), allocatable :: pvXIndexOnCell ! PV index on Cell in x direction
  real   , dimension(:,:,:,:), allocatable :: pvYIndexOnCell ! PV index on Cell in y direction
  
  type(mesh_info) :: PVmesh
  type(mesh_info) :: VIAmesh
  
  contains
  
  subroutine initMesh
    integer :: iPV, jPV, iCell, jCell, iPatch, iVar, iDOF, jDOF
    integer :: iPVs, iPVe, jPVs, jPVe
    real    :: rho
    real    :: sinx, cosx, secx, cscx, tanx, cotx, &
               siny, cosy, secy, cscy, tany, coty
    
    ! Allocate arrays in structures
    allocate( PVmesh%x        (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( PVmesh%y        (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( PVmesh%xi       (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( PVmesh%eta      (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( PVmesh%lon      (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( PVmesh%lat      (      ips:ipe, jps:jpe, ifs:ife) )
    
    allocate( PVmesh%sqrtG    (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( PVmesh%matrixG  (2, 2, ips:ipe, jps:jpe, ifs:ife) )
    allocate( PVmesh%matrixIG (2, 2, ips:ipe, jps:jpe, ifs:ife) )
    allocate( PVmesh%matrixA  (2, 2, ips:ipe, jps:jpe, ifs:ife) )
    allocate( PVmesh%matrixIA (2, 2, ips:ipe, jps:jpe, ifs:ife) )
    
    allocate( PVmesh%dsqrtGdx (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( PVmesh%dsqrtGdy (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( PVmesh%dG11dx   (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( PVmesh%dG11dy   (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( PVmesh%dG12dx   (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( PVmesh%dG12dy   (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( PVmesh%dG22dx   (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( PVmesh%dG22dy   (      ips:ipe, jps:jpe, ifs:ife) )
    
    allocate( PVmesh%sinx     (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( PVmesh%cosx     (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( PVmesh%tanx     (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( PVmesh%cotx     (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( PVmesh%secx     (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( PVmesh%cscx     (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( PVmesh%siny     (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( PVmesh%cosy     (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( PVmesh%tany     (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( PVmesh%coty     (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( PVmesh%secy     (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( PVmesh%cscy     (      ips:ipe, jps:jpe, ifs:ife) )
    
    allocate( VIAmesh%x   (ics:ice, jcs:jce, ifs:ife) )
    allocate( VIAmesh%y   (ics:ice, jcs:jce, ifs:ife) )
    allocate( VIAmesh%xi  (ics:ice, jcs:jce, ifs:ife) )
    allocate( VIAmesh%eta (ics:ice, jcs:jce, ifs:ife) )
    allocate( VIAmesh%lon (ics:ice, jcs:jce, ifs:ife) )
    allocate( VIAmesh%lat (ics:ice, jcs:jce, ifs:ife) )
    
    allocate( pvXIndexOnCell(DOF, ics:ice, jcs:jce, ifs:ife) )
    allocate( pvYIndexOnCell(DOF, ics:ice, jcs:jce, ifs:ife) )
    
    ! Calculate the pv index on cell
    do iPatch = ifs, ife
      do jCell = jcs, jce
        do iCell = ics, ice
          do iDOF= 1, DOF
            pvXIndexOnCell(iDOF, iCell, jCell, iPatch) = (iCell - 1)*DOF + (iDOF - iCell + 1)
          enddo
          
          do jDOF = 1, DOF
            pvYIndexOnCell(jDOF, iCell, jCell, iPatch) = (jCell - 1)*DOF + (jDOF - jCell + 1)
          enddo
        end do
      end do
    end do
    
    ! Calculate mesh infomation on PV
    do iPatch = ifs, ife
      do jPV = jps, jpe
        do iPV = ips, ipe
          PVmesh%x  (iPV, jPV, iPatch) = (iPV - 1) * dx/(DOF - 1) + x_min
          PVmesh%y  (iPV, jPV, iPatch) = (jPV - 1) * dy/(DOF - 1) + y_min
          
          PVmesh%xi (iPV, jPV, iPatch) = PVmesh%x(iPV, jPV, iPatch) * radius
          PVmesh%eta(iPV, jPV, iPatch) = PVmesh%y(iPV, jPV, iPatch) * radius
          
          call pointProjPlane2Sphere(PVmesh%lon(iPV, jPV, iPatch), PVmesh%lat(iPV, jPV, iPatch), &
                                     PVmesh%x  (iPV, jPV, iPatch), PVmesh%y  (iPV, jPV, iPatch), iPatch)
          
          PVmesh%sinx(iPV, jPV, iPatch) = sin(PVmesh%x(iPV, jPV, iPatch))
          PVmesh%cosx(iPV, jPV, iPatch) = cos(PVmesh%x(iPV, jPV, iPatch))
          PVmesh%tanx(iPV, jPV, iPatch) = tan(PVmesh%x(iPV, jPV, iPatch))
          PVmesh%cotx(iPV, jPV, iPatch) = 1. / PVmesh%tanx(iPV, jPV, iPatch)
          PVmesh%secx(iPV, jPV, iPatch) = 1. / PVmesh%cosx(iPV, jPV, iPatch)
          PVmesh%cscx(iPV, jPV, iPatch) = 1. / PVmesh%sinx(iPV, jPV, iPatch)
          PVmesh%siny(iPV, jPV, iPatch) = sin(PVmesh%y(iPV, jPV, iPatch))
          PVmesh%cosy(iPV, jPV, iPatch) = cos(PVmesh%y(iPV, jPV, iPatch))
          PVmesh%tany(iPV, jPV, iPatch) = tan(PVmesh%y(iPV, jPV, iPatch))
          PVmesh%coty(iPV, jPV, iPatch) = 1. / PVmesh%tany(iPV, jPV, iPatch)
          PVmesh%secy(iPV, jPV, iPatch) = 1. / PVmesh%cosy(iPV, jPV, iPatch)
          PVmesh%cscy(iPV, jPV, iPatch) = 1. / PVmesh%siny(iPV, jPV, iPatch)
          
          call calc_matrixG (PVmesh%matrixG (:, :, iPV, jPV, iPatch), PVmesh%x  (iPV, jPV, iPatch), PVmesh%y  (iPV, jPV, iPatch))
          call calc_matrixIG(PVmesh%matrixIG(:, :, iPV, jPV, iPatch), PVmesh%x  (iPV, jPV, iPatch), PVmesh%y  (iPV, jPV, iPatch))
          call calc_matrixA (PVmesh%matrixA (:, :, iPV, jPV, iPatch), PVmesh%lon(iPV, jPV, iPatch), PVmesh%lat(iPV, jPV, iPatch), iPatch)
          call calc_matrixIA(PVmesh%matrixIA(:, :, iPV, jPV, iPatch), PVmesh%lon(iPV, jPV, iPatch), PVmesh%lat(iPV, jPV, iPatch), iPatch)
          
          sinx = PVmesh%sinx(iPV, jPV, iPatch)
          cosx = PVmesh%cosx(iPV, jPV, iPatch)
          tanx = PVmesh%tanx(iPV, jPV, iPatch)
          cotx = PVmesh%cotx(iPV, jPV, iPatch)
          secx = PVmesh%secx(iPV, jPV, iPatch)
          cscx = PVmesh%cscx(iPV, jPV, iPatch)
          siny = PVmesh%siny(iPV, jPV, iPatch)
          cosy = PVmesh%cosy(iPV, jPV, iPatch)
          tany = PVmesh%tany(iPV, jPV, iPatch)
          coty = PVmesh%coty(iPV, jPV, iPatch)
          secy = PVmesh%secy(iPV, jPV, iPatch)
          cscy = PVmesh%cscy(iPV, jPV, iPatch)
          
          PVmesh%sqrtG(iPV, jPV, iPatch) = radius**2/((1 + tanx**2 + tany**2)**(3./2.) * cosx**2 * cosy**2)
          
          PVmesh%dsqrtGdx(iPV, jPV, iPatch) = radius * secy**2 * tanx * (-3*secx**4 + 2*secx**2 * (secy**2 + tanx**2)) / (secy**2 + tanx**2)**(5/2)
          PVmesh%dsqrtGdy(iPV, jPV, iPatch) = radius * secx**2 * tany * (   secy**4 - 2*secy**2 *  tanx**2           ) / (secy**2 + tanx**2)**(5/2)
          PVmesh%dG11dx  (iPV, jPV, iPatch) = 2. * sinx * cosx * tany**2 / radius**3
          PVmesh%dG12dx  (iPV, jPV, iPatch) = siny * (cosx**2 * secy + 0.5 * (3. + 2. * cosx**2 - 1.) * cosy * tanx**2 - sinx**2 * siny * tany) / radius**3
          PVmesh%dG22dx  (iPV, jPV, iPatch) = 2. * cosy**2 * secx**2 * tanx / radius**3
          PVmesh%dG11dy  (iPV, jPV, iPatch) = 2. * cosx**2 * secy**2 * tany / radius**3
          PVmesh%dG12dy  (iPV, jPV, iPatch) = sinx * (cosx * secy**2 + (2. * cosx**2 - 1.) * sinx * tanx) / radius**3
          PVmesh%dG22dy  (iPV, jPV, iPatch) = 2. * siny * cosy * tanx**2 / radius**3
          
          !print*,PVmesh%dsqrtGdx(iPV, jPV, iPatch)/PVmesh%sqrtG (iPV, jPV, iPatch)
          !print*,PVmesh%dG11dx  (iPV, jPV, iPatch),PVmesh%dG11dy(iPV, jPV, iPatch)
          !print*,PVmesh%dG12dx  (iPV, jPV, iPatch),PVmesh%dG12dy(iPV, jPV, iPatch)
          !print*,PVmesh%dG22dx  (iPV, jPV, iPatch),PVmesh%dG22dy(iPV, jPV, iPatch)
        end do
      end do
    end do
    
    ! Calculate mesh infomation on VIA
    do iPatch = ifs, ife
      do jCell = jcs, jce
        do iCell = ics, ice
          iPVs = pvXIndexOnCell(1  , iCell, jCell, iPatch)
          iPVe = pvXIndexOnCell(DOF, iCell, jCell, iPatch)
          jPVs = pvYIndexOnCell(1  , iCell, jCell, iPatch)
          jPVe = pvYIndexOnCell(DOF, iCell, jCell, iPatch)
          
          VIAmesh%x  (iCell, jCell, iPatch) = 0.5 * (PVmesh%x(iPVs, 1   , iPatch) + PVmesh%x(iPVe, 1   , iPatch))
          VIAmesh%y  (iCell, jCell, iPatch) = 0.5 * (PVmesh%y(1   , jPVs, iPatch) + PVmesh%y(1   , jPVe, iPatch))
          
          VIAmesh%xi (iCell, jCell, iPatch) = VIAmesh%x(iCell, jCell, iPatch) * radius
          VIAmesh%eta(iCell, jCell, iPatch) = VIAmesh%y(iCell, jCell, iPatch) * radius
          
          call pointProjPlane2Sphere(VIAmesh%lon(iCell, jCell, iPatch), VIAmesh%lat(iCell, jCell, iPatch), &
                                     VIAmesh%x  (iCell, jCell, iPatch), VIAmesh%y  (iCell, jCell, iPatch), iPatch)
        end do
      end do
    end do
    
  end subroutine initMesh
  

END MODULE mesh_mod

