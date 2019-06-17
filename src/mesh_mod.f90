MODULE mesh_mod
  use constants_mod
  use parameters_mod
  use projection_mod
  
  implicit none
  
  ! coordinate
  type mesh_info
    real, dimension(:,:,:    ), allocatable :: xP        ! central angle on x direction for points on each patch, unit: radian
    real, dimension(:,:,:    ), allocatable :: yP        ! central angle on y direction for points on each patch, unit: radian
    real, dimension(:,:,:    ), allocatable :: xC        ! central angle on x direction for cells on each patch, unit: radian
    real, dimension(:,:,:    ), allocatable :: yC        ! central angle on y direction for cells on each patch, unit: radian
    real, dimension(:,:,:    ), allocatable :: lonP      ! longitude on points
    real, dimension(:,:,:    ), allocatable :: latP      ! latitude on points
    real, dimension(:,:,:    ), allocatable :: lonC      ! longitude on cells
    real, dimension(:,:,:    ), allocatable :: latC      ! latitude on cells
    real, dimension(:,:,:    ), allocatable :: sqrtG    ! jacobian of Transformation, sqrt(G)
    real, dimension(:,:,:,:,:), allocatable :: matrixG  ! horizontal metric Tensor, which transform covariant vectors to contravariant vectors
    real, dimension(:,:,:,:,:), allocatable :: matrixIG ! horizontal metric Tensor, which transform contravariant vectors to covariant vectors
    real, dimension(:,:,:,:,:), allocatable :: matrixA  ! horizontal metric Tensor, which transform 
    real, dimension(:,:,:,:,:), allocatable :: matrixIA ! horizontal metric Tensor, which transform
    
    real, dimension(:,:,:    ), allocatable :: f       ! Coriolis parameter
    
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
    
    real, dimension(:,:,:    ), allocatable :: phi_s   ! surface geopotential height
    
    real, dimension(:,:      ), allocatable :: areaCell
    real                                    :: weightsOnPV(DOF,DOF)

  end type mesh_info
  
  ! PV index on Cell
  integer, dimension(:,:    ), allocatable :: pvIdx          ! PV index on Cell
  integer, dimension(:,:    ), allocatable :: pvIdy          ! PV index on Cell
  integer, dimension(:,:,:,:), allocatable :: pvXIndexOnCell ! PV index on Cell in x direction
  integer, dimension(:,:,:,:), allocatable :: pvYIndexOnCell ! PV index on Cell in y direction
  
  type(mesh_info) :: mesh
  
  contains
  
  subroutine initMesh
    integer :: iPV, jPV, iCell, jCell, iPatch, iVar, iDOF, jDOF
    integer :: iPVs, iPVe, jPVs, jPVe
    real    :: rho
    
    ! Allocate arrays in structures
    allocate( mesh%xP       (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( mesh%yP       (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( mesh%lonP     (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( mesh%latP     (      ips:ipe, jps:jpe, ifs:ife) )
    
    allocate( mesh%xC       (      ics:ice, jcs:jce, ifs:ife) )
    allocate( mesh%yC       (      ics:ice, jcs:jce, ifs:ife) )
    allocate( mesh%lonC     (      ics:ice, jcs:jce, ifs:ife) )
    allocate( mesh%latC     (      ics:ice, jcs:jce, ifs:ife) )
    
    allocate( mesh%sqrtG    (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( mesh%matrixG  (2, 2, ips:ipe, jps:jpe, ifs:ife) )
    allocate( mesh%matrixIG (2, 2, ips:ipe, jps:jpe, ifs:ife) )
    allocate( mesh%matrixA  (2, 2, ips:ipe, jps:jpe, ifs:ife) )
    allocate( mesh%matrixIA (2, 2, ips:ipe, jps:jpe, ifs:ife) )
    
    allocate( mesh%f        (      ips:ipe, jps:jpe, ifs:ife) )
    
    allocate( mesh%sinx     (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( mesh%cosx     (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( mesh%tanx     (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( mesh%cotx     (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( mesh%secx     (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( mesh%cscx     (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( mesh%siny     (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( mesh%cosy     (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( mesh%tany     (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( mesh%coty     (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( mesh%secy     (      ips:ipe, jps:jpe, ifs:ife) )
    allocate( mesh%cscy     (      ips:ipe, jps:jpe, ifs:ife) )
    
    allocate( mesh%phi_s    (      ips:ipe, jps:jpe, ifs:ife) )
    
    allocate( mesh%areaCell (Nx, Ny) )
    
    allocate( pvIdx         (DOF, ics:ice                  ) )
    allocate( pvIdy         (DOF, ics:ice                  ) )
    allocate( pvXIndexOnCell(DOF, ics:ice, jcs:jce, ifs:ife) )
    allocate( pvYIndexOnCell(DOF, ics:ice, jcs:jce, ifs:ife) )
    
    ! Calculate the pv index on cell
    do iCell = ics, ice
      do iDOF = 1, DOF
        pvIdx(iDOF,iCell) = (iCell - 1)*DOF + (iDOF - iCell + 1)
      enddo
    enddo
    
    do jCell = jcs, jce
      do jDOF = 1, DOF
        pvIdy(jDOF,jCell) = (jCell - 1)*DOF + (jDOF - jCell + 1)
      enddo
    enddo
    
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
          mesh%xP(iPV, jPV, iPatch) = (iPV - 1) * dx/(DOF - 1) + x_min
          mesh%yP(iPV, jPV, iPatch) = (jPV - 1) * dy/(DOF - 1) + y_min
          
          call pointProjPlane2Sphere(mesh%lonP(iPV, jPV, iPatch), mesh%latP(iPV, jPV, iPatch), &
                                     mesh%xP  (iPV, jPV, iPatch), mesh%yP  (iPV, jPV, iPatch), iPatch)
          
          mesh%sinx(iPV, jPV, iPatch) = sin(mesh%xP(iPV, jPV, iPatch))
          mesh%cosx(iPV, jPV, iPatch) = cos(mesh%xP(iPV, jPV, iPatch))
          mesh%tanx(iPV, jPV, iPatch) = tan(mesh%xP(iPV, jPV, iPatch))
          mesh%cotx(iPV, jPV, iPatch) = 1. / mesh%tanx(iPV, jPV, iPatch)
          mesh%secx(iPV, jPV, iPatch) = 1. / mesh%cosx(iPV, jPV, iPatch)
          mesh%cscx(iPV, jPV, iPatch) = 1. / mesh%sinx(iPV, jPV, iPatch)
          mesh%siny(iPV, jPV, iPatch) = sin(mesh%yP(iPV, jPV, iPatch))
          mesh%cosy(iPV, jPV, iPatch) = cos(mesh%yP(iPV, jPV, iPatch))
          mesh%tany(iPV, jPV, iPatch) = tan(mesh%yP(iPV, jPV, iPatch))
          mesh%coty(iPV, jPV, iPatch) = 1. / mesh%tany(iPV, jPV, iPatch)
          mesh%secy(iPV, jPV, iPatch) = 1. / mesh%cosy(iPV, jPV, iPatch)
          mesh%cscy(iPV, jPV, iPatch) = 1. / mesh%siny(iPV, jPV, iPatch)
          
          call calc_matrixG (mesh%matrixG (:, :, iPV, jPV, iPatch), mesh%xP  (iPV, jPV, iPatch), mesh%yP  (iPV, jPV, iPatch))
          call calc_matrixIG(mesh%matrixIG(:, :, iPV, jPV, iPatch), mesh%xP  (iPV, jPV, iPatch), mesh%yP  (iPV, jPV, iPatch))
          call calc_matrixA (mesh%matrixA (:, :, iPV, jPV, iPatch), mesh%lonP(iPV, jPV, iPatch), mesh%latP(iPV, jPV, iPatch), iPatch)
          call calc_matrixIA(mesh%matrixIA(:, :, iPV, jPV, iPatch), mesh%lonP(iPV, jPV, iPatch), mesh%latP(iPV, jPV, iPatch), iPatch)
          call calc_Jacobian(mesh%sqrtG   (      iPV, jPV, iPatch), mesh%xP  (iPV, jPV, iPatch), mesh%yP  (iPV, jPV, iPatch))
          
          mesh%f(iPV, jPV, iPatch) = 2. * Omega * sin(mesh%latP(iPV,jPV,iPatch))
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
          
          mesh%xC  (iCell, jCell, iPatch) = 0.5 * (mesh%xP(iPVs, 1   , iPatch) + mesh%xP(iPVe, 1   , iPatch))
          mesh%yC  (iCell, jCell, iPatch) = 0.5 * (mesh%yP(1   , jPVs, iPatch) + mesh%yP(1   , jPVe, iPatch))
          
          call pointProjPlane2Sphere(mesh%lonC(iCell, jCell, iPatch), mesh%latC(iCell, jCell, iPatch), &
                                     mesh%xC  (iCell, jCell, iPatch), mesh%yC  (iCell, jCell, iPatch), iPatch)
        end do
      end do
    end do
    
#ifdef MCV3
      ! Calculate weights of points in a cell For MCV3 only
      mesh%weightsOnPV(1,1) = 1.
      mesh%weightsOnPV(1,2) = 4.
      mesh%weightsOnPV(1,3) = 1.
      mesh%weightsOnPV(2,:) = 4. * mesh%weightsOnPV(1,:)
      mesh%weightsOnPV(3,:) = 1. * mesh%weightsOnPV(1,:)
      mesh%weightsOnPV      = mesh%weightsOnPV / 36.
#endif

#ifdef MCV4
      ! Calculate weights of points in a cell For MCV4 only
      mesh%weightsOnPV(1,1) = 1.
      mesh%weightsOnPV(1,2) = 3.
      mesh%weightsOnPV(1,3) = 3.
      mesh%weightsOnPV(1,4) = 1.
      mesh%weightsOnPV(2,:) = 3. * mesh%weightsOnPV(1,:)
      mesh%weightsOnPV(3,:) = 3. * mesh%weightsOnPV(1,:)
      mesh%weightsOnPV(4,:) = mesh%weightsOnPV(1,:)
      mesh%weightsOnPV      = mesh%weightsOnPV / 80.
#endif
    
#ifdef CUBE
    ! Calculate areaCell
    call EquiangularAllAreas(Nx, mesh%areaCell)
    mesh%areaCell = mesh%areaCell * radius **2
#endif

  end subroutine initMesh
  
#ifdef CUBE
  !------------------------------------------------------------------------------
  ! SUBROUTINE EquiangularAllAreas
  !
  ! Description:
  !   Compute the area of all cubed sphere grid cells, storing the results in
  !   a two dimensional array.
  !
  ! Parameters: 
  !   icube - Cell number (Nx or Ny) of the cubed sphere
  !   dA (OUT) - Output array containing the area of all cubed sphere grid cells
  !------------------------------------------------------------------------------
  SUBROUTINE EquiangularAllAreas(icube, dA)
    IMPLICIT NONE
    
    INTEGER,                         INTENT(IN)  :: icube
    REAL   , DIMENSION(icube,icube), INTENT(OUT) :: dA
    
    ! Local variables
    INTEGER                           :: k, k1, k2
    REAL                              :: a1, a2, a3, a4
    REAL , DIMENSION(icube+1,icube+1) :: ang
    REAL , DIMENSION(icube+1)         :: gp
    
    !#ifdef DBG 
    REAL    :: dbg !DBG
    !#endif
    
    ! Recall that we are using equi-angular spherical gridding
    !   Compute the angle between equiangular cubed sphere projection grid lines.
    DO k = 1, icube+1
      gp(k) = -0.25 * pi + (pi/DBLE(2*(icube))) * DBLE(k-1)
    ENDDO
    
    DO k2=1,icube+1
      DO k1=1,icube+1
        ang(k1,k2) = ACOS(-SIN(gp(k1)) * SIN(gp(k2)))
      ENDDO
    ENDDO
    
    DO k2=1,icube
      DO k1=1,icube
        a1 =      ang(k1  , k2  )
        a2 = pi - ang(k1+1, k2  )
        a3 = pi - ang(k1  , k2+1)
        a4 =      ang(k1+1, k2+1)      
        ! area = r*r*(-2*pi+sum(interior angles))
        DA(k1,k2) = -2.0*pi+a1+a2+a3+a4
      ENDDO
    ENDDO
    
    ! Only for debugging - test consistency
    dbg = 0.0
    DO k2=1,icube
      DO k1=1,icube
        dbg = dbg + DA(k1,k2)
      ENDDO
    ENDDO
    
    !print*,''
    !print*,'total area error     : ', dbg - 4. * pi / 6. !DBG
  END SUBROUTINE EquiangularAllAreas
#endif
END MODULE mesh_mod

