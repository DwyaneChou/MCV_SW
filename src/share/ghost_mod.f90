MODULE ghost_mod
  use constants_mod
  use parameters_mod
  use mesh_mod
  
  implicit none
  
  ! ghost location
  type ghostLocation
    integer, dimension(:    ), allocatable :: patchIndex ! The patch where the ghost cells locate
    integer, dimension(:    ), allocatable :: axis       ! The axis where the ghost cells locate, 1 for x axis, 2 for y axis
    real   , dimension(:,:,:), allocatable :: X          ! The x coordinate of ghost cell
    real   , dimension(:,:,:), allocatable :: Y          ! The y coordinate of ghost cell
    integer, dimension(:,:,:), allocatable :: XIndex     ! The ghost cell index on x direction
    integer, dimension(:,:,:), allocatable :: YIndex     ! The ghost cell index on y direction
  end type ghostLocation
  
  type(ghostLocation) :: ghostL ! ghost location on left edge
  type(ghostLocation) :: ghostR ! ghost location on right edge
  type(ghostLocation) :: ghostT ! ghost location on top edge
  type(ghostLocation) :: ghostB ! ghost location on bottom edge
  
  contains
  
  subroutine initGhost
    integer :: iPV, jPV, iPatch, iDOF
    real    :: lambda, theta

    allocate( ghostL%patchIndex(              ifs:ife) )
    allocate( ghostL%axis      (              ifs:ife) )
    allocate( ghostL%X         (DOF, ips:ipe, ifs:ife) )
    allocate( ghostL%Y         (DOF, ips:ipe, ifs:ife) )
    allocate( ghostL%XIndex    (DOF, ips:ipe, ifs:ife) )
    allocate( ghostL%YIndex    (DOF, ips:ipe, ifs:ife) )
    
    allocate( ghostR%patchIndex(              ifs:ife) )
    allocate( ghostR%axis      (              ifs:ife) )
    allocate( ghostR%X         (DOF, ips:ipe, ifs:ife) )
    allocate( ghostR%Y         (DOF, ips:ipe, ifs:ife) )
    allocate( ghostR%XIndex    (DOF, ips:ipe, ifs:ife) )
    allocate( ghostR%YIndex    (DOF, ips:ipe, ifs:ife) )

    allocate( ghostT%patchIndex(              ifs:ife) )
    allocate( ghostT%axis      (              ifs:ife) )
    allocate( ghostT%X         (DOF, ips:ipe, ifs:ife) )
    allocate( ghostT%Y         (DOF, ips:ipe, ifs:ife) )
    allocate( ghostT%XIndex    (DOF, ips:ipe, ifs:ife) )
    allocate( ghostT%YIndex    (DOF, ips:ipe, ifs:ife) )
    
    allocate( ghostB%patchIndex(              ifs:ife) )
    allocate( ghostB%axis      (              ifs:ife) )
    allocate( ghostB%X         (DOF, ips:ipe, ifs:ife) )
    allocate( ghostB%Y         (DOF, ips:ipe, ifs:ife) )
    allocate( ghostB%XIndex    (DOF, ips:ipe, ifs:ife) )
    allocate( ghostB%YIndex    (DOF, ips:ipe, ifs:ife) )
    
    ! Calculate ghost location info
    ! patch index
    ghostL%patchIndex(1) = 4
    ghostL%patchIndex(2) = 1
    ghostL%patchIndex(3) = 2
    ghostL%patchIndex(4) = 3
    ghostL%patchIndex(5) = 4
    ghostL%patchIndex(6) = 4
    
    ghostR%patchIndex(1) = 2
    ghostR%patchIndex(2) = 3
    ghostR%patchIndex(3) = 4
    ghostR%patchIndex(4) = 1
    ghostR%patchIndex(5) = 2
    ghostR%patchIndex(6) = 2
    
    ghostT%patchIndex(1) = 5
    ghostT%patchIndex(2) = 5
    ghostT%patchIndex(3) = 5
    ghostT%patchIndex(4) = 5
    ghostT%patchIndex(5) = 3
    ghostT%patchIndex(6) = 1
    
    ghostB%patchIndex(1) = 6
    ghostB%patchIndex(2) = 6
    ghostB%patchIndex(3) = 6
    ghostB%patchIndex(4) = 6
    ghostB%patchIndex(5) = 1
    ghostB%patchIndex(6) = 3
    
    ! axis
    ghostL%axis      (1) = 2
    ghostL%axis      (2) = 2
    ghostL%axis      (3) = 2
    ghostL%axis      (4) = 2
    ghostL%axis      (5) = 1
    ghostL%axis      (6) = 1
    
    ghostR%axis      (1) = 2
    ghostR%axis      (2) = 2
    ghostR%axis      (3) = 2
    ghostR%axis      (4) = 2
    ghostR%axis      (5) = 1
    ghostR%axis      (6) = 1
    
    ghostT%axis      (1) = 1
    ghostT%axis      (2) = 2
    ghostT%axis      (3) = 1
    ghostT%axis      (4) = 2
    ghostT%axis      (5) = 1
    ghostT%axis      (6) = 1
    
    ghostB%axis      (1) = 1
    ghostB%axis      (2) = 2
    ghostB%axis      (3) = 1
    ghostB%axis      (4) = 2
    ghostB%axis      (5) = 1
    ghostB%axis      (6) = 1
    
    ! calculate X, Y coordinate on ghost patch, here we assume ips==jps and ipe==jpe
    do iPatch = ifs, ife
      do iPV = 1, DOF - 1
        do jPV = 1, nPVx
          iDOF = iPV+ips-1
          call pointProjPlane2Sphere(lambda, theta, PVmesh%x(iDOF, jPV, iPatch), PVmesh%y(iDOF, jPV, iPatch), iPatch)
          call pointProjSphere2Plane(ghostL%X(iPV, jPV, iPatch), ghostL%Y(iPV, jPV, iPatch), lambda, theta, ghostL%patchIndex(iPatch))
          
          iDOF = iPV+nPVx
          call pointProjPlane2Sphere(lambda, theta, PVmesh%x(iDOF, jPV, iPatch), PVmesh%y(iDOF, jPV, iPatch), iPatch)
          call pointProjSphere2Plane(ghostR%X(iPV, jPV, iPatch), ghostR%Y(iPV, jPV, iPatch), lambda, theta, ghostR%patchIndex(iPatch))
          
          iDOF = iPV+nPVx
          call pointProjPlane2Sphere(lambda, theta, PVmesh%x(jPV, iDOF, iPatch), PVmesh%y(jPV, iDOF, iPatch), iPatch)
          call pointProjSphere2Plane(ghostT%X(iPV, jPV, iPatch), ghostT%Y(iPV, jPV, iPatch), lambda, theta, ghostT%patchIndex(iPatch))
          
          iDOF = iPV+ips-1
          call pointProjPlane2Sphere(lambda, theta, PVmesh%x(jPV, iDOF, iPatch), PVmesh%y(jPV, iDOF, iPatch), iPatch)
          call pointProjSphere2Plane(ghostB%X(iPV, jPV, iPatch), ghostB%Y(iPV, jPV, iPatch), lambda, theta, ghostB%patchIndex(iPatch))
        enddo
      enddo
    enddo
    
    !! Ghost Location check, the first value should be variable, and the 2nd one should be stationary
    !do iPatch = ifs, ife
    !  do iPV = 1, nPVx
    !    if(ghostB%axis(iPatch) == 1)then
    !      print*,iPatch,ghostB%X(3, iPV, iPatch), ghostB%Y(3, iPV, iPatch) 
    !    else
    !      print*,iPatch,ghostB%Y(3, iPV, iPatch), ghostB%X(3, iPV, iPatch)
    !    endif
    !  enddo
    !enddo
  end subroutine initGhost
  
END MODULE ghost_mod

