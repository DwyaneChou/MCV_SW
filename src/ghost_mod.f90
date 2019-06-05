MODULE ghost_mod
  use constants_mod
  use parameters_mod
  use projection_mod
  use mesh_mod
  use stat_mod
  
  implicit none
  
  integer nPVHalo
  
  ! ghost location
  type ghostLocation
    real   , dimension(:,:), allocatable :: X     ! X coordinate of ghost points
    real   , dimension(:,:), allocatable :: Y     ! Y coordinate of ghost points, not needed in interpolation
    real   , dimension(:,:), allocatable :: coefL ! Left coefficient of linear interpolation
    real   , dimension(:,:), allocatable :: coefR ! Right coefficient of linear interpolation
    integer, dimension(:,:), allocatable :: iref  ! Index of reference point
  end type ghostLocation
  
  type(ghostLocation) :: ghost
  contains
  
  subroutine initGhost
    integer :: i,j
    integer :: iPatch ! computational patch
    integer :: gpatch ! ghost patch
    integer :: ihalo  ! halo index
    real    :: lambda, theta
    real    :: coefL
    real    :: coefR
    integer :: iref
    
    real    :: X_RAW(ids:ide)
    
    nPVHalo = xhalo * (DOF - 1)
    
    allocate(ghost%X    (ids:ide,nPVHalo))
    allocate(ghost%Y    (ids:ide,nPVHalo))
    allocate(ghost%coefL(ids:ide,nPVHalo))
    allocate(ghost%coefR(ids:ide,nPVHalo))
    allocate(ghost%iref (ids:ide,nPVHalo))
    
    ! Calculate ghost points location
    iPatch = 1
    gpatch = 5
    do j = jde+1, jpe
      do i = ids, ide
        ihalo = j - jde
        call pointProjPlane2Sphere(lambda          , theta           , mesh%xP(i,j,iPatch), mesh%yP(i,j,iPatch), iPatch)
        call pointProjSphere2Plane(ghost%X(i,ihalo), ghost%Y(i,ihalo), lambda             , theta              , gPatch)
      enddo
    enddo
    
    X_RAW = mesh%xP(ids:ide,1,iPatch)
    
    ! Calculate reference points
    ! Now apply linear interpolation to obtain edge components
    DO j = 1, nPVHalo
      ! Reset the reference index
      iref = ids
      
      DO i = ids, ide
        !
        ! Find reference points
        !
        DO WHILE ((iref .NE. ide) .AND. (ghost%X(i,j) > X_RAW(iref)))
          iref = iref + 1
        ENDDO
        ghost%iref(i,j) = iref
    
        IF ((ghost%X(i,j) >= X_RAW(iref-1)) .AND. (ghost%X(i,j) <= X_RAW(iref  )))THEN
            
          coefR = (ghost%X(i,j) - X_RAW(iref-1)) / (X_RAW(iref) - X_RAW(iref-1))
          coefL = 1. - coefR
          
          ghost%coefL(i,j) = coefL
          ghost%coefR(i,j) = coefR
          
          !print*,i,iref,ghost%coefL(i,j),ghost%coefR(i,j)
          
          IF ((coefR < 0.0) .OR. (coefR > 1.)) THEN
            WRITE (*,*) 'FAIL in CubedSphereFillHalo_Linear'
            WRITE (*,*) 'a out of bounds'
            STOP
          ENDIF
        ENDIF
      ENDDO
    ENDDO
    
  end subroutine initGhost
  
  !------------------------------------------------------------------------------
  ! SUBROUTINE CubedSphereFillHalo
  !
  ! Description:
  !   Recompute the cubed sphere data storage array, with the addition of a
  !   halo region around the specified panel.
  !------------------------------------------------------------------------------
  SUBROUTINE CubedSphereFillHalo(field)

    IMPLICIT NONE

    REAL, DIMENSION(ips:ipe,jps:jpe,ifs:ife), INTENT(INOUT) :: field
    
    REAL, DIMENSION(ids:ide,jds:jde,ifs:ife) :: field_inner

    ! Local variables
    INTEGER :: i
    
    !zarg = 0.0 !DBG
    field_inner = field(ids:ide,jds:jde,ifs:ife)

    ! Equatorial panels
    !IF (np==1) THEN
       DO i=1,nPVHalo
          field(1-i    ,jds:jde,1) = field_inner(ide-i  ,jds:jde,4)  !exchange left
          field(ide+i  ,jds:jde,1) = field_inner(i+1    ,jds:jde,2)  !exchange right
          field(ids:ide,1-i    ,1) = field_inner(ids:ide,jde-i  ,6)  !exchange below
          field(ids:ide,jde+i  ,1) = field_inner(ids:ide,i+1    ,5)  !exchange over
       ENDDO
    !ELSE IF (np==2) THEN
       DO i=1,nPVHalo
          field(1-i    ,jds:jde,2) = field_inner(ide-i  ,jds:jde   ,1)  !exchange left
          field(ide+i  ,jds:jde,2) = field_inner(i+1    ,jds:jde   ,3)  !exchange right
          field(ids:ide,1-i    ,2) = field_inner(ide-i  ,jde:jds:-1,6)  !exchange below
          field(ids:ide,jde+i  ,2) = field_inner(ide-i  ,jds:jde   ,5)  !exchange over
       ENDDO
    !ELSE IF (np==3) THEN
       DO i=1,nPVHalo
          field(1-i    ,jds:jde,3) = field_inner(ide-i     ,jds:jde,2)  !exchange left
          field(ide+i  ,jds:jde,3) = field_inner(i+1       ,jds:jde,4)  !exchange right
          field(ids:ide,1-i    ,3) = field_inner(ide:ids:-1,i+1    ,6)  !exchange below
          field(ids:ide,jde+i  ,3) = field_inner(ide:ids:-1,jde-i  ,5)  !exchange over
       ENDDO
    !ELSE IF (np==4) THEN
       DO i=1,nPVHalo
          field(1-i    ,jds:jde,4) = field_inner(ide-i     ,jds:jde   ,3) !exchange left
          field(ide+i  ,jds:jde,4) = field_inner(i+1       ,jds:jde   ,1) !exchange right
          field(ids:ide,1-i    ,4) = field_inner(i+1       ,jds:jde   ,6) !exchange below
          field(ids:ide,jde+i  ,4) = field_inner(i+1       ,jde:jds:-1,5) !exchange over
       ENDDO
    ! Top panel
    !ELSE IF (np==5) THEN
       DO i=1,nPVHalo
          field(1-i    ,jds:jde,5) = field_inner(ide:ids:-1,jde-i,4) !exchange left
          field(ide+i  ,jds:jde,5) = field_inner(ids:ide   ,jde-i,2) !exchange right
          field(ids:ide,1-i    ,5) = field_inner(ids:ide   ,jde-i,1) !exchange below
          field(ids:ide,jde+i  ,5) = field_inner(ide:ids:-1,jde-i,3) !exchange over
       ENDDO
    ! Bottom panel
    !ELSE IF (np==6) THEN
       DO i=1,nPVHalo
          field(1-i    ,jds:jde,6) = field_inner(ids:ide   ,i+1  ,4) !exchange left
          field(ide+i  ,jds:jde,6) = field_inner(ide:ids:-1,i+1  ,2) !exchange right
          field(ids:ide,1-i    ,6) = field_inner(ide:ids:-1,i+1  ,3) !exchange below
          field(ids:ide,jde+i  ,6) = field_inner(ids:ide   ,i+1  ,1) !exchange over
       ENDDO
    !ELSE
    !   WRITE (*,*) 'Fatal error: In CubedSphereFillHalo'
    !   WRITE (*,*) 'Invalid panel id ', np
    !   STOP
    !ENDIF
  END SUBROUTINE CubedSphereFillHalo
  
  subroutine CubedSphereFillGhost(field)
    real, intent(inout) :: field(ips:ipe,jps:jpe,ifs:ife)
    
    real ghost_target(ips:ipe,jps:jpe,ifs:ife)
    
    integer i,j,iPatch
    
    ghost_target                          = FillValue
    ghost_target(ids:ide,jds:jde,ifs:ife) = field(ids:ide,jds:jde,ifs:ife)
    
    call CubedSphereFillHalo(field)
    
    do iPatch = ifs, ife
      do j = 1, nPVHalo
        ! Left boundary
        call linear_interp(ghost_target(ids-j,jds:jde,iPatch),field(ids-j,jds:jde,iPatch),ghost%iref(:,j),ghost%coefL(:,j),ghost%coefR(:,j))
        ! Rigth boundary
        call linear_interp(ghost_target(ide+j,jds:jde,iPatch),field(ide+j,jds:jde,iPatch),ghost%iref(:,j),ghost%coefL(:,j),ghost%coefR(:,j))
        ! Top boundary
        call linear_interp(ghost_target(ids:ide,jde+j,iPatch),field(ids:ide,jde+j,iPatch),ghost%iref(:,j),ghost%coefL(:,j),ghost%coefR(:,j))
        ! Bottom boundary
        call linear_interp(ghost_target(ids:ide,jds-j,iPatch),field(ids:ide,jds-j,iPatch),ghost%iref(:,j),ghost%coefL(:,j),ghost%coefR(:,j))
      enddo
    enddo
    
    field = ghost_target
  
  end subroutine CubedSphereFillGhost
  
  ! Fill up halo with ghost points  
  subroutine fill_ghost(stat)
    type(stat_field), intent(inout) :: stat

    call CubedSphereFillGhost(stat%zonal_wind     )
    call CubedSphereFillGhost(stat%meridional_wind)
    call CubedSphereFillGhost(stat%phi            )
  end subroutine fill_ghost
  
  subroutine linear_interp(dest,src,iref,coefL,coefR)
    real   , intent(out) :: dest (ids:ide)
    real   , intent(in ) :: src  (ids:ide)
    integer, intent(in ) :: iref (ids:ide)
    real   , intent(in ) :: coefL(ids:ide)
    real   , intent(in ) :: coefR(ids:ide)
    
    dest = coefL * src(iref-1) + coefR * src(iref)
    
  end subroutine linear_interp
  
END MODULE ghost_mod