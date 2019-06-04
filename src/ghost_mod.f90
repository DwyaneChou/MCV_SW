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
    do i = ids,ide
      print*,i,X_RAW(i),ghost%X(i,1),ghost%X(i,2)
    enddo
    
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
      stop
    ENDDO
    
  end subroutine initGhost
  
  !------------------------------------------------------------------------------
  ! SUBROUTINE CubedSphereFillHalo
  !
  ! Description:
  !   Recompute the cubed sphere data storage array, with the addition of a
  !   halo region around the specified panel.
  !
  ! Parameters:
  !   parg - Current panel values
  !   zarg (OUT) - Calculated panel values with halo/ghost region
  !   np - Panel number
  !   ncube - Dimension of the cubed sphere (# of grid lines)
  !   nhalo - Number of halo/ghost elements around each panel
  !------------------------------------------------------------------------------
  SUBROUTINE CubedSphereFillHalo(parg, zarg, np, ncube, nhalo)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: np, ncube,nhalo

    REAL, DIMENSION(ncube-1              , ncube-1              , 6), INTENT(IN ) :: parg
    REAL, DIMENSION(1-nhalo:ncube+nhalo-1, 1-nhalo:ncube+nhalo-1, 6), INTENT(OUT) :: zarg

    ! Local variables
    INTEGER                :: jh
    
    !zarg = 0.0 !DBG
    zarg(1:ncube-1          ,1:ncube-1          ,np) = parg(1:ncube-1,1:ncube-1,np)
    zarg(1-nhalo:0          ,1-nhalo:0          ,np) = 0.0
    zarg(1-nhalo:0          ,ncube:ncube+nhalo-1,np) = 0.0
    zarg(ncube:ncube+nhalo-1,1-nhalo:0          ,np) = 0.0
    zarg(ncube:ncube+nhalo-1,ncube:ncube+nhalo-1,np) = 0.0

    ! Equatorial panels
    IF (np==1) THEN
       DO jh=1,nhalo
          zarg(ncube+jh-1,1:ncube-1 ,1) = parg(jh       ,1:ncube-1 ,2)  !exchange right
          zarg(1-jh      ,1:ncube-1 ,1) = parg(ncube-jh ,1:ncube-1 ,4)  !exchange left
          zarg(1:ncube-1 ,1-jh      ,1) = parg(1:ncube-1,ncube-jh  ,6)  !exchange below
          zarg(1:ncube-1 ,ncube+jh-1,1) = parg(1:ncube-1,jh        ,5)  !exchange over
       ENDDO
    ELSE IF (np==2) THEN
       DO jh=1,nhalo
          zarg(1-jh      ,1:ncube-1 ,2) = parg(ncube-jh,1:ncube-1   ,1)  !exchange left
          zarg(ncube+jh-1,1:ncube-1 ,2) = parg(jh      ,1:ncube-1   ,3)  !exchange right
          zarg(1:ncube-1 ,1-jh      ,2) = parg(ncube-jh,ncube-1:1:-1,6)  !exchange below
          zarg(1:ncube-1 ,ncube+jh-1,2) = parg(ncube-jh,1:ncube-1   ,5)  !exchange over
       ENDDO
    ELSE IF (np==3) THEN
       DO jh=1,nhalo
          zarg(ncube+jh-1,1:ncube-1 ,3) = parg(jh          ,1:ncube-1,4)  !exchange right
          zarg(1-jh      ,1:ncube-1 ,3) = parg(ncube-jh    ,1:ncube-1,2)  !exchange left
          zarg(1:ncube-1 ,1-jh      ,3) = parg(ncube-1:1:-1,jh       ,6)  !exchange below
          zarg(1:ncube-1 ,ncube+jh-1,3) = parg(ncube-1:1:-1,ncube-jh ,5)  !exchange over
       ENDDO
    ELSE IF (np==4) THEN
       DO jh=1,nhalo
          zarg(1-jh      ,1:ncube-1 ,4) = parg(ncube-jh,1:ncube-1   ,3) !exchange left
          zarg(ncube+jh-1,1:ncube-1 ,4) = parg(jh      ,1:ncube-1   ,1) !exchange right
          zarg(1:ncube-1 ,1-jh      ,4) = parg(jh      ,1:ncube-1   ,6) !exchange below
          zarg(1:ncube-1 ,ncube+jh-1,4) = parg(jh      ,ncube-1:1:-1,5) !exchange over
       ENDDO
    ! Top panel
    ELSE IF (np==5) THEN
       DO jh=1,nhalo
          zarg(1-jh      ,1:ncube-1 ,5) = parg(ncube-1:1:-1,ncube-jh,4) !exchange left
          zarg(ncube+jh-1,1:ncube-1 ,5) = parg(1:ncube-1   ,ncube-jh,2) !exchange right
          zarg(1:ncube-1 ,1-jh      ,5) = parg(1:ncube-1   ,ncube-jh,1) !exchange below
          zarg(1:ncube-1 ,ncube+jh-1,5) = parg(ncube-1:1:-1,ncube-jh,3) !exchange over
       ENDDO
    ! Bottom panel
    ELSE IF (np==6) THEN
       DO jh=1,nhalo
          zarg(1-jh      ,1:ncube-1 ,6) = parg(1:ncube-1   ,jh      ,4) !exchange left
          zarg(ncube+jh-1,1:ncube-1 ,6) = parg(ncube-1:1:-1,jh      ,2) !exchange right
          zarg(1:ncube-1 ,1-jh      ,6) = parg(ncube-1:1:-1,jh      ,3) !exchange below
          zarg(1:ncube-1 ,ncube+jh-1,6) = parg(1:ncube-1   ,jh      ,1) !exchange over
       ENDDO
    ELSE
       WRITE (*,*) 'Fatal error: In CubedSphereFillHalo'
       WRITE (*,*) 'Invalid panel id ', np
       STOP
    ENDIF
  END SUBROUTINE CubedSphereFillHalo

  subroutine CubedSphereFillGhost
  
  end subroutine CubedSphereFillGhost
  
  ! Fill up halo with ghost points  
  subroutine fill_ghost(stat)
    type(stat_field), intent(inout) :: stat
    
    integer iPatch

    do iPatch = ifs, ife
      !call CubedSphereFillHalo_Linear_extended(stat%zonal_wind     (ids:ide,jds:jde,:), stat%zonal_wind     (:,:,iPatch), iPatch, ide+1, xhalo*(DOF-1))
      !call CubedSphereFillHalo_Linear_extended(stat%meridional_wind(ids:ide,jds:jde,:), stat%meridional_wind(:,:,iPatch), iPatch, ide+1, xhalo*(DOF-1))
      !call CubedSphereFillHalo_Linear_extended(stat%phi            (ids:ide,jds:jde,:), stat%phi            (:,:,iPatch), iPatch, ide+1, xhalo*(DOF-1))
    enddo
  end subroutine fill_ghost
  
END MODULE ghost_mod