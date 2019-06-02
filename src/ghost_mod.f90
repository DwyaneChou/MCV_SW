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
          call pointProjPlane2Sphere(lambda, theta, mesh%xP(iDOF, jPV, iPatch), mesh%yP(iDOF, jPV, iPatch), iPatch)
          call pointProjSphere2Plane(ghostL%X(iPV, jPV, iPatch), ghostL%Y(iPV, jPV, iPatch), lambda, theta, ghostL%patchIndex(iPatch))
          
          iDOF = iPV+nPVx
          call pointProjPlane2Sphere(lambda, theta, mesh%xP(iDOF, jPV, iPatch), mesh%yP(iDOF, jPV, iPatch), iPatch)
          call pointProjSphere2Plane(ghostR%X(iPV, jPV, iPatch), ghostR%Y(iPV, jPV, iPatch), lambda, theta, ghostR%patchIndex(iPatch))
          
          iDOF = iPV+nPVx
          call pointProjPlane2Sphere(lambda, theta, mesh%xP(jPV, iDOF, iPatch), mesh%yP(jPV, iDOF, iPatch), iPatch)
          call pointProjSphere2Plane(ghostT%X(iPV, jPV, iPatch), ghostT%Y(iPV, jPV, iPatch), lambda, theta, ghostT%patchIndex(iPatch))
          
          iDOF = iPV+ips-1
          call pointProjPlane2Sphere(lambda, theta, mesh%xP(jPV, iDOF, iPatch), mesh%yP(jPV, iDOF, iPatch), iPatch)
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

  !------------------------------------------------------------------------------
  ! SUBROUTINE CubedSphereFillHalo_Linear
  !
  ! Description:
  !   Recompute the cubed sphere data storage array, with the addition of a
  !   2-element halo region around the specified panel.  Use linear order
  !   interpolation to translate between panels.
  !
  ! Parameters:
  !   parg - Current panel values
  !   zarg (OUT) - Calculated panel values with halo/ghost region
  !   np - Panel number
  !   ncube - Dimension of the cubed sphere (# of grid lines)
  !------------------------------------------------------------------------------
  SUBROUTINE CubedSphereFillHalo_Linear(parg, zarg, np, ncube)

    IMPLICIT NONE

    INTEGER, PARAMETER :: nhalo = 2

    INTEGER, INTENT(IN) :: np, ncube

    REAL, DIMENSION(ncube-1              , ncube-1              , 6), INTENT(IN ) :: parg
    REAL, DIMENSION(1-nhalo:ncube+nhalo-1, 1-nhalo:ncube+nhalo-1, 6), INTENT(OUT) :: zarg

    ! Local variables
    INTEGER  :: ii, iref, jj, imin, imax
    REAL :: width, beta, a, newbeta

    REAL, DIMENSION(0:ncube, nhalo) :: prealpha
    REAL, DIMENSION(0:ncube, nhalo) :: newalpha

    REAL, DIMENSION(1-nhalo:ncube+nhalo-1, 1-nhalo:ncube+nhalo-1, 6) :: yarg

    ! Use 0.0 order interpolation to begin
    CALL CubedSphereFillHalo(parg, yarg, np, ncube, nhalo)

    zarg(:,:,np) = yarg(:,:,np)

    ! Calculate the overlapping alpha coordinates
    width = 0.5 * pi / DBLE(ncube-1) ! dx

    DO jj = 1, nhalo
      DO ii = 0, ncube
        prealpha(ii, jj) =  width * (DBLE(ii-1) + 0.5) - 0.25 * pi
        beta             = -width * (DBLE(jj-1) + 0.5) - 0.25 * pi

        CALL CubedSphereABPFromABP(prealpha(ii,jj), beta, 1, 5, &
                                   newalpha(ii,jj), newbeta)
      ENDDO
    ENDDO

    ! Now apply linear interpolation to obtain edge components
    DO jj = 1, nhalo
      ! Reset the reference index
      iref = 2

      ! Interpolation can be applied to more elements after first band
      IF (jj == 1) THEN
        imin = 1
        imax = ncube-1
      ELSE
        imin = 0
        imax = ncube
      ENDIF

      ! Apply linear interpolation
      DO ii = imin, imax
        DO WHILE ((iref .NE. ncube-1) .AND. &
                  (newalpha(ii,jj) > prealpha(iref,jj)))
          iref = iref + 1
        ENDDO

        IF ((newalpha(ii,jj)  >   prealpha(iref-1,jj)) .AND. &
            (newalpha(ii,jj) .LE. prealpha(iref  ,jj)))THEN
          a = (newalpha(ii  ,jj) - prealpha(iref-1,jj)) / &
              (prealpha(iref,jj) - prealpha(iref-1,jj))

          IF ((a < 0.0) .OR. (a > 1.)) THEN
            WRITE (*,*) 'FAIL in CubedSphereFillHalo_Linear'
            WRITE (*,*) 'a out of bounds'
            STOP
          ENDIF

          ! Bottom edge of panel
          zarg(ii, 1-jj, np) =                   &
            (1. - a) * yarg(iref-1, 1-jj, np) + &
                  a  * yarg(iref  , 1-jj, np)

          ! Left edge of panel
          zarg(1-jj, ii, np) =                   &
            (1. - a) * yarg(1-jj, iref-1, np) + &
                  a  * yarg(1-jj, iref  , np)

          ! Top edge of panel
          zarg(ii, ncube+jj-1, np) =                   &
            (1. - a) * yarg(iref-1, ncube+jj-1, np) + &
                  a  * yarg(iref  , ncube+jj-1, np)

          ! Right edge of panel
          zarg(ncube+jj-1, ii, np) =                   &
            (1. - a) * yarg(ncube+jj-1, iref-1, np) + &
                  a  * yarg(ncube+jj-1, iref  , np)

        ELSE
          WRITE (*,*) 'FAIL in CubedSphereFillHalo_Linear'
          WRITE (*,*) 'ii: ', ii, ' jj: ', jj
          WRITE (*,*) 'newalpha: ', newalpha(ii,jj)
          WRITE (*,*) 'prealpha: ', prealpha(iref-1,jj), '-', prealpha(iref,jj)
          STOP
        ENDIF
      ENDDO
    ENDDO

    ! Fill in corner bits
    zarg(0, 0, np) =                          &
      0.25 * (zarg( 1,0,np) + zarg(0, 1,np) + &
              zarg(-1,0,np) + zarg(0,-1,np))
    zarg(0, ncube, np) =                                 &
      0.25 * (zarg( 0,ncube-1,np) + zarg(0,ncube+1,np) + &
              zarg(-1,ncube  ,np) + zarg(1,ncube  ,np))
    zarg(ncube, 0, np) =                                 &
      0.25 * (zarg(ncube-1, 0,np) + zarg(ncube+1,0,np) + &
              zarg(ncube  ,-1,np) + zarg(ncube  ,1,np))
    zarg(ncube, ncube, np) =                                        &
      0.25 * (zarg(ncube-1,ncube  ,np) + zarg(ncube+1,ncube  ,np) + &
              zarg(ncube  ,ncube-1,np) + zarg(ncube  ,ncube+1,np))

  END SUBROUTINE CubedSphereFillHalo_Linear

  !--------------------------------------------------------------------------------
  ! SUBROUTINE CubedSphereFillHalo_Linear_extended
  !
  ! Same as CubedSphereFillHalo_Linear but it also fills the halo i<1 and i>ncube-1
  !
  !---------------------------------------------------------------------------------
  SUBROUTINE CubedSphereFillHalo_Linear_extended(parg, parg_halo, np, ncube, nhalo)

!    USE CubedSphereTrans  ! Cubed sphere transforms

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nhalo
    INTEGER, INTENT(IN) :: np, ncube

    REAL, DIMENSION(ncube-1, ncube-1, 6), INTENT(IN) :: parg

    REAL, DIMENSION(1-nhalo:ncube+nhalo-1, 1-nhalo:ncube+nhalo-1), INTENT(OUT) :: parg_halo

    REAL, DIMENSION(-nhalo:ncube+nhalo, -nhalo:ncube+nhalo):: zarg

    ! Local variables
    INTEGER  :: ii, iref, jj, imin, imax
    REAL :: width,beta, a, newbeta

    REAL, DIMENSION(-nhalo:ncube+nhalo, nhalo+1) :: prealpha !changed compared to non-extended
    REAL, DIMENSION(-nhalo:ncube+nhalo, nhalo+1) :: newalpha !changed compared to non-extended

    REAL, DIMENSION(-nhalo:ncube+nhalo, -nhalo:ncube+nhalo, 6) :: yarg

    ! Use 0.0 order interpolation to begin
    CALL CubedSphereFillHalo(parg, yarg, np, ncube, nhalo+1)

    zarg(:,:) = yarg(:,:,np)

    ! Calculate the overlapping alpha coordinates
    width = 0.5 * pi / DBLE(ncube-1) ! dx

    newalpha = -999999999999.0
    prealpha = -999999999999.0

    DO jj = 1, nhalo+1
!      DO ii = 0, ncube
      DO ii = 1-jj, ncube-1+jj !changed compared to non-extended
        prealpha(ii, jj) =  width * (DBLE(ii-1) + 0.5) - 0.25 * pi
        beta             = -width * (DBLE(jj-1) + 0.5) - 0.25 * pi

        CALL CubedSphereABPFromABP(prealpha(ii,jj), beta   , 1, 5, &
                                   newalpha(ii,jj), newbeta)
      ENDDO
    ENDDO

    ! Now apply linear interpolation to obtain edge components
    DO jj = 1, nhalo+1
      ! Reset the reference index
      iref = 3-jj 

      ! Interpolation can be applied to more elements after first band
      !
      imin = 2-jj       !changed compared to non-extended
      imax = ncube-2+jj !changed compared to non-extended
      !
      ! Apply linear interpolation
      !
      DO ii = imin, imax
        DO WHILE ((iref .NE. ncube-1) .AND. &
                  (newalpha(ii,jj) > prealpha(iref,jj)))
          iref = iref + 1
        ENDDO

        IF ((newalpha(ii,jj)  >   prealpha(iref-1,jj)) .AND.    &
            (newalpha(ii,jj) .LE. prealpha(iref  ,jj)))      &
        THEN
          a = (newalpha(ii,jj)   - prealpha(iref-1,jj)) / &
              (prealpha(iref,jj) - prealpha(iref-1,jj))

          IF ((a < 0.0) .OR. (a > 1.)) THEN
            WRITE (*,*) 'FAIL in CubedSphereFillHalo_Linear'
            WRITE (*,*) 'a out of bounds'
            STOP
          ENDIF

          ! Bottom edge of panel
          zarg(ii, 1-jj) =                   &
            (1. - a) * yarg(iref-1, 1-jj, np) + &
                  a  * yarg(iref  , 1-jj, np)
          
          ! Left edge of panel
          zarg(1-jj, ii) =                   &
            (1. - a) * yarg(1-jj, iref-1, np) + &
                  a  * yarg(1-jj, iref  , np)

          ! Top edge of panel
          zarg(ii, ncube+jj-1) =                   &
            (1. - a) * yarg(iref-1, ncube+jj-1, np) + &
                  a  * yarg(iref  , ncube+jj-1, np)

          ! Right edge of panel
          zarg(ncube+jj-1, ii) =                   &
            (1. - a) * yarg(ncube+jj-1, iref-1, np) + &
                  a  * yarg(ncube+jj-1, iref  , np)

        ELSE
          WRITE (*,*) 'FAIL in CubedSphereFillHalo_Linear'
          WRITE (*,*) 'ii: ', ii, ' jj: ', jj
          WRITE (*,*) 'newalpha: ', newalpha(ii,jj)
          WRITE (*,*) 'prealpha: ', prealpha(iref-1,jj), '-', prealpha(iref,jj)
          STOP
        ENDIF
      ENDDO
    ENDDO

    ! Fill in corner bits
    DO ii=0,nhalo-1
      !
      ! Diagonal lower left
      !
      zarg(0-ii, 0-ii) =                         &
           0.25 * (zarg( 1-ii,0-ii) + zarg(0-ii, 1-ii) + &
                   zarg(-1-ii,0-ii) + zarg(0-ii,-1-ii))
      !
      ! Diagonal upper left
      !
      zarg(0-ii, ncube+ii) =                                 &
      0.25 * (zarg(0-ii,ncube-1+ii) + zarg(0-ii,ncube+1+ii) + &
              zarg(-1-ii,ncube+ii ) + zarg(1-ii,ncube+ii  ))
      !
      ! Diagonal lower right
      !
      zarg(ncube+ii, 0-ii) =                                 &
           0.25 * (zarg(ncube-1+ii, 0-ii) + zarg(ncube+1+ii,0-ii) + &
                   zarg(ncube+ii  ,-1-ii) + zarg(ncube+ii  ,1-ii))
      !
      ! Diagonal upper right
      !
      zarg(ncube+ii, ncube+ii) =                                     &
           0.25 * (zarg(ncube-1+ii,ncube+ii  ) + zarg(ncube+1+ii,ncube+ii  ) + &
                   zarg(ncube+ii  ,ncube-1+ii) + zarg(ncube+ii  ,ncube+1+ii))
    END DO

    parg_halo = zarg(1-nhalo:ncube+nhalo-1, 1-nhalo:ncube+nhalo-1)
  END SUBROUTINE CubedSphereFillHalo_Linear_extended


  !------------------------------------------------------------------------------
  ! SUBROUTINE CubedSphereFillHalo_Cubic
  !
  ! Description:
  !   Recompute the cubed sphere data storage array, with the addition of a
  !   2-element halo region around the specified panel.  Use higher order 
  !   interpolation to translate between panels.
  !
  ! Parameters:
  !   parg - Current panel values
  !   zarg (OUT) - Calculated panel values with halo/ghost region
  !   np - Panel number
  !   ncube - Dimension of the cubed sphere (# of grid lines)
  !------------------------------------------------------------------------------
  SUBROUTINE CubedSphereFillHalo_Cubic(parg, zarg, np, ncube)

!    USE CubedSphereTrans  ! Cubed sphere transforms
!    USE MathUtils         ! Has function for 1D cubic interpolation

    IMPLICIT NONE

    INTEGER , PARAMETER :: nhalo = 2

    INTEGER , INTENT(IN) :: np, ncube

    REAL , DIMENSION(ncube-1, ncube-1, 6), INTENT(IN) :: parg

    REAL ,                                            &
         DIMENSION(1-nhalo:ncube+nhalo-1, 1-nhalo:ncube+nhalo-1, 6), &
         INTENT(OUT) :: zarg


    ! Local variables
    INTEGER  :: ii, iref, ibaseref, jj, imin, imax
    REAL :: width, beta, newbeta

    REAL, DIMENSION(0:ncube, nhalo) :: prealpha
    REAL, DIMENSION(0:ncube, nhalo) :: newalpha

    REAL, DIMENSION(1-nhalo:ncube+nhalo-1, 1-nhalo:ncube+nhalo-1, 6) :: yarg

    ! Use 0.0 order interpolation to begin
    CALL CubedSphereFillHalo(parg, yarg, np, ncube, nhalo)

    zarg(:,:,np) = yarg(:,:,np)

    ! Calculate the overlapping alpha coordinates
    width = 0.5 * pi / DBLE(ncube-1)

    DO jj = 1, nhalo
      DO ii = 0, ncube
        !
        ! alpha,beta for the cell center (extending the panel)
        !
        prealpha(ii, jj) = width * (DBLE(ii-1) + 0.5) - 0.25 * pi
        beta = - width * (DBLE(jj-1) + 0.5) - 0.25 * pi

        CALL CubedSphereABPFromABP(prealpha(ii,jj), beta, 1, 5, &
                                   newalpha(ii,jj), newbeta)
      ENDDO
    ENDDO

    ! Now apply cubic interpolation to obtain edge components
    DO jj = 1, nhalo
      ! Reset the reference index, which gives the element in newalpha that
      ! is closest to ii, looking towards larger values of alpha.
      iref = 2 

      ! Interpolation can be applied to more elements after first band
      !IF (jj == 1) THEN
      !  imin = 1
      !  imax = ncube-1
      !ELSE
        imin = 0
        imax = ncube
      !ENDIF

      ! Apply cubic interpolation
      DO ii = imin, imax
        DO WHILE ((iref .NE. ncube-1) .AND. &
                  (newalpha(ii,jj) > prealpha(iref,jj)))
          iref = iref + 1
        ENDDO

        ! Smallest index for cubic interpolation - apply special consideration
        IF (iref == 2) THEN
          ibaseref = iref-1

        ! Largest index for cubic interpolation - apply special consideration
        ELSEIF (iref == ncube-1) THEN
          ibaseref = iref-3

        ! Normal range
        ELSE
          ibaseref = iref-2
        ENDIF

        ! Bottom edge of panel
        zarg(ii, 1-jj, np) =                                &
          CUBIC_EQUISPACE_INTERP(                           &
            width, newalpha(ii,jj) - prealpha(ibaseref,jj), &
            yarg(ibaseref:ibaseref+3, 1-jj, np))

        ! Left edge of panel
        zarg(1-jj, ii, np) =                                &
          CUBIC_EQUISPACE_INTERP(                           &
            width, newalpha(ii,jj) - prealpha(ibaseref,jj), &
            yarg(1-jj, ibaseref:ibaseref+3, np))

        ! Top edge of panel
        zarg(ii, ncube+jj-1, np) =                          &
          CUBIC_EQUISPACE_INTERP(                           &
            width, newalpha(ii,jj) - prealpha(ibaseref,jj), &
            yarg(ibaseref:ibaseref+3, ncube+jj-1, np))

        ! Right edge of panel
        zarg(ncube+jj-1, ii, np) =                          &
          CUBIC_EQUISPACE_INTERP(                           &
            width, newalpha(ii,jj) - prealpha(ibaseref,jj), &
            yarg(ncube+jj-1, ibaseref:ibaseref+3, np))

      ENDDO
    ENDDO

    ! Fill in corner bits
    zarg(0, 0, np) =                         &
      0.25 * (zarg(1,0,np) + zarg(0,1,np) + &
               zarg(-1,0,np) + zarg(0,-1,np))
    zarg(0, ncube, np) =                                 &
      0.25 * (zarg(0,ncube-1,np) + zarg(0,ncube+1,np) + &
               zarg(-1,ncube,np)  + zarg(1,ncube,np))
    zarg(ncube, 0, np) =                                 &
      0.25 * (zarg(ncube-1,0,np) + zarg(ncube+1,0,np) + &
               zarg(ncube,-1,np)  + zarg(ncube,1,np))
    zarg(ncube, ncube, np) =                                     &
      0.25 * (zarg(ncube-1,ncube,np) + zarg(ncube+1,ncube,np) + &
               zarg(ncube,ncube-1,np) + zarg(ncube,ncube+1,np))

  END SUBROUTINE CubedSphereFillHalo_Cubic

  !------------------------------------------------------------------------------
  ! SUBROUTINE CubedSphereABPFromABP
  !
  ! Description:
  !   Determine the (alpha,beta,idest) coordinate of a source point on
  !   panel isource.
  !
  ! Parameters:
  !   alpha_in - Alpha coordinate in
  !   beta_in - Beta coordinate in
  !   isource - Source panel
  !   idest - Destination panel
  !   alpha_out (OUT) - Alpha coordinate out
  !   beta_out (OUT) - Beta coordiante out
  !------------------------------------------------------------------------------
  SUBROUTINE CubedSphereABPFromABP(alpha_in,  beta_in, isource, idest, &
                                   alpha_out, beta_out)

    IMPLICIT NONE

    REAL, INTENT(IN)  :: alpha_in, beta_in
    INTEGER , INTENT(IN)  :: isource, idest
    REAL, INTENT(OUT) :: alpha_out, beta_out

    ! Local variables
    REAL :: a1, b1
    REAL :: xx, yy, zz
    REAL :: sx, sy, sz

    ! Convert to relative Cartesian coordinates
    a1 = TAN(alpha_in)
    b1 = TAN(beta_in)

    sz = (1. + a1 * a1 + b1 * b1)**(-0.5)
    sx = sz * a1
    sy = sz * b1

    ! Convert to full Cartesian coordinates
    IF (isource == 6) THEN
      yy = sx; xx = -sy; zz = sz

    ELSEIF (isource == 5) THEN
      yy = sx; xx = sy; zz = -sz

    ELSEIF (isource == 1) THEN
      yy = sx; zz = sy; xx = sz

    ELSEIF (isource == 3) THEN
      yy = -sx; zz = sy; xx = -sz

    ELSEIF (isource == 2) THEN
      xx = -sx; zz = sy; yy = sz

    ELSEIF (isource == 4) THEN
      xx = sx; zz = sy; yy = -sz

    ELSE
      WRITE(*,*) 'Fatal Error: Source panel invalid in CubedSphereABPFromABP'
      WRITE(*,*) 'panel = ', isource
      STOP
    ENDIF

    ! Convert to relative Cartesian coordinates on destination panel
    IF (idest == 6) THEN
      sx = yy; sy = -xx; sz = zz

    ELSEIF (idest == 5) THEN
      sx = yy; sy = xx; sz = -zz

    ELSEIF (idest == 1) THEN
      sx = yy; sy = zz; sz = xx

    ELSEIF (idest == 3) THEN
      sx = -yy; sy = zz; sz = -xx

    ELSEIF (idest == 2) THEN
      sx = -xx; sy = zz; sz = yy

    ELSEIF (idest == 4) THEN
      sx = xx; sy = zz; sz = -yy

    ELSE
      WRITE(*,*) 'Fatal Error: Dest panel invalid in CubedSphereABPFromABP'
      WRITE(*,*) 'panel = ', idest
      STOP
    ENDIF
    IF (sz < 0) THEN
      WRITE(*,*) 'Fatal Error: In CubedSphereABPFromABP'
      WRITE(*,*) 'Invalid relative Z coordinate'
      STOP
    ENDIF

    ! Use panel information to calculate (alpha, beta) coords
    alpha_out = ATAN(sx / sz)
    beta_out = ATAN(sy / sz)

  END SUBROUTINE


  !------------------------------------------------------------------------------
  ! FUNCTION CUBIC_EQUISPACE_INTERP
  !
  ! Description:
  !   Apply cubic interpolation on the specified array of values, where all
  !   points are equally spaced.
  !
  ! Parameters:
  !   dx - Spacing of points
  !   x - X coordinate where interpolation is to be applied
  !   y - Array of 4 values = f(x + k * dx) where k = 0,1,2,3
!------------------------------------------------------------------------------
  FUNCTION CUBIC_EQUISPACE_INTERP(dx, x, y)
    
    IMPLICIT NONE
    
    REAL  :: CUBIC_EQUISPACE_INTERP
    REAL  :: dx, x
    REAL , DIMENSION(1:4) :: y
    
    CUBIC_EQUISPACE_INTERP =                                                   &
         (-y(1) / (6.0 * dx**3)) * (x - dx) * (x - 2.0 * dx) * (x - 3.0 * dx) + &
         ( y(2) / (2.0 * dx**3)) * (x) * (x - 2.0 * dx) * (x - 3.0 * dx) +      &
         (-y(3) / (2.0 * dx**3)) * (x) * (x - dx) * (x - 3.0 * dx) +            &
         ( y(4) / (6.0 * dx**3)) * (x) * (x - dx) * (x - 2.0 * dx)
    
  END FUNCTION CUBIC_EQUISPACE_INTERP
  
  subroutine convert_ghost_wind(u,v)
    real, intent(out) :: u(ips:ipe,jps:jpe,ifs:ife)
    real, intent(out) :: v(ips:ipe,jps:jpe,ifs:ife)
    
    real :: matrixG (2,2,ifs:ife)
    real :: matrixIG(2,2,ifs:ife)
    real :: matrixA (2,2,ifs:ife)
    real :: matrixIA(2,2,ifs:ife)
    
    integer i,j,iPatch
    integer gPatch ! ghost patch
    
    do iPatch = ifs, ife
      ! Left boundary
      if(iPatch==1) gpatch = 4
      if(iPatch==2) gpatch = 1
      if(iPatch==3) gpatch = 2
      if(iPatch==4) gpatch = 3
      if(iPatch==5) gpatch = 4
      if(iPatch==6) gpatch = 4
      
      do j = jds, jde
        do i = ips, ids-1
          matrixG  = mesh%matrixG (:,:,i,j,:)
          matrixIG = mesh%matrixIG(:,:,i,j,:)
          matrixA  = mesh%matrixA (:,:,i,j,:)
          matrixIA = mesh%matrixIA(:,:,i,j,:)
          
          call wind_convert_P2P(u(i,j,iPatch),v(i,j,iPatch),gPatch,u(i,j,iPatch),v(i,j,iPatch),iPatch,matrixG,matrixIG,matrixA,matrixIA)
        enddo
      enddo
      
      ! Right boundary
      if(iPatch==1) gpatch = 2
      if(iPatch==2) gpatch = 3
      if(iPatch==3) gpatch = 4
      if(iPatch==4) gpatch = 1
      if(iPatch==5) gpatch = 2
      if(iPatch==6) gpatch = 2
      
      do j = jds, jde
        do i = ide+1, ipe
          matrixG  = mesh%matrixG (:,:,i,j,:)
          matrixIG = mesh%matrixIG(:,:,i,j,:)
          matrixA  = mesh%matrixA (:,:,i,j,:)
          matrixIA = mesh%matrixIA(:,:,i,j,:)
          
          call wind_convert_P2P(u(i,j,iPatch),v(i,j,iPatch),gPatch,u(i,j,iPatch),v(i,j,iPatch),iPatch,matrixG,matrixIG,matrixA,matrixIA)
        enddo
      enddo
      
      ! Top boundary
      if(iPatch==1) gpatch = 5
      if(iPatch==2) gpatch = 5
      if(iPatch==3) gpatch = 5
      if(iPatch==4) gpatch = 5
      if(iPatch==5) gpatch = 3
      if(iPatch==6) gpatch = 1
      
      do j = jde+1, jpe
        do i = ids, ide
          matrixG  = mesh%matrixG (:,:,i,j,:)
          matrixIG = mesh%matrixIG(:,:,i,j,:)
          matrixA  = mesh%matrixA (:,:,i,j,:)
          matrixIA = mesh%matrixIA(:,:,i,j,:)
          
          call wind_convert_P2P(u(i,j,iPatch),v(i,j,iPatch),gPatch,u(i,j,iPatch),v(i,j,iPatch),iPatch,matrixG,matrixIG,matrixA,matrixIA)
        enddo
      enddo
      
      ! Bottom boundary
      if(iPatch==1) gpatch = 6
      if(iPatch==2) gpatch = 6
      if(iPatch==3) gpatch = 6
      if(iPatch==4) gpatch = 6
      if(iPatch==5) gpatch = 1
      if(iPatch==6) gpatch = 3
      
      do j = jps, jds-1
        do i = ids, ide
          matrixG  = mesh%matrixG (:,:,i,j,:)
          matrixIG = mesh%matrixIG(:,:,i,j,:)
          matrixA  = mesh%matrixA (:,:,i,j,:)
          matrixIA = mesh%matrixIA(:,:,i,j,:)
          
          call wind_convert_P2P(u(i,j,iPatch),v(i,j,iPatch),gPatch,u(i,j,iPatch),v(i,j,iPatch),iPatch,matrixG,matrixIG,matrixA,matrixIA)
        enddo
      enddo
      
    enddo
  
  end subroutine convert_ghost_wind

  ! convert vector from patch1 to patch2
  subroutine wind_convert_P2P(u1,v1,patch1,u2,v2,patch2,matrixG,matrixIG,matrixA,matrixIA)
    real   , intent(in ) :: u1
    real   , intent(in ) :: v1
    integer, intent(in ) :: patch1
    real   , intent(out) :: u2
    real   , intent(out) :: v2
    integer, intent(in ) :: patch2
    real   , intent(in ) :: matrixG (2,2,ifs:ife)
    real   , intent(in ) :: matrixIG(2,2,ifs:ife)
    real   , intent(in ) :: matrixA (2,2,ifs:ife)
    real   , intent(in ) :: matrixIA(2,2,ifs:ife)
    
    real contraU1
    real contraV1
    real contraU2
    real contraV2
    real u
    real v
    
    call cov2contrav            (contraU1 , contraV1, u1       , v1      , matrixIG(:,:,Patch1))
    call contravProjPlane2Sphere(u        , v       , contraU1 , contraV1, matrixA (:,:,Patch1))
    call contravProjSphere2Plane(contraU2 , contraV2, u        , v       , matrixIA(:,:,Patch2))
    call contrav2cov            (u2       , v2      , contraU2 , contraV2, matrixG (:,:,Patch2))
  end subroutine wind_convert_P2P
  
END MODULE ghost_mod

