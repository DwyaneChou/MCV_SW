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
    real   , dimension(:,:), allocatable :: coef  ! Position in cell
    integer, dimension(:,:), allocatable :: iref  ! Index of reference point
  end type ghostLocation
  
  type(ghostLocation) :: ghost
  contains
  
  subroutine initGhost
#ifdef CUBE
    integer :: i,j
    integer :: iPatch ! computational patch
    integer :: gpatch ! ghost patch
    integer :: ihalo  ! halo index
    real    :: lambda, theta
    real    :: coef
    integer :: iref
    
    real    :: X_RAW(Nx+1)
    
    nPVHalo = xhalo * (DOF - 1)
    
    allocate(ghost%X   (ids:ide,nPVHalo))
    allocate(ghost%Y   (ids:ide,nPVHalo))
    allocate(ghost%coef(ids:ide,nPVHalo))
    allocate(ghost%iref(ids:ide,nPVHalo))
    
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
    
    X_RAW = mesh%xP(pvIdx(1,1:Nx+1),1,iPatch)
    
    ! Calculate reference points
    ! Now apply linear interpolation to obtain edge components
    DO j = 1, nPVHalo
      ! Reset the reference index
      iref = 1
      
      DO i = ids, ide
        !
        ! Find reference points
        !
        DO WHILE (ghost%X(i,j) > X_RAW(iref))
          iref = iref + 1
        ENDDO
        
        IF ((ghost%X(i,j) >= X_RAW(iref-1)) .AND. (ghost%X(i,j) <= X_RAW(iref  )))THEN
            
          coef = ghost%X(i,j) - X_RAW(iref-1)
          
          ghost%coef(i,j) = coef
          ghost%iref(i,j) = iref - 1
          
        ENDIF
      ENDDO
    ENDDO
#endif
  end subroutine initGhost
  
#ifdef CUBE  
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
        call linear_interp(ghost_target(ids-j,jds:jde,iPatch),field(ids-j,jds:jde,iPatch),ghost%iref(:,j),ghost%coef(:,j))
        ! Rigth boundary
        call linear_interp(ghost_target(ide+j,jds:jde,iPatch),field(ide+j,jds:jde,iPatch),ghost%iref(:,j),ghost%coef(:,j))
        ! Top boundary
        call linear_interp(ghost_target(ids:ide,jde+j,iPatch),field(ids:ide,jde+j,iPatch),ghost%iref(:,j),ghost%coef(:,j))
        ! Bottom boundary
        call linear_interp(ghost_target(ids:ide,jds-j,iPatch),field(ids:ide,jds-j,iPatch),ghost%iref(:,j),ghost%coef(:,j))
      enddo
    enddo
    
    field = ghost_target
  
  end subroutine CubedSphereFillGhost
#endif  
  
  subroutine linear_interp(dest,src,iref,coef)
    real   , intent(out) :: dest(ids:ide)
    real   , intent(in ) :: src (ids:ide)
    integer, intent(in ) :: iref(ids:ide)
    real   , intent(in ) :: coef(ids:ide)
    
    integer i
    integer P1,P2,P3,P4
    real    q1,q2,q3,q4
    
#ifdef MCV3
      ! For MCV3 only
      do i = ids,ide
        P1 = pvIdx(1,iref(i))
        P2 = pvIdx(2,iref(i))
        P3 = pvIdx(3,iref(i))
        
        q1 = src(P1)
        q2 = src(P2)
        q3 = src(P3)
        
        dest(i) = q1 - ((3.*q1 - 4.*q2 + q3)*coef(i)) / dx + (2.*(q1 - 2.*q2 + q3)*coef(i)**2) / (dx**2)
      
      enddo
#endif

#ifdef MCV4
      ! For MCV4 only
      do i = ids,ide
        P1 = pvIdx(1,iref(i))
        P2 = pvIdx(2,iref(i))
        P3 = pvIdx(3,iref(i))
        P4 = pvIdx(4,iref(i))
        
        q1 = src(P1)
        q2 = src(P2)
        q3 = src(P3)
        q4 = src(P4)
        
        dest(i) = (2.*dx**3.*q1 + dx**2.*(-11.*q1 + 18.*q2 - 9.*q3 + 2.*q4) * coef(i) + 9.*dx*(2.*q1 - 5.*q2 + 4.*q3 - q4)*coef(i)**2. + 9.*(-q1 + 3.*q2 - 3.*q3 + q4)*coef(i)**3.) / (2.*dx**3.)
      
      enddo
#endif
    
  end subroutine linear_interp
  
#ifdef LONLAT
  subroutine lonlatFillGhost(field)
  real, intent(inout) :: field(ips:ipe,jps:jpe,ifs:ife)
  
  integer i,j
  integer iref, jref
  
  ! Zonal
  do j = jds, jde
    do i = ips, ids
      iref         = ide + i - 1
      field(i,j,:) = field(iref,j,:)
    enddo
  enddo
  
  do j = jds, jde
    do i = ide, ipe
      iref         = ids + i - ide
      field(i,j,:) = field(iref,j,:)
    enddo
  enddo
  
  ! Meridional, the points over the pole are selected on opposite location
  do i = ids, ide
    iref = i + (nPVx - 1) / 2 - int(2 * i / (nPVx - 1)) * (nPVx - 1)
    if(i == (nPVx - 1) / 2) iref = nPVx - 1
    if(i >=  nPVx - 1     ) iref = (nPVx - 1) / 2 + i - nPVx + 1
    
    do j = jps, jds
      jref         = jds - j + 1
      field(i,j,:) = field(iref,jref,:)
    enddo
    
    do j = jde, jpe
      jref         = jde - (j - jde)
      field(i,j,:) = field(iref,jref,:)
    enddo
  enddo
  
  end subroutine lonlatFillGhost
#endif  
  
  ! Fill up halo with ghost points  
  subroutine fill_ghost(stat)
    type(stat_field), intent(inout) :: stat
#ifdef CUBE
    call CubedSphereFillGhost(stat%zonal_wind     )
    call CubedSphereFillGhost(stat%meridional_wind)
    call CubedSphereFillGhost(stat%phi            )
#endif

#ifdef LONLAT
    call lonlatFillGhost(stat%u   )
    call lonlatFillGhost(stat%v   )
    call lonlatFillGhost(stat%phi )
    call lonlatFillGhost(stat%phiG)
#endif
  end subroutine fill_ghost
  
END MODULE ghost_mod