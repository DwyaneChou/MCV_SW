MODULE ghost_mod
  use constants_mod
  use parameters_mod
  use mesh_mod
  
  implicit none
  
  
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
  !   iPatch - Panel index
  !   ncube - Dimension of the cubed sphere (# of grid lines), ncube = nPVx
  !   nhalo - Number of halo/ghost elements around each panel
  !------------------------------------------------------------------------------
  SUBROUTINE CubedSphereFillHalo(parg, zarg, iPatch)

    IMPLICIT NONE

    REAL   , INTENT(IN ) :: parg(ids:ide,jds:jde,ifs:ife)
    REAL   , INTENT(OUT) :: zarg(ips:ipe,jps:jpe,ifs:ife)
    INTEGER, INTENT(IN ) :: iPatch

    ! Local variables
    INTEGER :: ihalo
    
    zarg = 0.
    zarg(ids:ide,jde:jde,iPatch) = parg(ids:ide,jde:jde,iPatch)

    ! Equatorial panels
    IF (iPatch==1) THEN
       DO ihalo=1,nhalo
          zarg(1-ihalo,jds:jde ,iPatch) = parg(ncube-ihalo ,1:ncube-1 ,4)  !exchange left
          zarg(ncube+ihalo-1,1:ncube-1 ,iPatch) = parg(ihalo       ,1:ncube-1 ,2)  !exchange right
          zarg(1:ncube-1 ,1-ihalo      ,iPatch) = parg(1:ncube-1,ncube-ihalo  ,6)  !exchange below
          zarg(1:ncube-1 ,ncube+ihalo-1,iPatch) = parg(1:ncube-1,ihalo        ,5)  !exchange over
       ENDDO
    ELSE IF (iPatch==2) THEN
       DO ihalo=1,nhalo
          zarg(1-ihalo      ,1:ncube-1 ,2) = parg(ncube-ihalo,1:ncube-1   ,1)  !exchange left
          zarg(ncube+ihalo-1,1:ncube-1 ,2) = parg(ihalo      ,1:ncube-1   ,3)  !exchange right
          zarg(1:ncube-1 ,1-ihalo      ,2) = parg(ncube-ihalo,ncube-1:1:-1,6)  !exchange below
          zarg(1:ncube-1 ,ncube+ihalo-1,2) = parg(ncube-ihalo,1:ncube-1   ,5)  !exchange over
       ENDDO
    ELSE IF (iPatch==3) THEN
       DO ihalo=1,nhalo
          zarg(1-ihalo      ,1:ncube-1 ,3) = parg(ncube-ihalo    ,1:ncube-1,2) !exchange left
          zarg(ncube+ihalo-1,1:ncube-1 ,3) = parg(ihalo          ,1:ncube-1,4) !exchange right
          zarg(1:ncube-1 ,1-ihalo      ,3) = parg(ncube-1:1:-1,ihalo       ,6) !exchange below
          zarg(1:ncube-1 ,ncube+ihalo-1,3) = parg(ncube-1:1:-1,ncube-ihalo ,5) !exchange over
       ENDDO
    ELSE IF (iPatch==4) THEN
       DO ihalo=1,nhalo
          zarg(1-ihalo      ,1:ncube-1 ,4) = parg(ncube-ihalo,1:ncube-1   ,3) !exchange left
          zarg(ncube+ihalo-1,1:ncube-1 ,4) = parg(ihalo      ,1:ncube-1   ,1) !exchange right
          zarg(1:ncube-1 ,1-ihalo      ,4) = parg(ihalo      ,1:ncube-1   ,6) !exchange below
          zarg(1:ncube-1 ,ncube+ihalo-1,4) = parg(ihalo      ,ncube-1:1:-1,5) !exchange over
       ENDDO
    ! Top panel
    ELSE IF (iPatch==5) THEN
       DO ihalo=1,nhalo
          zarg(1-ihalo      ,1:ncube-1 ,5) = parg(ncube-1:1:-1,ncube-ihalo,4) !exchange left
          zarg(ncube+ihalo-1,1:ncube-1 ,5) = parg(1:ncube-1   ,ncube-ihalo,2) !exchange right
          zarg(1:ncube-1 ,1-ihalo      ,5) = parg(1:ncube-1   ,ncube-ihalo,1) !exchange below
          zarg(1:ncube-1 ,ncube+ihalo-1,5) = parg(ncube-1:1:-1,ncube-ihalo,3) !exchange over
       ENDDO
    ! Bottom panel
    ELSE IF (iPatch==6) THEN
       DO ihalo=1,nhalo
          zarg(1-ihalo      ,1:ncube-1 ,6) = parg(1:ncube-1   ,ihalo      ,4) !exchange left
          zarg(ncube+ihalo-1,1:ncube-1 ,6) = parg(ncube-1:1:-1,ihalo      ,2) !exchange right
          zarg(1:ncube-1 ,1-ihalo      ,6) = parg(ncube-1:1:-1,ihalo      ,3) !exchange below
          zarg(1:ncube-1 ,ncube+ihalo-1,6) = parg(1:ncube-1   ,ihalo      ,1) !exchange over
       ENDDO
    ELSE
       WRITE (*,*) 'Fatal error: In CubedSphereFillHalo'
       WRITE (*,*) 'Invalid panel id ', iPatch
       STOP
    ENDIF
  END SUBROUTINE CubedSphereFillHalo
  
  ! Correct wind on ghost cells
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

