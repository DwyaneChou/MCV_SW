MODULE domain

  USE global ,only:nummoi,numvar,d2r,np

  type glob

    ! Namelist parameters
    real  :: dt     !  time step
    real  :: dx     !  grid-spacing in the x-direction
    real  :: dy     !  grid-spacing in the y-direction
    real  :: dz     !  grid-spacing in the z-direction

    real  :: x_min  !  start location of x-direction
    real  :: x_max  !  end location of x-direction
    real  :: y_min  !  start location of y-direction
    real  :: y_max  !  end location of y-direction
    real  :: z_min  !  start location of z-direction
    real  :: z_max  !  end location of z-direction

    integer  :: xhalo  !  halo number of x-diretion
    integer  :: yhalo  !  halo number of y-diretion
    integer  :: zhalo  !  halo number of z-diretion

    integer  :: xhalom  ! moisture halo number of x-diretion
    integer  :: yhalom  ! moisture halo number of y-diretion
    integer  :: zhalom  ! moisture halo number of z-diretion

    integer  :: max_step     ! the max-runing step of model
    integer  :: verticalstretching   ! The options of vertical stretching grid. 
    character(len=8)   :: gridsys    ! Grid system used: Cubed grid or Latlon grid system.

    !  other physical-porcess parameters
    !...

    ! Index parameter
    integer  :: ids    ! The starting index in the x-direction (Physical domain)
    integer  :: ide    ! The ending index in the x-direction  (Physical domain)
    integer  :: jds    ! The starting index in the y-direction  (Physical domain)
    integer  :: jde    ! The ending index in the y-direction  (Physical domain)
    integer  :: kds    ! The starting index in the z-direction  (Physical domain)
    integer  :: kde    ! The ending index in the z-direction  (Physical domain)

    integer  :: ics    ! The starting index in the x-direction (Physical cell/element domain)
    integer  :: ice    ! The ending index in the x-direction  (Physical cell/element domain)
    integer  :: jcs    ! The starting index in the y-direction  (Physical cell/element domain)
    integer  :: jce    ! The ending index in the y-direction  (Physical cell/element domain)
    integer  :: kcs    ! The starting index in the z-direction  (Physical cell/element domain)
    integer  :: kce    ! The ending index in the z-direction  (Physical cell/element domain)

    integer  :: ims    ! The starting index in the x-direction (Memory domain)
    integer  :: ime    ! The ending index in the x-direction  (Memory domain)
    integer  :: jms    ! The starting index in the y-direction  (Memory domain)
    integer  :: jme    ! The ending index in the y-direction  (Memory domain)
    integer  :: kms    ! The starting index in the z-direction  (Memory domain)
    integer  :: kme    ! The ending index in the z-direction  (Memory domain)

    integer  :: Nx     ! Element numbers in the x-direction
    integer  :: Ny     ! Element numbers in the y-direction
    integer  :: Nz     ! Element numbers in the z-direction
    integer  :: Nf     ! The number of cube faces

    integer  :: Npvx     ! Point-value numbers in the x-direction
    integer  :: Npvy     ! Point-value numbers in the y-direction
    integer  :: Npvz     ! Point-value numbers in the z-direction

    integer  :: Nll,Ntt  ! grid points in the lambda/theta direction.

    contains
      procedure,pass :: initialize   => InitParameter
      procedure,pass :: setIndex     => get_ijk_index
      procedure,pass :: setpvindex   => index_via2pv
      procedure,pass :: setviaindex  => index_pv2via
      procedure,pass :: setcubeindex => get_cube_index
  end type glob

  !/// MCV basic definiton ///
  type,extends(glob)  :: cell
    real,dimension(:,:,:,:),pointer :: pv   !  Point values.
    real,dimension(:,:,:,:),pointer :: via  !  Volume integrated average.
  end type cell

  type,extends(glob)  :: jame
    real,dimension(:,:),pointer       :: jac      !  Jacobian of transformation
    real,dimension(:,:,:),pointer     :: metric   !  Horizontal metrics
    real,dimension(:,:,:,:),pointer   :: jacv     !  Jacobian of vertical transformation
    real,dimension(:,:,:,:,:),pointer :: metricv  !  Vertical transformation metrics
    real,dimension(:,:,:,:),pointer   :: jab      !  Composite Jacobian of transformation
  end type jame

  type,extends(glob)  :: gloc            ! ghost cell location
    integer,dimension(:,:),pointer     :: bci
    real,dimension(:),pointer          :: bcr
  end type gloc

  type,extends(glob)  :: cord
    real,dimension(:),pointer   :: x 
    real,dimension(:),pointer   :: y 
    real,dimension(:),pointer   :: z 
    real,dimension(:),pointer   :: xc
    real,dimension(:),pointer   :: yc
    real,dimension(:),pointer   :: zc
  end type cord

  type,extends(glob)  :: gmap           ! geometric vertical mapping
    real,dimension(:,:,:),pointer     :: zs   ! Topography height
    real,dimension(:,:,:,:),pointer   :: zz   ! The geometry height/altitude containing the terrain.
  end type gmap

  type :: pv2d                   ! 2D array on six faces
    real,dimension(:,:,:),pointer   :: pv2d
  end type pv2d
  type :: pv3d                   ! 3D array on six faces
    real,dimension(:,:,:,:),pointer   :: pv3d
  end type pv3d

  ! Declare the model varialbes
  type(glob),pointer  :: gp      ! global parameters
  type(cord),pointer  :: coord   ! coordinate
  type(cell),dimension(:),pointer   :: q,qh       ! prognostic variables, hydrostatic model variables.
  type(cell),dimension(:),pointer   :: qmoist     ! prognostic moisture variables
  type(cell),dimension(:),pointer   :: qn,q0      ! storage of prognostic variables
  type(cell),dimension(:),pointer   :: rqsb       ! hydrostatic R(qh)
  type(cell),pointer  :: ph,dpx,dpy      ! hydrostatic pressure, its derivatives.
  type(gmap),pointer  :: gm              ! geometric mapping quantity.
  type(jame),pointer  :: jm              ! Jacobian and metrics
  type(gloc),pointer  :: gl              ! ghost cell location info. for cubed-sphere.

contains

  subroutine InitParameter(gp)

    implicit none

    class(glob) :: gp
    ! Local variables
    integer   :: ids,ide,jds,jde,kds,kde,   &
                 ims,ime,jms,jme,kms,kme,   &
                 ics,ice,jcs,jce,kcs,kce,   &
                 Nf
    character(len=256) :: namelist_file
    character(len=9)   :: gridsys
    integer   :: max_step,xhalo,yhalo,zhalo,verticalstretching,Nx,Ny,Nz,xhalom,yhalom,zhalom
    real      :: dt,dx,dy,dz,x_min,x_max,y_min,y_max,z_min,z_max
    namelist /namelist_dyn/ gridsys,verticalstretching,max_step,dt,Nx,Ny,Nz,dx,dy,dz, &
                            x_min,x_max,y_min,y_max,z_min,z_max,xhalo,yhalo,zhalo,xhalom,yhalom,zhalom
    logical            :: file_stat
    integer            :: iostat,fileunit
    ! ------------------------
    ! ///////////////////////
    !   Read namelist file
    ! //////////////////////
    ! &namelist_dyn
    ! default values
    gridsys ='CUBE'
    verticalstretching=0   ! default: equidistant
    max_step = 1     ! Max running steps
    nx    = 2
    ny    = 2
    nz    = 2

    dt    = 0.0      ! time step, unit: seconds
    dx    = 1.0      ! x-direction resolution, unit: degree
    dy    = 1.0      ! y-direction resolution, unit: degree
    dz    = 1.0      ! Maybe not use, unit: m
    xhalo = 1        ! ghost cell numbers in x-direction
    yhalo = 1        ! ghost cell numbers in y-direction
    zhalo = 1        ! ghost cell numbers in z-direction
    xhalom= 1        ! moistures ghost cell numbers in x-direction
    yhalom= 1        ! moistures ghost cell numbers in y-direction
    zhalom= 1        ! moistures ghost cell numbers in z-direction
    x_min =-45.      ! Unit: degree, starting location in the x-computational domain
    x_max = 45.      ! Unit: degree, ending location in the x-computational domain
    y_min =-45.      ! Unit: degree, starting location in the y-computational domain
    y_max = 45.      ! Unit: degree, ending location in the y-computational domain
    z_min = 0.       ! Bottom: unit: m, starting location in z-the computational domain
    z_max = 10000.   ! Top level, unit: m, ending location in z-the computational domain

    !---------------------------
    !   Read namelist file
    !---------------------------
    namelist_file="./namelist.atm"

    inquire(file=trim(namelist_file),exist=file_stat)
    if(file_stat)then
      fileunit = 10
      open(fileunit,file=trim(namelist_file),form="formatted",status='old')
    else
      write(6,*) "Failure in opening namelist file" 
      stop
    endif

    !---------------------------
    !   namelist_dyn
    !---------------------------
    rewind(fileunit)
    read  (fileunit,nml=namelist_dyn,iostat=iostat)

    if( iostat > 0 ) then
      print*, "Failure in reading namelist &namelist_dyn."
      stop
    end if
    gridsys=trim(gridsys)

    if(gridsys=='CUBE' .or. gridsys=='cube' )then
      write(6,*)
      write(6,*) "//////////////////////////////////////////////////"
      write(6,*) "  Current spherical system is cuded-sphere grid "
      write(6,*) "//////////////////////////////////////////////////"
      write(6,*)

      x_min = -45.
      x_max =  45.
      y_min = -45.
      y_max =  45

    elseif(gridsys=='LONLAT' .or. gridsys=='lonlat' )then  ! LATLON grid
      write(6,*)
      write(6,*) "///////////////////////////////////////////////////////"
      write(6,*) "  Current spherical system is latitude-longitude grid "
      write(6,*) "///////////////////////////////////////////////////////"
      write(6,*)

    elseif( gridsys=='CARTESIAN' .or. gridsys=='cartesian' )then
      write(6,*)
      write(6,*) "///////////////////////////////////////////////////////"
      write(6,*) "     Current grid system is Cartesian system "
      write(6,*) "///////////////////////////////////////////////////////"
      write(6,*)
    else
      stop " Give the correct name of gridsys in the namelist file!"
    endif

    ! check domain size
    if(AINT((x_max-x_min)/dx)/=nx)then
      stop 'Check the values of dx and nx which is not matched. Stop!'
    endif
    if(AINT((y_max-y_min)/dy)/=ny)then
      stop 'Check the values of dy and ny which is not matched. Stop!'
    endif

    ! /*****  set moisture halo *****/
    xhalom=xhalo
    yhalom=yhalo
    zhalom=zhalo

    write(6,*)
    write(6,*)'  /// Namelist information ///  '
    write(6,*)

    write(*,nml=namelist_dyn)

    ! Convert degree into radian 
    dx = dx*D2R
    dy = dy*D2R
    x_min = x_min*D2R
    x_max = x_max*D2R
    y_min = y_min*D2R
    y_max = y_max*D2R

    close(fileunit)

! //////////////////////////////////
!   set the element/cell numbers
! /////////////////////////////////

    !Nx = AINT((x_max-x_min)/dx)
    !Ny = AINT((y_max-y_min)/dy)
    !Nz = AINT((z_max-z_min)/dz)

! /////////////////////////////
!   assignment of parameters
! /////////////////////////////
    if(gridsys=='CUBE' .or. gridsys=='cube')then
      Nf   = 6        ! For cubed grid system.
    elseif(gridsys=='LONLAT' .or. gridsys=='lonlat' )then  ! LATLON grid
      Nf   = 1
    elseif( gridsys=='CARTESIAN' .or. gridsys=='cartesian' )then   ! Cartesian sys.
      Nf   = 1
    else
      stop " Please give the correct name of gridsys in the namelist file!"
    endif

    ids=1;ide=2*Nx+1
    jds=1;jde=2*Ny+1
    kds=1;kde=2*Nz+1

    ! parent objects
    gp%dx = dx
    gp%dy = dy
    gp%dz = dz
    gp%x_min = x_min
    gp%x_max = x_max 
    gp%y_min = y_min
    gp%y_max = y_max
    gp%z_min = z_min
    gp%z_max = z_max

    gp%xhalo = xhalo
    gp%yhalo = yhalo
    gp%zhalo = zhalo
      
    gp%xhalom= xhalom
    gp%yhalom= yhalom
    gp%zhalom= zhalom

    gp%max_step = max_step
    gp%verticalstretching = verticalstretching

    gp%gridsys = "cube"
    gp%dt = dt

    !
    ims=ids-2*xhalo
    ime=ide+2*xhalo
    jms=jds-2*yhalo
    jme=jde+2*yhalo
    kms=kds-2*zhalo
    kme=kde+2*zhalo

    ics=1 -xhalom
    ice=Nx+xhalom
    jcs=1 -yhalom
    jce=Ny+yhalom
    kcs=1 -zhalom
    kce=Nz+zhalom

    gp%ids = ids
    gp%ide = ide
    gp%jds = jds
    gp%jde = jde
    gp%kds = kds
    gp%kde = kde

    gp%ics = ics
    gp%ice = ice
    gp%jcs = jcs
    gp%jce = jce
    gp%kcs = kcs
    gp%kce = kce

    gp%ims = ims
    gp%ime = ime
    gp%jms = jms
    gp%jme = jme
    gp%kms = kms
    gp%kme = kme

    gp%Nx  = Nx
    gp%Ny  = Ny
    gp%Nz  = Nz
    gp%Nf  = Nf

    gp%Npvx= ide-ids+1
    gp%Npvy= jde-jds+1
    gp%Npvz= kde-kds+1

    if(gridsys=='CUBE' .or. gridsys=='cube')then
      gp%Nll = 2*Nx*4+0
      gp%Ntt = 2*Ny*2+0
    elseif( gridsys=='LONLAT' .or. gridsys=='lonlat' )then
      gp%Nll = 2*Nx+1
      gp%Ntt = 2*Ny+1
    elseif( gridsys=='CARTESIAN' .or. gridsys=='cartesian' )then
      gp%Nll = 2*Nx+1
      gp%Ntt = 2*Ny+1
    else
      stop " Give the correct name of gridsys in the namelist file!"
    endif

    select type (gp)
    type is (glob)
      ! do nothing
    class is (cell)
      print*,' Type cell array allocated by subroutine AllocateTypecell'
    class default
      stop 'initialize: unexpected type for gp object!'
    end select

    write(6,*)''
    write(6,*)' /// PV moment range information /// '
    write(6,*)''
    write(6,*)'  Physical domain range: '
    write(6,*)'    Start/ending index of x-direction -- ids, ide =',ids,ide
    write(6,*)'    Start/ending index of y-direction -- jds, jde =',jds,jde
    write(6,*)'    Start/ending index of z-direction -- kds, kde =',kds,kde
    write(6,*)'  Memory domain range: '
    write(6,*)'    Start/ending index of x-direction -- ims, ime =',ims,ime
    write(6,*)'    Start/ending index of y-direction -- jms, jme =',jms,jme
    write(6,*)'    Start/ending index of z-direction -- kms, kme =',kms,kme
    write(6,*)''

  end subroutine InitParameter

  subroutine get_ijk_index(gp,ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,ics,ice,jcs,jce,kcs,kce)

    implicit none

    class(glob),intent(in)     :: gp
    integer,intent(out)        :: ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme
    integer,optional,intent(out) :: ics,ice,jcs,jce,kcs,kce
    ! -----------------
    ids = gp%ids
    jds = gp%jds
    kds = gp%kds

    ide = gp%ide
    jde = gp%jde
    kde = gp%kde

    ims = gp%ims
    jms = gp%jms
    kms = gp%kms

    ime = gp%ime
    jme = gp%jme
    kme = gp%kme

    if(present(ics))ics = gp%ics
    if(present(ice))ice = gp%ice
    if(present(jcs))jcs = gp%jcs
    if(present(jce))jce = gp%jce
    if(present(kcs))kcs = gp%kcs
    if(present(kce))kce = gp%kce

  end subroutine get_ijk_index

  subroutine get_cube_index(gp,Nx,Ny,Nz,Nf)

    implicit none

    class(glob),intent(in)    :: gp
    integer,optional,intent(out)  :: Nx,Ny,Nz,Nf
    ! -----------------
    if(present(Nx))Nx = gp%Nx
    if(present(Ny))Ny = gp%Ny
    if(present(Nz))Nz = gp%Nz
    if(present(Nf))Nf = gp%Nf

  end subroutine get_cube_index

  subroutine index_pv2via(gp,ic,jc,m,n,i,j)

  ! i,j: pv point index
  ! ic,jc: cell index
  ! m,n :: entry index whithin the cell

    implicit none

    class(glob),intent(in)     :: gp
    integer,intent(in)         :: i,j
    integer,intent(inout)      :: ic,jc,m,n
  ! Local variables
    integer  :: ids,ide,jds,jde,kds,kde, &
                ims,ime,jms,jme,kms,kme
  ! -------------------
    call gp%setIndex( ids,ide,jds,jde,kds,kde,  &
                      ims,ime,jms,jme,kms,kme )

    if(mod(i,2)==0)then
      ic=i/2
    else
      ic=(i-1)/2+1
    endif
    if(i==ide)then
      ic=gp%Nx    ! i=ide, put it back one grid cell
    endif

    if(mod(j,2)==0)then
      jc=j/2
    else
      jc=(j-1)/2+1
    endif
    if(j==jde)then
      jc=gp%Ny   ! j=jde, put it back one grid cell
    endif

    ! entry indices in the cell
    m=i-(2*(ic-1)+1)+1
    n=j-(2*(jc-1)+1)+1
    if(i==ime)then
      m=3
      ic=(i-1)/2
    endif
    if(j==jme)then
      n=3
      jc=(j-1)/2
    endif

  end subroutine index_pv2via

  subroutine index_via2pv(gp,i,j,ic,jc)

  ! i,j: pv point index
  ! ic,jc: cell index

    implicit none

    class(glob),intent(in)     :: gp
    integer,intent(in)         :: ic,jc
    integer,intent(inout)      :: i,j
  ! Local variables

  ! -------------------
    i=2*(ic-1)+1
    j=2*(jc-1)+1

  end subroutine index_via2pv

  subroutine SetComputationalSpace

    implicit none

    integer  :: i,j,k,ii,jj,kk,n1,n2,n3
    integer  :: ids,ide,jds,jde,kds,kde, &
                ims,ime,jms,jme,kms,kme, &
                ics,ice,jcs,jce,kcs,kce
    real     :: xstart,ystart,zstart
    real     :: b,dzbar
    real     :: D1,D2,D3,D,avezeta,hs,dmumin,dmum,dmuu,aa,al,delta,summ
    real,dimension(:),pointer :: zbar
    ! ------------------------
    call gp%setIndex( ids,ide,jds,jde,kds,kde,  &
                      ims,ime,jms,jme,kms,kme,  & 
                      ics,ice,jcs,jce,kcs,kce ) 

    xstart=gp%x_min-gp%xhalo*gp%dx
    ystart=gp%y_min-gp%yhalo*gp%dy
    zstart=gp%z_min-gp%zhalo*gp%dz

    do i=ims,ime
      coord%x(i) = xstart + (i-ims)*gp%dx/2.
    enddo
    do i=ics,ice
      ii = 2*(i-1)+1
      coord%xc(i) = (coord%x(ii)+coord%x(ii+2))/2.
    enddo

    do j=jms,jme
      coord%y(j) = ystart + (j-jms)*gp%dx/2.
    enddo
    do j=jcs,jce
      jj = 2*(j-1)+1
      coord%yc(j) = (coord%y(jj)+coord%y(jj+2))/2.
    enddo

    ! z-gridspacing
    allocate(zbar(kms:kme))
    if(gp%verticalstretching==0)then ! equidistant grid spacing.
      do k=kms,kme
        coord%z(k) = zstart +(k-kms)*gp%dz/2.
      enddo
      do k=1,gp%Nz
        kk = 2*(k-1)+1
        coord%zc(k)= (coord%z(kk)+coord%z(kk+1))/2.
      enddo
    elseif(gp%verticalstretching==1)then
      b=10.;dzbar=1./float(gp%Nz)
      do k=kds,kde,2
        kk=(k-1)/2+1
        zbar(k)=float(kk-1)*dzbar
        zbar(k+1)=zbar(k)+dzbar/2.
      enddo

      do k=kds,kde,2
        coord%z(k)=gp%z_min+(gp%z_max-gp%z_min)*(sqrt(b*zbar(k)**2.+1.)-1.)/(sqrt(b+1.)-1.)
      enddo
      do k=1,gp%Nz
        kk=2*k
        coord%z(kk)=(coord%z(kk-1)+coord%z(kk+1))/2.
      enddo

      do k=kms,ids-1
        coord%z(k) = 2.*gp%z_min-coord%z(kds-k+1)
      enddo

      do k=kde+1,kme
        coord%z(k) = 2.*gp%z_max-coord%z(kde-(k-kde))
      enddo
    elseif(gp%verticalstretching==2)then
      b=10.;dzbar=1./float(gp%Nz)
      do k=kds,kde,2
        kk=(k-1)/2+1
        zbar(k)=float(kk-1)*dzbar
        zbar(k+1)=zbar(k)+dzbar/2.
      enddo
      do k=kds,kde,2
        coord%z(k)=gp%z_min+(gp%z_max-gp%z_min)*dzbar*float(k-1)/2.
      enddo

      do k=1,gp%Nz
        kk=2*k
        coord%z(kk)=(coord%z(kk-1)+coord%z(kk+1))/2.
      enddo
      do k=kms,ids-1
        coord%z(k) = 2.*gp%z_min-coord%z(kds-k+1)
      enddo

      do k=kde+1,kme
        coord%z(k) = 2.*gp%z_max-coord%z(kde-(k-kde))
      enddo

    elseif(gp%verticalstretching==3)then
      hs=0.
      avezeta=(gp%z_max-gp%z_min)/gp%Nz
      D =gp%z_max-gp%z_min
      D1=0.
      D2=D
      D3=D-D1-D2
      dmumin=377.
      n1=AINT(D1/dmumin)
      n3=0
      n2=gp%Nz-n1-n3
      dmum=(D1*dmumin - 2*D2*dmumin - D3*dmumin - dmumin**2*gp%Nz -  &
         Sqrt(4*D2*dmumin**2*(2*D1 - 2*dmumin*gp%Nz) +   &
           (-(D1*dmumin) + 2*D2*dmumin + D3*dmumin + dmumin**2*gp%Nz)**2))/  &
         (4.*(D1 - dmumin*gp%Nz))
      dmuu=2.*dmum-dmumin

      al=0.4
      aa=(1+n2)/2.
      summ=gp%z_min

      do k=1,n1  ! cell index
        delta=dmumin
        kk=2*(k-1)+1
        coord%z(kk) = summ
        coord%z(kk+1) = summ + delta/2.
        summ=summ+delta 
      enddo

      do k=n1+1,n2+n1
        delta=dmum+(dmumin-dmum)/tanh(2.*al)*tanh(2.*al/(1.-aa)*(k-aa))
        kk=2*(k-1)+1
        coord%z(kk) = summ
        coord%z(kk+1) = summ + delta/2.
        summ=summ+delta 
      enddo

      do k=n1+n2+1,n1+n2+n3
        delta=dmuu
        kk=2*(k-1)+1
        coord%z(kk) = summ
        coord%z(kk+1) = summ + delta/2.
        summ=summ+delta 
      enddo
      coord%z(kk+2) = gp%z_max
      
      do k=kms,ids-1
        coord%z(k) = 2.*gp%z_min-coord%z(kds-k+1)
      enddo

      do k=kde+1,kme
        coord%z(k) = 2.*gp%z_max-coord%z(kde-(k-kde))
      enddo
    else
      print*,'Verticalstretching /= 0 or 1 or 2. Stop!'
      stop
    endif

    write(6,*)''
    write(6,*)' /// PV moment vertical level information /// '
    write(6,*)''
    write(6,*)'  Computational heights from bottom to top'
    write(6,*)coord%z(kds:kde)
    write(6,*)''

    deallocate(zbar)

  end subroutine SetComputationalSpace

  subroutine copy(q,qn,nvar,cha)

    implicit none

    integer                                  ,intent(in )   :: nvar
    type(cell)      ,dimension(nvar)         ,intent(in )   :: qn
    type(cell)      ,dimension(nvar)         ,intent(out)   :: q
    character(len=3)                ,optional,intent(in )   :: cha
  ! Local variables
    integer  :: i,j,k,n,ivar
  ! -----------------------

!!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,n,ivar)
    do ivar=1,nvar
      do n=1,gp%Nf

        do k=gp%kms,gp%kme
          do j=gp%jms,gp%jme
            do i=gp%ims,gp%ime
              q(ivar)%pv(i,j,k,n)=qn(ivar)%pv(i,j,k,n)
            enddo
          enddo
        enddo

        if(present(cha))then
          if(cha=='via')then
            do k=gp%kcs,gp%kce
              do j=gp%jcs,gp%jce
                do i=gp%ics,gp%ice
                  q(ivar)%via(i,j,k,n)=qn(ivar)%via(i,j,k,n)
                enddo
              enddo
            enddo
          else
            write(6,*)"Error: please use the character- via when calling sub copy."
          endif
        endif

      enddo
    enddo
!!$OMP END PARALLEL DO

  end subroutine copy

  subroutine copyvia(qmoist,nmoist)

    implicit none

    integer,intent(in)    :: nmoist
    type(cell),dimension(nmoist),intent(inout)   :: qmoist
  ! Local variables
    integer   :: i,j,k,n,ii,jj,kk,ivar
    integer   :: ids,ide,jds,jde,kds,kde, &
                 ims,ime,jms,jme,kms,kme, &
                 ics,ice,jcs,jce,kcs,kce,Nx,Ny,Nz,Nf
  ! ------------------------
    call gp%SetIndex(ids,ide,jds,jde,kds,kde, &
                     ims,ime,jms,jme,kms,kme, &
                     ics,ice,jcs,jce,kcs,kce )
    call gp%SetCubeindex(Nf=Nf)

    do ivar=1,nmoist
      do n=1,Nf
        !do k=kcs,kce
          !kk=2*k
        do k=kds,kde
          kk=k
          do j=jcs,jce
            jj=2*j
            do i=ics,ice
              ii=2*i
              qmoist(ivar)%via(i,j,k,n) = qmoist(ivar)%pv(ii,jj,kk,n)
            enddo
          enddo
        enddo
      enddo
    enddo

  end subroutine copyvia

  subroutine zerocell(q,nvar,cha)

    implicit none

    integer,intent(in)   :: nvar
    character(len=3),optional,intent(in)    :: cha
    type(cell),dimension(nvar),intent(out)  :: q
  ! Local variables
    integer  :: i,j,k,n,ivar
  ! -----------------------

!!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,n,ivar)
    do ivar=1,nvar
      do n=1,gp%Nf

        do k=gp%kms,gp%kme
          do j=gp%jms,gp%jme
            do i=gp%ims,gp%ime
              q(ivar)%pv(i,j,k,n)=0.
            enddo
          enddo
        enddo

        if(present(cha))then
          if(cha=='via')then
            do k=gp%kcs,gp%kce
              do j=gp%jcs,gp%jce
                do i=gp%ics,gp%ice
                  q(ivar)%via(i,j,k,n)=0.
                enddo
              enddo
            enddo
          else
            write(6,*)"Error: please use the character- via when calling sub zerocell."
          endif
        endif

      enddo
    enddo
!!$OMP END PARALLEL DO

  end subroutine zerocell

  subroutine allocateTypeCell(Var,NumVar,via)

    implicit none

    integer,intent(in)    :: NumVar
    character(len=3),optional,intent(in)  :: via
    type(cell),dimension(:),pointer    :: Var
    ! Local variables
    integer   :: ivar,ids,ide,jds,jde,kds,kde, &
                 ims,ime,jms,jme,kms,kme
    ! --------------
    allocate(Var(NumVar)) 
    call gp%setIndex( ids,ide,jds,jde,kds,kde,  &
                      ims,ime,jms,jme,kms,kme ) 
    
    do ivar=1,NumVar
      allocate(Var(ivar)%pv (ims:ime,jms:jme,kms:kme,gp%Nf))
      Var(ivar)%pv=0.
    enddo
    if(present(via))then
      if( via=='via')then
        do ivar=1,NumVar
          allocate(Var(ivar)%via(ims:ime,jms:jme,kms:kme,gp%Nf))
          Var(ivar)%via=0.
        enddo
      endif
    else
      do ivar=1,NumVar
        nullify(Var(ivar)%via)
      enddo
    endif
    return

  end subroutine allocateTypeCell

  subroutine freeTypeCell(Var)

    implicit none

    type(cell),dimension(:),pointer    :: Var
    ! Local variables
    integer   :: ivar,dimsize
    ! --------------
    dimsize=size(Var)

    if(associated(Var))then
      do ivar=1,dimsize
        if(associated(Var(ivar)%pv))then
          deallocate(Var(ivar)%pv)
        endif
      enddo

      do ivar=1,dimsize
        if(associated(Var(ivar)%via))then
          deallocate(Var(ivar)%via)
        endif
      enddo

      deallocate(Var)
    endif
    return

  end subroutine freeTypeCell

  subroutine allocate_domain

    implicit none

    ! Local variables
    integer   :: ivar,ntol
    ! ----------------------------------
    allocate(gp)
    call gp%initialize()

    ! Model prognostic and reference-state variables 
    call allocateTypeCell(q     ,NumVar,via='via')     ! model dynamical variables
    call allocateTypeCell(qn    ,NumVar,via='via')     ! temporary variables storing the values at time level n.
    call allocateTypeCell(q0    ,NumVar,via='via')     ! storage of initial fields
    call allocateTypeCell(qh    ,NumVar,via='via')     ! hydrostatic variables
    call allocateTypeCell(rqsb  ,NumVar,via='via')     ! balance source term
    call allocateTypeCell(qmoist,NumMOI,via='via')     ! moisture variables

    allocate(ph,dpx,dpy)
    allocate(ph %pv(gp%ims:gp%ime,gp%jms:gp%jme,gp%kms:gp%kme,gp%Nf))
    allocate(dpx%pv(gp%ims:gp%ime,gp%jms:gp%jme,gp%kms:gp%kme,gp%Nf))
    allocate(dpy%pv(gp%ims:gp%ime,gp%jms:gp%jme,gp%kms:gp%kme,gp%Nf))

    ! Computational coordinates
    allocate(coord)
    allocate(coord%x (gp%ims:gp%ime))
    allocate(coord%y (gp%jms:gp%jme))
    allocate(coord%z (gp%kms:gp%kme))
    allocate(coord%xc(gp%ics:gp%ice))
    allocate(coord%yc(gp%jcs:gp%jce))
    allocate(coord%zc(gp%kcs:gp%kce))

    ! Jacobian and metric
    allocate(jm)
    allocate(jm%jac(gp%ims:gp%ime,gp%jms:gp%jme));jm%jac=0.
    allocate(jm%metric(gp%ims:gp%ime,gp%jms:gp%jme,6));jm%metric=0.
    allocate(jm%jacv(gp%ims:gp%ime,gp%jms:gp%jme,gp%kms:gp%kme,gp%Nf));jm%jacv=0.
    allocate(jm%metricv(gp%ims:gp%ime,gp%jms:gp%jme,gp%kms:gp%kme,gp%Nf,4));jm%metricv=0.
    allocate(jm%jab(gp%ims:gp%ime,gp%jms:gp%jme,gp%kms:gp%kme,gp%Nf));jm%jab=0.

    ! Ghost cell location information
    allocate(gl)
    ntol=(gp%npvx*2*gp%xhalo*2+gp%npvy*2*gp%yhalo*2)*6   ! Extension of two ghost cells each side. 
    allocate(gl%bci(10,ntol));gl%bci=0.
    allocate(gl%bcr(ntol));gl%bcr=0.

    ! Topography height and geomentric height
    allocate(gm)
    allocate(gm%zs(gp%ims:gp%ime,gp%jms:gp%jme,gp%Nf));gm%zs=0.
    allocate(gm%zz(gp%ims:gp%ime,gp%jms:gp%jme,gp%kms:gp%kme,gp%Nf));gm%zz=0.

  end subroutine allocate_domain

  subroutine solver_domain

    call allocate_domain
    call SetComputationalSpace

  end subroutine solver_domain

END MODULE domain

