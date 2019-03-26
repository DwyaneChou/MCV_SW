MODULE postprocess

  USE global     ,only: np,p0,rd,r2d,d2r,gamma,gra,cp,PI,ra
  USE domain     ,only: jame,gmap,gloc,cord,cell,pv2d,pv3d,coord,gp
  USE projection ,only: pprosp2p,pprop2sp,contravprop2sp
  !USE bc

  implicit none

contains

  subroutine grads(q,jm,gm,coord,nvar,filename,fileid)

    implicit none

    integer,intent(in)    :: nvar,fileid
    character(len=30),intent(in)  :: filename
    type(jame),intent(in) :: jm
    type(gmap),intent(in) :: gm
    type(cord),intent(in) :: coord
    type(cell),dimension(nvar),intent(in)  :: q
  ! Local variables  
    integer   :: i,j,k,n,ivar,Nll,Ntt,nlev
    integer   :: ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,Nf
    logical   :: lex
    type(pv2d),dimension(nvar)   :: q_2d
    type(pv3d),dimension(nvar)   :: q_3d
    type(pv2d)   :: pressure,zs,qt_2d
    type(pv3d)   :: geopotential
    real,dimension(:,:),pointer   :: qplot
    real,dimension(:),pointer     :: plevs,lambda,theta
    real      :: lon0,lat0,dlambda,dtheta,u,v,us,vs,lam,phi
  ! ---------------------------------------------
    nlev=26
    Nll=gp%Nll
    Ntt=gp%Ntt
    allocate(qplot(Nll,Ntt),plevs(nlev),lambda(Nll),theta(Ntt))

    call gp%SetIndex(ids,ide,jds,jde,kds,kde, &
                     ims,ime,jms,jme,kms,kme )
    call gp%SetCubeindex(Nf=Nf)

    do ivar=1,nvar
      allocate(q_2d(ivar)%pv2d(ims:ime,jms:jme,Nf));q_2d(ivar)%pv2d=0.
    enddo
    do ivar=1,nvar
      allocate(q_3d(ivar)%pv3d(ims:ime,jms:jme,1:nlev,Nf));q_3d(ivar)%pv3d=0.
    enddo

    allocate(pressure%pv2d(ims:ime,jms:jme,Nf),geopotential%pv3d(ims:ime,jms:jme,1:nlev,Nf),zs%pv2d(ims:ime,jms:jme,Nf))
    allocate(qt_2d%pv2d(ims:ime,jms:jme,Nf))


    do n=1,Nf
      do j=jms,jme
        do i=ims,ime
          zs%pv2d(i,j,n)=gm%zs(i,j,n)   ! topography

          !call pprop2sp(lam,phi,coord%x(i),coord%y(j),n)
          !u=jm%metricv(i,j,1,n,1)  !  G_v^{13}
          !v=jm%metricv(i,j,1,n,2)  !  G_v^{23}
          !call contravprop2sp(us,vs,u,v,n,lam,phi)
          !zs%pv2d(i,j,n)=vs   
        enddo
      enddo
    enddo

    !call result2d_isobar(q_2d,q,jm,nvar,850.E2)    ! primitive variables
    !call geopotentialheight(geopotential,q,jm,gm,nvar,500.E2)
    !call output_isobar_2d(qplot,q_2d,Nll,Ntt)

    plevs(1) =1000.E2
    plevs(2) =962.5E2
    plevs(3) =925.0E2
    plevs(4) =887.5E2
    plevs(5) =850.0E2
    plevs(6) =800.0E2
    plevs(7) =750.0E2
    plevs(8) =700.0E2
    plevs(9) =650.0E2
    plevs(10)=600.0E2
    plevs(11)=550.0E2
    plevs(12)=500.0E2
    plevs(13)=450.0E2
    plevs(14)=400.0E2
    plevs(15)=350.0E2
    plevs(16)=300.0E2
    plevs(17)=250.0E2
    plevs(18)=200.0E2
    plevs(19)=175.0E2
    plevs(20)=150.0E2
    plevs(21)=125.0E2
    plevs(22)=100.0E2
    plevs(23)=50.0E2
    plevs(24)=30.0E2
    plevs(25)=20.0E2
    plevs(26)=10.0E2



    inquire(file=trim(filename),exist=lex)
    if(lex)then
    !  open(fileid,file=trim(filename),form='unformatted',access='sequential',position='append')
    endif

    call result2d_isobar(q_3d,q,jm,gm,plevs,nvar,nlev)
    ! RHO
    do k=1,nlev
      do n=1,Nf
        do j=jds,jde
          do i=ids,ide
            qt_2d%pv2d(i,j,n)=q_3d(1)%pv3d(i,j,k,n)
          enddo
        enddo
      enddo
      call output_lonlat_2d(qplot,lambda,theta,qt_2d,Nll,Ntt)
      write(fileid)((real(qplot(i,j),4),i=1,Nll),j=1,Ntt)
    enddo
    ! U
    do k=1,nlev
      do n=1,Nf
        do j=jds,jde
          do i=ids,ide
            qt_2d%pv2d(i,j,n)=q_3d(2)%pv3d(i,j,k,n)
          enddo
        enddo
      enddo
      call output_lonlat_2d(qplot,lambda,theta,qt_2d,Nll,Ntt)
      write(fileid)((real(qplot(i,j),4),i=1,Nll),j=1,Ntt)
    enddo
    ! V
    do k=1,nlev
      do n=1,Nf
        do j=jds,jde
          do i=ids,ide
            qt_2d%pv2d(i,j,n)=q_3d(3)%pv3d(i,j,k,n)
          enddo
        enddo
      enddo
      call output_lonlat_2d(qplot,lambda,theta,qt_2d,Nll,Ntt)
      write(fileid)((real(qplot(i,j),4),i=1,Nll),j=1,Ntt)
    enddo

    ! W
    do k=1,nlev
      do n=1,Nf
        do j=jds,jde
          do i=ids,ide
            qt_2d%pv2d(i,j,n)=q_3d(4)%pv3d(i,j,k,n)
          enddo
        enddo
      enddo
      call output_lonlat_2d(qplot,lambda,theta,qt_2d,Nll,Ntt)
      write(fileid)((real(qplot(i,j),4),i=1,Nll),j=1,Ntt)
    enddo

    ! T
    do k=1,nlev
      do n=1,Nf
        do j=jds,jde
          do i=ids,ide
            qt_2d%pv2d(i,j,n)=q_3d(5)%pv3d(i,j,k,n)
          enddo
        enddo
      enddo
      call output_lonlat_2d(qplot,lambda,theta,qt_2d,Nll,Ntt)
      write(fileid)((real(qplot(i,j),4),i=1,Nll),j=1,Ntt)
    enddo

    ! Geopotential height
    call geopotentialheight(geopotential,q,jm,gm,plevs,nvar,nlev)
    do k=1,nlev
    !  call geopotentialheight(geopotential,q,jm,gm,nvar,plevs(k))
      do n=1,Nf
        do j=jds,jde
          do i=ids,ide
            qt_2d%pv2d(i,j,n)=geopotential%pv3d(i,j,k,n)
          enddo
        enddo
      enddo
      call output_lonlat_2d(qplot,lambda,theta,qt_2d,Nll,Ntt)
      write(fileid)((real(qplot(i,j),4),i=1,Nll),j=1,Ntt)
    enddo

    ! PS
    call surfacepresure(pressure,q,jm,nvar)
    !call diagnose_ps(pressure,q_3d(5),geopotential,log(plevs/100.),gm,nlev)

    call output_lonlat_2d(qplot,lambda,theta,pressure,Nll,Ntt)
    write(fileid)((real(qplot(i,j),4),i=1,Nll),j=1,Ntt)
    call output_lonlat_2d(qplot,lambda,theta,zs,Nll,Ntt)
    write(fileid)((real(qplot(i,j),4),i=1,Nll),j=1,Ntt)
    !close(fileid)

#ifdef CUBE
    dlambda= 2.*PI/Nll*R2D
    dtheta = PI/Ntt*R2D
    lon0   = dlambda/2.
    lat0   = PI/Ntt/2*R2D-90.
#endif
#ifdef LONLAT
    dlambda= gp%dx*R2D/2.
    dtheta = gp%dy*R2D/2.
    lon0   = gp%x_min*R2D
    lat0   = gp%y_min*R2D
#endif
#ifdef CARTESIAN
    dlambda= 90./(Nll-1)   ! any reasonable range in the lonlat grid, eg., 90.
    dtheta = 45./(Ntt-1)
    lon0   = 0.
    lat0   = 0.
#endif

    open(10,file='./post.ctl',status='unknown')
    write(10,*)'dset ^./grads.dat'
    write(10,*)'options sequential big_endian'
    write(10,*)'title MCV output'
    write(10,*)'undef -99999.9'
    write(10,*)'xdef ',Nll, 'linear',lon0, dlambda
    write(10,*)'ydef ',Ntt, 'linear',lat0, dtheta
    write(10,*)'zdef ',nlev, 'levels '
    do k=1,nlev
    write(10,*)plevs(k)/100.
    enddo
    write(10,*)'tdef 10000  linear  12z01may2017 720mn'
    write(10,*)'vars 8'
    write(10,*)'rho',nlev,' 0 density'
    write(10,*)'u',nlev,' 0 u wind'
    write(10,*)'v',nlev,' 0 v wind'
    write(10,*)'w',nlev,' 0 w wind'
    write(10,*)'t',nlev,' 0 temperature'
    write(10,*)'h',nlev,' 0 geopotential height'
    write(10,*)'ps 0 0 surface pressure'
    write(10,*)'zs 0 0 surface pressure'
    write(10,*)'endvars'
    close(10)

    do ivar=1,nvar
      deallocate(q_2d(ivar)%pv2d)
      deallocate(q_3d(ivar)%pv3d)
    enddo
    deallocate(pressure%pv2d,geopotential%pv3d,qt_2d%pv2d,zs%pv2d,qplot,plevs,lambda,theta)

  end subroutine grads

  subroutine gradsmodelvar(q,qh,qmoist,gl,jm,nvar,nmoist,fileid)

    implicit none

    integer,intent(in)    :: nvar,fileid,nmoist
    type(gloc),intent(in) :: gl
    type(jame),intent(in) :: jm
    type(cell),dimension(nvar),intent(in)    :: q,qh
    type(cell),dimension(nmoist),intent(in)  :: qmoist
  ! Local variables  
    integer   :: i,j,k,n,ivar,Nll,Ntt,nlev
    integer   :: ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,Nf
    type(pv3d),dimension(nvar)    :: q_3d,qh_3d,qm_3d
    type(pv2d)                    :: qt_2d,qht_2d
    real,dimension(:,:),pointer   :: qplot,qploth
    real,dimension(:),pointer     :: lambda,theta
  ! ---------------------------------------------
    nlev=gp%Npvz
    Nll=gp%Nll
    Ntt=gp%Ntt
    allocate(qplot(Nll,Ntt),qploth(Nll,Ntt))

    call gp%SetIndex(ids,ide,jds,jde,kds,kde, &
                     ims,ime,jms,jme,kms,kme )
    call gp%SetCubeindex(Nf=Nf)

    allocate(qt_2d%pv2d(ims:ime,jms:jme,Nf));qt_2d%pv2d=0.
    allocate(qht_2d%pv2d(ims:ime,jms:jme,Nf));qht_2d%pv2d=0.
    do ivar=1,nvar
      allocate(q_3d(ivar)%pv3d(ims:ime,jms:jme,1:nlev,Nf));q_3d(ivar)%pv3d=0.
      allocate(qh_3d(ivar)%pv3d(ims:ime,jms:jme,1:nlev,Nf));qh_3d(ivar)%pv3d=0.
      allocate(qm_3d(ivar)%pv3d(ims:ime,jms:jme,1:nlev,Nf));qm_3d(ivar)%pv3d=0.
    enddo
    allocate(lambda(Nll),theta(Ntt))

    !  Model primitive variable 
    call result2d_modelvar(q_3d,q,jm,nvar)
    call result2d_modelvar(qh_3d,qh,jm,nvar)    ! primitive variables

    ! rho
    do k=1,nlev
      do n=1,Nf
        do j=jds,jde
          do i=ids,ide
            qt_2d%pv2d(i,j,n)=q_3d(1)%pv3d(i,j,k,n)
          enddo
        enddo
      enddo
      call output_lonlat_2d(qplot,lambda,theta,qt_2d,Nll,Ntt)

#ifdef CARTESIAN
      do n=1,Nf
        do j=jds,jde
          do i=ids,ide
            qht_2d%pv2d(i,j,n)=qh_3d(1)%pv3d(i,j,k,n)
          enddo
        enddo
      enddo
      call output_lonlat_2d(qploth,lambda,theta,qht_2d,Nll,Ntt)

      write(fileid)((real(qplot(i,j)-qploth(i,j),4),i=1,Nll),j=1,Ntt)
#endif
#ifdef CUBE
      write(fileid)((real(qplot(i,j),4),i=1,Nll),j=1,Ntt)
#endif
#ifdef LONLAT
      write(fileid)((real(qplot(i,j),4),i=1,Nll),j=1,Ntt)
#endif
    enddo

    ! U
    do k=1,nlev
      do n=1,Nf
        do j=jds,jde
          do i=ids,ide
            qt_2d%pv2d(i,j,n)=q_3d(2)%pv3d(i,j,k,n)
          enddo
        enddo
      enddo
      call output_lonlat_2d(qplot,lambda,theta,qt_2d,Nll,Ntt)
      write(fileid)((real(qplot(i,j),4),i=1,Nll),j=1,Ntt)
    enddo

    ! V
    do k=1,nlev
      do n=1,Nf
        do j=jds,jde
          do i=ids,ide
            qt_2d%pv2d(i,j,n)=q_3d(3)%pv3d(i,j,k,n)
          enddo
        enddo
      enddo
      call output_lonlat_2d(qplot,lambda,theta,qt_2d,Nll,Ntt)
      write(fileid)((real(qplot(i,j),4),i=1,Nll),j=1,Ntt)
    enddo

    ! W
    do k=1,nlev
      do n=1,Nf
        do j=jds,jde
          do i=ids,ide
            qt_2d%pv2d(i,j,n)=q_3d(4)%pv3d(i,j,k,n)
          enddo
        enddo
      enddo
      call output_lonlat_2d(qplot,lambda,theta,qt_2d,Nll,Ntt)
      write(fileid)((real(qplot(i,j),4),i=1,Nll),j=1,Ntt)
    enddo

    ! T
    do k=1,nlev
      do n=1,Nf
        do j=jds,jde
          do i=ids,ide
            qt_2d%pv2d(i,j,n)=q_3d(5)%pv3d(i,j,k,n)
          enddo
        enddo
      enddo
      call output_lonlat_2d(qplot,lambda,theta,qt_2d,Nll,Ntt)

#ifdef CARTESIAN
      do n=1,Nf
        do j=jds,jde
          do i=ids,ide
            qht_2d%pv2d(i,j,n)=qh_3d(5)%pv3d(i,j,k,n)
          enddo
        enddo
      enddo
      call output_lonlat_2d(qploth,lambda,theta,qht_2d,Nll,Ntt)
      write(fileid)((real(qplot(i,j)-qploth(i,j),4),i=1,Nll),j=1,Ntt)
#endif
#ifdef CUBE
      write(fileid)((real(qplot(i,j),4),i=1,Nll),j=1,Ntt)
#endif
#ifdef LONLAT
      write(fileid)((real(qplot(i,j),4),i=1,Nll),j=1,Ntt)
#endif
    enddo

    !  Moist
    call result2d_modelmoist(qm_3d,qmoist,q,jm,nvar,nmoist)
    do k=1,nlev
      do n=1,Nf
        do j=jds,jde
          do i=ids,ide
            qt_2d%pv2d(i,j,n)=qm_3d(1)%pv3d(i,j,k,n)
          enddo
        enddo
      enddo
      call output_lonlat_2d(qplot,lambda,theta,qt_2d,Nll,Ntt)
      write(fileid)((real(qplot(i,j),4),i=1,Nll),j=1,Ntt)
    enddo

    do k=1,nlev
      do n=1,Nf
        do j=jds,jde
          do i=ids,ide
            qt_2d%pv2d(i,j,n)=qm_3d(2)%pv3d(i,j,k,n)
          enddo
        enddo
      enddo
      call output_lonlat_2d(qplot,lambda,theta,qt_2d,Nll,Ntt)
      write(fileid)((real(qplot(i,j),4),i=1,Nll),j=1,Ntt)
    enddo

    do k=1,nlev
      do n=1,Nf
        do j=jds,jde
          do i=ids,ide
            qt_2d%pv2d(i,j,n)=qm_3d(3)%pv3d(i,j,k,n)
          enddo
        enddo
      enddo
      call output_lonlat_2d(qplot,lambda,theta,qt_2d,Nll,Ntt)
      write(fileid)((real(qplot(i,j),4),i=1,Nll),j=1,Ntt)
    enddo

    do k=1,nlev
      do n=1,Nf
        do j=jds,jde
          do i=ids,ide
            qt_2d%pv2d(i,j,n)=qm_3d(4)%pv3d(i,j,k,n)
          enddo
        enddo
      enddo
      call output_lonlat_2d(qplot,lambda,theta,qt_2d,Nll,Ntt)
      write(fileid)((real(qplot(i,j),4),i=1,Nll),j=1,Ntt)
    enddo

    do k=1,nlev
      do n=1,Nf
        do j=jds,jde
          do i=ids,ide
            qt_2d%pv2d(i,j,n)=qm_3d(5)%pv3d(i,j,k,n)
          enddo
        enddo
      enddo
      call output_lonlat_2d(qplot,lambda,theta,qt_2d,Nll,Ntt)
      write(fileid)((real(qplot(i,j),4),i=1,Nll),j=1,Ntt)
    enddo

    open(10,file='./model.ctl',status='unknown')
    write(10,*)'dset ^./gradsmodelvar.dat'
    write(10,*)'options sequential big_endian'
    write(10,*)'title MCV output'
    write(10,*)'undef -99999.9'
    write(10,*)'xdef ',Nll, 'linear',PI/Nll*R2D, 2*PI/Nll*R2D
    write(10,*)'ydef ',Ntt, 'linear',PI/Ntt/2*R2D-90.,PI/Ntt*R2D
    write(10,*)'zdef ',nlev, 'linear 1 ', gp%Npvz
    write(10,*)'tdef 100  linear  12z01may2017 720mn'
    write(10,*)'vars 10'
    write(10,*)'rho  ',gp%Npvz,' 0 density'
    write(10,*)'u  ',gp%Npvz,' 0 u wind'
    write(10,*)'v  ',gp%Npvz,' 0 v wind'
    write(10,*)'w  ',gp%Npvz,' 0 w wind'
    write(10,*)'t  ',gp%Npvz,' 0 temperature'
    write(10,*)'q1  ',gp%Npvz,' 0 temperature'
    write(10,*)'q2  ',gp%Npvz,' 0 temperature'
    write(10,*)'q3  ',gp%Npvz,' 0 temperature'
    write(10,*)'q4  ',gp%Npvz,' 0 temperature'
    write(10,*)'q5  ',gp%Npvz,' 0 temperature'
    write(10,*)'endvars'
    close(10)

    do ivar=1,nvar
      deallocate(q_3d(ivar)%pv3d,qh_3d(ivar)%pv3d,qm_3d(ivar)%pv3d)
    enddo
    deallocate(qt_2d%pv2d)
    deallocate(qplot,qploth,lambda,theta)

  end subroutine gradsmodelvar

  subroutine output_patch(q,gl,jm,coord,nvar,npatch,start)

    implicit none

    integer,intent(in)    :: nvar,npatch,start
    type(gloc),intent(in) :: gl
    type(jame),intent(in) :: jm
    type(cord),intent(in) :: coord
    type(cell),dimension(nvar),intent(in)  :: q
  ! Local variables
    integer   :: i,j,k,ivar
    integer   :: ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,Nf
    real      :: lambda,theta
  ! --------------------------
    call gp%SetIndex(ids,ide,jds,jde,kds,kde, &
                     ims,ime,jms,jme,kms,kme )
    call gp%SetCubeindex(Nf=Nf)

    if (start == 0) then
    open(unit=101,file='Result_patch.plt',form='FORMATTED',status='UNKNOWN')
    else
    open(unit=101,file='Result_patch',form='FORMATTED',status='UNKNOWN',position='APPEND')
    end if
    write(101,*)'TITLE = ADVECTION CUBED-SPHERE'
    write (101,*) 'VARIABLES = "LAMBDA", "THETA", "RHO","RHOU","RHOV", "RHOW", "RHOT"'
    write (101,*) 'ZONE I=',ide-ids+1, ' J=',jde-jds+1,' DATAPACKING=POINT'

    !k=10
    k=5
    do j=jds,jde
      do i=ids,ide

        call pprop2sp(lambda,theta,coord%x(i),coord%y(j),npatch)
        !write(101,*)lambda*R2D,theta*R2D,                  &
        write(101,*)i,j,&
            q(1)%pv(i,j,k,npatch)/jm%jab(i,j,k,npatch),    &
            q(2)%pv(i,j,k,npatch)/jm%jab(i,j,k,npatch)*Ra, &
            q(3)%pv(i,j,k,npatch)/jm%jab(i,j,k,npatch)*Ra, &
            q(4)%pv(i,j,k,npatch)/jm%jab(i,j,k,npatch),    &
            q(5)%pv(i,j,k,npatch)/jm%jab(i,j,k,npatch)

      end do
    end do

    close(101)

  end subroutine output_patch

  subroutine output(q,gl,jm,nvar)

    implicit none

    integer,intent(in)    :: nvar
    type(gloc),intent(in) :: gl
    type(jame),intent(in) :: jm
    type(cell),dimension(nvar),intent(in)  :: q
  ! Local variables  
    integer   :: i,j,k,ivar
    integer   :: ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,Nf
    type(pv2d),dimension(nvar)   :: q_2d
    type(pv2d)   :: pressure,geopotential
  ! ---------------------------------------------
    call gp%SetIndex(ids,ide,jds,jde,kds,kde, &
                     ims,ime,jms,jme,kms,kme )
    call gp%SetCubeindex(Nf=Nf)

    do ivar=1,nvar
      allocate(q_2d(ivar)%pv2d(ims:ime,jms:jme,Nf));q_2d(ivar)%pv2d=0.
    enddo
    allocate(pressure%pv2d(ims:ime,jms:jme,Nf),geopotential%pv2d(ims:ime,jms:jme,Nf))

!    call result2d_isobar(q_2d,q,jm,nvar,850.E2)    ! primitive variables
    call surfacepresure(pressure,q,jm,nvar)
!    call geopotentialheight(geopotential,q,jm,gm,nvar,500.E2)

    ! output result
    call output_2d_modelvariables(q_2d,5,0)
    call output_2d_one(geopotential,1,0)
    call output_2d_one(pressure,2,0)

    ! output errors


    do ivar=1,nvar
      deallocate(q_2d(ivar)%pv2d)
    enddo
    deallocate(pressure%pv2d,geopotential%pv2d)

  end subroutine output


! -----------------------  LONLAT  -------------------------------
#ifdef LONLAT
  subroutine output_2d_modelvariables(q_2d,nvar,outtype)

    implicit none

    integer,intent(in)    :: nvar,outtype
    type(pv2d),dimension(nvar),intent(in) :: q_2d
  ! Local varaibles
    integer   :: i,j,m,n,Nll,Ntt,ii,jj,kk,Nx,Ny,ivar,ig,jg,Nf
    integer   :: ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme
    real      :: dlambda,dtheta,xc,yc,dx,dy,a
    real,dimension(:),pointer       :: lambda,theta
    real,dimension(:,:,:),pointer   :: qplot
    real,dimension(:,:),pointer     :: qplot_tmp
    type(pv2d),pointer              :: q_2d_tmp
  ! -----------------------------------
    call gp%SetIndex(ids,ide,jds,jde,kds,kde, &
                     ims,ime,jms,jme,kms,kme )
    call gp%SetCubeindex(Nf=Nf)

    a =Ra
    Nll=gp%Nll
    Ntt=gp%Ntt
    allocate(lambda(Nll),theta(Ntt),qplot(Nll,Ntt,nvar),qplot_tmp(Nll,Ntt))
    allocate(q_2d_tmp)
    allocate(q_2d_tmp%pv2d(ims:ime,jms:jme,gp%Nf))

    do ivar=1,nvar

      do n=1,Nf
        do j=jms,jme
          do i=ims,ime
            q_2d_tmp%pv2d(i,j,n)=q_2d(ivar)%pv2d(i,j,n)
          enddo
        enddo
      enddo

      call output_lonlat_2d(qplot_tmp,lambda,theta,q_2d_tmp,Nll,Ntt)
      qplot(:,:,ivar)=qplot_tmp(:,:)

    enddo


    if(outtype == 0) then
      open(unit=101,file='./output/result.plt',form='FORMATTED',status='UNKNOWN')
    else
      open(unit=101,file='./output/error.plt',form='FORMATTED',status='UNKNOWN')
    endif

    write(101,*)'MCV OUTPUT'
    write (101,*) 'VARIABLES = "LAMBDA", "THETA", "RHO", "U", "V", "W", "T"'
    write (101,*) 'ZONE I=',Nll, ' J=',Ntt,' DATAPACKING=POINT' 

    do j=1,Ntt
      do i=1,Nll

        write(101,*)lambda(i)/D2R,theta(j)/D2R, &
        qplot(i,j,1),qplot(i,j,2),qplot(i,j,3),qplot(i,j,4),qplot(i,j,5)

      end do
    end do

    close(101)

    deallocate(lambda,theta,qplot,qplot_tmp,q_2d_tmp)

  end subroutine output_2d_modelvariables

  subroutine output_2d_one(q_2d,outputtype,outtype)

    implicit none

    integer,intent(in)    :: outputtype,outtype
    type(pv2d),intent(in) :: q_2d
  ! Local varaibles
    integer   :: i,j,m,n,Nll,Ntt,ii,jj,kk,Nx,Ny,ivar,ig,jg
    integer   :: ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme
    real      :: dlambda,dtheta,xc,yc,dx,dy,a
    real,dimension(:),pointer       :: lambda,theta,temp,x,y
    real,dimension(:,:),pointer     :: qi
    real,dimension(:,:),pointer     :: qplot
  ! -----------------------------------
    call gp%SetIndex(ids,ide,jds,jde,kds,kde, &
                     ims,ime,jms,jme,kms,kme )

    a =Ra
    Nll=gp%Nll
    Ntt=gp%Ntt
    allocate(lambda(Nll),theta(Ntt),qplot(Nll,Ntt))

    call output_lonlat_2d(qplot,lambda,theta,q_2d,Nll,Ntt)

    if (outputtype == 1) then
      if ( outtype == 0 ) then
        open(unit=101,file='./output/height.plt',form='FORMATTED',status='UNKNOWN')
      else
        open(unit=101,file='./output/height_error.plt',form='FORMATTED',status='UNKNOWN')
      end if
      write(101,*)'MCV OUTPUT'
      write (101,*) 'VARIABLES = "LAMBDA", "THETA", "HEIGHT"'
      write (101,*) 'ZONE I=',Nll, ' J=',Ntt,' DATAPACKING=POINT' 
    else
      if ( outtype == 0 ) then
        open(unit=101,file='./output/pressure.plt',form='FORMATTED',status='UNKNOWN')
      else
        open(unit=101,file='./output/pressure_error.plt',form='FORMATTED',status='UNKNOWN')
      end if

      write(101,*)'MCV OUTPUT'
      write (101,*) 'VARIABLES = "LAMBDA", "THETA", "PS"'
      write (101,*) 'ZONE I=',Nll, ' J=',Ntt,' DATAPACKING=POINT' 

    end if

    do j=1,Ntt
      do i=1,Nll

        write(101,*)lambda(i)/D2R,theta(j)/D2R,qplot(i,j)

      end do
    end do

    close(101)

    deallocate(lambda,theta,qplot)

  end subroutine output_2d_one

  subroutine output_lonlat_2d(qplot,lambda,theta,q_2d,Nll,Ntt)

    implicit none

    integer,intent(in)           :: Nll,Ntt
    type(pv2d),intent(in)        :: q_2d
    real,dimension(Nll),intent(out)        :: lambda
    real,dimension(Ntt),intent(out)        :: theta
    real,dimension(Nll,Ntt),intent(out)    :: qplot
  ! Local varaibles
    integer   :: i,j,m,n,ii,jj,kk,Nx,Ny,ivar,ig,jg
    integer   :: ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme
    real      :: dlambda,dtheta,xc,yc,dx,dy,a,xstart,ystart
  ! -----------------------------------
    call gp%SetIndex(ids,ide,jds,jde,kds,kde, &
                     ims,ime,jms,jme,kms,kme )

    a =Ra
    dx=gp%dx
    dy=gp%dy
    Nx =gp%Nx
    Ny =gp%Ny

    dlambda= dx/2.
    dtheta = dy/2.
    xstart = gp%x_min
    ystart = gp%y_min

    do i=1,Nll
      lambda(i)= (i-1)*dlambda + xstart
    enddo
    do j=1,Ntt
      theta(j) = (j-1)*dtheta + ystart
    enddo

    kk=1
    do j=1,Ntt
      do i=1,Nll
        qplot(i,j)=q_2d%pv2d(i,j,kk)
      enddo
    enddo

  end subroutine output_lonlat_2d

#endif   ! End of LONLAT

! -----------------------  CUBE  -------------------------------

#ifdef CUBE
  subroutine output_2d_modelvariables(q_2d,nvar,outtype)

    implicit none

    integer,intent(in)    :: nvar,outtype
    type(pv2d),dimension(nvar),intent(in) :: q_2d
  ! Local varaibles
    integer   :: i,j,m,n,Nll,Ntt,ii,jj,kk,Nx,Ny,ivar,ig,jg,Nf
    integer   :: ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme
    real      :: dlambda,dtheta,xc,yc,dx,dy,a
    real,dimension(:),pointer       :: lambda,theta
    real,dimension(:,:,:),pointer   :: qplot
    real,dimension(:,:),pointer     :: qplot_tmp
    type(pv2d),pointer              :: q_2d_tmp
  ! -----------------------------------
    call gp%SetIndex(ids,ide,jds,jde,kds,kde, &
                     ims,ime,jms,jme,kms,kme )
    call gp%SetCubeindex(Nf=Nf)

    a =Ra
    Nll=gp%Nll
    Ntt=gp%Ntt
    allocate(lambda(Nll),theta(Ntt),qplot(Nll,Ntt,nvar),qplot_tmp(Nll,Ntt))
    allocate(q_2d_tmp)
    allocate(q_2d_tmp%pv2d(ims:ime,jms:jme,gp%Nf))

    do ivar=1,nvar

      do n=1,Nf
        do j=jms,jme
          do i=ims,ime
            q_2d_tmp%pv2d(i,j,n)=q_2d(ivar)%pv2d(i,j,n)
          enddo
        enddo
      enddo

      call output_lonlat_2d(qplot_tmp,lambda,theta,q_2d_tmp,Nll,Ntt)
      qplot(:,:,ivar)=qplot_tmp(:,:)

    enddo


    if(outtype == 0) then
      open(unit=101,file='./output/result.plt',form='FORMATTED',status='UNKNOWN')
    else
      open(unit=101,file='./output/error.plt',form='FORMATTED',status='UNKNOWN')
    endif

    write(101,*)'MCV OUTPUT'
    write (101,*) 'VARIABLES = "LAMBDA", "THETA", "RHO", "U", "V", "W", "T"'
    write (101,*) 'ZONE I=',Nll, ' J=',Ntt,' DATAPACKING=POINT' 

    do j=1,Ntt
      do i=1,Nll

        write(101,*)lambda(i)/D2R,theta(j)/D2R, &
        qplot(i,j,1),qplot(i,j,2),qplot(i,j,3),qplot(i,j,4),qplot(i,j,5)

      end do
    end do

    close(101)

    deallocate(lambda,theta,qplot,qplot_tmp,q_2d_tmp)

  end subroutine output_2d_modelvariables

  subroutine output_2d_one(q_2d,outputtype,outtype)

    implicit none

    integer,intent(in)    :: outputtype,outtype
    type(pv2d),intent(in) :: q_2d
  ! Local varaibles
    integer   :: i,j,m,n,Nll,Ntt,ii,jj,kk,Nx,Ny,ivar,ig,jg
    integer   :: ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme
    real      :: dlambda,dtheta,xc,yc,dx,dy,a
    real,dimension(:),pointer       :: lambda,theta,temp,x,y
    real,dimension(:,:),pointer     :: qi
    real,dimension(:,:),pointer     :: qplot
  ! -----------------------------------
    call gp%SetIndex(ids,ide,jds,jde,kds,kde, &
                     ims,ime,jms,jme,kms,kme )

    a =Ra
    Nll=gp%Nll
    Ntt=gp%Ntt
    allocate(lambda(Nll),theta(Ntt),qplot(Nll,Ntt))

    call output_lonlat_2d(qplot,lambda,theta,q_2d,Nll,Ntt)

    if (outputtype == 1) then
      if ( outtype == 0 ) then
        open(unit=101,file='./output/height.plt',form='FORMATTED',status='UNKNOWN')
      else
        open(unit=101,file='./output/height_error.plt',form='FORMATTED',status='UNKNOWN')
      end if
      write(101,*)'MCV OUTPUT'
      write (101,*) 'VARIABLES = "LAMBDA", "THETA", "HEIGHT"'
      write (101,*) 'ZONE I=',Nll, ' J=',Ntt,' DATAPACKING=POINT' 
    else
      if ( outtype == 0 ) then
        open(unit=101,file='./output/pressure.plt',form='FORMATTED',status='UNKNOWN')
      else
        open(unit=101,file='./output/pressure_error.plt',form='FORMATTED',status='UNKNOWN')
      end if

      write(101,*)'MCV OUTPUT'
      write (101,*) 'VARIABLES = "LAMBDA", "THETA", "PS"'
      write (101,*) 'ZONE I=',Nll, ' J=',Ntt,' DATAPACKING=POINT' 

    end if

    do j=1,Ntt
      do i=1,Nll

        write(101,*)lambda(i)/D2R,theta(j)/D2R,qplot(i,j)

      end do
    end do

    close(101)

    deallocate(lambda,theta,qplot)

  end subroutine output_2d_one

  subroutine output_lonlat_2d(qplot,lambda,theta,q_2d,Nll,Ntt)

    implicit none

    integer,intent(in)           :: Nll,Ntt
    type(pv2d),intent(in)        :: q_2d
    real,dimension(Nll),intent(out)        :: lambda
    real,dimension(Ntt),intent(out)        :: theta
    real,dimension(Nll,Ntt),intent(out)    :: qplot
  ! Local varaibles
    integer   :: i,j,m,n,ii,jj,kk,Nx,Ny,ivar,ig,jg
    integer   :: ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme
    real      :: dlambda,dtheta,xc,yc,dx,dy,a,x1,x2,y1,y2,q1,q2,q3,q4,qx1,qx2
    real,dimension(:),pointer       :: temp,x,y
    real,dimension(:,:),pointer     :: qi
  ! -----------------------------------
    call gp%SetIndex(ids,ide,jds,jde,kds,kde, &
                     ims,ime,jms,jme,kms,kme )

    a =Ra
    dx=gp%dx
    dy=gp%dy
    Nx =gp%Nx
    Ny =gp%Ny

    allocate(qi(NP,NP),temp(NP))
    allocate(x(NP),y(NP))

    dlambda=2.*PI/Nll
    dtheta =PI/Ntt

    do i=1,Nll
      !lambda(i)=-PI/4.+(i-1)*dlambda + dlambda/2.
      lambda(i)= (i-1)*dlambda + dlambda/2.
    enddo
    do j=1,Ntt
      theta(j)=-PI/2.+(j-1)*dtheta + dtheta/2.
    enddo

    do j=1,Ntt
      do i=1,Nll
        kk=1
        if (lambda(i) >    pi/4. .and. lambda(i)< 3.*pi/4.) kk=2
        if (lambda(i) > 3.*pi/4. .and. lambda(i)< 5.*pi/4.) kk=3
        if (lambda(i) > 5.*pi/4. .and. lambda(i)< 7.*pi/4.) kk=4

        call pprosp2p(xc,yc,lambda(i),theta(j),kk)

        if (yc < -pi/4.) then; kk=6; call pprosp2p(xc,yc,lambda(i),theta(j),kk); end if
        if (yc >  pi/4.) then; kk=5; call pprosp2p(xc,yc,lambda(i),theta(j),kk); end if

#ifdef xxx     
        ii=int((xc+pi/4.)/dx)+1; ii=max(min(ii,Nx),1)
        jj=int((yc+pi/4.)/dy)+1; jj=max(min(jj,Ny),1)

        call gp%setpvindex(ig,jg,ii,jj)
        do n=1,NP
          y(n)=coord%y(jg+(n-1))
          do m=1,NP
            x(m)=coord%x(ig+(m-1))
            qi(m,n)=q_2d%pv2d(ig+(m-1),jg+(n-1),kk)
          enddo
        enddo

        do n=1,NP
          temp(n)=interplotant_1d_post(qi(:,n),xc,x)
        enddo
        qplot(i,j)=interplotant_1d_post(temp,yc,y)
#else
        ii=int((xc+pi/4.)/(dx/2.))+1
        jj=int((yc+pi/4.)/(dx/2.))+1
        x1=coord%x(ii);x2=coord%x(ii+1)
        y1=coord%y(jj);y2=coord%y(jj+1)
        q1=q_2d%pv2d(ii  ,jj,kk)
        q2=q_2d%pv2d(ii+1,jj,kk)
        q3=q_2d%pv2d(ii  ,jj+1,kk)
        q4=q_2d%pv2d(ii+1,jj+1,kk)

        qx1=linear_1d(q1,q2,xc,x1,x2)
        qx2=linear_1d(q3,q4,xc,x1,x2)

        qplot(i,j)=linear_1d(qx1,qx2,yc,y1,y2)
#endif

      enddo
    enddo

    deallocate(qi,temp,x,y)

  end subroutine output_lonlat_2d

#endif    ! End of CUBE


! -----------------------  CARTESIAN  -------------------------------
#ifdef CARTESIAN

  subroutine output_2d_modelvariables(q_2d,nvar,outtype)

    implicit none

    integer,intent(in)    :: nvar,outtype
    type(pv2d),dimension(nvar),intent(in) :: q_2d
  ! Local varaibles
    integer   :: i,j,m,n,Nll,Ntt,ii,jj,kk,Nx,Ny,ivar,ig,jg,Nf
    integer   :: ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme
    real      :: dlambda,dtheta,xc,yc,dx,dy,a
    real,dimension(:),pointer       :: lambda,theta
    real,dimension(:,:,:),pointer   :: qplot
    real,dimension(:,:),pointer     :: qplot_tmp
    type(pv2d),pointer              :: q_2d_tmp
  ! -----------------------------------
    call gp%SetIndex(ids,ide,jds,jde,kds,kde, &
                     ims,ime,jms,jme,kms,kme )
    call gp%SetCubeindex(Nf=Nf)

    a =Ra
    Nll=gp%Nll
    Ntt=gp%Ntt
    allocate(lambda(Nll),theta(Ntt),qplot(Nll,Ntt,nvar),qplot_tmp(Nll,Ntt))
    allocate(q_2d_tmp)
    allocate(q_2d_tmp%pv2d(ims:ime,jms:jme,gp%Nf))

    do ivar=1,nvar

      do n=1,Nf
        do j=jms,jme
          do i=ims,ime
            q_2d_tmp%pv2d(i,j,n)=q_2d(ivar)%pv2d(i,j,n)
          enddo
        enddo
      enddo

      call output_lonlat_2d(qplot_tmp,lambda,theta,q_2d_tmp,Nll,Ntt)
      qplot(:,:,ivar)=qplot_tmp(:,:)

    enddo


    if(outtype == 0) then
      open(unit=101,file='./output/result.plt',form='FORMATTED',status='UNKNOWN')
    else
      open(unit=101,file='./output/error.plt',form='FORMATTED',status='UNKNOWN')
    endif

    write(101,*)'MCV OUTPUT'
    write (101,*) 'VARIABLES = "LAMBDA", "THETA", "RHO", "U", "V", "W", "T"'
    write (101,*) 'ZONE I=',Nll, ' J=',Ntt,' DATAPACKING=POINT' 

    do j=1,Ntt
      do i=1,Nll

        write(101,*)lambda(i)/D2R,theta(j)/D2R, &
        qplot(i,j,1),qplot(i,j,2),qplot(i,j,3),qplot(i,j,4),qplot(i,j,5)

      end do
    end do

    close(101)

    deallocate(lambda,theta,qplot,qplot_tmp,q_2d_tmp)

  end subroutine output_2d_modelvariables

  subroutine output_2d_one(q_2d,outputtype,outtype)

    implicit none

    integer,intent(in)    :: outputtype,outtype
    type(pv2d),intent(in) :: q_2d
  ! Local varaibles
    integer   :: i,j,m,n,Nll,Ntt,ii,jj,kk,Nx,Ny,ivar,ig,jg
    integer   :: ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme
    real      :: dlambda,dtheta,xc,yc,dx,dy,a
    real,dimension(:),pointer       :: lambda,theta,temp,x,y
    real,dimension(:,:),pointer     :: qi
    real,dimension(:,:),pointer     :: qplot
  ! -----------------------------------
    call gp%SetIndex(ids,ide,jds,jde,kds,kde, &
                     ims,ime,jms,jme,kms,kme )

    a =Ra
    Nll=gp%Nll
    Ntt=gp%Ntt
    allocate(lambda(Nll),theta(Ntt),qplot(Nll,Ntt))

    call output_lonlat_2d(qplot,lambda,theta,q_2d,Nll,Ntt)

    if (outputtype == 1) then
      if ( outtype == 0 ) then
        open(unit=101,file='./output/height.plt',form='FORMATTED',status='UNKNOWN')
      else
        open(unit=101,file='./output/height_error.plt',form='FORMATTED',status='UNKNOWN')
      end if
      write(101,*)'MCV OUTPUT'
      write (101,*) 'VARIABLES = "LAMBDA", "THETA", "HEIGHT"'
      write (101,*) 'ZONE I=',Nll, ' J=',Ntt,' DATAPACKING=POINT' 
    else
      if ( outtype == 0 ) then
        open(unit=101,file='./output/pressure.plt',form='FORMATTED',status='UNKNOWN')
      else
        open(unit=101,file='./output/pressure_error.plt',form='FORMATTED',status='UNKNOWN')
      end if

      write(101,*)'MCV OUTPUT'
      write (101,*) 'VARIABLES = "LAMBDA", "THETA", "PS"'
      write (101,*) 'ZONE I=',Nll, ' J=',Ntt,' DATAPACKING=POINT' 

    end if

    do j=1,Ntt
      do i=1,Nll

        write(101,*)lambda(i)/D2R,theta(j)/D2R,qplot(i,j)

      end do
    end do

    close(101)

    deallocate(lambda,theta,qplot)

  end subroutine output_2d_one

  subroutine output_lonlat_2d(qplot,lambda,theta,q_2d,Nll,Ntt)

    implicit none

    integer,intent(in)           :: Nll,Ntt
    type(pv2d),intent(in)        :: q_2d
    real,dimension(Nll),intent(out)        :: lambda
    real,dimension(Ntt),intent(out)        :: theta
    real,dimension(Nll,Ntt),intent(out)    :: qplot
  ! Local varaibles
    integer   :: i,j,m,n,ii,jj,kk,Nx,Ny,ivar,ig,jg
    integer   :: ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme
    real      :: dlambda,dtheta,xc,yc,dx,dy,a,xstart,ystart
  ! -----------------------------------
    call gp%SetIndex(ids,ide,jds,jde,kds,kde, &
                     ims,ime,jms,jme,kms,kme )

    a =Ra
    dx=gp%dx
    dy=gp%dy
    Nx =gp%Nx
    Ny =gp%Ny

    dlambda= dx/2.
    dtheta = dy/2.
    xstart = gp%x_min
    ystart = gp%y_min

    do i=1,Nll
      lambda(i)= (i-1)*dlambda + xstart
    enddo
    do j=1,Ntt
      theta(j) = (j-1)*dtheta + ystart
    enddo

    kk=1
    do j=1,Ntt
      do i=1,Nll

        !do ivar=1,5
        qplot(i,j)=q_2d%pv2d(i,j,kk)
        !enddo

      enddo
    enddo

  end subroutine output_lonlat_2d
#endif     ! End of CARTESIAN


! ########################################################################
!         
!                Common parts ready for outputting
!
! ########################################################################
  subroutine result2d_modelmoist(q_3d,qmoist,q,jm,nvar,nmoist)
!
!   primitive variables: (rho,u,v,w,t) at model level kk.
!
    implicit none

    integer,intent(in)    :: nmoist,nvar
    type(jame),intent(in) :: jm
    type(cell),dimension(nvar),intent(in)    :: q
    type(cell),dimension(nmoist),intent(in)  :: qmoist
    type(pv3d),dimension(nmoist),intent(out) :: q_3d
  ! Local variables
    integer :: i,j,k,n,ivar
    integer :: ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,Nf
  ! -----------------------------------------
    call gp%SetIndex(ids,ide,jds,jde,kds,kde, &
                     ims,ime,jms,jme,kms,kme )
    call gp%SetCubeindex(Nf=Nf)

    do n=1,Nf
      do k=kds,kde
        do j=jds,jde
          do i=ids,ide

            do ivar=1,nmoist
              q_3d(ivar)%pv3d(i,j,k,n)=qmoist(ivar)%pv(i,j,k,n)/q(1)%pv(i,j,k,n)
            enddo

          enddo
        enddo
      enddo
    enddo

  end subroutine result2d_modelmoist

  subroutine result2d_modelvar(q_3d,q,jm,nvar)
!
!   primitive variables: (rho,u,v,w,t) at model level kk.
!
    implicit none

    integer,intent(in)    :: nvar
    type(jame),intent(in) :: jm
    type(cell),dimension(nvar),intent(in)  :: q
    type(pv3d),dimension(nvar),intent(out) :: q_3d
  ! Local variables
    integer :: i,j,k,n,ivar
    integer :: ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,Nf
    real    :: rho,rhou,rhov,rhow,rhotheta,rhocontrau,rhocontrav,pl,pr,pressure
    real    :: rhol,rhowl,rhothetal,rhocontraul,rhocontravl
    real    :: rhor,rhowr,rhothetar,rhocontraur,rhocontravr
    real    :: dzz,l,lambda,theta
    real,dimension(:,:,:),allocatable  :: tec_var
  ! -----------------------------------------
    call gp%SetIndex(ids,ide,jds,jde,kds,kde, &
                     ims,ime,jms,jme,kms,kme )
    call gp%SetCubeindex(Nf=Nf)

    do n=1,Nf
      do j=jds,jde
        do i=ids,ide

          do k=kds,kde
            rhocontrau=q(2)%pv(i,j,k,n)/jm%jab(i,j,k,n)
            rhocontrav=q(3)%pv(i,j,k,n)/jm%jab(i,j,k,n)
            call pprop2sp(lambda,theta,coord%x(i),coord%y(j),n)
            call contravprop2sp(rhou,rhov,rhocontrau,rhocontrav,n,lambda,theta)
            rhow=q(4)%pv(i,j,k,n)/jm%jab(i,j,k,n)

            rho=q(1)%pv(i,j,k  ,n)/jm%jab(i,j,k,n)
            rhotheta=q(5)%pv(i,j,k  ,n)/jm%jab(i,j,k,n)
            pressure=(Rd*rhotheta/p0)**gamma*p0

            q_3d(1)%pv3d(i,j,k,n) = rho
            q_3d(2)%pv3d(i,j,k,n) = rhou/rho
            q_3d(3)%pv3d(i,j,k,n) = rhov/rho
            q_3d(4)%pv3d(i,j,k,n) = rhow/rho
            q_3d(5)%pv3d(i,j,k,n) = pressure/(Rd*rho)
          enddo

        enddo
      enddo
    enddo

  end subroutine result2d_modelvar

  subroutine result2d_isobar(q_3d,q,jm,gm,plevs,nvar,nlev)

    implicit none

    integer,intent(in)    :: nvar,nlev
    type(jame),intent(in) :: jm
    type(gmap),intent(in) :: gm
    real,dimension(nlev),intent(in)        :: plevs
    type(cell),dimension(nvar),intent(in)  :: q
    type(pv3d),dimension(nvar),intent(out) :: q_3d
  ! Local variables
    integer :: i,j,k,n,kk,LMAX,IMAX
    integer :: ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,Nf
    real    :: rho,rhotheta,contrau,contrav,th,pix
    real    :: u,v,w
    real    :: lambda,theta
    real,dimension(:),allocatable   :: PLN,Z,ZI,PILN,Z2,Z3
    real    :: t_top,q_top,z_top,logpw_top,ts, qs,zs,logpw1,pre
  ! -----------------------------------------
    call gp%SetIndex(ids,ide,jds,jde,kds,kde, &
                     ims,ime,jms,jme,kms,kme )
    call gp%SetCubeindex(Nf=Nf)

    LMAX=kde-kds+1;IMAX=nlev
    allocate(PLN(LMAX),Z(LMAX),Z2(LMAX),Z3(LMAX),PILN(IMAX),ZI(IMAX))
    PLN=0.;Z=0.;ZI=0.

    do n=1,Nf
      do j=jds,jde
        do i=ids,ide

          PILN=plevs          
          do k=kds,kde
            rho = q(1)%pv(i,j,k,n)/jm%jab(i,j,k,n)
            rhotheta=q(5)%pv(i,j,k,n)/jm%jab(i,j,k,n)
            PLN(k)=p0*(Rd*(rhotheta)/p0)**gamma
            th  = rhotheta/rho
            pix = (PLN(k)/P0)**(Rd/CP)
            Z(k)=th*pix
          enddo

!          PILN=log(PILN/100.)
!          PLN=log(PLN/100.)
!          call SPLINE(ZI,PILN,IMAX,Z,PLN,LMAX,1,GRA,Rd)
          call SPLINE(ZI,log(PILN/100.),IMAX,Z,log(PLN/100.),LMAX,1,GRA,Rd)

          do k=1,nlev
            q_3d(5)%pv3d(i,j,k,n)=ZI(k)
            q_3d(1)%pv3d(i,j,k,n)=PILN(k)/ZI(k)/Rd
          enddo

          do k=kds,kde
            rho = q(1)%pv(i,j,k,n)/jm%jab(i,j,k,n)
            rhotheta=q(5)%pv(i,j,k,n)/jm%jab(i,j,k,n)
            contrau=q(2)%pv(i,j,k,n)/jm%jab(i,j,k,n)/rho
            contrav=q(3)%pv(i,j,k,n)/jm%jab(i,j,k,n)/rho
            w=q(4)%pv(i,j,k,n)/jm%jab(i,j,k,n)/rho
            PLN(k)=p0*(Rd*(rhotheta)/p0)**gamma

            call pprop2sp(lambda,theta,coord%x(i),coord%y(j),n)
            call contravprop2sp(u,v,contrau,contrav,n,lambda,theta)

            Z (k)=u
            Z2(k)=v
            Z3(k)=w
          enddo

          call SPLINE(ZI,log(PILN/100.),IMAX,Z,log(PLN/100.),LMAX,2,GRA,Rd)
          do k=1,nlev
            q_3d(2)%pv3d(i,j,k,n)=ZI(k)
          enddo
          call SPLINE(ZI,log(PILN/100.),IMAX,Z2,log(PLN/100.),LMAX,2,GRA,Rd)
          do k=1,nlev
            q_3d(3)%pv3d(i,j,k,n)=ZI(k)
          enddo
          call SPLINE(ZI,log(PILN/100.),IMAX,Z3,log(PLN/100.),LMAX,2,GRA,Rd)
          do k=1,nlev
            q_3d(4)%pv3d(i,j,k,n)=ZI(k)
          enddo

        enddo
      enddo
    enddo

    deallocate(PLN,Z,ZI,PILN,Z2,Z3)

  end subroutine result2d_isobar

  subroutine geopotentialheight(q_3d,q,jm,gm,plevs,nvar,nlev)

    implicit none

    integer,intent(in)    :: nvar,nlev
    type(jame),intent(in) :: jm
    type(gmap),intent(in) :: gm
    real,dimension(nlev),intent(in)          :: plevs
    type(cell),dimension(nvar),intent(in)  :: q
    type(pv3d),intent(out) :: q_3d
  ! Local variables
    integer    :: i,j,k,n,kk,LMAX,IMAX
    integer    :: ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,Nf
    real       :: rhotheta,pl,pr,l,dzz,rho,pre
    real       :: t_top,q_top,z_top,logpw_top,ts, qs,zs,logpw1
    real,dimension(:),allocatable   :: PLN,Z,ZI,PILN
  ! -----------------------------
    call gp%SetIndex(ids,ide,jds,jde,kds,kde, &
                     ims,ime,jms,jme,kms,kme )
    call gp%SetCubeindex(Nf=Nf)

    LMAX=kde-kds+1;IMAX=nlev
    allocate(PLN(LMAX),Z(LMAX),PILN(IMAX),ZI(IMAX))
    PLN=0.;Z=0.;PILN=plevs;Z=0.
    do n=1,Nf
      do j=jds,jde
        do i=ids,ide
          
          do k=kds,kde
            rhotheta=q(5)%pv(i,j,k,n)/jm%jab(i,j,k,n)
            PLN(k)=p0*(Rd*(rhotheta)/p0)**gamma
            Z(k)=gm%zz(i,j,k,n)
          enddo

          call SPLINE(ZI,log(PILN/100.),IMAX,Z,log(PLN/100.),LMAX,1,GRA,Rd)

          !PILN=log(PILN/100.)
          !PLN=log(PLN/100.)
          !call SPLINE(ZI,PILN,IMAX,Z,PLN,LMAX,1,GRA,Rd)

          !k=kde
          !rho=q(1)%pv(i,j,k,n)/jm%jab(i,j,k,n)
          !rhotheta=q(5)%pv(i,j,k,n)/jm%jab(i,j,k,n)
          !pre=p0*(Rd*(rhotheta)/p0)**gamma
          !t_top=pre/rho/Rd
          !q_top=0.;z_top=gm%zz(i,j,k,n)
          !logpw_top=log(pre)

          !k=kds
          !rho=q(1)%pv(i,j,k,n)/jm%jab(i,j,k,n)
          !rhotheta=q(5)%pv(i,j,k,n)/jm%jab(i,j,k,n)
          !pre=p0*(Rd*(rhotheta)/p0)**gamma
          !ts=pre/rho/Rd
          !qs=0.;zs=gm%zz(i,j,k,n)
          !logpw1=log(pre)

          !call SPLINE(ZI,PILN,LMAX,Z,PLN,IMAX,                  &
          !            t_top, q_top, z_top, logpw_top,           &
          !            ts, qs, zs, logpw1, 1)

          do k=1,nlev
            q_3d%pv3d(i,j,k,n)=ZI(k)
          enddo

        enddo
      enddo
    enddo

    deallocate(PLN,Z,PILN,ZI)

  end subroutine geopotentialheight

  subroutine surfacepresure(q_2d,q,jm,nvar)

    implicit none

    integer,intent(in)    :: nvar
    type(jame),intent(in) :: jm
    type(cell),dimension(nvar),intent(in)  :: q
    type(pv2d),intent(out) :: q_2d
  ! Local variables
    integer    :: i,j,k,n,kk
    integer    :: ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,Nf
    real       :: rhotheta,pp
  ! -----------------------------
    call gp%SetIndex(ids,ide,jds,jde,kds,kde, &
                     ims,ime,jms,jme,kms,kme )
    call gp%SetCubeindex(Nf=Nf)

    do n=1,Nf
      do j=jds,jde
        do i=ids,ide
          rhotheta = q(5)%pv(i,j,1,n)/jm%jab(i,j,1,n)
          pp=p0*(Rd*(rhotheta)/p0)**(gamma)
          q_2d%pv2d(i,j,n) = pp/100.    ! convert it into hPa
        enddo
      enddo
    enddo

  end subroutine surfacepresure

  subroutine diagnose_ps(psfc,tp,zp,ALNP,gm,nzp)

    implicit none

    integer,intent(in)    :: nzp
    type(gmap),intent(in) :: gm
    real,dimension(nzp),intent(in)    :: ALNP   ! Ln(pressure)
    type(pv3d),intent(in)  :: zp, &   ! the geopotentional height at the isobar.
                              tp      ! the temperature.
    type(pv2d),intent(out) :: psfc
  ! Local variables
    integer    :: i,j,k,n,kp
    integer    :: ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,Nf
    real       :: dh,dphi,a,b,phi3,sq,alnps
  ! ---------------------------------
    call gp%SetIndex(ids,ide,jds,jde,kds,kde, &
                     ims,ime,jms,jme,kms,kme )
    call gp%SetCubeindex(Nf=Nf)

    psfc%pv2d=0.
    do n=1,Nf
      do j=jds,jde
        do i=ids,ide

          DO K=2,nzp
            KP=K
            IF(zp%pv3d(i,j,k,n).GE.gm%zs(i,j,n)) exit
          ENDDO 

          DH=ALNP(KP-1)-ALNP(KP)
          DPHI=GRA* (ZP%pv3d(i,j,KP-1,n)-ZP%pv3d(i,j,KP,n))
          A=-DPHI/DH
          B=Rd*(TP%pv3d(i,j,KP-1,n)-TP%pv3d(i,j,KP,n))/DH
          PHI3=0.5*GRA*(ZP%pv3d(i,j,KP-1,n)+ZP%pv3d(i,j,KP,n))+0.125*B*(DH**2)
          SQ=A**2-2.*B*(GRA*gm%zs(i,j,n) - PHI3)

          IF(SQ.GE.0.) THEN 
            ALNPS=0.5*(ALNP(KP-1)+ALNP(KP))-2.*(GRA*gm%zs(i,j,n)-PHI3)/(A+SQRT(SQ))
          ELSE
            SQ=A**2
          ENDIF

          psfc%pv2d(i,j,n) = exp(ALNPS)

        enddo
      enddo
    enddo

  end subroutine diagnose_ps

  real function interplotant_1d_post(q,xxx,xxxx)

    USE global, only: NP,Ra

    implicit none

    real,intent(in)  :: q(NP),xxx,xxxx(NP)
  ! Local varialbes
    real :: xx,x(NP),a
  ! -------------------------------
    a=Ra

    !xx=xxx/a; x=xxxx/a
    xx=xxx; x=xxxx

    interplotant_1d_post=&
                    q(1)*(xx-x(2))/(x(1)-x(2))*(xx-x(3))/(x(1)-x(3))+&
                    q(2)*(xx-x(1))/(x(2)-x(1))*(xx-x(3))/(x(2)-x(3))+&
                    q(3)*(xx-x(1))/(x(3)-x(1))*(xx-x(2))/(x(3)-x(2))
  
    return
  end function interplotant_1d_post

  real function linear_1d(q0,q1,x,x0,x1)

    implicit none 

    real,intent(in)  :: x0,x1,x
    real,intent(in)  :: q0,q1
  ! local variables

  ! ---------------------------
    linear_1d=(x-x1)/(x0-x1)*q0+(x-x0)/(x1-x0)*q1

  end function linear_1d

! --------------------------------------------------------------------------
! INTERPOLATION USING A CUBIC NONPERIODIC SPLINE
! ZI(I) : INTERPOLATION VALUE AT PILN(I)
! PILN(I) : LN(P) COORDINATE IN DECREASING ORDER
! LMAX : NUMBER OF INTERPOLATION POINT
! Z(I) : DATA AT PLN(I)
! PLN(I) : LN(P) COORDINATE IN DECREASING ORDER
! IMAX : NUMBER OF DATA POINT
! SM(I) : SECOND DERIVATIVEA AT DATA POINT
! ******* 1 INTERPOLATION FOR HEIGHT
! * IDX * 2 INTERPOLATION FOR WIND
! ******* 3 INTERPOLATION FOR TEMPERATURE
! CUBIC SPLINE FOR HEIGHT IS DIFFERENCIATED WITH RESPECT TO LN(P)
! ----------- HYDROSTATIC RELATION
!
  SUBROUTINE SPLINE(ZI,PILN,LMAX,Z,PLN,IMAX,IDX,grav,r_d)

    IMPLICIT NONE

    INTEGER , INTENT(IN ) :: LMAX, IMAX,IDX
    REAL, INTENT(IN ) :: grav,r_d
    REAL, DIMENSION(LMAX),INTENT(IN ) :: PILN
    REAL, DIMENSION(IMAX),INTENT(IN ) :: Z,PLN
    REAL, DIMENSION(LMAX),INTENT(OUT) :: ZI
  ! Local variables
    INTEGER :: i, k, l, im1,IB, klevel
    real :: GR, x
    real,dimension(200) :: SM,H,AL,AM,AP,C
  ! ------------------------------------------
    x=0.
    klevel = 0

    IM1=IMAX-1

    GR=-grav/r_d

    DO I=2,IMAX
       H(I)=PLN(I)-PLN(I-1)
    ENDDO

    DO I=2,IM1
       AL(I)=0.5*H(I+1)/(H(I)+H(I+1))
       AM(I)=0.5-AL(I)
    ENDDO

! END CONDITION FOR HEIGHT AND TEMPERATURE (LAPSE RATE IS CONSTANT)
! SM(1)=SM(2) ; SM(IMAX-1)=SM(IMAX)

    select case ( idx )
      case (1,3)
       AL(1)=-1.
       AM(IMAX)=-1.
! END CONDITION FOR WIND (WIND SHEAR IS CONSTANT)
! SM(1)=0. ; SM(IMAX)=0.
     case(2)
       AL(1)=0.
       AM(IMAX)=0.
    end select

    AL(IMAX)=0.
    DO I=2,IMAX
       AP(I)=1./(1.-AL(I-1)*AM(I))
       AL(I)=AL(I)*AP(I)
    ENDDO
!
    C(1)=0.
    C(IMAX)=0.
    DO I=2,IM1
      C(I)=3.*((Z(I+1)-Z(I))/H(I+1)-(Z(I)-Z(I-1))/H(I))/(H(I)+H(I+1))
    ENDDO

! FORWARD SUBSTITUTION
    DO I=2,IMAX
       C(I)=(C(I)-C(I-1)*AM(I))*AP(I)
    ENDDO
    SM(IMAX)=C(IMAX)

! BACKWARD SUBSTITUTION
    DO K=1,IM1
       I=IMAX-K
       SM(I)=C(I)-AL(I)*SM(I+1)
    ENDDO

! INTERPOLATION
    IB=2
    DO L=1,LMAX
       X=PILN(L)
       DO I=IB,IMAX
          IF( pln(i-1) >= x .and. x >= pln(i)) then
             klevel = i
          ELSEIF( x >= pln(1)) then
             klevel = 2
          ELSEIF( x <= pln(imax)) then
             klevel = imax
          ENDIF
       ENDDO

       select case ( idx )
         case (1,2)
           ZI(L)=(PLN(klevel)-X)/H(klevel)*(Z(klevel-1)-SM(klevel-1)/6.*(X-PLN(klevel-1)) &
               *(H(klevel)+PLN(klevel)-X))+(X-PLN(klevel-1))/H(klevel)*(Z(klevel)-SM(klevel)/ &
               6.*(PLN(klevel)-X)*(H(klevel)+X-PLN(klevel-1)))
! DIFFERENTIAL CALCULUS OF CUBIC SPLINE FOR HEIGHT
! GR=-G/R : COEFFICIENT IN HYDROSTATIC EQUATION
         case (3)
           ZI(L)=SM(klevel-1)*(-(PLN(klevel)-X)**2/(2.*H(klevel))+H(klevel)/6.) &
                +SM(klevel)*((X-PLN(klevel-1))**2/(2.*H(klevel))-H(klevel)/6.) &
                +(Z(klevel)-Z(klevel-1))/H(klevel)
           ZI(L)=ZI(L)*GR
       end select
    ENDDO

    RETURN
  END SUBROUTINE SPLINE

END MODULE postprocess
