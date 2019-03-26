MODULE diag

  USE global     ,only: np
  USE domain     ,only: jame,cord,cell,gp
  !USE projection
  !USE pprocess

  implicit none

contains

  real function totalmass(q,jm,coord,nvar)

    implicit none 

    integer,intent(in)    :: nvar
    type(jame),intent(in) :: jm
    type(cord),intent(in) :: coord
    type(cell),dimension(nvar),intent(in)  :: q
  ! Local variables
    integer   :: n,i,j,k,ii,jj,kk,inv,ip,jp,kp
    integer   :: ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,Nf,Nx,Ny,Nz
    real      :: mass,m0,dx,dy,dz,masslayer1,masslayer2,masslayer3
    real,dimension(Np,NP,NP)   :: qm
    real,dimension(NP,NP)      :: mkj
    real,dimension(NP)         :: mk
  !  --------------------

    call gp%SetIndex(ids,ide,jds,jde,kds,kde, &
                     ims,ime,jms,jme,kms,kme )
    call gp%SetCubeindex(Nx=Nx,Ny=Ny,Nz=Nz,Nf=Nf)

    !call output_patch(q,jm,coord,5,1,0)

    mass=0.
    inv=2
    do n=1,Nf
      do k=1,Nz
        do j=1,Ny
          do i=1,Nx
            
            do kp=1,NP
              kk= (k-1)*inv+1+(kp-1)
              do jp=1,NP
                jj = (j-1)*inv+1+(jp-1)
                do ip=1,NP
                  ii = (i-1)*inv+1+(ip-1)
                  qm(ip,jp,kp)=q(1)%pv(ii,jj,kk,n)
                enddo
              enddo
            enddo

            do kp=1,NP
              do jp=1,NP
                mkj(jp,kp)=1./6.*qm(1,jp,kp)+2./3.*qm(2,jp,kp)+1./6.*qm(3,jp,kp)
              enddo
            enddo
            
            do kp=1,NP
              mk(kp)=1./6.*mkj(1,kp)+2./3.*mkj(2,kp)+1./6.*mkj(3,kp)
            enddo
            
            ii=2*i;jj=2*j;kk=2*k
            dx=coord%x(ii+1)-coord%x(ii-1)
            dy=coord%y(jj+1)-coord%y(jj-1)
            dz=coord%z(kk+1)-coord%z(kk-1)
            mass=mass+(1./6.*mk(1)+2./3.*mk(2)+1./6.*mk(3))*dx*dy*dz

          enddo
        enddo
      enddo
    enddo

    totalmass=mass

  end function totalmass

  real function l2error(q,q0,jm,coord,nvar)

    implicit none 

    integer,intent(in)    :: nvar
    type(jame),intent(in) :: jm
    type(cord),intent(in) :: coord
    type(cell),dimension(nvar),intent(in)  :: q,q0
  ! Local variables
    integer   :: n,i,j,k,ii,jj,kk,inv,ip,jp,kp
    integer   :: ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,Nf,Nx,Ny,Nz
    real      :: error2,mass2,dx,dy,dz
    real,dimension(NP,NP,NP)   :: qm
    real,dimension(NP,NP)      :: mkj
    real,dimension(NP)         :: mk
  ! ----------------------
    call gp%SetIndex(ids,ide,jds,jde,kds,kde, &
                     ims,ime,jms,jme,kms,kme )
    call gp%SetCubeindex(Nx=Nx,Ny=Ny,Nz=Nz,Nf=Nf)

    inv=2
    do n=1,Nf
      do k=1,Nz
        do j=1,Ny
          do i=1,Nx
            ! q
            do kp=1,NP
              kk= (k-1)*inv+1+(kp-1)
              do jp=1,NP
                jj = (j-1)*inv+1+(jp-1)
                do ip=1,NP
                  ii = (i-1)*inv+1+(ip-1)
                  qm(ip,jp,kp)=q(1)%pv(ii,jj,kk,n)
                enddo
              enddo
            enddo

            do kp=1,NP
              do jp=1,NP
                mkj(jp,kp)=1./6.*qm(1,jp,kp)+2./3.*qm(2,jp,kp)+1./6.*qm(3,jp,kp)
              enddo
            enddo
            
            do kp=1,NP
              mk(kp)=1./6.*mkj(1,kp)+2./3.*mkj(2,kp)+1./6.*mkj(3,kp)
            enddo

            q(1)%via(i,j,k,n)=1./6.*mk(1)+2./3.*mk(2)+1./6.*mk(3)

            ! q0
            do kp=1,NP
              kk= (k-1)*inv+1+(kp-1)
              do jp=1,NP
                jj = (j-1)*inv+1+(jp-1)
                do ip=1,NP
                  ii = (i-1)*inv+1+(ip-1)
                  qm(ip,jp,kp)=q0(1)%pv(ii,jj,kk,n)
                enddo
              enddo
            enddo

            do kp=1,NP
              do jp=1,NP
                mkj(jp,kp)=1./6.*qm(1,jp,kp)+2./3.*qm(2,jp,kp)+1./6.*qm(3,jp,kp)
              enddo
            enddo
            
            do kp=1,NP
              mk(kp)=1./6.*mkj(1,kp)+2./3.*mkj(2,kp)+1./6.*mkj(3,kp)
            enddo

            q0(1)%via(i,j,k,n)=1./6.*mk(1)+2./3.*mk(2)+1./6.*mk(3)

          enddo
        enddo
      enddo
    enddo

    error2=0.
    mass2 =0.
    do n=1,Nf
      do k=1,Nz
        do j=1,Ny
          do i=1,Nx

            ii=2*i;jj=2*j;kk=2*k
            dx=coord%x(ii+1)-coord%x(ii-1)
            dy=coord%y(jj+1)-coord%y(jj-1)
            dz=coord%z(kk+1)-coord%z(kk-1)

            error2=error2+(q(1)%via(i,j,k,n)-q0(1)%via(i,j,k,n))**2*dx*dy*dz
            mass2 =mass2+q(1)%via(i,j,k,n)**2*dx*dy*dz

          enddo
        enddo
      enddo
    enddo

    l2error=error2/mass2

  end function l2error

END MODULE diag

