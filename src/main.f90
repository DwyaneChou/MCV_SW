PROGRAM main

  USE omp_lib       , only: omp_set_num_threads,OMP_GET_WTIME
  USE global        , only: numvar
  USE domain        , only: solver_domain,copy,gp,gl,gm,coord,jm,q,q0,qh,qmoist,qn,rqsb
  USE init          , only: solver_init
  USE debug         , only: resultmax,resultmin
  USE timemarching  , only: rk3_ex,rk2_im_v
  USE postprocess   , only: gradsmodelvar,grads

  implicit none

 ! Local variables
  integer   :: step,Nf,fileid,fileid1
  integer   :: ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme
  real      :: tend,tstart,time
  character(len=30)   :: filename
 !-------------------------------
  
  call omp_set_num_threads(6)

  call solver_domain
  call solver_init

  call copy(q0,q,5)

  fileid1=101
  open(fileid1,file='./gradsmodelvar.dat',form='unformatted',access='sequential')
  call gradsmodelvar(q,qh,qmoist,gl,jm,5,5,fileid1)

  fileid=100
  filename='./grads.dat'
  open(fileid,file=trim(filename),form='unformatted',access='sequential')
  call grads(q,jm,gm,coord,5,filename,fileid)

  call resultmax(q,gl,jm,coord,5)
  call resultmin(q,gl,jm,coord,5)

  call gp%SetIndex(ids,ide,jds,jde,kds,kde, &
                   ims,ime,jms,jme,kms,kme )
  call gp%SetCubeindex(Nf=Nf)

  ! prognostic variables - PV Moments
  do step=1,gp%max_step

    print*,'Time step =', step ;tstart=OMP_GET_WTIME()

    time=(step-1)*gp%dt

    call RK3_ex(q,qn,qh,rqsb,gl,jm,coord,gp%dt/2.,1,NumVar); call copy(q,qn,5)
    call RK2_im_v(q,qn,qh,gl,jm,coord,gp%dt,NumVar); call copy(q,qn,5)
    call RK3_ex(q,qn,qh,rqsb,gl,jm,coord,gp%dt/2.,1,NumVar); call copy(q,qn,5)


    tend=OMP_GET_WTIME()
    write(6,*)'CPU time of step',step,':',tend-tstart
    
    if(mod(step,192)==0)then
      call grads(q,jm,gm,coord,5,filename,fileid)
      call gradsmodelvar(q,qh,qmoist,gl,jm,5,5,fileid1)
      call resultmax(q,gl,jm,coord,5)
      call resultmin(q,gl,jm,coord,5)
    endif

  enddo

  close(fileid)
  stop 'main stop!'

END PROGRAM main

