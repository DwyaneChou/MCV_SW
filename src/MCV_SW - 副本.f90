    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Multi-moment Constrained finite Volume shallow water dynamic core !
    ! Written by Lilong Zhou, March 26, 2019                            !
    ! Numerical Weather Prediction Center of CMA                        !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    program MCV_SW
      use parameters_mod
      use stat_mod
      use tend_mod
      use mesh_mod
      use test_case_mod
      use time_scheme_mod
      use output_mod
      use spatial_operators_mod
      implicit none
      
      integer it
      integer :: old = 0
      integer :: new = 1
      
      integer :: timeStart,timeEnd
      
      !Timing start
      call SYSTEM_CLOCK(timeStart)
      
      call time_init
      call initParameters
      call initMesh
      call initStat
      call initTend
      call initTestCase
      
      call history_init
      
      ! time integration
      do it = 1,nsteps
        call RK4(stat(new),stat(old))
        
        if(mod(it*dt,float(history_interval))==0.and.(it*dt>=history_interval))then
          print*,it,'/',nsteps
          call history_write_stat(stat(new),it*dt)
        endif
        
        call switch_stat(stat(old), stat(new))
      enddo
      
      ! Timing end
      call SYSTEM_CLOCK(timeEnd)
      
      print*,'It took ',dble(timeEnd-timeStart)/10000.0,' seconds to run this program'
    end program MCV_SW
