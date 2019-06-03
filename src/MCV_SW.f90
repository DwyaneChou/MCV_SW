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
      
      integer :: output_idx, total_output_num
      
      integer :: timeStart,timeEnd
      
      !Timing start
      call SYSTEM_CLOCK(timeStart)
      
      call initParameters
      call initMesh
      call initStat
      call initTend
      call initTestCase
      
      call history_init(stat(old))
      call history_write_stat(stat(old),1)
      
      ! time integration
      output_idx       = 1
      total_output_num = nsteps * dt / history_interval + 1
      do it = 1,nsteps
        call RK4(stat(new),stat(old))
        
        if(mod(it*dt,float(history_interval))==0.and.(it*dt>=history_interval))then
          output_idx = output_idx + 1
          call history_write_stat(stat(new),output_idx)
          print*,'output index/total : ',output_idx,'/',total_output_num
        endif
        
        call switch_stat(stat(old), stat(new))
      enddo
      
      ! Timing end
      call SYSTEM_CLOCK(timeEnd)
      
      print*,'It took ',dble(timeEnd-timeStart)/1000.0,' seconds to run this program'
    end program MCV_SW
