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
      use ghost_mod
      use test_case_mod
      use time_scheme_mod
      use output_mod
      use spatial_operators_mod
      use diag_mod
      implicit none
      
      integer it
      
      integer :: old = 0
      integer :: new = 1
      
      real    :: total_mass0,total_energy0
      real    :: total_mass ,total_energy
      
      integer :: output_idx, total_output_num
      
      integer :: timeStart,timeEnd
      
      !Timing start
      call SYSTEM_CLOCK(timeStart)
      
      print*,'Start the MCV Shallow Water model'
      call initParameters
      call initMesh
      call initStat
      call initTend
      call initGhost
      call initTestCase
      
      print*,'The MCV order is               ', DOF
      print*,''
      print*,'Temporal integration scheme is '//trim(adjustl(integral_scheme))
      print*,''
      
      ! time integration
      output_idx       = 0
      total_output_num = nsteps * dt / history_interval
      call history_init      (stat(old))
      call history_write_stat(stat(old),1)
      call calc_total_mass   (total_mass0  ,stat(old))
      call calc_total_energy (total_energy0,stat(old))
      print*,'output index/total, MCR, ECR :',output_idx,'/',total_output_num,' ',0., 0.
      
      do it = 1,nsteps
        if(trim(adjustl(integral_scheme))=='RK3_TVD')call RK3_TVD(stat(new),stat(old))
        if(trim(adjustl(integral_scheme))=='RK4'    )call RK4    (stat(new),stat(old))
        
        if(mod(it*dt,float(history_interval))==0.and.(it*dt>=history_interval))then
          output_idx = output_idx + 1
          call addFillValue(stat(new))
          call history_write_stat(stat(new),output_idx)
          call calc_total_mass  (total_mass  ,stat(new))
          call calc_total_energy(total_energy,stat(old))
          print*,'output index/total, MCR, ECR :',output_idx,'/',total_output_num,' ',(total_mass-total_mass0)/total_mass0, (total_energy-total_energy0)/total_energy0
        endif
        
        call switch_stat(stat(old), stat(new))
      enddo
      
      ! Timing end
      call SYSTEM_CLOCK(timeEnd)
      
      print*,'It took ',dble(timeEnd-timeStart)/1000.0,' seconds to run this program'
    end program MCV_SW
