  module time_scheme_mod
    use stat_mod
    use tend_mod
    use parameters_mod
    use spatial_operators_mod
    use ghost_mod, only : CubedSphereFillHalo_Linear_extended
    implicit none
    
    contains
    
    subroutine RK4(stat_new,stat_old)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(in   ) :: stat_old

      integer, parameter :: old = 0
      integer, parameter :: new = 1
      
      integer, parameter :: one   = -1
      integer, parameter :: two   = -2
      integer, parameter :: three = -3
      integer, parameter :: four  = -4
      
      integer iPatch
      
      stat(one)%phi = stat_old%phi
      stat(one)%u   = stat_old%u
      stat(one)%v   = stat_old%v
      
      call spatial_operator(stat(one), tend(one))
      call update_stat(stat(two), stat(one), tend(one), 0.5 * dt)
      print*,ips,ipe
      do iPatch = ifs, ife
        call CubedSphereFillHalo_Linear_extended(stat(two)%u  (ids:ide,jds:jde,iPatch), stat(two)%u  (ips:ipe,jps:jpe,iPatch), iPatch, ide+1, xhalo*(DOF-1))
        call CubedSphereFillHalo_Linear_extended(stat(two)%v  (ids:ide,jds:jde,iPatch), stat(two)%v  (ips:ipe,jps:jpe,iPatch), iPatch, ide+1, xhalo*(DOF-1))
        call CubedSphereFillHalo_Linear_extended(stat(two)%phi(ids:ide,jds:jde,iPatch), stat(two)%phi(ips:ipe,jps:jpe,iPatch), iPatch, ide+1, xhalo*(DOF-1))
      enddo
      
      call spatial_operator(stat(two), tend(two))
      call update_stat(stat(three), stat(one), tend(two), 0.5 * dt)
      do iPatch = ifs, ife
        call CubedSphereFillHalo_Linear_extended(stat(three)%u  (ids:ide,jds:jde,iPatch), stat(three)%u  (ips:ipe,jps:jpe,iPatch), iPatch, ide+1, xhalo*(DOF-1))
        call CubedSphereFillHalo_Linear_extended(stat(three)%v  (ids:ide,jds:jde,iPatch), stat(three)%v  (ips:ipe,jps:jpe,iPatch), iPatch, ide+1, xhalo*(DOF-1))
        call CubedSphereFillHalo_Linear_extended(stat(three)%phi(ids:ide,jds:jde,iPatch), stat(three)%phi(ips:ipe,jps:jpe,iPatch), iPatch, ide+1, xhalo*(DOF-1))
      enddo
      
      call spatial_operator(stat(three), tend(three))
      call update_stat(stat(four), stat(one), tend(three), dt)
      do iPatch = ifs, ife
        call CubedSphereFillHalo_Linear_extended(stat(four)%u  (ids:ide,jds:jde,iPatch), stat(four)%u  (ips:ipe,jps:jpe,iPatch), iPatch, ide+1, xhalo*(DOF-1))
        call CubedSphereFillHalo_Linear_extended(stat(four)%v  (ids:ide,jds:jde,iPatch), stat(four)%v  (ips:ipe,jps:jpe,iPatch), iPatch, ide+1, xhalo*(DOF-1))
        call CubedSphereFillHalo_Linear_extended(stat(four)%phi(ids:ide,jds:jde,iPatch), stat(four)%phi(ips:ipe,jps:jpe,iPatch), iPatch, ide+1, xhalo*(DOF-1))
      enddo
      
      call spatial_operator(stat(four), tend(four))
      
      tend(new)%phi = (tend(one)%phi + 2.0 * tend(two)%phi + 2.0 * tend(three)%phi + tend(four)%phi) / 6.
      tend(new)%u   = (tend(one)%u   + 2.0 * tend(two)%u   + 2.0 * tend(three)%u   + tend(four)%u  ) / 6.
      tend(new)%v   = (tend(one)%v   + 2.0 * tend(two)%v   + 2.0 * tend(three)%v   + tend(four)%v  ) / 6.
      
      call update_stat(stat_new, stat_old, tend(new), dt)

    end subroutine RK4
    
    subroutine update_stat(stat_new, stat_old, tend, dt)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(in   ) :: stat_old
      type(tend_field), intent(in   ) :: tend
      
      real   , intent(in ) :: dt
      
      integer iPatch
      
      stat_new%phi(ids:ide,jds:jde,ifs:ife) = stat_old%phi(ids:ide,jds:jde,ifs:ife) + dt * tend%phi(ids:ide,jds:jde,ifs:ife)
      stat_new%u  (ids:ide,jds:jde,ifs:ife) = stat_old%u  (ids:ide,jds:jde,ifs:ife) + dt * tend%u  (ids:ide,jds:jde,ifs:ife)
      stat_new%v  (ids:ide,jds:jde,ifs:ife) = stat_old%v  (ids:ide,jds:jde,ifs:ife) + dt * tend%v  (ids:ide,jds:jde,ifs:ife)
      
      do iPatch = ifs, ife
        call CubedSphereFillHalo_Linear_extended(stat_new%u  (ids:ide,jds:jde,iPatch), stat_new%u  (ips:ipe,jps:jpe,iPatch), iPatch, ide+1, xhalo*(DOF-1))
        call CubedSphereFillHalo_Linear_extended(stat_new%v  (ids:ide,jds:jde,iPatch), stat_new%v  (ips:ipe,jps:jpe,iPatch), iPatch, ide+1, xhalo*(DOF-1))
        call CubedSphereFillHalo_Linear_extended(stat_new%phi(ids:ide,jds:jde,iPatch), stat_new%phi(ips:ipe,jps:jpe,iPatch), iPatch, ide+1, xhalo*(DOF-1))
      enddo
      
    end subroutine update_stat
    
    subroutine switch_stat(stat_old,stat_new)
      type(stat_field), intent(inout) :: stat_old
      type(stat_field), intent(in   ) :: stat_new
      
      stat_old%phi = stat_new%phi
      stat_old%u   = stat_new%u  
      stat_old%v   = stat_new%v
    end subroutine switch_stat
    
  end module time_scheme_mod