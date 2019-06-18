  module time_scheme_mod
    use stat_mod
    use tend_mod
    use mesh_mod
    use parameters_mod
    use spatial_operators_mod
    use ghost_mod
    use output_mod
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
      
      call copyStat(stat(one),stat_old)
      
      call spatial_operator (stat(one), tend(one))
      call update_stat      (stat(two), stat(one), tend(one), 0.5 * dt)
      
      call spatial_operator (stat(two), tend(two))
      call update_stat      (stat(three), stat(one), tend(two), 0.5 * dt)
      
      call spatial_operator (stat(three), tend(three))
      call update_stat      (stat(four), stat(one), tend(three), dt)
      
      call spatial_operator(stat(four), tend(four))
      
      tend(new)%phiG = (tend(one)%phiG + 2. * tend(two)%phiG + 2. * tend(three)%phiG + tend(four)%phiG) / 6.
      
      call update_stat      (stat_new, stat_old, tend(new), dt)

    end subroutine RK4
    
    subroutine RK3_TVD(stat_new,stat_old)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(in   ) :: stat_old

      integer, parameter :: old = 0
      integer, parameter :: new = 1
      
      integer, parameter :: one   = -1
      integer, parameter :: two   = -2
      integer, parameter :: three = -3
      integer, parameter :: four  = -4
      
      call copyStat(stat(one),stat_old)
      
      call spatial_operator (stat(one), tend(one))
      call update_stat      (stat(two), stat(one), tend(one), dt)
      
      call spatial_operator (stat(two), tend(two))
      call update_stat_RK3_TVD_1(stat(three), stat(one), stat(two), tend(two))
      
      call spatial_operator (stat(three), tend(three))
      call update_stat_RK3_TVD_2(stat_new, stat(one), stat(three), tend(three))

    end subroutine RK3_TVD
    
    subroutine update_stat(stat_new, stat_old, tend, inc_t)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(in   ) :: stat_old
      type(tend_field), intent(in   ) :: tend
      real            , intent(in   ) :: inc_t
      
      stat_new%phiG(ids:ide,jds:jde,ifs:ife) = stat_old%phiG(ids:ide,jds:jde,ifs:ife) + inc_t * tend%phiG(ids:ide,jds:jde,ifs:ife)
      
      stat_new%phi = stat_new%phiG / mesh%sqrtG
      
      call correct_bdy_ghost(stat_new)
      
      stat_new%phiG = stat_new%phi * mesh%sqrtG
      
      stat_new%u               = stat_old%u
      stat_new%v               = stat_old%v
      stat_new%contraU         = stat_old%contraU
      stat_new%contraV         = stat_old%contraV
      stat_new%zonal_wind      = stat_old%zonal_wind
      stat_new%meridional_wind = stat_old%meridional_wind
    end subroutine update_stat
    
    subroutine update_stat_RK3_TVD_1(stat_new, stat_old,stat1, tend)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(in   ) :: stat_old
      type(stat_field), intent(in   ) :: stat1
      type(tend_field), intent(in   ) :: tend
      
      stat_new%phiG(ids:ide,jds:jde,ifs:ife) = 0.75 * stat_old%phiG(ids:ide,jds:jde,ifs:ife) + 0.25 * stat1%phiG(ids:ide,jds:jde,ifs:ife) + 0.25 * dt * tend%phiG(ids:ide,jds:jde,ifs:ife)
      
      stat_new%phi = stat_new%phiG / mesh%sqrtG
      
      call correct_bdy_ghost(stat_new)
      
      stat_new%phiG = stat_new%phi * mesh%sqrtG
      
      stat_new%u               = stat_old%u
      stat_new%v               = stat_old%v
      stat_new%contraU         = stat_old%contraU
      stat_new%contraV         = stat_old%contraV
      stat_new%zonal_wind      = stat_old%zonal_wind
      stat_new%meridional_wind = stat_old%meridional_wind
    end subroutine update_stat_RK3_TVD_1
    
    subroutine update_stat_RK3_TVD_2(stat_new, stat_old,stat2, tend)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(in   ) :: stat_old
      type(stat_field), intent(in   ) :: stat2
      type(tend_field), intent(in   ) :: tend
      
      stat_new%phiG(ids:ide,jds:jde,ifs:ife) = stat_old%phiG(ids:ide,jds:jde,ifs:ife) / 3. + 2./3. * stat2%phiG(ids:ide,jds:jde,ifs:ife) + 2./3. * dt * tend%phiG(ids:ide,jds:jde,ifs:ife)
      
      stat_new%phi = stat_new%phiG / mesh%sqrtG
      
      call correct_bdy_ghost(stat_new)
      
      stat_new%phiG = stat_new%phi * mesh%sqrtG
      
      stat_new%u               = stat_old%u
      stat_new%v               = stat_old%v
      stat_new%contraU         = stat_old%contraU
      stat_new%contraV         = stat_old%contraV
      stat_new%zonal_wind      = stat_old%zonal_wind
      stat_new%meridional_wind = stat_old%meridional_wind
    end subroutine update_stat_RK3_TVD_2
    
    subroutine correct_bdy_ghost(stat)
      type(stat_field), intent(inout) :: stat
      
      ! Fill ghost band and unify values of common points on boundary
      call fill_ghost              (stat)
      call unify_bdy_stat          (stat)
      
    end subroutine correct_bdy_ghost
    
  end module time_scheme_mod