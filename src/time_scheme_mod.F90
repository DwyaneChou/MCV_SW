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
      
      stat(one)%phi             = stat_old%phi
      stat(one)%phiG            = stat_old%phiG
      stat(one)%u               = stat_old%u
      stat(one)%v               = stat_old%v
      stat(one)%contraU         = stat_old%contraU
      stat(one)%contraV         = stat_old%contraV
      stat(one)%zonal_wind      = stat_old%zonal_wind
      stat(one)%meridional_wind = stat_old%meridional_wind
      
      call spatial_operator (stat(one), tend(one))
      call update_stat      (stat(two), stat(one), tend(one), 0.5 * dt)
      call convert_wind_P2SP(stat(two))
      call fill_ghost       (stat(two))
      call unify_bdy_stat   (stat(two))
      call convert_wind_SP2P(stat(two))
      call convert_wind_cov2contrav(stat(two))
      
      call spatial_operator (stat(two), tend(two))
      call update_stat      (stat(three), stat(one), tend(two), 0.5 * dt)
      call convert_wind_P2SP(stat(three))
      call fill_ghost       (stat(three))
      call unify_bdy_stat   (stat(three))
      call convert_wind_SP2P(stat(three))
      call convert_wind_cov2contrav(stat(three))
      
      call spatial_operator (stat(three), tend(three))
      call update_stat      (stat(four), stat(one), tend(three), dt)
      call convert_wind_P2SP(stat(four))
      call fill_ghost       (stat(four))
      call unify_bdy_stat   (stat(four))
      call convert_wind_SP2P(stat(four))
      call convert_wind_cov2contrav(stat(four))
      
      call spatial_operator(stat(four), tend(four))
      
      tend(new)%phiG = (tend(one)%phiG + 2. * tend(two)%phiG + 2. * tend(three)%phiG + tend(four)%phiG) / 6.
      tend(new)%u    = (tend(one)%u    + 2. * tend(two)%u    + 2. * tend(three)%u    + tend(four)%u   ) / 6.
      tend(new)%v    = (tend(one)%v    + 2. * tend(two)%v    + 2. * tend(three)%v    + tend(four)%v   ) / 6.
      
      call update_stat      (stat_new, stat_old, tend(new), dt)
      call convert_wind_P2SP(stat_new)
      call fill_ghost       (stat_new)
      call unify_bdy_stat   (stat_new)
      call convert_wind_SP2P(stat_new)
      call convert_wind_cov2contrav(stat_new)

    end subroutine RK4
    
    subroutine update_stat(stat_new, stat_old, tend, inc_t)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(in   ) :: stat_old
      type(tend_field), intent(in   ) :: tend
      real            , intent(in   ) :: inc_t
      
      stat_new%phiG(ids:ide,jds:jde,ifs:ife) = stat_old%phiG(ids:ide,jds:jde,ifs:ife) + inc_t * tend%phiG(ids:ide,jds:jde,ifs:ife)
      stat_new%u   (ids:ide,jds:jde,ifs:ife) = stat_old%u   (ids:ide,jds:jde,ifs:ife) + inc_t * tend%u   (ids:ide,jds:jde,ifs:ife)
      stat_new%v   (ids:ide,jds:jde,ifs:ife) = stat_old%v   (ids:ide,jds:jde,ifs:ife) + inc_t * tend%v   (ids:ide,jds:jde,ifs:ife)
      
      stat_new%phi = stat_new%phiG / mesh%sqrtG
      
    end subroutine update_stat
    
    subroutine switch_stat(stat_old,stat_new)
      type(stat_field), intent(inout) :: stat_old
      type(stat_field), intent(in   ) :: stat_new
      
      stat_old%phiG = stat_new%phiG
      stat_old%phi  = stat_new%phi
      stat_old%u    = stat_new%u  
      stat_old%v    = stat_new%v
    end subroutine switch_stat
    
  end module time_scheme_mod