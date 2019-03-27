module timer_mod

  use hash_table_mod
  use log_mod

  implicit none

  private

  public timer_init
  public timer_start
  public timer_end
  public timer_duration

  type timer_type
    integer :: start_clock = 0.0d0
    integer :: end_clock   = 0.0d0
    real(8) :: duration      = 0.0d0
    integer :: count      = 0.0d0
  end type timer_type

  real(8) count_rate
  type(hash_table_type) timers

contains

  subroutine timer_init()

    call system_clock(count_rate=count_rate)
    count_rate = count_rate / 1.0d6
    timers = hash_table()

  end subroutine timer_init

  subroutine timer_start(name)

    character(*), intent(in) :: name

    type(timer_type), pointer :: timer

    if (timers%hashed(name)) then
      select type (value => timers%value(name))
      type is (timer_type)
        timer => value
        timer%duration = timer%duration + (timer%end_clock - timer%start_clock) / count_rate
        timer%count = timer%count + 1
        call system_clock(timer%start_clock)
      end select
    else
      allocate(timer)
      timer%count = timer%count + 1
      call system_clock(timer%start_clock)
      call timers%insert(name, timer)
    end if

  end subroutine timer_start

  subroutine timer_end(name)

    character(*), intent(in) :: name

    type(timer_type), pointer :: timer

    timer => get_timer(name)
    call system_clock(timer%end_clock)

    timer%duration = timer%duration + (timer%end_clock - timer%start_clock) / count_rate

  end subroutine timer_end

  real(8) function timer_duration(name) result(res)

    character(*), intent(in) :: name

    type(timer_type), pointer :: timer

    timer => get_timer(name)
    res = timer%duration

  end function timer_duration

  real(8) function timer_averaged_duration(name) result(res)

    character(*), intent(in) :: name

    type(timer_type), pointer :: timer

    timer => get_timer(name)
    res = timer%duration / timer%count

  end function timer_averaged_duration

  function get_timer(name) result(res)

    character(*), intent(in) :: name
    type(timer_type), pointer :: res

    if (timers%hashed(name)) then
      select type (value => timers%value(name))
      type is (timer_type)
        res => value
        return
      class default
        call log_error('Unknown timer ' // trim(name) // '!')
      end select
    end if

  end function get_timer

end module timer_mod