module random_number_mod

  implicit none

  private

  public random_number_init
  public random_number_get

  integer n
  integer, allocatable :: seed(:)

  interface random_number_get
    module procedure random_number_get_double
    module procedure random_number_get_integer
  end interface random_number_get

contains

  subroutine random_number_init

    integer clock, i

    call random_seed(size=n)
    allocate (seed(n))

    call system_clock(count=clock)
    seed = clock + 37*[(i - 1, i=1, n)]

    call random_seed(put=seed)

  end subroutine random_number_init

  subroutine random_number_get_double(a, b, r)

    real(8), intent(in) :: a, b
    real(8), intent(out) :: r

    call random_number(r)
    r = (r*(b - a)) + a

  end subroutine random_number_get_double

  subroutine random_number_get_integer(a, b, r)

    integer, intent(in) :: a, b
    integer, intent(out) :: r(:)

    integer dimSize(1), i
    real, allocatable :: rand(:)

    dimSize = shape(r)

    allocate (rand(dimSize(1)))

    call random_number(rand)
    do i = 1, dimSize(1)
      r(i) = int(rand(i)*(b + 1 - a)) + a
    end do

  end subroutine random_number_get_integer

end module random_number_mod
