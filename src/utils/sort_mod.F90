module sort_mod

  implicit none

  private

  public qsort

  interface qsort
    module procedure qsort_integer
    module procedure qsort_real
  end interface qsort

contains

  recursive subroutine qsort_integer(x, left_idx_, right_idx_)

    integer, intent(inout) :: x(:)
    integer, intent(in), optional :: left_idx_
    integer, intent(in), optional :: right_idx_

    integer part_idx, i
    integer left_idx, right_idx

    left_idx  = merge(left_idx_, 1, present(left_idx_))
    right_idx = merge(right_idx_, size(x), present(right_idx_))

    ! Array x should be already sorted.
    if (left_idx >= right_idx) return

#ifndef NDEBUG
    if (size(x) <= 20) then
      do i = 1, size(x)
        write(*, '(A)', advance='no') '-----'
      end do
      write(*, *)
      do i = 1, size(x)
        write(*, '(I5)', advance='no') i
      end do
      write(*, *)
      do i = 1, size(x)
        if (i >= left_idx .and. i <= right_idx) then
          write(*, '(I5)', advance='no') x(i)
        else
          write(*, '("     ")', advance='no')
        end if
      end do
      write(*, *)
    end if
#endif

    ! Partition the array so that the left elements are smaller than the right ones.
    part_idx = partition(x, left_idx, right_idx)

#ifndef NDEBUG
    if (size(x) <= 20) then
      do i = 1, size(x)
        if (i >= left_idx .and. i <= right_idx) then
          write(*, '(I5)', advance='no') x(i)
        else
          write(*, '("     ")', advance='no')
        end if
      end do
      write(*, *)
      do i = 1, part_idx - 1
        write(*, '("     ")', advance='no')
      end do
      write(*, '("    ^")')
    end if
#endif

    call qsort_integer(x, left_idx, part_idx - 1)
    call qsort_integer(x, part_idx + 1, right_idx)

  contains

    integer function partition(x, left_idx, right_idx) result(part_idx)

      integer, intent(inout) :: x(:)
      integer, intent(inout) :: left_idx
      integer, intent(inout) :: right_idx

      integer tmp
      integer pivot_idx, i, j

      pivot_idx = left_idx
      i = left_idx
      j = -1
      do i = left_idx + 1, right_idx
        if (x(i) < x(pivot_idx) .and. j /= -1) then
          ! Swap the small element with the first large element.
          tmp = x(i); x(i) = x(j); x(j) = tmp
          ! Shift the recorded first large element.
          j = j + 1
        else if (i /= pivot_idx .and. x(i) >= x(pivot_idx) .and. j == -1) then
          ! Record the first element that is larger than or equal to pivot, called first large element.
          j = i
        end if
      end do
      if (j == -1) j = right_idx + 1
      j = j - 1
      if (j /= pivot_idx) then
        ! Swap pivot with the last small element.
        tmp = x(pivot_idx); x(pivot_idx) = x(j); x(j) = tmp
        part_idx = j
      else
        part_idx = pivot_idx
      end if

    end function partition

  end subroutine qsort_integer

  recursive subroutine qsort_real(x, left_idx_, right_idx_)

    real, intent(inout) :: x(:)
    integer, intent(in), optional :: left_idx_
    integer, intent(in), optional :: right_idx_

    integer part_idx, i
    integer left_idx, right_idx

    left_idx  = merge(left_idx_, 1, present(left_idx_))
    right_idx = merge(right_idx_, size(x), present(right_idx_))

    ! Array x should be already sorted.
    if (left_idx >= right_idx) return

#ifndef NDEBUG
    if (size(x) <= 20) then
      do i = 1, size(x)
        write(*, '(A)', advance='no') '--------'
      end do
      write(*, *)
      do i = 1, size(x)
        write(*, '(I8)', advance='no') i
      end do
      write(*, *)
      do i = 1, size(x)
        if (i >= left_idx .and. i <= right_idx) then
          write(*, '(F8.2)', advance='no') x(i)
        else
          write(*, '("        ")', advance='no')
        end if
      end do
      write(*, *)
    end if
#endif

    ! Partition the array so that the left elements are smaller than the right ones.
    part_idx = partition(x, left_idx, right_idx)

#ifndef NDEBUG
    if (size(x) <= 20) then
      do i = 1, size(x)
        if (i >= left_idx .and. i <= right_idx) then
          write(*, '(F8.2)', advance='no') x(i)
        else
          write(*, '("        ")', advance='no')
        end if
      end do
      write(*, *)
      do i = 1, part_idx - 1
        write(*, '("        ")', advance='no')
      end do
      write(*, '("       ^")')
    end if
#endif

    call qsort_real(x, left_idx, part_idx)
    call qsort_real(x, part_idx + 1, right_idx)

  contains

    integer function partition(x, left_idx, right_idx) result(part_idx)

      real, intent(inout) :: x(:)
      integer, intent(inout) :: left_idx
      integer, intent(inout) :: right_idx

      integer tmp
      integer pivot_idx, i, j

      pivot_idx = left_idx
      i = left_idx + 1
      j = -1
      do i = left_idx + 1, right_idx
        if (x(i) < x(pivot_idx) .and. j /= -1) then
          ! Swap the small element with the first large element.
          tmp = x(i); x(i) = x(j); x(j) = tmp
          ! Shift the recorded first large element.
          j = j + 1
        else if (x(i) >= x(pivot_idx) .and. j == -1) then
          ! Record the first element that is larger than or equal to pivot, called first large element.
          j = i
        end if
      end do
      if (j == -1) j = right_idx + 1
      j = j - 1
      if (j /= pivot_idx) then
        ! Swap pivot with the last small element.
        tmp = x(pivot_idx); x(pivot_idx) = x(j); x(j) = tmp
        part_idx = j - 1
      else
        part_idx = pivot_idx
      end if

    end function partition

  end subroutine qsort_real

end module sort_mod