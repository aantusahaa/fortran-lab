! Estimate f(x) at x = 3.5 using Newton's backward difference interpolation  for the data points:
! (1,1), (2, 8), (3, 27), (4, 64)

program NewtonBackwardDifference
  implicit none
  integer :: n, i, j
  real :: p, h, u, term, result, x(10), y(10), table(10, 10)

  write (*, *) "Enter the number of data points:"
  read (*, *) n
  write (*, *) "Enter the data points (x, y):"
  read (*, *) (x(i), y(i), i = 1, n)
  
  do i = 1, n
    table(i, 1) = y(i)
  end do

  do j = 2, n
    do i = n, j, -1
      table(i, j) = table(i, j - 1) - table(i - 1, j - 1)
    end do
  end do

  write (*, *) "Backward Difference Table"
  write (*, *) repeat("-", n * 10)
  do i = 1, n
    do j = 1, i
      write(*, "(f10.4)", advance='no') table(i, j)
    end do
    write(*, *)
  end do
  write (*, *) repeat("-", n * 10)

  write(*, *) "Enter the interpolation point:"
  read (*, *) p

  h = x(2) - x(1)
  u = (p - x(n)) / h
  term = 1.0
  result = table(n, 1)

  do i = 2, n
    term = term * (u + (i - 2))
    result = result + (term * table(n, i) / factorial(i - 1))
  end do

  write (*, *) "Interpolated value at x =", p, "is", result

  contains
    integer function factorial(n)
      integer, intent(in) :: n
      integer :: i

      factorial = 1
      do i = 1, n
        factorial = factorial * i
      end do
    end function factorial
end program NewtonBackwardDifference