! Estimate f(x) at x = 2.5 using Lagrange interpolation  for the data points:
! (1,2), (3, 6), (4, 5)

program LagrangeInterpolation
  implicit none
  integer :: n, i, j
  real :: p, term, result, x(10), y(10)

  write (*, *) "Enter the number of data points:"
  read (*, *) n
  write (*, *) "Enter the data points (x, y):"
  read (*, *) (x(i), y(i), i = 1, n)

  write(*, *) "Enter the interpolation point:"
  read (*, *) p

  result = 0.0
  do i = 1, n
    term = y(i)
    do j = 1, n
      if (j /= i) then
        term = term * (p - x(j)) / (x(i) - x(j))
      end if
    end do
    result = result + term
  end do

  write (*, *) "Interpolated value at x =", p, "is", result
end program LagrangeInterpolation