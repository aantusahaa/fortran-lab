! Estimate f(x) at x = 33 using Stirling's central difference interpolation  for the data points:
! (25, 0.25), (30, 0.3), (35, 0.33), (40, 0.37), (45, 0.43)
! Reference: https://theengineeringmaths.com/wp-content/uploads/2017/11/interpolation-web.pdf

program StirlingCentralDifference
  implicit none
  integer :: i, j, k, l, m, n
  real :: p, h, u, d, et, ot, result, x(10), y(10), table(10, 10)

  write (*, *) "Enter the number of data points:"
  read (*, *) n
  write (*, *) "Enter the data points (x, y):"
  read (*, *) (x(i), y(i), i = 1, n)
  
  do i = 1, n
    table(i, 1) = y(i)
  end do

  do j = 2, n
    do i = 1, n - j + 1
      table(i, j) = table(i + 1, j - 1) - table(i, j - 1)
    end do
  end do

  write (*, *) "Central Difference Table"
  write (*, *) repeat("-", n * 10)
  do i = 1, n
    do j = 1, n - i + 1
      write(*, "(f10.4)", advance='no') table(i, j)
    end do
    write(*, *)
  end do
  write (*, *) repeat("-", n * 10)

  write(*, *) "Enter the interpolation point:"
  read (*, *) p

  h = x(2) - x(1)
  m = (n + 1) / 2
  u = (p - x(m)) / h
  result = table(m, 1)

  et = 1.0; k = 1;
  ot = 1.0; l = 1;

  do i = 2, n
    if (mod(i, 2) == 0) then
      if (k == 1) then 
        et = et * u
      else
        et = et * (u ** 2 - (k - 1) ** 2 )
      end if 
      d = (table(m, i) + table(m - 1, i)) / 2
      result = result + (et * d / factorial(i - 1))
      m = m - 1
      k = k + 1
    else
      ot = ot * (u ** 2 - (l - 1) ** 2)
      d = table(m, i)
      result = result + (ot * d / factorial(i - 1))
      l = l + 1
    end if
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
end program StirlingCentralDifference