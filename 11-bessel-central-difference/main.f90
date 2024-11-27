! Estimate f(x) at x = 0.273 using Bessel's central difference interpolation  for the data points:
! (0.0, 0.381300), (0.1, 0.285603), (0.2, 0.190092), (0.3, 0.096327), (0.4, 0.008268), (0.5, −0.067725)
! Reference: https://phys.libretexts.org/Bookshelves/Astronomy__Cosmology/Celestial_Mechanics_(Tatum)/01%3A_Numerical_Methods/1.10%3A_1.10-Besselian_Interpolation

program BesselCentralDifference
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
  result = (table(m, 1) + table(m + 1, 1)) / 2.0

  et = 1.0; k = 1;
  ot = 1.0; l = 1;

  do i = 2, n
    if (mod(i, 2) == 0) then
      if (k == 1) then 
        et = et * (u - 0.5)
      else
        et = et * (u + (k - 2)) * (u - (k - 1))
      end if 
      d = table(m, i)
      result = result + (et * d / factorial(i - 1))
      k = k + 1
    else
      ot = ot * (u + (l - 1)) * (u - l)
      d = (table(m, i) + table(m - 1, i)) / 2
      result = result + (ot * d / factorial(i - 1))
      m = m - 1
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
end program BesselCentralDifference