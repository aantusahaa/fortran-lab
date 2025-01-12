! Use the Weddle's Rule to approximate f(x) = x^2 integrating within [0, 1] with n = 6

program WeddleRule
  implicit none
  real :: a, b, h, sum, x
  integer :: i, n

  write(*, *) "Enter a, b and number of sub intervals:"
  read(*, *) a, b, n

  if (mod(n, 6) /= 0) stop "Sub interval should be multiple of 6"

  h = (b - a) / n
  sum = f(a) + f(b)

  do i = 1, n - 1
    x = a + i * h
    if (mod(i, 6) == 0) then
      sum = sum + (2 * f(x))
    else if (mod(i, 3) == 0) then
      sum = sum + (6 * f(x))
    else if (mod(i, 2) == 0) then 
      sum = sum + f(x)
    else 
      sum = sum + (5 * f(x))
    end if 
  end do

  sum = sum * h * 3 / 10

  write(*, *) "Approximation of the integral is", sum

  contains
    real function f(x) 
      real, intent(in) :: x
      f = x ** 2
    end function f
end program WeddleRule