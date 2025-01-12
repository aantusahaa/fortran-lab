! Use the Simpson's 3/8 Rule to approximate f(x) = x^2 integrating within [0, 1] with n = 6

program Simpsons38Rule
  implicit none
  real :: a, b, h, sum, x
  integer :: i, n

  write(*, *) "Enter a, b and number of sub intervals:"
  read(*, *) a, b, n

  if (mod(n, 3) /= 0) stop "Sub interval should be multiple of 3"

  h = (b - a) / n
  sum = f(a) + f(b)

  do i = 1, n - 1
    x = a + i * h
    if (mod(i, 3) == 0) then
      sum = sum + (2 * f(x))
    else
      sum = sum + (3 * f(x))
    end if 
  end do

  sum = sum * h * 3 / 8

  write(*, *) "Approximation of the integral is", sum

  contains
    real function f(x) 
      real, intent(in) :: x
      f = x ** 2
    end function f
end program Simpsons38Rule