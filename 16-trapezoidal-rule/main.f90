! Use the Trapezoidal Rule to approximate f(x) = x^2 integrating within [0, 1] with n = 4

program TrapezoidalRule
  implicit none
  real :: a, b, h, sum, x
  integer :: i, n

  write(*, *) "Enter a, b and number of sub intervals:"
  read(*, *) a, b, n

  h = (b - a) / n
  sum = f(a) + f(b)

  do i = 1, n - 1
    x = a + i * h
    sum = sum + (2 * f(x))
  end do

  sum = sum * h / 2

  write(*, *) "Approximation of the integral is", sum

  contains
    real function f(x) 
      real, intent(in) :: x
      f = x ** 2
    end function f
end program TrapezoidalRule