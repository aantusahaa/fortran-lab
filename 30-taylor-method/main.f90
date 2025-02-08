! Solve the initial value problem for the following differential equation using Taylor's method (2nd order):
! dy/dx = y - x^2 where y(0) = 0

program EulerMethod
  implicit none

  real :: a, b, h, x, y
  integer :: n, i
  character(30) :: separator = repeat("-", 30)

  write (*, *) "Enter a, b, h, y(a):"
  read (*, *) a, b, h, y

  write (*, *) !prints new line
  write (*, *) "Second-Order Taylor's Method"
  write (*, *) separator
  write (*, *) "Step |     x     |     y"
  write (*, *) separator

  n = int((b - a) / h) ! Calculate number of steps based on step size
  x = a
  
  do i = 1, n
    y = y + h * (f(x, y) + h / 2.0 * f_prime(x, y))
    x = a + i * h
    write (*, '(I4, 2F12.6)') i, x, y
  end do
  write (*, *) separator

  contains
    real function f(x, y)
      real, intent(in) :: x, y
      f = y - x**2
    end function f

    real function f_prime(x, y)
      real, intent(in) :: x, y
      f_prime = f(x,y) - 2*x
    end function f_prime
end program EulerMethod