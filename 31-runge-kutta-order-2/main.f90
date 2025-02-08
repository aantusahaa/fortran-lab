! Solve the initial value problem for the following differential equation using Runge-Kutta method (2nd order):
! dy/dx = x^2 + y^2 where y(0) = 0

program RungeKuttaMethod
  implicit none

  real :: a, b, h, x, y, k1, k2
  integer :: n, i
  character(30) :: separator = repeat("-", 30)

  write (*, *) "Enter a, b, h, y(a):"
  read (*, *) a, b, h, y

  write (*, *) !prints new line
  write (*, *) "Second-Order Runge-Kutta Method"
  write (*, *) separator
  write (*, *) "Step |     x     |     y"
  write (*, *) separator

  n = int((b - a) / h) ! Calculate number of steps based on step size
  x = a
  
  do i = 1, n
    k1 = h * f(x, y)
    x = a + i * h
    k2 = h * f(x, y + k1)
    y = y + 0.5 * (k1 + k2)
    write (*, '(I4, 2F12.6)') i, x, y
  end do
  write (*, *) separator

  contains
    real function f(x, y)
      real, intent(in) :: x, y
      f = x**2 + y**2
    end function f
end program RungeKuttaMethod