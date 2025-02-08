! Solve the initial value problem for the following differential equation using Runge-Kutta method (4th order):
! dy/dx = y - x^2 + 1 where y(0) = 0.5

program RungeKuttaMethod
  implicit none

  real :: a, b, h, x, y, k1, k2, k3, k4
  integer :: n, i
  character(30) :: separator = repeat("-", 30)

  write (*, *) "Enter a, b, h, y(a):"
  read (*, *) a, b, h, y

  write (*, *) !prints new line
  write (*, *) "Fourth-Order Runge-Kutta Method"
  write (*, *) separator
  write (*, *) "Step |     x     |     y"
  write (*, *) separator

  n = int((b - a) / h) ! Calculate number of steps based on step size
  x = a
  
  do i = 1, n
    k1 = h * f(x, y)
    k2 = h * f(x + 0.5 * h, y + 0.5 * k1)
    k3 = h * f(x + 0.5 * h, y + 0.5 * k2)
    x = a + i * h
    k4 = h * f(x, y + k3)
    y = y + (1 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
    write (*, '(I4, 2F12.6)') i, x, y
  end do
  write (*, *) separator

  contains
    real function f(x, y)
      real, intent(in) :: x, y
      f = y - x**2 + 1
    end function f
end program RungeKuttaMethod