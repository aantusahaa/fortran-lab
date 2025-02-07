! Solve the initial value problem for the following differential equation using Euler's method:
! dy/dx = x + y where y(0) = 1

program EulerMethod
  implicit none

  real :: a, b, h, x, y
  integer :: n, i
  character(45) :: separator = repeat("-", 45)

  write (*, *) "Enter a, b, h, y(a):"
  read (*, *) a, b, h, y

  write (*, *) separator
  write (*, *) "Euler's Method for dy/dx = x + y, y(0) = 1"
  write (*, *) separator
  write (*, *) ! prints a new line
  write (*, *) "Step |     x     |     y"
  write (*, *) separator

  n = int((b - a) / h) ! Calculate number of steps based on step size
  x = a
  
  do i = 1, n
    y = y + h * f(x, y)
    x = a + i * h
    write (*, '(I4, 2F12.6)') i, x, y
  end do
  write (*, *) separator

  contains
    real function f(x, y)
      real, intent(in) :: x, y
      f = x + y
    end function f
end program EulerMethod