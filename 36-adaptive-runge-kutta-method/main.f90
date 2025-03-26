! Solve the initial value problem for the following differential equation using Runge-Kutta method with adaptive step size:
! dy/dx = y.cos(x) where y(0) = 1, hmax = 2.0, hmin = 0.5, h = 0.1, tolerance = 0.001
! https://jonshiach.github.io/ODEs-book/_pages/2.5_Adaptive_step_size_control.html

program AdaptiveRungeKuttaMethod
  real :: x, x_end, y, yp, yp1, h, h_min, h_max, h_new, tol, error, k1, k2, k3, k4
  integer, parameter :: max_steps = 10000
  real :: x_points(max_steps), y_points(max_steps), h_points(max_steps)
  integer :: n, i
  character(45) :: separator = repeat("-", 45)

  write (*, *) "Enter a, b, y(a), h, hmin, hmax, tolerance"
  read (*, *) x, x_end, y, h, h_min, h_max, tol

  x_points(1) = x
  y_points(1) = y
  h_points(1) = h
  n = 1

  do while (x < x_end .and. n < max_steps)
    h = min(h, x_end - x)

    ! We will use the Bogacki-Shampine 2(3) embedded RK method
    ! Butcher Tableau:
    ! 0   | 
    ! 1/2 | 1/2
    ! 3/4 | 0    3/4
    ! 1   | 2/9  1/3 4/9
    !-----------------------
    !     | 2/9  1/3 4/9 0
    !     | 7/24 1/4 1/3 1/8

    k1 = f(x, y)
    k2 = f(x + h * (1.0 / 2.0), y + h * (1.0 / 2.0) * k1)
    k3 = f(x + h * (3.0 / 4.0), y + h * (3.0 / 4.0) * k2)
    k4 = f(x + h, y + h * (2.0 / 9.0 * k1 + 1.0 / 3.0 * k2 + 4.0 / 9.0 * k3))

    yp1 = y + h / 9.0 * (2 * k1 + 3 * k2 + 4 * k3)
    yp = y + h / 24.0 * (7 * k1 + 6 * k2 + 8 * k3 + 3 * k4)

    error = abs(yp1 - yp)
    
    if (error < tol) then
      ! Error acceptable, accept yp1
      y = yp1
      x = x + h
      n = n + 1
      x_points(n) = x
      y_points(n) = y
      h_points(n) = h
    end if

    h_new = 0.9 * (tol / error) ** (1.0/3.0)
    h = h * max(h_min, min(h_max, h_new)) ! Ensure new h doesn't go beyond h_max and h_min

    if (x >= x_end) exit

  end do

  if (n >= max_steps) then
    write(*, *) "Warning: Maximum number of steps reached. Solution may be incomplete."
  end if

  write (*, *) "Adaptive Runge-Kutta Method"
  write (*, *) separator
  write (*, *) "Step  |     x     |     y(x)    |   Step Size"
  write (*, *) separator
  
  do i = 1, n
    write (*, '(I4, 2X, F10.6, 2X, F13.8, 2X, F10.7)') i, x_points(i), y_points(i), h_points(i)
  end do
  write (*, *) separator

  contains
    real function f(x, y)
      real, intent(in) :: x, y
      f = y * cos(x)
    end function f
end program AdaptiveRungeKuttaMethod
