! Solve the initial value problem for the following differential equation using Euler's method with adaptive step size:
! dy/dx = y - x^2 + 1 where y(0) = 0.5, hmax = 1.0, hmin = 0.001, h = 0.1, tolerance = 0.001

program AdaptiveEulerMethod
  real :: x, x_end, y, y_one_step, y_half_step, y_two_step, h, h_min, h_max, h_new, tol, error
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
    y_one_step = y + h * f(x, y)

    y_half_step = y + (h / 2.0) * f(x, y)
    y_two_step = y_half_step + (h / 2.0) * f(x + h / 2.0, y_half_step)

    error = abs(y_two_step - y_one_step)
    
    if (error > tol) then 
      ! Error too large
      h_new = h * 0.8 ! Conservative reduction
      h = max(h_min, h_new) ! Ensure it doesn't fall below h_min
      cycle ! Recalculate same step with new h
    else
      ! Error acceptable, accept y_two_step
      y = y_two_step
      x = x + h
      n = n + 1
      x_points(n) = x
      y_points(n) = y
      h_points(n) = h

      h_new = h * 1.2 ! Conservative increase
      h = min(h_max, h_new) ! Ensure it doesn't go above h_max
      h = max(h, h_min) ! Ensure it doesn't fall below h_min
    end if

    if (x >= x_end) exit

  end do

  if (n >= max_steps) then
    write(*, *) "Warning: Maximum number of steps reached. Solution may be incomplete."
  end if

  write (*, *) "Adaptive Euler Method"
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
      f = y - x**2 + 1.0
    end function f
end program AdaptiveEulerMethod
