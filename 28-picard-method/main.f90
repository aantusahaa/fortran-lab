! Solve the initial value problem for the following differential equation using Picard’s method:
! dy/dx = y^2 where y(0) = 1

program PicardMethod
  implicit none

  real :: x, x_end, step, y_i
  integer :: n, max_iterations, i
  character(len = 40) :: separator = repeat("-", 40)

  x_end = 1.0    
  step = 0.1
  max_iterations = 3

  ! Initial condition: y(0) = 1
  ! For Picard's method, we start with an initial guess y_0(x) = y(0) = 1
  ! As we cannot perform general integration with Fortran we will calculate y_i manually up to y_3.

  write (*, *) "Picard's Method for dy/dx = y^2, y(0) = 1"
  write (*, *) ! prints a new line
  write (*, *) ! prints a new line
  write (*, *) "Iteration |   x   |   y_i(x)  "
  write (*, *) separator

  do n = 1, max_iterations ! Loop for Picard iterations

    do i = 0, int(x_end / step) ! Loop for x values
      x = real(i) * step

      if (n == 1) then
        ! First iteration (n=0): y_1(x) = y(0) + integral from 0 to x of (y_0(t))^2 dt
        ! Using y_0(t) = 1, we get y_1(x) = 1 + integral from 0 to x of (1)^2 dt = 1 + x
        y_i = 1 + x

      else if (n == 2) then
        ! Second iteration (n=1): y_2(x) = y(0) + integral from 0 to x of (y_1(t))^2 dt
        ! Using y_1(t) = 1 + t, we get y_2(x) = 1 + integral from 0 to x of (1+t)^2 dt
        ! y_2(x) = 1 + [t + t^2 + t^3/3]
        y_i = 1 + x + x**2 + (x**3) / 3.0

      else if (n == 3) then
        ! Third iteration (n=2): y_3(x) = y(0) + integral from 0 to x of (y_2(t))^2 dt
        ! Using y_2(t) = 1 + t + t^2 + t^3/3, we get y_3(x) = 1 + integral from 0 to x of (1 + t + t^2 + t^3/3)^2 dt
        ! y_3(x) = 1 + [t + t^2 + t^3 + (2/3)t^4 + (1/3)t^5 + (1/9)t^6 + (1/63)t^7]
        y_i = 1 + x + x**2 + x**3 + (2.0 / 3.0) * x**4 + (1.0 / 3.0) * x**5 + (1.0 / 9.0) * x**6 + (1.0 / 63.0) * x**7
      endif

      write(*, '(I9, F7.1, 2F15.8)') n, x, y_i
      
    end do ! End loop for x values
    write(*, *) separator
  end do ! End loop for Picard iterations


end program PicardMethod