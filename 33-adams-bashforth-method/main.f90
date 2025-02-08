! Solve the initial value problem for the following differential equation using Adams-Bashforth method:
! dy/dx = (e^-x) - y where y(0) = 1

program AdamsBashforthMethod
  implicit none

  real :: a, b, h, ya
  real, allocatable :: x(:), y(:)
  integer :: n, i, step
  character(30) :: separator = repeat("-", 30)

  write (*, *) "Enter a, b, h, y(a):"
  read (*, *) a, b, h, ya

  write (*, *) "Enter the Adam-Bashforth step to use (2 - 5):"
  read (*, *) step

  ! We will use precomputed formula for allowed steps
  if (step < 2 .or. step > 5) stop "Invalid step."

  write (*, *) !prints new line
  write (*, *) "Adams-Bashforth Method"
  write (*, *) separator
  write (*, *) "Step |     x     |     y"
  write (*, *) separator

  n = int((b - a) / h) ! Calculate number of steps based on step size
  allocate(x(0:n)); allocate(y(0:n)); ! Allocate required memory, the 0:n syntax assures the array index should start from 0 instead of 1
  x(0) = a
  y(0) = ya
  call calculateInitialApproximation(a, h, step - 1, x, y)

  ! Print the approximations from RK4
  do i = 0, step - 1
    write (*, '(I4, 2F12.6)') i, x(i), y(i)
  end do

  ! Calculate and print Adam-Bashforth approximations
  do i = step, n
    if (step == 2) then
      y(i) = y(i - 1) + (h / 2.0) * &
      (3 * f(x(i - 1), y(i - 1)) &
      - f(x(i - 2), y(i - 2)))
    else if (step == 3) then
      y(i) = y(i - 1) + (h / 12.0) * &
      ( 23 * f(x(i - 1), y(i - 1)) &
      - 16 * f(x(i - 2), y(i - 2)) &
      + 5 * f(x(i - 3), y(i - 3)))
    else if (step == 4) then 
      y(i) = y(i - 1) + (h / 24.0) * &
      ( 55 * f(x(i - 1), y(i - 1)) &
      - 59 * f(x(i - 2), y(i - 2)) &
      + 37 * f(x(i - 3), y(i - 3)) &
      - 9 * f(x(i - 4), y(i - 4)))
    else if (step == 5) then
      y(i) = y(i - 1) + (h / 720.0) * &
      ( 1901 * f(x(i - 1), y(i - 1)) &
      - 2774 * f(x(i - 2), y(i - 2)) &
      + 2616 * f(x(i - 3), y(i - 3)) &
      - 1274 * f(x(i - 4), y(i - 4)) &
      + 251 * f(x(i - 5), y(i - 5)))
    end if
    x(i) = a + i * h
    write (*, '(I4, 2F12.6)') i, x(i), y(i)
  end do
  write (*, *) separator

  contains
    real function f(x, y)
      real, intent(in) :: x, y
      f = exp(-x) - y
    end function f

    ! We will calculate the initial approximations with RK4
    subroutine calculateInitialApproximation(a, h, n, x, y)
      implicit none

      integer :: n
      real :: a, h, k1, k2, k3, k4, x(0:n), y(0:n)

      do i = 1, n
        k1 = h * f(x(i - 1), y(i - 1))
        k2 = h * f(x(i - 1) + 0.5 * h, y(i - 1) + 0.5 * k1)
        k3 = h * f(x(i - 1) + 0.5 * h, y(i - 1) + 0.5 * k2)
        x(i) = a + i * h
        k4 = h * f(x(i), y(i - 1) + k3)
        y(i) = y(i - 1) + (1 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
      end do
    end subroutine calculateInitialApproximation

end program AdamsBashforthMethod