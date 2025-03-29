! Solve the following nonlinear system of equations with Fixed Point iteration method.
! x^2 + y - 3 = 0
! x - y + 1 = 0

program FixedPointIteration
  implicit none

  integer :: i, max
  real :: tol, x, y, x_new, y_new

  write (*, *) "Enter the tolerance and max iteration"
  read (*, *) tol, max
  write (*, *) "Enter the initial guess"
  read (*, *) x, y

  write(*, *) ! Adds a new line
  write (*, *) "Iteration  table"
  write (*, *) repeat("-", 50)

  do i = 1, max
    x_new = sqrt(3 - y)
    y_new = x + 1

    write (*, *) i, x_new, y_new
    if (abs(x_new - x) < tol .and. abs(y_new - y) < tol) then
      write (*, *) repeat("-", 50)
      write (*, *) "Approximate solution: ", x_new, y_new
      stop
    end if
    x = x_new
    y = y_new
  end do

  write (*, *) repeat("-", 50)
  write (*, *) "Maximum number of iterations exceeded"
end program FixedPointIteration