! Solve the follow system of equation with Jacobi method. Use (0, 0, 0, 0) as initial guess and tolerance of 0.001
! 10w - x + 2y = 6
! - w + 11x - y + 3z = 25
! 2w - x + 10y - z = -11
! 3x - y + 8z = 15

program Jacobi
  implicit none

  real, allocatable :: A(:, :), b(:), xo(:), x(:)
  integer :: i, j, k, n, max
  real :: tol

  write (*, *) "Enter the number of equations, max iteration, tolerance:"
  read (*, *) n, max, tol

  allocate(A(n,n), b(n), xo(n), x(n)) ! Allocate necessary memory

  write (*, *) "Enter the coefficient matrix A:"
  read (*, *) ((A(i, j), j = 1, n), i = 1, n)
  write (*, *) "Enter the constant matrix b:"
  read (*, *) (b(i), i = 1, n)
  write (*, *) "Enter initial guess:"
  read (*, *) (x(i), i = 1, n)

  write(*, *) ! Adds a new line
  write (*, *) "Iteration  table"
  write (*, *) repeat("-", 100)

  k = 1
  do while (k <= max)
    do i = 1, n
      x(i) = (b(i) - sum(A(i, 1:i-1) * xo(1:i-1)) - sum(A(i, i+1:n) * xo(i+1:n))) / A(i, i)
    end do

    write (*, *) k, x
    if (norm2(x - xo) < tol) then
      write (*, *) repeat("-", 100)
      write (*, *) "Approximate solution: ", x
      stop
    end if
    k = k + 1
    xo = x
  end do
  write (*, *) repeat("-", 100)
  write (*, *) "Maximum number of iterations exceeded"
end program Jacobi
