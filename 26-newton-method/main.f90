! Solve the following nonlinear system of equations with Newton's method.
! x^2 + y - 3 = 0
! x - y + 1 = 0

program Newton
  implicit none

  integer :: i, max
  real :: tol, x(2), y(2)

  write (*, *) "Enter the tolerance and max iteration"
  read (*, *) tol, max
  write (*, *) "Enter the initial guess"
  read (*, *) x(1), x(2)

  write(*, *) ! Adds a new line
  write (*, *) "Iteration  table"
  write (*, *) repeat("-", 50)

  do i = 1, max
    y = - (matmul(inverse(J(x(1), x(2))), F(x(1), x(2))))
    x = x + y

    write (*, *) i, x
    if (norm2(y) < tol) then
      write (*, *) repeat("-", 50)
      write (*, *) "Approximate solution: ", x
      stop
    end if
  end do

  write (*, *) repeat("-", 50)
  write (*, *) "Maximum number of iterations exceeded"

  contains
    function F(x, y) result(val)
      real, intent(in) :: x, y
      real :: val(2)

      val = (/ x**2 + y - 3, x - y + 1 /)
    end function F

    function J(x, y) result (val)
      real, intent(in) :: x, y
      real :: val(2, 2)

      val(1, 1) = 2*x
      val(1, 2) = 1
      val(2, 1) = 1
      val(2, 2) = -1
    end function J

    function inverse(A) result (Ainv)
      real, intent(in) :: A(2, 2)
      real :: Ainv(2, 2), detinv

      ! Calculate the inverse determinant of the matrix
      detinv = 1/(A(1,1)*A(2,2) - A(1,2)*A(2,1))

      ! Calculate the inverse of the matrix
      Ainv(1,1) = +detinv * A(2,2)
      Ainv(2,1) = -detinv * A(2,1)
      Ainv(1,2) = -detinv * A(1,2)
      Ainv(2,2) = +detinv * A(1,1)
    end function inverse
end program Newton