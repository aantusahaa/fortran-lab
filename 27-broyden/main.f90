! Solve the following nonlinear system of equations with Newton's method.
! x^2 + y - 3 = 0
! x - y + 1 = 0

program Broyden
  implicit none

  real :: fR(2), jR(2, 2), s(2), x(2), z(2), tol
  integer :: i, max

  write (*, *) "Enter the tolerance and max iteration"
  read (*, *) tol, max
  write (*, *) "Enter the initial guess"
  read (*, *) x(1), x(2)

  write(*, *) ! Adds a new line
  write (*, *) "Iteration  table"
  write (*, *) repeat("-", 50)

  do i = 1, max
    fR = F(x)
    jR = J(x)
    s = matmul(inverse(jR), -fR)
    x = x + s
    z = F(x) - fR
    jR = jR + outer_product(z - matmul(jR, s), s) / dot_product(s, s)

    write (*, *) i, x
    if (norm2(s) < tol) then
      write (*, *) repeat("-", 50)
      write (*, *) "Approximate solution: ", x
      stop
    end if 
  end do

  write (*, *) repeat("-", 50)
  write (*, *) "Maximum number of iterations exceeded"

  contains
    function F(v) result(val)
      real, intent(in) :: v(2)
      real :: x, y
      real :: val(2)

      x = v(1)
      y = v(2)
      val = (/ x**2 + y - 3, x - y + 1 /)
    end function F

    function J(v) result (val)
      real, intent(in) :: v(2)
      real :: x, y
      real :: val(2, 2)

      x = v(1)
      y = v(2)
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

    function outer_product(u,v) result(p)
      real,  intent(in) :: u(:), v(:)
      real :: p(size(u),size(v))

      integer :: n, m

      n = size(u)
      m = size(v)
      p = spread(u,dim=2,ncopies=m) * spread(v,dim=1,ncopies=n)
    end function
end program Broyden