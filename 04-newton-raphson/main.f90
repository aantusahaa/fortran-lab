! Employ the Newton-Raphson method to locate the root of `f(x) = cos(x) - x` with an accuracy of 0.0001

program NewtonRaphson
  implicit none
  real :: p0, p, tol
  integer :: i, max_iteration

  write (*, *) "Enter p0, tolerance, max iteration: "
  read (*, *) p0, tol, max_iteration

  write(*, *)
  write (*, "(a10, 3a20)") "Iteration", "p0", "p", "Tolerance"
  write (*, *) repeat("-", 70)

  do i = 1, max_iteration
    if (f_prime(p0) == 0) stop "Please enter a different initial approximation"

    p = p0 - f(p0) / f_prime(p0)

    write (*, "(i10, 3f20.8)") i, p0, p, abs(p - p0)

    if (abs(p - p0) < tol) exit

    p0 = p

  end do

  write (*, *) repeat("-", 70)
  write (*, *) "Approximate root at x =", p

  contains
    real function f(x)
      real, intent(in) :: x

      f = cos(x) - x
    end function f

    real function f_prime(x)
      real, intent(in) :: x

      f_prime = - sin(x) - 1
    end function f_prime
end program NewtonRaphson