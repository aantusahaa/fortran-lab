! Employ the Bisection method to locate the root of `f(x) = cos(x) - x` in the interval [0,1] with an accuracy of 0.0001

program BisectionMethod
  implicit none
  real :: a, b, c, fa, fb, fc, tol
  integer :: i, max_iteration

  write (*, *) "Enter a, b, tolerance, max iteration: "
  read (*, *) a, b, tol, max_iteration

  fa = f(a); fb = f(b)

  if (fa * fb > 0.0) stop "No root in given interval"

  write(*, *)
  write (*, "(a10, 5a20)") "Iteration", "a", "b", "c", "f(c)", "Tolerance"
  write (*, *) repeat("-", 110)

  do i = 1, max_iteration
    c = a + (b - a) / 2.0
    fc = f(c)

    write (*, "(i10, 5f20.8)") i, a, b, c, fc, abs(b - a) / 2

    if (fc == 0.0 .or. abs(b - a) / 2 < tol) exit

    if (fa * fc > 0.0) then
      a = c; fa = fc
    else  
      b = c
    end if
  end do

  write (*, *) repeat("-", 110)
  write (*, *) "Approximate root at x =", c

  contains
    real function f(x)
      real, intent(in) :: x

      f = cos(x) - x
    end function f
end program BisectionMethod