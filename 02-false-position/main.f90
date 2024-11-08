! Employ the method of False Position to locate the root of `f(x) = cos(x) - x` in the interval [0,1] with an accuracy of 0.0001

program FalsePosition
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
    c = b - fb * (b - a) / (fb - fa)
    fc = f(c)

    write (*, "(i10, 5f20.8)") i, a, b, c, fc, abs(c-b)

    if (fc == 0.0 .or. abs(c-b) < tol) exit

    if (fb * fc < 0.0) then
      a = b; fa = fb
    end if
    b = c; fb = fc

  end do

  write (*, *) repeat("-", 110)
  write (*, *) "Approximate root at x =", c

  contains
    real function f(x)
      real, intent(in) :: x

      f = cos(x) - x
    end function f
end program FalsePosition