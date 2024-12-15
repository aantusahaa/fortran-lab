! Use the Forward Difference Formula to approximate f′(1) for f (x) = x^2 with h = 0.1.

program ForwardDifference
  implicit none
  real :: h, x, result

  write(*, *) "Enter approximation point (x) and step size (h):"
  read(*, *) x, h

  result = (f(x + h) - f(x)) / h

  write(*, *) "Derivative at x =", x, "is", result
  contains
    real function f(x) 
      real, intent(in) :: x
      f = x ** 2
    end function f
end program ForwardDifference