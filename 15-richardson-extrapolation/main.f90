! Use the Central Divided Difference Formula to approximate f′(1) for f (x) = e^x with h = 0.1 and refine the result with Richardson Extrapolation
! Reference: https://math.libretexts.org/Bookshelves/Calculus/CLP-2_Integral_Calculus_(Feldman_Rechnitzer_and_Yeager)/04%3A_Appendices/4.03%3A_C%3A_More_About_Numerical_Integration/4.3.01%3A_C.1%3A_Richardson_Extrapolation

program RichardsonExtrapolation
  implicit none
  real :: h, x, result, a_h, a_h2
  integer :: k = 2 ! k represents the order of error. As we are using cdd here the order of error is 2. For Simpson's rule it would have been 4.

  write(*, *) "Enter approximation point (x) and initial step size (h):"
  read(*, *) x, h

  a_h = cdd(x, h)
  a_h2 = cdd(x, h / 2)
  result = (((2 ** k) * a_h2) - a_h) / ((2 ** k) - 1)

  write(*, *) "Extrapolated value:", result

  contains
    real function f(x) 
      real, intent(in) :: x
      f = exp(x)
    end function f

    real function cdd(x, h)
      real, intent(in) :: x, h
      cdd = (f(x + h) - f(x - h)) / (2 * h)
    end function cdd
end program RichardsonExtrapolation