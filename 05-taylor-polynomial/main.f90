! Determine the third degree Taylor polynomial of e^x centered at x = 0 and approximate the value at x = 1

program TaylorPolynomial
  implicit none
  real :: x, sum
  integer :: n, i

  write (*, *) "Enter x, Degree n:"
  read (*, * ) x, n

  sum = 0.0

  do i = 0, n
    sum = sum + ((x ** i) / factorial(i)) ! e^x = 1 + x + (x^2 / 2!) + (x^3 / 3!) + ...
  end do

  write (*, *) "Taylor polynomial P(x) at x =", x, 'is', sum

  contains
  integer function factorial(n)
    integer, intent(in) :: n
    integer :: i

    factorial = 1
    do i = 1, n
      factorial = factorial * i
    end do
  end function factorial
end program TaylorPolynomial