! Solve the following boundary value problem with Linear Shooting Method with N = 20, M = 10, TOL = 10^-5
! y'' = (1 / 8)(32 + 2x^3 - yy'), for 1<=x<=3 with y(1) = 17 and y(3) = 43/3

program NonlinearShooting
  implicit none

  real :: a, b, alpha, beta, h, x, tk, k(4, 2), k_p(4, 2), u1, u2, tol
  integer :: n, i, j, m
  real, allocatable :: w(:, :)

  write(*, *) "a, b, α, β, n, m, tol"
  read (*, *) a, b, alpha, beta, n, m, tol

  allocate(w(2, 0:n))

  h = (b - a) / n
  tk = (beta - alpha) / (b - a)

  do j = 1, m
    w(1, 0) = alpha
    w(2, 0) = tk
    u1 = 0
    u2 = 1

    do i = 0, n - 1
      x = a + (i) * h
  
      ! We will use RK4 for the approximation
      k(1, 1) = h * w(2, i)
      k(1, 2) = h * f(x, w(1, i), w(2, i))
      k(2, 1) = h * (w(2, i) + 0.5 * k(1, 2))
      k(2, 2) = h * f(x + h / 2, w(1, i) + 0.5 * k(1, 1), w(2, i) + 0.5 * k(1, 2))
      k(3, 1) = h * (w(2, i) + 0.5 * k(2, 2))
      k(3, 2) = h * f(x + h / 2, w(1, i) + 0.5 * k(2, 1), w(2, i) + 0.5 * k(2, 2))
      k(4, 1) = h * (w(2, i) + 0.5 * k(3, 2))
      k(4, 2) = h * f(x + h, w(1, i) + k(3, 1), w(2, i) + k(3, 2))
      w(1, i + 1) = w(1, i) + (1.0 / 6.0) * (k(1,1) + 2 * k(2, 1) + 2 * k(3, 1) + k(4, 1))
      w(2, i + 1) = w(2, i) + (1.0 / 6.0) * (k(1,2) + 2 * k(2, 2) + 2 * k(3, 2) + k(4, 2))
  
      k_p(1, 1) = h * u2
      k_p(1, 2) = h * (fy(x, w(1, i), w(2, i)) * u1 + fy_p(x, w(1, i), w(2, i)) * u2)
      k_p(2, 1) = h * (u2 + 0.5 * k_p(1, 2))
      k_p(2, 2) = h * (fy(x + h / 2, w(1, i), w(2, i)) * (u1 + 0.5 * k_p(1, 1)) &
                   + fy_p(x + h / 2, w(1, i), w(2, i)) * (u2 + 0.5 * k_p(1, 2)))
      k_p(3, 1) = h * (u2 + 0.5 * k_p(2, 2))
      k_p(3, 2) = h * (fy(x + h / 2, w(1, i), w(2, i)) * (u1 + 0.5 * k_p(2, 1)) &
                   + fy_p(x + h / 2, w(1, i), w(2, i)) * (u2 + 0.5 * k_p(2, 2)))
      k_p(4, 1) = h * (u2 + k_p(3, 2))
      k_p(4, 2) = h * (fy(x + h, w(1, i), w(2, i)) * (u1 + 0.5 * k_p(3, 1)) &
                   + fy_p(x + h, w(1, i), w(2, i)) * (u2 + 0.5 * k_p(3, 2)))
      u1 = u1 + (1.0 / 6.0) * (k_p(1,1) + 2 * k_p(2, 1) + 2 * k_p(3, 1) + k_p(4, 1))
      u2 = u2 + (1.0 / 6.0) * (k_p(1,2) + 2 * k_p(2, 2) + 2 * k_p(3, 2) + k_p(4, 2))
    end do

    if (abs(w(1, n) - beta) <= tol) then
      write (*, *) "x               | w1             | w2"
      write (*, *) repeat("-", 50)
      do i = 0, n
        x = a + i * h
        write(*, *) x, w(1, i), w(2, i)
      end do
      stop
    end if

    tk = tk - (w(1, n) - beta) / u1 !Newton's method
  end do

  write(*, *) "Maximum number of iterations exceeded"

  contains
    real function f(x, y , y_p)
      real, intent(in) :: x, y, y_p
      f = (1.0 / 8.0) * (32 + 2 * x**3 - y * y_p)
    end function f

    real function fy(x, y , y_p)
      real, intent(in) :: x, y, y_p
      fy = (1.0 / 8.0) * (- y_p)
    end function fy

    real function fy_p(x, y , y_p)
      real, intent(in) :: x, y, y_p
      fy_p = (1.0 / 8.0) * (- y)
    end function fy_p

end program NonlinearShooting