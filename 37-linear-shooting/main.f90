! Solve the following boundary value problem with Linear Shooting Method with N = 10
! y'' = -(2/x)y' + (2/x^2)y + (sin(lnx)/x^2), for 1<=x<=2 with y(1) = 1 and y(2) = 2 

program LinearShooting
  implicit none

  real :: a, b, alpha, beta, h, x, k(4, 2), k_prime(4, 2), w, w1, w2
  integer :: n, i
  real, allocatable :: u(:, :), v(:, :)

  write(*, *) "a, b, α, β, n"
  read (*, *) a, b, alpha, beta, n

  allocate(u(2, 0:n), v(2, 0:n))

  h = (b - a) / N
  u(1, 0) = alpha
  u(2, 0) = 0
  v(1, 0) = 0
  v(2, 0) = 1

  do i = 0, n - 1
    x = a + i * h

    ! We will use RK4 for the approximation
    k(1, 1) = h * u(2, i)
    k(1, 2) = h * (p(x) * u(2, i) + q(x) * u(1, i) + r(x))
    k(2, 1) = h * (u(2, i) + 0.5 * k(1, 2))
    k(2, 2) = h * (p(x + h / 2) * (u(2, i) + 0.5 * k(1, 2)) + q(x + h / 2) * (u(1, i) + 0.5 * k(1, 1)) + r(x + h / 2))
    k(3, 1) = h * (u(2, i) + 0.5 * k(2, 2))
    k(3, 2) = h * (p(x + h / 2) * (u(2, i) + 0.5 * k(2, 2)) + q(x + h / 2) * (u(1, i) + 0.5 * k(2, 1)) + r(x + h / 2))
    k(4, 1) = h * (u(2, i) + 0.5 * k(3, 2))
    k(4, 2) = h * (p(x + h) * (u(2, i) + k(3, 2)) + q(x + h) * (u(1, i) + k(3, 1)) + r(x + h))
    u(1, i + 1) = u(1, i) + (1.0 / 6.0) * (k(1,1) + 2 * k(2, 1) + 2 * k(3, 1) + k(4, 1))
    u(2, i + 1) = u(2, i) + (1.0 / 6.0) * (k(1,2) + 2 * k(2, 2) + 2 * k(3, 2) + k(4, 2))

    k_prime(1, 1) = h * v(2, i)
    k_prime(1, 2) = h * (p(x) * v(2, i) + q(x) * v(1, i))
    k_prime(2, 1) = h * (v(2, i) + 0.5 * k_prime(1, 2))
    k_prime(2, 2) = h * (p(x + h / 2) * (v(2, i) + 0.5 * k_prime(1, 2)) + q(x + h / 2) * (v(1, i) + 0.5 * k_prime(1, 1)))
    k_prime(3, 1) = h * (v(2, i) + 0.5 * k_prime(2, 2))
    k_prime(3, 2) = h * (p(x + h / 2) * (v(2, i) + 0.5 * k_prime(2, 2)) + q(x + h / 2) * (v(1, i) + 0.5 * k_prime(2, 1)))
    k_prime(4, 1) = h * (v(2, i) + 0.5 * k_prime(3, 2))
    k_prime(4, 2) = h * (p(x + h) * (v(2, i) + k_prime(3, 2)) + q(x + h) * (v(1, i) + k_prime(3, 1)))
    v(1, i + 1) = v(1, i) + (1.0 / 6.0) * (k_prime(1,1) + 2 * k_prime(2, 1) + 2 * k_prime(3, 1) + k_prime(4, 1))
    v(2, i + 1) = v(2, i) + (1.0 / 6.0) * (k_prime(1,2) + 2 * k_prime(2, 2) + 2 * k_prime(3, 2) + k_prime(4, 2))
  end do 

  write (*, *) "x               | w1             | w2"
  write (*, *) repeat("-", 50)
  w = (beta - u(1, n)) / v(1, n)
  write(*, *) a, alpha, w

  do i = 1, n
    w1 = u(1, i) + w * v(1, i) ! Approximation of y
    w2 = u(2, i) + w * v(2, i) ! Approximation of y'
    x = a + i * h
    write(*, *) x, w1, w2
  end do 

  contains
    real function p(x)
      real, intent(in) :: x
      p = - 2 / x
    end function p

    real function q(x)
      real, intent(in) :: x
      q = 2 / x**2
    end function q

    real function r(x)
      real, intent(in) :: x
      r = sin(log(x)) / x**2
    end function r

end program LinearShooting