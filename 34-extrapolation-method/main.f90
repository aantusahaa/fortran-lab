! Solve the initial value problem for the following differential equation using Extrapolation method:
! dy/dx = 2x - y where y(0) = 1, hamx = 0.2, hmin = 0.01, tolerance = 10^-4

program ExtrapolationMethod
  implicit none

  real :: a, b, h_max, h_min, tol, ya, to, t, wo, w1, w2, w3, h, hk, v, q(30, 30), y(30)
  integer :: i, j, k, flag, nflag, nk(8) = (/2, 4, 6, 8, 12, 16, 24, 32/)
  character(30) :: separator = repeat("-", 30)

  write (*, *) "Enter a, b, y(a), h_max, h_min, tol:"
  read (*, *) a, b, ya, h_max, h_min, tol

  write (*, *) !prints new line
  write (*, *) "Extrapolation Method"
  write (*, *) separator
  write (*, *) "  x  |     y     |   h   |  k"
  write (*, *) separator

  to = a
  wo = ya
  h = h_max
  flag = 1 ! This flag will be used to exit the main loop at line 31

  do i = 1, 7
    do j = 1, i
      q(i, j) = (nk(i + 1) / nk(j))**2
    end do
  end do

  do while (flag == 1)
    k = 1; nflag = 0; ! When desired accuracy is achieved nflag is set to 1

    do while (k <= 8 .and. nflag == 0)
      hk = h / nk(k)
      t = to
      w2 = wo
      w3 = w2 + hk * f(t, w2) ! Euler's first step
      t = to + hk

      do j = 1, nk(k) - 1
        w1 = w2
        w2 = w3
        w3 = w1 + 2 * hk * f(t, w2) ! Midpoint method
        t = to + (j + 1) * hk
      end do

      y(k) = (w3 + w2 + hk * f(t, w3)) / 2.0 ! Endpoint correction

      if (k >= 2) then
        j = k; v = y(1)

        do while (j >= 2)
          y(j - 1) = y(j) + (y(j) - y(j - 1)) / (q(k - 1, j - 1) - 1) ! Extrapolation
          j = j - 1
        end do
        
        if (abs(y(1) - v) <= tol) nflag = 1 ! y1 is accepted as new w
      end if
      k = k + 1
    end do

    k = k - 1
    if (nflag == 0) then
      h = h / 2.0 ! New value for w rejected, decrease h
      if (h < h_min) then
        write (*, *) "h_min exceeded"
        flag = 0 ! True branch complete, go back to line 31
      end if
    else 
      wo = y(1) ! New value for w accepted
      to = to + h
      write (*, "(F5.3, F12.6, F8.3, I5)") to, wo, h, k
      
      if (to >= b) then
        flag = 0 ! Procedure completed successfully
      else if (to +  h > b) then
        h = b - to ! Terminate at t = b
      else if (k <= 3 .and. h < 0.5*h_max) then
        h = 2*h ! Increase the step size if possible
      end if
    end if 
  end do

  write (*, *) separator

  contains
    real function f(x, y)
      real, intent(in) :: x, y
      f = 2*x - y
    end function f
end program ExtrapolationMethod