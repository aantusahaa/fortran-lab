! Solve the following system of linear equation with Gaussian Elimination Method
! w - x + 2y - z = -8
! 2w - 2x + 3y - 3z = -20
! w + x + y = -2
! w - x + 4y + 3z = 4

program GaussianEliminationWithBackwardSubstitution
  implicit none

  integer :: n, i, j, k, p
  real, allocatable :: a(:,:), x(:)
  real :: m

  write (*, *) "Enter the number of equations:"
  read (*, *) n

  allocate(a(n, n + 1), x(n)) ! Allocate necessary memory

  write (*, *) "Enter the augmented matrix [A|b]:"
  read (*, *) ((A(i, j), j = 1, n + 1), i = 1, n)

  do i = 1, n - 1
    p = -1 ! Indicates we haven't found any pivot
    do k = i, n
      if (a(k, i) /= 0.0) then 
        p = k
        exit
      end if
    end do

    if (p == -1) stop "No unique solution exists" ! No pivot was found

    if (p /= i) then  ! Swap a(p) <-> a(i)
      do k = 1, n + 1
        m = a(p, k) ! Using m as temporary variable for swapping
        a(p, k) = a(i, k)
        a(i, k) = m
      end do
    end if

    do j = i + 1, n
      m = a(j, i) / a(i, i)
      a(j, i:n+1) = a(j, i:n+1) - m * a(i, i:n+1) ! a(j) - m * a(i) -> a(j)
    end do
  end do

  if (a(n, n) == 0) stop "No unique solution exists" ! Infinite solutions

  ! Backward Substitution
  x(n) = a(n, n + 1) / a(n, n)
  do i = n - 1, 1, -1
    x(i) = (a(i, n + 1) - sum(a(i, i+1:n) * x(i+1:n))) / a(i, i)
  end do

  write (*, *) "Solution Vector:", x

end program GaussianEliminationWithBackwardSubstitution 
