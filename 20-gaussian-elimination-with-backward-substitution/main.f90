! Solve the following system of linear equation with Gaussian Elimination Method
! w - x + 2y - z = -8
! 2w - 2x + 3y - 3z = -20
! w + x + y = -2
! w - x + 4y + 3z = 4

program GaussianEliminationWithBackwardSubstitution
  implicit none

  integer :: n, i, j, k, p
  real, allocatable :: A(:,:), x(:)
  real :: m

  write (*, *) "Enter the number of equations:"
  read (*, *) n

  allocate(A(n, n + 1), x(n)) ! Allocate necessary memory

  write (*, *) "Enter the augmented matrix [A|b]:"
  read (*, *) ((A(i, j), j = 1, n + 1), i = 1, n)

  do i = 1, n - 1
    p = -1 ! Indicates we haven't found any pivot
    do k = i, n
      if (A(k, i) /= 0.0) then 
        p = k
        exit
      end if
    end do

    if (p == -1) stop "No unique solution exists" ! No pivot was found

    if (p /= i) then  ! Swap A(p) <-> A(i)
      do k = 1, n + 1
        m = A(p, k) ! Using m as temporary variable for swapping
        A(p, k) = A(i, k)
        A(i, k) = m
      end do
    end if

    do j = i + 1, n
      m = A(j, i) / A(i, i)
      A(j, i:n+1) = A(j, i:n+1) - m * A(i, i:n+1) ! A(j) - m * A(i) -> A(j)
    end do
  end do

  if (A(n, n) == 0) stop "No unique solution exists" ! Infinite solutions

  ! Backward Substitution
  x(n) = A(n, n + 1) / A(n, n)
  do i = n - 1, 1, -1
    x(i) = (A(i, n + 1) - sum(A(i, i+1:n) * x(i+1:n))) / A(i, i)
  end do

  write (*, *) "Solution Vector:", x

end program GaussianEliminationWithBackwardSubstitution 
