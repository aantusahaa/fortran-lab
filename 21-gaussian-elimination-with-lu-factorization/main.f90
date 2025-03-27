! Solve the following system of linear equation with Gaussian Elimination Method, using LU factorization
! w + x + 3z = 8
! 2w + x - y + z = 7
! 3w - x - y + 2z = 14
! - w + 2x + 3y - z = -7

program GaussianEliminationWithLUDecomposition
  implicit none

  integer :: n, i, j
  real, allocatable :: L(:,:), U(:,:), b(:), y(:), x(:)
  real :: m

  write (*, *) "Enter the number of equations:"
  read (*, *) n

  allocate(L(n,n), U(n, n), b(n), y(n), x(n)) ! Allocate necessary memory

  write (*, *) "Enter the coefficient matrix A:"
  read (*, *) ((U(i, j), j = 1, n), i = 1, n)
  write (*, *) "Enter the constant matrix b:"
  read (*, *) (b(i), i = 1, n)

  ! Initialize L as an identity matrix of n dimension
  l = 0
  do i = 1, n
    L(i, i) = 1
  end do

  do i = 1, n - 1
    do j = i + 1, n
      m = U(j, i) / U(i, i)
      L(j, i) = m
      U(j, i:n) = U(j, i:n) - m * U(i, i:n) ! U(j) - m * U(i) -> U(j)
    end do
  end do


  write (*, *) "Upper triangular matrix"
  call print_matrix(U, n)
  write(*, *)
  write (*, *) "Lower triangular matrix"
  call print_matrix(L, n)

  ! Forward Substitution
  y(1) = b(1) / L(1, 1)
  do i = 2, n
    y(i) = (b(i) - sum(L(i, 1:i-1) * y(1:i-1))) / L(i, i)
  end do

  ! Backward Substitution
  x(n) = y(n) / U(n, n)
  do i = n - 1, 1, -1
    x(i) = (y(i) - sum(U(i, i+1:n) * x(i+1:n))) / U(i, i)
  end do

  write (*, *) "Solution Vector:", x

  contains
    subroutine print_matrix(matrix, n)
      real, dimension(n,n), intent(in) :: matrix
      integer :: i, j, n

      do i = 1, n
        do j = 1, n
            write(*, '(F8.2)', advance="no") matrix(i, j)
        end do
        write(*, *)
      end do
  end subroutine print_matrix

end program GaussianEliminationWithLUDecomposition 
