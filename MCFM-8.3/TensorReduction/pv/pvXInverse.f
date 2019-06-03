************************************************************************
* pvInverse computes the inverse of a matrix.
* Input:
*   A: n-by-n matrix A
*   n: dimension of A
* Output:
*   A: mangled LU decomposition of A
*   Ainv: inverse of A
*   perm: permutation vector

      subroutine pvXInverse(A, Ainv, n, perm)
      implicit none
      include 'types.f'
      include 'TRconstants.f'
      integer::n, perm(n)
      complex(dp):: A(n,n), Ainv(n,n)
      integer:: i, j

      call XLUDecomp(A, n, perm)
      do i = 1, n
        do j = 1, n
          Ainv(j,i) = czip
        enddo
        Ainv(i,i) = cone
        call XLUBackSubst(A, n, perm, Ainv(1,i))
      enddo
      end
