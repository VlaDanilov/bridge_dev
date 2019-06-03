************************************************************************
* Det computes the determinant of a matrix.
* Input:
*   A: n-by-n matrix A
*   n: dimension of A
* Output:
*   determinant of A
* Warning: A is overwritten

      function pvXDet(A, n)
      implicit none
      include 'types.f'
      include 'TRconstants.f'
      include 'pvNmax.f'
      complex(dp):: pvXDet
      integer:: n
      complex(dp):: A(n,n)
      integer:: i, perm(Nmax)

      call XLUDecomp(A, n, perm)
      pvXDet = cone
      do i = 1, n
        pvXDet = pvXDet*A(i,i)
        if( perm(i) .ne. i ) pvXDet = -pvXDet
      enddo
      end

