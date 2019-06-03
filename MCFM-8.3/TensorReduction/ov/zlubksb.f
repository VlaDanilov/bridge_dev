      subroutine zlubksb(a,n,indx,bin,b)
c--- Adapted from Numerical Recipes
C--- back substitution for real a and complex b
c--- extended so that original vector b is not destroyed
      implicit none
      include 'types.f'
      integer n,indx(n)
      real(dp):: a(n,n)
      complex(dp):: b(n),bin(n),sum
      integer i,ii,j,ll
      b=bin
      ii=0
      do i=1,n
         ll=indx(i)
         sum=b(ll)
         b(ll)=b(i)
         if (ii.ne.0)then
            do j=ii,i-1
               sum=sum-cmplx(a(i,j),kind=dp)*b(j)
            enddo               ! 11
         else if (sum.ne.cmplx(0._dp,0._dp,kind=dp)) then
            ii=i
         endif
         b(i)=sum
      enddo                     ! 12
      do i=n,1,-1
         sum=b(i)
         if(i.lt.n)then
            do j=i+1,n
               sum=sum-cmplx(a(i,j),kind=dp)*b(j)
            enddo               ! 13
         endif
         b(i)=sum/cmplx(a(i,i),kind=dp)
      enddo                     ! 14
      return
      end
