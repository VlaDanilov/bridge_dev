      subroutine zsvbksb(u,w,v,m,n,b,x) 
c--- Adapted from Numerical Recipes
C--- back substitution for real u,w,v and complex b
      implicit none
      include 'types.f'
      integer i,j,jj,m,n,nmax
      parameter (nmax=100) 
      real(dp):: u(m,n),w(n),v(n,n) 
      complex(dp):: b(m),x(n),s,tmp(nmax)

      do j=1,n 
        s=cmplx(0._dp,0._dp,kind=dp)
        if(w(j).ne.0d0)then 
          do i=1,m 
            s=s+u(i,j)*b(i) 
          enddo
          s=s/w(j) 
        endif 
        tmp(j)=s 
      enddo 

      do j=1,n 
        s=cmplx(0._dp,0._dp,kind=dp)
        do jj=1,n 
          s=s+v(j,jj)*tmp(jj) 
        enddo
        x(j)=s 
      enddo
      
      return 
      end 
