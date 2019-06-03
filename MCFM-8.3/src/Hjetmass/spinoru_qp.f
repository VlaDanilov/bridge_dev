      subroutine spinoru_qp(N,p,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'
      complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
      real*16 s(mxpart,mxpart)
      real*16 p(mxpart,4)
      integer N
      
      call spinoru_qp_s(N,p,za,zb,s)
      return
      end

      subroutine spinoru_qp_s(N,p,za,zb,s)
c---Calculate spinor products
c---extended to deal with negative energies ie with all momenta outgoing
c---Arbitrary conventions of Bern, Dixon, Kosower, Weinzierl, 
c---za(i,j)*zb(j,i)=s(i,j)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'
      complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
      real*16 s(mxpart,mxpart)
      real*16 p(mxpart,4), rt(mxpart)
      complex*32 c23(mxpart),f(mxpart)
      integer i,j,N
      
c---if one of the vectors happens to be zero this routine fails.
      do j=1,N
         za(j,j)=(0q0,0q0)
         zb(j,j)=za(j,j)

C-----positive energy case
         if (p(j,4) .gt. 0q0) then
            rt(j)=sqrt(p(j,4)+p(j,1))
            c23(j)=complex(p(j,3),-p(j,2))
            f(j)=(1q0,0q0)
         else
C-----negative energy case
            rt(j)=sqrt(-p(j,4)-p(j,1))
            c23(j)=complex(-p(j,3),p(j,2))
            f(j)=(0q0,1q0)
         endif
      enddo
      do i=2,N
         do j=1,i-1
         s(i,j)=2q0*(p(i,4)*p(j,4)-p(i,1)*p(j,1)
     &              -p(i,2)*p(j,2)-p(i,3)*p(j,3))
         za(i,j)=f(i)*f(j)
     &   *(c23(i)*complex(rt(j)/rt(i),0q0)-
     &    c23(j)*complex(rt(i)/rt(j),0q0))


         if (abs(s(i,j)).lt.1d-5) then
         zb(i,j)=-(f(i)*f(j))**2*conjg(za(i,j))
         else
         zb(i,j)=-complex(s(i,j),0q0)/za(i,j)
         endif

         za(j,i)=-za(i,j)
         zb(j,i)=-zb(i,j)
         s(j,i)=s(i,j)
         enddo
      enddo

      return
      end
