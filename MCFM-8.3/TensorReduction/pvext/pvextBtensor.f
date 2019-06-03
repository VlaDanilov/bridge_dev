      subroutine pvextBtensor(p1,m1s,m2s,FB0,FB1,FB2)
      implicit none
C     p1 is the external momenta
C     m1s,m2s are the squares of the internal masses
      include 'types.f'
      include 'pvBnames.f'
      include 'pvextBv.f'
      include 'TRydef.f'
      include 'TRmetric.f'
      complex(dp)::FB0(-2:0),FB1(y1max,-2:0),FB2(y2max,-2:0)
      real(dp)::p1(4),p1Dp1,m1s,m2s
      integer n1,n2,B0i,pvextBcache
      logical,save:: first=.true.
!$omp threadprivate(first)

      if (first) then
      first=.false.
      call pvarraysetup
      endif


      p1Dp1=p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2
      B0i=pvextBcache(p1Dp1,m1s,m2s)
      FB0(:)=Bv(B0i+bb0,:)

      do n1=1,4
      FB1(n1,:)=Bv(B0i+bb1,:)*p1(n1)
      enddo


      do n1=1,4
      do n2=n1,4 
      FB2(y2(n1,n2),:)=p1(n1)*p1(n2)*Bv(B0i+bb11,:)
     . +g(n1,n2)*Bv(B0i+bb00,:)
      enddo
      enddo

      return
      end

