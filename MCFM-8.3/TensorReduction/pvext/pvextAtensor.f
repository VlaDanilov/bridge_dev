      subroutine pvextAtensor(m1s,FA0,FA1,FA2)
      implicit none
      include 'types.f'
      include 'TRconstants.f'
      include 'pvAnames.f'
      include 'pvextAv.f'
      include 'TRydef.f'
      include 'TRmetric.f'
      complex(dp):: FA0(-2:0),FA1(y1max,-2:0),FA2(y2max,-2:0)
      real(dp)::m1s
      integer n1,n2,A0i,pvextAcache
      logical,save:: first=.true.
!$omp threadprivate(first)
      if (first) then
      first=.false.
      call pvarraysetup
      endif

      A0i=pvextAcache(m1s)

      FA0(:)=Av(A0i+aa0,:)

      do n1=1,4
      FA1(n1,:)=czip
      enddo

      do n1=1,4
      do n2=n1,4 
      FA2(y2(n1,n2),:)=g(n1,n2)*Av(A0i+aa00,:)
      enddo
      enddo

      return
      end

