      subroutine phase5Vgam(r,p1,p2,p3,p4,p5,p6,p7,wt)
c---- generate phase space for 2-->5 process with V+gamma+parton+parton final state
c---- r(mxdim),p1(4),p2(4) are inputs reversed in sign from physical values 
c---- phase space for -p1-p2 --> p3+p4+p5+p6
c---- with all 2 pi's (ie 1/(2*pi)^8)
c----
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'limits.f'
      include 'mxdim.f'
      include 'zerowidth.f'
      include 'breit.f'
      include 'ipsgen.f'
      real(dp):: r(mxdim)
      real(dp):: p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4)
      real(dp):: p12(4),p345(4),p34(4),p45(4),p67(4),s34min
      real(dp):: wt,wt12,wt345,wt34,wt45,wt67,wt0
      real(dp) :: mass, width
      
      integer:: j
      parameter(wt0=1._dp/twopi**3)

      wt=0._dp
      do j=1,4
      p12(j)=-p1(j)-p2(j)
      enddo
      s34min=wsqmin
      if (zerowidth) s34min=mass**2
      
      if (ipsgen == 1) then
        call phi1_2(r(1),r(2),r(3),r(4),p12,p345,p67,wt12,*99)
        call phi1_2m_bw(zip,r(5),r(6),r(7),s34min,p345,p5,p34,zmass,zwidth,wt345,*99)
        call phi3m0(r(8),r(11),p34,p3,p4,wt34,*99)
        wt=wt34
      else
        call phi1_2bw(r(1),r(2),r(3),r(4),p12,p345,p67,zmass,zwidth,wt12,*99)
        call phi1_2m_nobw(zip,r(5),r(6),r(7),0._dp,p345,p3,p45,wt345,*99)
        call phi3m0(r(8),r(11),p45,p4,p5,wt45,*99)
        wt=wt45
      endif
      
      call phi3m0(r(12),r(13),p67,p6,p7,wt67,*99)

      wt=wt0*wt12*wt345*wt67*wt

      return
 99   wt=0._dp
      return
      
      end

