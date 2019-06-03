      subroutine phase4Vgam(r,p1,p2,p3,p4,p5,p6,wt,*)
c---- generate phase space for 2-->4 process with V+gamma+parton final state
c---- r(mxdim),p1(4),p2(4) are inputs reversed in sign from physical values 
c---- phase space for -p1-p2 --> p3+p4+p5+p6
c---- with all 2 pi's (ie 1/(2*pi)^8)
c----
c---- At present hardcoded for V = Z !!!!!!!!!!!!!!!!!
c---- At present hardcoded for V = Z !!!!!!!!!!!!!!!!!
c---- At present hardcoded for V = Z !!!!!!!!!!!!!!!!!
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
      include 'kprocess.f'
      real(dp):: r(mxdim)
      real(dp):: p1(4),p2(4),p3(4),p4(4),p5(4),p6(4)
      real(dp):: p12(4),p345(4),p34(4),s34min
      real(dp):: wt,wt3456,wt345,wt34,wt0
      real(dp) :: mass, width
      
      integer:: j
      parameter(wt0=1._dp/twopi**2)

      if (kcase == kWgajet) then
          mass = wmass
          width = wwidth
      elseif (kcase == kZgajet) then
          mass = zmass
          width = zwidth
      else
          print *, "phase4Vgam.f: bad kcase"
          stop
      endif

      wt=0._dp
      do j=1,4
      p12(j)=-p1(j)-p2(j)
      enddo
      s34min=wsqmin
      if (zerowidth) s34min=mass**2

      call phi1_2m_nobw(zip,r(1),r(2),r(3),s34min,p12,p6,p345,wt3456,*99)

      call phi1_2m_bw(zip,r(4),r(5),r(6),s34min,p345,p5,p34,mass,width,wt345,*99)
      
!      if (p5(4) == p5(4)) then
!        continue
!      else
!        write(6,*) 'p345',p345
!        write(6,*) 'p5',p5
!        write(6,*) 'p34',p34
!        write(6,*) 'r(1:6)',r(1:6)
!        stop
!      endif
      
      call phi3m0(r(7),r(8),p34,p3,p4,wt34,*99)

      wt=wt0*wt3456*wt345*wt34

      return
 99   wt=0._dp
      return 1
      end

      subroutine phase4Vgam_fsr(r,p1,p2,p3,p4,p5,p6,wt,*)
c---- generate phase space for 2-->4 process with V+gamma+parton final state
c---- r(mxdim),p1(4),p2(4) are inputs reversed in sign from physical values 
c---- phase space for -p1-p2 --> p3+p4+p5+p6
c---- with all 2 pi's (ie 1/(2*pi)^8)
c----
c---- At present hardcoded for V = Z !!!!!!!!!!!!!!!!!
c---- At present hardcoded for V = Z !!!!!!!!!!!!!!!!!
c---- At present hardcoded for V = Z !!!!!!!!!!!!!!!!!
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
      include 'kprocess.f'
      real(dp):: r(mxdim)
      real(dp):: p1(4),p2(4),p3(4),p4(4),p5(4),p6(4)
      real(dp):: p12(4),p345(4),p45(4),s34min
      real(dp):: wt,wt3456,wt345,wt34,wt0
      real(dp) :: mass, width
      
      integer:: j
      parameter(wt0=1._dp/twopi**2)

      if (kcase == kWgajet) then
          mass = wmass
          width = wwidth
      elseif (kcase == kZgajet) then
          mass = zmass
          width = zwidth
      else
          print *, "phase4Vgam.f: bad kcase"
          stop
      endif

      wt=0._dp
      do j=1,4
      p12(j)=-p1(j)-p2(j)
      enddo
      s34min=wsqmin
      if (zerowidth) s34min=mass**2

      call phi1_2m_bw(zip,r(1),r(2),r(3),s34min,p12,p6,p345,mass,width,wt3456,*99)

      call phi1_2m_nobw(zip,r(4),r(5),r(6),0._dp,p345,p3,p45,wt345,*99)

      call phi3m0(r(7),r(8),p45,p4,p5,wt34,*99)

      wt=wt0*wt3456*wt345*wt34

      return
 99   wt=0._dp
      return 1
      end

