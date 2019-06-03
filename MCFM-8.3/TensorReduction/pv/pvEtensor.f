      subroutine pvEtensor(q1,q2,q3,q4,m1s,m2s,m3s,m4s,m5s,
     &  FE0,FE1,FE2,FE3,FE4,FE5)
C****NB The arguments of this routine are the momentum offsets "q",
C       in the loop not the external momenta ******
C****   m1s,m2s,m3s,m4s,m5s are the internal masses squared.
      implicit none
      include 'types.f'
      include 'TRconstants.f'
      include 'TRydef.f'
      real(dp):: q1(4),q2(4),q3(4),q4(4)
      real(dp):: q1Dq1,q2Dq2,q3Dq3,q4Dq4,f1,f2,f3,f4
      real(dp):: p1(4),p2(4),p3(4),p4(4),p5(4),p234(4)
      complex(dp):: FE0(-2:0),FE1(y1max,-2:0),FE2(y2max,-2:0),
     & FE3(y3max,-2:0),FE4(y4max,-2:0),FE5(y5max,-2:0),
     & FD01(-2:0),FD11(y1max,-2:0),FD21(y2max,-2:0),
     & FD31(y3max,-2:0),FD41(y4max,-2:0),FD51(y5max,-2:0),
     & FD61(y6max,-2:0),
     & FD02(-2:0),FD12(y1max,-2:0),FD22(y2max,-2:0),
     & FD32(y3max,-2:0),FD42(y4max,-2:0),FD52(y5max,-2:0),
     & FD62(y6max,-2:0),
     & FD03(-2:0),FD13(y1max,-2:0),FD23(y2max,-2:0),
     & FD33(y3max,-2:0),FD43(y4max,-2:0),FD53(y5max,-2:0),
     & FD63(y6max,-2:0),
     & FD04(-2:0),FD14(y1max,-2:0),FD24(y2max,-2:0),
     & FD34(y3max,-2:0),FD44(y4max,-2:0),FD54(y5max,-2:0),
     & FD64(y6max,-2:0),
     & FD05(-2:0),FD15(y1max,-2:0),FD25(y2max,-2:0),
     & FD35(y3max,-2:0),FD45(y4max,-2:0),FD55(y5max,-2:0),
     & FD65(y6max,-2:0),
     & FD15a(y1max,-2:0),FD25a(y2max,-2:0),
     & FD35a(y3max,-2:0),FD45a(y4max,-2:0),
     & Gram(4,4)
      real(dp):: p1Dp1,p2Dp2,p3Dp3,p4Dp4,p5Dp5,
     & m1s,m2s,m3s,m4s,m5s,s12,s23,s34,s45,s51,
     & p12(4),p23(4),p34(4),p45(4),p51(4),pvdot,v1(4),v2(4),v3(4),v4(4)
      integer:: n1,n2,n3,n4,n5,nn,ep
      logical:: pvGramsing,singmat
      common/singmat/singmat
      logical,save:: first=.true.
!$omp threadprivate(first,/singmat/)

      if (first) then
        first=.false.
        call pvarraysetup
      endif

      p1(:)=q1(:)
      p2(:)=q2(:)-q1(:)
      p3(:)=q3(:)-q2(:)
      p4(:)=q4(:)-q3(:)
      p5(:)=-q4(:)
     
      p12(:)=p1(:)+p2(:)
      p23(:)=p2(:)+p3(:)
      p34(:)=p3(:)+p4(:)
      p45(:)=p4(:)+p5(:)
      p51(:)=p5(:)+p1(:)
      p234(:)=p2(:)+p3(:)+p4(:)

      p1Dp1=p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2
      p2Dp2=p2(4)**2-p2(1)**2-p2(2)**2-p2(3)**2
      p3Dp3=p3(4)**2-p3(1)**2-p3(2)**2-p3(3)**2
      p4Dp4=p4(4)**2-p4(1)**2-p4(2)**2-p4(3)**2
      p5Dp5=p5(4)**2-p5(1)**2-p5(2)**2-p5(3)**2
      s12=p12(4)**2-p12(1)**2-p12(2)**2-p12(3)**2
      s23=p23(4)**2-p23(1)**2-p23(2)**2-p23(3)**2
      s34=p34(4)**2-p34(1)**2-p34(2)**2-p34(3)**2
      s45=p45(4)**2-p45(1)**2-p45(2)**2-p45(3)**2
      s51=p51(4)**2-p51(1)**2-p51(2)**2-p51(3)**2
      call pvE0scalar(FE0,p1Dp1,p2Dp2,p3Dp3,p4Dp4,p5Dp5,
     & s12,s23,s34,s45,s51,m1s,m2s,m3s,m4s,m5s)

      Gram(1,1)=cmplx(two*p1Dp1,0._dp,kind=dp)
      Gram(2,2)=cmplx(two*p2Dp2,0._dp,kind=dp)
      Gram(3,3)=cmplx(two*p3Dp3,0._dp,kind=dp)
      Gram(4,4)=cmplx(two*p4Dp4,0._dp,kind=dp)
      Gram(1,2)=cmplx(two*pvdot(p1,p2),0._dp,kind=dp)
      Gram(1,3)=cmplx(two*pvdot(p1,p3),0._dp,kind=dp)
      Gram(1,4)=cmplx(two*pvdot(p1,p4),0._dp,kind=dp)
      Gram(2,3)=cmplx(two*pvdot(p2,p3),0._dp,kind=dp)
      Gram(2,4)=cmplx(two*pvdot(p2,p4),0._dp,kind=dp)
      Gram(3,4)=cmplx(two*pvdot(p3,p4),0._dp,kind=dp)
      Gram(2,1)=Gram(1,2)
      Gram(3,1)=Gram(1,3)
      Gram(4,1)=Gram(1,4)
      Gram(3,2)=Gram(2,3)
      Gram(4,2)=Gram(2,4)
      Gram(4,3)=Gram(3,4)

      singmat=pvGramsing(Gram,4)

c      if (singmat) then
c      write(6,*) 'Etensor_alt:singmat=',singmat
c      call Etensor_alt(q1,q2,q3,q4,FE0,FE1,FE2,FE3,FE4,FE5)
c      return

c      else
      call pvvcalc(q1,q2,q3,q4,v1,v2,v3,v4)
      call pvDtensor(q2,q3,q4,m1s,m3s,m4s,m5s,
     & FD01,FD11,FD21,FD31,FD41,FD51,FD61)
      call pvDtensor(q1,q3,q4,m1s,m2s,m4s,m5s,
     & FD02,FD12,FD22,FD32,FD42,FD52,FD62)
      call pvDtensor(q1,q2,q4,m1s,m2s,m3s,m5s,
     & FD03,FD13,FD23,FD33,FD43,FD53,FD63)
      call pvDtensor(q1,q2,q3,m1s,m2s,m3s,m4s,
     & FD04,FD14,FD24,FD34,FD44,FD54,FD64)
      call pvDtensor(p2,p23,p234,m2s,m3s,m4s,m5s,
     & FD05,FD15,FD25,FD35,FD45,FD55,FD65)
     
      q1Dq1=pvdot(q1,q1)
      q2Dq2=pvdot(q2,q2)
      q3Dq3=pvdot(q3,q3)
      q4Dq4=pvdot(q4,q4)

      f1=m2s-m1s-q1Dq1
      f2=m3s-m1s-q2Dq2
      f3=m4s-m1s-q3Dq3
      f4=m5s-m1s-q4Dq4

      do ep=-2,0
      do n1=1,4
      FE1(n1,ep)=
     & +half*(FD01(ep)-FD05(ep)+f1*FE0(ep))*v1(n1)
     & +half*(FD02(ep)-FD05(ep)+f2*FE0(ep))*v2(n1)
     & +half*(FD03(ep)-FD05(ep)+f3*FE0(ep))*v3(n1)
     & +half*(FD04(ep)-FD05(ep)+f4*FE0(ep))*v4(n1)
      enddo
      enddo


      call pvswitch1(q1,FD05,FD15,FD15a)
      do ep=-2,0
      do n1=1,4
      do n2=n1,4
      FE2(y2(n1,n2),ep)=
     & +half*(FD11(n1,ep)-FD15a(n1,ep)+f1*FE1(n1,ep))*v1(n2)
     & +half*(FD12(n1,ep)-FD15a(n1,ep)+f2*FE1(n1,ep))*v2(n2)
     & +half*(FD13(n1,ep)-FD15a(n1,ep)+f3*FE1(n1,ep))*v3(n2)
     & +half*(FD14(n1,ep)-FD15a(n1,ep)+f4*FE1(n1,ep))*v4(n2)
      enddo
      enddo
      enddo

      call pvswitch2(q1,FD05,FD15,FD25,FD25a)
      do ep=-2,0
      do n1=1,4
      do n2=n1,4
      do n3=n2,4
      nn=y2(n1,n2)
      FE3(y3(n1,n2,n3),ep)=
     & +half*(FD21(nn,ep)-FD25a(nn,ep)+f1*FE2(nn,ep))*v1(n3)
     & +half*(FD22(nn,ep)-FD25a(nn,ep)+f2*FE2(nn,ep))*v2(n3)
     & +half*(FD23(nn,ep)-FD25a(nn,ep)+f3*FE2(nn,ep))*v3(n3)
     & +half*(FD24(nn,ep)-FD25a(nn,ep)+f4*FE2(nn,ep))*v4(n3)
      enddo
      enddo
      enddo
      enddo


      call pvswitch3(q1,FD05,FD15,FD25,FD35,FD35a)
      do ep=-2,0
      do n1=1,4
      do n2=n1,4
      do n3=n2,4
      do n4=n3,4
      nn=y3(n1,n2,n3)
      FE4(y4(n1,n2,n3,n4),ep)=
     & +half*(FD31(nn,ep)-FD35a(nn,ep)+f1*FE3(nn,ep))*v1(n4)
     & +half*(FD32(nn,ep)-FD35a(nn,ep)+f2*FE3(nn,ep))*v2(n4)
     & +half*(FD33(nn,ep)-FD35a(nn,ep)+f3*FE3(nn,ep))*v3(n4)
     & +half*(FD34(nn,ep)-FD35a(nn,ep)+f4*FE3(nn,ep))*v4(n4)
      enddo
      enddo
      enddo
      enddo
      enddo



      call pvswitch4(q1,FD05,FD15,FD25,FD35,FD45,FD45a)
      do ep=-2,0
      do n1=1,4
      do n2=n1,4
      do n3=n2,4
      do n4=n3,4
      do n5=n4,4
      nn=y4(n1,n2,n3,n4)
      FE5(y5(n1,n2,n3,n4,n5),ep)=
     & +half*(FD41(nn,ep)-FD45a(nn,ep)+f1*FE4(nn,ep))*v1(n5)
     & +half*(FD42(nn,ep)-FD45a(nn,ep)+f2*FE4(nn,ep))*v2(n5)
     & +half*(FD43(nn,ep)-FD45a(nn,ep)+f3*FE4(nn,ep))*v3(n5)
     & +half*(FD44(nn,ep)-FD45a(nn,ep)+f4*FE4(nn,ep))*v4(n5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo


c      call pvswitch5(q1,FD05,FD15,FD25,FD35,FD45,FD55,FD55a)
c      do ep=-2,0
c      do n1=1,4
c      do n2=n1,4
c      do n3=n2,4
c      do n4=n3,4
c      do n5=n4,4
c      do n6=n5,4
c      nn=y5(n1,n2,n3,n4,n5)
c      FE6(y6(n1,n2,n3,n4,n5,n6),ep)=
c     & +half*(FD51(nn,ep)-FD55a(nn,ep)+f1*FE5(nn,ep))*v1(n6)
c     & +half*(FD52(nn,ep)-FD55a(nn,ep)+f2*FE5(nn,ep))*v2(n6)
c     & +half*(FD53(nn,ep)-FD55a(nn,ep)+f3*FE5(nn,ep))*v3(n6)
c     & +half*(FD54(nn,ep)-FD55a(nn,ep)+f4*FE5(nn,ep))*v4(n6)
c      enddo
c      enddo
c      enddo
c      enddo
c      enddo
c      enddo
c      enddo

c      endif
      return
      end

