      subroutine pvextDtensor(p1,p2,p3,m1s,m2s,m3s,m4s,
     & FD0,FD1,FD2,FD3,FD4)
      implicit none
C     p1,p2,p3 are the external momenta
C     m1s,m2s,m3s,m4s are the squares of the masses in the propagators
      include 'types.f'
      include 'pvDnames.f'
      include 'pvextDv.f'
      include 'TRydef.f'
      include 'TRmetric.f'
      complex(dp):: FD0(-2:0),FD1(y1max,-2:0),FD2(y2max,-2:0),
     & FD3(y3max,-2:0),FD4(y4max,-2:0)
      real(dp):: p1Dp1,p2Dp2,p3Dp3,p4Dp4,
     & s12,s23,p1(4),p2(4),p3(4),p4(4),p12(4),p23(4),m1s,m2s,m3s,m4s
      real(dp):: pvSPK,pvSPKL,pvSPKK,pvSDDP,
     & pvSPKKK,pvSPPKK,pvSPPKL,pvSDDPP,pvSDDPK,pvSDDDD
      integer:: n1,n2,n3,n4,pvextDcache,D01
      logical,save:: first=.true.
!$omp threadprivate(first)

      if (first) then
      first=.false.
      call pvarraysetup
      endif
    
      call pvYcalc(p1,p2,p3,m1s,m2s,m3s,m4s)


C     p1,p2,p3,p4 are the external momenta
      p4(:)=-p1(:)-p2(:)-p3(:)
      p12(:)=p1(:)+p2(:)
      p23(:)=p2(:)+p3(:)

      p1Dp1=p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2
      p2Dp2=p2(4)**2-p2(1)**2-p2(2)**2-p2(3)**2
      p3Dp3=p3(4)**2-p3(1)**2-p3(2)**2-p3(3)**2
      p4Dp4=p4(4)**2-p4(1)**2-p4(2)**2-p4(3)**2
      s12=p12(4)**2-p12(1)**2-p12(2)**2-p12(3)**2
      s23=p23(4)**2-p23(1)**2-p23(2)**2-p23(3)**2


      D01=pvextDcache(p1Dp1,p2Dp2,p3Dp3,p4Dp4,s12,s23,m1s,m2s,m3s,m4s)
      
      FD0(:)=Dv(D01+dd0,:)

      do n1=1,4
      FD1(n1,:)=
     & +Dv(D01+dd1,:)*p1(n1)
     & +Dv(D01+dd2,:)*p2(n1)
     & +Dv(D01+dd3,:)*p3(n1)
      enddo


      do n1=1,4
      do n2=n1,4
      FD2(y2(n1,n2),:)=
     & +p1(n1)*p1(n2)*Dv(D01+dd11,:)
     & +p2(n1)*p2(n2)*Dv(D01+dd22,:)
     & +p3(n1)*p3(n2)*Dv(D01+dd33,:)
     & +pvSPK(n1,n2,p1,p2)*Dv(D01+dd12,:)
     & +pvSPK(n1,n2,p1,p3)*Dv(D01+dd13,:)
     & +pvSPK(n1,n2,p2,p3)*Dv(D01+dd23,:)
     & +g(n1,n2)*Dv(D01+dd00,:)
      enddo
      enddo


      do n1=1,4
      do n2=n1,4
      do n3=n2,4
      FD3(y3(n1,n2,n3),:)=
     & +p1(n1)*p1(n2)*p1(n3)*Dv(D01+dd111,:)
     & +p2(n1)*p2(n2)*p2(n3)*Dv(D01+dd222,:)
     & +p3(n1)*p3(n2)*p3(n3)*Dv(D01+dd333,:)
     & +pvSPKK(n1,n2,n3,p2,p1)*Dv(D01+dd112,:)
     & +pvSPKK(n1,n2,n3,p3,p1)*Dv(D01+dd113,:)
     & +pvSPKK(n1,n2,n3,p1,p2)*Dv(D01+dd122,:)
     & +pvSPKK(n1,n2,n3,p1,p3)*Dv(D01+dd133,:)
     & +pvSPKK(n1,n2,n3,p3,p2)*Dv(D01+dd223,:)
     & +pvSPKK(n1,n2,n3,p2,p3)*Dv(D01+dd233,:)
     & +pvSPKL(n1,n2,n3,p1,p2,p3)*Dv(D01+dd123,:)
     & +pvSDDP(n1,n2,n3,p1)*Dv(D01+dd001,:)
     & +pvSDDP(n1,n2,n3,p2)*Dv(D01+dd002,:)
     & +pvSDDP(n1,n2,n3,p3)*Dv(D01+dd003,:)
      enddo
      enddo
      enddo


      do n1=1,4
      do n2=n1,4
      do n3=n2,4
      do n4=n3,4
      FD4(y4(n1,n2,n3,n4),:)=
     & +p1(n1)*p1(n2)*p1(n3)*p1(n4)*Dv(D01+dd1111,:)
     & +p2(n1)*p2(n2)*p2(n3)*p2(n4)*Dv(D01+dd2222,:)
     & +p3(n1)*p3(n2)*p3(n3)*p3(n4)*Dv(D01+dd3333,:)
     & +pvSPKKK(n1,n2,n3,n4,p2,p1)*Dv(D01+dd1112,:)
     & +pvSPKKK(n1,n2,n3,n4,p3,p1)*Dv(D01+dd1113,:)
     & +pvSPKKK(n1,n2,n3,n4,p1,p2)*Dv(D01+dd1222,:)
     & +pvSPKKK(n1,n2,n3,n4,p3,p2)*Dv(D01+dd2223,:)
     & +pvSPKKK(n1,n2,n3,n4,p1,p3)*Dv(D01+dd1333,:)
     & +pvSPKKK(n1,n2,n3,n4,p2,p3)*Dv(D01+dd2333,:)

     & +pvSPPKK(n1,n2,n3,n4,p1,p2)*Dv(D01+dd1122,:)
     & +pvSPPKK(n1,n2,n3,n4,p1,p3)*Dv(D01+dd1133,:)
     & +pvSPPKK(n1,n2,n3,n4,p2,p3)*Dv(D01+dd2233,:)

     & +pvSPPKL(n1,n2,n3,n4,p1,p2,p3)*Dv(D01+dd1123,:)
     & +pvSPPKL(n1,n2,n3,n4,p2,p1,p3)*Dv(D01+dd1223,:)
     & +pvSPPKL(n1,n2,n3,n4,p3,p1,p2)*Dv(D01+dd1233,:)

     & +pvSDDPP(n1,n2,n3,n4,p1)*Dv(D01+dd0011,:)
     & +pvSDDPP(n1,n2,n3,n4,p2)*Dv(D01+dd0022,:)
     & +pvSDDPP(n1,n2,n3,n4,p3)*Dv(D01+dd0033,:)

     & +pvSDDPK(n1,n2,n3,n4,p1,p2)*Dv(D01+dd0012,:)
     & +pvSDDPK(n1,n2,n3,n4,p2,p3)*Dv(D01+dd0023,:)
     & +pvSDDPK(n1,n2,n3,n4,p1,p3)*Dv(D01+dd0013,:)

     & +pvSDDDD(n1,n2,n3,n4)*Dv(D01+dd0000,:)
      enddo
      enddo
      enddo
      enddo

      return
      end

