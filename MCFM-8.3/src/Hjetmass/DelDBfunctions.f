      subroutine DelDBfunctions(p1,p2,p3,msq,Bfxn)
      implicit none
C     p1,p2,p3 are the external momenta
C     msq are the squares of the masses in the propagators
      include 'types.f'
      include 'constants.f'
      include 'pvDnames.f'
      include 'pvCnames.f'
      include 'pvextDv.f'
      include 'pvextCv.f'
      complex(dp):: Bfxn(3)
      real(dp):: p1Dp1,p2Dp2,p3Dp3,p4Dp4,p2Dp3,p1Dp2,
     & s12,s13,s23,p1(4),p2(4),p3(4),p4(4),p12(4),p23(4),p13(4),msq
      integer:: pvextDcache,pvextCcache,D123,D231,D312,C1x23
      integer,parameter::a=1,b=2,c=3
      real(dp):: tiny=1.e-9_dp

C     p1,p2,p3,p4 are the external momenta
      p4(:)=-p1(:)-p2(:)-p3(:)
      p12(:)=p1(:)+p2(:)
      p13(:)=p1(:)+p3(:)
      p23(:)=p2(:)+p3(:)

      p1Dp1=p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2
      p2Dp2=p2(4)**2-p2(1)**2-p2(2)**2-p2(3)**2
      p3Dp3=p3(4)**2-p3(1)**2-p3(2)**2-p3(3)**2
      p4Dp4=p4(4)**2-p4(1)**2-p4(2)**2-p4(3)**2
!      if (abs(p1Dp1) < tiny) p1Dp1=zip
!      if (abs(p2Dp2) < tiny) p2Dp2=zip
!      if (abs(p3Dp3) < tiny) p3Dp3=zip
!      if (abs(p4Dp4) < tiny) p4Dp4=zip
      s12=p12(4)**2-p12(1)**2-p12(2)**2-p12(3)**2
      s13=p13(4)**2-p13(1)**2-p13(2)**2-p13(3)**2
      s23=p23(4)**2-p23(1)**2-p23(2)**2-p23(3)**2
      p2Dp3=half*(s23-p2Dp2-p3Dp3)
      p1Dp2=half*(s12-p1Dp1-p2Dp2)

      C1X23=pvextCcache(p1Dp1,s23,p4Dp4,msq,msq,msq)
      D123=pvextDcache(p1Dp1,p2Dp2,p3Dp3,p4Dp4,s12,s23,msq,msq,msq,msq)
      D231=pvextDcache(p2Dp2,p3Dp3,p1Dp1,p4Dp4,s23,s13,msq,msq,msq,msq)
      D312=pvextDcache(p3Dp3,p1Dp1,p2Dp2,p4Dp4,s13,s12,msq,msq,msq,msq)
      
      Bfxn(a)=half*p2Dp3*(Dv(D123+dd0,0)+Dv(D231+dd0,0)+Dv(D312+dd0,0))
     & -p1Dp2*(Dv(D123+dd0,0)+Dv(D123+dd1,0)
     &        +Dv(D231+dd3,0)+Dv(D312+dd2,0))
     & +Cv(C1x23+cc0,0)
     & +four*(Dv(D123+dd00,0)+Dv(D123+dd001,0)
     &       +Dv(D231+dd003,0)+Dv(D312+dd002,0))       

      Bfxn(b)=
     & Dv(D123+dd3,0)+Dv(D231+dd2,0)+Dv(D312+dd1,0)+Dv(D312+dd0,0)
     & +four*(
     & +Dv(D123+dd33,0)+Dv(D123+dd133,0)
     & +Dv(D231+dd23,0)+Dv(D231+dd223,0)
     & +Dv(D312+dd112,0)+two*Dv(D312+dd12,0)+Dv(D312+dd2,0))

      Bfxn(c)=-half*(Dv(D123+dd0,0)+Dv(D231+dd0,0)+Dv(D312+dd0,0))
     & +four*(
     & +Dv(D123+dd23,0)+Dv(D123+dd123,0)
     & +Dv(D231+dd23,0)+Dv(D231+dd123,0)
     & +Dv(D312+dd23,0)+Dv(D312+dd123,0))

      return
      end
