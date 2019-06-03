      subroutine DelDTfunctions(p1,p2,msq,FT)
      implicit none
C     p1,p2 are the external momenta
C     msq is the square of the masses in the propagators
      include 'types.f'
      include 'constants.f'
      include 'pvBnames.f'
      include 'pvCnames.f'
      include 'pvextCv.f'
      include 'pvextBv.f'
      complex(dp):: FL,FT
      real(dp):: p1Dp1,p2Dp2,p3Dp3,p1Dp2,
     & s12,s13,s23,p1(4),p2(4),p12(4),Gram,msq
      integer:: pvextBcache,pvextCcache,C1x2,B1,B2,B12
      real(dp):: tiny=1.e-9_dp

C     p1,p2 are the external momenta
      p12(:)=p1(:)+p2(:)

      p1Dp1=p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2
      p2Dp2=p2(4)**2-p2(1)**2-p2(2)**2-p2(3)**2
      p3Dp3=p12(4)**2-p12(1)**2-p12(2)**2-p12(3)**2
!      if (abs(p1Dp1) < tiny) p1Dp1=zip
!      if (abs(p2Dp2) < tiny) p2Dp2=zip
!      if (abs(p3Dp3) < tiny) p3Dp3=zip
      p1Dp2=half*(p3Dp3-p1Dp1-p2Dp2)
      Gram=p1Dp1*p2Dp2-p1Dp2**2
      C1X2=pvextCcache(p1Dp1,p2Dp2,p3Dp3,msq,msq,msq)
      B1=pvextBcache(p1Dp1,msq,msq)
      B2=pvextBcache(p2Dp2,msq,msq)
      B12=pvextBcache(p3Dp3,msq,msq)
      FL=-half/Gram*(
     & (two-three*p1Dp1*(p1Dp2+p2Dp2)/Gram)*(Bv(B1+bb0,0)-Bv(B12+bb0,0))
     &+(two-three*p2Dp2*(p1Dp2+p1Dp1)/Gram)*(Bv(B2+bb0,0)-Bv(B12+bb0,0)) ! changed factor p1Dp1 to p2Dp2
     &-(four*msq+p1Dp1+p2Dp2+p3Dp3-three*p1Dp1*p2Dp2*p3Dp3/Gram)
     & *Cv(C1x2+cc0,0)-two)
      FT=-half/Gram
     & *(p3Dp3*(Bv(B1+bb0,0)+Bv(B2+bb0,0)-two*Bv(B12+bb0,0)
     & -two*p1Dp2*Cv(C1x2+cc0,0))
     & +(p1Dp1-p2Dp2)*(Bv(B1+bb0,0)-Bv(B2+bb0,0)))-p1Dp2*FL
      return
      end
