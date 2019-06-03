C --- bx1,...,bx4 are box corrections for qa -> qa with different topologies
C --- other box corrections for quark-induced processes can be obtained via the 
C --- crossing relations (for both weak-QCD and pure QCD boxes)
C --- for pure QCD boxes, set Mv = 0

C --- bx1 = Re[(box)_{t-channel}*(born)_{s-channel}]
      subroutine dijet_bx1(bx,ss,tt,uu,cpl,Mv,ep)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'scale.f'
      real(dp):: bx,ss,tt,uu,cpl(1:2,0:1),Mv,xI2,xI3,xI4,
     .     afi,aff,vfi,vff,fac
      integer ep

      afi = cpl(1,0)
      aff = cpl(2,0)
      vfi = cpl(1,1)
      vff = cpl(2,1)

      fac = 2._dp*gsq**2*esq/(4._dp*pi)**2*aveqq

      bx = 
     .     + 8._dp*(-2._dp*uu*(aff*afi + vff*vfi)*xI2(ss,Mv**2,0._dp,musq,
     .     ep) + 2._dp*uu*(aff*afi + vff*vfi)*xI2(tt,0._dp,0._dp,musq,ep) 
     .     + (aff*afi*(Mv**2 - tt + uu)*(tt + uu) + (3._dp*tt**2 + uu**2 
     .     + Mv**2*(tt + uu))*vff*vfi)*xI3(0._dp,0._dp,ss,Mv**2,0._dp,0._dp,
     .     musq,ep) - tt*(Mv**2 + tt - uu)*(aff*afi + vff*vfi)*xI3(0._dp,
     .     0._dp,tt,0._dp,0._dp,0._dp,musq,ep) - tt*(Mv**2 + tt - uu)*(aff*
     .     afi + vff*vfi)*xI3(0._dp,tt,0._dp,Mv**2,0._dp,0._dp,musq,ep) 
     .     + (aff*afi*(Mv**2 - tt + uu)*(tt + uu) + (3._dp*tt**2 + uu**2 
     .     + Mv**2*(tt + uu))*vff*vfi)*xI3(ss,0._dp,0._dp,Mv**2,0._dp,0._dp,
     .     musq,ep) + tt*(-aff*afi*(Mv**4 + 2._dp*Mv**2*tt - tt**2 
     .     + uu**2) - (Mv**4 + 2._dp*Mv**2*tt + 3._dp*tt**2 + uu**2)*vff*
     .     vfi)*xI4(0._dp,0._dp,0._dp,0._dp,ss,tt,Mv**2,0._dp,0._dp,0._dp,musq,
     .     ep))

      bx = fac*bx/ss

      end subroutine dijet_bx1



C --- bx2 = Re[(box)_{u-channel}*(born)_{s-channel}]
      subroutine dijet_bx2(bx,ss,tt,uu,cpl,Mv,ep)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'scale.f'
      real(dp):: bx,ss,tt,uu,cpl(1:2,0:1),Mv,xI2,xI3,xI4,
     .     afi,aff,vfi,vff,fac
      integer ep

      afi = cpl(1,0)
      aff = cpl(2,0)
      vfi = cpl(1,1)
      vff = cpl(2,1)

      fac = 2._dp*gsq**2*esq/(4._dp*pi)**2*aveqq

      bx = 
     .     - 8._dp*(2._dp*tt*(aff*afi - vff*vfi)*xI2(ss,Mv**2,0._dp,musq,
     .     ep) - 2._dp*tt*(aff*afi - vff*vfi)*xI2(uu,0._dp,0._dp,musq,ep) 
     .     - (aff*afi*(Mv**2 + tt - uu)*(tt + uu) - (tt**2 + 3._dp*uu**2 
     .     + Mv**2*(tt + uu))*vff*vfi)*xI3(0._dp,0._dp,ss,Mv**2,0._dp,0._dp,
     .     musq,ep) + uu*(Mv**2 - tt + uu)*(aff*afi - vff*vfi)*xI3(0._dp,
     .     0._dp,uu,0._dp,0._dp,0._dp,musq,ep) + uu*(Mv**2 - tt + uu)*(aff*
     .     afi - vff*vfi)*xI3(0._dp,uu,0._dp,Mv**2,0._dp,0._dp,musq,ep) 
     .     - (aff*afi*(Mv**2 + tt - uu)*(tt + uu) - (tt**2 + 3._dp*uu**2 
     .     + Mv**2*(tt + uu))*vff*vfi)*xI3(ss,0._dp,0._dp,Mv**2,0._dp,0._dp,
     .     musq,ep) - uu*(-aff*afi*(Mv**4 + tt**2 + 2._dp*Mv**2*uu 
     .     - uu**2) + (Mv**4 + tt**2 + 2._dp*Mv**2*uu + 3._dp*uu**2)*vff*
     .     vfi)*xI4(0._dp,0._dp,0._dp,0._dp,ss,uu,Mv**2,0._dp,0._dp,0._dp,musq,
     .     ep))

      bx = fac*bx/ss

      end subroutine dijet_bx2



C --- bx3 = Re[(box)_{t-channel}*(born)_{t-channel}]
      subroutine dijet_bx3(bx,ss,tt,uu,cpl,Mv,ep)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'scale.f'
      real(dp):: bx,ss,tt,uu,cpl(1:2,0:1),Mv,xI2,xI3,xI4,
     .     afi,aff,vfi,vff,fac,cs
      integer ep

      afi = cpl(1,0)
      aff = cpl(2,0)
      vfi = cpl(1,1)
      vff = cpl(2,1)

c      fac = 2._dp*gsq**2*esq/(4._dp*pi)**2*aveqq
      cs = -2._dp/3._dp
      fac = -cs*gsq**2*esq/(4._dp*pi)**2*aveqq

      bx = 
     .     + 8._dp*(aff*afi + vff*vfi)*(2._dp*uu*xI2(ss,Mv**2,0._dp,musq,
     .     ep) - 2._dp*uu*xI2(tt,0._dp,0._dp,musq,ep) - (tt**2 + uu**2 
     .     + Mv**2*(tt + uu))*xI3(0._dp,0._dp,ss,Mv**2,0._dp,0._dp,musq,ep) 
     .     + tt*(Mv**2 + tt - uu)*xI3(0._dp,0._dp,tt,0._dp,0._dp,0._dp,musq,
     .     ep) + tt*(Mv**2 + tt - uu)*xI3(0._dp,tt,0._dp,Mv**2,0._dp,0._dp,
     .     musq,ep) - (tt**2 + uu**2 + Mv**2*(tt + uu))*xI3(ss,0._dp,0._dp
     .     ,Mv**2,0._dp,0._dp,musq,ep) + tt*(Mv**4 + 2._dp*Mv**2*tt + tt**2 
     .     + uu**2)*xI4(0._dp,0._dp,0._dp,0._dp,ss,tt,Mv**2,0._dp,0._dp,0._dp,
     .     musq,ep))


      bx = fac*bx/tt


      end subroutine dijet_bx3



C --- bx4 = Re[(box)_{u-channel}*(born)_{t-channel}]
      subroutine dijet_bx4(bx,ss,tt,uu,cpl,Mv,ep)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'scale.f'
      real(dp):: bx,ss,tt,uu,cpl(1:2,0:1),Mv,xI2,xI3,xI4,
     .     afi,aff,vfi,vff,fac,cs
      integer ep

      afi = cpl(1,0)
      aff = cpl(2,0)
      vfi = cpl(1,1)
      vff = cpl(2,1)

c      fac = 2._dp*gsq**2*esq/(4._dp*pi)**2*aveqq
      cs = -2._dp/3._dp
      fac = -cs*gsq**2*esq/(4._dp*pi)**2*aveqq

      bx = 
     .     - 16._dp*uu**2*(aff*afi + vff*vfi)*(-xI3(0._dp,0._dp,ss,Mv**2
     .     ,0._dp,0._dp,musq,ep) - xI3(ss,0._dp,0._dp,Mv**2,0._dp,0._dp,musq,
     .     ep) + uu*xI4(0._dp,0._dp,0._dp,0._dp,ss,uu,Mv**2,0._dp,0._dp,0._dp,
     .     musq,ep))

      bx = fac*bx/tt

      end subroutine dijet_bx4



C***********************************************************************C
C The following are the alternative subroutines for type I - IV boxes, 
C which do not call QCDLoop Lib as the four routines above.
C***********************************************************************C

      subroutine dijet_bx1_new(bx,ss,tt,uu,cpl,Mv)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      real(dp):: bx,ss,tt,uu,cpl(1:2,0:1),Mv,c0fin,xxI2,xxI3,xxI4,
     .     afi,aff,vfi,vff,fac

      afi = cpl(1,0)
      aff = cpl(2,0)
      vfi = cpl(1,1)
      vff = cpl(2,1)

      fac = 2._dp*gsq**2*esq/(4._dp*pi)**2*aveqq

       bx = 
     .     + 8._dp*(-tt*(Mv**2 + tt - uu)*(aff*afi + vff*vfi)*c0fin(tt,
     .     Mv**2) - 2._dp*uu*(aff*afi + vff*vfi)*xxI2(ss,Mv**2) 
     .     + 2._dp*uu*(aff*afi + vff*vfi)*xxI2(tt,0._dp) + 2._dp*(aff*afi*(
     .     Mv**2 - tt + uu)*(tt + uu) + (3._dp*tt**2 + uu**2 + Mv**2*(tt 
     .     + uu))*vff*vfi)*xxI3(ss,Mv**2) - tt*(Mv**2 + tt - uu)*
     .     (aff*afi + vff*vfi)*xxI3(tt,0._dp) + tt*(-aff*afi*(Mv**4 
     .     + 2._dp*Mv**2*tt - tt**2 + uu**2) - (Mv**4 + 2._dp*Mv**2*tt 
     .     + 3._dp*tt**2 + uu**2)*vff*vfi)*xxI4(ss,tt,Mv**2))

      bx = fac*bx/ss

      end subroutine dijet_bx1_new



      subroutine dijet_bx2_new(bx,ss,tt,uu,cpl,Mv)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      real(dp):: bx,ss,tt,uu,cpl(1:2,0:1),Mv,c0fin,xxI2,xxI3,xxI4,
     .     afi,aff,vfi,vff,fac

      afi = cpl(1,0)
      aff = cpl(2,0)
      vfi = cpl(1,1)
      vff = cpl(2,1)

      fac = 2._dp*gsq**2*esq/(4._dp*pi)**2*aveqq

      bx = 
     .     + 8._dp*(-uu*(Mv**2 - tt + uu)*(aff*afi - vff*vfi)*c0fin(uu,
     .     Mv**2) - 2._dp*tt*(aff*afi - vff*vfi)*xxI2(ss,Mv**2) 
     .     + 2._dp*tt*(aff*afi - vff*vfi)*xxI2(uu,0._dp) + 2._dp*(aff*afi*(
     .     Mv**2 + tt - uu)*(tt + uu) - (tt**2 + 3._dp*uu**2 + Mv**2*(tt 
     .     + uu))*vff*vfi)*xxI3(ss,Mv**2) - uu*(Mv**2 - tt + uu)*
     .     (aff*afi - vff*vfi)*xxI3(uu,0._dp) + uu*(-aff*afi*(Mv**4 
     .     + tt**2 + 2._dp*Mv**2*uu - uu**2) + (Mv**4 + tt**2 
     .     + 2._dp*Mv**2*uu + 3._dp*uu**2)*vff*vfi)*xxI4(ss,uu,Mv**2))

      bx = fac*bx/ss

      end subroutine dijet_bx2_new



      subroutine dijet_bx3_new(bx,ss,tt,uu,cpl,Mv)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      real(dp):: bx,ss,tt,uu,cpl(1:2,0:1),Mv,c0fin,xxI2,xxI3,xxI4,
     .     afi,aff,vfi,vff,fac,cs

      afi = cpl(1,0)
      aff = cpl(2,0)
      vfi = cpl(1,1)
      vff = cpl(2,1)

      cs = -2._dp/3._dp
C --- the overall '-' from the fermion-loop in contraction
      fac = -cs*gsq**2*esq/(4._dp*pi)**2*aveqq

      bx = 
     .     + 8._dp*(aff*afi + vff*vfi)*(tt*(Mv**2 + tt - uu)*c0fin(tt,
     .     Mv**2) + 2._dp*uu*xxI2(ss,Mv**2) - 2._dp*uu*xxI2(tt,0._dp) 
     .     - 2._dp*(tt**2 + uu**2 + Mv**2*(tt + uu))*xxI3(ss,Mv**2) 
     .     + tt*(Mv**2 + tt - uu)*xxI3(tt,0._dp) + tt*(Mv**4 
     .     + 2._dp*Mv**2*tt + tt**2 + uu**2)*xxI4(ss,tt,Mv**2))

      bx = fac*bx/tt

      end subroutine dijet_bx3_new



      subroutine dijet_bx4_new(bx,ss,tt,uu,cpl,Mv)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      real(dp):: bx,ss,tt,uu,cpl(1:2,0:1),Mv,c0fin,xxI2,xxI3,xxI4,
     .     afi,aff,vfi,vff,fac,cs

      afi = cpl(1,0)
      aff = cpl(2,0)
      vfi = cpl(1,1)
      vff = cpl(2,1)

      cs = -2._dp/3._dp
C --- the overall '-' from the fermion-loop in contraction
      fac = -cs*gsq**2*esq/(4._dp*pi)**2*aveqq

      bx = 
     .     -16._dp*uu**2*(aff*afi + vff*vfi)*(-2._dp*xxI3(ss,Mv**2) 
     .     + uu*xxI4(ss,uu,Mv**2))

      bx = fac*bx/tt

      end subroutine dijet_bx4_new



C --- Analytical expressions for some special loop integrals
      function xxI2(psq,msq)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'scale.f'
      include 'epinv.f'
      include 'cplx.h'
      real(dp):: psq,msq,xxI2

      xxI2 = epinv + 2._dp

C --- assume psq /= 0
      if((msq == 0._dp) .or. (msq == psq)) then
         xxI2 = xxI2 + real(log(cplx1(-musq/psq)),dp)
      else
         xxI2 = xxI2 
     .        + (msq - psq)/psq*real(log(cplx1((msq - psq)/msq)),dp)
     .        + real(log(cplx1(musq/msq)),dp)
      end if

c      write(6,*) "comparison: "
c      write(6,*) "1 --- ", xI2(psq,msq,0._dp,musq,0) 
c     .     + epinv*xI2(psq,msq,0._dp,musq,-1)
c      write(6,*) "2 --- ", xxI2

      end function xxI2




      function xxI3(psq,msq)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'scale.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'cplx.h'
      real(dp):: psq,msq,ddilog,xxI3

      if(msq == 0._dp) then
         xxI3 = epinv*epinv2 + real(epinv*log(cplx1(-musq/psq)) 
     .        + 0.5_dp*log(cplx1(-musq/psq))**2,dp)
      else if(msq == psq) then
         xxI3 = - 0.5_dp*(epinv*epinv2 + epinv*log(musq/msq) + half*log(musq/msq)**2) 
     .        - pi**2/12._dp
      else
         xxI3 = real(epinv*log(cplx1(msq/(msq - psq))) 
     .        + log(cplx1(msq/(msq - psq)))**2
     .        + log(cplx1(musq/msq))*log(cplx1(msq/(msq - psq))),dp)
     .        + ddilog(psq/msq)
      end if

      xxI3 = xxI3/psq

c      write(6,*) "comparison: "
c      write(6,*) "1 --- ", xI3(0._dp,0._dp,psq,msq,0._dp,0._dp,musq,0) 
c     .     + epinv*xI3(0._dp,0._dp,psq,msq,0._dp,0._dp,musq,-1)
c     .     + epinv*epinv2*xI3(0._dp,0._dp,psq,msq,0._dp,0._dp,musq,-2)
c      write(6,*) "2 --- ", xxI3

      end function xxI3



C --- the finit triangle integral and it goes divergent when msq = 0
      function c0fin(psq,msq)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'scale.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'cplx.h'
      real(dp):: psq,msq,ddilog,xxI3,c0fin

      if(psq == -msq) then
         c0fin = ddilog(1._dp)
      else if(msq == 0._dp) then
         c0fin = xxI3(psq,msq)*psq
      else
C --- following expression is exact when psq/msq > -1
C --- there is phase difference when psq/msq < -1. 
C --- However, we just need the real part.
         c0fin = real(
     .        + log(cplx1(-psq/msq))*log(cplx1(1._dp + psq/msq)),dp)
     .          + ddilog(-psq/msq)
      end if

      c0fin = c0fin/psq

c      write(6,*) "comparison: "
c      write(6,*) "1 --- ", xI3(0._dp,psq,0._dp,msq,0._dp,0._dp,musq,0) 
c     .     + epinv*xI3(0._dp,psq,0._dp,msq,0._dp,0._dp,musq,-1)
c     .     + epinv*epinv2*xI3(0._dp,psq,0._dp,msq,0._dp,0._dp,musq,-2)
c      write(6,*) "2 --- ", c0fin
c      print*, "psq, msq: ", psq, msq, xxI3(psq,msq)

      end function c0fin



      function xxI4(p12,p23,msq)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'scale.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'cplx.h'
      real(dp):: p12,p23,msq,ddilog,xxI4

      if(msq == 0._dp) then
         xxI4 = 4._dp*epinv*epinv2 
     .        + real(2._dp*epinv*(log(cplx1(-musq/p12)) 
     .        + log(cplx1(-musq/p23)))
     .        + log(cplx1(-musq/p12))**2
     .        + log(cplx1(-musq/p23))**2
     .        - log(cplx1(p12/p23))**2,dp)
     .        - pi**2
      else 
         xxI4 = epinv*epinv2 
     .        + real(
     .        + epinv*(
     .        + log(cplx1(-musq/p23)) 
     .        + 2._dp*log(cplx1(msq/(msq - p12)))
     .        )
     .        + 0.5_dp*log(cplx1(-musq/p23))**2
     .        - 0.5_dp*log(cplx1(-msq/p23))**2
     .        + 2._dp*log(cplx1(-musq/p23))*log(cplx1(msq
     .        /(msq-p12))),dp)
     .        - 4._dp*ddilog(1._dp - msq/(msq - p12))
     .        - ddilog(1._dp + msq/p23) - pi**2/6._dp
      end if

      xxI4 = xxI4/(p12 - msq)/p23

c      write(6,*) "comparison: "
c      write(6,*) "1 --- ", 
c     .     + xI4(0._dp,0._dp,0._dp,0._dp,p12,p23,msq,0._dp,0._dp,0._dp,musq,0) 
c     .     + epinv*xI4(0._dp,0._dp,0._dp,0._dp,p12,p23,msq,0._dp,0._dp,0._dp,musq,-1)
c     .     + epinv*epinv2*xI4(0._dp,0._dp,0._dp,0._dp,p12,p23,msq,0._dp,0._dp,0._dp,
c     .     musq,-2)
c      write(6,*) "1.5 --- ", 
c     .     + xI4(0._dp,0._dp,0._dp,0._dp,p12,p23,0._dp,0._dp,0._dp,msq,musq,0) 
c     .     + epinv*xI4(0._dp,0._dp,0._dp,0._dp,p12,p23,0._dp,0._dp,0._dp,msq,musq,-1)
c     .     + epinv*epinv2*xI4(0._dp,0._dp,0._dp,0._dp,p12,p23,0._dp,0._dp,0._dp,msq,
c     .     musq,-2)
c      write(6,*) "2 --- ", xxI4
c      write(6,*) "msq: ", msq

      end function xxI4
