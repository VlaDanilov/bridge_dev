! This file is part of SusHi.
! 
! It includes the real corrections to gluon fusion,
! which are called within sigma.f .
! 
C-{{{ coefficients

C-{{{ C1SM

      function C1SMh(s,t,u,m2)
      implicit none
      double complex C1SMh
      double precision s,t,u,m2
      double complex Integral3,Integral4
      C1SMh=(8.d0*(s + t + u) - 
     &  (-4.d0*m2 + s + t + u)*
     &   (2.d0*((t + u)*Integral3(s,s + t + u,m2)
     &        + (s + u)*Integral3(t,s + t + u,m2)
     &        + (s + t)*Integral3(u,s + t + u,m2))
     &   + s*u*Integral4(s,t,u,m2)
     &   + s*t*Integral4(s,u,t,m2)
     &     + t*u*Integral4(t,s,u,m2)))
      end

      function C1SMh_real(s,t,u,m2)
      implicit none
      double precision C1SMh_real
      double precision s,t,u,m2
      double complex C1SMh
      C1SMh_real = DReal(C1SMh(s,t,u,m2))
      end

      function C1SMh_abs(s,t,u,m2)
      implicit none
      double precision C1SMh_abs
      double precision s,t,u,m2
      double complex C1SMh
      C1SMh_abs = cdabs(C1SMh(s,t,u,m2))
      end

      function C1SMA(s,t,u,m2)
      implicit none
      double complex C1SMA
      double precision s,t,u,m2
      double complex Integral3,Integral4
      C1SMA = - (s + t + u) * (
     &     2 * ( (t + u) * Integral3(s,s+t+u,m2)
     &         + (s + u) * Integral3(t,s+t+u,m2)
     &         + (s + t) * Integral3(u,s+t+u,m2) )
     &   + s*u*Integral4(s,t,u,m2)
     &   + s*t*Integral4(s,u,t,m2)
     &   + t*u*Integral4(t,s,u,m2) )
      end

C-}}}
C-{{{ C1SUSY

      function C1SUSYh(s,t,u,m2)
      implicit none
      double complex C1SUSYh
      double precision s,t,u,m2
      double complex Integral3,Integral4
      C1SUSYh=-(2.d0*(s + t + u) + 
     &     m2*(2.d0*(
     &     (t + u)*Integral3(s,s + t + u,m2) + 
     &     (s + u)*Integral3(t,s + t + u,m2) + 
     &     (s + t)*Integral3(u,s + t + u,m2)) + 
     &     s*u*Integral4(s,t,u,m2) + 
     &     s*t*Integral4(s,u,t,m2) + 
     &     t*u*Integral4(t,s,u,m2)))
      end

C-}}}
C-{{{ C2SM

      function C2SMh(s,t,u,m2)
      implicit none
      double complex C2SMh
      double precision s,t,u,m2
      double complex Integral2,Integral3,Integral4
      C2SMh=(8.d0*(-((s*t)/(s + t)) + s**2/(s + u) + 
     &     t*u*(((2.d0*s + u)*Integral2(t,m2))/(s + u)**2 + 
     &     ((2.d0*s + t)*Integral2(u,m2))/(s + t)**2 - 
     &     ((2.d0*s + t)/(s + t)**2 + 
     &        (2.d0*s + u)/(s + u)**2)*
     &      Integral2(s + t + u,m2) + 
     &  ((t*Integral3(0.d0,t,m2) + 
     &       u*Integral3(0.d0,u,m2)))/(2.d0*s))) + 
     &  2.d0*((4.d0*m2 - s)*(t + u)*Integral3(s,s + t + u,m2) + 
     &  (-s**2 - 2.d0*(2.d0*m2 + t)*u - (2.d0*t*u**2)/s + 
     &     s*(4.d0*m2 + u - (8.d0*m2*u)/(s + u)))*
     &   Integral3(t,s + t + u,m2) + 
     &  ((-2.d0*t**2*u)/s - (-4.d0*m2+ s)*s + 
     &     t*(-4.d0*m2 - 2.d0*u + s - (8.d0*m2*s)/(s + t)))*
     &   Integral3(u,s + t + u,m2)) +  
     &  (4.d0*m2 - s)*s*u*Integral4(s,t,u,m2) + 
     &  (4.d0*m2 - s)*s*t*Integral4(s,u,t,m2) + 
     & t*u*(-12.d0*m2 + s - (4.d0*t*u)/s)*Integral4(t,s,u,m2))
      end

      function C2SMh_real(s,t,u,m2)
      implicit none
      double precision C2SMh_real
      double precision s,t,u,m2
      double complex C2SMh
      C2SMh_real = DReal(C2SMh(s,t,u,m2))
      end

      function C2SMh_abs(s,t,u,m2)
      implicit none
      double precision C2SMh_abs
      double precision s,t,u,m2
      double complex C2SMh
      C2SMh_abs = cdabs(C2SMh(s,t,u,m2))
      end

      function C2SMA(s,t,u,m2)
      implicit none
      double complex C2SMA
      double precision s,t,u,m2
      double complex Integral3,Integral4
      C2SMA = - s * (
     &     2 * ( (t + u) * Integral3(s,s+t+u,m2)
     &         + (s - u) * Integral3(t,s+t+u,m2)
     &         + (s - t) * Integral3(u,s+t+u,m2) )
     &   + s*u*Integral4(s,t,u,m2)
     &   + s*t*Integral4(s,u,t,m2)
     &   - t*u*Integral4(t,s,u,m2) )
      end

C-}}}
C-{{{ C2SUSY

      function C2SUSYh(s,t,u,m2)
      implicit none
      double complex C2SUSYh
      double precision s,t,u,m2
      double complex Integral2,Integral3,Integral4
      C2SUSYh=-(2.d0*(-((s*t)/(s + t)) + s**2/(s + u) + 
     &     t*u*(((2.d0*s + u)*Integral2(t,m2))/(s + u)**2 + 
     &     ((2.d0*s + t)*Integral2(u,m2))/(s + t)**2 - 
     &     ((2.d0*s + t)/(s + t)**2 + 
     &        (2.d0*s + u)/(s + u)**2)*
     &      Integral2(s + t + u,m2) + 
     &  ((t*Integral3(0.d0,t,m2) + 
     &       u*Integral3(0.d0,u,m2)))/(2.d0*s))) + 
     &  (4.d0*m2*(t + u)*Integral3(s,s + t + u,m2) + 
     &  (- 2.d0*t*u - (2.d0*t*u**2)/s + 
     &     s*4.d0*m2 - s*(8.d0*m2*u)/(s + u)-4.d0*m2*u)*
     &   Integral3(t,s + t + u,m2) + 
     &  ((-2.d0*t**2*u)/s + 4.d0*m2*s + 
     &     t*(-4.d0*m2 - 2.d0*u - (8.d0*m2*s)/(s + t)))*
     &   Integral3(u,s + t + u,m2))/2.d0 + 
     &  (4.d0*m2*s*u*Integral4(s,t,u,m2) + 
     &  4.d0*m2*s*t*Integral4(s,u,t,m2) + 
     & t*u*(-12.d0*m2 - (4.d0*t*u)/s)*Integral4(t,s,u,m2))/4.d0)
      end

C-}}}
C-{{{ C1eff

      function C1effh(s,t,u)
      implicit none
      double precision s,t,u,C1effh
      C1effh=4/3.d0*(s+t+u)**2
      end

      function C1effA(s,t,u)
      implicit none
      double precision s,t,u,C1effA
      C1effA=2*(s+t+u)**2
      end

C-}}}
C-{{{ C2eff

      function C2effh(s,t,u)
      implicit none
      double precision s,t,u,C2effh
      C2effh=4/3.d0*s**2
      end

      function C2effA(s,t,u)
      implicit none
      double precision s,t,u,C2effA
      C2effA=2*s**2
      end

C-}}}
C-{{{ ASM

      function ASMh(s,tpu,m2)
      implicit none
      double precision s,tpu,m2
      double complex ASMh
      double complex Integral2,Integral3
      ASMh=-3.d0*(2.d0*s*Integral2(s,m2)-2.d0*s*Integral2(s+tpu,m2)+
     &   (tpu)*(-2.d0 +
     &   (tpu - 4.d0*m2)*Integral3(s,s + tpu,m2)))/tpu**2
      end

      function ASMh_real(s,t,u,m2)
      implicit none
      double precision ASMh_real
      double precision s,t,u,m2
      double complex ASMh
      ASMh_real = DReal(ASMh(s,t+u,m2))
      end

      function ASMh_abs(s,t,u,m2)
      implicit none
      double precision ASMh_abs
      double precision s,t,u,m2
      double complex ASMh
      ASMh_abs = cdabs(ASMh(s,t+u,m2))
      end

      function ASMA(s,tpu,m2)
      implicit none
      double precision s,tpu,m2
      double complex ASMA
      double complex Integral2,Integral3
      ASMA=3.d0*Integral3(s,s + tpu,m2)
      end

C-}}}
C-{{{ ASUSY

      function ASUSYh(s,tpu,m2)
      implicit none
      double precision s,tpu,m2
      double complex ASUSYh
      double complex Integral2,Integral3 
      ASUSYh=3/4.d0*(2*s*Integral2(s,m2)-2*s*Integral2(s+tpu,m2)-
     &   (tpu)*(2 +
     &   4*m2*Integral3(s,s + tpu,m2)))/tpu**2
      end

C-}}}
C-{{{ ASMeff

      function ASMeffh()
      implicit none
      double precision ASMeffh
      ASMeffh=1.d0
      end

      function ASMeffA()
      implicit none
      double precision ASMeffA
      ASMeffA=-3/2.d0
      end

C-}}}
C-{{{ ASUSYeff

      function ASUSYeffh()
      implicit none
      double precision ASUSYeffh
      ASUSYeffh=1/8.d0
      end

C-}}}
C-{{{ ALOSM

      function ALOSMh(mh2,m2)
      implicit none
      double precision mh2,m2
      double complex ALOSMh
      double complex Integral3
      ALOSMh=3*(2+(4*m2-
     &mh2)*Integral3(0.d0,mh2,m2))/mh2
      end

      function ALOSMA(mh2,m2)
      implicit none
      double precision mh2,m2
      double complex ALOSMA
      double complex Integral3
      ALOSMA=-3*Integral3(0.d0,mh2,m2)
      end

C-}}}
C-{{{ ALOSMeff

      function ALOSMeffh()
      implicit none
      double precision ALOSMeffh
      ALOSMeffh=1.d0
      end

      function ALOSMeffA()
      implicit none
      double precision ALOSMeffA
      ALOSMeffA=3/2.d0
      end

C-}}}
C-{{{ ALOSUSY

      function ALOSUSYh(mh2,mb2)
      implicit none
      double complex ALOSUSYh
      double precision mh2,mb2
      double complex Integral3
      ALOSUSYh=-3/2.d0*(1+2*mb2*Integral3(0.d0,mh2,mb2))/mh2
      end

C-}}}
C-{{{ ALOSUSYeff

      function ALOSUSYeffh()
      implicit none
      double precision ALOSUSYeffh
      ALOSUSYeffh=1/8.d0
      end

C-}}}

C-}}}
C-{{{ squared amplidute for gg -> gh

      function AMPgg(s,t,u)
      implicit none
      double precision AMPgg,s,t,u
      double complex c1,c2,c3,c4
      double precision C1effh,C1effA,C2effh,C2effA,sAMPLO
      double complex C1SMh,C2SMh,C1SUSYh,C2SUSYh,AMPLOeff,
     &C1SMA,C2SMA

      include 'common-quark.f'
      include 'common-vars.f'

      c1 = 0.d0
      c2 = 0.d0
      c3 = 0.d0
      c4 = 0.d0

      if (pseudo.eq.1) then

C-{{{ pseudoscalar higgs

         if(gc.ne.0.d0) then
            c1 = c1 + mc2*gc*C1SMA(s,t,u,mc2)
            c2 = c2 + mc2*gc*C2SMA(s,t,u,mc2)
            c3 = c3 + mc2*gc*C2SMA(t,s,u,mc2)
            c4 = c4 + mc2*gc*C2SMA(u,s,t,mc2)
         endif

         if(gb.ne.0.d0) then
            c1 = c1 + mb2*gb*C1SMA(s,t,u,mb2)
            c2 = c2 + mb2*gb*C2SMA(s,t,u,mb2)
            c3 = c3 + mb2*gb*C2SMA(t,s,u,mb2)
            c4 = c4 + mb2*gb*C2SMA(u,s,t,mb2)
         endif

         if(gt.ne.0.d0) then
            if(htl) then
               c1 = c1 + gt*C1effA(s,t,u)
               c2 = c2 + gt*C2effA(s,t,u)
               c3 = c3 + gt*C2effA(t,s,u)
               c4 = c4 + gt*C2effA(u,s,t)
            else
               c1 = c1 + mt2*gt*C1SMA(s,t,u,mt2)
               c2 = c2 + mt2*gt*C2SMA(s,t,u,mt2)
               c3 = c3 + mt2*gt*C2SMA(t,s,u,mt2)
               c4 = c4 + mt2*gt*C2SMA(u,s,t,mt2)
            endif
         endif

         if(gbp.ne.0.d0) then
            c1 = c1 + mbp2*gbp*C1SMA(s,t,u,mbp2)
            c2 = c2 + mbp2*gbp*C2SMA(s,t,u,mbp2)
            c3 = c3 + mbp2*gbp*C2SMA(t,s,u,mbp2)
            c4 = c4 + mbp2*gbp*C2SMA(u,s,t,mbp2)
         endif

         if(gtp.ne.0.d0) then
            c1 = c1 + mtp2*gtp*C1SMA(s,t,u,mtp2)
            c2 = c2 + mtp2*gtp*C2SMA(s,t,u,mtp2)
            c3 = c3 + mtp2*gtp*C2SMA(t,s,u,mtp2)
            c4 = c4 + mtp2*gtp*C2SMA(u,s,t,mtp2)
         endif

C-}}}

      else

C-{{{ scalar higgs

         if(gc.ne.0.d0) then
            c1 = c1 + mc2*gc*C1SMh(s,t,u,mc2)
            c2 = c2 + mc2*gc*C2SMh(s,t,u,mc2)
            c3 = c3 + mc2*gc*C2SMh(t,s,u,mc2)
            c4 = c4 + mc2*gc*C2SMh(u,s,t,mc2)
         endif

         if(gb.ne.0.d0) then
            c1 = c1 + mb2*gb*C1SMh(s,t,u,mb2)
            c2 = c2 + mb2*gb*C2SMh(s,t,u,mb2)
            c3 = c3 + mb2*gb*C2SMh(t,s,u,mb2)
            c4 = c4 + mb2*gb*C2SMh(u,s,t,mb2)
         endif

         if(gt.ne.0.d0) then
            if(htl) then
               c1 = c1 + gt*C1effh(s,t,u)
               c2 = c2 + gt*C2effh(s,t,u)
               c3 = c3 + gt*C2effh(t,s,u)
               c4 = c4 + gt*C2effh(u,s,t)
            else
               c1 = c1 + mt2*gt*C1SMh(s,t,u,mt2)
               c2 = c2 + mt2*gt*C2SMh(s,t,u,mt2)
               c3 = c3 + mt2*gt*C2SMh(t,s,u,mt2)
               c4 = c4 + mt2*gt*C2SMh(u,s,t,mt2)
            endif
         endif

         if(gbp.ne.0.d0) then
            c1 = c1 + mbp2*gbp*C1SMh(s,t,u,mbp2)
            c2 = c2 + mbp2*gbp*C2SMh(s,t,u,mbp2)
            c3 = c3 + mbp2*gbp*C2SMh(t,s,u,mbp2)
            c4 = c4 + mbp2*gbp*C2SMh(u,s,t,mbp2)
         endif

         if(gtp.ne.0.d0) then
            c1 = c1 + mtp2*gtp*C1SMh(s,t,u,mtp2)
            c2 = c2 + mtp2*gtp*C2SMh(s,t,u,mtp2)
            c3 = c3 + mtp2*gtp*C2SMh(t,s,u,mtp2)
            c4 = c4 + mtp2*gtp*C2SMh(u,s,t,mtp2)
         endif

         if(gb1.ne.0.d0) then
            c1 = c1 + mbsb2*gb1*C1SUSYh(s,t,u,msb12)
     &              + mbsb2*gb2*C1SUSYh(s,t,u,msb22)
            c2 = c2 + mbsb2*gb1*C2SUSYh(s,t,u,msb12)
     &              + mbsb2*gb2*C2SUSYh(s,t,u,msb22)
            c3 = c3 + mbsb2*gb1*C2SUSYh(t,s,u,msb12)
     &              + mbsb2*gb2*C2SUSYh(t,s,u,msb22)
            c4 = c4 + mbsb2*gb1*C2SUSYh(u,s,t,msb12)
     &              + mbsb2*gb2*C2SUSYh(u,s,t,msb22)
         endif

         if(gt1.ne.0.d0) then
            c1 = c1 + mt2*(gt1*C1SUSYh(s,t,u,mst12)
     &                   + gt2*C1SUSYh(s,t,u,mst22))
            c2 = c2 + mt2*(gt1*C2SUSYh(s,t,u,mst12)
     &                   + gt2*C2SUSYh(s,t,u,mst22))
            c3 = c3 + mt2*(gt1*C2SUSYh(t,s,u,mst12)
     &                   + gt2*C2SUSYh(t,s,u,mst22))
            c4 = c4 + mt2*(gt1*C2SUSYh(u,s,t,mst12)
     &                   + gt2*C2SUSYh(u,s,t,mst22))
         endif

         if(gbp1.ne.0.d0) then
            c1 = c1 + mbp2*(gbp1*C1SUSYh(s,t,u,msbp12)
     &                    + gbp2*C1SUSYh(s,t,u,msbp22))
            c2 = c2 + mbp2*(gbp1*C2SUSYh(s,t,u,msbp12)
     &                    + gbp2*C2SUSYh(s,t,u,msbp22))
            c3 = c3 + mbp2*(gbp1*C2SUSYh(t,s,u,msbp12)
     &                    + gbp2*C2SUSYh(t,s,u,msbp22))        
            c4 = c4 + mbp2*(gbp1*C2SUSYh(u,s,t,msbp12)
     &                    + gbp2*C2SUSYh(u,s,t,msbp22))
         endif

         if(gtp1.ne.0.d0) then
            c1 = c1 + mtp2*(gtp1*C1SUSYh(s,t,u,mstp12)
     &                    + gtp2*C1SUSYh(s,t,u,mstp22))
            c2 = c2 + mtp2*(gtp1*C2SUSYh(s,t,u,mstp12)
     &                    + gtp2*C2SUSYh(s,t,u,mstp22))
            c3 = c3 + mtp2*(gtp1*C2SUSYh(t,s,u,mstp12)
     &                    + gtp2*C2SUSYh(t,s,u,mstp22))
            c4 = c4 + mtp2*(gtp1*C2SUSYh(u,s,t,mstp12)
     &                    + gtp2*C2SUSYh(u,s,t,mstp22))
         endif

C-}}}

      endif

      if(htl) then
         AMPgg=(cdabs(c1)**2+cdabs(c2)**2+cdabs(c3)**2+cdabs(c4)**2)
     &        /s/t/u*9/32.d0/(s+t+u)
      else
         AMPgg=(cdabs(c1)**2+cdabs(c2)**2+cdabs(c3)**2+cdabs(c4)**2)
     &        /s/t/u*9/32.d0/(s+t+u)
      endif
      end

C-}}}
C-{{{ squared amplidute for qg -> qh

      function AMPqg(s,t,u)
      implicit none
      double precision AMPqg,s,t,u
      double complex A
      double precision ASMeffh,ASMeffA,sAMPLO
      double complex ASMh,ASMA,ASUSYh,AMPLOeff

      include 'common-quark.f'
      include 'common-vars.f'      

      A = 0.d0

      if (pseudo.eq.1) then

C-{{{ pseudoscalar higgs

         if(gc.ne.0.d0) then
            A = A + mc2*gc*ASMA(t,s+u,mc2)
         endif
         if(gb.ne.0.d0) then
            A = A + mb2*gb*ASMA(t,s+u,mb2)
         endif
         if(gt.ne.0.d0) then
            if(htl) then
               A = A + gt*ASMeffA()
            else
               A = A + mt2*gt*ASMA(t,s+u,mt2)
            endif
         endif
         if(gbp.ne.0.d0) then
            A = A + mbp2*gbp*ASMA(t,s+u,mbp2)
         endif
         if(gtp.ne.0.d0) then
            A = A + mtp2*gtp*ASMA(t,s+u,mtp2)
         endif

C-}}}

      else

C-{{{ scalar higgs

         if(gc.ne.0.d0) then
            A = A + mc2*gc*ASMh(t,s+u,mc2)
         endif
         if(gb.ne.0.d0) then
            A = A + mb2*gb*ASMh(t,s+u,mb2)
         endif
         if(gt.ne.0.d0) then
            if(htl) then
               A = A + gt*ASMeffh()
            else
               A = A + mt2*gt*ASMh(t,s+u,mt2)
            endif
         endif
         if(gbp.ne.0.d0) then
            A = A + mbp2*gbp*ASMh(t,s+u,mbp2)
         endif
         if(gtp.ne.0.d0) then
            A = A + mtp2*gtp*ASMh(t,s+u,mtp2)
         endif
         if(gb1.ne.0.d0) then
            A = A + mbsb2*(gb1*ASUSYh(t,s+u,msb12)
     &                 + gb2*ASUSYh(t,s+u,msb22))
         endif
         if(gt1.ne.0.d0) then
            A = A + mt2*(gt1*ASUSYh(t,s+u,mst12)
     &                 + gt2*ASUSYh(t,s+u,mst22))
         endif
         if(gbp1.ne.0.d0) then
            A = A + mbp2*(gbp1*ASUSYh(t,s+u,msbp12)
     &                  + gbp2*ASUSYh(t,s+u,msbp22))
         endif
         if(gtp1.ne.0.d0) then
            A = A + mtp2*(gt1*ASUSYh(t,s+u,mstp12)
     &                  + gt2*ASUSYh(t,s+u,mstp22))
         endif

C-}}}

      endif

      if(htl) then
         AMPqg=-(s**2+u**2)/t*cdabs(A)**2/(s+t+u)
      else
         AMPqg=-(s**2+u**2)/t*cdabs(A)**2/(s+t+u)
      endif

      end

C-}}}
C-{{{ squared amplidute for qq -> gh

      function AMPqq(s,tpu)
      implicit none
      double precision AMPqq,s,tpu
      double complex A
      double precision ASMeffh,ASMeffA,sAMPLO
      double complex ASMh,ASMA,ASUSYh,AMPLOeff

      include 'common-quark.f'
      include 'common-vars.f'

      A = 0.d0

      if (pseudo.eq.1) then

C-{{{ pseudoscalar higgs

         if(gc.ne.0.d0) then
            A = A + mc2*gc*ASMA(s,tpu,mc2)
         endif
         if(gb.ne.0.d0) then
            A = A + mb2*gb*ASMA(s,tpu,mb2)
         endif
         if(gt.ne.0.d0) then
            if(htl) then
               A = A + gt*ASMeffA()
            else
               A = A + mt2*gt*ASMA(s,tpu,mt2)
            endif
         endif
         if(gbp.ne.0.d0) then
            A = A + mbp2*gbp*ASMA(s,tpu,mbp2)
         endif
         if(gtp.ne.0.d0) then
            A = A + mtp2*gtp*ASMA(s,tpu,mtp2)
         endif

C-}}}

      else

C-{{{ scalar higgs

         if(gc.ne.0.d0) then
            A = A + mc2*gc*ASMh(s,tpu,mc2)
         endif
         if(gb.ne.0.d0) then
            A = A + mb2*gb*ASMh(s,tpu,mb2)
         endif
         if(gt.ne.0.d0) then
            if(htl) then
               A = A + gt*ASMeffh()
            else
               A = A + mt2*gt*ASMh(s,tpu,mt2)
            endif
         endif
         if(gbp.ne.0.d0) then
            A = A + mbp2*gbp*ASMh(s,tpu,mbp2)
         endif
         if(gtp.ne.0.d0) then
            A = A + mtp2*gtp*ASMh(s,tpu,mtp2)
         endif
         if(gb1.ne.0.d0) then
            A = A + mbsb2*(gb1*ASUSYh(s,tpu,msb12)
     &                   + gb2*ASUSYh(s,tpu,msb22))
         endif
         if(gt1.ne.0.d0) then
            A = A + mt2*(gt1*ASUSYh(s,tpu,mst12)
     &                 + gt2*ASUSYh(s,tpu,mst22))
         endif
         if(gbp1.ne.0.d0) then
            A = A + mbp2*(gbp1*ASUSYh(s,tpu,msbp12)
     &                  + gbp2*ASUSYh(s,tpu,msbp22))
         endif
         if(gtp1.ne.0.d0) then
            A = A + mtp2*(gtp1*ASUSYh(s,tpu,mstp12)
     &                  + gtp2*ASUSYh(s,tpu,mstp22))
         endif

C-}}}

      endif

      if(htl) then
         AMPqq=cdabs(A)**2
      else
         AMPqq=cdabs(A)**2
      endif
      end

C-}}}
C-{{{ squared amplidute for gg -> h

      function sAMPLO(mhiggs2)
      implicit none
      double precision sAMPLO,mhiggs2
      double complex AMP
c      double precision ALOSMeffh
      double complex ALOSMh,ALOSMA,ALOSUSYh

      include 'common-quark.f'
      include 'common-vars.f'

      AMP = 0.d0

      if (pseudo.eq.1) then

C-{{{ pseudoscalar higgs

         if(gc.ne.0.d0) then
            AMP = AMP + mc2*gc*ALOSMA(mhiggs2,mc2)
         endif
         if(gb.ne.0.d0) then
            AMP = AMP + mb2*gb*ALOSMA(mhiggs2,mb2)
         endif
         if(gt.ne.0.d0) then
            AMP = AMP + mt2*gt*ALOSMA(mhiggs2,mt2)
     &                   !+ gt*ALOSMeffA(mhiggs2)
         endif
         if(gbp.ne.0.d0) then
            AMP = AMP + mbp2*gbp*ALOSMA(mhiggs2,mbp2)
         endif
         if(gtp.ne.0.d0) then
            AMP = AMP + mtp2*gtp*ALOSMA(mhiggs2,mtp2)
         endif

C-}}}

      else

C-{{{ scalar higgs

         if(gc.ne.0.d0) then
            AMP = AMP + mc2*gc*ALOSMh(mhiggs2,mc2)
         endif
         if(gb.ne.0.d0) then
            AMP = AMP + mb2*gb*ALOSMh(mhiggs2,mb2)
         endif
         if(gt.ne.0.d0) then
            AMP = AMP + mt2*gt*ALOSMh(mhiggs2,mt2)
     &                   !+ gt*ALOSMeffh(mhiggs2)
         endif
         if(gbp.ne.0.d0) then
            AMP = AMP + mbp2*gbp*ALOSMh(mhiggs2,mbp2)
         endif
         if(gtp.ne.0.d0) then
            AMP = AMP + mtp2*gtp*ALOSMh(mhiggs2,mtp2)
         endif
         if(gb1.ne.0.d0) then
            AMP = AMP + mbsb2*(gb1*ALOSUSYh(mhiggs2,msb12)
     &                       + gb2*ALOSUSYh(mhiggs2,msb22))
         endif
         if(gt1.ne.0.d0) then
            AMP = AMP + mt2*(gt1*ALOSUSYh(mhiggs2,mst12)
     &                     + gt2*ALOSUSYh(mhiggs2,mst22))
         endif
         if(gbp1.ne.0.d0) then
            AMP = AMP + mbp2*(gbp1*ALOSUSYh(mhiggs2,msbp12)
     &                      + gbp2*ALOSUSYh(mhiggs2,msbp22))
         endif
         if(gtp1.ne.0.d0) then
            AMP = AMP + mtp2*(gtp1*ALOSUSYh(mhiggs2,mstp12)
     &                      + gtp2*ALOSUSYh(mhiggs2,mstp22))
         endif

C-}}}

      endif

      sAMPLO=cdabs(AMP)**2

      end

C-}}}
C-{{{ squared amplidute for gg -> h

      function AMPLOpure(mhiggs2)
      implicit none
      double precision mhiggs2
      double complex AMP,AMPLOpure
c      double precision ALOSMeffh
      double complex ALOSMh,ALOSMA,ALOSUSYh

      include 'common-quark.f'
      include 'common-vars.f'

      AMP = 0.d0

      if (pseudo.eq.1) then

C-{{{ pseudoscalar higgs

         if(gc.ne.0.d0) then
            AMP = AMP + mc2*gc*ALOSMA(mhiggs2,mc2)
         endif
         if(gb.ne.0.d0) then
            AMP = AMP + mb2*gb*ALOSMA(mhiggs2,mb2)
         endif
         if(gt.ne.0.d0) then
            AMP = AMP + mt2*gt*ALOSMA(mhiggs2,mt2)
     &                   !+ gt*ALOSMeffA(mhiggs2)
         endif
         if(gbp.ne.0.d0) then
            AMP = AMP + mbp2*gbp*ALOSMA(mhiggs2,mbp2)
         endif
         if(gtp.ne.0.d0) then
            AMP = AMP + mtp2*gtp*ALOSMA(mhiggs2,mtp2)
         endif

C-}}}

      else

C-{{{ scalar higgs

         if(gc.ne.0.d0) then
            AMP = AMP + mc2*gc*ALOSMh(mhiggs2,mc2)
         endif
         if(gb.ne.0.d0) then
            AMP = AMP + mb2*gb*ALOSMh(mhiggs2,mb2)
         endif
         if(gt.ne.0.d0) then
            AMP = AMP + mt2*gt*ALOSMh(mhiggs2,mt2)
     &                   !+ gt*ALOSMeffh(mhiggs2)
         endif
         if(gbp.ne.0.d0) then
            AMP = AMP + mbp2*gbp*ALOSMh(mhiggs2,mbp2)
         endif
         if(gtp.ne.0.d0) then
            AMP = AMP + mtp2*gtp*ALOSMh(mhiggs2,mtp2)
         endif
         if(gb1.ne.0.d0) then
            AMP = AMP + mbsb2*(gb1*ALOSUSYh(mhiggs2,msb12)
     &                       + gb2*ALOSUSYh(mhiggs2,msb22))
         endif
         if(gt1.ne.0.d0) then
            AMP = AMP + mt2*(gt1*ALOSUSYh(mhiggs2,mst12)
     &                     + gt2*ALOSUSYh(mhiggs2,mst22))
         endif
         if(gbp1.ne.0.d0) then
            AMP = AMP + mbp2*(gbp1*ALOSUSYh(mhiggs2,msbp12)
     &                      + gbp2*ALOSUSYh(mhiggs2,msbp22))
         endif
         if(gtp1.ne.0.d0) then
            AMP = AMP + mtp2*(gtp1*ALOSUSYh(mhiggs2,mstp12)
     &                      + gtp2*ALOSUSYh(mhiggs2,mstp22))
         endif

C-}}}

      endif

      AMPLOpure=AMP

      end

C-}}}
C-{{{ amplidute for gg -> h (only b and sb)

      function AMPLOb(mhiggs2)
      implicit none
      double precision mhiggs2
      double complex AMPLOb,AMP
      double complex ALOSMh,ALOSMA,ALOSUSYh

      include 'common-quark.f'
      include 'common-vars.f'

      AMP = 0.d0

      if (pseudo.eq.1) then

         if(gb.ne.0.d0) then
            AMP = AMP + mb2*gb*ALOSMA(mhiggs2,mb2)
         endif

      else

         if(gb.ne.0.d0) then
            AMP = AMP + mb2*gb*ALOSMh(mhiggs2,mb2)
         endif
         if(gb1.ne.0.d0) then
            AMP = AMP + mbsb2*(gb1*ALOSUSYh(mhiggs2,msb12)
     &                       + gb2*ALOSUSYh(mhiggs2,msb22))
         endif

      endif

      AMPLOb = AMP

      end

C-}}}
C-{{{ amplidute for gg -> h (only c)

      function AMPLOc(mhiggs2)
      implicit none
      double precision mhiggs2
      double complex AMPLOc,AMP
      double complex ALOSMh,ALOSMA

      include 'common-quark.f'
      include 'common-vars.f'

      AMP = 0.d0

      if (pseudo.eq.1) then

         if(gc.ne.0.d0) then
            AMP = AMP + mc2*gc*ALOSMA(mhiggs2,mc2)
         endif

      else

         if(gc.ne.0.d0) then
            AMP = AMP + mc2*gc*ALOSMh(mhiggs2,mc2)
         endif

      endif

      AMPLOc = AMP

      end

C-}}}
C-{{{ amplidute for gg -> h eff

      function AMPLOeff(mhiggs2)
      implicit none
      double precision mhiggs2
      double complex AMPLOeff
      double precision ALOSUSYeffh,ALOSMeffh,ALOSMeffA
      double complex ALOSMh,ALOSMA

      include 'common-quark.f'
      include 'common-vars.f'

      if(pseudo.eq.1) then

         if(htl) then
            AMPLOeff = mc2 * gc * ALOSMA(mhiggs2,mc2)
     &            +    mb2 * gb * ALOSMA(mhiggs2,mb2)
     &            +          gt * ALOSMeffA()
     &            +    mbp2* gbp* ALOSMA(mhiggs2,mbp2)
     &            +    mtp2* gtp* ALOSMA(mhiggs2,mtp2)
         else
            AMPLOeff = mc2 * gc * ALOSMA(mhiggs2,mc2)
     &            +    mb2 * gb * ALOSMA(mhiggs2,mb2)
     &            +    mt2 * gt * ALOSMA(mhiggs2,mt2)
     &            +    mbp2* gbp* ALOSMA(mhiggs2,mbp2)
     &            +    mtp2* gtp* ALOSMA(mhiggs2,mtp2)
         endif

      else

         if(htl) then
            AMPLOeff = mc2 * gc * ALOSMh(mhiggs2,mc2)
     &            +    mb2 * gb * ALOSMh(mhiggs2,mb2)
     &            +          gt * ALOSMeffh()
     &            +    mbp2* gbp* ALOSMh(mhiggs2,mbp2)
     &            +    mtp2* gtp* ALOSMh(mhiggs2,mtp2)
     &            + gte + gtpe + gbpe
            if(gb1.ne.0.d0) then
               AMPLOeff = AMPLOeff +
     &              mbsb2 * (gb1/msb12+gb2/msb22)*ALOSUSYeffh()
            endif
         else
            AMPLOeff = mc2 * gc * ALOSMh(mhiggs2,mc2)
     &            +    mb2 * gb * ALOSMh(mhiggs2,mb2)
     &            +    mt2 * gt * ALOSMh(mhiggs2,mt2)
     &            +    mbp2* gbp* ALOSMh(mhiggs2,mbp2)
     &            +    mtp2* gtp* ALOSMh(mhiggs2,mtp2)
     &            + gte + gtpe + gbpe
            if(gb1.ne.0.d0) then
               AMPLOeff = AMPLOeff +
     &              mbsb2 *(gb1/msb12+gb2/msb22)*ALOSUSYeffh()
            endif
         endif
      endif
      end
C-}}}
