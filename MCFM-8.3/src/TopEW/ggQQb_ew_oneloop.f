      subroutine ggQQb_ew_oneloop(vew,s,beta,z)
C --- the analytical result is taken from [arXiv:hep-ph/0610335] 
C --- by J. H. Kuehn, A. Scharf and P. Uwer
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'scale.f'
      include 'anomcoup.f'
      real(dp):: vew,s,t,beta,z,vrts(5),slf(5),vrt(5),bx(5),alpha,
     .     genfacs,genfacself,genfacvert,genfacbox,genfacschi,sigma0,
     .     chifacs,chifac,gvt,gat,xI1,xI2,xI3,xI4,sw2,cw2,
     .     mw,mz,mh,born,aemmz,trih,trizx
      real(dp) :: db0
      integer ep
      common/em/aemmz

      ep = 0
c      musq = (2._dp*mt)**2
c      alpha = 1._dp/126.3_dp
c--- alpha -> standard MCFM value
      alpha=aemmz
      sw2 = xw


      mw = wmass
      mz = zmass
      mh = hmass

c      sw2 = 1._dp - mw**2/mz**2
      cw2 = 1._dp-sw2

      gvt = (0.5_dp-4._dp*sw2/3._dp)/2._dp/sqrt(sw2*cw2)
      gat = 0.5_dp/2._dp/sqrt(sw2*cw2)
      gw = 0.5_dp/sqrt(2._dp*sw2)

      genfacs = Nc*mt**2/s**2*z**2/(1._dp-beta**2*z**2)
      genfacself = (2._dp - Nc**2*(1._dp-beta*z))/Nc/(1._dp-beta**2*z**2)
      genfacvert = (2._dp - Nc**2*(1._dp-beta*z))/Nc/(1._dp+beta*z)**2
      genfacbox = (Nc**2*(1._dp-beta*z)-2._dp)/Nc/(1._dp-beta**2*z**2)
      genfacschi = Nc*mt**4/s**2*z**2/(1._dp-beta**2*z**2)
      chifacs = gat**2
      chifac = gat**2/mz**2
c --- sigma0 is jacobian additional to matrix element square
c --- sigma0 = gsq**2*beta/64._dp/pi/s/(Nc**2-1._dp)
      sigma0 = 1._dp
c      sigma0 = 1._dp/16._dp

      t = -s*(1._dp+beta*z)/2._dp+mt**2

C --- leading order differential cross section, use to normalize the result
C --- in the ref. (II.1)
      born = sigma0*(Nc**2*(1._dp+beta**2*z**2)-2._dp)/
     .     Nc/(1._dp-beta**2*z**2)**2*(1._dp-beta**4*z**4
     .     +2._dp*beta**2*(1._dp-beta**2)*(1._dp-z**2))
C --- the differential cross section is normalized when 
C --- sigma0 = (pi^2*alpha_s^2)/2 = gsq^2/32


C --- add triangle t,b contribution (II.16-17)
      trih = 
     .     + alpha/pi*sigma0*mt**2/mw**2/sw2*tevscale*
     .     beta**2/(1._dp - beta**2*z**2)/(s - mh**2)*
     .     ( + tevscale*
     .     ( + mt**2*(s - 4._dp*mt**2)*
     .     xI3(0._dp,0._dp,s,mt**2,mt**2,mt**2,musq,ep) - 2._dp*mt**2)
     .     + mb**2*(s - 4._dp*mb**2)*
     .     xI3(0._dp,0._dp,s,mb**2,mb**2,mb**2,musq,ep) - 2._dp*mb**2)

      
      trizx = 
     .     + 16._dp*alpha/pi*sigma0*gat*mt**2/mz**2
     .     /(1._dp-beta**2*z**2)*gat*(
     .     + mt**2*xI3(0._dp,0._dp,s,mt**2,mt**2,mt**2,musq,ep)
     .     - mb**2*xI3(0._dp,0._dp,s,mb**2,mb**2,mb**2,musq,ep) !used gab = - gat
     .     )


C --- In the ref. (II.18-22), and Appendix B (B.1-15).
      vrts(1) = 
     .     (-2._dp*alpha*genfacs*sigma0*
     .     (beta**2*s*(2._dp*(gat**2 + gvt**2)*mz**2 
     .     + (1._dp - beta**2)*(-3._dp*gat**2 + gvt**2)*s)*
     .     db0(mt**2,mt**2,mz**2) + 2._dp*(2._dp*(gat**2 + gvt**2)*mz**2 
     .     - beta**2*(-3._dp*gat**2 + gvt**2)*s)*( - xI2(mt**2,
     .     mt**2,mz**2,musq,ep) + xI2(s,mt**2,mt**2,musq,ep)) + (4._dp*
     .     (gat**2 + gvt**2)*mz**4 + 8._dp*beta**2*gat**2*mz**2*s
     .     - beta**2*(gat**2 + gvt**2 + beta**2*(-3._dp*gat**2 + gvt**2)
     .     )*s**2)*xI3(s,mt**2,mt**2,mt**2,mt**2,mz**2,musq,ep)))/pi  



      vrts(2) = 
     .     (-8._dp*alpha*genfacs*gw**2*sigma0*
     .     (-0.25_dp*beta**2*s*(4._dp*(mb**2 - mw**2) + (1._dp 
     .     - beta**2)*s)*db0(mt**2,mb**2,mw**2) + 0.5_dp*(4._dp*
     .     (-mb**2 + mw**2) + (1._dp + beta**2)*s)*
     .     (-xI2(mt**2,mb**2,mw**2,musq,ep) 
     .     + xI2(s,mb**2,mb**2,musq,ep)) 
     .     + 0.125_dp*(-4._dp*beta**2*s**2 + (4._dp*(-mb**2 + mw**2) 
     .     + (1._dp + beta**2)*s)**2)*
     .     xI3(s,mt**2,mt**2,mb**2,mb**2,mw**2,musq,ep)))/pi



      vrts(3) = 
     .     (-8._dp*alpha*chifacs*genfacschi*sigma0*
     .     (beta**2*s*db0(mt**2,mt**2,mz**2) + 2._dp*
     .     (-xI2(mt**2,mt**2,mz**2,musq,ep) 
     .     + xI2(s,mt**2,mt**2,musq,ep)) 
     .     + (2._dp*mz**2 + beta**2*s)*
     .     xI3(s,mt**2,mt**2,mt**2,mt**2,mz**2,musq,ep)))/pi



      vrts(4) = 
     .     (-8._dp*alpha*genfacs*gw**2*sigma0*
     .     (-0.03125_dp*(16._dp*beta**2*mb**4*s - 4._dp*beta**2*(1._dp 
     .     - beta**2)*mw**2*s**2 + beta**2*(1._dp - beta**2 )**2*s**3
     .     - 8._dp*mb**2*(2._dp*beta**2*mw**2*s + beta**2*(1._dp 
     .     - beta**2)*s**2))*db0(mt**2,mb**2,mw**2) + 0.0625_dp*
     .     (16._dp*mb**4 - 4._dp*(1._dp - beta**2)*mw**2*s 
     .     - (1._dp - beta**4)*s**2 + mb**2*(-16._dp*mw**2 
     .     + 8._dp*beta**2*s))*(xI2(mt**2,mb**2,mw**2,musq,ep) 
     .     - xI2(s,mb**2,mb**2,musq,ep)) + 0.015625_dp*(64._dp*mb**6 
     .     + 16._dp*(1._dp - beta**2)*mw**4*s + 8._dp*(1._dp 
     .     - beta**4)*mw**2*s**2 + (1._dp - beta**2)**3*s**3 
     .     - 16._dp*mb**4*(8._dp*mw**2 + (1._dp - beta**2)*s) 
     .     + 4._dp*mb**2*(16._dp*mw**4 - (1._dp - beta**2)**2*
     .     s**2))*xI3(s,mt**2,mt**2,mb**2,mb**2,mw**2,musq,ep)))
     .     /(mw**2*pi)



      vrts(5) = 
     .     (4._dp*alpha*genfacs*gw**2*mt**2*sigma0*
     .     (beta**2*s*(-mh**2 + (1._dp - beta**2)*s)*db0(mt**2
     .     ,mt**2,mh**2) + 2._dp*(mh**2 + beta**2*s)*(xI2(mt**2,mt**2,
     .     mh**2,musq,ep) - xI2(s,mt**2,mt**2,musq,ep)) + 
     .     (-2._dp*mh**4 - 3._dp*beta**2*mh**2*s + beta**2*(1._dp 
     .     - beta**2)*s**2)*xI3(s,mt**2,mt**2,mt**2,mt**2,mh**2
     .     ,musq,ep)))/(mw**2*pi)



      slf(1) = 
     .     (0.125_dp*alpha*genfacself*sigma0*
     .     ((2._dp*(2._dp*(gat**2 + gvt**2)*mz**2 + 
     .     (1._dp - beta**2)*(-3._dp*gat**2 + gvt**2)*s)*
     .     (1._dp - beta**4*z**4 + (1._dp - beta**2)*
     .     (2._dp*beta**2 - beta*z - 3._dp*beta**2*z**2))*
     .     db0(mt**2,mt**2,mz**2))/(1._dp + beta*z) + 
     .     (16._dp*(gat**2 + gvt**2)*(1._dp - beta**4*z**4 + 
     .     beta**2*(1._dp - beta**2)*(1._dp - 3._dp*z**2))*
     .     (-xI1(mt**2,musq,ep) + xI1(mz**2,musq,ep)))/ 
     .     ((1._dp - beta**2)*s*(1._dp + beta**2 + 2._dp*beta*z))  
     .     - (4._dp*((1._dp - beta**2)*(-3._dp*gat**2 + gvt**2)*
     .     (1._dp + beta**2 - 2._dp*beta**4 -
     .     beta**2*(2._dp - 3._dp*beta**2)*z**2 - beta**4*z**4) 
     .     + (2._dp*(gat**2 + gvt**2)*mz**2*(2._dp + 2._dp*beta**2 - 
     .     5._dp*beta**4 + 2._dp*beta**6 + beta**3*(3._dp - 2._dp*beta**2)*z 
     .     - 3._dp*beta**2*(2._dp - 3._dp*beta**2 + beta**4)*z**2 - 
     .     3._dp*beta**3*(1._dp - beta**2)*z**3 - beta**4*(2._dp 
     .     - beta**2)*z**4 - beta**5*z**5))/((1._dp 
     .     - beta**2)*s))*xI2(mt**2,mt**2,mz**2,musq,ep))/(1._dp 
     .     + beta*z)**2 + (4._dp*((2._dp*beta**2*(gat**2 + gvt**2)*mz**2*
     .     (1._dp - z**2)*(2._dp + beta**2 - 2._dp*beta**4 + 
     .     (3._dp*beta - 2._dp*beta**3)*z + beta**4*z**2 + beta**3*z**3))
     .     /s + gvt**2*(2._dp + 2._dp*beta**2 - 4._dp*beta**4 - 
     .     beta**6 + 2._dp*beta**8 + beta*(4._dp + 2._dp*beta**2 - 
     .     8._dp*beta**4 + 4._dp*beta**6)*z + beta**2*(-4._dp + 
     .     7._dp*beta**2 + beta**4 - 3._dp*beta**6)*z**2 - beta**3*
     .     (10._dp - 16._dp*beta**2 + 6._dp*beta**4)*z**3 + beta**4*(-5._dp 
     .     + 3._dp*beta**2 + beta**4)*z**4 - beta**5*(4._dp 
     .     - 2._dp*beta**2)*z**5 - beta**6*z**6) - gat**2*
     .     (2._dp + 2._dp*beta**2 - 8._dp*beta**4 - 3._dp*beta**6 
     .     + 6._dp*beta**8 + 2._dp*beta*(2._dp - beta**2 
     .     - 8._dp*beta**4 + 6._dp*beta**6)*z - beta**2*(4._dp 
     .     - 5._dp*beta**2 - 7._dp*beta**4 + 9._dp*beta**6)*z**2 
     .     - 6._dp*beta**3*(1._dp - 4._dp*beta**2 + 3._dp*beta**4)*z**3 
     .     + beta**4*(1._dp - 3._dp*beta**2 + 3._dp*beta**4)*z**4 
     .     - beta**5*(4._dp - 6._dp*beta**2)*z**5 + beta**6*z**6))*
     .     xI2(t,mt**2,mz**2,musq,ep))/((1._dp + beta*z)**2*(1._dp 
     .     + beta**2 + 2._dp*beta*z))))/pi



      slf(2) = 
     .     (alpha*genfacself*gw**2*sigma0*
     .     ((0.25_dp*(-4._dp*mb**2 + 4._dp*mw**2 - s + beta**2*s)*
     .     (1._dp - beta**4*z**4 + beta*(1._dp - beta**2)*
     .     (2._dp*beta - z - 3._dp*beta*z**2))*
     .     db0(mt**2,mb**2,mw**2))/(1._dp + beta*z) + (4._dp*(1._dp - 
     .     beta**4*z**4 + beta**2*(1._dp - beta**2)*(1._dp - 
     .     3._dp*z**2))*(-xI1(mb**2,musq,ep) + xI1(mw**2,musq,ep)))
     .     /((1._dp - beta**2)*s*(1._dp + beta**2 + 2._dp*beta*z)) + 
     .     (0.5_dp*((-4._dp*(mb**2 - mw**2)*(-2._dp - 2._dp*beta**2 + 
     .     5._dp*beta**4 - 2._dp*beta**6 - beta**3*(3._dp - 
     .     2._dp*beta**2)*z + 3._dp*beta**2*(2._dp - 3._dp*beta**2 + 
     .     beta**4)*z**2 + 3._dp*beta**3*(1._dp - beta**2)*z**3 + 
     .     beta**4*(2._dp - beta**2)*z**4 + beta**5*z**5))/(1._dp - 
     .     beta**2) - beta**2*s*(1._dp - z**2)*(2._dp + 
     .     beta**2 - 2._dp*beta**4 + beta*(3._dp - 2._dp*beta**2)*z + 
     .     beta**3*z**2*(beta + z)))*xI2(mt**2,mb**2,mw**2,musq,ep))/
     .     (s*(1._dp + beta*z)**2) + (0.5_dp*beta**2*(1._dp - z**2)*
     .     (2._dp + beta**2 - 2._dp*beta**4 + beta*(3._dp - 2._dp*beta**2)*z 
     .     + beta**3*z**2*(beta + z))*(4._dp*(-mb**2 + mw**2) + 
     .     s*(1._dp + beta**2 + 2._dp*beta*z))*xI2(t,mb**2,mw**2,musq,ep))
     .     /(s*(1._dp + beta*z)**2*(1._dp + beta**2 + 2._dp*beta*z))))/pi



      slf(3) = 
     .     (2._dp*alpha*chifac*genfacself*sigma0*
     .     ((0.125_dp*(1._dp - beta**2)*mz**2*s*(1._dp 
     .     - beta**4*z**4 + beta*(1._dp - beta**2)*
     .     (2._dp*beta - z - 3._dp*beta*z**2))*db0(mt**2,mt**2,
     .     mz**2))/(1._dp + beta*z) + (0.5_dp*(1._dp - beta**4*z**4 
     .     + (1._dp - beta**2)*(beta**2 - 3._dp*beta**2*z**2))*
     .     (-xI1(mt**2,musq,ep) + xI1(mz**2,musq,ep)))/(1._dp 
     .     + beta**2 + 2._dp*beta*z) - (0.25_dp*mz**2*(2._dp + 2._dp*beta**2 
     .     - 5._dp*beta**4 + 2._dp*beta**6 + beta**3*(3._dp 
     .     - 2._dp*beta**2)*z - 3._dp*beta**2*(2._dp - 3._dp*beta**2 
     .     + beta**4)*z**2 - 3._dp*beta**3*(1._dp - beta**2)*z**3 
     .     - beta**4*(2._dp - beta**2)*z**4 
     .     - beta**5*z**5)*xI2(mt**2,mt**2,mz**2,musq,ep))/
     .     (1._dp + beta*z)**2 + (0.125_dp*(1._dp - beta**2)*
     .     (s*(1._dp + beta**2 - beta**4 - 3._dp*beta**2*(1._dp 
     .     - beta**2)*z**2 - beta**4*z**4) + 
     .     (2._dp*beta**2*mz**2*(1._dp - z**2)*(2._dp + beta**2 
     .     - 2._dp*beta**4 + beta*(3._dp - 2._dp*beta**2)*z 
     .     + beta**3*z**2*(beta + z)))/(1._dp + beta*z)**2)*
     .     xI2(t,mt**2,mz**2,musq,ep))/
     .     (1._dp + beta**2 + 2._dp*beta*z)))/pi



      slf(4) = 
     .     (alpha*genfacself*gw**2*sigma0*
     .     ((-0.03125_dp*(16._dp*mb**2*(mb**2 - mw**2) 
     .     - 4._dp*(1._dp - beta**2)*(2._dp*mb**2 + mw**2)*s 
     .     + (1._dp - beta**2)**2*s**2)*(1._dp - beta**4*z**4 
     .     + beta*(1._dp - beta**2)*(2._dp*beta - z 
     .     - 3._dp*beta*z**2))*db0(mt**2,mb**2,mw**2))/(1._dp + beta*z) 
     .     + (0.5_dp*(4._dp*mb**2 + (1._dp - beta**2)*s)*(1._dp 
     .     - beta**4*z**4 + beta**2*(1._dp - beta**2)*(1._dp 
     .     - 3._dp*z**2))*(-xI1(mb**2,musq,ep) 
     .     + xI1(mw**2,musq,ep)))/((1._dp - beta**2)*s*
     .     (1._dp + beta**2 + 2._dp*beta*z)) - (0.0625_dp*(8._dp*(1._dp 
     .     - beta**2)**2*mb**2*s*(1._dp + beta**2 - 2._dp*beta**4 
     .     - beta**2*(2._dp - 3._dp*beta**2)*z**2 
     .     - beta**4*z**4) - 4._dp*(4._dp*mb**2*(-mb**2 + mw**2) 
     .     + (1._dp - beta**2)*mw**2*s)*(-2._dp - 2._dp*beta**2 
     .     + 5._dp*beta**4 - 2._dp*beta**6 - beta**3*(3._dp 
     .     - 2._dp*beta**2)*z + 3._dp*beta**2*(2._dp - 3._dp*beta**2 
     .     + beta**4)*z**2 + 3._dp*beta**3*(1._dp - beta**2)*z**3 
     .     + beta**4*(2._dp - beta**2)*z**4 + beta**5*z**5) 
     .     + beta**2*(1._dp - beta**2)**2*s**2*(1._dp - z**2)*
     .     (2._dp + beta**2 - 2._dp*beta**4 + beta*(3._dp - 2._dp*beta**2)*z 
     .     + beta**3*z**2*(beta + z)))*xI2(mt**2,mb**2,mw**2,musq,ep))/
     .     ((1._dp - beta**2)*s*(1._dp + beta*z)**2) + 
     .     (0.0625_dp*(-8._dp*mb**2*s*(-2._dp - 2._dp*beta**2 + 4._dp*beta**4 
     .     + beta**6 - 2._dp*beta**8 - 2._dp*beta*(2._dp + beta**2 
     .     - 4._dp*beta**4 + 2._dp*beta**6)*z + beta**2*(4._dp 
     .     - 7._dp*beta**2 - beta**4 + 3._dp*beta**6)*z**2 
     .     + 2._dp*beta**3*(5._dp - 8._dp*beta**2 + 3._dp*beta**4)*z**3 
     .     + beta**4*(5._dp - 3._dp*beta**2 - beta**4)*z**4 
     .     + 2._dp*beta**5*(2._dp - beta**2)*z**5 + beta**6*z**6) 
     .     - beta**2*(1._dp - z**2)*(2._dp + beta**2 
     .     - 2._dp*beta**4 + beta*(3._dp - 2._dp*beta**2)*z 
     .     + beta**3*z**2*(beta + z))*(16._dp*mb**2*(mb**2 - mw**2) 
     .     - 4._dp*(1._dp - beta**2)*mw**2*s - s**2*(1._dp 
     .     - beta**4 + 2._dp*beta*(1._dp - beta**2)*z)))*
     .     xI2(t,mb**2,mw**2,musq,ep))/(s*(1._dp + beta*z)**2*
     .     (1._dp + beta**2 + 2._dp*beta*z))))/(mw**2*pi)



      slf(5) = 
     .     (alpha*genfacself*gw**2*sigma0*
     .     ((-0.125_dp*(1._dp - beta**2)*s*(-mh**2 
     .     + (1._dp - beta**2)*s)*(1._dp - beta**4*z**4 
     .     + beta*(1._dp - beta**2)*(2._dp*beta - z 
     .     - 3._dp*beta*z**2))*db0(mt**2,mt**2,mh**2))/(1._dp + beta*z) 
     .     + (0.5_dp*(1._dp - beta**4*z**4 + beta**2*(1._dp 
     .     - beta**2)*(1._dp - 3._dp*z**2))*(xI1(mh**2,musq,ep) 
     .     - xI1(mt**2,musq,ep)))/(1._dp + beta**2 + 2._dp*beta*z) 
     .     + (0.25_dp*((1._dp - beta**2)**2*s*(1._dp + beta**2 
     .     - 2._dp*beta**4 - beta**2*(2._dp - 3._dp*beta**2)*z**2 
     .     - beta**4*z**4) + mh**2*(-2._dp - 2._dp*beta**2 
     .     + 5._dp*beta**4 - 2._dp*beta**6 - beta**3*(3._dp 
     .     - 2._dp*beta**2)*z + 3._dp*beta**2*(2._dp - 3._dp*beta**2 
     .     + beta**4)*z**2 + 3._dp*beta**3*(1._dp - beta**2)*z**3 
     .     + beta**4*(2._dp - beta**2)*z**4 + beta**5*z**5))*
     .     xI2(mt**2,mt**2,mh**2,musq,ep))/(1._dp + beta*z)**2 - 
     .     (0.125_dp*(1._dp - beta**2)*(s*(1._dp + beta**2 
     .     - 5._dp*beta**4 - 2._dp*beta**6 + 4._dp*beta**8 
     .     + 2._dp*beta*(1._dp - beta**2 - 5._dp*beta**4 
     .     + 4._dp*beta**6)*z - beta**2*(2._dp - 2._dp*beta**2 
     .     - 5._dp*beta**4 + 6._dp*beta**6)*z**2 - 2._dp*beta**3*(1._dp 
     .     - 7._dp*beta**2 + 6._dp*beta**4)*z**3 + beta**4*(2._dp 
     .     - 3._dp*beta**2 + 2._dp*beta**4)*z**4 - 2._dp*beta**5*(1._dp 
     .     - 2._dp*beta**2)*z**5 + beta**6*z**6) - 2._dp*beta**2*mh**2*
     .     (1._dp - z**2)*(2._dp + beta**2 - 2._dp*beta**4 + 
     .     beta*(3._dp - 2._dp*beta**2)*z + beta**3*z**2*(beta + z)))*
     .     xI2(t,mt**2,mh**2,musq,ep))/((1._dp + beta*z)**2*
     .     (1._dp + beta**2 + 2._dp*beta*z))))/(mw**2*pi)



      vrt(1) = 
     .     (0.125_dp*alpha*genfacvert*sigma0*
     .     (4._dp*(1._dp - beta**2)*(gat**2 + gvt**2) 
     .     - (4._dp*(2._dp*(gat**2 + gvt**2)*mz**2 + 
     .     (1._dp - beta**2)*(-3._dp*gat**2 + gvt**2)*s)*
     .     (1._dp - beta**4*z**4 + beta*(1._dp - beta**2)*
     .     (2._dp*beta - z - 3._dp*beta*z**2))*db0(mt**2,mt**2,mz**2))/
     .     (1._dp - beta*z) + (16._dp*(gat**2 + gvt**2)*(1._dp 
     .     + 2._dp*beta**2 - beta**4 + beta*(1._dp + 4._dp*beta**2 
     .     - 3._dp*beta**4)*z - 2._dp*beta**4*z**4*(1._dp + beta*z) 
     .     - 2._dp*(1._dp - beta**2)*(2._dp*beta**2*z**2 
     .     + 3._dp*beta**3*z**3))*(xI1(mt**2,musq,ep) 
     .     - xI1(mz**2,musq,ep)))/((1._dp - beta**2)*s*
     .     (1._dp - beta*z)*(1._dp + beta**2 + 2._dp*beta*z)) 
     .     + (4._dp*((2._dp*(gat**2 + gvt**2)*mz**2*(1._dp + 8._dp*beta**2 
     .     - 11._dp*beta**4 + 4._dp*beta**6 + beta*(1._dp + 4._dp*beta**2 
     .     - 3._dp*beta**4)*z - 2._dp*beta**2*(5._dp - 8._dp*beta**2 
     .     + 3._dp*beta**4)*z**2 - 6._dp*beta**3*(1._dp - beta**2)*
     .     z**3 - 2._dp*beta**4*(2._dp - beta**2)*z**4 
     .     - 2._dp*beta**5*z**5))/s + (1._dp - beta**2)**2*
     .     (gvt**2*(1._dp + beta**2 - 4._dp*beta**4 
     .     - beta*(1._dp - beta**2)*z - 2._dp*beta**2*(1._dp 
     .     - 3._dp*beta**2)*z**2 - 2._dp*beta**4*z**4) + gat**2*(1._dp 
     .     - 7._dp*beta**2 + 12._dp*beta**4 - beta*(1._dp 
     .     - beta**2)*z + 6._dp*beta**2*(1._dp - 3._dp*beta**2)*z**2 
     .     + 6._dp*beta**4*z**4)))*xI2(mt**2,mt**2,mz**2,musq,ep))/
     .     ((1._dp - beta**2)*(1._dp - beta**2*z**2)) + 
     .     (4._dp*((2._dp*(gat**2 + gvt**2)*mz**2*(1._dp - 4._dp*beta**2 
     .     - 3._dp*beta**4 + 4._dp*beta**6 + beta*(1._dp - 8._dp*beta**2 
     .     + 5._dp*beta**4)*z + 2._dp*beta**2*(1._dp + 2._dp*beta**2 
     .     - 3._dp*beta**4)*z**2 + 6._dp*beta**3*(1._dp - beta**2)*
     .     z**3 + 2._dp*beta**6*z**4 + 2._dp*beta**5*z**5))/s 
     .     - gvt**2*(2._dp + 3._dp*beta**2 - 6._dp*beta**4 
     .     - beta**6 + 4._dp*beta**8 + beta*(3._dp + 7._dp*beta**2 
     .     - 13._dp*beta**4 + 7._dp*beta**6)*z - beta**2*(7._dp 
     .     - 18._dp*beta**2 + 3._dp*beta**4 + 6._dp*beta**6)*z**2 
     .     - beta**3*(14._dp - 26._dp*beta**2 + 12._dp*beta**4)*z**3 
     .     - beta**4*(10._dp - 6._dp*beta**2 - 2._dp*beta**4)*z**4 
     .     - beta**5*(8._dp - 4._dp*beta**2)*z**5 
     .     - 2._dp*beta**6*z**6) + gat**2*(-2._dp + 5._dp*beta**2 
     .     - 10._dp*beta**4 - 7._dp*beta**6 + 12._dp*beta**8 
     .     - beta*(3._dp - 9._dp*beta**2 + 35._dp*beta**4 
     .     - 25._dp*beta**6)*z - beta**2*(1._dp - 6._dp*beta**2 
     .     - 11._dp*beta**4 + 18._dp*beta**6)*z**2 - 2._dp*beta**3*(1._dp 
     .     - 19._dp*beta**2 + 18._dp*beta**4)*z**3 + 2._dp*beta**4*(1._dp 
     .     - 3._dp*beta**2 + 3._dp*beta**4)*z**4 - 4._dp*beta**5*(2._dp 
     .     - 3._dp*beta**2)*z**5 + 2._dp*beta**6*z**6))*
     .     xI2(t,mt**2,mz**2,musq,ep))/((1._dp + beta**2 + 2._dp*beta*z)*
     .     (1._dp - beta**2*z**2)) + 2._dp*(1._dp - beta**2)**2*
     .     (gat**2 + gvt**2)*s*
     .     xI3(0._dp,mt**2,t,mt**2,mt**2,mz**2,musq,ep)))/pi



      vrt(2) = 
     .     (alpha*genfacvert*gw**2*sigma0*
     .     (1._dp - beta**2 + (0.5_dp*(4._dp*(mb**2 - mw**2) + 
     .     (1._dp - beta**2)*s)*(1._dp - beta**4*z**4 + 
     .     beta*(1._dp - beta**2)*(2._dp*beta - z 
     .     - 3._dp*beta*z**2))*db0(mt**2,mb**2,mw**2))/
     .     (1._dp - beta*z) + (4._dp*(1._dp + 2._dp*beta**2 
     .     - beta**4 + beta*(1._dp + 4._dp*beta**2 - 3._dp*beta**4)*z 
     .     - 2._dp*beta**4*z**4*(1._dp + beta*z) 
     .     - 2._dp*(1._dp - beta**2)*(2._dp*beta**2*z**2 
     .     + 3._dp*beta**3*z**3))*(xI1(mb**2,musq,ep) 
     .     - xI1(mw**2,musq,ep)))/((1._dp - beta**2)*s*
     .     (1._dp - beta*z)*(1._dp + beta**2 + 2._dp*beta*z)) + 
     .     (0.5_dp*((1._dp - beta**2)*s*(3._dp + 3._dp*beta**4 
     .     - 4._dp*beta**6 - beta*(1._dp - 8._dp*beta**2 
     .     + 5._dp*beta**4)*z - 6._dp*beta**2*(1._dp - beta**4)*z**2 
     .     - 6._dp*beta**3*(1._dp - beta**2)*z**3 - 2._dp*beta**6*z**4 
     .     - 2._dp*beta**5*z**5) + 4._dp*(-mb**2 + mw**2)*(1._dp 
     .     + 8._dp*beta**2 - 11._dp*beta**4 + 4._dp*beta**6 + beta*(1._dp 
     .     + 4._dp*beta**2 - 3._dp*beta**4)*z - 2._dp*beta**2*(5._dp 
     .     - 8._dp*beta**2 + 3._dp*beta**4)*z**2 - 6._dp*beta**3*(1._dp 
     .     - beta**2)*z**3 - 2._dp*beta**4*(2._dp 
     .     - beta**2)*z**4 - 2._dp*beta**5*z**5))*
     .     xI2(mt**2,mb**2,mw**2,musq,ep))/((1._dp - beta**2)*s*
     .     (1._dp - beta**2*z**2)) + (0.5_dp*(-s*(1._dp 
     .     + beta**2 + 2._dp*beta*z)*(3._dp + 3._dp*beta**4 - 4._dp*beta**6 
     .     - beta*(1._dp - 8._dp*beta**2 + 5._dp*beta**4)*z 
     .     - 6._dp*beta**2*(1._dp - beta**4)*z**2 - 6._dp*beta**3*
     .     (1._dp - beta**2)*z**3 - 2._dp*beta**6*z**4 
     .     - 2._dp*beta**5*z**5) + 4._dp*(-mb**2 + mw**2)*
     .     (1._dp - 4._dp*beta**2 - 3._dp*beta**4 + 4._dp*beta**6 
     .     + beta*(1._dp - 8._dp*beta**2 + 5._dp*beta**4)*z 
     .     + 2._dp*beta**2*(1._dp + 2._dp*beta**2 - 3._dp*beta**4)*z**2 
     .     + 6._dp*beta**3*(1._dp - beta**2)*z**3 
     .     + 2._dp*beta**6*z**4 + 2._dp*beta**5*z**5))*
     .     xI2(t,mb**2,mw**2,musq,ep))/(s*(1._dp + beta**2 
     .     + 2._dp*beta*z)*(1._dp - beta**2*z**2)) 
     .     + 2._dp*(1._dp - beta**2)*mb**2*
     .     xI3(0._dp,mt**2,t,mb**2,mb**2,mw**2,musq,ep)))/pi




      vrt(3) = 
     .     (2._dp*alpha*chifac*genfacvert*sigma0*
     .     (0.125_dp*(1._dp - beta**2)**2*s - (0.25_dp*(1._dp 
     .     - beta**2)*mz**2*s*(1._dp - beta**4*z**4 
     .     + beta*(1._dp - beta**2)*(2._dp*beta - z 
     .     - 3._dp*beta*z**2))*db0(mt**2,mt**2,mz**2))/
     .     (1._dp - beta*z) + (0.5_dp*(1._dp + 2._dp*beta**2 
     .     - beta**4 + beta*(1._dp + 4._dp*beta**2 - 3._dp*beta**4)*z 
     .     - 2._dp*beta**4*z**4*(1._dp + beta*z) - 2._dp*(1._dp 
     .     - beta**2)*(2._dp*beta**2*z**2 + 3._dp*beta**3*z**3))*
     .     (xI1(mt**2,musq,ep) - xI1(mz**2,musq,ep)))/
     .     ((1._dp - beta*z)*(1._dp + beta**2 + 2._dp*beta*z)) 
     .     + (0.125_dp*(-(1._dp - beta**2)**2*s*(1._dp 
     .     - beta*z)*(1._dp + beta**2 + 2._dp*beta*z) 
     .     + 2._dp*mz**2*(1._dp + 8._dp*beta**2 - 11._dp*beta**4 
     .     + 4._dp*beta**6 + beta*(1._dp + 4._dp*beta**2 - 3._dp*beta**4)*z 
     .     - 2._dp*beta**2*(5._dp - 8._dp*beta**2 + 3._dp*beta**4)*z**2 
     .     - 6._dp*beta**3*(1._dp - beta**2)*z**3 - 2._dp*beta**4*
     .     (2._dp - beta**2)*z**4 - 2._dp*beta**5*z**5))*
     .     xI2(mt**2,mt**2,mz**2,musq,ep))/(1._dp - beta**2*z**2) 
     .     + (0.125_dp*(1._dp - beta**2)*(2._dp*mz**2*(1._dp 
     .     - 4._dp*beta**2 - 3._dp*beta**4 + 4._dp*beta**6 + beta*(1._dp 
     .     - 8._dp*beta**2 + 5._dp*beta**4)*z + 2._dp*beta**2*(1._dp 
     .     + 2._dp*beta**2 - 3._dp*beta**4)*z**2 + 6._dp*beta**3*(1._dp 
     .     - beta**2)*z**3 + 2._dp*beta**6*z**4 + 2._dp*beta**5*z**5) 
     .     + beta*s*(-beta*(1._dp + beta**4) + (1._dp - 7._dp*beta**2 
     .     + beta**4 + beta**6)*z + beta*(3._dp - 12._dp*beta**2 
     .     + 7._dp*beta**4)*z**2 + 6._dp*beta**2*(1._dp - beta**2)*
     .     z**3 + 2._dp*beta**3*(4._dp - 3._dp*beta**2)*z**4 
     .     + 4._dp*beta**4*z**5 + 2._dp*beta**5*z**6))*
     .     xI2(t,mt**2,mz**2,musq,ep))/((1._dp + beta**2 + 2._dp*beta*z)*
     .     (1._dp - beta**2*z**2)) 
     .     - 0.0625_dp*(1._dp - beta**2)**2*s**2*
     .     (1._dp + beta**2 + 2._dp*beta*z)*
     .     xI3(0._dp,mt**2,t,mt**2,mt**2,mz**2,musq,ep)))/pi



      vrt(4) = 
     .     (alpha*genfacvert*gw**2*sigma0*
     .     (0.125_dp*(1._dp - beta**2)*(4._dp*mb**2 + (1._dp 
     .     - beta**2)*s) + (0.0625_dp*(16._dp*mb**2*(mb**2 
     .     - mw**2) - 4._dp*(1._dp - beta**2)*(2._dp*mb**2 
     .     + mw**2)*s + (1._dp - beta**2)**2*s**2)*(1._dp 
     .     - beta**4*z**4 + beta*(1._dp - beta**2)*
     .     (2._dp*beta - z - 3._dp*beta*z**2))*db0(mt**2,mb**2,
     .     mw**2))/(1._dp - beta*z) + (0.5_dp*(4._dp*mb**2 + 
     .     (1._dp - beta**2)*s)*(1._dp + 2._dp*beta**2 - beta**4 
     .     + beta*(1._dp + 4._dp*beta**2 - 3._dp*beta**4)*z - 2._dp*beta**4
     .     *z**4*(1._dp + beta*z) - 2._dp*(1._dp - beta**2)*
     .     (2._dp*beta**2*z**2 + 3._dp*beta**3*z**3))*(xI1(mb**2,musq,ep) 
     .     - xI1(mw**2,musq,ep)))/((1._dp - beta**2)*s*
     .     (1._dp - beta*z)*(1._dp + beta**2 + 2._dp*beta*z)) + 
     .     (0.0625_dp*(-8._dp*(1._dp - beta**2)**2*mb**2*s*(1._dp 
     .     - beta**2 + 4._dp*beta**4 + beta*(1._dp - beta**2)*z 
     .     - 6._dp*beta**4*z**2 + 2._dp*beta**4*z**4) + (-16._dp*mb**4 
     .     + 4._dp*mw**2*(4._dp*mb**2 + (1._dp - beta**2)*s))*
     .     (1._dp + 8._dp*beta**2 - 11._dp*beta**4 + 4._dp*beta**6 + 
     .     beta*(1._dp + 4._dp*beta**2 - 3._dp*beta**4)*z - 2._dp*beta**2*
     .     (5._dp - 8._dp*beta**2 + 3._dp*beta**4)*z**2 - 6._dp*beta**3*
     .     (1._dp - beta**2)*z**3 - 2._dp*beta**4*
     .     (2._dp - beta**2)*z**4 - 2._dp*beta**5*z**5) 
     .     - (1._dp - beta**2)**2*s**2*(1._dp - 4._dp*beta**2 
     .     - 3._dp*beta**4 + 4._dp*beta**6 + beta*(1._dp - 8._dp*beta**2 
     .     + 5._dp*beta**4)*z + 2._dp*beta**2*(1._dp + 2._dp*beta**2 
     .     - 3._dp*beta**4)*z**2 + 6._dp*beta**3*(1._dp 
     .     - beta**2)*z**3 + 2._dp*beta**6*z**4 
     .     + 2._dp*beta**5*z**5))*xI2(mt**2,mb**2,mw**2,musq,ep))/
     .     ((1._dp - beta**2)*s*(1._dp - beta**2*z**2)) 
     .     + (0.0625_dp*((-16._dp*mb**4 + 4._dp*mw**2*(4._dp*mb**2 
     .     + (1._dp - beta**2)*s))*(1._dp - 4._dp*beta**2 
     .     - 3._dp*beta**4 + 4._dp*beta**6 + beta*(1._dp - 8._dp*beta**2 
     .     + 5._dp*beta**4)*z + 2._dp*beta**2*(1._dp + 2._dp*beta**2 
     .     - 3._dp*beta**4)*z**2 + 6._dp*beta**3*(1._dp 
     .     - beta**2)*z**3 + 2._dp*beta**6*z**4 + 2._dp*beta**5*z**5) 
     .     + (1._dp - beta**2)*s**2*(1._dp + beta**2 + 2._dp*beta*z)*
     .     (1._dp - 4._dp*beta**2 - 3._dp*beta**4 + 4._dp*beta**6 
     .     + beta*(1._dp - 8._dp*beta**2 + 5._dp*beta**4)*z 
     .     + 2._dp*beta**2*(1._dp + 2._dp*beta**2 - 3._dp*beta**4)*z**2 
     .     + 6._dp*beta**3*(1._dp - beta**2)*z**3 + 2._dp*beta**6*z**4 
     .     + 2._dp*beta**5*z**5) + 8._dp*beta*mb**2*s*(-beta*(3._dp 
     .     - 4._dp*beta**2 - beta**4 + 4._dp*beta**6) + (1._dp 
     .     - 11._dp*beta**2 + 13._dp*beta**4 - 7._dp*beta**6)*z 
     .     + beta*(5._dp - 18._dp*beta**2 + 5._dp*beta**4 
     .     + 6._dp*beta**6)*z**2 + 2._dp*beta**2*(5._dp - 11._dp*beta**2 
     .     + 6._dp*beta**4)*z**3 + 2._dp*beta**3*(5._dp - 3._dp*beta**2 
     .     - beta**4)*z**4 + 4._dp*beta**4*(2._dp 
     .     - beta**2)*z**5 + 2._dp*beta**5*z**6))*
     .     xI2(t,mb**2,mw**2,musq,ep))/(s*(1._dp + beta**2 
     .     + 2._dp*beta*z)*(1._dp - beta**2*z**2)) 
     .     - 0.25_dp*(1._dp - beta**2)*mb**2*(-4._dp*mb**2 + 3._dp*s 
     .     + beta**2*s + 4._dp*beta*s*z)*
     .     xI3(0._dp,mt**2,t,mb**2,mb**2,mw**2,musq,ep)))/(mw**2*pi)



      vrt(5) = 
     .     (alpha*genfacvert*gw**2*sigma0*
     .     (0.125_dp*(1._dp - beta**2)**2*s + (0.25_dp*(1._dp 
     .     - beta**2)*s*(-mh**2 + (1._dp - beta**2)*s)*
     .     (1._dp - beta**4*z**4 + beta*(1._dp - beta**2)*
     .     (2._dp*beta - z - 3._dp*beta*z**2))*
     .     db0(mt**2,mt**2,mh**2))/(1._dp - beta*z) 
     .     + (0.5_dp*(1._dp + 2._dp*beta**2 - beta**4 
     .     + beta*(1._dp + 4._dp*beta**2 - 3._dp*beta**4)*z 
     .     - 2._dp*beta**4*z**4*(1._dp + beta*z) 
     .     - 2._dp*(1._dp - beta**2)*(2._dp*beta**2*z**2 
     .     + 3._dp*beta**3*z**3))*(-xI1(mh**2,musq,ep) 
     .     + xI1(mt**2,musq,ep)))/((1._dp - beta*z)*
     .     (1._dp + beta**2 + 2._dp*beta*z)) + (0.125_dp*(-(1._dp 
     .     - beta**2)**2*s*(1._dp + 5._dp*beta**2 - 8._dp*beta**4 
     .     + beta*(1._dp - beta**2)*z - 6._dp*beta**2*(1._dp 
     .     - 2._dp*beta**2)*z**2 - 4._dp*beta**4*z**4) 
     .     + 2._dp*mh**2*(1._dp + 8._dp*beta**2 - 11._dp*beta**4 
     .     + 4._dp*beta**6 + beta*(1._dp + 4._dp*beta**2 - 3._dp*beta**4)*z 
     .     - 2._dp*beta**2*(5._dp - 8._dp*beta**2 + 3._dp*beta**4)*z**2 
     .     - 6._dp*beta**3*(1._dp - beta**2)*z**3 - 2._dp*beta**4*
     .     (2._dp - beta**2)*z**4 - 2._dp*beta**5*z**5))*
     .     xI2(mt**2,mt**2,mh**2,musq,ep))/(1._dp - beta**2*z**2) 
     .     + (0.125_dp*(1._dp - beta**2)*(2._dp*mh**2*(1._dp 
     .     - 4._dp*beta**2 - 3._dp*beta**4 + 4._dp*beta**6 + beta*(1._dp 
     .     - 8._dp*beta**2 + 5._dp*beta**4)*z + 2._dp*beta**2*(1._dp 
     .     + 2._dp*beta**2 - 3._dp*beta**4)*z**2 + 6._dp*beta**3*(1._dp 
     .     - beta**2)*z**3 + 2._dp*beta**6*z**4 + 2._dp*beta**5*z**5) 
     .     + beta*s*(beta*(3._dp - 8._dp*beta**2 - 5._dp*beta**4 
     .     + 8._dp*beta**6) + (1._dp + beta**2 - 23._dp*beta**4 
     .     + 17._dp*beta**6)*z - beta*(1._dp - 11._dp*beta**4 
     .     + 12._dp*beta**6)*z**2 - 2._dp*beta**2*(1._dp - 13._dp*beta**2 
     .     + 12._dp*beta**4)*z**3 + 2._dp*beta**3*(2._dp - 3._dp*beta**2 
     .     + 2._dp*beta**4)*z**4 - 4._dp*beta**4*(1._dp 
     .     - 2._dp*beta**2)*z**5 + 2._dp*beta**5*z**6))*
     .     xI2(t,mt**2,mh**2,musq,ep))/((1._dp + beta**2 + 2._dp*beta*z)*
     .     (1._dp - beta**2*z**2)) 
     .     + 0.0625_dp*(1._dp - beta**2)**2*s**2*(3._dp 
     .     - beta**2 + 2._dp*beta*z)*
     .     xI3(0._dp,mt**2,t,mt**2,mt**2,mh**2,musq,ep)))/(mw**2*pi)

      bx(1) = 
     .     (0.125_dp*alpha*genfacbox*sigma0*
     .     (-2._dp*beta*(gat**2 + gvt**2)*z*(1._dp - z**2) + 
     .     (8._dp*beta*(gat**2 + gvt**2)*(1._dp - z**2)*
     .     (4._dp*beta - z - 2._dp*beta*z**2)*(xI1(mt**2,musq,ep) 
     .     - xI1(mz**2,musq,ep)))/(s*(1._dp + beta**2 
     .     + 2._dp*beta*z)) + (4._dp*((-(gat**2 + gvt**2)*mz**2*
     .     (-6._dp*beta**3 - (5._dp - 8._dp*beta**2 + 4._dp*beta**4)*z 
     .     - beta*(5._dp - 12._dp*beta**2)*z**2 + (3._dp 
     .     - 4._dp*beta**2 + 2._dp*beta**4)*z**3 + beta*(3._dp 
     .     - 4._dp*beta**2)*z**4))/(beta*s) - gat**2*(3._dp*beta**2 
     .     - 5._dp*beta**4 - 3._dp*beta*(3._dp - 5._dp*beta**2 
     .     + 2._dp*beta**4)*z - (1._dp + 8._dp*beta**2 
     .     - 9._dp*beta**4)*z**2 + beta*(4._dp - 7._dp*beta**2 
     .     + 3._dp*beta**4)*z**3 + beta**2*(5._dp - 3._dp*beta**2)*z**4) 
     .     + gvt**2*(5._dp*beta**2 - 3._dp*beta**4 + beta*(1._dp + beta**2 
     .     - 2._dp*beta**4)*z + (1._dp - 4._dp*beta**2 + 3._dp*beta**4)*z**2 
     .     - beta**3*(1._dp - beta**2)*z**3 - beta**2*
     .     (1._dp + beta**2)*z**4))*xI2(mt**2,mt**2,mz**2,musq,ep))/
     .     (1._dp + beta*z) - 2._dp*((2._dp*(gat**2 + gvt**2)*mz**2*z*(5._dp 
     .     - 4._dp*beta**2 - (3._dp - 2._dp*beta**2)*z**2))/(beta*s) 
     .     - gat**2*(2._dp - 3._dp*beta*(5._dp - 4._dp*beta**2)*z 
     .     - 2._dp*z**2 + 3._dp*beta*(3._dp - 2._dp*beta**2)*z**3) 
     .     - gvt**2*(2._dp + beta*(1._dp - 4._dp*beta**2)*z 
     .     - 2._dp*z**2 + beta*(1._dp + 2._dp*beta**2)*z**3))*
     .     xI2(s,mt**2,mt**2,musq,ep) + (4._dp*((2._dp*beta*(gat**2 
     .     + gvt**2)*mz**2*(beta + z)*(1._dp - 3._dp*beta**2 + beta*(1._dp 
     .     - 2._dp*beta**2)*z + 2._dp*beta**2*z**2 + beta**3*z**3))/s 
     .     - gvt**2*(1._dp + 8._dp*beta**2 - 3._dp*beta**6 
     .     + beta*(4._dp + 14._dp*beta**2 - 11._dp*beta**4 
     .     - 2._dp*beta**6)*z - beta**2*(2._dp + beta**2 
     .     + 3._dp*beta**4)*z**2 - beta**3*(11._dp - 6._dp*beta**2 
     .     - beta**4)*z**3 - 2._dp*beta**4*(1._dp 
     .     - beta**2)*z**4 - beta**5*z**5) - gat**2*
     .     (1._dp + 5._dp*beta**6 + beta*(4._dp - 10._dp*beta**2 
     .     + 5._dp*beta**4 + 6._dp*beta**6)*z + beta**2*(2._dp 
     .     - 17._dp*beta**2 + 9._dp*beta**4)*z**2 + beta**3*(1._dp 
     .     - 2._dp*beta**2 - 3._dp*beta**4)*z**3 + 6._dp*beta**4*(1._dp 
     .     - beta**2)*z**4 - beta**5*z**5))*
     .     xI2(t,mt**2,mz**2,musq,ep))/((1._dp + beta*z)*(1._dp 
     .     + beta**2 + 2._dp*beta*z)) + 2._dp*((8._dp*(gat**2 
     .     + gvt**2)*mz**4)/s + 4._dp*mz**2*(-gat**2*(1._dp 
     .     - 4._dp*beta**2 - beta*z) + gvt**2*(3._dp + beta*z)) 
     .     + s*(gat**2*(4._dp - 6._dp*beta**2 + 7._dp*beta**4 
     .     - 2._dp*beta*(2._dp - 3._dp*beta**2)*z + beta**2*z**2) 
     .     + gvt**2*(4._dp + 2._dp*beta**2 - beta**4 
     .     + 2._dp*beta*(2._dp - beta**2)*z + beta**2*z**2)))*
     .     xI3(0._dp,0._dp,s,mt**2,mt**2,mt**2,musq,ep) - 2._dp*((8._dp* 
     .     (gat**2 + gvt**2)*mz**4*(1._dp + beta*z))/s + 4._dp*mz**2*(1._dp 
     .     + beta*z)*(-gat**2 + 4._dp*beta**2*gat**2 + 3._dp*gvt**2 
     .     + beta*(gat**2 + gvt**2)*z) + s*(gat**2*(4._dp - 6._dp*beta**2 
     .     + 7._dp*beta**4 - beta*(1._dp - 2._dp*beta**2 
     .     - 6._dp*beta**4)*z - 3._dp*beta**2*(1._dp - 2._dp*beta**2)*z**2 
     .     + beta**3*z**3) + gvt**2*(4._dp + 2._dp*beta**2 - beta**4 
     .     + beta*(7._dp + 2._dp*beta**2 - 2._dp*beta**4)*z + beta**2*(5._dp 
     .     - 2._dp*beta**2)*z**2 + beta**3*z**3)))*xI3(0._dp,mt**2,t,mt**2
     .     ,mt**2,mz**2,musq,ep) + 2._dp*((-2._dp*(gat**2 + gvt**2)*
     .     mz**4*z*(5._dp - 8._dp*beta**2 - 3._dp*z**2 
     .     + 2._dp*beta**2*z**2))/(beta*s) + beta*s*(gvt**2*(beta*(4._dp 
     .     - beta**2) + 4._dp*(1._dp + beta**2 - beta**4)*z 
     .     - beta**3*z**2 + beta**2*z**3 + beta**4*z**3) 
     .     - gat**2*(-3._dp*beta**3 - 4._dp*(1._dp - 3._dp*beta**2 
     .     + 3._dp*beta**4)*z + beta*(4._dp - 3._dp*beta**2)*z**2 
     .     - beta**2*(5._dp - 3._dp*beta**2)*z**3)) + 2._dp*mz**2*
     .     (gvt**2*(1._dp + beta**2 + 4._dp*beta*z - (1._dp 
     .     - beta**2)*z**2 + 2._dp*beta*z**3) + gat**2*(1._dp 
     .     + beta**2 - 4._dp*beta*(3._dp - 4._dp*beta**2)*z - (1._dp 
     .     - beta**2)*z**2 + 2._dp*beta*(3._dp - 2._dp*beta**2)*z**3))
     .     )*xI3(s,mt**2,mt**2,mt**2,mt**2,mz**2,musq,ep) + ((16._dp*
     .     (gat**2 + gvt**2)*mz**6)/s + 8._dp*mz**4*(-gat**2*(1._dp 
     .     - 5._dp*beta**2 - 2._dp*beta*z) + gvt**2*(3._dp + beta**2 
     .     + 2._dp*beta*z)) + 2._dp*mz**2*s*(gvt**2*(4._dp + 9._dp*beta**2 
     .     - 2._dp*beta**4 + 10._dp*beta*z + beta**2*(2._dp + beta**2)*z**2
     .     ) + gat**2*(4._dp - 7._dp*beta**2 + 14._dp*beta**4 - 2._dp*beta*
     .     (3._dp - 8._dp*beta**2)*z + beta**2*(2._dp + beta**2)*z**2)) 
     .     + s**2*(-gat**2*(3._dp - 7._dp*beta**2 + 5._dp*beta**4 
     .     - 6._dp*beta**6 - beta*(3._dp - 8._dp*beta**2 
     .     + 12._dp*beta**4)*z + 3._dp*beta**2*(1._dp - beta**2 
     .     - beta**4)*z**2 - beta**3*z**3) + gvt**2*(1._dp 
     .     + 7._dp*beta**2 - beta**4 - 2._dp*beta**6 + beta*(3._dp 
     .     + 8._dp*beta**2 - 4._dp*beta**4)*z + beta**2*(1._dp 
     .     + 3._dp*beta**2 - beta**4)*z**2 + beta**3*z**3)))*
     .     xI4(0._dp,0._dp,mt**2,mt**2,s,t,mt**2,mt**2,mt**2,mz**2,
     .     musq,ep)))/pi



      bx(2) = 
     .     (0.5_dp*alpha*genfacbox*gw**2*sigma0*
     .     (-beta*z*(1._dp - z**2) + (4._dp*beta*(1._dp 
     .     - z**2)*(4._dp*beta - z - 2._dp*beta*z**2)*
     .     (xI1(mb**2,musq,ep) - xI1(mw**2,musq,ep)))/(s*(1._dp 
     .     + beta**2 + 2._dp*beta*z)) + (0.5_dp*(4._dp*(-mb**2 
     .     + mw**2)*(6._dp*beta**3 + (5._dp - 8._dp*beta**2 + 4._dp*beta**4
     .     )*z + beta*(5._dp - 12._dp*beta**2)*z**2 - (3._dp 
     .     - 4._dp*beta**2 + 2._dp*beta**4)*z**3 - beta*(3._dp 
     .     - 4._dp*beta**2)*z**4) + s*(2._dp*beta**3*(5._dp - beta**2) 
     .     - beta*(3._dp + 5._dp*beta**2)*z**4 + (1._dp 
     .     - beta**2)*((5._dp + 12._dp*beta**2 - 4._dp*beta**4)*z 
     .     + 9._dp*beta*z**2 - (3._dp + 4._dp*beta**2 
     .     - 2._dp*beta**4)*z**3)))*xI2(mt**2,mb**2,mw**2,musq,ep))/
     .     (beta*s*(1._dp + beta*z)) - (0.5_dp*(4._dp*(-mb**2 
     .     + mw**2)*z*(5._dp - 4._dp*beta**2 - (3._dp - 2._dp*beta**2
     .     )*z**2) - s*(4._dp*beta - (5._dp + 5._dp*beta**2 
     .     - 4._dp*beta**4)*z - 4._dp*beta*z**2 + (3._dp + 5._dp*beta**2 
     .     - 2._dp*beta**4)*z**3))*xI2(s,mb**2,mb**2,musq,ep))/(beta*s) + 
     .     ((4._dp*beta*(-mb**2 + mw**2)*(beta + z)*(1._dp 
     .     - 3._dp*beta**2 + beta*(1._dp - 2._dp*beta**2)*z 
     .     + 2._dp*beta**2*z**2 + beta**3*z**3) - s*(1._dp + beta**2 
     .     + 2._dp*beta*z)*(2._dp + 5._dp*beta**2 - beta**4 + beta*
     .     (3._dp - 6._dp*beta**2 + 2._dp*beta**4)*z - beta**2*(7._dp 
     .     - 2._dp*beta**2)*z**2 + beta**3*(2._dp - beta**2)*z**3 
     .     - beta**4*z**4))*xI2(t,mb**2,mw**2,musq,ep))/(s*(1._dp 
     .     + beta*z)*(1._dp + beta**2 + 2._dp*beta*z)) + (0.5_dp*(16._dp*
     .     (-mb**2 + mw**2)**2 - 8._dp*beta*mb**2*s*(2._dp*beta + z) 
     .     + 8._dp*mw**2*s*(2._dp + beta**2 + beta*z) + s**2*(7._dp 
     .     + 2._dp*beta**2 + beta**4 + 2._dp*beta*(1._dp + beta**2)*z 
     .     + 2._dp*beta**2*z**2))*
     .     xI3(0._dp,0._dp,s,mb**2,mb**2,mb**2,musq,ep))
     .     /s - (0.5_dp*(16._dp*(-mb**2 + mw**2)**2*(1._dp + beta*z) 
     .     - 8._dp*beta*mb**2*s*(beta + z)*(2._dp + beta*z) + 8._dp*mw**2*
     .     s*(1._dp + beta*z)*(2._dp + beta**2 + beta*z) + s**2*(1._dp 
     .     + beta*z)*(7._dp + 2._dp*beta**2 + beta**4 + 2._dp*beta*(1._dp 
     .     + beta**2)*z + 2._dp*beta**2*z**2))*xI3(0._dp,mt**2,t,mb**2
     .     ,mb**2,mw**2,musq,ep))/s + (0.125_dp*(-16._dp*(-mb**2 
     .     + mw**2)**2*z*(5._dp - 8._dp*beta**2 - (3._dp 
     .     - 2._dp*beta**2)*z**2) - 8._dp*mb**2*s*(2._dp*beta*(1._dp 
     .     + beta**2) - (5._dp - 8._dp*beta**2)*(1._dp + beta**2)*z 
     .     - 2._dp*beta*(1._dp - beta**2)*z**2 + (3._dp 
     .     - 2._dp*beta**2)*(1._dp + beta**2)*z**3) + 8._dp*mw**2*s*
     .     (2._dp*beta*(1._dp + beta**2) - (5._dp - 5._dp*beta**2 
     .     - 8._dp*beta**4)*z - 2._dp*beta*(1._dp - beta**2)*z**2 + 
     .     (3._dp + 3._dp*beta**2 - 2._dp*beta**4)*z**3) - s**2*
     .     (-4._dp*beta*(1._dp + 4._dp*beta**2 + beta**4) + (5._dp 
     .     - 30._dp*beta**2 + beta**4 - 8._dp*beta**6)*z + 4._dp*beta*
     .     (1._dp + 2._dp*beta**2 - beta**4)*z**2 - (3._dp 
     .     + 4._dp*beta**2 + 11._dp*beta**4 - 2._dp*beta**6)*z**3))*
     .     xI3(s,mt**2,mt**2,mb**2,mb**2,mw**2,musq,ep))/(beta*s) 
     .     + (0.125_dp*(64._dp*(-mb**2 + mw**2)**3 + s**3*(1._dp 
     .     + beta**2 + 2._dp*beta*z)*(7._dp + 2._dp*beta**2 + beta**4 
     .     + 2._dp*beta*(1._dp + beta**2)*z + 2._dp*beta**2*z**2) 
     .     + 8._dp*s*(10._dp*mw**4 + 2._dp*(-mb**2 + mw**2)**2*
     .     (3._dp*beta**2 + 4._dp*beta*z) - 4._dp*mb**2*mw**2*(3._dp 
     .     + beta**2*z**2) + 2._dp*mb**4*(1._dp + 2._dp*beta**2*z**2)) 
     .     + 4._dp*s**2*(mw**2*(11._dp + 8._dp*beta**2 + 3._dp*beta**4 
     .     + 4._dp*beta*(3._dp + 2._dp*beta**2)*z + 6._dp*beta**2*z**2) 
     .     - mb**2*(11._dp - 4._dp*beta**2 + 3._dp*beta**4 
     .     + 8._dp*beta*(1._dp + beta**2)*z + 2._dp*beta**2*
     .     (6._dp + beta**2)*z**2)))*
     .     xI4(0._dp,0._dp,mt**2,mt**2,s,t,mb**2,mb**2,mb**2,mw**2,
     .     musq,ep))/s))/pi


      bx(3) = 
     .     (0.5_dp*alpha*chifac*genfacbox*mt**2*sigma0*
     .     (-beta*z*(1._dp - z**2) + (4._dp*beta*(1._dp 
     .     - z**2)*(4._dp*beta - z - 2._dp*beta*z**2)*
     .     (xI1(mt**2,musq,ep) - xI1(mz**2,musq,ep)))/(s*(1._dp 
     .     + beta**2 + 2._dp*beta*z)) + ((-2._dp*beta*(1._dp 
     .     - beta**2)*s*(2._dp - beta**2 + beta*z - z**2 
     .     - beta*z**3) + 2._dp*mz**2*(6._dp*beta**3 + (5._dp 
     .     - 8._dp*beta**2 + 4._dp*beta**4)*z + beta*(5._dp 
     .     - 12._dp*beta**2)*z**2 - (3._dp - 4._dp*beta**2 
     .     + 2._dp*beta**4)*z**3 - beta*(3._dp - 4._dp*beta**2)*z**4)
     .     )*xI2(mt**2,mt**2,mz**2,musq,ep))/(beta*s*(1._dp + beta*z)) 
     .     - (1._dp*(-beta*s*(2._dp + beta*z)*(1._dp - z**2) 
     .     + 2._dp*mz**2*z*(5._dp - 4._dp*beta**2 - (3._dp 
     .     - 2._dp*beta**2)*z**2))*xI2(s,mt**2,mt**2,musq,ep))/(beta*s) 
     .     + ((4._dp*beta*mz**2*(beta + z)*(1._dp - 3._dp*beta**2 
     .     + beta*(1._dp - 2._dp*beta**2)*z + 2._dp*beta**2*z**2 
     .     + beta**3*z**3) + 2._dp*s*(1._dp - 4._dp*beta**2 + beta**6 
     .     + beta*(2._dp - 10._dp*beta**2 + 3._dp*beta**4)*z + beta**2*
     .     (3._dp - 5._dp*beta**2)*z**2 + 2._dp*beta**3*(3._dp 
     .     - beta**2)*z**3 + 4._dp*beta**4*z**4 + beta**5*z**5))*
     .     xI2(t,mt**2,mz**2,musq,ep))/(s*(1._dp + beta*z)*(1._dp 
     .     + beta**2 + 2._dp*beta*z)) + ((8._dp*mz**4 + 4._dp*beta*mz**2*s*
     .     (beta + z) + s**2*(2._dp - 2._dp*beta**2 + beta**4 
     .     + 2._dp*beta*z + beta**2*z**2))*
     .     xI3(0._dp,0._dp,s,mt**2,mt**2,mt**2,musq,ep))/s
     .     - (1._dp*(s**2*(2._dp - 2._dp*beta**2 + beta**4 
     .     + 3._dp*beta*z + 3._dp*beta**2*z**2 + beta**3*z**3) + 4._dp*
     .     (1._dp + beta*z)*(2._dp*mz**4 + beta*mz**2*s*(beta + z)))*
     .     xI3(0._dp,mt**2,t,mt**2,mt**2,mz**2,musq,ep))/s - 
     .     (1._dp*(-beta**2*s**2*(beta + z)*(1._dp + beta*z) 
     .     + 2._dp*mz**4*z*(5._dp - 8._dp*beta**2 - (3._dp 
     .     - 2._dp*beta**2)*z**2) - 2._dp*beta*mz**2*s*(1._dp + beta**2 
     .     - 2._dp*beta*(1._dp - 2._dp*beta**2)*z - (1._dp 
     .     - beta**2)*z**2 + beta*(1._dp - beta**2)*z**3))*
     .     xI3(s,mt**2,mt**2,mt**2,mt**2,mz**2,musq,ep))/(beta*s) + 
     .     (0.5_dp*(16._dp*mz**6 + 16._dp*beta*mz**4*s*(beta + z) 
     .     + beta*s**3*(beta + z)*(1._dp + beta*z)**2 + 2._dp*mz**2*s**2*
     .     (2._dp - beta**2 + 2._dp*beta**4 + 2._dp*beta*(1._dp 
     .     + 2._dp*beta**2)*z + beta**2*(2._dp + beta**2)*z**2))*
     .     xI4(0._dp,0._dp,mt**2,mt**2,s,t,mt**2,mt**2,mt**2,mz**2,musq,
     .     ep))/s))/pi



      bx(4) = 
     .     (0.25_dp*alpha*genfacbox*gw**2*sigma0*
     .     (-0.25_dp*beta*(4._dp*mb**2 + (1._dp - beta**2)*s)*z*(1._dp 
     .     - z**2) + (beta*(4._dp*mb**2 + (1._dp - beta**2)*s)*
     .     (1._dp - z**2)*(4._dp*beta - z - 2._dp*beta*z**2)*
     .     (xI1(mb**2,musq,ep) - xI1(mw**2,musq,ep)))/(s*(1._dp 
     .     + beta**2 + 2._dp*beta*z)) - (0.125_dp*((1._dp - beta**2
     .     )**2*s**2*(8._dp*beta - 2._dp*beta**3 - (5._dp 
     .     - 4._dp*beta**2 - 4._dp*beta**4)*z - 9._dp*beta*z**2 + (3._dp 
     .     - 4._dp*beta**2 - 2._dp*beta**4)*z**3 + 3._dp*beta*z**4) 
     .     + 16._dp*beta*(1._dp - beta**2)*mb**2*s*(2._dp 
     .     - 3._dp*beta**2 + beta*(3._dp - 2._dp*beta**2)*z 
     .     - (1._dp - 3._dp*beta**2)*z**2 - beta*(2._dp 
     .     - beta**2)*z**3 - beta**2*z**4) + (16._dp*mb**2*
     .     (mb**2 - mw**2) - 4._dp*(1._dp - beta**2)*mw**2*s)*
     .     (6._dp*beta**3 + (5._dp - 8._dp*beta**2 + 4._dp*beta**4)*z 
     .     + beta*(5._dp - 12._dp*beta**2)*z**2 - (3._dp 
     .     - 4._dp*beta**2 + 2._dp*beta**4)*z**3 - beta*(3._dp 
     .     - 4._dp*beta**2)*z**4))*xI2(mt**2,mb**2,mw**2,musq,ep))/
     .     (beta*s*(1._dp + beta*z)) + (0.125_dp*((16._dp*mb**2*(mb**2 
     .     - mw**2) - 4._dp*(1._dp - beta**2)*mw**2*s)*z*(5._dp 
     .     - 4._dp*beta**2 - (3._dp - 2._dp*beta**2)*z**2) 
     .     + 8._dp*beta*mb**2*s*(2._dp + beta*(5._dp - 4._dp*beta**2)*z 
     .     - 2._dp*z**2 - beta*(3._dp - 2._dp*beta**2)*z**3) + (1._dp 
     .     - beta**2)*s**2*(4._dp*beta - (5._dp - 3._dp*beta**2 
     .     - 4._dp*beta**4)*z - 4._dp*beta*z**2 + (3._dp - 3._dp*beta**2 
     .     - 2._dp*beta**4)*z**3))*xI2(s,mb**2,mb**2,musq,ep))/(beta*s) 
     .     - (0.25_dp*(beta*(16._dp*mb**2*(mb**2 - mw**2) 
     .     - 4._dp*(1._dp - beta**2)*mw**2*s)*(beta + z)*(1._dp 
     .     - 3._dp*beta**2 + beta*(1._dp - 2._dp*beta**2)*z 
     .     + 2._dp*beta**2*z**2 + beta**3*z**3) - (1._dp 
     .     - beta**2)*s**2*(1._dp + beta**2 + 2._dp*beta*z)*
     .     (2._dp - 5._dp*beta**2 + beta**4 + beta*(1._dp - 2._dp*beta**2 
     .     - 2._dp*beta**4)*z + beta**2*(3._dp - 2._dp*beta**2)*z**2 
     .     + beta**3*(2._dp + beta**2)*z**3 + beta**4*z**4) 
     .     - 8._dp*mb**2*s*(1._dp - 6._dp*beta**2 + 3._dp*beta**6 
     .     + beta*(2._dp - 16._dp*beta**2 + 7._dp*beta**4 + 2._dp*beta**6)*z 
     .     + beta**2*(4._dp - 9._dp*beta**2 + 3._dp*beta**4)*z**2 
     .     + beta**3*(9._dp - 4._dp*beta**2 - beta**4)*z**3 
     .     + 2._dp*beta**4*(3._dp - beta**2)*z**4 + beta**5*z**5))*
     .     xI2(t,mb**2,mw**2,musq,ep))/(s*(1._dp + beta*z)*(1._dp 
     .     + beta**2 + 2._dp*beta*z)) + (0.125_dp*(64._dp*mb**6 
     .     + 16._dp*mw**4*(4._dp*mb**2 + (1._dp - beta**2)*s) 
     .     + 8._dp*beta*(1._dp - beta**2)*mw**2*s**2*(beta + z) 
     .     + 32._dp*mb**2*mw**2*s*(1._dp + beta*z) + 16._dp*mb**4*
     .     (-8._dp*mw**2 + (1._dp - beta**2)*s - 2._dp*beta*s*z) 
     .     + 4._dp*mb**2*s**2*(3._dp - 2._dp*beta**2 + beta**4 
     .     + 4._dp*beta*(2._dp - beta**2)*z + 2._dp*beta**2*z**2) 
     .     + (1._dp - beta**2)*s**3*(3._dp - 2._dp*beta**2 + beta**4 
     .     + 2._dp*beta*(1._dp + beta**2)*z + 2._dp*beta**2*z**2))*
     .     xI3(0._dp,0._dp,s,mb**2,mb**2,mb**2,musq,ep))/s 
     .     - (0.125_dp*(64._dp*mb**6*(1._dp + beta*z) 
     .     + 16._dp*(1._dp - beta**2)*mw**4*s*(1._dp + beta*z) 
     .     + 8._dp*beta*(1._dp - beta**2)*mw**2*s**2*(beta + z)*
     .     (1._dp + beta*z) + (1._dp - beta**2)*s**3*(1._dp + beta*z)*
     .     (3._dp - 2._dp*beta**2 + beta**4 + 2._dp*beta*(1._dp + beta**2)*z 
     .     + 2._dp*beta**2*z**2) - 16._dp*mb**4*(8._dp*mw**2*(1._dp 
     .     + beta*z) - s*(1._dp - beta**2 - beta*(3._dp 
     .     - beta**2)*z - 2._dp*beta**2*z**2)) + mb**2*(4._dp*s**2*
     .     (3._dp - 2._dp*beta**2 + beta**4 + 9._dp*beta*z - 2._dp*beta**3*z 
     .     - beta**5*z + 10._dp*beta**2*z**2 - 4._dp*beta**4*z**2 
     .     + 2._dp*beta**3*z**3) + 32._dp*(1._dp + beta*z)*(2._dp*mw**4 
     .     + mw**2*s*(1._dp + beta*z))))*xI3(0._dp,mt**2,t,mb**2,mb**2
     .     ,mw**2,musq,ep))/s + (0.03125_dp*(-16._dp*(4._dp*mb**6 
     .     - 8._dp*mb**4*mw**2 + mw**4*(4._dp*mb**2 + (1._dp - beta**2
     .     )*s))*z*(5._dp - 8._dp*beta**2 - (3._dp - 2._dp*beta**2
     .     )*z**2) + 64._dp*beta*mb**2*mw**2*s*(1._dp + beta**2 
     .     + 2._dp*beta*z - (1._dp - beta**2)*z**2) 
     .     + 8._dp*(1._dp - beta**2)*mw**2*s**2*(2._dp*beta*(1._dp 
     .     + beta**2) - (5._dp - beta**2 - 8._dp*beta**4)*z 
     .     - 2._dp*beta*(1._dp - beta**2)*z**2 + (3._dp - beta**2 
     .     - 2._dp*beta**4)*z**3) - 16._dp*mb**4*s*(4._dp*beta*(1._dp 
     .     + beta**2) - (5._dp - 17._dp*beta**2 + 8._dp*beta**4)*z 
     .     - 4._dp*beta*(1._dp - beta**2)*z**2 + (3._dp - 9._dp*beta**2 
     .     + 2._dp*beta**4)*z**3) + 4._dp*mb**2*s**2*(8._dp*beta**3*(2._dp 
     .     - beta**2) + (5._dp - 2._dp*beta**2 + 21._dp*beta**4 
     .     - 8._dp*beta**6)*z + 8._dp*beta**3*(2._dp - beta**2)*z**2 
     .     - (3._dp - beta**4 - 2._dp*beta**6)*z**3) + (1._dp 
     .     - beta**2)*s**3*(4._dp*beta*(1._dp + beta**4) - (5._dp 
     .     - 22._dp*beta**2 + 9._dp*beta**4 - 8._dp*beta**6)*z 
     .     - 4._dp*beta*(1._dp - 2._dp*beta**2 - beta**4)*z**2 
     .     + (3._dp - 4._dp*beta**2 + 3._dp*beta**4 - 2._dp*beta**6)*z**3))*
     .     xI3(s,mt**2,mt**2,mb**2,mb**2,mw**2,musq,ep))/(beta*s) 
     .     - (0.03125_dp*(256._dp*mb**8 - 64._dp*(1._dp - beta**2)*
     .     mw**6*s - 16._dp*(1._dp - beta**2)*mw**4*s**2*(1._dp 
     .     + 3._dp*beta**2 + 4._dp*beta*z) - (1._dp - beta**2)*
     .     s**4*(1._dp + beta**2 + 2._dp*beta*z)*(3._dp - 2._dp*beta**2 
     .     + beta**4 + 2._dp*beta*(1._dp + beta**2)*z + 2._dp*beta**2*z**2) 
     .     - 4._dp*(1._dp - beta**2)*mw**2*s**3*(3._dp + 3._dp*beta**4 
     .     + 4._dp*beta*(1._dp + 2._dp*beta**2)*z + 6._dp*beta**2*z**2) 
     .     - 128._dp*mb**6*(6._dp*mw**2 + beta*s*z*(2._dp + beta*z)) 
     .     + 32._dp*mb**4*(24._dp*mw**4 + 2._dp*mw**2*s*(3._dp + beta**2 
     .     + 8._dp*beta*z + 2._dp*beta**2*z**2) - s**2*(1._dp 
     .     - 3._dp*beta**2 + beta**4 - 6._dp*beta*z + 2._dp*beta**3*z 
     .     - 5._dp*beta**2*z**2 + 2._dp*beta**4*z**2)) - 8._dp*mb**2*
     .     (32._dp*mw**6 + 16._dp*mw**4*s*(1._dp + beta**2 + 2._dp*beta*z) 
     .     + 2._dp*mw**2*s**2*(5._dp + beta**4 + 12._dp*beta*z 
     .     + 2._dp*beta**2*(2._dp + beta**2)*z**2) + s**3*(-2._dp 
     .     + 6._dp*beta**2 - 2._dp*beta**4 + 4._dp*beta*z + 4._dp*beta**3*z 
     .     - 2._dp*beta**5*z + 9._dp*beta**2*z**2 - 4._dp*beta**4*z**2 
     .     + beta**6*z**2 + 2._dp*beta**3*z**3)))*xI4(0._dp,0._dp,mt**2,
     .     mt**2,s,t,mb**2,mb**2,mb**2,mw**2,musq,ep))/s))/(mw**2*pi)



      bx(5) = 
     .     (-0.25_dp*alpha*genfacbox*gw**2*mt**2*sigma0*
     .     (beta*z*(1._dp - z**2) - (4._dp*beta*(1._dp - z**2)*
     .     (4._dp*beta - z - 2._dp*beta*z**2)*(-xI1(mh**2,musq,
     .     ep) + xI1(mt**2,musq,ep)))/(s*(1._dp + beta**2 + 2._dp*beta*z)) 
     .     + ((2._dp*beta*(1._dp - beta**2)*s*(2._dp + 3._dp*beta**2 
     .     - beta*(3._dp - 4._dp*beta**2)*z - (1._dp 
     .     + 6._dp*beta**2)*z**2 + beta*(1._dp - 2._dp*beta**2)*z**3 
     .     + 2._dp*beta**2*z**4) + 2._dp*mh**2*(-6._dp*beta**3 - (5._dp 
     .     - 8._dp*beta**2 + 4._dp*beta**4)*z - beta*(5._dp 
     .     - 12._dp*beta**2)*z**2 + (3._dp - 4._dp*beta**2 + 2._dp*beta**4
     .     )*z**3 + beta*(3._dp - 4._dp*beta**2)*z**4))*xI2(mt**2,mt**2,
     .     mh**2,musq,ep))/(beta*s*(1._dp + beta*z)) + ((2._dp*mh**2*z*
     .     (5._dp - 4._dp*beta**2 - (3._dp - 2._dp*beta**2)*z**2) 
     .     - beta*s*(2._dp - beta*(7._dp - 8._dp*beta**2)*z 
     .     - 2._dp*z**2 + beta*(3._dp - 4._dp*beta**2)*z**3))*
     .     xI2(s,mt**2,mt**2,musq,ep))/(beta*s) + ((-4._dp*beta*mh**2*
     .     (beta + z)*(1._dp - 3._dp*beta**2 + beta*(1._dp - 2._dp*beta**2
     .     )*z + 2._dp*beta**2*z**2 + beta**3*z**3) - 2._dp*s*(1._dp 
     .     - 3._dp*beta**6 + beta*(2._dp + 2._dp*beta**2 - 5._dp*beta**4 
     .     - 4._dp*beta**6)*z + beta**2*(1._dp + 3._dp*beta**2 
     .     - 6._dp*beta**4)*z**2 + 2._dp*beta**5*(1._dp + beta**2)*z**3 
     .     + 4._dp*beta**6*z**4 + beta**5*z**5))*xI2(t,mt**2,mh**2,musq,
     .     ep))/(s*(1._dp + beta*z)*(1._dp + beta**2 + 2._dp*beta*z)) - 
     .     (1._dp*(8._dp*mh**4 - 4._dp*mh**2*s*(2._dp - 3._dp*beta**2 
     .     - beta*z) + s**2*(6._dp - 10._dp*beta**2 + 5._dp*beta**4 
     .     - 2._dp*beta*(1._dp - 2._dp*beta**2)*z + beta**2*z**2))*
     .     xI3(0._dp,0._dp,s,mt**2,mt**2,mt**2,musq,ep))/s + ((s**2*(6._dp 
     .     - 10._dp*beta**2 + 5._dp*beta**4 + beta*(3._dp - 4._dp*beta**2 
     .     + 4._dp*beta**4)*z - beta**2*(1._dp - 4._dp*beta**2)*z**2 
     .     + beta**3*z**3) + (1._dp + beta*z)*(8._dp*mh**4 - 4._dp*mh**2*s*
     .     (2._dp - 3._dp*beta**2 - beta*z)))*xI3(0._dp,mt**2,t,mt**2
     .     ,mt**2,mh**2,musq,ep))/s + ((2._dp*mh**4*z*(5._dp 
     .     - 8._dp*beta**2 - (3._dp - 2._dp*beta**2)*z**2) 
     .     - 2._dp*beta*mh**2*s*(1._dp + beta**2 - 2._dp*beta*(5._dp 
     .     - 6._dp*beta**2)*z - (1._dp - beta**2)*z**2 
     .     + 3._dp*beta*(1._dp - beta**2)*z**3) + beta**2*s**2*
     .     (beta - 2._dp*beta**3 - (7._dp - 13._dp*beta**2 
     .     + 8._dp*beta**4)*z + beta*(1._dp - 2._dp*beta**2)*z**2 
     .     - 2._dp*beta**2*(1._dp - beta**2)*z**3))*
     .     xI3(s,mt**2,mt**2,mt**2,mt**2,mh**2,musq,ep))/(beta*s) 
     .     + (0.5_dp*(-16._dp*mh**6 + 16._dp*mh**4*s*(1._dp - 2._dp*beta**2 
     .     - beta*z) - 2._dp*mh**2*s**2*(6._dp - 13._dp*beta**2 
     .     + 10._dp*beta**4 - 6._dp*beta*(1._dp - 2._dp*beta**2)*z 
     .     + beta**2*(2._dp + beta**2)*z**2) - s**3*(2._dp + beta**2 
     .     - 6._dp*beta**4 + 4._dp*beta**6 + beta*(5._dp - 10._dp*beta**2 
     .     + 8._dp*beta**4)*z + beta**4*(1._dp + 2._dp*beta**2)*z**2 
     .     + beta**3*z**3))*xI4(0._dp,0._dp,mt**2,mt**2,s,t,mt**2,mt**2,
     .     mt**2,mh**2,musq,ep))/s))/(mw**2*pi)

c--- use the value of "tevscale" (passed via anomcoup.f)
c--- as an anomalous top Yukawa coupling:
c--- g(top Y) = (tevscale) x g(SM, top Y)
c--- it therefore affects Higgs diagrams as the square
      vrts(5)=vrts(5)*tevscale**2
      slf(5) =slf(5) *tevscale**2
      vrt(5) =vrt(5) *tevscale**2
      bx(5)  =bx(5)  *tevscale**2

      vew = 
     .     + (vrts(1) + vrts(2) + vrts(3) + vrts(4) + vrts(5))/2._dp
     .     + slf(1) + slf(2) + slf(3) + slf(4) + slf(5)
     .     + vrt(1) + vrt(2) + vrt(3) + vrt(4) + vrt(5)
     .     + bx(1) + bx(2) + bx(3) + bx(4) + bx(5)

      vew = vew + (trih + trizx)/2._dp

      vew = vew/born


      end subroutine ggQQb_ew_oneloop

