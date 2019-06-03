      subroutine z_ew_box(corr,ss,tt,uu)
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'constants.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'masses.f'
      include 'scale.f'
      include 'zcouple.f'
      real(dp):: corr(nf,8),ss,tt,uu,gv1(nf),ga1(nf),gv2,ga2,
     .     xI2,xI3,xI4,fac,Mvsq,mz,mw,aemmz,propa,propz,D0t,D0u,Qq(nf),
     .     Ql,bxzt(nf),bxzu(nf),bxwt(nf),bxwu(nf),bxzat(nf),bxzau(nf),
     .     bxwat(nf),bxwau(nf)
      integer ep
      common/em/aemmz

C -- no poles in Box contribution, finite part needed only
      ep = 0
      mz = zmass
      mw = wmass

      gv1(1:5) = (l(1:5) + r(1:5))/2._dp
      ga1(1:5) = (l(1:5) - r(1:5))/2._dp

      gv2 = (l1 + r1)/2._dp
      ga2 = (l1 - r1)/2._dp
      gw = 1._dp/sqrt(2._dp*xw)
      Qq(1:nf) = Q(1:nf)
      Ql = q1

C -- prefactor first '2' comes from 2*dM*M, coupling 'el**6', spin & color ave
      fac = 2._dp*aemmz**3*4._dp*pi/xn/4._dp
C -- propz for Z and propa for photon in Born
      propz = (ss-mz**2)/((ss-mz**2)**2+mz**2*zwidth**2)
      propa = 1._dp/ss

C -- set gauge boson mass to that of W
      Mvsq = mw**2

C -- light initial quark no bottom with W in box, 
C -- internal fermions remain massless
      bxwt(1:4) = 
     .     + 4._dp*(ga1(1:4) + gv1(1:4))*(ga2 + gv2)*gw**4*(-2._dp*uu*
     .     xI2(ss,Mvsq,Mvsq,musq,ep) + 2._dp*uu*xI2(tt,zip,zip,musq,
     .     ep) + (tt**2 + uu**2 + 2._dp*Mvsq*(tt + uu))*xI3(zip,zip,
     .     ss,Mvsq,zip,Mvsq,musq,ep) - tt*(2._dp*Mvsq + tt - uu)*
     .     xI3(zip,zip,tt,zip,Mvsq,zip,musq,ep) - tt*(2._dp*Mvsq 
     .     + tt - uu)*xI3(zip,tt,zip,Mvsq,zip,zip,musq,ep) + (tt**2 
     .     + uu**2 + 2._dp*Mvsq*(tt + uu))*xI3(ss,zip,zip,Mvsq,
     .     Mvsq,zip,musq,ep) - (4._dp*Mvsq*tt**2 + 2._dp*Mvsq**2*(2._dp*tt 
     .     + uu) + tt*(tt**2 + uu**2))*xI4(zip,zip,zip,zip,ss,tt,
     .     Mvsq,zip,Mvsq,zip,musq,ep))

      bxwu = 
     .     + 8._dp*(ga1 + gv1)*(ga2 + gv2)*gw**4*uu**2*(-xI3(zip,zip,
     .     ss,Mvsq,zip,Mvsq,musq,ep) - xI3(ss,zip,zip,Mvsq,
     .     Mvsq,zip,musq,ep) + uu*xI4(zip,zip,zip,zip,ss,uu,
     .     Mvsq,zip,Mvsq,zip,musq,ep))

      bxwat(1:4) = 
     .     + 4._dp*gw**4*Ql*Qq(1:4)*(-2._dp*uu*xI2(ss,Mvsq,Mvsq,musq,ep) 
     .     + 2._dp*uu*xI2(tt,zip,zip,musq,ep) + (tt**2 + uu**2 
     .     + 2._dp*Mvsq*(tt + uu))*xI3(zip,zip,ss,Mvsq,zip,Mvsq,
     .     musq,ep) - tt*(2._dp*Mvsq + tt - uu)*xI3(zip,zip,tt,zip,
     .     Mvsq,zip,musq,ep) - tt*(2._dp*Mvsq + tt - uu)*
     .     xI3(zip,tt,zip,Mvsq,zip,zip,musq,ep) + (tt**2 + uu**2 
     .     + 2._dp*Mvsq*(tt + uu))*xI3(ss,zip,zip,Mvsq,Mvsq,zip,
     .     musq,ep) - (4._dp*Mvsq*tt**2 + 2._dp*Mvsq**2*(2._dp*tt + uu) 
     .     + tt*(tt**2 + uu**2))*xI4(zip,zip,zip,zip,ss,tt,
     .     Mvsq,zip,Mvsq,zip,musq,ep))


      bxwau = 
     .     + 8._dp*gw**4*Ql*Qq*uu**2*(-xI3(zip,zip,ss,Mvsq,zip,Mvsq,
     .     musq,ep) - xI3(ss,zip,zip,Mvsq,Mvsq,zip,musq,ep) 
     .     + uu*xI4(zip,zip,zip,zip,ss,uu,Mvsq,zip,Mvsq,zip,
     .     musq,ep))

C -- when initial quark is bottom, W correction gives virtual top in box
      bxwt(5) = 
     .     + 4._dp*(ga1(5) + gv1(5))*(ga2 + gv2)*gw**4*uu**2*((2._dp*
     .     xI2(ss,Mvsq,Mvsq,musq,ep))/(ss + tt) - (2._dp*xI2(tt,
     .     Mt**2,zip,musq,ep))/(ss + tt) + ((Mt**2*ss - 2._dp*Mvsq*ss 
     .     + ss**2 + 2._dp*ss*tt + 2._dp*tt**2)*xI3(zip,zip,ss,Mvsq,
     .     Mt**2,Mvsq,musq,ep))/(ss + tt)**2 - (tt*(-Mt**2 + 2._dp*Mvsq 
     .     + ss + 2._dp*tt)*xI3(zip,zip,tt,Mt**2,Mvsq,zip,musq,ep)) 
     .     /(ss + tt)**2 - (tt*(-Mt**2 + 2._dp*Mvsq + ss + 2._dp*tt)*
     .     xI3(zip,tt,zip,Mvsq,Mt**2,zip,musq,ep))/(ss + tt)**2 
     .     + ((Mt**2*ss - 2._dp*Mvsq*ss + ss**2 + 2._dp*ss*tt 
     .     + 2._dp*tt**2)*xI3(ss,zip,zip,Mvsq,Mvsq,zip,musq,ep))/(ss 
     .     + tt)**2 - 2._dp*tt*xI4(zip,zip,zip,zip,ss,tt,Mvsq,zip,
     .     Mvsq,zip,musq,ep) + ((Mt**4*ss + 2._dp*Mvsq**2*(ss - tt) 
     .     - 4._dp*Mvsq*tt**2 + ss*tt*(ss + 2._dp*tt) + Mt**2*(ss**2 
     .     - 2._dp*Mvsq*(ss - tt) + ss*tt + 2._dp*tt**2))*
     .     xI4(zip,zip,zip,zip,ss,tt,Mvsq,Mt**2,Mvsq,zip,
     .     musq,ep))/(ss + tt)**2)


      bxwat(5) = 
     .     + 4._dp*gw**4*Ql*Qq(5)*uu**2*((2._dp*xI2(ss,Mvsq,Mvsq,musq,
     .     ep))/(ss + tt) - (2._dp*xI2(tt,Mt**2,zip,musq,ep))/(ss + tt) 
     .     + ((Mt**2*ss - 2._dp*Mvsq*ss + ss**2 + 2._dp*ss*tt 
     .     + 2._dp*tt**2)*xI3(zip,zip,ss,Mvsq,Mt**2,Mvsq,musq,ep))
     .     /(ss + tt)**2 - (tt*(-Mt**2 + 2._dp*Mvsq + ss + 2._dp*tt)*
     .     xI3(zip,zip,tt,Mt**2,Mvsq,zip,musq,ep))/(ss + tt)**2 
     .     - (tt*(-Mt**2 + 2._dp*Mvsq + ss + 2._dp*tt)*xI3(zip,tt,zip,
     .     Mvsq,Mt**2,zip,musq,ep))/(ss + tt)**2 + ((Mt**2*ss 
     .     - 2._dp*Mvsq*ss + ss**2 + 2._dp*ss*tt + 2._dp*tt**2)*
     .     xI3(ss,zip,zip,Mvsq,Mvsq,zip,musq,ep))/(ss + tt)**2 
     .     - 2._dp*tt*xI4(zip,zip,zip,zip,ss,tt,Mvsq,zip,
     .     Mvsq,zip,musq,ep) + ((Mt**4*ss + 2._dp*Mvsq**2*(ss - tt) 
     .     - 4._dp*Mvsq*tt**2 + ss*tt*(ss + 2._dp*tt) + Mt**2*(ss**2 
     .     - 2._dp*Mvsq*(ss - tt) + ss*tt + 2._dp*tt**2))*
     .     xI4(zip,zip,zip,zip,ss,tt,Mvsq,Mt**2,Mvsq,zip,musq,
     .     ep))/(ss + tt)**2)
     

C -- set gauge boson mass to that of Z
      Mvsq = mz**2

      bxzt = 
     .     + 8._dp*(-2._dp*(ga1*ga2 + gv1*gv2)*(8._dp*ga1*ga2*gv1*gv2 
     .     + gv1**2*(3._dp*ga2**2 + gv2**2) + ga1**2*(ga2**2 
     .     + 3._dp*gv2**2))*uu*xI2(ss,Mvsq,Mvsq,musq,ep) 
     .     + 2._dp*(ga1*ga2 + gv1*gv2)*(8._dp*ga1*ga2*gv1*gv2 
     .     + gv1**2*(3._dp*ga2**2 + gv2**2) + ga1**2*(ga2**2 
     .     + 3._dp*gv2**2))*uu*xI2(tt,zip,zip,musq,ep) 
     .     + (ga1**3*ga2*(ga2**2 + 3._dp*gv2**2)*(2._dp*Mvsq - tt + uu)*
     .     (tt + uu) + 3._dp*ga1*ga2*gv1**2*(ga2**2 + 3._dp*gv2**2)*
     .     (2._dp*Mvsq - tt + uu)*(tt + uu) + 3._dp*ga1**2*gv1*gv2*
     .     (3._dp*ga2**2 + gv2**2)*(3._dp*tt**2 + uu**2 + 2._dp*Mvsq*(tt 
     .     + uu)) + gv1**3*gv2*(3._dp*ga2**2 + gv2**2)*(3._dp*tt**2 
     .     + uu**2 + 2._dp*Mvsq*(tt + uu)))*xI3(zip,zip,ss,Mvsq,zip,
     .     Mvsq,musq,ep) - (ga1*ga2 + gv1*gv2)*(8._dp*ga1*ga2*gv1*gv2 
     .     + gv1**2*(3._dp*ga2**2 + gv2**2) + ga1**2*(ga2**2 
     .     + 3._dp*gv2**2))*tt*(2._dp*Mvsq + tt - uu)*
     .     xI3(zip,zip,tt,zip,Mvsq,zip,musq,ep) - (ga1*ga2 
     .     + gv1*gv2)*(8._dp*ga1*ga2*gv1*gv2 + gv1**2*(3._dp*ga2**2 
     .     + gv2**2) + ga1**2*(ga2**2 + 3._dp*gv2**2))*tt*(2._dp*Mvsq 
     .     + tt - uu)*xI3(zip,tt,zip,Mvsq,zip,zip,musq,ep) 
     .     + (ga1**3*ga2*(ga2**2 + 3._dp*gv2**2)*(2._dp*Mvsq - tt + uu)*
     .     (tt + uu) + 3._dp*ga1*ga2*gv1**2*(ga2**2 + 3._dp*gv2**2)*
     .     (2._dp*Mvsq - tt + uu)*(tt + uu) + 3._dp*ga1**2*gv1*gv2*
     .     (3._dp*ga2**2 + gv2**2)*(3._dp*tt**2 + uu**2 + 2._dp*Mvsq*(tt 
     .     + uu)) + gv1**3*gv2*(3._dp*ga2**2 + gv2**2)*(3._dp*tt**2 
     .     + uu**2 + 2._dp*Mvsq*(tt + uu)))*xI3(ss,zip,zip,Mvsq,
     .     Mvsq,zip,musq,ep) - (ga1**3*ga2*(ga2**2 + 3._dp*gv2**2)*
     .     (4._dp*Mvsq*tt**2 - tt**3 + tt*uu**2 + 2._dp*Mvsq**2*(2._dp*tt 
     .     + uu)) + 3._dp*ga1*ga2*gv1**2*(ga2**2 + 3._dp*gv2**2)*
     .     (4._dp*Mvsq*tt**2 - tt**3 + tt*uu**2 + 2._dp*Mvsq**2*(2._dp*tt 
     .     + uu)) + 3._dp*ga1**2*gv1*gv2*(3._dp*ga2**2 + gv2**2)*
     .     (4._dp*Mvsq*tt**2 + 2._dp*Mvsq**2*(2._dp*tt + uu) 
     .     + tt*(3._dp*tt**2 + uu**2)) + gv1**3*gv2*(3._dp*ga2**2 
     .     + gv2**2)*(4._dp*Mvsq*tt**2 + 2._dp*Mvsq**2*(2._dp*tt + uu) 
     .     + tt*(3._dp*tt**2 + uu**2)))*xI4(zip,zip,zip,zip,ss,tt,
     .     Mvsq,zip,Mvsq,zip,musq,ep))

      bxzu = 
     .     + 8._dp*(-2._dp*(ga1*ga2 - gv1*gv2)*(-8._dp*ga1*ga2*gv1*gv2 
     .     + gv1**2*(3._dp*ga2**2 + gv2**2) + ga1**2*(ga2**2 
     .     + 3._dp*gv2**2))*tt*xI2(ss,Mvsq,Mvsq,musq,ep) 
     .     + 2._dp*(ga1*ga2 - gv1*gv2)*(-8._dp*ga1*ga2*gv1*gv2 
     .     + gv1**2*(3._dp*ga2**2 + gv2**2) + ga1**2*(ga2**2 
     .     + 3._dp*gv2**2))*tt*xI2(uu,zip,zip,musq,ep) 
     .     + (ga1**3*ga2*(ga2**2 + 3._dp*gv2**2)*(2._dp*Mvsq + tt - uu)*
     .     (tt + uu) + 3._dp*ga1*ga2*gv1**2*(ga2**2 + 3._dp*gv2**2)*
     .     (2._dp*Mvsq + tt - uu)*(tt + uu) - 3._dp*ga1**2*gv1*gv2*
     .     (3._dp*ga2**2 + gv2**2)*(tt**2 + 3._dp*uu**2 + 2._dp*Mvsq*(tt 
     .     + uu)) - gv1**3*gv2*(3._dp*ga2**2 + gv2**2)*(tt**2 
     .     + 3._dp*uu**2 + 2._dp*Mvsq*(tt + uu)))*xI3(zip,zip,ss,
     .     Mvsq,zip,Mvsq,musq,ep) + (ga1*ga2 - gv1*gv2)*
     .     (-8._dp*ga1*ga2*gv1*gv2 + gv1**2*(3._dp*ga2**2 + gv2**2) 
     .     + ga1**2*(ga2**2 + 3._dp*gv2**2))*(-2._dp*Mvsq + tt - uu)*uu*
     .     xI3(zip,zip,uu,zip,Mvsq,zip,musq,ep) + (ga1*ga2 
     .     - gv1*gv2)*(-8._dp*ga1*ga2*gv1*gv2 + gv1**2*(3._dp*ga2**2 
     .     + gv2**2) + ga1**2*(ga2**2 + 3._dp*gv2**2))*(-2._dp*Mvsq + tt 
     .     - uu)*uu*xI3(zip,uu,zip,Mvsq,zip,zip,musq,ep) 
     .     + (ga1**3*ga2*(ga2**2 + 3._dp*gv2**2)*(2._dp*Mvsq + tt - uu)*
     .     (tt + uu) + 3._dp*ga1*ga2*gv1**2*(ga2**2 + 3._dp*gv2**2)*
     .     (2._dp*Mvsq + tt - uu)*(tt + uu) - 3._dp*ga1**2*gv1*gv2*
     .     (3._dp*ga2**2 + gv2**2)*(tt**2 + 3._dp*uu**2 + 2._dp*Mvsq*(tt 
     .     + uu)) - gv1**3*gv2*(3._dp*ga2**2 + gv2**2)*(tt**2 
     .     + 3._dp*uu**2 + 2._dp*Mvsq*(tt + uu)))*xI3(ss,zip,zip,Mvsq,
     .     Mvsq,zip,musq,ep) - (ga1**3*ga2*(ga2**2 + 3._dp*gv2**2)*
     .     (4._dp*Mvsq*uu**2 + 2._dp*Mvsq**2*(tt + 2._dp*uu) + uu*(tt**2 
     .     - uu**2)) + 3._dp*ga1*ga2*gv1**2*(ga2**2 + 3._dp*gv2**2)*
     .     (4._dp*Mvsq*uu**2 + 2._dp*Mvsq**2*(tt + 2._dp*uu) + uu*(tt**2 
     .     - uu**2)) - 3._dp*ga1**2*gv1*gv2*(3._dp*ga2**2 + gv2**2)*
     .     (4._dp*Mvsq*uu**2 + 2._dp*Mvsq**2*(tt + 2._dp*uu) + uu*(tt**2 
     .     + 3._dp*uu**2)) - gv1**3*gv2*(3._dp*ga2**2 + gv2**2)*
     .     (4._dp*Mvsq*uu**2 + 2._dp*Mvsq**2*(tt + 2._dp*uu) + uu*(tt**2 
     ,     + 3._dp*uu**2)))*xI4(zip,zip,zip,zip,ss,uu,Mvsq,zip,
     .     Mvsq,zip,musq,ep))

      bxzat = 
     .     + 8._dp*Ql*Qq*(-2._dp*(4._dp*ga1*ga2*gv1*gv2 + ga1**2*(ga2**2 
     .     + gv2**2) + gv1**2*(ga2**2 + gv2**2))*uu*xI2(ss,Mvsq,Mvsq,
     .     musq,ep) + 2._dp*(4._dp*ga1*ga2*gv1*gv2 + ga1**2*(ga2**2 
     .     + gv2**2) + gv1**2*(ga2**2 + gv2**2))*uu*xI2(tt,zip,zip,
     .     musq,ep) + (4._dp*ga1*ga2*gv1*gv2*(2._dp*Mvsq - tt + uu)*(tt 
     .     + uu) + ga1**2*(ga2**2 + gv2**2)*(3._dp*tt**2 + uu**2 
     .     + 2._dp*Mvsq*(tt + uu)) + gv1**2*(ga2**2 + gv2**2)*
     .     (3._dp*tt**2 + uu**2 + 2._dp*Mvsq*(tt + uu)))*xI3(zip,zip,
     .     ss,Mvsq,zip,Mvsq,musq,ep) - (4._dp*ga1*ga2*gv1*gv2 
     .     + ga1**2*(ga2**2 + gv2**2) + gv1**2*(ga2**2 + gv2**2))*tt*
     .     (2._dp*Mvsq + tt - uu)*xI3(zip,zip,tt,zip,Mvsq,zip,
     .     musq,ep) - (4._dp*ga1*ga2*gv1*gv2 + ga1**2*(ga2**2 + gv2**2) 
     .     + gv1**2*(ga2**2 + gv2**2))*tt*(2._dp*Mvsq + tt - uu)*
     .     xI3(zip,tt,zip,Mvsq,zip,zip,musq,ep) 
     .     + (4._dp*ga1*ga2*gv1*gv2*(2._dp*Mvsq - tt + uu)*(tt + uu) 
     .     + ga1**2*(ga2**2 + gv2**2)*(3._dp*tt**2 + uu**2 
     .     + 2._dp*Mvsq*(tt + uu)) + gv1**2*(ga2**2 + gv2**2)*
     .     (3._dp*tt**2 + uu**2 + 2._dp*Mvsq*(tt + uu)))*
     .     xI3(ss,zip,zip,Mvsq,Mvsq,zip,musq,ep) 
     .     - (4._dp*ga1*ga2*gv1*gv2*(4._dp*Mvsq*tt**2 - tt**3 + tt*uu**2 
     .     + 2._dp*Mvsq**2*(2._dp*tt + uu)) + ga1**2*(ga2**2 + gv2**2)*
     .     (4._dp*Mvsq*tt**2 + 2._dp*Mvsq**2*(2._dp*tt + uu) 
     .     + tt*(3._dp*tt**2 + uu**2)) + gv1**2*(ga2**2 + gv2**2)*
     .     (4._dp*Mvsq*tt**2 + 2._dp*Mvsq**2*(2._dp*tt + uu) 
     .     + tt*(3._dp*tt**2 + uu**2)))*xI4(zip,zip,zip,zip,ss,tt,
     .     Mvsq,zip,Mvsq,zip,musq,ep))

      bxzau = 
     .     + 8._dp*Ql*Qq*(2._dp*(-4._dp*ga1*ga2*gv1*gv2 + ga1**2*(ga2**2 
     .     + gv2**2) + gv1**2*(ga2**2 + gv2**2))*tt*xI2(ss,Mvsq,Mvsq,
     .     musq,ep) - 2._dp*(-4._dp*ga1*ga2*gv1*gv2 + ga1**2*(ga2**2 
     .     + gv2**2) + gv1**2*(ga2**2 + gv2**2))*tt*xI2(uu,zip,zip,
     .     musq,ep) - (-4._dp*ga1*ga2*gv1*gv2*(2._dp*Mvsq + tt - uu)*
     .     (tt + uu) + ga1**2*(ga2**2 + gv2**2)*(tt**2 + 3._dp*uu**2 
     .     + 2._dp*Mvsq*(tt + uu)) + gv1**2*(ga2**2 + gv2**2)*(tt**2 
     .     + 3._dp*uu**2 + 2._dp*Mvsq*(tt + uu)))*xI3(zip,zip,ss,
     .     Mvsq,zip,Mvsq,musq,ep) + (-4._dp*ga1*ga2*gv1*gv2 
     .     + ga1**2*(ga2**2 + gv2**2) + gv1**2*(ga2**2 + gv2**2))*uu*
     .     (2._dp*Mvsq - tt + uu)*xI3(zip,zip,uu,zip,Mvsq,zip,
     .     musq,ep) + (-4._dp*ga1*ga2*gv1*gv2 + ga1**2*(ga2**2 + gv2**2) 
     .     + gv1**2*(ga2**2 + gv2**2))*uu*(2._dp*Mvsq - tt + uu)*
     .     xI3(zip,uu,zip,Mvsq,zip,zip,musq,ep) 
     .     - (-4._dp*ga1*ga2*gv1*gv2*(2._dp*Mvsq + tt - uu)*(tt + uu) 
     .     + ga1**2*(ga2**2 + gv2**2)*(tt**2 + 3._dp*uu**2 
     .     + 2._dp*Mvsq*(tt + uu)) + gv1**2*(ga2**2 + gv2**2)*(tt**2 
     .     + 3._dp*uu**2 + 2._dp*Mvsq*(tt + uu)))*xI3(ss,zip,zip,
     .     Mvsq,Mvsq,zip,musq,ep) + (-4._dp*ga1*ga2*gv1*gv2*
     .     (4._dp*Mvsq*uu**2 + 2._dp*Mvsq**2*(tt + 2._dp*uu) + uu*(tt**2 
     .     - uu**2)) + ga1**2*(ga2**2 + gv2**2)*(4._dp*Mvsq*uu**2 
     .     + 2._dp*Mvsq**2*(tt + 2._dp*uu) + uu*(tt**2 + 3._dp*uu**2)) 
     .     + gv1**2*(ga2**2 + gv2**2)*(4._dp*Mvsq*uu**2 + 2._dp*Mvsq**2*(tt 
     .     + 2._dp*uu) + uu*(tt**2 + 3._dp*uu**2)))*
     .     xI4(zip,zip,zip,zip,ss,uu,Mvsq,zip,Mvsq,zip,musq,ep))

C -- for the use of comparison individually
      corr(:,1) = fac*propz*bxzt(:)
      corr(:,2) = fac*propz*bxzu(:)
      corr(:,3) = fac*propz*bxwt(:)
      corr(:,4) = fac*propz*bxwu(:)
      corr(:,5) = fac*propa*bxzat(:)
      corr(:,6) = fac*propa*bxzau(:)
      corr(:,7) = fac*propa*bxwat(:)
      corr(:,8) = fac*propa*bxwau(:)

      end subroutine z_ew_box
