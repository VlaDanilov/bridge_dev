      function amp_wqqagg(st,i1,i2,i3,i4,i5,i6,i7,xa,xb)
      implicit none
      include 'types.f'
      complex(dp):: amp_wqqagg

*******************************************************************
*     R. Mondini, January 2017                                    *
*     0 -> q(p1)+qb(p2)+l(p3)+vb(p4)+gam(p5)+g(p6)+g(p7)          *
*******************************************************************
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      integer:: i1,i2,i3,i4,i5,i6,i7,j1,j2,j3,j4,j5,j6,j7
      integer:: j,k,st
      real(dp):: s345,s134,s267,s167,s256,s234,s157
      complex(dp):: xa(mxpart,mxpart),xb(mxpart,mxpart),t1ab,t2ab,t11aa,t11bb,t22aa,t22bb,t31aa,t32aa,t23bb
c---- statement functions
c---- <j1|j2|j3]
      t1ab(j1,j2,j3)=xa(j1,j2)*xb(j2,j3)
c---- <j1|j2+j3|j4]
      t2ab(j1,j2,j3,j4)=xa(j1,j2)*xb(j2,j4)+xa(j1,j3)*xb(j3,j4)
c---- <j1|j2|j3|j4>
      t11aa(j1,j2,j3,j4)=xa(j1,j2)*xb(j2,j3)*xa(j3,j4)
c---- [j1|j2|j3|j4]
      t11bb(j1,j2,j3,j4)=xb(j1,j2)*xa(j2,j3)*xb(j3,j4)
c---- <j1|j2+j3|j4+j5|j6>
      t22aa(j1,j2,j3,j4,j5,j6)=t11aa(j1,j2,j4,j6)+t11aa(j1,j2,j5,j6)
     & +t11aa(j1,j3,j4,j6)+t11aa(j1,j3,j5,j6)
c---- [j1|j2+j3|j4+j5|j6]
      t22bb(j1,j2,j3,j4,j5,j6)=t11bb(j1,j2,j4,j6)+t11bb(j1,j2,j5,j6)
     & +t11bb(j1,j3,j4,j6)+t11bb(j1,j3,j5,j6)
c---- <j1|j2+j3+j4|j5|j6>
      t31aa(j1,j2,j3,j4,j5,j6)=t11aa(j1,j2,j5,j6)+t11aa(j1,j3,j5,j6)
     & +t11aa(j1,j4,j5,j6)
c---- <j1|j2+j3+j4|j5+j6|j7>
      t32aa(j1,j2,j3,j4,j5,j6,j7)=t11aa(j1,j2,j5,j7)+t11aa(j1,j3,j5,j7)
     & +t11aa(j1,j4,j5,j7)+t11aa(j1,j2,j6,j7)+t11aa(j1,j3,j6,j7)
     & +t11aa(j1,j4,j6,j7)
c---- [j1|j2+j3|j4+j5+j6|j7]
      t23bb(j1,j2,j3,j4,j5,j6,j7)=t11bb(j1,j2,j4,j7)+t11bb(j1,j2,j5,j7)
     & +t11bb(j1,j2,j6,j7)+t11bb(j1,j3,j4,j7)+t11bb(j1,j3,j5,j7)
     & +t11bb(j1,j3,j6,j7)

      s134=s(i1,i3)+s(i1,i4)+s(i3,i4)
      s234=s(i2,i3)+s(i2,i4)+s(i3,i4)
      s345=s(i3,i4)+s(i3,i5)+s(i4,i5)
      s267=s(i2,i6)+s(i2,i7)+s(i6,i7)
      s167=s(i1,i6)+s(i1,i7)+s(i6,i7)
      s157=s(i1,i5)+s(i1,i7)+s(i5,i7)
      s256=s(i2,i5)+s(i2,i6)+s(i5,i6)

c---- helicity stamps (st)
c---- Qd = photon on leg 2 = [a] in Dixon,Signer [hep-ph/9803250] (e.g. see (4.9)-(4.12) for NLO)
c---- Qu = photon on leg 1 = [b] in Dixon,Signer [hep-ph/9803250] (e.g. see (4.9)-(4.12) for NLO)
c---- Ql = photon emitted by lepton (amplitudes include part of the Wga vertex diagrams)

c---- 'q-,qb+,l-,vb+,ga+,g+,g+,Qd'=1
c---- 'q-,qb+,l-,vb+,ga+,g+,g-,Qd'=2
c---- 'q-,qb+,l-,vb+,ga+,g-,g+,Qd'=3
c---- 'q-,qb+,l-,vb+,ga-,g+,g+,Qd'=4

c---- 'q-,qb+,l-,vb+,ga+,g+,g+,Qu'=5
c---- 'q-,qb+,l-,vb+,ga+,g+,g-,Qu'=6
c---- 'q-,qb+,l-,vb+,ga+,g-,g+,Qu'=7
c---- 'q-,qb+,l-,vb+,ga-,g+,g+,Qu'=8

c---- 'q-,qb+,l-,vb+,ga+,g+,g+,Ql'=9
c---- 'q-,qb+,l-,vb+,ga+,g+,g-,Ql'=10
c---- 'q-,qb+,l-,vb+,ga+,g-,g+,Ql'=11
c---- 'q-,qb+,l-,vb+,ga-,g+,g+,Ql'=12
c---- 'q-,qb+,l-,vb+,ga-,g-,g-,Ql'=13
c---- 'q-,qb+,l-,vb+,ga-,g-,g+,Ql'=14
c---- 'q-,qb+,l-,vb+,ga-,g+,g-,Ql'=15
c---- 'q-,qb+,l-,vb+,ga+,g-,g-,Ql'=16

c---- A7[a](i1q-,i2qb+,i3l-,i4vb+,i5ga+,i6g+,i7g+)
      if (st==1) then 
      amp_wqqagg = (xa(i1,i3)**2*t2ab(i2,i3,i4,i5))/(xa(i1,i7)*
     & xa(i2,i5)*xa(i2,i6)*xa(i3,i4)*xa(i6,i7)*t2ab(i5,i3,i4,i5))

c---- A7[a](i1q-,i2qb+,i3l-,i4vb+,i5ga+,i6g+,i7g-)
      elseif (st==2) then
      amp_wqqagg = -((xa(i1,i3)*xa(i2,i7)*t2ab(i7,i1,i3,i4)*
     & t2ab(i7,i2,i6,i5))/(s134*s267*xa(i2,i5)*xa(i2,i6)*xa(i3,i4)*
     & xb(i4,i3)*xa(i6,i7))-(xa(i2,i7)*t32aa(i7,i2,i6,i7,i4,i5,i3)*
     & (xa(i1,i3)*xb(i5,i1)*xa(i2,i7)*xb(i4,i3)-xa(i1,i7)*xb(i5,i4)*
     & t2ab(i2,i6,i7,i1)+xb(i4,i3)*xa(i3,i7)*t2ab(i2,i6,i7,i5)))/(s267*
     & xa(i2,i5)*xa(i2,i6)*xa(i3,i4)*xb(i4,i3)*xa(i6,i7)*t2ab(i5,i3,i4,i5)*
     & t2ab(i2,i6,i7,i1)))-(xb(i6,i1)*t2ab(i3,i1,i7,i6)**2*
     & t2ab(i2,i3,i4,i5))/(s167*xb(i7,i1)*xa(i2,i5)*xa(i3,i4)*xb(i7,i6)*
     & t2ab(i5,i3,i4,i5)*t2ab(i2,i6,i7,i1))

c---- A7[a](i1q-,i2qb+,i3l-,i4vb+,i5ga+,i6g-,i7g+)
      elseif (st==3) then
      amp_wqqagg = (xa(i1,i3)**2*t2ab(i6,i2,i5,i7)**3)/(s134*s256*
     & xa(i2,i5)*xa(i3,i4)*t2ab(i5,i2,i6,i7)*t32aa(i6,i2,i5,i6,i3,i4,i1))+
     & (xa(i1,i3)**2*xb(i7,i2)**3*t22bb(i5,i3,i4,i2,i6,i7))/(s267*
     & xb(i6,i2)*xa(i3,i4)*xb(i7,i6)*t2ab(i5,i3,i4,i5)*t2ab(i1,i6,i7,i2)*
     & t2ab(i5,i2,i6,i7))+((xa(i1,i6)**3*t2ab(i3,i4,i5,i2)*
     & t23bb(i4,i2,i5,i1,i6,i7,i5))/(s167*xa(i1,i7)*xa(i2,i5)*xa(i3,i4)*
     & xb(i4,i3)*xa(i6,i7)*t2ab(i5,i3,i4,i5)*t2ab(i1,i6,i7,i2))-(xa(i1,i3)*
     & xa(i1,i6)**3*xb(i5,i2)*t2ab(i6,i2,i5,i4))/(xa(i1,i7)*xa(i2,i5)*
     & xa(i3,i4)*xb(i4,i3)*xa(i6,i7)*t2ab(i1,i6,i7,i2)*
     & t22aa(i1,i3,i4,i2,i5,i6)))

c---- A7[a](i1q-,i2qb+,i3l-,i4vb+,i5ga-,i6g+,i7g+)
      elseif (st==4) then
      amp_wqqagg = -((s267*xa(i1,i3)*t2ab(i5,i1,i3,i4))/(s134*xa(i2,i6)*
     & xa(i3,i4)*xb(i4,i3)*xa(i6,i7)*t2ab(i7,i2,i6,i5))-(t2ab(i1,i3,i5,i4)*
     & (xa(i1,i5)*t22aa(i7,i2,i6,i1,i5,i3)-s267*xa(i1,i3)*xa(i5,i7)))/
     & (xa(i1,i7)*xa(i2,i6)*xa(i3,i4)*xb(i4,i3)*xa(i6,i7)*t2ab(i5,i3,i4,i5)*
     & t2ab(i7,i2,i6,i5)))+(xa(i1,i3)**2*xb(i6,i2)**2)/(s256*xa(i1,i7)*
     & xb(i5,i2)*xa(i3,i4)*t2ab(i7,i2,i6,i5))

c---- A7[b](i1q-,i2qb+,i3l-,i4vb+,i5ga+,i6g+,i7g+)
      elseif (st==5) then
      amp_wqqagg = -((xa(i1,i3)**2*t2ab(i1,i3,i4,i5))/(xa(i1,i5)*
     & xa(i1,i7)*xa(i2,i6)*xa(i3,i4)*xa(i6,i7)*t2ab(i5,i3,i4,i5)))

c---- A7[b](i1q-,i2qb+,i3l-,i4vb+,i5ga+,i6g+,i7g-)
      elseif (st==6) then
      amp_wqqagg = -((xb(i6,i1)*t2ab(i3,i1,i7,i6)**2*
     & t22bb(i5,i3,i4,i1,i7,i6))/(s167*xb(i7,i1)*xa(i3,i4)*xb(i7,i6)*
     & t2ab(i5,i3,i4,i5)*t2ab(i2,i6,i7,i1)*t2ab(i5,i1,i7,i6)))
     & -((xa(i2,i7)*(s267*xa(i3,i7)-xa(i1,i3)*t2ab(i7,i2,i6,i1))*
     & (s267*xa(i1,i7)*xb(i5,i4)-xa(i1,i3)*xb(i4,i3)*t2ab(i7,i2,i6,i5)))/
     & (s267*xa(i1,i5)*xa(i2,i6)*xa(i3,i4)*xb(i4,i3)*xa(i6,i7)*
     & t2ab(i5,i3,i4,i5)*t2ab(i2,i6,i7,i1))-(xa(i1,i7)*xa(i2,i7)*
     & t2ab(i7,i2,i6,i4)*(xa(i2,i7)*t2ab(i3,i2,i4,i5)+xa(i2,i3)*
     & t1ab(i7,i6,i5)))/(xa(i1,i5)*xa(i2,i6)*xa(i3,i4)*xb(i4,i3)*
     & xa(i6,i7)*t2ab(i2,i6,i7,i1)*t22aa(i2,i3,i4,i1,i5,i7)))-(xa(i1,i7)**2*
     & t2ab(i7,i1,i5,i6)*t2ab(i3,i2,i4,i6)**2)/(s157*s234*xa(i1,i5)*
     & xa(i3,i4)*t2ab(i5,i1,i7,i6)*t22aa(i2,i3,i4,i1,i5,i7))

c---- A7[b](i1q-,i2qb+,i3l-,i4vb+,i5ga+,i6g-,i7g+)
      elseif (st==7) then
      amp_wqqagg = ((xa(i1,i6)**3*t2ab(i3,i4,i5,i2)*(xb(i5,i2)*
     & t1ab(i1,i3,i4)+xb(i5,i4)*t2ab(i1,i6,i7,i2)))/(s167*xa(i1,i5)*xa(i1,i7)*
     & xa(i3,i4)*xb(i4,i3)*xa(i6,i7)*t2ab(i5,i3,i4,i5)*t2ab(i1,i6,i7,i2))-
     & (xa(i1,i6)**3*xb(i4,i2)*t2ab(i3,i2,i4,i5))/(s167*s234*xa(i1,i5)*xa(i1,i7)*
     & xa(i3,i4)*xb(i4,i3)*xa(i6,i7)))+(xa(i1,i3)**2*xb(i7,i2)**3*t2ab(i1,i3,i4,i5))/
     & (s267*xa(i1,i5)*xb(i6,i2)*xa(i3,i4)*xb(i7,i6)*t2ab(i5,i3,i4,i5)*t2ab(i1,i6,i7,i2))

c---- A7[b](i1q-,i2qb+,i3l-,i4vb+,i5ga-,i6g+,i7g+)
      elseif (st==8) then
      amp_wqqagg = -((xa(i1,i5)*t2ab(i7,i2,i6,i4)*t22aa(i2,i6,i7,i1,i5,i3))/
     & (xb(i5,i1)*xa(i1,i7)*xa(i2,i6)*xa(i3,i4)*xb(i4,i3)*xa(i6,i7)*
     & t22aa(i2,i3,i4,i1,i5,i7))-(t2ab(i1,i3,i5,i4)*t32aa(i5,i2,i6,i7,i1,i5,i3))/
     & (xb(i5,i1)*xa(i1,i7)*xa(i2,i6)*xa(i3,i4)*xb(i4,i3)*xa(i6,i7)*t2ab(i5,i3,i4,i5)))
     & +(xa(i1,i5)**2*t2ab(i3,i2,i4,i6)**2)/(s157*s234*xa(i1,i7)*xa(i3,i4)*
     & t22aa(i7,i1,i5,i3,i4,i2))

c---- A7[l](i1q-,i2qb+,i3l-,i4vb+,i5ga+,i6g+,i7g+)
      elseif (st==9) then
      amp_wqqagg = (xa(i1,i3)**2*xb(i5,i4))/(xa(i1,i7)*
     & xa(i2,i6)*xa(i3,i5)*xa(i6,i7)*t2ab(i5,i3,i4,i5))

c---- A7[l](i1q-,i2qb+,i3l-,i4vb+,i5ga+,i6g+,i7g-)
      elseif (st==10) then
      amp_wqqagg = -((xb(i6,i1)*xb(i5,i4)*t2ab(i3,i1,i7,i6)**2)/
     & (s167*xb(i7,i1)*xa(i3,i5)*xb(i7,i6)*t2ab(i5,i3,i4,i5)*t2ab(i2,i6,i7,i1)))
     & -(xa(i2,i7)*xb(i5,i4)*t22aa(i7,i2,i6,i4,i5,i3)**2)/(s267*s345*
     & xa(i2,i6)*xa(i3,i5)*xa(i6,i7)*t2ab(i5,i3,i4,i5)*t2ab(i2,i6,i7,i1))

c---- A7[l](i1q-,i2qb+,i3l-,i4vb+,i5ga+,i6g-,i7g+)
      elseif (st==11) then
      amp_wqqagg = -(xa(i1,i6)**3*xb(i5,i4)*t2ab(i3,i4,i5,i2)**2)/
     & (s167*s345*xa(i1,i7)*xa(i3,i5)*xa(i6,i7)*t2ab(i5,i3,i4,i5)*
     & t2ab(i1,i6,i7,i2))-(xa(i1,i3)**2*xb(i7,i2)**3*xb(i5,i4))/
     & (s267*xb(i6,i2)*xa(i3,i5)*xb(i7,i6)*t2ab(i5,i3,i4,i5)*
     & t2ab(i1,i6,i7,i2))

c---- A7[l](i1q-,i2qb+,i3l-,i4vb+,i5ga-,i6g+,i7g+)
      elseif (st==12) then
      amp_wqqagg = -((xa(i4,i5)*t2ab(i1,i3,i5,i4)**2)/
     & (s345*xa(i1,i7)*xa(i2,i6)*xb(i5,i3)*xa(i6,i7)*t2ab(i5,i3,i4,i5)))

c---- A7[l](i1q-,i2qb+,i3l-,i4vb+,i5ga-,i6g-,i7g-)
      elseif (st==13) then
      amp_wqqagg = (xb(i4,i2)**2*xa(i4,i5))/(xb(i7,i1)*
     & xb(i6,i2)*xb(i5,i3)*xb(i7,i6)*t2ab(i5,i3,i4,i5))

c---- A7[l](i1q-,i2qb+,i3l-,i4vb+,i5ga-,i6g-,i7g+)
      elseif (st==14) then
      amp_wqqagg = -(xa(i1,i6)**3*xb(i4,i2)**2*xa(i4,i5))/
     & (s167*xa(i1,i7)*xb(i5,i3)*xa(i6,i7)*t2ab(i5,i3,i4,i5)*
     & t2ab(i1,i6,i7,i2))-(xb(i7,i2)**3*xa(i4,i5)*
     & t2ab(i1,i3,i5,i4)**2)/(s267*s345*xb(i6,i2)*xb(i5,i3)*
     & xb(i7,i6)*t2ab(i5,i3,i4,i5)*t2ab(i1,i6,i7,i2))

c---- A7[l](i1q-,i2qb+,i3l-,i4vb+,i5ga-,i6g+,i7g-)
      elseif (st==15) then
      amp_wqqagg = -(xb(i6,i1)*xa(i4,i5)*t22bb(i4,i3,i5,i1,i7,i6)**2)/
     & (s167*s345*xb(i7,i1)*xb(i5,i3)*xb(i7,i6)*t2ab(i5,i3,i4,i5)*
     & t2ab(i2,i6,i7,i1))-(xa(i2,i7)*xa(i4,i5)*t2ab(i7,i2,i6,i4)**2)/
     & (s267*xa(i2,i6)*xb(i5,i3)*xa(i6,i7)*t2ab(i5,i3,i4,i5)*
     & t2ab(i2,i6,i7,i1))

c---- A7[l](i1q-,i2qb+,i3l-,i4vb+,i5ga+,i6g-,i7g-)
      elseif (st==16) then
      amp_wqqagg = (-xb(i5,i4)*t2ab(i3,i4,i5,i2)**2)/
     & (s345*xb(i7,i1)*xb(i6,i2)*xa(i3,i5)*xb(i7,i6)*t2ab(i5,i3,i4,i5))

      endif

      return 
      end

