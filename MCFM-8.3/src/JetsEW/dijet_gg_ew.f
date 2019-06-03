C --- virtual weak correction to gluon fusion to qqb
C --- matrix element w/ input arguments ss,tt,uu
      subroutine dijet_gg_ew(msq,ss,tt,uu)
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'constants.f'
      include 'qcdcouple.f'
      real(dp):: msq(nf),msq1(nf),msq2(nf),ss,tt,uu,z

C --- constraint ss + tt + uu = 0
      z = -2._dp*tt/ss - 1._dp
      call twoj_gg_ew(msq1,ss,z)
      call twoj_gg_ew(msq2,ss,-z)

      msq(1:nf) = msq1(1:nf) + msq2(1:nf)
C --- take off the initial state average
      msq = msq/avegg
c     msq = 0._dp
      
      return
      end





      subroutine twoj_gg_ew(ggqqb,s,z)
C --- The analytical result is taken from the literature
C --- arXiv:0909.0059 by Kuehn et al
C --- Published in Phys.Rev. D82 (2010) 013007
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'constants.f'
      include 'masses.f'
      include 'scale.f'
      include 'ewcharge.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'cplx.h'
      real(dp):: ggqqb(nf),z,s,t,u,mz,mw,y1,y2,xx,yy,f,
     .     xI1,xI2,xI3,xI4,ddilog,sz,sw,spX,vz,vw,vp,bz,bw,bp,
     .     gvq(nf),gaq(nf),T3(nf),facbox,facself,sw2,cw2,
     .     sw5,sp5,bw5,bp5,facbx
      integer ep

C --- finite result so that no coefficients of poles need be checked (zeros)
C --- namely, it is sufficient to consider only ep = 0
      ep = 0
      mz = zmass
      mw = wmass

C --- mixing angles, for on-shell, use sw2 = 1 - mw^2/mz^2
      sw2 = xw
      cw2 = 1._dp - sw2

C --- 3rd isospin to the up- and down-type quarks
      T3 = (/-0.5_dp,0.5_dp,-0.5_dp,0.5_dp,-0.5_dp/)

C --- vector and axial Z couplings to quarks
      gvq(1:5) = (T3(1:5)-2._dp*sw2*Q(1:5))/2._dp/sqrt(sw2*cw2)
      gaq(1:5) = T3(1:5)/2._dp/sqrt(sw2*cw2)
C --- W coulpling to quarks
      gw = 0.5_dp/sqrt(2._dp*sw2)

C --- Mandelstam variables t and u, written in terms of s and z = cos(theta)
C --- here they use '1+z' to t and '1-z' to u; while in the ttb calculation, 
C --- they adoptted the opposite definition. see in 'qqb_ew.f'. 
      t = -s/2._dp*(1._dp+z)
      u = -s/2._dp*(1._dp-z)

C --- switch to ss,tt,uu language since it is easier to use crossing relation
C --- in the other routine
c      z = -2._dp*t/s - 1._dp
      
C --- y1 and y2 enter function f(y1,y2) see (III.37) and (III.38)
C --- only for the case envolves final state b-quark
C --- should be carefull to deal with yy when crossing rel t -> s
c      if(.false.) then
      xx = t+mw**2-mt**2
      yy = real(sqrt(cplx1(1._dp - 4._dp*mw**2*t/xx**2)),dp)
      y1 = xx/2._dp/t*(1._dp + yy)
      y2 = xx/2._dp/t*(1._dp - yy)
c      end if

C --- initializing each msq
      sz = 0._dp
      sw = 0._dp
      spX = 0._dp
      vz = 0._dp
      vw = 0._dp
      vp = 0._dp
      bz = 0._dp
      bw = 0._dp
      bp = 0._dp

      facbox = (2._dp-Nc**2+Nc**2*z)/(Nc*(-1._dp+Nc**2)*s*(1._dp+z))
      facbox = facbox*as**2*esq/4._dp
      facself = facbox*(1._dp+z**2)/(1._dp+z)

C --- for final state quarks u,d,s,c, only Z and W corrections are non-zero
      sz = facself*(
     .     + 4._dp*mz**2 + 3._dp*s*(1._dp + z) - (2._dp*(2._dp*mz**2
     .     + s*(1._dp + z))**2*real(log(cplx1(1._dp + (0.5_dp*s*(1._dp 
     .     + z))/mz**2))))/(s*(1._dp + z))
     .     )

C --- W correction is similar with Z correction. Results of sw and bw are 
C --- obtained by doing mz -> mw; (gvq^2 + gaq^2) -> 2*gw^2 in sz and bz, 
C --- respectively, where '2' is the factor in front of 'facself' and 'facbox'. 
      sw = 2._dp*facself*(
     .     + 4._dp*mw**2 + 3._dp*s*(1._dp + z) - (2._dp*(2._dp*mw**2
     .     + s*(1._dp + z))**2*real(log(cplx1(1._dp + (0.5_dp*s*(1._dp 
     .     + z))/mw**2))))/(s*(1._dp + z))
     .     )

      bz = facbox*(
     .     real(
     .     + 2._dp*(4._dp*(mz**2 + s) 
     .     + (4._dp*mz**2 + s)*z + (2._dp*mz**2 + 5._dp*s)*z**2 
     .     - (1._dp + z)*(-2._dp*s + (2._dp*mz**2 + 3._dp*s)*z)*log( 
     .     cplx1(s/mz**2)) - (2._dp*(mz**2 + s)**2*(1._dp + z)*(4._dp 
     .     - z + z**2)*(ddilog((-s)/mz**2) 
     .     + log(cplx1(s/mz**2))*log(cplx1(1._dp + s/mz**2))))/
     .     (s*(1._dp - z)) - ((4._dp*mz**2 + 6._dp*s - (-2._dp*mz**2 + s)*z 
     .     + s*z**2)*(2._dp*mz**2 + s*(1._dp + z))*
     .     log(cplx1(1._dp - (t)/mz**2)))/(s*(1._dp + z)) 
     .     + (2._dp*(8._dp*mz**4 + 12._dp*mz**2*s + 5._dp*s**2 
     .     + 2._dp*s*(2._dp*mz**2 + s)*z + s**2*z**2)*(ddilog(t/mz**2)
     .     + log(cplx1(s/mz**2))*log(cplx1(1._dp - (t)/mz**2))))/
     .     (s*(1._dp - z)))
     .     ,dp)
     .     )

      bw = 2._dp*facbox*(
     .     real(
     .     + 2._dp*(4._dp*(mw**2 + s) 
     .     + (4._dp*mw**2 + s)*z + (2._dp*mw**2 + 5._dp*s)*z**2 
     .     - (1._dp + z)*(-2._dp*s + (2._dp*mw**2 + 3._dp*s)*z)*log( 
     .     cplx1(s/mw**2)) - (2._dp*(mw**2 + s)**2*(1._dp + z)*(4._dp 
     .     - z + z**2)*(ddilog((-s)/mw**2) 
     .     + log(cplx1(s/mw**2))*log(cplx1(1._dp + s/mw**2))))/
     .     (s*(1._dp - z)) - ((4._dp*mw**2 + 6._dp*s - (-2._dp*mw**2 + s)*z 
     .     + s*z**2)*(2._dp*mw**2 + s*(1._dp + z))*
     .     log(cplx1(1._dp - (t)/mw**2)))/(s*(1._dp + z)) 
     .     + (2._dp*(8._dp*mw**4 + 12._dp*mw**2*s + 5._dp*s**2 
     .     + 2._dp*s*(2._dp*mw**2 + s)*z + s**2*z**2)*(ddilog(t/mw**2)
     .     + log(cplx1(s/mw**2))*log(cplx1(1._dp - (t)/mw**2))))/
     .     (s*(1._dp - z)))
     .     ,dp)
     .     )


c--- BEGIN versions from Jia using scalar integrals (robust under crossing)
c--- BEGIN versions from Jia using scalar integrals (robust under crossing)
c--- BEGIN versions from Jia using scalar integrals (robust under crossing)

      facbx = 2._dp*(2._dp+Nc**2*(-1._dp+z))/Nc/(Nc**2-1)**2*as**2*esq
      bz = facbx*(
     .     + 2._dp*z - (2._dp*(1._dp + z + 2._dp*z**2)*xI1(mz**2,musq,ep))/
     .     (s*(1._dp + z)) - (2._dp*(3._dp*mz**2*(1._dp + z)
     .     + 4._dp*s*(1._dp + z**2))*xI2(0._dp,mz**2,0._dp,musq,ep))/
     .     (s*(1._dp + z)) + (2._dp*(2._dp*mz**2*z + s*(-2._dp
     .     + 3._dp*z))*xI2(s,0._dp,0._dp,musq,ep))/s
     .     + (2._dp*(2._dp*mz**2*(2._dp + z) + s*(6._dp - 1._dp*z
     .     + z**2))*xI2(t,mz**2,0._dp,musq,ep))/(s*(1._dp + z))
     .     + (2._dp*(8._dp*mz**4 + 4._dp*mz**2*s*(3._dp + z) + s**2*(5._dp
     .     + 2._dp*z + z**2))*xI3(0._dp,0._dp,s,0._dp,0._dp,0._dp,musq,ep))/
     .     (s*(-1._dp + z)*(1._dp + z)) - (1._dp*(8._dp*mz**4
     .     + 4._dp*mz**2*s*(3._dp + z) + s**2*(5._dp + 2._dp*z
     .     + z**2))*xI3(0._dp,0._dp,t,mz**2,0._dp,0._dp,musq,ep))/
     .     (s*(-1._dp + z)) + (2._dp*(2._dp*mz**4*z*(3._dp + z**2)
     .     + 4._dp*mz**2*s*(1._dp + 2._dp*z + z**3) + s**2*(3._dp + 4._dp*z
     .     - 1._dp*z**2 + 2._dp*z**3))*xI3(0._dp,s,0._dp,mz**2,0._dp,0._dp,
     .     musq,ep))/(s*(-1._dp + z)*(1._dp + z)) - (1._dp*(8._dp*mz**4
     .     + 4._dp*mz**2*s*(3._dp + z) + s**2*(5._dp + 2._dp*z
     .     + z**2))*xI3(t,0._dp,0._dp,mz**2,0._dp,0._dp,musq,ep))/
     .     (s*(-1._dp + z)) + ((2._dp*mz**2 + s + s*z)*(8._dp*mz**4
     .     + 4._dp*mz**2*s*(3._dp + z) + s**2*(5._dp + 2._dp*z
     .     + z**2))*xI4(0._dp,0._dp,0._dp,0._dp,t,s,mz**2,0._dp,0._dp,0._dp,
     .     musq,ep))/(s*(-1._dp + z)*(1._dp + z))
     .     ) 

      bw = 2._dp*facbx*(
     .     + 2._dp*z - (2._dp*(1._dp + z + 2._dp*z**2)*xI1(mw**2,musq,ep))/
     .     (s*(1._dp + z)) - (2._dp*(3._dp*mw**2*(1._dp + z)
     .     + 4._dp*s*(1._dp + z**2))*xI2(0._dp,mw**2,0._dp,musq,ep))/
     .     (s*(1._dp + z)) + (2._dp*(2._dp*mw**2*z + s*(-2._dp
     .     + 3._dp*z))*xI2(s,0._dp,0._dp,musq,ep))/s
     .     + (2._dp*(2._dp*mw**2*(2._dp + z) + s*(6._dp - 1._dp*z
     .     + z**2))*xI2(t,mw**2,0._dp,musq,ep))/(s*(1._dp + z))
     .     + (2._dp*(8._dp*mw**4 + 4._dp*mw**2*s*(3._dp + z) + s**2*(5._dp
     .     + 2._dp*z + z**2))*xI3(0._dp,0._dp,s,0._dp,0._dp,0._dp,musq,ep))/
     .     (s*(-1._dp + z)*(1._dp + z)) - (1._dp*(8._dp*mw**4
     .     + 4._dp*mw**2*s*(3._dp + z) + s**2*(5._dp + 2._dp*z
     .     + z**2))*xI3(0._dp,0._dp,t,mw**2,0._dp,0._dp,musq,ep))/
     .     (s*(-1._dp + z)) + (2._dp*(2._dp*mw**4*z*(3._dp + z**2)
     .     + 4._dp*mw**2*s*(1._dp + 2._dp*z + z**3) + s**2*(3._dp + 4._dp*z
     .     - 1._dp*z**2 + 2._dp*z**3))*xI3(0._dp,s,0._dp,mw**2,0._dp,0._dp,
     .     musq,ep))/(s*(-1._dp + z)*(1._dp + z)) - (1._dp*(8._dp*mw**4
     .     + 4._dp*mw**2*s*(3._dp + z) + s**2*(5._dp + 2._dp*z
     .     + z**2))*xI3(t,0._dp,0._dp,mw**2,0._dp,0._dp,musq,ep))/
     .     (s*(-1._dp + z)) + ((2._dp*mw**2 + s + s*z)*(8._dp*mw**4
     .     + 4._dp*mw**2*s*(3._dp + z) + s**2*(5._dp + 2._dp*z
     .     + z**2))*xI4(0._dp,0._dp,0._dp,0._dp,t,s,mw**2,0._dp,0._dp,0._dp,
     .     musq,ep))/(s*(-1._dp + z)*(1._dp + z))
     .     ) 
     
c--- END versions from Jia using scalar integrals (robust under crossing)
c--- END versions from Jia using scalar integrals (robust under crossing)
c--- END versions from Jia using scalar integrals (robust under crossing)


C --- when it is a b-quark

      sw5 = facself*(
     .     + 4._dp*((0.5_dp*(-4._dp*(mt**2 - mw**2)**2 + (mt**2 
     .     - 3._dp*mw**2)*s*(1._dp + z)))/(mt**2 - mw**2) + (-2._dp*mt**2 
     .     + 2._dp*mw**2 + s*(1._dp + z))*f(y1,y2) + (mw**2*(-2._dp*(mt**2 
     .     - mw**2)**2 + (2._dp*mt**2 - mw**2)*s*(1._dp + z))*
     .     real(log(cplx1(mt**2/mw**2)),dp))/(mt**2 - mw**2)**2)
     .     )

      bw5 = facbox*(
     .     + 4._dp*(s*z*(1._dp + z) + (2._dp*(2._dp*(-mt**2 + mw**2 
     .     + s) - 2._dp*(mt**2 - mw**2)*z + (-mt**2 + mw**2 
     .     + 2._dp*s)*z**2)*(-xI1(mt**2,musq,ep) + xI1(mw**2,musq,ep)))
     .     /(mt**2 - mw**2) - (1._dp + z)*(2._dp*s - (-2._dp*mt**2 
     .     + 2._dp*mw**2 + 3._dp*s)*z)*xI2(s,mt**2,mt**2,musq,ep) 
     .     + (-4._dp*mt**2 + 4._dp*mw**2 + 6._dp*s - (2._dp*mt**2 
     .     - 2._dp*mw**2 + s)*z + s*z**2)*xI2(t,mt**2,mw**2,musq,ep) 
     .     - ((s*(-4._dp*mt**2 + 4._dp*mw**2 + 3._dp*s) + (6._dp*(mt**2 
     .     - mw**2)**2 + 2._dp*s*(-3._dp*mt**2 + 4._dp*mw**2 + 2._dp*s))*z 
     .     - s**2*z**2 + 2._dp*((mt**2 - mw**2)**2 + s*(-mt**2 
     .     + 2._dp*mw**2 + s))*z**3)*xI3(0._dp,0._dp,s,mt**2,mw**2,
     .     mt**2,musq,ep))/(1._dp - z) - ((8._dp*(mt**2 - mw**2)**2 
     .     + 4._dp*(-2._dp*mt**2 + 3._dp*mw**2)*s + 5._dp*s**2 
     .     + 2._dp*s*(-2._dp*mt**2 + 2._dp*mw**2 + s)*z + s**2*z**2)*
     .     (xI3(0._dp,0._dp,s,mt**2,mt**2,mt**2,musq,ep) - 
     .     (1._dp + z)*xI3(0._dp,0._dp,t,mt**2,mt**2,mw**2,musq,ep)))
     .     /(1._dp - z) - (0.5_dp*(-16._dp*(mt**2 - mw**2)**3
     .     + 16._dp*(mt**4 - 3._dp*mt**2*mw**2 + 2._dp*mw**4)*s
     .     - 10._dp*mt**2*s**2 + 22._dp*mw**2*s**2
     .     + 5._dp*s**3 + s*(16._dp*(mt**2 - mw**2)**2 
     .     - 16._dp*mt**2*s + 20._dp*mw**2*s + 7._dp*s**2)*z 
     .     + s*(8._dp*mt**2*(mt**2 - mw**2) - 14._dp*mt**2*s 
     .     + 6._dp*mw**2*s + 3._dp*s**2)*z**2 + s**3*z**3)*
     .     xI4(0._dp,0._dp,0._dp,s + t + u,t,s,mt**2,mt**2,mw**2,
     .     mt**2,musq,ep))/(1._dp - z))
     .     )

         bp5 = facbox*(
     .     (2._dp*mt**2*(s*z*(1._dp + z) + 2._dp*(z**2 + 2._dp*(1._dp 
     .     + z))*(xI1(mt**2,musq,ep) - xI1(mw**2,musq,ep)) 
     .     - (1._dp + z)*(2._dp*s + (2._dp*mt**2 - 2._dp*mw**2 + s)*z)*
     .     xI2(s,mt**2,mt**2,musq,ep) + (2._dp + z)*(-2._dp*(mt**2 
     .     - mw**2) + s*(1._dp + z))*xI2(t,mt**2,mw**2,musq,ep) 
     .     - ((s*(-4._dp*(mt**2 - mw**2) + s) + 2._dp*(3._dp*(mt**2 
     .     - mw**2)**2 - mt**2*s + s**2 + 2._dp*s*mw**2)*z + s**2*z**2 
     .     + 2._dp*((mt**2 - mw**2)**2 + mt**2*s)*z**3)*
     .     xI3(0._dp,0._dp,s,mt**2,mw**2,mt**2,musq,ep))/
     .     (1._dp - z) - ((8._dp*(mt**2 
     .     - mw**2)**2 + 4._dp*mw**2*s + s**2 + 2._dp*s*(-2._dp*mt**2 
     .     + 2._dp*mw**2 + s)*z + s**2*z**2)*(xI3(0._dp,0._dp,s,mt**2,
     .     mt**2,mt**2,musq,ep) - (1._dp + z)*
     .     xI3(0._dp,0._dp,t,mt**2,mt**2,mw**2,musq,ep)))/(1._dp - z) 
     .     - (0.5_dp*(16._dp*(-mt**2 
     .     + mw**2)**3 + 16._dp*mw**2*(-mt**2 + mw**2)*s 
     .     - 2._dp*mt**2*s**2 + 6._dp*mw**2*s**2 + s**3  
     .     + s*(16._dp*(mt**2 - mw**2)**2 - 8._dp*mt**2*s 
     .     + 12._dp*mw**2*s + 3._dp*s**2)*z 
     .     + s*(8._dp*mt**2*(mt**2 - mw**2) + 3._dp*s*(-2._dp*mt**2 
     .     + 2._dp*mw**2 + s))*z**2 + s**3*z**3)*xI4(0._dp,0._dp,0._dp,
     .     s + t + u,t,s,mt**2,mt**2,mw**2,mt**2,musq,ep))/
     .     (1._dp - z)))/mw**2
     .     )

C --- see (III.39), (III.40)
      spX = 0._dp
      vz = -2._dp*sz
      vw = -2._dp*sw
      vp = -2._dp*spX

      ggqqb(1:4) = (sz+vz+bz)*(gvq(1:4)**2+gaq(1:4)**2)
     .     +gw**2*(sw+spX+vw+vp+bw+bp)


      sp5 = mt**2/2._dp/mw**2*sw5 ! for final states b-quark only
      vw = -2._dp*sw5
      vp = -2._dp*sp5
      ggqqb(5) = (sz+vz+bz)*(gvq(5)**2+gaq(5)**2)
     .     +gw**2*(sw5+sp5+vw+vp+bw5+bp5)
C --- no b final state contribution
      ggqqb(5) = 0._dp

      end subroutine twoj_gg_ew


C --- function specified from two point function B0, see (III.37)
      function f(y1,y2)
      implicit none
      include 'types.f'
      include 'cplx.h'
      real(dp):: y1,y2,f

      f = real(- y1*log(cplx1(y1/(y1-1._dp))) 
     .     - y2*log(cplx1(y2/(y2-1._dp))),dp)

      end function f
