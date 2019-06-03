
      ! UV,IR renormalized in MSBAR scheme
      ! using 1004.3653 and 1309.3245
      subroutine quark_formfactor(s12, Fout)
      implicit none
      include 'types.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scale.f'
      include 'constants.f'
      include 'zeta.f'
      include 'nf.f'

      real(dp), intent(in) :: s12
      complex(dp), intent(out) :: Fout(0:2)

      complex(dp) :: F0q, F1q, F2q
      complex(dp) :: Lnrat, xl12
      complex(dp) :: oneloop, twoloop
        
      xl12=Lnrat(-s12,musq)

      oneloop =
     &  + cf * (  - 4 + zeta2/2 + 3._dp/2._dp*xl12 - 1._dp/2.
     &    _dp*xl12**2 )

      twoloop =
     &  + cf*nf * ( 4085._dp/1296._dp + 1._dp/18._dp*zeta3 + 23._dp/216.
     &    _dp*pi**2 + 209._dp/108._dp*xl12 + 1._dp/18._dp*xl12*pi**2 + 
     &    19._dp/36._dp*xl12**2 + 1._dp/18._dp*xl12**3 )
      twoloop = twoloop + cf**2 * ( 255._dp/32._dp - 15._dp/2._dp*zeta3
     &     + 7._dp/8._dp*pi**2 - 83._dp/1440._dp*pi**4 + 45._dp/8._dp*
     &    xl12 - 6*xl12*zeta3 + 3._dp/8._dp*xl12*pi**2 + 25._dp/8._dp*
     &    xl12**2 - 1._dp/24._dp*xl12**2*pi**2 + 3._dp/4._dp*xl12**3 + 
     &    1._dp/8._dp*xl12**4 )
      twoloop = twoloop + ca*cf * (  - 51157._dp/2592._dp + 313._dp/36._
     &    dp*zeta3 - 337._dp/432._dp*pi**2 + 11._dp/180._dp*pi**4 - 
     &    2545._dp/216._dp*xl12 + 13._dp/2._dp*xl12*zeta3 - 11._dp/36._d
     &    p*xl12*pi**2 - 233._dp/72._dp*xl12**2 + 1._dp/12._dp*xl12**2*
     &    pi**2 - 11._dp/36._dp*xl12**3 )

      Fout(0) = 1
      Fout(1) = oneloop
      Fout(2) = twoloop

      end subroutine

      function wzgamma_amp_RR(i5,i6,i1,i3,i2,za,zb,
     &                        c_alpha,c_beta,c_gamma)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'nf.f'

      complex(dp) :: wzgamma_amp_RR
      integer, intent(in) :: i5,i6,i1,i3,i2
      complex(dp), intent(in) :: c_alpha, c_beta, c_gamma
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

      wzgamma_amp_RR = -2*sqrt(2._dp)*(
     &     za(i2,i5)*za(i1,i2)*zb(i1,i6)/za(i1,i3)/za(i2,i3) * c_alpha
     &   - za(i2,i5)*zb(i3,i6)/za(i1,i3) * c_beta
     &   + za(i1,i5)*zb(i1,i3)*zb(i3,i6)/za(i1,i3)/zb(i2,i3) * c_gamma
     &  )

      end function 
      

      subroutine omega_wzgamma(mp, i5,i6,i1,i3,i2, maxLoops,
     &                            c_alphai, c_betai, c_gammai)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f' ! im, v, xnsq, cf
      include 'nf.f' 
      include 'ewcharge.f' ! q
      include 'zeta.f'
      include 'tiny.f'

      real(dp), intent(in) :: mp(mxpart,4)
      integer, intent(in) :: i5,i6,i1,i3,i2
      integer, intent(in) :: maxLoops
      complex(dp), intent(out) :: c_alphai(0:2,4)
      complex(dp), intent(out) :: c_betai(0:2,4)
      complex(dp), intent(out) :: c_gammai(0:2,4)

      real(dp) :: xu,xv,OneMu,OneMv,OneMuMv,uPv
      real(dp) :: s12,s13,s23,s123;

      integer, parameter :: nw=4

      ! for two-dimensional harmonic polylogarithms
      real(dp) :: HZ1,HZ2,HZ3,HZ4,GYZ1,GYZ2,GYZ3,GYZ4
      dimension HZ1(0:1),HZ2(0:1,0:1),HZ3(0:1,0:1,0:1),
     &          HZ4(0:1,0:1,0:1,0:1)
      dimension GYZ1(0:3),GYZ2(0:3,0:3),GYZ3(0:3,0:3,0:3),
     &          GYZ4(0:3,0:3,0:3,0:3)

      real(dp) :: dotvec
      real(dp) :: Li2

      complex(dp) :: w(1101)
      complex(dp) :: alphaG1WL0, alphaG2WL0, alphaG3WL0
      complex(dp) :: alphaG1WL1, alphaG2WL1, alphaG3WL1
      complex(dp) :: alphaG1WL2, alphaG2WL2, alphaG3WL2, alphaG4WL2
      complex(dp) :: betaG1WL0, betaG2WL0, betaG3WL0
      complex(dp) :: betaG1WL1, betaG2WL1, betaG3WL1
      complex(dp) :: betaG1WL2, betaG2WL2, betaG3WL2, betaG4WL2
      complex(dp) :: gammaG1WL0, gammaG2WL0, gammaG3WL0
      complex(dp) :: gammaG1WL1, gammaG2WL1, gammaG3WL1
      complex(dp) :: gammaG1WL2, gammaG2WL2, gammaG3WL2, gammaG4WL2

      s12 = dotvec(mp(i1,:)+mp(i2,:),mp(i1,:)+mp(i2,:))
      s13 = dotvec(mp(i1,:)+mp(i3,:),mp(i1,:)+mp(i3,:))
      s23 = dotvec(mp(i2,:)+mp(i3,:),mp(i2,:)+mp(i3,:))
      s123 = dotvec(mp(i1,:)+mp(i2,:)+mp(i3,:),mp(i1,:)+mp(i2,:)+mp(i3,:))

      xu = - s13/s12
      xv = s123/s12
      OneMu = 1-xu
      OneMv = 1-xv
      OneMuMv = 1-xu-xv
      uPv = xu+xv

      if (xu < 0._dp .or. xu > 1-xv .or. xv < 0._dp .or. xv > 1) then
        write (*,*) "broken point! xu, xv", xu, xv
        if ( abs(xv) < 100*tiny .or. abs(xu) < 100*tiny ) then
          write (*,*) "trying to fix"
          xv = abs(xv)
          xu = abs(xu)
        else
          ! for now no further fixes, let's hope tdhpl is fine with it
        endif
      endif

      if (maxLoops == 2) then
        call tdhpl(xu,xv,nw,GYZ1,GYZ2,GYZ3,GYZ4,HZ1,HZ2,HZ3,HZ4)
        include 'omega_w_2loop.f'
      elseif (maxLoops == 1) then
        ! this does not contain the two-loop coefficients
        include 'omega_w_1loop.f'
      elseif (maxLoops == 0) then
        include 'omega_w_0loop.f'
      else
        call abort
      endif

      c_alphai(0,1) = alphaG1WL0
      c_alphai(0,2) = alphaG2WL0
      c_alphai(0,3) = alphaG3WL0
      c_alphai(0,4) = 0

      c_alphai(1,1) = alphaG1WL1
      c_alphai(1,2) = alphaG2WL1
      c_alphai(1,3) = alphaG3WL1
      c_alphai(1,4) = 0

      c_alphai(2,1) = alphaG1WL2
      c_alphai(2,2) = alphaG2WL2
      c_alphai(2,3) = alphaG3WL2
      c_alphai(2,4) = alphaG4WL2

      c_betai(0,1) = betaG1WL0
      c_betai(0,2) = betaG2WL0
      c_betai(0,3) = betaG3WL0
      c_betai(0,4) = 0

      c_betai(1,1) = betaG1WL1
      c_betai(1,2) = betaG2WL1
      c_betai(1,3) = betaG3WL1
      c_betai(1,4) = 0

      c_betai(2,1) = betaG1WL2
      c_betai(2,2) = betaG2WL2
      c_betai(2,3) = betaG3WL2
      c_betai(2,4) = betaG4WL2

      c_gammai(0,1) = gammaG1WL0
      c_gammai(0,2) = gammaG2WL0
      c_gammai(0,3) = gammaG3WL0
      c_gammai(0,4) = 0

      c_gammai(1,1) = gammaG1WL1
      c_gammai(1,2) = gammaG2WL1
      c_gammai(1,3) = gammaG3WL1
      c_gammai(1,4) = 0

      c_gammai(2,1) = gammaG1WL2
      c_gammai(2,2) = gammaG2WL2
      c_gammai(2,3) = gammaG3WL2
      c_gammai(2,4) = gammaG4WL2

c     ! various checks on HPL's and coefficients

c     write (*,*) "c_gammai(1,1)", c_gammai(1,1),
c    & xv*(1-xu-xv)/2d0/(1-xv)**2 *  HZ1(0) + (1-xu-xv)/2d0/(1-xv)
c     
c     write (*,*) "c_alphai(1,1)" , c_alphai(1,1),
c    & zeta2 + GYZ1(3)*im*pi - GYZ1(3)*HZ1(0) + GYZ2(3,0) + GYZ1(0)*
c    & HZ1(0) + HZ1(0)*im*pi - HZ2(0,0) + HZ2(1,0)

c     write (*,*) "c_alphai(1,1) dilog", c_alphai(1,1),
c    &  + zeta2 + Log(1 + xu*xv**(-1))*Pi*im + Log(xu)*Log(1 + xu*
c    &    xv**(-1)) + Log(xu)*Log(xv) + Log(xv)*Pi*im - Log(xv)*Log(1
c    &     - xv) - Log(xv)*Log(1 + xu*xv**(-1)) - 1.D0/2.D0*Log(xv)**2
c    &     + ddilog( - xu*xv**(-1)) - ddilog(xv)


c     write (*,*) "c_alphai(1,2)", c_alphai(1,2),
c    &  - 7.D0/2.D0 - 1.D0/2.D0*xv*OneMu**(-1) + Pi*im*xu**(-1) - 2*Pi*
c    &    im*xu**(-1)*xv + Pi*im*xu**(-1)*xv**2 + 1.D0/2.D0*Pi*im - Pi*
c    &    im*xv*OneMu**(-1) - 1.D0/2.D0*Pi*im*xv**2*OneMu**(-2) + Pi*im
c    &    *xv**2*OneMu**(-1) + HZ1(0)*xu**(-1)*xv - HZ1(0)*xu**(-1)*
c    &    xv**2 + HZ1(0)*xv*OneMu**(-1) + 1.D0/2.D0*HZ1(0)*xv**2*
c    &    OneMu**(-2) - HZ1(0)*xv**2*OneMu**(-1) - HZ1(0)*GYZ1(1)*
c    &    xu**(-2) + 2*HZ1(0)*GYZ1(1)*xu**(-2)*xv - HZ1(0)*GYZ1(1)*
c    &    xu**(-2)*xv**2 + HZ1(0)*GYZ1(2)*xu**(-2) - 2*HZ1(0)*GYZ1(2)*
c    &    xu**(-2)*xv + HZ1(0)*GYZ1(2)*xu**(-2)*xv**2 - 1.D0/2.D0*HZ1(1
c    &    ) - HZ1(1)*xu**(-1) + 2*HZ1(1)*xu**(-1)*xv - HZ1(1)*xu**(-1)*
c    &    xv**2 + HZ1(1)*xv*OneMu**(-1) + 1.D0/2.D0*HZ1(1)*xv**2*
c    &    OneMu**(-2) - HZ1(1)*xv**2*OneMu**(-1) - HZ1(1)*GYZ1(1)*
c    &    xu**(-2) + 2*HZ1(1)*GYZ1(1)*xu**(-2)*xv - HZ1(1)*GYZ1(1)*
c    &    xu**(-2)*xv**2 + GYZ1(1)*Pi*im*xu**(-2) - 2*GYZ1(1)*Pi*im*
c    &    xu**(-2)*xv
c    &    + GYZ1(1)*Pi*im*xu**(-2)*xv**2 + 1.D
c    &    0/2.D0*GYZ1(2) + GYZ1(2)*xu**(-1) - 2*GYZ1(2)*xu**(-1)*xv + 
c    &    GYZ1(2)*xu**(-1)*xv**2 - GYZ1(2)*xv*OneMu**(-1) - 1.D0/2.D0*
c    &    GYZ1(2)*xv**2*OneMu**(-2) + GYZ1(2)*xv**2*OneMu**(-1) + GYZ2(
c    &    1,2)*xu**(-2) - 2*GYZ2(1,2)*xu**(-2)*xv + GYZ2(1,2)*xu**(-2)*
c    &    xv**2

c     write (*,*) "c_alphai(1,2) dilog", c_alphai(1,2),
c    &  - 7.D0/2.D0 - 1.D0/2.D0*xv*OneMu**(-1) + Pi*im*xu**(-1) - 2*Pi*
c    &    im*xu**(-1)*xv + Pi*im*xu**(-1)*xv**2 + 1.D0/2.D0*Pi*im - Pi*
c    &    im*xv*OneMu**(-1) - 1.D0/2.D0*Pi*im*xv**2*OneMu**(-2) + Pi*im
c    &    *xv**2*OneMu**(-1) + 1.D0/2.D0*Log(1 - xv) + Log(1 - xv)*
c    &    xu**(-1) - 2*Log(1 - xv)*xu**(-1)*xv + Log(1 - xv)*xu**(-1)*
c    &    xv**2 - Log(1 - xv)*xv*OneMu**(-1) - 1.D0/2.D0*Log(1 - xv)*
c    &    xv**2*OneMu**(-2) + Log(1 - xv)*xv**2*OneMu**(-1) + Log(1 - 
c    &    xu)*Pi*im*xu**(-2) - 2*Log(1 - xu)*Pi*im*xu**(-2)*xv + Log(1
c    &     - xu)*Pi*im*xu**(-2)*xv**2 + 1.D0/2.D0*Log(1 - xu)**2*
c    &    xu**(-2) - Log(1 - xu)**2*xu**(-2)*xv + 1.D0/2.D0*Log(1 - xu)
c    &    **2*xu**(-2)*xv**2 + 1.D0/2.D0*Log(1 - 1/(1 - xv)*xu) + Log(1
c    &     - 1/(1 - xv)*xu)*xu**(-1) - 2*Log(1 - 1/(1 - xv)*xu)*
c    &    xu**(-1)*xv + Log(1 - 1/(1 - xv)*xu)*xu**(-1)*xv**2 - Log(1
c    &     - 1/(1 - xv)*xu)*xv*OneMu**(-1) - 1.D0/2.D0*Log(1 - 1/(1 - 
c    &    xv)*xu)*xv**2*OneMu**(-2)
c    &     + Log(1 - 1/(1 - xv)*xu)*xv**2*OneMu**(-1) + Log(xv)*
c    &    xu**(-1)*xv - Log(xv)*xu**(-1)*xv**2 + Log(xv)*xv*OneMu**(-1)
c    &     + 1.D0/2.D0*Log(xv)*xv**2*OneMu**(-2) - Log(xv)*xv**2*
c    &    OneMu**(-1) - Log(xv)*Log(1 - xu)*xu**(-2) + 2*Log(xv)*Log(1
c    &     - xu)*xu**(-2)*xv - Log(xv)*Log(1 - xu)*xu**(-2)*xv**2 + 
c    &    Log(xv)*Log(1 - 1/(1 - xv)*xu)*xu**(-2) - 2*Log(xv)*Log(1 - 
c    &    1/(1 - xv)*xu)*xu**(-2)*xv + Log(xv)*Log(1 - 1/(1 - xv)*xu)*
c    &    xu**(-2)*xv**2 + ddilog(1/(1 - xu)*xv)*xu**(-2) - 2*ddilog(
c    &    1/(1 - xu)*xv)*xu**(-2)*xv + ddilog(1/(1 - xu)*xv)*xu**(-2)*
c    &    xv**2 - ddilog(xv)*xu**(-2) + 2*ddilog(xv)*xu**(-2)*xv - 
c    &    ddilog(xv)*xu**(-2)*xv**2


c     write (*,*) "c_alphai(1,3)", c_alphai(1,3),
c    &  - 4d0*(1-xu-xv)/(1-xv)


c     write (*,*) "c_betai(1,1)", c_betai(1,1),
c    &  - 9.D0/2.D0 + zeta2 + 1.D0/2.D0*xu*OneMv**(-1) + 3.D0/2.D0*HZ1(
c    & 0) - 3.D0/2.D0*HZ1(0)*OneMv**(-1) + 1.D0/2.D0*HZ1(0)*xu*
c    & OneMv**(-2) - 1.D0/2.D0*HZ1(0)*xu*OneMv**(-1) + HZ1(0)*Pi*im + 
c    & HZ1(0)*GYZ1(0) - HZ1(0)*GYZ1(3) - HZ2(0,0) + HZ2(1,0) + GYZ1(3)*
c    & Pi*im + GYZ2(3,0)

c     write (*,*) "c_betai(1,2)", c_betai(1,2),
c    & 1.D0/2.D0 + 1.D0/2.D0*xu*OneMv**(-1) + Pi*im*xu**(-1) - Pi*im*
c    & xu**(-1)*xv + 1.D0/2.D0*Pi*im + 1.D0/2.D0*Pi*im*xv*OneMu**(-1)
c    &  - 3.D0/2.D0*HZ1(0) + HZ1(0)*xu**(-1)*xv + 3.D0/2.D0*HZ1(0)*
c    & OneMv**(-1) - 1.D0/2.D0*HZ1(0)*xv*OneMu**(-1) + 1.D0/2.D0*HZ1(0)
c    & *xu*OneMv**(-2) - 1.D0/2.D0*HZ1(0)*xu*OneMv**(-1) - HZ1(0)*GYZ1(
c    & 1)*xu**(-2) + HZ1(0)*GYZ1(1)*xu**(-2)*xv - HZ1(0)*GYZ1(1)*
c    & xu**(-1)*xv + HZ1(0)*GYZ1(2)*xu**(-2) - HZ1(0)*GYZ1(2)*xu**(-2)*
c    & xv + HZ1(0)*GYZ1(2)*xu**(-1)*xv - 1.D0/2.D0*HZ1(1) - HZ1(1)*
c    & xu**(-1) + HZ1(1)*xu**(-1)*xv - 1.D0/2.D0*HZ1(1)*xv*OneMu**(-1)
c    &  - HZ1(1)*GYZ1(1)*xu**(-2) + HZ1(1)*GYZ1(1)*xu**(-2)*xv - HZ1(1)
c    & *GYZ1(1)*xu**(-1)*xv + GYZ1(1)*Pi*im*xu**(-2) - GYZ1(1)*Pi*im*
c    & xu**(-2)*xv + GYZ1(1)*Pi*im*xu**(-1)*xv + 1.D0/2.D0*GYZ1(2) + 
c    & GYZ1(2)*xu**(-1) - GYZ1(2)*xu**(-1)*xv + 1.D0/2.D0*GYZ1(2)*xv*
c    & OneMu**(-1) + GYZ2(1,2)*xu**(-2) - GYZ2(1,2)*xu**(-2)*xv + GYZ2(
c    & 1,2)*xu**(-1)*xv

c     write (*,*) "c_betai(1,2) dilog", c_betai(1,2),
c    &  + 1.D0/2.D0 + 1.D0/2.D0*xu*OneMv**(-1) + Pi*im*xu**(-1) - Pi*im
c    &    *xu**(-1)*xv + 1.D0/2.D0*Pi*im + 1.D0/2.D0*Pi*im*xv*
c    &    OneMu**(-1) + 1.D0/2.D0*Log(1 - xv) + Log(1 - xv)*xu**(-1) - 
c    &    Log(1 - xv)*xu**(-1)*xv + 1.D0/2.D0*Log(1 - xv)*xv*
c    &    OneMu**(-1) + Log(1 - xu)*Pi*im*xu**(-2) - Log(1 - xu)*Pi*im*
c    &    xu**(-2)*xv + Log(1 - xu)*Pi*im*xu**(-1)*xv + 1.D0/2.D0*Log(1
c    &     - xu)**2*xu**(-2) - 1.D0/2.D0*Log(1 - xu)**2*xu**(-2)*xv + 1.
c    &    D0/2.D0*Log(1 - xu)**2*xu**(-1)*xv + 1.D0/2.D0*Log(1 - 1/(1
c    &     - xv)*xu) + Log(1 - 1/(1 - xv)*xu)*xu**(-1) - Log(1 - 1/(1
c    &     - xv)*xu)*xu**(-1)*xv + 1.D0/2.D0*Log(1 - 1/(1 - xv)*xu)*xv*
c    &    OneMu**(-1) - 3.D0/2.D0*Log(xv) + Log(xv)*xu**(-1)*xv + 3.D0/
c    &    2.D0*Log(xv)*OneMv**(-1) - 1.D0/2.D0*Log(xv)*xv*OneMu**(-1)
c    &     + 1.D0/2.D0*Log(xv)*xu*OneMv**(-2) - 1.D0/2.D0*Log(xv)*xu*
c    &    OneMv**(-1) - Log(xv)*Log(1 - xu)*xu**(-2) + Log(xv)*Log(1 - 
c    &    xu)*xu**(-2)*xv
c    &     - Log(xv)*Log(1 - xu)*xu**(-1)*xv + Log(xv)*Log(1 - 
c    &    1/(1 - xv)*xu)*xu**(-2) - Log(xv)*Log(1 - 1/(1 - xv)*xu)*
c    &    xu**(-2)*xv + Log(xv)*Log(1 - 1/(1 - xv)*xu)*xu**(-1)*xv + 
c    &    ddilog(1/(1 - xu)*xv)*xu**(-2) - ddilog(1/(1 - xu)*xv)*
c    &    xu**(-2)*xv + ddilog(1/(1 - xu)*xv)*xu**(-1)*xv - ddilog(xv)*
c    &    xu**(-2) + ddilog(xv)*xu**(-2)*xv - ddilog(xv)*xu**(-1)*xv




c     write (*,*) "c_betai(1,3)", c_betai(1,3),
c    &   4d0*xu/(1-xv)


c     ! checks on HPL's and 2dHPL's
c     write (*,*) "v = ", xv
c     write (*,*) "H(0,v)  = ", HZ1(0), log(xv)
c     write (*,*) "H(1,v)  = ", HZ1(1), -log(1-xv)

c     write (*,*) "H(0,0,v) = ", HZ2(0,0),
c    & 0.5d0 * log(xv)**2
c     write (*,*) "H(1,0,v) = ", HZ2(1,0),
c    &    -log(xv)*log(1-xv)-ddilog(xv)

c     write (*,*) "G(0,u) = ", GYZ1(0), log(xu)
c     write (*,*) "G(1,xu) = ", GYZ1(1), log(1-xu)
c     write (*,*) "G(1-v,u) = ", GYZ1(2), Log(1-xu/(1-xv))
c     write (*,*) "G(-v,u) = ", GYZ1(3), Log((xu+xv)/xv)

c     write (*,*) "G(1,1-v,u) = ", GYZ2(1,2),
c    & 1/2d0*Log(1-xu)**2 -Log(1-xu)*Log(1-xv)
c    & + ddilog(xv/(1-xu)) - ddilog(xv);


c     write (*,*) "G(-v,0,u) = ", GYZ2(3,0),
c    &  Log(xu)*Log((xu+xv)/xv) + ddilog(-xu/xv);

      end subroutine
