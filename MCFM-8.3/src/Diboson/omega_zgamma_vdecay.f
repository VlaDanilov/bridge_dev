      !!
      !! only for Z decay, omega_3 coefficients for W are not included!
      !!
      subroutine omega_zgamma_vdecay(mp, i5,i6,i1,i3,i2, maxLoops,
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

      real(dp) :: y,z,OneMy,OneMz,yPz,OneMyMz;
      real(dp) :: s12,s13,s23,s123;

      integer :: nw ! maximum weight for 2D HPL's

      ! for two-dimensional harmonic polylogarithms
      real(dp) :: HZ1,HZ2,HZ3,HZ4,GYZ1,GYZ2,GYZ3,GYZ4
      dimension HZ1(0:1),HZ2(0:1,0:1),HZ3(0:1,0:1,0:1),
     &          HZ4(0:1,0:1,0:1,0:1)
      dimension GYZ1(0:3),GYZ2(0:3,0:3),GYZ3(0:3,0:3,0:3),
     &          GYZ4(0:3,0:3,0:3,0:3)

      real(dp) :: dotvec

      complex(dp) :: w(772)
      complex(dp) :: alphaG1WL0, alphaG2WL0
      complex(dp) :: alphaG1WL1, alphaG2WL1 
      complex(dp) :: alphaG1WL2, alphaG2WL2, alphaG4WL2
      complex(dp) :: betaG1WL0, betaG2WL0
      complex(dp) :: betaG1WL1, betaG2WL1
      complex(dp) :: betaG1WL2, betaG2WL2, betaG4WL2
      complex(dp) :: gammaG1WL0, gammaG2WL0
      complex(dp) :: gammaG1WL1, gammaG2WL1
      complex(dp) :: gammaG1WL2, gammaG2WL2, gammaG4WL2

      s12 = dotvec(mp(i1,:)+mp(i2,:),mp(i1,:)+mp(i2,:))
      s13 = dotvec(mp(i1,:)+mp(i3,:),mp(i1,:)+mp(i3,:))
      s23 = dotvec(mp(i2,:)+mp(i3,:),mp(i2,:)+mp(i3,:))
      s123 = dotvec(mp(i1,:)+mp(i2,:)+mp(i3,:),mp(i1,:)+mp(i2,:)+mp(i3,:))

      y = s13/s123
      z = s23/s123
      OneMy = 1-y
      OneMz = 1-z
      yPz = y+z
      OneMyMz = 1-y-z

      if (y < 0._dp .or. y > 1-z .or. z < 0._dp .or. z > 1) then
        write (*,*) "broken point! y,z", y,z
        if ( abs(y) < 100*tiny .or. abs(z) < 100*tiny ) then
          write (*,*) "trying to fix"
          y = abs(y)
          z = abs(z)
        else
          ! for now no further fixes, let's hope tdhpl is fine with it
        endif
      endif

      w = 0._dp

      if (maxLoops == 2) then
        nw = 4
        call tdhpl(y,z,nw,GYZ1,GYZ2,GYZ3,GYZ4,HZ1,HZ2,HZ3,HZ4)
        include 'omega_w_2loop_vdecay.f'
      elseif (maxLoops == 1) then
        nw = 2
        call tdhpl(y,z,nw,GYZ1,GYZ2,GYZ3,GYZ4,HZ1,HZ2,HZ3,HZ4)
        include 'omega_w_1loop_vdecay.f'
      elseif (maxLoops == 0) then
        include 'omega_w_0loop_vdecay.f'
      else
        call abort
      endif

      c_alphai(0,1) = alphaG1WL0
      c_alphai(0,2) = alphaG2WL0
      c_alphai(0,3) = 0 ! not included
      c_alphai(0,4) = 0

      c_alphai(1,1) = alphaG1WL1
      c_alphai(1,2) = alphaG2WL1
      c_alphai(1,3) = 0 !not included
      c_alphai(1,4) = 0

      c_alphai(2,1) = alphaG1WL2
      c_alphai(2,2) = alphaG2WL2
      c_alphai(2,3) = 0 ! not included
      c_alphai(2,4) = alphaG4WL2

      c_betai(0,1) = betaG1WL0
      c_betai(0,2) = betaG2WL0
      c_betai(0,3) = 0 ! not included
      c_betai(0,4) = 0

      c_betai(1,1) = betaG1WL1
      c_betai(1,2) = betaG2WL1
      c_betai(1,3) = 0 ! not included
      c_betai(1,4) = 0

      c_betai(2,1) = betaG1WL2
      c_betai(2,2) = betaG2WL2
      c_betai(2,3) = 0 ! not included
      c_betai(2,4) = betaG4WL2

      c_gammai(0,1) = gammaG1WL0
      c_gammai(0,2) = gammaG2WL0
      c_gammai(0,3) = 0 ! not included
      c_gammai(0,4) = 0

      c_gammai(1,1) = gammaG1WL1
      c_gammai(1,2) = gammaG2WL1
      c_gammai(1,3) = 0 ! not included
      c_gammai(1,4) = 0

      c_gammai(2,1) = gammaG1WL2
      c_gammai(2,2) = gammaG2WL2
      c_gammai(2,3) = 0 ! not included
      c_gammai(2,4) = gammaG4WL2

      end subroutine
