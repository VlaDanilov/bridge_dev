      function omega_I1(sij)
      implicit none
      include 'types.f'
      complex(dp):: omega_I1
      
      include 'epinv.f'
      include 'epinv2.f'
      include 'scale.f'
      include 'constants.f' ! cf
      real(dp):: sij
      complex(dp):: Lnrat,xl12, cln
        
      xl12=cln(-sij/musq,-1._dp)

      omega_I1=-epinv*epinv2+epinv*(-1.5_dp+xl12)
     &   -0.5_dp*xl12**2+1.5_dp*xl12 + 0.5_dp
      ! +  0.5_dp DRED conversion factor
      omega_I1 = omega_I1 * cf

      return
      end

      ! MSBAR IR+UV evolution for two initial state quarks i1,i2
      ! assuming mu0sq as reference scale squared
      subroutine evolve_msbar_ir_qq(i5,i6,i1,i3,i2,za,zb,amps,mu0sq)
          implicit none
          include 'types.f'
          include 'constants.f'
          include 'zeta.f'
          include 'mxpart.f'
          include 'nf.f'
          include 'scale.f'

          integer, intent(in) :: i5,i6,i1,i3,i2
          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          complex(dp), intent(inout) :: amps(0:2)
          real(dp), intent(in) :: mu0sq

          complex(dp) :: ampsmu0(0:2)
          real(dp) :: s12, xi, xi0

          s12 = real(za(i1,i2)*zb(i2,i1))
          xi = log(musq/s12)
          xi0 = log(mu0sq/s12)

          ampsmu0(:) = amps(:)

          amps(0) = ampsmu0(0)
          amps(1) = ampsmu0(1) + CF*ampsmu0(0)*(
     &          (-3._dp/2 - im*pi)*xi - xi**2/2 +
     &          (+3._dp/2 + im*pi)*xi0 + xi0**2/2 )
          ! numerically tested against mathematica notebook
          amps(2) = ampsmu0(2) + 
     &      ampsmu0(1) * (
     &        CA*(11*xi/6 - 11*xi0/6) + NF*(-xi/3 + xi0/3)
     &       + CF*( (-3._dp/2-im*pi)*xi - xi**2/2
     &             +(+3._dp/2+im*pi)*xi0 +  xi0**2/2)) 
     &     + ampsmu0(0) * (
     &     (CF*NF*((19._dp/36 + (Im/6)*Pi)*xi**2 + xi**3/18 + 
     &              (-65._dp/108 - ((5*Im)/9)*Pi - Pi**2/12)*xi0 + 
     &              (-1._dp/36 + (Im/6)*Pi)*xi0**2 + xi0**3/9
     &               + xi*(65._dp/108 + ((5*Im)/9)*Pi + Pi**2/12 + 
     &               (-1._dp/2 - (Im/3)*Pi)*xi0 - xi0**2/6))
     &       + CF**2*((3._dp/4 + (Im/2)*Pi)*xi**3 + xi**4/8 + 
     &               (9._dp/8 + ((3*Im)/2)*Pi - Pi**2/2)*xi0**2 +
     &                 (3._dp/4 + (Im/2)*Pi)*xi0**3 + xi0**4/8 + 
     &                  xi**2*(9._dp/8 + ((3*Im)/2)*Pi - Pi**2/2 +
     &                 (-3._dp/4 - (Im/2)*Pi)*xi0 - xi0**2/4) + 
     &                 xi*(-3._dp/8 + Pi**2/2 + (-9._dp/4 - (3*Im)*Pi
     &                  + Pi**2)*xi0 + (-3._dp/4 - (Im/2)*Pi)*xi0**2
     &                  - 6*zeta3) + 
     &                 xi0*(3._dp/8 - Pi**2/2 + 6*zeta3))
     &        + CA*CF*((-233._dp/72 - ((11*Im)/12)*Pi + Pi**2/12)*xi**2
     &                  - (11*xi**3)/36 + 
     &                  (35._dp/72 - ((11*Im)/12)*Pi - Pi**2/12)*xi0**2
     &                    - (11*xi0**3)/18 + 
     &                 xi0*(961._dp/216 + ((67*Im)/18)*Pi + (11*Pi**2)/24
     &                    - (Im/6)*Pi**3 - (13*zeta3)/2) + 
     &                 xi*(-961._dp/216 - ((67*Im)/18)*Pi - (11*Pi**2)/24
     &                 + (Im/6)*Pi**3 + (11._dp/4 + ((11*Im)/6)*Pi)*xi0 + 
     &                   (11*xi0**2)/12 + (13*zeta3)/2)))
     &      )

      end subroutine

      ! see 1309.3245
      subroutine zgamma_catani_to_msbar(i5,i6,i1,i3,i2,za,zb,amps,musq)
          implicit none
          include 'types.f'
          include 'constants.f'
          include 'nf.f'
          include 'mxpart.f'
          include 'zeta.f'

          integer, intent(in) :: i5,i6,i1,i3,i2
          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          complex(dp), intent(inout) :: amps(0:2)
          real(dp), intent(in) :: musq

          complex(dp) :: ampsfin(0:2)
          real(dp) :: s12, s123
          complex(dp) :: Lnrat,xl12

          complex(dp) :: C0,C1
          real(dp) :: gamma0cusp, gamma1cusp
          real(dp) :: gamma0q
          real(dp) :: beta0
          complex(dp) :: GammaZero, GammaPrimeZero

          s12 = real(za(i1,i2)*zb(i2,i1))
          s123 = real(za(i1,i2)*zb(i2,i1) + za(i1,i3)*zb(i3,i1) +
     &                             za(i2,i3)*zb(i3,i2))
          xl12=Lnrat(musq,-s12)
          ampsfin(:) = amps(:)

          ! all specific to two partons, being quarks

          gamma0cusp = 4
          gamma1cusp = (268._dp/9 - 4*pisq/3)*cA - 80*tr*nf/9
          gamma0q = -3*cf
          GammaZero = 2*(-cf/2*gamma0cusp*xl12) + 2*gamma0q
          GammaPrimeZero = -gamma0cusp*(2*cf)

          beta0 = 11*cA/3 - 4*tr*nf/3

          C0 = 2*(-cf/16*(gamma0cusp*xl12**2 - 4*gamma0q/cf*xl12))
     &          - pisq/96*GammaPrimeZero
          C1 = 2*(-cf/48*(gamma0cusp*xl12**3 - 6*gamma0q/cf*xl12**2) )
     &          - pisq/48*GammaZero - zeta3/24*GammaPrimeZero

          amps(1) = ampsfin(1) + C0*ampsfin(0)
          amps(2) = ampsfin(2) + C0*ampsfin(1) + ampsfin(0)*(
     &      C0**2/2 + gamma1cusp/8*(C0 + pisq/128*GammaPrimeZero)
     &      + beta0/2*(C1 + pisq/32*GammaZero
     &                    + 7*zeta3/96*GammaPrimeZero)
     &     )

      end subroutine
