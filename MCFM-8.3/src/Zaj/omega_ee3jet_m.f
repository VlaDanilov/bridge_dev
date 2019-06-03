      ! gamma*(p4) - >  q(p1) +  qbar(p2) + g(p3)
      ! based on hep-ph/0206067 and analytic continuations hep-ph/0207020
      ! analytic continuations to regions (2a)+, (3a)+, (4a)+
      module omega_ee3jet_m
        implicit none
        include 'types.f'
        include 'nf.f'
        include 'constants.f'

        real(dp), parameter, private :: beta0 = (11*CA-4*TR*Nf)/6

      abstract interface 
        function anomAmp(i1,i2,i3,i4,i5,i6,za,zb,flip,c_alpha,c_beta,c_gamma)
          implicit none
          include 'types.f'
          include 'mxpart.f'

          complex(dp) :: anomAmp(0:1)
          integer, intent(in) :: i1,i2,i3,i4,i5,i6
          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          logical, intent(in) :: flip
          complex(dp), intent(in) :: c_alpha(0:1), c_beta(0:1), c_gamma(0:1)
        end function
      end interface

      contains

        ! for normal amplitudes the delta contribution is always zero,
        ! I include it nevertheless as a small check
        function amp_zaj_fsr_pp(i1,i2,i3,i4,i5,i6, za,zb,
     &                c_alphai, c_betai, c_gammai)
          include 'types.f'
          include 'mxpart.f'
          include 'scale.f'

          complex(dp) :: amp_zaj_fsr_pp(0:1)
          integer, intent(in) :: i1,i2,i3,i4,i5,i6
          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          complex(dp), intent(in) :: c_alphai(0:1), c_betai(0:1), c_gammai(0:1)

          complex(dp) :: c_deltai(0:1)

          real(dp) :: s,t
          s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
          t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

          c_deltai(:) = (c_alphai(:) - c_betai(:) - c_gammai(:))
     &                      * s(i1,i2)/t(i1,i2,i3)/2._dp

          amp_zaj_fsr_pp(:) = 
     -  (c_alphai(:)*(-(s(i1,i2)*s(i4,i5)*s(i4,i6)*za(i2,i5)*zb(i3,i2)*
     -          zb(i4,i1)) + 
     -       s(i1,i2)*s(i4,i6)*za(i2,i5)*za(i5,i6)*zb(i3,i2)*zb(i5,i4)*
     -        zb(i6,i1)))/
     -   (s(i2,i3)*s(i4,i5)*s(i4,i6)*za(i1,i3)*za(i4,i6)*zb(i2,i1)) + 
     -  (c_betai(:)*(s(i2,i3)*s(i4,i5)*s(i4,i6)*za(i2,i5)*zb(i2,i1)*
     -        zb(i4,i3) - s(i2,i3)*s(i4,i6)*za(i2,i5)*za(i5,i6)*
     -        zb(i2,i1)*zb(i5,i4)*zb(i6,i3)))/
     -   (s(i2,i3)*s(i4,i5)*s(i4,i6)*za(i1,i3)*za(i4,i6)*zb(i2,i1)) + 
     -  (c_deltai(:)*(s(i2,i3)*s(i4,i5)*s(i4,i6)*za(i1,i5)*zb(i3,i1)*
     -        zb(i4,i1) + s(i2,i3)*s(i4,i5)*s(i4,i6)*za(i2,i5)*
     -        zb(i3,i1)*zb(i4,i2) + 
     -       s(i2,i3)*s(i4,i5)*s(i4,i6)*za(i3,i5)*zb(i3,i1)*zb(i4,i3) - 
     -       s(i2,i3)*s(i4,i6)*za(i1,i5)*za(i5,i6)*zb(i3,i1)*zb(i5,i4)*
     -        zb(i6,i1) - s(i2,i3)*s(i4,i6)*za(i2,i5)*za(i5,i6)*
     -        zb(i3,i1)*zb(i5,i4)*zb(i6,i2) - 
     -       s(i2,i3)*s(i4,i6)*za(i3,i5)*za(i5,i6)*zb(i3,i1)*zb(i5,i4)*
     -        zb(i6,i3)))/
     -   (s(i2,i3)*s(i4,i5)*s(i4,i6)*za(i1,i3)*za(i4,i6)*zb(i2,i1)) + 
     -  (c_gammai(:)*(-(s(i4,i5)*s(i4,i6)*za(i1,i5)*za(i2,i3)*zb(i2,i1)*
     -          zb(i3,i1)*zb(i4,i3)) + 
     -       s(i4,i6)*za(i1,i5)*za(i2,i3)*za(i5,i6)*zb(i2,i1)*zb(i3,i1)*
     -        zb(i5,i4)*zb(i6,i3)))/
     -   (s(i2,i3)*s(i4,i5)*s(i4,i6)*za(i1,i3)*za(i4,i6)*zb(i2,i1))

          ! scale evolution
          amp_zaj_fsr_pp(1) = amp_zaj_fsr_pp(1) + 
     &        beta0/2*log(musq/t(i1,i2,i3)) * amp_zaj_fsr_pp(0)

        end function

        function amp_zaj_fsr_pm(i1,i2,i3,i4,i5,i6, za,zb,
     &                c_alphai, c_betai, c_gammai)
          include 'types.f'
          include 'mxpart.f'
          include 'scale.f'

          complex(dp) :: amp_zaj_fsr_pm(0:1)
          integer, intent(in) :: i1,i2,i3,i4,i5,i6
          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          complex(dp), intent(in) :: c_alphai(0:1), c_betai(0:1), c_gammai(0:1)

          complex(dp) :: c_deltai(0:1)

          real(dp) :: s,t
          s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
          t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

          c_deltai(:) = (c_alphai(:) - c_betai(:) - c_gammai(:))
     &                      * s(i1,i2)/t(i1,i2,i3)/2._dp

          amp_zaj_fsr_pm(:) = 
     -  -((c_alphai(:)*za(i1,i2)*zb(i6,i1)*
     -       (-(za(i2,i4)*zb(i6,i4)) - za(i2,i5)*zb(i6,i5)))/
     -     (za(i1,i3)*za(i2,i3)*zb(i5,i4)*zb(i6,i4))) - 
     -  (c_betai(:)*(s(i2,i3)*za(i2,i4)*zb(i2,i1)*zb(i6,i3)*zb(i6,i4) + 
     -       s(i2,i3)*za(i2,i5)*zb(i2,i1)*zb(i6,i3)*zb(i6,i5)))/
     -   (s(i2,i3)*za(i1,i3)*zb(i2,i1)*zb(i5,i4)*zb(i6,i4)) - 
     -  (c_deltai(:)*(s(i2,i3)*za(i1,i4)*zb(i3,i1)*zb(i6,i1)*zb(i6,i4) + 
     -       s(i2,i3)*za(i2,i4)*zb(i3,i1)*zb(i6,i2)*zb(i6,i4) + 
     -       s(i2,i3)*za(i3,i4)*zb(i3,i1)*zb(i6,i3)*zb(i6,i4) + 
     -       s(i2,i3)*za(i1,i5)*zb(i3,i1)*zb(i6,i1)*zb(i6,i5) + 
     -       s(i2,i3)*za(i2,i5)*zb(i3,i1)*zb(i6,i2)*zb(i6,i5) + 
     -       s(i2,i3)*za(i3,i5)*zb(i3,i1)*zb(i6,i3)*zb(i6,i5)))/
     -   (s(i2,i3)*za(i1,i3)*zb(i2,i1)*zb(i5,i4)*zb(i6,i4)) - 
     -  (c_gammai(:)*(-(za(i1,i4)*za(i2,i3)*zb(i2,i1)*zb(i3,i1)*zb(i6,i3)*
     -          zb(i6,i4)) - 
     -       za(i1,i5)*za(i2,i3)*zb(i2,i1)*zb(i3,i1)*zb(i6,i3)*zb(i6,i5)
     -       ))/(s(i2,i3)*za(i1,i3)*zb(i2,i1)*zb(i5,i4)*zb(i6,i4))
          
          ! scale evolution
          amp_zaj_fsr_pm(1) = amp_zaj_fsr_pm(1) + 
     &        beta0/2*log(musq/t(i1,i2,i3)) * amp_zaj_fsr_pm(0)

        end function

        function amp_zaj_anomZZ_pp(i1,i2,i3,i4,i5,i6, za,zb, flip,
     &                c_alphai, c_betai, c_gammai)
          include 'types.f'
          include 'mxpart.f'
          include 'scale.f'
          include 'anomcoup.f'

          complex(dp) :: amp_zaj_anomZZ_pp(0:1)
          integer, intent(in) :: i1,i2,i3,i4,i5,i6
          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          logical, intent(in) :: flip
          complex(dp), intent(in) :: c_alphai(0:1), c_betai(0:1), c_gammai(0:1)

          complex(dp) :: c_deltai(0:1)
          complex(dp) :: ha, hb
          real(dp) :: hSign

          real(dp) :: s,t
          s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
          t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

          c_deltai(:) = (c_alphai(:) - c_betai(:) - c_gammai(:))
     &                      * s(i1,i2)/t(i1,i2,i3)/2._dp

          hsign = 1._dp
          if (flip) then
            hSign = -1._dp
          endif

          ha = hSign*im*hitZ(1) + hitZ(3)
          hb = hSign*im*hitZ(2) + hitZ(4)

          amp_zaj_anomZZ_pp(:) = 
     -  c_alphai(:)*((ha*(s(i4,i5) + s(i4,i6))*
     -        (2*s(i1,i2)*za(i2,i6)*za(i4,i5)*zb(i3,i2)*zb(i4,i1) + 
     -          2*s(i1,i2)*za(i2,i4)*za(i5,i6)*zb(i3,i2)*zb(i4,i1))*
     -        zb(i6,i4))/(4._dp*s(i2,i3)*za(i1,i3)*za(i4,i6)*zb(i2,i1)) + 
     -     (hb*(s(i4,i5) + s(i4,i6))*
     -        (-(s(i1,i2)*s(i4,i5)*za(i2,i6)*za(i4,i5)*zb(i3,i2)*
     -             zb(i4,i1)) - 
     -          s(i1,i2)*s(i4,i6)*za(i2,i6)*za(i4,i5)*zb(i3,i2)*
     -           zb(i4,i1) - 
     -          s(i1,i2)*s(i4,i5)*za(i2,i4)*za(i5,i6)*zb(i3,i2)*
     -           zb(i4,i1))*zb(i6,i4))/
     -      (4._dp*s(i2,i3)*za(i1,i3)*za(i4,i6)*zb(i2,i1))) + 
     -  c_betai(:)*((ha*(s(i4,i5) + s(i4,i6))*
     -        (-2*s(i2,i3)*za(i2,i6)*za(i4,i5)*zb(i2,i1)*zb(i4,i3) - 
     -          2*s(i2,i3)*za(i2,i4)*za(i5,i6)*zb(i2,i1)*zb(i4,i3))*
     -        zb(i6,i4))/(4._dp*s(i2,i3)*za(i1,i3)*za(i4,i6)*zb(i2,i1)) + 
     -     (hb*(s(i4,i5) + s(i4,i6))*
     -        (s(i2,i3)*s(i4,i5)*za(i2,i6)*za(i4,i5)*zb(i2,i1)*
     -           zb(i4,i3) + 
     -          s(i2,i3)*s(i4,i6)*za(i2,i6)*za(i4,i5)*zb(i2,i1)*
     -           zb(i4,i3) + 
     -          s(i2,i3)*s(i4,i5)*za(i2,i4)*za(i5,i6)*zb(i2,i1)*
     -           zb(i4,i3))*zb(i6,i4))/
     -      (4._dp*s(i2,i3)*za(i1,i3)*za(i4,i6)*zb(i2,i1))) + 
     -  c_deltai(:)*((ha*(s(i4,i5) + s(i4,i6))*
     -        (-2*s(i1,i4)*s(i2,i3)*za(i5,i6)*zb(i3,i1) - 
     -          2*s(i2,i3)*s(i2,i4)*za(i5,i6)*zb(i3,i1) - 
     -          2*s(i2,i3)*s(i3,i4)*za(i5,i6)*zb(i3,i1) - 
     -          2*s(i2,i3)*za(i1,i6)*za(i4,i5)*zb(i3,i1)*zb(i4,i1) - 
     -          2*s(i2,i3)*za(i2,i6)*za(i4,i5)*zb(i3,i1)*zb(i4,i2) - 
     -          2*s(i2,i3)*za(i3,i6)*za(i4,i5)*zb(i3,i1)*zb(i4,i3))*
     -        zb(i6,i4))/(4._dp*s(i2,i3)*za(i1,i3)*za(i4,i6)*zb(i2,i1)) + 
     -     (hb*(s(i4,i5) + s(i4,i6))*
     -        (s(i1,i4)*s(i2,i3)*s(i4,i5)*za(i5,i6)*zb(i3,i1) + 
     -          s(i2,i3)*s(i2,i4)*s(i4,i5)*za(i5,i6)*zb(i3,i1) + 
     -          s(i2,i3)*s(i3,i4)*s(i4,i5)*za(i5,i6)*zb(i3,i1) + 
     -          s(i2,i3)*s(i4,i5)*za(i1,i6)*za(i4,i5)*zb(i3,i1)*
     -           zb(i4,i1) + 
     -          s(i2,i3)*s(i4,i6)*za(i1,i6)*za(i4,i5)*zb(i3,i1)*
     -           zb(i4,i1) + 
     -          s(i2,i3)*s(i4,i5)*za(i2,i6)*za(i4,i5)*zb(i3,i1)*
     -           zb(i4,i2) + 
     -          s(i2,i3)*s(i4,i6)*za(i2,i6)*za(i4,i5)*zb(i3,i1)*
     -           zb(i4,i2) + 
     -          s(i2,i3)*s(i4,i5)*za(i3,i6)*za(i4,i5)*zb(i3,i1)*
     -           zb(i4,i3) + 
     -          s(i2,i3)*s(i4,i6)*za(i3,i6)*za(i4,i5)*zb(i3,i1)*
     -           zb(i4,i3))*zb(i6,i4))/
     -      (4._dp*s(i2,i3)*za(i1,i3)*za(i4,i6)*zb(i2,i1))) + 
     -  c_gammai(:)*((ha*(s(i4,i5) + s(i4,i6))*
     -        (2*za(i1,i6)*za(i2,i3)*za(i4,i5)*zb(i2,i1)*zb(i3,i1)*
     -           zb(i4,i3) + 
     -          2*za(i1,i4)*za(i2,i3)*za(i5,i6)*zb(i2,i1)*zb(i3,i1)*
     -           zb(i4,i3))*zb(i6,i4))/
     -      (4._dp*s(i2,i3)*za(i1,i3)*za(i4,i6)*zb(i2,i1)) + 
     -     (hb*(s(i4,i5) + s(i4,i6))*
     -        (-(s(i4,i5)*za(i1,i6)*za(i2,i3)*za(i4,i5)*zb(i2,i1)*
     -             zb(i3,i1)*zb(i4,i3)) - 
     -          s(i4,i6)*za(i1,i6)*za(i2,i3)*za(i4,i5)*zb(i2,i1)*
     -           zb(i3,i1)*zb(i4,i3) - 
     -          s(i4,i5)*za(i1,i4)*za(i2,i3)*za(i5,i6)*zb(i2,i1)*
     -           zb(i3,i1)*zb(i4,i3))*zb(i6,i4))/
     -      (4._dp*s(i2,i3)*za(i1,i3)*za(i4,i6)*zb(i2,i1)))

          amp_zaj_anomZZ_pp(:) = -amp_zaj_anomZZ_pp(:)

          ! scale evolution
          amp_zaj_anomZZ_pp(1) = amp_zaj_anomZZ_pp(1) + 
     &        beta0/2*log(musq/t(i1,i2,i3)) * amp_zaj_anomZZ_pp(0)
        end function

        function amp_zaj_anomZZ_pm(i1,i2,i3,i4,i5,i6, za,zb, flip,
     &                c_alphai, c_betai, c_gammai)
          include 'types.f'
          include 'mxpart.f'
          include 'scale.f'
          include 'anomcoup.f'

          complex(dp) :: amp_zaj_anomZZ_pm(0:1)
          integer, intent(in) :: i1,i2,i3,i4,i5,i6
          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          logical, intent(in) :: flip
          complex(dp), intent(in) :: c_alphai(0:1), c_betai(0:1), c_gammai(0:1)

          complex(dp) :: c_deltai(0:1)
          complex(dp) :: ha, hb, hSign

          real(dp) :: s,t
          s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
          t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

          c_deltai(:) = (c_alphai(:) - c_betai(:) - c_gammai(:))
     &                      * s(i1,i2)/t(i1,i2,i3)/2._dp

          hsign = 1._dp
          if (flip) then
            hSign = -1._dp
          endif

          ha = -hSign*im*hitZ(1) + hitZ(3)
          hb = -hSign*im*hitZ(2) + hitZ(4)


          amp_zaj_anomZZ_pm(:) = 
     -  c_deltai(:)*((ha*(s(i4,i5) + s(i4,i6))*za(i4,i5)*
     -        (-2*s(i2,i3)*za(i1,i4)*zb(i3,i1)*zb(i6,i1) - 
     -          2*s(i2,i3)*za(i2,i4)*zb(i3,i1)*zb(i6,i2) - 
     -          2*s(i2,i3)*za(i3,i4)*zb(i3,i1)*zb(i6,i3)))/
     -      (4._dp*s(i2,i3)*za(i1,i3)*zb(i2,i1)) + 
     -     (hb*(s(i4,i5) + s(i4,i6))*za(i4,i5)*
     -        (s(i2,i3)*s(i4,i5)*za(i1,i4)*zb(i3,i1)*zb(i6,i1) + 
     -          s(i2,i3)*s(i4,i6)*za(i1,i4)*zb(i3,i1)*zb(i6,i1) + 
     -          s(i2,i3)*s(i4,i5)*za(i2,i4)*zb(i3,i1)*zb(i6,i2) + 
     -          s(i2,i3)*s(i4,i6)*za(i2,i4)*zb(i3,i1)*zb(i6,i2) + 
     -          s(i2,i3)*s(i4,i5)*za(i3,i4)*zb(i3,i1)*zb(i6,i3) + 
     -          s(i2,i3)*s(i4,i6)*za(i3,i4)*zb(i3,i1)*zb(i6,i3) + 
     -          s(i1,i4)*s(i2,i3)*za(i4,i5)*zb(i3,i1)*zb(i6,i5) + 
     -          s(i2,i3)*s(i2,i4)*za(i4,i5)*zb(i3,i1)*zb(i6,i5) + 
     -          s(i2,i3)*s(i3,i4)*za(i4,i5)*zb(i3,i1)*zb(i6,i5)))/
     -      (4._dp*s(i2,i3)*za(i1,i3)*zb(i2,i1))) + 
     -  c_alphai(:)*((ha*(s(i4,i5) + s(i4,i6))*za(i1,i2)*za(i2,i4)*
     -        za(i4,i5)*zb(i6,i1))/(2._dp*za(i1,i3)*za(i2,i3)) + 
     -     (hb*(s(i4,i5) + s(i4,i6))*za(i4,i5)*
     -        (-(s(i1,i2)*s(i4,i5)*za(i2,i4)*zb(i3,i2)*zb(i6,i1)) - 
     -          s(i1,i2)*s(i4,i6)*za(i2,i4)*zb(i3,i2)*zb(i6,i1) - 
     -          s(i1,i2)*za(i2,i4)*za(i4,i5)*zb(i3,i2)*zb(i4,i1)*
     -           zb(i6,i5)))/(4._dp*s(i2,i3)*za(i1,i3)*zb(i2,i1))) + 
     -  c_betai(:)*(-(ha*(s(i4,i5) + s(i4,i6))*za(i2,i4)*za(i4,i5)*
     -         zb(i6,i3))/(2._dp*za(i1,i3)) + 
     -     (hb*(s(i4,i5) + s(i4,i6))*za(i4,i5)*
     -        (s(i2,i3)*s(i4,i5)*za(i2,i4)*zb(i2,i1)*zb(i6,i3) + 
     -          s(i2,i3)*s(i4,i6)*za(i2,i4)*zb(i2,i1)*zb(i6,i3) + 
     -          s(i2,i3)*za(i2,i4)*za(i4,i5)*zb(i2,i1)*zb(i4,i3)*
     -           zb(i6,i5)))/(4._dp*s(i2,i3)*za(i1,i3)*zb(i2,i1))) + 
     -  c_gammai(:)*((ha*(s(i4,i5) + s(i4,i6))*za(i1,i4)*za(i4,i5)*
     -        zb(i3,i1)*zb(i6,i3))/(2._dp*za(i1,i3)*zb(i3,i2)) + 
     -     (hb*(s(i4,i5) + s(i4,i6))*za(i4,i5)*
     -        (-(s(i4,i5)*za(i1,i4)*za(i2,i3)*zb(i2,i1)*zb(i3,i1)*
     -             zb(i6,i3)) - 
     -          s(i4,i6)*za(i1,i4)*za(i2,i3)*zb(i2,i1)*zb(i3,i1)*
     -           zb(i6,i3) - 
     -          za(i1,i4)*za(i2,i3)*za(i4,i5)*zb(i2,i1)*zb(i3,i1)*
     -           zb(i4,i3)*zb(i6,i5)))/(4._dp*s(i2,i3)*za(i1,i3)*zb(i2,i1))
     -     )

          ! scale evolution
          amp_zaj_anomZZ_pm(1) = amp_zaj_anomZZ_pm(1) + 
     &        beta0/2*log(musq/t(i1,i2,i3)) * amp_zaj_anomZZ_pm(0)
        end function

        function amp_zaj_anomZa_pp(i1,i2,i3,i4,i5,i6, za,zb, flip,
     &                c_alphai, c_betai, c_gammai)
          include 'types.f'
          include 'mxpart.f'
          include 'scale.f'
          include 'anomcoup.f'

          complex(dp) :: amp_zaj_anomZa_pp(0:1)
          integer, intent(in) :: i1,i2,i3,i4,i5,i6
          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          logical, intent(in) :: flip
          complex(dp), intent(in) :: c_alphai(0:1), c_betai(0:1), c_gammai(0:1)

          complex(dp) :: c_deltai(0:1)
          complex(dp) :: ha, hb
          real(dp) :: hSign

          real(dp) :: s,t
          s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
          t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

          c_deltai(:) = (c_alphai(:) - c_betai(:) - c_gammai(:))
     &                      * s(i1,i2)/t(i1,i2,i3)/2._dp

          hsign = 1._dp
          if (flip) then
            hSign = -1._dp
          endif

          ha = hSign*im*hitGam(1) + hitGam(3)
          hb = hSign*im*hitGam(2) + hitGam(4)

          amp_zaj_anomZa_pp(:) = 
     -  (c_deltai(:)*s(i5,i6)*(2*ha*s(i1,i4)*s(i2,i3)*za(i5,i6)*
     -        zb(i3,i1) + 2*ha*s(i2,i3)*s(i2,i4)*za(i5,i6)*
     -        zb(i3,i1) + 2*ha*s(i2,i3)*s(i3,i4)*za(i5,i6)*
     -        zb(i3,i1) - hb*s(i1,i5)*s(i2,i3)*s(i4,i6)*za(i5,i6)*
     -        zb(i3,i1) - hb*s(i1,i6)*s(i2,i3)*s(i4,i6)*za(i5,i6)*
     -        zb(i3,i1) - hb*s(i2,i3)*s(i2,i5)*s(i4,i6)*za(i5,i6)*
     -        zb(i3,i1) - hb*s(i2,i3)*s(i2,i6)*s(i4,i6)*za(i5,i6)*
     -        zb(i3,i1) - hb*s(i2,i3)*s(i3,i5)*s(i4,i6)*za(i5,i6)*
     -        zb(i3,i1) - hb*s(i2,i3)*s(i3,i6)*s(i4,i6)*za(i5,i6)*
     -        zb(i3,i1) + 2*ha*s(i2,i3)*za(i1,i6)*za(i4,i5)*
     -        zb(i3,i1)*zb(i4,i1) + 
     -       2*ha*s(i2,i3)*za(i2,i6)*za(i4,i5)*zb(i3,i1)*
     -        zb(i4,i2) + 2*ha*s(i2,i3)*za(i3,i6)*za(i4,i5)*
     -        zb(i3,i1)*zb(i4,i3))*zb(i6,i4))/
     -   (4._dp*s(i2,i3)*za(i1,i3)*za(i4,i6)*zb(i2,i1)) + 
     -  (c_alphai(:)*s(i5,i6)*za(i1,i2)*
     -     (-2*ha*za(i2,i4)*za(i5,i6)*zb(i4,i1) + 
     -       hb*s(i4,i6)*za(i2,i5)*za(i5,i6)*zb(i5,i1) + 
     -       za(i2,i6)*(-2*ha*za(i4,i5)*zb(i4,i1) + 
     -          hb*s(i4,i6)*za(i5,i6)*zb(i6,i1)))*zb(i6,i4))/
     -   (4._dp*za(i1,i3)*za(i2,i3)*za(i4,i6)) + 
     -  (c_betai(:)*s(i5,i6)*(2*ha*s(i2,i3)*za(i2,i6)*za(i4,i5)*zb(i2,i1)*
     -        zb(i4,i3) + 2*ha*s(i2,i3)*za(i2,i4)*za(i5,i6)*
     -        zb(i2,i1)*zb(i4,i3) - 
     -       hb*s(i2,i3)*s(i4,i6)*za(i2,i5)*za(i5,i6)*zb(i2,i1)*
     -        zb(i5,i3) - hb*s(i2,i3)*s(i4,i6)*za(i2,i6)*za(i5,i6)*
     -        zb(i2,i1)*zb(i6,i3))*zb(i6,i4))/
     -   (4._dp*s(i2,i3)*za(i1,i3)*za(i4,i6)*zb(i2,i1)) + 
     -  (c_gammai(:)*s(i5,i6)*(-2*ha*za(i1,i6)*za(i2,i3)*za(i4,i5)*
     -        zb(i2,i1)*zb(i3,i1)*zb(i4,i3) - 
     -       2*ha*za(i1,i4)*za(i2,i3)*za(i5,i6)*zb(i2,i1)*zb(i3,i1)*
     -        zb(i4,i3) + hb*s(i4,i6)*za(i1,i5)*za(i2,i3)*za(i5,i6)*
     -        zb(i2,i1)*zb(i3,i1)*zb(i5,i3) + 
     -       hb*s(i4,i6)*za(i1,i6)*za(i2,i3)*za(i5,i6)*zb(i2,i1)*
     -        zb(i3,i1)*zb(i6,i3))*zb(i6,i4))/
     -   (4._dp*s(i2,i3)*za(i1,i3)*za(i4,i6)*zb(i2,i1))

          amp_zaj_anomZa_pp(:) = -amp_zaj_anomZa_pp(:)

          ! scale evolution
          amp_zaj_anomZa_pp(1) = amp_zaj_anomZa_pp(1) + 
     &        beta0/2*log(musq/t(i1,i2,i3)) * amp_zaj_anomZa_pp(0)
        end function

        function amp_zaj_anomZa_pm(i1,i2,i3,i4,i5,i6, za,zb, flip,
     &                c_alphai, c_betai, c_gammai)
          include 'types.f'
          include 'mxpart.f'
          include 'scale.f'
          include 'anomcoup.f'

          complex(dp) :: amp_zaj_anomZa_pm(0:1)
          integer, intent(in) :: i1,i2,i3,i4,i5,i6
          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          logical, intent(in) :: flip
          complex(dp), intent(in) :: c_alphai(0:1), c_betai(0:1), c_gammai(0:1)

          complex(dp) :: c_deltai(0:1)
          complex(dp) :: ha, hb
          real(dp) :: hSign

          real(dp) :: s,t
          s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
          t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

          c_deltai(:) = (c_alphai(:) - c_betai(:) - c_gammai(:))
     &                      * s(i1,i2)/t(i1,i2,i3)/2._dp

          hsign = 1._dp
          if (flip) then
            hSign = -1._dp
          endif

          ha = -hSign*im*hitGam(1) + hitGam(3)
          hb = -hSign*im*hitGam(2) + hitGam(4)

          amp_zaj_anomZa_pm(:) = 
     -  -(c_deltai(:)*s(i5,i6)*za(i4,i5)*
     -      (-2*ha*s(i2,i3)*za(i1,i4)*zb(i3,i1)*zb(i6,i1) - 
     -        2*ha*s(i2,i3)*za(i2,i4)*zb(i3,i1)*zb(i6,i2) - 
     -        2*ha*s(i2,i3)*za(i3,i4)*zb(i3,i1)*zb(i6,i3) - 
     -        hb*s(i1,i5)*s(i2,i3)*za(i4,i5)*zb(i3,i1)*zb(i6,i5) - 
     -        hb*s(i1,i6)*s(i2,i3)*za(i4,i5)*zb(i3,i1)*zb(i6,i5) - 
     -        hb*s(i2,i3)*s(i2,i5)*za(i4,i5)*zb(i3,i1)*zb(i6,i5) - 
     -        hb*s(i2,i3)*s(i2,i6)*za(i4,i5)*zb(i3,i1)*zb(i6,i5) - 
     -        hb*s(i2,i3)*s(i3,i5)*za(i4,i5)*zb(i3,i1)*zb(i6,i5) - 
     -        hb*s(i2,i3)*s(i3,i6)*za(i4,i5)*zb(i3,i1)*zb(i6,i5)))/
     -   (4._dp*s(i2,i3)*za(i1,i3)*zb(i2,i1)) - 
     -  (c_alphai(:)*s(i5,i6)*za(i1,i2)*za(i4,i5)*
     -     (2*ha*za(i2,i4)*zb(i6,i1) + 
     -       hb*za(i4,i5)*
     -        (za(i2,i5)*zb(i5,i1) + za(i2,i6)*zb(i6,i1))*zb(i6,i5)))/
     -   (4._dp*za(i1,i3)*za(i2,i3)) - 
     -  (c_betai(:)*s(i5,i6)*za(i4,i5)*
     -     (-2*ha*s(i2,i3)*za(i2,i4)*zb(i2,i1)*zb(i6,i3) - 
     -       hb*s(i2,i3)*za(i2,i5)*za(i4,i5)*zb(i2,i1)*zb(i5,i3)*
     -        zb(i6,i5) - hb*s(i2,i3)*za(i2,i6)*za(i4,i5)*zb(i2,i1)*
     -        zb(i6,i3)*zb(i6,i5)))/(4._dp*s(i2,i3)*za(i1,i3)*zb(i2,i1)) - 
     -  (c_gammai(:)*s(i5,i6)*za(i4,i5)*
     -     (2*ha*za(i1,i4)*za(i2,i3)*zb(i2,i1)*zb(i3,i1)*
     -        zb(i6,i3) + hb*za(i1,i5)*za(i2,i3)*za(i4,i5)*
     -        zb(i2,i1)*zb(i3,i1)*zb(i5,i3)*zb(i6,i5) + 
     -       hb*za(i1,i6)*za(i2,i3)*za(i4,i5)*zb(i2,i1)*zb(i3,i1)*
     -        zb(i6,i3)*zb(i6,i5)))/(4._dp*s(i2,i3)*za(i1,i3)*zb(i2,i1))

          ! scale evolution
          amp_zaj_anomZa_pm(1) = amp_zaj_anomZa_pm(1) + 
     &        beta0/2*log(musq/t(i1,i2,i3)) * amp_zaj_anomZa_pm(0)
        end function


        function amp_zaj_anomaZ_pp(i1,i2,i3,i4,i5,i6, za,zb, flip,
     &                c_alphai, c_betai, c_gammai)
          include 'types.f'
          include 'mxpart.f'
          include 'scale.f'
          include 'anomcoup.f'

          complex(dp) :: amp_zaj_anomaZ_pp(0:1)
          integer, intent(in) :: i1,i2,i3,i4,i5,i6
          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          logical, intent(in) :: flip
          complex(dp), intent(in) :: c_alphai(0:1), c_betai(0:1), c_gammai(0:1)

          complex(dp) :: c_deltai(0:1)
          complex(dp) :: ha, hb
          real(dp) :: hSign

          real(dp) :: s,t
          s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
          t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

          c_deltai(:) = (c_alphai(:) - c_betai(:) - c_gammai(:))
     &                      * s(i1,i2)/t(i1,i2,i3)/2._dp

          hsign = 1._dp
          if (flip) then
            hSign = -1._dp
          endif

          ha = hSign*im*hitGam(1) + hitGam(3)
          hb = hSign*im*hitGam(2) + hitGam(4)

          amp_zaj_anomaZ_pp(:) = 
     -  (c_alphai(:)*(s(i4,i5) + s(i4,i6) + s(i5,i6))*za(i1,i2)*
     -     (-((-2*ha + hb*s(i4,i5) + hb*s(i4,i6))*za(i2,i6)*
     -          za(i4,i5)) + 
     -       (2*ha - hb*s(i4,i5))*za(i2,i4)*za(i5,i6))*
     -     zb(i4,i1)*zb(i6,i4))/(4._dp*za(i1,i3)*za(i2,i3)*za(i4,i6)) - 
     -  (c_betai(:)*(s(i4,i5) + s(i4,i6) + s(i5,i6))*
     -     (-(s(i2,i3)*(-2*ha + hb*s(i4,i5) + hb*s(i4,i6))*
     -          za(i2,i6)*za(i4,i5)*zb(i2,i1)*zb(i4,i3)) + 
     -       s(i2,i3)*(2*ha - hb*s(i4,i5))*za(i2,i4)*za(i5,i6)*
     -        zb(i2,i1)*zb(i4,i3))*zb(i6,i4))/
     -   (4._dp*s(i2,i3)*za(i1,i3)*za(i4,i6)*zb(i2,i1)) - 
     -  (c_gammai(:)*(s(i4,i5) + s(i4,i6) + s(i5,i6))*
     -     ((-2*ha + hb*s(i4,i5) + hb*s(i4,i6))*za(i1,i6)*
     -        za(i2,i3)*za(i4,i5)*zb(i2,i1)*zb(i3,i1)*zb(i4,i3) - 
     -       (2*ha - hb*s(i4,i5))*za(i1,i4)*za(i2,i3)*za(i5,i6)*
     -        zb(i2,i1)*zb(i3,i1)*zb(i4,i3))*zb(i6,i4))/
     -   (4._dp*s(i2,i3)*za(i1,i3)*za(i4,i6)*zb(i2,i1)) - 
     -  (c_deltai(:)*(s(i4,i5) + s(i4,i6) + s(i5,i6))*
     -     (s(i1,i4)*s(i2,i3)*(2*ha - hb*s(i4,i5))*za(i5,i6)*
     -        zb(i3,i1) + s(i2,i3)*s(i2,i4)*
     -        (2*ha - hb*s(i4,i5))*za(i5,i6)*zb(i3,i1) - 
     -       s(i2,i3)*(-2*ha + hb*s(i4,i5) + hb*s(i4,i6))*
     -        za(i1,i6)*za(i4,i5)*zb(i3,i1)*zb(i4,i1) - 
     -       s(i2,i3)*(-2*ha + hb*s(i4,i5) + hb*s(i4,i6))*
     -        za(i2,i6)*za(i4,i5)*zb(i3,i1)*zb(i4,i2) - 
     -       s(i2,i3)*zb(i3,i1)*zb(i4,i3)*
     -        (hb*s(i4,i5)*za(i3,i6)*za(i4,i5) - 
     -          2*ha*za(i3,i4)*za(i5,i6) - 
     -          za(i4,i5)*(-((-2*ha + hb*s(i4,i6))*za(i3,i6)) - 
     -             hb*za(i3,i4)*za(i5,i6)*zb(i5,i4))))*zb(i6,i4))/
     -   (4._dp*s(i2,i3)*za(i1,i3)*za(i4,i6)*zb(i2,i1))

          amp_zaj_anomaZ_pp(:) = -amp_zaj_anomaZ_pp(:)

          ! scale evolution
          amp_zaj_anomaZ_pp(1) = amp_zaj_anomaZ_pp(1) + 
     &        beta0/2*log(musq/t(i1,i2,i3)) * amp_zaj_anomaZ_pp(0)
        end function

        function amp_zaj_anomaZ_pm(i1,i2,i3,i4,i5,i6, za,zb, flip,
     &                c_alphai, c_betai, c_gammai)
          include 'types.f'
          include 'mxpart.f'
          include 'scale.f'
          include 'anomcoup.f'

          complex(dp) :: amp_zaj_anomaZ_pm(0:1)
          integer, intent(in) :: i1,i2,i3,i4,i5,i6
          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          logical, intent(in) :: flip
          complex(dp), intent(in) :: c_alphai(0:1), c_betai(0:1), c_gammai(0:1)

          complex(dp) :: c_deltai(0:1)
          complex(dp) :: ha, hb
          real(dp) :: hSign

          real(dp) :: s,t
          s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
          t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

          c_deltai(:) = (c_alphai(:) - c_betai(:) - c_gammai(:))
     &                      * s(i1,i2)/t(i1,i2,i3)/2._dp

          hsign = 1._dp
          if (flip) then
            hSign = -1._dp
          endif

          ha = -hSign*im*hitGam(1) + hitGam(3)
          hb = -hSign*im*hitGam(2) + hitGam(4)

          amp_zaj_anomaZ_pm(:) = 
     -  (c_alphai(:)*(s(i4,i5) + s(i4,i6) + s(i5,i6))*za(i1,i2)*za(i2,i4)*
     -     za(i4,i5)*((2*ha - hb*s(i4,i5) - hb*s(i4,i6))*
     -        zb(i6,i1) - hb*za(i4,i5)*zb(i4,i1)*zb(i6,i5)))/
     -   (4._dp*za(i1,i3)*za(i2,i3)) + 
     -  (c_betai(:)*(s(i4,i5) + s(i4,i6) + s(i5,i6))*za(i4,i5)*
     -     (hb*s(i2,i3)*s(i4,i5)*za(i2,i4)*zb(i2,i1)*zb(i6,i3) + 
     -       s(i2,i3)*(-2*ha + hb*s(i4,i6))*za(i2,i4)*zb(i2,i1)*
     -        zb(i6,i3) + hb*s(i2,i3)*za(i2,i4)*za(i4,i5)*zb(i2,i1)*
     -        zb(i4,i3)*zb(i6,i5)))/(4._dp*s(i2,i3)*za(i1,i3)*zb(i2,i1)) + 
     -  (c_gammai(:)*(s(i4,i5) + s(i4,i6) + s(i5,i6))*za(i4,i5)*
     -     (-(hb*s(i4,i5)*za(i1,i4)*za(i2,i3)*zb(i2,i1)*zb(i3,i1)*
     -          zb(i6,i3)) + 
     -       (2*ha - hb*s(i4,i6))*za(i1,i4)*za(i2,i3)*zb(i2,i1)*
     -        zb(i3,i1)*zb(i6,i3) - 
     -       hb*za(i1,i4)*za(i2,i3)*za(i4,i5)*zb(i2,i1)*zb(i3,i1)*
     -        zb(i4,i3)*zb(i6,i5)))/(4._dp*s(i2,i3)*za(i1,i3)*zb(i2,i1)) + 
     -  (c_deltai(:)*(s(i4,i5) + s(i4,i6) + s(i5,i6))*za(i4,i5)*
     -     (s(i2,i3)*(-2*ha + hb*s(i4,i5) + hb*s(i4,i6))*
     -        za(i1,i4)*zb(i3,i1)*zb(i6,i1) + 
     -       s(i2,i3)*(-2*ha + hb*s(i4,i5) + hb*s(i4,i6))*
     -        za(i2,i4)*zb(i3,i1)*zb(i6,i2) + 
     -       hb*s(i1,i4)*s(i2,i3)*za(i4,i5)*zb(i3,i1)*zb(i6,i5) + 
     -       hb*s(i2,i3)*s(i2,i4)*za(i4,i5)*zb(i3,i1)*zb(i6,i5) - 
     -       s(i2,i3)*za(i3,i4)*zb(i3,i1)*
     -        (-((-2*ha + hb*s(i4,i5) + hb*s(i4,i6))*
     -             zb(i6,i3)) - hb*za(i4,i5)*zb(i4,i3)*zb(i6,i5))))/
     -   (4._dp*s(i2,i3)*za(i1,i3)*zb(i2,i1))

          ! scale evolution
          amp_zaj_anomaZ_pm(1) = amp_zaj_anomaZ_pm(1) + 
     &        beta0/2*log(musq/t(i1,i2,i3)) * amp_zaj_anomaZ_pm(0)
        end function


        ! since amp_zaj_fsr_* returns the UV+IR finite amplitude, this routine
        ! provides the terms to undo the subtractions for debugging comparison
        ! with the old MCFM routines:
        ! (I1 +  beta0/2/ep - CA/6/2)
        ! The latter term, associated with a finite UV renormalization term,
        ! is added manually in the old routines, and not included in the
        ! amplitudes itself.
        function mcfm_subtraction(i1,i2,i3,i4,i5,i6, za,zb)
          include 'types.f'
          include 'mxpart.f'
          include 'epinv.f'
          include 'constants.f'
          include 'nf.f'
          include 'scale.f'

          complex(dp) :: mcfm_subtraction
          integer, intent(in) :: i1,i2,i3,i4,i5,i6
          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

          complex(dp) :: lnrat

          real(dp) :: s,t
          s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
          t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

          mcfm_subtraction = 0._dp

          mcfm_subtraction = mcfm_subtraction + epinv**2*(1/(2._dp*Nc) - Nc)

          mcfm_subtraction = mcfm_subtraction + epinv*(
     -  (3 - 3*Nc**2 + 2*lnrat(musq,-s(i1,i2)) - 
     -    2*Nc**2*lnrat(musq,-s(i1,i3)) - 2*Nc**2*lnrat(musq,-s(i2,i3)))
     -   /(4._dp*Nc)
     -    )

          mcfm_subtraction = mcfm_subtraction + (
     -  (-2*CA*Nc + 18*lnrat(musq,-s(i1,i2)) + 
     -    6*lnrat(musq,-s(i1,i2))**2 - 
     -    6*beta0*Nc*lnrat(musq,-s(i1,i3)) - 
     -    9*Nc**2*lnrat(musq,-s(i1,i3)) - 
     -    6*Nc**2*lnrat(musq,-s(i1,i3))**2 - 
     -    6*beta0*Nc*lnrat(musq,-s(i2,i3)) - 
     -    9*Nc**2*lnrat(musq,-s(i2,i3)) - 
     -    6*Nc**2*lnrat(musq,-s(i2,i3))**2)/(24._dp*Nc)
     -    )

          ! tH-V to dred scheme conversion
          mcfm_subtraction = mcfm_subtraction + (CF + CA/6)/2

        end function

        subroutine omega_ee3jet(p, i1,i2,i3, c_alphai, c_betai, c_gammai)
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'nf.f'
          include 'zeta.f'
          include 'tiny.f'

          real(dp), intent(in) :: p(mxpart,4)
          integer, intent(in) :: i1,i2,i3
          ! these arrays could be extended to 0:2 to include the two loop results
          complex(dp), intent(out) :: c_alphai(0:1), c_betai(0:1), c_gammai(0:1)

          real(dp) :: dotvec

          real(dp) :: s,t
          s(i1,i2) = dotvec(p(i1,:) + p(i2,:), p(i1,:) + p(i2,:))
          t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

          real(dp) :: s12, s13, s23, s123
          complex(dp) :: w(36)
          real(dp) :: y,z,yPz,oneMz
          real(dp) :: n
          real(dp) :: li2, xu, xv

          s12 = s(i1,i2)
          s13 = s(i1,i3)
          s23 = s(i2,i3)
          s123 = s12 + s13 + s23

          n = Nc

          c_alphai(:) = 0._dp
          c_betai(:) = 0._dp
          c_gammai(:) = 0._dp

          y = s13/s123
          z = s23/s123
          yPz = (s123-s12)/s123
          oneMz = (s123-s23)/s123

          if( i3 == 3 ) then
            ! region (2a)+
            xu = -s13/s12
            xv = s123/s12
          elseif ( i2 == 3) then
            ! region (3a)+
            xu = -s23/s13
            xv = s123/s13
          elseif ( i1 == 3) then
            ! region (4a)+
            xu = -s13/s23
            xv = s123/s23
          else
            call abort
          endif

          if (xu < 0._dp .or. xu > 1-xv .or. xv < 0._dp .or. xv > 1) then
            if ( abs(xv) < 100*tiny .or. abs(xu) < 100*tiny ) then
              xv = abs(xv)
              xu = abs(xu)
            elseif (xu+xv > 1._dp + 100*tiny) then
              ! no fixes necessary
              write (*,*) "broken point!", "u, v, u+v", xu, xv, xu+xv
            else
              ! for now no further fixes until problems arise
              write (*,*) "broken point!", "u, v, u+v", xu, xv, xu+xv
            endif
          endif

          if( i3 == 3 ) then
            ! region (2a)+
            include 'omega_ee3jet_2a_inc.f'
          elseif ( i2 == 3) then
            ! region (3a)+
            include 'omega_ee3jet_3a_inc.f'
          elseif ( i1 == 3) then
            ! region (4a)+
            include 'omega_ee3jet_4a_inc.f'
          else
            call abort
          endif


        end subroutine

      ! This calculates all crossings just as zaj_crossings_l in zaj_utils.f,
      ! but for the newly calculated amplitudes based on the ee3jet code.
      ! It is modified to compute the tensor coefficients alpha,beta,gamma
      ! for each crossing of the quarks, since it only depends on i1,i2,i3
      subroutine zaj_crossings_l_new(p, i1,i2,i3,i4,i5,i6,za,zb,amps)
        use helicities 
        implicit none
        include 'types.f'
        include 'constants.f'
        include 'mxpart.f'

        real(dp), intent(in) :: p(mxpart,4)
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        ! helicities for particles q, l, g, gam
        complex(dp), intent(out) :: amps(0:1, 2,2,2,2)

        complex(dp) :: mcfmSub
        complex(dp) :: c_alphai_12(0:1), c_betai_12(0:1), c_gammai_12(0:1)
        complex(dp) :: c_alphai_21(0:1), c_betai_21(0:1), c_gammai_21(0:1)

        mcfmSub = mcfm_subtraction(i1,i2,i3,i4,i5,i6,za,zb)
        call omega_ee3jet(p, i1,i2,i3, c_alphai_12, c_betai_12, c_gammai_12)
        call omega_ee3jet(p, i2,i1,i3, c_alphai_21, c_betai_21, c_gammai_21)

        ! first block q=+, l=+
        amps(:, helPlus, helPlus, helPlus, helPlus) = 
     &    amp_zaj_fsr_pp(i1,i2,i3,i4,i5,i6,za,zb,c_alphai_12,c_betai_12,c_gammai_12)
        amps(:, helPlus, helPlus, helPlus, helMinus) =
     &    amp_zaj_fsr_pm(i1,i2,i3,i4,i5,i6,za,zb,c_alphai_12,c_betai_12,c_gammai_12)
        amps(:, helPlus, helPlus, helMinus, helPlus) =
     &    amp_zaj_fsr_pm(i2,i1,i3,i4,i6,i5,zb,za,c_alphai_21,c_betai_21,c_gammai_21)
        amps(:, helPlus, helPlus, helMinus, helMinus) =
     &    amp_zaj_fsr_pp(i2,i1,i3,i4,i6,i5,zb,za,c_alphai_21,c_betai_21,c_gammai_21)

c       ! second block, q=+, l=-, in principle with a minus sign
        amps(:,helPlus, helMinus, helPlus, helPlus) = 
     &    amp_zaj_fsr_pp(i1,i2,i3,i4,i6,i5,za,zb,c_alphai_12,c_betai_12,c_gammai_12)
        amps(:,helPlus, helMinus, helPlus, helMinus) =
     &    amp_zaj_fsr_pm(i1,i2,i3,i4,i6,i5,za,zb,c_alphai_12,c_betai_12,c_gammai_12)
        amps(:,helPlus, helMinus, helMinus, helPlus) =
     &    amp_zaj_fsr_pm(i2,i1,i3,i4,i5,i6,zb,za,c_alphai_21,c_betai_21,c_gammai_21)
        amps(:,helPlus, helMinus, helMinus, helMinus) = 
     &    amp_zaj_fsr_pp(i2,i1,i3,i4,i5,i6,zb,za,c_alphai_21,c_betai_21,c_gammai_21)

c       ! third block q=-, l=-
        amps(:,helMinus, helMinus, helPlus, helPlus) = 
     &    amp_zaj_fsr_pp(i2,i1,i3,i4,i6,i5,za,zb,c_alphai_21,c_betai_21,c_gammai_21)
        amps(:,helMinus, helMinus, helPlus, helMinus) =
     &    amp_zaj_fsr_pm(i2,i1,i3,i4,i6,i5,za,zb,c_alphai_21,c_betai_21,c_gammai_21)
        amps(:,helMinus, helMinus, helMinus, helPlus) =
     &    amp_zaj_fsr_pm(i1,i2,i3,i4,i5,i6,zb,za,c_alphai_12,c_betai_12,c_gammai_12)
        amps(:,helMinus, helMinus, helMinus, helMinus) =
     &    amp_zaj_fsr_pp(i1,i2,i3,i4,i5,i6,zb,za,c_alphai_12,c_betai_12,c_gammai_12)

c       ! fourth block q=-, l=+, in principle with a minus sign
        amps(:,helMinus, helPlus, helPlus, helPlus) =
     &    amp_zaj_fsr_pp(i2,i1,i3,i4,i5,i6,za,zb,c_alphai_21,c_betai_21,c_gammai_21)
        amps(:,helMinus, helPlus, helPlus, helMinus) = 
     &    amp_zaj_fsr_pm(i2,i1,i3,i4,i5,i6,za,zb,c_alphai_21,c_betai_21,c_gammai_21)
        amps(:,helMinus, helPlus, helMinus, helPlus) =
     &    amp_zaj_fsr_pm(i1,i2,i3,i4,i6,i5,zb,za,c_alphai_12,c_betai_12,c_gammai_12)
        amps(:,helMinus, helPlus, helMinus, helMinus) =
     &    amp_zaj_fsr_pp(i1,i2,i3,i4,i6,i5,zb,za,c_alphai_12,c_betai_12,c_gammai_12)

        ! conversion to IR+UV unsubtracted "scheme" in MCFM
        ! including conversion from tH-V to dred, such that the returned amp
        ! matches exactly the zaj_virtamp_l1 result
        amps(1, :,:,:,:) = (amps(1,:,:,:,:) +  mcfmSub*amps(0, :,:,:,:))*2

        ! I believe that the old mcfm loop amplitudes are off by a factor of two.
        ! This factor is removed again at the end of the qqb_zaj_v and
        ! qqb_zaj_v_new routines. For the sake of comparing with the old amplitudes
        ! I will just replicate this factor of two here.

        ! In principle we could set the scheme to tH-V, don't add mcfmSub
        ! and define a qqb_zaj_z_new which sets Q1 and Q2 to zero.
     
      end subroutine

      ! and the same again for the anomalous coupling amplitudes which require
      ! a manual sign flip on za <-> zb
      subroutine zaj_crossings_anom_new(p, i1,i2,i3,i4,i5,i6,za,zb, ampPP, ampPM, amps)
        use helicities 
        implicit none
        include 'types.f'
        include 'constants.f'
        include 'mxpart.f'

        real(dp), intent(in) :: p(mxpart,4)
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        ! helicities for particles q, l, g, gam
        complex(dp), intent(out) :: amps(0:1, 2,2,2,2)

        procedure (anomAmp) :: ampPP
        procedure (anomAmp) :: ampPM

        complex(dp) :: mcfmSub
        complex(dp) :: c_alphai_12(0:1), c_betai_12(0:1), c_gammai_12(0:1)
        complex(dp) :: c_alphai_21(0:1), c_betai_21(0:1), c_gammai_21(0:1)


        mcfmSub = mcfm_subtraction(i1,i2,i3,i4,i5,i6,za,zb)
        call omega_ee3jet(p, i1,i2,i3, c_alphai_12, c_betai_12, c_gammai_12)
        call omega_ee3jet(p, i2,i1,i3, c_alphai_21, c_betai_21, c_gammai_21)

        ! first block q=+, l=+
        amps(:, helPlus, helPlus, helPlus, helPlus) = 
     &    ampPP(i1,i2,i3,i4,i5,i6,za,zb,.false.,c_alphai_12,c_betai_12,c_gammai_12)
        amps(:, helPlus, helPlus, helPlus, helMinus) =
     &    ampPM(i1,i2,i3,i4,i5,i6,za,zb,.false.,c_alphai_12,c_betai_12,c_gammai_12)
        amps(:, helPlus, helPlus, helMinus, helPlus) =
     &    ampPM(i2,i1,i3,i4,i6,i5,zb,za,.true.,c_alphai_21,c_betai_21,c_gammai_21)
        amps(:, helPlus, helPlus, helMinus, helMinus) =
     &    ampPP(i2,i1,i3,i4,i6,i5,zb,za,.true.,c_alphai_21,c_betai_21,c_gammai_21)

c       ! second block, q=+, l=-, in principle with a minus sign
        amps(:,helPlus, helMinus, helPlus, helPlus) = 
     &    ampPP(i1,i2,i3,i4,i6,i5,za,zb,.false.,c_alphai_12,c_betai_12,c_gammai_12)
        amps(:,helPlus, helMinus, helPlus, helMinus) =
     &    ampPM(i1,i2,i3,i4,i6,i5,za,zb,.false.,c_alphai_12,c_betai_12,c_gammai_12)
        amps(:,helPlus, helMinus, helMinus, helPlus) =
     &    ampPM(i2,i1,i3,i4,i5,i6,zb,za,.true.,c_alphai_21,c_betai_21,c_gammai_21)
        amps(:,helPlus, helMinus, helMinus, helMinus) = 
     &    ampPP(i2,i1,i3,i4,i5,i6,zb,za,.true.,c_alphai_21,c_betai_21,c_gammai_21)

c       ! third block q=-, l=-
        amps(:,helMinus, helMinus, helPlus, helPlus) = 
     &    ampPP(i2,i1,i3,i4,i6,i5,za,zb,.false.,c_alphai_21,c_betai_21,c_gammai_21)
        amps(:,helMinus, helMinus, helPlus, helMinus) =
     &    ampPM(i2,i1,i3,i4,i6,i5,za,zb,.false.,c_alphai_21,c_betai_21,c_gammai_21)
        amps(:,helMinus, helMinus, helMinus, helPlus) =
     &    ampPM(i1,i2,i3,i4,i5,i6,zb,za,.true.,c_alphai_12,c_betai_12,c_gammai_12)
        amps(:,helMinus, helMinus, helMinus, helMinus) =
     &    ampPP(i1,i2,i3,i4,i5,i6,zb,za,.true.,c_alphai_12,c_betai_12,c_gammai_12)

c       ! fourth block q=-, l=+, in principle with a minus sign
        amps(:,helMinus, helPlus, helPlus, helPlus) =
     &    ampPP(i2,i1,i3,i4,i5,i6,za,zb,.false.,c_alphai_21,c_betai_21,c_gammai_21)
        amps(:,helMinus, helPlus, helPlus, helMinus) = 
     &    ampPM(i2,i1,i3,i4,i5,i6,za,zb,.false.,c_alphai_21,c_betai_21,c_gammai_21)
        amps(:,helMinus, helPlus, helMinus, helPlus) =
     &    ampPM(i1,i2,i3,i4,i6,i5,zb,za,.true.,c_alphai_12,c_betai_12,c_gammai_12)
        amps(:,helMinus, helPlus, helMinus, helMinus) =
     &    ampPP(i1,i2,i3,i4,i6,i5,zb,za,.true.,c_alphai_12,c_betai_12,c_gammai_12)

        amps(1, :,:,:,:) = (amps(1,:,:,:,:) +  mcfmSub*amps(0, :,:,:,:))*2
     
      end subroutine


      end module
