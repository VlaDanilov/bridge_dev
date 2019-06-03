      module singletop2_realamps_m

      public :: singletop2_amp_tree, singletop2_real, singletop2_hamp_tree
      public :: singletop2_real_light, singletop2_real_heavy
      public :: singletop2_nonrestree_c10, singletop2_nonrestree_c4

      private :: ubtdg_h_full, ubtdg_l_new

      ! calling it "protos" because it has been compared pointwise
      ! to the Protos code http://jaguilar.web.cern.ch/jaguilar/protos/
      public :: singletop2_amp_protos

      private

      contains

      ! non-resonant tree amplitude for C_4 operator, for renormalization of non-resonant
      ! C4 and C7 contributions
      subroutine singletop2_nonrestree_c4(za,zb, ju,jb,jn,je,jc,jd, amps, sc)
          use types
          use anomcoup_tbW
          use eftcouple
          use constants
          implicit none

          include 'nf.f'
          include 'mxpart.f'

          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          integer, intent(in) :: ju,jb,jn,je,jc,jd
          real(dp), intent(in) :: sc

          ! 1: PM, 2: MP
          complex(dp) :: amps(2)

          integer :: j,k,l
          real(dp) :: s,t
          s(j,k) = real(za(j,k)*zb(k,j))
          t(j,k,l) = s(j,k) + s(j,l) + s(k,l)

          complex(dp) :: prop16, prop34, prop25
          real(dp) :: s12, s15, s16, s25, s56, s14, s36
          real(dp) :: s26, s45, s46, s24, s13, s23, s34, s35

          prop16 = 1._dp / (s(ju,jd) - wmass**2*sc)
          prop34 = 1._dp / (s(jn,je) - wmass**2*sc + im*wmass*wwidth*sc)
          prop25 = 1._dp / (s(jb,jc) - zmass**2*sc)

          s12 = s(ju,jb)
          s13 = s(ju,jn)
          s14 = s(ju,je)
          s15 = s(ju,jc)
          s16 = s(ju,jd)
          s23 = s(jb,jn)
          s24 = s(jb,je)
          s25 = s(jb,jc)
          s26 = s(jb,jd)
          s34 = s(jn,je)
          s35 = s(jn,jc)
          s36 = s(jn,jd)
          s45 = s(je,jc)
          s46 = s(je,jd)
          s24 = s(jb,je)
          s56 = s(jc,jd)

          amps(1) = 
     -  (gw**4*prop34*zb(jc,jb)*(((gb**2 + 3*gw**2)*prop25*
     -         ((-2*s13 + s25 + s26 + 2*(s45 + s46) + s56)*za(jn,jd) - 
     -           2*za(jb,jd)*za(ju,jn)*zb(jb,ju) + 2*za(jb,jn)*za(je,jd)*zb(je,jb))*zb(je,ju))/
     -       (s25 + s26 + s56) + (2*gb**2*
     -         ((-2*s13 + s25 + s26 + 2*(s45 + s46) + s56)*za(jn,jd) - 
     -           2*za(jb,jd)*za(ju,jn)*zb(jb,ju) + 2*za(jb,jn)*za(je,jd)*zb(je,jb))*zb(je,ju))/
     -       (s25*(s25 + s26 + s56)) + 6*gw**2*prop16*prop25*
     -       (2*za(jb,jd)*zb(jb,ju)*(za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb)) + 
     -         (2*s12 + s16 + s25 + 2*s26 - s34)*za(jn,jd)*zb(je,ju) + 
     -         2*za(jb,jn)*zb(je,jb)*(za(je,jd)*zb(je,ju) + za(jn,jd)*zb(jn,ju))) + 
     -      (6*gb**2*prop16*(2*za(jb,jd)*zb(jb,ju)*(za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb)) + 
     -           (2*s12 + s16 + s25 + 2*s26 - s34)*za(jn,jd)*zb(je,ju) + 
     -           2*za(jb,jn)*zb(je,jb)*(za(je,jd)*zb(je,ju) + za(jn,jd)*zb(jn,ju))))/s25 + 
     -      ((gb**2 - 3*gw**2)*prop25*(2*za(jb,jd)*zb(jb,ju)*
     -            (za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb) + za(ju,jn)*zb(je,ju)) + 
     -           za(jn,jd)*((2*s12 + 2*s13 + s16 + s26 - s34 - 2*(s45 + s46) - s56)*zb(je,ju) + 
     -              2*za(jb,jn)*zb(je,jb)*zb(jn,ju))))/(s16 + s26 - s34 + s56) - 
     -      (4*gb**2*(2*za(jb,jd)*zb(jb,ju)*
     -            (za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb) + za(ju,jn)*zb(je,ju)) + 
     -           za(jn,jd)*((2*s12 + 2*s13 + s16 + s26 - s34 - 2*(s45 + s46) - s56)*zb(je,ju) + 
     -              2*za(jb,jn)*zb(je,jb)*zb(jn,ju))))/(s25*(s16 + s26 - s34 + s56))))/
     -  (3._dp*ecossin**2)



          amps(2) =
     -  (gw**4*prop34*za(jb,jc)*(2*gb**2*(s16 + s26 - s34 + s56)*
     -       ((2*s13 - s25 - s26 - 2*(s45 + s46) - s56)*za(jn,jd) + 
     -         2*za(jb,jd)*za(ju,jn)*zb(jb,ju) - 2*za(jb,jn)*za(je,jd)*zb(je,jb))*zb(je,ju) + 
     -      (gb**2 + 3*gw**2)*prop25*s25*(s16 + s26 - s34 + s56)*
     -       ((2*s13 - s25 - s26 - 2*(s45 + s46) - s56)*za(jn,jd) + 
     -         2*za(jb,jd)*za(ju,jn)*zb(jb,ju) - 2*za(jb,jn)*za(je,jd)*zb(je,jb))*zb(je,ju) - 
     -      6*gb**2*prop16*(s25 + s26 + s56)*(s16 + s26 - s34 + s56)*
     -       (2*za(jb,jd)*zb(jb,ju)*(za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb)) + 
     -         (2*s12 + s16 + s25 + 2*s26 - s34)*za(jn,jd)*zb(je,ju) + 
     -         2*za(jb,jn)*zb(je,jb)*(za(je,jd)*zb(je,ju) + za(jn,jd)*zb(jn,ju))) - 
     -      6*gw**2*prop16*prop25*s25*(s25 + s26 + s56)*(s16 + s26 - s34 + s56)*
     -       (2*za(jb,jd)*zb(jb,ju)*(za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb)) + 
     -         (2*s12 + s16 + s25 + 2*s26 - s34)*za(jn,jd)*zb(je,ju) + 
     -         2*za(jb,jn)*zb(je,jb)*(za(je,jd)*zb(je,ju) + za(jn,jd)*zb(jn,ju))) + 
     -      4*gb**2*(s25 + s26 + s56)*(2*za(jb,jd)*zb(jb,ju)*
     -          (za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb) + za(ju,jn)*zb(je,ju)) + 
     -         za(jn,jd)*((2*s12 + 2*s13 + s16 + s26 - s34 - 2*(s45 + s46) - s56)*zb(je,ju) + 
     -            2*za(jb,jn)*zb(je,jb)*zb(jn,ju))) - 
     -      (gb**2 - 3*gw**2)*prop25*s25*(s25 + s26 + s56)*
     -       (2*za(jb,jd)*zb(jb,ju)*(za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb) + 
     -            za(ju,jn)*zb(je,ju)) + 
     -         za(jn,jd)*((2*s12 + 2*s13 + s16 + s26 - s34 - 2*(s45 + s46) - s56)*zb(je,ju) + 
     -            2*za(jb,jn)*zb(je,jb)*zb(jn,ju)))))/
     -  (3._dp*ecossin**2*s25*(s25 + s26 + s56)*(s16 + s26 - s34 + s56))


      end subroutine 

      ! non-resonant tree amplitude for C_dB operator, for renormalization of non-resonant
      ! C7 contributions
      subroutine singletop2_nonrestree_c10(za,zb, ju,jb,jn,je,jc,jd, amps, sc)
          use types
          use anomcoup_tbW
          use eftcouple
          use constants
          implicit none

          include 'nf.f'
          include 'mxpart.f'

          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          integer, intent(in) :: ju,jb,jn,je,jc,jd
          real(dp), intent(in) :: sc

          ! 1: PM, 2: MP
          complex(dp) :: amps(2)

          integer :: j,k,l
          real(dp) :: s,t
          s(j,k) = real(za(j,k)*zb(k,j))
          t(j,k,l) = s(j,k) + s(j,l) + s(k,l)

          complex(dp) :: prop16, prop34, prop25
          real(dp) :: s12, s15, s16, s25, s56, s14, s36
          real(dp) :: s26, s45, s46, s24, s13, s23, s34, s35

          prop16 = 1._dp / (s(ju,jd) - wmass**2*sc)
          prop34 = 1._dp / (s(jn,je) - wmass**2*sc + im*wmass*wwidth*sc)
          prop25 = 1._dp / (s(jb,jc) - zmass**2*sc)

          s12 = s(ju,jb)
          s13 = s(ju,jn)
          s14 = s(ju,je)
          s15 = s(ju,jc)
          s16 = s(ju,jd)
          s23 = s(jb,jn)
          s24 = s(jb,je)
          s25 = s(jb,jc)
          s26 = s(jb,jd)
          s34 = s(jn,je)
          s35 = s(jn,jc)
          s36 = s(jn,jd)
          s45 = s(je,jc)
          s46 = s(je,jd)
          s24 = s(jb,je)
          s56 = s(jc,jd)

          amps(1) = 
     &  (gb**2*gw**2*prop34*zb(jc,jb)*
     -    (2*gw**2*(s16 + s26 - s34 + s56)*
     -       ((2*s13 - s25 - s26 - 2*(s45 + s46) - s56)*za(jn,jd) + 
     -         2*za(jb,jd)*za(ju,jn)*zb(jb,ju) - 2*za(jb,jn)*za(je,jd)*zb(je,jb))*
     -       zb(je,ju) - (gb**2 + 3*gw**2)*prop25*s25*(s16 + s26 - s34 + s56)*
     -       ((2*s13 - s25 - s26 - 2*(s45 + s46) - s56)*za(jn,jd) + 
     -         2*za(jb,jd)*za(ju,jn)*zb(jb,ju) - 2*za(jb,jn)*za(je,jd)*zb(je,jb))*
     -       zb(je,ju) - 6*gw**2*prop16*(s25 + s26 + s56)*(s16 + s26 - s34 + s56)*
     -       (2*za(jb,jd)*zb(jb,ju)*(za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb)) + 
     -         (2*s12 + s16 + s25 + 2*s26 - s34)*za(jn,jd)*zb(je,ju) + 
     -         2*za(jb,jn)*zb(je,jb)*(za(je,jd)*zb(je,ju) + za(jn,jd)*zb(jn,ju))) + 
     -      6*gw**2*prop16*prop25*s25*(s25 + s26 + s56)*(s16 + s26 - s34 + s56)*
     -       (2*za(jb,jd)*zb(jb,ju)*(za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb)) + 
     -         (2*s12 + s16 + s25 + 2*s26 - s34)*za(jn,jd)*zb(je,ju) + 
     -         2*za(jb,jn)*zb(je,jb)*(za(je,jd)*zb(je,ju) + za(jn,jd)*zb(jn,ju))) + 
     -      4*gw**2*(s25 + s26 + s56)*
     -       (2*za(jb,jd)*zb(jb,ju)*(za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb) + 
     -            za(ju,jn)*zb(je,ju)) + 
     -         za(jn,jd)*((2*s12 + 2*s13 + s16 + s26 - s34 - 2*(s45 + s46) - s56)*
     -             zb(je,ju) + 2*za(jb,jn)*zb(je,jb)*zb(jn,ju))) + 
     -      (gb**2 - 3*gw**2)*prop25*s25*(s25 + s26 + s56)*
     -       (2*za(jb,jd)*zb(jb,ju)*(za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb) + 
     -            za(ju,jn)*zb(je,ju)) + 
     -         za(jn,jd)*((2*s12 + 2*s13 + s16 + s26 - s34 - 2*(s45 + s46) - s56)*
     -             zb(je,ju) + 2*za(jb,jn)*zb(je,jb)*zb(jn,ju)))))/
     -  (3._dp*ecossin**2*s25*(s25 + s26 + s56)*(s16 + s26 - s34 + s56))



          amps(2) =
     &  (gb**2*gw**2*prop34*za(jb,jc)*
     -    (-(((gb**2 + 3*gw**2)*prop25*
     -           ((-2*s13 + s25 + s26 + 2*(s45 + s46) + s56)*za(jn,jd) - 
     -             2*za(jb,jd)*za(ju,jn)*zb(jb,ju) + 2*za(jb,jn)*za(je,jd)*zb(je,jb))*
     -           zb(je,ju))/(s25 + s26 + s56)) + 
     -      (2*gw**2*((-2*s13 + s25 + s26 + 2*(s45 + s46) + s56)*za(jn,jd) - 
     -           2*za(jb,jd)*za(ju,jn)*zb(jb,ju) + 2*za(jb,jn)*za(je,jd)*zb(je,jb))*
     -         zb(je,ju))/(s25*(s25 + s26 + s56)) - 
     -      6*gw**2*prop16*prop25*(2*za(jb,jd)*zb(jb,ju)*
     -          (za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb)) + 
     -         (2*s12 + s16 + s25 + 2*s26 - s34)*za(jn,jd)*zb(je,ju) + 
     -         2*za(jb,jn)*zb(je,jb)*(za(je,jd)*zb(je,ju) + za(jn,jd)*zb(jn,ju))) + 
     -      (6*gw**2*prop16*(2*za(jb,jd)*zb(jb,ju)*
     -            (za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb)) + 
     -           (2*s12 + s16 + s25 + 2*s26 - s34)*za(jn,jd)*zb(je,ju) + 
     -           2*za(jb,jn)*zb(je,jb)*(za(je,jd)*zb(je,ju) + za(jn,jd)*zb(jn,ju))))/s25
     -       - ((gb**2 - 3*gw**2)*prop25*
     -         (2*za(jb,jd)*zb(jb,ju)*
     -            (za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb) + za(ju,jn)*zb(je,ju)) + 
     -           za(jn,jd)*((2*s12 + 2*s13 + s16 + s26 - s34 - 2*(s45 + s46) - s56)*
     -               zb(je,ju) + 2*za(jb,jn)*zb(je,jb)*zb(jn,ju))))/
     -       (s16 + s26 - s34 + s56) - 
     -      (4*gw**2*(2*za(jb,jd)*zb(jb,ju)*
     -            (za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb) + za(ju,jn)*zb(je,ju)) + 
     -           za(jn,jd)*((2*s12 + 2*s13 + s16 + s26 - s34 - 2*(s45 + s46) - s56)*
     -               zb(je,ju) + 2*za(jb,jn)*zb(je,jb)*zb(jn,ju))))/
     -       (s25*(s16 + s26 - s34 + s56))))/(3._dp*ecossin**2)


      end subroutine 

      ! scale parameter sc just for rescaling 'included' masses
      subroutine singletop2_hamp_tree(za,zb, ju,jb,jn,je,jc,jd, mtsq, ampsSM, ampsAnomL2, ampsAnomNonresL2, sc)
        use types
        use anomcoup_tbW
        use eftcouple
        use constants
        implicit none

        include 'nf.f'
        include 'mxpart.f'

        complex(dp), intent(out) :: ampsSM(2) ! first: ----, second: ++-- (only non-resonant)
        complex(dp), intent(out) :: ampsAnomL2(5), ampsAnomNonresL2(5)
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        integer, intent(in) :: ju,jb,jn,je,jc,jd
        complex(dp), intent(in) :: mtsq
        real(dp), intent(in) :: sc

        complex(dp) :: ampsNonresSM(2)

        real(dp):: st,prop
        complex(dp) :: prop16, prop34, prop126, prop25

        integer :: j,k,l
        real(dp) :: s,t
        s(j,k) = real(za(j,k)*zb(k,j))
        t(j,k,l) = s(j,k) + s(j,l) + s(k,l)

        real(dp) :: s126, s12, s15, s16, s25, s56, s14, s36
        real(dp) :: s26, s45, s46, s24, s13, s23, s34, s35

        complex(dp) :: anomc1_sc, anomc2_sc, anomc3_sc, anomc4_sc
        complex(dp) :: anomc6_sc, anomc7_sc, anomc8_sc, anomc9_sc

        anomc1_sc = anomc1 * sc
        anomc2_sc = anomc2 * sc
        anomc3_sc = anomc3 * sqrt(sc)
        anomc4_sc = anomc4 * sqrt(sc)
        anomc6_sc = anomc6 * sqrt(sc)
        anomc7_sc = anomc7 * sqrt(sc)
        anomc8_sc = anomc8
        anomc9_sc = anomc9

        prop16 = 1._dp / (s(ju,jd)-wmass**2*sc)
        prop34 = 1._dp / (s(jn,je) - wmass**2*sc + im*wmass*wwidth*sc)
        prop126 = 1._dp / (t(ju,jb,jd) - mtsq)
        prop25 = 1._dp / (s(jb,jc) - zmass**2*sc)

        s12 = s(ju,jb)
        s13 = s(ju,jn)
        s14 = s(ju,je)
        s15 = s(ju,jc)
        s16 = s(ju,jd)
        s23 = s(jb,jn)
        s24 = s(jb,je)
        s25 = s(jb,jc)
        s26 = s(jb,jd)
        s34 = s(jn,je)
        s35 = s(jn,jc)
        s36 = s(jn,jd)
        s45 = s(je,jc)
        s46 = s(je,jd)
        s24 = s(jb,je)
        s56 = s(jc,jd)
        s126 = t(ju,jb,jd)

        ampsSM = 0._dp
        ampsAnomL2 = 0._dp

        ampsNonresSM = 0._dp
        ampsAnomNonresL2 = 0._dp

        if (enable_resonant) then
            ampsSM(1) = -(prop126*prop16*prop34*(-2*za(jc,jd)*za(jn,jc)*zb(jb,ju)*zb(jc,je) -
     -           2*za(jn,jc)*za(jn,jd)*zb(je,ju)*zb(jn,jb) + 2*za(jn,jc)*za(jn,jd)*zb(je,jb)*zb(jn,ju)))/2._dp

        if (enable_lambda2) then
            ampsAnomL2(1) = 2*prop126*prop34*za(jn,jc)*(2*im*Sqrt(mtsq)*prop16*imag(anomc3_sc)*zb(jb,ju)*
     -     (-(za(jc,jd)*zb(jc,je)) + za(jb,jd)*zb(je,jb)) + 
     -    prop16*real(anomc1_sc)*(za(jc,jd)*zb(jb,ju)*zb(jc,je) + 
     -       za(jn,jd)*(zb(je,ju)*zb(jn,jb) - zb(je,jb)*zb(jn,ju))) + 
     -    2*(real(anomc8_sc)*(za(jc,jd)*zb(jb,ju)*zb(jc,je) + 
     -          za(jn,jd)*(zb(je,ju)*zb(jn,jb) - zb(je,jb)*zb(jn,ju))) + 
     -       Sqrt(mtsq)*prop16*real(anomc3_sc)*
     -        (za(jc,jd)*zb(jb,ju)*zb(jc,je) - za(jb,jd)*zb(jb,ju)*zb(je,jb) + 
     -          2*za(jn,jd)*(zb(je,ju)*zb(jn,jb) - zb(je,jb)*zb(jn,ju)))))
        endif

        if (enable_lambda4) then
            ampsAnomL2(2) = prop16*prop34*(-(im*Sqrt(mtsq)*prop126*imag(anomc2_sc)*za(jn,jd)*zb(jb,ju)*zb(jc,je)) +
     -    Sqrt(mtsq)*prop126*real(anomc2_sc)*za(jn,jd)*zb(jb,ju)*zb(jc,je) +
     -    2*im*imag(anomc4_sc)*za(jn,jd)*(2*zb(jb,ju)*zb(jc,je) - zb(jc,jb)*zb(je,ju)) +
     -    2*real(anomc4_sc)*za(jn,jd)*(-2*zb(jb,ju)*zb(jc,je) + zb(jc,jb)*zb(je,ju)) +
     -    4*prop126*real(anomc4_sc)*(-(za(jb,jd)*zb(jb,ju)*
     -          (za(jb,jn)*zb(jc,jb)*zb(je,jb) + za(jn,jd)*zb(jd,jc)*zb(je,jb) +
     -            za(ju,jn)*(zb(jc,ju)*zb(je,jb) + zb(jc,jb)*zb(je,ju)))) +
     -       zb(je,ju)*(-(za(ju,jd)*za(ju,jn)*zb(jb,ju)*zb(jc,ju)) +
     -          za(jn,jd)*(za(ju,jb)*zb(jb,ju)*zb(jc,jb) +
     -             za(ju,jd)*(-(zb(jc,ju)*zb(jd,jb)) + zb(jc,jb)*zb(jd,ju)) -
     -             za(ju,jn)*(zb(jb,ju)*zb(jc,jn) + zb(jc,ju)*zb(jn,jb) - zb(jc,jb)*zb(jn,ju))))) +
     -    4*im*prop126*imag(anomc4_sc)*(za(jb,jd)*zb(jb,ju)*
     -        (za(jb,jn)*zb(jc,jb)*zb(je,jb) + za(jn,jd)*zb(jd,jc)*zb(je,jb) +
     -          za(ju,jn)*(zb(jc,ju)*zb(je,jb) + zb(jc,jb)*zb(je,ju))) +
     -       zb(je,ju)*(za(ju,jd)*za(ju,jn)*zb(jb,ju)*zb(jc,ju) +
     -          za(jn,jd)*(-(za(ju,jb)*zb(jb,ju)*zb(jc,jb)) +
     -             za(ju,jd)*(zb(jc,ju)*zb(jd,jb) - zb(jc,jb)*zb(jd,ju)) +
     -             za(ju,jn)*(zb(jb,ju)*zb(jc,jn) + zb(jc,ju)*zb(jn,jb) - zb(jc,jb)*zb(jn,ju))))))

            ampsAnomL2(3) = prop16*prop34*(-(Sqrt(mtsq)*prop126*(im*imag(anomc2_sc)
     &          + real(anomc2_sc))*za(jb,jd)*za(jn,jc)*zb(je,ju)) -
     -    2*im*imag(anomc4_sc)*(za(jb,jc)*(2*prop126*za(jb,jd)*zb(jb,ju)*
     -           (za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb)) + za(jn,jd)*zb(je,ju)) +
     -       2*za(jb,jd)*((za(jn,jc)*(-1 + prop126*za(ju,jb)*zb(jb,ju) + prop126*za(jb,jd)*zb(jd,jb)) +
     -             prop126*(za(jn,jd)*za(ju,jc) - za(jc,jd)*za(ju,jn))*zb(jd,ju))*zb(je,ju) +
     -          prop126*za(jb,jn)*zb(jb,ju)*(za(jc,jd)*zb(jd,je) + za(ju,jc)*zb(je,ju)))) -
     -    2*real(anomc4_sc)*(za(jb,jc)*(2*prop126*za(jb,jd)*zb(jb,ju)*
     -           (za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb)) + za(jn,jd)*zb(je,ju)) +
     -       2*za(jb,jd)*((za(jn,jc)*(-1 + prop126*za(ju,jb)*zb(jb,ju) + prop126*za(jb,jd)*zb(jd,jb)) +
     -             prop126*(za(jn,jd)*za(ju,jc) - za(jc,jd)*za(ju,jn))*zb(jd,ju))*zb(je,ju) +
     -          prop126*za(jb,jn)*zb(jb,ju)*(za(jc,jd)*zb(jd,je) + za(ju,jc)*zb(je,ju)))))

            ampsAnomL2(4) = -2*Sqrt(mtsq)*prop126*prop34*real(anomc9_sc)*za(jn,jc)*za(ju,jb)*zb(jd,je)


        endif ! enable_lambda4

            ampsSM(:) = ampsSM(:) * gw**4
            ampsAnomL2(:) = ampsAnomL2(:) * gw**4

        endif ! enable_resonant

        ! could simplify things in sum, but we want to keep things hackable (for now)
        if (enable_nonresonant) then
            ! WWG vertex contribution
            ampsNonresSM(1) = ampsNonresSM(1) + (2*gb**2*gw**4*prop16*prop34*za(jn,jc)*
     -    ((s12 - s126 + s14 + s16 - s23 - s25)*(s13 + s23 + s34)*za(jb,jd)*zb(jb,ju)*zb(je,jb) - 
     -      s12*(s12 - s126 + s14 + s16 - s23 - s25)*za(jn,jd)*zb(je,jb)*zb(jn,ju) + 
     -      (s126*s13 - s13*s16 + (s126 + s13 - s16)*s23 + s23**2 + (s13 + s23)*s25 + 
     -         (s126 - s16 + s23 + s25)*s34 + s12*s35)*
     -       (za(jc,jd)*zb(jb,ju)*zb(jc,je) + 
     -         za(jn,jd)*(zb(je,ju)*zb(jn,jb) - zb(je,jb)*zb(jn,ju)))))/
     -  (3._dp*ecossin**2*s12*s25*s35)

            ampsNonresSM(2) = ampsNonresSM(2) + 
     -              (2*gb**2*gw**4*prop16*prop34*
     -    (za(jn,jd)*za(ju,jb)*zb(jc,ju)*zb(je,ju) + 
     -      za(jb,jd)*(za(jn,jc)*zb(jc,je)*zb(jc,ju) + za(jb,jn)*zb(jc,ju)*zb(je,jb) + 
     -         za(jn,jd)*zb(jd,jc)*zb(je,ju)) - 
     -      za(jb,jn)*zb(jc,je)*(za(je,jd)*zb(je,ju) + za(jn,jd)*zb(jn,ju))))/
     -  (3._dp*ecossin**2*s25)

c           ! WWZ vertex contribution
            ampsNonresSM(1) = ampsNonresSM(1) +
     -  (gw**4*(gb**2 + 3*gw**2)*prop16*prop25*prop34*za(jn,jc)*
     -    ((s12 - s126 + s14 + s16 - s23 - s25)*(s13 + s23 + s34)*za(jb,jd)*zb(jb,ju)*zb(je,jb) - 
     -      s12*(s12 - s126 + s14 + s16 - s23 - s25)*za(jn,jd)*zb(je,jb)*zb(jn,ju) + 
     -      (s126*s13 - s13*s16 + (s126 + s13 - s16)*s23 + s23**2 + (s13 + s23)*s25 + 
     -         (s126 - s16 + s23 + s25)*s34 + s12*s35)*
     -       (za(jc,jd)*zb(jb,ju)*zb(jc,je) + 
     -         za(jn,jd)*(zb(je,ju)*zb(jn,jb) - zb(je,jb)*zb(jn,ju)))))/(3._dp*ecossin**2*s12*s35)

            ampsNonresSM(2) = ampsNonresSM(2) + 
     -  (-2*gb**2*gw**4*prop16*prop25*prop34*
     -    (za(jn,jd)*za(ju,jb)*zb(jc,ju)*zb(je,ju) + 
     -      za(jb,jd)*(za(jn,jc)*zb(jc,je)*zb(jc,ju) + za(jb,jn)*zb(jc,ju)*zb(je,jb) + 
     -         za(jn,jd)*zb(jd,jc)*zb(je,ju)) - 
     -      za(jb,jn)*zb(jc,je)*(za(je,jd)*zb(je,ju) + za(jn,jd)*zb(jn,ju))))/
     -  (3._dp*ecossin**2)

c           ! GR, G/photon exchange with W decay after t-channel exchange (Right side)
            ampsNonresSM(1) = ampsNonresSM(1) + 
     -  (4*gb**2*gw**4*prop34*za(jn,jc)*
     -    ((s13 + s23 + s34 + s35)*za(jc,jd)*zb(jb,ju)*zb(jc,je) + 
     -      (s12 - s126 + s13 + s14 + s16 - s25 + s34 + s35)*za(jn,jd)*
     -       (zb(je,ju)*zb(jn,jb) - zb(je,jb)*zb(jn,ju))))/
     -  (9._dp*ecossin**2*s25*(s12 + s15 + s25)*s35)

            ampsNonresSM(2) = ampsNonresSM(2) +
     -  (4*gb**2*gw**4*prop34*(za(jb,jd)*
     -       ((s126 + s14 - s16 + s23 + s25 + s35)*zb(jc,ju)*
     -          (za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb)) + 
     -         ((s126 + s13 + s14 - s16 + 2*s23 + s25 + s34 + 2*s35)*za(ju,jn)*
     -             zb(jc,ju) + (s12 + s13 + s14 + s16)*za(jn,jd)*zb(jd,jc))*zb(je,ju)) + 
     -      za(jn,jd)*((s12 - s126 + s16)*za(ju,jb)*zb(jc,ju)*zb(je,ju) + 
     -         (s126 - s16 + s23 + s25 + s35)*za(jb,jn)*
     -          (zb(jc,jn)*zb(je,ju) - zb(jc,je)*zb(jn,ju)))))/
     -  (9._dp*ecossin**2*s25*(s12 + s15 + s25)*(s126 - s16 + s23 + s25 + s35))

c           ! GL, G/photon exchange with W decay before t-channel exchange (Left side)
            ampsNonresSM(1) = ampsNonresSM(1) +
     -  (2*gb**2*gw**4*prop34*za(jn,jc)*
     -    (-((s12 - s126 + s16 - s23 - s25)*(s13 + s23 + s34)*za(jc,jd)*zb(jb,ju)*
     -         zb(jc,je)) + (s12 - s126 + s14 + s16 - s23 - s25)*(s13 + s23 + s34)*
     -       za(jb,jd)*zb(jb,ju)*zb(je,jb) + 
     -      za(jn,jd)*((-s12**2 + s12*(s126 - s13 - s14 - s16 + s25 - s34) + 
     -            (s126 - s16 + s23 + s25)*(s13 + s23 + s34))*zb(je,ju)*zb(jn,jb) + 
     -         (s12 - s126 + s16 - s23 - s25)*(s13 + s23 + s34)*zb(je,jb)*zb(jn,ju))))/
     -  (9._dp*ecossin**2*s12*s25*s35*(s25 + s26 + s56))

            ampsNonresSM(2) = ampsNonresSM(2) +
     -  (-2*gb**2*gw**4*prop34*(((s126 - s16 + s23 + s25 + s35)*za(jb,jn)*
     -          (za(je,jd)*zb(jc,je) + za(jn,jd)*zb(jc,jn)) - 
     -         (-s12 + 2*s126 - 2*s16 + s23 + s25 + s35)*za(jn,jd)*za(ju,jb)*zb(jc,ju))*
     -       zb(je,ju) + za(jb,jd)*(s14*zb(jc,ju)*
     -          (za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb)) + 
     -         ((s126 + s13 + s14 - s16 + 2*s23 + s25 + s34 + 2*s35)*za(ju,jn)*
     -             zb(jc,ju) + (s12 - s126 + s13 + s14 + 2*s16 - s23 - s25 - s35)*
     -             za(jn,jd)*zb(jd,jc))*zb(je,ju))))/
     -  (9._dp*ecossin**2*s25*(s126 - s16 + s23 + s25 + s35)*(s25 + s26 + s56))

c           ! ZR, Z exchange with W decay after t-channel exchange (Right side)
            ampsNonresSM(1) = ampsNonresSM(1) +
     -  (gw**2*(gb**2 - 3*gw**2)*(gb**2 + 3*gw**2)*prop25*prop34*za(jn,jc)*
     -    (-((s13 + s23 + s34 + s35)*za(jc,jd)*zb(jb,ju)*zb(jc,je)) + 
     -      (s12 - s126 + s13 + s14 + s16 - s25 + s34 + s35)*za(jn,jd)*
     -       (-(zb(je,ju)*zb(jn,jb)) + zb(je,jb)*zb(jn,ju))))/
     -  (18._dp*ecossin**2*(s12 + s15 + s25)*s35)

            ampsNonresSM(2) = ampsNonresSM(2) +
     -  (gb**2*gw**2*(gb**2 - 3*gw**2)*prop25*prop34*
     -    (za(jb,jd)*((s126 + s14 - s16 + s23 + s25 + s35)*zb(jc,ju)*
     -          (za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb)) + 
     -         ((s126 + s13 + s14 - s16 + 2*s23 + s25 + s34 + 2*s35)*za(ju,jn)*
     -             zb(jc,ju) + (s12 + s13 + s14 + s16)*za(jn,jd)*zb(jd,jc))*zb(je,ju)) + 
     -      za(jn,jd)*((s12 - s126 + s16)*za(ju,jb)*zb(jc,ju)*zb(je,ju) + 
     -         (s126 - s16 + s23 + s25 + s35)*za(jb,jn)*
     -          (zb(jc,jn)*zb(je,ju) - zb(jc,je)*zb(jn,ju)))))/
     -  (9._dp*ecossin**2*(s12 + s15 + s25)*(s126 - s16 + s23 + s25 + s35))

c           ! ZL, Z exchange with W decay before t-channel exchange (Left side)
            ampsNonresSM(1) = ampsNonresSM(1) + 
     -  (gw**2*(gb**2 + 3*gw**2)**2*prop25*prop34*za(jn,jc)*
     -    (-((s12 - s126 + s16 - s23 - s25)*(s13 + s23 + s34)*za(jc,jd)*zb(jb,ju)*
     -         zb(jc,je)) + (s12 - s126 + s14 + s16 - s23 - s25)*(s13 + s23 + s34)*
     -       za(jb,jd)*zb(jb,ju)*zb(je,jb) + 
     -      za(jn,jd)*((-s12**2 + s12*(s126 - s13 - s14 - s16 + s25 - s34) + 
     -            (s126 - s16 + s23 + s25)*(s13 + s23 + s34))*zb(je,ju)*zb(jn,jb) + 
     -         (s12 - s126 + s16 - s23 - s25)*(s13 + s23 + s34)*zb(je,jb)*zb(jn,ju))))/
     -  (18._dp*ecossin**2*s12*s35*(s25 + s26 + s56))

            ampsNonresSM(2) = ampsNonresSM(2) + 
     -  (gb**2*gw**2*(gb**2 + 3*gw**2)*prop25*prop34*
     -    (((s126 - s16 + s23 + s25 + s35)*za(jb,jn)*
     -          (za(je,jd)*zb(jc,je) + za(jn,jd)*zb(jc,jn)) - 
     -         (-s12 + 2*s126 - 2*s16 + s23 + s25 + s35)*za(jn,jd)*za(ju,jb)*zb(jc,ju))*
     -       zb(je,ju) + za(jb,jd)*(s14*zb(jc,ju)*
     -          (za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb)) + 
     -         ((s126 + s13 + s14 - s16 + 2*s23 + s25 + s34 + 2*s35)*za(ju,jn)*
     -             zb(jc,ju) + (s12 - s126 + s13 + s14 + 2*s16 - s23 - s25 - s35)*
     -             za(jn,jd)*zb(jd,jc))*zb(je,ju))))/
     -  (9._dp*ecossin**2*(s126 - s16 + s23 + s25 + s35)*(s25 + s26 + s56))

            ampsSM(:) = ampsSM(:) + ampsNonresSM(:)

            if (enable_lambda2) then
                ampsAnomL2(1) = ampsAnomL2(1) + real(anomc1_sc)*ampsNonresSM(1)
                ampsAnomL2(5) = ampsAnomL2(5) + real(anomc1_sc)*ampsNonresSM(2)
            endif


            ! these non resonant EFT contributions don't interfere with the SM
            ! helicity amplitudes
            if (enable_lambda4) then

                ! WWG PM
                ampsAnomNonresL2(2) = ampsAnomNonresL2(2) + 
     -  (-2*gb**2*gw**4*prop16*prop34*(im*imag(anomc4_sc) - real(anomc4_sc))*zb(jc,jb)*
     -    (2*za(jb,jd)*zb(jb,ju)*(za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb)) + 
     -      2*za(jb,jn)*za(je,jd)*zb(je,jb)*zb(je,ju) + 
     -      za(jn,jd)*((2*s126 - s16 + s25 - s34)*zb(je,ju) + 
     -         2*za(jb,jn)*zb(je,jb)*zb(jn,ju))))/(ecossin**2*s25)
                ! WWG MP
                ampsAnomNonresL2(3) = ampsAnomNonresL2(3) +
     -  (-2*gb**2*gw**4*prop16*prop34*(im*imag(anomc4_sc) + real(anomc4_sc))*za(jb,jc)*
     -    (2*za(jb,jd)*zb(jb,ju)*(za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb)) + 
     -      2*za(jb,jn)*za(je,jd)*zb(je,jb)*zb(je,ju) + 
     -      za(jn,jd)*((2*s126 - s16 + s25 - s34)*zb(je,ju) + 
     -         2*za(jb,jn)*zb(je,jb)*zb(jn,ju))))/(ecossin**2*s25)

                ! WWZ PM
                ampsAnomNonresL2(2) = ampsAnomNonresL2(2) +
     -  (-2*gw**6*prop16*prop25*prop34*(im*imag(anomc4_sc) - real(anomc4_sc))*zb(jc,jb)*
     -    (2*za(jb,jd)*zb(jb,ju)*(za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb)) + 
     -      2*za(jb,jn)*za(je,jd)*zb(je,jb)*zb(je,ju) + 
     -      za(jn,jd)*((2*s126 - s16 + s25 - s34)*zb(je,ju) + 
     -         2*za(jb,jn)*zb(je,jb)*zb(jn,ju))))/ecossin**2
                ! WWZ MP
                ampsAnomNonresL2(3) = ampsAnomNonresL2(3) +
     -  (-2*gw**6*prop16*prop25*prop34*(im*imag(anomc4_sc) + real(anomc4_sc))*za(jb,jc)*
     -    (2*za(jb,jd)*zb(jb,ju)*(za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb)) + 
     -      2*za(jb,jn)*za(je,jd)*zb(je,jb)*zb(je,ju) + 
     -      za(jn,jd)*((2*s126 - s16 + s25 - s34)*zb(je,ju) + 
     -         2*za(jb,jn)*zb(je,jb)*zb(jn,ju))))/ecossin**2

                ! GR PM
                ampsAnomNonresL2(2) = ampsAnomNonresL2(2) +
     -  (-4*gb**2*gw**4*prop34*(im*imag(anomc4_sc) - real(anomc4_sc))*zb(jc,jb)*
     -    (2*za(jb,jd)*zb(jb,ju)*(za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb) + 
     -         za(ju,jn)*zb(je,ju)) + 
     -      za(jn,jd)*((s13 + s14 + s16 - 2*s23 - s25)*zb(je,ju) + 
     -         2*za(jb,jn)*zb(je,jb)*zb(jn,ju))))/(3.*ecossin**2*s25*(s12 + s15 + s25))
                ! GR MP
                ampsAnomNonresL2(3) = ampsAnomNonresL2(3) +
     -  (-4*gb**2*gw**4*prop34*(im*imag(anomc4_sc) + real(anomc4_sc))*za(jb,jc)*
     -    (2*za(jb,jd)*zb(jb,ju)*(za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb) + 
     -         za(ju,jn)*zb(je,ju)) + 
     -      za(jn,jd)*((s13 + s14 + s16 - 2*s23 - s25)*zb(je,ju) + 
     -         2*za(jb,jn)*zb(je,jb)*zb(jn,ju))))/(3.*ecossin**2*s25*(s12 + s15 + s25))

                ! GL PM
                ampsAnomNonresL2(2) = ampsAnomNonresL2(2) +
     -  (-2*gb**2*gw**4*prop34*(im*imag(anomc4_sc) - real(anomc4_sc))*zb(jc,jb)*
     -    ((2*s126 - s13 - s14 - 2*s16 + 2*(s23 + s25) - s34)*za(jn,jd) - 
     -      2*za(jb,jd)*za(ju,jn)*zb(jb,ju) + 2*za(jb,jn)*za(je,jd)*zb(je,jb))*zb(je,ju))/
     -  (3.*ecossin**2*s25*(s25 + s26 + s56))
                ! GL MP
                ampsAnomNonresL2(3) = ampsAnomNonresL2(3) +
     -  (-2*gb**2*gw**4*prop34*(im*imag(anomc4_sc) + real(anomc4_sc))*za(jb,jc)*
     -    ((2*s126 - s13 - s14 - 2*s16 + 2*(s23 + s25) - s34)*za(jn,jd) - 
     -      2*za(jb,jd)*za(ju,jn)*zb(jb,ju) + 2*za(jb,jn)*za(je,jd)*zb(je,jb))*zb(je,ju))/
     -  (3.*ecossin**2*s25*(s25 + s26 + s56))

                ! ZR PM
                ampsAnomNonresL2(2) = ampsAnomNonresL2(2) +
     -  (gw**4*(gb**2 - 3*gw**2)*prop25*prop34*(im*imag(anomc4_sc) - real(anomc4_sc))*zb(jc,jb)*
     -    (2*za(jb,jd)*zb(jb,ju)*(za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb) + 
     -         za(ju,jn)*zb(je,ju)) + 
     -      za(jn,jd)*((s13 + s14 + s16 - 2*s23 - s25)*zb(je,ju) + 
     -         2*za(jb,jn)*zb(je,jb)*zb(jn,ju))))/(3.*ecossin**2*(s12 + s15 + s25))
                ! ZR MP
                ampsAnomNonresL2(3) = ampsAnomNonresL2(3) +
     -  (gw**4*(gb**2 - 3*gw**2)*prop25*prop34*(im*imag(anomc4_sc) + real(anomc4_sc))*za(jb,jc)*
     -    (2*za(jb,jd)*zb(jb,ju)*(za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb) + 
     -         za(ju,jn)*zb(je,ju)) + 
     -      za(jn,jd)*((s13 + s14 + s16 - 2*s23 - s25)*zb(je,ju) + 
     -         2*za(jb,jn)*zb(je,jb)*zb(jn,ju))))/(3.*ecossin**2*(s12 + s15 + s25))

                ! ZL PM
                ampsAnomNonresL2(2) = ampsAnomNonresL2(2)
     -  -(gw**4*(gb**2 + 3*gw**2)*prop25*prop34*(im*imag(anomc4_sc) - real(anomc4_sc))*
     -     zb(jc,jb)*((2*s126 - s13 - s14 - 2*s16 + 2*(s23 + s25) - s34)*za(jn,jd) - 
     -       2*za(jb,jd)*za(ju,jn)*zb(jb,ju) + 2*za(jb,jn)*za(je,jd)*zb(je,jb))*zb(je,ju))
     -   /(3.*ecossin**2*(s25 + s26 + s56))
                ! ZL MP
                ampsAnomNonresL2(3) = ampsAnomNonresL2(3)
     -  -(gw**4*(gb**2 + 3*gw**2)*prop25*prop34*(im*imag(anomc4_sc) + real(anomc4_sc))*
     -     za(jb,jc)*((2*s126 - s13 - s14 - 2*s16 + 2*(s23 + s25) - s34)*za(jn,jd) - 
     -       2*za(jb,jd)*za(ju,jn)*zb(jb,ju) + 2*za(jb,jn)*za(je,jd)*zb(je,jb))*zb(je,ju))
     -   /(3.*ecossin**2*(s25 + s26 + s56))

                ampsAnomL2(:) = ampsAnomL2(:) + ampsAnomNonresL2(:)

            endif


        endif ! enable_nonresonant

      end

      function singletop2_amp_tree(p, ju,jb,jn,je,jc,jd)
        use types
        use anomcoup_tbW
        use debugtools_m
        implicit none
        real(dp) :: singletop2_amp_tree

        include 'constants.f'
        include 'nf.f'
        include 'mxpart.f'
        include 'cplx.h'
        include 'zprods_com.f'
        include 'masses.f'
        real(dp), intent(in) :: p(mxpart,4)
        integer, intent(in) :: ju,jb,jn,je,jc,jd
        real(dp):: st,prop
        complex(dp) :: prop16, prop34, prop126

        complex(dp) :: ampsSM(2)
        complex(dp) :: ampsAnomL2(5)
        complex(dp) :: ampsAnomNonresL2(5)

        complex(dp) :: mtsq

        integer :: j,k,l
        real(dp) :: s,t
        s(j,k) = real(za(j,k)*zb(k,j))
        t(j,k,l) = s(j,k) + s(j,l) + s(k,l)

        mtsq = mt**2 - im*mt*twidth

        call singletop2_hamp_tree(za,zb, ju,jb,jn,je,jc,jd, mtsq, ampsSM, ampsAnomL2, ampsAnomNonresL2, 1d0)

        singletop2_amp_tree = 0._dp

        if (.not. disable_sm) then
            singletop2_amp_tree = sum(abs(ampsSM)**2)
        endif

        if (enable_lambda2) then
            singletop2_amp_tree = singletop2_amp_tree 
     &              + 2*real(ampsSM(1)*conjg(ampsAnomL2(1)))/lambda**2
     &              + 2*real(ampsSM(2)*conjg(ampsAnomL2(5)))/lambda**2
        endif

        if (enable_lambda4) then
            singletop2_amp_tree = singletop2_amp_tree + sum(abs(ampsAnomL2)**2)/lambda**4
        endif

      end

      subroutine singletop2_real_light(p,msq)
        use types
        implicit none
        include 'nf.f'
        include 'mxpart.f'

        real(dp), intent(in) :: p(mxpart,4)
        real(dp), intent(out) :: msq(-nf:nf, -nf:nf)

        call singletop2_real(p,msq,.true.,.false.)
      end subroutine

      subroutine singletop2_real_heavy(p,msq)
        use types
        implicit none
        include 'nf.f'
        include 'mxpart.f'

        real(dp), intent(in) :: p(mxpart,4)
        real(dp), intent(out) :: msq(-nf:nf, -nf:nf)

        call singletop2_real(p,msq,.false.,.true.)
      end subroutine

      ! largely taken from qqb_tbb_g
      subroutine singletop2_real(p,msq,light,heavy)
        use types
        use singletop2_scale_m
        implicit none

        include 'constants.f'
        include 'nf.f'
        include 'mxpart.f'
        include 'ewcouple.f'
        include 'qcdcouple.f'
        include 'ckm.f'
        include 'nwz.f'
        include 'zprods_com.f'

        real(dp), intent(in) :: p(mxpart,4)
        real(dp), intent(out) :: msq(-nf:nf, -nf:nf)
        logical, intent(in) :: light, heavy

        real(dp) :: ubtdg_l, ubtdg_h

        real(dp) :: fac, qqb,qbq,qg,gq,qbg,gqb
        real(dp) :: ub,bu,ubarb,bubar,bg,gb


        integer :: j,k

        call spinoru(7,p,za,zb)
        fac=2._dp*cf*gw**8*xn**2*4._dp*pi

        ub = 0._dp
        bu = 0._dp
        ubarb = 0._dp
        bubar = 0._dp
        bg = 0._dp
        gb = 0._dp

        qqb = 0._dp
        qbq = 0._dp
        qg = 0._dp
        gq = 0._dp
        qbg = 0._dp
        gqb = 0._dp

c---      Additional radiation terms below that correspond for example to
c---      diagrams of the form g+b-->W(s+c)+t(->W+b)

        if (nwz == +1) then
          if (light) then
            ub = ub + aveqq*fac*ubtdg_l_new(1,2,3,4,5,6,7,p) * as_light_beam1
            bu = bu + aveqq*fac*ubtdg_l_new(2,1,3,4,5,6,7,p) * as_light_beam2
            ubarb = ubarb + aveqq*fac*ubtdg_l_new(6,2,3,4,5,1,7,p) * as_light_beam1
            bubar = bubar + aveqq*fac*ubtdg_l_new(6,1,3,4,5,2,7,p) * as_light_beam2

            bg = bg + aveqg*fac*ubtdg_l_new(7,1,3,4,5,6,2,p) * as_light_beam2
            gb = gb + aveqg*fac*ubtdg_l_new(7,2,3,4,5,6,1,p) * as_light_beam1

            ! qg = qg + aveqg*fac*ubtdg_l_new(1,7,3,4,5,6,2,p)
            ! gq = gq + aveqg*fac*ubtdg_l_new(2,7,3,4,5,6,1,p)
            ! qbg = qbg + aveqg*fac*ubtdg_l_new(6,7,3,4,5,1,2,p)
            ! gqb = qgb + aveqg*fac*ubtdg_l_new(6,7,3,4,5,2,1,p)
          endif

          if (heavy) then
            ub = ub + aveqq*fac*ubtdg_h_full(1,2,3,4,5,6,7,p) * as_heavy_beam2
            bu = bu + aveqq*fac*ubtdg_h_full(2,1,3,4,5,6,7,p) * as_heavy_beam1
            ubarb = ubarb + aveqq*fac*ubtdg_h_full(6,2,3,4,5,1,7,p) * as_heavy_beam2
            bubar = bubar + aveqq*fac*ubtdg_h_full(6,1,3,4,5,2,7,p) * as_heavy_beam1

            ! bg = bg + aveqg*fac*ubtdg_h_full(7,1,3,4,5,6,2,p)
            ! gb = gb + +aveqg*fac*ubtdg_h_full(7,2,3,4,5,6,1,p)
            qg = qg + aveqg*fac*ubtdg_h_full(1,7,3,4,5,6,2,p) * as_heavy_beam2
            gq = gq + aveqg*fac*ubtdg_h_full(2,7,3,4,5,6,1,p) * as_heavy_beam1
            qbg = qbg + aveqg*fac*ubtdg_h_full(6,7,3,4,5,1,2,p) * as_heavy_beam2
            gqb = gqb + aveqg*fac*ubtdg_h_full(6,7,3,4,5,2,1,p) * as_heavy_beam1
          endif

        elseif (nwz == -1) then
          if (light) then
            ub = ub + aveqq*fac*ubtdg_l_new(6,2,4,3,5,1,7,p) * as_light_beam1
            bu = bu + aveqq*fac*ubtdg_l_new(6,1,4,3,5,2,7,p) * as_light_beam2
            ubarb = ubarb + aveqq*fac*ubtdg_l_new(1,2,4,3,5,6,7,p) * as_light_beam1
            bubar = bubar + aveqq*fac*ubtdg_l_new(2,1,4,3,5,6,7,p) * as_light_beam2

            bg = bg + aveqg*fac*(ubtdg_l_new(6,1,4,3,5,7,2,p)) * as_light_beam2
            gb = gb + aveqg*fac*(ubtdg_l_new(6,2,4,3,5,7,1,p)) * as_light_beam1
          endif

          if (heavy) then
            ub = ub + aveqq*fac*ubtdg_h_full(6,2,4,3,5,1,7,p) * as_heavy_beam2
            bu = bu + aveqq*fac*ubtdg_h_full(6,1,4,3,5,2,7,p) * as_heavy_beam1
            ubarb = ubarb + aveqq*fac*ubtdg_h_full(1,2,4,3,5,6,7,p) * as_heavy_beam2
            bubar = bubar + aveqq*fac*ubtdg_h_full(2,1,4,3,5,6,7,p) * as_heavy_beam1

            qg = qg + aveqg*fac*ubtdg_h_full(6,7,4,3,5,1,2,p) * as_heavy_beam2
            gq = gq + aveqg*fac*ubtdg_h_full(6,7,4,3,5,2,1,p) * as_heavy_beam1
            qbg = qbg + aveqg*fac*ubtdg_h_full(1,7,4,3,5,6,2,p) * as_heavy_beam2
            gqb = gqb + aveqg*fac*ubtdg_h_full(2,7,4,3,5,6,1,p) * as_heavy_beam1
          endif
        endif

      msq(:,:)=0._dp

      if (nwz == +1) then
      do j=-nf,nf
      do k=-nf,nf
c--- Q-Qbar
      if     ((j > 0) .and. (k < 0)) then
        if (j == 5) then
          msq(j,k)=(Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)
     &             +Vsq(+4,k)+Vsq(+5,k))*bubar
        else
          msq(j,k)=Vsq(j,k)*qqb
        endif
c--- Qbar-Q
      elseif ((j < 0) .and. (k > 0)) then
        if (k == 5) then
          msq(j,k)=(Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)
     &             +Vsq(j,+4)+Vsq(j,+5))*ubarb
        else
          msq(j,k)=Vsq(j,k)*qbq
        endif
c--- Q-Q
      elseif ((j == 5) .and. (k > 0)) then
        msq(j,k)=(Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)
     &           +Vsq(-4,k)+Vsq(-5,k))*bu
      elseif ((j > 0) .and. (k == 5)) then
        msq(j,k)=(Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)
     &           +Vsq(j,-4)+Vsq(j,-5))*ub
c--- g-Q
      elseif ((j == 0) .and. (k > 0)) then
        if (k == 5) then
          msq(j,k)=2._dp*gb
        else
          msq(j,k)=(Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)
     &             +Vsq(-4,k)+Vsq(-5,k))*gq
        endif
c--- g-Qbar
      elseif ((j == 0) .and. (k < 0)) then
        msq(j,k)=(Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)
     &           +Vsq(+4,k)+Vsq(+5,k))*gqb
c--- Q-g
      elseif ((j > 0) .and. (k == 0)) then
        if (j == 5) then
          msq(j,k)=2._dp*bg
        else
          msq(j,k)=(Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)
     &             +Vsq(j,-4)+Vsq(j,-5))*qg
        endif
c--- Qbar-g
      elseif ((j < 0) .and. (k == 0)) then
        msq(j,k)=(Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)
     &           +Vsq(j,+4)+Vsq(j,+5))*qbg
      endif

      enddo
      enddo
      elseif (nwz == -1) then

      do j=-nf,nf
      do k=-nf,nf

c--- Q-Qbar
      if     ((j > 0) .and. (k < 0)) then
        if (k == -5) then
          msq(j,k)=(Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)
     &             +Vsq(j,-4)+Vsq(j,-5))*ub
        else
          msq(j,k)=Vsq(j,k)*qqb
        endif
c--- Qbar-Q
      elseif ((j < 0) .and. (k > 0)) then
        if (j == -5) then
          msq(j,k)=(Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)
     &             +Vsq(-4,k)+Vsq(-5,k))*bu
        else
          msq(j,k)=Vsq(j,k)*qbq
        endif
c--- Qbar-Qbar
      elseif ((j == -5) .and. (k < 0)) then
        msq(j,k)=(Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)
     &           +Vsq(+4,k)+Vsq(+5,k))*bubar
      elseif ((j < 0) .and. (k == -5)) then
        msq(j,k)=(Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)
     &           +Vsq(j,+4)+Vsq(j,+5))*ubarb
c--- g-Q
      elseif ((j == 0) .and. (k > 0)) then
          msq(j,k)=(Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)
     &             +Vsq(-4,k)+Vsq(-5,k))*gq
c--- g-Qbar
      elseif ((j == 0) .and. (k < 0)) then
        if (k == -5) then
          msq(j,k)=2._dp*gb
        else
        msq(j,k)=(Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)
     &           +Vsq(+4,k)+Vsq(+5,k))*gqb
        endif
c--- Q-g
      elseif ((j > 0) .and. (k == 0)) then
          msq(j,k)=(Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)
     &             +Vsq(j,-4)+Vsq(j,-5))*qg
c--- Qbar-g
      elseif ((j < 0) .and. (k == 0)) then
        if (j == -5) then
          msq(j,k)=2._dp*bg
        else
        msq(j,k)=(Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)
     &           +Vsq(j,+4)+Vsq(j,+5))*qbg
        endif
      endif

      enddo
      enddo

      endif

      end subroutine

      function ubtdg_h_full(ju,jb,jn,je,jc,jd,jg,p)
        use types
        use anomcoup_tbW
        use singletop2_realamps_nonres_m
        implicit none

        include 'nf.f'
        include 'mxpart.f'
        include 'zprods_com.f'

        real(dp) :: ubtdg_h_full

        integer, intent(in) :: ju,jb,jn,je,jc,jd,jg
        real(dp), intent(in) :: p(mxpart,4)

        ! helicity amplitude labeling:
        ! first index gluon (1 = Plus, 2 = Minus)
        ! second index: (5,2,6,1) helicities in (1 = ----, 2=+---, 3=-+--, 4=-+++, 5=++--), so that 1 = resonant SM

        complex(dp) :: amps(2,5)
        complex(dp) :: ampsNonres(2,5)
        complex(dp) :: ampsAnomL2(2,5)

        amps(:,:) = 0._dp
        ampsNonres(:,:) = 0._dp
        ampsAnomL2(:,:) = 0._dp

        if (enable_resonant) then
            amps(2,1) = amps(2,1) + streal_heavyResonant_MMMM_M_SM(ju,jb,jn,je,jc,jd,jg,za,zb) 
            amps(1,1) = amps(1,1) + streal_heavyResonant_MMMM_P_SM(ju,jb,jn,je,jc,jd,jg,za,zb)
        endif

        if (enable_nonresonant) then
            ampsNonres(1,1) = ampsNonres(1,1) + streal_heavyWWG_MMMM_P_SM(ju,jb,jn,je,jc,jd,jg,za,zb)
            ampsNonres(2,1) = ampsNonres(2,1) + streal_heavyWWG_MMMM_M_SM(ju,jb,jn,je,jc,jd,jg,za,zb)

            ampsNonres(1,1) = ampsNonres(1,1) + streal_heavyWWZ_MMMM_P_SM(ju,jb,jn,je,jc,jd,jg,za,zb)
            ampsNonres(2,1) = ampsNonres(2,1) + streal_heavyWWZ_MMMM_M_SM(ju,jb,jn,je,jc,jd,jg,za,zb)

            ampsNonres(1,1) = ampsNonres(1,1) + streal_heavyGR_MMMM_P_SM(ju,jb,jn,je,jc,jd,jg,za,zb)
            ampsNonres(2,1) = ampsNonres(2,1) + streal_heavyGR_MMMM_M_SM(ju,jb,jn,je,jc,jd,jg,za,zb)

            ampsNonres(1,1) = ampsNonres(1,1) + streal_heavyGL_MMMM_P_SM(ju,jb,jn,je,jc,jd,jg,za,zb)
            ampsNonres(2,1) = ampsNonres(2,1) + streal_heavyGL_MMMM_M_SM(ju,jb,jn,je,jc,jd,jg,za,zb)

            ampsNonres(1,1) = ampsNonres(1,1) + streal_heavyZR_MMMM_P_SM(ju,jb,jn,je,jc,jd,jg,za,zb)
            ampsNonres(2,1) = ampsNonres(2,1) + streal_heavyZR_MMMM_M_SM(ju,jb,jn,je,jc,jd,jg,za,zb)

            ampsNonres(1,1) = ampsNonres(1,1) + streal_heavyZL_MMMM_P_SM(ju,jb,jn,je,jc,jd,jg,za,zb)
            ampsNonres(2,1) = ampsNonres(2,1) + streal_heavyZL_MMMM_M_SM(ju,jb,jn,je,jc,jd,jg,za,zb)



            ampsNonres(1,5) = ampsNonres(1,5) + streal_heavyWWG_PPMM_P_SM(ju,jb,jn,je,jc,jd,jg,za,zb)
            ampsNonres(2,5) = ampsNonres(2,5) + streal_heavyWWG_PPMM_M_SM(ju,jb,jn,je,jc,jd,jg,za,zb)

            ampsNonres(1,5) = ampsNonres(1,5) + streal_heavyWWZ_PPMM_P_SM(ju,jb,jn,je,jc,jd,jg,za,zb)
            ampsNonres(2,5) = ampsNonres(2,5) + streal_heavyWWZ_PPMM_M_SM(ju,jb,jn,je,jc,jd,jg,za,zb)

            ampsNonres(1,5) = ampsNonres(1,5) + streal_heavyGR_PPMM_P_SM(ju,jb,jn,je,jc,jd,jg,za,zb)
            ampsNonres(2,5) = ampsNonres(2,5) + streal_heavyGR_PPMM_M_SM(ju,jb,jn,je,jc,jd,jg,za,zb)

            ampsNonres(1,5) = ampsNonres(1,5) + streal_heavyGL_PPMM_P_SM(ju,jb,jn,je,jc,jd,jg,za,zb)
            ampsNonres(2,5) = ampsNonres(2,5) + streal_heavyGL_PPMM_M_SM(ju,jb,jn,je,jc,jd,jg,za,zb)

            ampsNonres(1,5) = ampsNonres(1,5) + streal_heavyZR_PPMM_P_SM(ju,jb,jn,je,jc,jd,jg,za,zb)
            ampsNonres(2,5) = ampsNonres(2,5) + streal_heavyZR_PPMM_M_SM(ju,jb,jn,je,jc,jd,jg,za,zb)

            ampsNonres(1,5) = ampsNonres(1,5) + streal_heavyZL_PPMM_P_SM(ju,jb,jn,je,jc,jd,jg,za,zb)
            ampsNonres(2,5) = ampsNonres(2,5) + streal_heavyZL_PPMM_M_SM(ju,jb,jn,je,jc,jd,jg,za,zb)

            amps(:,:) = amps(:,:) + ampsNonres(:,:)
        endif

        ubtdg_h_full = 0._dp

        if (.not. disable_sm) then
            ubtdg_h_full = sum(abs(amps)**2)
        endif

        if (enable_lambda2) then
            if (enable_resonant) then
                ampsAnomL2(2,1) = ampsAnomL2(2,1) + streal_heavyResonant_MMMM_M_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
                ampsAnomL2(1,1) = ampsAnomL2(1,1) + streal_heavyResonant_MMMM_P_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
            endif

            if (enable_nonresonant) then
                ampsAnomL2(:,:) = ampsAnomL2(:,:) + ampsNonres(:,:)*real(anomc1)
            endif

            ubtdg_h_full = ubtdg_h_full + sum(2*real(amps*conjg(ampsAnomL2)))/lambda**2
        endif

        if (enable_lambda4) then
            if (enable_resonant) then
                ampsAnomL2(2,2) = ampsAnomL2(2,2) + streal_heavyResonant_PMMM_M_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
                ampsAnomL2(1,2) = ampsAnomL2(1,2) + streal_heavyResonant_PMMM_P_L2(ju,jb,jn,je,jc,jd,jg,za,zb)

                ampsAnomL2(2,3) = ampsAnomL2(2,3) + streal_heavyResonant_MPMM_M_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
                ampsAnomL2(1,3) = ampsAnomL2(1,3) + streal_heavyResonant_MPMM_P_L2(ju,jb,jn,je,jc,jd,jg,za,zb)

                ! only from 4R
                ampsAnomL2(2,4) = ampsAnomL2(2,4) + streal_heavyResonant_MPPP_M_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
                ampsAnomL2(1,4) = ampsAnomL2(1,4) + streal_heavyResonant_MPPP_P_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
            endif

            if (enable_nonresonant) then
                ampsAnomL2(1,2) = ampsAnomL2(1,2) + streal_heavyWWG_PMMM_P_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
                ampsAnomL2(2,2) = ampsAnomL2(2,2) + streal_heavyWWG_PMMM_M_L2(ju,jb,jn,je,jc,jd,jg,za,zb)

                ampsAnomL2(1,2) = ampsAnomL2(1,2) + streal_heavyWWZ_PMMM_P_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
                ampsAnomL2(2,2) = ampsAnomL2(2,2) + streal_heavyWWZ_PMMM_M_L2(ju,jb,jn,je,jc,jd,jg,za,zb)

                ampsAnomL2(1,2) = ampsAnomL2(1,2) + streal_heavyGR_PMMM_P_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
                ampsAnomL2(2,2) = ampsAnomL2(2,2) + streal_heavyGR_PMMM_M_L2(ju,jb,jn,je,jc,jd,jg,za,zb)

                ampsAnomL2(1,2) = ampsAnomL2(1,2) + streal_heavyGL_PMMM_P_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
                ampsAnomL2(2,2) = ampsAnomL2(2,2) + streal_heavyGL_PMMM_M_L2(ju,jb,jn,je,jc,jd,jg,za,zb)

                ampsAnomL2(1,2) = ampsAnomL2(1,2) + streal_heavyZR_PMMM_P_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
                ampsAnomL2(2,2) = ampsAnomL2(2,2) + streal_heavyZR_PMMM_M_L2(ju,jb,jn,je,jc,jd,jg,za,zb)

                ampsAnomL2(1,2) = ampsAnomL2(1,2) + streal_heavyZL_PMMM_P_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
                ampsAnomL2(2,2) = ampsAnomL2(2,2) + streal_heavyZL_PMMM_M_L2(ju,jb,jn,je,jc,jd,jg,za,zb)



                ampsAnomL2(1,3) = ampsAnomL2(1,3) + streal_heavyWWG_MPMM_P_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
                ampsAnomL2(2,3) = ampsAnomL2(2,3) + streal_heavyWWG_MPMM_M_L2(ju,jb,jn,je,jc,jd,jg,za,zb)

                ampsAnomL2(1,3) = ampsAnomL2(1,3) + streal_heavyWWZ_MPMM_P_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
                ampsAnomL2(2,3) = ampsAnomL2(2,3) + streal_heavyWWZ_MPMM_M_L2(ju,jb,jn,je,jc,jd,jg,za,zb)

                ampsAnomL2(1,3) = ampsAnomL2(1,3) + streal_heavyGR_MPMM_P_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
                ampsAnomL2(2,3) = ampsAnomL2(2,3) + streal_heavyGR_MPMM_M_L2(ju,jb,jn,je,jc,jd,jg,za,zb)

                ampsAnomL2(1,3) = ampsAnomL2(1,3) + streal_heavyGL_MPMM_P_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
                ampsAnomL2(2,3) = ampsAnomL2(2,3) + streal_heavyGL_MPMM_M_L2(ju,jb,jn,je,jc,jd,jg,za,zb)

                ampsAnomL2(1,3) = ampsAnomL2(1,3) + streal_heavyZR_MPMM_P_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
                ampsAnomL2(2,3) = ampsAnomL2(2,3) + streal_heavyZR_MPMM_M_L2(ju,jb,jn,je,jc,jd,jg,za,zb)

                ampsAnomL2(1,3) = ampsAnomL2(1,3) + streal_heavyZL_MPMM_P_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
                ampsAnomL2(2,3) = ampsAnomL2(2,3) + streal_heavyZL_MPMM_M_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
            endif

            ubtdg_h_full = ubtdg_h_full + sum(abs(ampsAnomL2)**2)/lambda**4
        endif

      end function

      function ubtdg_l_new(ju,jb,jn,je,jc,jd,jg,p)
        use types
        use anomcoup_tbW
        use singletop2_realamps_nonres_m
        implicit none
        include 'nf.f'
        include 'mxpart.f'
        include 'zprods_com.f'

        real(dp):: ubtdg_l_new

        integer, intent(in) :: ju,jb,jn,je,jc,jd,jg
        real(dp), intent(in) :: p(mxpart,4)

        complex(dp) :: amps(2,5)
        complex(dp) :: ampsNonres(2,5)
        complex(dp) :: ampsAnomL2(2,5)

        amps(:,:) = 0._dp
        ampsNonres(:,:) = 0._dp
        ampsAnomL2(:,:) = 0._dp

        if (enable_resonant) then
            amps(1,1) = amps(1,1) + streal_lightResonant_MMMM_P_SM(ju,jb,jn,je,jc,jd,jg,za,zb)
            amps(2,1) = amps(2,1) + streal_lightResonant_MMMM_M_SM(ju,jb,jn,je,jc,jd,jg,za,zb)
        endif

        if (enable_nonresonant) then
            ampsNonres(1,1) = ampsNonres(1,1) + streal_lightWWG_MMMM_P_SM(ju,jb,jn,je,jc,jd,jg,za,zb)
            ampsNonres(2,1) = ampsNonres(2,1) + streal_lightWWG_MMMM_M_SM(ju,jb,jn,je,jc,jd,jg,za,zb)

            ampsNonres(1,1) = ampsNonres(1,1) + streal_lightWWZ_MMMM_P_SM(ju,jb,jn,je,jc,jd,jg,za,zb)
            ampsNonres(2,1) = ampsNonres(2,1) + streal_lightWWZ_MMMM_M_SM(ju,jb,jn,je,jc,jd,jg,za,zb)

            ampsNonres(1,1) = ampsNonres(1,1) + streal_lightGR_MMMM_P_SM(ju,jb,jn,je,jc,jd,jg,za,zb)
            ampsNonres(2,1) = ampsNonres(2,1) + streal_lightGR_MMMM_M_SM(ju,jb,jn,je,jc,jd,jg,za,zb)

            ampsNonres(1,1) = ampsNonres(1,1) + streal_lightGL_MMMM_P_SM(ju,jb,jn,je,jc,jd,jg,za,zb)
            ampsNonres(2,1) = ampsNonres(2,1) + streal_lightGL_MMMM_M_SM(ju,jb,jn,je,jc,jd,jg,za,zb)

            ampsNonres(1,1) = ampsNonres(1,1) + streal_lightZR_MMMM_P_SM(ju,jb,jn,je,jc,jd,jg,za,zb)
            ampsNonres(2,1) = ampsNonres(2,1) + streal_lightZR_MMMM_M_SM(ju,jb,jn,je,jc,jd,jg,za,zb)

            ampsNonres(1,1) = ampsNonres(1,1) + streal_lightZL_MMMM_P_SM(ju,jb,jn,je,jc,jd,jg,za,zb)
            ampsNonres(2,1) = ampsNonres(2,1) + streal_lightZL_MMMM_M_SM(ju,jb,jn,je,jc,jd,jg,za,zb)



            ampsNonres(1,5) = ampsNonres(1,5) + streal_lightWWG_PPMM_P_SM(ju,jb,jn,je,jc,jd,jg,za,zb)
            ampsNonres(2,5) = ampsNonres(2,5) + streal_lightWWG_PPMM_M_SM(ju,jb,jn,je,jc,jd,jg,za,zb)

            ampsNonres(1,5) = ampsNonres(1,5) + streal_lightWWZ_PPMM_P_SM(ju,jb,jn,je,jc,jd,jg,za,zb)
            ampsNonres(2,5) = ampsNonres(2,5) + streal_lightWWZ_PPMM_M_SM(ju,jb,jn,je,jc,jd,jg,za,zb)

            ampsNonres(1,5) = ampsNonres(1,5) + streal_lightGR_PPMM_P_SM(ju,jb,jn,je,jc,jd,jg,za,zb)
            ampsNonres(2,5) = ampsNonres(2,5) + streal_lightGR_PPMM_M_SM(ju,jb,jn,je,jc,jd,jg,za,zb)

            ampsNonres(1,5) = ampsNonres(1,5) + streal_lightGL_PPMM_P_SM(ju,jb,jn,je,jc,jd,jg,za,zb)
            ampsNonres(2,5) = ampsNonres(2,5) + streal_lightGL_PPMM_M_SM(ju,jb,jn,je,jc,jd,jg,za,zb)

            ampsNonres(1,5) = ampsNonres(1,5) + streal_lightZR_PPMM_P_SM(ju,jb,jn,je,jc,jd,jg,za,zb)
            ampsNonres(2,5) = ampsNonres(2,5) + streal_lightZR_PPMM_M_SM(ju,jb,jn,je,jc,jd,jg,za,zb)

            ampsNonres(1,5) = ampsNonres(1,5) + streal_lightZL_PPMM_P_SM(ju,jb,jn,je,jc,jd,jg,za,zb)
            ampsNonres(2,5) = ampsNonres(2,5) + streal_lightZL_PPMM_M_SM(ju,jb,jn,je,jc,jd,jg,za,zb)

            amps(:,:) = amps(:,:) + ampsNonres(:,:)
        endif

        ubtdg_l_new = 0._dp

        if (.not. disable_sm) then
            ubtdg_l_new = sum(abs(amps)**2)
        endif

        if (enable_lambda2) then
            if (enable_resonant) then
                ampsAnomL2(1,1) = streal_lightResonant_MMMM_P_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
                ampsAnomL2(2,1) = streal_lightResonant_MMMM_M_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
            endif

            if (enable_nonresonant) then
                ampsAnomL2(:,:) = ampsAnomL2(:,:) + ampsNonres(:,:)*real(anomc1)
            endif

            ubtdg_l_new = ubtdg_l_new + sum(2*real(amps*conjg(ampsAnomL2)))/lambda**2
        endif

        if (enable_lambda4) then
            if (enable_resonant) then
                ! enters as square
                ampsAnomL2(1,2) = streal_lightResonant_PMMM_P_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
                ampsAnomL2(2,2) = streal_lightResonant_PMMM_M_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
                ampsAnomL2(1,3) = streal_lightResonant_MPMM_P_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
                ampsAnomL2(2,3) = streal_lightResonant_MPMM_M_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
                ampsAnomL2(1,4) = streal_lightResonant_MPPP_P_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
                ampsAnomL2(2,4) = streal_lightResonant_MPPP_M_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
            endif

            if (enable_nonresonant) then
                ampsAnomL2(1,2) = ampsAnomL2(1,2) + streal_lightWWG_PMMM_P_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
                ampsAnomL2(2,2) = ampsAnomL2(2,2) + streal_lightWWG_PMMM_M_L2(ju,jb,jn,je,jc,jd,jg,za,zb)

                ampsAnomL2(1,2) = ampsAnomL2(1,2) + streal_lightWWZ_PMMM_P_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
                ampsAnomL2(2,2) = ampsAnomL2(2,2) + streal_lightWWZ_PMMM_M_L2(ju,jb,jn,je,jc,jd,jg,za,zb)

                ampsAnomL2(1,2) = ampsAnomL2(1,2) + streal_lightGR_PMMM_P_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
                ampsAnomL2(2,2) = ampsAnomL2(2,2) + streal_lightGR_PMMM_M_L2(ju,jb,jn,je,jc,jd,jg,za,zb)

                ampsAnomL2(1,2) = ampsAnomL2(1,2) + streal_lightGL_PMMM_P_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
                ampsAnomL2(2,2) = ampsAnomL2(2,2) + streal_lightGL_PMMM_M_L2(ju,jb,jn,je,jc,jd,jg,za,zb)

                ampsAnomL2(1,2) = ampsAnomL2(1,2) + streal_lightZR_PMMM_P_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
                ampsAnomL2(2,2) = ampsAnomL2(2,2) + streal_lightZR_PMMM_M_L2(ju,jb,jn,je,jc,jd,jg,za,zb)

                ampsAnomL2(1,2) = ampsAnomL2(1,2) + streal_lightZL_PMMM_P_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
                ampsAnomL2(2,2) = ampsAnomL2(2,2) + streal_lightZL_PMMM_M_L2(ju,jb,jn,je,jc,jd,jg,za,zb)



                ampsAnomL2(1,3) = ampsAnomL2(1,3) + streal_lightWWG_MPMM_P_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
                ampsAnomL2(2,3) = ampsAnomL2(2,3) + streal_lightWWG_MPMM_M_L2(ju,jb,jn,je,jc,jd,jg,za,zb)

                ampsAnomL2(1,3) = ampsAnomL2(1,3) + streal_lightWWZ_MPMM_P_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
                ampsAnomL2(2,3) = ampsAnomL2(2,3) + streal_lightWWZ_MPMM_M_L2(ju,jb,jn,je,jc,jd,jg,za,zb)

                ampsAnomL2(1,3) = ampsAnomL2(1,3) + streal_lightGR_MPMM_P_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
                ampsAnomL2(2,3) = ampsAnomL2(2,3) + streal_lightGR_MPMM_M_L2(ju,jb,jn,je,jc,jd,jg,za,zb)

                ampsAnomL2(1,3) = ampsAnomL2(1,3) + streal_lightGL_MPMM_P_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
                ampsAnomL2(2,3) = ampsAnomL2(2,3) + streal_lightGL_MPMM_M_L2(ju,jb,jn,je,jc,jd,jg,za,zb)

                ampsAnomL2(1,3) = ampsAnomL2(1,3) + streal_lightZR_MPMM_P_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
                ampsAnomL2(2,3) = ampsAnomL2(2,3) + streal_lightZR_MPMM_M_L2(ju,jb,jn,je,jc,jd,jg,za,zb)

                ampsAnomL2(1,3) = ampsAnomL2(1,3) + streal_lightZL_MPMM_P_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
                ampsAnomL2(2,3) = ampsAnomL2(2,3) + streal_lightZL_MPMM_M_L2(ju,jb,jn,je,jc,jd,jg,za,zb)
            endif

            ubtdg_l_new = ubtdg_l_new + sum(abs(ampsAnomL2)**2)/lambda**4
        endif

      end

      function singletop2_amp_protos(p, ju,jb,jn,je,jc,jd)
        use types
        use anomcoup_tbW
        use eftcouple
        implicit none
        real(dp) :: singletop2_amp_protos

        include 'constants.f'
        include 'nf.f'
        include 'mxpart.f'
        include 'zprods_com.f'
        real(dp), intent(in) :: p(mxpart,4)
        integer, intent(in) :: ju,jb,jn,je,jc,jd
        real(dp):: st,prop
        complex(dp) :: prop16, prop34, prop126

        complex(dp) :: ampSM
        complex(dp) :: ampsAnomL2(5), ampsAnomL4(5)

        real(dp) :: s12, s15, s16, s25, s56, s14, s36, s126
        real(dp) :: s26, s45, s46, s24, s13, s23, s34, s35

        complex(dp) :: mtsq

        real(dp) :: tagm, tagp

        integer :: j,k,l
        real(dp) :: s,t
        s(j,k) = real(za(j,k)*zb(k,j))
        t(j,k,l) = s(j,k) + s(j,l) + s(k,l)

        mtsq = mt**2 - im*mt*twidth

        prop16 = 1._dp / (s(ju,jd)-wmass**2)
        prop34 = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
        prop126 = 1._dp / (t(ju,jb,jd) - mtsq)

        s12 = s(ju,jb)
        s13 = s(ju,jn)
        s14 = s(ju,je)
        s15 = s(ju,jc)
        s16 = s(ju,jd)
        s23 = s(jb,jn)
        s24 = s(jb,je)
        s25 = s(jb,jc)
        s26 = s(jb,jd)
        s34 = s(jn,je)
        s35 = s(jn,jc)
        s36 = s(jn,jd)
        s45 = s(je,jc)
        s46 = s(je,jd)
        s24 = s(jb,je)
        s56 = s(jc,jd)
        s126 = t(ju,jb,jd)

        ampsAnomL2(:) = 0._dp
        ampsAnomL4(:) = 0._dp

        tagm = 1d0
        tagp = 1d0

        ampSM = gw**4*prop126*prop16*prop34*za(jn,jc)*
     -  (za(jc,jd)*zb(jb,ju)*zb(jc,je) + 
     -    za(jn,jd)*(zb(je,ju)*zb(jn,jb) - zb(je,jb)*zb(jn,ju)))

        ! --
        ampsAnomL2(1) = 2*gw**4*prop126*prop34*tagm*za(jn,jc)*
     -  (2*im*Sqrt(mt**2)*prop16*tagp*imag(anomc3)*zb(jb,ju)*
     -     (-(za(jc,jd)*zb(jc,je)) + za(jb,jd)*zb(je,jb)) + 
     -    prop16*tagp*real(anomc1)*(za(jc,jd)*zb(jb,ju)*zb(jc,je) + 
     -       za(jn,jd)*(zb(je,ju)*zb(jn,jb) - zb(je,jb)*zb(jn,ju))) + 
     -    2*(real(anomc8)*(za(jc,jd)*zb(jb,ju)*zb(jc,je) + 
     -          za(jn,jd)*(zb(je,ju)*zb(jn,jb) - zb(je,jb)*zb(jn,ju))) + 
     -       Sqrt(mt**2)*prop16*tagp*real(anomc3)*
     -        (za(jc,jd)*zb(jb,ju)*zb(jc,je) - za(jb,jd)*zb(jb,ju)*zb(je,jb) + 
     -          2*za(jn,jd)*(zb(je,ju)*zb(jn,jb) - zb(je,jb)*zb(jn,ju)))))

        ! --
        ampsAnomL4(1) = gw**4*prop126*prop34*tagm*za(jn,jc)*
     -  (-4*im*Sqrt(mt**2)*prop16*tagp*imag(anomc3)*real(anomc1)*zb(jb,ju)*
     -     (za(jc,jd)*zb(jc,je) - za(jb,jd)*zb(je,jb)) + 
     -    16*im*Sqrt(mt**2)*imag(anomc3)*real(anomc8)*za(jn,jd)*
     -     (zb(je,ju)*zb(jn,jb) - zb(je,jb)*zb(jn,ju)) + 
     -    16*Sqrt(mt**2)*real(anomc3)*real(anomc8)*za(jn,jd)*
     -     (zb(je,ju)*zb(jn,jb) - zb(je,jb)*zb(jn,ju)) + 
     -    prop16*tagp*real(anomc1)**2*
     -     (za(jc,jd)*zb(jb,ju)*zb(jc,je) + 
     -       za(jn,jd)*(zb(je,ju)*zb(jn,jb) - zb(je,jb)*zb(jn,ju))) + 
     -    4*real(anomc1)*real(anomc8)*
     -     (za(jc,jd)*zb(jb,ju)*zb(jc,je) + 
     -       za(jn,jd)*(zb(je,ju)*zb(jn,jb) - zb(je,jb)*zb(jn,ju))) + 
     -    4*Sqrt(mt**2)*prop16*tagp*real(anomc1)*real(anomc3)*
     -     (za(jc,jd)*zb(jb,ju)*zb(jc,je) - za(jb,jd)*zb(jb,ju)*zb(je,jb) + 
     -       2*za(jn,jd)*(zb(je,ju)*zb(jn,jb) - zb(je,jb)*zb(jn,ju))) + 
     -    16*prop16*tagp*imag(anomc3)**2*
     -     (-(s16*za(jc,jd)*zb(jb,ju)*zb(jc,je)) + 
     -       (-s12 + s45 + s46 + s56)*za(jb,jd)*zb(jb,ju)*zb(je,jb) - 
     -       s12*za(jn,jd)*zb(je,jb)*zb(jn,ju) + 
     -       (s12 + s13 + s16 - s45 - s46 - s56)*
     -        (za(jc,jd)*zb(jb,ju)*zb(jc,je) + 
     -          za(jn,jd)*(zb(je,ju)*zb(jn,jb) - zb(je,jb)*zb(jn,ju)))) + 
     -    16*prop16*tagp*real(anomc3)**2*
     -     (-(s16*za(jc,jd)*zb(jb,ju)*zb(jc,je)) + 
     -       (-s12 + s45 + s46 + s56)*za(jb,jd)*zb(jb,ju)*zb(je,jb) - 
     -       s12*za(jn,jd)*zb(je,jb)*zb(jn,ju) + 
     -       (s12 + s13 + s16 - s45 - s46 - s56)*
     -        (za(jc,jd)*zb(jb,ju)*zb(jc,je) + 
     -          za(jn,jd)*(zb(je,ju)*zb(jn,jb) - zb(je,jb)*zb(jn,ju)))))

        ! ++
        ampsAnomL4(5) = (gw**4*prop126*prop16*prop34*tagm*tagp*
     -    (8*im*mt*(s126 - s16 + s23 + s25 + s35)*imag(anomc4)*real(anomc2)*
     -       ((s12 - s126 + s14 + s16)*za(jb,jd)*za(jn,jc)*zb(jc,je)*zb(jc,ju) + 
     -         (s12 - s126 + s14 + s16)*za(jb,jd)*za(jb,jn)*zb(jc,ju)*zb(je,jb) + 
     -         (s12 - s126 + s16)*za(jb,jn)*za(je,jd)*zb(jc,je)*zb(je,ju) + 
     -         (s12 - s126 + s16)*za(jb,jn)*za(jn,jd)*zb(jc,jn)*zb(je,ju) + 
     -         (s12 - s126 + s13 + s14 + s16 + s23 + s34 + s35)*za(jb,jd)*za(ju,jn)*
     -          zb(jc,ju)*zb(je,ju) + 
     -         (s12 + s13 + s14 + s16)*za(jb,jd)*za(jn,jd)*zb(jd,jc)*zb(je,ju)) - 
     -      8*im*mt*(s126 - s16 + s23 + s25 + s35)*imag(anomc2)*real(anomc4)*
     -       ((s12 - s126 + s14 + s16)*za(jb,jd)*za(jn,jc)*zb(jc,je)*zb(jc,ju) + 
     -         (s12 - s126 + s14 + s16)*za(jb,jd)*za(jb,jn)*zb(jc,ju)*zb(je,jb) + 
     -         (s12 - s126 + s16)*za(jb,jn)*za(je,jd)*zb(jc,je)*zb(je,ju) + 
     -         (s12 - s126 + s16)*za(jb,jn)*za(jn,jd)*zb(jc,jn)*zb(je,ju) + 
     -         (s12 - s126 + s13 + s14 + s16 + s23 + s34 + s35)*za(jb,jd)*za(ju,jn)*
     -          zb(jc,ju)*zb(je,ju) + 
     -         (s12 + s13 + s14 + s16)*za(jb,jd)*za(jn,jd)*zb(jd,jc)*zb(je,ju)) - 
     -      8*mt*(s126 - s16 + s23 + s25 + s35)*imag(anomc2)*imag(anomc4)*
     -       ((s12 - s126 + s16)*za(jb,jn)*(za(je,jd)*zb(jc,je) + za(jn,jd)*zb(jc,jn))*
     -          zb(je,ju) + za(jb,jd)*
     -          ((-s12 + s126 + s14 - s16)*za(jn,jc)*zb(jc,je)*zb(jc,ju) + 
     -            (-s12 + s126 + s14 - s16)*za(jb,jn)*zb(jc,ju)*zb(je,jb) + 
     -            ((-s12 + s126 + s13 + s14 - s16 + s23 + s34 + s35)*za(ju,jn)*
     -                zb(jc,ju) + (-s12 + 2*s126 + s13 + s14 - s16)*za(jn,jd)*zb(jd,jc))*
     -             zb(je,ju))) - 8*mt*(s126 - s16 + s23 + s25 + s35)*real(anomc2)*
     -       real(anomc4)*((s12 - s126 + s16)*za(jb,jn)*
     -          (za(je,jd)*zb(jc,je) + za(jn,jd)*zb(jc,jn))*zb(je,ju) + 
     -         za(jb,jd)*((-s12 + s126 + s14 - s16)*za(jn,jc)*zb(jc,je)*zb(jc,ju) + 
     -            (-s12 + s126 + s14 - s16)*za(jb,jn)*zb(jc,ju)*zb(je,jb) + 
     -            ((-s12 + s126 + s13 + s14 - s16 + s23 + s34 + s35)*za(ju,jn)*
     -                zb(jc,ju) + (-s12 + 2*s126 + s13 + s14 - s16)*za(jn,jd)*zb(jd,jc))*
     -             zb(je,ju))) - 2*(s12 - s126 + s16)*imag(anomc2)**2*
     -       (((s126 - s16 + s23 + s25 + s35)*za(jb,jn)*
     -             (za(je,jd)*zb(jc,je) + za(jn,jd)*zb(jc,jn)) + 
     -            (s12 - 2*s126 + 2*s16 - s23 - s25 - s35)*za(jn,jd)*za(ju,jb)*zb(jc,ju))*
     -          zb(je,ju) + za(jb,jd)*
     -          (-((s126 - s14 - s16 + s23 + s25 + s35)*za(jn,jc)*zb(jc,je)*zb(jc,ju)) + 
     -            s14*za(jb,jn)*zb(jc,ju)*zb(je,jb) + 
     -            ((s13 + s14 + s23 + s34 + s35)*za(ju,jn)*zb(jc,ju) + 
     -               (s12 - s126 + s13 + s14 + 2*s16 - s23 - s25 - s35)*za(jn,jd)*
     -                zb(jd,jc))*zb(je,ju))) - 
     -      2*(s12 - s126 + s16)*real(anomc2)**2*
     -       (((s126 - s16 + s23 + s25 + s35)*za(jb,jn)*
     -             (za(je,jd)*zb(jc,je) + za(jn,jd)*zb(jc,jn)) + 
     -            (s12 - 2*s126 + 2*s16 - s23 - s25 - s35)*za(jn,jd)*za(ju,jb)*zb(jc,ju))*
     -          zb(je,ju) + za(jb,jd)*
     -          (-((s126 - s14 - s16 + s23 + s25 + s35)*za(jn,jc)*zb(jc,je)*zb(jc,ju)) + 
     -            s14*za(jb,jn)*zb(jc,ju)*zb(je,jb) + 
     -            ((s13 + s14 + s23 + s34 + s35)*za(ju,jn)*zb(jc,ju) + 
     -               (s12 - s126 + s13 + s14 + 2*s16 - s23 - s25 - s35)*za(jn,jd)*
     -                zb(jd,jc))*zb(je,ju))) + 
     -      16*imag(anomc4)**2*(-2*(s12 - s126 + s16)*
     -          ((s126 - s16 + s23 + s25 + s35)*za(jb,jn)*
     -             ((-2*s12 + s126 - 2*s14 - 2*s16 + s23 + s25 + s35)*za(je,jd)*
     -                zb(jc,je) + (-s12 - 2*s14 - s16 + s23 + s25 + s35)*za(jn,jd)*
     -                zb(jc,jn)) + (s12**2 - 2*s12*(s126 - s16) + 
     -               (s126 - s16)*(2*s126 - 2*s16 + s23 + s25 + s35))*za(jn,jd)*za(ju,jb)*
     -             zb(jc,ju))*zb(je,ju) - 
     -         2*za(jb,jd)*((s126**3 + s12**2*s14 + 2*s14**2*s16 + 2*s14*s16**2 - 
     -               s16**3 - 2*s14**2*s23 - 2*s14*s16*s23 + 2*s16**2*s23 + s14*s23**2 - 
     -               s16*s23**2 - 2*s14**2*s25 - 2*s14*s16*s25 + 2*s16**2*s25 + 
     -               2*s14*s23*s25 - 2*s16*s23*s25 + s14*s25**2 - s16*s25**2 - 
     -               2*s14**2*s35 - 2*s14*s16*s35 + s16**2*s35 + 2*s14*s23*s35 - 
     -               s16*s23*s35 + 2*s14*s25*s35 - s16*s25*s35 + s14*s35**2 + 
     -               s126**2*(s14 - 3*s16 + 2*s23 + 2*s25 + s35) + 
     -               s126*(-2*s14**2 + 3*s16**2 + (s23 + s25)*(s23 + s25 + s35) + 
     -                  s14*(-3*s16 + s23 + s25 + s35) - 2*s16*(2*s23 + 2*s25 + s35)) - 
     -               s12*(s126**2 - 2*s14*s16 + 
     -                  (s16 - s23 - s25)*(s16 - s23 - s25 - s35) + 
     -                  s126*(2*s14 - 2*s16 + 2*s23 + 2*s25 + s35)))*za(jn,jc)*zb(jc,je)*
     -             zb(jc,ju) + ((2*s126**3 + 2*s13*s14*s16 + 2*s14**2*s16 + 
     -                  2*s13*s16**2 + 4*s14*s16**2 - 2*s13*s14*s23 - 2*s14**2*s23 - 
     -                  2*s13*s16*s23 - 2*s14*s16*s23 + 2*s16**2*s23 + s13*s23**2 - 
     -                  s14*s23**2 - 2*s16*s23**2 + s23**3 - 2*s13*s14*s25 - 
     -                  2*s14**2*s25 - 2*s13*s16*s25 - 4*s14*s16*s25 + 2*s13*s23*s25 - 
     -                  2*s16*s23*s25 + 2*s23**2*s25 + s13*s25**2 + s14*s25**2 + 
     -                  s23*s25**2 + 2*s14*s16*s34 + 2*s16**2*s34 - 2*s14*s23*s34 - 
     -                  2*s16*s23*s34 + s23**2*s34 - 2*s14*s25*s34 - 2*s16*s25*s34 + 
     -                  2*s23*s25*s34 + s25**2*s34 - 2*s13*s14*s35 - 2*s14**2*s35 - 
     -                  2*s13*s16*s35 - 2*s14*s16*s35 + s16**2*s35 + 2*s13*s23*s35 - 
     -                  2*s14*s23*s35 - 3*s16*s23*s35 + 3*s23**2*s35 + 2*s13*s25*s35 - 
     -                  s16*s25*s35 + 4*s23*s25*s35 + s25**2*s35 - 2*s14*s34*s35 - 
     -                  2*s16*s34*s35 + 2*s23*s34*s35 + 2*s25*s34*s35 + s13*s35**2 - 
     -                  s14*s35**2 - s16*s35**2 + 3*s23*s35**2 + 2*s25*s35**2 + 
     -                  s34*s35**2 + s35**3 + s12**2*(s13 + s14 + s23 + s34 + s35) + 
     -                  s126**2*(s13 + 3*s14 - 4*s16 + 3*s23 + 2*s25 + s34 + 2*s35) + 
     -                  s126*(-2*s14**2 + 2*s16**2 - 5*s16*s23 + s23**2 - 2*s16*s25 + 
     -                     s23*s25 - 3*s16*s34 + s23*s34 + s25*s34 - 3*s16*s35 + 
     -                     s23*s35 + s34*s35 + s13*(-2*s14 - 3*s16 + s23 + s25 + s35) + 
     -                     s14*(-7*s16 + s23 + 3*s25 - 2*s34 + s35)) + 
     -                  s12*(-2*s126**2 + 2*s13*s16 + 4*s14*s16 - 2*s14*s23 + 2*s16*s23 - 
     -                     2*s14*s25 + 2*s16*s34 - 2*s14*s35 + s16*s35 + s23*s35 + 
     -                     s25*s35 + s35**2 - 
     -                     s126*(2*s13 + 4*s14 - 2*s16 + 4*s23 + 2*s25 + 2*s34 + 3*s35)))*
     -                za(ju,jn)*zb(jc,ju) + 
     -               (s12**3 + s126**3 + 2*s13*s14*s16 + 2*s14**2*s16 + 2*s13*s16**2 + 
     -                  4*s14*s16**2 + 2*s16**3 + s12**2*(-2*s126 + s13 + s14 + 3*s16) - 
     -                  2*s13*s14*s23 - 2*s14**2*s23 - 2*s13*s16*s23 - 4*s14*s16*s23 - 
     -                  s16**2*s23 + s13*s23**2 + s14*s23**2 - 2*s13*s14*s25 - 
     -                  2*s14**2*s25 - 2*s13*s16*s25 - 4*s14*s16*s25 - s16**2*s25 + 
     -                  2*s13*s23*s25 + 2*s14*s23*s25 + s13*s25**2 + s14*s25**2 - 
     -                  2*s13*s14*s35 - 2*s14**2*s35 - 2*s13*s16*s35 - 4*s14*s16*s35 - 
     -                  2*s16**2*s35 + 2*s13*s23*s35 + 2*s14*s23*s35 + s16*s23*s35 + 
     -                  2*s13*s25*s35 + 2*s14*s25*s35 + s16*s25*s35 + s13*s35**2 + 
     -                  s14*s35**2 + s16*s35**2 + 
     -                  s126**2*(s13 + s14 - s16 + 2*s23 + 2*s25 + s35) - 
     -                  s12*(-2*s13*s16 - 4*s14*s16 - 4*s16**2 + 2*s14*s23 + s16*s23 + 
     -                     2*s14*s25 + s16*s25 + 
     -                     s126*(2*s13 + 4*s14 + 4*s16 + s23 + s25) + 2*s14*s35 + 
     -                     2*s16*s35 - s23*s35 - s25*s35 - s35**2) + 
     -                  s126*(-2*s14**2 - 2*s16**2 - 2*s16*s23 + s23**2 - 2*s16*s25 + 
     -                     2*s23*s25 + s25**2 + s23*s35 + s25*s35 + 
     -                     s14*(-5*s16 + s23 + s25 + s35) + 
     -                     s13*(-2*s14 - 3*s16 + s23 + s25 + s35)))*za(jn,jd)*zb(jd,jc))*
     -             zb(je,ju) + za(jb,jn)*
     -             ((s126**3 + s12**2*s14 + 2*s14**2*s16 + 4*s14*s16**2 - 2*s14**2*s23 - 
     -                  4*s14*s16*s23 + s14*s23**2 - 2*s14**2*s25 - 4*s14*s16*s25 + 
     -                  2*s14*s23*s25 + s14*s25**2 + 
     -                  s126**2*(3*s14 - 2*s16 + s23 + s25) - 2*s14**2*s35 - 
     -                  4*s14*s16*s35 - s16**2*s35 + 2*s14*s23*s35 + s16*s23*s35 + 
     -                  2*s14*s25*s35 + s16*s25*s35 + s14*s35**2 + s16*s35**2 - 
     -                  s126*(2*s14**2 + 7*s14*s16 - s16**2 + s16*(s23 + s25 - s35) - 
     -                     3*s14*(s23 + s25 + s35) + s35*(s23 + s25 + s35)) - 
     -                  s12*(s126**2 + s126*(4*s14 - s16 + s23 + s25) + 
     -                     (s16 - s23 - s25 - s35)*s35 + 2*s14*(-2*s16 + s23 + s25 + s35))
     -                  )*zb(jc,ju)*zb(je,jb) + 
     -               (s12 - s126 + s16)*(s126 - s16 + s23 + s25 + s35)*za(ju,jn)*
     -                zb(jb,ju)*zb(jc,jn)*zb(je,ju)))) + 
     -      16*real(anomc4)**2*(-2*(s12 - s126 + s16)*
     -          ((s126 - s16 + s23 + s25 + s35)*za(jb,jn)*
     -             ((-2*s12 + s126 - 2*s14 - 2*s16 + s23 + s25 + s35)*za(je,jd)*
     -                zb(jc,je) + (-s12 - 2*s14 - s16 + s23 + s25 + s35)*za(jn,jd)*
     -                zb(jc,jn)) + (s12**2 - 2*s12*(s126 - s16) + 
     -               (s126 - s16)*(2*s126 - 2*s16 + s23 + s25 + s35))*za(jn,jd)*za(ju,jb)*
     -             zb(jc,ju))*zb(je,ju) - 
     -         2*za(jb,jd)*((s126**3 + s12**2*s14 + 2*s14**2*s16 + 2*s14*s16**2 - 
     -               s16**3 - 2*s14**2*s23 - 2*s14*s16*s23 + 2*s16**2*s23 + s14*s23**2 - 
     -               s16*s23**2 - 2*s14**2*s25 - 2*s14*s16*s25 + 2*s16**2*s25 + 
     -               2*s14*s23*s25 - 2*s16*s23*s25 + s14*s25**2 - s16*s25**2 - 
     -               2*s14**2*s35 - 2*s14*s16*s35 + s16**2*s35 + 2*s14*s23*s35 - 
     -               s16*s23*s35 + 2*s14*s25*s35 - s16*s25*s35 + s14*s35**2 + 
     -               s126**2*(s14 - 3*s16 + 2*s23 + 2*s25 + s35) + 
     -               s126*(-2*s14**2 + 3*s16**2 + (s23 + s25)*(s23 + s25 + s35) + 
     -                  s14*(-3*s16 + s23 + s25 + s35) - 2*s16*(2*s23 + 2*s25 + s35)) - 
     -               s12*(s126**2 - 2*s14*s16 + 
     -                  (s16 - s23 - s25)*(s16 - s23 - s25 - s35) + 
     -                  s126*(2*s14 - 2*s16 + 2*s23 + 2*s25 + s35)))*za(jn,jc)*zb(jc,je)*
     -             zb(jc,ju) + ((2*s126**3 + 2*s13*s14*s16 + 2*s14**2*s16 + 
     -                  2*s13*s16**2 + 4*s14*s16**2 - 2*s13*s14*s23 - 2*s14**2*s23 - 
     -                  2*s13*s16*s23 - 2*s14*s16*s23 + 2*s16**2*s23 + s13*s23**2 - 
     -                  s14*s23**2 - 2*s16*s23**2 + s23**3 - 2*s13*s14*s25 - 
     -                  2*s14**2*s25 - 2*s13*s16*s25 - 4*s14*s16*s25 + 2*s13*s23*s25 - 
     -                  2*s16*s23*s25 + 2*s23**2*s25 + s13*s25**2 + s14*s25**2 + 
     -                  s23*s25**2 + 2*s14*s16*s34 + 2*s16**2*s34 - 2*s14*s23*s34 - 
     -                  2*s16*s23*s34 + s23**2*s34 - 2*s14*s25*s34 - 2*s16*s25*s34 + 
     -                  2*s23*s25*s34 + s25**2*s34 - 2*s13*s14*s35 - 2*s14**2*s35 - 
     -                  2*s13*s16*s35 - 2*s14*s16*s35 + s16**2*s35 + 2*s13*s23*s35 - 
     -                  2*s14*s23*s35 - 3*s16*s23*s35 + 3*s23**2*s35 + 2*s13*s25*s35 - 
     -                  s16*s25*s35 + 4*s23*s25*s35 + s25**2*s35 - 2*s14*s34*s35 - 
     -                  2*s16*s34*s35 + 2*s23*s34*s35 + 2*s25*s34*s35 + s13*s35**2 - 
     -                  s14*s35**2 - s16*s35**2 + 3*s23*s35**2 + 2*s25*s35**2 + 
     -                  s34*s35**2 + s35**3 + s12**2*(s13 + s14 + s23 + s34 + s35) + 
     -                  s126**2*(s13 + 3*s14 - 4*s16 + 3*s23 + 2*s25 + s34 + 2*s35) + 
     -                  s126*(-2*s14**2 + 2*s16**2 - 5*s16*s23 + s23**2 - 2*s16*s25 + 
     -                     s23*s25 - 3*s16*s34 + s23*s34 + s25*s34 - 3*s16*s35 + 
     -                     s23*s35 + s34*s35 + s13*(-2*s14 - 3*s16 + s23 + s25 + s35) + 
     -                     s14*(-7*s16 + s23 + 3*s25 - 2*s34 + s35)) + 
     -                  s12*(-2*s126**2 + 2*s13*s16 + 4*s14*s16 - 2*s14*s23 + 2*s16*s23 - 
     -                     2*s14*s25 + 2*s16*s34 - 2*s14*s35 + s16*s35 + s23*s35 + 
     -                     s25*s35 + s35**2 - 
     -                     s126*(2*s13 + 4*s14 - 2*s16 + 4*s23 + 2*s25 + 2*s34 + 3*s35)))*
     -                za(ju,jn)*zb(jc,ju) + 
     -               (s12**3 + s126**3 + 2*s13*s14*s16 + 2*s14**2*s16 + 2*s13*s16**2 + 
     -                  4*s14*s16**2 + 2*s16**3 + s12**2*(-2*s126 + s13 + s14 + 3*s16) - 
     -                  2*s13*s14*s23 - 2*s14**2*s23 - 2*s13*s16*s23 - 4*s14*s16*s23 - 
     -                  s16**2*s23 + s13*s23**2 + s14*s23**2 - 2*s13*s14*s25 - 
     -                  2*s14**2*s25 - 2*s13*s16*s25 - 4*s14*s16*s25 - s16**2*s25 + 
     -                  2*s13*s23*s25 + 2*s14*s23*s25 + s13*s25**2 + s14*s25**2 - 
     -                  2*s13*s14*s35 - 2*s14**2*s35 - 2*s13*s16*s35 - 4*s14*s16*s35 - 
     -                  2*s16**2*s35 + 2*s13*s23*s35 + 2*s14*s23*s35 + s16*s23*s35 + 
     -                  2*s13*s25*s35 + 2*s14*s25*s35 + s16*s25*s35 + s13*s35**2 + 
     -                  s14*s35**2 + s16*s35**2 + 
     -                  s126**2*(s13 + s14 - s16 + 2*s23 + 2*s25 + s35) - 
     -                  s12*(-2*s13*s16 - 4*s14*s16 - 4*s16**2 + 2*s14*s23 + s16*s23 + 
     -                     2*s14*s25 + s16*s25 + 
     -                     s126*(2*s13 + 4*s14 + 4*s16 + s23 + s25) + 2*s14*s35 + 
     -                     2*s16*s35 - s23*s35 - s25*s35 - s35**2) + 
     -                  s126*(-2*s14**2 - 2*s16**2 - 2*s16*s23 + s23**2 - 2*s16*s25 + 
     -                     2*s23*s25 + s25**2 + s23*s35 + s25*s35 + 
     -                     s14*(-5*s16 + s23 + s25 + s35) + 
     -                     s13*(-2*s14 - 3*s16 + s23 + s25 + s35)))*za(jn,jd)*zb(jd,jc))*
     -             zb(je,ju) + za(jb,jn)*
     -             ((s126**3 + s12**2*s14 + 2*s14**2*s16 + 4*s14*s16**2 - 2*s14**2*s23 - 
     -                  4*s14*s16*s23 + s14*s23**2 - 2*s14**2*s25 - 4*s14*s16*s25 + 
     -                  2*s14*s23*s25 + s14*s25**2 + 
     -                  s126**2*(3*s14 - 2*s16 + s23 + s25) - 2*s14**2*s35 - 
     -                  4*s14*s16*s35 - s16**2*s35 + 2*s14*s23*s35 + s16*s23*s35 + 
     -                  2*s14*s25*s35 + s16*s25*s35 + s14*s35**2 + s16*s35**2 - 
     -                  s126*(2*s14**2 + 7*s14*s16 - s16**2 + s16*(s23 + s25 - s35) - 
     -                     3*s14*(s23 + s25 + s35) + s35*(s23 + s25 + s35)) - 
     -                  s12*(s126**2 + s126*(4*s14 - s16 + s23 + s25) + 
     -                     (s16 - s23 - s25 - s35)*s35 + 2*s14*(-2*s16 + s23 + s25 + s35))
     -                  )*zb(jc,ju)*zb(je,jb) + 
     -               (s12 - s126 + s16)*(s126 - s16 + s23 + s25 + s35)*za(ju,jn)*
     -                zb(jb,ju)*zb(jc,jn)*zb(je,ju))))))/
     -  (2._dp*(s12 - s126 + s16)*(s126 - s16 + s23 + s25 + s35))

        ! PM
        ampsAnomL2(2) = gw**4*prop126*prop16*prop34*tagm*tagp*
     -  (-(im*mt*imag(anomc2)*za(jn,jd)*zb(jb,ju)*zb(jc,je)) + 
     -    mt*real(anomc2)*za(jn,jd)*zb(jb,ju)*zb(jc,je) + 
     -    4*(im*imag(anomc4) - real(anomc4))*
     -     (za(jb,jd)*zb(jb,ju)*(za(jb,jn)*zb(jc,jb)*zb(je,jb) + 
     -          za(jn,jd)*zb(jd,jc)*zb(je,jb) + 
     -          za(ju,jn)*(zb(jc,ju)*zb(je,jb) + zb(jc,jb)*zb(je,ju))) + 
     -       zb(je,ju)*(za(ju,jd)*za(ju,jn)*zb(jb,ju)*zb(jc,ju) + 
     -          za(jn,jd)*(-((s12 + s13 + s16)*zb(jc,jb)) + 
     -             za(ju,jd)*zb(jc,ju)*zb(jd,jb) + 
     -             za(ju,jn)*(zb(jb,ju)*zb(jc,jn) + zb(jc,ju)*zb(jn,jb))))))

        ampsAnomL4(2) = gw**4*prop126*prop34*tagm*(-(im*mt*prop16*tagp*imag(anomc2)*real(anomc1)*
     -       za(jn,jd)*zb(jb,ju)*zb(jc,je)) + 
     -    mt*prop16*tagp*real(anomc1)*real(anomc2)*za(jn,jd)*zb(jb,ju)*zb(jc,je) - 
     -    4*im*mt*imag(anomc2)*real(anomc8)*za(jn,jd)*zb(jb,ju)*zb(jc,je) + 
     -    4*mt*real(anomc2)*real(anomc8)*za(jn,jd)*zb(jb,ju)*zb(jc,je) - 
     -    4*prop16*tagp*imag(anomc2)*imag(anomc3)*
     -     (za(jn,jd)*(zb(jb,ju)*((s126 - s16)*zb(jc,je) + 
     -             za(jb,jd)*zb(jd,jc)*zb(je,jb)) + 
     -          (s16*zb(jc,jb) - za(ju,jd)*zb(jc,ju)*zb(jd,jb))*zb(je,ju)) + 
     -       zb(jb,ju)*(za(jb,jd)*(za(jn,jc)*zb(jc,jb)*zb(jc,je) + 
     -             (za(jb,jn)*zb(jc,jb) + za(ju,jn)*zb(jc,ju))*zb(je,jb)) - 
     -          za(ju,jd)*zb(jc,ju)*(za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb) + 
     -             za(ju,jn)*zb(je,ju)))) - 
     -    4*im*prop16*tagp*imag(anomc3)*real(anomc2)*
     -     (za(jn,jd)*(zb(jb,ju)*((s126 - s16)*zb(jc,je) + 
     -             za(jb,jd)*zb(jd,jc)*zb(je,jb)) + 
     -          (s16*zb(jc,jb) - za(ju,jd)*zb(jc,ju)*zb(jd,jb))*zb(je,ju)) + 
     -       zb(jb,ju)*(za(jb,jd)*(za(jn,jc)*zb(jc,jb)*zb(jc,je) + 
     -             (za(jb,jn)*zb(jc,jb) + za(ju,jn)*zb(jc,ju))*zb(je,jb)) - 
     -          za(ju,jd)*zb(jc,ju)*(za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb) + 
     -             za(ju,jn)*zb(je,ju)))) - 
     -    4*im*prop16*tagp*imag(anomc2)*real(anomc3)*
     -     (za(jn,jd)*(zb(jb,ju)*((s126 - s16)*zb(jc,je) + 
     -             za(jb,jd)*zb(jd,jc)*zb(je,jb)) + 
     -          (s16*zb(jc,jb) - za(ju,jd)*zb(jc,ju)*zb(jd,jb))*zb(je,ju)) + 
     -       zb(jb,ju)*(za(jb,jd)*(za(jn,jc)*zb(jc,jb)*zb(jc,je) + 
     -             (za(jb,jn)*zb(jc,jb) + za(ju,jn)*zb(jc,ju))*zb(je,jb)) - 
     -          za(ju,jd)*zb(jc,ju)*(za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb) + 
     -             za(ju,jn)*zb(je,ju)))) + 
     -    4*prop16*tagp*real(anomc2)*real(anomc3)*
     -     (za(jn,jd)*(zb(jb,ju)*((s126 - s16)*zb(jc,je) + 
     -             za(jb,jd)*zb(jd,jc)*zb(je,jb)) + 
     -          (s16*zb(jc,jb) - za(ju,jd)*zb(jc,ju)*zb(jd,jb))*zb(je,ju)) + 
     -       zb(jb,ju)*(za(jb,jd)*(za(jn,jc)*zb(jc,jb)*zb(jc,je) + 
     -             (za(jb,jn)*zb(jc,jb) + za(ju,jn)*zb(jc,ju))*zb(je,jb)) - 
     -          za(ju,jd)*zb(jc,ju)*(za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb) + 
     -             za(ju,jn)*zb(je,ju)))) - 
     -    16*mt*prop16*tagp*imag(anomc3)*imag(anomc4)*zb(je,ju)*
     -     (-(za(ju,jn)*zb(jb,ju)*(za(jb,jd)*zb(jc,jb) + za(ju,jd)*zb(jc,ju))) + 
     -       za(jn,jd)*((s12 + s13 + s16)*zb(jc,jb) - za(ju,jd)*zb(jc,ju)*zb(jd,jb) - 
     -          za(ju,jn)*(zb(jb,ju)*zb(jc,jn) + zb(jc,ju)*zb(jn,jb)))) - 
     -    16*im*mt*prop16*tagp*imag(anomc4)*real(anomc3)*zb(je,ju)*
     -     (-(za(ju,jn)*zb(jb,ju)*(za(jb,jd)*zb(jc,jb) + za(ju,jd)*zb(jc,ju))) + 
     -       za(jn,jd)*((s12 + s13 + s16)*zb(jc,jb) - za(ju,jd)*zb(jc,ju)*zb(jd,jb) - 
     -          za(ju,jn)*(zb(jb,ju)*zb(jc,jn) + zb(jc,ju)*zb(jn,jb)))) - 
     -    16*im*mt*prop16*tagp*imag(anomc3)*real(anomc4)*zb(je,ju)*
     -     (-(za(ju,jn)*zb(jb,ju)*(za(jb,jd)*zb(jc,jb) + za(ju,jd)*zb(jc,ju))) + 
     -       za(jn,jd)*((s12 + s13 + s16)*zb(jc,jb) - za(ju,jd)*zb(jc,ju)*zb(jd,jb) - 
     -          za(ju,jn)*(zb(jb,ju)*zb(jc,jn) + zb(jc,ju)*zb(jn,jb)))) + 
     -    16*mt*prop16*tagp*real(anomc3)*real(anomc4)*zb(je,ju)*
     -     (-(za(ju,jn)*zb(jb,ju)*(za(jb,jd)*zb(jc,jb) + za(ju,jd)*zb(jc,ju))) + 
     -       za(jn,jd)*((s12 + s13 + s16)*zb(jc,jb) - za(ju,jd)*zb(jc,ju)*zb(jd,jb) - 
     -          za(ju,jn)*(zb(jb,ju)*zb(jc,jn) + zb(jc,ju)*zb(jn,jb)))) + 
     -    2*prop16*tagp*real(anomc1)*real(anomc4)*
     -     (-2*za(jb,jd)*zb(jb,ju)*(za(jb,jn)*zb(jc,jb)*zb(je,jb) + 
     -          za(jn,jd)*zb(jd,jc)*zb(je,jb) + 
     -          za(ju,jn)*(zb(jc,ju)*zb(je,jb) + zb(jc,jb)*zb(je,ju))) + 
     -       2*zb(je,ju)*(-(za(ju,jd)*za(ju,jn)*zb(jb,ju)*zb(jc,ju)) + 
     -          za(jn,jd)*((s12 + s13 + s16)*zb(jc,jb) - za(ju,jd)*zb(jc,ju)*zb(jd,jb) - 
     -             za(ju,jn)*(zb(jb,ju)*zb(jc,jn) + zb(jc,ju)*zb(jn,jb))))) + 
     -    8*real(anomc4)*real(anomc8)*
     -     (-2*za(jb,jd)*zb(jb,ju)*(za(jb,jn)*zb(jc,jb)*zb(je,jb) + 
     -          za(jn,jd)*zb(jd,jc)*zb(je,jb) + 
     -          za(ju,jn)*(zb(jc,ju)*zb(je,jb) + zb(jc,jb)*zb(je,ju))) + 
     -       2*zb(je,ju)*(-(za(ju,jd)*za(ju,jn)*zb(jb,ju)*zb(jc,ju)) + 
     -          za(jn,jd)*((s12 + s13 + s16)*zb(jc,jb) - za(ju,jd)*zb(jc,ju)*zb(jd,jb) - 
     -             za(ju,jn)*(zb(jb,ju)*zb(jc,jn) + zb(jc,ju)*zb(jn,jb))))) + 
     -    4*im*prop16*tagp*imag(anomc4)*real(anomc1)*
     -     (za(jb,jd)*zb(jb,ju)*(za(jb,jn)*zb(jc,jb)*zb(je,jb) + 
     -          za(jn,jd)*zb(jd,jc)*zb(je,jb) + 
     -          za(ju,jn)*(zb(jc,ju)*zb(je,jb) + zb(jc,jb)*zb(je,ju))) + 
     -       zb(je,ju)*(za(ju,jd)*za(ju,jn)*zb(jb,ju)*zb(jc,ju) + 
     -          za(jn,jd)*(-((s12 + s13 + s16)*zb(jc,jb)) + 
     -             za(ju,jd)*zb(jc,ju)*zb(jd,jb) + 
     -             za(ju,jn)*(zb(jb,ju)*zb(jc,jn) + zb(jc,ju)*zb(jn,jb))))) + 
     -    16*im*imag(anomc4)*real(anomc8)*
     -     (za(jb,jd)*zb(jb,ju)*(za(jb,jn)*zb(jc,jb)*zb(je,jb) + 
     -          za(jn,jd)*zb(jd,jc)*zb(je,jb) + 
     -          za(ju,jn)*(zb(jc,ju)*zb(je,jb) + zb(jc,jb)*zb(je,ju))) + 
     -       zb(je,ju)*(za(ju,jd)*za(ju,jn)*zb(jb,ju)*zb(jc,ju) + 
     -          za(jn,jd)*(-((s12 + s13 + s16)*zb(jc,jb)) + 
     -             za(ju,jd)*zb(jc,ju)*zb(jd,jb) + 
     -             za(ju,jn)*(zb(jb,ju)*zb(jc,jn) + zb(jc,ju)*zb(jn,jb))))))

        ! MP
        ampsAnomL2(3) = -(gw**4*prop126*prop16*prop34*tagm*tagp*za(jb,jd)*
     -    (mt*(im*imag(anomc2) + real(anomc2))*za(jn,jc)*zb(je,ju) + 
     -      4*im*imag(anomc4)*(za(jb,jc)*zb(jb,ju)*
     -          (za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb)) + 
     -         ((s126 - s16)*za(jn,jc) + 
     -            (za(jn,jd)*za(ju,jc) - za(jc,jd)*za(ju,jn))*zb(jd,ju))*zb(je,ju) + 
     -         za(jb,jn)*zb(jb,ju)*(za(jc,jd)*zb(jd,je) + za(ju,jc)*zb(je,ju))) + 
     -      4*real(anomc4)*(za(jb,jc)*zb(jb,ju)*
     -          (za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb)) + 
     -         ((s126 - s16)*za(jn,jc) + 
     -            (za(jn,jd)*za(ju,jc) - za(jc,jd)*za(ju,jn))*zb(jd,ju))*zb(je,ju) + 
     -         za(jb,jn)*zb(jb,ju)*(za(jc,jd)*zb(jd,je) + za(ju,jc)*zb(je,ju)))))

        ampsAnomL4(3) = -(gw**4*prop126*prop16*prop34*tagm*tagp*za(jb,jd)*
     -    (real(anomc1)*(mt*(im*imag(anomc2) + real(anomc2))*za(jn,jc)*zb(je,ju) + 
     -         4*real(anomc4)*(za(jb,jc)*zb(jb,ju)*
     -             (za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb)) + 
     -            ((s126 - s16)*za(jn,jc) + 
     -               (za(jn,jd)*za(ju,jc) - za(jc,jd)*za(ju,jn))*zb(jd,ju))*zb(je,ju) + 
     -            za(jb,jn)*zb(jb,ju)*(za(jc,jd)*zb(jd,je) + za(ju,jc)*zb(je,ju)))) + 
     -      4*imag(anomc4)*(im*real(anomc1)*
     -          (za(jb,jc)*zb(jb,ju)*(za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb)) + 
     -            ((s126 - s16)*za(jn,jc) + 
     -               (za(jn,jd)*za(ju,jc) - za(jc,jd)*za(ju,jn))*zb(jd,ju))*zb(je,ju) + 
     -            za(jb,jn)*zb(jb,ju)*(za(jc,jd)*zb(jd,je) + za(ju,jc)*zb(je,ju))) + 
     -         4*mt*(imag(anomc3) - im*real(anomc3))*za(jn,jc)*
     -          (za(jn,jd)*zb(jd,jn)*zb(je,ju) + 
     -            (za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb) + za(ju,jn)*zb(je,ju))*
     -             zb(jn,ju))) + 4*za(jn,jc)*
     -       (imag(anomc2)*(imag(anomc3) - im*real(anomc3))*
     -          (za(jn,jd)*zb(jd,jn)*zb(je,ju) + 
     -            (za(jn,jc)*zb(jc,je) + za(ju,jn)*zb(je,ju))*zb(jn,ju) + 
     -            za(jb,jn)*(-(zb(jb,ju)*zb(je,jn)) + zb(je,jb)*zb(jn,ju))) + 
     -         (im*imag(anomc3) + real(anomc3))*
     -          (-4*mt*real(anomc4)*(za(jn,jd)*zb(jd,jn)*zb(je,ju) + 
     -               (za(jn,jc)*zb(jc,je) + za(jb,jn)*zb(je,jb) + za(ju,jn)*zb(je,ju))*
     -                zb(jn,ju)) - real(anomc2)*
     -             (za(jn,jd)*zb(jd,jn)*zb(je,ju) + 
     -               (za(jn,jc)*zb(jc,je) + za(ju,jn)*zb(je,ju))*zb(jn,ju) + 
     -               za(jb,jn)*(-(zb(jb,ju)*zb(je,jn)) + zb(je,jb)*zb(jn,ju)))))))

        ! MPPP
        ampsAnomL2(4) = -2*gw**4*mt*prop126*prop34*tagm*real(anomc9)*za(jn,jc)*za(ju,jb)*zb(jd,je)
        
        ampsAnomL4(4) = -2*gw**4*prop126*prop34*tagm*real(anomc9)*
     -  (mt*real(anomc1)*za(jn,jc)*za(ju,jb)*zb(jd,je) - 
     -    4*(im*imag(anomc3) + real(anomc3))*
     -     (s13*za(jb,jc)*za(ju,jn)*zb(jd,je) - 
     -       za(jn,jc)*za(ju,jb)*za(ju,jn)*zb(jd,jn)*zb(je,ju) - 
     -       za(jb,jn)*(za(jn,jc)*za(ju,jb)*zb(jd,jb)*zb(je,jn) + 
     -          za(ju,jc)*za(ju,jn)*zb(jd,je)*zb(jn,ju))))

        singletop2_amp_protos = 0._dp

        if (disable_sm) ampSM = 0._dp

        ! -- amplitude, interferes with SM
        singletop2_amp_protos = singletop2_amp_protos + abs(
     &            ampSM
     &          + ampsAnomL2(1)/lambda**2
     &          + ampsAnomL4(1)/lambda**4
     &      )**2

        ! other amplitudes 1/Lambda^4 up to 1/Lambda^8
        singletop2_amp_protos = singletop2_amp_protos + sum(abs(
     &            ampsAnomL2(2:5)/lambda**2
     &          + ampsANomL4(2:5)/lambda**4
     &      )**2)

      end

      end module
