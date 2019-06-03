      module singletop2_virtamps_m
        use ddmodule
        use types

      public :: singletop2_amp_virt
      private

      contains

      function singletop2_amp_virt(p,ju,jb,jn,je,jc,jd, musqLight, musqHeavy, light, heavy)
        use types
        implicit none
        real(dp):: singletop2_amp_virt

        include 'constants.f'
        include 'nf.f'
        include 'mxpart.f'
        include 'zprods_com.f'
        include 'epinv.f'
        include 'epinv2.f'
        include 'badpoint.f'

        real(dp), intent(in) :: p(mxpart,4)
        integer, intent(in) :: ju,jd,jn,je,jc,jb
        real(dp), intent(in) :: musqLight, musqHeavy
        logical, intent(in) :: light,heavy

        ! This parameter effectively determines the minimum number of
        ! significant correct digits.
        real(dp), parameter :: cutoff = 1d-4
        real(dp), parameter :: sc = 1.44_dp

        real(dp) :: dpresult, dpresult_sc
        real(dp) :: ddresult, ddresult_sc

        ! the additional sc parameter is just for rescaling masses and widths, nothing else!
        dpresult_sc = singletop2_amp_virt_dp(p*sqrt(sc),za*sqrt(sc),zb*sqrt(sc),
     &          ju,jb,jn,je,jc,jd, musqLight*sc, musqHeavy*sc, light, heavy, sc)

        dpresult = singletop2_amp_virt_dp(p,za,zb,ju,jb,jn,je,jc,jd, musqLight, musqHeavy, light, heavy, 1._dp)

        ! first cutoff determines when to use dd-precision to check result
        if (abs(sc**2 * dpresult_sc / dpresult - 1._dp) > cutoff) then
            ddresult =  singletop2_amp_virt_dd(p,za,zb,ju,jb,jn,je,jc,jd,
     &              musqLight, musqHeavy, light, heavy, 1._dp)
            ddresult_sc =  singletop2_amp_virt_dd(p*sqrt(sc),za*sqrt(sc),zb*sqrt(sc),
     &          ju,jb,jn,je,jc,jd, musqLight*sc, musqHeavy*sc, light, heavy, sc)

            ! hard cutoff when to discard result
            if (abs(sc**2 * ddresult_sc / ddresult - 1._dp) > cutoff) then
c               write (*,*) "dpresult / ddresult", dpresult / ddresult,
c    &              abs(sc**2 * dpresult_sc / dpresult - 1._dp), abs(sc**2 * ddresult_sc / ddresult - 1._dp)
                badpoint = .true.
                return
            else
                singletop2_amp_virt = ddresult
                return
            endif
        endif

        singletop2_amp_virt = dpresult


      end function

      function singletop2_amp_virt_dp(p,za,zb,ju,jb,jn,je,jc,jd, musqLight, musqHeavy, light, heavy, sc)
        use types
        use mod_qcdloop_c
        use eftcouple
        use singletop2_realamps_m
        use singletop2_ints_new_m
        use singletop2_ints_nonres_m
        use anomcoup_tbW
      implicit none
      real(dp):: singletop2_amp_virt_dp

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'epinv.f'
      include 'epinv2.f'

      real(dp), intent(in) :: p(mxpart,4)
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
      integer, intent(in) :: ju,jd,jn,je,jc,jb
      real(dp), intent(in) :: musqLight, musqHeavy
      logical, intent(in) :: light,heavy
      real(dp), intent(in) :: sc

      complex(dp) :: mtsq

      complex(dp) :: ampSMtree(2)
      complex(dp) :: ampsAnomL2tree(5)
      complex(dp) :: ampsAnomNonresL2tree(5)

      complex(dp) :: ampsSMLight(2), ampsSMHeavy(2), ampsSMResHeavy(2)
      complex(dp) :: ampsNonresSMLight(2), ampsNonresSMHeavy(2)

      complex(dp) :: ampsAnomLight(5), ampsAnomHeavy(5)

      complex(dp) :: ampNonresTreeC4(2)
      complex(dp) :: ampNonresTreeC10(2)

      complex(dp) :: res

      complex(dp) :: anomc1_sc, anomc2_sc, anomc3_sc, anomc4_sc
      complex(dp) :: anomc6_sc, anomc7_sc, anomc8_sc, anomc9_sc

      real(dp) :: s126, s12, s15, s16, s25, s56, s14, s36
      real(dp) :: s26, s45, s46, s24, s13, s23, s34, s35
      complex(dp) :: prop16, prop34, prop25
      integer :: i,j,k,l
      real(dp) :: s,t
      s(j,k) = real(za(j,k)*zb(k,j))
      t(j,k,l) = s(j,k) + s(j,l) + s(k,l)

      anomc1_sc = anomc1 * sc
      anomc2_sc = anomc2 * sc
      anomc3_sc = anomc3 * sqrt(sc)
      anomc4_sc = anomc4 * sqrt(sc)
      anomc6_sc = anomc6 * sqrt(sc)
      anomc7_sc = anomc7 * sqrt(sc)
      anomc8_sc = anomc8
      anomc9_sc = anomc9

      ampsAnomLight = 0._dp
      ampsAnomHeavy = 0._dp
      ampsSMLight = 0._dp
      ampsSMHeavy = 0._dp
      ampsSMResHeavy = 0._dp
      ampSMtree = 0._dp
      ampsAnomL2tree = 0._dp

      ampsNonresSMLight = 0._dp
      ampsNonresSMHeavy = 0._dp

      ampNonresTreeC4 = 0._dp
      ampNonresTreeC10 = 0._dp

      mtsq = (mt**2 - im*mt*twidth)*sc

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

      prop16 = 1._dp / (s(ju,jd) - wmass**2 * sc)
      prop34 = 1._dp / (s(jn,je) - wmass**2*sc + im*wmass*wwidth*sc)
      prop25 = 1._dp / (s(jb,jc) - zmass**2*sc)

      ! resonant
      call initamp(musqHeavy, mtsq, za,zb,ju,jb,jn,je,jc,jd, epinv, epinv*epinv2)
      ! non resonant 
      call initamp_nonres(musqHeavy, musqLight, za,zb,ju,jb,jn,je,jc,jd, epinv, epinv*epinv2, sc)

      call singletop2_hamp_tree(za,zb, ju,jb,jn,je,jc,jd, mtsq, ampSMtree, ampsAnomL2tree, ampsAnomNonresL2tree, sc)

      if (light) then

          if (enable_resonant) then
              ampsSMLight(1) = (mytree*prop16*prop34*(8 + epinv*(3 + 2*epinv2)
     &          + (3 + 2*epinv)*Log(-(MusqLight/s16)) + 
     &          Log(-(MusqLight/s16))**2))/(-mtsq + s126)
          endif

          if (enable_nonresonant) then
              ! MM contributions,
              res = ampNonresonantLightFullMM()
            ampsNonresSMLight(1) = ampsNonresSMLight(1) + res

              ! PP contributions
              res = ampNonresonantLightFullPP()
            ampsNonresSMLight(2) = ampsNonresSMLight(2) + res

            ampsSMLight(1:2) = ampsSMLight(1:2) + ampsNonresSMLight(1:2)
          endif

    
          if (enable_lambda2) then
              if (enable_resonant) then
                  ampsAnomLight(1) = ampsAnomLight(1) + (2*prop34*(8 + epinv*(3 + 2*epinv2) + (3 + 2*epinv)*Log(-(MusqLight/s16)) + 
     -            Log(-(MusqLight/s16))**2)*(im*sqrt(mtsq)*prop16*(struc4Step5 + struc9Step5)*imag(anomc3_sc) - 
     -            mytree*prop16*real(anomc1_sc) - 4*sqrt(mtsq)*mytree*prop16*real(anomc3_sc) + 
     -            sqrt(mtsq)*prop16*struc4Step5*real(anomc3_sc) - sqrt(mtsq)*prop16*struc9Step5*real(anomc3_sc) - 
     -            2*mytree*real(anomc8_sc)))/(mtsq - s126)
              endif

              if (enable_nonresonant) then
                ampsAnomLight(1) = ampsAnomLight(1) + ampsNonresSMLight(1)*real(anomc1_sc)
                ampsAnomLight(5) = ampsAnomLight(5) + ampsNonresSMLight(2)*real(anomc1_sc)
              endif

          endif


          if (enable_lambda4) then
            if (enable_resonant) then
             ! PM
             ampsAnomLight(2) = (prop16*prop34*(8 + epinv*(3 + 2*epinv2) + 
     -      Log(-(MusqLight/s16))*(3 + 2*epinv + Log(-(MusqLight/s16))))*
     -    (im*Sqrt(mtsq)*struc59PM*imag(anomc2_sc) + 
     -      4*im*(2*struc14PM + 2*struc16PM - 2*struc17PM + 
     -         (s126 + 2*(s13 + s34 + s35))*struc1PM - s34*struc59PM + 
     -         mtsq*(-struc1PM + struc59PM) + 
     -         2*(struc43PM + struc56PM + 2*(struc5PM + struc8PM)))*imag(anomc4_sc) - 
     -      Sqrt(mtsq)*struc59PM*real(anomc2_sc) - 
     -      4*(2*struc14PM + 2*struc16PM - 2*struc17PM + 
     -         (s126 + 2*(s13 + s34 + s35))*struc1PM - s34*struc59PM + 
     -         mtsq*(-struc1PM + struc59PM) + 
     -         2*(struc43PM + struc56PM + 2*(struc5PM + struc8PM)))*real(anomc4_sc)))/
     -  (4*(mtsq - s126))

             ! MP
             ampsAnomLight(3) = -(prop16*prop34*(8 + epinv*(3 + 2*epinv2) + 
     -       Log(-(MusqLight/s16))*(3 + 2*epinv + Log(-(MusqLight/s16))))*
     -     (im*Sqrt(mtsq)*struc59MP*imag(anomc2_sc) - 
     -       2*im*(-4*struc17MP + 2*(mtsq - s126 + 2*(s13 + s35))*struc1MP + 
     -          (-2*mtsq + s126 + s34 - s35)*struc59MP + 
     -          2*(2*struc39MP + 2*struc41MP + struc42MP + 2*struc43MP + struc58MP - 
     -             2*struc8MP))*imag(anomc4_sc) + Sqrt(mtsq)*struc59MP*real(anomc2_sc) - 
     -       2*(-4*struc17MP + 2*(mtsq - s126 + 2*(s13 + s35))*struc1MP + 
     -          (-2*mtsq + s126 + s34 - s35)*struc59MP + 
     -          2*(2*struc39MP + 2*struc41MP + struc42MP + 2*struc43MP + struc58MP - 
     -             2*struc8MP))*real(anomc4_sc)))/(4.*(mtsq - s126))

              ampsAnomLight(4) = (-2*Sqrt(mtsq)*prop34*(8 + epinv*(3 + 2*epinv2) + 
     -      (3 + 2*epinv)*Log(-(MusqLight/s16)) + Log(-(MusqLight/s16))**2)*real(anomc9_sc)*
     -    za(jn,jc)*za(ju,jb)*zb(jd,je))/(mtsq - s126)
          endif

          endif

          if (enable_nonresonant) then

              if( real(anomc4_sc) /= 0._dp ) then
                  res = ampNonresonantLightFullReC4PM()
                  ampsAnomLight(2) = ampsAnomLight(2) + res*real(anomc4_sc)

                  res = ampNonresonantLightFullReC4MP()
                  ampsAnomLight(3) = ampsAnomLight(3) + res*real(anomc4_sc)
              endif

              if( aimag(anomc4_sc) /= 0._dp ) then
                  res = ampNonresonantLightFullImC4PM()
                  ampsAnomLight(2) = ampsAnomLight(2) + res*aimag(anomc4_sc)

                  res = ampNonresonantLightFullImC4MP()
                  ampsAnomLight(3) = ampsAnomLight(3) + res*aimag(anomc4_sc)
              endif

          endif ! enable_nonresonant
      endif

      ! in the following: 
      ! no factor of two, since we calculate in 1/(4pi), but mcfm convention is 1/(2pi)

      singletop2_amp_virt_dp = 0._dp

c assemble light amplitudes
      if (light) then
          if (.not. disable_sm) then
              singletop2_amp_virt_dp = singletop2_amp_virt_dp + sum(real(ampSMtree(:)*conjg(ampsSMLight(:))))
          endif

          if (enable_lambda2) then
              singletop2_amp_virt_dp = singletop2_amp_virt_dp 
c                       ---- contribution
     &                  + real(ampSMtree(1)*conjg(ampsAnomLight(1)))/lambda**2/sc
     &                  + real(ampsAnomL2tree(1)*conjg(ampsSMLight(1)))/lambda**2/sc
c                       ++-- contribution, for non-resonant pieces
c                           which is 5 in anom array and 2 in SM array
     &                  + real(ampSMtree(2)*conjg(ampsAnomLight(5)))/lambda**2/sc
     &                  + real(ampsAnomL2tree(5)*conjg(ampsSMLight(2)))/lambda**2/sc
          endif

          if (enable_lambda4) then
              singletop2_amp_virt_dp = singletop2_amp_virt_dp
     &                  + sum(real(ampsAnomL2tree*conjg(ampsAnomLight)))/lambda**4/sc/sc
          endif
      endif

      if (heavy) then

          if (enable_resonant) then
              res = ampHeavyMM()
              ampsSMResHeavy(1) = (res 
c mass renormalization
     &          + (2*mtsq*mytree*(4 + 3*epinv + 3*Log(MusqHeavy/mtsq)))/(mtsq - s126)**2
     &              )*prop16*prop34
          endif

          ampsSMHeavy(1) = ampsSMHeavy(1) + ampsSMResHeavy(1)

          ! non-resonant contributions
          if (enable_nonresonant) then
              res = ampNonresonantHeavyMM()
              ampsNonresSMHeavy(1) = ampsNonresSMHeavy(1) + res

              res = ampNonresonantHeavyPP()
              ampsNonresSMHeavy(2) = ampsNonresSMHeavy(2) + res
              
              ampsSMHeavy(1:2) = ampsSMHeavy(1:2) + ampsNonresSMHeavy(1:2)

              if (real(anomc1_sc) /= 0._dp) then
                  ampsAnomHeavy(1) = ampsAnomHeavy(1) + ampsNonresSMHeavy(1)*real(anomc1_sc)
                  ampsAnomHeavy(5) = ampsAnomHeavy(5) + ampsNonresSMHeavy(2)*real(anomc1_sc)
              endif

          endif

          if (.not. disable_sm) then
              singletop2_amp_virt_dp = singletop2_amp_virt_dp + sum(real(ampSMtree(:)*conjg(ampsSMHeavy(:))))
          endif

          if (enable_nonresonant) then
              ! treeamps for renormalization
              if (real(anomc4_sc) /= 0._dp .or. aimag(anomc4_sc) /= 0._dp .or.
     &            real(anomc7_sc) /= 0._dp .or. aimag(anomc7_sc) /= 0._dp) then
                 call singletop2_nonrestree_c4(za,zb, ju,jb,jn,je,jc,jd, ampNonresTreeC4, sc)
              endif

              if (real(anomc7_sc) /= 0._dp .or. aimag(anomc7_sc) /= 0._dp) then
                 call singletop2_nonrestree_c10(za,zb, ju,jb,jn,je,jc,jd, ampNonresTreeC10, sc)
              endif

              if (real(anomc4_sc) /= 0._dp) then
                  res = ampNonresonantHeavyReC4PM()
                  ampsAnomHeavy(2) = ampsAnomHeavy(2) + res*real(anomc4_sc)
                  ampsAnomHeavy(2) = ampsAnomHeavy(2) + (ampNonresTreeC4(1)/gw**4)*real(anomc4_sc)*epinv

                  res = ampNonresonantHeavyReC4MP()
                  ampsAnomHeavy(3) = ampsAnomHeavy(3) + res*real(anomc4_sc)
                  ampsAnomHeavy(3) = ampsAnomHeavy(3) + (ampNonresTreeC4(2)/gw**4)*real(anomc4_sc)*epinv
              endif

              if (aimag(anomc4_sc) /= 0._dp) then
                  res = ampNonresonantHeavyImC4PM()
                  ampsAnomHeavy(2) = ampsAnomHeavy(2) + res*aimag(anomc4_sc)
                  ampsAnomHeavy(2) = ampsAnomHeavy(2) + (ampNonresTreeC4(1)/gw**4)*(-im)*aimag(anomc4_sc)*epinv

                  res = ampNonresonantHeavyImC4MP()
                  ampsAnomHeavy(3) = ampsAnomHeavy(3) + res*aimag(anomc4_sc)
                  ampsAnomHeavy(3) = ampsAnomHeavy(3) + (ampNonresTreeC4(2)/gw**4)*(im)*aimag(anomc4_sc)*epinv

              endif

              if (real(anomc7_sc) /= 0._dp) then
                  res = ampNonresonantHeavyReC7PM()
                  ampsAnomHeavy(2) = ampsAnomHeavy(2) + res*real(anomc7_sc)
                  ampsAnomHeavy(2) = ampsAnomHeavy(2) + epinv*real(anomc7_sc)*
     &                 + (ampNonresTreeC4(1) - 1._dp/3._dp*ampNonresTreeC10(1))/gw**4


                  res = ampNonresonantHeavyReC7MP()
                  ampsAnomHeavy(3) = ampsAnomHeavy(3) + res*real(anomc7_sc)
                  ampsAnomHeavy(3) = ampsAnomHeavy(3) + epinv*real(anomc7_sc)*
     &                 + (ampNonresTreeC4(2)- 1._dp/3._dp*ampNonresTreeC10(2))/gw**4
              endif

              if (aimag(anomc7_sc) /= 0._dp) then
                  res = ampNonresonantHeavyImC7PM()
                  ampsAnomHeavy(2) = ampsAnomHeavy(2) + res*aimag(anomc7_sc)
                  ampsAnomHeavy(2) = ampsAnomHeavy(2) + epinv*(-im)*aimag(anomc7_sc)*
     &                 + (ampNonresTreeC4(1) - 1._dp/3._dp*ampNonresTreeC10(1))/gw**4

                  res = ampNonresonantHeavyImC7MP()
                  ampsAnomHeavy(3) = ampsAnomHeavy(3) + res*aimag(anomc7_sc)
                  ampsAnomHeavy(3) = ampsAnomHeavy(3) + epinv*(im)*aimag(anomc7_sc)*
     &                 + (ampNonresTreeC4(2)- 1._dp/3._dp*ampNonresTreeC10(2))/gw**4
              endif
          endif

          if (enable_resonant) then

          ! assemble EFT amplitudes
          ! reuse standard model amplitude
          if (real(anomc1_sc) /= 0._dp) then
              ampsAnomHeavy(1) = ampsAnomHeavy(1) + 2*ampsSMResHeavy(1)*real(anomc1_sc)
          endif

          if (real(anomc8_sc) /= 0._dp) then
              res = ampHeavyMMrealc8()
              ! this has different normalization!
              ampsAnomHeavy(1) = ampsAnomHeavy(1) + (res + 
c mass renormalization
     &           (8*mtsq*mytree*(4 + 3*epinv + 3*Log(MusqHeavy/mtsq)))/(mtsq - s126)**2
     &        )*real(anomc8_sc)*prop34
          endif

          if (real(anomc3_sc) /= 0._dp) then
              res = ampHeavyMMrealc3()
              ampsAnomHeavy(1) = ampsAnomHeavy(1) + (res +
c mass renormalization
     &          (2*Sqrt(mtsq)*(mtsq + s126)*(4*mytree - struc4Step5 + struc9Step5)*
     &          (4 + 3*epinv + 3*Log(MusqHeavy/mtsq)))/(mtsq - s126)**2 +
c c3 renormalization
     &          (2*epinv*Sqrt(mtsq)*(4*mytree - struc4Step5 + struc9Step5))/(mtsq - s126)
     &        )*real(anomc3_sc)*prop16*prop34
          endif

          if (aimag(anomc3_sc) /= 0._dp) then
              res = ampHeavyMMimagc3()
              ampsAnomHeavy(1) = ampsAnomHeavy(1) + (res +
c mass renormalization
     &          (-2*im*Sqrt(mtsq)*(mtsq + s126)*(struc4Step5 + struc9Step5)*
     &            (4 + 3*epinv + 3*Log(MusqHeavy/mtsq)))/(mtsq - s126)**2 + 
c c3 renormalization
     &          (-2*epinv*im*Sqrt(mtsq)*(struc4Step5 + struc9Step5))/(mtsq - s126)
     &        )*aimag(anomc3_sc)*prop16*prop34
          endif

          if (real(anomc6_sc) /= 0._dp) then
              res = ampHeavyMMrealc6()
              ampsAnomHeavy(1) = ampsAnomHeavy(1) + (res +
c mass renormalization
     &          (8*mtsq**1.5_dp*mytree*(1 + 3*epinv + 3*Log(MusqHeavy/mtsq)))/(mtsq - s126)**2 +
c c3 renormalization
     &          (2*epinv*Sqrt(mtsq)*(4*mytree - struc4Step5 + struc9Step5))/(mtsq - s126)
     &        )*real(anomc6_sc)*prop16*prop34
          endif

          if (aimag(anomc6_sc) /= 0._dp) then
              res = ampHeavyMMimagc6()
              ampsAnomHeavy(1) = ampsAnomHeavy(1) + (res +
     &      (-2*epinv*im*Sqrt(mtsq)*(struc4Step5 + struc9Step5))/(mtsq - s126)
c c3 renormalization
     &        )*aimag(anomc6_sc)*prop16*prop34
          endif


          if (real(anomc9_sc) /= 0._dp) then
              res = ampHeavyMPPPrealc9()
              ampsAnomHeavy(4) = ampsAnomHeavy(4) + (res +
c mass renormalization
     &        (2*Sqrt(mtsq)*(mtsq + s126)*(4 + 3*epinv + 3*Log(MusqHeavy/mtsq))*za(jn,jc)*
     &    za(ju,jb)*zb(jd,je))/(mtsq - s126)**2
     &        )*real(anomc9_sc)*prop34
          endif

c PM and MP below


          if (real(anomc2_sc) /= 0._dp) then
c PM
              res = ampHeavyPMrealc2()
              ampsAnomHeavy(2) = ampsAnomHeavy(2) + (res +
c mass renormalization
     &          (Sqrt(mtsq)*(mtsq + s126)*struc59PM*
     &    (4 + 3*epinv + 3*Log(MusqHeavy/mtsq)))/(4*(mtsq - s126)**2)
     &          )*real(anomc2_sc)*prop16*prop34

c MP
              res = ampHeavyMPrealc2()
              ampsAnomHeavy(3) = ampsAnomHeavy(3) + (res +
c mass renormalization
     &        (Sqrt(mtsq)*(mtsq + s126)*struc59MP*
     &    (4 + 3*epinv + 3*Log(MusqHeavy/mtsq)))/(4*(mtsq - s126)**2)
     &          )*real(anomc2_sc)*prop16*prop34
          endif

          if (aimag(anomc2_sc) /= 0._dp) then
c PM
              res = ampHeavyPMimagc2()
              ampsAnomHeavy(2) = ampsAnomHeavy(2) + (res
c mass renormalization
     &          -(im*Sqrt(mtsq)*(mtsq + s126)*struc59PM*
     &     (4 + 3*epinv + 3*Log(MusqHeavy/mtsq)))/(4*(mtsq - s126)**2)
     &          )*aimag(anomc2_sc)*prop16*prop34

c MP
              res = ampHeavyMPimagc2()
              ampsAnomHeavy(3) = ampsAnomHeavy(3) + (res +
c mass renormalization
     &        (im*Sqrt(mtsq)*(mtsq + s126)*struc59MP*
     &    (4 + 3*epinv + 3*Log(MusqHeavy/mtsq)))/(4*(mtsq - s126)**2)
     &          )*aimag(anomc2_sc)*prop16*prop34
          endif

          if (real(anomc4_sc) /= 0._dp) then
c PM
              res = ampHeavyPMrealc4()
              ampsAnomHeavy(2) = ampsAnomHeavy(2) + (res +
c mass renormalization
     &        (-4*mtsq*(struc11PM + struc13PM + struc17PM - 
     &      (s12 + s13 + s16)*struc1PM + struc27PM + 2*struc2PM + struc37PM + 
     &      2*struc3PM + struc55PM)*(4 + 3*epinv + 3*Log(MusqHeavy/mtsq)))/
     &  (mtsq - s126)**2
c c4 renormalization
     &        -((epinv*(2*struc11PM + 2*struc13PM + 2*struc17PM - 
     &        (2*s12 + s126 + 2*(s13 + s16))*struc1PM + 
     &        2*(struc27PM + 2*struc2PM + struc37PM + 2*struc3PM + struc55PM) + 
     &        mtsq*(struc1PM - struc59PM) + s126*struc59PM))/(mtsq - s126))
     &          )*real(anomc4_sc)*prop16*prop34

c MP
              res = ampHeavyMPrealc4()
              ampsAnomHeavy(3) = ampsAnomHeavy(3) + (res +
c mass renormalization
     &        (2*mtsq*(2*struc13MP + 2*struc27MP + 4*struc3MP + 2*struc53MP + 
     &      2*struc55MP + s126*struc59MP - s16*struc59MP + 4*struc8MP)*
     &    (4 + 3*epinv + 3*Log(MusqHeavy/mtsq)))/(mtsq - s126)**2 +
c c4 renormalization
     &        (epinv*(2*struc13MP + s126*struc1MP - s16*struc59MP + 
     &      mtsq*(-struc1MP + struc59MP) + 
     &      2*(struc27MP + 2*struc3MP + struc53MP + struc55MP + 2*struc8MP)))/
     &  (mtsq - s126)
     &          )*real(anomc4_sc)*prop16*prop34

          endif

          if (aimag(anomc4_sc) /= 0._dp) then
c PM
              res = ampHeavyPMimagc4()
              ampsAnomHeavy(2) = ampsAnomHeavy(2) + (res +
c mass renormalization
     &        (4*im*mtsq*(struc11PM + struc13PM + struc17PM - 
     &      (s12 + s13 + s16)*struc1PM + struc27PM + 2*struc2PM + struc37PM + 
     &      2*struc3PM + struc55PM)*(4 + 3*epinv + 3*Log(MusqHeavy/mtsq)))/
     &  (mtsq - s126)**2
c c4 renormalization
     &      +  (epinv*im*(2*struc11PM + 2*struc13PM + 2*struc17PM - 
     &      (2*s12 + s126 + 2*(s13 + s16))*struc1PM + 
     &      2*(struc27PM + 2*struc2PM + struc37PM + 2*struc3PM + struc55PM) + 
     &      mtsq*(struc1PM - struc59PM) + s126*struc59PM))/(mtsq - s126)
     &          )*aimag(anomc4_sc)*prop16*prop34

c MP
              res = ampHeavyMPimagc4()
              ampsAnomHeavy(3) = ampsAnomHeavy(3) + (res +
c mass renormalization
     &        (2*im*mtsq*(2*struc13MP + 2*struc27MP + 4*struc3MP + 2*struc53MP + 
     &      2*struc55MP + s126*struc59MP - s16*struc59MP + 4*struc8MP)*
     &    (4 + 3*epinv + 3*Log(MusqHeavy/mtsq)))/(mtsq - s126)**2 +
c c4 renormalization
     &        (epinv*im*(2*struc13MP + s126*struc1MP - s16*struc59MP + 
     &      mtsq*(-struc1MP + struc59MP) + 
     &      2*(struc27MP + 2*struc3MP + struc53MP + struc55MP + 2*struc8MP)))/
     &  (mtsq - s126)
     &          )*aimag(anomc4_sc)*prop16*prop34
          endif

          if (real(anomc7_sc) /= 0._dp) then
c PM
              res = ampHeavyPMrealc7()
              ampsAnomHeavy(2) = ampsAnomHeavy(2) + (res
c c4 renormalization
     &         -((epinv*(2*struc11PM + 2*struc13PM + 2*struc17PM + mtsq*struc1PM - 
     &        2*s12*struc1PM - s126*struc1PM - 2*s13*struc1PM - 
     &        2*s16*struc1PM + 2*struc27PM + 4*struc2PM + 2*struc37PM + 
     &        4*struc3PM + 2*struc55PM - mtsq*struc59PM + s126*struc59PM))/
     &    (mtsq - s126))
     &          )*real(anomc7_sc)*prop16*prop34

c MP
              res = ampHeavyMPrealc7()
              ampsAnomHeavy(3) = ampsAnomHeavy(3) + (res +
c c4 renormalization
     &        (epinv*(2*struc13MP + s126*struc1MP + 2*struc27MP + 4*struc3MP + 
     &      2*struc53MP + 2*struc55MP - s16*struc59MP + 
     &      mtsq*(-struc1MP + struc59MP) + 4*struc8MP))/(mtsq - s126)
     &          )*real(anomc7_sc)*prop16*prop34
          endif

          if (aimag(anomc7_sc) /= 0._dp) then
c PM
              res = ampHeavyPMimagc7()
              ampsAnomHeavy(2) = ampsAnomHeavy(2) + (res +
c c4 renormalization
     &        (epinv*im*(2*struc11PM + 2*struc13PM + 2*struc17PM + mtsq*struc1PM - 
     &      2*s12*struc1PM - s126*struc1PM - 2*s13*struc1PM - 2*s16*struc1PM + 
     &      2*struc27PM + 4*struc2PM + 2*struc37PM + 4*struc3PM + 2*struc55PM - 
     &      mtsq*struc59PM + s126*struc59PM))/(mtsq - s126)
     &          )*aimag(anomc7_sc)*prop16*prop34

c MP
              res = ampHeavyMPimagc7()
              ampsAnomHeavy(3) = ampsAnomHeavy(3) + (res +
c c4 renormalization
     &        (epinv*im*(2*struc13MP + s126*struc1MP + 2*struc27MP + 4*struc3MP + 
     &      2*struc53MP + 2*struc55MP - s16*struc59MP + 
     &      mtsq*(-struc1MP + struc59MP) + 4*struc8MP))/(mtsq - s126)
     &          )*aimag(anomc7_sc)*prop16*prop34
          endif

          endif

c HEAVY ASSEMBLY
          if (enable_lambda2) then
              singletop2_amp_virt_dp = singletop2_amp_virt_dp
c                       ----
     &                  + real(ampSMtree(1)*conjg(ampsAnomHeavy(1)))/lambda**2/sc
     &                  + real(ampsAnomL2tree(1)*conjg(ampsSMHeavy(1)))/lambda**2/sc
c                       ++--, for non-resonant pieces
     &                  + real(ampSMtree(2)*conjg(ampsAnomHeavy(5)))/lambda**2/sc
     &                  + real(ampsAnomL2tree(5)*conjg(ampsSMHeavy(2)))/lambda**2/sc
          endif

          if (enable_lambda4) then
              singletop2_amp_virt_dp = singletop2_amp_virt_dp
     &                  + sum(real(ampsAnomL2tree*conjg(ampsAnomHeavy)))/lambda**4/sc/sc

          endif
      endif


      end

      ! Copy of dp function, except for additional initialization of spinor structures
      ! This routine is a mixture of dp and dd, where just the main matrix elements are in dd
      function singletop2_amp_virt_dd(p,za,zb,ju,jb,jn,je,jc,jd, musqLight, musqHeavy, light, heavy, sc)
        use ddmodule 
        use types
        use mod_qcdloop_c
        use eftcouple
        use singletop2_realamps_m
        use singletop2_ints_dd_m
        use singletop2_ints_nonres_dd_m
        use anomcoup_tbW
      implicit none
      real(dp):: singletop2_amp_virt_dd

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'epinv.f'
      include 'epinv2.f'

      real(dp), intent(in) :: p(mxpart,4)
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
      integer, intent(in) :: ju,jd,jn,je,jc,jb
      real(dp), intent(in) :: musqLight, musqHeavy
      logical, intent(in) :: light,heavy
      real(dp), intent(in) :: sc

      complex(dp) :: mtsq

      complex(dp) :: ampSMtree(2)
      complex(dp) :: ampsAnomL2tree(5)
      complex(dp) :: ampsAnomNonresL2tree(5)

      complex(dp) :: ampsSMLight(2), ampsSMHeavy(2), ampsSMResHeavy(2)
      complex(dp) :: ampsNonresSMLight(2), ampsNonresSMHeavy(2)

      complex(dp) :: ampsAnomLight(5), ampsAnomHeavy(5)

      complex(dp) :: ampNonresTreeC4(2)
      complex(dp) :: ampNonresTreeC10(2)

      complex(dp) :: anomc1_sc, anomc2_sc, anomc3_sc, anomc4_sc
      complex(dp) :: anomc6_sc, anomc7_sc, anomc8_sc, anomc9_sc

      complex(dp) :: res

      real(dp) :: s126, s12, s15, s16, s25, s56, s14, s36
      real(dp) :: s26, s45, s46, s24, s13, s23, s34, s35
      complex(dp) :: prop16, prop34, prop25
      integer :: i,j,k,l
      real(dp) :: s,t
      s(j,k) = real(za(j,k)*zb(k,j))
      t(j,k,l) = s(j,k) + s(j,l) + s(k,l)

      complex(dp) :: mytree, struc11step4, struc4step5, struc8step5, struc9step5

      complex(dp) :: struc10MP, struc13MP, struc14MP, struc15MP, struc17MP, struc18MP, struc19MP
      complex(dp) :: struc1MP, struc22MP, struc24MP, struc27MP, struc29MP, struc2MP
      complex(dp) :: struc31MP, struc33MP, struc34MP, struc35MP, struc36MP, struc37MP
      complex(dp) :: struc38MP, struc39MP, struc3MP, struc42MP, struc43MP, struc44MP
      complex(dp) :: struc45MP, struc47MP, struc50MP, struc51MP, struc53MP, struc54MP
      complex(dp) :: struc55MP, struc56MP, struc57MP, struc58MP, struc59MP, struc5MP
      complex(dp) :: struc6MP, struc7MP, struc8MP, struc9MP, struc41MP

      complex(dp) :: struc10PM, struc11PM, struc12PM, struc13PM, struc14PM, struc15PM, struc16PM
      complex(dp) :: struc17PM, struc18PM, struc19PM, struc1PM, struc21PM, struc22PM
      complex(dp) :: struc23PM, struc24PM, struc25PM, struc26PM, struc27PM, struc28PM
      complex(dp) :: struc29PM, struc2PM, struc31PM, struc32PM, struc33PM, struc34PM
      complex(dp) :: struc35PM, struc36PM, struc37PM, struc38PM, struc3PM, struc42PM
      complex(dp) :: struc43PM, struc44PM, struc45PM, struc47PM, struc48PM, struc49PM
      complex(dp) :: struc50PM, struc51PM, struc52PM, struc55PM, struc56PM, struc57PM
      complex(dp) :: struc59PM, struc5PM, struc6PM, struc7PM, struc8PM, struc9PM

      anomc1_sc = anomc1 * sc
      anomc2_sc = anomc2 * sc
      anomc3_sc = anomc3 * sqrt(sc)
      anomc4_sc = anomc4 * sqrt(sc)
      anomc6_sc = anomc6 * sqrt(sc)
      anomc7_sc = anomc7 * sqrt(sc)
      anomc8_sc = anomc8
      anomc9_sc = anomc9

      mytree = za(jc,jn)*zb(ju,jb)*(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd))
      struc9Step5 = 2._dp*za(jb,jd)*za(jn,jc)*zb(jb,ju)*zb(je,jb)
      struc4Step5 = -2._dp*za(jn,jc)*za(jc,jd)*zb(jb,ju)*zb(jc,je)
      struc8Step5 = 2._dp*za(jn,jc)*za(jn,jd)*zb(jn,ju)*zb(je,jb)
      struc11step4 = -(za(jc,jd)*za(jn,jc)*za(jn,jd)*zb(jc,je)*zb(jd,jb)*zb(jn,ju))


      struc10PM = -(za(je,jd)*za(jn,jc)*zb(jc,jb)*zb(jc,je)*zb(je,ju))
      struc11PM = -2*za(ju,jd)*za(ju,jn)*zb(jb,ju)*zb(jc,ju)*zb(je,ju)
      struc12PM = -2*za(jb,jn)*za(ju,jd)*zb(jb,ju)*zb(jc,ju)*zb(je,jb)
      struc13PM = -2*za(jb,jd)*za(ju,jn)*zb(jb,ju)*zb(jc,ju)*zb(je,jb)
      struc14PM = -2*za(jn,jd)*za(ju,jn)*zb(jc,ju)*zb(je,jb)*zb(jn,ju)
      struc15PM = -2*za(je,jd)*za(ju,jn)*zb(jc,ju)*zb(je,jb)*zb(je,ju)
      struc16PM = -2*za(jn,jc)*za(ju,jd)*zb(jb,ju)*zb(jc,je)*zb(jc,ju)
      struc17PM = -2*za(jn,jd)*za(ju,jn)*zb(jc,ju)*zb(je,ju)*zb(jn,jb)
      struc18PM = -(za(jb,jd)*za(ju,jn)**2*zb(jb,ju)*zb(jc,ju)*zb(je,ju)*zb(jn,jb))
      struc19PM = -(za(jb,jd)*za(jb,jn)*za(ju,jn)*zb(jb,ju)*zb(jc,ju)*zb(je,jb)*zb(jn,jb))
      struc1PM = -2*za(jn,jd)*zb(jc,jb)*zb(je,ju)
      struc21PM = -(za(jb,jn)*za(jn,jd)*za(ju,jn)*zb(jc,ju)*zb(je,jb)*zb(jn,jb)*zb(jn,ju))
      struc22PM = -(za(je,jd)*za(ju,jn)**2*zb(jc,ju)*zb(je,ju)**2*zb(jn,jb))
      struc23PM = -(za(jb,jn)*za(je,jd)*za(ju,jn)*zb(jc,ju)*zb(je,jb)*zb(je,ju)*zb(jn,jb))
      struc24PM = -(za(jb,jd)*za(jn,jc)*za(ju,jn)*zb(jb,ju)*zb(jc,je)*zb(jc,ju)*zb(jn,jb))
      struc25PM = -(za(jn,jc)*za(jn,jd)*za(ju,jn)*zb(jc,je)*zb(jc,ju)*zb(jn,jb)*zb(jn,ju))
      struc26PM = -(za(je,jd)*za(jn,jc)*za(ju,jn)*zb(jc,je)*zb(jc,ju)*zb(je,ju)*zb(jn,jb))
      struc27PM = -2*za(jn,jd)*za(ju,jd)*zb(jc,ju)*zb(jd,jb)*zb(je,ju)
      struc28PM = -(za(jb,jd)*za(ju,jd)*za(ju,jn)*zb(jb,ju)*zb(jc,ju)*zb(jd,jb)*zb(je,ju))
      struc29PM = -(za(jb,jd)*za(jb,jn)*za(ju,jd)*zb(jb,ju)*zb(jc,ju)*zb(jd,jb)*zb(je,jb))
      struc2PM = -(za(jb,jd)*za(ju,jn)*zb(jb,ju)*zb(jc,jb)*zb(je,ju))
      struc31PM = -(za(jb,jn)*za(jn,jd)*za(ju,jd)*zb(jc,ju)*zb(jd,jb)*zb(je,jb)*zb(jn,ju))
      struc32PM = -(za(je,jd)*za(ju,jd)*za(ju,jn)*zb(jc,ju)*zb(jd,jb)*zb(je,ju)**2)
      struc33PM = -(za(jb,jn)*za(je,jd)*za(ju,jd)*zb(jc,ju)*zb(jd,jb)*zb(je,jb)*zb(je,ju))
      struc34PM = -(za(jb,jd)*za(jn,jc)*za(ju,jd)*zb(jb,ju)*zb(jc,je)*zb(jc,ju)*zb(jd,jb))
      struc35PM = -(za(jn,jc)*za(jn,jd)*za(ju,jd)*zb(jc,je)*zb(jc,ju)*zb(jd,jb)*zb(jn,ju))
      struc36PM = -(za(je,jd)*za(jn,jc)*za(ju,jd)*zb(jc,je)*zb(jc,ju)*zb(jd,jb)*zb(je,ju))
      struc37PM = -2*za(jn,jd)*za(ju,jn)*zb(jb,ju)*zb(jc,jn)*zb(je,ju)
      struc38PM = -2*za(jb,jn)*za(jn,jd)*zb(jb,ju)*zb(jc,jn)*zb(je,jb)
      struc3PM = -(za(jb,jd)*za(jb,jn)*zb(jb,ju)*zb(jc,jb)*zb(je,jb))
      struc42PM = -2*za(jn,jc)*za(jn,jd)*zb(jb,ju)*zb(jc,je)*zb(jc,jn)
      struc43PM = -2*za(jn,jd)**2*zb(jc,jn)*zb(jd,jb)*zb(je,ju)
      struc44PM = -(za(jb,jd)*za(jn,jd)*za(ju,jn)*zb(jb,ju)*zb(jc,jn)*zb(jd,jb)*zb(je,ju))
      struc45PM = -(za(jb,jd)*za(jb,jn)*za(jn,jd)*zb(jb,ju)*zb(jc,jn)*zb(jd,jb)*zb(je,jb))
      struc47PM = -(za(jb,jn)*za(jn,jd)**2*zb(jc,jn)*zb(jd,jb)*zb(je,jb)*zb(jn,ju))
      struc48PM = -(za(je,jd)*za(jn,jd)*za(ju,jn)*zb(jc,jn)*zb(jd,jb)*zb(je,ju)**2)
      struc49PM = -(za(jb,jn)*za(je,jd)*za(jn,jd)*zb(jc,jn)*zb(jd,jb)*zb(je,jb)*zb(je,ju))
      struc50PM = -(za(jb,jd)*za(jn,jc)*za(jn,jd)*zb(jb,ju)*zb(jc,je)*zb(jc,jn)*zb(jd,jb))
      struc51PM = -(za(jn,jc)*za(jn,jd)**2*zb(jc,je)*zb(jc,jn)*zb(jd,jb)*zb(jn,ju))
      struc52PM = -(za(je,jd)*za(jn,jc)*za(jn,jd)*zb(jc,je)*zb(jc,jn)*zb(jd,jb)*zb(je,ju))
      struc55PM = -2*za(jb,jd)*za(jn,jd)*zb(jb,ju)*zb(jd,jc)*zb(je,jb)
      struc56PM = -2*za(jn,jd)**2*zb(jd,jc)*zb(je,jb)*zb(jn,ju)
      struc57PM = -2*za(je,jd)*za(jn,jd)*zb(jd,jc)*zb(je,jb)*zb(je,ju)
      struc59PM = -4*za(jn,jd)*zb(jb,ju)*zb(jc,je)
      struc5PM = -(za(jb,jn)*za(jn,jd)*zb(jc,jb)*zb(je,jb)*zb(jn,ju))
      struc6PM = -(za(je,jd)*za(ju,jn)*zb(jc,jb)*zb(je,ju)**2)
      struc7PM = -(za(jb,jn)*za(je,jd)*zb(jc,jb)*zb(je,jb)*zb(je,ju))
      struc8PM = -(za(jb,jd)*za(jn,jc)*zb(jb,ju)*zb(jc,jb)*zb(jc,je))
      struc9PM = -(za(jn,jc)*za(jn,jd)*zb(jc,jb)*zb(jc,je)*zb(jn,ju))


      struc10MP = za(jb,jc)*za(je,jd)*za(jn,jc)*zb(jc,je)*zb(je,ju)
      struc13MP = 2*za(jb,jd)*za(jb,jn)*za(ju,jc)*zb(jb,ju)*zb(je,ju)
      struc14MP = 2*za(jb,jn)*za(jn,jd)*za(ju,jc)*zb(je,ju)*zb(jn,ju)
      struc15MP = 2*za(jb,jn)*za(je,jd)*za(ju,jc)*zb(je,ju)**2
      struc17MP = 2*za(jb,jn)*za(jn,jd)*za(ju,jc)*zb(je,ju)*zb(jn,ju)
      struc18MP = za(jb,jd)*za(jb,jn)*za(ju,jc)*za(ju,jn)*zb(jb,ju)*zb(je,ju)*zb(jn,ju)
      struc19MP = za(jb,jd)*za(jb,jn)**2*za(ju,jc)*zb(jb,ju)*zb(je,jb)*zb(jn,ju)
      struc1MP = 2*za(jb,jc)*za(jn,jd)*zb(je,ju)
      struc22MP = za(jb,jn)*za(je,jd)*za(ju,jc)*za(ju,jn)*zb(je,ju)**2*zb(jn,ju)
      struc24MP = za(jb,jd)*za(jb,jn)*za(jn,jc)*za(ju,jc)*zb(jb,ju)*zb(jc,je)*zb(jn,ju)
      struc27MP = 2*za(jb,jd)*za(jn,jd)*za(ju,jc)*zb(jd,ju)*zb(je,ju)
      struc29MP = za(jb,jd)**2*za(jb,jn)*za(ju,jc)*zb(jb,ju)*zb(jd,ju)*zb(je,jb)
      struc2MP = za(jb,jc)*za(jb,jd)*za(ju,jn)*zb(jb,ju)*zb(je,ju)
      struc31MP = za(jb,jd)*za(jb,jn)*za(jn,jd)*za(ju,jc)*zb(jd,ju)*zb(je,jb)*zb(jn,ju)
      struc33MP = za(jb,jd)*za(jb,jn)*za(je,jd)*za(ju,jc)*zb(jd,ju)*zb(je,jb)*zb(je,ju)
      struc34MP = za(jb,jd)**2*za(jn,jc)*za(ju,jc)*zb(jb,ju)*zb(jc,je)*zb(jd,ju)
      struc35MP = za(jb,jd)*za(jn,jc)*za(jn,jd)*za(ju,jc)*zb(jc,je)*zb(jd,ju)*zb(jn,ju)
      struc36MP = za(jb,jd)*za(je,jd)*za(jn,jc)*za(ju,jc)*zb(jc,je)*zb(jd,ju)*zb(je,ju)
      struc37MP = 2*za(jb,jd)*za(jn,jc)*za(ju,jn)*zb(je,ju)*zb(jn,ju)
      struc38MP = 2*za(jb,jd)*za(jb,jn)*za(jn,jc)*zb(je,jb)*zb(jn,ju)
      struc39MP = 2*za(jb,jd)*za(jb,jn)*za(jn,jc)*zb(jb,ju)*zb(je,jn)
      struc3MP = za(jb,jc)*za(jb,jd)*za(jb,jn)*zb(jb,ju)*zb(je,jb)
      struc42MP = 2*za(jb,jd)*za(jn,jc)**2*zb(jc,je)*zb(jn,ju)
      struc43MP = 2*za(jb,jd)*za(jn,jc)*za(jn,jd)*zb(jd,jn)*zb(je,ju)
      struc44MP = za(jb,jd)**2*za(jn,jc)*za(ju,jn)*zb(jb,ju)*zb(jd,jn)*zb(je,ju)
      struc45MP = za(jb,jd)**2*za(jb,jn)*za(jn,jc)*zb(jb,ju)*zb(jd,jn)*zb(je,jb)
      struc47MP = za(jb,jd)*za(jb,jn)*za(jn,jc)*za(jn,jd)*zb(jd,jn)*zb(je,jb)*zb(jn,ju)
      struc50MP = za(jb,jd)**2*za(jn,jc)**2*zb(jb,ju)*zb(jc,je)*zb(jd,jn)
      struc51MP = za(jb,jd)*za(jn,jc)**2*za(jn,jd)*zb(jc,je)*zb(jd,jn)*zb(jn,ju)
      struc53MP = -2*za(jb,jd)*za(jc,jd)*za(ju,jn)*zb(jd,ju)*zb(je,ju)
      struc54MP = -2*za(jb,jd)*za(jb,jn)*za(jc,jd)*zb(jd,ju)*zb(je,jb)
      struc55MP = 2*za(jb,jd)*za(jb,jn)*za(jc,jd)*zb(jb,ju)*zb(jd,je)
      struc56MP = 2*za(jb,jn)*za(jc,jd)*za(jn,jd)*zb(jd,je)*zb(jn,ju)
      struc57MP = 2*za(jb,jn)*za(jc,jd)*za(je,jd)*zb(jd,je)*zb(je,ju)
      struc58MP = -2*za(jb,jd)*za(jc,jd)*za(jn,jc)*zb(jc,je)*zb(jd,ju)
      struc59MP = 4*za(jb,jd)*za(jn,jc)*zb(je,ju)
      struc5MP = za(jb,jc)*za(jb,jn)*za(jn,jd)*zb(je,jb)*zb(jn,ju)
      struc6MP = za(jb,jc)*za(je,jd)*za(ju,jn)*zb(je,ju)**2
      struc7MP = za(jb,jc)*za(jb,jn)*za(je,jd)*zb(je,jb)*zb(je,ju)
      struc8MP = za(jb,jc)*za(jb,jd)*za(jn,jc)*zb(jb,ju)*zb(jc,je)
      struc9MP = za(jb,jc)*za(jn,jc)*za(jn,jd)*zb(jc,je)*zb(jn,ju)
      struc41MP = 2*za(jb,jn)*za(je,jd)*za(jn,jc)*zb(je,jn)*zb(je,ju)


      ampsAnomLight = 0._dp
      ampsAnomHeavy = 0._dp
      ampsSMLight = 0._dp
      ampsSMHeavy = 0._dp
      ampsSMResHeavy = 0._dp
      ampSMtree = 0._dp
      ampsAnomL2tree = 0._dp

      ampsNonresSMLight = 0._dp
      ampsNonresSMHeavy = 0._dp

      ampNonresTreeC4 = 0._dp
      ampNonresTreeC10 = 0._dp

      mtsq = (mt**2 - im*mt*twidth)*sc

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

      prop16 = 1._dp / (s(ju,jd) - wmass**2 * sc)
      prop34 = 1._dp / (s(jn,je) - wmass**2*sc + im*wmass*wwidth*sc)
      prop25 = 1._dp / (s(jb,jc) - zmass**2*sc)

      ! resonant
      call initamp(musqHeavy, mtsq, za,zb,ju,jb,jn,je,jc,jd, epinv, epinv*epinv2, sc)
      ! non resonant 
      call initamp_nonres(musqHeavy, musqLight, za,zb,ju,jb,jn,je,jc,jd, epinv, epinv*epinv2, sc)

      call singletop2_hamp_tree(za,zb, ju,jb,jn,je,jc,jd, mtsq, ampSMtree, ampsAnomL2tree, ampsAnomNonresL2tree, sc)

      if (light) then

          if (enable_resonant) then
              ampsSMLight(1) = (cmplx(mytree)*prop16*prop34*(8 + epinv*(3 + 2*epinv2)
     &          + (3 + 2*epinv)*Log(-(MusqLight/s16)) + 
     &          Log(-(MusqLight/s16))**2))/(-mtsq + s126)
          endif

          if (enable_nonresonant) then
              ! MM contributions,
              res = ampNonresonantLightFullMM()
            ampsNonresSMLight(1) = ampsNonresSMLight(1) + res

              ! PP contributions
              res = ampNonresonantLightFullPP()
            ampsNonresSMLight(2) = ampsNonresSMLight(2) + res

            ampsSMLight(1:2) = ampsSMLight(1:2) + ampsNonresSMLight(1:2)
          endif

    
          if (enable_lambda2) then
              if (enable_resonant) then
                  ampsAnomLight(1) = ampsAnomLight(1) + (2*prop34*(8 + epinv*(3 + 2*epinv2) + (3 + 2*epinv)*Log(-(MusqLight/s16)) + 
     -            Log(-(MusqLight/s16))**2)*(im*sqrt(mtsq)*prop16*(struc4Step5 + struc9Step5)*imag(anomc3_sc) - 
     -            mytree*prop16*real(anomc1_sc) - 4*sqrt(mtsq)*mytree*prop16*real(anomc3_sc) + 
     -            sqrt(mtsq)*prop16*struc4Step5*real(anomc3_sc) - sqrt(mtsq)*prop16*struc9Step5*real(anomc3_sc) - 
     -            2*mytree*real(anomc8_sc)))/(mtsq - s126)
              endif

              if (enable_nonresonant) then
                ampsAnomLight(1) = ampsAnomLight(1) + ampsNonresSMLight(1)*real(anomc1_sc)
                ampsAnomLight(5) = ampsAnomLight(5) + ampsNonresSMLight(2)*real(anomc1_sc)
              endif

          endif


          if (enable_lambda4) then
            if (enable_resonant) then
             ! PM
             ampsAnomLight(2) = (prop16*prop34*(8 + epinv*(3 + 2*epinv2) + 
     -      Log(-(MusqLight/s16))*(3 + 2*epinv + Log(-(MusqLight/s16))))*
     -    (im*Sqrt(mtsq)*struc59PM*imag(anomc2_sc) + 
     -      4*im*(2*struc14PM + 2*struc16PM - 2*struc17PM + 
     -         (s126 + 2*(s13 + s34 + s35))*struc1PM - s34*struc59PM + 
     -         mtsq*(-struc1PM + struc59PM) + 
     -         2*(struc43PM + struc56PM + 2*(struc5PM + struc8PM)))*imag(anomc4_sc) - 
     -      Sqrt(mtsq)*struc59PM*real(anomc2_sc) - 
     -      4*(2*struc14PM + 2*struc16PM - 2*struc17PM + 
     -         (s126 + 2*(s13 + s34 + s35))*struc1PM - s34*struc59PM + 
     -         mtsq*(-struc1PM + struc59PM) + 
     -         2*(struc43PM + struc56PM + 2*(struc5PM + struc8PM)))*real(anomc4_sc)))/
     -  (4*(mtsq - s126))

             ! MP
             ampsAnomLight(3) = -(prop16*prop34*(8 + epinv*(3 + 2*epinv2) + 
     -       Log(-(MusqLight/s16))*(3 + 2*epinv + Log(-(MusqLight/s16))))*
     -     (im*Sqrt(mtsq)*struc59MP*imag(anomc2_sc) - 
     -       2*im*(-4*struc17MP + 2*(mtsq - s126 + 2*(s13 + s35))*struc1MP + 
     -          (-2*mtsq + s126 + s34 - s35)*struc59MP + 
     -          2*(2*struc39MP + 2*struc41MP + struc42MP + 2*struc43MP + struc58MP - 
     -             2*struc8MP))*imag(anomc4_sc) + Sqrt(mtsq)*struc59MP*real(anomc2_sc) - 
     -       2*(-4*struc17MP + 2*(mtsq - s126 + 2*(s13 + s35))*struc1MP + 
     -          (-2*mtsq + s126 + s34 - s35)*struc59MP + 
     -          2*(2*struc39MP + 2*struc41MP + struc42MP + 2*struc43MP + struc58MP - 
     -             2*struc8MP))*real(anomc4_sc)))/(4.*(mtsq - s126))

              ampsAnomLight(4) = (-2*Sqrt(mtsq)*prop34*(8 + epinv*(3 + 2*epinv2) + 
     -      (3 + 2*epinv)*Log(-(MusqLight/s16)) + Log(-(MusqLight/s16))**2)*real(anomc9_sc)*
     -    za(jn,jc)*za(ju,jb)*zb(jd,je))/(mtsq - s126)
          endif

          endif

          if (enable_nonresonant) then

              if( real(anomc4_sc) /= 0._dp ) then
                  res = ampNonresonantLightFullReC4PM()
                  ampsAnomLight(2) = ampsAnomLight(2) + res*real(anomc4_sc)

                  res = ampNonresonantLightFullReC4MP()
                  ampsAnomLight(3) = ampsAnomLight(3) + res*real(anomc4_sc)
              endif

              if( aimag(anomc4_sc) /= 0._dp ) then
                  res = ampNonresonantLightFullImC4PM()
                  ampsAnomLight(2) = ampsAnomLight(2) + res*aimag(anomc4_sc)

                  res = ampNonresonantLightFullImC4MP()
                  ampsAnomLight(3) = ampsAnomLight(3) + res*aimag(anomc4_sc)
              endif

          endif ! enable_nonresonant
      endif

      ! in the following: 
      ! no factor of two, since we calculate in 1/(4pi), but mcfm convention is 1/(2pi)

      singletop2_amp_virt_dd = 0._dp

c assemble light amplitudes
      if (light) then
          if (.not. disable_sm) then
              singletop2_amp_virt_dd = singletop2_amp_virt_dd + sum(real(ampSMtree(:)*conjg(ampsSMLight(:))))
          endif

          if (enable_lambda2) then
              singletop2_amp_virt_dd = singletop2_amp_virt_dd 
c                       ---- contribution
     &                  + real(ampSMtree(1)*conjg(ampsAnomLight(1)))/lambda**2/sc
     &                  + real(ampsAnomL2tree(1)*conjg(ampsSMLight(1)))/lambda**2/sc
c                       ++-- contribution, for non-resonant pieces
c                           which is 5 in anom array and 2 in SM array
     &                  + real(ampSMtree(2)*conjg(ampsAnomLight(5)))/lambda**2/sc
     &                  + real(ampsAnomL2tree(5)*conjg(ampsSMLight(2)))/lambda**2/sc
          endif

          if (enable_lambda4) then
              singletop2_amp_virt_dd = singletop2_amp_virt_dd
     &                  + sum(real(ampsAnomL2tree*conjg(ampsAnomLight)))/lambda**4/sc/sc
          endif
      endif

      if (heavy) then

          if (enable_resonant) then
              res = ampHeavyMM()
              ampsSMResHeavy(1) = (res 
c mass renormalization
     &          + (2*mtsq*mytree*(4 + 3*epinv + 3*Log(MusqHeavy/mtsq)))/(mtsq - s126)**2
     &              )*prop16*prop34
          endif

          ampsSMHeavy(1) = ampsSMHeavy(1) + ampsSMResHeavy(1)

          ! non-resonant contributions
          if (enable_nonresonant) then
              res = ampNonresonantHeavyMM()
              ampsNonresSMHeavy(1) = ampsNonresSMHeavy(1) + res

              res = ampNonresonantHeavyPP()
              ampsNonresSMHeavy(2) = ampsNonresSMHeavy(2) + res
              
              ampsSMHeavy(1:2) = ampsSMHeavy(1:2) + ampsNonresSMHeavy(1:2)

              if (real(anomc1_sc) /= 0._dp) then
                  ampsAnomHeavy(1) = ampsAnomHeavy(1) + ampsNonresSMHeavy(1)*real(anomc1_sc)
                  ampsAnomHeavy(5) = ampsAnomHeavy(5) + ampsNonresSMHeavy(2)*real(anomc1_sc)
              endif

          endif

          if (.not. disable_sm) then
              singletop2_amp_virt_dd = singletop2_amp_virt_dd + sum(real(ampSMtree(:)*conjg(ampsSMHeavy(:))))
          endif

          if (enable_nonresonant) then
              ! treeamps for renormalization
              if (real(anomc4_sc) /= 0._dp .or. aimag(anomc4_sc) /= 0._dp .or.
     &            real(anomc7_sc) /= 0._dp .or. aimag(anomc7_sc) /= 0._dp) then
                 call singletop2_nonrestree_c4(za,zb, ju,jb,jn,je,jc,jd, ampNonresTreeC4, sc)
              endif

              if (real(anomc7_sc) /= 0._dp .or. aimag(anomc7_sc) /= 0._dp) then
                 call singletop2_nonrestree_c10(za,zb, ju,jb,jn,je,jc,jd, ampNonresTreeC10, sc)
              endif

              if (real(anomc4_sc) /= 0._dp) then
                  res = ampNonresonantHeavyReC4PM()
                  ampsAnomHeavy(2) = ampsAnomHeavy(2) + res*real(anomc4_sc)
                  ampsAnomHeavy(2) = ampsAnomHeavy(2) + (ampNonresTreeC4(1)/gw**4)*real(anomc4_sc)*epinv

                  res = ampNonresonantHeavyReC4MP()
                  ampsAnomHeavy(3) = ampsAnomHeavy(3) + res*real(anomc4_sc)
                  ampsAnomHeavy(3) = ampsAnomHeavy(3) + (ampNonresTreeC4(2)/gw**4)*real(anomc4_sc)*epinv
              endif

              if (aimag(anomc4_sc) /= 0._dp) then
                  res = ampNonresonantHeavyImC4PM()
                  ampsAnomHeavy(2) = ampsAnomHeavy(2) + res*aimag(anomc4_sc)
                  ampsAnomHeavy(2) = ampsAnomHeavy(2) + (ampNonresTreeC4(1)/gw**4)*(-im)*aimag(anomc4_sc)*epinv

                  res = ampNonresonantHeavyImC4MP()
                  ampsAnomHeavy(3) = ampsAnomHeavy(3) + res*aimag(anomc4_sc)
                  ampsAnomHeavy(3) = ampsAnomHeavy(3) + (ampNonresTreeC4(2)/gw**4)*(im)*aimag(anomc4_sc)*epinv

              endif

              if (real(anomc7_sc) /= 0._dp) then
                  res = ampNonresonantHeavyReC7PM()
                  ampsAnomHeavy(2) = ampsAnomHeavy(2) + res*real(anomc7_sc)
                  ampsAnomHeavy(2) = ampsAnomHeavy(2) + epinv*real(anomc7_sc)*
     &                 + (ampNonresTreeC4(1) - 1._dp/3._dp*ampNonresTreeC10(1))/gw**4


                  res = ampNonresonantHeavyReC7MP()
                  ampsAnomHeavy(3) = ampsAnomHeavy(3) + res*real(anomc7_sc)
                  ampsAnomHeavy(3) = ampsAnomHeavy(3) + epinv*real(anomc7_sc)*
     &                 + (ampNonresTreeC4(2)- 1._dp/3._dp*ampNonresTreeC10(2))/gw**4
              endif

              if (aimag(anomc7_sc) /= 0._dp) then
                  res = ampNonresonantHeavyImC7PM()
                  ampsAnomHeavy(2) = ampsAnomHeavy(2) + res*aimag(anomc7_sc)
                  ampsAnomHeavy(2) = ampsAnomHeavy(2) + epinv*(-im)*aimag(anomc7_sc)*
     &                 + (ampNonresTreeC4(1) - 1._dp/3._dp*ampNonresTreeC10(1))/gw**4

                  res = ampNonresonantHeavyImC7MP()
                  ampsAnomHeavy(3) = ampsAnomHeavy(3) + res*aimag(anomc7_sc)
                  ampsAnomHeavy(3) = ampsAnomHeavy(3) + epinv*(im)*aimag(anomc7_sc)*
     &                 + (ampNonresTreeC4(2)- 1._dp/3._dp*ampNonresTreeC10(2))/gw**4
              endif
          endif

          if (enable_resonant) then

          ! assemble EFT amplitudes
          ! reuse standard model amplitude
          if (real(anomc1_sc) /= 0._dp) then
              ampsAnomHeavy(1) = ampsAnomHeavy(1) + 2*ampsSMResHeavy(1)*real(anomc1_sc)
          endif

          if (real(anomc8_sc) /= 0._dp) then
              res = ampHeavyMMrealc8()
              ! this has different normalization!
              ampsAnomHeavy(1) = ampsAnomHeavy(1) + (res + 
c mass renormalization
     &           (8*mtsq*mytree*(4 + 3*epinv + 3*Log(MusqHeavy/mtsq)))/(mtsq - s126)**2
     &        )*real(anomc8_sc)*prop34
          endif

          if (real(anomc3_sc) /= 0._dp) then
              res = ampHeavyMMrealc3()
              ampsAnomHeavy(1) = ampsAnomHeavy(1) + (res +
c mass renormalization
     &          (2*Sqrt(mtsq)*(mtsq + s126)*(4*mytree - struc4Step5 + struc9Step5)*
     &          (4 + 3*epinv + 3*Log(MusqHeavy/mtsq)))/(mtsq - s126)**2 +
c c3 renormalization
     &          (2*epinv*Sqrt(mtsq)*(4*mytree - struc4Step5 + struc9Step5))/(mtsq - s126)
     &        )*real(anomc3_sc)*prop16*prop34
          endif

          if (aimag(anomc3_sc) /= 0._dp) then
              res = ampHeavyMMimagc3()
              ampsAnomHeavy(1) = ampsAnomHeavy(1) + (res +
c mass renormalization
     &          (-2*im*Sqrt(mtsq)*(mtsq + s126)*(struc4Step5 + struc9Step5)*
     &            (4 + 3*epinv + 3*Log(MusqHeavy/mtsq)))/(mtsq - s126)**2 + 
c c3 renormalization
     &          (-2*epinv*im*Sqrt(mtsq)*(struc4Step5 + struc9Step5))/(mtsq - s126)
     &        )*aimag(anomc3_sc)*prop16*prop34
          endif

          if (real(anomc6_sc) /= 0._dp) then
              res = ampHeavyMMrealc6()
              ampsAnomHeavy(1) = ampsAnomHeavy(1) + (res +
c mass renormalization
     &          (8*mtsq**1.5_dp*mytree*(1 + 3*epinv + 3*Log(MusqHeavy/mtsq)))/(mtsq - s126)**2 +
c c3 renormalization
     &          (2*epinv*Sqrt(mtsq)*(4*mytree - struc4Step5 + struc9Step5))/(mtsq - s126)
     &        )*real(anomc6_sc)*prop16*prop34
          endif

          if (aimag(anomc6_sc) /= 0._dp) then
              res = ampHeavyMMimagc6()
              ampsAnomHeavy(1) = ampsAnomHeavy(1) + (res +
     &      (-2*epinv*im*Sqrt(mtsq)*(struc4Step5 + struc9Step5))/(mtsq - s126)
c c3 renormalization
     &        )*aimag(anomc6_sc)*prop16*prop34
          endif


          if (real(anomc9_sc) /= 0._dp) then
              res = ampHeavyMPPPrealc9()
              ampsAnomHeavy(4) = ampsAnomHeavy(4) + (res +
c mass renormalization
     &        (2*Sqrt(mtsq)*(mtsq + s126)*(4 + 3*epinv + 3*Log(MusqHeavy/mtsq))*za(jn,jc)*
     &    za(ju,jb)*zb(jd,je))/(mtsq - s126)**2
     &        )*real(anomc9_sc)*prop34
          endif

c PM and MP below


          if (real(anomc2_sc) /= 0._dp) then
c PM
              res = ampHeavyPMrealc2()
              ampsAnomHeavy(2) = ampsAnomHeavy(2) + (res +
c mass renormalization
     &          (Sqrt(mtsq)*(mtsq + s126)*struc59PM*
     &    (4 + 3*epinv + 3*Log(MusqHeavy/mtsq)))/(4*(mtsq - s126)**2)
     &          )*real(anomc2_sc)*prop16*prop34

c MP
              res = ampHeavyMPrealc2()
              ampsAnomHeavy(3) = ampsAnomHeavy(3) + (res +
c mass renormalization
     &        (Sqrt(mtsq)*(mtsq + s126)*struc59MP*
     &    (4 + 3*epinv + 3*Log(MusqHeavy/mtsq)))/(4*(mtsq - s126)**2)
     &          )*real(anomc2_sc)*prop16*prop34
          endif

          if (aimag(anomc2_sc) /= 0._dp) then
c PM
              res = ampHeavyPMimagc2()
              ampsAnomHeavy(2) = ampsAnomHeavy(2) + (res
c mass renormalization
     &          -(im*Sqrt(mtsq)*(mtsq + s126)*struc59PM*
     &     (4 + 3*epinv + 3*Log(MusqHeavy/mtsq)))/(4*(mtsq - s126)**2)
     &          )*aimag(anomc2_sc)*prop16*prop34

c MP
              res = ampHeavyMPimagc2()
              ampsAnomHeavy(3) = ampsAnomHeavy(3) + (res +
c mass renormalization
     &        (im*Sqrt(mtsq)*(mtsq + s126)*struc59MP*
     &    (4 + 3*epinv + 3*Log(MusqHeavy/mtsq)))/(4*(mtsq - s126)**2)
     &          )*aimag(anomc2_sc)*prop16*prop34
          endif

          if (real(anomc4_sc) /= 0._dp) then
c PM
              res = ampHeavyPMrealc4()
              ampsAnomHeavy(2) = ampsAnomHeavy(2) + (res +
c mass renormalization
     &        (-4*mtsq*(struc11PM + struc13PM + struc17PM - 
     &      (s12 + s13 + s16)*struc1PM + struc27PM + 2*struc2PM + struc37PM + 
     &      2*struc3PM + struc55PM)*(4 + 3*epinv + 3*Log(MusqHeavy/mtsq)))/
     &  (mtsq - s126)**2
c c4 renormalization
     &        -((epinv*(2*struc11PM + 2*struc13PM + 2*struc17PM - 
     &        (2*s12 + s126 + 2*(s13 + s16))*struc1PM + 
     &        2*(struc27PM + 2*struc2PM + struc37PM + 2*struc3PM + struc55PM) + 
     &        mtsq*(struc1PM - struc59PM) + s126*struc59PM))/(mtsq - s126))
     &          )*real(anomc4_sc)*prop16*prop34

c MP
              res = ampHeavyMPrealc4()
              ampsAnomHeavy(3) = ampsAnomHeavy(3) + (res +
c mass renormalization
     &        (2*mtsq*(2*struc13MP + 2*struc27MP + 4*struc3MP + 2*struc53MP + 
     &      2*struc55MP + s126*struc59MP - s16*struc59MP + 4*struc8MP)*
     &    (4 + 3*epinv + 3*Log(MusqHeavy/mtsq)))/(mtsq - s126)**2 +
c c4 renormalization
     &        (epinv*(2*struc13MP + s126*struc1MP - s16*struc59MP + 
     &      mtsq*(-struc1MP + struc59MP) + 
     &      2*(struc27MP + 2*struc3MP + struc53MP + struc55MP + 2*struc8MP)))/
     &  (mtsq - s126)
     &          )*real(anomc4_sc)*prop16*prop34

          endif

          if (aimag(anomc4_sc) /= 0._dp) then
c PM
              res = ampHeavyPMimagc4()
              ampsAnomHeavy(2) = ampsAnomHeavy(2) + (res +
c mass renormalization
     &        (4*im*mtsq*(struc11PM + struc13PM + struc17PM - 
     &      (s12 + s13 + s16)*struc1PM + struc27PM + 2*struc2PM + struc37PM + 
     &      2*struc3PM + struc55PM)*(4 + 3*epinv + 3*Log(MusqHeavy/mtsq)))/
     &  (mtsq - s126)**2
c c4 renormalization
     &      +  (epinv*im*(2*struc11PM + 2*struc13PM + 2*struc17PM - 
     &      (2*s12 + s126 + 2*(s13 + s16))*struc1PM + 
     &      2*(struc27PM + 2*struc2PM + struc37PM + 2*struc3PM + struc55PM) + 
     &      mtsq*(struc1PM - struc59PM) + s126*struc59PM))/(mtsq - s126)
     &          )*aimag(anomc4_sc)*prop16*prop34

c MP
              res = ampHeavyMPimagc4()
              ampsAnomHeavy(3) = ampsAnomHeavy(3) + (res +
c mass renormalization
     &        (2*im*mtsq*(2*struc13MP + 2*struc27MP + 4*struc3MP + 2*struc53MP + 
     &      2*struc55MP + s126*struc59MP - s16*struc59MP + 4*struc8MP)*
     &    (4 + 3*epinv + 3*Log(MusqHeavy/mtsq)))/(mtsq - s126)**2 +
c c4 renormalization
     &        (epinv*im*(2*struc13MP + s126*struc1MP - s16*struc59MP + 
     &      mtsq*(-struc1MP + struc59MP) + 
     &      2*(struc27MP + 2*struc3MP + struc53MP + struc55MP + 2*struc8MP)))/
     &  (mtsq - s126)
     &          )*aimag(anomc4_sc)*prop16*prop34
          endif

          if (real(anomc7_sc) /= 0._dp) then
c PM
              res = ampHeavyPMrealc7()
              ampsAnomHeavy(2) = ampsAnomHeavy(2) + (res
c c4 renormalization
     &         -((epinv*(2*struc11PM + 2*struc13PM + 2*struc17PM + mtsq*struc1PM - 
     &        2*s12*struc1PM - s126*struc1PM - 2*s13*struc1PM - 
     &        2*s16*struc1PM + 2*struc27PM + 4*struc2PM + 2*struc37PM + 
     &        4*struc3PM + 2*struc55PM - mtsq*struc59PM + s126*struc59PM))/
     &    (mtsq - s126))
     &          )*real(anomc7_sc)*prop16*prop34

c MP
              res = ampHeavyMPrealc7()
              ampsAnomHeavy(3) = ampsAnomHeavy(3) + (res +
c c4 renormalization
     &        (epinv*(2*struc13MP + s126*struc1MP + 2*struc27MP + 4*struc3MP + 
     &      2*struc53MP + 2*struc55MP - s16*struc59MP + 
     &      mtsq*(-struc1MP + struc59MP) + 4*struc8MP))/(mtsq - s126)
     &          )*real(anomc7_sc)*prop16*prop34
          endif

          if (aimag(anomc7_sc) /= 0._dp) then
c PM
              res = ampHeavyPMimagc7()
              ampsAnomHeavy(2) = ampsAnomHeavy(2) + (res +
c c4 renormalization
     &        (epinv*im*(2*struc11PM + 2*struc13PM + 2*struc17PM + mtsq*struc1PM - 
     &      2*s12*struc1PM - s126*struc1PM - 2*s13*struc1PM - 2*s16*struc1PM + 
     &      2*struc27PM + 4*struc2PM + 2*struc37PM + 4*struc3PM + 2*struc55PM - 
     &      mtsq*struc59PM + s126*struc59PM))/(mtsq - s126)
     &          )*aimag(anomc7_sc)*prop16*prop34

c MP
              res = ampHeavyMPimagc7()
              ampsAnomHeavy(3) = ampsAnomHeavy(3) + (res +
c c4 renormalization
     &        (epinv*im*(2*struc13MP + s126*struc1MP + 2*struc27MP + 4*struc3MP + 
     &      2*struc53MP + 2*struc55MP - s16*struc59MP + 
     &      mtsq*(-struc1MP + struc59MP) + 4*struc8MP))/(mtsq - s126)
     &          )*aimag(anomc7_sc)*prop16*prop34
          endif

          endif

c HEAVY ASSEMBLY
          if (enable_lambda2) then
              singletop2_amp_virt_dd = singletop2_amp_virt_dd
c                       ----
     &                  + real(ampSMtree(1)*conjg(ampsAnomHeavy(1)))/lambda**2/sc
     &                  + real(ampsAnomL2tree(1)*conjg(ampsSMHeavy(1)))/lambda**2/sc
c                       ++--, for non-resonant pieces
     &                  + real(ampSMtree(2)*conjg(ampsAnomHeavy(5)))/lambda**2/sc
     &                  + real(ampsAnomL2tree(5)*conjg(ampsSMHeavy(2)))/lambda**2/sc
          endif

          if (enable_lambda4) then
              singletop2_amp_virt_dd = singletop2_amp_virt_dd
     &                  + sum(real(ampsAnomL2tree*conjg(ampsAnomHeavy)))/lambda**4/sc/sc

          endif
      endif

      end

      end module
