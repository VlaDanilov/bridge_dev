      subroutine qqb_zaj(p_mcfm,msq)
          use helicities 
          use VVconfig_m
          use zaj_treeamps_m, only : zaj_tree_isr_pp, zaj_tree_isr_pm,
     &            zaj_tree_isr_mp, zaj_tree_fsr_pp, zaj_tree_fsr_pm,
     &            zaj_tree_anomZZ_pp, zaj_tree_anomZZ_pm,
     &            zaj_tree_anomZZ_mp, zaj_tree_anomZa_pp,
     &            zaj_tree_anomZa_pm, zaj_tree_anomZa_mp,
     &            zaj_tree_anomaZ_pp, zaj_tree_anomaZ_pm,
     &            zaj_tree_anomaZ_mp
          use omega_ee3jet_m, only : amp_zaj_fsr_pp, amp_zaj_fsr_pm,
     &            amp_zaj_anomZZ_pp, amp_zaj_anomZZ_pm,
     &            amp_zaj_anomZa_pp, amp_zaj_anomZa_pm,
     &            amp_zaj_anomaZ_pp, amp_zaj_anomaZ_pm,
     &            zaj_crossings_l_new,
     &            zaj_crossings_anom_new

          implicit none
          include 'types.f'
          include 'constants.f'
          include 'nf.f'
          include 'mxpart.f'
          include 'anomcoup.f' ! logical :: anomtgc

          real(dp), intent(in) :: p_mcfm(mxpart,4)
          real(dp), intent(out) :: msq(-nf:nf, -nf:nf)

          complex(dp) :: amps_q(2,2,2,2), amps_l(2,2,2,2), ampsDress(2,2,2,2)
          !complex(dp) :: amps_l_new(0:1, 2,2,2,2)
          !complex(dp) :: amps_anomZZ(2,2,2,2), amps_anomZa(2,2,2,2)
          !complex(dp) :: amps_anomaZ(2,2,2,2)
          complex(dp) :: amps_anomZZ_new(0:1, 2,2,2,2)
          complex(dp) :: amps_anomZa_new(0:1, 2,2,2,2)
          complex(dp) :: amps_anomaZ_new(0:1, 2,2,2,2)
          complex(dp) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          real(dp) :: p(mxpart,4)
          integer :: qi

          integer (kind(decayElAntiEl)) :: vDecay

          integer :: i1,i2,i3
          real(dp) :: s,t
          s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
          t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

          msq(:,:) = 0._dp
          amps_q = 0._dp
          amps_l = 0._dp
          ampsDress = 0._dp

          !amps_anomZZ = 0._dp
          !amps_anomZa = 0._dp
          !amps_anomaZ = 0._dp

          amps_anomZZ_new = 0._dp
          amps_anomZa_new = 0._dp
          amps_anomaZ_new = 0._dp

          vDecay = decayChannel()

          p(1,:) = p_mcfm(2,:) ! quark
          p(2,:) = p_mcfm(1,:) ! antiquark
          p(3,:) = p_mcfm(6,:) ! gluon
          p(4,:) = p_mcfm(5,:) ! photon
          p(5,:) = p_mcfm(4,:) ! antilepton
          p(6,:) = p_mcfm(3,:) ! lepton

          call spinoru_s(6,p,za,zb)

          call zaj_crossings(1,2,3,4,5,6,za,zb,
     &       zaj_tree_isr_pp, zaj_tree_isr_pm, zaj_tree_isr_mp, amps_q)
          if (vDecay == decayElAntiEl .or. vDecay == decayQuarks) then
            call zaj_crossings_l(1,2,3,4,5,6,za,zb,
     &         zaj_tree_fsr_pp, zaj_tree_fsr_pm, amps_l)
            amps_l = amps_l*t(4,5,6)
          endif

          !call zaj_crossings_l_new(p, 1,2,3,4,5,6,za,zb,amps_l_new)

          if (anomtgc) then
c           call zaj_crossings_anom(1,2,3,4,5,6,za,zb,
c    &         zaj_tree_anomZZ_pp, zaj_tree_anomZZ_pm, zaj_tree_anomZZ_mp,
c    &         amps_anomZZ)
c           call zaj_crossings_anom(1,2,3,4,5,6,za,zb,
c    &         zaj_tree_anomZa_pp, zaj_tree_anomZa_pm, zaj_tree_anomZa_mp,
c    &         amps_anomZa)
c           call zaj_crossings_anom(1,2,3,4,5,6,za,zb,
c    &         zaj_tree_anomaZ_pp, zaj_tree_anomaZ_pm, zaj_tree_anomaZ_mp,
c    &         amps_anomaZ)

            call zaj_crossings_anom_new(p, 1,2,3,4,5,6,za,zb,
     &          amp_zaj_anomZZ_pp, amp_zaj_anomZZ_pm, amps_anomZZ_new)
            call zaj_crossings_anom_new(p, 1,2,3,4,5,6,za,zb,
     &          amp_zaj_anomZa_pp, amp_zaj_anomZa_pm, amps_anomZa_new)
            call zaj_crossings_anom_new(p, 1,2,3,4,5,6,za,zb,
     &          amp_zaj_anomaZ_pp, amp_zaj_anomaZ_pm, amps_anomaZ_new)

          endif
          do qi=1,2
            call zaj_treeamps_dress(qi,za,zb,amps_q,amps_l,
     &              amps_anomZZ_new(0,:,:,:,:), amps_anomZa_new(0,:,:,:,:),
     &              amps_anomaZ_new(0,:,:,:,:), ampsDress)
            msq(qi,-qi) = aveqq*sum(abs(ampsDress)**2)
          enddo

          call zaj_crossings(2,1,3,4,5,6,za,zb,
     &       zaj_tree_isr_pp, zaj_tree_isr_pm, zaj_tree_isr_mp, amps_q)
          if (vDecay == decayElAntiEl .or. vDecay == decayQuarks) then
            call zaj_crossings_l(2,1,3,4,5,6,za,zb,
     &         zaj_tree_fsr_pp, zaj_tree_fsr_pm, amps_l)
            amps_l = amps_l*t(4,5,6)
          endif

          if (anomtgc) then
c           call zaj_crossings_anom(2,1,3,4,5,6,za,zb,
c    &         zaj_tree_anomZZ_pp, zaj_tree_anomZZ_pm, zaj_tree_anomZZ_mp,
c    &         amps_anomZZ)
c           call zaj_crossings_anom(2,1,3,4,5,6,za,zb,
c    &         zaj_tree_anomZa_pp, zaj_tree_anomZa_pm, zaj_tree_anomZa_mp,
c    &         amps_anomZa)
c           call zaj_crossings_anom(2,1,3,4,5,6,za,zb,
c    &         zaj_tree_anomaZ_pp, zaj_tree_anomaZ_pm, zaj_tree_anomaZ_mp,
c    &         amps_anomaZ)

            call zaj_crossings_anom_new(p, 2,1,3,4,5,6,za,zb,
     &          amp_zaj_anomZZ_pp, amp_zaj_anomZZ_pm, amps_anomZZ_new)
            call zaj_crossings_anom_new(p, 2,1,3,4,5,6,za,zb,
     &          amp_zaj_anomZa_pp, amp_zaj_anomZa_pm, amps_anomZa_new)
            call zaj_crossings_anom_new(p, 2,1,3,4,5,6,za,zb,
     &          amp_zaj_anomaZ_pp, amp_zaj_anomaZ_pm, amps_anomaZ_new)

          endif
          do qi=1,2
            call zaj_treeamps_dress(qi,za,zb,amps_q,amps_l,
     &              amps_anomZZ_new(0,:,:,:,:), amps_anomZa_new(0,:,:,:,:),
     &              amps_anomaZ_new(0,:,:,:,:), ampsDress)
            msq(-qi,qi) = aveqq*sum(abs(ampsDress)**2)
          enddo

          call zaj_crossings(3,2,1,4,5,6,za,zb,
     &       zaj_tree_isr_pp, zaj_tree_isr_pm, zaj_tree_isr_mp, amps_q)
          if (vDecay == decayElAntiEl .or. vDecay == decayQuarks) then
            call zaj_crossings_l(3,2,1,4,5,6,za,zb,
     &         zaj_tree_fsr_pp, zaj_tree_fsr_pm, amps_l)
            amps_l = amps_l*t(4,5,6)
          endif

          if (anomtgc) then
c           call zaj_crossings_anom(3,2,1,4,5,6,za,zb,
c    &         zaj_tree_anomZZ_pp, zaj_tree_anomZZ_pm, zaj_tree_anomZZ_mp,
c    &         amps_anomZZ)
c           call zaj_crossings_anom(3,2,1,4,5,6,za,zb,
c    &         zaj_tree_anomZa_pp, zaj_tree_anomZa_pm, zaj_tree_anomZa_mp,
c    &         amps_anomZa)
c           call zaj_crossings_anom(3,2,1,4,5,6,za,zb,
c    &         zaj_tree_anomaZ_pp, zaj_tree_anomaZ_pm, zaj_tree_anomaZ_mp,
c    &         amps_anomaZ)

            call zaj_crossings_anom_new(p, 3,2,1,4,5,6,za,zb,
     &          amp_zaj_anomZZ_pp, amp_zaj_anomZZ_pm, amps_anomZZ_new)
            call zaj_crossings_anom_new(p, 3,2,1,4,5,6,za,zb,
     &          amp_zaj_anomZa_pp, amp_zaj_anomZa_pm, amps_anomZa_new)
            call zaj_crossings_anom_new(p, 3,2,1,4,5,6,za,zb,
     &          amp_zaj_anomaZ_pp, amp_zaj_anomaZ_pm, amps_anomaZ_new)

          endif
          do qi=1,2
            call zaj_treeamps_dress(qi,za,zb,amps_q,amps_l,
     &              amps_anomZZ_new(0,:,:,:,:), amps_anomZa_new(0,:,:,:,:),
     &              amps_anomaZ_new(0,:,:,:,:), ampsDress)
            msq(qi,0) = aveqg*sum(abs(ampsDress)**2)
          enddo

          call zaj_crossings(2,3,1,4,5,6,za,zb,
     &       zaj_tree_isr_pp, zaj_tree_isr_pm, zaj_tree_isr_mp, amps_q)
          if (vDecay == decayElAntiEl .or. vDecay == decayQuarks) then
            call zaj_crossings_l(2,3,1,4,5,6,za,zb,
     &         zaj_tree_fsr_pp, zaj_tree_fsr_pm, amps_l)
            amps_l = amps_l*t(4,5,6)
          endif

          if (anomtgc) then
c           call zaj_crossings_anom(2,3,1,4,5,6,za,zb,
c    &         zaj_tree_anomZZ_pp, zaj_tree_anomZZ_pm, zaj_tree_anomZZ_mp,
c    &         amps_anomZZ)
c           call zaj_crossings_anom(2,3,1,4,5,6,za,zb,
c    &         zaj_tree_anomZa_pp, zaj_tree_anomZa_pm, zaj_tree_anomZa_mp,
c    &         amps_anomZa)
c           call zaj_crossings_anom(2,3,1,4,5,6,za,zb,
c    &         zaj_tree_anomaZ_pp, zaj_tree_anomaZ_pm, zaj_tree_anomaZ_mp,
c    &         amps_anomaZ)

            call zaj_crossings_anom_new(p, 2,3,1,4,5,6,za,zb,
     &          amp_zaj_anomZZ_pp, amp_zaj_anomZZ_pm, amps_anomZZ_new)
            call zaj_crossings_anom_new(p, 2,3,1,4,5,6,za,zb,
     &          amp_zaj_anomZa_pp, amp_zaj_anomZa_pm, amps_anomZa_new)
            call zaj_crossings_anom_new(p, 2,3,1,4,5,6,za,zb,
     &          amp_zaj_anomaZ_pp, amp_zaj_anomaZ_pm, amps_anomaZ_new)

          endif
          do qi=1,2
            call zaj_treeamps_dress(qi,za,zb,amps_q,amps_l,
     &              amps_anomZZ_new(0,:,:,:,:), amps_anomZa_new(0,:,:,:,:),
     &              amps_anomaZ_new(0,:,:,:,:), ampsDress)
            msq(-qi,0) = aveqg*sum(abs(ampsDress)**2)
          enddo

          call zaj_crossings(3,1,2,4,5,6,za,zb,
     &       zaj_tree_isr_pp, zaj_tree_isr_pm, zaj_tree_isr_mp, amps_q)
          if (vDecay == decayElAntiEl .or. vDecay == decayQuarks) then
            call zaj_crossings_l(3,1,2,4,5,6,za,zb,
     &         zaj_tree_fsr_pp, zaj_tree_fsr_pm, amps_l)
            amps_l = amps_l*t(4,5,6)
          endif

          if (anomtgc) then
c           call zaj_crossings_anom(3,1,2,4,5,6,za,zb,
c    &         zaj_tree_anomZZ_pp, zaj_tree_anomZZ_pm, zaj_tree_anomZZ_mp,
c    &         amps_anomZZ)
c           call zaj_crossings_anom(3,1,2,4,5,6,za,zb,
c    &         zaj_tree_anomZa_pp, zaj_tree_anomZa_pm, zaj_tree_anomZa_mp,
c    &         amps_anomZa)
c           call zaj_crossings_anom(3,1,2,4,5,6,za,zb,
c    &         zaj_tree_anomaZ_pp, zaj_tree_anomaZ_pm, zaj_tree_anomaZ_mp,
c    &         amps_anomaZ)

            call zaj_crossings_anom_new(p, 3,1,2,4,5,6,za,zb,
     &          amp_zaj_anomZZ_pp, amp_zaj_anomZZ_pm, amps_anomZZ_new)
            call zaj_crossings_anom_new(p, 3,1,2,4,5,6,za,zb,
     &          amp_zaj_anomZa_pp, amp_zaj_anomZa_pm, amps_anomZa_new)
            call zaj_crossings_anom_new(p, 3,1,2,4,5,6,za,zb,
     &          amp_zaj_anomaZ_pp, amp_zaj_anomaZ_pm, amps_anomaZ_new)

          endif
          do qi=1,2
            call zaj_treeamps_dress(qi,za,zb,amps_q,amps_l,
     &              amps_anomZZ_new(0,:,:,:,:), amps_anomZa_new(0,:,:,:,:),
     &              amps_anomaZ_new(0,:,:,:,:), ampsDress)
            msq(0,qi) = aveqg*sum(abs(ampsDress)**2)
          enddo

          call zaj_crossings(1,3,2,4,5,6,za,zb,
     &       zaj_tree_isr_pp, zaj_tree_isr_pm, zaj_tree_isr_mp, amps_q)
          if (vDecay == decayElAntiEl .or. vDecay == decayQuarks) then
            call zaj_crossings_l(1,3,2,4,5,6,za,zb,
     &         zaj_tree_fsr_pp, zaj_tree_fsr_pm, amps_l)
            amps_l = amps_l*t(4,5,6)
          endif

          if (anomtgc) then
c           call zaj_crossings_anom(1,3,2,4,5,6,za,zb,
c    &         zaj_tree_anomZZ_pp, zaj_tree_anomZZ_pm, zaj_tree_anomZZ_mp,
c    &         amps_anomZZ)
c           call zaj_crossings_anom(1,3,2,4,5,6,za,zb,
c    &         zaj_tree_anomZa_pp, zaj_tree_anomZa_pm, zaj_tree_anomZa_mp,
c    &         amps_anomZa)
c           call zaj_crossings_anom(1,3,2,4,5,6,za,zb,
c    &         zaj_tree_anomaZ_pp, zaj_tree_anomaZ_pm, zaj_tree_anomaZ_mp,
c    &         amps_anomaZ)

            call zaj_crossings_anom_new(p, 1,3,2,4,5,6,za,zb,
     &          amp_zaj_anomZZ_pp, amp_zaj_anomZZ_pm, amps_anomZZ_new)
            call zaj_crossings_anom_new(p, 1,3,2,4,5,6,za,zb,
     &          amp_zaj_anomZa_pp, amp_zaj_anomZa_pm, amps_anomZa_new)
            call zaj_crossings_anom_new(p, 1,3,2,4,5,6,za,zb,
     &          amp_zaj_anomaZ_pp, amp_zaj_anomaZ_pm, amps_anomaZ_new)

          endif
          do qi=1,2
            call zaj_treeamps_dress(qi,za,zb,amps_q,amps_l,
     &              amps_anomZZ_new(0,:,:,:,:), amps_anomZa_new(0,:,:,:,:),
     &              amps_anomaZ_new(0,:,:,:,:), ampsDress)
            msq(0,-qi) = aveqg*sum(abs(ampsDress)**2)
          enddo

          ! other qq channels
          msq(-3,3) = msq(-1,1)
          msq(-4,4) = msq(-2,2)
          msq(-5,5) = msq(-1,1)

          msq(3,-3) = msq(1,-1)
          msq(4,-4) = msq(2,-2)
          msq(5,-5) = msq(1,-1)

          ! other qg channels
          msq(-3,0) = msq(-1,0)
          msq(-4,0) = msq(-2,0)
          msq(-5,0) = msq(-1,0)

          msq(3,0) = msq(1,0)
          msq(4,0) = msq(2,0)
          msq(5,0) = msq(1,0)

          ! other gq channels
          msq(0,3) = msq(0,1)
          msq(0,4) = msq(0,2)
          msq(0,5) = msq(0,1)

          msq(0,-3) = msq(0,-1)
          msq(0,-4) = msq(0,-2)
          msq(0,-5) = msq(0,-1)

          msq(:,:) = msq(:,:) * 8._dp

      end subroutine

c     dresses raw helicity amps obtained from zaj_treeamps with
c       coupling factors, Z propagator, given quark flavor qi
      subroutine zaj_treeamps_dress(qi,za,zb,amps_q,amps_l,
     &             amps_anomZZ, amps_anomZa, amps_anomaZ, ampsDress)
          use helicities 

          implicit none
          include 'types.f'
          include 'constants.f'
          include 'mxpart.f'
          include 'nf.f'
          include 'zcouple.f'
          include 'ewcouple.f'
          include 'ewcharge.f'
          include 'qcdcouple.f'
          include 'masses.f'
          include 'anomcoup.f' ! anomtgc

          integer, intent(in) :: qi
          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          complex(dp), intent(in) :: amps_q(2,2,2,2),amps_l(2,2,2,2)
          complex(dp), intent(in) :: amps_anomZZ(2,2,2,2)
          complex(dp), intent(in) :: amps_anomZa(2,2,2,2)
          complex(dp), intent(in) :: amps_anomaZ(2,2,2,2)
          complex(dp), intent(out) :: ampsDress(2,2,2,2)

          integer :: h1,h2,h3,h4
          real(dp) :: vq, vl, Ql
          complex(dp) :: propzQQ, propzQL
          complex(dp) :: propZ_56, propZ_456, propA_56, propA_456

          integer :: i1,i2,i3
          real(dp) :: s,t
          s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
          t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

          propzQQ = s(5,6)/(s(5,6) - zmass**2 + im*zmass*zwidth)
          propzQL = t(4,5,6)/(t(4,5,6) - zmass**2 + im*zmass*zwidth)
          ! used for old anomZZ amplitudes
          !cprop = (t(4,5,6)-s(5,6))/(t(4,5,6) - zmass**2 + im*zmass*zwidth)

          propZ_56 = 1/(s(5,6) - zmass**2 + im*zmass*zwidth)
          propZ_456 = 1/(t(4,5,6) - zmass**2 + im*zmass*zwidth)
          propA_56 = 1/s(5,6)
          propA_456 = 1/t(4,5,6)

          ampsDress = 0._dp

          ! (h4,h3,h2,h1) = helicities of (q,l,g,gam)
          do h1=1,2
            do h2=1,2
              do h3=1,2
                if (h3 == helLeft) then
                  vl = l1
                elseif (h3 == helRight) then
                  vl = r1
                endif

                do h4=1,2
                  if (h4 == helLeft) then
                    vq = l(qi)
                  elseif (h4 == helRight) then
                    vq = r(qi)
                  endif

                  if (h3 == h4) then
                    Ql = -abs(q1)
                  else
                    Ql = +abs(q1)
                  endif

                  ampsDress(h4,h3,h2,h1) = twort2*sqrt(esq)**3*sqrt(gsq)*(
c XXX comment next two lines
     &              amps_q(h4,h3,h2,h1) * Q(qi)*(q1*Q(qi) + vl*vq*propzQQ)
     &            + amps_l(h4,h3,h2,h1) * Ql*(q1*Q(qi)*propA_456 + vl*vq*propZ_456)
     &               )

                  if (anomtgc) then
                    ampsDress(h4,h3,h2,h1) = ampsDress(h4,h3,h2,h1) +
     &                  twort2*sqrt(esq)**3*sqrt(gsq)*(
     &              + amps_anomZZ(h4,h3,h2,h1) * vl*vq * propZ_456 * propZ_56
     &              + amps_anomZa(h4,h3,h2,h1) * q1*vq * propZ_456 * propA_56
     &              + amps_anomaZ(h4,h3,h2,h1) * Q(qi)*vl * propA_456 * propZ_56
     &               )
                  endif

                enddo
              enddo
            enddo
          enddo

      end subroutine

