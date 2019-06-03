      subroutine qqb_zaj_vdecay(p_mcfm,msq)
        implicit none
        include 'types.f'
        include 'constants.f'
        include 'nf.f'
        include 'mxpart.f'

        real(dp), intent(in) :: p_mcfm(mxpart,4)
        real(dp), intent(out) :: msq(-nf:nf, -nf:nf)

        call qqb_zaj_vdecay_i(1,2,3,4,5,6,p_mcfm,msq)

      end subroutine

      subroutine qqb_zaj_vdecay_cross(p_mcfm,msq)
        implicit none
        include 'types.f'
        include 'constants.f'
        include 'nf.f'
        include 'mxpart.f'

        real(dp), intent(in) :: p_mcfm(mxpart,4)
        real(dp), intent(out) :: msq(-nf:nf, -nf:nf)

        call qqb_zaj_vdecay_i(1,2,3,6,5,4,p_mcfm,msq)

      end subroutine

      subroutine qqb_zaj_vdecay_i(i1,i2,i3,i4,i5,i6,p_mcfm,msq)
          use helicities 
          use VVconfig_m
          use zaj_treeamps_m, only : zaj_tree_isr_pp, zaj_tree_isr_pm,
     &            zaj_tree_isr_mp, zaj_tree_fsr_pp, zaj_tree_fsr_pm
          use omega_ee3jet_m, only : amp_zaj_fsr_pp, amp_zaj_fsr_pm,
     &            zaj_crossings_l_new

          implicit none
          include 'types.f'
          include 'constants.f'
          include 'nf.f'
          include 'mxpart.f'

          integer, intent(in) :: i1,i2,i3,i4,i5,i6
          real(dp), intent(in) :: p_mcfm(mxpart,4)
          real(dp), intent(out) :: msq(-nf:nf, -nf:nf)

          complex(dp) :: amps_q(2,2,2,2), amps_l(2,2,2,2), ampsDress(2,2,2,2,2)
          complex(dp) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          real(dp) :: p(mxpart,4)
          integer :: qi

          real(dp) :: s,t
          s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
          t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

          msq(:,:) = 0._dp
          amps_q = 0._dp
          amps_l = 0._dp
          ampsDress = 0._dp

          p(1,:) = p_mcfm(i2,:) ! quark
          p(2,:) = p_mcfm(i1,:) ! antiquark
          p(3,:) = p_mcfm(i6,:) ! gluon
          p(4,:) = p_mcfm(i5,:) ! photon
          p(5,:) = p_mcfm(i4,:) ! antilepton
          p(6,:) = p_mcfm(i3,:) ! lepton

          call spinoru_s(6,p,za,zb)

c q qb
          call zaj_crossings(6,5,3,4,2,1,za,zb,
     &       zaj_tree_isr_pp, zaj_tree_isr_pm, zaj_tree_isr_mp, amps_q)
          call zaj_crossings_l(6,5,3,4,2,1,za,zb,
     &       zaj_tree_fsr_pp, zaj_tree_fsr_pm, amps_l)
          amps_l = amps_l*t(4,1,2)

          do qi=1,2
            call zaj_treeamps_dress_vdecay(qi,6,5,3,4,2,1,za,zb,amps_q,amps_l,ampsDress)
            msq(qi,-qi) =
     &        + n_light_down*aveqq*sum(abs(ampsDress(1,:,:,:,:))**2)
     &        + n_light_up*aveqq*sum(abs(ampsDress(2,:,:,:,:))**2)
          enddo


c qb q
          call zaj_crossings(6,5,3,4,1,2,za,zb,
     &       zaj_tree_isr_pp, zaj_tree_isr_pm, zaj_tree_isr_mp, amps_q)
          call zaj_crossings_l(6,5,3,4,1,2,za,zb,
     &       zaj_tree_fsr_pp, zaj_tree_fsr_pm, amps_l)
          amps_l = amps_l*t(4,1,2)

          do qi=1,2
            call zaj_treeamps_dress_vdecay(qi,6,5,3,4,1,2,za,zb,amps_q,amps_l,ampsDress)
            msq(-qi,qi) =
     &        + n_light_down*aveqq*sum(abs(ampsDress(1,:,:,:,:))**2)
     &        + n_light_up*aveqq*sum(abs(ampsDress(2,:,:,:,:))**2)
          enddo


          ! other qq channels
          msq(-3,3) = msq(-1,1)
          msq(-4,4) = msq(-2,2)
          msq(-5,5) = msq(-1,1)

          msq(3,-3) = msq(1,-1)
          msq(4,-4) = msq(2,-2)
          msq(5,-5) = msq(1,-1)

          msq(:,:) = msq(:,:) * 8._dp * cA

      end subroutine

c     dresses raw helicity amps obtained from zaj_treeamps with
c       coupling factors, Z propagator, given quark flavor qi
      subroutine zaj_treeamps_dress_vdecay(qi,i1,i2,i3,i4,i5,i6,za,zb,amps_q,amps_l,ampsDress)
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

          integer, intent(in) :: qi
          integer, intent(in) :: i1,i2,i3,i4,i5,i6
          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          complex(dp), intent(in) :: amps_q(2,2,2,2),amps_l(2,2,2,2)
          complex(dp), intent(out) :: ampsDress(2,2,2,2,2)

          integer :: h1,h2,h3,h4
          real(dp) :: vq, vl, Ql
          complex(dp) :: propzQQ, propzQL
          complex(dp) :: propZ_56, propZ_456, propA_56, propA_456

          integer :: flavor
          real(dp) :: s,t
          s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
          t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

          propzQQ = s(i5,i6)/(s(i5,i6) - zmass**2 + im*zmass*zwidth)
          propzQL = t(i4,i5,i6)/(t(i4,i5,i6) - zmass**2 + im*zmass*zwidth)

          propZ_56 = 1/(s(i5,i6) - zmass**2 + im*zmass*zwidth)
          propZ_456 = 1/(t(i4,i5,i6) - zmass**2 + im*zmass*zwidth)
          propA_56 = 1/s(i5,i6)
          propA_456 = 1/t(i4,i5,i6)

          ampsDress = 0._dp
          vl = 0._dp
          vq = 0._dp
          ql = 0._dp

          ! (h4,h3,h2,h1) = helicities of (q,l,g,gam)
          do h1=1,2
            do h2=1,2
              do h3=1,2
                do h4=1,2
                  do flavor=1,2
                    if (h3 == helLeft) then
                      vq = l(qi)
                    elseif (h3 == helRight) then
                      vq = r(qi)
                    endif

                    if (h4 == helLeft) then
                      vl = l(flavor)
                    elseif (h4 == helRight) then
                      vl = r(flavor)
                    endif

                    if (h3 == h4) then
                      Ql = +Q(flavor)
                    else
                      Ql = -Q(flavor)
                    endif

                    ampsDress(flavor,h4,h3,h2,h1) = twort2*sqrt(esq)**3*sqrt(gsq)*(
     &              + amps_q(h4,h3,h2,h1) * Ql*(Q(flavor)*Q(qi) + vl*vq*propzQQ)
     &              + amps_l(h4,h3,h2,h1) * Q(qi)*(Q(flavor)*Q(qi)*propA_456 + vl*vq*propZ_456)
     &                 )

                  enddo
                enddo
              enddo
            enddo
          enddo

      end subroutine

