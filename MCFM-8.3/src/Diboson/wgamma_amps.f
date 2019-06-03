      subroutine wgam_mat(p,msq)
          implicit none
          include 'types.f'
          include 'mxpart.f'
          include 'nf.f'
          include 'constants.f'
          include 'ewcouple.f'
          include 'qcdcouple.f'
          include 'nwz.f'
          include 'ckm.f'

          real(dp), intent(in) :: p(mxpart,4)
          real(dp), intent(out) :: msq(-nf:nf,-nf:nf,0:2)

          complex(dp) :: ampsW(2,2,2,0:2)
          real(dp) :: qqb(0:2), qbq(0:2)
          integer :: j,k,h12,h56,h3
          real(dp) :: facLO

          call wgamma_amps(p,3,4,1,5,2,ampsW)

          facLO = aveqq*gwsq**2*esq*xn/4

          if (nwz == -1) then
            qqb(0) = ( abs(ampsW(1,2,1,0))**2 + abs(ampsW(1,2,2,0))**2 )
            qqb(1) = ( 2*real(conjg(ampsW(1,2,1,0))*ampsW(1,2,1,1))
     &               + 2*real(conjg(ampsW(1,2,2,0))*ampsW(1,2,2,1)))
            qqb(2) = ( abs(ampsW(1,2,1,1))**2 + abs(ampsW(1,2,2,1))**2
     &               + 2*real(conjg(ampsW(1,2,1,0))*ampsW(1,2,1,2))
     &               + 2*real(conjg(ampsW(1,2,2,0))*ampsW(1,2,2,2)))

            qbq(0) = ( abs(ampsW(2,2,1,0))**2 + abs(ampsW(2,2,2,0))**2)
            qbq(1) = ( 2*real(conjg(ampsW(2,2,1,0))*ampsW(2,2,1,1))
     &               + 2*real(conjg(ampsW(2,2,2,0))*ampsW(2,2,2,1)))
            qbq(2) = ( abs(ampsW(2,2,1,1))**2 + abs(ampsW(2,2,2,1))**2
     &               + 2*real(conjg(ampsW(2,2,1,0))*ampsW(2,2,1,2))
     &               + 2*real(conjg(ampsW(2,2,2,0))*ampsW(2,2,2,2)))
          elseif (nwz == +1) then
            qqb(0) = ( abs(ampsW(1,1,1,0))**2 + abs(ampsW(1,1,2,0))**2)
            qqb(1) = ( 2*real(conjg(ampsW(1,1,1,0))*ampsW(1,1,1,1))
     &               + 2*real(conjg(ampsW(1,1,2,0))*ampsW(1,1,2,1)))
            qqb(2) = ( abs(ampsW(1,1,1,1))**2 + abs(ampsW(1,1,2,1))**2
     &               + 2*real(conjg(ampsW(1,1,1,0))*ampsW(1,1,1,2))
     &               + 2*real(conjg(ampsW(1,1,2,0))*ampsW(1,1,2,2)))

            qbq(0) = ( abs(ampsW(2,1,1,0))**2 + abs(ampsW(2,1,2,0))**2)
            qbq(1) = ( 2*real(conjg(ampsW(2,1,1,0))*ampsW(2,1,1,1))
     &               + 2*real(conjg(ampsW(2,1,2,0))*ampsW(2,1,2,1)))
            qbq(2) = ( abs(ampsW(2,1,1,1))**2 + abs(ampsW(2,1,2,1))**2
     &               + 2*real(conjg(ampsW(2,1,1,0))*ampsW(2,1,1,2))
     &               + 2*real(conjg(ampsW(2,1,2,0))*ampsW(2,1,2,2)))
          endif

          msq(:,:,:) = 0

          do j=-nf,nf
            do k=-nf,nf
              if (j>0 .and. k<0) then
                  msq(j,k,:) = Vsq(j,k)*qqb(:)
              elseif (j<0 .and. k>0) then
                  msq(j,k,:) = Vsq(j,k)*qbq(:)
              endif
            enddo
          enddo

          msq(:,:,0) = msq(:,:,0) * facLO
          msq(:,:,1) = msq(:,:,1) * facLO * ason2pi
          msq(:,:,2) = msq(:,:,2) * facLO * ason2pi**2

      end subroutine

      subroutine wgamma_amps(p,i5,i6,i1,i3,i2, ampsW)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'ewcharge.f'
      include 'masses.f' ! wmass, wwidth
      include 'constants.f' ! im
      include 'epinv.f'
      include 'epinv2.f'
      include 'scale.f'

      real(dp), intent(in) :: p(mxpart,4)
      integer, intent(in) :: i5,i6,i1,i3,i2
      real(dp) :: mp(mxpart,4)
      complex(dp), intent(out) :: ampsW(2,2,2,0:2)

      complex(dp) :: ampsWf(2,2,2,0:2)
      complex(dp) :: amp_rr_W(0:2)

      complex(dp) :: wzgamma_amp_rr
      complex(dp) :: za(mxpart,mxpart), zb(mxpart,mxpart)

      complex(dp) :: c_alphai(0:2,4)
      complex(dp) :: c_betai(0:2,4)
      complex(dp) :: c_gammai(0:2,4)

      complex(dp) :: c_alpha, c_beta, c_gamma
      complex(dp) :: c_alpha_cc, c_beta_cc, c_gamma_cc

      complex(dp) :: c_alpha_1l, c_beta_1l, c_gamma_1l
      complex(dp) :: c_alpha_cc_1l, c_beta_cc_1l, c_gamma_cc_1l

      real(dp) :: eqpr, eq
      complex(dp) :: amp_rr(0:2)

      complex(dp) :: qbqamp(2,2,2)
      complex(dp) :: qbqamp_nlo(2,2,2)

      real(dp) :: dotvec
      real(dp) :: s123, s34, s56
      real(dp) :: s12, s13, s23
      complex(dp) :: agamtree, agamvirt
      complex(dp) :: DW, DWf
      integer i,j,k,h12,h34

      ! eqpr -  eq = 1
      eq = -2._dp/3
      eqpr = 1._dp/3

      mp(:,:) = p(:,:)
      mp(1,:) = -p(i2,:)
      mp(2,:) = -p(i1,:)
      mp(3,:) = -p(i3,:)
      mp(5,:) = p(i5,:)
      mp(6,:) = p(i6,:)

      call spinoru_s(6,mp,za,zb)

      ampsW(:,:,:,:) = 0._dp
      ampsWf(:,:,:,:) = 0._dp

      s56 = dotvec(mp(5,:)+mp(6,:),mp(5,:)+mp(6,:))
      s12 = dotvec(mp(1,:)+mp(2,:),mp(1,:)+mp(2,:))

      DW = s56 - wmass**2 + im*wwidth*wmass
      DWf = s12 - wmass**2 + im*wwidth*wmass

      ! first block, 56132 ordering
      call omega_wzgamma(mp,5,6,1,3,2,c_alphai,c_betai,c_gammai)

      call wgamma_amp_rr(5,6,1,3,2,za,zb,c_alphai,c_betai,c_gammai,
     &                           eq,eqpr,amp_rr_W)
      ampsW(2,2,2,:) = amp_rr_W(:)/Dw

      call wgamma_ampf_rr(5,6,1,3,2,za,zb,2,amp_rr_W)
      ampsWf(2,2,2,:) = amp_rr_W(:)/Dwf
      call wgamma_ampf_rr(5,6,1,3,2,za,zb,1,amp_rr_W)
      ampsWf(2,2,1,:) = amp_rr_W(:)/Dwf

      call wgamma_amp_rr(5,6,1,3,2,za,zb,c_alphai,c_betai,c_gammai,
     &                             -eqpr,-eq,amp_rr_W)
      ampsW(2,1,1,:) = -amp_rr_W(:)/Dw

      ! second block, 65231 ordering
      call omega_wzgamma(mp,6,5,2,3,1,c_alphai,c_betai,c_gammai)

      call wgamma_amp_rr(6,5,2,3,1,zb,za,c_alphai,c_betai,c_gammai,
     &                             -eqpr,-eq,amp_rr_W)
      ampsW(2,2,1,:) = -amp_rr_W(:)/Dw

      call wgamma_amp_rr(6,5,2,3,1,zb,za,c_alphai,c_betai,c_gammai,
     &                             eq,eqpr,amp_rr_W)
      ampsW(2,1,2,:) = amp_rr_W(:)/Dw

      call wgamma_ampf_rr(6,5,2,3,1,zb,za,2,amp_rr_W)
      ampsWf(2,1,2,:) = amp_rr_W(:)/Dwf
      call wgamma_ampf_rr(6,5,2,3,1,zb,za,1,amp_rr_W)
      ampsWf(2,1,1,:) = amp_rr_W(:)/Dwf


      ! third block, 56231 ordering
      call omega_wzgamma(mp,5,6,2,3,1,c_alphai,c_betai,c_gammai)

      call wgamma_amp_rr(5,6,2,3,1,za,zb,c_alphai,c_betai,c_gammai,
     &                             eq,eqpr,amp_rr_W)
      ampsW(1,2,2,:) = amp_rr_W(:)/Dw

      call wgamma_ampf_rr(5,6,2,3,1,za,zb,2,amp_rr_W)
      ampsWf(1,2,2,:) = amp_rr_W(:)/Dwf
      call wgamma_ampf_rr(5,6,2,3,1,za,zb,1,amp_rr_W)
      ampsWf(1,2,1,:) = amp_rr_W(:)/Dwf

      call wgamma_amp_rr(5,6,2,3,1,za,zb,c_alphai,c_betai,c_gammai,
     &                             -eqpr,-eq,amp_rr_W)
      ampsW(1,1,1,:) = -amp_rr_W(:)/Dw

      ! fourth and last block, 65132 ordering
      call omega_wzgamma(mp,6,5,1,3,2,c_alphai,c_betai,c_gammai)

      call wgamma_amp_rr(6,5,1,3,2,zb,za,c_alphai,c_betai,c_gammai,
     &                             -eqpr,-eq,amp_rr_W)
      ampsW(1,2,1,:) = -amp_rr_W(:)/Dw

      call wgamma_amp_rr(6,5,1,3,2,zb,za,c_alphai,c_betai,c_gammai,
     &                             eq,eqpr,amp_rr_W)
      ampsW(1,1,2,:) = amp_rr_W(:)/Dw

      call wgamma_ampf_rr(6,5,1,3,2,zb,za,2,amp_rr_W)
      ampsWf(1,1,2,:) = amp_rr_W(:)/Dwf
      call wgamma_ampf_rr(6,5,1,3,2,zb,za,1,amp_rr_W)
      ampsWf(1,1,1,:) = amp_rr_W(:)/Dwf

      ampsW(:,:,:,:) = ampsW(:,:,:,:) + ampsWf(:,:,:,:)

      return

c--- debug below

      write (*,*) "----------------------"

      s12 = dotvec(mp(1,:)+mp(2,:),mp(1,:)+mp(2,:))
      s13 = dotvec(mp(1,:)+mp(3,:),mp(1,:)+mp(3,:))
      s23 = dotvec(mp(2,:)+mp(3,:),mp(2,:)+mp(3,:))
      s123 = dotvec(mp(1,:)+mp(2,:)+mp(3,:),mp(1,:)+mp(2,:)+mp(3,:))
      s34 =  dotvec(mp(5,:)+mp(6,:),mp(5,:)+mp(6,:))

      call spinoru(5,p,za,zb)

      !born amps
      qbqamp(2,2,1) = agamtree(1,2,3,4,5,za,zb,-1)
      qbqamp(2,2,2) = agamtree(1,2,3,4,5,za,zb,+1)

      qbqamp(1,2,1) = agamtree(2,1,3,4,5,za,zb,-1)
      qbqamp(1,2,2) = agamtree(2,1,3,4,5,za,zb,+1)

      qbqamp(2,1,1) = agamtree(2,1,4,3,5,zb,za,-1)
      qbqamp(2,1,2) = agamtree(2,1,4,3,5,zb,za,+1)

      qbqamp(1,1,1) = agamtree(1,2,4,3,5,zb,za,-1)
      qbqamp(1,1,2) = agamtree(1,2,4,3,5,zb,za,+1)

      !loop amps
      qbqamp_nlo(2,2,1) = agamvirt(1,2,3,4,5,za,zb,-1)
      qbqamp_nlo(2,2,2) = agamvirt(1,2,3,4,5,za,zb,+1)

      qbqamp_nlo(1,2,1) = agamvirt(2,1,3,4,5,za,zb,-1)
      qbqamp_nlo(1,2,2) = agamvirt(2,1,3,4,5,za,zb,+1)

      qbqamp_nlo(2,1,1) = agamvirt(2,1,4,3,5,zb,za,-1)
      qbqamp_nlo(2,1,2) = agamvirt(2,1,4,3,5,zb,za,+1)

      qbqamp_nlo(1,1,1) = agamvirt(1,2,4,3,5,zb,za,-1)
      qbqamp_nlo(1,1,2) = agamvirt(1,2,4,3,5,zb,za,+1)


      write (*,*) "ratios born"
      do h12=1,2
        do h34=1,2
            write (*,*) "h12,h34,h5", h12,h34,
     &        2*sqrt(2._dp)*qbqamp(h12,h34,1) / ampsWf(h12,h34,1,0),
     &        2*sqrt(2._dp)*qbqamp(h12,h34,2) / ampsWf(h12,h34,2,0)
        enddo
      enddo

      write (*,*) "ratios loop"
      do h12=1,2
        do h34=1,2
            write (*,*) "h12,h34,h5", h12,h34,
     & 2*cf*sqrt(2._dp)*qbqamp_nlo(h12,h34,1) / ampsWf(h12,h34,1,1),
     & 2*cf*sqrt(2._dp)*qbqamp_nlo(h12,h34,2) / ampsWf(h12,h34,2,1)
        enddo
      enddo


      pause

      end subroutine

      subroutine wgamma_ampf_RR(i5,i6,i1,i3,i2,za,zb,hgamma,amps)
         use VVconfig_m
        implicit none
        include 'types.f'
        include 'mxpart.f'
        include 'constants.f'
        include 'scale.f'
        include 'nf.f'

        integer, intent(in) :: i5,i6,i1,i3,i2
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        integer, intent(in) :: hgamma
        complex(dp), intent(out) :: amps(0:2)

        real(dp) :: s12, s56
        complex(dp) :: Fq(0:2), agamtree
        integer, parameter :: scheme = schemeMSBAR

        s12 = real(za(i1,i2)*zb(i2,i1))
        s56 = real(za(i5,i6)*zb(i6,i5))
        call quark_formfactor(s12, Fq)

        if (hgamma == 1) then
          agamtree=-zb(i1,i6)**2/(s12-s56)*(
     &      Qu*(-za(i6,i3)/(zb(i2,i1)*zb(i5,i3)))
     &     +Qd*(+za(i6,i3)/(zb(i2,i1)*zb(i5,i3))))
        elseif (hgamma == 2) then
          agamtree=-za(i2,i5)**2/(s12-s56)*(
     &      Qd*(+zb(i6,i3)/(za(i1,i2)*za(i5,i3)))
     &     +Qu*(-zb(i6,i3)/(za(i1,i2)*za(i5,i3))))
        endif

        amps(0) = 2*sqrt(2._dp)*Fq(0)*agamtree*s12
        amps(1) = Fq(1) * amps(0)
        amps(2) = Fq(2) * amps(0)

        if (scheme == schemeMCFM) then
            amps(1) = amps(1) + cf*(-pisq/12 + 0.5)*amps(0)
        elseif (scheme == schemeMSBAR) then
            ! already prepared for MSBAR scheme
        else
            print *, "undefined scheme for final state radiation"
            call exit(1)
        endif

      end subroutine

      subroutine wgamma_amp_RR(i5,i6,i1,i3,i2,za,zb,
     &   c_alphai, c_betai, c_gammai, eq,eqpr,amps)
        use VVconfig_m
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'nf.f'

      integer, intent(in) :: i5,i6,i1,i3,i2
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
      complex(dp), intent(in) :: c_alphai(0:2,4),
     &                           c_betai(0:2,4), c_gammai(0:2,4)
      real(dp), intent(in) :: eq,eqpr
c     integer (kind(bosZ)), intent(in) :: bos
      complex(dp), intent(out) :: amps(0:2)

      complex(dp) :: c_alpha(0:2)
      complex(dp) :: c_beta(0:2)
      complex(dp) :: c_gamma(0:2)
      real(dp) :: s12, s123, dotvec
      complex(dp) :: wzgamma_amp_RR, omega_I1

      integer, parameter :: scheme = schemeMSBAR

      s12 = real(za(i1,i2)*zb(i2,i1))
      s123 = real(za(i1,i2)*zb(i2,i1) + za(i1,i3)*zb(i3,i1) +
     &                             za(i2,i3)*zb(i3,i2))

      s12 = real(za(i1,i2)*zb(i2,i1))

      c_alpha(:) = eqpr*c_alphai(:,1) + eq*c_alphai(:,2) + c_alphai(:,3)
      c_beta(:) = eqpr*c_betai(:,1) + eq*c_betai(:,2) + c_betai(:,3)
      c_gamma(:) = eqpr*c_gammai(:,1) + eq*c_gammai(:,2) + c_gammai(:,3)

      amps(0) = wzgamma_amp_RR(i5,i6,i1,i3,i2,za,zb,
     &             c_alpha(0),c_beta(0),c_gamma(0))

      amps(1) = wzgamma_amp_RR(i5,i6,i1,i3,i2,za,zb,
     &             c_alpha(1),c_beta(1),c_gamma(1))

      amps(2) = wzgamma_amp_RR(i5,i6,i1,i3,i2,za,zb,
     &             c_alpha(2),c_beta(2),c_gamma(2))

      if (scheme == schemeMCFM) then
          amps(1) = amps(1) + omega_I1(s12)*amps(0)
      elseif (scheme == schemeMSBAR) then
          call zgamma_catani_to_msbar(i5,i6,i1,i3,i2,za,zb,amps,s123)
          call evolve_msbar_IR_qq(i5,i6,i1,i3,i2,za,zb,amps,s123)
      else
          print *, "scheme unsupported"
          call exit(1)
      endif

      end subroutine

