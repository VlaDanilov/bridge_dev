      subroutine zgam_mat(p,msq)
        use VVconfig_m, only: schemeMCFM, zgam_scheme, decayChannel,
     &                    decayQuarks, schemeMSBAR
          implicit none
          include 'types.f'
          include 'mxpart.f'
          include 'nf.f'
          include 'constants.f'
          include 'ewcouple.f'
          include 'qcdcouple.f'

          real(dp), intent(in) :: p(mxpart,4)
          real(dp), intent(out) :: msq(-nf:nf,-nf:nf,0:2)

          complex(dp) :: ampsZ(-nf:nf,-nf:nf,2,2,2,0:2)
          complex(dp) :: ampsZ_du(-nf:nf,-nf:nf,1:2,2,2,2,0:2)
          integer :: j, h12,h56,h3
          real(dp) :: facLO

          zgam_scheme = schemeMSBAR

          facLO = aveqq*esq**3*xn

          ! call with maxLoops=2
          if (decayChannel() == decayQuarks) then
            call zgamma_amps_vdecay(p, 2, ampsZ_du)
            facLO = facLO * cA
          else
            call zgamma_amps(p,4,3,1,5,2, 2, ampsZ)
          endif

          msq(:,:,:) = 0

          if (decayChannel() == decayQuarks) then
            do j=1,2
              ! fill LO
              msq(-j,j,0) = 
     &              + n_light_down*sum(abs(ampsZ_du(-j,j,1,:,:,:,0))**2)
     &              + n_light_up*sum(abs(ampsZ_du(-j,j,2,:,:,:,0))**2)
              msq(j,-j,0) = 
     &              + n_light_down*sum(abs(ampsZ_du(j,-j,1,:,:,:,0))**2)
     &              + n_light_up*sum(abs(ampsZ_du(j,-j,2,:,:,:,0))**2)

              ! fill NLO
              msq(-j,j,1) = 
     &         + n_light_down*sum(2*real(conjg(ampsZ_du(-j,j,1,:,:,:,0))*ampsZ_du(-j,j,1,:,:,:,1)))
     &         + n_light_up*sum(2*real(conjg(ampsZ_du(-j,j,2,:,:,:,0))*ampsZ_du(-j,j,2,:,:,:,1)))
              msq(j,-j,1) = 
     &         + n_light_down*sum(2*real(conjg(ampsZ_du(j,-j,1,:,:,:,0))*ampsZ_du(j,-j,1,:,:,:,1)))
     &         + n_light_up*sum(2*real(conjg(ampsZ_du(j,-j,2,:,:,:,0))*ampsZ_du(j,-j,2,:,:,:,1)))

              ! fill NNLO
              msq(-j,j,2) =
     &          + n_light_down*sum(2*real(
     &              conjg(ampsZ_du(-j,j,1,:,:,:,0))*ampsZ_du(-j,j,1,:,:,:,2))
     &          + n_light_down*abs(ampsZ_du(-j,j,1,:,:,:,1))**2)
     &          + n_light_up*sum(2*real(
     &              conjg(ampsZ_du(-j,j,2,:,:,:,0))*ampsZ_du(-j,j,2,:,:,:,2))
     &          + n_light_up*abs(ampsZ_du(-j,j,2,:,:,:,1))**2)

              msq(j,-j,2) = 
     &          + n_light_down*sum(2*real(
     &              conjg(ampsZ_du(j,-j,1,:,:,:,0))*ampsZ_du(j,-j,1,:,:,:,2))
     &          + n_light_down*abs(ampsZ_du(j,-j,1,:,:,:,1))**2)
     &          + n_light_up*sum(2*real(
     &              conjg(ampsZ_du(j,-j,2,:,:,:,0))*ampsZ_du(j,-j,2,:,:,:,2))
     &          + n_light_up*abs(ampsZ_du(j,-j,2,:,:,:,1))**2)
            enddo
          else
            do j=1,2
              ! fill LO
              msq(-j,j,0) = sum(abs(ampsZ(-j,j,:,:,:,0))**2)
              msq(j,-j,0) = sum(abs(ampsZ(j,-j,:,:,:,0))**2)

              ! fill NLO
              msq(-j,j,1) = sum(2*real(
     &                conjg(ampsZ(-j,j,:,:,:,0))*ampsZ(-j,j,:,:,:,1)))
              msq(j,-j,1) = sum(2*real(
     &                conjg(ampsZ(j,-j,:,:,:,0))*ampsZ(j,-j,:,:,:,1)))

              ! fill NNLO
              msq(-j,j,2) = sum(2*real(
     &       conjg(ampsZ(-j,j,:,:,:,0))*ampsZ(-j,j,:,:,:,2))
     &          + abs(ampsZ(-j,j,:,:,:,1))**2)

              msq(j,-j,2) = sum(2*real(
     &       conjg(ampsZ(j,-j,:,:,:,0))*ampsZ(j,-j,:,:,:,2))
     &          + abs(ampsZ(j,-j,:,:,:,1))**2)
            enddo
          endif

          msq(-3,3,:) = msq(-1,1,:)
          msq(-4,4,:) = msq(-2,2,:)
          msq(-5,5,:) = msq(-1,1,:)

          msq(3,-3,:) = msq(1,-1,:)
          msq(4,-4,:) = msq(2,-2,:)
          msq(5,-5,:) = msq(1,-1,:)

          msq(:,:,0) = msq(:,:,0) * facLO
          msq(:,:,1) = msq(:,:,1) * facLO * ason2pi
          msq(:,:,2) = msq(:,:,2) * facLO * ason2pi**2

      end subroutine

      subroutine qqb_zgam_v_new(p,msq,msq0)
        use VVconfig_m, only: schemeMCFM, zgam_scheme, decayChannel,
     &                    decayQuarks
          implicit none
          include 'types.f'
          include 'mxpart.f'
          include 'nf.f'
          include 'scheme.f'
          include 'constants.f'
          include 'ewcouple.f'
          include 'qcdcouple.f'

          real(dp), intent(in) :: p(mxpart,4)
          real(dp), intent(out) :: msq(-nf:nf,-nf:nf)
          real(dp), intent(out) :: msq0(-nf:nf,-nf:nf)

          complex(dp) :: ampsZ(-nf:nf,-nf:nf,2,2,2,0:2)
          complex(dp) :: ampsZ_du(-nf:nf,-nf:nf,1:2,2,2,2,0:2)
          integer :: j, h12,h56,h3
          real(dp) :: fac

          zgam_scheme = schemeMCFM
          scheme = 'dred'

          fac = aveqq*esq**3*cA

          !call with maxLoops = 1
          if (decayChannel() == decayQuarks) then
            call zgamma_amps_vdecay(p, 1, ampsZ_du)
            fac = fac * cA
          else
            call zgamma_amps(p,4,3,1,5,2, 1,ampsZ)
          endif

          msq(:,:) = 0._dp
          msq0(:,:) = 0._dp

          if (decayChannel() == decayQuarks) then
            do j=1,2
              msq(-j,j) = 
     &          + n_light_down*2._dp*sum(real(conjg(ampsZ_du(-j,j,1,:,:,:,0))*ampsZ_du(-j,j,1,:,:,:,1)))
     &          + n_light_up*2._dp*sum(real(conjg(ampsZ_du(-j,j,2,:,:,:,0))*ampsZ_du(-j,j,2,:,:,:,1)))
              msq(j,-j) = 
     &          + n_light_down*2._dp*sum(real(conjg(ampsZ_du(j,-j,1,:,:,:,0))*ampsZ_du(j,-j,1,:,:,:,1)))
     &          + n_light_up*2._dp*sum(real(conjg(ampsZ_du(j,-j,2,:,:,:,0))*ampsZ_du(j,-j,2,:,:,:,1)))
              msq0(-j,j) = 
     &          + n_light_down*sum(abs(ampsZ_du(-j,j,1,:,:,:,0))**2)
     &          + n_light_up*sum(abs(ampsZ_du(-j,j,2,:,:,:,0))**2)
              msq0(j,-j) = 
     &          + n_light_down*sum(abs(ampsZ_du(j,-j,1,:,:,:,0))**2)
     &          + n_light_up*sum(abs(ampsZ_du(j,-j,2,:,:,:,0))**2)
            enddo
          else
            do j=1,2
              msq(-j,j) = 2._dp*sum(real(conjg(ampsZ(-j,j,:,:,:,0))*ampsZ(-j,j,:,:,:,1)))
              msq(j,-j) = 2._dp*sum(real(conjg(ampsZ(j,-j,:,:,:,0))*ampsZ(j,-j,:,:,:,1)))
              msq0(-j,j) = sum(abs(ampsZ(-j,j,:,:,:,0))**2)
              msq0(j,-j) = sum(abs(ampsZ(j,-j,:,:,:,0))**2)
            enddo
          endif

          msq(-3,3) = msq(-1,1)
          msq(-4,4) = msq(-2,2)
          msq(-5,5) = msq(-1,1)

          msq(3,-3) = msq(1,-1)
          msq(4,-4) = msq(2,-2)
          msq(5,-5) = msq(1,-1)

          msq0(-3,3) = msq0(-1,1)
          msq0(-4,4) = msq0(-2,2)
          msq0(-5,5) = msq0(-1,1)

          msq0(3,-3) = msq0(1,-1)
          msq0(4,-4) = msq0(2,-2)
          msq0(5,-5) = msq0(1,-1)

          msq(:,:) = msq(:,:) * fac * ason2pi
          msq0(:,:) = msq0(:,:) * fac

      end subroutine

      subroutine qqb_zgam_new(p,msq)
            use VVconfig_m, only: decayChannel, decayQuarks
          implicit none
          include 'types.f'
          include 'mxpart.f'
          include 'nf.f'
          include 'scheme.f'
          include 'constants.f'
          include 'ewcouple.f'
          include 'ewcharge.f'
          include 'zcouple.f'

          real(dp), intent(in) :: p(mxpart,4)
          real(dp), intent(out) :: msq(-nf:nf,-nf:nf)

          ! results for lepton decay
          complex(dp) :: ampsZ(-nf:nf,-nf:nf,2,2,2,0:2)
          ! for quark decay
          complex(dp) :: ampsZ_du(-nf:nf,-nf:nf,1:2,2,2,2,0:2)
          integer :: j, h12,h56,h3
          real(dp) :: fac

          fac = aveqq*esq**3*cA

          ! call with maxLoops = 0
          if (decayChannel() == decayQuarks) then
            call zgamma_amps_vdecay(p, 0, ampsZ_du)
            fac = fac * cA
          else
            call zgamma_amps(p,4,3,1,5,2, 0, ampsZ)
          endif


c         !! DEBUG CHECKS
c         ! Both routines, zgamma_amps and zgamma_amps_vdecay must at
c         ! tree-level return the same results if the couplings q1,r1
c         ! and l1 which are used in zgamma_amps are adjusted accordingly
c         q1=Q(2)
c         r1=r(2)
c         l1=l(2)
c         call zgamma_amps_vdecay(p, 0, ampsZ_du)
c         call zgamma_amps(p,4,3,1,5,2, 0, ampsZ)
c         write (*,*) "qbar q, up-type decay"
c         do h12=1,2; do h56=1,2; do h3=1,2
c         write (*,*) h12,h56,h3, ampsZ_du(-1,1, 2,h12,h56,h3,0)/ampsZ(-1,1,h12,h56,h3,0)
c         enddo; enddo; enddo
c         write (*,*) "q qbar, up-type decay"
c         do h12=1,2; do h56=1,2; do h3=1,2
c         write (*,*) h12,h56,h3, ampsZ_du(1,-1, 2,h12,h56,h3,0)/ampsZ(1,-1,h12,h56,h3,0)
c         enddo; enddo; enddo

c         q1=Q(1)
c         r1=r(1)
c         l1=l(1)
c         call zgamma_amps_vdecay(p, 0, ampsZ_du)
c         call zgamma_amps(p,4,3,1,5,2, 0, ampsZ)
c         write (*,*) "qbar q, down-type decay"
c         do h12=1,2; do h56=1,2; do h3=1,2
c         write (*,*) h12,h56,h3, ampsZ_du(-1,1, 1,h12,h56,h3,0)/ampsZ(-1,1,h12,h56,h3,0)
c         enddo; enddo; enddo
c         write (*,*) "q qbar, down-type decay"
c         do h12=1,2; do h56=1,2; do h3=1,2
c         write (*,*) h12,h56,h3, ampsZ_du(1,-1, 1,h12,h56,h3,0)/ampsZ(1,-1,h12,h56,h3,0)
c         enddo; enddo; enddo

c         pause

c         !!! END DEBUG CHECKS

          msq(:,:) = 0

          if (decayChannel() == decayQuarks) then
            do j=1,2
              msq(-j,j) =
     &                + n_light_down*sum(abs(ampsZ_du(-j,j,1,:,:,:,0))**2)
     &                + n_light_up*sum(abs(ampsZ_du(-j,j,2,:,:,:,0))**2)
              msq(j,-j) =
     &                + n_light_down*sum(abs(ampsZ_du(j,-j,1,:,:,:,0))**2)
     &                + n_light_up*sum(abs(ampsZ_du(j,-j,2,:,:,:,0))**2)
            enddo
          else
            do j=1,2
              msq(-j,j) = sum(abs(ampsZ(-j,j,:,:,:,0))**2)
              msq(j,-j) = sum(abs(ampsZ(j,-j,:,:,:,0))**2)
            enddo
          endif


          msq(-3,3) = msq(-1,1)
          msq(-4,4) = msq(-2,2)
          msq(-5,5) = msq(-1,1)

          msq(3,-3) = msq(1,-1)
          msq(4,-4) = msq(2,-2)
          msq(5,-5) = msq(1,-1)

          msq(:,:) = msq(:,:) * fac

      end subroutine

      subroutine zgamma_amps(p,i5,i6,i1,i3,i2,maxLoops,ampsZ)
          use VVconfig_m
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
      include 'anomcoup.f' ! logical :: anomtgc
      include 'zerowidth.f'
      include 'zcouple.f'

      real(dp), intent(in) :: p(mxpart,4)
      integer, intent(in) :: i5,i6,i1,i3,i2
      integer, intent(in) :: maxLoops
      complex(dp), intent(out) :: ampsZ(-nf:nf,-nf:nf,2,2,2,0:2)
      complex(dp) :: ampsZf(-nf:nf,-nf:nf,2,2,2,0:2)
      complex(dp) :: ampsAnom(-nf:nf,-nf:nf,2,2,2,0:2)
      real(dp) :: mp(mxpart,4)

      complex(dp) :: za(mxpart,mxpart), zb(mxpart,mxpart)

      complex(dp) :: amp_rr_Z(0:2), amp_rr_Gamma(0:2)
      integer (kind(decayElAntiEl)) :: vDecay

      ! anomalous coupling contributions
      complex(dp) :: amp_ZZ(0:2), amp_aZ(0:2), amp_Za(0:2)

      real(dp) :: dotvec
      real(dp) :: s123
      real(dp) :: s12, s13, s23, s56

      complex(dp) :: c_alphai(0:2,4)
      complex(dp) :: c_betai(0:2,4)
      complex(dp) :: c_gammai(0:2,4)

      complex(dp) :: qbqamp(2,2,2,2)
      complex(dp) :: qbqamp_nlo(2,2,2,2)
      integer j,h12,h34
      complex(dp) :: DZ, DGamma
      complex(dp) :: DZf, DGammaf

      mp(:,:) = p(:,:)
      mp(1,:) = -p(i1,:) ! 1
      mp(2,:) = -p(i2,:) ! 2
      mp(3,:) = -p(i3,:) ! 5
      mp(5,:) = p(i5,:) ! 4
      mp(6,:) = p(i6,:) ! 3

      call spinoru_s(6,mp,za,zb)

      ampsZ(:,:,:,:,:,:) = 0._dp 
      ampsZf(:,:,:,:,:,:) = 0._dp

      s56 = dotvec(mp(5,:)+mp(6,:),mp(5,:)+mp(6,:))
      s12 = dotvec(mp(1,:)+mp(2,:),mp(1,:)+mp(2,:))

      DZ = s56 - zmass**2 + im*zwidth*zmass
      DGamma = s56

      DZf = s12 - zmass**2 + im*zwidth*zmass
      DGammaf = s12

      vDecay = decayChannel()

      ! to make comparison with qqb_zgam easier I use the same global
      ! and relative signs for the amplitudes with different helicities
      ! here whether or not they are strictly consistent among each
      ! other

      !first block, 56132 ordering
      call omega_wzgamma(mp,5,6,1,3,2,maxLoops,c_alphai,c_betai,c_gammai)
      do j=1,2

      call zgamma_amp_rr(5,6,1,3,2,za,zb,c_alphai,c_betai,c_gammai,
     &                      j,helRight,bosZ,amp_rr_Z)
      call zgamma_amp_rr(5,6,1,3,2,za,zb,c_alphai,c_betai,c_gammai,
     &                      j,helRight,bosGamma,amp_rr_Gamma)
      ampsZ(-j,j, 2,2,2,:) =
     &    amp_rr_Z(:) * r(j) * r1 / DZ
     &   +amp_rr_Gamma(:) * Q(j) * q1 / DGamma

      if (vDecay == decayElAntiEl .or. vDecay == decayQuarks) then
      call zgamma_ampf_rr(5,6,1,3,2,za,zb,amp_rr_Z)
      ampsZf(-j,j, 2,2,2,:) = 
     &    amp_rr_Z(:) * r(j) * r1 / DZf
     &   +amp_rr_Z(:) * Q(j) * q1 / DGammaf
      ampsZf(j,-j, 1,2,2,:) = 
     &    amp_rr_Z(:) * l(j) * r1 / DZf
     &   +amp_rr_Z(:) * Q(j) * q1 / DGammaf
      endif

      call zgamma_amp_rr(5,6,1,3,2,zb,za,c_alphai,c_betai,c_gammai,
     &                      j,helLeft,bosZ,amp_rr_Z)
      call zgamma_amp_rr(5,6,1,3,2,zb,za,c_alphai,c_betai,c_gammai,
     &                      j,helLeft,bosGamma,amp_rr_Gamma)
      ampsZ(-j,j, 1,1,1,:) =
     &    -amp_rr_Z(:) * l(j) * l1 / DZ
     &    -amp_rr_Gamma(:) * Q(j) * q1 / DGamma

      if (vDecay == decayElAntiEl .or. vDecay == decayQuarks) then
      call zgamma_ampf_rr(5,6,1,3,2,zb,za,amp_rr_Z)
      ampsZf(-j,j, 1,1,1,:) =
     &    -amp_rr_Z(:) * l(j) * l1 / DZf
     &    -amp_rr_Z(:) * Q(j) * q1 / DGammaf
      ampsZf(j,-j, 2,1,1,:) =
     &    -amp_rr_Z(:) * r(j) * l1 / DZf
     &    -amp_rr_Z(:) * Q(j) * q1 / DGammaf
      endif

      call zgamma_amp_rr(5,6,1,3,2,za,zb,c_alphai,c_betai,c_gammai,
     &                      -j,helLeft,bosZ,amp_rr_Z)
      call zgamma_amp_rr(5,6,1,3,2,za,zb,c_alphai,c_betai,c_gammai,
     &                      -j,helLeft,bosGamma,amp_rr_Gamma)
      ampsZ(j,-j, 1,2,2,:) =
     &     amp_rr_Z(:) * l(j) * r1 / DZ
     &    +amp_rr_Gamma(:) * Q(j) * q1 / DGamma

      call zgamma_amp_rr(5,6,1,3,2,zb,za,c_alphai,c_betai,c_gammai,
     &                      -j,helRight,bosZ,amp_rr_Z)
      call zgamma_amp_rr(5,6,1,3,2,zb,za,c_alphai,c_betai,c_gammai,
     &                      -j,helRight,bosGamma,amp_rr_Gamma)
      ampsZ(j,-j, 2,1,1,:) =
     &    -amp_rr_Z(:) * r(j) * l1 / DZ
     &    -amp_rr_Gamma(:) * Q(j) * q1 / DGamma

      enddo

      !second block, 65231 ordering
      call omega_wzgamma(mp,6,5,2,3,1,maxLoops,c_alphai,c_betai,c_gammai)
      do j=1,2

      call zgamma_amp_rr(6,5,2,3,1,zb,za,c_alphai,c_betai,c_gammai,
     &                      j,helRight,bosZ,amp_rr_Z)
      call zgamma_amp_rr(6,5,2,3,1,zb,za,c_alphai,c_betai,c_gammai,
     &                      j,helRight,bosGamma,amp_rr_Gamma)
      ampsZ(-j,j, 2,2,1,:) =
     &     amp_rr_Z(:) * r(j) * r1 / DZ
     &    +amp_rr_Gamma(:) * Q(j) * q1 / DGamma

      if (vDecay == decayElAntiEl .or. vDecay == decayQuarks) then
      call zgamma_ampf_rr(6,5,2,3,1,zb,za,amp_rr_Z)
      ampsZf(-j,j, 2,2,1,:) =
     &    +amp_rr_Z(:) * r(j) * r1 / DZf
     &    +amp_rr_Z(:) * Q(j) * q1 / DGammaf
      ampsZf(j,-j, 1,2,1,:) =
     &     amp_rr_Z(:) * l(j) * r1 / DZf
     &    +amp_rr_Z(:) * Q(j) * q1 / DGammaf
      endif

      call zgamma_amp_rr(6,5,2,3,1,za,zb,c_alphai,c_betai,c_gammai,
     &                      j,helLeft,bosZ,amp_rr_Z)
      call zgamma_amp_rr(6,5,2,3,1,za,zb,c_alphai,c_betai,c_gammai,
     &                      j,helLeft,bosGamma,amp_rr_Gamma)
      ampsZ(-j,j, 1,1,2,:) =
     &    -amp_rr_Z(:) * l(j) * l1 / DZ
     &    -amp_rr_Gamma(:) * Q(j) * q1 / DGamma
  
      if (vDecay == decayElAntiEl .or. vDecay == decayQuarks) then
      call zgamma_ampf_rr(6,5,2,3,1,za,zb,amp_rr_Z)
      ampsZf(-j,j, 1,1,2,:) =
     &    -amp_rr_Z(:) * l(j) * l1 / DZf
     &    -amp_rr_Z(:) * Q(j) * q1 / DGammaf
      ampsZf(j,-j, 2,1,2,:) =
     &    -amp_rr_Z(:) * r(j) * l1 / DZf
     &    -amp_rr_Z(:) * Q(j) * q1 / DGammaf
      endif

      call zgamma_amp_rr(6,5,2,3,1,zb,za,c_alphai,c_betai,c_gammai,
     &                      -j,helLeft,bosZ,amp_rr_Z)
      call zgamma_amp_rr(6,5,2,3,1,zb,za,c_alphai,c_betai,c_gammai,
     &                      -j,helLeft,bosGamma,amp_rr_Gamma)
      ampsZ(j,-j, 1,2,1,:) =
     &     amp_rr_Z(:) * l(j) * r1 / DZ
     &    +amp_rr_Gamma(:) * Q(j) * q1 / DGamma

      call zgamma_amp_rr(6,5,2,3,1,za,zb,c_alphai,c_betai,c_gammai,
     &                      -j,helRight,bosZ,amp_rr_Z)
      call zgamma_amp_rr(6,5,2,3,1,za,zb,c_alphai,c_betai,c_gammai,
     &                      -j,helRight,bosGamma,amp_rr_Gamma)
      ampsZ(j,-j, 2,1,2,:) =
     &    -amp_rr_Z(:) * r(j) * l1 / DZ
     &    -amp_rr_Gamma(:) * Q(j) * q1 / DGamma

      enddo

      !third block, 56231 ordering
      call omega_wzgamma(mp,5,6,2,3,1,maxLoops,c_alphai,c_betai,c_gammai)
      do j=1,2

      call zgamma_amp_rr(5,6,2,3,1,za,zb,c_alphai,c_betai,c_gammai,
     &                      j,helLeft,bosZ,amp_rr_Z)
      call zgamma_amp_rr(5,6,2,3,1,za,zb,c_alphai,c_betai,c_gammai,
     &                      j,helLeft,bosGamma,amp_rr_Gamma)
      ampsZ(-j,j, 1,2,2,:) =
     &    -amp_rr_Z(:) * l(j) * r1 / DZ
     &    -amp_rr_Gamma(:) * Q(j) * q1 / DGamma

      if (vDecay == decayElAntiEl .or. vDecay == decayQuarks) then
      call zgamma_ampf_rr(5,6,2,3,1,za,zb,amp_rr_Z)
      ampsZf(-j,j, 1,2,2,:) =
     &     amp_rr_Z(:) * l(j) * r1 / DZf
     &    +amp_rr_Z(:) * Q(j) * q1 / DGammaf
      ampsZf(j,-j, 2,2,2,:) =
     &     amp_rr_Z(:) * r(j) * r1 / DZf
     &    +amp_rr_Z(:) * Q(j) * q1 / DGammaf
      endif

      call zgamma_amp_rr(5,6,2,3,1,zb,za,c_alphai,c_betai,c_gammai,
     &                      j,helRight,bosZ,amp_rr_Z)
      call zgamma_amp_rr(5,6,2,3,1,zb,za,c_alphai,c_betai,c_gammai,
     &                      j,helRight,bosGamma,amp_rr_Gamma)
      ampsZ(-j,j, 2,1,1,:) =
     &     amp_rr_Z(:) * r(j) * l1 / DZ
     &    +amp_rr_Gamma(:) * Q(j) * q1 / DGamma

      if (vDecay == decayElAntiEl .or. vDecay == decayQuarks) then
      call zgamma_ampf_rr(5,6,2,3,1,zb,za,amp_rr_Z)
      ampsZf(-j,j, 2,1,1,:) =
     &    -amp_rr_Z(:) * r(j) * l1 / DZf
     &    -amp_rr_Z(:) * Q(j) * q1 / DGammaf
      ampsZf(j,-j, 1,1,1,:) =
     &    -amp_rr_Z(:) * l(j) * l1 / DZf
     &    -amp_rr_Z(:) * Q(j) * q1 / DGammaf
      endif

      call zgamma_amp_rr(5,6,2,3,1,za,zb,c_alphai,c_betai,c_gammai,
     &                      -j,helRight,bosZ,amp_rr_Z)
      call zgamma_amp_rr(5,6,2,3,1,za,zb,c_alphai,c_betai,c_gammai,
     &                      -j,helRight,bosGamma,amp_rr_Gamma)
      ampsZ(j,-j, 2,2,2,:) =
     &    -amp_rr_Z(:) * r(j) * r1 / DZ
     &    -amp_rr_Gamma(:) * Q(j) * q1 / DGamma

      call zgamma_amp_rr(5,6,2,3,1,zb,za,c_alphai,c_betai,c_gammai,
     &                      -j,helLeft,bosZ,amp_rr_Z)
      call zgamma_amp_rr(5,6,2,3,1,zb,za,c_alphai,c_betai,c_gammai,
     &                      -j,helLeft,bosGamma,amp_rr_Gamma)
      ampsZ(j,-j, 1,1,1,:) =
     &     amp_rr_Z(:) * l(j) * l1 / DZ
     &    +amp_rr_Gamma(:) * Q(j) * q1 / DGamma


      enddo

      !fourth and last block, 65132 ordering
      call omega_wzgamma(mp,6,5,1,3,2,maxLoops,c_alphai,c_betai,c_gammai)
      do j=1,2

      call zgamma_amp_rr(6,5,1,3,2,zb,za,c_alphai,c_betai,c_gammai,
     &                      j,helLeft,bosZ,amp_rr_Z)
      call zgamma_amp_rr(6,5,1,3,2,zb,za,c_alphai,c_betai,c_gammai,
     &                      j,helLeft,bosGamma,amp_rr_Gamma)
      ampsZ(-j,j, 1,2,1,:) =
     &    -amp_rr_Z(:) * l(j) * r1 / DZ
     &    -amp_rr_Gamma(:) * Q(j) * q1 / DGamma

      if (vDecay == decayElAntiEl .or. vDecay == decayQuarks) then
      call zgamma_ampf_rr(6,5,1,3,2,zb,za,amp_rr_Z)
      ampsZf(-j,j, 1,2,1,:) =
     &     +amp_rr_Z(:) * l(j) * r1 / DZf
     &    +amp_rr_Z(:) * Q(j) * q1 / DGammaf
      ampsZf(j,-j, 2,2,1,:) =
     &     amp_rr_Z(:) * r(j) * r1 / DZf
     &    +amp_rr_Z(:) * Q(j) * q1 / DGammaf
      endif

      call zgamma_amp_rr(6,5,1,3,2,za,zb,c_alphai,c_betai,c_gammai,
     &                      j,helRight,bosZ,amp_rr_Z)
      call zgamma_amp_rr(6,5,1,3,2,za,zb,c_alphai,c_betai,c_gammai,
     &                      j,helRight,bosGamma,amp_rr_Gamma)
      ampsZ(-j,j, 2,1,2,:) =
     &     amp_rr_Z(:) * r(j) * l1 / DZ
     &    +amp_rr_Gamma(:) * Q(j) * q1 / DGamma

      if (vDecay == decayElAntiEl .or. vDecay == decayQuarks) then
      call zgamma_ampf_rr(6,5,1,3,2,za,zb,amp_rr_Z)
      ampsZf(-j,j, 2,1,2,:) =
     &     -amp_rr_Z(:) * r(j) * l1 / DZf
     &    -amp_rr_Z(:) * Q(j) * q1 / DGammaf
      ampsZf(j,-j, 1,1,2,:) =
     &    -amp_rr_Z(:) * l(j) * l1 / DZf
     &    -amp_rr_Z(:) * Q(j) * q1 / DGammaf
      endif

      call zgamma_amp_rr(6,5,1,3,2,zb,za,c_alphai,c_betai,c_gammai,
     &                      -j,helRight,bosZ,amp_rr_Z)
      call zgamma_amp_rr(6,5,1,3,2,zb,za,c_alphai,c_betai,c_gammai,
     &                      -j,helRight,bosGamma,amp_rr_Gamma)
      ampsZ(j,-j, 2,2,1,:) =
     &    -amp_rr_Z(:) * r(j) * r1 / DZ
     &    -amp_rr_Gamma(:) * Q(j) * q1 / DGamma

      call zgamma_amp_rr(6,5,1,3,2,za,zb,c_alphai,c_betai,c_gammai,
     &                      -j,helLeft,bosZ,amp_rr_Z)
      call zgamma_amp_rr(6,5,1,3,2,za,zb,c_alphai,c_betai,c_gammai,
     &                      -j,helLeft,bosGamma,amp_rr_Gamma)
      ampsZ(j,-j, 1,1,2,:) =
     &     amp_rr_Z(:) * l(j) * l1 / DZ
     &    +amp_rr_Gamma(:) * Q(j) * q1 / DGamma

      enddo

c =========================== 
c === anomalous couplings === 
c =========================== 

      ampsAnom = 0._dp

      if (anomtgc) then

      call zgamma_amp_anom(5,6,1,3,2,za,zb,.false.,amp_ZZ,amp_aZ,amp_Za)
      do j=1,2
        ampsAnom(-j,j, 2,2,2,:) = 
     &       amp_ZZ(:) / DZ
     &          * r(j)
     &          * r1
     &      +amp_Za(:) / DZf
     &          * r(j)
     &          * q1
     &      +amp_aZ(:) / DZ
     &          * Q(j)
     &          * r1
        ampsAnom(j,-j, 1,2,2,:) = 
     &        amp_aZ(:) / DZ
     &          * (-Q(j))
     &          * r1
     &       +amp_ZZ(:) / DZ
     &          * (-l(j))
     &          * r1
     &       +amp_Za(:) / DZf
     &          * (-l(j))
     &          * q1
      enddo

      call zgamma_amp_anom(6,5,2,3,1,zb,za,.true.,amp_ZZ,amp_aZ,amp_Za)
      do j=1,2
        ampsAnom(-j,j, 2,2,1,:) = 
     &        amp_ZZ(:) / DZ
     &          * r(j)
     &          * r1
     &       +amp_aZ(:) / DZ
     &          * Q(j)
     &          * r1
     &       +amp_Za(:) / DZf
     &          * r(j)
     &          * q1

        ampsAnom(j,-j, 1,2,1,:) = 
     &        amp_aZ(:) / DZ
     &          * (-Q(j))
     &          * r1
     &       +amp_ZZ(:) / DZ
     &          * (-l(j))
     &          * r1
     &       +amp_Za(:) / DZf
     &          * (-l(j))
     &          * q1

      enddo

      call zgamma_amp_anom(6,5,1,3,2,za,zb,.false.,amp_ZZ,amp_aZ,amp_Za)
      do j=1,2
        ampsAnom(-j,j, 2,1,2,:) = 
     &        amp_ZZ(:) / DZ
     &          * r(j)
     &          * l1
     &       +amp_aZ(:) / DZ
     &          * Q(j)
     &          * l1
     &       +amp_Za(:) / DZf
     &          * r(j)
     &          * q1

        ampsAnom(j,-j, 1,1,2,:) = 
     &        amp_ZZ(:) / DZ
     &          * (-l(j))
     &          * l1
     &       +amp_aZ(:) / DZ
     &          * (-Q(j))
     &          * l1
     &       +amp_Za(:) / DZf
     &          * (-l(j))
     &          * q1
      enddo

      call zgamma_amp_anom(5,6,2,3,1,zb,za,.true.,amp_ZZ,amp_aZ,amp_Za)
      do j=1,2
        ampsAnom(-j,j, 2,1,1,:) = 
     &        amp_ZZ(:) / DZ
     &          * r(j)
     &          * l1
     &       +amp_aZ(:) / DZ
     &          * Q(j)
     &          * l1
     &       +amp_Za(:) / DZf
     &          * r(j)
     &          * q1

        ampsAnom(j,-j, 1,1,1,:) = 
     &        amp_ZZ(:) / DZ
     &          * (-l(j))
     &          * l1
     &       +amp_aZ(:) / DZ
     &          * (-Q(j))
     &          * l1
     &       +amp_Za(:) / DZf
     &          * (-l(j))
     &          * q1
      enddo

      call zgamma_amp_anom(5,6,2,3,1,za,zb,.false.,amp_ZZ,amp_aZ,amp_Za)
      do j=1,2
        ampsAnom(-j,j, 1,2,2,:) = 
     &       -amp_ZZ(:) / DZ
     &          * l(j)
     &          * r1
     &       -amp_aZ(:) / DZ
     &          * Q(j)
     &          * r1
     &       -amp_Za(:) / DZf
     &          * l(j)
     &          * q1

        ampsAnom(j,-j, 2,2,2,:) = 
     &       -amp_ZZ(:) / DZ
     &          * (-r(j))
     &          * r1
     &       -amp_aZ(:) / DZ
     &          * (-Q(j))
     &          * r1
     &       -amp_Za(:) / DZf
     &          * (-r(j))
     &          * q1
      enddo

      call zgamma_amp_anom(6,5,1,3,2,zb,za,.true.,amp_ZZ,amp_aZ,amp_Za)
      do j=1,2
        ampsAnom(-j,j, 1,2,1,:) = 
     &       -amp_ZZ(:) / DZ
     &          * l(j)
     &          * r1
     &       -amp_aZ(:) / DZ
     &          * Q(j)
     &          * r1
     &       -amp_Za(:) / DZf
     &          * l(j)
     &          * q1

        ampsAnom(j,-j, 2,2,1,:) = 
     &       -amp_ZZ(:) / DZ
     &          * (-r(j))
     &          * r1
     &       -amp_aZ(:) / DZ
     &          * (-Q(j))
     &          * r1
     &       -amp_Za(:) / DZf
     &          * (-r(j))
     &          * q1
      enddo

      call zgamma_amp_anom(6,5,2,3,1,za,zb,.false.,amp_ZZ,amp_aZ,amp_Za)
      do j=1,2
        ampsAnom(-j,j, 1,1,2,:) = 
     &       -amp_ZZ(:) / DZ
     &          * l(j)
     &          * l1
     &       -amp_aZ(:) / DZ
     &          * Q(j)
     &          * l1
     &       -amp_Za(:) / DZf
     &          * l(j)
     &          * q1

        ampsAnom(j,-j, 2,1,2,:) = 
     &       -amp_ZZ(:) / DZ
     &          * (-r(j))
     &          * l1
     &       -amp_aZ(:) / DZ
     &          * (-Q(j))
     &          * l1
     &       -amp_Za(:) / DZf
     &          * (-r(j))
     &          * q1
      enddo

      call zgamma_amp_anom(5,6,1,3,2,zb,za,.true.,amp_ZZ,amp_aZ,amp_Za)
      do j=1,2
        ampsAnom(-j,j, 1,1,1,:) = 
     &       -amp_ZZ(:) / DZ
     &          * l(j)
     &          * l1
     &       -amp_aZ(:) / DZ
     &          * Q(j)
     &          * l1
     &       -amp_Za(:) / DZf
     &          * l(j)
     &          * q1

        ampsAnom(j,-j, 2,1,1,:) = 
     &       -amp_ZZ(:) / DZ
     &          * (-r(j))
     &          * l1
     &       -amp_aZ(:) / DZ
     &          * (-Q(j))
     &          * l1
     &       -amp_Za(:) / DZf
     &          * (-r(j))
     &          * q1
      enddo

      endif

      if (zerowidth) ampsZf = 0._dp

      ampsZ = ampsZ + ampsZf + ampsAnom

      return

      ! ======================================= debug below

      j = 2
      print *, "quark flavor: j=", j

      if (1 == 2) then

      write (*,*) "tree level qbar q"
      write (*,*) "11", ampsZ(-j,j, 1,1,1,0), ampsZ(-j,j, 1,1,2,0)
      write (*,*) "12", ampsZ(-j,j, 1,2,1,0), ampsZ(-j,j, 1,2,2,0)
      write (*,*) "21", ampsZ(-j,j, 2,1,1,0), ampsZ(-j,j, 2,1,2,0)
      write (*,*) "22", ampsZ(-j,j, 2,2,1,0), ampsZ(-j,j, 2,2,2,0)

      write (*,*)
      write (*,*) "1-loop qbar q"
      write (*,*) "11", ampsZ(-j,j, 1,1,1,1), ampsZ(-j,j, 1,1,2,1)
      write (*,*) "12", ampsZ(-j,j, 1,2,1,1), ampsZ(-j,j, 1,2,2,1)
      write (*,*) "21", ampsZ(-j,j, 2,1,1,1), ampsZ(-j,j, 2,1,2,1)
      write (*,*) "22", ampsZ(-j,j, 2,2,1,1), ampsZ(-j,j, 2,2,2,1)

      write (*,*) "tree level q qbar"
      write (*,*) "11", ampsZ(j,-j, 1,1,1,0), ampsZ(j,-j, 1,1,2,0)
      write (*,*) "12", ampsZ(j,-j, 1,2,1,0), ampsZ(j,-j, 1,2,2,0)
      write (*,*) "21", ampsZ(j,-j, 2,1,1,0), ampsZ(j,-j, 2,1,2,0)
      write (*,*) "22", ampsZ(j,-j, 2,2,1,0), ampsZ(j,-j, 2,2,2,0)

      write (*,*)
      write (*,*) "1-loop q qbar"
      write (*,*) "11", ampsZ(j,-j, 1,1,1,1), ampsZ(j,-j, 1,1,2,1)
      write (*,*) "12", ampsZ(j,-j, 1,2,1,1), ampsZ(j,-j, 1,2,2,1)
      write (*,*) "21", ampsZ(j,-j, 2,1,1,1), ampsZ(j,-j, 2,1,2,1)
      write (*,*) "22", ampsZ(j,-j, 2,2,1,1), ampsZ(j,-j, 2,2,2,1)

      endif

      if (1==1) then

      write (*,*) "final"
      write (*,*) "11", ampsZf(-j, j, 1,1,1,0), ampsZf(-j, j, 1,1,2,0)
      write (*,*) "12", ampsZf(-j, j, 1,2,1,0), ampsZf(-j, j, 1,2,2,0)
      write (*,*) "21", ampsZf(-j, j, 2,1,1,0), ampsZf(-j, j, 2,1,2,0)
      write (*,*) "22", ampsZf(-j, j, 2,2,1,0), ampsZf(-j, j, 2,2,2,0)

      write (*,*) "final loop"
      write (*,*) "11", ampsZf(-j, j, 1,1,1,1), ampsZf(-j, j, 1,1,2,1)
      write (*,*) "12", ampsZf(-j, j, 1,2,1,1), ampsZf(-j, j, 1,2,2,1)
      write (*,*) "21", ampsZf(-j, j, 2,1,1,1), ampsZf(-j, j, 2,1,2,1)
      write (*,*) "22", ampsZf(-j, j, 2,2,1,1), ampsZf(-j, j, 2,2,2,1)

      write (*,*) "final exc"
      write (*,*) "11", ampsZf( j,-j, 1,1,1,0), ampsZf( j,-j, 1,1,2,0)
      write (*,*) "12", ampsZf( j,-j, 1,2,1,0), ampsZf( j,-j, 1,2,2,0)
      write (*,*) "21", ampsZf( j,-j, 2,1,1,0), ampsZf( j,-j, 2,1,2,0)
      write (*,*) "22", ampsZf( j,-j, 2,2,1,0), ampsZf( j,-j, 2,2,2,0)

      write (*,*) "final loop exc"
      write (*,*) "11", ampsZf( j,-j, 1,1,1,1), ampsZf( j,-j, 1,1,2,1)
      write (*,*) "12", ampsZf( j,-j, 1,2,1,1), ampsZf( j,-j, 1,2,2,1)
      write (*,*) "21", ampsZf( j,-j, 2,1,1,1), ampsZf( j,-j, 2,1,2,1)
      write (*,*) "22", ampsZf( j,-j, 2,2,1,1), ampsZf( j,-j, 2,2,2,1)

      endif


      s12 = dotvec(mp(1,:)+mp(2,:),mp(1,:)+mp(2,:))
      s13 = dotvec(mp(1,:)+mp(3,:),mp(1,:)+mp(3,:))
      s23 = dotvec(mp(2,:)+mp(3,:),mp(2,:)+mp(3,:))
      s123 = dotvec(mp(1,:)+mp(2,:)+mp(3,:),mp(1,:)+mp(2,:)+mp(3,:))

      call spinoru(5,p,za,zb)

      qbqamp(:,:,:,:) = czip
      qbqamp_nlo(:,:,:,:) = czip
      call zvirtamps(1,2,3,4,5,za,zb,qbqamp_nlo,qbqamp)

      qbqamp(:,:,:,:) = qbqamp(:,:,:,:) * 2._dp * sqrt(2._dp)
      qbqamp_nlo(:,:,:,:) = qbqamp_nlo(:,:,:,:) * 2._dp * sqrt(2._dp)

      write (*,*) "mcfm trees" 
        do h12=1,2
            do h34=1,2
                    write (*,*) "j,h12,h34,h5",h12,h34,
     &                    qbqamp(j,h12,h34,1) / ampsZ(-j,j,h12,h34,1,0),
     &                    qbqamp(j,h12,h34,2) / ampsZ(-j,j,h12,h34,2,0)
            enddo
        enddo

      write (*,*) "mcfm loops" 
        do h12=1,2
            do h34=1,2
                    write (*,*) "j,h12,h34,h5",h12,h34,
     &           qbqamp_nlo(j,h12,h34,1)*cf / ampsZ(-j,j,h12,h34,1,1),
     &           qbqamp_nlo(j,h12,h34,2)*cf / ampsZ(-j,j,h12,h34,2,1)
            enddo
        enddo

      qbqamp(:,:,:,:) = czip
      qbqamp_nlo(:,:,:,:) = czip
      call zvirtamps(2,1,3,4,5,za,zb,qbqamp_nlo,qbqamp)

      qbqamp(:,:,:,:) = qbqamp(:,:,:,:) * 2._dp * sqrt(2._dp)
      qbqamp_nlo(:,:,:,:) = qbqamp_nlo(:,:,:,:) * 2._dp * sqrt(2._dp)

      write (*,*) "mcfm trees exc" 
        do h12=1,2
            do h34=1,2
                    write (*,*) "j,h12,h34,h5",h12,h34,
     &                    qbqamp(j,h12,h34,1) / ampsZ(j,-j,h12,h34,1,0),
     &                    qbqamp(j,h12,h34,2) / ampsZ(j,-j,h12,h34,2,0)
            enddo
        enddo

      write (*,*) "mcfm loops exc" 
        do h12=1,2
            do h34=1,2
                    write (*,*) "j,h12,h34,h5",h12,h34,
     &          qbqamp_nlo(j,h12,h34,1)*cf / ampsZ(j,-j,h12,h34,1,1),
     &          qbqamp_nlo(j,h12,h34,2)*cf / ampsZ(j,-j,h12,h34,2,1)
            enddo
        enddo

      pause

      end subroutine

      subroutine zgamma_amp_qff(i5,i6,i1,i3,i2,za,zb,tree,amps)
        use VVconfig_m
        implicit none
        include 'types.f'
        include 'mxpart.f'
        include 'constants.f'
        include 'scale.f'
        include 'nf.f'

        integer, intent(in) :: i5,i6,i1,i3,i2
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: tree
        complex(dp), intent(out) :: amps(0:2)

        real(dp) :: s12
        complex(dp) :: Fq(0:2)

        s12 = real(za(i1,i2)*zb(i2,i1))
        call quark_formfactor(s12, Fq)

        amps(0) = tree * (-Fq(0) * sqrt(2._dp) * 2._dp)
        amps(1) = Fq(1) * amps(0)
        amps(2) = Fq(2) * amps(0)
 
        if (zgam_scheme == schemeMCFM) then
            amps(1) = amps(1) + cf*(-pisq/12 + 0.5)*amps(0)
        elseif (zgam_scheme == schemeMSBAR) then
            ! already prepared for MSBAR scheme
        else
            print *, "undefined scheme in zgamma_amp_qff"
            call exit(1)
        endif

      end subroutine

      subroutine zgamma_amp_anom(i5,i6,i1,i3,i2,za,zb,flip,
     &              amps_ZZ, amps_aZ, amps_Za)
        implicit none
        include 'types.f'
        include 'mxpart.f'
        include 'nf.f'
        include 'anomcoup.f'
        include 'masses.f'
        include 'constants.f'

        integer, intent(in) :: i5,i6,i1,i3,i2
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        logical, intent(in) :: flip
        complex(dp), intent(out) :: amps_ZZ(0:2), amps_aZ(0:2),
     &                              amps_Za(0:2)
        complex(dp) :: amps(0:2)
        complex(dp) :: tree_ZZ, tree_aZ, tree_Za
        complex(dp) :: cprop
        real(dp) :: hSign, s12, s56

        s12 = real(za(i1,i2)*zb(i2,i1))
        s56 = real(za(i5,i6)*zb(i6,i5))

        ! we need to flip the sign of h1tZ,h2tZ,h1tgam,h2tgam
        ! when za<>zb
        hSign = 1._dp
        if (flip) then
          hSign = -1._dp
        endif

        cprop = (s12-s56)/(s12-zmass**2 + im*zmass*zwidth )

        tree_ZZ = 1._dp/4._dp*cprop
     & *((hSign*im*hitZ(1)+hitZ(3))*2._dp*za(i2,i5)*zb(i1,i3)*zb(i6,i3)
     &  +(hSign*im*hitZ(2)+hitZ(4))*za(i2,i1)*za(i5,i3)*zb(i3,i6)*zb(i1,i3)**2)
        tree_aZ = 1._dp/4._dp
     & *((hSign*im*hitgam(1)+hitgam(3))*2._dp*za(i2,i5)*zb(i1,i3)*zb(i6,i3)
     &  +(hSign*im*hitgam(2)+hitgam(4))*za(i2,i1)*za(i5,i3)*zb(i3,i6)*zb(i1,i3)**2)
        tree_Za = 1._dp/4._dp
     & *((hSign*im*hitgam(1)+hitgam(3))*2._dp*za(i5,i2)*zb(i6,i3)*zb(i1,i3)
     &  -(hSign*im*hitgam(1)+hitgam(4))*za(i5,i6)*za(i2,i3)*zb(i3,i1)*zb(i6,i3)**2)

        call zgamma_amp_qff(i5,i6,i1,i3,i2,za,zb,(1._dp,0._dp),amps)

        amps_ZZ = tree_ZZ * amps
        amps_aZ = tree_aZ * amps
        amps_Za = tree_Za * amps

      end subroutine

      subroutine zgamma_ampf_RR(i5,i6,i1,i3,i2,za,zb,amps)
          use VVconfig_m
        implicit none
        include 'types.f'
        include 'mxpart.f'
        include 'nf.f'
        include 'zcouple.f'

        integer, intent(in) :: i5,i6,i1,i3,i2
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(out) :: amps(0:2)
        complex(dp) :: tree

        ! the coupling is in principle q1, but
        ! when switching to a quark decay by replacing q1,l1,r1, with an
        ! augmented color factor this would break
        tree = -q1*za(i2,i5)**2*zb(i2,i1)/(za(i5,i3)*za(i6,i3))

        call zgamma_amp_qff(i5,i6,i1,i3,i2,za,zb,tree,amps)

      end subroutine

      subroutine zgamma_amp_RR(i5,i6,i1,i3,i2,za,zb,
     &   c_alphai, c_betai, c_gammai, qf,helq, bos,amps)
        use VVconfig_m
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'scale.f'
      include 'constants.f'
      include 'ewcharge.f'
      include 'zcouple.f'

      integer, intent(in) :: i5,i6,i1,i3,i2
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
      complex(dp), intent(in) :: c_alphai(0:2,4)
      complex(dp), intent(in) :: c_betai(0:2,4)
      complex(dp), intent(in) :: c_gammai(0:2,4)
      integer, intent(in) :: qf
      integer (kind(helLeft)), intent(in) :: helq
      integer (kind(bosZ)), intent(in) :: bos
      complex(dp), intent(out) :: amps(0:2)

      complex(dp) :: c_alpha(0:2)
      complex(dp) :: c_beta(0:2)
      complex(dp) :: c_gamma(0:2)
      real(dp) :: s12, s123
      complex(dp) :: wzgamma_amp_RR, omega_I1
      real(dp) :: eq,vl,vr
      integer :: j
      real(dp) :: NFV

      eq = Q(qf)

      s12 = real(za(i1,i2)*zb(i2,i1))
      s123 = real(za(i1,i2)*zb(i2,i1) + za(i1,i3)*zb(i3,i1) +
     &                             za(i2,i3)*zb(i3,i2))

      NFV = 0
      if (bos == bosZ) then
        do j=1,nf
            NFV = NFV - Q(j)*(l(j) + r(j))
        enddo

        if (helQ == helLeft ) then
          vl = l(abs(qf))*sign(1._dp,-real(qf,dp))
          NFV = NFV/2/vl
        elseif (helQ == helRight) then
          vr = r(abs(qf))*sign(1._dp,-real(qf,dp))
          NFV = NFV/2/vr
        endif
      elseif (bos == bosGamma) then
        NFV = sum(Q(1:nf)**2)/eq
      endif

      c_alpha(:) = eq*(c_alphai(:,1) + c_alphai(:,2))
     &   + c_alphai(:,4) * NFV
      c_beta(:) = eq*(c_betai(:,1) + c_betai(:,2))
     &   + c_betai(:,4) * NFV
      c_gamma(:) = eq*(c_gammai(:,1) + c_gammai(:,2))
     &   + c_gammai(:,4) * NFV

      amps(0) = wzgamma_amp_RR(i5,i6,i1,i3,i2,za,zb,
     &             c_alpha(0),c_beta(0),c_gamma(0))

      amps(1) = wzgamma_amp_RR(i5,i6,i1,i3,i2,za,zb,
     &             c_alpha(1),c_beta(1),c_gamma(1))

      amps(2) = wzgamma_amp_RR(i5,i6,i1,i3,i2,za,zb,
     &             c_alpha(2),c_beta(2),c_gamma(2))


      ! IR scheme and scale handling
      if (zgam_scheme == schemeMCFM) then
          ! in case of schemeMCFM, we assume that only an NLO cross
          ! section is calculated, amps(2) is left untouched
          amps(1) = amps(1) + omega_I1(s12)*amps(0)
      elseif (zgam_scheme == schemeMSBAR) then
          ! perform scheme change to msbar at renormalization point s123
          call zgamma_catani_to_msbar(i5,i6,i1,i3,i2,za,zb,amps,s123)
          !MSBAR IR+UV evolution to scale musq from s123
          call evolve_msbar_IR_qq(i5,i6,i1,i3,i2,za,zb,amps,s123)
      else
          print *, "scheme unsupported"
          call exit(1)
      endif

      end subroutine

