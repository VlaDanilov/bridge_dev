      subroutine zgamma_amps_vdecay(p,maxLoops,ampsZ)
          use VVconfig_m
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'ewcharge.f'
      include 'masses.f' ! wmass, wwidth
      include 'anomcoup.f' ! anomtgc
      include 'constants.f' ! im
      include 'epinv.f'
      include 'epinv2.f'
      include 'scale.f'
      include 'zerowidth.f'
      include 'zcouple.f'

      real(dp), intent(in) :: p(mxpart,4)
      !integer, intent(in) :: i5,i6,i1,i3,i2
      integer, intent(in) :: maxLoops
      ! labels: initial partons; 1: down-type decay, 2: up-type decay
      !         helicities h12,h34,h5 (crossed), loop order
      complex(dp), intent(out) :: ampsZ(-nf:nf,-nf:nf,1:2,2,2,2,0:2)
      complex(dp) :: ampsZf(-nf:nf,-nf:nf,1:2,2,2,2,0:2)

      real(dp) :: mp(mxpart,4)

      complex(dp) :: za(mxpart,mxpart), zb(mxpart,mxpart)

      complex(dp) :: amp_rr_Z(0:2), amp_rr_Gamma(0:2)
      integer (kind(decayElAntiEl)) :: vDecay
      complex(dp) :: wzgamma_amp_rr

      real(dp) :: dotvec
      real(dp) :: s123, s12, s56
      integer :: j,flavor

      complex(dp) :: c_alphai(0:2,4)
      complex(dp) :: c_betai(0:2,4)
      complex(dp) :: c_gammai(0:2,4)

      complex(dp) :: DZ, DGamma
      complex(dp) :: DZf, DGammaf

      vDecay = decayChannel()
      if (vDecay /= decayQuarks) call abort
      if (anomtgc) call abort

      ! these are crossed kinematics, so amp and ampf have switched
      ! meanings in this routine
      mp(:,:) = p(:,:)
      mp(1,:) = -p(3,:)
      mp(2,:) = -p(4,:)
      mp(3,:) = -p(5,:)
      mp(5,:) = p(2,:)
      mp(6,:) = p(1,:)

      call spinoru_s(6,mp,za,zb)

      ampsZ = 0._dp 
      ampsZf = 0._dp

      s56 = dotvec(mp(5,:)+mp(6,:),mp(5,:)+mp(6,:))
      s12 = dotvec(mp(1,:)+mp(2,:),mp(1,:)+mp(2,:))

      DZ = s56 - zmass**2 + im*zwidth*zmass
      DGamma = s56

      DZf = s12 - zmass**2 + im*zwidth*zmass
      DGammaf = s12

      ! +++
      do j=1,2
      call zgamma_ampf_rr_vdecay(5,6,1,3,2,j,za,zb,amp_rr_Z)
      do flavor=1,2
      ampsZf(-j,j, flavor, helRight,helRight,helRight,:) = 
     &    amp_rr_Z(:) * r(j) * r(flavor) / DZf
     &   +amp_rr_Z(:) * Q(j) * q(flavor) / DGammaf

      ampsZf(j,-j, flavor, helLeft,helRight,helRight,:) = 
     &   -amp_rr_Z(:) * l(j) * r(flavor) / DZf
     &   -amp_rr_Z(:) * Q(j) * q(flavor) / DGammaf
      enddo; enddo

      ! ---
      do j=1,2
      call zgamma_ampf_rr_vdecay(5,6,1,3,2,j,zb,za,amp_rr_Z)
      do flavor=1,2
      ampsZf(-j,j, flavor, helLeft,helLeft,helLeft,:) = 
     &   -amp_rr_Z(:) * l(j) * l(flavor) / DZf
     &   -amp_rr_Z(:) * Q(j) * q(flavor) / DGammaf

      ampsZf(j,-j, flavor, helRight,helLeft,helLeft,:) = 
     &   +amp_rr_Z(:) * r(j) * l(flavor) / DZf
     &   +amp_rr_Z(:) * Q(j) * q(flavor) / DGammaf
      enddo; enddo

      ! -++
      do j=1,2
      call zgamma_ampf_rr_vdecay(6,5,1,3,2,j,za,zb,amp_rr_Z)
      do flavor=1,2
      ampsZf(-j,j, flavor, helLeft,helRight,helRight,:) = 
     &   -amp_rr_Z(:) * l(j) * r(flavor) / DZf
     &   -amp_rr_Z(:) * Q(j) * q(flavor) / DGammaf

      ampsZf(j,-j, flavor, helRight,helRight,helRight,:) = 
     &   +amp_rr_Z(:) * r(j) * r(flavor) / DZf
     &   +amp_rr_Z(:) * Q(j) * q(flavor) / DGammaf
      enddo; enddo

      ! +--
      do j=1,2
      call zgamma_ampf_rr_vdecay(6,5,1,3,2,j,zb,za,amp_rr_Z)
      do flavor=1,2
      ampsZf(-j,j, flavor, helRight,helLeft,helLeft,:) = 
     &   +amp_rr_Z(:) * r(j) * l(flavor) / DZf
     &   +amp_rr_Z(:) * Q(j) * q(flavor) / DGammaf

      ampsZf(j,-j, flavor, helLeft,helLeft,helLeft,:) = 
     &   -amp_rr_Z(:) * l(j) * l(flavor) / DZf
     &   -amp_rr_Z(:) * Q(j) * q(flavor) / DGammaf
      enddo; enddo

      ! --+
      do j=1,2
      do flavor=1,2
      call zgamma_ampf_rr_vdecay(6,5,2,3,1,j,za,zb,amp_rr_Z)
      ampsZf(-j,j, flavor, helLeft,helLeft,helRight,:) = 
     &   -amp_rr_Z(:) * l(j) * l(flavor) / DZf
     &   -amp_rr_Z(:) * Q(j) * q(flavor) / DGammaf

      ampsZf(j,-j, flavor, helRight,helLeft,helRight,:) = 
     &   +amp_rr_Z(:) * r(j) * l(flavor) / DZf
     &   +amp_rr_Z(:) * Q(j) * q(flavor) / DGammaf
      enddo; enddo

      ! +-+
      do j=1,2
      call zgamma_ampf_rr_vdecay(5,6,2,3,1,j,za,zb,amp_rr_Z)
      do flavor=1,2
      ampsZf(-j,j, flavor, helRight,helLeft,helRight,:) = 
     &    amp_rr_Z(:) * r(j) * l(flavor) / DZf
     &   +amp_rr_Z(:) * Q(j) * q(flavor) / DGammaf

      ampsZf(j,-j, flavor, helLeft,helLeft,helRight,:) = 
     &   -amp_rr_Z(:) * l(j) * l(flavor) / DZf
     &   -amp_rr_Z(:) * Q(j) * q(flavor) / DGammaf
      enddo; enddo

      ! ++-
      do j=1,2
      call zgamma_ampf_rr_vdecay(6,5,2,3,1,j,zb,za,amp_rr_Z)
      do flavor=1,2
      ampsZf(-j,j, flavor, helRight,helRight,helLeft,:) = 
     &    amp_rr_Z(:) * r(j) * r(flavor) / DZf
     &   +amp_rr_Z(:) * Q(j) * q(flavor) / DGammaf

      ampsZf(j,-j, flavor, helLeft,helRight,helLeft,:) = 
     &   -amp_rr_Z(:) * l(j) * r(flavor) / DZf
     &   -amp_rr_Z(:) * Q(j) * q(flavor) / DGammaf
      enddo; enddo


      ! -+-
      do j=1,2
      call zgamma_ampf_rr_vdecay(5,6,2,3,1,j,zb,za,amp_rr_Z)
      do flavor=1,2
      ampsZf(-j,j, flavor, helLeft,helRight,helLeft,:) = 
     &   -amp_rr_Z(:) * l(j) * r(flavor) / DZf
     &   -amp_rr_Z(:) * Q(j) * q(flavor) / DGammaf

      ampsZf(j,-j, flavor, helRight,helRight,helLeft,:) = 
     &   +amp_rr_Z(:) * r(j) * r(flavor) / DZf
     &   +amp_rr_Z(:) * Q(j) * q(flavor) / DGammaf
      enddo; enddo


      call omega_zgamma_vdecay(mp,5,6,1,3,2,maxLoops,c_alphai,c_betai,c_gammai)
      do flavor=1,2

      call zgamma_amp_rr_vdecay(5,6,1,3,2,za,zb,c_alphai,c_betai,c_gammai,
     &                      flavor,helRight,bosZ,amp_rr_Z)
      call zgamma_amp_rr_vdecay(5,6,1,3,2,za,zb,c_alphai,c_betai,c_gammai,
     &                      flavor,helRight,bosGamma,amp_rr_Gamma)

      do j=1,2
      ampsZ(-j,j, flavor, 2,2,2,:) =
     &   +amp_rr_Z(:) * r(j) * r(flavor) / DZ
     &   +amp_rr_Gamma(:) * Q(j) * q(flavor) / DGamma

      ampsZ(j,-j, flavor, 1,2,2,:) =
     &   +amp_rr_Z(:) * l(j) * r(flavor) / DZ
     &   +amp_rr_Gamma(:) * Q(j) * q(flavor) / DGamma
      enddo
      
      call zgamma_amp_rr_vdecay(5,6,1,3,2,zb,za,c_alphai,c_betai,c_gammai,
     &                      flavor,helRight,bosZ,amp_rr_Z)
      call zgamma_amp_rr_vdecay(5,6,1,3,2,zb,za,c_alphai,c_betai,c_gammai,
     &                      flavor,helRight,bosGamma,amp_rr_Gamma)

      do j=1,2
      ampsZ(-j,j, flavor, 1,1,1,:) =
     &   -amp_rr_Z(:) * l(j) * l(flavor) / DZ
     &   -amp_rr_Gamma(:) * Q(j) * q(flavor) / DGamma

      ampsZ(j,-j, flavor, 2,1,1,:) =
     &   -amp_rr_Z(:) * r(j) * l(flavor) / DZ
     &   -amp_rr_Gamma(:) * Q(j) * q(flavor) / DGamma
      enddo
      
      enddo

      call omega_zgamma_vdecay(mp,5,6,2,3,1,maxLoops,c_alphai,c_betai,c_gammai)
      do flavor=1,2

      call zgamma_amp_rr_vdecay(5,6,2,3,1,za,zb,c_alphai,c_betai,c_gammai,
     &                      flavor,helRight,bosZ,amp_rr_Z)
      call zgamma_amp_rr_vdecay(5,6,2,3,1,za,zb,c_alphai,c_betai,c_gammai,
     &                      flavor,helRight,bosGamma,amp_rr_Gamma)

      do j=1,2
      ampsZ(-j,j, flavor, 2,1,2,:) =
     &   -amp_rr_Z(:) * r(j) * l(flavor) / DZ
     &   -amp_rr_Gamma(:) * Q(j) * q(flavor) / DGamma

      ampsZ(j,-j, flavor, 1,1,2,:) =
     &   -amp_rr_Z(:) * l(j) * l(flavor) / DZ
     &   -amp_rr_Gamma(:) * Q(j) * q(flavor) / DGamma
      enddo


      call zgamma_amp_rr_vdecay(5,6,2,3,1,zb,za,c_alphai,c_betai,c_gammai,
     &                      flavor,helRight,bosZ,amp_rr_Z)
      call zgamma_amp_rr_vdecay(5,6,2,3,1,zb,za,c_alphai,c_betai,c_gammai,
     &                      flavor,helRight,bosGamma,amp_rr_Gamma)

      do j=1,2
      ampsZ(-j,j, flavor, 1,2,1,:) =
     &   +amp_rr_Z(:) * l(j) * r(flavor) / DZ
     &   +amp_rr_Gamma(:) * Q(j) * q(flavor) / DGamma

      ampsZ(j,-j, flavor, 2,2,1,:) =
     &   +amp_rr_Z(:) * r(j) * r(flavor) / DZ
     &   +amp_rr_Gamma(:) * Q(j) * q(flavor) / DGamma
      enddo

      enddo


      call omega_zgamma_vdecay(mp,6,5,1,3,2,maxLoops,c_alphai,c_betai,c_gammai)
      do flavor=1,2

      call zgamma_amp_rr_vdecay(6,5,1,3,2,za,zb,c_alphai,c_betai,c_gammai,
     &                      flavor,helRight,bosZ,amp_rr_Z)
      call zgamma_amp_rr_vdecay(6,5,1,3,2,za,zb,c_alphai,c_betai,c_gammai,
     &                      flavor,helRight,bosGamma,amp_rr_Gamma)

      do j=1,2
      ampsZ(-j,j, flavor, 1,2,2,:) =
     &   +amp_rr_Z(:) * l(j) * r(flavor) / DZ
     &   +amp_rr_Gamma(:) * Q(j) * q(flavor) / DGamma

      ampsZ(j,-j, flavor, 2,2,2,:) =
     &   +amp_rr_Z(:) * r(j) * r(flavor) / DZ
     &   +amp_rr_Gamma(:) * Q(j) * q(flavor) / DGamma
      enddo
      
      call zgamma_amp_rr_vdecay(6,5,1,3,2,zb,za,c_alphai,c_betai,c_gammai,
     &                      flavor,helRight,bosZ,amp_rr_Z)
      call zgamma_amp_rr_vdecay(6,5,1,3,2,zb,za,c_alphai,c_betai,c_gammai,
     &                      flavor,helRight,bosGamma,amp_rr_Gamma)

      do j=1,2
      ampsZ(-j,j, flavor, 2,1,1,:) =
     &   -amp_rr_Z(:) * r(j) * l(flavor) / DZ
     &   -amp_rr_Gamma(:) * Q(j) * q(flavor) / DGamma

      ampsZ(j,-j, flavor, 1,1,1,:) =
     &   -amp_rr_Z(:) * l(j) * l(flavor) / DZ
     &   -amp_rr_Gamma(:) * Q(j) * q(flavor) / DGamma
      enddo
      
      enddo


      call omega_zgamma_vdecay(mp,6,5,2,3,1,maxLoops,c_alphai,c_betai,c_gammai)
      do flavor=1,2

      call zgamma_amp_rr_vdecay(6,5,2,3,1,za,zb,c_alphai,c_betai,c_gammai,
     &                      flavor,helRight,bosZ,amp_rr_Z)
      call zgamma_amp_rr_vdecay(6,5,2,3,1,za,zb,c_alphai,c_betai,c_gammai,
     &                      flavor,helRight,bosGamma,amp_rr_Gamma)

      do j=1,2
      ampsZ(-j,j, flavor, 1,1,2,:) =
     &   -amp_rr_Z(:) * l(j) * l(flavor) / DZ
     &   -amp_rr_Gamma(:) * Q(j) * q(flavor) / DGamma

      ampsZ(j,-j, flavor, 2,1,2,:) =
     &   -amp_rr_Z(:) * r(j) * l(flavor) / DZ
     &   -amp_rr_Gamma(:) * Q(j) * q(flavor) / DGamma
      enddo
      
      call zgamma_amp_rr_vdecay(6,5,2,3,1,zb,za,c_alphai,c_betai,c_gammai,
     &                      flavor,helRight,bosZ,amp_rr_Z)
      call zgamma_amp_rr_vdecay(6,5,2,3,1,zb,za,c_alphai,c_betai,c_gammai,
     &                      flavor,helRight,bosGamma,amp_rr_Gamma)

      do j=1,2
      ampsZ(-j,j, flavor, 2,2,1,:) =
     &   +amp_rr_Z(:) * r(j) * r(flavor) / DZ
     &   +amp_rr_Gamma(:) * Q(j) * q(flavor) / DGamma

      ampsZ(j,-j, flavor, 1,2,1,:) =
     &   +amp_rr_Z(:) * l(j) * r(flavor) / DZ
     &   +amp_rr_Gamma(:) * Q(j) * q(flavor) / DGamma
      enddo
      
      enddo

      if (zerowidth) ampsZ = 0._dp

      ampsZ = ampsZ + ampsZf

      end subroutine

      ! this is copied from zgamma_ampf_RR with added qf argument
      ! to be used for crossed amplitude
      subroutine zgamma_ampf_RR_vdecay(i5,i6,i1,i3,i2,qf,za,zb,amps)
          use VVconfig_m
        implicit none
        include 'types.f'
        include 'mxpart.f'
        include 'nf.f'
        include 'zcouple.f'
        include 'ewcharge.f'

        integer, intent(in) :: i5,i6,i1,i3,i2
        integer, intent(in) :: qf
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(out) :: amps(0:2)
        complex(dp) :: tree

        tree = -Q(qf)*za(i2,i5)**2*zb(i2,i1)/(za(i5,i3)*za(i6,i3))

        call zgamma_amp_qff(i5,i6,i1,i3,i2,za,zb,tree,amps)

      end subroutine

      ! this is copied from zgamma_amp_RR with fixed eq = q1
      ! to be used for crossed amplitude
      subroutine zgamma_amp_RR_vdecay(i5,i6,i1,i3,i2,za,zb,
     &   c_alphai, c_betai, c_gammai, flavor, helq, bos,amps)
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
      integer, intent(in) :: flavor
      integer (kind(helLeft)), intent(in) :: helq
      integer (kind(bosZ)), intent(in) :: bos
      complex(dp), intent(out) :: amps(0:2)

      complex(dp) :: c_alpha(0:2)
      complex(dp) :: c_beta(0:2)
      complex(dp) :: c_gamma(0:2)
      real(dp) :: s12, s123
      complex(dp) :: wzgamma_amp_RR, omega_I1
      real(dp) :: vl,vr
      integer :: j
      real(dp) :: NFV
      real(dp) :: eq

      s12 = real(za(i1,i2)*zb(i2,i1))
      s123 = real(za(i1,i2)*zb(i2,i1) + za(i1,i3)*zb(i3,i1) +
     &                             za(i2,i3)*zb(i3,i2))

      eq = Q(flavor)

      NFV = 0
      if (bos == bosZ) then
        do j=1,nf
            NFV = NFV - Q(j)*(l(j) + r(j))
        enddo

        if (helQ == helLeft ) then
          ! I just hate signs, are these really OK?
          ! How can we check this?
          vl = l(flavor)
          NFV = NFV/2/vl
        elseif (helQ == helRight) then
          vr = r(flavor)
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

