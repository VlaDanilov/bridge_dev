      module hjetmass_highpt
          use debugtools_m
          implicit none

          public hjetmass_highpt_lo, hjetmass_highpt_amps, I1g, I1q

          private

      contains

      subroutine hjetmass_highpt_amps(p,ampsqq,ampsgq,ampsqg,ampsgg)
          use hjetmass_hel
          implicit none
          include 'types.f'
          include 'mxpart.f'
          include 'nf.f'
          include 'constants.f'
          include 'masses.f'
          include 'hdecaymode.f'
          include 'qcdcouple.f'
          include 'ewcouple.f'
          include 'scale.f'
          include 'zprods_decl.f'

          real(dp), intent(in) :: p(mxpart,4)
          complex(dp), intent(out) :: ampsqq(2, 2,2), ampsgq(2, 2,2)
          complex(dp), intent(out) :: ampsqg(2, 2,2), ampsgg(2, 2,2,2)

          integer, parameter :: iglue = 5

          complex(dp) :: omega2l(2,3)

          call spinoru_dp(5,p,za,zb)

          ampsqq = 0._dp
          ampsgq = 0._dp
          ampsqg = 0._dp
          ampsgg = 0._dp

          ! 125
          call hhpt_omega_2l(p,1,2,iglue,omega2l)

          ampsgg(:, 1,1,1) = hhpt_ggg_ppp(za,zb, 1,2,iglue)*omega2l(:,1)
          ampsgg(:, 2,2,2) = hhpt_ggg_ppp(zb,za, 1,2,iglue)*omega2l(:,1)

          ampsgg(:, 1,2,1) = hhpt_ggg_pmp(za,zb, 1,2,iglue)*omega2l(:,2)
          ampsgg(:, 2,1,2) = hhpt_ggg_pmp(zb,za, 1,2,iglue)*omega2l(:,2)

          ! -++
          ampsqq(:, 1,1) = hhpt_qqg_mpp(za,zb, 1,2,iglue)*omega2l(:,3)
          ! +s--
          ampsqq(:, 2,2) = hhpt_qqg_mpp(zb,za, 1,2,iglue)*omega2l(:,3)

          ! 215
          call hhpt_omega_2l(p,2,1,iglue,omega2l)

          ampsgg(:, 1,2,2) = hhpt_ggg_pmp(zb,za, 2,1,iglue)*omega2l(:,2)
          ampsgg(:, 2,1,1) = hhpt_ggg_pmp(za,zb, 2,1,iglue)*omega2l(:,2)

          ! +-+
          ampsqq(:, 2,1) = hhpt_qqg_mpp(za,zb, 2,1,iglue)*omega2l(:,3)
          ! -+-
          ampsqq(:, 1,2) = hhpt_qqg_mpp(zb,za, 2,1,iglue)*omega2l(:,3)

          ! 152
          call hhpt_omega_2l(p,1,iglue,2,omega2l)

          ampsgg(:, 1,1,2) = hhpt_ggg_pmp(za,zb, 1,iglue,2)*omega2l(:,2)
          ampsgg(:, 2,2,1) = hhpt_ggg_pmp(zb,za, 1,iglue,2)*omega2l(:,2)

          ! -++
          ampsqg(:, 1,1) = hhpt_qqg_mpp(za,zb, 1,iglue,2)*omega2l(:,3)
          ! +--
          ampsqg(:, 2,2) = hhpt_qqg_mpp(zb,za, 1,iglue,2)*omega2l(:,3)

          ! 512
          call hhpt_omega_2l(p,iglue,1,2,omega2l)
          ! +-+
          ampsqg(:, 2,1) = hhpt_qqg_mpp(za,zb, iglue,1,2)*omega2l(:,3)
          ! -+-
          ampsqg(:, 1,2) = hhpt_qqg_mpp(zb,za, iglue,1,2)*omega2l(:,3)


          ! 521
          call hhpt_omega_2l(p,iglue,2,1,omega2l)
          ! -++
          ampsgq(:, 1,1) = hhpt_qqg_mpp(za,zb, iglue,2,1)*omega2l(:,3)
          ! +--
          ampsgq(:, 2,2) = hhpt_qqg_mpp(zb,za, iglue,2,1)*omega2l(:,3)

          ! 251
          call hhpt_omega_2l(p,2,iglue,1,omega2l)
          ! +-+
          ampsgq(:, 2,1) = hhpt_qqg_mpp(za,zb, 2,iglue,1)*omega2l(:,3)
          ! -+-
          ampsgq(:, 1,2) = hhpt_qqg_mpp(zb,za, 2,iglue,1)*omega2l(:,3)

      end subroutine
      
      subroutine hjetmass_highpt_lo(p,qq,gq,qg,gg)
          use hjetmass_hel
          implicit none
          include 'types.f'
          include 'mxpart.f'
          include 'nf.f'
          include 'constants.f'
          include 'masses.f'
          include 'hdecaymode.f'
          include 'qcdcouple.f'
          include 'ewcouple.f'
          include 'scale.f'
          include 'zprods_decl.f'
          include 'asymptotic.f'

          real(dp), intent(in) :: p(mxpart,4)
          real(dp), intent(out) :: qq,gq,qg,gg

          integer, parameter :: iglue = 5

          complex(dp) :: omega(2,3)
          complex(dp) :: ampqqg(2, 2,2), ampggg(2, 2,2,2)
          !complex(dp) :: ampex(2,2,2), ampqqg_ex(2,2)

          call spinoru_dp(5,p,za,zb)

          ampqqg = 0._dp
          ampggg = 0._dp

          ! 125
          call hhpt_omega(p,1,2,iglue,omega)

          ampggg(:,1,1,1) = hhpt_ggg_ppp(za,zb, 1,2,iglue)*omega(:,1)
          ampggg(:,2,2,2) = hhpt_ggg_ppp(zb,za, 1,2,iglue)*omega(:,1)

          !ampex(1,1,1) = hjetmass_ggg_ppp(za,zb,1,2,iglue)
          !ampex(2,2,2) = hjetmass_ggg_ppp(zb,za,1,2,iglue)

          ampggg(:,1,2,1) = hhpt_ggg_pmp(za,zb, 1,2,iglue)*omega(:,2)
          ampggg(:,2,1,2) = hhpt_ggg_pmp(zb,za, 1,2,iglue)*omega(:,2)

          !ampex(1,2,1) = hjetmass_ggg_pmp(za,zb,1,2,iglue)
          !ampex(2,1,2) = hjetmass_ggg_pmp(zb,za,1,2,iglue)

          ! -++
          ampqqg(:,1,1) = hhpt_qqg_mpp(za,zb, 1,2,iglue)*omega(:,3)
          ! +--
          ampqqg(:,2,2) = hhpt_qqg_mpp(zb,za, 1,2,iglue)*omega(:,3)

          ! 215
          call hhpt_omega(p,2,1,iglue,omega)

          ampggg(:,1,2,2) = hhpt_ggg_pmp(zb,za, 2,1,iglue)*omega(:,2)
          ampggg(:,2,1,1) = hhpt_ggg_pmp(za,zb, 2,1,iglue)*omega(:,2)

          !ampex(1,2,2) = hjetmass_ggg_pmp(zb,za,2,1,iglue)
          !ampex(2,1,1) = hjetmass_ggg_pmp(za,zb,2,1,iglue)

          ! +-+
          ampqqg(:,2,1) = hhpt_qqg_mpp(za,zb, 2,1,iglue)*omega(:,3)
          ! -+-
          ampqqg(:,1,2) = hhpt_qqg_mpp(zb,za, 2,1,iglue)*omega(:,3)

          !ampqqg_ex(2,2) = hjetmass_qqg_mpm(za,zb,1,2,iglue)
          !ampqqg_ex(1,1) = hjetmass_qqg_mpm(zb,za,1,2,iglue)
          !ampqqg_ex(2,1) = hjetmass_qqg_mpp(za,zb,1,2,iglue)
          !ampqqg_ex(1,2) = hjetmass_qqg_mpp(zb,za,1,2,iglue)

          if (mtex >= 0) then
            qq = sum(abs(ampqqg(1,:,:))**2)
          endif

          if (mtex >= 2) then
            qq = qq + sum(ampqqg(1,:,:)*conjg(ampqqg(2,:,:)))*mt**2
            qq = qq + sum(ampqqg(2,:,:)*conjg(ampqqg(1,:,:)))*mt**2
          endif

          if (mtex >= 4) then
            qq = qq + sum(ampqqg(2,:,:)*conjg(ampqqg(2,:,:)))*mt**4
          endif

          ampqqg = 0._dp

          ! 152
          call hhpt_omega(p,1,iglue,2,omega)

          ampggg(:,1,1,2) = hhpt_ggg_pmp(za,zb, 1,iglue,2)*omega(:,2)
          ampggg(:,2,2,1) = hhpt_ggg_pmp(zb,za, 1,iglue,2)*omega(:,2)

          !ampex(1,1,2) = hjetmass_ggg_pmp(za,zb,1,iglue,2)
          !ampex(2,2,1) = hjetmass_ggg_pmp(zb,za,1,iglue,2)

          ! -++
          ampqqg(:,1,1) = hhpt_qqg_mpp(za,zb, 1,iglue,2)*omega(:,3)
          ! +--
          ampqqg(:,2,2) = hhpt_qqg_mpp(zb,za, 1,iglue,2)*omega(:,3)

          ! 512
          call hhpt_omega(p,iglue,1,2,omega)
          ! +-+
          ampqqg(:,2,1) = hhpt_qqg_mpp(za,zb, iglue,1,2)*omega(:,3)
          ! -+-
          ampqqg(:,1,2) = hhpt_qqg_mpp(zb,za, iglue,1,2)*omega(:,3)

          if (mtex >= 0) then
            qg = sum(abs(ampqqg(1,:,:))**2)
          endif

          if (mtex >= 2) then
            qg = qg + sum(ampqqg(1,:,:)*conjg(ampqqg(2,:,:)))*mt**2
            qg = qg + sum(ampqqg(2,:,:)*conjg(ampqqg(1,:,:)))*mt**2
          endif

          if (mtex >= 4) then
            qg = qg + sum(ampqqg(2,:,:)*conjg(ampqqg(2,:,:)))*mt**4
          endif

          ampqqg = 0._dp


          ! 521
          call hhpt_omega(p,iglue,2,1,omega)
          ! -++
          ampqqg(:,1,1) = hhpt_qqg_mpp(za,zb, iglue,2,1)*omega(:,3)
          ! +--
          ampqqg(:,2,2) = hhpt_qqg_mpp(zb,za, iglue,2,1)*omega(:,3)

          ! 251
          call hhpt_omega(p,2,iglue,1,omega)
          ! +-+
          ampqqg(:,2,1) = hhpt_qqg_mpp(za,zb, 2,iglue,1)*omega(:,3)
          ! -+-
          ampqqg(:,1,2) = hhpt_qqg_mpp(zb,za, 2,iglue,1)*omega(:,3)

          if (mtex >= 0) then
            gq = sum(abs(ampqqg(1,:,:))**2)
          endif

          if (mtex >= 2) then
            gq = gq + sum(ampqqg(1,:,:)*conjg(ampqqg(2,:,:)))*mt**2
            gq = gq + sum(ampqqg(2,:,:)*conjg(ampqqg(1,:,:)))*mt**2
          endif

          if (mtex >= 4) then
            gq = gq + sum(ampqqg(2,:,:)*conjg(ampqqg(2,:,:)))*mt**4
          endif

          if (mtex >= 0) then
            gg = sum(abs(ampggg(1,:,:,:))**2)
          endif

          if (mtex >= 2) then
            gg = gg + sum(ampggg(1,:,:,:)*conjg(ampggg(2,:,:,:)))*mt**2
            gg = gg + sum(ampggg(2,:,:,:)*conjg(ampggg(1,:,:,:)))*mt**2
          endif

          if (mtex >= 4) then
            gg = gg + sum(ampggg(2,:,:,:)*conjg(ampggg(2,:,:,:)))*mt**4
          endif

      end subroutine

      pure function hhpt_qqg_mpp(za,zb,i1,i2,i3)
          implicit none
          include 'types.f'
          include 'mxpart.f'

          complex(dp) :: hhpt_qqg_mpp

          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          integer, intent(in) :: i1,i2,i3

          real(dp) :: s

          s = real(za(i1,i2)*zb(i2,i1))

          hhpt_qqg_mpp = 1._dp/sqrt(2._dp) * zb(i2,i3)**2/zb(i1,i2)/s

      end function

      pure function hhpt_ggg_ppp(za,zb,i1,i2,i3)
          implicit none
          include 'types.f'
          include 'mxpart.f'

          complex(dp) :: hhpt_ggg_ppp

          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          integer, intent(in) :: i1,i2,i3

          real(dp) :: s
          s = real(za(i1,i2)*zb(i2,i1))

          hhpt_ggg_ppp = s/sqrt(2._dp)/za(i1,i2)/za(i2,i3)/za(i3,i1)
      end function

      pure function hhpt_ggg_pmp(za,zb,i1,i2,i3)
          implicit none
          include 'types.f'
          include 'mxpart.f'

          complex(dp) :: hhpt_ggg_pmp
          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          integer, intent(in) :: i1,i2,i3

          real(dp) :: s
          s = real(za(i1,i2)*zb(i2,i1))

          hhpt_ggg_pmp = zb(i1,i3)**3/sqrt(2._dp)/zb(i1,i2)/zb(i3,i2)/s

      end function

      subroutine hhpt_omega(p,i1,i2,i3,omega)
      !! this subroutine is for the one loop amps only
      !! for different expansion orders
          implicit none
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'masses.f'
          include 'qcdcouple.f'
          include 'ewcouple.f'

          real(dp), intent(in) :: p(mxpart,4)
          integer, intent(in) :: i1,i2,i3
          complex(dp), intent(out) :: omega(2,3)

          complex(dp) :: omega2l(2,3)

          complex(dp) :: w(40)
          complex(dp) :: omega_ggg_ppp_1l_2a_mt0
          complex(dp) :: omega_ggg_pmp_1l_2a_mt0
          complex(dp) :: omega_qag_mpp_1l_2a_mt0
          complex(dp) :: omega_ggg_ppp_1l_3a_mt0
          complex(dp) :: omega_ggg_pmp_1l_3a_mt0
          complex(dp) :: omega_qag_mpp_1l_3a_mt0
          complex(dp) :: omega_ggg_ppp_1l_4a_mt0
          complex(dp) :: omega_ggg_pmp_1l_4a_mt0
          complex(dp) :: omega_qag_mpp_1l_4a_mt0

          complex(dp) :: omega_ggg_ppp_1l_2a_mt2
          complex(dp) :: omega_ggg_pmp_1l_2a_mt2
          complex(dp) :: omega_qag_mpp_1l_2a_mt2
          complex(dp) :: omega_ggg_ppp_1l_3a_mt2
          complex(dp) :: omega_ggg_pmp_1l_3a_mt2
          complex(dp) :: omega_qag_mpp_1l_3a_mt2
          complex(dp) :: omega_ggg_ppp_1l_4a_mt2
          complex(dp) :: omega_ggg_pmp_1l_4a_mt2
          complex(dp) :: omega_qag_mpp_1l_4a_mt2

          integer, parameter :: n1=0, n2=1, nw = 2

          complex(dp):: Hc1(n1:n2),Hc2(n1:n2,n1:n2),Hc3(n1:n2,n1:n2,n1:n2), 
     &              Hc4(n1:n2,n1:n2,n1:n2,n1:n2) 
          real(dp):: Hr1(n1:n2),Hr2(n1:n2,n1:n2),Hr3(n1:n2,n1:n2,n1:n2), 
     &              Hr4(n1:n2,n1:n2,n1:n2,n1:n2) 
          real(dp):: Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     &              Hi4(n1:n2,n1:n2,n1:n2,n1:n2)

          real(dp) :: s,t,u,u2,u3,u4
          real(dp) :: mt2
          real(dp) :: dotvec

          mt2 = mt**2

          w = 0._dp
          omega = 0._dp

          s = dotvec(p(i1,:)+p(i2,:), p(i1,:)+p(i2,:))
          t = dotvec(p(i1,:)+p(i3,:), p(i1,:)+p(i3,:))
          u = dotvec(p(i2,:)+p(i3,:), p(i2,:)+p(i3,:))

          u2 = -u/s
          u3 = -s/t
          u4 = -s/u

          if (s > 0._dp) then
            call hplog(u2,nw,Hc1,Hc2,Hc3,Hc4,Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)

            Hr1(1) = -Hr1(1)
            Hr2(0,1) = -Hr2(0,1)
            Hr2(1,0) = -Hr2(1,0)

            include 'hhpt_omega_1l_2a.f'
            omega(1,1) = omega_ggg_ppp_1l_2a_mt0
            omega(1,2) = omega_ggg_pmp_1l_2a_mt0
            omega(1,3) = omega_qag_mpp_1l_2a_mt0

            call hhpt_omega_2l(p,i1,i2,i3,omega2l)

            omega(2,1) = omega_ggg_ppp_1l_2a_mt2
            omega(2,2) = omega_ggg_pmp_1l_2a_mt2
            omega(2,3) = omega_qag_mpp_1l_2a_mt2
          else if (t > 0._dp) then
            call hplog(u3,nw,Hc1,Hc2,Hc3,Hc4,Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)

            Hr1(1) = -Hr1(1)
            Hr2(0,1) = -Hr2(0,1)
            Hr2(1,0) = -Hr2(1,0)

            include 'hhpt_omega_1l_3a.f'
            omega(1,1) = omega_ggg_ppp_1l_3a_mt0
            omega(1,2) = omega_ggg_pmp_1l_3a_mt0
            omega(1,3) = omega_qag_mpp_1l_3a_mt0

            omega(2,1) = omega_ggg_ppp_1l_3a_mt2
            omega(2,2) = omega_ggg_pmp_1l_3a_mt2
            omega(2,3) = omega_qag_mpp_1l_3a_mt2
          else if (u > 0._dp) then
            call hplog(u4,nw,Hc1,Hc2,Hc3,Hc4,Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)

            Hr1(1) = -Hr1(1)
            Hr2(0,1) = -Hr2(0,1)
            Hr2(1,0) = -Hr2(1,0)

            include 'hhpt_omega_1l_4a.f'
            omega(1,1) = omega_ggg_ppp_1l_4a_mt0
            omega(1,2) = omega_ggg_pmp_1l_4a_mt0
            omega(1,3) = omega_qag_mpp_1l_4a_mt0

            omega(2,1) = omega_ggg_ppp_1l_4a_mt2
            omega(2,2) = omega_ggg_pmp_1l_4a_mt2
            omega(2,3) = omega_qag_mpp_1l_4a_mt2
          else
              call abort
          endif

          omega = omega * mt**2 / sqrt(vevsq) * sqrt(as**3/pi)

      end subroutine

      subroutine hhpt_omega_2l(p,i1,i2,i3,omega)
      !! this subroutine is for the one and two loop results combined
      !! using HPLs
          implicit none
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'masses.f'
          include 'qcdcouple.f'
          include 'ewcouple.f'
          include 'zeta.f'
          include 'nf.f' ! for Nf in amplitude
          include 'scale.f'

          real(dp), intent(in) :: p(mxpart,4)
          integer, intent(in) :: i1,i2,i3
          complex(dp), intent(out) :: omega(2,3)

          complex(dp) :: w(382)
          complex(dp) :: omega_ggg_ppp_1l_2a
          complex(dp) :: omega_ggg_pmp_1l_2a
          complex(dp) :: omega_qag_mpp_1l_2a
          complex(dp) :: omega_ggg_ppp_1l_3a
          complex(dp) :: omega_ggg_pmp_1l_3a
          complex(dp) :: omega_qag_mpp_1l_3a
          complex(dp) :: omega_ggg_ppp_1l_4a
          complex(dp) :: omega_ggg_pmp_1l_4a
          complex(dp) :: omega_qag_mpp_1l_4a

          complex(dp) :: omega_ggg_ppp_2l_2a
          complex(dp) :: omega_ggg_pmp_2l_2a
          complex(dp) :: omega_qag_mpp_2l_2a
          complex(dp) :: omega_ggg_ppp_2l_3a
          complex(dp) :: omega_ggg_pmp_2l_3a
          complex(dp) :: omega_qag_mpp_2l_3a
          complex(dp) :: omega_ggg_ppp_2l_4a
          complex(dp) :: omega_ggg_pmp_2l_4a
          complex(dp) :: omega_qag_mpp_2l_4a

          real(dp) :: s,t,u,u2,u3,u4
          real(dp) :: mt2
          real(dp) :: dotvec
          real(dp) :: beta0

          integer, parameter :: n1=0, n2=1, nw = 4
          
          complex(dp) :: li2

          complex(dp):: Hc1(n1:n2),Hc2(n1:n2,n1:n2),Hc3(n1:n2,n1:n2,n1:n2), 
     &              Hc4(n1:n2,n1:n2,n1:n2,n1:n2) 
          real(dp):: Hr1(n1:n2),Hr2(n1:n2,n1:n2),Hr3(n1:n2,n1:n2,n1:n2), 
     &              Hr4(n1:n2,n1:n2,n1:n2,n1:n2) 
          real(dp):: Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     &              Hi4(n1:n2,n1:n2,n1:n2,n1:n2)

          mt2 = mt**2

          w = 0._dp
          omega = 0._dp

          s = dotvec(p(i1,:)+p(i2,:), p(i1,:)+p(i2,:))
          t = dotvec(p(i1,:)+p(i3,:), p(i1,:)+p(i3,:))
          u = dotvec(p(i2,:)+p(i3,:), p(i2,:)+p(i3,:))

          u2 = -u/s
          u3 = -s/t
          u4 = -s/u

          beta0 = 11d0/6d0*CA - 2d0/3d0*TR*Nf

          if (s > 0._dp) then
            call hplog(u2,nw,Hc1,Hc2,Hc3,Hc4,Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)

            ! H to G conversion
            Hr1(1) = -Hr1(1)
            Hr2(0,1) = -Hr2(0,1)
            Hr2(1,0) = -Hr2(1,0)

            Hr3(0,0,1) = -Hr3(0,0,1)
            Hr3(0,1,0) = -Hr3(0,1,0)
            Hr3(1,0,0) = -Hr3(1,0,0)
            Hr3(1,1,1) = -Hr3(1,1,1)

            Hr4(0,0,0,1) = -Hr4(0,0,0,1)
            Hr4(0,0,1,0) = -Hr4(0,0,1,0)
            Hr4(0,1,0,0) = -Hr4(0,1,0,0)
            Hr4(1,0,0,0) = -Hr4(1,0,0,0)

            Hr4(0,1,1,1) = -Hr4(0,1,1,1)
            Hr4(1,0,1,1) = -Hr4(1,0,1,1)
            Hr4(1,1,0,1) = -Hr4(1,1,0,1)
            Hr4(1,1,1,0) = -Hr4(1,1,1,0)

            include 'hhpt_omega_2l_2a.f'

            omega(1, 1) = omega_ggg_ppp_1l_2a
            omega(1, 2) = omega_ggg_pmp_1l_2a
            omega(1, 3) = omega_qag_mpp_1l_2a

            omega(2, 1) = omega_ggg_ppp_2l_2a
            omega(2, 2) = omega_ggg_pmp_2l_2a
            omega(2, 3) = omega_qag_mpp_2l_2a

            omega(2,:) = omega(2,:) + omega(1,:) * 3d0/2d0 * beta0 * log(musq/s)

          else if (t > 0._dp) then
            call hplog(u3,nw,Hc1,Hc2,Hc3,Hc4,Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)
            Hr1(1) = -Hr1(1)
            Hr2(0,1) = -Hr2(0,1)
            Hr2(1,0) = -Hr2(1,0)

            Hr3(0,0,1) = -Hr3(0,0,1)
            Hr3(0,1,0) = -Hr3(0,1,0)
            Hr3(1,0,0) = -Hr3(1,0,0)
            Hr3(1,1,1) = -Hr3(1,1,1)

            Hr4(0,0,0,1) = -Hr4(0,0,0,1)
            Hr4(0,0,1,0) = -Hr4(0,0,1,0)
            Hr4(0,1,0,0) = -Hr4(0,1,0,0)
            Hr4(1,0,0,0) = -Hr4(1,0,0,0)

            Hr4(0,1,1,1) = -Hr4(0,1,1,1)
            Hr4(1,0,1,1) = -Hr4(1,0,1,1)
            Hr4(1,1,0,1) = -Hr4(1,1,0,1)
            Hr4(1,1,1,0) = -Hr4(1,1,1,0)

            include 'hhpt_omega_2l_3a.f'

            omega(1, 1) = omega_ggg_ppp_1l_3a
            omega(1, 2) = omega_ggg_pmp_1l_3a
            omega(1, 3) = omega_qag_mpp_1l_3a

            omega(2, 1) = omega_ggg_ppp_2l_3a
            omega(2, 2) = omega_ggg_pmp_2l_3a
            omega(2, 3) = omega_qag_mpp_2l_3a

            omega(2,:) = omega(2,:) + omega(1,:) * 3d0/2d0 * beta0 * log(musq/t)

          else if (u > 0._dp) then
            call hplog(u4,nw,Hc1,Hc2,Hc3,Hc4,Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)
            Hr1(1) = -Hr1(1)
            Hr2(0,1) = -Hr2(0,1)
            Hr2(1,0) = -Hr2(1,0)

            Hr3(0,0,1) = -Hr3(0,0,1)
            Hr3(0,1,0) = -Hr3(0,1,0)
            Hr3(1,0,0) = -Hr3(1,0,0)
            Hr3(1,1,1) = -Hr3(1,1,1)

            Hr4(0,0,0,1) = -Hr4(0,0,0,1)
            Hr4(0,0,1,0) = -Hr4(0,0,1,0)
            Hr4(0,1,0,0) = -Hr4(0,1,0,0)
            Hr4(1,0,0,0) = -Hr4(1,0,0,0)

            Hr4(0,1,1,1) = -Hr4(0,1,1,1)
            Hr4(1,0,1,1) = -Hr4(1,0,1,1)
            Hr4(1,1,0,1) = -Hr4(1,1,0,1)
            Hr4(1,1,1,0) = -Hr4(1,1,1,0)
            include 'hhpt_omega_2l_4a.f'

            omega(1, 1) = omega_ggg_ppp_1l_4a
            omega(1, 2) = omega_ggg_pmp_1l_4a
            omega(1, 3) = omega_qag_mpp_1l_4a

            omega(2, 1) = omega_ggg_ppp_2l_4a
            omega(2, 2) = omega_ggg_pmp_2l_4a
            omega(2, 3) = omega_qag_mpp_2l_4a

            omega(2,:) = omega(2,:) + omega(1,:) * 3d0/2d0 * beta0 * log(musq/u)
          else
              call abort
          endif

          ! Amps defined in terms of strong coupling that runs with Nf+1
          ! active flavors. Convert back to Nf here.
          omega(2,:) = omega(2,:) + 1/2._dp * log(musq/mt2) * omega(1,:)

          ! dress with couplings
          omega = omega * mt**2 / sqrt(vevsq) * sqrt(as**3/pi)
          omega(2,:) = omega(2,:) * as/2._dp/pi


      end subroutine

      function I1g(s,t,u,musq)
          implicit none
          include 'types.f'
          include 'nf.f'
          include 'constants.f'
          include 'epinv.f'
          include 'epinv2.f'

          complex(dp) :: I1g
          real(dp), intent(in) :: s,t,u, musq

          real(dp) :: beta0

          beta0 = 11d0/6d0*CA - 2d0/3d0*TR*Nf

          if (s>0) then
      I1g =
     &  - 3.D0/2.D0*epinv*epinv2*CA - 3.D0/2.D0*beta0*epinv - 1.D0/2.D0*Pi*
     & im*epinv*CA - 1.D0/2.D0*Pi*im*beta0 + 1.D0/8.D0*Pi**2*CA - 1.D0/
     & 4.D0*Pi**2*im**2*CA + 1.D0/2.D0*Log( - u*musq**(-1))*epinv*CA + 
     & 1.D0/2.D0*Log( - u*musq**(-1))*beta0 - 1.D0/4.D0*Log( - u*
     & musq**(-1))**2*CA + 1.D0/2.D0*Log( - t*musq**(-1))*epinv*CA + 1.D
     & 0/2.D0*Log( - t*musq**(-1))*beta0 - 1.D0/4.D0*Log( - t*
     & musq**(-1))**2*CA + 1.D0/2.D0*Log(s*musq**(-1))*epinv*CA + 1.D0/
     & 2.D0*Log(s*musq**(-1))*beta0 + 1.D0/2.D0*Log(s*musq**(-1))*Pi*im
     & *CA - 1.D0/4.D0*Log(s*musq**(-1))**2*CA
           else if (t>0) then
      I1g =
     &  - 3.D0/2.D0*epinv*epinv2*CA - 3.D0/2.D0*beta0*epinv - 1.D0/2.D0*Pi*
     & im*epinv*CA - 1.D0/2.D0*Pi*im*beta0 + 1.D0/8.D0*Pi**2*CA - 1.D0/
     & 4.D0*Pi**2*im**2*CA + 1.D0/2.D0*Log( - u*musq**(-1))*epinv*CA + 
     & 1.D0/2.D0*Log( - u*musq**(-1))*beta0 - 1.D0/4.D0*Log( - u*
     & musq**(-1))**2*CA + 1.D0/2.D0*Log(t*musq**(-1))*epinv*CA + 1.D0/
     & 2.D0*Log(t*musq**(-1))*beta0 + 1.D0/2.D0*Log(t*musq**(-1))*Pi*im
     & *CA - 1.D0/4.D0*Log(t*musq**(-1))**2*CA + 1.D0/2.D0*Log( - s*
     & musq**(-1))*epinv*CA + 1.D0/2.D0*Log( - s*musq**(-1))*beta0 - 1.D
     & 0/4.D0*Log( - s*musq**(-1))**2*CA
           else if (u>0) then
      I1g =
     &  - 3.D0/2.D0*epinv*epinv2*CA - 3.D0/2.D0*beta0*epinv - 1.D0/2.D0*Pi*
     & im*epinv*CA - 1.D0/2.D0*Pi*im*beta0 + 1.D0/8.D0*Pi**2*CA - 1.D0/
     & 4.D0*Pi**2*im**2*CA + 1.D0/2.D0*Log(u*musq**(-1))*epinv*CA + 1.D0
     & /2.D0*Log(u*musq**(-1))*beta0 + 1.D0/2.D0*Log(u*musq**(-1))*Pi*
     & im*CA - 1.D0/4.D0*Log(u*musq**(-1))**2*CA + 1.D0/2.D0*Log( - t*
     & musq**(-1))*epinv*CA + 1.D0/2.D0*Log( - t*musq**(-1))*beta0 - 1.D
     & 0/4.D0*Log( - t*musq**(-1))**2*CA + 1.D0/2.D0*Log( - s*
     & musq**(-1))*epinv*CA + 1.D0/2.D0*Log( - s*musq**(-1))*beta0 - 1.D
     & 0/4.D0*Log( - s*musq**(-1))**2*CA
           else
               call abort
           endif


      end function

      function I1q(s,t,u,musq)
          implicit none
          include 'types.f'
          include 'nf.f'
          include 'constants.f'
          include 'epinv.f'
          include 'epinv2.f'

          complex(dp) :: I1q
          real(dp), intent(in) :: s,t,u, musq

          real(dp) :: beta0

          beta0 = 11d0/6d0*CA - 2d0/3d0*TR*Nf

          if (s>0) then
      I1q =
     & 3.D0/4.D0*epinv*CA**(-1) - 3.D0/4.D0*epinv*CA + 1.D0/2.D0*
     & epinv*epinv2*CA**(-1) - epinv*epinv2*CA - 1.D0/2.D0*beta0*epinv + 3.D0/4.
     & D0*Pi*im*CA**(-1) + 1.D0/2.D0*Pi*im*epinv*CA**(-1) - 1.D0/24.D0*
     & Pi**2*CA**(-1) + 1.D0/12.D0*Pi**2*CA + 1.D0/4.D0*Pi**2*im**2*
     & CA**(-1) + 3.D0/8.D0*Log( - u*musq**(-1))*CA + 1.D0/2.D0*Log( - 
     & u*musq**(-1))*epinv*CA + 1.D0/4.D0*Log( - u*musq**(-1))*beta0 - 
     & 1.D0/4.D0*Log( - u*musq**(-1))**2*CA + 3.D0/8.D0*Log( - t*
     & musq**(-1))*CA + 1.D0/2.D0*Log( - t*musq**(-1))*epinv*CA + 1.D0/
     & 4.D0*Log( - t*musq**(-1))*beta0 - 1.D0/4.D0*Log( - t*musq**(-1))
     & **2*CA - 3.D0/4.D0*Log(s*musq**(-1))*CA**(-1) - 1.D0/2.D0*Log(s*
     & musq**(-1))*epinv*CA**(-1) - 1.D0/2.D0*Log(s*musq**(-1))*Pi*im*
     & CA**(-1) + 1.D0/4.D0*Log(s*musq**(-1))**2*CA**(-1)
          else if (t>0) then
      I1q =
     & 3.D0/4.D0*epinv*CA**(-1) - 3.D0/4.D0*epinv*CA + 1.D0/2.D0*
     & epinv*epinv2*CA**(-1) - epinv*epinv2*CA - 1.D0/2.D0*beta0*epinv - 3.D0/8.
     & D0*Pi*im*CA - 1.D0/2.D0*Pi*im*epinv*CA - 1.D0/4.D0*Pi*im*beta0
     &  - 1.D0/24.D0*Pi**2*CA**(-1) + 1.D0/12.D0*Pi**2*CA - 1.D0/4.D0*
     & Pi**2*im**2*CA + 3.D0/8.D0*Log( - u*musq**(-1))*CA + 1.D0/2.D0*
     & Log( - u*musq**(-1))*epinv*CA + 1.D0/4.D0*Log( - u*musq**(-1))*
     & beta0 - 1.D0/4.D0*Log( - u*musq**(-1))**2*CA + 3.D0/8.D0*Log(t*
     & musq**(-1))*CA + 1.D0/2.D0*Log(t*musq**(-1))*epinv*CA + 1.D0/4.D0
     & *Log(t*musq**(-1))*beta0 + 1.D0/2.D0*Log(t*musq**(-1))*Pi*im*CA
     &  - 1.D0/4.D0*Log(t*musq**(-1))**2*CA - 3.D0/4.D0*Log( - s*
     & musq**(-1))*CA**(-1) - 1.D0/2.D0*Log( - s*musq**(-1))*epinv*
     & CA**(-1) + 1.D0/4.D0*Log( - s*musq**(-1))**2*CA**(-1)
          else if (u>0) then
      I1q =
     & 3.D0/4.D0*epinv*CA**(-1) - 3.D0/4.D0*epinv*CA + 1.D0/2.D0*
     & epinv*epinv2*CA**(-1) - epinv*epinv2*CA - 1.D0/2.D0*beta0*epinv - 3.D0/8.
     & D0*Pi*im*CA - 1.D0/2.D0*Pi*im*epinv*CA - 1.D0/4.D0*Pi*im*beta0
     &  - 1.D0/24.D0*Pi**2*CA**(-1) + 1.D0/12.D0*Pi**2*CA - 1.D0/4.D0*
     & Pi**2*im**2*CA + 3.D0/8.D0*Log(u*musq**(-1))*CA + 1.D0/2.D0*Log(
     & u*musq**(-1))*epinv*CA + 1.D0/4.D0*Log(u*musq**(-1))*beta0 + 1.D0
     & /2.D0*Log(u*musq**(-1))*Pi*im*CA - 1.D0/4.D0*Log(u*musq**(-1))**
     & 2*CA + 3.D0/8.D0*Log( - t*musq**(-1))*CA + 1.D0/2.D0*Log( - t*
     & musq**(-1))*epinv*CA + 1.D0/4.D0*Log( - t*musq**(-1))*beta0 - 1.D
     & 0/4.D0*Log( - t*musq**(-1))**2*CA - 3.D0/4.D0*Log( - s*
     & musq**(-1))*CA**(-1) - 1.D0/2.D0*Log( - s*musq**(-1))*epinv*
     & CA**(-1) + 1.D0/4.D0*Log( - s*musq**(-1))**2*CA**(-1)
          else
             call abort
          endif

      end function

      end module
