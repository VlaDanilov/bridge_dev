      subroutine hjetmass_v(p,msq)
        use hjetmass_hel
        use hjetmass_highpt
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'constants.f'
      include 'masses.f'
      include 'hdecaymode.f'
      include 'sprods_com.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'scheme.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scale.f'
      include 'asymptotic.f'
      include 'zprods_decl.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),gg,qg,gq,hdecay
      double precision gg_new, gq_old, gg_mix
      real(dp) :: pt, pttwo
      double precision gq_new, gq_mix
      double precision s34, sman, tman, uman
      double precision msqgamgam
      double precision virtme_gg, virtme_gq, qa, aq
      double precision virtme_gq_test
      double precision virtme_gq_new, virtme_gg_multi
      double precision mesq_gg, mesq_gq
      double precision msqlo(-nf:nf, -nf:nf)
      double precision msqv(-nf:nf, - nf:nf)
      double precision virtgg, virtqa, virtaq, virtqg, virtgq, Asq, fac
      double precision gg_insert, gq_insert, dummy
      double precision x1000, x100, x10, ep2, ep1, fin, epsave
      double precision ep2_mcfm, ep1_mcfm, fin_mcfm
      double precision insarray(-nf:nf,-nf:nf)
      double precision insarray_new(-nf:nf,-nf:nf)
      double precision insarray_fixedscale(-nf:nf, -nf:nf)

      real(dp) :: virt_gg, virt_gq, virt_qg, virt_qq
      real(dp) :: virt_gg_hhpt, virt_gq_hhpt, virt_qg_hhpt, virt_qq_hhpt

      integer, parameter :: iglue = 5

      complex(dp) :: ampsqq_mtex(2, 2,2)
      complex(dp) :: ampsqq_hhpt(2, 2,2)
      complex(dp) :: ampsgq_hhpt(2, 2,2)
      complex(dp) :: ampsqg_hhpt(2, 2,2)
      complex(dp) :: ampsqq_lo(2,2)
      complex(dp) :: ampsgq_lo(2,2)
      complex(dp) :: ampsqg_lo(2,2)

      complex(dp) :: ampsgg_mtex(2, 2,2,2)
      complex(dp) :: ampsgg_hhpt(2, 2,2,2)
      complex(dp) :: ampsgg_lo(2,2,2)

      complex(dp) :: omega2l(2,3)

      complex(dp) :: kfac_qq(2,2), kfac_gq(2,2), kfac_qg(2,2)
      complex(dp) :: kfac_gg(2,2,2)

      real(dp) :: dotvec
      real(dp) :: epinv_sav,epinv2_sav,musq_sav

      s34=(p(3,4)+p(4,4))**2
     & -(p(3,1)+p(4,1))**2-(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2

      ampsqq_mtex = 0._dp
      ampsqq_hhpt = 0._dp
      ampsqq_lo = 0._dp

      ampsgg_mtex = 0._dp
      ampsgg_hhpt = 0._dp
      ampsgg_lo = 0._dp

      omega2l = 0._dp

      scheme = 'tH-V'

C   Deal with Higgs decay
      if (hdecaymode == 'tlta') then
          call htautaudecay(p,3,4,hdecay)
      elseif (hdecaymode == 'bqba') then
          call hbbdecay(p,3,4,hdecay)
      elseif (hdecaymode == 'gaga') then
          hdecay=msqgamgam(hmass)
      else
      write(6,*) 'Unimplemented process in hjetmass'
      stop
      endif

      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)

c     !! old 1/mt expansion amps
c     sman = s(1,2)
c     tman = s(1,5)
c     uman = s(2,5)

c     gg = virtme_gg_multi(sman,tman,uman,mtex) * hdecay
c     gq = virtme_gq_new(uman,tman,sman,mtex) * hdecay
c     qg = virtme_gq_new(tman,uman,sman,mtex) * hdecay
c     aq = -8d0/3d0 * virtme_gq_new(sman,uman,tman,mtex) * hdecay
c     qa = -8d0/3d0 * virtme_gq_new(sman,tman,uman,mtex) * hdecay
c     !!

c     !! new helicity amplitudes

        call spinoru_dp(5,p,za,zb)

      sman = dotvec(p(1,:)+p(2,:), p(1,:)+p(2,:))
      tman = dotvec(p(1,:)+p(5,:), p(1,:)+p(5,:))
      uman = dotvec(p(2,:)+p(5,:), p(2,:)+p(5,:))

      if (mtex == 100) then

        call hjetmass_highpt_amps(p,ampsqq_hhpt,ampsgq_hhpt,ampsqg_hhpt,ampsgg_hhpt)

        ampsqq_lo(2,2) = -im*hjetmass_qqg_mpm(za,zb,1,2,iglue)
        ampsqq_lo(1,1) = -im*hjetmass_qqg_mpm(zb,za,1,2,iglue)
        ampsqq_lo(2,1) = im*hjetmass_qqg_mpp(za,zb,1,2,iglue)
        ampsqq_lo(1,2) = im*hjetmass_qqg_mpp(zb,za,1,2,iglue)

        ampsqq_hhpt(2,:,:) = ampsqq_hhpt(2,:,:) + I1q(sman,tman,uman,musq)*ampsqq_lo(:,:)*(as/2d0/pi)

        virt_qq_hhpt = sum(2._dp*real(conjg(ampsqq_lo(:,:))*ampsqq_hhpt(2,:,:)))*hdecay/36._dp*v*tr
c       write (*,*) "pttwo", pttwo(3,4,p)
c       write (*,*) "amp comparison qq"
c       write (*,*) ampsqq_hhpt(1,:,:)/ampsqq_lo(:,:)


        ampsqg_lo(2,2) = -im*hjetmass_qqg_mpm(za,zb,1,iglue,2)
        ampsqg_lo(1,1) = -im*hjetmass_qqg_mpm(zb,za,1,iglue,2)
        ampsqg_lo(2,1) = im*hjetmass_qqg_mpp(za,zb,1,iglue,2)
        ampsqg_lo(1,2) = im*hjetmass_qqg_mpp(zb,za,1,iglue,2)

        ampsqg_hhpt(2,:,:) = ampsqg_hhpt(2,:,:) + I1q(tman,uman,sman,musq)*ampsqg_lo(:,:)*(as/2d0/pi)

        virt_qg_hhpt = sum(2._dp*real(conjg(ampsqg_lo(:,:))*ampsqg_hhpt(2,:,:)))*hdecay/96._dp*v*tr

        !write (*,*) "pttwo", pttwo(3,4,p)
        !write (*,*) "amp comparison qg"
        !write (*,*) ampsqg_hhpt(1,:,:)/ampsqg_lo(:,:)


        ampsgq_lo(2,2) = -im*hjetmass_qqg_mpm(za,zb,iglue,2,1)
        ampsgq_lo(1,1) = -im*hjetmass_qqg_mpm(zb,za,iglue,2,1)
        ampsgq_lo(2,1) = im*hjetmass_qqg_mpp(za,zb,iglue,2,1)
        ampsgq_lo(1,2) = im*hjetmass_qqg_mpp(zb,za,iglue,2,1)

        ampsgq_hhpt(2,:,:) = ampsgq_hhpt(2,:,:) + I1q(uman,tman,sman,musq)*ampsgq_lo(:,:)*(as/2d0/pi)

        virt_gq_hhpt = sum(2._dp*real(conjg(ampsgq_lo(:,:))*ampsgq_hhpt(2,:,:)))*hdecay/96._dp*v*tr

        !write (*,*) "pttwo", pttwo(3,4,p)
        !write (*,*) "amp comparison gq"
        !write (*,*) ampsgq_hhpt(1,:,:)/ampsgq_lo(:,:)

        ampsgg_lo(1,1,1) = hjetmass_ggg_ppp(za,zb,1,2,iglue)
        ampsgg_lo(2,2,2) = hjetmass_ggg_ppp(zb,za,1,2,iglue)
        ampsgg_lo(1,1,2) = -hjetmass_ggg_pmp(za,zb,1,iglue,2)
        ampsgg_lo(2,2,1) = -hjetmass_ggg_pmp(zb,za,1,iglue,2)
        ampsgg_lo(1,2,1) = -hjetmass_ggg_pmp(za,zb,1,2,iglue)
        ampsgg_lo(2,1,2) = -hjetmass_ggg_pmp(zb,za,1,2,iglue)
        ampsgg_lo(1,2,2) = -hjetmass_ggg_pmp(zb,za,2,1,iglue)
        ampsgg_lo(2,1,1) = -hjetmass_ggg_pmp(za,zb,2,1,iglue)

        !write (*,*) "pttwo", pttwo(3,4,p)
        !write (*,*) "amp comparison"
        !write (*,*) ampsgg_hhpt(1,:,:,:)/ampsgg_lo(:,:,:)

        ampsgg_hhpt(2,:,:,:) = ampsgg_hhpt(2,:,:,:) + I1g(sman,tman,uman,musq)*ampsgg_lo(:,:,:)*(as/2d0/pi)

        virt_gg_hhpt = sum(2._dp*real(conjg(ampsgg_lo(:,:,:))*ampsgg_hhpt(2,:,:,:)))*hdecay/256._dp*v*ca

        qa = virt_qq_hhpt
        aq = virt_qq_hhpt
        qg = virt_qg_hhpt
        gq = virt_gq_hhpt
        gg = virt_gg_hhpt
      else

        call hjetmass_qqg_mpm_2l_mtex(za,zb,1,2,iglue,ampsqq_mtex(:,2,2)) 
        ampsqq_lo(2,2) = hjetmass_qqg_mpm(za,zb,1,2,iglue)
        call hjetmass_qqg_mpm_2l_mtex(zb,za,1,2,iglue,ampsqq_mtex(:,1,1)) 
        ampsqq_lo(1,1) = hjetmass_qqg_mpm(zb,za,1,2,iglue)
        call hjetmass_qqg_mpp_2l_mtex(za,zb,1,2,iglue,ampsqq_mtex(:,2,1)) 
        ampsqq_lo(2,1) = hjetmass_qqg_mpp(za,zb,1,2,iglue)
        call hjetmass_qqg_mpp_2l_mtex(zb,za,1,2,iglue,ampsqq_mtex(:,1,2)) 
        ampsqq_lo(1,2) = hjetmass_qqg_mpp(zb,za,1,2,iglue)

        if (mtex == 0) then
            kfac_qq = ampsqq_lo(:,:)/ampsqq_mtex(1,:,:)
        else
            kfac_qq = 1._dp
        endif

        virt_qq = sum(2._dp*real(conjg(ampsqq_lo(:,:))*kfac_qq*ampsqq_mtex(2,:,:)))*hdecay/36._dp*v*tr

        call hjetmass_qqg_mpm_2l_mtex(za,zb,1,iglue,2,ampsqq_mtex(:,2,2)) 
        ampsqq_lo(2,2) = hjetmass_qqg_mpm(za,zb,1,iglue,2)
        call hjetmass_qqg_mpm_2l_mtex(zb,za,1,iglue,2,ampsqq_mtex(:,1,1)) 
        ampsqq_lo(1,1) = hjetmass_qqg_mpm(zb,za,1,iglue,2)
        call hjetmass_qqg_mpp_2l_mtex(za,zb,1,iglue,2,ampsqq_mtex(:,2,1)) 
        ampsqq_lo(2,1) = hjetmass_qqg_mpp(za,zb,1,iglue,2)
        call hjetmass_qqg_mpp_2l_mtex(zb,za,1,iglue,2,ampsqq_mtex(:,1,2)) 
        ampsqq_lo(1,2) = hjetmass_qqg_mpp(zb,za,1,iglue,2)


        if (mtex == 0) then
            kfac_qg = ampsqq_lo(:,:)/ampsqq_mtex(1,:,:)
        else
            kfac_qg = 1._dp
        endif

        virt_qg = sum(2._dp*real(conjg(ampsqq_lo(:,:))*kfac_qg*ampsqq_mtex(2,:,:)))*hdecay/96._dp*v*tr

        call hjetmass_qqg_mpm_2l_mtex(za,zb,iglue,2,1,ampsqq_mtex(:,2,2)) 
        ampsqq_lo(2,2) = hjetmass_qqg_mpm(za,zb,iglue,2,1)
        call hjetmass_qqg_mpm_2l_mtex(zb,za,iglue,2,1,ampsqq_mtex(:,1,1)) 
        ampsqq_lo(1,1) = hjetmass_qqg_mpm(zb,za,iglue,2,1)
        call hjetmass_qqg_mpp_2l_mtex(za,zb,iglue,2,1,ampsqq_mtex(:,2,1)) 
        ampsqq_lo(2,1) = hjetmass_qqg_mpp(za,zb,iglue,2,1)
        call hjetmass_qqg_mpp_2l_mtex(zb,za,iglue,2,1,ampsqq_mtex(:,1,2)) 
        ampsqq_lo(1,2) = hjetmass_qqg_mpp(zb,za,iglue,2,1)

        if (mtex == 0) then
            kfac_gq = ampsqq_lo(:,:)/ampsqq_mtex(1,:,:)
        else
            kfac_gq = 1._dp
        endif

        virt_gq = sum(2._dp*real(conjg(ampsqq_lo(:,:))*kfac_gq*ampsqq_mtex(2,:,:)))*hdecay/96._dp*v*tr

        call hjetmass_ggg_ppp_2l_mtex(za,zb,1,2,iglue,ampsgg_mtex(:,1,1,1))
        ampsgg_lo(1,1,1) = hjetmass_ggg_ppp(za,zb,1,2,iglue)
        call hjetmass_ggg_ppp_2l_mtex(zb,za,1,2,iglue,ampsgg_mtex(:,2,2,2))
        ampsgg_lo(2,2,2) = hjetmass_ggg_ppp(zb,za,1,2,iglue)

        call hjetmass_ggg_ppm_2l_mtex(za,zb,1,2,iglue,ampsgg_mtex(:,1,1,2))
        ampsgg_lo(1,1,2) = hjetmass_ggg_pmp(za,zb,1,iglue,2)
        call hjetmass_ggg_ppm_2l_mtex(zb,za,1,2,iglue,ampsgg_mtex(:,2,2,1))
        ampsgg_lo(2,2,1) = hjetmass_ggg_pmp(zb,za,1,iglue,2)

        call hjetmass_ggg_pmp_2l_mtex(za,zb,1,2,iglue,ampsgg_mtex(:,1,2,1))
        ampsgg_lo(1,2,1) = hjetmass_ggg_pmp(za,zb,1,2,iglue)
        call hjetmass_ggg_pmp_2l_mtex(zb,za,1,2,iglue,ampsgg_mtex(:,2,1,2))
        ampsgg_lo(2,1,2) = hjetmass_ggg_pmp(zb,za,1,2,iglue)

        call hjetmass_ggg_pmp_2l_mtex(zb,za,2,1,iglue,ampsgg_mtex(:,1,2,2))
        ampsgg_lo(1,2,2) = hjetmass_ggg_pmp(zb,za,2,1,iglue)
        call hjetmass_ggg_pmp_2l_mtex(za,zb,2,1,iglue,ampsgg_mtex(:,2,1,1))
        ampsgg_lo(2,1,1) = hjetmass_ggg_pmp(za,zb,2,1,iglue)

        if (mtex == 0) then
            kfac_gg = ampsgg_lo(:,:,:)/ampsgg_mtex(1,:,:,:)
        else
            kfac_gg = 1._dp
        endif

        virt_gg = sum(2._dp*real(conjg(ampsgg_lo(:,:,:))*kfac_gg*ampsgg_mtex(2,:,:,:)))*hdecay/256._dp*v*ca

        qa = -virt_qq
        aq = -virt_qq
        qg = -virt_qg
        gq = -virt_gq
        gg = virt_gg
      endif


      insarray_new = 0d0
      if (mtex /= 100) then
        call insert_exact_new(p,insarray_new)
      endif

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0 

      if ((j.eq.0).and.(k.eq.0)) then
        msq(j,k) = (gg - insarray_new(0,0))
      endif
      if ((j.gt.0).and.(k.eq.-j)) then
        msq(j,k) = (qa - insarray_new(j,k))
      endif
      if ((j.lt.0).and.(k.eq.-j)) then
        msq(j,k) = (aq - insarray_new(j,k))
      endif
      if ((j.eq.0).and.(k.ne.0)) then
        msq(j,k) = (gq - insarray_new(j,k))
      endif
      if ((j.ne.0).and.(k.eq.0)) then
        msq(j,k) = (qg - insarray_new(j,k))
      endif

      enddo
      enddo

      end

      subroutine insert_exact_new(p, msqlo)
        implicit none
        include 'types.f'
        include 'mxpart.f'
        include 'nf.f'
        include 'constants.f'
        include 'masses.f'
        include 'scale.f'
        include 'epinv.f'
        include 'epinv2.f'
        include 'qcdcouple.f'
        include 'ewcouple.f'
        include 'sprods_com.f'
        double precision, intent(in) :: p(mxpart,4)
        double precision, intent(out) :: msqlo(-nf:nf, -nf:nf)
        double precision sman, tman, uman
        double complex insert_gg
        double complex cdlogwrap
        double precision kg, mh, nl, pcut, insgq_new, bf
        double complex LogAbsMus, LogAbsMut, LogAbsMuu
        integer j,k

        sman = s(1,2)
        tman = s(1,5)
        uman = s(2,5)

        mh = hmass

        call hjetmass(p,msqlo)

        nl = 5d0
        pcut = 1d0
        bf = (11d0/6d0*ca - 2d0/3d0*tr*nl)
        kg = ca*(67d0/18d0 - pisqo6) - 10d0/9d0 * tr*nl
     &          - ca*dlog(pcut)**2
     &          + (11d0/6d0*ca - 2d0/3d0*tr*nl)*(pcut - 1 - dlog(pcut))

        LogAbsMus = cdlogwrap(-musq/sman)
        LogAbsMut = cdlogwrap(-musq/tman)
        LogAbsMuu = cdlogwrap(-musq/uman)

      insert_gg =
     & 3.D0*Kg + 3.D0*bf+ 3.D0*bf*epinv + bf*LogAbsMus + 
     & bf*LogAbsMut + bf*LogAbsMuu - cA*pi**2 + 3.D0*cA*epinv2*epinv
     &  + cA*LogAbsMus*epinv + 1.D0/2.D0*cA*LogAbsMus**2 + cA*LogAbsMut
     & *epinv + 1.D0/2.D0*cA*LogAbsMut**2 + cA*LogAbsMuu*epinv + 1.D0/2.
     & D0*cA*LogAbsMuu**2

      do j=-nf,nf
      do k=-nf,nf

      if ((j.eq.0).and.(k.eq.0)) then
        msqlo(j,k) = msqlo(j,k)*dreal(insert_gg)
      endif
      if ((j.gt.0).and.(k.eq.-j)) then
        msqlo(j,k) = msqlo(j,k) * insgq_new(sman,tman,uman)
      endif
      if ((j.lt.0).and.(k.eq.-j)) then
        msqlo(j,k) = msqlo(j,k) * insgq_new(sman,uman,tman)
      endif
      if ((j.eq.0).and.(k.ne.0)) then
        msqlo(j,k) = msqlo(j,k) * insgq_new(uman,tman,sman)
      endif
      if ((j.ne.0).and.(k.eq.0)) then
        msqlo(j,k) = msqlo(j,k) * insgq_new(tman,uman,sman)
      endif

      ! convert prefactor from lo to virt
      msqlo(j,k) = msqlo(j,k) * as/2/pi

      enddo
      enddo
    
      end subroutine insert_exact_new

      double precision function insgq_new(sman, tman, uman)
        implicit none
        include 'types.f'
        include 'constants.f'
        include 'masses.f'
        include 'scale.f'
        include 'epinv.f'
        include 'qcdcouple.f'
        include 'ewcouple.f'
        double precision kg, kq
        double precision, intent(in) :: sman, tman, uman
        double precision kquark, kgluon, pcut, nl

        double complex LogMums, LogMumt, LogMumu
        double complex insop
        double complex cdlogwrap


        nl = 5d0
        pcut = 1d0
        kg = ca*(67d0/18d0 - pisqo6) - 10d0/9d0 * tr*nl
     &          - ca*dlog(pcut)**2
     &          + (11d0/6d0*ca - 2d0/3d0*tr*nl)*(pcut - 1 - dlog(pcut))
        kq = cf*(7d0/2d0 - pisqo6)
     &          - cf*dlog(pcut)**2
     &          + 3d0/2d0*cf*(pcut - 1 - dlog(pcut))

        kquark = kq
        kgluon = kg

        LogMumt = cdlogwrap(-musq/tman)
        LogMumu = cdlogwrap(-musq/uman)
        LogMums = cdlogwrap(-musq/sman)

      insop =
     &  + epinv * ( 3.D0*CF + 2.D0*CF*LogMums - 2.D0/3.D0*nl*tr + 11.D0/
     &    6.D0*cA + cA*LogMumt + cA*LogMumu - cA*LogMums )
      insop = insop + epinv**2 * ( 2.D0*CF + cA )
      insop = insop + 2.D0*Kquark + Kgluon + 3.D0*CF + 3.D0*CF*LogMums
     &     + CF*LogMums**2 - 2.D0/3.D0*CF*pi**2 - 2.D0/3.D0*nl*tr - 1.D0
     &    /3.D0*nl*tr*LogMumt - 1.D0/3.D0*nl*tr*LogMumu + 11.D0/6.D0*cA
     &     + 5.D0/3.D0*cA*LogMumt + 1.D0/2.D0*cA*LogMumt**2 + 5.D0/3.D0
     &    *cA*LogMumu + 1.D0/2.D0*cA*LogMumu**2 - 3.D0/2.D0*cA*LogMums
     &     - 1.D0/2.D0*cA*LogMums**2 - 1.D0/3.D0*cA*pi**2

      insgq_new = dreal(insop)

      end function insgq_new


