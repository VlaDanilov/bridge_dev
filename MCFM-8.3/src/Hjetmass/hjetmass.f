      subroutine hjetmass(p,msq)
          use hjetmass_hel
          use hjetmass_highpt
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'constants.f'
      include 'masses.f'
      include 'hdecaymode.f'
      double precision s(mxpart,mxpart)
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'asymptotic.f'
      include 'first.f'
      include 'zprods_decl.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),gg,qg,gq,qq,hdecay
      real(dp) :: hhpt_qq, hhpt_gq, hhpt_qg, hhpt_gg
      double precision msq_tmp2(-nf:nf, -nf:nf)
      double precision ehsvm3,ehsvm4,s34
      double precision msqgamgam
      double precision gg_htl, gq_htl, qg_htl, qq_htl
      double precision dot

      double precision sman, tman, uman

      ! sushi amplitude defs
      double precision cex_prefactor, aex_prefactor
      double precision c1ex, c2ex, c3ex, c4ex
      double precision c1exSq, c2exSq, c3exSq, c4exSq
      double precision aexSq, fac, Asq
      double precision asmh_real, asmh_abs
      double precision c1smh_real, c2smh_real
      double precision c1smh_abs, c2smh_abs
      double complex c1smh, c2smh, asmh
      double precision sushi_gg, sushi_gq, sushi_qg, sushi_qq
      ! end sushi amplitude defs
      integer, parameter :: iglue = 5
      real(dp) :: pttwo, dotvec, etaraptwo, yraptwo
      real(dp) :: mesq_gg,mesq_gq

      complex(dp) :: ampggg(2,2,2)
      real(dp) :: ampgggsq
      complex(dp) :: ampqqg(2,2)
      real(dp) :: ampqqgsq
      complex(dp) :: ampsgg_mtex(3, 2,2,2)
      complex(dp) :: ampsqq_mtex(3, 2,2)

      call dotem(iglue,p,s)

      s34=(p(3,4)+p(4,4))**2
     & -(p(3,1)+p(4,1))**2-(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2

C   Deal with Higgs decay
      if (hdecaymode == 'tlta') then
          call htautaudecay(p,3,4,hdecay)
      elseif (hdecaymode == 'bqba') then
          call hbbdecay(p,3,4,hdecay)
      elseif (hdecaymode == 'gaga') then
          hdecay=msqgamgam(sqrt(s34))
      else
      write(6,*) 'Unimplemented process in hjetmass'
      stop
      endif

      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)


      call spinoru_dp(5,p,za,zb)

c ########### mt expansion amplitudes ##########

c     call hjetmass_ggg_ppp_1l_mtex(za,zb,1,2,iglue,ampsgg_mtex(:,1,1,1))
c     call hjetmass_ggg_ppp_1l_mtex(zb,za,1,2,iglue,ampsgg_mtex(:,2,2,2))

c     call hjetmass_ggg_ppm_1l_mtex(za,zb,1,2,iglue,ampsgg_mtex(:,1,1,2))
c     call hjetmass_ggg_ppm_1l_mtex(zb,za,1,2,iglue,ampsgg_mtex(:,2,2,1))

c     call hjetmass_ggg_pmp_1l_mtex(za,zb,1,2,iglue,ampsgg_mtex(:,1,2,1))
c     call hjetmass_ggg_pmp_1l_mtex(zb,za,1,2,iglue,ampsgg_mtex(:,2,1,2))

c     call hjetmass_ggg_pmp_1l_mtex(zb,za,2,1,iglue,ampsgg_mtex(:,1,2,2))
c     call hjetmass_ggg_pmp_1l_mtex(za,zb,2,1,iglue,ampsgg_mtex(:,2,1,1))

c     if (mtex >= 0) then
c       gg = sum(abs(ampsgg_mtex(1,:,:,:))**2)*hdecay/256d0*v*ca
c     endif
c     
c     if (mtex >= 2) then
c       gg = gg + sum(ampsgg_mtex(1,:,:,:)*conjg(ampsgg_mtex(2,:,:,:)))*hdecay/256d0*v*ca/mt**2
c       gg = gg + sum(conjg(ampsgg_mtex(1,:,:,:))*ampsgg_mtex(2,:,:,:))*hdecay/256d0*v*ca/mt**2
c     endif

c     if (mtex >= 4) then
c       gg = gg + sum(ampsgg_mtex(1,:,:,:)*conjg(ampsgg_mtex(3,:,:,:)))*hdecay/256d0*v*ca/mt**4
c       gg = gg + sum(conjg(ampsgg_mtex(1,:,:,:))*ampsgg_mtex(3,:,:,:))*hdecay/256d0*v*ca/mt**4
c       gg = gg + sum(conjg(ampsgg_mtex(2,:,:,:))*ampsgg_mtex(2,:,:,:))*hdecay/256d0*v*ca/mt**4
c     endif

c     call hjetmass_qqg_mpm_1l_mtex(za,zb,1,2,iglue,ampsqq_mtex(:,2,2)) 
c     call hjetmass_qqg_mpm_1l_mtex(zb,za,1,2,iglue,ampsqq_mtex(:,1,1)) 
c     call hjetmass_qqg_mpp_1l_mtex(za,zb,1,2,iglue,ampsqq_mtex(:,2,1)) 
c     call hjetmass_qqg_mpp_1l_mtex(zb,za,1,2,iglue,ampsqq_mtex(:,1,2)) 

c     if (mtex >= 0) then
c       qq = sum(abs(ampsqq_mtex(1,:,:))**2)*hdecay/36d0*v*tr
c     endif

c     if (mtex >= 2) then
c       qq = qq + sum(ampsqq_mtex(1,:,:)*conjg(ampsqq_mtex(2,:,:)))*hdecay/36d0*v*tr/mt**2
c       qq = qq + sum(ampsqq_mtex(2,:,:)*conjg(ampsqq_mtex(1,:,:)))*hdecay/36d0*v*tr/mt**2
c     endif

c     if (mtex >= 4) then
c       qq = qq + sum(ampsqq_mtex(1,:,:)*conjg(ampsqq_mtex(3,:,:)))*hdecay/36d0*v*tr/mt**4
c       qq = qq + sum(ampsqq_mtex(3,:,:)*conjg(ampsqq_mtex(1,:,:)))*hdecay/36d0*v*tr/mt**4
c       qq = qq + sum(ampsqq_mtex(2,:,:)*conjg(ampsqq_mtex(2,:,:)))*hdecay/36d0*v*tr/mt**4
c     endif

c     call hjetmass_qqg_mpm_1l_mtex(za,zb,1,iglue,2,ampsqq_mtex(:,2,2)) 
c     call hjetmass_qqg_mpm_1l_mtex(zb,za,1,iglue,2,ampsqq_mtex(:,1,1)) 
c     call hjetmass_qqg_mpp_1l_mtex(za,zb,1,iglue,2,ampsqq_mtex(:,2,1)) 
c     call hjetmass_qqg_mpp_1l_mtex(zb,za,1,iglue,2,ampsqq_mtex(:,1,2)) 

c     if (mtex >= 0) then
c       qg = sum(abs(ampsqq_mtex(1,:,:))**2)*hdecay/96d0*v*tr
c     endif

c     if (mtex >= 2) then
c       qg = qg + sum(ampsqq_mtex(1,:,:)*conjg(ampsqq_mtex(2,:,:)))*hdecay/96d0*v*tr/mt**2
c       qg = qg + sum(ampsqq_mtex(2,:,:)*conjg(ampsqq_mtex(1,:,:)))*hdecay/96d0*v*tr/mt**2
c     endif

c     if (mtex >= 4) then
c       qg = qg + sum(ampsqq_mtex(1,:,:)*conjg(ampsqq_mtex(3,:,:)))*hdecay/96d0*v*tr/mt**4
c       qg = qg + sum(ampsqq_mtex(3,:,:)*conjg(ampsqq_mtex(1,:,:)))*hdecay/96d0*v*tr/mt**4
c       qg = qg + sum(ampsqq_mtex(2,:,:)*conjg(ampsqq_mtex(2,:,:)))*hdecay/96d0*v*tr/mt**4
c     endif

c     call hjetmass_qqg_mpm_1l_mtex(za,zb,iglue,2,1,ampsqq_mtex(:,2,2)) 
c     call hjetmass_qqg_mpm_1l_mtex(zb,za,iglue,2,1,ampsqq_mtex(:,1,1)) 
c     call hjetmass_qqg_mpp_1l_mtex(za,zb,iglue,2,1,ampsqq_mtex(:,2,1)) 
c     call hjetmass_qqg_mpp_1l_mtex(zb,za,iglue,2,1,ampsqq_mtex(:,1,2)) 

c     if (mtex >= 0) then
c       gq = sum(abs(ampsqq_mtex(1,:,:))**2)*hdecay/96d0*v*tr
c     endif

c     if (mtex >= 2) then
c       gq = gq + sum(ampsqq_mtex(1,:,:)*conjg(ampsqq_mtex(2,:,:)))*hdecay/96d0*v*tr/mt**2
c       gq = gq + sum(ampsqq_mtex(2,:,:)*conjg(ampsqq_mtex(1,:,:)))*hdecay/96d0*v*tr/mt**2
c     endif

c     if (mtex >= 4) then
c       gq = gq + sum(ampsqq_mtex(1,:,:)*conjg(ampsqq_mtex(3,:,:)))*hdecay/96d0*v*tr/mt**4
c       gq = gq + sum(ampsqq_mtex(3,:,:)*conjg(ampsqq_mtex(1,:,:)))*hdecay/96d0*v*tr/mt**4
c       gq = gq + sum(ampsqq_mtex(2,:,:)*conjg(ampsqq_mtex(2,:,:)))*hdecay/96d0*v*tr/mt**4
c     endif

c ########### mt exact amplitudes ##########
      ampggg = 0._dp

      ampggg(1,1,1) = hjetmass_ggg_ppp(za,zb,1,2,iglue)
      ampggg(2,2,2) = hjetmass_ggg_ppp(zb,za,1,2,iglue)

      ampggg(1,1,2) = hjetmass_ggg_pmp(za,zb,1,iglue,2)
      ampggg(2,2,1) = hjetmass_ggg_pmp(zb,za,1,iglue,2)

      ampggg(1,2,1) = hjetmass_ggg_pmp(za,zb,1,2,iglue)
      ampggg(2,1,2) = hjetmass_ggg_pmp(zb,za,1,2,iglue)

      ampggg(1,2,2) = hjetmass_ggg_pmp(zb,za,2,1,iglue)
      ampggg(2,1,1) = hjetmass_ggg_pmp(za,zb,2,1,iglue)
      gg = sum(abs(ampggg)**2)*hdecay/256d0*v*ca

      ampqqg = 0._dp
      ampqqg(2,2) = hjetmass_qqg_mpm(za,zb,1,2,iglue)
      ampqqg(1,1) = hjetmass_qqg_mpm(zb,za,1,2,iglue)
      ampqqg(2,1) = hjetmass_qqg_mpp(za,zb,1,2,iglue)
      ampqqg(1,2) = hjetmass_qqg_mpp(zb,za,1,2,iglue)
      qq = sum(abs(ampqqg)**2)*hdecay/36._dp * v*tr

      ampqqg = 0._dp
      ampqqg(2,2) = hjetmass_qqg_mpm(za,zb,1,iglue,2)
      ampqqg(1,1) = hjetmass_qqg_mpm(zb,za,1,iglue,2)
      ampqqg(2,1) = hjetmass_qqg_mpp(za,zb,1,iglue,2)
      ampqqg(1,2) = hjetmass_qqg_mpp(zb,za,1,iglue,2)
      qg = sum(abs(ampqqg)**2)*hdecay/96._dp * v*tr

      ampqqg = 0._dp
      ampqqg(2,2) = hjetmass_qqg_mpm(za,zb,iglue,2,1)
      ampqqg(1,1) = hjetmass_qqg_mpm(zb,za,iglue,2,1)
      ampqqg(2,1) = hjetmass_qqg_mpp(za,zb,iglue,2,1)
      ampqqg(1,2) = hjetmass_qqg_mpp(zb,za,iglue,2,1)
      gq = sum(abs(ampqqg)**2)*hdecay/96._dp * v*tr



c ########### high pt amplitudes ##########
c     call hjetmass_highpt_lo(p,hhpt_qq,hhpt_gq,hhpt_qg,hhpt_gg)

c     hhpt_qq = hhpt_qq * hdecay * v*tr/36._dp
c     hhpt_gq = hhpt_gq * hdecay * v*tr/96._dp
c     hhpt_qg = hhpt_qg * hdecay * v*tr/96._dp
c     hhpt_gg = hhpt_gg * hdecay * v*ca/256._dp

c     qq = hhpt_qq
c     gq = hhpt_gq
c     qg = hhpt_qg
c     gg = hhpt_gg


c ########### old mt expansion amplitudes ##########
c     sman = s(1,2)
c     tman = s(1,iglue)
c     uman = s(2,iglue)

c     write (*,*) gg/(mesq_gg(s(1,2),s(1,5),s(2,5),mtex) * hdecay)
c     write (*,*) qg/(mesq_gq(s(1,2),s(2,5),s(1,5),mtex) * hdecay)
c     write (*,*) gq/(mesq_gq(s(1,2),s(1,5),s(2,5),mtex) * hdecay)
c     write (*,*) qq/(-8d0/3d0*mesq_gq(s(5,2),s(5,1),s(2,1),mtex) * hdecay)
c     write (*,*) ""


c ########### old mt exact amplitudes ##########

c     !! amplitude from SusHi
c     cex_prefactor = sqrt(gsq)**3 / 4d0 / pi**2 / sqrt(vevsq)

c     if (first .eqv. .true.) then
c       call sushi_bernini(18)
c     else
c       first =  .false. 
c     endif

c     c1exSq = (c1smh_abs(sman,tman,uman,mt**2) *
c    & cex_prefactor * mt**2 / 16d0)**2

c     c2exSq = (c2smh_abs(sman,tman,uman,mt**2) *
c    & cex_prefactor * mt**2 / 16d0)**2

c     c3exSq = (c2smh_abs(tman,sman,uman,mt**2) *
c    & cex_prefactor * mt**2 / 16d0)**2

c     c4exSq = (c2smh_abs(uman,sman,tman,mt**2) *
c    & cex_prefactor * mt**2 / 16d0)**2

c     sushi_gg = (c1exSq + c2exSq + c3exSq + c4exSq)*
c    &              v*ca*16/sman/tman/uman
c     gg =  sushi_gg*hdecay/256d0

c     write (*,*) "ratio: ", ampgggsq/gg
c     write (*,*) ""
c     pause

c     aex_prefactor = sqrt(gsq)**3 / 32d0 / pi**2 / sqrt(vevsq)
c     aexSq = (asmh_abs(tman,sman,uman,mt**2) / tman * aex_prefactor *
c    &             mt**2 * 8d0/3d0)**2
c     sushi_qg = aexSq*8d0/2d0 * (-uman**2 * tman - sman**2 * tman)
c     qg =  sushi_qg*hdecay/96d0

c     aexSq = (asmh_abs(uman,sman,tman,mt**2) / uman * aex_prefactor *
c    &             mt**2 * 8d0/3d0)**2
c     sushi_gq = aexSq*8d0/2d0 * (-tman**2 * uman - sman**2 * uman)
c     gq = sushi_gq*hdecay/96d0

c     aexSq = (asmh_abs(sman,uman,tman,mt**2) / sman * aex_prefactor *
c    &             mt**2 * 8d0/3d0)**2
c     sushi_qq = aexSq*8d0/2d0 * (-tman**2 * sman - uman**2 * sman)
c     qq = sushi_qq*hdecay/(-36d0)

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0

      if ((j.eq. 0) .or. (k.eq.0)) then
           if ((j.eq. 0) .and. (k.eq.0)) then
                msq(j,k)=gg
           elseif ((j.eq.0).and.(k.ne.0)) then
                msq(j,k)=gq
           elseif ((j.ne.0).and.(k.eq.0)) then
                msq(j,k)=qg
           endif
      elseif ((j.eq.-k).and. (j.ne.0)) then
             msq(j,k)=qq
      endif

      enddo
      enddo

      end

c     double precision function mesq_gg(sman,tman,uman,mtex)
c       implicit none
c       double precision, intent(in) :: sman,tman,uman
c       integer, intent(in) :: mtex
c       include 'types.f'
c       include 'constants.f'
c       include 'qcdcouple.f'
c       include 'ewcouple.f'
c       include 'masses.f'

c       double precision m1
c       m1 = mt

c       if (mtex >= 0) then
c         mesq_gg =
c    &  + 2.D0*sman**(-1)*tman**(-1)*uman**3 + 4.D0*sman**(-1)*uman**2
c    &     + 6.D0*sman**(-1)*tman*uman + 4.D0*sman**(-1)*tman**2 + 2.D0
c    &    *sman**(-1)*tman**3*uman**(-1) + 4.D0*tman**(-1)*uman**2 + 12.
c    &    D0*uman + 12.D0*tman + 4.D0*tman**2*uman**(-1) + 6.D0*sman*
c    &    tman**(-1)*uman + 12.D0*sman + 6.D0*sman*tman*uman**(-1) + 4.D
c    &    0*sman**2*tman**(-1) + 4.D0*sman**2*uman**(-1) + 2.D0*sman**3
c    &    *tman**(-1)*uman**(-1)
c       endif

c       if (mtex >= 2) then
c        mesq_gg = mesq_gg
c    &  + 7.D0/30.D0*M1**(-2)*sman**(-1)*tman**(-1)*uman**4 + 7.D0/10.D0
c    &    *M1**(-2)*sman**(-1)*uman**3 + 7.D0/6.D0*M1**(-2)*sman**(-1)*
c    &    tman*uman**2 + 7.D0/6.D0*M1**(-2)*sman**(-1)*tman**2*uman + 7.
c    &    D0/10.D0*M1**(-2)*sman**(-1)*tman**3 + 7.D0/30.D0*M1**(-2)*
c    &    sman**(-1)*tman**4*uman**(-1) + 7.D0/10.D0*M1**(-2)*
c    &    tman**(-1)*uman**3 + 67.D0/30.D0*M1**(-2)*uman**2 + 33.D0/10.D
c    &    0*M1**(-2)*tman*uman + 67.D0/30.D0*M1**(-2)*tman**2 + 7.D0/10.
c    &    D0*M1**(-2)*tman**3*uman**(-1) + 7.D0/6.D0*M1**(-2)*sman*
c    &    tman**(-1)*uman**2 + 33.D0/10.D0*M1**(-2)*sman*uman + 33.D0/
c    &    10.D0*M1**(-2)*sman*tman + 7.D0/6.D0*M1**(-2)*sman*tman**2*
c    &    uman**(-1) + 7.D0/6.D0*M1**(-2)*sman**2*tman**(-1)*uman + 67.D
c    &    0/30.D0*M1**(-2)*sman**2 + 7.D0/6.D0*M1**(-2)*sman**2*tman*
c    &    uman**(-1) + 7.D0/10.D0*M1**(-2)*sman**3*tman**(-1) + 7.D0/10.
c    &    D0*M1**(-2)*sman**3*uman**(-1) + 7.D0/30.D0*M1**(-2)*sman**4*
c    &    tman**(-1)*uman**(-1)
c       endif

c       if (mtex >= 4) then
c        mesq_gg = mesq_gg
c    &  + 1543.D0/50400.D0*M1**(-4)*sman**(-1)*tman**(-1)*uman**5 + 
c    &    1543.D0/12600.D0*M1**(-4)*sman**(-1)*uman**4 + 1543.D0/6300.D0
c    &    *M1**(-4)*sman**(-1)*tman*uman**3 + 1543.D0/5040.D0*M1**(-4)*
c    &    sman**(-1)*tman**2*uman**2 + 1543.D0/6300.D0*M1**(-4)*
c    &    sman**(-1)*tman**3*uman + 1543.D0/12600.D0*M1**(-4)*
c    &    sman**(-1)*tman**4 + 1543.D0/50400.D0*M1**(-4)*sman**(-1)*
c    &    tman**5*uman**(-1) + 1543.D0/12600.D0*M1**(-4)*tman**(-1)*
c    &    uman**4 + 12197.D0/25200.D0*M1**(-4)*uman**3 + 3739.D0/4200.D0
c    &    *M1**(-4)*tman*uman**2 + 3739.D0/4200.D0*M1**(-4)*tman**2*
c    &    uman + 12197.D0/25200.D0*M1**(-4)*tman**3 + 1543.D0/12600.D0*


c    &    M1**(-4)*tman**4*uman**(-1) + 1543.D0/6300.D0*M1**(-4)*sman*
c    &    tman**(-1)*uman**3 + 3739.D0/4200.D0*M1**(-4)*sman*uman**2 + 
c    &    1059.D0/800.D0*M1**(-4)*sman*tman*uman + 3739.D0/4200.D0*
c    &    M1**(-4)*sman*tman**2 + 1543.D0/6300.D0*M1**(-4)*sman*tman**3
c    &    *uman**(-1)
c         mesq_gg = mesq_gg +
c    &    1543.D0/5040.D0*M1**(-4)*sman**2*tman**(-1)*
c    &    uman**2 + 3739.D0/4200.D0*M1**(-4)*sman**2*uman 
c    &    + 3739.D0/4200.D0


c    &    *M1**(-4)*sman**2*tman + 1543.D0/5040.D0*M1**(-4)*sman**2*
c    &    tman**2*uman**(-1) + 1543.D0/6300.D0*M1**(-4)*sman**3*
c    &    tman**(-1)*uman + 12197.D0/25200.D0*M1**(-4)*sman**3 + 1543.D0
c    &    /6300.D0*M1**(-4)*sman**3*tman*uman**(-1) + 1543.D0/12600.D0*
c    &    M1**(-4)*sman**4*tman**(-1) + 1543.D0/12600.D0*M1**(-4)*
c    &    sman**4*uman**(-1) + 1543.D0/50400.D0*M1**(-4)*sman**5*
c    &    tman**(-1)*uman**(-1)
c       endif

c       mesq_gg = mesq_gg/256d0*gsq*24d0*(as/3d0/pi)**2/vevsq

c     end function

c     double precision function mesq_gq(sman,tman,uman,mtex)
c       implicit none
c       double precision, intent(in) :: sman,tman,uman
c       integer, intent(in) :: mtex
c       include 'types.f'
c       include 'constants.f'
c       include 'qcdcouple.f'
c       include 'ewcouple.f'
c       include 'masses.f'

c       double precision m1
c       m1 = mt

c       if (mtex >= 0) then
c       mesq_gq = - tman**2*uman**(-1) - sman**2*uman**(-1)
c       endif

c       if (mtex >= 2) then
c         mesq_gq =  mesq_gq
c    &  + M1**(-2) * (  - 3./10.*tman**2 - 7./60.*tman**3*uman**(-1) - 
c    &    7./60.*sman*tman**2*uman**(-1) - 3./10.*sman**2 - 7./60.*
c    &    sman**2*tman*uman**(-1) - 7./60.*sman**3*uman**(-1) )
c       endif


c       if (mtex >= 4) then
c       mesq_gq =  mesq_gq
c    &  + M1**(-4) * (  - 223./2800.*tman**2*uman - 169./2800.*tman**3
c    &     - 1543./100800.*tman**4*uman**(-1) - 169./2800.*sman*tman**2
c    &     - 1543./50400.*sman*tman**3*uman**(-1) - 223./2800.*sman**2*
c    &    uman - 169./2800.*sman**2*tman - 1543./50400.*sman**2*tman**2
c    &    *uman**(-1) - 169./2800.*sman**3 - 1543./50400.*sman**3*tman*
c    &    uman**(-1) - 1543./100800.*sman**4*uman**(-1) )
c       endif
c       
c       mesq_gq = mesq_gq/4d0/3d0/8d0*gsq*4*(as/3d0/pi)**2/vevsq

c     end function

