      subroutine hjetmass_r(p_mcfm,msq)
          use hjetmass_r_amps_qp_mod
          use hjetmass_r_amps_dp_mod
      implicit none
c---Matrix element squared averaged over initial colors and spins
c
c     g(-p1)+g(-p2) -->  H(p3)+g(p_iglue1=5)+g(p_iglue2=6) 
      include 'types.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      include 'msq_struc.f'
      include 'nflav.f'
      include 'hdecaymode.f'
      include 'scale.f'
      include 'badpoint.f'
      include 'noglue.f'

      real(dp), intent(in) :: p_mcfm(mxpart,4)
      real(dp), intent(out) ::  msq(-nf:nf,-nf:nf)
      
      real(dp) :: hdecay,s34

      integer, parameter :: iglue1=5, iglue2=6

      real(dp) :: Asq,fac
      real(qp) :: p_qp(mxpart,4)
      real(dp):: p_dp(mxpart,4)
   
      real(dp) :: massvec, msqgamgam

      real(dp) :: Hqagg,Haqgg,Hgqqg,Hgaag,Hqgqg,Hagag
      real(dp) :: Hqgqg_ab,Hqgqg_ba,Hqgqg_sym
      real(dp) :: Hgqqg_ab,Hgqqg_ba,Hgqqg_sym
      real(dp) :: Hagag_ab,Hagag_ba,Hagag_sym
      real(dp) :: Hgaag_ab,Hgaag_ba,Hgaag_sym
      real(dp) :: Hqagg_ab,Hqagg_ba,Hqagg_sym
      real(dp) :: Haqgg_ab,Haqgg_ba,Haqgg_sym

      real(dp) :: Hgggg,Hgggg_1256,Hgggg_1265,Hgggg_1625
      real(dp) :: Hggqa, Hggqa_ab,Hggqa_ba,Hggqa_sym

      real(dp) :: Hqqqq, Hqqqq_a,Hqqqq_b,Hqqqq_i
      real(dp) :: Hqaqa, Hqaqa_a,Hqaqa_b,Hqaqa_i
      real(dp) :: Haqaq, Haqaq_a,Haqaq_b,Haqaq_i
      real(dp) :: Hqaaq, Hqaaq_a,Hqaaq_b,Hqaaq_i

      real(dp) :: Hqrqr,Habab,Haaaa,Hqarb,Hqbqb,Haqbr,Hbqbq

      logical ggchan, gqchan, qqchan

      logical :: passed_scaletest
      logical :: ret
      real(dp) :: dum1,dum2,dum3
      integer :: j,k

      integer i1,i2,i3
      real(dp) :: t
      t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

      ggchan = .true.
      gqchan = .true.
      qqchan = .true.

      if (noglue .eqv. .true.) then
        ggchan = .false.
        gqchan = .false.
        qqchan = .true.
      end if

      if (gqonly .eqv. .true.) then
        ggchan = .false.
        gqchan = .true.
        qqchan = .false.
      end if

      if (ggonly .eqv. .true.) then
        ggchan = .true.
        gqchan = .false.
        qqchan = .false.
      end if

      Hqagg = 0d0; Haqgg = 0d0; Hgqqg = 0d0
      Hgaag = 0d0; Hqgqg = 0d0; Hagag = 0d0
      Hqgqg_ab = 0d0; Hqgqg_ba = 0d0; Hqgqg_sym = 0d0
      Hgqqg_ab = 0d0; Hgqqg_ba = 0d0; Hgqqg_sym = 0d0
      Hagag_ab = 0d0; Hagag_ba = 0d0; Hagag_sym = 0d0
      Hgaag_ab = 0d0; Hgaag_ba = 0d0; Hgaag_sym = 0d0
      Hqagg_ab = 0d0; Hqagg_ba = 0d0; Hqagg_sym = 0d0
      Haqgg_ab = 0d0; Haqgg_ba = 0d0; Haqgg_sym = 0d0

      Hgggg = 0d0; Hgggg_1256 = 0d0
      Hgggg_1265 = 0d0; Hgggg_1625 = 0d0
      Hggqa = 0d0; Hggqa_ab = 0d0; Hggqa_ba = 0d0; Hggqa_sym = 0d0

      Hqqqq = 0d0; Hqqqq_a = 0d0; Hqqqq_b = 0d0; Hqqqq_i = 0d0
      Hqaqa = 0d0; Hqaqa_a = 0d0; Hqaqa_b = 0d0; Hqaqa_i = 0d0
      Haqaq = 0d0; Haqaq_a = 0d0; Haqaq_b = 0d0; Haqaq_i = 0d0
      Hqaaq = 0d0; Hqaaq_a = 0d0; Hqaaq_b = 0d0; Hqaaq_i = 0d0

      Hqrqr = 0d0; Habab = 0d0; Haaaa = 0d0; Hqarb = 0d0
      Hqbqb = 0d0; Haqbr = 0d0; Hbqbq = 0d0

      !call pvextclearcache

C   Deal with Higgs decay
      s34 = massvec(p_mcfm(3,:)+p_mcfm(4,:))

      if (hdecaymode == 'tlta') then
          call htautaudecay(p_mcfm,3,4,hdecay)
      elseif (hdecaymode == 'bqba') then
          call hbbdecay(p_mcfm,3,4,hdecay)
      elseif (hdecaymode == 'gaga') then
          hdecay=msqgamgam(sqrt(s34))
      else
        write(6,*) 'Unimplemented process in hjetmass_r'
        call abort
      endif
      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)
C   Done with decay

      Asq=(as/(3d0*pi))**2/vevsq
      fac=gsq**2*Asq *hdecay

      p_qp(:,:) = 0d0
      p_qp(1,:) = p_mcfm(1,:)
      p_qp(2,:) = p_mcfm(2,:)
      p_qp(3,:) = p_mcfm(iglue1,:)
      p_qp(4,:) = p_mcfm(iglue2,:)

      p_dp(:,:) = 0d0
      p_dp(1,:) = p_mcfm(1,:)
      p_dp(2,:) = p_mcfm(2,:)
      p_dp(3,:) = p_mcfm(iglue1,:)
      p_dp(4,:) = p_mcfm(iglue2,:)

! ========== 
! === qqchan
! ========== 

      if (qqchan) then
! qarb
      ret = passed_scaletest(hjetmass_HqqbQQb_nonid_dp,
     &                       hjetmass_HqqbQQb_nonid,
     &              p_dp, p_qp, 4,3,1,2,
     &              Hqarb, dum1, dum2, dum3, "Hqarb")

      if (ret .eqv. .false.) then
        badpoint = .true.
        return
      endif

! qbqb
      ret = passed_scaletest(hjetmass_HqqbQQb_nonid_dp,
     &                       hjetmass_HqqbQQb_nonid,
     &              p_dp, p_qp, 4,2,1,3,
     &              Hqbqb, dum1, dum2, dum3, "Hqbqb")

      if (ret .eqv. .false.) then
        badpoint = .true.
        return
      endif

! qrqr
      ret = passed_scaletest(hjetmass_HqqbQQb_nonid_dp,
     &                       hjetmass_HqqbQQb_nonid,
     &              p_dp, p_qp, 2,4,1,3,
     &              Hqrqr, dum1, dum2, dum3, "Hqrqr")

      if (ret .eqv. .false.) then
        badpoint = .true.
        return
      endif

! qqqq
      ret = passed_scaletest(hjetmass_HqqbQQb_ident_dp,
     &                       hjetmass_HqqbQQb_ident,
     &              p_dp, p_qp, 2,4,1,3,
     &              Hqqqq, Hqqqq_a, Hqqqq_b, Hqqqq_i, "Hqqqq")

      if (ret .eqv. .false.) then
        badpoint = .true.
        return
      endif

! qaqa
      ret = passed_scaletest(hjetmass_hqqbqqb_ident_dp,
     &                       hjetmass_hqqbqqb_ident,
     &              p_dp, p_qp, 4,2,1,3,
     &              Hqaqa, Hqaqa_a, Hqaqa_b, Hqaqa_i, "Hqaqa")

      if (ret .eqv. .false.) then
        badpoint = .true.
        return
      endif

      endif ! qqchan


! ========== 
! === gqchan
! ========== 

      if (gqchan) then
! qagg
      ret = passed_scaletest(hjetmass_Hqqbgg_dp, hjetmass_Hqqbgg,
     &              p_dp, p_qp, 1,2,3,4,
     &              Hqagg, Hqagg_ab, Hqagg_ba, Hqagg_sym, "Hqagg")

      if (ret .eqv. .false.) then
        badpoint = .true.
        return
      endif

! aqgg
      ret = passed_scaletest(hjetmass_Hqqbgg_dp, hjetmass_Hqqbgg,
     &              p_dp, p_qp, 2,1,3,4,
     &              Haqgg, Haqgg_ab, Haqgg_ba, Haqgg_sym, "Haqgg")

      if (ret .eqv. .false.) then
        badpoint = .true.
        return
      endif

! qgqg
      ret = passed_scaletest(hjetmass_Hqqbgg_dp, hjetmass_Hqqbgg,
     &              p_dp, p_qp, 1,3,2,4,
     &              Hqgqg, Hqgqg_ab, Hqgqg_ba, Hqgqg_sym, "Hqgqg")

      if (ret .eqv. .false.) then
        badpoint = .true.
        return
      endif

! agag
      ret = passed_scaletest(hjetmass_Hqqbgg_dp, hjetmass_Hqqbgg,
     &              p_dp, p_qp, 3,1,2,4,
     &              Hagag, Hagag_ab, Hagag_ba, Hagag_sym, "Hagag")

      if (ret .eqv. .false.) then
        badpoint = .true.
        return
      endif

! gqqg
      ret = passed_scaletest(hjetmass_Hqqbgg_dp, hjetmass_Hqqbgg,
     &              p_dp, p_qp, 2,3,1,4,
     &              Hgqqg, Hgqqg_ab, Hgqqg_ba, Hgqqg_sym, "Hgqqg")

      if (ret .eqv. .false.) then
        badpoint = .true.
        return
      endif

! gaag
      ret = passed_scaletest(hjetmass_Hqqbgg_dp, hjetmass_Hqqbgg,
     &              p_dp, p_qp, 3,2,1,4,
     &              Hgaag, Hgaag_ab, Hgaag_ba, Hgaag_sym, "Hgaag")

      if (ret .eqv. .false.) then
        badpoint = .true.
        return
      endif

      end if ! gqchan

! ========== 
! === ggchan
! ========== 

      if (ggchan) then
! ggqa
      ret = passed_scaletest(hjetmass_Hqqbgg_dp, hjetmass_Hqqbgg,
     &              p_dp, p_qp, 4,3,1,2,
     &              Hggqa, Hggqa_ab, Hggqa_ba, Hggqa_sym, "Hggqa")

      if (ret .eqv. .false.) then
        badpoint = .true.
        return
      endif

! gggg
      ! 1,2,3,4 passed as dummy arguments
      ret = passed_scaletest(hjetmass_Hgggg_dp, hjetmass_Hgggg,
     &          p_dp, p_qp, 1,2,3,4,
     &          Hgggg, Hgggg_1256, Hgggg_1265, Hgggg_1625, "Hgggg")

      if (ret .eqv. .false.) then
        badpoint = .true.
        return
      endif

      end if ! ggchan

! set other identical channels 
      Haqbr=Hqarb
      
      Haqaq=Hqaqa
      Haqaq_a=Hqaqa_a
      Haqaq_b=Hqaqa_b
      Haqaq_i=Hqaqa_i
      Hbqbq=Hqbqb

      Habab=Hqrqr
      Haaaa=Hqqqq

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      msq_struc(:,j,k)=0d0

      if ((j.gt.0).and.(k.gt.0)) then 
        if (j.eq.k) then
          msq(j,k)=0.5d0*aveqq*fac*Hqqqq
          msq_struc(iqq_a,j,k)=0.5d0*aveqq*fac*Hqqqq_a
          msq_struc(iqq_b,j,k)=0.5d0*aveqq*fac*Hqqqq_b
          msq_struc(iqq_i,j,k)=0.5d0*aveqq*fac*Hqqqq_i
        else
          msq(j,k)=aveqq*fac*Hqrqr
          msq_struc(iqq_a,j,k)=msq(j,k)
          msq_struc(iqq_b,j,k)=0d0
          msq_struc(iqq_i,j,k)=0d0
        endif
      endif
      
      if ((j.lt.0).and.(k.lt.0)) then 
        if (j.eq.k) then
          msq(j,k)=0.5d0*aveqq*fac*Haaaa
        else
          msq(j,k)=aveqq*fac*Habab
          msq_struc(iqq_a,j,k)=msq(j,k)
          msq_struc(iqq_b,j,k)=0d0
          msq_struc(iqq_i,j,k)=0d0
        endif
      endif

      if ((j.gt.0).and.(k.lt.0)) then
        if (j.eq.-k) then
          msq(j,k)=aveqq*fac*(0.5d0*Hqagg+Hqaqa+dfloat(nflav-1)*Hqarb)
          msq_struc(iqr,j,k)=aveqq*fac*dfloat(nflav-1)*Hqarb
          msq_struc(iqq_a,j,k)=aveqq*fac*Hqaqa_a
          msq_struc(iqq_b,j,k)=aveqq*fac*Hqaqa_b
          msq_struc(iqq_i,j,k)=aveqq*fac*Hqaqa_i
          msq_struc(igg_ab,j,k)=aveqq*fac*0.5d0*Hqagg_ab
          msq_struc(igg_ba,j,k)=aveqq*fac*0.5d0*Hqagg_ba
          msq_struc(igg_sym,j,k)=aveqq*fac*0.5d0*Hqagg_sym
        else
          msq(j,k)=aveqq*fac*Hqbqb
          msq_struc(iqq_a,j,k)=msq(j,k)
          msq_struc(iqq_b,j,k)=0d0
          msq_struc(iqq_i,j,k)=0d0
        endif
      endif

      if ((j.lt.0).and.(k.gt.0)) then
        if (j.eq.-k) then
          msq(j,k)=aveqq*fac*(0.5d0*Haqgg+Haqaq+dfloat(nflav-1)*Haqbr)
          msq_struc(iqr,j,k)=aveqq*fac*dfloat(nflav-1)*Haqbr
          msq_struc(iqq_a,j,k)=aveqq*fac*Haqaq_a
          msq_struc(iqq_b,j,k)=aveqq*fac*Haqaq_b
          msq_struc(iqq_i,j,k)=aveqq*fac*Haqaq_i
          msq_struc(igg_ab,j,k)=aveqq*fac*0.5d0*Haqgg_ab
          msq_struc(igg_ba,j,k)=aveqq*fac*0.5d0*Haqgg_ba
          msq_struc(igg_sym,j,k)=aveqq*fac*0.5d0*Haqgg_sym
        else
          msq(j,k)=aveqq*fac*Hbqbq
          msq_struc(iqq_a,j,k)=msq(j,k)
          msq_struc(iqq_b,j,k)=0d0
          msq_struc(iqq_i,j,k)=0d0
        endif
      endif

      if ((j.gt.0).and.(k.eq.0)) then
        msq(j,0)=aveqg*fac*Hqgqg
        msq_struc(igg_ab,j,0)=aveqg*fac*Hqgqg_ab
        msq_struc(igg_ba,j,0)=aveqg*fac*Hqgqg_ba
        msq_struc(igg_sym,j,0)=aveqg*fac*Hqgqg_sym
      endif
      
      if ((j.lt.0).and.(k.eq.0)) then
        msq(j,0)=aveqg*fac*Hagag
        msq_struc(igg_ab,j,0)=aveqg*fac*Hagag_ab
        msq_struc(igg_ba,j,0)=aveqg*fac*Hagag_ba
        msq_struc(igg_sym,j,0)=aveqg*fac*Hagag_sym
      endif

      if ((j.eq.0).and.(k.gt.0)) then
        msq(0,k)=aveqg*fac*Hgqqg
        msq_struc(igg_ab,0,k)=aveqg*fac*Hgqqg_ab
        msq_struc(igg_ba,0,k)=aveqg*fac*Hgqqg_ba
        msq_struc(igg_sym,0,k)=aveqg*fac*Hgqqg_sym
      endif

      if ((j.eq.0).and.(k.lt.0)) then
        msq(0,k)=aveqg*fac*Hgaag
        msq_struc(igg_ab,0,k)=aveqg*fac*Hgaag_ab
        msq_struc(igg_ba,0,k)=aveqg*fac*Hgaag_ba
        msq_struc(igg_sym,0,k)=aveqg*fac*Hgaag_sym
      endif

      if ((j.eq.0).and.(k.eq.0)) then
        msq(0,0)=avegg*fac*(0.5d0*Hgggg+dfloat(nflav)*Hggqa)
        msq_struc(igg_ab,0,0)=avegg*fac*dfloat(nflav)*Hggqa_ab
        msq_struc(igg_ba,0,0)=avegg*fac*dfloat(nflav)*Hggqa_ba
        msq_struc(igg_sym,0,0)=avegg*fac*dfloat(nflav)*Hggqa_sym
        msq_struc(igggg_a,0,0)=avegg*fac*0.5d0*Hgggg_1256
        msq_struc(igggg_b,0,0)=avegg*fac*0.5d0*Hgggg_1625
        msq_struc(igggg_c,0,0)=avegg*fac*0.5d0*Hgggg_1265
      endif
      
      enddo
      enddo

c--- subtraction matrix elements use qa->aq; calculate this and
c--- artificially store it in msq_struc(iqr,0,0), which is not
c--- used for anything else

      if (qqchan) then

      ret = passed_scaletest(hjetmass_HqqbQQb_ident_dp,
     &                       hjetmass_HqqbQQb_ident,
     &          p_dp, p_qp, 3,2,1,4,
     &          Hqaaq, Hqaaq_a, Hqaaq_b, Hqaaq_i, "Hqaaq")

      if (ret .eqv. .false.) then
        badpoint = .true.
        return
      endif

      end if ! qqchan

      msq_struc(iqr,0,0)=aveqq*fac*Hqaaq

      end subroutine

      logical function passed_scaletest(sub_msq_dp, sub_msq_qp,
     &          p_dp, p_qp, i1,i2,i3,i4,a,b,c,d,str)
          implicit none
          include 'types.f'
          include 'mxpart.f'
          include 'masses.f'
          include 'scale.f'
          include 'debug.f'
          include 'verbose.f'

          real(dp), intent(in) :: p_dp(mxpart,4)
          real(qp), intent(in) :: p_qp(mxpart,4)
          integer, intent(in) :: i1,i2,i3,i4
          real(dp), intent(out) :: a,b,c,d
          character(*), intent(in) :: str

          real(dp) :: a_dp,b_dp,c_dp,d_dp
          real(dp) :: a_dp_scale,b_dp_scale,c_dp_scale,d_dp_scale

          real(qp) :: a_qp,b_qp,c_qp,d_qp
          real(qp) :: a_qp_scale,b_qp_scale,c_qp_scale,d_qp_scale

          real(dp) :: p_dp_scale(mxpart,4)
          real(qp) :: p_qp_scale(mxpart,4)

          ! factor to use for scaling test
          real(dp), parameter :: scaletest = 1.2_dp
          ! If the ratio of the scaled results minus one differs by more than
          ! scaletest_cutoff, we either continue with quad precision
          ! if enable_qp is .true. or discard the point.
          ! If the quad precision scale test fails the point is
          ! discarded.
          real(dp), parameter :: scaletest_cutoff = 1.e-6_dp
          real(dp) :: scaletest_dynamic

          real(dp) :: ratio_dp
          real(qp) :: ratio_qp

          logical, parameter :: enable_qp = .true.
          ! Set debugprint to 1 for showing discarded points.
          ! Set it to 2 to also show qp fixed points.
          ! Currently it's set to 2 if debug from the input.DAT
          ! is set to .true. and 1 if verbose is .true.
          integer :: debugprint

          interface
            subroutine sub_msq_dp(i1,i2,i3,i4,p_dp,a,b,c,d)
                implicit none
                include 'types.f'
                include 'mxpart.f'
                integer, intent(in) :: i1,i2,i3,i4
                real(dp), intent(in) :: p_dp(mxpart,4)
                real(dp), intent(out) :: a,b,c,d
            end subroutine

            subroutine sub_msq_qp(i1,i2,i3,i4,p_qp,a,b,c,d)
                implicit none
                include 'types.f'
                include 'mxpart.f'
                integer, intent(in) :: i1,i2,i3,i4
                real(qp), intent(in) :: p_qp(mxpart,4)
                real(qp), intent(out) :: a,b,c,d
            end subroutine
          end interface

        debugprint = 0
        if (verbose) debugprint = 1
        if (debug) debugprint = 2

        p_dp_scale = p_dp*scaletest
        p_qp_scale = p_qp*scaletest

        call sub_msq_dp(i1,i2,i3,i4, p_dp,a_dp,b_dp,c_dp,d_dp)

        mt = mt*scaletest
        hmass = hmass*scaletest
        musq = musq*scaletest**2

        call sub_msq_dp(i1,i2,i3,i4,p_dp_scale,
     &          a_dp_scale,b_dp_scale,c_dp_scale,d_dp_scale)

        mt = mt/scaletest
        hmass = hmass/scaletest
        musq = musq/scaletest**2

        ratio_dp = a_dp_scale/a_dp

        scaletest_dynamic = min(1.e-6_dp, 10._dp**(-abs(log(a_dp)/log(10._dp))))

        if (abs(ratio_dp - 1._dp) > scaletest_dynamic
     &              .or. isnan(ratio_dp)) then

            if (enable_qp) then

              call sub_msq_qp(i1,i2,i3,i4,p_qp,a_qp,b_qp,c_qp,d_qp)

              mt = mt*scaletest
              hmass = hmass*scaletest
              musq = musq*scaletest**2

              call sub_msq_qp(i1,i2,i3,i4,p_qp_scale,
     &              a_qp_scale,b_qp_scale,c_qp_scale,d_qp_scale)

              mt = mt/scaletest
              hmass = hmass/scaletest
              musq = musq/scaletest**2

              ratio_qp = a_qp_scale/a_qp

              if(abs(ratio_qp - 1._dp) > scaletest_dynamic .or.
     &                                              isnan(ratio_qp)) then
                if (debugprint > 0) write (*,*) str,
     &              " discarding broken point, scaling qp,dp:", real(ratio_qp,dp), ratio_dp, scaletest_dynamic
                passed_scaletest = .false.
                return
              else
                if (debugprint > 1) write (*,*) str,
     &              " fixed broken point, scaling qp,dp:", real(ratio_qp,dp), ratio_dp
                a = real(a_qp,dp)
                b = real(b_qp,dp)
                c = real(c_qp,dp)
                d = real(d_qp,dp)
                passed_scaletest = .true.
                return
              endif

            else
                if (debugprint > 0) write (*,*) str,
     &              " discarding broken point, scaling dp:", ratio_dp, scaletest_dynamic
                passed_scaletest = .false.
                return
            endif
        else
            a = a_dp
            b = b_dp
            c = c_dp
            d = d_dp
            passed_scaletest = .true.
            return
        endif
      end function
      
