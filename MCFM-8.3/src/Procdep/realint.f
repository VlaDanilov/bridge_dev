      function realint(vector,wgt)
          use ieee_arithmetic
          use VVconfig_m
          use singletop2_m, only : singletop2_gs, singletop2_tree, singletop2_gs_light, singletop2_gs_heavy
          use singletop2_realamps_m, only : singletop2_real, singletop2_real_light, singletop2_real_heavy
          use singletop2_scale_m
          use bbfrac_m, only : bbfrac
      implicit none
      include 'types.f'
      real(dp):: realint
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'debug.f'
      include 'noglue.f'
      include 'nflav.f'
      include 'vegas_common.f'
      include 'ptilde.f'
      include 'phasemin.f'
      include 'new_pspace.f'
      include 'npart.f'
      include 'scale.f'
      include 'facscale.f'
      include 'dynamicscale.f'
      include 'realwt.f'
      include 'efficiency.f'
      include 'maxwt.f'
      include 'kprocess.f'
      include 'PDFerrors.f'
      include 'masses.f'
      include 'wts_bypart.f'
      include 'dipolescale.f'
      include 'stopscales.f'
      include 'flags.f'
      include 'frag.f'
      include 'ipsgen.f'
      include 'decay1q2a.f'
      include 'outputoptions.f'
      include 'breit.f'
      include 'dm_params.f' 
      include 'outputflags.f'
      include 'runstring.f'
      include 'bypart.f'
      include 'energy.f'
      include 'incldip.f'
      include 'nproc.f'
      include 'initialscales.f'
      include 'badpoint.f'
      include 'nqcdjets.f'
      include 'taucut.f'
      include 'qcdcouple.f'
      include 'scalevar.f'
      include 'couple.f'
      include 'nlooprun.f'
      include 'cutoff.f'
      include 'ewcorr.f'

c--- APPLgrid - enable grids
      include 'APPLinclude.f'
c      include 'qcdcouple.f'
c      include 'nlooprun.f'
c      include 'couple.f' 
      double precision psCR, currason2pi
c--- APPLgrid - end

c---- SSbegin
      include 'reweight.f'
c---- SSend
cz
cz //
      integer:: ih1,ih2,j,k,nd,nmax,nmin,nvec,ii,t
      integer itrial
      real(dp):: alphas,msqtrial,
     & fx1up(-nf:nf),fx2up(-nf:nf),fx1dn(-nf:nf),fx2dn(-nf:nf),
     & dipfx1up(0:maxd,-nf:nf),dipfx2up(0:maxd,-nf:nf),
     & dipfx1dn(0:maxd,-nf:nf),dipfx2dn(0:maxd,-nf:nf),
     & xmsqvar(2,0:maxd),asorig,scaleup,scaledn
      real(dp):: vector(mxdim),W,val,val2,valsum,xint,ptmp
      real(dp):: fx1(-nf:nf),fx2(-nf:nf),
     & dipfx1(0:maxd,-nf:nf),dipfx2(0:maxd,-nf:nf),
     & fx1_H(-nf:nf),fx2_H(-nf:nf),fx1_L(-nf:nf),fx2_L(-nf:nf)
      real(dp):: p(mxpart,4),pjet(mxpart,4),p1ext(4),p2ext(4)
      real(dp):: pswt
      real(dp):: s(mxpart,mxpart),wgt,msq(-nf:nf,-nf:nf)
      ! to distinguish for singletop2 between light and heavy line emissions
      real(dp) :: msq_light(-nf:nf,-nf:nf), msq_heavy(-nf:nf,-nf:nf)
      real(dp) :: msqc_light(maxd,-nf:nf,-nf:nf), msqc_heavy(maxd,-nf:nf,-nf:nf)

      real(dp):: msqa(-nf:nf,-nf:nf)
      real(dp):: msqc(maxd,-nf:nf,-nf:nf),xmsq(0:maxd),msqc_new(maxd,-nf:nf,-nf:nf),bit1,bit2
      real(dp):: flux,BrnRat,xreal,xreal2
      real(dp):: xx1,xx2,q(mxpart,4)
      real(dp):: m3,m4,m5,R,Rbbmin,dot
      real(dp):: xmsq_bypart(0:maxd,-1:1,-1:1),xmsqjk,
     & plo(mxpart,4),pswtdip,psave(mxpart,4)
      integer:: sgnj,sgnk
      common/xreal/xreal,xreal2
      common/Rbbmin/Rbbmin
      logical:: bin,failed
      logical:: includedipole,includereal
      real(dp):: QandGint
      external qqb_w2jet_g,qqb_w2jet_gs,qqb_z2jet_g,qqb_z2jet_gs,
     & qqb_w2jet,qqb_w1jet_gs,qqb_z2jet,qqb_z1jet_gs,qqb_Hg_g,qqb_Hg_gs,
     & qqb_hww_g,qqb_hww_gs,qqb_zbb_g,qqb_zbb_gs,
     & qqb_wh_ww,qqb_wh_ww_gs,
     & qqb_wh_zz,qqb_wh_zz_gs,
     & qqb_wbb_g,qqb_wbb_gs,
     & qqb_dirgam_g,qqb_dirgam_gs,qqb_hflgam_g,qqb_hflgam_gs,
     & qqb_trigam_g,qqb_trigam_gs,qqb_gmgmjt_g,qqb_gmgmjt_gs,
     & qqb_w_g,qqb_w_gs,qqb_z1jet,qqb_z_gs,qqb_ww_g,qqb_ww_gs,
     & qqb_wz_g,qqb_wz_gs,qqb_zz_g,qqb_zz_gs,qqb_wgam_g,qqb_wgam_gs,
     & qqb_QQb_g,qqb_QQb_gs,
     & VV_Hqq_g,VV_Hqq_gs,VV_HWW_g,VV_HWW_gs,
     & gg_Hg,gg_H_gs,
     & gg_HWWgg,gg_HWWg_gs,gg_HZZgg,gg_HZZg_gs,
     & gg_Hgg,gg_Hg_gs,
     & gQ_zQ_g,gQ_zQ_gs,qqb_tbb_g,qqb_tbb_gs,
     & qqb_w_tndk_g,qqb_w_tndk_gs,
     & qqb_w_twdk_g,qqb_w_twdk_gs,qqb_w_twdk_gdk,qqb_w_twdk_gsdk,
     & qqb_zbjet_g,qqb_zbjet_gs,qqb_w_cjet_g,qqb_w_cjet_gs,
     & qqb_wbfromc_g,qqb_wbfromc_gs,
     & gg_hggg,gg_hgg_gs,qq_tchan_htqg,qq_tchan_htq_gs,
     & qg_tbq_g,qg_tbq_gs,qq_tbg_g,qq_tbg_gs,epem3j_g,epem3j_gs,
     & qq_tchan_ztqg_dk,qq_tchan_ztq_dk_gs,
     & qq_tchan_ztqg,qq_tchan_ztq_gs,
     & qqb_QQbdk_g,qqb_QQbdk_gs,qqb_gamgam_g,qqb_gamgam_gs,
     & qqb_zaj_gs,qqb_zaj_g,
     & qqb_gmgmjt_g_mad,
     & qqb_zaa_g,qqb_zaa_gs,
     & qqb_waa_g,qqb_waa_gs,
     & qqb_WH1jet_g,qqb_WH1jet_gs,
     & qqb_ZH1jet_g,qqb_ZH1jet_gs,
     & qqb_waa_g_mad,qqb_tottth_g,qqb_tottth_g_mad
      common/density/ih1,ih2
      common/bin/bin
      common/Pext/p1ext,p2ext
!$omp threadprivate(/Pext/)
      common/nmax/nmax
      common/BrnRat/BrnRat
      common/nmin/nmin

      real(dp):: realeventp(mxpart,4)
      common/realeventp/realeventp

      external :: qqb_zgam_new, qqb_zgam


      if (1 .eq. 2) then
      vector(           1 )=  1.13928677147470225e-003_dp
      vector(           2 )=  1.13877556155377317e-003_dp
      vector(           3 )=  7.53534215141593577e-002_dp
      vector(           4 )=  1.86355430277158583e-002_dp
      vector(           5 )=  0.55134166051107103_dp
      vector(           6 )=  0.75627564371201084_dp
      vector(           7 )=  0.14147198578671921_dp
      vector(           8 )=  0.64150193209252071_dp
      vector(           9 )=  0.66515300845468606_dp
      vector(          10 )=  0.47746039559714804_dp
      vector(          11 )=  0.57882453871580786_dp
      vector(          12 )=  0.18107661738880171_dp
      vector(          13 )=  0.46006342202118478_dp
      vector(          14 )=  0.52616037406314187_dp
      vector(          15 )=  0.50150127313169712_dp
      endif

      QandGflag=.false.
      p(:,:)=0._dp
!$omp atomic
      ntotshot=ntotshot+1
      pswt=0._dp
      realint=0._dp      
c--- ensure isolation code does not think this is fragmentation piece
      z_frag=0._dp

      W=sqrts**2
      
c---- set default reweight to 1 (hence no reweighting)
      reweight = 1.0_dp

c      if (first) then
c         write(6,*)
c         write(6,*) 'nmin=',nmin,',nmax=',nmax
c         write(6,*)
c         first=.false.
c      endif

      p(:,:)=0._dp
      pjet(:,:)=0._dp

c-- note: new_pspace now signifies multi-channel integration
      if (new_pspace) then
        call gen_lops(vector,plo,pswt,*999)
!        call writeout(plo)
        npart=npart+1
        call multichan(vector(ndim-2),vector(ndim-1),vector(ndim),
     &                 vector(ndim+2),plo,p,pswtdip,*999)

        pswt=pswt*pswtdip
      else
        call gen_realps(vector,p,pswt,*999)
      endif
      psave(:,:)=p(:,:)

      if (1 .eq. 2) then
      p(  1,  1) =    0.00000000000000_dp
      p(  1,  2) =    0.00000000000000_dp
      p(  1,  3) =   -474.572402233294_dp
      p(  1,  4) =   -474.572402233294_dp
      p(  2,  1) =    0.00000000000000_dp
      p(  2,  2) =    0.00000000000000_dp
      p(  2,  3) =    16631.8709665902_dp
      p(  2,  4) =   -16631.8709665902_dp
      p(  3,  1) =    203.771496931320_dp
      p(  3,  2) =    316.299769747859_dp
      p(  3,  3) =   -541.279649719025_dp
      p(  3,  4) =    659.207997701128_dp
      p(  4,  1) =   -207.497710831028_dp
      p(  4,  2) =   -299.958895141825_dp
      p(  4,  3) =   -1261.60604651191_dp
      p(  4,  4) =    1313.27210169836_dp
      p(  5,  1) =  -0.389411990227555e-01_dp
      p(  5,  2) =   0.737966704545456e-02_dp
      p(  5,  3) =   -1578.53713848121_dp
      p(  5,  4) =    1578.53713903073_dp
      p(  6,  1) =   0.479675234604849e-02_dp
      p(  6,  2) =   0.237889496590332e-01_dp
      p(  6,  3) =   -13165.4699384615_dp
      p(  6,  4) =    13165.4699384465_dp
      p(  7,  1) =    3.76035834638489_dp
      p(  7,  2) =   -16.3720432227384_dp
      p(  7,  3) =    389.594208816684_dp
      p(  7,  4) =    389.956191946837_dp

      p(  1,  1) =    0.00000000000000_dp
      p(  1,  2) =    0.00000000000000_dp
      p(  1,  3) =   -3328.06591634386_dp
      p(  1,  4) =   -3328.06591634386_dp
      p(  2,  1) =    0.00000000000000_dp
      p(  2,  2) =    0.00000000000000_dp
      p(  2,  3) =    2849.35582416532_dp
      p(  2,  4) =   -2849.35582416532_dp
      p(  3,  1) =   -227.427440223096_dp
      p(  3,  2) =    187.136720354957_dp
      p(  3,  3) =    49.3278427388677_dp
      p(  3,  4) =    298.629848591713_dp
      p(  4,  1) =   -52.0615550785897_dp
      p(  4,  2) =   -375.182416823698_dp
      p(  4,  3) =   -665.903451353387_dp
      p(  4,  4) =    766.095826685034_dp
      p(  5,  1) =    7.24633951429031_dp
      p(  5,  2) =   -30.7841063522462_dp
      p(  5,  3) =    58.4179192751088_dp
      p(  5,  4) =    66.4290895067507_dp
      p(  6,  1) =    279.488971781183_dp
      p(  6,  2) =    188.045675240935_dp
      p(  6,  3) =   -1990.46526590120_dp
      p(  6,  4) =    2018.76874259616_dp
      p(  7,  1) =   -7.24631599378821_dp
      p(  7,  2) =    30.7841275800517_dp
      p(  7,  3) =    3027.33304741915_dp
      p(  7,  4) =    3027.49823312953_dp
      endif

c--- sanity check for NaN
      do j=1,npart+2
      if (p(j,4) /= p(j,4)) then
!        call writeout(p)
!        do k=1,ndim+2
!        write(6,*) '     vector(',k,')=',vector(k),'._dp'
!        enddo
!        write(6,*) 'bad PS point above'
!        write(6,*)
!          write(6,*) 'discarding NaN in phase space point'
        goto 999
      endif
      enddo

      nvec=npart+2
      call dotem(nvec,p,s)
      
c----calculate the x's for the incoming partons from generated momenta
      xx1=two*(p(1,4)*p2ext(4)-p(1,3)*p2ext(3))/W
      xx2=two*(p(2,4)*p1ext(4)-p(2,3)*p1ext(3))/W

c      if (debug) write(*,*) 'Reconstructed x1,x2 ',xx1,xx2
      
      if ((xx1 >  1._dp) .or. (xx2 >  1._dp)
     &.or.(xx1 < xmin) .or. (xx2 < xmin)) then
         goto 999
      endif

c--- (moved to includedipole) impose cuts on final state
c      call masscuts(p,*999)

      if (usescet) then
c----reject event if any tau is too small
        call smalltau(p,npart,*999)
      else
c----reject event if any s(i,j) is too small
        call smalls(s,npart,*999)
      endif
     
c--- extra cut to divide WQj/ZQj regions
      if ( (nproc == 312) .or. (nproc == 317)
     & .or.(nproc == 322) .or. (nproc == 327)
     & .or.(nproc == 342) .or. (nproc == 352)) then
        if (R(p,5,6) < Rbbmin) goto 999
      endif

c--- see whether this point will pass cuts - if it will not, do not
c--- bother calculating the matrix elements for it, instead set to zero
      includereal=includedipole(0,p)
      incldip(0)=includereal

c     it is always a good idea to initialize all elements as zero
      msq(:,:) = 0
      msqc(:,:,:) = 0
 
      if (includereal .eqv. .false.) then
        do j=-nf,nf
        do k=-nf,nf
          msqLH(j,k)=0._dp   ! for stop+b process
          msqHL(j,k)=0._dp   ! for stop+b process
        enddo
        enddo
      endif
      
      do j=1,mxpart
      do k=1,4
        realeventp(j,k)=p(j,k)
      enddo
      enddo
      
c--- test to see whether we need Gflag and Qflag together
      if ( ((kcase==kW_2jet) .or. (kcase==kZ_2jet))
     &.and. (Qflag) .and. (Gflag) ) then
        QandGflag=.true.
        QandGint=0._dp
c--- first pass: Gflag
        Gflag=.true.
        Qflag=.false.
      endif
      
c--- restart from here when calculating with Qflag and Gflag
c--- (W+2 jet and Z+2 jet processes only)
   44 continue   
      
      if (dynamicscale) then
        call scaleset(initscale,initfacscale,p)
        dipscale(0)=facscale
      endif

      if (kcase==ktopanom .and. use_DDIS) then
        ! we store an array for dipole modified factorization scales, dipole pdfs and pdfs here, which need to be reset
        call singletop2_scale_reset()
      endif
      
      if ((doscalevar) .and. (foundpow .eqv. .false.)) then
        itrial=1
      endif
   66 continue

c---- generate collinear points that satisfy the jet cuts (for checking)
c      call singgen(p,s,*998)
            
c--- Calculate the required matrix elements
      if     (kcase==kW_only) then
c        call singcheck(qqb_w_g,qqb_w_gs,p)         ! Checked 11/30/01
        if (includereal) call qqb_w_g(p,msq)
        call qqb_w_gs(p,msqc)
      elseif (kcase==kW_1jet) then
c        call singcheck(qqb_w2jet,qqb_w1jet_gs,p)   ! Checked 11/16/01
        if (includereal) call qqb_w2jet(p,msq)
        call qqb_w1jet_gs(p,msqc)  
      elseif (kcase==kWgamma) then
!         if(includereal) call singcheck(qqb_wgam_g,qqb_wgam_gs,p) ! Checked 08/27/02
        if (includereal) call qqb_wgam_g(p,msq)
        call qqb_wgam_gs(p,msqc)  
      elseif (kcase==kWbfrmc) then
        if (includereal) call qqb_wbfromc_g(p,msq)
        call qqb_wbfromc_gs(p,msqc)  
      elseif (kcase==kW_cjet) then
c        call singcheck(qqb_w_cjet_g,qqb_w_cjet_gs,p) ! Checked 15/05/07
        if (includereal) call qqb_w_cjet_g(p,msq)
        call qqb_w_cjet_gs(p,msqc)
      elseif (kcase==kZgamma) then
        if (decayChannel() == decayQuarks) then
          !if(includereal) call singcheck(qqb_zaj_vdecay,qqb_zgam_gs,p)
          if (includereal) call qqb_zaj_vdecay(p,msq)
          call qqb_zgam_gs(p,msqc)
        else
          !if(includereal) call singcheck(qqb_zaj,qqb_zgam_gs,p)
          call set_anomcoup(p)
          if (includereal) call qqb_zaj(p,msq)
          call qqb_zgam_gs(p,msqc)
        endif
      elseif (kcase==kZ_2gam) then
c         if(includereal) call singcheck(qqb_zaa_g,qqb_zaa_gs,p)
        if (includereal) call qqb_zaa_g(p,msq)
        call qqb_zaa_gs(p,msqc)
      elseif (kcase==kW_2gam) then
c        if(includereal) call singcheck(qqb_waa_g,qqb_waa_gs,p)
c        call compare_madgraph(p,qqb_waa_g,qqb_waa_g_mad)
         stop
c        if (includereal) call qqb_waa_g(p,msq)      
c        call qqb_waa_gs(p,msqc)
      elseif (kcase==kZgajet) then
          if (decayChannel() == decayQuarks) then
c             if(includereal) call singcheck(qqb_zaj_g_vdecay,qqb_zaj_gs_vdecay,p)
              if (includereal) call qqb_zaj_g_vdecay(p,msq)
              call qqb_zaj_gs_vdecay(p,msqc)
          else
              !if(includereal) call singcheck(qqb_zaj_g,qqb_zaj_gs,p)
              call set_anomcoup(p)
              if (includereal) call qqb_zaj_g(p,msq)
              call qqb_zaj_gs(p,msqc)
          endif
      elseif (kcase==kWbbmas) then
c        call singcheck(qqb_wbbm_g,qqb_wbbm_gs,p)     ! Checked 10/21/10
      if (includereal) call qqb_wbbm_g(p,msq)      
        call qqb_wbbm_gs(p,msqc)      
      elseif (kcase==kWttmas) then
c        call singcheck(qqb_wbbm_g,qqb_wbbm_gs,p)     ! Checked 10/21/10
      if (includereal) call qqb_wbbm_g(p,msq)      
        call qqb_wbbm_gs(p,msqc)      
      elseif (kcase==kWbbbar) then
c        call singcheck(qqb_wbb_g,qqb_wbb_gs,p)     ! Checked 11/30/01
        if (includereal) call qqb_wbb_g(p,msq)      
        call qqb_wbb_gs(p,msqc)      
      elseif (kcase==kW_2jet) then  
c        if (includereal) call singcheck(qqb_w2jet_g,qqb_w2jet_gs_new,p) ! Re-checked June 09
        if (includereal)  call qqb_w2jet_g(p,msq)
        call qqb_w2jet_gs_new(p,msqc)
      elseif (kcase==kZ_only) then
c        call singcheck(qqb_z1jet,qqb_z_gs,p)         ! Checked 11/30/01
        if (includereal) call qqb_z1jet(p,msq)      
        call qqb_z_gs(p,msqc)     
      elseif (kcase==kZ_1jet) then
c        call singcheck(qqb_z2jet,qqb_z1jet_gs,p)   ! Checked 11/16/01
        if (includereal) call qqb_z2jet(p,msq)      
        call qqb_z1jet_gs(p,msqc)  
      elseif (kcase==kZ_2jet) then
c        if (includereal) call singcheck(qqb_z2jet_g,qqb_z2jet_gs_new,p) ! Re-checked June 09
        if (includereal) call qqb_z2jet_g(p,msq)
        call qqb_z2jet_gs_new(p,msqc)
c        call qqb_z2jet_gs(p,msqc) 
c        write(6,*) 'Gflag,Qflag',Gflag,Qflag
c        do j=-nf,nf
c        do k=-nf,nf
c        bit1=0._dp
c        bit2=0._dp
c        bit1=sum(msqc(1:ndmax,j,k))
c        bit2=sum(msqc_new(1:ndmax,j,k))
c        if ((abs(bit1) .gt. 1.e-100_dp) .or. (abs(bit2) .gt. 1.e-100_dp)) then
c          if (abs(bit1/bit2-1._dp) .gt. 1.e-12_dp) then
c            write(6,*) 'j,k,bit1,bit1/bit2',j,k,bit1,bit2,bit1/bit2-1._dp
c          endif
c        endif
c        enddo
c        enddo
c        pause
      elseif (kcase==kZbbbar) then
c        call singcheck(qqb_zbb_g,qqb_zbb_gs,p)     ! Checked 11/30/01
        if (includereal) call qqb_zbb_g(p,msq)
        call qqb_zbb_gs(p,msqc) 
      elseif (kcase==kWWqqbr) then
c        call singcheck(qqb_ww_g,qqb_ww_gs,p)       ! Checked 11/30/01
        if (includereal) call qqb_ww_g(p,msq)
        call qqb_ww_gs(p,msqc)      
      elseif (kcase==kWWqqdk) then
        if (includereal) call dkqqb_ww_g(p,msq)
        call dkqqb_ww_gs(p,msqc)      
      elseif (kcase==kWZbbar) then
c        call singcheck(qqb_wz_g,qqb_wz_gs,p)       ! Checked 12/05/01
        if (includereal) call qqb_wz_g(p,msq)
        call qqb_wz_gs(p,msqc)      
      elseif (kcase==kZZlept) then
c        call singcheck(qqb_zz_g,qqb_zz_gs,p)       ! Checked 12/05/01
        if (includereal) call qqb_zz_g(p,msq)
        call qqb_zz_gs(p,msqc)      
      elseif (kcase==kWHbbar) then
c        call singcheck(qqb_wh_g,qqb_wh_gs,p)
        if (includereal) call qqb_wh_g(p,msq)
        call qqb_wh_gs(p,msqc)
      elseif (kcase==kWHbbdk) then
        if (includereal) call dkqqb_wh_g(p,msq)
        call dkqqb_wh_gs(p,msqc)
      elseif (kcase==kWH1jet) then
c     call singcheck(qqb_WH1jet_g,qqb_WH1jet_gs,p)
         if (toponly) then
            msq=zip
            msqc=zip
         else
            if (includereal) call qqb_WH1jet_g(p,msq)
            call qqb_WH1jet_gs(p,msqc)
         endif
      elseif (kcase==kWH__WW) then
c        call singcheck(qqb_wh_ww_g,qqb_wh_ww_gs,p)
        if (includereal) call qqb_wh_ww_g(p,msq)
        call qqb_wh_ww_gs(p,msqc)
      elseif (kcase==kWH__ZZ) then
c        call singcheck(qqb_wh_zz_g,qqb_wh_zz_gs,p)
        if (includereal) call qqb_wh_zz_g(p,msq)
        call qqb_wh_zz_gs(p,msqc)
      elseif (kcase==kWHgaga) then
c        call singcheck(qqb_wh_gaga_g,qqb_wh_gaga_gs,p)
        if (includereal) call qqb_wh_gaga_g(p,msq)
        call qqb_wh_gaga_gs(p,msqc) 
      elseif (kcase==kZHbbar) then
c        call singcheck(qqb_zh_g,qqb_zh_gs,p)
        if (includereal) call qqb_zh_g(p,msq)
        call qqb_zh_gs(p,msqc)     
      elseif (kcase==kZHbbdk) then
         if (includereal) call dkqqb_zh_g(p,msq)
         call dkqqb_zh_gs(p,msqc)
         
      elseif (kcase==kZHgaga) then
c        call singcheck(qqb_zh_gaga_g,qqb_zh_gaga_gs,p)
        if (includereal) call qqb_zh_gaga_g(p,msq)
        call qqb_zh_gaga_gs(p,msqc)     
      elseif (kcase==kZH__WW) then
c        call singcheck(qqb_zh_ww_g,qqb_zh_ww_gs,p)
        if (includereal) call qqb_zh_ww_g(p,msq)
        call qqb_zh_ww_gs(p,msqc)     
      elseif (kcase==kZH__ZZ) then
c        call singcheck(qqb_zh_zz_g,qqb_zh_zz_gs,p)
        if (includereal) call qqb_zh_zz_g(p,msq)
        call qqb_zh_zz_gs(p,msqc)     
      elseif (kcase==kZH1jet) then
c        call singcheck(qqb_ZH1jet_g,qqb_ZH1jet_gs,p)
c        if (includereal) call qqb_ZHgg_mad(p,msqa)
c         do j=-nf,nf
c         do k=-nf,nf
c         write(6,*) 'j,k,msqa',j,k,msq(j,k),msqa(j,k)
c         enddo
c         enddo
c         pause
        if (toponly) then
          msq=zip
          msqc=zip
        else
          if (includereal) call qqb_ZH1jet_g(p,msq)
          call qqb_ZH1jet_gs(p,msqc)
        endif
      elseif (kcase==ktwo_ew) then
c        if (includereal) call singcheck(qqb_twojet_mix_g,qqb_twojet_mix_gs,p) 
        if (  (-s(1,3) < cutoff) .or. (-s(2,3) < cutoff)
     &   .or. (-s(1,4) < cutoff) .or. (-s(2,4) < cutoff)
     &   .or. (-s(1,5) < cutoff) .or. (-s(2,5) < cutoff)
     &   .or. ( s(3,4) < cutoff) .or. ( s(3,5) < cutoff)
     &   .or. ( s(4,5) < cutoff)) goto 999
        if (includereal) call qqb_twojet_mix_g(p,msq)
        call qqb_twojet_mix_gs(p,msqc)
        wt_noew=zip ! non-EW calculation has no real contribution
      elseif (kcase==kdirgam) then
!        if (includereal) call singcheck(qqb_dirgam_g,qqb_dirgam_gs,p) 
        if (includereal) call qqb_dirgam_g(p,msq)
        call qqb_dirgam_gs(p,msqc)
      elseif (kcase==khflgam) then
c        if (includereal) call singcheck(qqb_hflgam_g,qqb_hflgam_gs,p)
        if (includereal) call qqb_hflgam_g(p,msq)
        call qqb_hflgam_gs(p,msqc)
      elseif (kcase==kgamgam) then
c        if (includereal) call singcheck(qqb_gamgam_g,qqb_gamgam_gs,p)
        if (includereal) call qqb_gamgam_g(p,msq)
        call qqb_gamgam_gs(p,msqc)
      elseif (kcase==kgg2gam) then
c        if (includereal) call singcheck(gg_2gam_g,gg_2gam_gs,p)
        if (includereal) call gg_2gam_g(p,msq)
        call gg_2gam_gs(p,msqc)
      elseif (kcase==kgmgmjt) then
c        if (includereal) call singcheck(qqb_gmgmjt_g,qqb_gmgmjt_gs,p)
        if (includereal) call qqb_gmgmjt_g(p,msq)
        call qqb_gmgmjt_gs(p,msqc)
      elseif (kcase==ktrigam) then
c        if (includereal) call singcheck(qqb_trigam_g,qqb_trigam_gs,p)
        if (includereal) call qqb_trigam_g(p,msq)
        call qqb_trigam_gs(p,msqc)
      elseif (kcase==kfourga) then 
c        if (includereal) call singcheck(qqb_fourgam_g,qqb_fourgam_gs,p)       
        if (includereal) call qqb_fourgam_g(p,msq)
        call qqb_fourgam_gs(p,msqc)
       elseif (kcase==kggfus0) then
c         call singcheck(gg_hg,gg_h_gs,p)       ! Checked 28/02/03
         if (includereal) call gg_hg(p,msq)
         call gg_h_gs(p,msqc)
       elseif (kcase==kHigaga) then
c         call singcheck(gg_hgamgamg,gg_hgamgam_gs,p)
         if (includereal) call gg_hgamgamg(p,msq)
         call gg_hgamgam_gs(p,msqc)
       elseif (kcase==kHi_Zga) then
c         call singcheck(gg_hzgamg,gg_hzgam_gs,p)
         if (includereal) call gg_hzgamg(p,msq)
         call gg_hzgam_gs(p,msqc)
      elseif ((kcase==kHWW_4l) .or. (kcase==kHWW2lq)) then
c        call singcheck(qqb_hww_g,qqb_hww_gs,p)
        if (includereal) call qqb_hww_g(p,msq)      
        call qqb_hww_gs(p,msqc)      
      elseif (kcase==kHWWdkW) then
        if (includereal) call dkqqb_hww_g(p,msq)      
        call dkqqb_hww_gs(p,msqc)      
      elseif (kcase==kHWWdkW) then
        if (includereal) call dkqqb_hww_g(p,msq)      
        call dkqqb_hww_gs(p,msqc)      
      elseif (kcase==kHZZ_4l) then
c        call singcheck(qqb_hzz_g,qqb_hzz_gs,p)
        if (includereal) call qqb_hzz_g(p,msq)      
        call qqb_hzz_gs(p,msqc)  
      elseif (kcase==kH_1jet) then
c        call singcheck(qqb_Hg_g,qqb_Hg_gs,p)       ! Checked 19/02/02
        if (includereal) call qqb_Hg_g(p,msq)  
        call qqb_Hg_gs(p,msqc) 
      elseif ((kcase==ktt_bbl) .or. (kcase==ktt_bbh)) then
c         call singcheck(qqb_QQbdk_g,qqb_QQbdk_gs,p) ! Checked 15/8/08
        if (includereal) call qqb_QQbdk_g(p,msq)  
        call qqb_QQbdk_gs(p,msqc) 
      elseif ((kcase==ktt_ldk) .or. (kcase==ktt_hdk)) then
        if     (decay1q2a == 1) then
c------ radiation in the decay of the top quark
        if (includereal) call dk1qqb_QQb_g(p,msq)
          call dk1qqb_QQb_gs(p,msqc)
      elseif (decay1q2a == 2) then 
c------ radiation in the decay of the anti-top quark
        if (includereal) call dk2qqb_QQb_g(p,msq)
          call dk2qqb_QQb_gs(p,msqc)
      else
        write(6,*) 'Problem generating PS, decay1q2a=',decay1q2a
        stop
      endif
      elseif (kcase==ktt_udk) then
        if     (decay1q2a == 1) then
c------ radiation in the decay of the top quark
        if (includereal) call dk1uqqb_QQb_g(p,msq)
          call dk1uqqb_QQb_gs(p,msqc)
      elseif (decay1q2a == 2) then 
c------ radiation in the decay of the anti-top quark
        if (includereal) call dk2uqqb_QQb_g(p,msq)
          call dk2uqqb_QQb_gs(p,msqc)
      else
        write(6,*) 'Problem generating PS, decay1q2a=',decay1q2a
        stop
      endif
      elseif (kcase==ktthWdk) then
        if     (decay1q2a == 1) then
c------ radiation in the decay of the top quark
        if (includereal) call dkW1qqb_QQb_g(p,msq)
          call dkW1qqb_QQb_gs(p,msqc)
      elseif (decay1q2a == 2) then 
c------ radiation in the decay of the anti-top quark
        if (includereal) call dkW2qqb_QQb_g(p,msq)
          call dkW2qqb_QQb_gs(p,msqc)
      else
        write(6,*) 'Problem generating PS, decay1q2a=',decay1q2a
        stop
      endif
      elseif (kcase==ktt_bbu) then
c         call singcheck(qqb_QQbdku_g,qqb_QQbdku_gs,p) !
c         pause
        if (includereal) call qqb_QQbdku_g(p,msq)  
        call qqb_QQbdku_gs(p,msqc) 
      elseif (kcase==ktt_mix) then
         if ((s(3,5) < cutoff) .or. (s(4,5) < cutoff)) goto 999
         if (abs(two*dot(p,3,4)+two*mt**2-zmass**2) .lt. 0.1_dp) then
           msq(:,:)=zip
           msqc(:,:,:)=zip
         else
c           call singcheck(qqb_QQb_mix_g,qqb_QQb_mix_gs,p)
           if (includereal) call qqb_QQb_mix_g(p,msq)
           call qqb_QQb_mix_gs(p,msqc)
           wt_noew=zip ! non-EW calculation has no real contribution
         endif
      elseif ((kcase==ktt_tot) .or. (kcase==kcc_tot)
     &   .or. (kcase==kbb_tot)) then
c        call singcheck(qqb_QQb_g,qqb_QQb_gs,p)
        if (includereal) call qqb_QQb_g(p,msq)
        call qqb_QQb_gs(p,msqc)
      elseif (kcase==ktopanom) then
        bbfrac = 0._dp

        msq_light = 0._dp
        msq_heavy = 0._dp
        msqc_light = 0._dp
        msqc_heavy = 0._dp
        ! this is required even when use_DDIS is .false.
        call singletop2_scale_setup(p)
        if (includereal) then
          ! setup once for real emission matrix elements
          call singletop2_real(p,msq_light,.true.,.false.)
          call singletop2_real(p,msq_heavy,.false.,.true.)
          msq = msq_light + msq_heavy
        endif
        ! scale setup for use_DDIS case is called within _gs and dipolesub
        call singletop2_gs(p,msqc_light,.true.,.false.)
        call singletop2_gs(p,msqc_heavy,.false.,.true.)
        msqc = msqc_light + msqc_heavy

        ! reset scale setup before singcheck, which calls the emission routine first
c       call singletop2_scale_setup(p)
c       call singcheck(singletop2_real_heavy, singletop2_gs_heavy, p)
      elseif (kcase==kbq_tpq) then
        bbfrac = 0._dp
c        call singcheck(qqb_tbb_g,qqb_tbb_gs,p)
        if (includereal) call qqb_tbb_g(p,msq)
        call qqb_tbb_gs(p,msqc)

      elseif (kcase==kt_bbar) then
c        call singcheck(qqb_tbb_g,qqb_tbb_gs,p)
        if (includereal) call qqb_tbbdk_g(p,msq)
        call qqb_tbbdk_gs(p,msqc)
      elseif (kcase==kttdkay) then
c        call singcheck(qqb_tbb_g_new,qqb_tbb_gs,p)
        if (includereal) call bq_tpq_gdk(p,msq)
        call bq_tpq_gsdk(p,msqc)
      elseif (kcase==ktdecay) then
c        call singcheck(qqb_tbb_g_new,qqb_tbb_gs,p)
      if (includereal) call dkqqb_tbbdk_g(p,msq)
        call dkqqb_tbbdk_gs(p,msqc)       
       elseif (kcase==kW_tndk) then
c        call singcheck(qqb_w_tndk_g,qqb_w_tndk_gs,p)      ! Checked 12/3/04
        if (includereal) call qqb_w_tndk_g(p,msq)
        call qqb_w_tndk_gs(p,msqc)
      elseif (kcase==kW_twdk) then
c        call singcheck(qqb_w_twdk_g,qqb_w_twdk_gs,p)            ! Checked 2/4/05
        if (includereal) call qqb_w_twdk_g(p,msq)
        call qqb_w_twdk_gs(p,msqc)
      elseif (kcase==kWtdkay) then
        if (includereal) call dkqqb_w_twdk_g(p,msq)
        call dkqqb_w_twdk_gs(p,msqc)
      elseif ( (kcase==kqq_ttw)) then 
        if (includereal) call qqb_ttw_g(p,msq)
        call qqb_ttw_gs(p,msqc)
      elseif ( (kcase==kttwldk)) then
        if     (decay1q2a == 1) then
c------ radiation in the decay of the top quark
        if (includereal) call dk1qqb_ttw_g(p,msq)
          call dk1qqb_ttw_gs(p,msqc)
      elseif (decay1q2a == 2) then 
c------ radiation in the decay of the anti-top quark
        if (includereal) call dk2qqb_ttw_g(p,msq)
          call dk2qqb_ttw_gs(p,msqc)
      else
        write(6,*) 'Problem generating PS, decay1q2a=',decay1q2a
        stop
      endif
      elseif (kcase==kggfus1) then
!        if (includereal) call singcheck(gg_hgg,gg_hg_gs,p)
        if (includereal) call gg_hgg(p,msq)
        call gg_hg_gs(p,msqc)
      elseif (kcase==kHi_Zaj) then
!        if (includereal) call singcheck(gg_hgg_zgam,gg_hg_zgam_gs,p)
        if (includereal) call gg_hgg_zgam(p,msq)
        call gg_hg_zgam_gs(p,msqc)
      elseif (kcase==khjetma) then
        badpoint = .false.

        if (includereal) call hjetmass_r(p,msq)
        if (badpoint .eqv. .false.) then
          call hjetmass_gs(p,msqc)
        else
            msq = 0._dp
            msqc = 0._dp
        endif
      elseif (kcase==kHgagaj) then
c        call singcheck(gg_hgagagg,gg_hgagag_gs,p)
        if (includereal) call gg_hgagagg(p,msq)
        call gg_hgagag_gs(p,msqc)
      elseif (kcase==kHWWjet) then
c        call singcheck(gg_hWWgg,gg_hWWg_gs,p)
        if (includereal) call gg_hWWgg(p,msq)
        call gg_hWWg_gs(p,msqc)
      elseif (kcase==kHWW2jt) then
c        call singcheck(gg_hWWgg,gg_hWWg_gs,p)
        if (includereal) call gg_hWWggg(p,msq)
        call gg_hWWgg_gs(p,msqc)
      elseif (kcase==kHZZjet) then
c        call singcheck(gg_hZZgg,gg_hZZg_gs,p)
        if (includereal) call gg_hZZgg(p,msq)
        call gg_hZZg_gs(p,msqc)
      elseif (kcase==kHZZ2jt) then
c        call singcheck(gg_hZZgg,gg_hZZg_gs,p)
        if (includereal) call gg_hZZggg(p,msq)
        call gg_hZZgg_gs(p,msqc)
c      write(6,*) msq
c      pause
      elseif (kcase==kqq_Hqq) then
c        call singcheck(VV_Hqq_g,VV_Hqq_gs,p)
        if (includereal) call VV_Hqq_g(p,msq)
        call VV_Hqq_gs(p,msqc)
      elseif (kcase==kqq_Hgg) then
c        call singcheck(VV_Hgaga_g,VV_Hgaga_gs,p)
        if (includereal) call VV_Hgaga_g(p,msq)
        call VV_Hgaga_gs(p,msqc)
      elseif (kcase==kqq_HWW) then
c        call singcheck(VV_HWW_g,VV_HWW_gs,p)
       if (includereal) call VV_HWW_g(p,msq)
        call VV_HWW_gs(p,msqc)
      elseif (kcase==kqq_HZZ) then
c        call singcheck(VV_HZZ_g,VV_HZZ_gs,p)
        if (includereal) call VV_HZZ_g(p,msq)
        call VV_HZZ_gs(p,msqc)
      elseif (kcase==kggfus2) then
c        call singcheck(gg_hggg,gg_hgg_gs,p)            ! Checked 10/29/09
        if (includereal) call gg_hggg(p,msq)
        call gg_hgg_gs(p,msqc)
      elseif (kcase==kgagajj) then
c        call singcheck(gg_hggg,gg_hgg_gs,p)
        if (includereal) call gg_hggg(p,msq)
        call gg_hgg_gs(p,msqc)
      elseif (kcase==kqg_tbq) then
c        call singcheck(qg_tbq_g,qg_tbq_gs,p)       ! Checked 2/4/08
       if (includereal) call qg_tbq_g(p,msq)
       call qg_tbq_gs(p,msqc)
      elseif (kcase==k4ftwdk) then
c        call singcheck(qg_tbq_g,qg_tbq_gs,p)
       if (includereal) call qg_tbqdk_g(p,msq)
       call qg_tbqdk_gs(p,msqc)
      elseif (kcase==kdk_4ft) then
c       if (includereal) call dkqg_tbqdk_g_old(p,msq)  
       if (includereal) call dkqg_tbqdk_g(p,msq)  
       call dkqg_tbqdk_gs(p,msqc)
      elseif (kcase==kqq_tbg) then
c        call singcheck(qq_tbg_g,qq_tbg_gs,p)       ! Checked 8/9/08
       if (includereal) call qq_tbg_g(p,msq)
       call qq_tbg_gs(p,msqc)
      elseif (kcase==kepem3j) then
c        call singcheck(epem3j_g,epem3j_gs,p)       ! Checked 17/11/08
       if (includereal) call epem3j_g(p,msq)
       call epem3j_gs(p,msqc)
      elseif (kcase==kgQ__ZQ) then
c        call singcheck(gQ_zQ_g,gQ_zQ_gs,p)
        if (includereal) call gQ_zQ_g(p,msq)
        call gQ_zQ_gs(p,msqc)
      elseif (kcase==kZ_bjet) then
c        call singcheck(qqb_zbjet_g,qqb_zbjet_gs,p)      ! Checked 07/18/05
        if (includereal) call qqb_zbjet_g(p,msq)
        call qqb_zbjet_gs(p,msqc)
      elseif (kcase==kW_bjet) then
c        call singcheck(qqb_wbjet_g,qqb_wbjet_gs,p) ! Rechecked 14/3/08
        if (includereal) call qqb_wbjet_g(p,msq)
        call qqb_wbjet_gs(p,msqc)
      elseif (kcase==kWcsbar) then
c        call singcheck(qqb_w_g,qqb_w_gs,p)         ! Checked 11/30/01
        if (includereal) call qqb_w_g(p,msq)      
        call qqb_w_gs(p,msqc)     
      elseif (kcase==kWcs_ms) then
        if (includereal) call qqb_w_cjet(p,msq)
        ndmax=0
      elseif (kcase==kZ_tjet) then
c        if (includereal) call singcheck(qq_tchan_ztqg,qq_tchan_ztq_gs,p)
        if (includereal) call qq_tchan_ztqg(p,msq)
        call qq_tchan_ztq_gs(p,msqc)
      elseif (kcase==kZ_tdkj) then
c        if (includereal) call singcheck(qq_tchan_ztqg_dk,
c     &                                  qq_tchan_ztq_dk_gs,p)
        if (includereal) call qq_tchan_ztqg_dk(p,msq)
        call qq_tchan_ztq_dk_gs(p,msqc)
      elseif (kcase==kH_tjet) then
c        if (includereal) call singcheck(qq_tchan_htqg,qq_tchan_htq_gs,p)
        if (includereal) call qq_tchan_htqg(p,msq)
        call qq_tchan_htq_gs(p,msqc)
      elseif (kcase==kH_tdkj) then
        if (includereal) call qq_tchan_htqg_dk(p,msq)
        call qq_tchan_htq_dk_gs(p,msqc)
      elseif (kcase==ktottth) then
c        if (includereal) call qqb_tottth_g(p,msq)
c        call qqb_tottth_gs(p,msqsc)
c        call singcheck(qqb_tottth_g,qqb_tottth_gs,p)
c        call compare_madgraph(p,qqb_tottth_g,qqb_tottth_g_mad)
      elseif (kcase==kdm_jet) then 
         if(includereal) call qqb_dm_monojet_g(p,msq)       
         call qqb_dm_monojet_gs(p,msqc)
      elseif (kcase==kdm_gam) then
!         if(includereal) call singcheck(qqb_dm_monophot_g
!     &        ,qqb_dm_monophot_gs,p)
         if(includereal) call qqb_dm_monophot_g(p,msq) 
         call qqb_dm_monophot_gs(p,msqc)

      endif
      
! code to find power of alpha-s for scale variation
      if ((doscalevar) .and. (foundpow .eqv. .false.)) then
        if (itrial == 1) then
          msqtrial=maxval(msq)
          if (msqtrial == 0) goto 999
          if (dynamicscale) call scaleset(initscale,initfacscale,p)
          as=as*two
          ason2pi=ason2pi*two
          ason4pi=ason4pi*two
          gsq=gsq*two
          itrial=itrial+1
          goto 66
        endif
        msqtrial=maxval(msq)/msqtrial
        alphaspow=-1
        if (abs(msqtrial-one) < 1.e-8) alphaspow=0
        if (abs(msqtrial-two) < 1.e-8) alphaspow=1
        if (abs(msqtrial-four) < 1.e-8) alphaspow=2
        if (abs(msqtrial-eight) < 1.e-8) alphaspow=3
        if (abs(msqtrial-16._dp) < 1.e-8) alphaspow=4
        if (abs(msqtrial-32._dp) < 1.e-8) alphaspow=5
        if (alphaspow == -1) then
          write(6,*) 'Unable to determine power of alpha-s for scale variation'
          write(6,*) ' msqtrial = ',msqtrial
          stop
        endif
        as=as/two
        ason2pi=ason2pi/two
        ason4pi=ason4pi/two
        gsq=gsq/two
        foundpow=.true.
!        write(6,*) 'Found alpha-s power: ',alphaspow
        goto 66
      endif

      xmsq = 0._dp
      xmsq_bypart = 0._dp
      bbfrac = 0._dp
 
c--- APPLgrid - initialize array
      if (creategrid.and.bin) then
       do nd=0,ndmax
         do j=-nf,nf
            do k=-nf,nf
               weightr(nd,j,k) = 0d0
            enddo
         enddo
       enddo  
      weightfactor = 1d0
      endif
c--- APPLgrid - end

      currentPDF=0
            
      flux=fbGeV2/(two*xx1*xx2*W)
c--- for mlm study, divide by (Ecm)**2=W
c      if (runstring(1:3) == 'mlm') then
c      flux=flux/W
c      endif

c--- initialize a PDF set here, if calculating errors
  777 continue

      xmsq = 0._dp
      bbfrac = 0._dp

c--- calculate PDF's

      if (PDFerrors) then
!$omp critical(critPDFerrors)
        if (dynamicscale) then
          do nd=ndmax,0,-1  ! so that fx1,fx2 correct for real kinematics
            if (dipscale(nd) < 1.e-8_dp) then        
c--- in case dipole is not used, set up dummy value of scale for safety
c--- and set all PDF entries to zero
              dipscale(nd)=dipscale(0)
              fx1(:)=0._dp
              fx2(:)=0._dp
            else
              call InitPDF(currentPDF)
              call fdist(ih1,xx1,dipscale(nd),fx1)
              call fdist(ih2,xx2,dipscale(nd),fx2)
              dipfx1(nd,:)=fx1(:)
              dipfx2(nd,:)=fx2(:)
            endif
          enddo
          if ((kcase==kqg_tbq) .or. (kcase==k4ftwdk)) then
            fx1_H(:)=fx1(:)
            fx1_L(:)=fx1(:)
            fx2_H(:)=fx2(:)
            fx2_L(:)=fx2(:)
          endif
        else
          if ((kcase==kqg_tbq) .or. (kcase==k4ftwdk)) then
c--- for single top + b, make sure to use two different scales
            call InitPDF(currentPDF)
            call fdist(ih1,xx1,facscale_H,fx1_H)
            call fdist(ih2,xx2,facscale_H,fx2_H)
            call fdist(ih1,xx1,facscale_L,fx1_L)
            call fdist(ih2,xx2,facscale_L,fx2_L)
            do j=-nf,nf
              if (j == 0) then  ! heavy quark line has gluon init. state
                fx1(j)=fx1_H(j)
                fx2(j)=fx2_H(j)
              else
                fx1(j)=fx1_L(j)
                fx2(j)=fx2_L(j)
              endif
            enddo
          else
            call InitPDF(currentPDF)
            call fdist(ih1,xx1,facscale,fx1)
            call fdist(ih2,xx2,facscale,fx2)
          endif
        endif
!$omp end critical(critPDFerrors)

      else ! case of no PDFerrors
        if (dynamicscale) then
          do nd=ndmax,0,-1  ! so that fx1,fx2 correct for real kinematics
            if (dipscale(nd) < 1.e-8_dp) then        
c--- in case dipole is not used, set up dummy value of scale for safety
c--- and set all PDF entries to zero
              dipscale(nd)=dipscale(0)
              fx1(:)=0._dp
              fx2(:)=0._dp
            else
              if (doscalevar) then
                call fdist(ih1,xx1,dipscale(nd)*two,fx1)
                call fdist(ih2,xx2,dipscale(nd)*two,fx2)
                dipfx1up(nd,:)=fx1(:)
                dipfx2up(nd,:)=fx2(:)
                call fdist(ih1,xx1,dipscale(nd)/two,fx1)
                call fdist(ih2,xx2,dipscale(nd)/two,fx2)
                dipfx1dn(nd,:)=fx1(:)
                dipfx2dn(nd,:)=fx2(:)
                xmsqvar(:,nd)=zip
              endif
              if (kcase == ktopanom .and. use_DDIS) then
                ! first part: pdf's for dipole contributions
                if (nd > 0) then
                    if (singletop2_dipscale(nd,1) /= 0._dp .and. singletop2_dipscale(nd,2) /= 0._dp) then
                        call fdist(ih1,xx1,singletop2_dipscale(nd,1),singletop2_dipole_pdfs(nd,1,:))
                        call fdist(ih2,xx2,singletop2_dipscale(nd,2),singletop2_dipole_pdfs(nd,2,:))
                    else
                        if (sum(abs(msqc(nd,:,:))) > 0._dp) then
                            write (*,*) "CRITICAL WARNING: dipole contribution ", sum(abs(msqc(nd,:,:))), " but scale zero"
                        endif
                    endif
                endif
              else ! usual case
                call fdist(ih1,xx1,dipscale(nd),fx1)
                call fdist(ih2,xx2,dipscale(nd),fx2)
                dipfx1(nd,:)=fx1(:)
                dipfx2(nd,:)=fx2(:)
              endif
            endif
          enddo
          if (kcase == ktopanom .and. use_DDIS) then
            ! scales were previously set by some dipole routine, reset with unmodified momenta
            call singletop2_scale_setup(p)
            ! in addition we need all the PDF sets for the real emission matrix elements
            call fdist(ih1,xx1,facscale_beam1_islight_onlight, singletop2_pdfs(i_beam1_islight_onlight,:))
            call fdist(ih2,xx2,facscale_beam2_isheavy_onlight, singletop2_pdfs(i_beam2_isheavy_onlight,:))

            call fdist(ih1,xx1,facscale_beam1_isheavy_onlight, singletop2_pdfs(i_beam1_isheavy_onlight,:))
            call fdist(ih2,xx2,facscale_beam2_islight_onlight, singletop2_pdfs(i_beam2_islight_onlight,:))

            call fdist(ih1,xx1,facscale_beam1_islight_onheavy, singletop2_pdfs(i_beam1_islight_onheavy,:))
            call fdist(ih2,xx2,facscale_beam2_isheavy_onheavy, singletop2_pdfs(i_beam2_isheavy_onheavy,:))

            call fdist(ih1,xx1,facscale_beam1_isheavy_onheavy, singletop2_pdfs(i_beam1_isheavy_onheavy,:))
            call fdist(ih2,xx2,facscale_beam2_islight_onheavy, singletop2_pdfs(i_beam2_islight_onheavy,:))
          endif
          if ((kcase==kqg_tbq) .or. (kcase==k4ftwdk)) then
            fx1_H(:)=fx1(:)
            fx1_L(:)=fx1(:)
            fx2_H(:)=fx2(:)
            fx2_L(:)=fx2(:)
          endif
        else ! no dynamic scale
          if ((kcase==kqg_tbq) .or. (kcase==k4ftwdk)) then
c--- for single top + b, make sure to use two different scales
            call fdist(ih1,xx1,facscale_H,fx1_H)
            call fdist(ih2,xx2,facscale_H,fx2_H)
            call fdist(ih1,xx1,facscale_L,fx1_L)
            call fdist(ih2,xx2,facscale_L,fx2_L)
            do j=-nf,nf
              if (j == 0) then  ! heavy quark line has gluon init. state
                fx1(j)=fx1_H(j)
                fx2(j)=fx2_H(j)
              else
                fx1(j)=fx1_L(j)
                fx2(j)=fx2_L(j)
              endif
            enddo
          else
c--- usual case            
            call fdist(ih1,xx1,facscale,fx1)
            call fdist(ih2,xx2,facscale,fx2)
            if (doscalevar) then
              call fdist(ih1,xx1,facscale*two,fx1up)
              call fdist(ih2,xx2,facscale*two,fx2up)
              call fdist(ih1,xx1,facscale/two,fx1dn)
              call fdist(ih2,xx2,facscale/two,fx2dn)
              xmsqvar(:,:)=zip
            endif
          endif
        endif

      endif

      do j=-nflav,nflav
      do k=-nflav,nflav

c     ! only gluon quark, this selects only the bbfrac piece
c     if (.not. ((j==0 .and. abs(k) /= 5) .or. (abs(j) /= 5 .and. k==0))) cycle

      if (ggonly) then
      if ((j.ne.0) .or. (k.ne.0)) cycle
      endif

      if (gqonly) then
      if (((j==0).and.(k==0)) .or. ((j.ne.0).and.(k.ne.0))) cycle
      endif      
      
      if (noglue) then 
      if ((j==0) .or. (k==0)) cycle
      endif

      if (omitgg) then 
      if ((j==0) .and. (k==0)) cycle
      endif

      if ((kcase==kWcsbar).and.(j .ne. 4).and.(k .ne. 4)) cycle

      if     (j > 0) then
        sgnj=+1
      elseif (j < 0) then
        sgnj=-1
      else
        sgnj=0
      endif
      if     (k > 0) then
        sgnk=+1
      elseif (k < 0) then
        sgnk=-1
      else
        sgnk=0
      endif

c     write (*,*) "scale checks"
c     write (*,*) "dipole 1", singletop2_dipscale(1,1), singletop2_dipscale(1,2)
c     write (*,*) "dipole 2", singletop2_dipscale(2,1), singletop2_dipscale(2,2)
c     write (*,*) "dipole 3", singletop2_dipscale(3,1), singletop2_dipscale(3,2)
c     write (*,*) "dipole 4", singletop2_dipscale(4,1), singletop2_dipscale(4,2)
c     write (*,*) "dipscales", dipscale(1), dipscale(2)
c     write (*,*) "dipolesum", sum(msqc(1,:,:)), sum(msqc(2,:,:)), sum(msqc(3,:,:)), sum(msqc(4,:,:))
c     call singletop2_scale_setup(p)
c     write (*,*) "msq ", renscale_beam1_islight_onlight, renscale_beam2_isheavy_onlight
c     write (*,*) ""

c--- for single top + b, make sure to use two different scales
      if ((kcase==kqg_tbq) .or. (kcase==k4ftwdk)) then
        xmsqjk=fx1_L(j)*fx2_H(k)*msqLH(j,k)
     &        +fx1_H(j)*fx2_L(k)*msqHL(j,k)
      else
        xmsqjk = 0._dp

        if (kcase == ktopanom .and. use_DDIS) then
          if (abs(k) == 5) then ! qb, gb
            xmsqjk = xmsqjk + singletop2_pdfs(i_beam1_islight_onlight,j) *
     &              singletop2_pdfs(i_beam2_isheavy_onlight,k) * msq_light(j,k)
            xmsqjk = xmsqjk + singletop2_pdfs(i_beam1_islight_onheavy,j) *
     &              singletop2_pdfs(i_beam2_isheavy_onheavy,k) * msq_heavy(j,k)
          elseif (abs(j) == 5) then ! bq, bg
            xmsqjk = xmsqjk + singletop2_pdfs(i_beam1_isheavy_onlight,j) *
     &              singletop2_pdfs(i_beam2_islight_onlight,k) * msq_light(j,k)
            xmsqjk = xmsqjk + singletop2_pdfs(i_beam1_isheavy_onheavy,j) *
     &              singletop2_pdfs(i_beam2_islight_onheavy,k) * msq_heavy(j,k)
          elseif (abs(j) /= 5 .and. k == 0) then ! this is qg
            xmsqjk = xmsqjk + singletop2_pdfs(i_beam1_islight_onheavy,j) *
     &              singletop2_pdfs(i_beam2_isheavy_onheavy,k) * msq_heavy(j,k)
          elseif (j == 0 .and. abs(k) /= 5) then ! this is crossed qg
            xmsqjk = xmsqjk + singletop2_pdfs(i_beam1_isheavy_onheavy,j) *
     &              singletop2_pdfs(i_beam2_islight_onheavy,k) * msq_heavy(j,k)
          else ! no other cases for our process
            if (msq(j,k) /= 0._dp) then
              write (*,*) "CRITICAL WARNING: unexpected channel contribution", j,k, msq(j,k)
            endif
            xmsqjk = 0._dp
          endif
        else ! usual case
          xmsqjk=fx1(j)*fx2(k)*msq(j,k)
        endif

      xmsq(0)=xmsq(0)+xmsqjk

      if (doscalevar) then
        if (dynamicscale) then
          xmsqvar(1,0)=xmsqvar(1,0)+dipfx1up(0,j)*dipfx2up(0,k)*msq(j,k)
          xmsqvar(2,0)=xmsqvar(2,0)+dipfx1dn(0,j)*dipfx2dn(0,k)*msq(j,k)
        else
          xmsqvar(1,0)=xmsqvar(1,0)+fx1up(j)*fx2up(k)*msq(j,k)
          xmsqvar(2,0)=xmsqvar(2,0)+fx1dn(j)*fx2dn(k)*msq(j,k)
        endif
      endif

      ! ktopanom: two b's, real emission contribution
      if (nproc==169 .or. nproc==164 .or. nproc==161 .or. nproc==166) then
        if ((j == 0 .and. abs(k) /= 5) .or. (abs(j) /= 5 .and. k ==0)) then
            bbfrac(0) = bbfrac(0) + xmsqjk
        endif
      endif

         if (currentPDF == 0) then
           xmsq_bypart(0,sgnj,sgnk)=xmsq_bypart(0,sgnj,sgnk)+xmsqjk
c--- APPLgrid - save weight nd = 0 
           if (creategrid.and.bin) then
              weightr(0,j,k) = weightr(0,j,k) + msq(j,k)
           endif
c--- APPLgrid - end
         endif

         ! add dipole contributions
         do nd=1,ndmax
         if (dynamicscale) then
             if (kcase == ktopanom .and. use_DDIS) then
                 xmsqjk = singletop2_dipole_pdfs(nd,1,j) * singletop2_dipole_pdfs(nd,2,k) * (-msqc(nd,j,k))
c                write (*,*) j,k,singletop2_dipole_pdfs(nd,1,j),singletop2_dipole_pdfs(nd,2,k),(-msqc(nd,j,k))
             else
               xmsqjk=dipfx1(nd,j)*dipfx2(nd,k)*(-msqc(nd,j,k))
               if (doscalevar) then
                 xmsqvar(1,nd)=xmsqvar(1,nd)+dipfx1up(nd,j)*dipfx2up(nd,k)*(-msqc(nd,j,k))
                 xmsqvar(2,nd)=xmsqvar(2,nd)+dipfx1dn(nd,j)*dipfx2dn(nd,k)*(-msqc(nd,j,k))
               endif
             endif
         else ! no dynamicscale
             xmsqjk=fx1(j)*fx2(k)*(-msqc(nd,j,k))
             if (doscalevar) then
               xmsqvar(1,nd)=xmsqvar(1,nd)+fx1up(j)*fx2up(k)*(-msqc(nd,j,k))
               xmsqvar(2,nd)=xmsqvar(2,nd)+fx1dn(j)*fx2dn(k)*(-msqc(nd,j,k))
             endif
         endif ! dynamicscale

      ! ktopanom: two b's, dipole contribution
      if (nproc==169 .or. nproc==164 .or. nproc==161 .or. nproc==166) then
c       if (nd == 13 .or. nd == 14) the
          if ((j == 0 .and. abs(k) /= 5) .or. (abs(j) /= 5 .and. k ==0)) then
              bbfrac(nd) = bbfrac(nd) + xmsqjk
          endif
c       endif
      endif

         xmsq(nd)=xmsq(nd)+xmsqjk
         if (currentPDF == 0) then
           xmsq_bypart(nd,sgnj,sgnk)=xmsq_bypart(nd,sgnj,sgnk)+xmsqjk
c--- APPLgrid - save weight nd = 1,ndmax 
           if (creategrid.and.bin) then
              weightr(nd,j,k) = weightr(nd,j,k) + (-msqc(nd,j,k))
           endif
c--- APPLgrid - end
         endif
         enddo
         
      endif

      enddo
      enddo

c--- loop over all PDF error sets, if necessary
      if (PDFerrors) then
        do nd=0,ndmax
          PDFxsec_nd(currentPDF,nd)=xmsq(nd)
        enddo
        currentPDF=currentPDF+1
        if (currentPDF <= maxPDFsets) goto 777
c--- reset xmsq to the central PDF values
        do nd=0,ndmax
          xmsq(nd)=PDFxsec_nd(0,nd)
        enddo
      endif

      realint=0._dp
      xint=0._dp

      valsum=0._dp ! running total of weights at this point

c--- zero out temporary histograms
c      if (bin) call zerorealhistos
      if (bin) call smartzero

c---trial with weight of real alone
c---first set up all dipole contributions
c---this is the value of integral including subtractions
      do nd=0,ndmax

        ! bbfrac, compute actual fraction
        if (kcase == ktopanom .or. kcase == kbq_tpq) then
          if (xmsq(nd) /= 0._dp) then
              bbfrac(nd) = bbfrac(nd)/xmsq(nd)
          endif
        endif

        xmsq(nd)=xmsq(nd)*flux*pswt/BrnRat
c        if (creatent) then
          wt_gg=xmsq_bypart(nd,0,0)*wgt*flux*pswt/BrnRat/real(itmx,dp)
          wt_gq=(xmsq_bypart(nd,+1,0)+xmsq_bypart(nd,-1,0)
     &          +xmsq_bypart(nd,0,+1)+xmsq_bypart(nd,0,-1)
     &          )*wgt*flux*pswt/BrnRat/real(itmx,dp)
          wt_qq=(xmsq_bypart(nd,+1,+1)+xmsq_bypart(nd,-1,-1)
     &          )*wgt*flux*pswt/BrnRat/real(itmx,dp)
          wt_qqb=(xmsq_bypart(nd,+1,-1)+xmsq_bypart(nd,-1,+1)
     &          )*wgt*flux*pswt/BrnRat/real(itmx,dp)
c        endif
        failed=.false.
        
        if (nd == 0) then
c---if there's no real contribution, record the event as failing to pass cuts
          if (xmsq(nd) == 0._dp) then
             failed=.true.
             goto 996
          endif
        else
c--- if this dipole has no contribution, go to end of loop
          if (xmsq(nd) == 0._dp) goto 997         
c---check whether each counter-event passes the cuts
          do j=1,mxpart
          do k=1,4
          q(j,k)=ptilde(nd,j,k)
          enddo
          enddo

          if (incldip(nd)) incldip(nd)=includedipole(nd,q)
          if (incldip(nd) .eqv. .false.) failed=.true.
          
c          if (incldip(nd)) then
c            failed=.false.
c            call smalltau(q,3,*118)
c            goto 119
c 118        failed=.true.
c 119        continue
c          endif
        endif

 996    if (failed) then
          if (nd == 0) then
!$omp atomic
            ncutzero=ncutzero+1
!$omp atomic
            ntotzero=ntotzero+1
          endif
          call dotem(nvec,p,s)
          xmsq(nd)=0._dp
          goto 997
        endif

c---if it does, add to total
!        xint=xint+xmsq(nd)
c---  SSbegin
        xint=xint+xmsq(nd)*reweight
c---  SSend

        val=wgt*flux*pswt/BrnRat
        do j=-1,1
        do k=-1,1
!$omp atomic
          lord_bypart(j,k)=lord_bypart(j,k)+
     &         val*xmsq_bypart(nd,j,k)
        enddo
        enddo

        val=xmsq(nd)*wgt
        val2=val**2

      valsum=valsum+val

c--- update PDF errors
        if (PDFerrors) then
          do currentPDF=0,maxPDFsets
          PDFwgt(currentPDF)=
     &       flux*pswt*PDFxsec_nd(currentPDF,nd)/BrnRat*wgt/itmx
!$omp atomic
          PDFxsec(currentPDF)=PDFxsec(currentPDF)
     &       +PDFwgt(currentPDF)
          enddo
        endif

c--- catch NaN before it enters histograms
        if (ieee_is_nan(val)) then
          write(6,*) 'discarding point with weight NaN: pswt=',pswt
          goto 999
        endif 

        if (.not. ieee_is_finite(val)) then
          write(6,*) 'discarding point with infinite weight: pswt=',pswt
          goto 999
        endif

                
! compute weights for scale variation
      if (doscalevar) then
        if (abs(xmsq(nd)) > zip) then
          if (dynamicscale) then
            asorig=alphas(dipscale(nd)*initscale/initfacscale,amz,nlooprun)
            scaleup=dipscale(nd)*initscale/initfacscale*two
            scaledn=dipscale(nd)*initscale/initfacscale/two
          else
            asorig=as
            scaleup=scale*two
            scaledn=scale/two
          endif
          scalereweight(1)=(alphas(scaleup,amz,nlooprun)/asorig)**alphaspow
          scalereweight(2)=(alphas(scaledn,amz,nlooprun)/asorig)**alphaspow
          scalereweight(1)=scalereweight(1)*xmsqvar(1,nd)*flux*pswt/BrnRat/xmsq(nd)
          scalereweight(2)=scalereweight(2)*xmsqvar(2,nd)*flux*pswt/BrnRat/xmsq(nd)
          if (maxscalevar == 6) then
            scalereweight(3)=scalereweight(1)*xmsq(nd)/xmsqvar(1,nd)/flux/pswt*BrnRat
            scalereweight(4)=scalereweight(2)*xmsq(nd)/xmsqvar(2,nd)/flux/pswt*BrnRat
            scalereweight(5)=xmsqvar(1,nd)*flux*pswt/BrnRat/xmsq(nd)
            scalereweight(6)=xmsqvar(2,nd)*flux*pswt/BrnRat/xmsq(nd)
          endif
        else
          scalereweight(:)=zip
        endif
      endif

c---if we're binning, add to histo too
        if (bin) then
c--- APPLgrid - writing out the common block
           if (creategrid) then

              if (dynamicscale) then
                currason2pi = alphas(dipscale(nd),amz,nlooprun)/twopi
              else
                currason2pi = ason2pi
              endif

              psCR = (1d0/currason2pi)**(lowest_order+1)
              do j=-nflav,nflav
                 do k=-nflav,nflav
                    weightr(nd,j,k)=weightr(nd,j,k)*psCR
                 enddo
              enddo
              contrib = 200
              dipole  = nd
              weightfactor = wgt*flux*pswt/BrnRat/dfloat(itmx)
c              weightfactor = weightfactor * eight*pisq/gsq
              ag_xx1       = xx1
              ag_xx2       = xx2
              if (dynamicscale) then
                ag_scale     = dipscale(nd)
              else
                ag_scale     = facscale
              endif
              refwt        = val/dfloat(itmx)
              refwt2       = val2/dfloat(itmx)
C               print*,"*************************************************"
C               print*, "meWeightFactor = ", weightfactor," D =
C               ",dipole,
C      *       " me(0, 2,-1) = " ,  weightr(0, 2 ,-1) ," ", msq(2,-1),
C      *       " me(0, -1,2) = " ,  weightr(0, -1 ,2) ," ", msq(-1,2),
C      *       " me(1, 2,-1) = " ,  weightr(1, 2 ,-1) ," ",
C      -msqc(1,2,-1),
C      *       " me(1, -1,2) = " ,  weightr(1, -1 ,2) ," ",
C      -msqc(1,-1,2)
C               print*, " x1 = ",xx1," x2 = ",xx2," sca = ",facscale
C               print *, "rewt = ", refwt
C               print*,"*************************************************"
C               flush(6)
           endif
c--- APPLgrid - end
C          print *, "In realint pjet = ", pjet    
          call getptildejet(nd,pjet)
          call dotem(nvec,pjet,s)
          call nplotter(pjet,val,val2,nd)
c--- POWHEG-style output if requested
          if (writepwg) then
            if (nd == 0) then 
!$omp critical(pwhgplot)
              call pwhgplotter(p,pjet,val,nd)
!$omp end critical(pwhgplot)
            else
             do j=1,mxpart
              do k=1,4
              q(j,k)=ptilde(nd,j,k)
              enddo
              enddo
!$omp critical(pwhgplot)
              call pwhgplotter(q,pjet,val,nd)
!$omp end critical(pwhgplot)
            endif
          endif
        endif
c---otherwise, skip contribution
 997    continue
      enddo

c 998  continue

c--- add temporary histograms to cumulative totals
      if (bin) then
         call smartadd(wgt) 
      endif
c--- update the maximum weight so far, if necessary
      if (abs(valsum) > wtmax) then
!$omp critical(MaxWgt)
        wtmax=abs(valsum)
!$omp end critical(MaxWgt)
      endif
      
      if (realwt) then
      realint=xmsq(0)
      else
      realint=xint
      endif
      xreal=xreal+xint*wgt/real(itmx,dp)
      xreal2=xreal2+(xint*wgt)**2/real(itmx,dp)
      
c--- handle special case of Qflag and Gflag
      if (QandGflag) then
        QandGint=QandGint+realint
        if ((Gflag) .and. (.not.(Qflag))) then
c--- go back for second pass (Qflag), calling includedipole first to reset
c--- the values of "jets" and "jetlabel"
          includereal=includedipole(0,p)
          Qflag=.true.
          Gflag=.false.
          goto 44
        else
c--- return both to .true. and assign value to realint (to return to VEGAS)
          Qflag=.true.
          Gflag=.true.
          realint=QandGint
        endif
      endif

      return

 999  realint=0._dp
!$omp atomic
      ntotzero=ntotzero+1
c--- safety catch
      if (QandGflag) then
        Qflag=.true.
        Gflag=.true.
      endif      
      return
      end

