      module hjetmass_r_amps_dp_mod

      private
      public :: hjetmass_HqqbQQb_ident_dp, hjetmass_HqqbQQB_nonid_dp,
     &   hjetmass_Hgggg_dp, hjetmass_Hqqbgg_dp

      contains
      
      complex*16 function
     &    hjetmass_HqqbQQb_amp_pmpm_dp(i1,i2,i3,i4,p,za,zb,
     &                              triints, bubints)
        implicit none
        include 'types.f'
        include 'mxpart.f'
        complex*16 check, rat
        integer i1,i2,i3,i4
        include 'constants.f'
        include 'masses.f'
        complex*16 za(mxpart,mxpart), zb(mxpart,mxpart)
        real*8 p(mxpart,4)

        complex*16 triints(1)
        complex*16 bubints(3)

        complex*16 hjetmass_qqbQQb_triangle_pmpm_s34_mhsq_s12_rat_dp
        complex*16 hjetmass_qqbQQb_triangle_pmpm_s34_mhsq_s12_dp
        complex*16 hjetmass_qqbQQb_bubble_pmpm_s34_dp
        complex*16 hjetmass_qqbQQb_bubble_pmpm_s12_dp
        complex*16 hjetmass_qqbQQb_bubble_pmpm_mhsq_dp

        check = 0d0
        rat = 0d0
      check = check + hjetmass_qqbqqb_triangle_pmpm_s34_mhsq_s12_dp
     & (i1,i2,i3,i4,za,zb,mt)
     & * triints(1)
      check = check + hjetmass_qqbqqb_bubble_pmpm_mhsq_dp
     & (i1,i2,i3,i4,za,zb)
     & * bubints(1)
      check = check + hjetmass_qqbqqb_bubble_pmpm_s34_dp
     & (i1,i2,i3,i4,za,zb)
     & * bubints(2)
      check = check + hjetmass_qqbqqb_bubble_pmpm_s12_dp
     & (i1,i2,i3,i4,za,zb)
     & * bubints(3)

      rat = hjetmass_qqbQQb_triangle_pmpm_s34_mhsq_s12_rat_dp
     &   (i1,i2,i3,i4,za,zb) / 2d0

      hjetmass_HqqbQQb_amp_pmpm_dp = -(check + rat)
      return

      end function hjetmass_HqqbQQb_amp_pmpm_dp

      complex*16 function
     &    hjetmass_HqqbQQb_amp_pmmp_dp(i1,i2,i3,i4,p,za,zb,
     &                              triints, bubints)
        implicit none
        include 'types.f'
        include 'mxpart.f'
        complex*16 check, rat
        integer i1,i2,i3,i4
        include 'constants.f'
        include 'masses.f'
        complex*16 za(mxpart,mxpart), zb(mxpart,mxpart)
        real*8 p(mxpart,4)

        complex*16 triints(1)
        complex*16 bubints(3)

        complex*16 hjetmass_qqbQQb_triangle_pmmp_s34_mhsq_s12_rat_dp
        complex*16 hjetmass_qqbQQb_triangle_pmmp_s34_mhsq_s12_dp
        complex*16 hjetmass_qqbQQb_bubble_pmmp_s34_dp
        complex*16 hjetmass_qqbQQb_bubble_pmmp_s12_dp
        complex*16 hjetmass_qqbQQb_bubble_pmmp_mhsq_dp

        check = 0d0
        rat = 0d0

      check = check + hjetmass_qqbqqb_triangle_pmmp_s34_mhsq_s12_dp
     & (i1,i2,i3,i4,za,zb,mt)
     & * triints(1)
      check = check + hjetmass_qqbqqb_bubble_pmmp_mhsq_dp
     & (i1,i2,i3,i4,za,zb)
     & * bubints(1)
      check = check + hjetmass_qqbqqb_bubble_pmmp_s34_dp
     & (i1,i2,i3,i4,za,zb)
     & * bubints(2)
      check = check + hjetmass_qqbqqb_bubble_pmmp_s12_dp
     & (i1,i2,i3,i4,za,zb)
     & * bubints(3)

      rat = hjetmass_qqbQQb_triangle_pmmp_s34_mhsq_s12_rat_dp
     &   (i1,i2,i3,i4,za,zb) / 2d0

      hjetmass_HqqbQQb_amp_pmmp_dp = (check + rat)
      return

      end function hjetmass_HqqbQQb_amp_pmmp_dp
      
      complex*16 function
     &    hjetmass_Hqqbgg_amp_pmpp_dp(i1,i2,i3,i4,p,za,zb,flip,
     &                              boxints, triints, bubints)
        implicit none
        include 'types.f'
        include 'mxpart.f'
        complex*16 check, rat, tmpcast
        integer i1,i2,i3,i4
        include 'constants.f'
        include 'masses.f'
        complex*16 za(mxpart,mxpart), zb(mxpart,mxpart)
        real*8 p(mxpart,4)
        logical flip

        complex*16 boxints(3)
        complex*16 triints(6)
        complex*16 bubints(5)

        complex*16 hjetmass_qqbgg_box_pmpp_0_0_s12_mhsq_s34_s124_dp
        complex*16 hjetmass_qqbgg_box_pmpp_0_0_s12_mhsq_s34_s123_dp
        complex*16 hjetmass_qqbgg_box_pmpp_0_s12_0_mhsq_s123_s124_dp

        complex*16 hjetmass_qqbgg_triangle_pmpp_s34_mhsq_s12_dp
        complex*16 hjetmass_qqbgg_triangle_pmpp_s124_mhsq_0_dp
        complex*16 hjetmass_qqbgg_triangle_pmpp_s123_mhsq_0_dp
        complex*16 hjetmass_qqbgg_triangle_pmpp_s12_s124_0_dp
        complex*16 hjetmass_qqbgg_triangle_pmpp_s12_s123_0_dp

        complex*16 hjetmass_qqbgg_triangle_pmpp_s34_mhsq_s12_rat_dp
        complex*16 hjetmass_qqbgg_triangle_pmpp_s124_mhsq_0_rat_dp
        complex*16 hjetmass_qqbgg_triangle_pmpp_s123_mhsq_0_rat_dp

        complex*16 hjetmass_qqbgg_bubble_pmpp_s12_dp
        complex*16 hjetmass_qqbgg_bubble_pmpp_s123_dp
        complex*16 hjetmass_qqbgg_bubble_pmpp_s124_dp
        complex*16 hjetmass_qqbgg_bubble_pmpp_mhsq_dp

        check = (0d0,0d0)
        rat = (0d0,0d0)

        check = check +
     &  hjetmass_qqbgg_box_pmpp_0_0_s12_mhsq_s34_s124_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(2)
        check = check +
     &  hjetmass_qqbgg_box_pmpp_0_0_s12_mhsq_s34_s123_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(1)
        check = check +
     &  hjetmass_qqbgg_box_pmpp_0_s12_0_mhsq_s123_s124_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(3)
        check = check +
     &  hjetmass_qqbgg_triangle_pmpp_s34_mhsq_s12_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(6)

        tmpcast =
     &  hjetmass_qqbgg_triangle_pmpp_s34_mhsq_s12_rat_dp
     &  (i1,i2,i3,i4,za,zb)/2d0
        rat = rat - tmpcast

        check = check +
     &  hjetmass_qqbgg_triangle_pmpp_s124_mhsq_0_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(2)

        tmpcast =
     &  hjetmass_qqbgg_triangle_pmpp_s124_mhsq_0_rat_dp
     &  (i1,i2,i3,i4,za,zb)/2d0
        rat = rat - tmpcast

        check = check +
     &  hjetmass_qqbgg_triangle_pmpp_s12_s124_0_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(4)

        check = check +
     &  hjetmass_qqbgg_triangle_pmpp_s123_mhsq_0_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(1)

        tmpcast =
     &  hjetmass_qqbgg_triangle_pmpp_s123_mhsq_0_rat_dp
     &  (i1,i2,i3,i4,za,zb)/2d0
        rat = rat - tmpcast

        check = check +
     &  hjetmass_qqbgg_triangle_pmpp_s12_s123_0_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(3)

        tmpcast =
     &  hjetmass_qqbgg_bubble_pmpp_mhsq_dp
     &  (i1,i2,i3,i4,za,zb,mt,p,flip)
        check = check + tmpcast
     &  * bubints(1)

        tmpcast =
     &  hjetmass_qqbgg_bubble_pmpp_s12_dp
     &  (i1,i2,i3,i4,za,zb,mt)
        check = check + tmpcast
     &  * bubints(4)

        tmpcast =
     &  hjetmass_qqbgg_bubble_pmpp_s124_dp
     &  (i1,i2,i3,i4,za,zb,mt,p,flip)
        check = check + tmpcast
     &  * bubints(3)

        tmpcast =
     &  hjetmass_qqbgg_bubble_pmpp_s123_dp
     &  (i1,i2,i3,i4,za,zb,mt,p,flip)
        check = check + tmpcast
     &  * bubints(2)

        !hjetmass_Hqqbgg_amp_pmpp = -(0d0,6d0)*mt**2*(check + rat)
        hjetmass_Hqqbgg_amp_pmpp_dp = (check + rat)
        return

      end function

      complex*16 function
     &    hjetmass_Hqqbgg_amp_pmpm_dp(i1,i2,i3,i4,p,za,zb,flip,
     &                              boxints, triints, bubints)
        implicit none
        include 'types.f'
        include 'mxpart.f'
        complex*16 check, rat, tmpcast
        integer i1,i2,i3,i4
        include 'constants.f'
        include 'masses.f'
        complex*16 za(mxpart,mxpart), zb(mxpart,mxpart)
        real*8 p(mxpart,4)
        logical flip

        complex*16 boxints(3)
        complex*16 triints(6)
        complex*16 bubints(5)

        complex*16 hjetmass_qqbgg_box_pmpm_0_0_s12_mhsq_s34_s124_dp
        complex*16 hjetmass_qqbgg_box_pmpm_0_0_s12_mhsq_s34_s123_dp
        complex*16 hjetmass_qqbgg_box_pmpm_0_s12_0_mhsq_s123_s124_dp
  
        complex*16 hjetmass_qqbgg_triangle_pmpm_s34_mhsq_s12_dp
        complex*16 hjetmass_qqbgg_triangle_pmpm_s124_mhsq_0_dp
        complex*16 hjetmass_qqbgg_triangle_pmpm_s123_mhsq_0_dp
        complex*16 hjetmass_qqbgg_triangle_pmpm_s12_s124_0_dp
        complex*16 hjetmass_qqbgg_triangle_pmpm_s34_0_0_dp
        complex*16 hjetmass_qqbgg_triangle_pmpm_s12_s123_0_dp

        complex*16 hjetmass_qqbgg_triangle_pmpm_s34_mhsq_s12_rat_dp
        complex*16 hjetmass_qqbgg_triangle_pmpm_s124_mhsq_0_rat_dp
        complex*16 hjetmass_qqbgg_triangle_pmpm_s123_mhsq_0_rat_dp

        complex*16 hjetmass_qqbgg_bubble_pmpm_s12_dp
        complex*16 hjetmass_qqbgg_bubble_pmpm_s34_dp
        complex*16 hjetmass_qqbgg_bubble_pmpm_s123_dp
        complex*16 hjetmass_qqbgg_bubble_pmpm_s124_dp
        complex*16 hjetmass_qqbgg_bubble_pmpm_mhsq_dp

        check = (0d0,0d0)
        rat = (0d0,0d0)

        check = check +
     &  hjetmass_qqbgg_box_pmpm_0_0_s12_mhsq_s34_s124_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(2)
        check = check +
     &  hjetmass_qqbgg_box_pmpm_0_0_s12_mhsq_s34_s123_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(1)
        check = check +
     &  hjetmass_qqbgg_box_pmpm_0_s12_0_mhsq_s123_s124_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(3)
        check = check +
     &  hjetmass_qqbgg_triangle_pmpm_s34_mhsq_s12_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(6)

        tmpcast = 
     &  hjetmass_qqbgg_triangle_pmpm_s34_mhsq_s12_rat_dp
     &  (i1,i2,i3,i4,za,zb)/2d0
        rat = rat - tmpcast

        check = check +
     &  hjetmass_qqbgg_triangle_pmpm_s124_mhsq_0_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(2)

        tmpcast =
     &  hjetmass_qqbgg_triangle_pmpm_s124_mhsq_0_rat_dp
     &  (i1,i2,i3,i4,za,zb)/2d0
        rat = rat - tmpcast

        check = check +
     &  hjetmass_qqbgg_triangle_pmpm_s12_s124_0_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(4)

        check = check +
     &  hjetmass_qqbgg_triangle_pmpm_s34_0_0_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(5)
        check = check +
     &  hjetmass_qqbgg_triangle_pmpm_s123_mhsq_0_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(1)

        tmpcast =
     &  hjetmass_qqbgg_triangle_pmpm_s123_mhsq_0_rat_dp
     &  (i1,i2,i3,i4,za,zb)/2d0
        rat = rat - tmpcast

        check = check +
     &  hjetmass_qqbgg_triangle_pmpm_s12_s123_0_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(3)

        tmpcast =
     &  hjetmass_qqbgg_bubble_pmpm_mhsq_dp
     &  (i1,i2,i3,i4,za,zb,mt,p,flip)
        check = check + tmpcast
     &  * bubints(1)

        tmpcast =
     &  hjetmass_qqbgg_bubble_pmpm_s34_dp
     &  (i1,i2,i3,i4,za,zb,mt)
        check = check + tmpcast
     &  * bubints(5)

        tmpcast =
     &  hjetmass_qqbgg_bubble_pmpm_s12_dp
     &  (i1,i2,i3,i4,za,zb,mt)
        check = check + tmpcast
     &  * bubints(4)

        tmpcast =
     &  hjetmass_qqbgg_bubble_pmpm_s124_dp
     &  (i1,i2,i3,i4,za,zb,mt,p,flip)
        check = check + tmpcast
     &  * bubints(3)

        tmpcast =
     &  hjetmass_qqbgg_bubble_pmpm_s123_dp
     &  (i1,i2,i3,i4,za,zb,mt,p,flip)
        check = check + tmpcast
     &  * bubints(2)

        !hjetmass_Hqqbgg_amp_pmpm = -(0d0,6d0)*mt**2*(check + rat)
        hjetmass_Hqqbgg_amp_pmpm_dp = (check + rat)
        return

      end function

      complex*16 function
     &    hjetmass_Hqqbgg_amp_pmmp_dp(i1,i2,i3,i4,p,za,zb,flip,
     &                              boxints, triints, bubints)
        implicit none
        include 'types.f'
        include 'mxpart.f'
        complex*16 check, rat, tmpcast
        integer i1,i2,i3,i4
        include 'constants.f'
        include 'masses.f'
        complex*16 za(mxpart,mxpart), zb(mxpart,mxpart)
        real*8 p(mxpart,4)
        logical flip

        complex*16 boxints(3)
        complex*16 triints(6)
        complex*16 bubints(5)

        complex*16 hjetmass_qqbgg_box_pmmp_0_0_s12_mhsq_s34_s124_dp
        complex*16 hjetmass_qqbgg_box_pmmp_0_0_s12_mhsq_s34_s123_dp
        complex*16 hjetmass_qqbgg_box_pmmp_0_s12_0_mhsq_s123_s124_dp

        complex*16 hjetmass_qqbgg_triangle_pmmp_s34_mhsq_s12_dp
        complex*16 hjetmass_qqbgg_triangle_pmmp_s124_mhsq_0_dp
        complex*16 hjetmass_qqbgg_triangle_pmmp_s123_mhsq_0_dp
        complex*16 hjetmass_qqbgg_triangle_pmmp_s12_s124_0_dp
        complex*16 hjetmass_qqbgg_triangle_pmmp_s34_0_0_dp
        complex*16 hjetmass_qqbgg_triangle_pmmp_s12_s123_0_dp

        complex*16 hjetmass_qqbgg_triangle_pmmp_s34_mhsq_s12_rat_dp
        complex*16 hjetmass_qqbgg_triangle_pmmp_s124_mhsq_0_rat_dp
        complex*16 hjetmass_qqbgg_triangle_pmmp_s123_mhsq_0_rat_dp

        complex*16 hjetmass_qqbgg_bubble_pmmp_s12_dp
        complex*16 hjetmass_qqbgg_bubble_pmmp_s34_dp
        complex*16 hjetmass_qqbgg_bubble_pmmp_s123_dp
        complex*16 hjetmass_qqbgg_bubble_pmmp_s124_dp
        complex*16 hjetmass_qqbgg_bubble_pmmp_mhsq_dp

        check = (0d0,0d0)
        rat = (0d0,0d0)

        check = check +
     &  hjetmass_qqbgg_box_pmmp_0_0_s12_mhsq_s34_s124_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(2)
        check = check +
     &  hjetmass_qqbgg_box_pmmp_0_0_s12_mhsq_s34_s123_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(1)
        check = check +
     &  hjetmass_qqbgg_box_pmmp_0_s12_0_mhsq_s123_s124_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(3)
        check = check +
     &  hjetmass_qqbgg_triangle_pmmp_s34_mhsq_s12_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(6)

        tmpcast =
     &  hjetmass_qqbgg_triangle_pmmp_s34_mhsq_s12_rat_dp
     &  (i1,i2,i3,i4,za,zb)/2d0
        rat = rat - tmpcast

        check = check +
     &  hjetmass_qqbgg_triangle_pmmp_s124_mhsq_0_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(2)

        tmpcast =
     &  hjetmass_qqbgg_triangle_pmmp_s124_mhsq_0_rat_dp
     &  (i1,i2,i3,i4,za,zb)/2d0
        rat = rat - tmpcast

        check = check +
     &  hjetmass_qqbgg_triangle_pmmp_s12_s124_0_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(4)

        check = check +
     &  hjetmass_qqbgg_triangle_pmmp_s34_0_0_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(5)
        check = check +
     &  hjetmass_qqbgg_triangle_pmmp_s123_mhsq_0_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(1)

        tmpcast =
     &  hjetmass_qqbgg_triangle_pmmp_s123_mhsq_0_rat_dp
     &  (i1,i2,i3,i4,za,zb)/2d0
        rat = rat - tmpcast

        check = check +
     &  hjetmass_qqbgg_triangle_pmmp_s12_s123_0_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(3)

        tmpcast =
     &  hjetmass_qqbgg_bubble_pmmp_mhsq_dp
     &  (i1,i2,i3,i4,za,zb,mt,p,flip)
        check = check + tmpcast
     &  * bubints(1)

        tmpcast =
     &  hjetmass_qqbgg_bubble_pmmp_s34_dp
     &  (i1,i2,i3,i4,za,zb,mt)
        check = check + tmpcast
     &  * bubints(5)

        tmpcast =
     &  hjetmass_qqbgg_bubble_pmmp_s12_dp
     &  (i1,i2,i3,i4,za,zb,mt)
        check = check + tmpcast
     &  * bubints(4)

        tmpcast =
     &  hjetmass_qqbgg_bubble_pmmp_s124_dp
     &  (i1,i2,i3,i4,za,zb,mt,p,flip)
        check = check + tmpcast
     &  * bubints(3)

        tmpcast =
     &  hjetmass_qqbgg_bubble_pmmp_s123_dp
     &  (i1,i2,i3,i4,za,zb,mt,p,flip)
        check = check + tmpcast
     &  * bubints(2)

      !hjetmass_Hqqbgg_amp_pmmp = -(0d0,6d0)*mt**2*(check + rat)
      hjetmass_Hqqbgg_amp_pmmp_dp = (check + rat)

      end function

      complex*16 function
     &    hjetmass_Hgggg_amp_pppp_dp(i1,i2,i3,i4,za,zb,
     &                              boxints, triints)
        implicit none
        include 'types.f'
        include 'mxpart.f'
        complex*16 check, rat
        integer i1,i2,i3,i4
        include 'constants.f'
        include 'masses.f'
        complex*16 za(mxpart,mxpart), zb(mxpart,mxpart)

        complex*16 boxints(16)
        complex*16 triints(18)

        complex*16 hjetmass_triangle_pppp_s234_0_mhsq_dp
        complex*16 hjetmass_triangle_pppp_s134_0_mhsq_dp
        complex*16 hjetmass_triangle_pppp_s124_0_mhsq_dp
        complex*16 hjetmass_triangle_pppp_s123_0_mhsq_dp
        complex*16 hjetmass_triangle_pppp_s234_0_mhsq_rat_dp
        complex*16 hjetmass_triangle_pppp_s134_0_mhsq_rat_dp
        complex*16 hjetmass_triangle_pppp_s124_0_mhsq_rat_dp
        complex*16 hjetmass_triangle_pppp_s123_0_mhsq_rat_dp

        complex*16 hjetmass_box_pppp_0_0_s12_mhsq_s34_s124_dp
        complex*16 hjetmass_box_pppp_0_0_s12_mhsq_s34_s123_dp
        complex*16 hjetmass_box_pppp_0_s12_0_mhsq_s123_s124_dp
        complex*16 hjetmass_box_pppp_0_0_s14_mhsq_s23_s134_dp
        complex*16 hjetmass_box_pppp_0_0_s14_mhsq_s23_s124_dp
        complex*16 hjetmass_box_pppp_0_s14_0_mhsq_s124_s134_dp
        complex*16 hjetmass_box_pppp_0_0_s23_mhsq_s14_s234_dp
        complex*16 hjetmass_box_pppp_0_0_s23_mhsq_s14_s123_dp
        complex*16 hjetmass_box_pppp_0_s23_0_mhsq_s123_s234_dp
        complex*16 hjetmass_box_pppp_0_0_s34_mhsq_s12_s234_dp
        complex*16 hjetmass_box_pppp_0_0_s34_mhsq_s12_s134_dp
        complex*16 hjetmass_box_pppp_0_s34_0_mhsq_s134_s234_dp
        complex*16 hjetmass_box_pppp_0_0_0_s124_s14_s12_dp
        complex*16 hjetmass_box_pppp_0_0_0_s123_s12_s23_dp
        complex*16 hjetmass_box_pppp_0_0_0_s134_s34_s14_dp
        complex*16 hjetmass_box_pppp_0_0_0_s234_s23_s34_dp

        check = (0d0,0d0)
        rat = (0d0,0d0)

        ! pppp
        check = check +
     &  hjetmass_box_pppp_0_0_s12_mhsq_s34_s124_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(6)
        check = check +
     &  hjetmass_box_pppp_0_0_s12_mhsq_s34_s123_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(5)
        check = check +
     &  hjetmass_box_pppp_0_s12_0_mhsq_s123_s124_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(13)
        check = check +
     &  hjetmass_box_pppp_0_0_s14_mhsq_s23_s134_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(8)
        check = check +
     &  hjetmass_box_pppp_0_0_s14_mhsq_s23_s124_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(7)
        check = check +
     &  hjetmass_box_pppp_0_s14_0_mhsq_s124_s134_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(14)
        check = check +
     &  hjetmass_box_pppp_0_0_s23_mhsq_s14_s234_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(10)
        check = check +
     &  hjetmass_box_pppp_0_0_s23_mhsq_s14_s123_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(9)
        check = check +
     &  hjetmass_box_pppp_0_s23_0_mhsq_s123_s234_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(15)
        check = check +
     &  hjetmass_box_pppp_0_0_s34_mhsq_s12_s234_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(12)
        check = check +
     &  hjetmass_box_pppp_0_0_s34_mhsq_s12_s134_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(11)
        check = check +
     &  hjetmass_box_pppp_0_s34_0_mhsq_s134_s234_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(16)
        check = check +
     &  hjetmass_box_pppp_0_0_0_s124_s14_s12_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(2)
        check = check +
     &  hjetmass_box_pppp_0_0_0_s123_s12_s23_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(1)
        check = check +
     &  hjetmass_box_pppp_0_0_0_s134_s34_s14_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(3)
        check = check +
     &  hjetmass_box_pppp_0_0_0_s234_s23_s34_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(4)
        check = check +
     &  hjetmass_triangle_pppp_s234_0_mhsq_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(11)
        rat = rat + 
     &  hjetmass_triangle_pppp_s234_0_mhsq_rat_dp
     &  (i1,i2,i3,i4,za,zb)/2d0
        check = check +
     &  hjetmass_triangle_pppp_s134_0_mhsq_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(6)
        rat = rat + 
     &  hjetmass_triangle_pppp_s134_0_mhsq_rat_dp
     &  (i1,i2,i3,i4,za,zb)/2d0
        check = check +
     &  hjetmass_triangle_pppp_s124_0_mhsq_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(3)
        rat = rat + 
     &  hjetmass_triangle_pppp_s124_0_mhsq_rat_dp
     &  (i1,i2,i3,i4,za,zb)/2d0
        check = check +
     &  hjetmass_triangle_pppp_s123_0_mhsq_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(2)
        rat = rat + 
     &  hjetmass_triangle_pppp_s123_0_mhsq_rat_dp
     &  (i1,i2,i3,i4,za,zb)/2d0

        hjetmass_Hgggg_amp_pppp_dp = (check + rat)
        return

      end function

      complex*16 function
     &    hjetmass_Hgggg_amp_pppm_dp(i1,i2,i3,i4,p,za,zb,flip,
     &                              boxints, triints, bubints)
        implicit none
        include 'types.f'
        include 'mxpart.f'
        complex*16 check, rat, tmpcast
        integer i1,i2,i3,i4
        include 'constants.f'
        include 'masses.f'
        complex*16 za(mxpart,mxpart), zb(mxpart,mxpart)
        logical flip
        real*8 p(mxpart,4)

        complex*16 boxints(16)
        complex*16 triints(18)
        complex*16 bubints(9)

        complex*16 hjetmass_box_pppm_0_0_s34_mhsq_s12_s234_dp
        complex*16 hjetmass_box_pppm_0_0_s34_mhsq_s12_s134_dp
        complex*16 hjetmass_box_pppm_0_0_s23_mhsq_s14_s123_dp
        complex*16 hjetmass_box_pppm_0_0_s23_mhsq_s14_s234_dp
        complex*16 hjetmass_box_pppm_0_0_s14_mhsq_s23_s124_dp
        complex*16 hjetmass_box_pppm_0_0_s14_mhsq_s23_s134_dp
        complex*16 hjetmass_box_pppm_0_0_s12_mhsq_s34_s123_dp
        complex*16 hjetmass_box_pppm_0_0_s12_mhsq_s34_s124_dp
        complex*16 hjetmass_box_pppm_0_s23_0_mhsq_s123_s234_dp
        complex*16 hjetmass_box_pppm_0_s34_0_mhsq_s134_s234_dp
        complex*16 hjetmass_box_pppm_0_s14_0_mhsq_s124_s134_dp
        complex*16 hjetmass_box_pppm_0_s12_0_mhsq_s123_s124_dp
        complex*16 hjetmass_box_pppm_0_0_0_s123_s12_s23_dp
        complex*16 hjetmass_box_pppm_0_0_0_s134_s34_s14_dp
        complex*16 hjetmass_box_pppm_0_0_0_s234_s23_s34_dp
        complex*16 hjetmass_box_pppm_0_0_0_s124_s14_s12_dp

        complex*16 hjetmass_triangle_pppm_s123_0_mhsq_dp
        complex*16 hjetmass_triangle_pppm_s124_0_mhsq_dp
        complex*16 hjetmass_triangle_pppm_s134_0_mhsq_dp
        complex*16 hjetmass_triangle_pppm_s234_0_mhsq_dp
        complex*16 hjetmass_triangle_pppm_s23_mhsq_s14_dp
        complex*16 hjetmass_triangle_pppm_s34_mhsq_s12_dp
        complex*16 hjetmass_triangle_pppm_s34_0_0_dp
        complex*16 hjetmass_triangle_pppm_s14_0_0_dp
        complex*16 hjetmass_triangle_pppm_s14_s134_0_dp
        complex*16 hjetmass_triangle_pppm_s14_s124_0_dp
        complex*16 hjetmass_triangle_pppm_s34_s234_0_dp
        complex*16 hjetmass_triangle_pppm_s34_s134_0_dp

        complex*16 hjetmass_triangle_pppm_s123_0_mhsq_rat_dp
        complex*16 hjetmass_triangle_pppm_s124_0_mhsq_rat_dp
        complex*16 hjetmass_triangle_pppm_s134_0_mhsq_rat_dp
        complex*16 hjetmass_triangle_pppm_s234_0_mhsq_rat_dp
        complex*16 hjetmass_triangle_pppm_s23_mhsq_s14_rat_dp
        complex*16 hjetmass_triangle_pppm_s34_mhsq_s12_rat_dp

        complex*16 hjetmass_bubble_pppm_mhsq_dp
        complex*16 hjetmass_bubble_pppm_s124_dp
        complex*16 hjetmass_bubble_pppm_s134_dp
        complex*16 hjetmass_bubble_pppm_s14_dp
        complex*16 hjetmass_bubble_pppm_s234_dp
        complex*16 hjetmass_bubble_pppm_s34_dp

        check = (0d0,0d0)
        rat = (0d0,0d0)

        ! pppm 
        check = check +
     &  hjetmass_box_pppm_0_0_s12_mhsq_s34_s124_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(6)
        check = check +
     &  hjetmass_box_pppm_0_0_s12_mhsq_s34_s123_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(5)
        check = check +
     &  hjetmass_box_pppm_0_s12_0_mhsq_s123_s124_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(13)
        check = check +
     &  hjetmass_box_pppm_0_0_s14_mhsq_s23_s134_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(8)
        check = check +
     &  hjetmass_box_pppm_0_0_s14_mhsq_s23_s124_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(7)
        check = check +
     &  hjetmass_box_pppm_0_s14_0_mhsq_s124_s134_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(14)
        check = check +
     &  hjetmass_box_pppm_0_0_s23_mhsq_s14_s234_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(10)
        check = check +
     &  hjetmass_box_pppm_0_0_s23_mhsq_s14_s123_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(9)
        check = check +
     &  hjetmass_box_pppm_0_s23_0_mhsq_s123_s234_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(15)
        check = check +
     &  hjetmass_box_pppm_0_0_s34_mhsq_s12_s234_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(12)
        check = check +
     &  hjetmass_box_pppm_0_0_s34_mhsq_s12_s134_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(11)
        check = check +
     &  hjetmass_box_pppm_0_s34_0_mhsq_s134_s234_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(16)
        check = check +
     &  hjetmass_box_pppm_0_0_0_s124_s14_s12_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(2)
        check = check +
     &  hjetmass_box_pppm_0_0_0_s123_s12_s23_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(1)
        check = check +
     &  hjetmass_box_pppm_0_0_0_s134_s34_s14_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(3)
        check = check +
     &  hjetmass_box_pppm_0_0_0_s234_s23_s34_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(4)
        check = check +
     &  hjetmass_triangle_pppm_s234_0_mhsq_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(11)
        rat = rat + 
     &  hjetmass_triangle_pppm_s234_0_mhsq_rat_dp
     &  (i1,i2,i3,i4,za,zb)/2d0
        check = check +
     &  hjetmass_triangle_pppm_s134_0_mhsq_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(6)
        rat = rat + 
     &  hjetmass_triangle_pppm_s134_0_mhsq_rat_dp
     &  (i1,i2,i3,i4,za,zb)/2d0
        check = check +
     &  hjetmass_triangle_pppm_s124_0_mhsq_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(3)
        rat = rat + 
     &  hjetmass_triangle_pppm_s124_0_mhsq_rat_dp
     &  (i1,i2,i3,i4,za,zb)/2d0
        check = check +
     &  hjetmass_triangle_pppm_s123_0_mhsq_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(2)
        rat = rat + 
     &  hjetmass_triangle_pppm_s123_0_mhsq_rat_dp
     &  (i1,i2,i3,i4,za,zb)/2d0
        check = check +
     &  hjetmass_triangle_pppm_s34_mhsq_s12_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(16)
        rat = rat + 
     &  hjetmass_triangle_pppm_s34_mhsq_s12_rat_dp
     &  (i1,i2,i3,i4,za,zb)/2d0
        check = check +
     &  hjetmass_triangle_pppm_s23_mhsq_s14_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(12)
        rat = rat + 
     &  hjetmass_triangle_pppm_s23_mhsq_s14_rat_dp
     &  (i1,i2,i3,i4,za,zb)/2d0
        check = check +
     &  hjetmass_triangle_pppm_s34_0_0_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(15)
        check = check +
     &  hjetmass_triangle_pppm_s14_s134_0_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(9)
        check = check +
     &  hjetmass_triangle_pppm_s14_s124_0_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(8)
        check = check +
     &  hjetmass_triangle_pppm_s14_0_0_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(7)
        check = check +
     &  hjetmass_triangle_pppm_s34_s234_0_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(18)
        check = check +
     &  hjetmass_triangle_pppm_s34_s134_0_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(17)

        tmpcast =
     &  hjetmass_bubble_pppm_s234_dp
     &  (i1,i2,i3,i4,za,zb,mt,p,flip)
        check = check + tmpcast
     &  * bubints(7)

        tmpcast =
     &  hjetmass_bubble_pppm_mhsq_dp
     &  (i1,i2,i3,i4,za,zb,mt,p,flip)
        check = check + tmpcast
     &  * bubints(1)

        tmpcast =
     &  hjetmass_bubble_pppm_s134_dp
     &  (i1,i2,i3,i4,za,zb,mt,p,flip)
        check = check + tmpcast
     &  * bubints(5)

        tmpcast =
     &  hjetmass_bubble_pppm_s124_dp
     &  (i1,i2,i3,i4,za,zb,mt,p,flip)
        check = check + tmpcast
     &  * bubints(3)

        tmpcast =
     &  hjetmass_bubble_pppm_s34_dp
     &  (i1,i2,i3,i4,za,zb,mt)
        check = check + tmpcast
     &  * bubints(9)

        tmpcast = 
     &  hjetmass_bubble_pppm_s14_dp
     &  (i1,i2,i3,i4,za,zb,mt)
        check = check + tmpcast
     &  * bubints(6)

        !hjetmass_Hgggg_amp_pppm = -(0d0,6d0)*mt**2*(check + rat)
        hjetmass_Hgggg_amp_pppm_dp = (check + rat)
        return

      end function

      complex*16 function
     &    hjetmass_Hgggg_amp_ppmm_dp(i1,i2,i3,i4,p,za,zb,flip,
     &                              boxints, triints, bubints)
        implicit none
        include 'types.f'
        include 'mxpart.f'
        complex*16 check, rat, tmpcast
        integer i1,i2,i3,i4
        include 'constants.f'
        include 'masses.f'
        complex*16 za(mxpart,mxpart), zb(mxpart,mxpart)
        logical flip
        real*8 p(mxpart,4)

        complex*16 boxints(16)
        complex*16 triints(18)
        complex*16 bubints(9)

        complex*16 hjetmass_triangle_ppmm_s234_0_mhsq_dp
        complex*16 hjetmass_triangle_ppmm_s134_0_mhsq_dp
        complex*16 hjetmass_triangle_ppmm_s124_0_mhsq_dp
        complex*16 hjetmass_triangle_ppmm_s123_0_mhsq_dp
        complex*16 hjetmass_triangle_ppmm_s23_mhsq_s14_dp
        complex*16 hjetmass_triangle_ppmm_s23_0_0_dp
        complex*16 hjetmass_triangle_ppmm_s14_0_0_dp
        complex*16 hjetmass_triangle_ppmm_s14_s134_0_dp
        complex*16 hjetmass_triangle_ppmm_s14_s124_0_dp
        complex*16 hjetmass_triangle_ppmm_s23_s234_0_dp
        complex*16 hjetmass_triangle_ppmm_s23_s123_0_dp

        complex*16 hjetmass_triangle_ppmm_s234_0_mhsq_rat_dp
        complex*16 hjetmass_triangle_ppmm_s134_0_mhsq_rat_dp
        complex*16 hjetmass_triangle_ppmm_s124_0_mhsq_rat_dp
        complex*16 hjetmass_triangle_ppmm_s123_0_mhsq_rat_dp
        complex*16 hjetmass_triangle_ppmm_s23_mhsq_s14_rat_dp

        complex*16 hjetmass_box_ppmm_0_0_s12_mhsq_s34_s124_dp
        complex*16 hjetmass_box_ppmm_0_0_s12_mhsq_s34_s123_dp
        complex*16 hjetmass_box_ppmm_0_s12_0_mhsq_s123_s124_dp
        complex*16 hjetmass_box_ppmm_0_0_s14_mhsq_s23_s134_dp
        complex*16 hjetmass_box_ppmm_0_0_s14_mhsq_s23_s124_dp
        complex*16 hjetmass_box_ppmm_0_s14_0_mhsq_s124_s134_dp
        complex*16 hjetmass_box_ppmm_0_0_s23_mhsq_s14_s234_dp
        complex*16 hjetmass_box_ppmm_0_0_s23_mhsq_s14_s123_dp
        complex*16 hjetmass_box_ppmm_0_s23_0_mhsq_s123_s234_dp
        complex*16 hjetmass_box_ppmm_0_0_s34_mhsq_s12_s234_dp
        complex*16 hjetmass_box_ppmm_0_0_s34_mhsq_s12_s134_dp
        complex*16 hjetmass_box_ppmm_0_s34_0_mhsq_s134_s234_dp
        complex*16 hjetmass_box_ppmm_0_0_0_s124_s14_s12_dp
        complex*16 hjetmass_box_ppmm_0_0_0_s123_s12_s23_dp
        complex*16 hjetmass_box_ppmm_0_0_0_s134_s34_s14_dp
        complex*16 hjetmass_box_ppmm_0_0_0_s234_s23_s34_dp

        complex*16 hjetmass_bubble_ppmm_mhsq_dp
        complex*16 hjetmass_bubble_ppmm_s123_dp
        complex*16 hjetmass_bubble_ppmm_s124_dp
        complex*16 hjetmass_bubble_ppmm_s134_dp
        complex*16 hjetmass_bubble_ppmm_s14_dp
        complex*16 hjetmass_bubble_ppmm_s234_dp
        complex*16 hjetmass_bubble_ppmm_s23_dp

        check = (0d0,0d0)
        rat = (0d0,0d0)
        tmpcast = (0d0,0d0)

        ! ppmm
        check = check +
     &  hjetmass_box_ppmm_0_0_s12_mhsq_s34_s124_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(6)
        check = check +
     &  hjetmass_box_ppmm_0_0_s12_mhsq_s34_s123_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(5)
        check = check +
     &  hjetmass_box_ppmm_0_s12_0_mhsq_s123_s124_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(13)
        check = check +
     &  hjetmass_box_ppmm_0_0_s14_mhsq_s23_s134_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(8)
        check = check +
     &  hjetmass_box_ppmm_0_0_s14_mhsq_s23_s124_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(7)
        check = check +
     &  hjetmass_box_ppmm_0_s14_0_mhsq_s124_s134_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(14)
        check = check +
     &  hjetmass_box_ppmm_0_0_s23_mhsq_s14_s234_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(10)
        check = check +
     &  hjetmass_box_ppmm_0_0_s23_mhsq_s14_s123_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(9)
        check = check +
     &  hjetmass_box_ppmm_0_s23_0_mhsq_s123_s234_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(15)
        check = check +
     &  hjetmass_box_ppmm_0_0_s34_mhsq_s12_s234_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(12)
        check = check +
     &  hjetmass_box_ppmm_0_0_s34_mhsq_s12_s134_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(11)
        check = check +
     &  hjetmass_box_ppmm_0_s34_0_mhsq_s134_s234_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(16)
        check = check +
     &  hjetmass_box_ppmm_0_0_0_s124_s14_s12_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(2)
        check = check +
     &  hjetmass_box_ppmm_0_0_0_s123_s12_s23_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(1)
        check = check +
     &  hjetmass_box_ppmm_0_0_0_s134_s34_s14_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(3)
        check = check +
     &  hjetmass_box_ppmm_0_0_0_s234_s23_s34_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(4)
        check = check +
     &  hjetmass_triangle_ppmm_s234_0_mhsq_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(11)
        rat = rat + 
     &  hjetmass_triangle_ppmm_s234_0_mhsq_rat_dp
     &  (i1,i2,i3,i4,za,zb)/2d0
        check = check +
     &  hjetmass_triangle_ppmm_s134_0_mhsq_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(6)
        rat = rat + 
     &  hjetmass_triangle_ppmm_s134_0_mhsq_rat_dp
     &  (i1,i2,i3,i4,za,zb)/2d0
        check = check +
     &  hjetmass_triangle_ppmm_s124_0_mhsq_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(3)
        rat = rat + 
     &  hjetmass_triangle_ppmm_s124_0_mhsq_rat_dp
     &  (i1,i2,i3,i4,za,zb)/2d0
        check = check +
     &  hjetmass_triangle_ppmm_s123_0_mhsq_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(2)
        rat = rat + 
     &  hjetmass_triangle_ppmm_s123_0_mhsq_rat_dp
     &  (i1,i2,i3,i4,za,zb)/2d0 
        check = check +
     &  hjetmass_triangle_ppmm_s23_mhsq_s14_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(12)
        rat = rat + 
     &  hjetmass_triangle_ppmm_s23_mhsq_s14_rat_dp
     &  (i1,i2,i3,i4,za,zb)/2d0
        check = check +
     &  hjetmass_triangle_ppmm_s14_s134_0_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(9)
        check = check +
     &  hjetmass_triangle_ppmm_s23_0_0_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(10)
        check = check +
     &  hjetmass_triangle_ppmm_s14_s124_0_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(8)
        check = check +
     &  hjetmass_triangle_ppmm_s23_s234_0_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(14)
        check = check +
     &  hjetmass_triangle_ppmm_s14_0_0_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(7)
        check = check +
     &  hjetmass_triangle_ppmm_s23_s123_0_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(13)

        tmpcast = 
     &  hjetmass_bubble_ppmm_s234_dp
     &  (i1,i2,i3,i4,za,zb,mt,p,flip)
        check = check + tmpcast
     &  * bubints(7)

        tmpcast = 
     &  hjetmass_bubble_ppmm_mhsq_dp
     &  (i1,i2,i3,i4,za,zb,mt,p,flip)
        check = check + tmpcast
     &  * bubints(1)

        tmpcast =
     &  hjetmass_bubble_ppmm_s134_dp
     &  (i1,i2,i3,i4,za,zb,mt,p,flip)
        check = check + tmpcast
     &  * bubints(5)

        tmpcast =
     &  hjetmass_bubble_ppmm_s124_dp
     &  (i1,i2,i3,i4,za,zb,mt,p,flip)
        check = check + tmpcast
     &  * bubints(3)

        tmpcast =
     &  hjetmass_bubble_ppmm_s123_dp
     &  (i1,i2,i3,i4,za,zb,mt,p,flip)
        check = check + tmpcast
     &  * bubints(2)

        tmpcast =
     &  hjetmass_bubble_ppmm_s23_dp
     &  (i1,i2,i3,i4,za,zb,mt)
        check = check + tmpcast
     &  * bubints(8)

        tmpcast =
     &  hjetmass_bubble_ppmm_s14_dp
     &  (i1,i2,i3,i4,za,zb,mt)
        check = check + tmpcast
     &  * bubints(6)

        hjetmass_Hgggg_amp_ppmm_dp = (check + rat)
        return

      end function

      complex*16 function
     &    hjetmass_Hgggg_amp_pmpm_dp(i1,i2,i3,i4,p,za,zb,flip,
     &                              boxints, triints, bubints)
        implicit none
        include 'types.f'
        include 'mxpart.f'
        complex*16 check, rat, tmpcast
        integer i1,i2,i3,i4
        include 'constants.f'
        include 'masses.f'
        complex*16 za(mxpart,mxpart), zb(mxpart,mxpart)
        logical flip
        real*8 p(mxpart,4)

        complex*16 boxints(16)
        complex*16 triints(18)
        complex*16 bubints(9)

        complex*16 hjetmass_triangle_pmpm_s234_0_mhsq_dp
        complex*16 hjetmass_triangle_pmpm_s134_0_mhsq_dp
        complex*16 hjetmass_triangle_pmpm_s124_0_mhsq_dp
        complex*16 hjetmass_triangle_pmpm_s123_0_mhsq_dp
        complex*16 hjetmass_triangle_pmpm_s23_mhsq_s14_dp
        complex*16 hjetmass_triangle_pmpm_s34_mhsq_s12_dp
        complex*16 hjetmass_triangle_pmpm_s34_0_0_dp
        complex*16 hjetmass_triangle_pmpm_s23_0_0_dp
        complex*16 hjetmass_triangle_pmpm_s14_0_0_dp
        complex*16 hjetmass_triangle_pmpm_s12_0_0_dp
        complex*16 hjetmass_triangle_pmpm_s12_s124_0_dp
        complex*16 hjetmass_triangle_pmpm_s12_s123_0_dp
        complex*16 hjetmass_triangle_pmpm_s14_s134_0_dp
        complex*16 hjetmass_triangle_pmpm_s14_s124_0_dp
        complex*16 hjetmass_triangle_pmpm_s23_s234_0_dp
        complex*16 hjetmass_triangle_pmpm_s23_s123_0_dp
        complex*16 hjetmass_triangle_pmpm_s34_s234_0_dp
        complex*16 hjetmass_triangle_pmpm_s34_s134_0_dp

        complex*16 hjetmass_triangle_pmpm_s234_0_mhsq_rat_dp
        complex*16 hjetmass_triangle_pmpm_s134_0_mhsq_rat_dp
        complex*16 hjetmass_triangle_pmpm_s124_0_mhsq_rat_dp
        complex*16 hjetmass_triangle_pmpm_s123_0_mhsq_rat_dp
        complex*16 hjetmass_triangle_pmpm_s23_mhsq_s14_rat_dp
        complex*16 hjetmass_triangle_pmpm_s34_mhsq_s12_rat_dp

        complex*16 hjetmass_box_pmpm_0_0_s12_mhsq_s34_s124_dp
        complex*16 hjetmass_box_pmpm_0_0_s12_mhsq_s34_s123_dp
        complex*16 hjetmass_box_pmpm_0_s12_0_mhsq_s123_s124_dp
        complex*16 hjetmass_box_pmpm_0_0_s14_mhsq_s23_s134_dp
        complex*16 hjetmass_box_pmpm_0_0_s14_mhsq_s23_s124_dp
        complex*16 hjetmass_box_pmpm_0_s14_0_mhsq_s124_s134_dp
        complex*16 hjetmass_box_pmpm_0_0_s23_mhsq_s14_s234_dp
        complex*16 hjetmass_box_pmpm_0_0_s23_mhsq_s14_s123_dp
        complex*16 hjetmass_box_pmpm_0_s23_0_mhsq_s123_s234_dp
        complex*16 hjetmass_box_pmpm_0_0_s34_mhsq_s12_s234_dp
        complex*16 hjetmass_box_pmpm_0_0_s34_mhsq_s12_s134_dp
        complex*16 hjetmass_box_pmpm_0_s34_0_mhsq_s134_s234_dp
        complex*16 hjetmass_box_pmpm_0_0_0_s124_s14_s12_dp
        complex*16 hjetmass_box_pmpm_0_0_0_s123_s12_s23_dp
        complex*16 hjetmass_box_pmpm_0_0_0_s134_s34_s14_dp
        complex*16 hjetmass_box_pmpm_0_0_0_s234_s23_s34_dp

        complex*16 hjetmass_bubble_pmpm_mhsq_dp
        complex*16 hjetmass_bubble_pmpm_s123_dp
        complex*16 hjetmass_bubble_pmpm_s124_dp
        complex*16 hjetmass_bubble_pmpm_s12_dp
        complex*16 hjetmass_bubble_pmpm_s134_dp
        complex*16 hjetmass_bubble_pmpm_s14_dp
        complex*16 hjetmass_bubble_pmpm_s234_dp
        complex*16 hjetmass_bubble_pmpm_s23_dp
        complex*16 hjetmass_bubble_pmpm_s34_dp

        check = (0d0,0d0)
        rat = (0d0,0d0)
        tmpcast = (0d0,0d0)

        ! pmpm
        check = check +
     &  hjetmass_box_pmpm_0_0_s12_mhsq_s34_s124_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(6)
        check = check +
     &  hjetmass_box_pmpm_0_0_s12_mhsq_s34_s123_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(5)
        check = check +
     &  hjetmass_box_pmpm_0_s12_0_mhsq_s123_s124_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(13)
        check = check +
     &  hjetmass_box_pmpm_0_0_s14_mhsq_s23_s134_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(8)
        check = check +
     &  hjetmass_box_pmpm_0_0_s14_mhsq_s23_s124_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(7)
        check = check +
     &  hjetmass_box_pmpm_0_s14_0_mhsq_s124_s134_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(14)
        check = check +
     &  hjetmass_box_pmpm_0_0_s23_mhsq_s14_s234_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(10)
        check = check +
     &  hjetmass_box_pmpm_0_0_s23_mhsq_s14_s123_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(9)
        check = check +
     &  hjetmass_box_pmpm_0_s23_0_mhsq_s123_s234_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(15)
        check = check +
     &  hjetmass_box_pmpm_0_0_s34_mhsq_s12_s234_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(12)
        check = check +
     &  hjetmass_box_pmpm_0_0_s34_mhsq_s12_s134_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(11)
        check = check +
     &  hjetmass_box_pmpm_0_s34_0_mhsq_s134_s234_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(16)
        check = check +
     &  hjetmass_box_pmpm_0_0_0_s124_s14_s12_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(2)
        check = check +
     &  hjetmass_box_pmpm_0_0_0_s123_s12_s23_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(1)
        check = check +
     &  hjetmass_box_pmpm_0_0_0_s134_s34_s14_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(3)
        check = check +
     &  hjetmass_box_pmpm_0_0_0_s234_s23_s34_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * boxints(4)
        check = check +
     &  hjetmass_triangle_pmpm_s234_0_mhsq_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(11)
        rat = rat + 
     &  hjetmass_triangle_pmpm_s234_0_mhsq_rat_dp
     &  (i1,i2,i3,i4,za,zb)/2d0
        check = check +
     &  hjetmass_triangle_pmpm_s134_0_mhsq_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(6)
        rat = rat + 
     &  hjetmass_triangle_pmpm_s134_0_mhsq_rat_dp
     &  (i1,i2,i3,i4,za,zb)/2d0
        check = check +
     &  hjetmass_triangle_pmpm_s124_0_mhsq_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(3)
        rat = rat + 
     &  hjetmass_triangle_pmpm_s124_0_mhsq_rat_dp
     &  (i1,i2,i3,i4,za,zb)/2d0
        check = check +
     &  hjetmass_triangle_pmpm_s123_0_mhsq_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(2)
        rat = rat + 
     &  hjetmass_triangle_pmpm_s123_0_mhsq_rat_dp
     &  (i1,i2,i3,i4,za,zb)/2d0
        check = check +
     &  hjetmass_triangle_pmpm_s34_mhsq_s12_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(16)
        rat = rat + 
     &  hjetmass_triangle_pmpm_s34_mhsq_s12_rat_dp
     &  (i1,i2,i3,i4,za,zb)/2d0
        check = check +
     &  hjetmass_triangle_pmpm_s23_mhsq_s14_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(12)
        rat = rat + 
     &  hjetmass_triangle_pmpm_s23_mhsq_s14_rat_dp
     &  (i1,i2,i3,i4,za,zb)/2d0
        check = check +
     &  hjetmass_triangle_pmpm_s12_s124_0_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(5)
        check = check +
     &  hjetmass_triangle_pmpm_s34_0_0_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(15)
        check = check +
     &  hjetmass_triangle_pmpm_s12_s123_0_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(4)
        check = check +
     &  hjetmass_triangle_pmpm_s14_s134_0_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(9)
        check = check +
     &  hjetmass_triangle_pmpm_s23_0_0_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(10)
        check = check +
     &  hjetmass_triangle_pmpm_s14_s124_0_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(8)
        check = check +
     &  hjetmass_triangle_pmpm_s23_s234_0_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(14)
        check = check +
     &  hjetmass_triangle_pmpm_s14_0_0_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(7)
        check = check +
     &  hjetmass_triangle_pmpm_s23_s123_0_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(13)
        check = check +
     &  hjetmass_triangle_pmpm_s34_s234_0_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(18)
        check = check +
     &  hjetmass_triangle_pmpm_s12_0_0_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(1)
        check = check +
     &  hjetmass_triangle_pmpm_s34_s134_0_dp
     &  (i1,i2,i3,i4,za,zb,mt)
     &  * triints(17)

        tmpcast = 
     &  hjetmass_bubble_pmpm_s234_dp
     &  (i1,i2,i3,i4,za,zb,mt,p,flip)
        check = check + tmpcast
     &  * bubints(7)

        tmpcast =
     &  hjetmass_bubble_pmpm_mhsq_dp
     &  (i1,i2,i3,i4,za,zb,mt,p,flip)
        check = check + tmpcast
     &  * bubints(1)

        tmpcast =
     &  hjetmass_bubble_pmpm_s134_dp
     &  (i1,i2,i3,i4,za,zb,mt,p,flip)
        check = check + tmpcast
     &  * bubints(5)

        tmpcast =
     &  hjetmass_bubble_pmpm_s124_dp
     &  (i1,i2,i3,i4,za,zb,mt,p,flip)
        check = check + tmpcast
     &  * bubints(3)

        tmpcast =
     &  hjetmass_bubble_pmpm_s123_dp
     &  (i1,i2,i3,i4,za,zb,mt,p,flip)
        check = check + tmpcast
     &  * bubints(2)

        tmpcast =
     &  hjetmass_bubble_pmpm_s34_dp
     &  (i1,i2,i3,i4,za,zb,mt)
        check = check + tmpcast
     &  * bubints(9)

        tmpcast =
     &  hjetmass_bubble_pmpm_s12_dp
     &  (i1,i2,i3,i4,za,zb,mt)
        check = check + tmpcast
     &  * bubints(4)

        tmpcast =
     &  hjetmass_bubble_pmpm_s23_dp
     &  (i1,i2,i3,i4,za,zb,mt)
        check = check + tmpcast
     &  * bubints(8)

        tmpcast =
     &  hjetmass_bubble_pmpm_s14_dp
     &  (i1,i2,i3,i4,za,zb,mt)
        check = check + tmpcast
     &  * bubints(6)

        hjetmass_Hgggg_amp_pmpm_dp = (check + rat)
        return

      end function

      ! follows structure of ggHg/HQAggnew
      subroutine hjetmass_Hqqbgg_dp(i1,i2,i3,i4,p,
     & ampsq,ampsq_ab,ampsq_ba,ampsq_sym)
        implicit none
        include 'types.f'
        include 'mxpart.f'
        include 'constants.f'
        include 'masses.f'

        integer, intent(in) :: i1,i2,i3,i4
        real(dp), intent(in) :: p(mxpart,4)
        real(dp), intent(out) :: ampsq, ampsq_ba, ampsq_ab, ampsq_sym

        integer j1,j2,j3
        complex*16 amp(2,2,2,2)

        ampsq = 0d0
        ampsq_ba = 0d0
        ampsq_ab = 0d0
        ampsq_sym = 0d0

        call hjetmass_Hqqbgg_amp_dp(p,amp,i1,i2,i3,i4)
        do j1=1,2; do j2=1,2; do j3=1,2
            ampsq_ab = ampsq_ab + cf*xn**2/2d0*abs(amp(1,j1,j2,j3))**2
            ampsq_ba = ampsq_ba + cf*xn**2/2d0*abs(amp(2,j1,j2,j3))**2
            ampsq_sym = ampsq_sym -
     &                cf/2d0*abs(amp(1,j1,j2,j3)+amp(2,j1,j2,j3))**2
        enddo; enddo; enddo

        ampsq = ampsq_ab + ampsq_ba + ampsq_sym

        ampsq =  ampsq*(6d0*mt**2)**2
        ampsq_ab =  ampsq_ab*(6d0*mt**2)**2
        ampsq_ba =  ampsq_ba*(6d0*mt**2)**2
        ampsq_sym =  ampsq_sym*(6d0*mt**2)**2

      end subroutine
      
      subroutine hjetmass_HqqbQQb_amp_dp(p,amp,i1,i2,i3,i4)
        implicit none
        include 'types.f'
        include 'mxpart.f'
        include 'constants.f'
        include 'masses.f'
        include 'zprods_decl.f'
        include 'HqqbQQblabels.f'
        real*8 s(mxpart,mxpart)
        include 'scale.f'
        include 'qcdcouple.f'
        include 'ewcouple.f'
        real*8 p(mxpart,4)
        integer i1,i2,i3,i4

        complex*16 amp(2,2)

        complex*16 triints(1)
        complex*16 bubints(3)
        za(:,:) = (0d0,0d0)
        zb(:,:) = (0d0,0d0)
        s(:,:) = (0d0,0d0)
        amp(:,:) = (0d0,0d0)

          call spinoru_dp_s(4,p,za,zb,s)

          call hjetmass_qqbqqb_triints_init_dp
     &                (i1,i2,i3,i4,triints,0,s)
          call hjetmass_qqbqqb_bubints_init_dp
     &                (i1,i2,i3,i4,bubints,0,s)

          amp(2,2) = hjetmass_HqqbQQb_amp_pmpm_dp
     &      (i1,i2,i3,i4,p,za,zb,triints,bubints)
          amp(1,1) = hjetmass_HqqbQQb_amp_pmpm_dp
     &      (i1,i2,i3,i4,p,zb,za,triints,bubints)

          amp(2,1) = hjetmass_HqqbQQb_amp_pmmp_dp
     &      (i1,i2,i3,i4,p,za,zb,triints,bubints)
          amp(1,2) = hjetmass_HqqbQQb_amp_pmmp_dp
     &      (i1,i2,i3,i4,p,zb,za,triints,bubints)
         
      end subroutine hjetmass_HqqbQQb_amp_dp

      ! follows structure of ggHg/Ampsq_AQaq.f
      ! dum1,dum2,dum3 introduced for consistent interface, unused
      subroutine hjetmass_HqqbQQb_nonid_dp(i1,i2,i3,i4,p,ampsq,
     &                  dum1,dum2,dum3)
        implicit none
        include 'types.f'
        include 'mxpart.f'
        include 'constants.f'
        include 'masses.f'

        integer, intent(in) :: i1,i2,i3,i4
        real(dp), intent(in) :: p(mxpart,4)
        real(dp), intent(out) :: ampsq,dum1,dum2,dum3

        complex(dp) :: amp(2,2)
        integer h1,h2

        call hjetmass_HqqbQQb_amp_dp(p,amp,i1,i2,i3,i4)

        ampsq = 0d0
        do h1=1,2; do h2=1,2
          ampsq = ampsq + cdabs(amp(h1,h2))**2
        enddo; enddo

        ampsq = ampsq*V/4d0*(6d0*mt**2)**2

        return
      end subroutine hjetmass_HqqbQQb_nonid_dp

      subroutine hjetmass_HqqbQQb_ident_dp(i1,i2,i3,i4,p,
     &  ampsq, ampsq_a, ampsq_b, ampsq_i)
        implicit none
        include 'types.f'
        include 'mxpart.f'
        include 'constants.f'
        include 'masses.f'

        integer, intent(in) :: i1,i2,i3,i4
        real(dp), intent(in) :: p(mxpart,4)
        real(dp), intent(out) :: ampsq, ampsq_a, ampsq_b, ampsq_i

        complex(dp) :: amp_a(2,2), amp_b(2,2)
        integer h1,h2

        call hjetmass_HqqbQQb_amp_dp(p,amp_a,i1,i2,i3,i4)
        call hjetmass_HqqbQQb_amp_dp(p,amp_b,i1,i4,i3,i2)

        ampsq = 0d0
        ampsq_a = 0d0
        ampsq_b = 0d0
        ampsq_i = 0d0

        do h1=1,2; do h2=1,2
          ampsq_a = ampsq_a + cdabs(amp_a(h1,h2))**2
          ampsq_b = ampsq_b + cdabs(amp_b(h1,h2))**2
          if (h1 .eq. h2) then
            ampsq_i = ampsq_i +
     &           2d0/xn*dble(amp_a(h1,h2)*dconjg(amp_b(h1,h2)))
          end if
        enddo; enddo

        ampsq_a = ampsq_a*V/4d0
        ampsq_b = ampsq_b*V/4d0
        ampsq_i = ampsq_i*V/4d0
        ampsq = ampsq_a + ampsq_b + ampsq_i

        ampsq_a = ampsq_a*(6d0*mt**2)**2
        ampsq_b = ampsq_b*(6d0*mt**2)**2
        ampsq_i = ampsq_i*(6d0*mt**2)**2
        ampsq   = ampsq*(6d0*mt**2)**2

        return
      end subroutine hjetmass_HqqbQQb_ident_dp


      subroutine hjetmass_Hqqbgg_amp_dp(p,amp,i1,i2,i3,i4)
        implicit none
        include 'types.f'
        include 'mxpart.f'
        include 'constants.f'
        include 'masses.f'
        complex*16 za(mxpart,mxpart), zb(mxpart,mxpart)
        real*8 s(mxpart,mxpart)
        include 'scale.f'
        include 'qcdcouple.f'
        include 'ewcouple.f'
        real*8 p(mxpart,4)
        integer i1,i2,i3,i4

        ! first index for two permutations
        ! second index for first quark helicity
        ! third, fourth indices: gluon helicities
        ! index 1 for minus, 2 for plus helicity
        complex*16 amp(2,2,2,2)

        complex*16 boxints(3)
        complex*16 triints(6)
        complex*16 bubints(5)

        integer cp(2,4)

        call spinoru_dp_s(4,p,za,zb,s)

        cp(1,:) = (/ i1,i2,i4,i3 /)
        cp(2,:) = (/ i1,i2,i3,i4 /)

        ! first ordering (ab)
            
          call hjetmass_qqbgg_boxints_init_dp
     &                (cp(1,1),cp(1,2),cp(1,3),cp(1,4),boxints,0,s)
          call hjetmass_qqbgg_triints_init_dp
     &                (cp(1,1),cp(1,2),cp(1,3),cp(1,4),triints,0,s)
          call hjetmass_qqbgg_bubints_init_dp
     &                (cp(1,1),cp(1,2),cp(1,3),cp(1,4),bubints,0,s)

          amp(1,2,2,2) = -hjetmass_Hqqbgg_amp_pmpp_dp
     &      (cp(1,1),cp(1,2),cp(1,3),cp(1,4),
     &       p,za,zb,.false.,boxints,triints,bubints)
          amp(1,1,1,1) = hjetmass_Hqqbgg_amp_pmpp_dp
     &      (cp(1,1),cp(1,2),cp(1,3),cp(1,4),
     &       p,zb,za,.true.,boxints,triints,bubints)

          amp(1,2,2,1) = -hjetmass_Hqqbgg_amp_pmpm_dp
     &      (cp(1,1),cp(1,2),cp(1,3),cp(1,4),
     &       p,za,zb,.false.,boxints,triints,bubints)
          amp(1,1,1,2) = hjetmass_Hqqbgg_amp_pmpm_dp
     &      (cp(1,1),cp(1,2),cp(1,3),cp(1,4),
     &       p,zb,za,.true.,boxints,triints,bubints)

          amp(1,2,1,2) = -hjetmass_Hqqbgg_amp_pmmp_dp
     &      (cp(1,1),cp(1,2),cp(1,3),cp(1,4),
     &       p,za,zb,.false.,boxints,triints,bubints)
          amp(1,1,2,1) = hjetmass_Hqqbgg_amp_pmmp_dp
     &      (cp(1,1),cp(1,2),cp(1,3),cp(1,4),
     &       p,zb,za,.true.,boxints,triints,bubints)


          ! maybe calculate pmmm directly to avoid additional qcdloop calls
          call hjetmass_qqbgg_boxints_init_dp
     &                (cp(1,2),cp(1,1),cp(1,4),cp(1,3),boxints,0,s)
          call hjetmass_qqbgg_triints_init_dp
     &                (cp(1,2),cp(1,1),cp(1,4),cp(1,3),triints,0,s)
          call hjetmass_qqbgg_bubints_init_dp
     &                (cp(1,2),cp(1,1),cp(1,4),cp(1,3),bubints,0,s)

          amp(1,2,1,1) = hjetmass_Hqqbgg_amp_pmpp_dp
     &      (cp(1,2),cp(1,1),cp(1,4),cp(1,3),
     &       p,zb,za,.true.,boxints,triints,bubints)
          amp(1,1,2,2) = -hjetmass_Hqqbgg_amp_pmpp_dp
     &      (cp(1,2),cp(1,1),cp(1,4),cp(1,3),
     &       p,za,zb,.false.,boxints,triints,bubints)

        ! second color ordering (ba)

          call hjetmass_qqbgg_boxints_init_dp
     &                (cp(2,1),cp(2,2),cp(2,3),cp(2,4),boxints,0,s)
          call hjetmass_qqbgg_triints_init_dp
     &                (cp(2,1),cp(2,2),cp(2,3),cp(2,4),triints,0,s)
          call hjetmass_qqbgg_bubints_init_dp
     &                (cp(2,1),cp(2,2),cp(2,3),cp(2,4),bubints,0,s)

          amp(2,2,2,2) = -hjetmass_Hqqbgg_amp_pmpp_dp
     &      (cp(2,1),cp(2,2),cp(2,3),cp(2,4),
     &       p,za,zb,.false.,boxints,triints,bubints)
          amp(2,1,1,1) = hjetmass_Hqqbgg_amp_pmpp_dp
     &      (cp(2,1),cp(2,2),cp(2,3),cp(2,4),
     &       p,zb,za,.true.,boxints,triints,bubints)

          amp(2,2,2,1) = -hjetmass_Hqqbgg_amp_pmmp_dp
     &      (cp(2,1),cp(2,2),cp(2,3),cp(2,4),
     &       p,za,zb,.false.,boxints,triints,bubints)
          amp(2,1,1,2) = hjetmass_Hqqbgg_amp_pmmp_dp
     &      (cp(2,1),cp(2,2),cp(2,3),cp(2,4),
     &       p,zb,za,.true.,boxints,triints,bubints)

          amp(2,2,1,2) = -hjetmass_Hqqbgg_amp_pmpm_dp
     &      (cp(2,1),cp(2,2),cp(2,3),cp(2,4),
     &       p,za,zb,.false.,boxints,triints,bubints)
          amp(2,1,2,1) = hjetmass_Hqqbgg_amp_pmpm_dp
     &      (cp(2,1),cp(2,2),cp(2,3),cp(2,4),
     &       p,zb,za,.true.,boxints,triints,bubints)


          ! maybe calculate pmmm directly to avoid additional qcdloop calls
          call hjetmass_qqbgg_boxints_init_dp
     &                (cp(2,2),cp(2,1),cp(2,4),cp(2,3),boxints,0,s)
          call hjetmass_qqbgg_triints_init_dp
     &                (cp(2,2),cp(2,1),cp(2,4),cp(2,3),triints,0,s)
          call hjetmass_qqbgg_bubints_init_dp
     &                (cp(2,2),cp(2,1),cp(2,4),cp(2,3),bubints,0,s)

          amp(2,2,1,1) = hjetmass_Hqqbgg_amp_pmpp_dp
     &      (cp(2,2),cp(2,1),cp(2,4),cp(2,3),
     &       p,zb,za,.true.,boxints,triints,bubints)
          amp(2,1,2,2) = -hjetmass_Hqqbgg_amp_pmpp_dp
     &      (cp(2,2),cp(2,1),cp(2,4),cp(2,3),
     &       p,za,zb,.false.,boxints,triints,bubints)

      end subroutine

      subroutine hjetmass_Hgggg_amp_dp(p,amp)
        implicit none
        include 'types.f'
        include 'mxpart.f'
        include 'constants.f'
        include 'masses.f'
        complex*16 za(mxpart,mxpart), zb(mxpart,mxpart)
        real*8 s(mxpart,mxpart)
        include 'scale.f'
        include 'qcdcouple.f'
        include 'ewcouple.f'
        real*8 p(mxpart,4)
        integer j

        ! first index for three permutations
        ! other indices: 1 for minus, 2 for plus helicity
        complex*16 amp(3,2,2,2,2)

        complex*16 boxints(16)
        complex*16 triints(18)
        complex*16 bubints(9)
        integer cp(3,4)

        call spinoru_dp_s(4,p,za,zb,s)

        ! summation over three non-cyclic permutations
        cp(1,:) = (/ 1,2,3,4 /)
        cp(2,:) = (/ 1,2,4,3 /)
        cp(3,:) = (/ 1,4,2,3 /)

        do j=1,3
          call hjetmass_boxints_init_dp
     &              (cp(j,1),cp(j,2),cp(j,3),cp(j,4),boxints,0,s)
          call hjetmass_triints_init_dp
     &              (cp(j,1),cp(j,2),cp(j,3),cp(j,4),triints,0,s)
          call hjetmass_bubints_init_dp
     &              (cp(j,1),cp(j,2),cp(j,3),cp(j,4),bubints,0,s)

          amp(j,2,2,2,2) = hjetmass_Hgggg_amp_pppp_dp
     &      (cp(j,1),cp(j,2),cp(j,3),cp(j,4),
     &       za,zb,boxints,triints)
          amp(j,1,1,1,1) = hjetmass_Hgggg_amp_pppp_dp
     &      (cp(j,1),cp(j,2),cp(j,3),cp(j,4),
     &       zb,za,boxints,triints)

          amp(j,2,2,2,1) = hjetmass_Hgggg_amp_pppm_dp
     &      (cp(j,1),cp(j,2),cp(j,3),cp(j,4),
     &       p,za,zb,.false.,boxints,triints,bubints)
          amp(j,1,1,1,2) = hjetmass_Hgggg_amp_pppm_dp
     &      (cp(j,1),cp(j,2),cp(j,3),cp(j,4),
     &       p,zb,za,.true.,boxints,triints,bubints)

          amp(j,2,2,1,1) = hjetmass_Hgggg_amp_ppmm_dp
     &      (cp(j,1),cp(j,2),cp(j,3),cp(j,4),
     &       p,za,zb,.false.,boxints,triints,bubints)
          amp(j,1,1,2,2) = hjetmass_Hgggg_amp_ppmm_dp
     &      (cp(j,1),cp(j,2),cp(j,3),cp(j,4),
     &       p,zb,za,.true.,boxints,triints,bubints)

          amp(j,2,1,2,1) = hjetmass_Hgggg_amp_pmpm_dp
     &      (cp(j,1),cp(j,2),cp(j,3),cp(j,4),
     &       p,za,zb,.false.,boxints,triints,bubints)
          amp(j,1,2,1,2) = hjetmass_Hgggg_amp_pmpm_dp
     &      (cp(j,1),cp(j,2),cp(j,3),cp(j,4),
     &       p,zb,za,.true.,boxints,triints,bubints)

          call hjetmass_boxints_init_dp
     &              (cp(j,4),cp(j,1),cp(j,2),cp(j,3),boxints,0,s)
          call hjetmass_triints_init_dp
     &              (cp(j,4),cp(j,1),cp(j,2),cp(j,3),triints,0,s)
          call hjetmass_bubints_init_dp
     &              (cp(j,4),cp(j,1),cp(j,2),cp(j,3),bubints,0,s)
          amp(j,2,1,1,2) = hjetmass_Hgggg_amp_ppmm_dp
     &      (cp(j,4),cp(j,1),cp(j,2),cp(j,3),
     &       p,za,zb,.false.,boxints,triints,bubints)
          amp(j,1,2,2,1) = hjetmass_Hgggg_amp_ppmm_dp
     &      (cp(j,4),cp(j,1),cp(j,2),cp(j,3),
     &       p,zb,za,.true.,boxints,triints,bubints)

          amp(j,2,2,1,2) = hjetmass_Hgggg_amp_pppm_dp
     &      (cp(j,4),cp(j,1),cp(j,2),cp(j,3),
     &       p,za,zb,.false.,boxints,triints,bubints)
          amp(j,1,1,2,1) = hjetmass_Hgggg_amp_pppm_dp
     &      (cp(j,4),cp(j,1),cp(j,2),cp(j,3),
     &       p,zb,za,.true.,boxints,triints,bubints)

          call hjetmass_boxints_init_dp
     &              (cp(j,3),cp(j,4),cp(j,1),cp(j,2),boxints,0,s)
          call hjetmass_triints_init_dp
     &              (cp(j,3),cp(j,4),cp(j,1),cp(j,2),triints,0,s)
          call hjetmass_bubints_init_dp
     &              (cp(j,3),cp(j,4),cp(j,1),cp(j,2),bubints,0,s)
          amp(j,2,1,2,2) = hjetmass_Hgggg_amp_pppm_dp
     &      (cp(j,3),cp(j,4),cp(j,1),cp(j,2),
     &       p,za,zb,.false.,boxints,triints,bubints)
          amp(j,1,2,1,1) = hjetmass_Hgggg_amp_pppm_dp
     &      (cp(j,3),cp(j,4),cp(j,1),cp(j,2),
     &       p,zb,za,.true.,boxints,triints,bubints)

          call hjetmass_boxints_init_dp
     &              (cp(j,2),cp(j,3),cp(j,4),cp(j,1),boxints,0,s)
          call hjetmass_triints_init_dp
     &              (cp(j,2),cp(j,3),cp(j,4),cp(j,1),triints,0,s)
          call hjetmass_bubints_init_dp
     &              (cp(j,2),cp(j,3),cp(j,4),cp(j,1),bubints,0,s)
          amp(j,1,2,2,2) = hjetmass_Hgggg_amp_pppm_dp
     &      (cp(j,2),cp(j,3),cp(j,4),cp(j,1),
     &       p,za,zb,.false.,boxints,triints,bubints)
          amp(j,2,1,1,1) = hjetmass_Hgggg_amp_pppm_dp
     &      (cp(j,2),cp(j,3),cp(j,4),cp(j,1),
     &       p,zb,za,.true.,boxints,triints,bubints)

        enddo

        return

        end subroutine

      ! i1,i2,i3,i4  dummy arguments for a consistent interface
      subroutine hjetmass_Hgggg_dp(i1,i2,i3,i4,p,
     & hjetmass_Hgggg_full, hjetmass_Hgggg_1234, hjetmass_Hgggg_1243,
     & hjetmass_Hgggg_1423)
        implicit none
        include 'types.f'
        include 'mxpart.f'
        include 'constants.f'
        include 'masses.f'
        complex*16 amp(3,2,2,2,2)
        integer h1,h2,h3,h4

        integer, intent(in) :: i1,i2,i3,i4
        real(dp), intent(in) :: p(mxpart,4)
        real(dp), intent(out) :: hjetmass_Hgggg_full,
     &   hjetmass_Hgggg_1234, hjetmass_Hgggg_1243, hjetmass_Hgggg_1423

      call hjetmass_Hgggg_amp_dp(p,amp)

      hjetmass_Hgggg_1234 = 0d0
      hjetmass_Hgggg_1243 = 0d0
      hjetmass_Hgggg_1423 = 0d0
      hjetmass_Hgggg_full = 0d0

      do h1=1,2; do h2=1,2; do h3=1,2; do h4=1,2
        hjetmass_Hgggg_1234 = hjetmass_Hgggg_1234 +
     &                  xn**2*V/2d0*abs(amp(1,h1,h2,h3,h4))**2
        hjetmass_Hgggg_1243 = hjetmass_Hgggg_1243 +
     &                  xn**2*V/2d0*abs(amp(2,h1,h2,h3,h4))**2
        hjetmass_Hgggg_1423 = hjetmass_Hgggg_1423 +
     &                  xn**2*V/2d0*abs(amp(3,h1,h2,h3,h4))**2
      enddo; enddo; enddo; enddo


       hjetmass_Hgggg_full = hjetmass_Hgggg_1234 + hjetmass_Hgggg_1243
     &      + hjetmass_Hgggg_1423

       hjetmass_Hgggg_full =  hjetmass_Hgggg_full*(6*mt**2)**2
       hjetmass_Hgggg_1234 =  hjetmass_Hgggg_1234*(6*mt**2)**2
       hjetmass_Hgggg_1243 =  hjetmass_Hgggg_1243*(6*mt**2)**2
       hjetmass_Hgggg_1423 =  hjetmass_Hgggg_1423*(6*mt**2)**2
        return
      end subroutine

      end module

