
      double complex function hjetmass_box_pppp_0_s12_0_mhsq_s123_s124_dp 
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double precision mt
          double complex hjetmass_box_twomass_easy_pppp_dp

      hjetmass_box_pppp_0_s12_0_mhsq_s123_s124_dp =
     &  hjetmass_box_twomass_easy_pppp_dp(i3,i2,i1,i4,za,zb,mt)
      return

      end function
