
      double complex function hjetmass_box_pppp_0_0_s23_mhsq_s14_s234_dp 
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double precision mt
          double complex hjetmass_box_twomass_hard_pppp_dp

      hjetmass_box_pppp_0_0_s23_mhsq_s14_s234_dp =
     &  hjetmass_box_twomass_hard_pppp_dp(i1,i4,i3,i2,za,zb,mt)
      return

      end function
