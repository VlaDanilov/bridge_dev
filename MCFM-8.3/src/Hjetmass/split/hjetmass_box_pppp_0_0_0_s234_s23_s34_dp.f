
      double complex function hjetmass_box_pppp_0_0_0_s234_s23_s34_dp 
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double precision mt
          double complex hjetmass_box_onemass_pppp_dp

      hjetmass_box_pppp_0_0_0_s234_s23_s34_dp =
     &  hjetmass_box_onemass_pppp_dp(i2,i3,i4,i1,za,zb,mt)
      return

      end function
