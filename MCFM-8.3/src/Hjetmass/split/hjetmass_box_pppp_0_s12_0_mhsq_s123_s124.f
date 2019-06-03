
      complex*32 function hjetmass_box_pppp_0_s12_0_mhsq_s123_s124
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          double precision mt
          complex*32 hjetmass_box_twomass_easy_pppp

      hjetmass_box_pppp_0_s12_0_mhsq_s123_s124 =
     &  hjetmass_box_twomass_easy_pppp(i3,i2,i1,i4,za,zb,mt)
      return

      end function
