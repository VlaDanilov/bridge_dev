
      complex*32 function hjetmass_box_pppp_0_0_0_s123_s12_s23
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          double precision mt
          complex*32 hjetmass_box_onemass_pppp

      hjetmass_box_pppp_0_0_0_s123_s12_s23 =
     &  hjetmass_box_onemass_pppp(i1,i2,i3,i4,za,zb,mt)
      return

      end function
