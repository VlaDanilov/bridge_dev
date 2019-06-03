
      double complex function hjetmass_qqbgg_bubble_pmmp_mhsq_dp 
     &     (i1,i2,i3,i4,za,zb,mt,p,flip)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double precision mt
          double precision p(mxpart,4)
          double complex hjetmass_qqbgg_bubble_pmmp_s124_dp
          double complex hjetmass_qqbgg_bubble_pmmp_s123_dp
          double complex hjetmass_qqbgg_bubble_pmmp_s34_dp
          double complex hjetmass_qqbgg_bubble_pmmp_s12_dp
          logical flip



      hjetmass_qqbgg_bubble_pmmp_mhsq_dp = 
     & -hjetmass_qqbgg_bubble_pmmp_s124_dp(i1,i2,i3,i4,za,zb,mt,p,flip)
     & -hjetmass_qqbgg_bubble_pmmp_s123_dp(i1,i2,i3,i4,za,zb,mt,p,flip)
     & -hjetmass_qqbgg_bubble_pmmp_s34_dp(i1,i2,i3,i4,za,zb,mt)
     & -hjetmass_qqbgg_bubble_pmmp_s12_dp(i1,i2,i3,i4,za,zb,mt)
      return

      end function
