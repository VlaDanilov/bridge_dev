
      double complex function hjetmass_bubble_pmpm_mhsq_dp 
     &     (i1,i2,i3,i4,za,zb,mt,p,flip)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double precision mt
          double precision p(mxpart,4)
          double complex hjetmass_bubble_pmpm_s234_dp
          double complex hjetmass_bubble_pmpm_s134_dp
          double complex hjetmass_bubble_pmpm_s124_dp
          double complex hjetmass_bubble_pmpm_s123_dp
          double complex hjetmass_bubble_pmpm_s34_dp
          double complex hjetmass_bubble_pmpm_s12_dp
          double complex hjetmass_bubble_pmpm_s23_dp
          double complex hjetmass_bubble_pmpm_s14_dp
          logical flip



      hjetmass_bubble_pmpm_mhsq_dp = 
     & -hjetmass_bubble_pmpm_s234_dp(i1,i2,i3,i4,za,zb,mt,p,flip)
     & -hjetmass_bubble_pmpm_s134_dp(i1,i2,i3,i4,za,zb,mt,p,flip)
     & -hjetmass_bubble_pmpm_s124_dp(i1,i2,i3,i4,za,zb,mt,p,flip)
     & -hjetmass_bubble_pmpm_s123_dp(i1,i2,i3,i4,za,zb,mt,p,flip)
     & -hjetmass_bubble_pmpm_s34_dp(i1,i2,i3,i4,za,zb,mt)
     & -hjetmass_bubble_pmpm_s12_dp(i1,i2,i3,i4,za,zb,mt)
     & -hjetmass_bubble_pmpm_s23_dp(i1,i2,i3,i4,za,zb,mt)
     & -hjetmass_bubble_pmpm_s14_dp(i1,i2,i3,i4,za,zb,mt)
      return

      end function
