
      complex*32 function hjetmass_bubble_pmpm_mhsq
     &     (i1,i2,i3,i4,za,zb,mt,p,flip)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          double precision mt
          real*16 p(mxpart,4)
          complex*32 hjetmass_bubble_pmpm_s234
          complex*32 hjetmass_bubble_pmpm_s134
          complex*32 hjetmass_bubble_pmpm_s124
          complex*32 hjetmass_bubble_pmpm_s123
          complex*32 hjetmass_bubble_pmpm_s34
          complex*32 hjetmass_bubble_pmpm_s12
          complex*32 hjetmass_bubble_pmpm_s23
          complex*32 hjetmass_bubble_pmpm_s14
          logical flip



      hjetmass_bubble_pmpm_mhsq = 
     & -hjetmass_bubble_pmpm_s234(i1,i2,i3,i4,za,zb,mt,p,flip)
     & -hjetmass_bubble_pmpm_s134(i1,i2,i3,i4,za,zb,mt,p,flip)
     & -hjetmass_bubble_pmpm_s124(i1,i2,i3,i4,za,zb,mt,p,flip)
     & -hjetmass_bubble_pmpm_s123(i1,i2,i3,i4,za,zb,mt,p,flip)
     & -hjetmass_bubble_pmpm_s34(i1,i2,i3,i4,za,zb,mt)
     & -hjetmass_bubble_pmpm_s12(i1,i2,i3,i4,za,zb,mt)
     & -hjetmass_bubble_pmpm_s23(i1,i2,i3,i4,za,zb,mt)
     & -hjetmass_bubble_pmpm_s14(i1,i2,i3,i4,za,zb,mt)
      return

      end function
