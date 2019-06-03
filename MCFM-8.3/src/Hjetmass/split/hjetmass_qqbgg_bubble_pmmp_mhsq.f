
      complex*32 function hjetmass_qqbgg_bubble_pmmp_mhsq
     &     (i1,i2,i3,i4,za,zb,mt,p,flip)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          double precision mt
          real*16 p(mxpart,4)
          complex*32 hjetmass_qqbgg_bubble_pmmp_s124
          complex*32 hjetmass_qqbgg_bubble_pmmp_s123
          complex*32 hjetmass_qqbgg_bubble_pmmp_s34
          complex*32 hjetmass_qqbgg_bubble_pmmp_s12
          logical flip



      hjetmass_qqbgg_bubble_pmmp_mhsq = 
     & -hjetmass_qqbgg_bubble_pmmp_s124(i1,i2,i3,i4,za,zb,mt,p,flip)
     & -hjetmass_qqbgg_bubble_pmmp_s123(i1,i2,i3,i4,za,zb,mt,p,flip)
     & -hjetmass_qqbgg_bubble_pmmp_s34(i1,i2,i3,i4,za,zb,mt)
     & -hjetmass_qqbgg_bubble_pmmp_s12(i1,i2,i3,i4,za,zb,mt)
      return

      end function
