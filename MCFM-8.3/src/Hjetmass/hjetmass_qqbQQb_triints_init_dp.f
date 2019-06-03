      subroutine hjetmass_qqbQQb_triints_init_dp
     & (i1,i2,i3,i4,triints,ord,s)
        use mod_qcdloop_c
          implicit none
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'scale.f'
          include 'masses.f'
          integer i1,i2,i3,i4
          integer ord
          complex*16 triints(1)
          real*8 t, s(mxpart,mxpart)
          real*8 mt2, mhsq

        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        mt2 = mt**2
        mhsq = s(i1,i2)+s(i1,i3)+s(i1,i4)+s(i2,i3)+s(i2,i4)+s(i3,i4)

        triints(:) = (0d0,0d0)

c   "( s34,mhsq, s12;mtsq,mtsq,mtsq)"
        triints(1) = 
     &  qlI3(s(i3,i4), mhsq, s(i1,i2), mt2, mt2, mt2, musq, ord)

      end subroutine
