      subroutine hjetmass_qqbQQb_triints_init
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
          complex*32 triints(1)
          real*16 t, s(mxpart,mxpart)
          real*16 mt2, mhsq, qmusq

        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        mt2 = mt**2
        mhsq = s(i1,i2)+s(i1,i3)+s(i1,i4)+s(i2,i3)+s(i2,i4)+s(i3,i4)
        qmusq = musq

        triints(:) = (0d0,0d0)

c   "( s34,mhsq, s12;mtsq,mtsq,mtsq)"
       triints(1) = qlI3q(s(i3,i4), mhsq, s(i1,i2), mt2, mt2, mt2, qmusq, ord)

      end subroutine
