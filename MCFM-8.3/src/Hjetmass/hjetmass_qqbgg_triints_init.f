      subroutine hjetmass_qqbgg_triints_init(i1,i2,i3,i4,triints,ord,s)
        use mod_qcdloop_c
          implicit none
          include 'types.f'
          include 'mxpart.f'

          include 'constants.f'
          real*16 s(mxpart,mxpart)
          include 'scale.f'
          include 'masses.f'
          integer i1,i2,i3,i4
          integer ord
          complex*32 triints(6)
          real*16 t
          real*16 mt2, mhsq, qmusq

        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        mt2 = mt**2
        mhsq = s(i1,i2)+s(i1,i3)+s(i1,i4)+s(i2,i3)+s(i2,i4)+s(i3,i4)
        qmusq = musq
      

        triints(:) = (0q0,0q0)

c   "(s123,mhsq,   0;mtsq,mtsq,mtsq)",
        triints(1) = qlI3q(t(i1,i2,i3), mhsq, 0q0, mt2, mt2, mt2, qmusq, ord)

c   "(s124,mhsq,   0;mtsq,mtsq,mtsq)",
        triints(2) = qlI3q(t(i1,i2,i4), mhsq, 0q0, mt2, mt2, mt2, qmusq, ord)

c   "( s12,s123,   0;mtsq,mtsq,mtsq)",
        triints(3) = qlI3q(s(i1,i2), t(i1,i2,i3), 0q0, mt2, mt2, mt2, qmusq, ord)

c   "( s12,s124,   0;mtsq,mtsq,mtsq)",
        triints(4) = qlI3q(s(i1,i2), t(i1,i2,i4), 0q0, mt2, mt2, mt2, qmusq, ord)

c   "( s34,   0,   0;mtsq,mtsq,mtsq)",
        triints(5) = qlI3q(s(i3,i4), 0q0, 0q0, mt2, mt2, mt2, qmusq, ord)

c   "( s34,mhsq, s12;mtsq,mtsq,mtsq)"
        triints(6) = qlI3q(s(i3,i4), mhsq, s(i1,i2), mt2, mt2, mt2, qmusq, ord)

      end subroutine
