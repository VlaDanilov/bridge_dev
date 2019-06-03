      subroutine hjetmass_triints_init(i1,i2,i3,i4,triints,ord,s)
        use mod_qcdloop_c
          implicit none
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'scale.f'
          include 'masses.f'
          integer i1,i2,i3,i4
          integer ord
          complex*32 triints(18)
          real*16 t, s(mxpart,mxpart)
          real*16 mt2, mhsq, qmusq

        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        mt2 = mt**2
        mhsq = s(i1,i2)+s(i1,i3)+s(i1,i4)+s(i2,i3)+s(i2,i4)+s(i3,i4)
        qmusq = musq

        triints(:) = (0q0,0q0)

!         ( s12,   0,   0;mtsq,mtsq,mtsq)
        triints(1) = qlI3q(s(i1,i2), 0q0, 0q0, mt2, mt2, mt2, qmusq, ord)

!         (s123,   0,mhsq;mtsq,mtsq,mtsq)
        triints(2) = qlI3q(t(i1,i2,i3), 0q0, mhsq, mt2, mt2, mt2, qmusq, ord)

!         (s124,   0,mhsq;mtsq,mtsq,mtsq)
        triints(3) = qlI3q(t(i1,i2,i4), 0q0, mhsq, mt2, mt2, mt2, qmusq, ord)

!         ( s12,s123,   0;mtsq,mtsq,mtsq)
        triints(4) = qlI3q(s(i1,i2), t(i1,i2,i3), 0q0, mt2, mt2, mt2, qmusq, ord)

!         ( s12,s124,   0;mtsq,mtsq,mtsq)
        triints(5) = qlI3q(s(i1,i2), t(i1,i2,i4), 0q0, mt2, mt2, mt2, qmusq, ord)

!         (s134,   0,mhsq;mtsq,mtsq,mtsq)
        triints(6) = qlI3q(t(i1,i3,i4), 0q0, mhsq, mt2, mt2, mt2, qmusq, ord)

!         ( s14,   0,   0;mtsq,mtsq,mtsq)
        triints(7) = qlI3q(s(i1,i4), 0q0, 0q0, mt2, mt2, mt2, qmusq, ord)

!         ( s14,s124,   0;mtsq,mtsq,mtsq)
        triints(8) = qlI3q(s(i1,i4), t(i1,i2,i4), 0q0, mt2, mt2, mt2, qmusq, ord)

!         ( s14,s134,   0;mtsq,mtsq,mtsq)
        triints(9) = qlI3q(s(i1,i4), t(i1,i3,i4), 0q0, mt2, mt2, mt2, qmusq, ord)

!         ( s23,   0,   0;mtsq,mtsq,mtsq)
        triints(10) = qlI3q(s(i2,i3), 0q0, 0q0, mt2, mt2, mt2, qmusq, ord)

!         (s234,   0,mhsq;mtsq,mtsq,mtsq)
        triints(11) = qlI3q(t(i2,i3,i4), 0q0, mhsq, mt2, mt2, mt2, qmusq, ord)

!         ( s23,mhsq, s14;mtsq,mtsq,mtsq)
        triints(12) = qlI3q(s(i2,i3), mhsq, s(i1,i4), mt2, mt2, mt2, qmusq, ord)

!         ( s23,s123,   0;mtsq,mtsq,mtsq)
        triints(13) = qlI3q(s(i2,i3), t(i1,i2,i3), 0q0, mt2, mt2, mt2, qmusq, ord)

!         ( s23,s234,   0;mtsq,mtsq,mtsq)
        triints(14) = qlI3q(s(i2,i3), t(i2,i3,i4), 0q0, mt2, mt2, mt2, qmusq, ord)

!         ( s34,   0,   0;mtsq,mtsq,mtsq)
        triints(15) = qlI3q(s(i3,i4), 0q0, 0q0, mt2, mt2, mt2, qmusq, ord)

!         ( s34,mhsq, s12;mtsq,mtsq,mtsq)
        triints(16) = qlI3q(s(i3,i4), mhsq, s(i1,i2), mt2, mt2, mt2, qmusq, ord)

!         ( s34,s134,   0;mtsq,mtsq,mtsq)
        triints(17) = qlI3q(s(i3,i4), t(i1,i3,i4), 0q0, mt2, mt2, mt2, qmusq, ord)

!         ( s34,s234,   0;mtsq,mtsq,mtsq)
        triints(18) = qlI3q(s(i3,i4), t(i2,i3,i4), 0q0, mt2, mt2, mt2, qmusq, ord)

      end subroutine
